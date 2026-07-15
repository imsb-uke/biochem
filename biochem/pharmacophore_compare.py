import os
import shutil
import subprocess

import pandas as pd


## ========= helper functions =============
## NOTE: these wrap CDPL's Python bindings (`CDPKit`). They were written
## without a live CDPKit install to test against, since it is not available
## on this machine. `psdcreate` (CLI) and the `.pml` pharmacophore I/O are
## well-documented and low-risk. The `.psd` screening-database read/write
## calls in `_build_psd_database`'s fallback branch and in
## `_align_library_to_query` (Pharm.PSDScreeningDBCreator,
## Pharm.PSDScreeningDBAccessor, and their methods) are the least certain —
## validate these first against the actual installed CDPKit version.

def _prepare_molecule(mol):
    """Run the standard CDPL property-perception sequence needed before
    pharmacophore generation."""
    from CDPL import Chem

    Chem.perceiveComponents(mol, False)
    Chem.perceiveSSSR(mol, False)
    Chem.setRingFlags(mol, False)
    Chem.calcImplicitHydrogenCounts(mol, False)
    Chem.perceiveHybridizationStates(mol, False)
    Chem.setAromaticityFlags(mol, False)


def _build_query_pharmacophore(query_sdf: str, query_pml_path: str):
    """Generate a CDPL pharmacophore from the query ligand's 3D structure and save it as .pml."""
    from CDPL import Chem, Pharm

    reader = Chem.MoleculeReader(query_sdf)
    mol = Chem.BasicMolecule()
    if not reader.read(mol):
        raise ValueError(f"Could not read a molecule from {query_sdf}")

    _prepare_molecule(mol)

    pharm_gen = Pharm.DefaultPharmacophoreGenerator()
    pharmacophore = Pharm.BasicPharmacophore()
    pharm_gen.generate(mol, pharmacophore)

    writer = Pharm.FilePMLFeatureContainerWriter(query_pml_path)
    writer.write(pharmacophore)
    writer.close()

    return pharmacophore


def _build_psd_database(candidate_sdf: str, psd_path: str, psdcreate_bin: str | None):
    """Build a CDPL pharmacophore screening database (.psd) from a multi-molecule SDF,
    preferring the psdcreate CLI and falling back to the Python CDPL API."""
    if psdcreate_bin and shutil.which(psdcreate_bin):
        cmd = [psdcreate_bin, "-i", candidate_sdf, "-o", psd_path]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"psdcreate failed:\n{result.stderr}")
        return

    # Fallback: build the .psd directly through the Python CDPL API
    from CDPL import Chem, Pharm

    reader = Chem.MoleculeReader(candidate_sdf)
    db_creator = Pharm.PSDScreeningDBCreator(psd_path, Pharm.ScreeningDBCreator.CREATE)

    mol = Chem.BasicMolecule()
    while reader.read(mol):
        _prepare_molecule(mol)
        db_creator.process(mol)
        mol = Chem.BasicMolecule()

    db_creator.close()


def _align_library_to_query(query_pharmacophore, psd_path: str, default_score: float):
    """Align every pharmacophore stored in the .psd database against the query
    pharmacophore and max-pool the fit score per candidate molecule name."""
    from CDPL import Chem, Pharm

    db_accessor = Pharm.PSDScreeningDBAccessor(psd_path)
    alignment = Pharm.PharmacophoreAlignment(True)
    alignment.addFeatures(query_pharmacophore, True)
    score_func = Pharm.PharmacophoreFitScore()

    best_scores = {}

    num_entries = db_accessor.getNumMolecules()
    for i in range(num_entries):
        mol = Chem.BasicMolecule()
        db_accessor.getMolecule(i, mol)
        candidate_name = Chem.getName(mol) or f"candidate_{i}"

        candidate_pharmacophore = Pharm.BasicPharmacophore()
        db_accessor.getFeatures(i, candidate_pharmacophore)

        alignment.clearEntities(False)
        alignment.addFeatures(candidate_pharmacophore, False)

        best_score = default_score
        while alignment.nextAlignment():
            transform = alignment.getTransform()
            score = score_func(query_pharmacophore, candidate_pharmacophore, transform)
            best_score = max(best_score, score)

        best_scores[candidate_name] = max(best_scores.get(candidate_name, default_score), best_score)

    return best_scores


## ========= Main function ==============

def screen_library_pharmacophore(query_sdf: str,
                                  candidate_ligand_df_file: str | None,
                                  file_dir: str,
                                  project_name: str = 'my_screening',
                                  psdcreate_bin: str | None = 'psdcreate',
                                  psd_path: str | None = None,
                                  default_score: float = 0.0,
                                 ) -> dict:
    """Screen a library of prepared candidate ligands (from `prepare_ligand`) against
    one known-active query molecule using CDPKit/CDPL pharmacophore alignment.
    Returns candidates ranked by continuous CDPL fit score to the query pharmacophore."""

    os.makedirs(file_dir, exist_ok=True)

    query_pml_path = os.path.join(file_dir, f"{project_name}_query.pml")
    scores_csv_path = os.path.join(file_dir, f"{project_name}_scores.csv")

    query_pharmacophore = _build_query_pharmacophore(query_sdf, query_pml_path)

    if psd_path and os.path.exists(psd_path):
        message_source = f"reused existing screening database at {psd_path}"
    else:
        if not candidate_ligand_df_file:
            raise ValueError(
                "candidate_ligand_df_file is required when psd_path is not "
                "provided or does not exist yet."
            )

        df_ligand = pd.read_csv(candidate_ligand_df_file)
        candidate_sdf_files = df_ligand['ligand_sdf_file'].tolist()

        merged_sdf_path = os.path.join(file_dir, f"{project_name}_candidates_merged.sdf")
        with open(merged_sdf_path, 'w') as outfile:
            for sdf_file in candidate_sdf_files:
                with open(sdf_file, 'r') as infile:
                    outfile.write(infile.read())

        psd_path = psd_path or os.path.join(file_dir, f"{project_name}_candidates.psd")
        _build_psd_database(merged_sdf_path, psd_path, psdcreate_bin)
        message_source = f"built a new screening database at {psd_path}"

    scores = _align_library_to_query(query_pharmacophore, psd_path, default_score)

    ranked = sorted(scores.items(), key=lambda kv: kv[1], reverse=True)
    df_scores = pd.DataFrame(ranked, columns=['candidate_name', 'score'])
    df_scores.to_csv(scores_csv_path, index=False)

    return {
        'message': f"Screened {len(scores)} candidates against the query pharmacophore ({message_source}).",
        'scores': dict(ranked),
        'ranked_hits': ranked,
        'saved_file': scores_csv_path,
        'psd_file': psd_path,
        'query_pharmacophore_file': query_pml_path,
    }
