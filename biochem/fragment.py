import string
from functools import lru_cache

from rdkit import Chem


FUNCTIONAL_GROUP_SMARTS = {
    "Carboxylic acid": "[CX3](=O)[OX2H1]",
    "Ester": "[#6][CX3](=[OX1])[OX2H0;$([OX2H0][#6])]",
    "Amide": "[NX3][CX3](=[OX1])[#6]",
    "Aldehyde": "[CX3H1](=O)[#6]",
    "Ketone": "[#6][CX3](=O)[#6]",
    "Ether": "[OX2H0;!$([OX2H0]C=O)]([#6])[#6]",
    "Alcohol": "[#6;!$([CX3]=[OX1])][OX2H]",
    "Amine": "[NX3;!$(NC=[OX1]);!$(N=*);!$([N-])]",
    "Nitrile": "[NX1]#[CX2]",
    "Nitro": "[$([NX3](=O)=O),$([NX3+](=O)[O-])]",
    "Benzene ring": "c1ccccc1",
    "Alkene": "[CX3]=[CX3]",
    "Alkyne": "[CX2]#[CX2]",
    "Acyl halide": "[CX3](=[OX1])[F,Cl,Br,I]",
    "Sulfonamide": "[#16X4](=[OX1])(=[OX1])[NX3]",
    "Sulfonic acid": "[#16X4](=[OX1])(=[OX1])[OX2H1]",
    "Thiol": "[SX2H]",
    "Carbamate": "[NX3][CX3](=[OX1])[OX2H0]",
    "Urea": "[NX3][CX3](=[OX1])[NX3]",
    "Halide": "[F,Cl,Br,I]",
}


@lru_cache(maxsize=1)
def _functional_group_mols() -> dict:
    """Compile the curated functional-group SMARTS patterns into RDKit Mols, cached after first call"""
    mols = {}
    for name, smarts in FUNCTIONAL_GROUP_SMARTS.items():
        mol = Chem.MolFromSmarts(smarts)
        if mol is not None:
            mols[name] = mol
    return mols


def _consume_ring_closures(smiles: str, i: int) -> tuple:
    """Consume ring-closure digits (or %nn) immediately following an atom; returns (positions, next_index)"""
    positions = []
    while i < len(smiles):
        if smiles[i] == "%" and smiles[i + 1:i + 3].isdigit():
            positions.extend([i, i + 1, i + 2])
            i += 3
        elif smiles[i].isdigit():
            positions.append(i)
            i += 1
        else:
            break
    return positions, i


def map_atoms_to_smiles_chars(smiles: str) -> dict:
    """Map each RDKit atom index to the character positions (symbol + ring-closure digits) it occupies in the SMILES string"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}

    atom_to_chars = {}
    atom_idx = 0
    i = 0
    while i < len(smiles):
        char = smiles[i]

        if char == "[":
            end = smiles.find("]", i)
            if end == -1:
                i += 1
                continue
            chars = list(range(i, end + 1))
            i = end + 1
            ring_chars, i = _consume_ring_closures(smiles, i)
            if atom_idx < mol.GetNumAtoms():
                atom_to_chars[atom_idx] = chars + ring_chars
                atom_idx += 1

        elif smiles[i:i + 2] in ("Cl", "Br"):
            chars = [i, i + 1]
            i += 2
            ring_chars, i = _consume_ring_closures(smiles, i)
            if atom_idx < mol.GetNumAtoms():
                atom_to_chars[atom_idx] = chars + ring_chars
                atom_idx += 1

        elif char in "BCNOPSFIHbcnops":
            chars = [i]
            i += 1
            ring_chars, i = _consume_ring_closures(smiles, i)
            if atom_idx < mol.GetNumAtoms():
                atom_to_chars[atom_idx] = chars + ring_chars
                atom_idx += 1

        else:
            i += 1

    return atom_to_chars


def get_functional_groups(smiles: str) -> list:
    """Find common functional groups in a SMILES string, with the atom and character positions each occupies"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    atom_to_chars = map_atoms_to_smiles_chars(smiles)

    matches = []
    for name, fg_mol in _functional_group_mols().items():
        for match in mol.GetSubstructMatches(fg_mol):
            char_positions = sorted({
                pos
                for atom_idx in match
                for pos in atom_to_chars.get(atom_idx, [])
            })
            if char_positions:
                span = (char_positions[0], char_positions[-1])
                matches.append({
                    "name": name,
                    "atom_indices": match,
                    "char_positions": char_positions,
                    "span": span,
                    "substring": smiles[span[0]:span[1] + 1],
                })

    return matches


def annotate_functional_groups(smiles: str) -> str:
    """Render a SMILES string with a letter under each character marking which functional group it belongs to, plus a legend"""
    matches = get_functional_groups(smiles)
    matches.sort(key=lambda m: m["span"][0])

    alphabet = string.ascii_uppercase + string.ascii_lowercase + string.digits
    letters = [alphabet[i % len(alphabet)] for i in range(len(matches))]

    labels = [" "] * len(smiles)
    for match, letter in sorted(
        zip(matches, letters),
        key=lambda pair: pair[0]["span"][1] - pair[0]["span"][0],
        reverse=True,
    ):
        start, end = match["span"]
        for pos in range(start, end + 1):
            labels[pos] = letter

    legend = [
        f"{match['substring']} : {match['name']} -> {letter}"
        for match, letter in zip(matches, letters)
    ]

    return "\n".join([smiles, "".join(labels), *legend])


if __name__ == "__main__":
    print(annotate_functional_groups("CC(=O)OC1=CC=CC=C1C(=O)O"))
