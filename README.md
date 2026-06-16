# biochem

Cheminformatics and structural biology tools for BioChemAIgent and Drug Discovery Platform. 

## Install

```bash
micromamba env create -f environment.yml -y
micromamba env create -f environment.yml -y --prefix ~/env/bcai
pip install -e .
```

## Key tools

- `render_structures` — render PDB/CIF/SDF/MOL2 files to interactive 3D HTML (py3Dmol). Exposes `style_rules`, `surface_rules`, `label_rules`, `chain_color_map`, `interaction_tables`, box constraint overlay.
- `get_protein_ligand_interaction` / `interaction_plot` — compute and visualise protein-ligand interaction networks.
- `run_af3` / `run_esm3` — AlphaFold3 and ESM3 structure prediction.
- `molecular_docking` tools — protein/ligand preparation, VINA/SMINA/GNINA/DiffDock docking.
- `protein_ligand_basics` — SMILES→3D, protonation, PDB extraction, ADMET prediction.

## Notes

- `render_structures.style_rules` and `surface_rules` use `TypedDict` (`StyleRule`, `SurfaceRule`) so FastMCP generates a schema with named properties. This is required for LLMs to correctly populate nested selection/style dicts — bare `Dict[str, Any]` causes models to emit empty `{}` objects.
- Cartoon style without an explicit `color` renders white on the default white background. Always include `"color": "spectrum"` or another colour when specifying cartoon style rules.

## License

This software is licensed under the GNU Affero General Public License v3.0 (AGPL-3.0).

You are free to use, modify, and distribute this software for non-commercial purposes,
provided that any derivative works are also licensed under AGPL-3.0 and made publicly
available.

Commercial use — including use in proprietary pipelines, SaaS products, or any
revenue-generating activity — requires a separate commercial license.
Contact: yousefi.bme@gmail.com, sbonn@uke.de

SPDX-License-Identifier: AGPL-3.0-only
