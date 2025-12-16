from __future__ import annotations

from .config import DockingConfig
from .scoring.terms import make_default_vina_terms
from .scoring.weighted_terms import WeightedTerms
from .scoring.precalculate import Precalculate, eval_interacting_pairs
from .pdbqt import parse_complex


def _infer_lig_name(ligand_path: Path) -> str:
    """
    Mirror lig_name extraction logic in main.cpp:
    it takes the directory name between first and second '/' of ligand path.

    For example:
        ../coreset_pdbqt/1bzc/1bzc_ligand.pdbqt  ->  '1bzc'
    """
    # Use POSIX-style parts to be robust across platforms
    parts = ligand_path.as_posix().split("/")
    if len(parts) >= 3:
        return parts[-2]
    # Fallback: strip extension from filename
    name = ligand_path.stem
    return name.split("_")[0] if "_" in name else name


def run_docking(cfg: DockingConfig) -> None:
    """
    High-level orchestration function, Python counterpart of:

        main() -> parse options -> main_procedure() -> do_search()

    At this stage this function is a scaffold: it prepares the same inputs
    (lig_name, search box, file paths) and will gradually be extended to call
    pure-Python implementations of:
      - PDBQT parsing and model construction
      - scoring (everything / weighted_terms / precalculate / igrid)
      - multi-objective search (LSHADE + NSGA2 / MODEA)
    """
    lig_name = _infer_lig_name(cfg.ligand)

    # 1) parse receptor & ligand PDBQT into a Python Model
    model = parse_complex(cfg.receptor, cfg.ligand)

    # 2) Construct a minimal vina-like scoring function with default weights
    terms = make_default_vina_terms()
    # weights follow everything.cpp defaults: [-0.035579, -0.005156, 0.840245, -0.035069, -0.587439, tors_div_weight]
    weights = [
        -0.035579,  # gauss1
        -0.005156,  # gauss2
        0.840245,   # repulsion
        -0.035069,  # hydrophobic
        -0.587439,  # non_dir_h_bond
    ]
    wt = WeightedTerms(terms=terms, weights=weights)
    prec = Precalculate(wt=wt)

    # 3) Evaluate a simple intramolecular energy for the ligand as a sanity check
    if model.ligands:
        intra_pairs = model.build_intra_ligand_pairs()
        e_intra = eval_interacting_pairs(prec, intra_pairs, model.atoms)
    else:
        e_intra = 0.0

    # 4) Evaluate a crude ligand-receptor interaction energy (all ligand atoms vs all grid atoms)
    e_inter = 0.0
    for idx, lig_atom in model.iter_ligand_atoms():
        for rec_atom in model.grid_atoms:
            e_inter += prec.eval_pair(lig_atom, rec_atom)

    # 4) TODO: build grids and run LSHADE+NSGA2 search.
    #    For now we stop after scoring sanity check.
    raise NotImplementedError(
        f"Python scoring skeleton ready. Ligand '{lig_name}' has "
        f"{model.num_movable_atoms} movable atoms; "
        f"intra-ligand rough energy = {e_intra:.3f}; "
        f"ligand-receptor rough energy = {e_inter:.3f}. "
        f"Next steps: implement full interacting_pairs, grids, and LSHADE+NSGA2 search."
    )


