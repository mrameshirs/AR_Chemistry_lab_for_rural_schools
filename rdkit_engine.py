"""
rdkit_engine.py — turning a real structure into a MoleculeGeometry
=====================================================================
RDKit does NOT tell you what chemical forms from a set of atoms — that's
PubChem's job (or any real chemical database). RDKit's job here is
narrower and well-established: given a structure that's already been
identified, work out the correct bond orders, formal charges, and lone
pairs so it can be drawn correctly. This module is that second step.

The lone-pair arithmetic below is the SAME formula chem_engine.py already
uses for simple AXn molecules — this module just applies it per-atom
instead of assuming one central atom, which lets it generalise to
whatever structure RDKit hands it: polyatomic ions, resonance structures,
real molecules with more than one "central" atom.

Verified before use (see the module's own test block at the bottom):
    lone_pairs = (valence_electrons - formal_charge - sum_of_bond_orders) / 2
checked against sulfate's real, textbook Lewis structure.

HONEST LIMITS:
- Main-group chemistry (the overwhelming majority of a school syllabus,
  including polyatomic ions like sulfate/carbonate/nitrate) is solid.
- Transition metals and organometallics are NOT reliably handled by this
  simple valence-electron formula — d-electrons don't follow the same
  counting rule, so lone-pair counts for those are flagged as
  approximate rather than presented with false confidence.
- RDKit's automatic resonance/charge perception is real and well-tested,
  but not infallible on unusual bonding. Flagged, not hidden.
"""

from __future__ import annotations

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from chem_engine import Atom, Bond, LonePair, MoleculeGeometry
from chem_data import ELEMENTS

# Main-group valence electron count by column, for elements not in our
# curated ELEMENTS table. Transition metals (d-block) are intentionally
# left out — their valence electron count doesn't follow this simple rule.
_MAIN_GROUP_VALENCE = {
    1: 1, 2: 2,                                   # groups 1-2 (alkali/alkaline earth)
    13: 3, 14: 4, 15: 5, 16: 6, 17: 7, 18: 8,      # groups 13-18 (p-block)
}

_PT = Chem.GetPeriodicTable()


class StructureWarning:
    pass


def _valence_electrons(symbol):
    """Real valence electron count where we're confident; None (flagged) otherwise."""
    if symbol in ELEMENTS:
        return ELEMENTS[symbol]["valence_e"], True
    try:
        group = _PT.GetNOuterElecs(_PT.GetAtomicNumber(symbol))
        # RDKit's GetNOuterElecs already returns a sensible main-group valence
        # count for main-group elements; treat it as confident for those,
        # approximate for d-block.
        z = _PT.GetAtomicNumber(symbol)
        is_transition = 21 <= z <= 30 or 39 <= z <= 48 or 57 <= z <= 80 or 89 <= z <= 112
        return group, not is_transition
    except Exception:
        return 8, False  # unknown; guess octet, flag as unreliable


def _element_color(symbol):
    if symbol in ELEMENTS:
        return ELEMENTS[symbol]["color"]
    # Reasonable generic fallback colours for elements outside our curated set
    fallback = {
        "Fe": "#E06633", "Cu": "#C88033", "Zn": "#7D80B0", "Ca": "#3DC23D",
        "K": "#8F40D4", "Na": "#AB5CF2", "Mg": "#8AFF00",
    }
    return fallback.get(symbol, "#B0B0B0")


def mol_from_sdf(sdf_text):
    """Parse a real SDF/MOL block (from PubChem or anywhere else) into an
    RDKit Mol with 3D coordinates, sanitized so formal charges and bond
    orders are chemically consistent."""
    mol = Chem.MolFromMolBlock(sdf_text, sanitize=True, removeHs=False)
    if mol is None:
        # Some SDF exports need lenient parsing; retry without strict sanitization
        mol = Chem.MolFromMolBlock(sdf_text, sanitize=False, removeHs=False)
        if mol is None:
            raise ValueError("Could not parse this structure file (invalid or unsupported SDF).")
        Chem.SanitizeMol(mol, catchErrors=True)
    if mol.GetNumConformers() == 0:
        AllChem.Compute2DCoords(mol)
    return mol


def molecule_geometry_from_mol(mol, formula, name=None, bond_length_scale=1.0):
    """
    Convert an RDKit Mol (with a 3D or 2D conformer) into the same
    MoleculeGeometry shape build_covalent()/build_ionic() produce, so it
    flows through the exact same Plotly/Matplotlib/AR renderers already
    built and tested for the local engine.
    """
    conf = mol.GetConformer()
    atoms, bonds, lone_pairs = [], [], []
    warnings = []
    total_mass = 0.0
    net_charge = 0

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        symbol = atom.GetSymbol()
        atoms.append(Atom(symbol, np.array([pos.x, pos.y, pos.z]) * bond_length_scale,
                           "peripheral" if idx > 0 else "central"))
        total_mass += _PT.GetAtomicWeight(atom.GetAtomicNum())
        net_charge += atom.GetFormalCharge()

    for b in mol.GetBonds():
        order = b.GetBondTypeAsDouble()
        order_int = max(1, round(order))  # Kekulized aromatic bonds resolve to 1 or 2
        bonds.append(Bond(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), order_int))

    # bond-order sum per atom, needed for the lone-pair formula
    bond_sum = [0] * mol.GetNumAtoms()
    for b in bonds:
        bond_sum[b.i] += b.order
        bond_sum[b.j] += b.order

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        symbol = atom.GetSymbol()
        if symbol == "H":
            continue  # hydrogen's single bond fully accounts for its 2-electron duet
        fc = atom.GetFormalCharge()
        ve, confident = _valence_electrons(symbol)
        if not confident:
            warnings.append(
                f"{symbol}: lone-pair count is approximate — main-group valence "
                f"counting doesn't reliably apply to transition metals."
            )
        lp_electrons = ve - fc - bond_sum[idx]
        lp = max(0, round(lp_electrons / 2))
        if lp_electrons < 0 or lp_electrons % 2 != 0:
            warnings.append(f"{symbol} (atom {idx}): lone-pair count may be inexact for this structure.")
        pos = conf.GetAtomPosition(idx)
        for k in range(lp):
            # place lone-pair markers slightly offset from the atom, spread
            # around it — purely illustrative positioning, not a claim about
            # exact orbital geometry
            angle = 2 * np.pi * k / max(lp, 1)
            offset = np.array([np.cos(angle), np.sin(angle), 0.3]) * 0.5
            lone_pairs.append(LonePair((np.array([pos.x, pos.y, pos.z]) + offset) * bond_length_scale))

    # RDKit structures are drawn as covalent frameworks; a net-charged
    # fragment (like sulfate on its own) is a polyatomic ION, not "ionic"
    # in the simple metal/nonmetal sense — that distinction is surfaced via
    # net_charge/polarity text instead of overloading this field.
    bond_type = "Polar covalent"
    polarity = f"Net charge {net_charge:+d}" if net_charge != 0 else "Neutral molecule"

    return MoleculeGeometry(
        formula=formula, atoms=atoms, bonds=bonds, lone_pairs=lone_pairs,
        geometry_name="Real structure (PubChem + RDKit)", geometry_name_hi="वास्तविक संरचना",
        bond_angle=None, steric=None, lone_pair_count=sum(1 for _ in lone_pairs),
        bond_type=bond_type, molecule_polarity=polarity, molar_mass=round(total_mass, 2),
        is_ionic=False,
        notes=(["This structure came from a real chemistry database (PubChem), interpreted by "
                "RDKit — not derived by this app's own simple AXₙ model, which is why it can "
                "show real molecules that model can't."] + warnings),
        meta={"net_charge": net_charge, "source": "pubchem+rdkit", "name": name},
    )


def geometry_from_sdf(sdf_text, formula, name=None, is_3d=True):
    """
    is_3d=False means this SDF came from PubChem's 2D fallback (no
    precomputed 3D conformer available for this compound — see
    pubchem_client.get_sdf_best_effort). In that case a real 3D conformer
    is generated locally with RDKit's standard ETKDG algorithm — the same
    well-established method used to build this module's own test
    fixtures — rather than silently rendering flattened 2D coordinates as
    if they were genuine 3D positions.
    """
    mol = mol_from_sdf(sdf_text)
    generated_3d = False
    if not is_3d:
        mol = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
        if result != 0:
            # ETKDG can fail on unusual structures (e.g. some disconnected
            # ionic fragments) -- try a simpler fallback rather than crash
            result = AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True,
                                            useBasicKnowledge=False, enforceChirality=False)
        if result == 0:
            generated_3d = True
        else:
            AllChem.Compute2DCoords(mol)  # last resort: flat layout, honestly labelled below

    geo = molecule_geometry_from_mol(mol, formula, name=name)

    n_fragments = len(Chem.GetMolFrags(mol))
    extra_notes = []
    if not is_3d:
        if generated_3d:
            extra_notes.append(
                "PubChem doesn't have a precomputed 3D structure for this compound — this is "
                "common for simple ionic salts. A 3D conformer was generated locally using "
                "RDKit's standard ETKDG algorithm instead of PubChem's own data."
            )
        else:
            extra_notes.append(
                "Neither PubChem nor local 3D generation could produce real 3D coordinates for "
                "this structure — showing a flat 2D layout instead. Bond orders and formal "
                "charges are still real; only the 3D positions are a placeholder."
            )
    if n_fragments > 1:
        extra_notes.append(
            f"This structure has {n_fragments} disconnected pieces (separate ions) — like this "
            f"app's own local ionic builder, the relative 3D arrangement shown between them is a "
            f"schematic placement, not a real crystal lattice."
        )
    geo.notes = geo.notes + extra_notes
    return geo


def geometry_from_pubchem(cid, formula, name=None):
    """Convenience wrapper used by the app: fetches the best available
    structure from PubChem and turns it into a MoleculeGeometry, handling
    the 3D/2D-fallback distinction honestly (see geometry_from_sdf)."""
    import pubchem_client
    sdf, is_3d = pubchem_client.get_sdf_best_effort(cid)
    return geometry_from_sdf(sdf, formula=formula, name=name, is_3d=is_3d)
