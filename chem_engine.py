"""
chem_engine.py — the actual chemistry
======================================
Every molecule shown by BondVision AR — whether picked from the curated
library or built atom-by-atom — is turned into 3D coordinates by the SAME
functions here. That keeps the two features consistent with each other by
construction, instead of risking two copies of the geometry logic quietly
disagreeing.

Scope, stated plainly: this models simple AXn molecules (one central atom,
n identical peripheral atoms, single bonds unless a curated bond_order says
otherwise). Real molecules with mixed substituents or branching chains
(ethanol, benzene, DDT, ...) are genuinely out of scope for this generator —
see the Molecule Library for real, more complex molecules described in text.
"""

from __future__ import annotations  # lets `float | None` work on Python 3.9/3.10 too

import math
from dataclasses import dataclass, field

import numpy as np

from chem_data import (
    ELEMENTS, GEOMETRY_INFO, EXACT_ANGLES, ION_CHARGES, GREEN_FACTS,
    GREEN_FALLBACK_NOTE, ROMAN, CENTRAL_ALLOWED_N, PERIPHERAL_ELEMENTS,
    IONIC_METALS, IONIC_ANIONS, CURATED_MOLECULES, COMPLEX_MOLECULES,
    COMPLEX_BY_FORMULA, ANION_STEM,
)

BOND_EN_CUTOFFS = (0.4, 1.7)  # commonly used classroom thresholds (Pauling scale)


@dataclass
class Atom:
    element: str
    pos: np.ndarray
    role: str  # "central" | "peripheral"


@dataclass
class Bond:
    i: int  # index into atoms list
    j: int
    order: int


@dataclass
class LonePair:
    pos: np.ndarray


@dataclass
class MoleculeGeometry:
    formula: str
    atoms: list
    bonds: list
    lone_pairs: list
    geometry_name: str
    geometry_name_hi: str
    bond_angle: float | None
    steric: int | None
    lone_pair_count: int
    bond_type: str          # "Nonpolar covalent" | "Polar covalent" | "Ionic"
    molecule_polarity: str  # "Polar" | "Nonpolar" | "Ionic — not applicable"
    molar_mass: float
    is_ionic: bool = False
    notes: list = field(default_factory=list)
    meta: dict = field(default_factory=dict)


# =============================================================================
# GEOMETRY VECTOR SETS
# =============================================================================
def _norm(v):
    v = np.array(v, dtype=float)
    return v / np.linalg.norm(v)


def _tetrahedral_vectors():
    raw = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
    return [_norm(v) for v in raw]


def _trigonal_vectors():
    return [_norm((math.cos(math.radians(a)), math.sin(math.radians(a)), 0))
            for a in (90, 210, 330)]


def _linear_vectors():
    return [np.array([0, 0, 1.0]), np.array([0, 0, -1.0])]


def _trigonal_bipyramidal_vectors():
    eq = [_norm((math.cos(math.radians(a)), math.sin(math.radians(a)), 0))
          for a in (0, 120, 240)]
    ax = [np.array([0, 0, 1.0]), np.array([0, 0, -1.0])]
    return eq + ax


def _octahedral_vectors():
    return [np.array(v, dtype=float) for v in
            [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]]


def _bent_vectors(angle_deg):
    """2 vectors in the xy-plane, symmetric about +x, exactly angle_deg apart."""
    half = math.radians(angle_deg / 2)
    return [np.array([math.cos(half), math.sin(half), 0]),
            np.array([math.cos(half), -math.sin(half), 0])]


def _tripod_vectors(angle_deg):
    """
    3 vectors arranged like a tripod (equal polar angle from +z, 120 deg
    apart in azimuth) such that the angle between any two is angle_deg.
    Solved from the spherical dot-product identity:
        cos(angle) = cos^2(phi) - 0.5 * sin^2(phi)
    """
    c = math.cos(math.radians(angle_deg))
    cos2phi = (c + 0.5) / 1.5
    cos2phi = max(0.0, min(1.0, cos2phi))
    cosphi = math.sqrt(cos2phi)
    sinphi = math.sqrt(max(0.0, 1 - cos2phi))
    vecs = []
    for k in range(3):
        theta = math.radians(120 * k)
        vecs.append(np.array([sinphi * math.cos(theta), sinphi * math.sin(theta), cosphi]))
    return vecs


def geometry_vectors(steric, lone_pairs, exact_angle=None):
    """
    Return (bonding_vectors, lone_pair_vectors, ideal_angle_used) for a given
    electron-domain count. exact_angle overrides the idealised value with a
    real, textbook-confirmed bond angle for a specific famous molecule.
    """
    key = (steric, lone_pairs)
    if key == (2, 0):
        return _linear_vectors(), [], 180.0
    if key == (3, 0):
        return _trigonal_vectors(), [], 120.0
    if key == (4, 0):
        return _tetrahedral_vectors(), [], 109.5
    if key == (4, 1):
        angle = exact_angle or 109.5
        return _tripod_vectors(angle), [_tetrahedral_vectors()[3]], angle
    if key == (4, 2):
        angle = exact_angle or 104.5
        lp_dirs = [v for v in _tetrahedral_vectors()][2:]  # 2 "spare" directions
        return _bent_vectors(angle), lp_dirs, angle
    if key == (3, 1):
        # Bent, derived from trigonal-planar electron geometry with 1 lone
        # pair (e.g. ozone, sulfur dioxide — both resonance structures).
        angle = exact_angle or 120.0
        return _bent_vectors(angle), [np.array([-1.0, 0.0, 0.0])], angle
    if key == (5, 0):
        return _trigonal_bipyramidal_vectors(), [], 90.0
    if key == (6, 0):
        return _octahedral_vectors(), [], 90.0
    raise ValueError(f"Unsupported electron geometry: steric={steric}, lone_pairs={lone_pairs}")


# =============================================================================
# BOND CLASSIFICATION
# =============================================================================
def classify_bond(en_a, en_b, force_covalent=False):
    """
    The 0.4 / 1.7 electronegativity-difference cutoff is a common classroom
    rule of thumb, and it is genuinely imperfect at the edges. Two textbook
    exceptions: HF (EN difference 1.78) and BF3 (1.94) both sit past the
    "ionic" line, yet both are real, well-known MOLECULAR covalent
    compounds — HF is a covalent gas/liquid, and BF3 is a covalent gas used
    as a Lewis-acid catalyst. Neither hydrogen nor boron forms true ionic
    bonds with a nonmetal this way.

    More generally: anything built via build_covalent() below was already
    chosen to be modelled as a discrete molecule with a VSEPR shape — and a
    VSEPR shape is only a meaningful concept for a covalent molecule in the
    first place. So build_covalent always passes force_covalent=True, and
    "Ionic" is reserved for the separate Ionic Compound builder, which
    models real metal–nonmetal lattices (NaCl, MgCl2, ...) instead.
    """
    d = abs(en_a - en_b)
    if force_covalent:
        return ("Nonpolar covalent" if d <= BOND_EN_CUTOFFS[0] else "Polar covalent"), d
    if d > BOND_EN_CUTOFFS[1]:
        return "Ionic", d
    if d > BOND_EN_CUTOFFS[0]:
        return "Polar covalent", d
    return "Nonpolar covalent", d


def molecule_polarity(bond_type, lone_pairs, n_peripheral, same_peripheral=True):
    """Overall molecular polarity from symmetry (identical terminal atoms only)."""
    if bond_type == "Ionic":
        return "Ionic — not a discrete polar/nonpolar molecule"
    if n_peripheral == 1:
        # heteronuclear diatomic is polar unless EN difference is ~0
        return "Polar" if bond_type == "Polar covalent" else "Nonpolar"
    if lone_pairs > 0:
        return "Polar"        # asymmetric electron distribution
    if same_peripheral:
        return "Nonpolar"     # symmetric arrangement, dipoles cancel
    return "Polar"


# =============================================================================
# BUILD A COVALENT AXn MOLECULE
# =============================================================================
def build_covalent(central, peripheral, n, bond_order=1, bond_length=1.35,
                    steric_override=None, lone_pairs_override=None,
                    exact_angle=None, formula=None):
    c = ELEMENTS[central]
    p = ELEMENTS[peripheral]

    if n == 1:
        # Diatomic: no VSEPR angle to speak of.
        pos_c = np.array([0.0, 0.0, bond_length / 2])
        pos_p = np.array([0.0, 0.0, -bond_length / 2])
        atoms = [Atom(central, pos_c, "central"), Atom(peripheral, pos_p, "peripheral")]
        bonds = [Bond(0, 1, bond_order)]
        bond_t, _ = classify_bond(c["en_neg"], p["en_neg"], force_covalent=True)
        polarity = molecule_polarity(bond_t, 0, 1)
        mass = c["mass"] + p["mass"]
        return MoleculeGeometry(
            formula=formula or f"{central}{peripheral}", atoms=atoms, bonds=bonds,
            lone_pairs=[], geometry_name="Diatomic (linear)", geometry_name_hi="द्विपरमाण्विक (रैखिक)",
            bond_angle=None, steric=None, lone_pair_count=0,
            bond_type=bond_t, molecule_polarity=polarity, molar_mass=round(mass, 2),
        )

    if steric_override is not None:
        steric = steric_override
        lone_pairs = lone_pairs_override or 0
    else:
        lone_pairs = (c["valence_e"] - n * bond_order) / 2
        if lone_pairs != int(lone_pairs) or lone_pairs < 0:
            raise ValueError(f"{central} cannot cleanly form {n} bonds of order {bond_order}")
        lone_pairs = int(lone_pairs)
        steric = n + lone_pairs

    bond_vecs, lp_vecs, angle_used = geometry_vectors(steric, lone_pairs, exact_angle)
    geom = GEOMETRY_INFO.get((steric, lone_pairs), {"name": "Custom", "name_hi": "कस्टम"})

    atoms = [Atom(central, np.array([0.0, 0.0, 0.0]), "central")]
    bonds = []
    for k in range(n):
        pos = bond_vecs[k] * bond_length
        atoms.append(Atom(peripheral, pos, "peripheral"))
        bonds.append(Bond(0, len(atoms) - 1, bond_order))

    lone_pairs_out = [LonePair(v * bond_length * 0.55) for v in lp_vecs]

    bond_t, _ = classify_bond(c["en_neg"], p["en_neg"], force_covalent=True)
    polarity = molecule_polarity(bond_t, lone_pairs, n, same_peripheral=True)
    mass = c["mass"] + n * p["mass"]

    return MoleculeGeometry(
        formula=formula or f"{central}{peripheral}{n if n > 1 else ''}",
        atoms=atoms, bonds=bonds, lone_pairs=lone_pairs_out,
        geometry_name=geom["name"], geometry_name_hi=geom["name_hi"],
        bond_angle=angle_used, steric=steric, lone_pair_count=lone_pairs,
        bond_type=bond_t, molecule_polarity=polarity, molar_mass=round(mass, 2),
    )


# =============================================================================
# BUILD AN IONIC COMPOUND  (charge-balanced formula unit, lattice-style view)
# =============================================================================
def ionic_formula(cation, anion, cation_charge=None, anion_charge=None):
    cation_charge = cation_charge or ION_CHARGES[cation][0]
    anion_charge = anion_charge or ION_CHARGES[anion][0]
    qc, qa = cation_charge, abs(anion_charge)
    g = math.gcd(qc, qa)
    n_cation, n_anion = qa // g, qc // g
    return n_cation, n_anion


def build_ionic(cation, anion, cation_charge=None, spacing=1.9):
    cation_charge = cation_charge or ION_CHARGES[cation][0]
    n_cat, n_an = ionic_formula(cation, anion, cation_charge=cation_charge)
    c, a = ELEMENTS[cation], ELEMENTS[anion]

    total = n_cat + n_an
    atoms = []
    # Simple alternating row — a schematic "slice" of a real 3D ionic
    # lattice, not a literal unit cell.
    xs = np.linspace(-(total - 1) / 2, (total - 1) / 2, total) * spacing
    seq = []
    ci, ai = 0, 0
    for i in range(total):
        if i % 2 == 0 and ci < n_cat:
            seq.append(cation); ci += 1
        elif ai < n_an:
            seq.append(anion); ai += 1
        else:
            seq.append(cation); ci += 1

    for x, sym in zip(xs, seq):
        atoms.append(Atom(sym, np.array([x, 0.0, 0.0]), "central" if sym == cation else "peripheral"))

    formula_str = f"{cation}{'' if n_cat == 1 else n_cat}{anion}{'' if n_an == 1 else n_an}"
    mass = n_cat * c["mass"] + n_an * a["mass"]

    chem_name = None
    if len(ION_CHARGES[cation]) > 1:
        stem = ANION_STEM.get(anion, ELEMENTS[anion]["name"].lower())
        chem_name = f"{ELEMENTS[cation]['name']}({ROMAN[cation_charge]}) {stem}"

    return MoleculeGeometry(
        formula=formula_str, atoms=atoms, bonds=[], lone_pairs=[],
        geometry_name="Ionic lattice (schematic slice)", geometry_name_hi="आयनिक जालक",
        bond_angle=None, steric=None, lone_pair_count=0,
        bond_type="Ionic", molecule_polarity="Ionic — not a discrete polar/nonpolar molecule",
        molar_mass=round(mass, 2), is_ionic=True,
        meta={"cation": cation, "anion": anion, "cation_charge": cation_charge,
              "n_cation": n_cat, "n_anion": n_an, "chem_name": chem_name},
        notes=["Ionic compounds form repeating 3D crystal lattices, not discrete "
               "molecules with a VSEPR shape — this view is a simplified schematic slice."],
    )


# =============================================================================
# GREEN CHEMISTRY SCORE
# =============================================================================
_LEVEL_SCORE = {"Low": 90, "Medium": 60, "High": 30, "Very High": 15}
_CLIMATE_SCORE = {"None": 95, "Minor": 70, "Greenhouse gas": 35, "Ozone-depleting": 15}


def green_score(formula, has_halogen=False, is_acid_base=False):
    facts = GREEN_FACTS.get(formula)
    if facts:
        parts = {
            "Persistence": _LEVEL_SCORE[facts["persistence"]],
            "Toxicity": _LEVEL_SCORE[facts["toxicity"]],
            "Climate impact": _CLIMATE_SCORE[facts["climate"]],
        }
        # The score is the WEAKEST factor, not the average. Averaging would
        # let one severe hazard get diluted by two unrelated good scores —
        # e.g. fluorine gas (F2) is famously one of the most hazardous
        # substances known, but is neither persistent nor a greenhouse gas;
        # averaging would still call it "generally low concern", which is
        # exactly the kind of wrong reassurance a safety-facing score must
        # not give. A chain is only as strong as its weakest link.
        weakest = min(parts, key=parts.get)
        score = parts[weakest]
        return {
            "score": score, "grade": _grade(score), "verified": True,
            "persistence": facts["persistence"], "toxicity": facts["toxicity"],
            "climate": facts["climate"], "note": facts["note"],
            "weakest_factor": weakest,
        }
    # Fallback structural heuristic for anything not in the curated table.
    score = 70
    if has_halogen:
        score -= 25
    if is_acid_base:
        score -= 15
    score = max(10, min(95, score))
    return {
        "score": score, "grade": _grade(score), "verified": False,
        "persistence": "Unknown", "toxicity": "Unknown", "climate": "Unknown",
        "note": GREEN_FALLBACK_NOTE,
    }


def _grade(score):
    if score >= 80:
        return "A+ — environmentally favourable"
    if score >= 65:
        return "A — generally low concern"
    if score >= 45:
        return "B — handle with caution"
    return "C — significant environmental or safety concern"


# =============================================================================
# ATOM WORKSPACE MATCHER
# =============================================================================
# The core of "tap atoms, get the real chemistry automatically". Given what
# a student has tapped into the workspace, this returns the best HONEST
# answer available — never a fabricated structure.
from collections import Counter


def _curated_atom_counts(m):
    counts = Counter()
    if m.get("is_ionic"):
        n_cat, n_an = ionic_formula(m["cation"], m["anion"])
        counts[m["cation"]] += n_cat
        counts[m["anion"]] += n_an
    else:
        counts[m["central"]] += 1
        counts[m["peripheral"]] += m["n"]
    return dict(counts)


def _build_from_curated(m):
    if m.get("is_ionic"):
        return build_ionic(m["cation"], m["anion"])
    kwargs = dict(central=m["central"], peripheral=m["peripheral"], n=m["n"],
                  bond_order=m["bond_order"], formula=m["formula"])
    if "steric_override" in m:
        kwargs["steric_override"] = m["steric_override"]
        kwargs["lone_pairs_override"] = m["lone_pairs_override"]
    if m["formula"] in EXACT_ANGLES:
        kwargs["exact_angle"] = EXACT_ANGLES[m["formula"]]
    return build_covalent(**kwargs)


def match_atoms(atom_counts):
    """
    atom_counts: {element_symbol: count} of whatever is currently in the
    workspace basket. Returns a dict describing the best real answer:

      kind == "empty"            nothing tapped yet
      kind == "curated"          exact match to a library molecule with a
                                  derivable 3D shape -> includes "geo"
      kind == "complex"          exact match to a real molecule this tool
                                  can't derive a structure for (glucose,
                                  benzene, ...) -> includes "molecule" (fact
                                  card only, no "geo")
      kind == "generic_ionic"    a metal + nonmetal ratio that correctly
                                  charge-balances, even though it isn't one
                                  of the specifically curated examples
      kind == "generic_covalent" a central/peripheral ratio the VSEPR engine
                                  can validly build, even though it isn't
                                  specifically curated
      kind == "no_match"         a combination with no verified real
                                  structure in this tool -> includes a
                                  computed total mass so it's still useful
    """
    counts = {k: v for k, v in atom_counts.items() if v > 0}
    if not counts:
        return {"kind": "empty"}

    for m in CURATED_MOLECULES:
        if _curated_atom_counts(m) == counts:
            return {"kind": "curated", "molecule": m, "geo": _build_from_curated(m)}

    for m in COMPLEX_MOLECULES:
        if m["atoms"] == counts:
            return {"kind": "complex", "molecule": m}

    if len(counts) == 2:
        (e1, n1), (e2, n2) = list(counts.items())

        for cat, cn, an, an_n in [(e1, n1, e2, n2), (e2, n2, e1, n1)]:
            if cat in IONIC_METALS and an in IONIC_ANIONS:
                for charge in ION_CHARGES[cat]:
                    n_cat, n_an = ionic_formula(cat, an, cation_charge=charge)
                    if n_cat == cn and n_an == an_n:
                        return {"kind": "generic_ionic", "geo": build_ionic(cat, an, cation_charge=charge)}

        for cen, cen_n, per, per_n in [(e1, n1, e2, n2), (e2, n2, e1, n1)]:
            if (cen_n == 1 and cen in CENTRAL_ALLOWED_N and per in PERIPHERAL_ELEMENTS
                    and per_n in CENTRAL_ALLOWED_N[cen]):
                try:
                    geo = build_covalent(central=cen, peripheral=per, n=per_n, bond_order=1)
                    if geo.bond_type != "Ionic":
                        return {"kind": "generic_covalent", "geo": geo}
                except Exception:
                    pass

    total_mass = round(sum(ELEMENTS[e]["mass"] * n for e, n in counts.items()), 2)
    return {"kind": "no_match", "counts": counts, "total_mass": total_mass}


def hill_formula(counts):
    """
    Standard Hill notation: Carbon first (if present), then Hydrogen, then
    every other element alphabetically — the convention chemical databases
    (including PubChem's formula search) expect. Counts of 1 are omitted.
    """
    counts = {k: v for k, v in counts.items() if v > 0}
    order = []
    if "C" in counts:
        order.append("C")
        if "H" in counts:
            order.append("H")
    remaining = sorted(e for e in counts if e not in order)
    order += remaining
    return "".join(f"{e}{counts[e] if counts[e] != 1 else ''}" for e in order)
