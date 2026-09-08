"""
pubchem_client.py — real chemical lookup via PubChem's PUG-REST API
=====================================================================
This is what makes the Atom Workspace able to recognise combinations far
beyond the local curated/generic engine — Mg + S + 4 O correctly resolves
to real magnesium sulfate records, not a guess.

Every URL pattern here was checked against PubChem's own documentation
and a live worked example (aspirin, CID 2244) before being used, not
constructed from memory. See the comments at each function.

HONEST LIMITS, stated up front:
- This needs internet. It is the one part of BondVision AR that does.
- PubChem asks that callers not exceed 5 requests/second and 400/minute —
  fine for one student clicking a button, never for a hot loop.
- A common name is genuinely ambiguous. "Magnesium sulfate" alone matches
  at least 8 different real PubChem records (anhydrous, monohydrate,
  dihydrate, ... heptahydrate/Epsom salt, nonahydrate) because hydration
  state makes a different compound. This client surfaces the candidates
  and lets the caller choose — it never silently guesses one.
"""

import time

import requests

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
TIMEOUT = 15  # seconds; PubChem's own stated limit is 30s per request
_last_call = [0.0]


def _throttle(min_interval=0.25):
    """Never exceed ~4 requests/second — PubChem's own stated limit is 5/sec."""
    elapsed = time.time() - _last_call[0]
    if elapsed < min_interval:
        time.sleep(min_interval - elapsed)
    _last_call[0] = time.time()


def _get(url, parse_json=True):
    _throttle()
    try:
        r = requests.get(url, timeout=TIMEOUT)
    except requests.exceptions.RequestException as exc:
        raise PubChemError(f"Could not reach PubChem — check your internet connection. ({exc})")
    if r.status_code == 404:
        raise PubChemNotFound("No PubChem record found for that search.")
    if r.status_code >= 400:
        raise PubChemError(f"PubChem returned HTTP {r.status_code}: {r.text[:200]}")
    return r.json() if parse_json else r.text


class PubChemError(Exception):
    pass


class PubChemNotFound(PubChemError):
    pass


def _resolve_cids(data, url_for_polling_base, max_attempts=15, poll_interval=2.0):
    """
    Some PUG-REST searches (formula search, similarity search — anything
    requiring a scan across PubChem's ~110M compounds rather than a direct
    lookup) don't return results immediately. Instead of {"IdentifierList":
    {"CID": [...]}}, the first response is {"Waiting": {"ListKey": "..."}},
    and the real result has to be polled from a separate endpoint until it's
    ready. This was missed in the first version of this client, which is
    exactly why a real compound (e.g. forsterite/Mg2SiO4, a very common
    mineral) could come back as a false "not found" — the job key was
    silently discarded rather than followed up on. Confirmed against
    PubChem's own documented example (async similarity search) before this
    fix, not guessed.
    """
    if "IdentifierList" in data:
        return data["IdentifierList"]["CID"]

    list_key = data.get("Waiting", {}).get("ListKey")
    if not list_key:
        return []

    poll_url = f"{PUG}/compound/listkey/{list_key}/cids/JSON"
    for _ in range(max_attempts):
        time.sleep(poll_interval)
        polled = _get(poll_url)
        if "IdentifierList" in polled:
            return polled["IdentifierList"]["CID"]
        if "Waiting" not in polled:
            break  # unexpected shape; stop polling rather than loop forever
    return []


def search_by_formula(formula, max_results=8):
    """
    Look up candidate compounds by molecular formula.
    Confirmed real endpoint pattern (PubChem PUG-REST docs):
        {PUG}/compound/formula/{formula}/cids/JSON?MaxRecords=N
    This search is asynchronous (see _resolve_cids above) — the first
    response is a job key, not the answer.
    """
    url = f"{PUG}/compound/formula/{formula}/cids/JSON?MaxRecords={max_results}"
    data = _get(url)
    cids = _resolve_cids(data, url)
    if not cids:
        raise PubChemNotFound(f"No compound found with formula {formula}.")
    return cids[:max_results]


def search_by_name(name):
    """
    Confirmed real endpoint pattern, matches the worked aspirin example:
        {PUG}/compound/name/{name}/cids/JSON
    Name search is documented as an immediate lookup (unlike formula
    search), but this still goes through _resolve_cids defensively — no
    cost to it, and it means a future PubChem change to this endpoint's
    behaviour wouldn't silently reintroduce the same bug.
    """
    url = f"{PUG}/compound/name/{requests.utils.quote(name)}/cids/JSON"
    data = _get(url)
    cids = _resolve_cids(data, url)
    if not cids:
        raise PubChemNotFound(f"No compound found with name '{name}'.")
    return cids


def get_properties(cids):
    """
    Confirmed real endpoint pattern (verified live against CID 2244 / aspirin,
    which returned MolecularFormula 'C9H8O4' exactly as PubChem's own docs
    show):
        {PUG}/compound/cid/{cids}/property/{props}/JSON
    Returns a list of dicts, one per CID, with CID/MolecularFormula/
    MolecularWeight/IUPACName/CanonicalSMILES.
    """
    ids = ",".join(str(c) for c in cids)
    props = "MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES"
    url = f"{PUG}/compound/cid/{ids}/property/{props}/JSON"
    data = _get(url)
    return data.get("PropertyTable", {}).get("Properties", [])


def get_sdf_3d(cid):
    """
    Confirmed real endpoint pattern from PubChem/IUPAC documentation:
        {PUG}/compound/cid/{cid}/record/SDF?record_type=3d
    Returns the raw SDF (structure-data file) text, containing real atom
    positions and bond orders for a computationally-generated 3D conformer
    (PubChem is explicit that this is computed, not experimentally
    measured — worth saying plainly to a student, not just to me).

    NOTE: not every PubChem compound has a precomputed 3D conformer —
    this is a real, documented PubChem behaviour (confirmed against a
    real user report of the exact same symptom before assuming it was a
    bug in this client). Simple ionic salts are disproportionately likely
    to be missing one, since 3D conformer generation applies most reliably
    to single connected covalent structures. Raises PubChemNotFound if
    there's genuinely no 3D record — callers should fall back to
    get_sdf_2d() + local conformer generation, which is exactly what
    rdkit_engine.geometry_from_cid() below does.
    """
    url = f"{PUG}/compound/cid/{cid}/record/SDF?record_type=3d"
    return _get(url, parse_json=False)


def get_sdf_2d(cid):
    """
    The default (no record_type parameter) SDF endpoint. Far more
    universally available than the 3D one — 2D layout is computed for
    essentially every compound in PubChem, including simple salts that
    lack a stored 3D conformer.
    """
    url = f"{PUG}/compound/cid/{cid}/record/SDF"
    return _get(url, parse_json=False)


def get_sdf_best_effort(cid):
    """
    Try the real PubChem 3D conformer first; if PubChem doesn't have one
    for this compound (a real, common case, not an error), fall back to
    its 2D structure. Returns (sdf_text, is_3d) so the caller can be
    honest with the user about which one they got — and, if it's 2D,
    that a 3D conformer needs to be generated locally rather than
    silently presenting flattened 2D coordinates as if they were real 3D.
    """
    try:
        return get_sdf_3d(cid), True
    except PubChemNotFound:
        return get_sdf_2d(cid), False


def lookup_formula_with_names(formula, max_results=6):
    """
    Convenience wrapper: formula -> [(cid, name, molecular_formula), ...]
    so the app can show real candidate names for disambiguation, exactly
    because a bare formula is often genuinely ambiguous (see module note).
    """
    cids = search_by_formula(formula, max_results=max_results)
    props = get_properties(cids)
    out = []
    for p in props:
        name = p.get("IUPACName") or p.get("MolecularFormula", "?")
        out.append({
            "cid": p["CID"],
            "name": name,
            "formula": p.get("MolecularFormula", "?"),
            "mass": p.get("MolecularWeight"),
            "smiles": p.get("CanonicalSMILES"),
        })
    return out
