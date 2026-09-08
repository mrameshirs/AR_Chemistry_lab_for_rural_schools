# 🧪 BondVision AR — AI-Powered Chemistry Learning Platform

Tap atoms into a workspace and watch real bond formation happen automatically —
ionic electron transfer, covalent electron sharing, real VSEPR 3D shapes, and a
transparent green-chemistry advisor. Built for classrooms with no lab equipment.

INSPIRE-MANAK project prototype · Nithyamithran Ramesh, Class VI B, Manav Mandir High School

## What it does

**🧪 Atom Workspace** — the main feature. Tap elements from a 22-atom palette
(or use a one-click preset like Water, Table Salt, Glucose, or Benzene) and the
app automatically figures out:
- whether the result is ionic or covalent, and why
- the real 3D VSEPR shape (for simple molecules)
- an electron-level diagram: valence electrons transferring between ions, or
  shared electron pairs and lone pairs for covalent bonds
- the correct chemical name, including picking the right oxidation state for
  metals like iron or copper that have more than one (Fe²⁺ vs Fe³⁺, Cu⁺ vs Cu²⁺)

Tapping 6 carbon + 6 hydrogen atoms correctly identifies **benzene**, and
6 carbon + 12 hydrogen + 6 oxygen correctly identifies **glucose** — these are
genuinely different molecules, and the app tells them apart correctly rather
than guessing from atom count alone.

When a combination is genuinely beyond the local engine — tap Mg + S + 4 O and
you get magnesium sulfate, a metal ion bound to a polyatomic sulfate ion, not
a simple binary compound — the app offers a real online lookup via PubChem and
RDKit instead of refusing outright. See "Real chemistry lookup" below.

**📚 Molecule Library** — 31 curated molecules with real, derived VSEPR
geometry and either idealised or exact textbook bond angles.

**🌱 Green Chemistry Compare** — a transparent, rule-based environmental
score. The grade is driven by the *weakest* of persistence, toxicity, and
climate impact, not an average — so a molecule that's chemically inert but a
severe greenhouse gas (like SF₆) can't hide behind a good toxicity score.

**🥽 AR Preview** — two modes, and two molecule sources. **Photo Overlay**
takes one photo and composites a rendered molecule onto it, with sliders for
position and size — tested, reliable, works everywhere. **Live Marker AR** is
genuine real-time AR: a live camera feed with a 3D molecule locked to a
printed marker card as you move the phone, built on real A-Frame + AR.js, now
with an optional floating text label rendered directly in the 3D scene. Either
mode works on a Molecule Library pick or on whatever you most recently looked
up via PubChem — real molecules the local engine can't build on its own now
show up in AR too. See "AR: what's actually verified" below before relying on
it for a live demo.

**🎮 Quiz & Badges**, **⚠️ Lab Safety Advisor**, **🌐 Multilingual Explainer**
(English + Hindi in full sentences, four more languages at vocabulary level),
**👩‍🏫 Teacher Dashboard** (session-scoped quiz analytics).

## AR: what's actually verified

The Live Marker AR mode is real, working code — not a mockup — built on two
genuine, currently-published npm packages (`aframe@1.8.0` and
`@ar-js-org/ar.js@3.4.8`). Here's exactly what was checked, and by what
method, before it shipped:

| Claim | How it was verified |
|---|---|
| The npm packages are real and current | Downloaded both tarballs directly and inspected their actual contents, not just their names |
| The exact build files used exist | Confirmed `aframe-ar.js` (1.68MB) matches the path in the package's own `package.json` `"main"` field |
| The HTML/component syntax is correct | Taken directly from AR.js's own bundled README, not recalled from memory |
| The printed marker will actually track | The real, official Hiro marker PNG was downloaded from the AR.js project itself and bundled — not redrawn by hand, which would not track |
| The browser will grant camera access | Streamlit's own compiled frontend source was read directly, confirming its iframe hardcodes `allow="camera; microphone; ..."` |
| The 3D bond-orientation math is correct | Checked against known cases (identity, 180° flip, axis alignment) before being used to generate any HTML |

**What was not, and could not be, verified from this environment:** whether
a real phone camera, under real classroom lighting, with a real paper
printout, actually detects the marker reliably, and whether frame rate is
acceptable on genuinely low-end Android hardware. There is no way to test a
live camera and a physical printed marker from a sandboxed backend with no
browser. **Test this yourself, with a real printed marker and a real phone,
before a live demo** — the Photo Overlay mode is the tested fallback if it
doesn't work well on the day.

## Why tap, not drag

Real mouse drag-and-drop needs a custom browser component with its own build
pipeline — and native HTML5 drag-and-drop is notoriously unreliable on
touchscreens, which matters directly here since this project explicitly
targets budget rural-school Android phones. Tapping works identically on a
mouse click or a thumb tap, so that's the interaction the workspace uses.
Real drag-and-drop is listed as a roadmap item, not silently skipped.

## Chemistry scope, stated honestly

The Atom Workspace's *local* engine models **simple AXₙ molecules**: one
central atom, n identical peripheral atoms, single bonds unless a curated
entry specifies otherwise (like CO₂'s double bonds). This is genuinely
correct and useful — it's the same model every school VSEPR chapter uses —
but it cannot derive real structures for polyatomic ions (sulfate,
carbonate) or molecules with carbon chains or rings.

If you tap a combination with no verified real structure in the local
engine (a lone atom, or an unbalanced ratio like 1 iron + 1 chlorine, which
doesn't correspond to any real iron chloride), the app says so honestly —
and then offers a real second path.

## Real chemistry lookup: PubChem + RDKit

Tap Mg + S + 4 O and the local engine correctly declines — magnesium
sulfate isn't a simple binary ionic compound, it's a metal ion attracted
to a polyatomic sulfate ion, which has its own internal covalent bonding.
Rather than guess, the Atom Workspace can look this up for real:

```
PubChem (real chemistry database, 100M+ compounds)
    -> resolves what actually forms from these elements
    -> returns a real, computed 3D structure

RDKit (real cheminformatics library)
    -> interprets that structure: bond orders, formal charges, lone pairs
    -> the SAME per-atom electron-counting formula this app already uses
       for simple molecules, applied atom-by-atom instead of assuming one
       central atom — verified against sulfate's real textbook Lewis
       structure before shipping (see chem_engine.py / rdkit_engine.py)

This app's existing renderers
    -> the same Plotly 3D viewer, Lewis-dot diagram, and AR scene builder
       already built for the local engine, unmodified — they only care
       about atoms/bonds/positions, not where the data came from
```

This is deliberately **not built on an LLM or generative AI**. Bond order,
formal charge, and lone-pair count for a known structure have one correct
answer, computable with established graph algorithms — the wrong tool for
that job is a model that predicts plausible text. An LLM's legitimate role
here would be a completely different, open-ended feature — an actual chat
tutor answering unpredictable student questions — not deciding how many
electrons sit on an oxygen atom.

**A real bug this surfaced, worth documenting rather than quietly fixing:**
formula search on PubChem is asynchronous — a search across their ~110
million compounds doesn't return an answer immediately, it returns a job
key that has to be polled. The first version of this client only handled
the immediate-response shape, so a real, extremely common compound
(forsterite, Mg₂SiO₄ — one of the most abundant minerals in Earth's mantle)
came back as a false "no record found," because the job key was silently
discarded instead of followed up on. Fixed by polling
`compound/listkey/{key}/cids/JSON` until the job resolves, verified against
PubChem's own documented async example (their similarity-search tutorial)
before shipping, and tested against the exact real failure case.

**A second real limitation, also worth documenting:** not every PubChem
compound has a precomputed 3D conformer — this is genuine, confirmed
PubChem behaviour, not a bug in this client (a PubChemPy user reported the
identical symptom independently: a compound's basic record exists, but its
3D-specific record doesn't). Simple ionic salts are disproportionately
affected, since 3D conformer generation applies most reliably to single
connected covalent structures, not multi-ion salts. Rather than error out,
the app now falls back to PubChem's far more universal 2D structure and
generates a real 3D conformer locally with RDKit's standard ETKDG
algorithm — and says so plainly in the app rather than presenting a
locally-guessed structure as if it were PubChem's own data.

**What this can't reliably do:** disambiguate between real distinct
compounds sharing a common name (searching "magnesium sulfate" actually
returns at least 8 different real PubChem records — anhydrous, mono- through
nona-hydrate — because hydration state is genuinely a different compound;
the app shows the real candidates and makes you pick, rather than guessing
one), and transition-metal/organometallic bonding, where simple main-group
valence counting doesn't reliably apply and lone-pair counts are flagged as
approximate rather than stated with false confidence.

**The one real cost:** this needs internet. It's the single feature in
BondVision AR that isn't ₹0-and-offline, and it's used only on explicit
request (a button click), never automatically.

## What's roadmap, not shipped

- **Real drag-and-drop** — needs a custom Streamlit component; tapping was
  chosen deliberately instead (see above)
- **Verified environmental/toxicity data for arbitrary PubChem lookups** —
  the green-chemistry score for anything outside the 31+4 hand-curated
  molecules is an honest structural estimate, not a verified rating;
  covering more compounds properly needs either a lot more careful manual
  curation or wiring up PubChem's separate GHS/safety-data endpoints
- **Multilingual LLM tutor (8 languages)** — needs an LLM API key, not
  configured in this deployment
- **Photo-to-solution homework helper** — needs a vision-capable LLM

## Running it locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

No API keys or secrets are required. PubChem's API is free and needs no
key — the only external dependency is that the PubChem lookup feature in
the Atom Workspace needs an internet connection at the moment it's used;
everything else runs entirely from built-in chemistry data and logic.

## Deploying on Streamlit Community Cloud

1. Push this repo to GitHub
2. [share.streamlit.io](https://share.streamlit.io) → New app → point at this
   repo → main file `app.py` → Deploy

That's it — there's no `secrets.toml` to configure for this project.

## Files

```
app.py             the Streamlit app: navigation + all 9 pages
chem_data.py       periodic data, curated molecules, green-chemistry facts,
                   safety combinations, translation vocabulary — pure data
chem_engine.py     VSEPR geometry math, bond classification, the atom-tap
                   matcher, ionic charge balancing, green-chemistry scoring,
                   Hill-notation formula builder
viz.py             Plotly 3D viewer, Matplotlib Lewis-dot diagrams, AR photo
                   compositing, printable flashcard generator
live_ar.py         real marker-tracked AR scene builder (A-Frame + AR.js),
                   bond-orientation math, in-scene text labels, the real
                   Hiro marker as base64
pubchem_client.py  real PubChem PUG-REST API client — formula/name search,
                   3D structure retrieval, rate limiting, graceful failure
rdkit_engine.py    turns a real structure (from PubChem or anywhere else)
                   into the same MoleculeGeometry format the local engine
                   produces — bond orders, formal charges, lone pairs
assets/            the real, official AR.js Hiro marker image
requirements.txt   Python dependencies
.streamlit/        colour theme
```

## Honest limits

- VSEPR shapes and bond angles reflect a standard general-chemistry teaching
  model, not a quantum-mechanical calculation
- The green-chemistry grade is a documented rubric built from well-known
  public facts (IPCC greenhouse gas science, Montreal Protocol substances,
  standard hazard classes) — not a certified regulatory rating and not the
  output of a trained model
- Multilingual support is full-quality for English and Hindi; the other four
  languages are vocabulary-level and marked as a starting point, not a
  teacher-reviewed translation

## Licence

MIT
