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

**📚 Molecule Library** — 31 curated molecules with real, derived VSEPR
geometry and either idealised or exact textbook bond angles.

**🌱 Green Chemistry Compare** — a transparent, rule-based environmental
score. The grade is driven by the *weakest* of persistence, toxicity, and
climate impact, not an average — so a molecule that's chemically inert but a
severe greenhouse gas (like SF₆) can't hide behind a good toxicity score.

**🥽 AR Preview** — takes one photo and composites a rendered molecule onto
it, with sliders for position and size.

**🎮 Quiz & Badges**, **⚠️ Lab Safety Advisor**, **🌐 Multilingual Explainer**
(English + Hindi in full sentences, four more languages at vocabulary level),
**👩‍🏫 Teacher Dashboard** (session-scoped quiz analytics).

## Why tap, not drag

Real mouse drag-and-drop needs a custom browser component with its own build
pipeline — and native HTML5 drag-and-drop is notoriously unreliable on
touchscreens, which matters directly here since this project explicitly
targets budget rural-school Android phones. Tapping works identically on a
mouse click or a thumb tap, so that's the interaction the workspace uses.
Real drag-and-drop is listed as a roadmap item, not silently skipped.

## Chemistry scope, stated honestly

The Atom Workspace and Molecule Builder engine model **simple AXₙ
molecules**: one central atom, n identical peripheral atoms, single bonds
unless a curated entry specifies otherwise (like CO₂'s double bonds). This is
a genuinely correct and useful model — it's the same one every school VSEPR
chapter uses — but it cannot derive real structures for molecules with
carbon chains or rings (ethanol, acetic acid, benzene, glucose). Those four
appear as **fact cards**: real formula, real molar mass, real environmental
facts, with an explicit note about why no 3D structure is shown rather than
a guessed one that would look plausible but be wrong.

If you tap a combination with no verified real structure in this tool (a
lone atom, or an unbalanced ratio like 1 iron + 1 chlorine, which doesn't
correspond to any real iron chloride), the app says so honestly rather than
fabricating an answer.

## What's roadmap, not shipped

- **Real drag-and-drop** — needs a custom Streamlit component; tapping was
  chosen deliberately instead (see above)
- **Live ArUco marker-tracked AR** — needs a native camera app or
  `streamlit-webrtc` with a reliable TURN server; a browser tab can only grab
  one still photo, not a live frame-by-frame video feed
- **A trained graph neural network on EPA/PubChem toxicity data** — needs a
  real labelled dataset and training pipeline; the green-chemistry score here
  is a transparent, documented rule-based rubric instead
- **Multilingual LLM tutor (8 languages)** — needs an LLM API key, not
  configured in this deployment
- **Photo-to-solution homework helper** — needs a vision-capable LLM

## Running it locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

No API keys or secrets are required — everything here runs from built-in
chemistry data and logic.

## Deploying on Streamlit Community Cloud

1. Push this repo to GitHub
2. [share.streamlit.io](https://share.streamlit.io) → New app → point at this
   repo → main file `app.py` → Deploy

That's it — there's no `secrets.toml` to configure for this project.

## Files

```
app.py           the Streamlit app: navigation + all 9 pages
chem_data.py     periodic data, curated molecules, green-chemistry facts,
                 safety combinations, translation vocabulary — pure data
chem_engine.py   VSEPR geometry math, bond classification, the atom-tap
                 matcher, ionic charge balancing, green-chemistry scoring
viz.py           Plotly 3D viewer, Matplotlib Lewis-dot diagrams, AR photo
                 compositing, printable flashcard generator
requirements.txt Python dependencies
.streamlit/      colour theme
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
