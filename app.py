"""
BondVision AR — AI-powered chemistry learning platform for rural schools
=========================================================================
INSPIRE-MANAK project prototype.
Nithyamithran Ramesh, Class VI B, Manav Mandir High School.

Honesty note for anyone reading this file: the original project brief also
described live ArUco-marker AR, a GNN trained on EPA/PubChem data, and an
LLM explaining things in 8 languages. None of those are shipped here — they
need a native camera pipeline, a labelled toxicology dataset, and an API key
respectively, none of which exist in this environment. What ships instead is
real, working, and clearly labelled: actual VSEPR chemistry, a transparent
rule-based green-chemistry rubric built from well-established environmental
facts, and a photo-overlay AR *preview*. See the Home page and README for
what's real today versus roadmap.
"""

import base64
import io
import random

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from chem_data import (
    ELEMENTS, CENTRAL_ALLOWED_N, PERIPHERAL_ELEMENTS, IONIC_METALS, IONIC_ANIONS,
    CURATED_MOLECULES, MOL_BY_FORMULA, GREEN_FACTS, SAFETY_COMBOS, VOCAB,
    EXACT_ANGLES, GEOMETRY_INFO, PALETTE_ORDER, COMPLEX_BY_FORMULA, ION_CHARGES,
)
from chem_engine import build_covalent, build_ionic, green_score, match_atoms, hill_formula
from viz import (
    plotly_molecule_figure, matplotlib_sticker, composite_ar, molecule_card,
    lewis_dot_covalent, lewis_dot_ionic,
)
from live_ar import molecule_to_ar_html, hiro_marker_base64
import pubchem_client
from rdkit_engine import geometry_from_sdf

st.set_page_config(page_title="BondVision AR", page_icon="🧪", layout="wide")

# =============================================================================
# STYLE
# =============================================================================
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Baloo+2:wght@600;700;800&family=Mukta:wght@400;600;700&display=swap');
html, body, [class*="css"] { font-family: 'Mukta', sans-serif; }
h1, h2, h3, h4 { font-family: 'Baloo 2', sans-serif !important; }

.hero {
    background: linear-gradient(120deg, #0E7C6B 0%, #11998E 45%, #2C7A4B 100%);
    padding: 26px 30px; border-radius: 22px; color: #fff; margin-bottom: 16px;
    box-shadow: 0 10px 28px rgba(14,124,107,.28);
}
.hero h1 { margin: 0; font-size: 2.3rem; font-weight: 800; }
.hero p { margin: 6px 0 0 0; font-weight: 600; opacity: .95; }
.tag { display:inline-block; background: rgba(255,255,255,.22); padding: 4px 13px;
       border-radius: 999px; font-size: .8rem; font-weight: 700; margin: 6px 6px 0 0; }

.tile { border-radius: 16px; padding: 14px 16px; color:#fff; height:100%; }
.tile b { display:block; font-family:'Baloo 2',sans-serif; font-size:1.7rem; }
.tile span { font-size:.82rem; font-weight:600; opacity:.95; }
.t-green {background: linear-gradient(135deg,#11998E,#38EF7D);}
.t-blue  {background: linear-gradient(135deg,#2193B0,#6DD5ED);}
.t-purple{background: linear-gradient(135deg,#8E2DE2,#4A00E0);}
.t-orange{background: linear-gradient(135deg,#F7971E,#FFD200); color:#3B2A00;}
.t-red   {background: linear-gradient(135deg,#EB3349,#F45C43);}

.grade-box { border-radius: 16px; padding: 18px 22px; color:#fff; margin: 8px 0 14px 0; }
.grade-box h3 { margin:0; font-size:1.5rem; }
.g-A2 {background: linear-gradient(120deg,#11998E,#38EF7D);}
.g-A  {background: linear-gradient(120deg,#2193B0,#6DD5ED);}
.g-B  {background: linear-gradient(120deg,#F2994A,#F2C94C); color:#3B2A00;}
.g-C  {background: linear-gradient(120deg,#EB3349,#F45C43);}

.safety-card {background:#FFF3E0; border:2px solid #F7971E; border-radius:14px; padding:14px 18px; margin-bottom:10px;}
.roadmap {background:#EFEFEF; border-left:4px solid #8E2DE2; border-radius:8px; padding:10px 16px; font-size:.88rem; color:#444;}
.small-note {font-size:.82rem; color:#6E7C76;}
</style>
""", unsafe_allow_html=True)

GRADE_CLASS = {"A+": "g-A2", "A": "g-A", "B": "g-B", "C": "g-C"}


def grade_css(grade_text):
    return GRADE_CLASS.get(grade_text.split(" ")[0], "g-B")


# session state defaults
if "quiz_score" not in st.session_state:
    st.session_state.quiz_score = 0
if "quiz_total" not in st.session_state:
    st.session_state.quiz_total = 0
if "quiz_log" not in st.session_state:
    st.session_state.quiz_log = []
if "current_geo" not in st.session_state:
    st.session_state.current_geo = None
    st.session_state.current_label = None


def set_current(geo, label):
    st.session_state.current_geo = geo
    st.session_state.current_label = label


# =============================================================================
# SIDEBAR NAV
# =============================================================================
PAGES = [
    "🏠 Home",
    "🧪 Atom Workspace",
    "📚 Molecule Library",
    "🌱 Green Chemistry Compare",
    "🥽 AR Preview",
    "🎮 Quiz & Badges",
    "⚠️ Lab Safety Advisor",
    "🌐 Multilingual Explainer",
    "👩‍🏫 Teacher Dashboard",
]
page = st.sidebar.radio("Navigate", PAGES, label_visibility="collapsed")
st.sidebar.markdown("---")
st.sidebar.caption(
    "Runs entirely in this browser tab — no camera image or score is uploaded "
    "anywhere except the single AR photo you choose to capture, which stays "
    "in your session only."
)

# =============================================================================
# HOME
# =============================================================================
if page == "🏠 Home":
    st.markdown("""
    <div class="hero">
      <h1>🧪 BondVision AR</h1>
      <p>3D chemistry, AR-style visualisation, and an AI-flavoured green chemistry advisor — for classrooms with no lab equipment.</p>
      <div>
        <span class="tag">Real VSEPR chemistry</span>
        <span class="tag">Transparent green score</span>
        <span class="tag">Works on any laptop</span>
        <span class="tag">₹0 to run</span>
        <span class="tag">INSPIRE-MANAK prototype</span>
      </div>
    </div>
    """, unsafe_allow_html=True)

    c1, c2, c3, c4 = st.columns(4)
    for col, cls, big, small in [
        (c1, "t-green", "31+4", "real molecules recognised"),
        (c2, "t-blue", "22", "elements on the tap palette"),
        (c3, "t-purple", "6", "safety warnings covered"),
        (c4, "t-orange", "0 ₹", "cost to run in a browser"),
    ]:
        col.markdown(f'<div class="tile {cls}"><b>{big}</b><span>{small}</span></div>', unsafe_allow_html=True)

    st.write("")
    st.markdown("#### What's real in this prototype — and what's roadmap")
    left, right = st.columns(2)
    with left:
        st.success(
            "**Working today**\n\n"
            "- Tap-to-build atom workspace with automatic bond formation\n"
            "- 31-molecule curated library + 4 complex-organic fact cards\n"
            "- **Real chemical lookup via PubChem + RDKit** — combinations the local "
            "engine can't build (polyatomic ions like sulfate, real molecules like "
            "aspirin) are looked up in a real 100M+ compound database and correctly "
            "interpreted for bonding and electron structure (needs internet)\n"
            "- Electron-level views: ionic transfer diagrams and Lewis dot structures, "
            "now including real formal-charge/lone-pair chemistry for looked-up ions\n"
            "- Transparent, rule-based green-chemistry scoring\n"
            "- Photo-overlay AR preview (tested, reliable fallback)\n"
            "- **Live marker-tracked AR** (real camera + real-time 3D, built on A-Frame/AR.js — "
            "now also works on real PubChem-sourced structures, with an optional in-scene "
            "floating label; shipped, but needs your own device test before a live demo)\n"
            "- Quiz with badges, session-scoped leaderboard\n"
            "- Lab safety advisor (real, well-known hazard pairs)\n"
            "- English + Hindi explanations, 4 more languages at vocabulary level"
        )
    with right:
        st.markdown(
            '<div class="roadmap"><b>Roadmap — needs infrastructure this prototype doesn\'t have</b><br><br>'
            "• Real drag-and-drop (mouse) — needs a custom browser component; "
            "tapping was chosen instead since it also works on touchscreens<br>"
            "• Verified environmental/toxicity data for PubChem-looked-up compounds "
            "beyond the 31+4 hand-curated ones — needs either a lot more careful manual "
            "curation or a second data source (e.g. PubChem's separate GHS/safety endpoints)<br>"
            "• Multilingual LLM tutor (8 languages) — needs an LLM API key; "
            "not configured in this deployment<br>"
            "• Photo-to-solution homework helper — needs a vision-capable LLM"
            "</div>",
            unsafe_allow_html=True,
        )

    st.write("")
    st.markdown("#### How a student would use it")
    st.markdown(
        "1. **Tap atoms** into the Atom Workspace, or try a famous combination with one click\n"
        "2. Watch the real bond type, electron transfer, and molecular geometry appear automatically\n"
        "3. Check its **Green Chemistry** grade and the actual environmental facts behind it\n"
        "4. Try the **AR Preview** to see it hover over a photo of your desk\n"
        "5. Take the **Quiz** and earn badges\n"
    )

# =============================================================================
# ATOM WORKSPACE
# =============================================================================
elif page == "🧪 Atom Workspace":
    st.markdown("## 🧪 Atom Workspace")
    st.caption(
        "Tap atoms into the workspace below. As soon as the counts match a real molecule, "
        "BondVision AR figures out the bonding automatically — ionic electron transfer, "
        "covalent electron sharing, the lot — and shows you exactly how it happened."
    )
    st.markdown(
        '<div class="roadmap">A note on the interaction: true drag-and-drop needs a custom '
        "browser component this tool doesn't include, and it's also notoriously unreliable on "
        "touchscreens — which matters, since this project targets budget rural-school phones. "
        "Tapping works identically on a mouse or a thumb, so that's what drives the workspace "
        "below.</div>", unsafe_allow_html=True,
    )

    if "basket" not in st.session_state:
        st.session_state.basket = {}
    if "basket_history" not in st.session_state:
        st.session_state.basket_history = []

    def add_atom(elem):
        st.session_state.basket[elem] = st.session_state.basket.get(elem, 0) + 1
        st.session_state.basket_history.append(elem)

    def undo_last():
        if st.session_state.basket_history:
            last = st.session_state.basket_history.pop()
            st.session_state.basket[last] -= 1
            if st.session_state.basket[last] <= 0:
                del st.session_state.basket[last]

    def clear_basket():
        st.session_state.basket = {}
        st.session_state.basket_history = []

    st.markdown("#### Try a famous combination")
    p1, p2, p3, p4, p5 = st.columns(5)
    if p1.button("💧 Water\n(2H + 1O)", width='stretch'):
        clear_basket(); add_atom("H"); add_atom("H"); add_atom("O")
    if p2.button("🧂 Table salt\n(1Na + 1Cl)", width='stretch'):
        clear_basket(); add_atom("Na"); add_atom("Cl")
    if p3.button("🫧 Carbon dioxide\n(1C + 2O)", width='stretch'):
        clear_basket(); add_atom("C"); add_atom("O"); add_atom("O")
    if p4.button("🍯 Glucose\n(6C+12H+6O)", width='stretch'):
        clear_basket()
        for _ in range(6): add_atom("C")
        for _ in range(12): add_atom("H")
        for _ in range(6): add_atom("O")
    if p5.button("⚪ Benzene\n(6C + 6H)", width='stretch'):
        clear_basket()
        for _ in range(6): add_atom("C")
        for _ in range(6): add_atom("H")

    st.markdown("#### Or tap atoms one at a time")
    cols = st.columns(7)
    for i, elem in enumerate(PALETTE_ORDER):
        info = ELEMENTS[elem]
        count = st.session_state.basket.get(elem, 0)
        label = f"{elem}{' · ' + str(count) if count else ''}"
        with cols[i % 7]:
            st.markdown(
                f'<div style="width:14px;height:14px;border-radius:4px;background:{info["color"]};'
                f'border:1px solid #333;margin:2px auto 0;"></div>', unsafe_allow_html=True,
            )
            if st.button(label, key=f"add_{elem}", width='stretch', help=info["name"]):
                add_atom(elem)

    bc1, bc2, bc3 = st.columns([2, 1, 1])
    basket_text = ", ".join(f"{e}×{n}" for e, n in st.session_state.basket.items()) or "empty"
    bc1.markdown(f"**In the workspace:** {basket_text}")
    if bc2.button("↩️ Undo last", width='stretch'):
        undo_last()
    if bc3.button("🗑️ Clear all", width='stretch'):
        clear_basket()

    st.divider()

    result = match_atoms(st.session_state.basket)

    if result["kind"] == "empty":
        st.info("Tap some atoms above, or try one of the famous combinations, to see what forms.")

    elif result["kind"] == "no_match":
        formula_str = hill_formula(st.session_state.basket)
        st.warning(
            f"**No verified real structure for this exact combination in the local engine** — "
            f"that's not the same as saying it can't exist, just that BondVision AR's own AXₙ "
            f"model won't guess at a structure it can't confirm. Total mass of what's in the "
            f"workspace: **{result['total_mass']} g/mol**."
        )
        st.caption(
            "For a *local* covalent molecule, one element usually needs to be at a count of "
            "exactly 1 (the central atom); for a *local* ionic compound, the ratio needs to "
            "exactly charge-balance. But real chemistry is bigger than that — try the real "
            "lookup below."
        )

        st.markdown("---")
        st.markdown(f"### 🔍 Look up **{formula_str}** on PubChem (needs internet)")
        st.caption(
            "PubChem is the U.S. National Institutes of Health's free public chemistry database — "
            "100+ million real compounds. This calls their live API, so it needs an internet "
            "connection and won't work if the classroom's connection is down; that's the one part "
            "of BondVision AR that needs internet."
        )

        if st.button("🌐 Search PubChem for this formula", type="primary"):
            with st.spinner(f"Searching PubChem for {formula_str}…"):
                try:
                    candidates = pubchem_client.lookup_formula_with_names(formula_str, max_results=6)
                    st.session_state.pubchem_candidates = candidates
                    st.session_state.pubchem_error = None
                except pubchem_client.PubChemNotFound:
                    st.session_state.pubchem_candidates = []
                    st.session_state.pubchem_error = (
                        f"PubChem has no record matching the exact formula {formula_str}. "
                        f"This can be genuinely correct — not every element combination forms "
                        f"a real, stable compound."
                    )
                except pubchem_client.PubChemError as e:
                    st.session_state.pubchem_candidates = []
                    st.session_state.pubchem_error = str(e)

        if st.session_state.get("pubchem_error"):
            st.error(st.session_state.pubchem_error)

        candidates = st.session_state.get("pubchem_candidates") or []
        if candidates:
            if len(candidates) > 1:
                st.info(
                    f"**{len(candidates)} real compounds share this exact formula** — this is "
                    f"genuinely common (different hydration states, isomers, or charge states "
                    f"are all distinct real compounds). Pick which one:"
                )
            labels = [f"{c['name']} (CID {c['cid']}, {c['mass']} g/mol)" for c in candidates]
            pick_idx = st.radio("Real PubChem matches", range(len(candidates)),
                                 format_func=lambda i: labels[i], label_visibility="collapsed")
            chosen = candidates[pick_idx]

            if st.button(f"📥 Load 3D structure for {chosen['name']}"):
                with st.spinner("Fetching real 3D structure from PubChem…"):
                    try:
                        sdf = pubchem_client.get_sdf_3d(chosen["cid"])
                        geo = geometry_from_sdf(sdf, formula=chosen["formula"], name=chosen["name"])
                        st.session_state.pubchem_geo = geo
                        st.session_state.pubchem_geo_label = chosen["name"]
                        set_current(geo, chosen["name"])
                    except Exception as e:
                        st.error(f"Could not load or parse the 3D structure: {e}")

        if st.session_state.get("pubchem_geo") is not None:
            geo = st.session_state.pubchem_geo
            label = st.session_state.pubchem_geo_label
            st.success(f"**{label}** — loaded from a real PubChem 3D structure, interpreted by RDKit")
            for note in geo.notes:
                st.info(note)

            tab3d, tab_electrons = st.tabs(["3D shape (from PubChem)", "⚡ Electrons & bonding (via RDKit)"])
            with tab3d:
                colL, colR = st.columns([3, 2])
                with colL:
                    st.plotly_chart(plotly_molecule_figure(geo), width='stretch')
                with colR:
                    st.markdown(f"### {geo.formula}")
                    st.metric("Net charge", f"{geo.meta['net_charge']:+d}" if geo.meta["net_charge"] else "0 (neutral)")
                    st.metric("Molar mass", f"{geo.molar_mass} g/mol")
                    st.caption("3D coordinates are PubChem's computed conformer, not an experimental "
                               "measurement — PubChem is explicit about this distinction too.")
            with tab_electrons:
                fig = lewis_dot_covalent(geo)
                st.pyplot(fig)
                plt.close(fig)
                st.caption("Bond orders and formal charges computed by RDKit from the real PubChem "
                           "structure — not this app's own simplified AXₙ model.")

            gscore = green_score(geo.formula)
            st.markdown(
                f'<div class="grade-box {grade_css(gscore["grade"])}"><h3>🌱 {gscore["grade"]}</h3>'
                f'<p>{gscore["note"]}</p></div>', unsafe_allow_html=True)
            if not gscore["verified"]:
                st.caption("⚠️ " + gscore["note"] + " (This compound isn't in BondVision AR's "
                           "hand-curated environmental facts table, so this is a structural estimate only.)")

            st.info("💡 Head to the **AR Preview** page and choose \"Use my last lookup\" to view "
                    "this real structure in AR.")

    elif result["kind"] == "complex":
        m = result["molecule"]
        mass = round(sum(ELEMENTS[e]["mass"] * n for e, n in m["atoms"].items()), 2)
        set_current(None, m["formula"])
        st.success(f"**Match found: {m['name']} ({m['formula']})**  ·  {m['name_hi']}")
        st.info(
            f"**Why there's no 3D shape shown:** {m['structure_note']} BondVision AR's builder "
            f"only handles a single central atom with identical surrounding atoms — genuinely "
            f"real, but genuinely limited. What follows is real reference data instead of a "
            f"guessed structure."
        )
        c1, c2 = st.columns(2)
        c1.metric("Molar mass", f"{mass} g/mol")
        gscore = green_score(m["formula"])
        c2.metric("Green chemistry grade", gscore["grade"].split(" — ")[0])
        st.markdown(
            f'<div class="grade-box {grade_css(gscore["grade"])}"><h3>🌱 {gscore["grade"]}</h3>'
            f'<p>{gscore["note"]}</p></div>', unsafe_allow_html=True)

    else:  # curated / generic_ionic / generic_covalent — a real buildable structure
        geo = result["geo"]
        set_current(geo, geo.formula)
        tag = {"curated": "📚 In the curated library", "generic_ionic": "🧮 Auto-derived (ionic)",
               "generic_covalent": "🧮 Auto-derived (covalent)"}[result["kind"]]
        st.success(f"**Match found: {geo.meta.get('chem_name') or geo.formula}**  ·  {tag}")

        tab3d, tab_electrons = st.tabs(["3D shape (VSEPR)", "⚡ Electrons & bond formation"])
        with tab3d:
            colL, colR = st.columns([3, 2])
            with colL:
                st.plotly_chart(plotly_molecule_figure(geo), width='stretch')
                for note in geo.notes:
                    st.info(note)
            with colR:
                st.markdown(f"### {geo.formula}")
                if not geo.is_ionic:
                    st.metric("Molecular geometry", geo.geometry_name)
                    if geo.bond_angle:
                        st.metric("Bond angle", f"{geo.bond_angle:.1f}°" +
                                  (" (idealised)" if geo.formula not in EXACT_ANGLES else " (real value)"))
                st.metric("Bond type", geo.bond_type)
                st.metric("Polarity", geo.molecule_polarity)
                st.metric("Molar mass", f"{geo.molar_mass} g/mol")

        with tab_electrons:
            if geo.is_ionic:
                meta = geo.meta
                fig = lewis_dot_ionic(meta["cation"], meta["anion"], meta["cation_charge"],
                                       meta["n_cation"], meta["n_anion"])
                st.pyplot(fig)
                plt.close(fig)
                st.caption(
                    f"{ELEMENTS[meta['cation']]['name']} gives up {meta['cation_charge']} "
                    f"electron(s) per atom; {ELEMENTS[meta['anion']]['name']} needs "
                    f"{abs(ION_CHARGES[meta['anion']][0])} to complete its outer shell. "
                    f"That's exactly why the formula is {geo.formula} and not some other ratio."
                )
            else:
                fig = lewis_dot_covalent(geo)
                st.pyplot(fig)
                plt.close(fig)
                st.caption(
                    "Each dot-pair on a bond line is one shared electron pair; the gold dot-pairs "
                    "are lone pairs not involved in any bond."
                )

        st.write("")
        gscore = green_score(geo.formula, has_halogen=any(
            a.element in ("F", "Cl", "Br", "I") for a in geo.atoms))
        st.markdown(
            f'<div class="grade-box {grade_css(gscore["grade"])}"><h3>🌱 {gscore["grade"]}</h3>'
            f'<p>{gscore["note"]}</p></div>', unsafe_allow_html=True)
        if not gscore["verified"]:
            st.caption("⚠️ " + gscore["note"])

# =============================================================================
# MOLECULE LIBRARY
# =============================================================================
elif page == "📚 Molecule Library":
    st.markdown("## 📚 Molecule Library")
    st.caption("31 molecules with genuine derived 3D shapes and bond angles — some idealised, some the exact textbook value. "
               "Complex organics like glucose and benzene appear in the Atom Workspace as fact cards instead.")

    formulas = [m["formula"] for m in CURATED_MOLECULES]
    pick = st.selectbox("Pick a molecule", formulas,
                         format_func=lambda f: f"{f} — {MOL_BY_FORMULA[f]['name']}")
    m = MOL_BY_FORMULA[pick]

    if m.get("is_ionic"):
        geo = build_ionic(m["cation"], m["anion"])
    else:
        kwargs = dict(central=m["central"], peripheral=m["peripheral"], n=m["n"],
                      bond_order=m["bond_order"], formula=m["formula"])
        if "steric_override" in m:
            kwargs["steric_override"] = m["steric_override"]
            kwargs["lone_pairs_override"] = m["lone_pairs_override"]
        if m["formula"] in EXACT_ANGLES:
            kwargs["exact_angle"] = EXACT_ANGLES[m["formula"]]
        geo = build_covalent(**kwargs)

    set_current(geo, pick)

    colL, colR = st.columns([3, 2])
    with colL:
        st.plotly_chart(plotly_molecule_figure(geo), width='stretch')
        if m.get("note_resonance"):
            st.info("Real bonding here is a resonance hybrid — the true structure is an average "
                    "of contributing forms, not a single fixed double/single bond pattern.")
        if m.get("extra_note"):
            st.info(m["extra_note"])
    with colR:
        st.markdown(f"### {m['name']}  ·  {m['name_hi']}")
        st.metric("Geometry", geo.geometry_name)
        if geo.bond_angle:
            exact = "real value" if pick in EXACT_ANGLES else "idealised"
            st.metric("Bond angle", f"{geo.bond_angle:.1f}° ({exact})")
        st.metric("Bond type", geo.bond_type)
        st.metric("Polarity", geo.molecule_polarity)
        st.metric("Molar mass", f"{geo.molar_mass} g/mol")

        card_buf = molecule_card(geo, m["name"], pick,
                                  green_score(pick)["grade"], lang_line=m["name_hi"])
        st.download_button("⬇️ Download printable flashcard", card_buf.getvalue(),
                            file_name=f"{pick}_card.png", mime="image/png",
                            width='stretch')

    st.write("")
    gscore = green_score(pick)
    st.markdown(
        f'<div class="grade-box {grade_css(gscore["grade"])}"><h3>🌱 {gscore["grade"]}</h3>'
        f'<p>{gscore["note"]}</p></div>', unsafe_allow_html=True)
    gc1, gc2, gc3 = st.columns(3)
    gc1.metric("Persistence", gscore["persistence"])
    gc2.metric("Toxicity", gscore["toxicity"])
    gc3.metric("Climate impact", gscore["climate"])

# =============================================================================
# GREEN CHEMISTRY COMPARE
# =============================================================================
elif page == "🌱 Green Chemistry Compare":
    st.markdown("## 🌱 Green Chemistry Compare")
    st.caption(
        "Pick two molecules and compare them side by side. The grade is driven by the "
        "**weakest** of persistence, toxicity, and climate impact — not an average — so a "
        "single serious hazard can't be hidden behind two good scores."
    )
    with st.expander("How is the score actually calculated?"):
        st.write(
            "Each of the three factors (persistence, toxicity, climate impact) is rated Low/Medium/"
            "High/Very High from well-established environmental chemistry — greenhouse gas science, "
            "Montreal Protocol ozone-depleting substances, and standard hazard classes. Each rating "
            "maps to a number, and the **overall score is the lowest of the three** — so a molecule "
            "that is chemically inert but a severe greenhouse gas (like SF₆) still gets flagged "
            "appropriately, instead of its good toxicity score covering for it."
        )
        st.caption(
            "This is a teaching rubric built from public, well-known facts — not a certified "
            "regulatory hazard rating (GHS/EPA/ECHA) and not the output of a trained model."
        )

    formulas = [m["formula"] for m in CURATED_MOLECULES]
    c1, c2 = st.columns(2)
    f1 = c1.selectbox("Molecule A", formulas, index=formulas.index("SF6"),
                       format_func=lambda f: f"{f} — {MOL_BY_FORMULA[f]['name']}")
    f2 = c2.selectbox("Molecule B", formulas, index=formulas.index("H2O"),
                       format_func=lambda f: f"{f} — {MOL_BY_FORMULA[f]['name']}")

    for col, f in [(c1, f1), (c2, f2)]:
        g = green_score(f)
        with col:
            st.markdown(f"#### {f} — {MOL_BY_FORMULA[f]['name']}")
            st.markdown(f'<div class="grade-box {grade_css(g["grade"])}"><h3>{g["grade"]}</h3>'
                        f'<p>{g["note"]}</p></div>', unsafe_allow_html=True)
            st.metric("Persistence", g["persistence"])
            st.metric("Toxicity", g["toxicity"])
            st.metric("Climate impact", g["climate"])

# =============================================================================
# AR PREVIEW
# =============================================================================
elif page == "🥽 AR Preview":
    st.markdown("## 🥽 AR Preview")

    source = st.radio(
        "Molecule source",
        ["📚 Molecule Library (31 curated shapes)",
         "🔍 My last PubChem lookup" + ("" if st.session_state.get("pubchem_geo") else " (nothing looked up yet)")],
        horizontal=False,
    )

    geo = None
    label = None
    if source.startswith("📚"):
        formulas = [m["formula"] for m in CURATED_MOLECULES]
        pick = st.selectbox("Molecule to preview", formulas,
                             format_func=lambda f: f"{f} — {MOL_BY_FORMULA[f]['name']}")
        m = MOL_BY_FORMULA[pick]
        if m.get("is_ionic"):
            st.info("Ionic lattices don't have a single discrete 3D shape to preview — pick a covalent molecule instead.")
        else:
            kwargs = dict(central=m["central"], peripheral=m["peripheral"], n=m["n"],
                          bond_order=m["bond_order"], formula=m["formula"])
            if "steric_override" in m:
                kwargs["steric_override"] = m["steric_override"]; kwargs["lone_pairs_override"] = m["lone_pairs_override"]
            if m["formula"] in EXACT_ANGLES:
                kwargs["exact_angle"] = EXACT_ANGLES[m["formula"]]
            geo = build_covalent(**kwargs)
            label = pick
    else:
        if st.session_state.get("pubchem_geo") is None:
            st.info("Nothing looked up yet — go to the **Atom Workspace**, tap a combination with no "
                    "local match, and use the PubChem lookup there first.")
        else:
            geo = st.session_state.pubchem_geo
            label = st.session_state.pubchem_geo_label
            st.success(f"Using **{label}** ({geo.formula}) — a real structure from PubChem, interpreted by RDKit. "
                       f"This is new: AR now works for real molecules the local engine can't build on its own, "
                       f"like polyatomic ions.")

    if geo is not None:
        ar_mode = st.radio(
            "AR mode",
            ["📸 Photo Overlay (tested, works everywhere)",
             "🎥 Live Marker AR (real-time, experimental — test on your own device first)"],
            horizontal=False,
        )

        if ar_mode.startswith("📸"):
            st.warning(
                "**What this does, in plain terms:** takes one still photo and pastes a rendered "
                "picture of the molecule on top of it, at a position you control with sliders. It "
                "does not track anything or update live — this is the reliable fallback, tested "
                "and known to work in every browser."
            )
            photo = st.camera_input("Point the camera at your desk and take a photo")
            c1, c2, c3 = st.columns(3)
            scale = c1.slider("Size", 0.15, 0.8, 0.4)
            x_frac = c2.slider("Left ↔ Right", 0.0, 1.0, 0.5)
            y_frac = c3.slider("Up ↔ Down", 0.0, 1.0, 0.5)

            if photo is not None:
                sticker = matplotlib_sticker(geo)
                result = composite_ar(photo.getvalue(), sticker, scale=scale, x_frac=x_frac, y_frac=y_frac)
                st.image(result, caption=f"{label}", width='stretch')
                buf = io.BytesIO(); result.save(buf, format="PNG")
                st.download_button("⬇️ Download this image", buf.getvalue(),
                                    file_name=f"{label}_ar_preview.png", mime="image/png")
            else:
                st.info("Take a photo above to see the molecule composited onto it.")

        else:
            st.error(
                "**This is real, live, marker-tracked AR** — a genuine camera feed with the "
                "molecule locked to a printed marker in real time, built on A-Frame + AR.js. "
                "It has been checked as far as possible without an actual browser: the real "
                "npm packages were downloaded and inspected, Streamlit's own compiled source "
                "was read to confirm its iframe grants camera permission, and the geometry math "
                "was verified against known cases. **What hasn't been tested is a real phone "
                "camera pointed at a real printed marker under real lighting — please try it "
                "yourself before showing it to a judge.**"
            )
            st.markdown("**Step 1 — print the marker** (this is the real, official AR.js Hiro marker):")
            marker_bytes = base64.b64decode(hiro_marker_base64())
            dl1, dl2 = st.columns([1, 3])
            dl1.download_button("⬇️ Download marker to print", marker_bytes,
                                 file_name="hiro_marker.png", mime="image/png")
            dl2.caption("Print at a decent size (at least 5-6 cm across) on plain paper, on a flat surface, "
                        "in good even lighting.")

            ar_scale = st.slider("Molecule size on the marker", 0.05, 0.6, 0.22, step=0.01)
            show_label = st.checkbox("Show a floating label in the AR scene itself (new)", value=True)
            label_text = None
            if show_label:
                gscore = green_score(geo.formula)
                label_text = f"{geo.formula} - {gscore['grade'].split(' — ')[0]}"

            st.markdown("**Step 2 — allow camera access below, then point it at the printed marker:**")
            ar_html = molecule_to_ar_html(geo, scale=ar_scale, label_text=label_text)
            try:
                st.iframe(ar_html, height=480)
            except AttributeError:
                components.html(ar_html, height=480, scrolling=False)

# =============================================================================
# QUIZ & BADGES
# =============================================================================
elif page == "🎮 Quiz & Badges":
    st.markdown("## 🎮 Quiz & Badges")
    st.caption("Scores reset when this browser tab's session ends — there's no shared server-side leaderboard yet.")

    def make_question():
        m = random.choice(CURATED_MOLECULES)
        if m.get("is_ionic"):
            geo = build_ionic(m["cation"], m["anion"])
        else:
            kwargs = dict(central=m["central"], peripheral=m["peripheral"], n=m["n"], bond_order=m["bond_order"], formula=m["formula"])
            if "steric_override" in m:
                kwargs["steric_override"] = m["steric_override"]; kwargs["lone_pairs_override"] = m["lone_pairs_override"]
            if m["formula"] in EXACT_ANGLES:
                kwargs["exact_angle"] = EXACT_ANGLES[m["formula"]]
            geo = build_covalent(**kwargs)

        kind = random.choice(["geometry", "polarity", "green"])
        if kind == "geometry":
            correct = geo.geometry_name
            pool = {v["name"] for v in GEOMETRY_INFO.values()} | {"Diatomic (linear)", "Ionic lattice (schematic slice)"}
            wrong = random.sample(list(pool - {correct}), k=min(3, len(pool) - 1))
            options = wrong + [correct]
            random.shuffle(options)
            q = f"What is the molecular geometry of **{m['formula']}** ({m['name']})?"
        elif kind == "polarity":
            correct = geo.molecule_polarity
            options = list({correct, "Polar", "Nonpolar", "Ionic — not a discrete polar/nonpolar molecule"})
            random.shuffle(options)
            q = f"Is **{m['formula']}** ({m['name']}) polar, nonpolar, or ionic overall?"
        else:
            g = green_score(m["formula"])
            correct = g["grade"]
            options = ["A+ — environmentally favourable", "A — generally low concern",
                       "B — handle with caution", "C — significant environmental or safety concern"]
            random.shuffle(options)
            q = f"Based on its real environmental profile, what green-chemistry grade does **{m['formula']}** ({m['name']}) get?"
        return q, options, correct, m["formula"]

    if "quiz_q" not in st.session_state:
        st.session_state.quiz_q = make_question()

    q, options, correct, formula = st.session_state.quiz_q
    st.markdown(f"#### {q}")
    choice = st.radio("Choose one", options, key=f"quiz_choice_{st.session_state.quiz_total}")

    c1, c2 = st.columns([1, 3])
    if c1.button("Submit", type="primary"):
        st.session_state.quiz_total += 1
        correct_now = (choice == correct)
        if correct_now:
            st.session_state.quiz_score += 1
            st.success(f"Correct! {formula} → {correct}")
        else:
            st.error(f"Not quite. {formula} → **{correct}**")
        st.session_state.quiz_log.append({"formula": formula, "correct": correct_now})
        st.session_state.quiz_q = make_question()
        st.rerun()

    st.write("")
    acc = (st.session_state.quiz_score / st.session_state.quiz_total * 100) if st.session_state.quiz_total else 0
    m1, m2, m3 = st.columns(3)
    m1.metric("Score", f"{st.session_state.quiz_score} / {st.session_state.quiz_total}")
    m2.metric("Accuracy", f"{acc:.0f}%")
    badge = ("🏆 Bond Master" if st.session_state.quiz_score >= 10 else
             "🥈 Green Chemist" if st.session_state.quiz_score >= 5 else
             "🥉 Getting started" if st.session_state.quiz_score >= 1 else "—")
    m3.metric("Badge", badge)

# =============================================================================
# LAB SAFETY ADVISOR
# =============================================================================
elif page == "⚠️ Lab Safety Advisor":
    st.markdown("## ⚠️ Virtual Lab Safety Advisor")
    st.caption(
        "These are well-known hazard warnings — the same kind printed on household bleach "
        "bottles and taught in every school chemistry lab. This page identifies hazards; "
        "it does not give procedures for making anything dangerous."
    )
    for combo in SAFETY_COMBOS:
        st.markdown(
            f'<div class="safety-card"><b>🚫 {combo["a"]} + {combo["b"]}</b><br>'
            f'<b>Danger:</b> {combo["danger"]}<br>'
            f'<span class="small-note">{combo["why"]}</span></div>',
            unsafe_allow_html=True,
        )
    st.markdown("#### General lab safety checklist")
    st.markdown(
        "- Always add **acid to water**, never water to acid\n"
        "- Wear goggles and gloves for anything corrosive\n"
        "- Never mix cleaning products without checking the label first\n"
        "- Work in a ventilated space when a reaction produces gas\n"
        "- Know where the eyewash station and fire extinguisher are before you start"
    )

# =============================================================================
# MULTILINGUAL EXPLAINER
# =============================================================================
elif page == "🌐 Multilingual Explainer":
    st.markdown("## 🌐 Multilingual Explainer")
    st.caption(
        "English and Hindi use full, natural sentence templates. Marathi, Tamil, Telugu and "
        "Bengali currently give a reviewed **vocabulary table** rather than a composed sentence — "
        "we'd rather ship an honest word-level translation than an untested full sentence."
    )

    formulas = [m["formula"] for m in CURATED_MOLECULES]
    pick = st.selectbox("Molecule", formulas, format_func=lambda f: f"{f} — {MOL_BY_FORMULA[f]['name']}")
    m = MOL_BY_FORMULA[pick]

    if m.get("is_ionic"):
        geo = build_ionic(m["cation"], m["anion"])
    else:
        kwargs = dict(central=m["central"], peripheral=m["peripheral"], n=m["n"], bond_order=m["bond_order"], formula=m["formula"])
        if "steric_override" in m:
            kwargs["steric_override"] = m["steric_override"]; kwargs["lone_pairs_override"] = m["lone_pairs_override"]
        if m["formula"] in EXACT_ANGLES:
            kwargs["exact_angle"] = EXACT_ANGLES[m["formula"]]
        geo = build_covalent(**kwargs)

    lang = st.selectbox("Language", ["English", "हिंदी Hindi", "मराठी Marathi", "தமிழ் Tamil",
                                      "తెలుగు Telugu", "বাংলা Bengali"])

    if lang == "English":
        polarity_word = "ionic" if geo.is_ionic else geo.molecule_polarity.lower()
        st.markdown(f"**{m['name']} ({pick})** has a **{geo.geometry_name.lower()}** shape "
                    f"and is **{polarity_word}**.")
    elif lang == "हिंदी Hindi":
        pol = VOCAB["hi"]["ionic"] if geo.is_ionic else (
            VOCAB["hi"]["polar"] if "Polar" in geo.molecule_polarity else VOCAB["hi"]["nonpolar"])
        st.markdown(f"**{m['name_hi']} ({pick})** का आकार **{geo.geometry_name_hi}** है "
                    f"और यह **{pol}** है।")
        if st.button("🔊 Play in Hindi"):
            try:
                from gtts import gTTS
                text = f"{m['name_hi']} का आकार {geo.geometry_name_hi} है और यह {pol} है।"
                buf = io.BytesIO()
                gTTS(text=text, lang="hi").write_to_fp(buf)
                st.audio(buf.getvalue(), format="audio/mp3")
            except Exception as e:
                st.warning(f"Voice playback needs internet access to Google's TTS service and isn't "
                           f"available right now ({e}). The written explanation above is still correct.")
    else:
        code = {"मराठी Marathi": "mr", "தமிழ் Tamil": "ta", "తెలుగు Telugu": "te", "বাংলা Bengali": "bn"}[lang]
        v = VOCAB[code]
        st.info(v["note"])
        st.markdown(f"**{pick}** — {geo.geometry_name}")
        cols = st.columns(3)
        cols[0].metric("Polar", v["polar"])
        cols[1].metric("Nonpolar", v["nonpolar"])
        cols[2].metric("Ionic", v["ionic"])

# =============================================================================
# TEACHER DASHBOARD
# =============================================================================
elif page == "👩‍🏫 Teacher Dashboard":
    st.markdown("## 👩‍🏫 Teacher Dashboard")
    st.caption(
        "Session-only analytics — this resets when the browser tab closes. A real multi-student "
        "dashboard would need a small database, which is listed as roadmap."
    )

    if not st.session_state.quiz_log:
        st.info("No quiz attempts yet in this session. Have a student try the Quiz page first.")
    else:
        df = pd.DataFrame(st.session_state.quiz_log)
        c1, c2, c3 = st.columns(3)
        c1.metric("Attempts", len(df))
        c2.metric("Correct", int(df["correct"].sum()))
        c3.metric("Accuracy", f"{df['correct'].mean() * 100:.0f}%")

        st.markdown("#### Most-missed molecules")
        missed = df[~df["correct"]]["formula"].value_counts()
        if len(missed):
            st.bar_chart(missed)
        else:
            st.success("No misses yet — a clean sheet!")

        st.markdown("#### Full attempt log")
        st.dataframe(df, width='stretch', hide_index=True)
        st.download_button("⬇️ Export as CSV", df.to_csv(index=False).encode(),
                            file_name="bondvision_quiz_log.csv", mime="text/csv")

# =============================================================================
# FOOTER
# =============================================================================
st.divider()
f1, f2 = st.columns([2, 3])
with f1:
    st.markdown(
        "**Nithyamithran Ramesh**  \nClass VI B · Manav Mandir High School  \n"
        "INSPIRE-MANAK project prototype"
    )
with f2:
    st.caption(
        "BondVision AR is a teaching tool, not a certified safety, toxicology, or regulatory "
        "reference. VSEPR shapes and bond angles reflect standard general-chemistry models; "
        "green-chemistry grades are a transparent rule-based rubric built from well-known "
        "public facts, not a trained scientific model or an official hazard rating."
    )
