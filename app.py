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

import io
import random

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

from chem_data import (
    ELEMENTS, CENTRAL_ALLOWED_N, PERIPHERAL_ELEMENTS, IONIC_METALS, IONIC_ANIONS,
    CURATED_MOLECULES, MOL_BY_FORMULA, GREEN_FACTS, SAFETY_COMBOS, VOCAB,
    EXACT_ANGLES, GEOMETRY_INFO, PALETTE_ORDER, COMPLEX_BY_FORMULA, ION_CHARGES,
)
from chem_engine import build_covalent, build_ionic, green_score, match_atoms
from viz import (
    plotly_molecule_figure, matplotlib_sticker, composite_ar, molecule_card,
    lewis_dot_covalent, lewis_dot_ionic,
)

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
            "- Electron-level views: ionic transfer diagrams and Lewis dot structures\n"
            "- Transparent, rule-based green-chemistry scoring\n"
            "- Photo-overlay AR preview\n"
            "- Quiz with badges, session-scoped leaderboard\n"
            "- Lab safety advisor (real, well-known hazard pairs)\n"
            "- English + Hindi explanations, 4 more languages at vocabulary level"
        )
    with right:
        st.markdown(
            '<div class="roadmap"><b>Roadmap — needs infrastructure this prototype doesn\'t have</b><br><br>'
            "• Real drag-and-drop (mouse) — needs a custom browser component; "
            "tapping was chosen instead since it also works on touchscreens<br>"
            "• Live ArUco marker-tracked AR — needs a native camera app or "
            "streamlit-webrtc with a reliable TURN server, not just a browser tab<br>"
            "• A trained graph neural network on EPA/PubChem toxicity data — "
            "needs a real labelled dataset and training pipeline<br>"
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
        st.warning(
            f"**No verified real structure for this exact combination** in BondVision AR's "
            f"library — that's not the same as saying it can't exist, just that this tool won't "
            f"guess at a structure it can't confirm. Total mass of what's in the workspace: "
            f"**{result['total_mass']} g/mol**."
        )
        st.caption(
            "Try adjusting the counts — for a covalent molecule, one element usually needs to be "
            "at a count of exactly 1 (the central atom); for an ionic compound, the ratio needs to "
            "exactly charge-balance (e.g. 1 Mg to 2 Cl, not 1 to 1)."
        )

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
    st.markdown("## 🥽 AR Concept Preview")
    st.warning(
        "**What this actually does, in plain terms:** this takes one still photo from your "
        "camera and pastes a rendered picture of the molecule on top of it, at a position you "
        "control with sliders. It does **not** track a printed marker or follow the camera in "
        "real time — a browser tab can only grab a single photo, not a live video feed for "
        "Python to process frame by frame. Real marker-tracked AR (detecting a printed card and "
        "locking the molecule to it as you move the phone) needs a native camera app, which is "
        "listed as roadmap on the Home page."
    )

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

        photo = st.camera_input("Point the camera at your desk and take a photo")
        c1, c2, c3 = st.columns(3)
        scale = c1.slider("Size", 0.15, 0.8, 0.4)
        x_frac = c2.slider("Left ↔ Right", 0.0, 1.0, 0.5)
        y_frac = c3.slider("Up ↔ Down", 0.0, 1.0, 0.5)

        if photo is not None:
            sticker = matplotlib_sticker(geo)
            result = composite_ar(photo.getvalue(), sticker, scale=scale, x_frac=x_frac, y_frac=y_frac)
            st.image(result, caption=f"{pick} — {m['name']}", width='stretch')
            buf = io.BytesIO(); result.save(buf, format="PNG")
            st.download_button("⬇️ Download this image", buf.getvalue(),
                                file_name=f"{pick}_ar_preview.png", mime="image/png")
        else:
            st.info("Take a photo above to see the molecule composited onto it.")

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
