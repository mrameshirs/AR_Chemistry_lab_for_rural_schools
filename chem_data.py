"""
chem_data.py — reference data for BondVision AR
================================================
Every number here is a standard, textbook value (Pauling electronegativity,
IUPAC atomic weight, covalent radius). Nothing here is invented.

Environmental / safety notes are well-established facts (IPCC greenhouse gas
classifications, Montreal Protocol substances, standard lab safety rules).
Where a claim is genuinely uncertain or a rounded figure, the note says so
in plain words rather than stating a fake precise number.
"""

# =============================================================================
# ELEMENTS
# =============================================================================
# en / hi   : display names
# Z         : atomic number
# en_neg    : Pauling electronegativity
# mass      : standard atomic weight (g/mol)
# radius_pm : covalent radius, picometres (relative sizing only, schematic)
# color     : CPK-style display colour (visualisation convention, not a fact)
# valence_e : main-group valence electron count
ELEMENTS = {
    "H":  {"name": "Hydrogen",  "name_hi": "हाइड्रोजन",   "Z": 1,  "en_neg": 2.20, "mass": 1.008,   "radius_pm": 31,  "color": "#F2F2F2", "valence_e": 1},
    "Li": {"name": "Lithium",   "name_hi": "लिथियम",      "Z": 3,  "en_neg": 0.98, "mass": 6.94,    "radius_pm": 128, "color": "#CC80FF", "valence_e": 1},
    "Be": {"name": "Beryllium", "name_hi": "बेरिलियम",    "Z": 4,  "en_neg": 1.57, "mass": 9.012,   "radius_pm": 96,  "color": "#C2FF00", "valence_e": 2},
    "B":  {"name": "Boron",     "name_hi": "बोरॉन",       "Z": 5,  "en_neg": 2.04, "mass": 10.811,  "radius_pm": 84,  "color": "#FFB5B5", "valence_e": 3},
    "C":  {"name": "Carbon",    "name_hi": "कार्बन",      "Z": 6,  "en_neg": 2.55, "mass": 12.011,  "radius_pm": 76,  "color": "#404040", "valence_e": 4},
    "N":  {"name": "Nitrogen",  "name_hi": "नाइट्रोजन",   "Z": 7,  "en_neg": 3.04, "mass": 14.007,  "radius_pm": 71,  "color": "#3050F8", "valence_e": 5},
    "O":  {"name": "Oxygen",    "name_hi": "ऑक्सीजन",    "Z": 8,  "en_neg": 3.44, "mass": 15.999,  "radius_pm": 66,  "color": "#FF0D0D", "valence_e": 6},
    "F":  {"name": "Fluorine",  "name_hi": "फ्लोरीन",     "Z": 9,  "en_neg": 3.98, "mass": 18.998,  "radius_pm": 57,  "color": "#8FE04A", "valence_e": 7},
    "Na": {"name": "Sodium",    "name_hi": "सोडियम",      "Z": 11, "en_neg": 0.93, "mass": 22.990,  "radius_pm": 166, "color": "#AB5CF2", "valence_e": 1},
    "Mg": {"name": "Magnesium", "name_hi": "मैग्नीशियम",  "Z": 12, "en_neg": 1.31, "mass": 24.305,  "radius_pm": 141, "color": "#8AFF00", "valence_e": 2},
    "Al": {"name": "Aluminium", "name_hi": "एल्युमिनियम", "Z": 13, "en_neg": 1.61, "mass": 26.982,  "radius_pm": 121, "color": "#BFA6A6", "valence_e": 3},
    "Si": {"name": "Silicon",   "name_hi": "सिलिकॉन",     "Z": 14, "en_neg": 1.90, "mass": 28.085,  "radius_pm": 111, "color": "#E8C89C", "valence_e": 4},
    "P":  {"name": "Phosphorus","name_hi": "फॉस्फोरस",    "Z": 15, "en_neg": 2.19, "mass": 30.974,  "radius_pm": 107, "color": "#FF8000", "valence_e": 5},
    "S":  {"name": "Sulfur",    "name_hi": "सल्फर",       "Z": 16, "en_neg": 2.58, "mass": 32.06,   "radius_pm": 105, "color": "#F2E017", "valence_e": 6},
    "Cl": {"name": "Chlorine",  "name_hi": "क्लोरीन",     "Z": 17, "en_neg": 3.16, "mass": 35.45,   "radius_pm": 102, "color": "#3DDC3D", "valence_e": 7},
    "K":  {"name": "Potassium", "name_hi": "पोटैशियम",    "Z": 19, "en_neg": 0.82, "mass": 39.098,  "radius_pm": 203, "color": "#8F40D4", "valence_e": 1},
    "Ca": {"name": "Calcium",   "name_hi": "कैल्शियम",    "Z": 20, "en_neg": 1.00, "mass": 40.078,  "radius_pm": 176, "color": "#3DC23D", "valence_e": 2},
    "Fe": {"name": "Iron",      "name_hi": "लोहा",        "Z": 26, "en_neg": 1.83, "mass": 55.845,  "radius_pm": 132, "color": "#E06633", "valence_e": 2},
    "Cu": {"name": "Copper",    "name_hi": "तांबा",       "Z": 29, "en_neg": 1.90, "mass": 63.546,  "radius_pm": 132, "color": "#C88033", "valence_e": 2},
    "Zn": {"name": "Zinc",      "name_hi": "जस्ता",       "Z": 30, "en_neg": 1.65, "mass": 65.38,   "radius_pm": 122, "color": "#7D80B0", "valence_e": 2},
    "Br": {"name": "Bromine",   "name_hi": "ब्रोमीन",     "Z": 35, "en_neg": 2.96, "mass": 79.904,  "radius_pm": 120, "color": "#A62929", "valence_e": 7},
    "I":  {"name": "Iodine",    "name_hi": "आयोडीन",      "Z": 53, "en_neg": 2.66, "mass": 126.904, "radius_pm": 139, "color": "#940094", "valence_e": 7},
}

# Elements the Builder allows as a CENTRAL atom, with the peripheral-atom
# counts that are chemically sensible for a single-bond-only model.
# (Values follow standard valence: B/Al are electron-deficient (0 lone
# pairs); P and S may expand their octet, which real PCl5 and SF6 use.)
CENTRAL_ALLOWED_N = {
    "B": [3], "Al": [3],
    "C": [4], "Si": [4],
    "N": [3], "P": [3, 5],
    "O": [2], "S": [2, 6],
}

# Elements the Builder allows as a PERIPHERAL atom (each forms one single bond).
PERIPHERAL_ELEMENTS = ["H", "F", "Cl", "Br", "I"]

# Common oxidation states for the Ionic Compound path. Most main-group metals
# have one dependable charge; the transition metals genuinely have more than
# one real, common charge, which is worth surfacing rather than hiding —
# picking which one is itself a real chemistry decision a student can make.
ION_CHARGES = {
    "Li": [1], "Na": [1], "K": [1],
    "Be": [2], "Mg": [2], "Ca": [2], "Zn": [2],
    "Al": [3],
    "Fe": [2, 3],   # iron(II) / iron(III) — both real and common
    "Cu": [1, 2],   # copper(I) / copper(II) — both real and common
    "F": [-1], "Cl": [-1], "Br": [-1], "I": [-1], "O": [-2],
}
IONIC_METALS = ["Li", "Na", "K", "Be", "Mg", "Ca", "Al", "Fe", "Cu", "Zn"]
IONIC_ANIONS = ["F", "Cl", "Br", "I", "O"]
ANION_STEM = {"F": "fluoride", "Cl": "chloride", "Br": "bromide", "I": "iodide", "O": "oxide"}

ROMAN = {1: "I", 2: "II", 3: "III"}

# =============================================================================
# ELECTRON-DOMAIN GEOMETRY LABELS  (steric number, lone pairs) -> info
# =============================================================================
GEOMETRY_INFO = {
    (2, 0): {"name": "Linear",              "name_hi": "रैखिक",              "ideal_angle": 180.0},
    (3, 1): {"name": "Bent (angular)",      "name_hi": "मुड़ा हुआ (कोणीय)",    "ideal_angle": 120.0},
    (3, 0): {"name": "Trigonal planar",     "name_hi": "त्रिकोणीय समतलीय",     "ideal_angle": 120.0},
    (4, 0): {"name": "Tetrahedral",         "name_hi": "चतुष्फलकीय",          "ideal_angle": 109.5},
    (4, 1): {"name": "Trigonal pyramidal",  "name_hi": "त्रिकोणीय पिरामिडी",   "ideal_angle": 109.5},
    (4, 2): {"name": "Bent (angular)",      "name_hi": "मुड़ा हुआ (कोणीय)",    "ideal_angle": 109.5},
    (5, 0): {"name": "Trigonal bipyramidal","name_hi": "त्रिकोणीय द्विपिरामिडी","ideal_angle": 90.0},
    (6, 0): {"name": "Octahedral",          "name_hi": "अष्टफलकीय",           "ideal_angle": 90.0},
}

# Real, well-known exact bond angles for a few famous molecules, used instead
# of the generic idealised angle when we have solid textbook figures.
EXACT_ANGLES = {
    "H2O": 104.5,
    "NH3": 107.0,
    "H2S": 92.0,
    "PCl3": 100.0,
    "O3": 117.0,
    "SO2": 119.0,
}

# =============================================================================
# CURATED MOLECULE LIBRARY
# =============================================================================
# central / peripheral / n / bond_order describe the STRUCTURE.
# Geometry, polarity and bond type are derived at run time by chem_engine,
# so the same maths that powers the Builder also powers this library —
# one formula, not two copies that could quietly disagree.
#
# steric_override / lone_pairs_override are used ONLY for the two resonance
# molecules (O3, SO2), whose true bonding is a hybrid that a simple integer
# bond order can't represent cleanly.
CURATED_MOLECULES = [
    # --- elemental diatomics ---
    {"formula": "H2",  "name": "Hydrogen",        "name_hi": "हाइड्रोजन गैस", "central": "H", "peripheral": "H", "n": 1, "bond_order": 1},
    {"formula": "N2",  "name": "Nitrogen",         "name_hi": "नाइट्रोजन गैस", "central": "N", "peripheral": "N", "n": 1, "bond_order": 3},
    {"formula": "O2",  "name": "Oxygen",           "name_hi": "ऑक्सीजन गैस", "central": "O", "peripheral": "O", "n": 1, "bond_order": 2},
    {"formula": "F2",  "name": "Fluorine",         "name_hi": "फ्लोरीन गैस", "central": "F", "peripheral": "F", "n": 1, "bond_order": 1},
    {"formula": "Cl2", "name": "Chlorine",         "name_hi": "क्लोरीन गैस", "central": "Cl", "peripheral": "Cl", "n": 1, "bond_order": 1},
    {"formula": "Br2", "name": "Bromine",          "name_hi": "ब्रोमीन",     "central": "Br", "peripheral": "Br", "n": 1, "bond_order": 1},
    {"formula": "I2",  "name": "Iodine",           "name_hi": "आयोडीन",      "central": "I", "peripheral": "I", "n": 1, "bond_order": 1},
    # --- hydrogen halides ---
    {"formula": "HCl", "name": "Hydrogen chloride","name_hi": "हाइड्रोजन क्लोराइड", "central": "H", "peripheral": "Cl", "n": 1, "bond_order": 1},
    {"formula": "HF",  "name": "Hydrogen fluoride","name_hi": "हाइड्रोजन फ्लोराइड", "central": "H", "peripheral": "F", "n": 1, "bond_order": 1},
    {"formula": "HBr", "name": "Hydrogen bromide", "name_hi": "हाइड्रोजन ब्रोमाइड", "central": "H", "peripheral": "Br", "n": 1, "bond_order": 1},
    {"formula": "HI",  "name": "Hydrogen iodide",  "name_hi": "हाइड्रोजन आयोडाइड", "central": "H", "peripheral": "I", "n": 1, "bond_order": 1},
    # --- small oxide / carbon gases ---
    {"formula": "CO",  "name": "Carbon monoxide",  "name_hi": "कार्बन मोनोऑक्साइड", "central": "C", "peripheral": "O", "n": 1, "bond_order": 3},
    {"formula": "NO",  "name": "Nitric oxide",     "name_hi": "नाइट्रिक ऑक्साइड", "central": "N", "peripheral": "O", "n": 1, "bond_order": 2},
    {"formula": "CO2", "name": "Carbon dioxide",   "name_hi": "कार्बन डाइऑक्साइड", "central": "C", "peripheral": "O", "n": 2, "bond_order": 2},
    # --- resonance molecules (hardcoded geometry, see EXACT_ANGLES) ---
    {"formula": "O3",  "name": "Ozone",  "name_hi": "ओज़ोन", "central": "O", "peripheral": "O", "n": 2, "bond_order": 1,
     "steric_override": 3, "lone_pairs_override": 1, "note_resonance": True},
    {"formula": "SO2", "name": "Sulfur dioxide", "name_hi": "सल्फर डाइऑक्साइड", "central": "S", "peripheral": "O", "n": 2, "bond_order": 1,
     "steric_override": 3, "lone_pairs_override": 1, "note_resonance": True},
    # --- classic VSEPR teaching set ---
    {"formula": "H2O",   "name": "Water",              "name_hi": "जल",              "central": "O",  "peripheral": "H",  "n": 2, "bond_order": 1},
    {"formula": "NH3",   "name": "Ammonia",             "name_hi": "अमोनिया",         "central": "N",  "peripheral": "H",  "n": 3, "bond_order": 1},
    {"formula": "CH4",   "name": "Methane",             "name_hi": "मीथेन",           "central": "C",  "peripheral": "H",  "n": 4, "bond_order": 1},
    {"formula": "BF3",   "name": "Boron trifluoride",   "name_hi": "बोरॉन ट्राइफ्लोराइड", "central": "B", "peripheral": "F", "n": 3, "bond_order": 1},
    {"formula": "AlCl3", "name": "Aluminium chloride",  "name_hi": "एल्युमिनियम क्लोराइड", "central": "Al", "peripheral": "Cl", "n": 3, "bond_order": 1,
     "extra_note": "Shown as the simple monomer; in real vapour, AlCl3 mostly exists as the dimer Al2Cl6."},
    {"formula": "SiH4",  "name": "Silane",              "name_hi": "साइलेन",          "central": "Si", "peripheral": "H",  "n": 4, "bond_order": 1},
    {"formula": "SiCl4", "name": "Silicon tetrachloride","name_hi": "सिलिकॉन टेट्राक्लोराइड", "central": "Si", "peripheral": "Cl", "n": 4, "bond_order": 1},
    {"formula": "PCl3",  "name": "Phosphorus trichloride","name_hi": "फॉस्फोरस ट्राइक्लोराइड", "central": "P", "peripheral": "Cl", "n": 3, "bond_order": 1},
    {"formula": "PCl5",  "name": "Phosphorus pentachloride","name_hi": "फॉस्फोरस पेंटाक्लोराइड", "central": "P", "peripheral": "Cl", "n": 5, "bond_order": 1},
    {"formula": "H2S",   "name": "Hydrogen sulfide",    "name_hi": "हाइड्रोजन सल्फाइड", "central": "S", "peripheral": "H", "n": 2, "bond_order": 1},
    {"formula": "SF6",   "name": "Sulfur hexafluoride", "name_hi": "सल्फर हेक्साफ्लोराइड", "central": "S", "peripheral": "F", "n": 6, "bond_order": 1},
    {"formula": "CCl4",  "name": "Carbon tetrachloride","name_hi": "कार्बन टेट्राक्लोराइड", "central": "C", "peripheral": "Cl", "n": 4, "bond_order": 1},
    # --- ionic compounds ---
    {"formula": "NaCl", "name": "Sodium chloride", "name_hi": "सोडियम क्लोराइड (नमक)", "is_ionic": True, "cation": "Na", "anion": "Cl"},
    {"formula": "MgCl2","name": "Magnesium chloride","name_hi": "मैग्नीशियम क्लोराइड", "is_ionic": True, "cation": "Mg", "anion": "Cl"},
    {"formula": "MgO",  "name": "Magnesium oxide",  "name_hi": "मैग्नीशियम ऑक्साइड", "is_ionic": True, "cation": "Mg", "anion": "O"},
]

# =============================================================================
# GREEN CHEMISTRY FACTS
# =============================================================================
# Levels are one of: "Low", "Medium", "High", "Very High" for persistence
# and toxicity; climate is one of: "None", "Minor", "Greenhouse gas",
# "Ozone-depleting". These qualitative calls follow well-established,
# widely taught environmental chemistry (IPCC greenhouse-gas science,
# Montreal Protocol substances, standard hazard classes) rather than a
# single disputed numeric figure.
GREEN_FACTS = {
    "H2":  {"persistence": "Low", "toxicity": "Low", "climate": "None",
            "note": "Burns to form only water. Central to India's National Green Hydrogen Mission as a clean fuel; the main hazard is flammability, not toxicity."},
    "N2":  {"persistence": "Low", "toxicity": "Low", "climate": "None",
            "note": "Makes up about 78% of the atmosphere and is chemically inert under normal conditions — not a pollutant."},
    "O2":  {"persistence": "Low", "toxicity": "Low", "climate": "None",
            "note": "Essential for respiration and combustion."},
    "F2":  {"persistence": "Low", "toxicity": "Very High", "climate": "None",
            "note": "One of the most reactive and corrosive gases known; extremely hazardous even in small amounts."},
    "Cl2": {"persistence": "Low", "toxicity": "Very High", "climate": "None",
            "note": "Toxic, corrosive gas. Historically used as a chemical weapon in World War I — a stark reminder of why lab safety rules exist."},
    "Br2": {"persistence": "Low", "toxicity": "High", "climate": "None",
            "note": "Corrosive liquid that gives off toxic, irritating vapour."},
    "I2":  {"persistence": "Low", "toxicity": "Medium", "climate": "None",
            "note": "Far less volatile and hazardous than chlorine or bromine; dilute iodine solutions are used as a mild antiseptic."},
    "HCl": {"persistence": "Low", "toxicity": "High", "climate": "None",
            "note": "Corrosive acidic gas; reacts with moisture in the lungs and eyes."},
    "HF":  {"persistence": "Low", "toxicity": "Very High", "climate": "None",
            "note": "Extremely dangerous — it can cause deep tissue damage that isn't immediately painful. Handled only with special training."},
    "HBr": {"persistence": "Low", "toxicity": "High", "climate": "None", "note": "Corrosive acidic gas, similar hazards to hydrogen chloride."},
    "HI":  {"persistence": "Low", "toxicity": "High", "climate": "None", "note": "Corrosive acidic gas, similar hazards to hydrogen chloride."},
    "CO":  {"persistence": "Medium", "toxicity": "Very High", "climate": "Minor",
            "note": "Colourless and odourless but binds to haemoglobin far more strongly than oxygen — a classic 'silent' poisoning hazard from incomplete combustion."},
    "NO":  {"persistence": "Low", "toxicity": "High", "climate": "Minor",
            "note": "Reactive and short-lived in air, quickly turning into other nitrogen oxides. In tiny amounts inside the body, it also acts as a signalling molecule — a 1998 Nobel Prize was awarded for that discovery."},
    "CO2": {"persistence": "High", "toxicity": "Low", "climate": "Greenhouse gas",
            "note": "Not directly toxic at normal air concentrations, but its long atmospheric lifetime makes it the main driver of human-caused climate change."},
    "O3":  {"persistence": "Medium", "toxicity": "Medium", "climate": "Minor",
            "note": "A genuine two-faced molecule: high in the stratosphere it shields us from UV radiation, but at ground level it's a lung-irritating air pollutant."},
    "SO2": {"persistence": "Medium", "toxicity": "High", "climate": "Minor",
            "note": "Major contributor to acid rain and a respiratory irritant; a classic industrial and vehicle-exhaust pollutant."},
    "H2O": {"persistence": "Low", "toxicity": "Low", "climate": "None",
            "note": "Essential for life. Technically a greenhouse gas too, but it's part of the natural water cycle rather than a pollutant of concern."},
    "NH3": {"persistence": "Low", "toxicity": "High", "climate": "Minor",
            "note": "Corrosive and sharp-smelling at high concentration, but breaks down fairly quickly and doesn't build up in the environment. Vital to the Haber process that makes most of the world's fertiliser."},
    "BF3": {"persistence": "Low", "toxicity": "High", "climate": "None", "note": "Toxic, corrosive gas used industrially as a catalyst; reacts readily with moisture."},
    "AlCl3": {"persistence": "Low", "toxicity": "High", "climate": "None",
              "note": "Reacts violently with water, releasing corrosive HCl fumes. Used industrially as a Friedel–Crafts catalyst."},
    "SiH4": {"persistence": "Low", "toxicity": "Medium", "climate": "None",
             "note": "Pyrophoric — it can ignite spontaneously in air. Widely used in semiconductor manufacturing under strict controls."},
    "SiCl4": {"persistence": "Low", "toxicity": "High", "climate": "None",
              "note": "Fumes strongly in moist air, forming corrosive HCl. Used in making optical fibre and semiconductor-grade silicon."},
    "PCl3": {"persistence": "Low", "toxicity": "High", "climate": "None", "note": "Reacts violently with water; corrosive and toxic."},
    "PCl5": {"persistence": "Low", "toxicity": "High", "climate": "None", "note": "Reacts violently with water; corrosive and toxic."},
    "H2S": {"persistence": "Low", "toxicity": "Very High", "climate": "None",
            "note": "Smells of rotten eggs at low concentration, but dangerously the sense of smell fails at higher, more lethal concentrations. Occurs naturally from volcanic activity and decaying organic matter."},
    "SF6": {"persistence": "Very High", "toxicity": "Low", "climate": "Greenhouse gas",
            "note": "A striking contradiction: chemically almost inert and non-toxic, yet one of the most potent greenhouse gases known, with an atmospheric lifetime measured in thousands of years. Used as an insulator in electrical switchgear; several countries are now phasing it down."},
    "CCl4": {"persistence": "High", "toxicity": "High", "climate": "Ozone-depleting",
             "note": "Liver-toxic and classified as an ozone-depleting substance; its production is now restricted under the Montreal Protocol. Formerly used as a dry-cleaning solvent and fire-extinguisher fluid."},
    "NaCl": {"persistence": "Low", "toxicity": "Low", "climate": "None", "note": "Ordinary table salt — naturally abundant and essential in small dietary amounts."},
    "MgCl2": {"persistence": "Low", "toxicity": "Low", "climate": "None", "note": "Used as a de-icing agent and in food processing (as nigari, a tofu coagulant)."},
    "MgO": {"persistence": "Low", "toxicity": "Low", "climate": "None", "note": "Used in antacids and as a heat-resistant refractory material."},
    "C6H6": {"persistence": "Medium", "toxicity": "Very High", "climate": "None",
             "note": "A known human carcinogen (IARC Group 1) with strict occupational exposure limits. Once a common industrial solvent, now heavily regulated."},
    "C6H12O6": {"persistence": "Low", "toxicity": "Low", "climate": "None",
                "note": "A simple sugar — readily biodegradable and non-toxic. The basic energy currency of nearly all living cells."},
    "C2H6O": {"persistence": "Low", "toxicity": "Medium", "climate": "None",
              "note": "Biodegradable and used as a renewable biofuel/green solvent when produced from biomass — a genuine green-chemistry success story. Still an intoxicant in the quantities people drink, so 'green' doesn't mean 'harmless to consume.'"},
    "C2H4O2": {"persistence": "Low", "toxicity": "Medium", "climate": "None",
               "note": "Readily biodegradable — it's the acid in ordinary vinegar. Corrosive in concentrated (glacial) form, but safe in the dilute household concentration."},
}

GREEN_FALLBACK_NOTE = (
    "No curated environmental data for this exact compound. This is a structural "
    "estimate from general trends only — halogen content, and whether the compound "
    "is a strong acid/base — not a verified toxicology or climate rating."
)

# =============================================================================
# LAB SAFETY — well known, label-level hazard combinations
# =============================================================================
SAFETY_COMBOS = [
    {"a": "Bleach (sodium hypochlorite)", "b": "Ammonia cleaners",
     "danger": "Produces toxic chloramine gases.",
     "why": "This is exactly why bleach and ammonia-based cleaners carry a printed warning never to be mixed — a genuine household hazard, not just a lab one."},
    {"a": "Bleach (sodium hypochlorite)", "b": "Vinegar or other acids",
     "danger": "Releases toxic chlorine gas.",
     "why": "Acid reacts with hypochlorite to free chlorine gas, the same gas historically used as a chemical weapon."},
    {"a": "Concentrated sulfuric acid", "b": "Water (added in the wrong order)",
     "danger": "Violent, spattering exothermic reaction.",
     "why": "Always add acid to water, slowly and with stirring — never water to acid. This is one of the first rules taught in any chemistry lab."},
    {"a": "Strong acid", "b": "Strong base",
     "danger": "Fast, highly exothermic neutralisation.",
     "why": "Mixing concentrated strong acids and bases releases a lot of heat quickly and can spatter — always dilute first and add slowly."},
    {"a": "Alkali metal (e.g. sodium)", "b": "Water",
     "danger": "Violent reaction releasing flammable hydrogen gas and heat.",
     "why": "A classic, dramatic demonstration reaction — impressive, but only ever done by a trained teacher behind a safety screen, in small amounts."},
    {"a": "Hydrogen peroxide", "b": "Vinegar",
     "danger": "Forms peracetic acid, a skin and eye irritant.",
     "why": "Both are common household items, but combining them is not advised without proper training."},
]

# =============================================================================
# MULTILINGUAL VOCABULARY (template-based, not free machine translation)
# =============================================================================
# Full, natural sentence templates exist for English and Hindi, which are the
# two languages we're confident reads naturally. For the others we provide a
# reviewed VOCABULARY table in a simple label format rather than a composed
# sentence — safer than guessing at grammar we can't verify, still genuinely
# useful for a classroom, and clearly marked as a starting point.
VOCAB = {
    "en": {
        "polar": "polar", "nonpolar": "nonpolar", "ionic": "ionic",
        "single": "single bond", "double": "double bond", "triple": "triple bond",
        "geometry": {k: v["name"] for k, v in GEOMETRY_INFO.items()},
    },
    "hi": {
        "polar": "ध्रुवीय", "nonpolar": "अध्रुवीय", "ionic": "आयनिक",
        "single": "एकल बंध", "double": "द्वि-बंध", "triple": "त्रि-बंध",
        "geometry": {k: v["name_hi"] for k, v in GEOMETRY_INFO.items()},
    },
    "mr": {  # Marathi — vocabulary table only
        "polar": "ध्रुवीय", "nonpolar": "अध्रुवीय", "ionic": "आयनिक",
        "linear": "रेषीय", "trigonal_planar": "त्रिकोणी सपाट", "tetrahedral": "चतुष्फलकीय",
        "bent": "वाकडे", "note": "शब्दसंग्रह संदर्भासाठी; कृपया शिक्षकांकडून पडताळणी करा.",
    },
    "ta": {  # Tamil
        "polar": "முனைவு", "nonpolar": "முனைவற்ற", "ionic": "அயனி",
        "linear": "நேரியல்", "trigonal_planar": "முக்கோண தட்டையான", "tetrahedral": "நான்முக வடிவம்",
        "bent": "வளைந்த", "note": "இது ஒரு சொல்லடைவு குறிப்பு மட்டுமே; ஆசிரியரிடம் சரிபார்க்கவும்.",
    },
    "te": {  # Telugu
        "polar": "ధ్రువ", "nonpolar": "అధ్రువ", "ionic": "అయానిక్",
        "linear": "రేఖీయ", "trigonal_planar": "త్రికోణ చదునైన", "tetrahedral": "చతుర్ముఖ",
        "bent": "వంగిన", "note": "ఇది కేవలం పదకోశ సూచన మాత్రమే; ఉపాధ్యాయుని ద్వారా నిర్ధారించుకోండి.",
    },
    "bn": {  # Bengali
        "polar": "মেরু", "nonpolar": "অমেরু", "ionic": "আয়নিক",
        "linear": "সরলরৈখিক", "trigonal_planar": "ত্রিকৌণিক সমতল", "tetrahedral": "চতুস্তলকীয়",
        "bent": "বাঁকা", "note": "এটি শুধুমাত্র একটি শব্দকোষ নির্দেশিকা; শিক্ষকের কাছে যাচাই করে নিন।",
    },
}

MOL_BY_FORMULA = {m["formula"]: m for m in CURATED_MOLECULES}

# =============================================================================
# COMPLEX ORGANIC MOLECULES — fact cards, not derived structures
# =============================================================================
# These have more than one central atom (carbon chains, rings) — genuinely
# outside what a single-central-atom AXn generator can build correctly. Real
# structural chemistry (ring closure, chain branching, isomer choice) is a
# harder problem than valence-electron counting. Rather than fabricate a
# structure that LOOKS right but isn't, these are shown as real reference
# facts only: formula, molar mass, and true environmental/safety notes.
COMPLEX_MOLECULES = [
    {"formula": "C6H6", "name": "Benzene", "name_hi": "बेंज़ीन", "atoms": {"C": 6, "H": 6},
     "structure_note": "A six-carbon ring with delocalised (shared across the whole ring) "
                        "bonding electrons — a genuinely different kind of bonding from the "
                        "simple single/double bonds this tool's builder can construct."},
    {"formula": "C6H12O6", "name": "Glucose", "name_hi": "ग्लूकोज़", "atoms": {"C": 6, "H": 12, "O": 6},
     "structure_note": "A ring of 5 carbons and 1 oxygen, with -OH groups branching off — real "
                        "structural chemistry (which carbon bonds to which, and the ring closure) "
                        "that a simple central-atom-plus-identical-peripherals model can't derive."},
    {"formula": "C2H6O", "name": "Ethanol", "name_hi": "एथेनॉल", "atoms": {"C": 2, "H": 6, "O": 1},
     "structure_note": "Two carbons joined to each other, one carrying an -OH group — a carbon "
                        "chain, not a single central atom with identical peripherals."},
    {"formula": "C2H4O2", "name": "Acetic acid", "name_hi": "ऐसिटिक अम्ल", "atoms": {"C": 2, "H": 4, "O": 2},
     "structure_note": "A two-carbon chain ending in a -COOH acid group — again a chain, not a "
                        "single central atom."},
]
COMPLEX_BY_FORMULA = {m["formula"]: m for m in COMPLEX_MOLECULES}

# Display order for the tap-to-build atom palette — metals first, then
# nonmetals, roughly by how often a student will actually reach for them.
PALETTE_ORDER = ["H", "C", "N", "O", "Na", "Cl", "Mg", "Al", "Ca", "K",
                  "Fe", "Cu", "Zn", "Si", "P", "S", "B", "F", "Br", "I", "Li", "Be"]
