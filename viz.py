"""
viz.py — turning a MoleculeGeometry into pictures
==================================================
Two renderers, for two different jobs:

- plotly_molecule_figure(): the interactive 3D viewer (rotate/zoom in the
  browser) used by the Builder and Library pages.
- matplotlib_sticker(): a small static PNG "sticker" of the same molecule,
  used only to composite onto a photo for the AR Preview page. Kept on a
  separate, lighter-weight library (matplotlib, not a headless-browser
  Plotly image export) so that feature doesn't add a heavy dependency.
"""

import io

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from PIL import Image

from chem_data import ELEMENTS, ION_CHARGES

BOND_COLOR = "#7A7A7A"
LONE_PAIR_COLOR = "#F2E017"


def _sphere_size(elem, base=16, scale=0.10):
    return base + ELEMENTS[elem]["radius_pm"] * scale


def plotly_molecule_figure(geo, show_lone_pairs=True, height=460):
    fig = go.Figure()

    # bonds — offset parallel lines represent double/triple bonds
    for b in geo.bonds:
        p1, p2 = geo.atoms[b.i].pos, geo.atoms[b.j].pos
        axis = p2 - p1
        length = np.linalg.norm(axis)
        if length < 1e-9:
            continue
        # a vector perpendicular to the bond axis, used to offset multiple lines
        arbitrary = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 * length else np.array([0.0, 1.0, 0.0])
        perp = np.cross(axis, arbitrary)
        perp = perp / (np.linalg.norm(perp) + 1e-9) * 0.09

        offsets = {1: [0.0], 2: [-1.0, 1.0], 3: [-1.2, 0.0, 1.2]}.get(b.order, [0.0])
        for off in offsets:
            o = perp * off
            fig.add_trace(go.Scatter3d(
                x=[p1[0] + o[0], p2[0] + o[0]], y=[p1[1] + o[1], p2[1] + o[1]], z=[p1[2] + o[2], p2[2] + o[2]],
                mode="lines", line=dict(color=BOND_COLOR, width=9), hoverinfo="skip", showlegend=False,
            ))

    # atoms, grouped by element so the legend is compact and hover text is per-element
    seen = set()
    for a in geo.atoms:
        info = ELEMENTS[a.element]
        fig.add_trace(go.Scatter3d(
            x=[a.pos[0]], y=[a.pos[1]], z=[a.pos[2]], mode="markers",
            marker=dict(size=_sphere_size(a.element), color=info["color"],
                        line=dict(color="#1a1a1a", width=1)),
            name=f"{a.element} ({info['name']})",
            legendgroup=a.element, showlegend=a.element not in seen,
            hovertext=f"{info['name']} ({a.element})<br>Electronegativity: {info['en_neg']}",
            hoverinfo="text",
        ))
        seen.add(a.element)

    if show_lone_pairs:
        lp_shown = False
        for lp in geo.lone_pairs:
            fig.add_trace(go.Scatter3d(
                x=[lp.pos[0]], y=[lp.pos[1]], z=[lp.pos[2]], mode="markers",
                marker=dict(size=10, color=LONE_PAIR_COLOR, opacity=0.55, symbol="circle"),
                name="Lone pair", legendgroup="lp", showlegend=not lp_shown,
                hovertext="Lone pair (not a bond)", hoverinfo="text",
            ))
            lp_shown = True

    fig.update_layout(
        height=height,
        margin=dict(l=0, r=0, t=10, b=0),
        scene=dict(
            xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False),
            aspectmode="data",
        ),
        legend=dict(orientation="h", yanchor="bottom", y=0.01, x=0.01),
        paper_bgcolor="rgba(0,0,0,0)",
    )
    return fig


def matplotlib_sticker(geo, figsize=2.6, dpi=150, elev=18, azim=35):
    """A small transparent-background PNG for compositing onto a photo (AR Preview)."""
    fig = plt.figure(figsize=(figsize, figsize), dpi=dpi)
    ax = fig.add_subplot(111, projection="3d")
    ax.set_axis_off()
    try:
        ax.set_box_aspect((1, 1, 1))
    except Exception:
        pass  # older matplotlib without set_box_aspect

    for b in geo.bonds:
        p1, p2 = geo.atoms[b.i].pos, geo.atoms[b.j].pos
        for order_shift in range(b.order):
            shift = (order_shift - (b.order - 1) / 2) * 0.10
            ax.plot([p1[0], p2[0]], [p1[1] + shift, p2[1] + shift], [p1[2], p2[2]],
                    color="#666666", linewidth=3, solid_capstyle="round")

    for a in geo.atoms:
        info = ELEMENTS[a.element]
        ax.scatter(*a.pos, s=_sphere_size(a.element) * 9, color=info["color"],
                   edgecolors="#222222", linewidths=0.6, depthshade=True)

    ax.view_init(elev=elev, azim=azim)
    lim = max(1.6, *[np.max(np.abs(a.pos)) for a in geo.atoms]) * 1.15
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)

    buf = io.BytesIO()
    fig.savefig(buf, format="png", transparent=True, bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert("RGBA")


def composite_ar(photo_bytes, sticker_img, scale=0.42, x_frac=0.5, y_frac=0.5):
    """Paste the molecule sticker onto a captured photo. Returns a PIL Image."""
    base = Image.open(io.BytesIO(photo_bytes)).convert("RGBA")
    bw, bh = base.size
    target_w = int(bw * scale)
    ratio = target_w / sticker_img.width
    target_h = int(sticker_img.height * ratio)
    sticker = sticker_img.resize((target_w, target_h), Image.LANCZOS)

    px = int(bw * x_frac - target_w / 2)
    py = int(bh * y_frac - target_h / 2)
    out = base.copy()
    out.alpha_composite(sticker, dest=(px, py))
    return out.convert("RGB")


def lewis_dot_covalent(geo, figsize=5.2, dpi=150):
    """
    2D Lewis-structure view of a covalent MoleculeGeometry: atom symbols,
    bond lines, and dots for both shared (bonding) and lone-pair electrons.
    Reuses the SAME 3D coordinates as the VSEPR viewer, just projected to
    (x, y) — Lewis structures are conventionally schematic 2D diagrams
    anyway, a deliberately different (and complementary) view from the true
    3D shape, exactly as real chemistry teaching presents both side by side.
    """
    fig, ax = plt.subplots(figsize=(figsize, figsize), dpi=dpi)
    ax.set_aspect("equal")
    ax.axis("off")

    pts2d = {i: (a.pos[0], a.pos[1]) for i, a in enumerate(geo.atoms)}

    for b in geo.bonds:
        x1, y1 = pts2d[b.i]
        x2, y2 = pts2d[b.j]
        dx, dy = x2 - x1, y2 - y1
        length = (dx ** 2 + dy ** 2) ** 0.5 or 1e-9
        perp = (-dy / length, dx / length)
        offsets = {1: [0.0], 2: [-0.09, 0.09], 3: [-0.14, 0.0, 0.14]}.get(b.order, [0.0])
        for off in offsets:
            ox, oy = perp[0] * off, perp[1] * off
            ax.plot([x1 + ox, x2 + ox], [y1 + oy, y2 + oy], color="#555555", linewidth=2, zorder=1)
            # a dot-pair at the midpoint of each line = one shared electron pair
            mx, my = (x1 + x2) / 2 + ox, (y1 + y2) / 2 + oy
            dot_perp = (perp[0] * 0.045, perp[1] * 0.045)
            for s in (-1, 1):
                ax.plot(mx + dot_perp[0] * s, my + dot_perp[1] * s, "o", color="#222222",
                        markersize=5, zorder=3)

    for lp in geo.lone_pairs:
        x, y = lp.pos[0], lp.pos[1]
        length = (x ** 2 + y ** 2) ** 0.5 or 1e-9
        perp = (-y / length, x / length)
        for s in (-1, 1):
            ax.plot(x + perp[0] * 0.06 * s, y + perp[1] * 0.06 * s, "o", color="#B8860B",
                    markersize=6, zorder=3)

    for i, a in enumerate(geo.atoms):
        x, y = pts2d[i]
        info = ELEMENTS[a.element]
        ax.scatter([x], [y], s=1500, color=info["color"], edgecolors="#222", linewidths=1.2, zorder=2)
        text_color = "#FFFFFF" if a.element == "C" else "#000000"
        ax.text(x, y, a.element, ha="center", va="center", fontsize=13, fontweight="bold",
                color=text_color, zorder=4)

    lim = max(1.4, max(abs(p[0]) for p in pts2d.values()), max(abs(p[1]) for p in pts2d.values())) * 1.35
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)

    legend_y = -lim * 0.94
    ax.plot(-lim * 0.5, legend_y, "o", color="#222222", markersize=5)
    ax.text(-lim * 0.5 + lim * 0.08, legend_y, "shared pair", va="center", fontsize=8, color="#444")
    ax.plot(lim * 0.05, legend_y, "o", color="#B8860B", markersize=6)
    ax.text(lim * 0.05 + lim * 0.08, legend_y, "lone pair", va="center", fontsize=8, color="#444")

    return fig


def lewis_dot_ionic(cation, anion, cation_charge, n_cation, n_anion, figsize=5.6, dpi=150):
    """
    Electron-transfer diagram for an ionic compound: cations on the left,
    anions on the right, with arrows showing exactly which valence electrons
    move where. The total electrons leaving the cations always equals the
    total electrons the anions need — that's what "charge balanced" means —
    so the diagram is drawn from that same accounting, not just decoration.
    """
    fig, ax = plt.subplots(figsize=(figsize, figsize * 0.62), dpi=dpi)
    ax.axis("off")

    c_info, a_info = ELEMENTS[cation], ELEMENTS[anion]
    anion_charge = ION_CHARGES[anion][0]           # anions here have exactly one real charge
    anion_need = abs(anion_charge)                 # electrons this anion needs to complete its octet

    def superscript_charge(n, sign):
        digit = {1: "", 2: "\u00b2", 3: "\u00b3"}.get(n, str(n))
        return digit + ("\u207a" if sign > 0 else "\u207b")

    # queue of (source_cation_index) electron units leaving the cations
    electron_queue = []
    for ci in range(n_cation):
        electron_queue += [ci] * cation_charge

    cat_x, an_x = -1.6, 1.6
    cat_ys = np.linspace((n_cation - 1) / 2, -(n_cation - 1) / 2, n_cation) * 1.1
    an_ys = np.linspace((n_anion - 1) / 2, -(n_anion - 1) / 2, n_anion) * 1.1

    def draw_atom(x, y, elem, dots, charge_label=None):
        info = ELEMENTS[elem]
        ax.scatter([x], [y], s=1700, color=info["color"], edgecolors="#222", linewidths=1.3, zorder=2)
        ax.text(x, y, elem, ha="center", va="center", fontsize=14, fontweight="bold",
                color=("#FFFFFF" if elem == "C" else "#000000"), zorder=4)
        # dots ringed around the atom
        for k in range(dots):
            ang = 2 * np.pi * k / max(dots, 1) + np.pi / 6
            dx, dy = 0.22 * np.cos(ang), 0.22 * np.sin(ang)
            ax.plot(x + dx, y + dy, "o", color="#222222", markersize=5, zorder=3)
        if charge_label:
            ax.text(x, y + 0.32, charge_label, ha="center", va="bottom", fontsize=11,
                    color="#B22222", fontweight="bold", zorder=5)

    used = 0
    for ci, y in enumerate(cat_ys):
        draw_atom(cat_x, y, cation, dots=0, charge_label=f"{cation}{superscript_charge(cation_charge, 1)}")
    for ai, y in enumerate(an_ys):
        draw_atom(an_x, y, anion, dots=8, charge_label=f"{anion}{superscript_charge(anion_need, -1)}")

    # arrows: distribute queue electrons to anions in order, drawing one
    # curved arrow per electron actually transferred
    q = list(electron_queue)
    need_per_anion = anion_need
    for ai, ay in enumerate(an_ys):
        for _ in range(need_per_anion):
            if not q:
                break
            ci = q.pop(0)
            cy = cat_ys[ci]
            ax.annotate("", xy=(an_x - 0.35, ay), xytext=(cat_x + 0.35, cy),
                        arrowprops=dict(arrowstyle="-|>", color="#1B6B93", lw=1.4,
                                        connectionstyle=f"arc3,rad={0.15 if cy != ay else 0.0}"),
                        zorder=1)

    ax.set_xlim(-2.6, 2.6)
    ax.set_ylim(-max(1.6, len(cat_ys), len(an_ys)), max(1.6, len(cat_ys), len(an_ys)))
    ax.text(0, -max(1.6, len(cat_ys), len(an_ys)) * 0.94,
            "Arrows show real valence electrons moving from metal to nonmetal, "
            "until every atom reaches a full outer shell.",
            ha="center", fontsize=8.5, color="#555")
    return fig


def molecule_card(geo, title, formula, grade_text, lang_line=None, figsize=(4, 5), dpi=170):
    """
    A simple printable 'flashcard' PNG: structure + name + green-chemistry grade.

    lang_line is intentionally NOT used for non-Latin scripts here: Matplotlib's
    default font (DejaVu Sans) has no Devanagari/Tamil/etc. glyphs, and we can't
    guarantee a compatible font is installed on the deployment host. Rather than
    ship a downloadable card with visible empty boxes where Hindi text should
    be, we keep non-English text on the Streamlit page itself (rendered by the
    browser, which handles this correctly) and keep the PNG card to scripts
    Matplotlib can actually draw.
    """
    lang_line = lang_line if (lang_line and lang_line.isascii()) else None
    fig = plt.figure(figsize=figsize, dpi=dpi)
    fig.patch.set_facecolor("#FFFDF5")

    ax3d = fig.add_axes([0.08, 0.42, 0.84, 0.5], projection="3d")
    ax3d.set_axis_off()
    try:
        ax3d.set_box_aspect((1, 1, 1))
    except Exception:
        pass
    for b in geo.bonds:
        p1, p2 = geo.atoms[b.i].pos, geo.atoms[b.j].pos
        for order_shift in range(b.order):
            shift = (order_shift - (b.order - 1) / 2) * 0.10
            ax3d.plot([p1[0], p2[0]], [p1[1] + shift, p2[1] + shift], [p1[2], p2[2]],
                      color="#666666", linewidth=3)
    for a in geo.atoms:
        info = ELEMENTS[a.element]
        ax3d.scatter(*a.pos, s=_sphere_size(a.element) * 8, color=info["color"], edgecolors="#222", linewidths=0.6)
    ax3d.view_init(elev=18, azim=35)
    lim = max(1.6, *[np.max(np.abs(a.pos)) for a in geo.atoms]) * 1.15
    ax3d.set_xlim(-lim, lim); ax3d.set_ylim(-lim, lim); ax3d.set_zlim(-lim, lim)

    ax_txt = fig.add_axes([0, 0, 1, 0.4])
    ax_txt.set_axis_off()
    ax_txt.text(0.5, 0.86, title, ha="center", va="top", fontsize=15, fontweight="bold", color="#12241C")
    ax_txt.text(0.5, 0.68, formula, ha="center", va="top", fontsize=12, color="#5C6B63")
    ax_txt.text(0.5, 0.48, geo.geometry_name, ha="center", va="top", fontsize=11, color="#0E7C6B")
    ax_txt.text(0.5, 0.30, grade_text, ha="center", va="top", fontsize=10.5, color="#3B2A00", wrap=True)
    if lang_line:
        ax_txt.text(0.5, 0.10, lang_line, ha="center", va="top", fontsize=10, color="#5C6B63")

    for spine_ax in [fig.add_axes([0, 0, 1, 1])]:
        spine_ax.set_axis_off()
        spine_ax.add_patch(plt.Rectangle((0.02, 0.02), 0.96, 0.96, fill=False,
                                          edgecolor="#F0A227", linewidth=3, linestyle="dashed",
                                          transform=spine_ax.transAxes))

    buf = io.BytesIO()
    fig.savefig(buf, format="png", facecolor=fig.get_facecolor())
    plt.close(fig)
    buf.seek(0)
    return buf
