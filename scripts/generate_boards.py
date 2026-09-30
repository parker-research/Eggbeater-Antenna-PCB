#!/usr/bin/env python3
"""Parametric generator for the Eggbeater Antenna PCBs (KiCad).

Produces two KiCad projects:

  kicad/main_board/      - round 2-layer board that mounts inside a 2in PVC cap,
                           carries the four No.8 loop bolts and two slots that
                           receive the phasing board tabs.
  kicad/phasing_board/   - rectangular 4-layer board: edge-launch feed connector,
                           RG-316 / BN-43-3312 common-mode choke (balun), and a
                           100 ohm balanced, shielded (broadside-coupled stripline)
                           quarter-wave phasing line feeding the two tabs.

The board files are written in the KiCad 9 s-expression format and then upgraded
to the installed KiCad version with `kicad-cli pcb upgrade` (see Makefile /
scripts/build.sh).  Only the Python standard library is required.

Usage:  python3 scripts/generate_boards.py [--freq-mhz 436.5]
"""

import argparse
import json
import math
import os
import uuid

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
UUID_NS = uuid.UUID("5b1d8c6e-6c55-4f7a-9d38-2f2f4b1f2a10")

C0 = 299_792_458.0

# ---------------------------------------------------------------------------
# Stackup / transmission line parameters (see scripts/stripline_solver.py)
# JLCPCB JLC04161H-7628, 1.6 mm, 4 layer:
#   F.Cu 35um | 7628 prepreg 0.2104mm er4.4 | In1 15.2um | core 1.065mm er4.6 |
#   In2 15.2um | 7628 prepreg 0.2104mm er4.4 | B.Cu 35um
# Broadside-coupled pair: P on In1.Cu, N on In2.Cu, stacked, GND on F/B + fences.
# 2D field solver result: w = 0.25 mm -> Zdiff = 100.0 ohm, eeff = 4.476
# ---------------------------------------------------------------------------
TRACE_W = 0.25
E_EFF = 4.476
INNER_GND_CLEARANCE = 1.0      # gap from stripline trace edge to inner-layer GND pour
LANE = 3.5                     # centre-to-centre spacing between adjacent line sections
FENCE_OFFSET = LANE / 2        # via fence offset from line centre
FENCE_PITCH = 2.0
VIA_D, VIA_DRILL = 0.6, 0.3

# Mechanical interface between the boards
TAB_W = 7.0                    # phasing board tab width
TAB_L = 3.2                    # tab length = 1.6 mm main board + 1.6 mm protrusion
TAB_Y = 5.0                    # tab centres at +-TAB_Y
TAB_PAD_PITCH = 2.2
TAB_PAD_W = 1.3
SLOT_W = 1.8                   # main board slot width  (1.6 mm board + clearance)
SLOT_L = 8.4                   # main board slot length (tab + corner-radius relief)


def uid(*key):
    return str(uuid.uuid5(UUID_NS, "/".join(str(k) for k in key)))


def f(v):
    s = f"{v:.4f}".rstrip("0").rstrip(".")
    return "0" if s in ("-0", "") else s


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def chamfer(points, c=1.0):
    """Replace every interior corner with a 45-degree chamfer of leg length c."""
    out = [points[0]]
    for i in range(1, len(points) - 1):
        p0, p1, p2 = points[i - 1], points[i], points[i + 1]
        d1 = math.dist(p0, p1)
        d2 = math.dist(p1, p2)
        u1 = ((p1[0] - p0[0]) / d1, (p1[1] - p0[1]) / d1)
        u2 = ((p2[0] - p1[0]) / d2, (p2[1] - p1[1]) / d2)
        if abs(u1[0] * u2[1] - u1[1] * u2[0]) < 1e-9:
            out.append(p1)
            continue
        cc = min(c, d1 / 2.01, d2 / 2.01)
        out.append((p1[0] - u1[0] * cc, p1[1] - u1[1] * cc))
        out.append((p1[0] + u2[0] * cc, p1[1] + u2[1] * cc))
    out.append(points[-1])
    return out


def polylen(pts):
    return sum(math.dist(a, b) for a, b in zip(pts, pts[1:]))


def seg_dist(p, a, b):
    ax, ay = a
    bx, by = b
    dx, dy = bx - ax, by - ay
    L2 = dx * dx + dy * dy
    t = 0 if L2 == 0 else max(0, min(1, ((p[0] - ax) * dx + (p[1] - ay) * dy) / L2))
    return math.dist(p, (ax + t * dx, ay + t * dy))


def rect_dist(p, r):
    """Distance from point to axis aligned rectangle r=(x0,y0,x1,y1) (0 inside)."""
    x0, y0, x1, y1 = r
    dx = max(x0 - p[0], 0, p[0] - x1)
    dy = max(y0 - p[1], 0, p[1] - y1)
    return math.hypot(dx, dy)


def point_in_poly(p, poly):
    x, y = p
    inside = False
    n = len(poly)
    for i in range(n):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % n]
        if (y0 > y) != (y1 > y):
            xi = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
            if xi > x:
                inside = not inside
    return inside


def poly_edge_dist(p, poly):
    return min(seg_dist(p, poly[i], poly[(i + 1) % len(poly)]) for i in range(len(poly)))


def offset_polyline_samples(pts, off, pitch):
    """Sample points at +-off from a polyline, every `pitch` mm."""
    out = []
    for a, b in zip(pts, pts[1:]):
        L = math.dist(a, b)
        if L < 1e-6:
            continue
        ux, uy = (b[0] - a[0]) / L, (b[1] - a[1]) / L
        nx, ny = -uy, ux
        n = max(1, int(L // pitch))
        for k in range(n + 1):
            t = k * L / n
            cx, cy = a[0] + ux * t, a[1] + uy * t
            for s in (1, -1):
                out.append((cx + s * nx * off, cy + s * ny * off))
    return out


# ---------------------------------------------------------------------------
# S-expression writer (KiCad 9 board format)
# ---------------------------------------------------------------------------

LAYERS_4 = """  (layers
    (0 "F.Cu" signal)
    (4 "In1.Cu" signal)
    (6 "In2.Cu" signal)
    (2 "B.Cu" signal)
    (9 "F.Adhes" user "F.Adhesive")
    (11 "B.Adhes" user "B.Adhesive")
    (13 "F.Paste" user)
    (15 "B.Paste" user)
    (5 "F.SilkS" user "F.Silkscreen")
    (7 "B.SilkS" user "B.Silkscreen")
    (1 "F.Mask" user)
    (3 "B.Mask" user)
    (17 "Dwgs.User" user "User.Drawings")
    (19 "Cmts.User" user "User.Comments")
    (21 "Eco1.User" user "User.Eco1")
    (23 "Eco2.User" user "User.Eco2")
    (25 "Edge.Cuts" user)
    (27 "Margin" user)
    (31 "F.CrtYd" user "F.Courtyard")
    (29 "B.CrtYd" user "B.Courtyard")
    (35 "F.Fab" user)
    (33 "B.Fab" user)
  )"""

LAYERS_2 = LAYERS_4.replace('    (4 "In1.Cu" signal)\n    (6 "In2.Cu" signal)\n', "")

STACKUP_4 = """    (stackup
      (layer "F.SilkS" (type "Top Silk Screen"))
      (layer "F.Paste" (type "Top Solder Paste"))
      (layer "F.Mask" (type "Top Solder Mask") (thickness 0.01))
      (layer "F.Cu" (type "copper") (thickness 0.035))
      (layer "dielectric 1" (type "prepreg") (thickness 0.2104) (material "7628") (epsilon_r 4.4) (loss_tangent 0.02))
      (layer "In1.Cu" (type "copper") (thickness 0.0152))
      (layer "dielectric 2" (type "core") (thickness 1.065) (material "FR4") (epsilon_r 4.6) (loss_tangent 0.02))
      (layer "In2.Cu" (type "copper") (thickness 0.0152))
      (layer "dielectric 3" (type "prepreg") (thickness 0.2104) (material "7628") (epsilon_r 4.4) (loss_tangent 0.02))
      (layer "B.Cu" (type "copper") (thickness 0.035))
      (layer "B.Mask" (type "Bottom Solder Mask") (thickness 0.01))
      (layer "B.Paste" (type "Bottom Solder Paste"))
      (layer "B.SilkS" (type "Bottom Silk Screen"))
      (copper_finish "HAL lead-free")
      (dielectric_constraints no)
    )"""

STACKUP_2 = """    (stackup
      (layer "F.SilkS" (type "Top Silk Screen"))
      (layer "F.Paste" (type "Top Solder Paste"))
      (layer "F.Mask" (type "Top Solder Mask") (thickness 0.01))
      (layer "F.Cu" (type "copper") (thickness 0.035))
      (layer "dielectric 1" (type "core") (thickness 1.51) (material "FR4") (epsilon_r 4.5) (loss_tangent 0.02))
      (layer "B.Cu" (type "copper") (thickness 0.035))
      (layer "B.Mask" (type "Bottom Solder Mask") (thickness 0.01))
      (layer "B.Paste" (type "Bottom Solder Paste"))
      (layer "B.SilkS" (type "Bottom Silk Screen"))
      (copper_finish "HAL lead-free")
      (dielectric_constraints no)
    )"""


class Board:
    def __init__(self, name, layers4, origin):
        self.name = name
        self.layers4 = layers4
        self.ox, self.oy = origin
        self.nets = [""]
        self.items = []
        self.n = 0

    def P(self, x, y):
        return self.ox + x, self.oy + y

    def net(self, name):
        if name is None:
            return 0
        if name not in self.nets:
            self.nets.append(name)
        return self.nets.index(name)

    def _id(self, kind):
        self.n += 1
        return uid(self.name, kind, self.n)

    # --- graphics -------------------------------------------------------
    def line(self, a, b, layer, w=0.15):
        (x0, y0), (x1, y1) = self.P(*a), self.P(*b)
        self.items.append(
            f'  (gr_line (start {f(x0)} {f(y0)}) (end {f(x1)} {f(y1)}) '
            f'(stroke (width {f(w)}) (type solid)) (layer "{layer}") (uuid "{self._id("gl")}"))')

    def arc(self, a, m, b, layer, w=0.15):
        (x0, y0), (xm, ym), (x1, y1) = self.P(*a), self.P(*m), self.P(*b)
        self.items.append(
            f'  (gr_arc (start {f(x0)} {f(y0)}) (mid {f(xm)} {f(ym)}) (end {f(x1)} {f(y1)}) '
            f'(stroke (width {f(w)}) (type solid)) (layer "{layer}") (uuid "{self._id("ga")}"))')

    def circle(self, c, r, layer, w=0.15, fill=False):
        (x0, y0) = self.P(*c)
        self.items.append(
            f'  (gr_circle (center {f(x0)} {f(y0)}) (end {f(x0 + r)} {f(y0)}) '
            f'(stroke (width {f(w)}) (type solid)) (fill {"yes" if fill else "no"}) '
            f'(layer "{layer}") (uuid "{self._id("gc")}"))')

    def poly_outline(self, pts, layer="Edge.Cuts", w=0.1):
        for a, b in zip(pts, pts[1:] + pts[:1]):
            self.line(a, b, layer, w)

    def rect(self, x0, y0, x1, y1, layer, w=0.15):
        self.poly_outline([(x0, y0), (x1, y0), (x1, y1), (x0, y1)], layer, w)

    def text(self, s, at, layer, size=1.0, thick=None, angle=0, justify=None, bold=False):
        x, y = self.P(*at)
        thick = thick if thick is not None else max(0.15, round(size * 0.15, 3))
        mirror = layer.startswith("B.")
        just = []
        if justify:
            just.append(justify)
        if mirror:
            just.append("mirror")
        js = f" (justify {' '.join(just)})" if just else ""
        b = " (bold yes)" if bold else ""
        s = s.replace('"', '\\"')
        self.items.append(
            f'  (gr_text "{s}" (at {f(x)} {f(y)} {f(angle)}) (layer "{layer}") (uuid "{self._id("gt")}")\n'
            f'    (effects (font (size {f(size)} {f(size)}) (thickness {f(thick)}){b}){js}))')

    def arrow(self, tail, head, layer, w=0.2, hl=1.0):
        self.line(tail, head, layer, w)
        ang = math.atan2(head[1] - tail[1], head[0] - tail[0])
        for s in (1, -1):
            a = ang + math.pi + s * math.radians(30)
            self.line(head, (head[0] + hl * math.cos(a), head[1] + hl * math.sin(a)), layer, w)

    # --- copper ---------------------------------------------------------
    def track(self, pts, width, layer, net):
        n = self.net(net)
        for a, b in zip(pts, pts[1:]):
            (x0, y0), (x1, y1) = self.P(*a), self.P(*b)
            self.items.append(
                f'  (segment (start {f(x0)} {f(y0)}) (end {f(x1)} {f(y1)}) (width {f(width)}) '
                f'(layer "{layer}") (net {n}) (uuid "{self._id("seg")}"))')

    def via(self, p, net, d=VIA_D, drill=VIA_DRILL):
        x, y = self.P(*p)
        self.items.append(
            f'  (via (at {f(x)} {f(y)}) (size {f(d)}) (drill {f(drill)}) (layers "F.Cu" "B.Cu") '
            f'(net {self.net(net)}) (uuid "{self._id("via")}"))')

    def zone(self, layer, net, pts, clearance, name=None, priority=0):
        n = self.net(net)
        ps = " ".join(f"(xy {f(self.P(*p)[0])} {f(self.P(*p)[1])})" for p in pts)
        nm = f' (name "{name}")' if name else ""
        self.items.append(
            f'  (zone (net {n}) (net_name "{net}") (layer "{layer}") (uuid "{self._id("zone")}"){nm} (hatch edge 0.5)\n'
            f'    (priority {priority})\n'
            f'    (connect_pads yes (clearance {f(clearance)}))\n'
            f'    (min_thickness 0.2) (filled_areas_thickness no)\n'
            f'    (fill yes (thermal_gap 0.5) (thermal_bridge_width 0.5) (island_removal_mode 1) (island_area_min 1))\n'
            f'    (polygon (pts {ps})))')

    def footprint(self, ref, value, at, pads, graphics=(), descr="", smd=False):
        """pads: list of dicts(num, kind, shape, at, size, layers, net, drill)."""
        x, y = self.P(*at)
        out = [f'  (footprint "EggbeaterAntenna:{value}" (layer "F.Cu") (uuid "{self._id("fp")}")',
               f'    (at {f(x)} {f(y)})',
               f'    (descr "{descr}")',
               f'    (property "Reference" "{ref}" (at 0 0 0) (layer "F.Fab") (hide yes) (uuid "{self._id("fpr")}")'
               f' (effects (font (size 1 1) (thickness 0.15))))',
               f'    (property "Value" "{value}" (at 0 0 0) (layer "F.Fab") (hide yes) (uuid "{self._id("fpv")}")'
               f' (effects (font (size 1 1) (thickness 0.15))))',
               f'    (attr {"smd" if smd else "through_hole"} exclude_from_pos_files exclude_from_bom)']
        for g in graphics:
            out.append("    " + g)
        for p in pads:
            px, py = p["at"]
            sx, sy = p["size"]
            layers = " ".join(f'"{l}"' for l in p["layers"])
            drill = f' (drill {f(p["drill"])})' if p.get("drill") else ""
            net = ""
            if p.get("net"):
                net = f' (net {self.net(p["net"])} "{p["net"]}")'
            extra = ""
            if p["kind"] == "smd" and p.get("zone_connect") is not None:
                extra = f' (zone_connect {p["zone_connect"]})'
            out.append(
                f'    (pad "{p["num"]}" {p["kind"]} {p["shape"]} (at {f(px)} {f(py)}) (size {f(sx)} {f(sy)})'
                f'{drill} (layers {layers}){net}{extra} (uuid "{self._id("pad")}"))')
        out.append("  )")
        self.items.append("\n".join(out))

    # --- output ---------------------------------------------------------
    def write(self, path, title, rev, comments=()):
        nets = "\n".join(f'  (net {i} "{n}")' for i, n in enumerate(self.nets))
        cm = "\n".join(f'      (comment {i + 1} "{c}")' for i, c in enumerate(comments))
        hdr = f"""(kicad_pcb
  (version 20241229)
  (generator "eggbeater_generate_boards")
  (generator_version "9.0")
  (general (thickness 1.6) (legacy_teardrops no))
  (paper "A4")
  (title_block
    (title "{title}")
    (rev "{rev}")
    (company "Eggbeater-Antenna-PCB (open hardware)")
{cm}
  )
{LAYERS_4 if self.layers4 else LAYERS_2}
  (setup
{STACKUP_4 if self.layers4 else STACKUP_2}
    (pad_to_mask_clearance 0)
    (allow_soldermask_bridges_in_footprints no)
    (tenting front back)
    (pcbplotparams
      (layerselection 0x00000000_00000000_55555555_5755f5ff)
      (plot_on_all_layers_selection 0x00000000_00000000_00000000_00000000)
      (usegerberextensions no)
      (usegerberattributes yes)
      (usegerberadvancedattributes yes)
      (creategerberjobfile yes)
      (outputformat 1)
      (mirror no)
      (drillshape 1)
      (scaleselection 1)
      (outputdirectory "")
    )
  )
{nets}
"""
        with open(path, "w") as fh:
            fh.write(hdr)
            fh.write("\n".join(self.items))
            fh.write("\n)\n")


def write_project(path, min_clear=0.15):
    pro = {
        "board": {
            "design_settings": {
                "defaults": {
                    "board_outline_line_width": 0.1,
                    "copper_line_width": 0.2,
                    "copper_text_size_h": 1.5, "copper_text_size_v": 1.5,
                    "copper_text_thickness": 0.3,
                    "silk_line_width": 0.15,
                    "silk_text_size_h": 1.0, "silk_text_size_v": 1.0,
                    "silk_text_thickness": 0.15,
                },
                # JLCPCB standard capabilities (with margin)
                "rules": {
                    "min_clearance": min_clear,
                    "min_track_width": 0.127,
                    "min_via_annular_width": 0.13,
                    "min_via_diameter": 0.5,
                    "min_through_hole_diameter": 0.3,
                    "min_hole_to_hole": 0.5,
                    "min_hole_clearance": 0.25,
                    "min_copper_edge_clearance": 0.3,
                    "min_silk_clearance": 0.0,
                    "min_text_height": 0.8,
                    "min_text_thickness": 0.15,
                    "solder_mask_to_copper_clearance": 0.0,
                },
                "rule_severities": {
                    "silk_over_copper": "ignore",
                    "silk_overlap": "ignore",
                    "silk_edge_clearance": "ignore",
                    "lib_footprint_issues": "ignore",
                    "lib_footprint_mismatch": "ignore",
                    "missing_courtyard": "ignore",
                    "footprint_type_mismatch": "ignore",
                    "text_height": "warning",
                },
            }
        },
        "net_settings": {
            "classes": [{
                "name": "Default", "clearance": 0.2, "track_width": 0.25,
                "via_diameter": VIA_D, "via_drill": VIA_DRILL,
                "microvia_diameter": 0.3, "microvia_drill": 0.1,
                "diff_pair_width": 0.2, "diff_pair_gap": 0.25, "diff_pair_via_gap": 0.25,
                "wire_width": 6, "bus_width": 12, "line_style": 0, "pcb_color": "rgba(0, 0, 0, 0.000)",
                "schematic_color": "rgba(0, 0, 0, 0.000)", "priority": 2147483647,
            }],
            "meta": {"version": 3},
        },
        "meta": {"filename": os.path.basename(path), "version": 1},
    }
    with open(path, "w") as fh:
        json.dump(pro, fh, indent=2)


# ---------------------------------------------------------------------------
# Phasing board
# ---------------------------------------------------------------------------

def pad(num, at, size, layers, net=None, kind="smd", shape="rect", drill=None):
    return dict(num=num, at=at, size=size, layers=layers, net=net, kind=kind, shape=shape, drill=drill)


def build_phasing_board(freq_mhz, outdir):
    dL = C0 / (4 * freq_mhz * 1e6 * math.sqrt(E_EFF)) * 1000  # mm
    W2 = 15.0                       # half width (board is 30 mm wide)
    LB = 99.0 - TAB_L               # body length; body + tabs = 99 mm (JLC <100 mm pricing)
    b = Board(f"phasing_{freq_mhz}", True, (50.0, 100.0))

    # ---- outline -----------------------------------------------------------
    NOTCH_W, NOTCH_D = 1.0, 0.8     # shoulder relief so the tab seats flush
    ch = 1.0                        # 45 deg chamfer on the outer corners
    outline = [(0, -W2 + ch), (ch, -W2), (LB - ch, -W2), (LB, -W2 + ch)]
    for yt in (-TAB_Y, TAB_Y):
        y0, y1 = yt - TAB_W / 2, yt + TAB_W / 2
        outline += [(LB, y0 - NOTCH_W), (LB - NOTCH_D, y0 - NOTCH_W), (LB - NOTCH_D, y0),
                    (LB + TAB_L, y0), (LB + TAB_L, y1), (LB - NOTCH_D, y1),
                    (LB - NOTCH_D, y1 + NOTCH_W), (LB, y1 + NOTCH_W)]
    outline += [(LB, W2 - ch), (LB - ch, W2), (ch, W2), (0, W2 - ch)]
    b.poly_outline(outline)

    # ---- feed connector (edge launch SMA / BNC / F) ------------------------
    # centre pin lands on B.Cu; ground legs clamp both faces.
    j_pads = [
        pad("1", (2.9, 0), (5.2, 1.2), ["B.Cu", "B.Mask"], "RF_IN"),
        pad("2", (3.3, 0), (6.0, 13.0), ["F.Cu", "F.Mask"], "GND"),
        pad("2", (3.3, -4.05), (6.0, 4.9), ["B.Cu", "B.Mask"], "GND"),
        pad("2", (3.3, 4.05), (6.0, 4.9), ["B.Cu", "B.Mask"], "GND"),
    ]
    b.footprint("J1", "EdgeLaunch_SMA_BNC_F_1.6mm", (0, 0), j_pads, smd=True,
                descr="Universal edge-launch pads for 1.6mm PCB: SMA, BNC or Type-F edge mount jack")
    rf_rects_B = [(0.3, -0.6, 5.5, 0.6)]
    gnd_pad_rects = [(0.3, -6.5, 6.3, 6.5)]

    # ---- CMCC (balun): RG-316 U-turn through BN-43-3312 --------------------
    HX = 14.0
    H_IN, H_OUT = (HX, 5.5), (HX, -5.5)
    HOLE = 1.8

    def shield_rect(h):   # F.Cu, braid pigtail lies here (toward the feed connector)
        return (h[0] - 6.4, h[1] - 1.5, h[0] - 1.4, h[1] + 1.5)

    def center_rect(h):   # B.Cu, centre conductor is bent over onto this
        return (h[0] - 4.4, h[1] - 0.8, h[0] - 1.3, h[1] + 0.8)

    def rpad(num, r, layers, net):
        return pad(num, ((r[0] + r[2]) / 2 - HX, (r[1] + r[3]) / 2), (r[2] - r[0], r[3] - r[1]), layers, net)

    core = (HX + 6.0, -9.7, HX + 6.0 + 25.4, 9.7)   # BN-43-3312 25.4 x 19.4 x 9.5 mm
    cm_pads = [
        rpad("1", center_rect(H_IN), ["B.Cu", "B.Mask"], "RF_IN"),
        rpad("2", shield_rect(H_IN), ["F.Cu", "F.Mask"], "GND"),
        rpad("3", center_rect(H_OUT), ["B.Cu", "B.Mask"], "LINE_N"),
        rpad("4", shield_rect(H_OUT), ["F.Cu", "F.Mask"], "LINE_P"),
        pad("", (0, H_IN[1]), (HOLE, HOLE), ["*.Cu", "*.Mask"], None, "np_thru_hole", "circle", HOLE),
        pad("", (0, H_OUT[1]), (HOLE, HOLE), ["*.Cu", "*.Mask"], None, "np_thru_hole", "circle", HOLE),
    ]
    b.footprint("L1", "CMCC_RG316_BN-43-3312", (HX, 0), cm_pads,
                descr="Common mode current choke: RG-316 U-turn through BN-43-3312, NPTH feed-through holes")

    # outer-layer copper that must be kept free of GND stitching vias
    V_P, V_N = (10.8, -9.0), (10.8, -10.4)
    J = (12.0, -9.7)
    b.track([(10.8, -7.2), V_P], 0.8, "F.Cu", "LINE_P")
    b.track([(9.9, -5.5), (9.9, -10.4), V_N], 0.4, "B.Cu", "LINE_N")
    b.track([(5.4, 0), (7.6, 0), (9.6, 2.0), (9.6, 4.9), (10.4, 5.5)], 0.36, "B.Cu", "RF_IN")
    b.via(V_P, "LINE_P")
    b.via(V_N, "LINE_N")
    nongnd_outer = [shield_rect(H_OUT), center_rect(H_OUT), center_rect(H_IN)] + rf_rects_B
    nongnd_tracks = [((10.8, -7.2), V_P, 0.8), ((9.9, -5.5), (9.9, -10.4), 0.4), ((9.9, -10.4), V_N, 0.4),
                     ((5.4, 0), (7.6, 0), 0.36), ((7.6, 0), (9.6, 2.0), 0.36), ((9.6, 2.0), (9.6, 4.9), 0.36)]
    sig_vias = [V_P, V_N]
    gnd_pad_rects += [shield_rect(H_IN)]

    # ---- tabs --------------------------------------------------------------
    XE = LB - 4.4                    # where the stacked pair fans out to the tab vias
    tab_ends = {}
    for name, yt in (("TAB1", -TAB_Y), ("TAB2", TAB_Y)):
        x0, x1 = LB - 2.2, LB + TAB_L - 0.3
        cx, L = (x0 + x1) / 2, x1 - x0
        pads_ = []
        for k, dy in enumerate((-TAB_PAD_PITCH, 0, TAB_PAD_PITCH)):
            netF = "LINE_P" if dy == 0 else "GND"
            netB = "LINE_N" if dy == 0 else "GND"
            pads_.append(pad(str(k + 1), (cx - LB, dy), (L, TAB_PAD_W), ["F.Cu", "F.Mask"], netF))
            pads_.append(pad(str(k + 4), (cx - LB, dy), (L, TAB_PAD_W), ["B.Cu", "B.Mask"], netB))
            r = (x0, yt + dy - TAB_PAD_W / 2, x1, yt + dy + TAB_PAD_W / 2)
            if dy == 0:
                nongnd_outer.append(r)
            else:
                gnd_pad_rects.append(r)
        b.footprint("P1" if name == "TAB1" else "P2", f"PhasingTab_{name}", (LB, yt), pads_, smd=True,
                    descr="Edge tab: 3 pads per face, middle = balanced line, outer = shield")
        tvp, tvn = (LB - 3.3, yt - 0.8), (LB - 3.3, yt + 0.8)
        b.via(tvp, "LINE_P")
        b.via(tvn, "LINE_N")
        b.track([tvp, (LB - 1.6, yt - 0.3)], 0.5, "F.Cu", "LINE_P")
        b.track([tvn, (LB - 1.6, yt + 0.3)], 0.5, "B.Cu", "LINE_N")
        nongnd_tracks += [(tvp, (LB - 1.6, yt - 0.3), 0.5), (tvn, (LB - 1.6, yt + 0.3), 0.5)]
        sig_vias += [tvp, tvn]
        tab_ends[name] = ((XE, yt), tvp, tvn)

    # ---- phasing line routing --------------------------------------------
    Y1 = -11.5                       # lane of the direct path
    Y2 = 11.5                        # lane of the delayed path
    YT_MIN = Y1 + LANE + 1.0         # deepest a meander U may reach
    XU0, XU_MAX = 22.0, XE - 6.5 - LANE - 1.0

    path1 = [J, (J[0] + 1.8, Y1), (XE - (tab_ends["TAB1"][0][1] - Y1), Y1), tab_ends["TAB1"][0]]
    path1 = chamfer(path1)

    def make_path2(n_u, yt):
        pts = [J, (13.8, -7.9), (16.0, -7.9), (17.5, -6.4), (17.5, Y2 - 1.5), (19.0, Y2)]
        for i in range(n_u):
            x = XU0 + 2 * LANE * i
            pts += [(x, Y2), (x, yt), (x + LANE, yt), (x + LANE, Y2)]
        pts += [(XE - (Y2 - tab_ends["TAB2"][0][1]), Y2), tab_ends["TAB2"][0]]
        return chamfer(pts)

    L1 = polylen(path1)
    n_u = 0
    while True:
        base = polylen(make_path2(n_u, Y2 - 2.5)) if n_u else polylen(make_path2(0, 0))
        need = L1 + dL
        if n_u and polylen(make_path2(n_u, YT_MIN)) >= need:
            lo, hi = YT_MIN, Y2 - 2.5
            for _ in range(60):
                mid = (lo + hi) / 2
                if polylen(make_path2(n_u, mid)) > need:
                    lo = mid
                else:
                    hi = mid
            yt = (lo + hi) / 2
            break
        if n_u == 0 and base >= need:
            raise SystemExit("delayed path already too long; frequency too high for this layout")
        n_u += 1
        if XU0 + 2 * LANE * (n_u - 1) + LANE > XU_MAX:
            raise SystemExit(f"{freq_mhz} MHz needs more meander than fits on a {LB} mm board")
    path2 = make_path2(n_u, yt)
    L2 = polylen(path2)

    # stacked broadside pair: identical geometry on In1 (P) and In2 (N)
    for pth in (path1, path2):
        b.track(pth, TRACE_W, "In1.Cu", "LINE_P")
        b.track(pth, TRACE_W, "In2.Cu", "LINE_N")
    b.track([V_P, J], TRACE_W, "In1.Cu", "LINE_P")
    b.track([V_N, J], TRACE_W, "In2.Cu", "LINE_N")
    for name in ("TAB1", "TAB2"):
        e, tvp, tvn = tab_ends[name]
        b.track([e, tvp], TRACE_W, "In1.Cu", "LINE_P")
        b.track([e, tvn], TRACE_W, "In2.Cu", "LINE_N")
    inner_segs = []
    for pth in (path1, path2, [V_P, J], [V_N, J]) + tuple([tab_ends[n][0], tab_ends[n][k]] for n in ("TAB1", "TAB2") for k in (1, 2)):
        inner_segs += list(zip(pth, pth[1:]))

    # ---- ground stitching / fences ----------------------------------------
    cands = []
    for pth in (path1, path2):
        cands += offset_polyline_samples(pth, FENCE_OFFSET, FENCE_PITCH)
    # perimeter rows + field grid
    for yy in (-W2 + 0.9, W2 - 0.9):
        cands += [(x, yy) for x in frange(1.0, LB - 0.9, 2.0)]
    cands += [(0.9 + 0.0, y) for y in frange(-W2 + 1, W2 - 1, 2.0)]
    cands += [(x, y) for x in frange(1.5, LB - 1.0, 2.5) for y in frange(-W2 + 2.4, W2 - 2.4, 2.5)]
    cands += [(7.0, y) for y in (-6.0, -3.5, 3.5, 6.0)]  # tie the connector pads together

    holes = [H_IN, H_OUT]
    vias = []
    for c in cands:
        if not point_in_poly(c, outline) or poly_edge_dist(c, outline) < 0.9:
            continue
        if any(seg_dist(c, a, bb) < TRACE_W / 2 + INNER_GND_CLEARANCE + 0.3 for a, bb in inner_segs):
            continue
        if any(math.dist(c, h) < HOLE / 2 + 0.3 + 0.6 for h in holes):
            continue
        if any(math.dist(c, v) < VIA_D + 0.5 for v in sig_vias):
            continue
        if any(rect_dist(c, r) < VIA_D / 2 + 0.35 for r in nongnd_outer):
            continue
        if any(seg_dist(c, a, bb) < w / 2 + VIA_D / 2 + 0.35 for a, bb, w in nongnd_tracks):
            continue
        if any(rect_dist(c, r) < VIA_D / 2 + 0.2 for r in gnd_pad_rects):
            continue
        if any(math.dist(c, v) < 1.3 for v in vias):
            continue
        vias.append(c)
    for v in vias:
        b.via(v, "GND")

    # ---- ground pours on all four layers -----------------------------------
    for layer, clr in (("F.Cu", 0.3), ("In1.Cu", INNER_GND_CLEARANCE), ("In2.Cu", INNER_GND_CLEARANCE), ("B.Cu", 0.3)):
        b.zone(layer, "GND", outline, clr)

    # ---- silkscreen --------------------------------------------------------
    for layer in ("F.SilkS", "B.SilkS"):
        # polarization edge labels (read with the tab end on the right)
        b.text("RHCP", (LB - 17, -W2 + 2.0), layer, 1.5, bold=True)
        b.arrow((LB - 13.8, -W2 + 2.0), (LB - 8.0, -W2 + 2.0), layer, 0.25)
        b.text("LHCP", (LB - 17, W2 - 2.0), layer, 1.5, angle=180, bold=True)
        b.arrow((LB - 13.8, W2 - 2.0), (LB - 8.0, W2 - 2.0), layer, 0.25)
    # F side: ferrite + shields
    b.rect(*core, "F.SilkS", 0.15)
    b.text("BN-43-3312", ((core[0] + core[2]) / 2, -3.6), "F.SilkS", 1.0)
    b.text("RG-316 x1 turn", ((core[0] + core[2]) / 2, -1.6), "F.SilkS", 1.0)
    b.text("zip tie", ((core[0] + core[2]) / 2, 1.6), "F.SilkS", 0.9)
    b.line(((core[0] + core[2]) / 2 - 3, 3.0), ((core[0] + core[2]) / 2 + 3, 3.0), "F.SilkS", 0.15)
    for x in ((core[0] + core[2]) / 2 - 2.5, (core[0] + core[2]) / 2 + 2.5):
        b.line((x, -W2 + 0.3), (x, -W2 + 1.3), "F.SilkS", 0.2)
        b.line((x, W2 - 1.3), (x, W2 - 0.3), "F.SilkS", 0.2)
    b.text("SHIELD", (HX - 3.9, 0), "F.SilkS", 0.9, angle=90)
    b.text("IN", (HX - 3.9, 8.3), "F.SilkS", 0.9)
    b.text("OUT", (8.2, -8.4), "F.SilkS", 0.9)
    b.text("CMCC", (HX + 2.3, 0), "F.SilkS", 0.9, angle=90)
    b.text(f"Eggbeater Phasing Board  {freq_mhz:g} MHz", (66, -5.6), "F.SilkS", 1.3, bold=True)
    b.text("100R balanced shielded 1/4 wave phasing line", (66, -3.4), "F.SilkS", 0.9)
    b.text("github.com/parker-research/Eggbeater-Antenna-PCB", (66, 3.0), "F.SilkS", 0.8)
    b.text("rev A  (inspired by Halibut Electronics EggNOGS)", (66, 4.8), "F.SilkS", 0.8)
    # B side: centre pads
    b.text("CENTER", (HX - 2.85, 2.9), "B.SilkS", 0.8)
    b.text("CENTER", (HX - 2.85, -2.9), "B.SilkS", 0.8)
    b.text("feed: SMA / BNC / F", (9.0, -12.9), "B.SilkS", 0.8)
    b.text(f"{freq_mhz:g} MHz  |  4L 1.6mm JLC04161H-7628", (66, -3.0), "B.SilkS", 1.0)
    b.text("In1/In2: 0.25mm broadside pair", (66, -1.0), "B.SilkS", 0.8)

    os.makedirs(outdir, exist_ok=True)
    pcb = os.path.join(outdir, "phasing_board.kicad_pcb")
    b.write(pcb, f"Eggbeater Phasing Board {freq_mhz:g} MHz", "A",
            comments=[f"4-layer 1.6mm, JLC04161H-7628. Line: 0.25mm broadside pair In1/In2 = 100R diff, eeff {E_EFF}",
                      f"Direct path {L1:.2f} mm, delayed path {L2:.2f} mm, delta {L2 - L1:.2f} mm (target {dL:.2f} mm = 90 deg)"])
    write_project(os.path.join(outdir, "phasing_board.kicad_pro"))
    return dict(freq=freq_mhz, L1=L1, L2=L2, dL=dL, n_u=n_u, yt=yt, vias=len(vias))


def frange(a, b, s):
    out = []
    x = a
    while x <= b + 1e-9:
        out.append(round(x, 4))
        x += s
    return out


# ---------------------------------------------------------------------------
# Main board
# ---------------------------------------------------------------------------

def build_main_board(outdir):
    R = 29.0              # 58 mm diameter: fits inside a 2in (60.3mm) Sch40 PVC cap
    RB = 23.0             # bolt circle radius
    BOLT_DRILL, BOLT_PAD = 4.3, 9.0   # No.8 (4.17mm) screw clearance
    b = Board("main", False, (100.0, 100.0))
    b.circle((0, 0), R, "Edge.Cuts", 0.1)

    # slots (stadium shaped, non plated)
    for yc in (-TAB_Y, TAB_Y):
        r = SLOT_W / 2
        y0, y1 = yc - SLOT_L / 2 + r, yc + SLOT_L / 2 - r
        b.line((-r, y0), (-r, y1), "Edge.Cuts", 0.1)
        b.line((r, y0), (r, y1), "Edge.Cuts", 0.1)
        b.arc((-r, y0), (0, y0 - r), (r, y0), "Edge.Cuts", 0.1)
        b.arc((-r, y1), (0, y1 + r), (r, y1), "Edge.Cuts", 0.1)

    # bolts; world view from the UP side, KiCad Y is down.
    k = RB / math.sqrt(2)
    bolts = {"UL": (-k, -k), "UR": (k, -k), "LL": (-k, k), "LR": (k, k)}
    # loop A (fed by TAB1 / direct path when RHCP selected): UL + LR
    # loop B: UR + LL
    bolt_net = {"UL": "LOOP_A1", "LR": "LOOP_A2", "LL": "LOOP_B1", "UR": "LOOP_B2"}
    for i, (nm, p) in enumerate(bolts.items()):
        b.footprint(f"H{i + 1}", "MountingHole_4.3mm_No8_Pad9mm", p,
                    [pad("1", (0, 0), (BOLT_PAD, BOLT_PAD), ["*.Cu", "*.Mask"], bolt_net[nm],
                         "thru_hole", "circle", BOLT_DRILL)],
                    descr="No.8-32 loop terminal bolt, plated")

    # slot pads: 3 per side, plated through
    PX = SLOT_W / 2 + 0.3 + 1.0      # pad centre (pads are 2.0 long in X)
    slot_nets = {-TAB_Y: ("LOOP_A1", "LOOP_A2"), TAB_Y: ("LOOP_B1", "LOOP_B2")}
    for j, yc in enumerate((-TAB_Y, TAB_Y)):
        pads_ = []
        nL, nR = slot_nets[yc]
        for i, dy in enumerate((-TAB_PAD_PITCH, 0, TAB_PAD_PITCH)):
            pads_.append(pad(str(i + 1), (-PX, dy), (2.0, TAB_PAD_W), ["*.Cu", "*.Mask"],
                             nL if dy == 0 else None, "thru_hole", "rect", 0.5))
            pads_.append(pad(str(i + 4), (PX, dy), (2.0, TAB_PAD_W), ["*.Cu", "*.Mask"],
                             nR if dy == 0 else None, "thru_hole", "rect", 0.5))
        b.footprint(f"S{j + 1}", "PhasingBoard_Slot_1.8x8.4mm", (0, yc), pads_,
                    descr="Slot + solder pads for phasing board tab")

    # traces: loop A on F.Cu, loop B on B.Cu (they cross)
    TW = 1.5
    b.track([(-PX, -TAB_Y), (-5.5, -TAB_Y), bolts["UL"]], TW, "F.Cu", "LOOP_A1")
    b.track([(PX, -TAB_Y), (5.5, -TAB_Y), bolts["LR"]], TW, "F.Cu", "LOOP_A2")
    b.track([(-PX, TAB_Y), (-5.5, TAB_Y), bolts["LL"]], TW, "B.Cu", "LOOP_B1")
    b.track([(PX, TAB_Y), (5.5, TAB_Y), bolts["UR"]], TW, "B.Cu", "LOOP_B2")

    # --- silkscreen, UP side (F) ---
    F, B = "F.SilkS", "B.SilkS"
    b.text("Eggbeater Main Board", (0, -24.0), F, 1.4, bold=True)
    b.text("THIS SIDE UP", (0, -21.6), F, 1.6, bold=True)
    b.text("Phasing board on OTHER side", (0, -12.7), F, 0.9)
    b.text("solder tabs here", (0, -11.0), F, 0.9)
    for yc in (-TAB_Y, TAB_Y):
        b.rect(-PX - 1.25, yc - 3.55, PX + 1.25, yc + 3.55, F, 0.15)
    for nm, (bx, by) in bolts.items():
        b.text("A" if bolt_net[nm].startswith("LOOP_A") else "B", (bx, by + (7.0 if by < 0 else -7.0)), F, 2.0, bold=True)
    b.text("Loop A: bolts A-A   Loop B: bolts B-B", (0, 12.0), F, 0.9)
    b.text("2in PVC cap  |  No.8-32 hardware", (0, 22.0), F, 1.0)
    b.text("rev A", (0, 24.0), F, 1.0)

    # --- silkscreen, DOWN side (B) ---
    b.text("THIS SIDE DOWN", (0, 21.6), B, 1.6, bold=True)
    b.text("insert phasing board here", (0, 12.7), B, 1.0)
    b.text("(solder here too)", (0, 14.5), B, 1.0)
    # polarization selector: the phasing board edge that sits here is the polarization you get
    b.rect(-SLOT_W / 2, -17.0, SLOT_W / 2, -12.0, B, 0.15)
    b.arrow((-5.0, -14.5), (-1.6, -14.5), B, 0.25)
    b.arrow((5.0, -14.5), (1.6, -14.5), B, 0.25)
    b.text("SELECT", (0, -22.4), B, 1.4, bold=True)
    b.text("POLARIZATION", (0, -20.3), B, 1.4, bold=True)
    b.text("edge label here = result", (0, -18.2), B, 0.8)

    os.makedirs(outdir, exist_ok=True)
    b.write(os.path.join(outdir, "main_board.kicad_pcb"), "Eggbeater Main Board", "A",
            comments=["2-layer 1.6mm FR4. Mates with phasing board at 90 degrees.",
                      "Loop A = bolts UL+LR (F.Cu), Loop B = bolts UR+LL (B.Cu)"])
    write_project(os.path.join(outdir, "main_board.kicad_pro"))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--freq-mhz", type=float, default=436.5)
    ap.add_argument("--outdir", default=None, help="phasing board output directory")
    a = ap.parse_args()
    od = a.outdir or os.path.join(ROOT, "kicad", "phasing_board")
    info = build_phasing_board(a.freq_mhz, od)
    build_main_board(os.path.join(ROOT, "kicad", "main_board"))
    print(json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in info.items()}, indent=1))
