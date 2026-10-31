"""Juris C-75A — CAM export.
Writes cnc/<part>/ with STEP + STL in the part's machining frame, a setup DXF (stock, part silhouette,
origin, registration pins), and runs a two-sided 3-axis reachability check on every machined part.
Units mm. Run after build.py:  python cnc.py
"""
import json, math, os, sys
import numpy as np, trimesh, ezdxf
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union
import manifold3d as m3
import build as B

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cnc')
os.makedirs(OUT, exist_ok=True)
MARGIN = 25.0        # stock margin around the part for clamps and pins
PIN_D = 6.35         # 1/4" registration pins
SKIP_STEP = '--no-step' in sys.argv


def tm_of(man):
    mh = man.simplify(0.005).to_mesh()
    tm = trimesh.Trimesh(np.array(mh.vert_properties)[:, :3].astype(float), np.array(mh.tri_verts), process=True)
    return tm


# ----------------------------------------------------------------------------- frames
def neck_to_machine():
    """4x4: body -> neck machining frame. X = u (toward headstock), Y = -v, Z = w + fb_t (Z=0 on the
    fingerboard gluing plane, neck carve toward -Z)."""
    o = B.n2b(0, 0, 0); eu = B.n2b(1, 0, 0) - o; ev = B.n2b(0, 1, 0) - o; ew = B.n2b(0, 0, 1) - o
    Rm = np.array([eu, -ev, ew])               # rows: machine axes in body coords (proper rotation)
    T = np.eye(4); T[:3, :3] = Rm; T[:3, 3] = -Rm @ o + np.array([0, 0, B.P['fb_t']])
    return T


def flip_x():
    T = np.eye(4); T[1, 1] = -1; T[2, 2] = -1; return T


def place(tm, T=None, z0='top'):
    """apply T, then centre XY on the bbox and put Z=0 at the top (or bottom) of the part"""
    tm = tm.copy()
    if T is not None: tm.apply_transform(T)
    lo, hi = tm.bounds
    off = np.array([-(lo[0] + hi[0]) / 2, -(lo[1] + hi[1]) / 2, -(hi[2] if z0 == 'top' else lo[2])])
    tm.apply_translation(off)
    return tm, off


# ----------------------------------------------------------------------------- reachability
def reach_check(tm, step=0.5):
    """Cast vertical rays on a grid. A part is machinable from the top and bottom (two setups, 3-axis)
    when every vertical line crosses its surface at most twice. Returns the area (mm^2) where more
    than two crossings occur."""
    lo, hi = tm.bounds
    xs = np.arange(lo[0] + step / 2, hi[0], step) + 0.0123
    ys = np.arange(lo[1] + step / 2, hi[1], step) + 0.0171
    X, Y = np.meshgrid(xs, ys)
    orig = np.c_[X.ravel(), Y.ravel(), np.full(X.size, hi[2] + 10)]
    dirs = np.tile([0, 0, -1.0], (len(orig), 1))
    loc, idx_ray, _ = tm.ray.intersects_location(orig, dirs, multiple_hits=True)
    # merge duplicate hits on shared edges
    cnt = np.zeros(len(orig), int)
    order = np.lexsort((loc[:, 2], idx_ray))
    loc, idx_ray = loc[order], idx_ray[order]
    keep = np.ones(len(idx_ray), bool)
    keep[1:] = ~((idx_ray[1:] == idx_ray[:-1]) & (np.abs(loc[1:, 2] - loc[:-1, 2]) < 1e-3))
    np.add.at(cnt, idx_ray[keep], 1)
    bad = cnt > 2
    pts = orig[bad, :2]
    return float(bad.sum() * step * step), pts, int((cnt > 0).sum())


# ----------------------------------------------------------------------------- STEP
def write_step(tm, path, unify=True):
    """triangle mesh -> sewn B-rep solid with coplanar triangles merged into true planar faces"""
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakePolygon, BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid
    from OCP.gp import gp_Pnt
    from OCP.TopoDS import TopoDS
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopAbs import TopAbs_SHELL
    from OCP.ShapeFix import ShapeFix_Shape
    from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain
    from OCP.STEPControl import STEPControl_Writer, STEPControl_AsIs
    from OCP.Interface import Interface_Static
    tm = tm.copy(); tm.merge_vertices(digits_vertex=4)
    tm.update_faces(tm.nondegenerate_faces(height=1e-5)); tm.remove_unreferenced_vertices()
    sew = BRepBuilderAPI_Sewing(1e-3)
    V = tm.vertices
    for f in tm.faces:
        poly = BRepBuilderAPI_MakePolygon(gp_Pnt(*V[f[0]]), gp_Pnt(*V[f[1]]), gp_Pnt(*V[f[2]]), True)
        mf = BRepBuilderAPI_MakeFace(poly.Wire(), True)
        if mf.IsDone(): sew.Add(mf.Face())
    sew.Perform()
    ex = TopExp_Explorer(sew.SewedShape(), TopAbs_SHELL)
    shells = []
    while ex.More(): shells.append(ex.Current()); ex.Next()
    assert len(shells) == 1 and sew.NbFreeEdges() == 0, (len(shells), sew.NbFreeEdges())
    solid = BRepBuilderAPI_MakeSolid(TopoDS.Shell(shells[0])).Solid()
    fx = ShapeFix_Shape(solid); fx.Perform()
    if unify:
        u = ShapeUpgrade_UnifySameDomain(fx.Shape(), True, True, True)
        u.SetLinearTolerance(1e-4); u.SetAngularTolerance(1e-4); u.Build()
        fx2 = ShapeFix_Shape(u.Shape()); fx2.Perform(); out = fx2.Shape()
    else:
        out = fx.Shape()
    from OCP.TopAbs import TopAbs_SOLID
    if not TopExp_Explorer(out, TopAbs_SOLID).More():          # fix returned a bare shell: wrap it again
        ex2 = TopExp_Explorer(out, TopAbs_SHELL)
        out = BRepBuilderAPI_MakeSolid(TopoDS.Shell(ex2.Current())).Solid()
    Interface_Static.SetCVal_s('write.step.unit', 'MM')
    w = STEPControl_Writer(); w.Transfer(out, STEPControl_AsIs); w.Write(path)


def check_step(path):
    from OCP.STEPControl import STEPControl_Reader
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.GProp import GProp_GProps
    from OCP.BRepGProp import BRepGProp
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopAbs import TopAbs_SOLID, TopAbs_FACE
    r = STEPControl_Reader(); r.ReadFile(path); r.TransferRoots(); s = r.OneShape()
    g = GProp_GProps(); BRepGProp.VolumeProperties_s(s, g)
    def n(k):
        ex = TopExp_Explorer(s, k); c = 0
        while ex.More(): c += 1; ex.Next()
        return c
    return dict(valid=bool(BRepCheck_Analyzer(s).IsValid()), solids=n(TopAbs_SOLID), faces=n(TopAbs_FACE), volume_cm3=round(g.Mass() / 1000, 2))


# ----------------------------------------------------------------------------- setup DXF
def silhouette(tm):
    polys = [Polygon(t[:, :2]) for t in tm.triangles if abs(Polygon(t[:, :2]).area) > 1e-6]
    sil = unary_union(polys).buffer(0.01).buffer(-0.01)
    return sil


def setup_dxf(path, tm, title, extra=None):
    lo, hi = tm.bounds
    L, W = hi[0] - lo[0] + 2 * MARGIN, hi[1] - lo[1] + 2 * MARGIN
    doc = ezdxf.new('R2010'); doc.units = ezdxf.units.MM; ms = doc.modelspace()
    for n, c in [('STOCK', 8), ('PART_OUTLINE', 7), ('PINS', 1), ('ORIGIN', 3), ('NOTES', 2)] + [(k, 5) for k in (extra or {})]:
        doc.layers.add(n, color=c)
    ms.add_lwpolyline([(-L / 2, -W / 2), (L / 2, -W / 2), (L / 2, W / 2), (-L / 2, W / 2)], close=True, dxfattribs={'layer': 'STOCK'})
    sil = silhouette(tm)
    for g in (sil.geoms if isinstance(sil, MultiPolygon) else [sil]):
        ms.add_lwpolyline(list(g.exterior.coords), close=True, dxfattribs={'layer': 'PART_OUTLINE'})
        for h in g.interiors:
            ms.add_lwpolyline(list(h.coords), close=True, dxfattribs={'layer': 'PART_OUTLINE'})
    for x in (-L / 2 + MARGIN / 2, L / 2 - MARGIN / 2):
        ms.add_circle((x, 0), PIN_D / 2, dxfattribs={'layer': 'PINS'})
    ms.add_line((-10, 0), (10, 0), dxfattribs={'layer': 'ORIGIN'}); ms.add_line((0, -10), (0, 10), dxfattribs={'layer': 'ORIGIN'})
    ms.add_text(title, height=5, dxfattribs={'layer': 'NOTES'}).set_placement((-L / 2, W / 2 + 4))
    for lay, geoms in (extra or {}).items():
        for g in geoms:
            ms.add_lwpolyline([tuple(map(float, p)) for p in g], close=True, dxfattribs={'layer': lay})
    doc.saveas(path)
    return L, W


# ----------------------------------------------------------------------------- parts
def main():
    B.MORTISE = B.mortise()
    lower, upper = B.build_rim_halves()
    TN = neck_to_machine()
    parts = {}
    # name: (manifold, transform body->machine before placing, z0, two_sided, notes)
    parts['top'] = (B.build_top(chamfer=False), None, 'top', True,
                    'Setup A outside (facets, skirt) from the top; setup B: flip about X (pin line), inside (beam, bars, plate) from the top. Soundholes through.')
    parts['back'] = (B.build_back(chamfer=False), flip_x(), 'top', True,
                     'Exported inside-up. Setup A inside (flat glue skirt and arch underside); setup B: flip about X, outside facets.')
    parts['rim_lower'] = (lower, None, 'top', False,
                          'One setup from the split face: pocket the recess to the ledge, cut through inside the ledge, drill dowels, profile.')
    parts['rim_upper'] = (upper, flip_x(), 'top', True,
                          'Exported split-face-up. Setup A same as the lower half. Setup B: flip about X and cut the open-top neck mortise from the outer face (or cut it after the halves are glued).')
    parts['neck'] = (B.build_neck(), TN, 'top', True,
                     'Neck frame: X along the neck, Z=0 at the fingerboard gluing plane. Setup A fingerboard face, truss channel, headstock face; setup B: flip about X, back carve, heel, volute, headstock back. Tuner holes: separate op with the headstock face held flat.')
    parts['fingerboard'] = (B.build_fingerboard(), TN, 'bottom', False,
                            'Flat side down. Radius the top (12"), cut 0.6 mm fret slots, profile the taper.')
    # the bridge is purchased (flat-bottomed base + compensated archtop top), so it is not machined

    report = {}
    for name, (man, T, z0, two, note) in parts.items():
        d = os.path.join(OUT, name); os.makedirs(d, exist_ok=True)
        tm = tm_of(man)
        tm, off = place(tm, T, z0)
        tm.export(os.path.join(d, f'{name}.stl'))
        if not SKIP_STEP:
            print('step', name, len(tm.faces), flush=True)
            write_step(trimesh.load(os.path.join(d, f'{name}.stl')), os.path.join(d, f'{name}.step'))
            stepchk = check_step(os.path.join(d, f'{name}.step'))
            if not (stepchk['valid'] and stepchk['solids'] == 1):     # freeform parts: keep the planar triangles as faces
                write_step(trimesh.load(os.path.join(d, f'{name}.stl')), os.path.join(d, f'{name}.step'), unify=False)
                stepchk = check_step(os.path.join(d, f'{name}.step'))
        else:
            stepchk = None
        area, badpts, hit = reach_check(tm, 0.5 if name not in ('fingerboard', 'bridge') else 0.25)
        lo, hi = tm.bounds
        size = (hi - lo)
        extra = None
        if name.startswith('rim'):
            through, recess = B.rim_regions()
            fy = -1 if name == 'rim_upper' else 1
            tr = lambda pts: [(x + off[0], fy * y + off[1]) for x, y in pts]
            circ = lambda c, r: [(c[0] + r * math.cos(t), c[1] + r * math.sin(t)) for t in np.linspace(0, 2 * math.pi, 64)[:-1]]
            extra = {'CUT_THROUGH_inside_ledge': [tr(p) for p in through.to_polygons()],
                     f'POCKET_recess_{B.P["rim_h"] / 2 - B.P["lining_h"]:.2f}mm_deep': [tr(p) for p in recess.to_polygons()],
                     'OUTER_PROFILE': [tr(B.OUTL)],
                     f'DOWEL_{B.P["dowel_d"]}mm_x{B.P["dowel_depth"]}mm_deep': [tr(circ(c, B.P['dowel_d'] / 2)) for c in B.DOWELS]}
            if name == 'rim_upper':
                extra[f'MORTISE_{B.P["tenon_depth"]}mm_deep_from_outer_face'] = [tr(B.tenon_poly(B.P['tenon_fit']))]
        if name in ('top', 'back'):
            tr = (lambda pts: [(x + off[0], y + off[1]) for x, y in pts]) if name == 'top' else (lambda pts: [(x + off[0], -y + off[1]) for x, y in pts])
            extra = {'OUTER_PROFILE': [tr(B.OUTL)]}
            if name == 'top':
                extra['SOUNDHOLES_through'] = [tr(p) for p in B.SLOTS]
        L, W = setup_dxf(os.path.join(d, f'{name}_setup.dxf'), tm, f'JURIS C-75A  {name.upper()}  origin: XY centre, Z0 = part {z0}', extra)
        report[name] = dict(part_size_mm=[round(float(x), 1) for x in size], stock_mm=[round(L), round(W), round(float(size[2]) + 2.0, 1)],
                            stock_in=[round(L / 25.4, 2), round(W / 25.4, 2), round((float(size[2]) + 2.0) / 25.4, 2)],
                            watertight=bool(tm.is_watertight), triangles=int(len(tm.faces)),
                            undercut_area_mm2=round(area, 1), undercut_sample=badpts[:5].round(1).tolist(),
                            step=stepchk, stl_volume_cm3=round(tm.volume / 1000, 2), setups=2 if two else 1, origin_offset=off.round(3).tolist(), notes=note)
        print(name, json.dumps({k: report[name][k] for k in ('part_size_mm', 'stock_in', 'watertight', 'undercut_area_mm2')}), flush=True)
    json.dump(report, open(os.path.join(OUT, 'cnc_report.json'), 'w'), indent=1)
    flat_parts()
    return report


# ----------------------------------------------------------------------------- 2D parts and positions
def flat_parts():
    d = os.path.join(OUT, 'flat'); os.makedirs(d, exist_ok=True)
    def doc_with(layers):
        doc = ezdxf.new('R2010'); doc.units = ezdxf.units.MM
        for n, c in layers: doc.layers.add(n, color=c)
        return doc, doc.modelspace()
    circ = lambda ms, c, dia, lay: ms.add_circle(tuple(map(float, c)), dia / 2, dxfattribs={'layer': lay})
    pl = lambda ms, pts, lay: ms.add_lwpolyline([tuple(map(float, p)) for p in pts], close=True, dxfattribs={'layer': lay})

    # faceplate (headstock-face coordinates: x across, y from the nut back face toward the tip)
    doc, ms = doc_with([('OUTLINE', 7), ('TUNER_HOLES', 1), ('TRUSS_ACCESS', 2), ('TRUSS_COVER', 3)])
    pl(ms, B.HEAD_POLY, 'OUTLINE')
    for v, s, _ in B.TUNER_LIST: circ(ms, (v, s), B.P['tuner_hole'], 'TUNER_HOLES')
    pl(ms, B.open_poly(B.TRUSS_ACCESS, 2.0), 'TRUSS_ACCESS')
    pl(ms, B.TRUSS_COVER, 'TRUSS_COVER')
    for s_ in (8.0, 38.0): circ(ms, (0, s_), 2.0, 'TRUSS_COVER')
    doc.saveas(os.path.join(d, 'faceplate_and_truss_cover.dxf'))

    # gold foil bracket flat pattern (strip along y from the fingerboard tab; shelf crosses at the far end)
    gi = json.load(open(os.path.join(B.OUT, 'report.json'))).get('goldfoil')
    if gi:
        fp = gi['flat_pattern']; g = B.GOLDFOIL; bb = B.GF_BRACKET
        w2, y1 = fp['width'] / 2, fp['tab']; y2 = y1 + fp['drop']; y3 = y2 + fp['arm']; sw = fp['shelf_W']; hs = fp['shelf_span'] / 2
        doc, ms = doc_with([('OUTLINE', 7), ('BEND_LINES', 1), ('HOLES', 3), ('NOTES', 8)])
        pl(ms, [(-w2, 0), (w2, 0), (w2, y3 - sw), (hs, y3 - sw), (hs, y3), (-hs, y3), (-hs, y3 - sw), (-w2, y3 - sw)], 'OUTLINE')
        for yb in (y1, y2):
            ms.add_line((-w2, yb), (w2, yb), dxfattribs={'layer': 'BEND_LINES'})
        for x in bb['fb_screws']: circ(ms, (x, y1 / 2), 3.0, 'HOLES')
        for sx in (-1, 1): circ(ms, (sx * g['screw_x'], y3 - sw / 2), bb['screw_d'], 'HOLES')
        ms.add_text(f"{fp['t']} mm brass or aluminium. Bend 90 deg down at the first line, 90 deg back at the second. "
                    "Lengths are to the inside of the bends; add your bend allowance. Gold foil sizes are placeholders.",
                    dxfattribs={'layer': 'NOTES', 'height': 2.5}).set_placement((-hs, y3 + 6))
        doc.saveas(os.path.join(d, 'goldfoil_bracket_flat.dxf'))

    # rim control plate (plate coordinates: x along the rim face, y up; origin at the plate centre)
    doc, ms = doc_with([('OUTLINE', 7), ('POT_AND_JACK_HOLES', 1), ('SCREW_HOLES', 3), ('RIM_WINDOW_behind', 8)])
    L, H, T = B.CTRL_PLATE
    pl(ms, [(-L / 2, -H / 2), (L / 2, -H / 2), (L / 2, H / 2), (-L / 2, H / 2)], 'OUTLINE')
    for name, x, dia in B.CTRL_HOLES: circ(ms, (x, 0), dia, 'POT_AND_JACK_HOLES')
    dx, dz, dl, dh = B.CTRL_DOOR
    doc.layers.add('BATTERY_DOOR_OPENING', color=2)
    pl(ms, B.open_poly(np.array([(dx - dl / 2, dz - dh / 2), (dx + dl / 2, dz - dh / 2), (dx + dl / 2, dz + dh / 2), (dx - dl / 2, dz + dh / 2)]), 2.0), 'BATTERY_DOOR_OPENING')
    for se, sz in B.CTRL_SCREWS: circ(ms, (se, sz), 2.8, 'SCREW_HOLES')
    W, Hw = B.CTRL_WIN
    pl(ms, [(-W / 2, -Hw / 2), (W / 2, -Hw / 2), (W / 2, Hw / 2), (-W / 2, Hw / 2)], 'RIM_WINDOW_behind')
    doc.saveas(os.path.join(d, f'control_plate_{T:.0f}mm.dxf'))

    # positions table
    fb_end = B.n2b(B.U_FB_END, 0, 0)
    pos = dict(
        fret_positions_from_nut_mm=[round(B.fret_s(n), 2) for n in range(1, B.P['frets'] + 1)],
        side_dots_between_frets=[3, 5, 7, 9, 12, 15, 17, 19, 21], side_dot_dia_mm=2.0,
        tuner_holes_face_coords_mm=[[round(v, 1), round(s, 1)] for v, s, _ in B.TUNER_LIST], tuner_hole_dia_mm=B.P['tuner_hole'],
        control_plate=dict(center_body_xy=[round(float(c), 1) for c in B.CTRL_C], z=round(B.Z_SPLIT, 1),
                           outward_normal=[round(float(c), 4) for c in B.CTRL_N], along_rim=[round(float(c), 4) for c in B.CTRL_E],
                           holes={n: dict(along_mm=x, dia_mm=dd) for n, x, dd in B.CTRL_HOLES}, battery_door=dict(zip(('along_mm', 'height_mm', 'length_mm', 'tall_mm'), B.CTRL_DOOR)), screws_mm=B.CTRL_SCREWS,
                           note='the window behind it is machined; drill the four 2 mm pilot holes square to the rim face'),
        pickup_wire_hole=dict(center_body_xy=list(B.WIRE_HOLE[:2]), dia_mm=B.WIRE_HOLE[2], note='through the top, under the fingerboard extension (machined with the top)'),
        endpin=dict(center_body_xy=[0.0, 0.0], z=round(B.Z_BACK_IN + B.P['rim_h'] / 2 - 4, 1), dia_mm=12.0),
        dowels_body_xy=B.DOWELS, dowel=f"{B.P['dowel_d']} mm x {2 * B.P['dowel_depth']} mm",
        truss_rod=f"{B.P['truss_len'] / 25.4:.0f} in dual-action, channel {B.P['truss_w']} x {B.P['truss_d']} mm",
        bridge_saddle=dict(scale_mm=B.P['scale'], compensation_mm=dict(treble=1.5, bass=3.5), extra_height_mm=B.BRIDGE_FIT),
        body_origin='body frame: x across (treble +), y from the tail toward the neck, z up from the back edge',
    )
    json.dump(pos, open(os.path.join(OUT, 'positions.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
