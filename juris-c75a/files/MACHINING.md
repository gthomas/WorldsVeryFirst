# Juris C-75A — machining guide

This folder holds the CAM-ready files for every machined part. It does not hold toolpaths or G-code: those depend on your machine, controller and CAM program, so you generate them from these models. Every part has been checked to be cuttable on a **3-axis router in one or two setups** (top, then flipped), with no undercuts.

## Assumptions

- **Machine:** 3-axis CNC router with at least 610 × 460 mm (24 × 18 in) of travel and about 100 mm (4 in) of Z clearance.
- **Smallest tool for pockets and cut-outs:** 1/4 in (6.35 mm) end mill. Every inside corner on a pocket or cut-out is at least 3.5 mm radius, so a 1/4 in tool reaches it.
- **Units:** mm throughout.
- **Origin for every model:** X/Y at the centre of the part's bounding box. Z0 is at the top of the part as exported, or at the bottom for the fingerboard and bridge. The setup DXF for each part draws the stock, the part outline and the origin.
- **Two-sided parts:** flip about the X axis, along the line through the two registration pins. The pin holes (1/4 in, 6.35 mm) are drawn on each setup DXF in the waste at both ends of the stock.

## Stock to buy

Standard sizes from a luthier supply. Finished part sizes are in the Parts table below.

| Qty | Item | Typical size | Used for | Notes |
|---|---|---|---|---|
| 2 | Mahogany body blank, 2-piece | 21 × 14 × 1¾–2 in | One rim half each | The half is 18.05 × 13.6 × 1.12 in finished. A 14 in blank leaves only about 0.2 in per side, so hold it by the ends (the pin holes in the setup DXF are at the ends), or buy oversized 15 in blanks. The centre seam is fine: it runs through the neck and tail blocks. |
| 1 | Mahogany neck blank | 30 × 3 × 1 in | Neck | Needs three glue-ups before machining: a 13° scarf-jointed headstock, a stacked heel (two 1 in layers glued under the heel end, about 3 in of total depth there), and headstock ears to reach 3.5 in width. Or buy a one-piece blank of at least 26 × 3.5 × 3 in and skip the glue-ups. |
| 1 | Quartersawn Sitka spruce archtop top wedges, bookmatched | 2 halves, each about 22 × 8 in, at least 1 in thick at the centre joint | Top | Joined it is about 20 × 16 in. The finished top is 0.95 in thick at the plateau, which sits on the centre line. Adirondack (red) spruce is the upgrade for more volume and stiffness, with no design changes. |
| 1 | Quartersawn bigleaf maple back set, bookmatched | 2 halves, each about 22 × 8 in, about 1/2 in thick | Back | Finished 0.36 in at most. Plain or flamed makes no acoustic difference. A stiff maple back reflects energy into the top and machines crisply. |
| 1 | Ebony fingerboard blank | 20½ × 2¾ × 3/8 in | Fingerboard | Finished 18 × 2.22 × 0.29 in. |
| 1 | Ebony bridge blank | 5 × 1 × 1½ in | Bridge | Finished 3.8 × 0.85 × 1.3 in. |
| 1 | Ebony or black fibre headstock veneer | 8 × 4 × 1/16 in | Faceplate | Finished head is 7.8 × 3.5 in. |
| 1 | Black pickguard sheet (acrylic or celluloid) | 8 × 5 × 0.09 in | Pickguard | |
| 1 | Brass sheet | 2 × 1 in, 1.5 mm | Fingerboard plate | |
| 4 | Hardwood dowels | 6 × 20 mm | Rim halves | |
| 1 | Dual-action truss rod | 15 in | Neck | Check the channel size against your rod. |

## Tools

| Tool | Used for |
|---|---|
| 1/2 in or 3/8 in flat end mill | Roughing the top, back, neck and rim pockets |
| 1/4 in (6.35 mm) flat end mill | Pockets, profiles, soundholes, mortise, truss rod channel (also its width) |
| 1/4 in ball-nose | 3D finishing of the top facets, back, neck carve, heel, volute, bridge |
| 1/8 in ball-nose | Finishing the volute and the heel-to-shaft blend |
| 0.6 mm (0.023 in) end mill or a fret-slot saw | Fret slots |
| 90° V-bit | Edge chamfer on the assembled body |
| 6 mm drill / end mill | Rim alignment dowels |

Feeds and speeds depend on your machine; use the CAM program's defaults for hardwood (mahogany, maple, ebony) and softwood (spruce) as a starting point. Run the grain lengthwise on the top and back, use compression or downcut end mills for profiles and a fresh ball-nose for the facets, and finish the spruce facets with a light 0.1–0.2 mm pass cutting with the grain where your CAM allows.

## Parts

| Part | Folder | Stock (L × W × T) | Setups | What gets cut |
|---|---|---|---|---|
| Top | `top/` | Bookmatched quartersawn Sitka spruce, 16 × 20 × 1-1/8 in | 2 | A: outside facets and flat skirt. B (flipped): inside of the plate, ring beam, two tone bars, upper-bout cross bar. Soundholes cut through. |
| Back | `back/` | Bookmatched bigleaf maple, 16 × 20 × 1/2 in | 2 | Exported inside-up. A: inside (flat glue skirt, arch underside). B: outside facets. |
| Rim, lower half | `rim_lower/` | Mahogany body blank, 21 × 14 × 1¾ in, planed to 1.12 in | 1 | From the split face: pocket the recess 20.55 mm deep, cut through inside the glue ledge, drill 4 dowel holes 10 mm deep, profile the outside. |
| Rim, upper half | `rim_upper/` | Same | 2 | Same as the lower half, then flip and cut the open-top neck mortise 26 mm deep from the outer face. You can also cut the mortise after the two halves are glued. |
| Neck | `neck/` | Glued-up neck blank (see Stock to buy), at least 26 × 3.5 × 3 in overall | 2 | Neck frame: X along the neck, Z0 at the fingerboard gluing plane. A: gluing plane, truss rod channel, headstock face. B (flipped): back carve, heel, tenon, volute, headstock back. |
| Fingerboard | `fingerboard/` | Ebony, 20½ × 2¾ × 3/8 in | 1 | Flat side down: 12 in radius top, 22 fret slots, taper profile. |
| Bridge | `bridge/` | Ebony, 5 × 1 × 1½ in | 2 | Base down: compensated, radiused saddle crown and sides. Flip: relief under the centre. |

The `flat/` folder has 2D cut files:
- `faceplate_and_truss_cover.dxf`: 1.6 mm ebony or black fibre headstock overlay with tuner holes and truss rod access, plus the truss rod cover and its screw holes.
- `pickguard_2.3mm.dxf`: flat pattern for the floating guard, laid out in the guard's own plane. It includes the two 9.5 mm pot holes, the three Lollar tab holes and the mounting holes.
- `fingerboard_plate_1.5mm_brass.dxf`: the small plate joining the guard to the underside of the fingerboard extension.

## Build order

1. **Rim:** glue the two halves on four 6 × 20 mm dowels. Check that the faces are flat and the ledges line up.
2. **Back:** glue the back onto the lower ledge.
3. **Truss rod:** fit a 15 in dual-action rod into the neck channel. It adjusts at the headstock, under the cover.
4. **Neck:** glue the neck tenon into the open-top mortise. Below the rim top the heel face is flat on the body's end face; the tenon floor sets the height.
5. **Top:** glue the top onto the upper ledge and neck block. It also covers the tenon.
6. **Chamfer:** cut the 45° edge chamfer around the top and back edges with a V-bit, 3.5 mm deep. It fades to nothing over the last 28 mm before the neck's end face; blend that taper by hand or with a ramped toolpath.
7. **Finish the top:** mask a 1/4 in (6.35 mm) band around the top's edge, including the chamfer, and dye the rest of the top black. The mask line is the `TOP_DYE_MASK` layer in `model/juris_c75a_templates.dxf`. Pull the mask, then shellac over everything; the undyed band reads as a faux binding.
8. **Fingerboard:** glue it on, fret it, and drill the side dots by hand.
9. **Hand-drilled holes:** drill the side jack and endpin square to the rim. Positions are in `positions.json`.
10. **Tuner holes:** drill them with the headstock face held flat, using a drill press with a 13° wedge or a separate CNC setup. They are 10 mm; change `tuner_hole` in `build.py` if your tuners are different.
11. **Bridge:** sand the saddle crown down to set the action. It is left 1.5 mm tall on purpose. Then file the string notches.

## Check against your hardware before cutting

These dimensions come from typical parts. Measure yours and change the parameter in `build.py`, then run `python build.py && python cnc.py` to regenerate everything.

| Part | Assumed | Parameter |
|---|---|---|
| Tuners | 10 mm post holes | `tuner_hole` |
| Truss rod | 15 in long, needs a 1/4 × 3/8 in channel | `truss_len`, `truss_w`, `truss_d` |
| Pickup | Lollar gold foil, floating mount (drawing dimensions) | `build_pickup` |
| Pots | 3/8 in bushings (9.5 mm holes) | pickguard DXF |
| Strings | 13–56 set for the action and bridge height | `action12` |

## Known limits

- **Outline accuracy:** the body outline was traced from the C-75 PDF at about 0.13 mm per pixel. If you have the original CAD outline, replace the traced edges in `outline.py`.
- **Curved surfaces:** the STEP models describe the neck carve, volute and fingerboard radius as many small flat faces, about 0.005 mm from the true shape. That's fine for 3D finishing, but they won't behave as single smooth surfaces in CAD.
- **Undercut check:** casts vertical rays on a 0.5 mm grid (0.25 mm for the fingerboard and bridge). One cell on the neck reads 0.2 mm²; recasting that spot gives a clean result, so it's numerical noise.
- **Hardware not modelled for machining:** the tailpiece, pickup, knobs and brackets are purchased or fabricated.
