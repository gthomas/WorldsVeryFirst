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
| 1 | Quartersawn Sitka spruce archtop top wedges, bookmatched | 2 halves, each about 22 × 8 in, at least 1 in thick at the centre joint | Top | Joined it is about 20 × 16 in. The finished top is 0.75 in thick at the plateau, which sits on the centre line. Adirondack (red) spruce is the upgrade for more volume and stiffness, with no design changes. |
| 1 | Quartersawn bigleaf maple back set, bookmatched | 2 halves, each about 22 × 8 in, about 3/4 in thick | Back | Finished 0.6 in at most. Plain or flamed makes no acoustic difference. A stiff maple back reflects energy into the top and machines crisply. |
| 1 | Richlite guitar fingerboard blank (Black Diamond), unslotted | 20½ × 2¾ × 3/8 in or larger | Fingerboard | Finished 18 × 2.22 × 0.29 in. Barely moves with humidity, so fret ends stay put in dry weather. Machines like a hard hardwood; cut fret slots with a sharp saw or slotting cutter and press the frets in (a drop of thin glue in each slot is a good idea). |
| 1 | Maple headstock veneer | 8 × 4 × 1/16 in | Faceplate | Finished head is 7.8 × 3.5 in. Finished to match the top (dye and shellac along with it). |
| 1 | Black fibre or aluminium sheet | 5 × 2½ in, 2 mm | Rim control plate | Holds the stacked pot, the stereo jack and the battery door. |
| 1 | Brass or steel strap | 4 × ½ in, 1.5 mm | Pickup bracket | Bent to hang the EMG 91 from the underside of the fingerboard extension, unless you use a fingerboard-mount bracket for it. |
| 1 | Stacked (concentric) pot, 25k, with a concentric knob | 3/8 in bushing | Volume and tone | EMG actives use 25k pots. |
| 1 | End-loading 9 V sleeve | 32 × 22 mm, 54 mm deep | Battery | Fold from 1.5 mm ABS or aluminium, or 3D-print it; it fixes to the back of the plate around the door opening. The door is cut from the plate offcut and held by a quarter-turn latch. |
| 4 | Hardwood dowels | 6 × 20 mm | Rim halves | |
| 1 | Dual-action truss rod | 15 in | Neck | Check the channel size against your rod. |

## Hardware to buy

| Item | Choice | Fit notes |
|---|---|---|
| Tuners | Genuine Grover Deluxe SD90, 3+3, nickel, safety post (133N), as on the C-75 | 15:1, 6 mm posts. 11/32 in post holes with 8.8 mm press-in bushings; one back screw per tuner. |
| Bridge base | Gretsch 008-0636-200 rosewood bridge base, flat bottom (for flat-top guitars), as on the C-75 | Sits flat on the bridge plateau with no fitting; thumbwheel posts at 74 mm. Measured base: 135 × 14.95 mm, 11.5–12 mm tall. The plateau is 141 mm (5.55 in) wide across the strings, so the base sits fully on it with 3 mm to spare each side, and its feet land over the 5 × 6 mm ring beam. Strength check with this base: about 0.62 mm sink under the feet, 0.87 mm at most, stresses at least 4× under clear-spruce strength. Expect it to settle a little more over the first months; take that up with the thumbwheels (checked without the soundhole lips, which only add stiffness). |
| Bridge top | A standard wooden compensated archtop bridge top (ebony or rosewood) with post holes at 74 mm | Must match the base's 74 mm post spacing. |
| Tailpiece | A trapeze tailpiece that screws to the tail block, for example the StewMac trapeze (5.1 in long, 50 mm string spread) | The flat at the tail is 56 mm (2.2 in) wide, so a standard 51 mm (2 in) mounting plate sits on it with about 2.5 mm to spare each side. Its screws go into the tail block (68 mm wide). See the break angle note below. |
| Pickup | EMG 91 | Comes with its mounting bracket, screws, output cable and battery cable. Use EMG's bracket from the fingerboard extension; the strap bracket in the model is a stand-in. |
| Pot | EMG 25KX2 stacked volume/tone, 25k, solderless (part 5129) | Comes with a split stacked knob. Check the bushing against the 3/8 in hole. |
| Output jack | Stereo 1/4 in jack with 9 V battery connector, e.g. EMG's prewired jack, or a Switchcraft stereo jack | 3/8 in hole; the 2 mm plate suits a short-thread jack. |
| Battery | 9 V alkaline or lithium | End-loads into the sleeve behind the door. |

**Break angle and tailpiece length.** Because the plateau is high, a shorter trapeze gives more break angle, not less. With the string anchor 5.5 mm above the top: 50 mm from the tail edge gives about 13°, 100 mm about 16°, and a 5 in (125 mm) trapeze about 21.5°. Over about 18° the bridge sees noticeably more down-pressure (about a third more at 22°). The top has the margin for it, but the extra drive also changes the sound, so a longer trapeze, with its anchor 2–3 in from the tail, is the conservative choice.

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

Feeds and speeds depend on your machine; use the CAM program's defaults for hardwood (mahogany, maple, Richlite) and softwood (spruce) as a starting point. Run the grain lengthwise on the top and back, use compression or downcut end mills for profiles and a fresh ball-nose for the facets, and finish the spruce facets with a light 0.1–0.2 mm pass cutting with the grain where your CAM allows.

## Parts

| Part | Folder | Stock (L × W × T) | Setups | What gets cut |
|---|---|---|---|---|
| Top | `top/` | Bookmatched quartersawn Sitka spruce, 16 × 20 × 1-1/8 in | 2 | A: outside facets and flat skirt. B (flipped): inside of the plate, ring beam, two tone bars, upper-bout cross bar, and the 4 × 8 mm lip around each soundhole (it joins the ring beam, so there's no slot too narrow for the 1/4 in cutter). Soundholes, the 6 mm pickup wire hole and the tenon notch at the neck end cut through. |
| Back | `back/` | Bookmatched bigleaf maple, 16 × 20 × 3/4 in | 2 | Exported inside-up. A: inside (flat glue border, underside of the ridge). B: outside facets, ridge and flat border. All faces are planes, so a flat end mill with a fine stepover finishes them. |
| Rim, lower half | `rim_lower/` | Mahogany body blank, 21 × 14 × 1¾ in, planed to 1.12 in | 1 | From the split face: pocket the recess 20.55 mm deep, cut through inside the glue ledge, cut the control window notch 15 mm deep through the treble lower-bout wall, drill 4 dowel holes 10 mm deep, profile the outside. |
| Rim, upper half | `rim_upper/` | Same | 2 | Same as the lower half (including the control window notch), then flip and cut the open-top neck mortise 26 mm deep from the outer face. You can also cut the mortise after the two halves are glued. |
| Neck | `neck/` | Glued-up neck blank (see Stock to buy), at least 26 × 3.5 × 3 in overall | 2 | Neck frame: X along the neck, Z0 at the fingerboard gluing plane. A: gluing plane (including the tongue under the fingerboard extension), truss rod channel, headstock face. B (flipped): back carve, heel, tongue underside, tenon, volute, headstock back. |
| Fingerboard | `fingerboard/` | Richlite, 20½ × 2¾ × 3/8 in | 1 | Flat side down: 12 in radius top, 22 fret slots, taper profile. |

The `flat/` folder has 2D cut files:
- `faceplate_and_truss_cover.dxf`: 1.6 mm maple headstock overlay (finished to match the top) with tuner holes and truss rod access, plus the truss rod cover and its screw holes.
- `control_plate_2mm.dxf`: the rim access plate, 100 × 44 mm, with 9.5 mm (3/8 in) holes for the stereo jack and the stacked pot, the 32 × 22 mm battery door opening, four screw holes, and the outline of the window it covers.
- `goldfoil_bracket_flat.dxf`: flat pattern for the optional gold foil bracket (see below): outline, two bend lines, two holes for screws into the underside of the fingerboard end, and two for the pickup's mounting screws. Sizes are placeholders until the pickup is measured.

## Build order

1. **Rim:** glue the two halves on four 6 × 20 mm dowels. Check that the faces are flat and the ledges line up.
2. **Back:** glue the back onto the lower ledge.
3. **Truss rod:** fit a 15 in dual-action rod into the neck channel. It adjusts at the headstock, under the cover.
4. **Top:** glue the top onto the upper ledge and neck block. It has a notch at the neck end, open to the end face, directly over the mortise.
5. **Neck:** set the neck from above, after the body is closed. The tenon drops through the notch in the top into the mortise, and the tongue lands on the top over the neck block; glue all three faces (tenon, tongue on top, heel on the end face). The tenon floor sets the height, so you can still check the neck angle with the body closed.
6. **Chamfer:** cut the 45° edge chamfer around the top and back edges with a V-bit, 3.5 mm deep. It fades to nothing over the last 28 mm before the neck's end face; blend that taper by hand or with a ramped toolpath.
7. **Wiring:** feed the EMG 91 lead through the 6 mm hole in the top under the fingerboard extension. Fix the battery sleeve to the back of the control plate around the door opening, with the battery clip at its inner end. Mount the stacked pot and the stereo jack on the plate, wire everything through the window (EMG's solderless harness works), and screw the plate on. The EMG is active, so there is no string or tailpiece ground; don't add one. To change the battery, open the door and slide the 9 V out; the plate stays put.
8. **Finish the top:** mask a 1/4 in (6.35 mm) band around the top's edge, including the chamfer, and dye the rest of the top black. The mask line is the `TOP_DYE_MASK` layer in `model/juris_c75a_templates.dxf`. Pull the mask, then shellac over everything; the undyed band reads as a faux binding. Finish the maple headstock face the same way (dye, then shellac), so the top and headstock match; mask the truss rod cover recess and tuner holes.
9. **Fingerboard:** glue it on, fret it, and drill the side dots by hand.
10. **Hand-drilled holes:** drill the endpin and the control plate's four pilot holes square to the rim. Positions are in `positions.json`.
11. **Tuner holes:** drill them with the headstock face held flat, using a drill press with a 13° wedge or a separate CNC setup. They are 11/32 in (8.73 mm) for the Grover 133's press-in bushings; check the fit in scrap first. Each tuner also takes one small screw on the back of the headstock; mark it from the tuner itself.
12. **Bridge:** set the flat-bottomed base on the plateau (it needs no fitting) and the compensated top on its thumbwheels. Put the treble end about 1.5 mm and the bass end about 3.5 mm behind the scale line, then set intonation by sliding the bridge. The saddle sits about 31.6 mm (1.24 in) above the plateau, so check that your base and top reach that with adjustment to spare.

## Check against your hardware before cutting

These dimensions come from typical parts. Measure yours and change the parameter in `build.py`, then run `python build.py && python cnc.py` to regenerate everything.

| Part | Assumed | Parameter |
|---|---|---|
| Tuners | Grover 133 Vintage Deluxe (safety post), 11/32 in holes, 8.8 mm press-in bushings | `tuner_hole` |
| Truss rod | 15 in long, needs a 1/4 × 3/8 in channel | `truss_len`, `truss_w`, `truss_d` |
| Bridge base | Gretsch 008-0636-200, measured 135 × 14.95 × 11.75 mm, posts at 74 mm; plateau 141 mm wide | `BASE`, `plat_half_w` |
| Tailpiece | Standard trapeze, 51 mm (2 in) mounting plate | `TAIL_HALF` in `outline.py` (flat is 2 × 28 mm) |
| Pickup | EMG 91, 86.3 × 38.1 × 10.9 mm (EMG sheet, 0.43 in max height), on a fingerboard-extension bracket | `EMG91`, `build_pickup_bracket` |
| Battery | 9 V, 48.5 × 26.5 × 17.5 mm, end-loaded through a door in the plate | `BATTERY`, `CTRL_DOOR`, `SLEEVE_DEPTH` |
| Output jack | Stereo (EMG switches the battery through the jack), 3/8 in hole | `CTRL_HOLES` |
| Pots and jack | Stacked pot and stereo jack, 3/8 in bushings (9.5 mm holes) on a 2 mm plate | `CTRL_HOLES`, `CTRL_PLATE` |
| Strings | 13–56 set for the action and bridge height | `action12` |

## Optional pickup: Novak DeArmond-style gold foil

The model and viewer carry a second pickup option: a Curtis Novak DeArmond-style gold foil in the mustache casing, hung on a bent bracket in place of the EMG (choose it under Pickup in the viewer).

- **Bracket:** one 1.5 mm brass or aluminium strap. A tab screws to the underside of the fingerboard end, it drops about 5 mm just past the end of the board, and a shelf runs across under the casing's ears. The pickup screws down to the shelf through its ears over springs, which sets the height (about 3 mm of spring travel). The shelf clears the top by about 7 mm. Flat pattern: `flat/goldfoil_bracket_flat.dxf`.
- **Sizes are placeholders.** I couldn't find published dimensions for Novak's mustache casing, so the model uses typical DeArmond mustache proportions: 124 mm tip to tip, a 78 × 34 mm housing 9 mm tall, and ear screws at ±56 mm. Measure the pickup when it arrives and change `GOLDFOIL` in `build.py`, then rebuild.
- **Wiring:** it's passive. Use a 500k stacked volume/tone pot instead of the EMG's 25k, run a string ground to the tailpiece, and ground the casing. The battery door and sleeve aren't needed; the control plate can be cut without the door, or keep it to swap back to the EMG later.
- **Lead:** it runs through the same 6 mm hole under the fingerboard end.

## Voicing the top

The top is graduated light: the plateau is 3.5 mm, the falls are 2.8 mm at the plateau rising to 3.2 mm at the skirt, the ring beam is 5 × 6 mm, the tone bars 6 × 5 mm and the cross bar 6 × 5 mm. It weighs about 150 g. Thinner than this gains nothing: the arched shape sets most of the stiffness, so the strength check shows the same flex under the bridge for anything lighter, with less safety margin. The arch was lowered from 14 to 11 mm for the same reason (the neck sits 3 mm lower to keep the bridge height); that made the top about 56% more flexible under the bridge. Going to 10 mm adds almost nothing and costs strength.

- **Creases:** each fold between facets has a concave side, inside or out. Finish all the facets with the 1/4 in ball-nose and leave its 3 mm radius in those inside corners; don't sharpen them by hand. The convex ridges can stay crisp, with just a light sanding break so they don't chip. With that radius the creases carry no more stress than the open facets.
- **Tap tuning (optional):** after machining and before gluing, hold the top at a node and tap the plateau. The main tap tone should sit between two semitones and away from the air resonance (about 127 Hz). If it lands right on a note you'd play a lot, take a few tenths of a millimetre off the tone bars by hand, which lowers it.
- **Neck dead spots:** `neck_modes.py` checks the neck, fingerboard and headstock with the tuners on, held rigid at the body. The mode most likely to deaden a note is the second up-and-down bend, near 190 Hz: around open G, the D string 4th–6th frets and the A string 9th–11th. In the real guitar it will land somewhat lower, so it can't be pinned to one note. If a note there sounds short once the guitar is strung up, try a small clamp-on mass on the headstock: a few grams moves the mode off the note. Lighter tuners raise it, heavier ones lower it. Carbon rods beside the truss rod were modelled and moved this mode by under 1%, so they're left out.
- **Humidity:** like any solid spruce top, it cracks along the grain if it dries out. Keep it between 40 and 55% relative humidity. The faceted shape gives a little more room for seasonal movement than a flat top.

## Known limits

- **Outline accuracy:** the body outline was traced from the C-75 PDF at about 0.13 mm per pixel. If you have the original CAD outline, replace the traced edges in `outline.py`.
- **Curved surfaces:** the STEP models describe the neck carve, volute and fingerboard radius as many small flat faces, about 0.005 mm from the true shape. That's fine for 3D finishing, but they won't behave as single smooth surfaces in CAD.
- **Undercut check:** casts vertical rays on a 0.5 mm grid (0.25 mm for the fingerboard and bridge). Four cells on the neck read 1 mm² in total, all where rays graze the curved sides; recasting those spots gives clean results, so it's numerical noise.
- **Hardware not modelled for machining:** the tailpiece, bridge, pickup, knobs, pickup bracket and control plate hardware are purchased or fabricated. The bridge in the 3D model is a stand-in at the right height and position.
