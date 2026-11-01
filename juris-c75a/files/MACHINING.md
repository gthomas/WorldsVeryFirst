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
| 1 | Mahogany body blank, 2-piece, standard thickness | 21 × 14 × 1¾ in (milled from 8/4) | Rim | One ring, cut whole from a standard blank: no stacking, glue-up or dowels. Finished 18.05 × 13.6 × 1.75 in, so a 14 in blank leaves about 0.2 in per side; hold it by the ends (the pin holes in the setup DXF are at the ends), or buy an oversized 15 in blank. The centre seam is fine: it runs through the neck and tail blocks. |
| 1 | Mahogany neck blank | 30 × 3 × 1 in | Neck | Needs three glue-ups before machining: a 13° scarf-jointed headstock, a stacked heel (two 1 in layers glued under the heel end, about 3 in of total depth there), and headstock ears to reach 3.5 in width. Or buy a one-piece blank of at least 26 × 3.5 × 3 in and skip the glue-ups. |
| 1 | Quartersawn (rift and quartered) Port Orford cedar, flat stock | 2 pieces, each 22 × 7½ in, at least 7/8 in thick after flattening | Top | Joined it is about 22 × 15 in. The finished top needs 0.76 in of thickness across most of its width (the plateau and the soundhole lips sit at different heights), so use flat boards, not tapered archtop wedges. A slip-match from one board is fine under the black dye. Sitka spruce works with no design changes (about 10 g lighter, a little weaker in bending). Wear a respirator and use extraction: POC dust and oil are respiratory irritants. |
| 1 | Quartersawn bigleaf maple back set, bookmatched | 2 halves, each about 22 × 8 in, about 1-1/8 in thick (5/4) | Back | Finished 1.0 in at most (the 22 mm arch). Plain or flamed makes no acoustic difference. A stiff maple back reflects energy into the top and machines crisply. |
| 1 | Richlite guitar fingerboard blank (Black Diamond), unslotted | 20 × 2¾ × 5/16 in (Richlite's standard blank) | Fingerboard | Finished 18 × 2.22 × 0.29 in. Barely moves with humidity, so fret ends stay put in dry weather. Machines like a hard hardwood; cut fret slots with a sharp saw or slotting cutter and press the frets in (a drop of thin glue in each slot is a good idea). |
| 1 | Maple headstock veneer | 8 × 4 × 1/16 in | Faceplate | Resaw it from the back stock's offcut. Finished head is 7.5 × 3.6 in (the built C-75's headstock shape), a hair bigger than the common 190 × 90 mm veneers. Finished to match the top (dye and shellac along with it). |
| 1 | Black fibre or aluminium sheet | 1½ × 1½ in, 2 mm | Rim control plate | Holds the stacked pot only. |
| 1 | ¼ in MDF or ply | 6 × 3 in | Router template for the control face | Cut from `rim_control_face_template.dxf`; rout with a pattern bit. |
| 1 | Brass or steel strap | 4 × ½ in, 1.5 mm | Pickup bracket | Bent to hang the EMG 91 from the underside of the fingerboard extension, unless you use a fingerboard-mount bracket for it. |
| 1 | Stacked (concentric) pot, 25k, with a concentric knob | 3/8 in bushing | Volume and tone | EMG actives use 25k pots. |
| 1 | Dual-action truss rod | 15 in | Neck | Check the channel size against your rod. |

## Lumber: sources and cost

Prices are as listed on 1 November 2026, before tax and shipping. Most are one-off pieces or depend on stock, so check before ordering.

| Part | Buy | Source | Price |
|---|---|---|---|
| Top | Port Orford cedar, 4/4 Select & Better, rift and quartered, one board at least 7½ in wide × 48 in (about 2.7 bd ft). Cut it in two for a slip-matched pair. Make sure it flattens to at least 7/8 in. If not, buy 8/4 and resaw it for a true bookmatch (about 5.3 bd ft, about $213). | [Hearne Hardwoods](https://www.hearnehardwoods.com/cedar-port-orford-lumber/), $40/bd ft (call for stock: 610-932-7400) | about $107 |
| Back | Bigleaf maple, 5/4 Select & Better, 8 in × 48 in (about 3.3 bd ft). Ask for quartered stock. The listed grain is mixed, and the cheapest figure grade is light quilt. | [Hearne Hardwoods](https://www.hearnehardwoods.com/maple-big-leaf-lumber/), from $15/bd ft | about $50 |
| Rim | Honduran mahogany 2-piece body blank, 20½ × 14 × 1⅞ in. Face it to 1¾ in. | [StewMac](https://www.stewmac.com/tonewoods/wood-parts/honduran-mahogany-body-blank-2-piece---505) (WoodStax piece 505; similar blanks run about $105–130) | $114.73 |
| Neck | Honduran mahogany neck blank, 1 × 3 × 30 in, plus a 1 × 4 × 36 in blank for the heel stack and headstock ears | [StewMac](https://www.stewmac.com/tonewoods/wood-parts/honduran-mahogany-neck-blanks), $35.49 + $35.49 | $70.98 |
| Fingerboard | Richlite Black Diamond guitar fretboard blank, 20 × 2.75 × 0.3125 in, unslotted | [Richlite](https://www.richlite.com/products/guitar-fretboard-blank-black-diamond) $20 (sold out); [Solo Music Gear](https://www.solomusicgear.com/product/richlite-guitar-fretboard-blank-20-unslotted-black-diamond/) CA$23.99 (in stock) | about $20 |
| Headstock veneer | From the back's offcut | — | $0 |
| **Total** | | | **about $363** (about $469 with 8/4 cedar resawn for a true bookmatch) |

Other options looked at:
- [Allied Lutherie's Port Orford cedar archtop set](https://alliedlutherie.com/products/port-orford-cedar-archtop), $300, sold out. Tapered archtop wedges don't suit this top.
- [NW Timber's musical-grade POC](https://nwtimber.com/rare-wide-boards/port-orford-white-cedar): 10/4 × 9⅝ × 40½ in at $697, or 1/4 in bookmatches (too thin).
- [Maderas Barber's European curly maple archtop backs](https://maderasbarber.com/en/archtop-guitar-backs/898-a-curly-maple-backs-guitar-archtop.html), $130, tapered 25 → 8 mm. These suit a carved back better than this machined one.

## Hardware to buy

| Item | Choice | Fit notes |
|---|---|---|
| Tuners | Genuine Grover Deluxe SD90, 3+3, nickel, safety post (133N), as on the C-75 | 15:1, 6 mm posts. 11/32 in post holes with 8.8 mm press-in bushings; one back screw per tuner. |
| Bridge base | Gretsch 008-0636-200 rosewood bridge base, flat bottom (for flat-top guitars), as on the C-75 | Sits flat on the bridge plateau with no fitting; thumbwheel posts at 74 mm. Measured base: 135 × 14.95 mm, 11.5–12 mm tall. The plateau is 141 mm (5.55 in) wide across the strings, so the base sits fully on it with 3 mm to spare each side, and its feet land over the 5 × 6 mm ring beam. Strength check with this base (Port Orford cedar top): about 0.51 mm sink under the feet, 0.72 mm at most, stresses at least 4× under clear-wood strength. Expect it to settle a little more over the first months; take that up with the thumbwheels (checked without the soundhole lips, which only add stiffness). |
| Bridge top | A standard wooden compensated archtop bridge top (ebony or rosewood) with post holes at 74 mm | Must match the base's 74 mm post spacing. |
| Tailpiece | A trapeze tailpiece that screws to the tail block, for example the StewMac trapeze (5.1 in long, 50 mm string spread) | The flat at the tail is 56 mm (2.2 in) wide, so a standard 51 mm (2 in) mounting plate sits on it with about 2.5 mm to spare each side. Its screws go into the tail block (68 mm wide). See the break angle note below. |
| Pickup | EMG 91 | Comes with its mounting bracket, screws, output cable and battery cable. Use EMG's bracket from the fingerboard extension; the strap bracket in the model is a stand-in. |
| Pot | EMG 25KX2 stacked volume/tone, 25k, solderless (part 5129) | Comes with a split stacked knob. Check the bushing against the 3/8 in hole. |
| Battery box | Gotoh BB-04 single 9 V (or StewMac #3577, same format) | 70 × 29 mm flange, four screws on 61.5 × 20.7 mm centres. The body drops through a 57 × 25 mm opening (1/4 in bit corners) and hangs about 13 mm inside the wall. Check the part before routing. |
| Output jack | Switchcraft 152B stereo long-thread jack (EMG switches the battery through the stereo jack) | 12 mm hole through the 20 mm wall in the treble tail facet; the 152B takes panels up to 1¼ in. |
| Battery | 9 V alkaline or lithium | Drops into the BB-04 under its hinged lid. |

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

Feeds and speeds depend on your machine; use the CAM program's defaults for hardwood (mahogany, maple, Richlite) and softwood (Port Orford cedar) as a starting point. Run the grain lengthwise on the top and back, use compression or downcut end mills for profiles and a fresh ball-nose for the facets, and finish the cedar facets with a light 0.1–0.2 mm pass cutting with the grain where your CAM allows.

## Parts

| Part | Folder | Stock (L × W × T) | Setups | What gets cut |
|---|---|---|---|---|
| Top | `top/` | Bookmatched quartersawn Port Orford cedar, 16 × 20 × 1-1/8 in | 2 | A: outside facets and flat skirt; the 25 mm crease blends in the valleys come from the ball-nose finishing pass (keep the stepover fine along them). B (flipped): inside of the plate, with the same blends under each ridge, ring beam, two tone bars, upper-bout cross bar, and the 4 × 8 mm lip around each soundhole (it joins the ring beam, so there's no slot too narrow for the 1/4 in cutter). Soundholes, the 6 mm pickup wire hole and the tenon notch at the neck end cut through. |
| Back | `back/` | Bookmatched bigleaf maple, 16 × 20 × 1-1/8 in | 2 | Exported inside-up. A: inside (flat glue border, underside of the ridge). B: outside facets, ridge and flat border. All faces are planes, so a flat end mill with a fine stepover finishes them. |
| Rim | `rim/` | Mahogany body blank, 21 × 14 × 1¾ in (standard) | 1 | Everything from the top face: pocket the recess inside the 8 mm wall down to the back's glue ledge (12 × 8 mm), cut through inside the ledge, cut the open-top neck mortise 26 mm deep, profile the outside. The top glues to the full 8 mm wall. Afterwards, rout the pot window and the BB-04 opening through the treble lower-bout wall from outside with `rim_control_face_template.dxf`, and drill the 12 mm jack hole in the treble tail facet (see `positions.json`). |
| Neck | `neck/` | Glued-up neck blank (see Stock to buy), at least 26 × 3.5 × 3 in overall | 2 | Neck frame: X along the neck, Z0 at the fingerboard gluing plane. A: gluing plane (including the tongue under the fingerboard extension), truss rod channel, headstock face. B (flipped): back carve, heel, tongue underside, tenon, volute, headstock back. |
| Fingerboard | `fingerboard/` | Richlite, 20½ × 2¾ × 3/8 in | 1 | Flat side down: 12 in radius top, 22 fret slots, taper profile. |

The `flat/` folder has 2D cut files:
- `faceplate_and_truss_cover.dxf`: 1.6 mm maple headstock overlay (finished to match the top) with tuner holes and truss rod access, plus the truss rod cover and its screw holes.
- `control_plate_2mm.dxf`: the pot plate, 36 × 36 mm with 4 mm corners, a 9.5 mm (3/8 in) hole for the stacked pot, four screw holes, and the 26 mm window it covers.
- `rim_control_face_template.dxf`: full-size template for the flat treble lower-bout face. It has the BB-04 opening, the pot window, the plate and flange outlines, the pilot holes, and the ends of the flat face for registration.
- `goldfoil_bracket_flat.dxf`: flat pattern for the optional gold foil bracket (see below): outline, two bend lines, two holes for screws into the underside of the fingerboard end, and two for the pickup's mounting screws. Sizes are placeholders until the pickup is measured.

## Build order

1. **Rim:** one piece, so there is nothing to glue. Rout the control-plate window through the treble lower-bout wall from outside with a template, centred at the height in `positions.json`.
2. **Back:** glue the back onto the ledge at the bottom of the rim.
3. **Truss rod:** fit a 15 in dual-action rod into the neck channel. It adjusts at the headstock, under the cover.
4. **Top:** glue the top onto the top of the rim wall (an 8 mm glue land) and the neck block. It has a notch at the neck end, open to the end face, directly over the mortise.
5. **Neck:** set the neck from above, after the body is closed. The tenon drops through the notch in the top into the mortise, and the tongue lands on the top over the neck block; glue all three faces (tenon, tongue on top, heel on the end face). The tenon floor sets the height, so you can still check the neck angle with the body closed.
6. **Chamfer:** cut the 45° edge chamfer around the top and back edges with a V-bit, 3.5 mm deep. It fades to nothing over the last 28 mm before the neck's end face; blend that taper by hand or with a ramped toolpath.
7. **Wiring:** feed the EMG 91 lead through the 6 mm hole in the top under the fingerboard extension. Mount the stacked pot on the plate and the 152B in the tail facet. Run the BB-04's clip lead in through its opening and wire everything through the pot window (EMG's solderless harness works), then screw the plate and the battery box on. The EMG is active, so there is no string or tailpiece ground; don't add one. To change the battery, flip open the BB-04 lid; no screws.
8. **Finish the top:** mask a 1/4 in (6.35 mm) band around the top's edge, including the chamfer, and dye the rest of the top black. The mask line is the `TOP_DYE_MASK` layer in `model/juris_c75a_templates.dxf`. Pull the mask, then shellac over everything; the undyed band reads as a faux binding. Finish the maple headstock face the same way (dye, then shellac), so the top and headstock match; mask the truss rod cover recess and tuner holes.
9. **Fingerboard:** glue it on, fret it, and drill the side dots by hand.
10. **Hand-drilled holes:** rout the control face's two openings with the template, drill the 12 mm jack hole in the treble tail facet, and drill the plate's and the battery box's pilot holes square to the rim. Positions are in `positions.json`.
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
| Battery box | Gotoh BB-04, 70 × 29 mm flange, 57 × 25 mm opening | `BB`, `BATTERY` |
| Output jack | Switchcraft 152B stereo, 12 mm hole in the treble tail facet | `JACK_HOLE`, `jack_frame` |
| Pot | Stacked pot, 3/8 in bushing (9.5 mm hole) on a 36 mm 2 mm plate over a 26 mm window | `CTRL_HOLES`, `CTRL_PLATE`, `CTRL_WIN_D` |
| Strings | 13–56 set for the action and bridge height | `action12` |

## Optional pickup: Novak DeArmond-style gold foil

The model and viewer carry a second pickup option: a Curtis Novak DeArmond-style gold foil in the mustache casing, hung on a bent bracket in place of the EMG (choose it under Pickup in the viewer).

- **Bracket:** one 1.5 mm brass or aluminium strap. A tab screws to the underside of the fingerboard end, it drops about 5 mm just past the end of the board, and a shelf runs across under the casing's ears. The pickup screws down to the shelf through its ears over springs, which sets the height (about 3 mm of spring travel). The shelf clears the top by about 7 mm. Flat pattern: `flat/goldfoil_bracket_flat.dxf`.
- **Sizes are placeholders.** I couldn't find published dimensions for Novak's mustache casing, so the model uses typical DeArmond mustache proportions: 124 mm tip to tip, a 78 × 34 mm housing 9 mm tall, and ear screws at ±56 mm. Measure the pickup when it arrives and change `GOLDFOIL` in `build.py`, then rebuild.
- **Wiring:** it's passive. Use a 500k stacked volume/tone pot instead of the EMG's 25k, run a string ground to the tailpiece, and ground the casing. The battery box isn't needed. Leave its opening uncut, or fit the box anyway so you can swap back to the EMG later.
- **Lead:** it runs through the same 6 mm hole under the fingerboard end.

## Voicing the top

The top is graduated light: the plateau is 3.5 mm, the falls are 2.8 mm at the plateau rising to 3.2 mm at the skirt, the ring beam is 5 × 6 mm, the tone bars 6 × 5 mm and the cross bar 6 × 5 mm. It weighs about 150 g. Thinner than this gains nothing: the arched shape sets most of the stiffness, so the strength check shows the same flex under the bridge for anything lighter, with less safety margin. The arch was lowered from 14 to 11 mm for the same reason (the neck sits 3 mm lower to keep the bridge height); that made the top about 56% more flexible under the bridge. Going to 10 mm adds almost nothing and costs strength.

- **Creases:** each fold between facets has a concave side, inside or out. Finish all the facets with the 1/4 in ball-nose and leave its 3 mm radius in those inside corners; don't sharpen them by hand. The convex ridges can stay crisp, with just a light sanding break so they don't chip. With that radius the creases carry no more stress than the open facets.
- **Tap tuning (optional):** after machining and before gluing, hold the top at a node and tap the plateau. The main tap tone should sit between two semitones and away from the air resonance (about 127 Hz). If it lands right on a note you'd play a lot, take a few tenths of a millimetre off the tone bars by hand, which lowers it.
- **Neck dead spots:** `neck_modes.py` checks the neck, fingerboard and headstock with the tuners on, held rigid at the body. The mode most likely to deaden a note is the second up-and-down bend, near 200 Hz with the current headstock: around G3–G#3, so open G or its 1st fret, the D string 5th–6th frets and the A string 10th–11th. In the real guitar it will land somewhat lower, so it can't be pinned to one note. If a note there sounds short once the guitar is strung up, try a small clamp-on mass on the headstock: a few grams moves the mode off the note. Lighter tuners raise it, heavier ones lower it. Carbon rods beside the truss rod were modelled and moved this mode by under 1%, so they're left out.
- **Humidity:** like any solid softwood top, it cracks along the grain if it dries out. Keep it between 40 and 55% relative humidity. The faceted shape gives a little more room for seasonal movement than a flat top.

## Known limits

- **Outline accuracy:** the body outline was traced from the C-75 PDF at about 0.13 mm per pixel. If you have the original CAD outline, replace the traced edges in `outline.py`.
- **Curved surfaces:** the STEP models describe the neck carve, volute and fingerboard radius as many small flat faces, about 0.005 mm from the true shape. That's fine for 3D finishing, but they won't behave as single smooth surfaces in CAD.
- **Undercut check:** casts vertical rays on a 0.5 mm grid (0.25 mm for the fingerboard and bridge). Four cells on the neck read 1 mm² in total, all where rays graze the curved sides; recasting those spots gives clean results, so it's numerical noise.
- **Hardware not modelled for machining:** the tailpiece, bridge, pickup, knobs, pickup bracket and control plate hardware are purchased or fabricated. The bridge in the 3D model is a stand-in at the right height and position.
