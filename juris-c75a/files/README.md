# Juris C-75A: faceted acoustic archtop on the C-75 outline

Parametric model. Units are mm. Axes: x+ toward the treble side, y+ toward the headstock (y = 0 at the tail), z+ up (z = 0 at the back edge).

## Files
- `juris_c75a.glb`: full assembly (open it in Blender, Fusion or any glTF viewer)
- `top.stl`, `back.stl`, `rim.stl`, `neck.stl`, `fingerboard.stl`, `bridge.stl`: watertight solids for CAM
- `juris_c75a_templates.dxf`: 2D templates on separate layers. OUTLINE, RIM_INNER (7 mm wall), LINING_INNER (12 mm ledge), TOP_PLATEAU, TOP_BEAM, TOP_RIB, TOP_FOLDS, SOUNDHOLES, BACK_PLATEAU, BACK_FOLDS, CENTERLINE, BRIDGE_LINE
- `report.json`: derived numbers (bridge height, clearances, weight estimate, air volume)
- `fea.py` and `fea.json`: static finite-element check of the top under string load (quadratic tets, orthotropic quartersawn spruce, top fixed where it is glued to the rim and neck block)
- `build.py` and `outline.py`: the generator. Edit `P` in build.py and run `python build.py` (needs manifold3d, trimesh, mapbox_earcut, ezdxf, numpy)

## How the geometry is built
- **Outline:** traced from Juris_C_75_Body_2D.pdf and calibrated on its 18.000 and 13.500 dimensions. The lower bout is symmetrised, and the corners are filleted (12 to 19 mm, sized for a 12 mm end mill).
- **Top:** made of planes only, symmetric about the centreline. A flat hexagonal bridge plateau (100 x 155 mm) sits 14 mm above the edge inside an eight-sided arch zone inset 22 mm from the outline; 14 triangular facets run from the plateau to that octagon, and outside it the top is a flat skirt at edge level that sits on the rim and glue ledges. Falls taper 3.8 to 3.2 mm, the plateau is 5 mm thick, a 5 x 6 mm ring beam runs under the plateau edge, and two 6 x 9 mm tone bars with 35 mm scalloped ends run lengthwise under the plateau's long sides, and a 6 x 7 mm cross bar runs from the neck block's lower face out to the bass upper bout, tapering over its last 35 mm and stopping 2 mm short of the glue ledge (all on DXF layer TOP_BARS). The two soundholes are mirror images in the upper side facets, outboard of the bars; the treble one ends under the pickguard so the pickup lead drops into the body.
- **Back:** quartersawn bigleaf maple (top is quartersawn Sitka spruce). The same scheme with a larger plateau and a 6 mm arch, 3.2 mm thick throughout, unbraced.
- **Rim:** routed from one 12/4 mahogany blank, finished 57.1 mm (2.25") tall for a 2.5" body at the edge. It has 7 mm walls, 12 x 8 mm glue ledges routed top and bottom, solid neck and tail blocks, a neck mortise, a 1/2" side jack bore and an endpin bore.
- **Edges:** 3.5 mm 45-degree chamfer around the top and back edges, cut after glue-up through the plate and into the rim.
- **Heel and volute:** the '59 round carve runs all the way into the heel, which sweeps down to finish flush with the back edge on the body face. The heel is 53 mm wide, 1.5 mm inside the body's sides, and lands on a square, unchamfered end face. Behind the nut the round neck section morphs into the square headstock section with no step; the volute loft forms the headstock itself up to the first tuners, so there is no seam at the base. The edge chamfer fades out over the last 28 mm before the neck end face.
- **Neck:** set neck at 5 degrees with a mortise and tenon. The tenon sits below the top, and the neck clears the body at about fret 15. '59 roundback profile, 12" radius, 22 frets. 13-degree 3x3 Juris headstock: sides flaring straight from the nut width (no wings), a pointed gable top, flat face with square edges for standard tuners, ebony faceplate, pointed nickel truss rod cover, posts following the flare. A lofted volute (about 4.5 mm proud) blends the neck into the headstock behind the nut.
- **Pickguard:** floats parallel to the strings and never touches the top. The pickup hangs from it by its mounting tab. It is held at two points: a strap screwed into the treble side of the rim, and a 1.5 mm plate under the end of the fingerboard extension, screwed up into the ebony and through the guard (the guard and the board's underside are in the same plane there).
