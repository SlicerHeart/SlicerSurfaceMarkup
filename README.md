# SurfaceMarkup extension for 3D Slicer

Adds support for grid surfaces such as NURBS or Bézier surfaces.
- New markup type called GridSurface that shows up in the Markups module
- Placement by the three initial control points defining three of the four corners thus the surface plane. The fourth corner is added symmetrically.
- The grid control points are added automatically to uniformly fill the defined surface area
- The surface can be edited by moving the control points
- Output model node can be associated for the surface patch the display options of which can be edited in the Models module

Example NURBS surface:
![NURBS](https://raw.githubusercontent.com/SlicerHeart/SlicerSurfaceMarkup/master/Screenshots/NURBS.png?raw=true)

NURBS surface wrapped around (i.e. cylinder-like) shaped like a heart valve:
![NURBS Valve](https://github.com/SlicerHeart/SlicerSurfaceMarkup/raw/master/Screenshots/NURBS_WrapAround_Valve.png)

Example Bézier surface:
![Bezier](https://raw.githubusercontent.com/SlicerHeart/SlicerSurfaceMarkup/master/Screenshots/Bezier.png?raw=true)

## Funding sources:

- The Cora Topolewski Fund at the Children's Hospital of Philadelphia (CHOP)

- CHOP Frontier Grant (Pediatric Valve Center)

- National Heart, Blood, and Lung Institute (NHLBI) (R01 HL153166)
