# GridSurface extension for 3D Slicer

Adds support for grid surfaces such as NURBS or Bézier surfaces.
- New markup type called GridSurface that shows up in the Markups module
- Placement by the three initial control points defining three of the four corners thus the surface plane. The fourth corner is added symmetrically.
- The grid control points are added automatically to uniformly fill the defined surface area
- The surface can be edited by moving the control points

![NURBS](https://raw.githubusercontent.com/cpinter/SlicerGridSurface/master/Screenshots/NURBS.png?raw=true)

![Bezier](https://raw.githubusercontent.com/cpinter/SlicerGridSurface/master/Screenshots/Bezier.png?raw=true)
