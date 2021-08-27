cmake_minimum_required(VERSION 3.13.4)

project(GridSurface)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://www.slicer.org/wiki/Documentation/Nightly/Extensions/GridSurface")
set(EXTENSION_CATEGORY "Liver")
set(EXTENSION_CONTRIBUTORS "Csaba Pinter (Pixel Medical / Ebatinca), Andras Lasso (PerkLab, Queen's), Rafael Palomar (Oslo University Hospital / NTNU)")
set(EXTENSION_DESCRIPTION "3D Slicer extension adding support for grid surfaces")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/GridSurface.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/GridSurface/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(GridSurfaceMarkups)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
