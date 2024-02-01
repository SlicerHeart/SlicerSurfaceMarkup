#ifndef GENERIC_SURFACE_SOURCE_H
#define GENERIC_SURFACE_SOURCE_H

#include "vtkBezierSurfaceSource.h"
#include "vtkNURBSSurfaceSource.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkSmartPointer.h"
#include "vtkPoints.h"

class GenericSurfaceSource {
public:
    // Enum to represent the type of surface source
    enum SurfaceType { Bezier, Nurbs, Custom };

    // Constructor and Destructor
    GenericSurfaceSource();
    ~GenericSurfaceSource();

    // Methods to set the surface source
    void SetBezier(vtkBezierSurfaceSource* source);
    void SetNurbs(vtkNURBSSurfaceSource* source);
    void SetCustom(vtkPolyDataAlgorithm* source);

    // Method to perform action based on the current type
    void PerformActionBasedOnType(vtkPoints* points);

private:
    // Current type of the surface source
    SurfaceType currentType;

    // Smart pointer to the current surface source
    vtkSmartPointer<vtkPolyDataAlgorithm> currentSource;

    // Preventing copying and assignment
    GenericSurfaceSource(const GenericSurfaceSource&) = delete;
    GenericSurfaceSource& operator=(const GenericSurfaceSource&) = delete;
};

#endif // GENERIC_SURFACE_SOURCE_H



