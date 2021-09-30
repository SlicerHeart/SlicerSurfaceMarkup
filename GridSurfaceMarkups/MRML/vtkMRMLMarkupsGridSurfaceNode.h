/*==============================================================================

 Distributed under the OSI-approved BSD 3-Clause License.

  Copyright (c) Oslo University Hospital. All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.

  * Neither the name of Oslo University Hospital nor the names
    of Contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  This file was originally developed by Rafael Palomar (The Intervention Centre,
  Oslo University Hospital) and was supported by The Research Council of Norway
  through the ALive project (grant nr. 311393).

==============================================================================*/

#ifndef __vtkMRMLMarkupsGridSurfaceNode_h_
#define __vtkMRMLMarkupsGridSurfaceNode_h_

#include "vtkSlicerGridSurfaceMarkupsModuleMRMLExport.h"

// MRML includes
#include <vtkMRMLMarkupsNode.h>
#include <vtkMRMLModelNode.h>

//VTK includes
#include <vtkWeakPointer.h>

//-----------------------------------------------------------------------------
class VTK_SLICER_GRIDSURFACEMARKUPS_MODULE_MRML_EXPORT vtkMRMLMarkupsGridSurfaceNode
: public vtkMRMLMarkupsNode
{
public:
  static vtkMRMLMarkupsGridSurfaceNode* New();
  vtkTypeMacro(vtkMRMLMarkupsGridSurfaceNode, vtkMRMLMarkupsNode);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  //--------------------------------------------------------------------------------
  // MRMLNode methods
  //--------------------------------------------------------------------------------
  const char* GetIcon() override {return ":/Icons/MarkupsGeneric.png";}
  const char* GetAddIcon() override {return ":/Icons/MarkupsGenericMouseModePlace.png";}
  const char* GetPlaceAddIcon() override {return ":/Icons/MarkupsGenericMouseModePlaceAdd.png";}

  vtkMRMLNode* CreateNodeInstance() override;

  /// Get node XML tag name (like Volume, Model)
  ///
  const char* GetNodeTagName() override {return "MarkupsGridSurface";}

  /// Get markup name
  const char* GetMarkupType() override {return "GridSurface";}

  /// Get markup short name
  const char* GetDefaultNodeNamePrefix() override {return "GS";}

  /// \sa vtkMRMLNode::CopyContent
  vtkMRMLCopyContentDefaultMacro(vtkMRMLMarkupsGridSurfaceNode);

  void ProcessMRMLEvents(vtkObject* caller, unsigned long event, void* callData) override;

public:
  // GridSurface type enum defines the calculation to create the surface from the control points
  enum
    {
    GridSurfaceTypeBezier,
    //GridSurfaceTypeThinPlate,
    GridSurfaceType_Last
    };

  /// GridSurfaceType represents the method that is used to calculate the size of the ROI.
  /// BOX ROI does not require control points to define a region, while the size of a BOUNDING_BOX ROI will be defined by the control points.
  vtkGetMacro(GridSurfaceType, int);
  void SetGridSurfaceType(int gridSurfaceType);
  static const char* GetGridSurfaceTypeAsString(int gridSurfaceType);
  static int GetGridSurfaceTypeFromString(const char* gridSurfaceType);

  ///@{
  /// Number of control points on each side of the grid
  vtkGetVector2Macro(GridResolution, int);
  void SetGridResolution(const int gridResolution[2]);
  void SetGridResolution(int x, int y);
  ///@}

  //TODO:
  void UpdateGridSurfaceFromControlPoints();
  //TOOD:
  void UpdateControlPointsFromGridSurface();

protected:
  int GridSurfaceType{vtkMRMLMarkupsGridSurfaceNode::GridSurfaceTypeBezier};

  bool IsUpdatingControlPointsFromGridSurface{false};
  bool IsUpdatingGridSurfaceFromControlPoints{false};

  int GridResolution[2] { 4, 4 };

  /// In order to be able to resample to a new grid resolution after a change.
  int PreviousGridResolution[2] { 0, 0 };

protected:
  void FillControlPointGridFromCorners(double position_0[3], double position_1[3], double position_2[3], vtkPoints* controlPoints);
  void ResampleToNewGridResolution();

protected:
  vtkMRMLMarkupsGridSurfaceNode();
  ~vtkMRMLMarkupsGridSurfaceNode() override = default;

private:
  vtkWeakPointer<vtkMRMLModelNode> Target;

private:
  vtkMRMLMarkupsGridSurfaceNode(const vtkMRMLMarkupsGridSurfaceNode&);
  void operator=(const vtkMRMLMarkupsGridSurfaceNode&);

  friend class vtkSlicerGridSurfaceRepresentation2D;
  friend class vtkSlicerGridSurfaceRepresentation3D;
};

#endif //__vtkMRMLMarkupsGridSurfaceNode_h_
