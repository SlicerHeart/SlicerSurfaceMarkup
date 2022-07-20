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

#include "vtkSlicerGridSurfaceRepresentation3D.h"

#include "vtkMRMLMarkupsGridSurfaceNode.h"

// MRML includes
#include <qMRMLThreeDWidget.h>
#include <vtkMRMLDisplayableManagerGroup.h>
#include <vtkMRMLModelDisplayableManager.h>

// Slicer includes
#include <qSlicerApplication.h>
#include <qSlicerLayoutManager.h>

// VTK includes
#include <vtkActor.h>
#include <vtkCollection.h>
#include <vtkNew.h>
#include <vtkPlaneSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyLine.h>
#include <vtkProperty.h>

//------------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerGridSurfaceRepresentation3D);

//------------------------------------------------------------------------------
vtkSlicerGridSurfaceRepresentation3D::vtkSlicerGridSurfaceRepresentation3D()
  : Superclass()
{
  this->ControlPolygonTubeFilter->SetInputData(this->ControlPolygonPolyData.GetPointer());
  this->ControlPolygonTubeFilter->SetRadius(0.25);
  this->ControlPolygonTubeFilter->SetNumberOfSides(20);

  this->ControlPolygonMapper->SetInputConnection(this->ControlPolygonTubeFilter->GetOutputPort());

  this->ControlPolygonActor->SetMapper(this->ControlPolygonMapper);

  this->GridSurfaceActor->SetProperty(this->GridSurfaceProperty);
  this->ControlPolygonActor->SetProperty(this->ControlPolygonProperty);

  this->UpdateInterpolatorConnection();
}

//------------------------------------------------------------------------------
vtkSlicerGridSurfaceRepresentation3D::~vtkSlicerGridSurfaceRepresentation3D() = default;

//----------------------------------------------------------------------
void vtkSlicerGridSurfaceRepresentation3D::UpdateFromMRML(vtkMRMLNode* caller, unsigned long event, void* callData /*=nullptr*/)
{
  this->Superclass::UpdateFromMRML(caller, event, callData);

  auto gridSurfaceMarkupsNode = vtkMRMLMarkupsGridSurfaceNode::SafeDownCast(this->GetMarkupsNode());
  if (!gridSurfaceMarkupsNode || !this->IsDisplayable())
  {
    this->VisibilityOff();
    return;
  }

  this->UpdateGridSurface(gridSurfaceMarkupsNode);
  this->UpdateControlPolygon(gridSurfaceMarkupsNode);

  double diameter = (this->MarkupsDisplayNode->GetCurveLineSizeMode() == vtkMRMLMarkupsDisplayNode::UseLineDiameter ?
    this->MarkupsDisplayNode->GetLineDiameter() : this->ControlPointSize * this->MarkupsDisplayNode->GetLineThickness());
  this->ControlPolygonTubeFilter->SetRadius(diameter * 0.5);

  // Set same visibility options for the surface patch and the control polygon as the control points
  int controlPointType = Active;
  if (this->MarkupsDisplayNode->GetActiveComponentType() != vtkMRMLMarkupsDisplayNode::ComponentLine)
  {
    controlPointType = this->GetAllControlPointsSelected() ? Selected : Unselected;
  }
  this->GridSurfaceProperty->DeepCopy(this->GetControlPointsPipeline(controlPointType)->Property);
  this->ControlPolygonProperty->DeepCopy(this->GetControlPointsPipeline(controlPointType)->Property);

  // Use fill and outline visibility/opacity for the surface patch and the control polygon, respectively
  this->GridSurfaceActor->SetVisibility(this->MarkupsDisplayNode->GetFillVisibility());
  this->GridSurfaceProperty->SetOpacity(this->MarkupsDisplayNode->GetFillOpacity());
  this->ControlPolygonActor->SetVisibility(this->MarkupsDisplayNode->GetOutlineVisibility());
  this->ControlPolygonProperty->SetOpacity(this->MarkupsDisplayNode->GetOutlineOpacity());

  this->NeedToRenderOn();
}

//----------------------------------------------------------------------
void vtkSlicerGridSurfaceRepresentation3D::GetActors(vtkPropCollection *pc)
{
  this->Superclass::GetActors(pc);
  this->GridSurfaceActor->GetActors(pc);
  this->ControlPolygonActor->GetActors(pc);
}

//----------------------------------------------------------------------
void vtkSlicerGridSurfaceRepresentation3D::ReleaseGraphicsResources(vtkWindow *win)
{
  this->Superclass::ReleaseGraphicsResources(win);
  this->GridSurfaceActor->ReleaseGraphicsResources(win);
  this->ControlPolygonActor->ReleaseGraphicsResources(win);
}

//----------------------------------------------------------------------
int vtkSlicerGridSurfaceRepresentation3D::RenderOverlay(vtkViewport *viewport)
{
  int count = this->Superclass::RenderOverlay(viewport);
  if (this->GridSurfaceActor->GetVisibility())
  {
    count += this->GridSurfaceActor->RenderOverlay(viewport);
    count += this->ControlPolygonActor->RenderOverlay(viewport);
  }
  return count;
}

//-----------------------------------------------------------------------------
int vtkSlicerGridSurfaceRepresentation3D::RenderOpaqueGeometry(vtkViewport *viewport)
{
  int count=0;
  count = this->Superclass::RenderOpaqueGeometry(viewport);
  if (this->GridSurfaceActor->GetVisibility())
  {
    count += this->GridSurfaceActor->RenderOpaqueGeometry(viewport);
  }

  if (this->ControlPolygonActor->GetVisibility())
  {
    double diameter = (this->MarkupsDisplayNode->GetCurveLineSizeMode() == vtkMRMLMarkupsDisplayNode::UseLineDiameter ?
      this->MarkupsDisplayNode->GetLineDiameter() : this->ControlPointSize * this->MarkupsDisplayNode->GetLineThickness());
    this->ControlPolygonTubeFilter->SetRadius(diameter * 0.5);
    count += this->ControlPolygonActor->RenderOpaqueGeometry(viewport);
  }

  return count;
}

//-----------------------------------------------------------------------------
int vtkSlicerGridSurfaceRepresentation3D::RenderTranslucentPolygonalGeometry(
  vtkViewport *viewport)
{
  int count = this->Superclass::RenderTranslucentPolygonalGeometry(viewport);
  if (this->GridSurfaceActor->GetVisibility())
  {
    // The internal actor needs to share property keys.
    // This ensures the mapper state is consistent and allows depth peeling to work as expected.
    this->GridSurfaceActor->SetPropertyKeys(this->GetPropertyKeys());
    count += this->GridSurfaceActor->RenderTranslucentPolygonalGeometry(viewport);
  }

  if (this->ControlPolygonActor->GetVisibility())
  {
    // The internal actor needs to share property keys.
    // This ensures the mapper state is consistent and allows depth peeling to work as expected.
    this->ControlPolygonActor->SetPropertyKeys(this->GetPropertyKeys());
    count += this->ControlPolygonActor->RenderTranslucentPolygonalGeometry(viewport);
  }

  return count;
}

//-----------------------------------------------------------------------------
vtkTypeBool vtkSlicerGridSurfaceRepresentation3D::HasTranslucentPolygonalGeometry()
{
  if (this->Superclass::HasTranslucentPolygonalGeometry())
  {
    return true;
  }
  if (this->GridSurfaceActor->GetVisibility() && this->GridSurfaceActor->HasTranslucentPolygonalGeometry())
  {
    return true;
  }
  if (this->ControlPolygonActor->GetVisibility() && this->ControlPolygonActor->HasTranslucentPolygonalGeometry())
  {
    return true;
  }
  return false;
}

//----------------------------------------------------------------------
double* vtkSlicerGridSurfaceRepresentation3D::GetBounds()
{
  vtkBoundingBox boundingBox;
  const std::vector<vtkProp*> actors({this->GridSurfaceActor, this->ControlPolygonActor});
  this->AddActorsBounds(boundingBox, actors, Superclass::GetBounds());
  boundingBox.GetBounds(this->Bounds);
  return this->Bounds;
}

//----------------------------------------------------------------------
// void vtkSlicerGridSurfaceRepresentation3D::CanInteract(
//   vtkMRMLInteractionEventData* interactionEventData,
//   int &foundComponentType, int &foundComponentIndex, double &closestDistance2)
// {
//   foundComponentType = vtkMRMLMarkupsDisplayNode::ComponentNone;
//   vtkMRMLMarkupsNode* markupsNode = this->GetMarkupsNode();
//   if ( !markupsNode || markupsNode->GetLocked() || markupsNode->GetNumberOfDefinedControlPoints(true) < 1
//     || !interactionEventData )
//     {
//     return;
//     }
//   Superclass::CanInteract(interactionEventData, foundComponentType, foundComponentIndex, closestDistance2);
//   if (foundComponentType != vtkMRMLMarkupsDisplayNode::ComponentNone)
//     {
//     // if mouse is near a control point then select that (ignore the line)
//     return;
//     }

//   this->CanInteractWithBezierSurface(interactionEventData, foundComponentType, foundComponentIndex, closestDistance2);
// }

//-----------------------------------------------------------------------------
void vtkSlicerGridSurfaceRepresentation3D::PrintSelf(ostream& os, vtkIndent indent)
{
  //Superclass typedef defined in vtkTypeMacro() found in vtkSetGet.h
  this->Superclass::PrintSelf(os, indent);

  if (this->GridSurfaceActor)
  {
    os << indent << "GridSurface Visibility: " << this->GridSurfaceActor->GetVisibility() << "\n";
  }
  else
  {
    os << indent << "GridSurface Visibility: (none)\n";
  }

  if (this->ControlPolygonActor)
  {
    os << indent << "ControlPolygon Visibility: " << this->ControlPolygonActor->GetVisibility() << "\n";
  }
  else
  {
    os << indent << "ControlPolygon Visibility: (none)\n";
  }
}

//-----------------------------------------------------------------------------
void vtkSlicerGridSurfaceRepresentation3D::UpdateInterpolatorConnection()
{
  if (!this->GridSurfaceControlPointSet->GetNumberOfPoints())
  {
    this->GridSurfaceNormals->SetInputConnection(nullptr);
    this->GridSurfaceMapper->SetInputConnection(nullptr);
    this->GridSurfaceActor->SetMapper(nullptr);
    return;
  }

  vtkMRMLMarkupsGridSurfaceNode* gridSurfaceNode = vtkMRMLMarkupsGridSurfaceNode::SafeDownCast(this->MarkupsNode);
  if (!gridSurfaceNode)
  {
    return;
  }

  switch (gridSurfaceNode->GetGridSurfaceType())
  {
  case vtkMRMLMarkupsGridSurfaceNode::NURBS:
  {
    vtkMRMLModelNode* outputSurfaceModelNode = gridSurfaceNode->GetOutputSurfaceModelNode();
    if (!outputSurfaceModelNode)
    {
      this->GridSurfaceNormals->SetInputConnection(this->NurbsSurfaceSource->GetOutputPort());
      this->GridSurfaceMapper->SetInputConnection(this->GridSurfaceNormals->GetOutputPort());
      this->GridSurfaceActor->SetMapper(this->GridSurfaceMapper);
    }
    else
    {
      this->GridSurfaceNormals->SetInputConnection(this->NurbsSurfaceSource->GetOutputPort());
      this->GridSurfaceMapper->SetInputConnection(nullptr);
      this->GridSurfaceActor->SetMapper(nullptr);
      // Connect NURBS surface output to the surface model node
      outputSurfaceModelNode->SetPolyDataConnection(this->GridSurfaceNormals->GetOutputPort());
    }
    break;
  }
  case vtkMRMLMarkupsGridSurfaceNode::Bezier:
  {
    this->GridSurfaceNormals->SetInputConnection(this->BezierSurfaceSource->GetOutputPort());
    this->GridSurfaceMapper->SetInputConnection(this->GridSurfaceNormals->GetOutputPort());
    this->GridSurfaceActor->SetMapper(this->GridSurfaceMapper);
    break;
  }
  default:
    vtkErrorMacro("UpdateInterpolatorConnection: Invalid interpolator type");
  }
}

//-----------------------------------------------------------------------------
 void vtkSlicerGridSurfaceRepresentation3D::UpdateInteractionPipeline()
 {
   vtkMRMLMarkupsGridSurfaceNode* gridSurfaceNode = vtkMRMLMarkupsGridSurfaceNode::SafeDownCast(this->MarkupsNode);

   if (gridSurfaceNode)
   {
     int gridResolution[2] ={0};
     gridSurfaceNode->GetGridResolution(gridResolution);
     if (gridSurfaceNode->GetNumberOfDefinedControlPoints(true) < gridResolution[0] * gridResolution[1])
     {
       this->InteractionPipeline->Actor->SetVisibility(false);
       return;
     }
   }
   // Final visibility handled by superclass in vtkSlicerMarkupsWidgetRepresentation
   Superclass::UpdateInteractionPipeline();
 }

//-----------------------------------------------------------------------------
void vtkSlicerGridSurfaceRepresentation3D::UpdateGridSurface(vtkMRMLMarkupsGridSurfaceNode* node)
{
  if (!node)
  {
    return;
  }

  int gridResolution[2] = {0};
  node->GetGridResolution(gridResolution);
  if ( node->GetNumberOfDefinedControlPoints() < gridResolution[0] * gridResolution[1]
    && this->GridSurfaceControlPointSet->GetNumberOfPoints() == 0 )
  {
    // Do not calculate surface while plane control points are being defined
    return;
  }

  // Change internal grid if grid resolution has changed
  if ( gridResolution[0] != this->NurbsSurfaceSource->GetInputResolution()[0]
    || gridResolution[1] != this->NurbsSurfaceSource->GetInputResolution()[1]
    || this->GridSurfaceControlPointSet->GetNumberOfPoints() == 0 ) // If grid was just defined
  {
    this->InitializeGridSurfaceControlPoints(gridResolution[0], gridResolution[1]);

    this->BezierSurfaceSource->SetNumberOfControlPoints(gridResolution[0], gridResolution[1]);
    this->BezierSurfaceSource->SetResolution(gridResolution[0] * 3, gridResolution[1] * 3);

    this->NurbsSurfaceSource->SetInputResolution(gridResolution[0], gridResolution[1]);
    // Ensure appropriate interpolation degrees for every grid resolution
    this->NurbsSurfaceSource->SetInterpolationDegrees(
      gridResolution[0] == 3 ? 2 : 3, gridResolution[1] == 3 ? 2 : 3);

    if (node->GetNumberOfControlPoints() != gridResolution[0] * gridResolution[1])
    {
      // Redefine control points to new resolution if needed
      node->ResampleToNewGridResolution();
    }
  }

  // Set supported parameters from MRML node
  this->NurbsSurfaceSource->SetExpansionFactor(node->GetExpansionFactor());
  this->NurbsSurfaceSource->SetWrapAround(node->GetWrapAround());
  this->NurbsSurfaceSource->SetDelta(1.0 / node->GetSamplingResolution());
  this->NurbsSurfaceSource->SetIterativeParameterSpaceCalculation(node->GetIterativeParameterSpaceCalculation());
  this->NurbsSurfaceSource->SetGenerateQuadMesh(node->GetGenerateQuadMesh());

  // Set markup control points to the surface source
  //TODO: Revisit if this is needed or only in case of the Bezier source
  vtkNew<vtkPoints> controlPoints;
  controlPoints->SetNumberOfPoints(gridResolution[0] * gridResolution[1]);
  if (node->GetNumberOfControlPoints() == gridResolution[0] * gridResolution[1])
  {
    for (int i=0; i<gridResolution[0]; ++i)
    {
      for (int j=0; j<gridResolution[1]; ++j)
      {
        // Transpose point array for Bezier source
        double point[3] = {0.0};
        node->GetNthControlPointPosition(j*gridResolution[0] + i, point);
        controlPoints->SetPoint(i*gridResolution[1] + j,
          static_cast<float>(point[0]),
          static_cast<float>(point[1]),
          static_cast<float>(point[2]));
      }
    }
    this->GridSurfaceControlPointSet->SetPoints(controlPoints);

    // Set control points to surface source
    switch (node->GetGridSurfaceType())
    {
    case vtkMRMLMarkupsGridSurfaceNode::NURBS:
      this->NurbsSurfaceSource->SetInputData(0, this->GridSurfaceControlPointSet);
      break;
    case vtkMRMLMarkupsGridSurfaceNode::Bezier:
      this->BezierSurfaceSource->SetControlPoints(controlPoints);
      break;
    default:
      vtkErrorMacro("UpdateGridSurface: Invalid interpolator type");
    }
  }

  this->UpdateInterpolatorConnection();
}

//-----------------------------------------------------------------------------
void vtkSlicerGridSurfaceRepresentation3D::UpdateControlPolygon(vtkMRMLMarkupsGridSurfaceNode* node)
{
  int gridResolution[2] = {0};
  node->GetGridResolution(gridResolution);
  if (node->GetNumberOfControlPoints() == gridResolution[0] * gridResolution[1])
  {
    // Generate topology
    vtkSmartPointer<vtkCellArray> planeCells = vtkSmartPointer<vtkCellArray>::New();
    for (int i=0; i<gridResolution[0]-1; ++i)
    {
      for (int j=0; j<gridResolution[1]-1; ++j)
      {
        vtkNew<vtkPolyLine> polyLine;
        polyLine->GetPointIds()->SetNumberOfIds(5);
        polyLine->GetPointIds()->SetId(0, i*gridResolution[1] + j);
        polyLine->GetPointIds()->SetId(1, i*gridResolution[1] + j+1);
        polyLine->GetPointIds()->SetId(2, (i+1)*gridResolution[1] + j+1);
        polyLine->GetPointIds()->SetId(3, (i+1)*gridResolution[1] + j);
        polyLine->GetPointIds()->SetId(4, i*gridResolution[1] + j);
        planeCells->InsertNextCell(polyLine);
      }
    }
    // Draw closing lines when wrap around is enabled
    if (node->GetWrapAround() == vtkMRMLMarkupsGridSurfaceNode::AlongU)
    {
      for (int j=0; j<gridResolution[1]; ++j)
      {
        vtkNew<vtkPolyLine> polyLine;
        polyLine->GetPointIds()->SetNumberOfIds(2);
        polyLine->GetPointIds()->SetId(0, j);
        polyLine->GetPointIds()->SetId(1, (gridResolution[0]-1)*gridResolution[1]+j);
        planeCells->InsertNextCell(polyLine);
      }
    }
    else if (node->GetWrapAround() == vtkMRMLMarkupsGridSurfaceNode::AlongV)
    {
      for (int i=0; i<gridResolution[0]; ++i)
      {
        vtkNew<vtkPolyLine> polyLine;
        polyLine->GetPointIds()->SetNumberOfIds(2);
        polyLine->GetPointIds()->SetId(0, i*gridResolution[1]);
        polyLine->GetPointIds()->SetId(1, (i+1)*gridResolution[1]-1);
        planeCells->InsertNextCell(polyLine);
      }
    }

    this->ControlPolygonPolyData->SetPoints(this->GridSurfaceControlPointSet->GetPoints());
    this->ControlPolygonPolyData->SetLines(planeCells);
  }
}

//-----------------------------------------------------------------------------
void vtkSlicerGridSurfaceRepresentation3D::InitializeGridSurfaceControlPoints(int resX, int resY)
{
  // Fill the initial point array of the Bezier surface with the appropriate number of points
  auto planeSource = vtkSmartPointer<vtkPlaneSource>::New();
  planeSource->SetResolution(resX - 1, resY - 1);
  planeSource->Update();

  this->GridSurfaceControlPointSet->SetPoints(planeSource->GetOutput()->GetPoints());
}
