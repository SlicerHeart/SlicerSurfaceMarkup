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
#include "vtkBezierSurfaceSource.h"

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
  this->BezierSurfaceSource = vtkSmartPointer<vtkBezierSurfaceSource>::New();

  this->BezierSurfaceNormals = vtkSmartPointer<vtkPolyDataNormals>::New();
  this->BezierSurfaceNormals->SetInputConnection(this->BezierSurfaceSource->GetOutputPort());

  this->InitializeBezierSurfaceControlPoints(4, 4);

  this->BezierSurfaceMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  this->BezierSurfaceMapper->SetInputConnection(this->BezierSurfaceNormals->GetOutputPort());
  this->BezierSurfaceActor = vtkSmartPointer<vtkActor>::New();
  this->BezierSurfaceActor->SetMapper(this->BezierSurfaceMapper);

  this->ControlPolygonPolyData = vtkSmartPointer<vtkPolyData>::New();
  this->ControlPolygonTubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
  this->ControlPolygonTubeFilter->SetInputData(this->ControlPolygonPolyData.GetPointer());
  this->ControlPolygonTubeFilter->SetRadius(1);
  this->ControlPolygonTubeFilter->SetNumberOfSides(20);

  this->ControlPolygonMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  this->ControlPolygonMapper->SetInputConnection(this->ControlPolygonTubeFilter->GetOutputPort());

  this->ControlPolygonActor = vtkSmartPointer<vtkActor>::New();
  this->ControlPolygonActor->SetMapper(this->ControlPolygonMapper);
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

  this->UpdateBezierSurface(gridSurfaceMarkupsNode);
  this->UpdateControlPolygon(gridSurfaceMarkupsNode);

  double diameter = (this->MarkupsDisplayNode->GetCurveLineSizeMode() == vtkMRMLMarkupsDisplayNode::UseLineDiameter ?
    this->MarkupsDisplayNode->GetLineDiameter() : this->ControlPointSize * this->MarkupsDisplayNode->GetLineThickness());
  this->ControlPolygonTubeFilter->SetRadius(diameter * 0.5);

  int controlPointType = Active;
  if (this->MarkupsDisplayNode->GetActiveComponentType() != vtkMRMLMarkupsDisplayNode::ComponentLine)
  {
    controlPointType = this->GetAllControlPointsSelected() ? Selected : Unselected;
  }
  this->ControlPolygonActor->SetProperty(this->GetControlPointsPipeline(controlPointType)->Property);

  this->NeedToRenderOn();
}

//----------------------------------------------------------------------
void vtkSlicerGridSurfaceRepresentation3D::GetActors(vtkPropCollection *pc)
{
  this->Superclass::GetActors(pc);
  this->BezierSurfaceActor->GetActors(pc);
  this->ControlPolygonActor->GetActors(pc);
}

//----------------------------------------------------------------------
void vtkSlicerGridSurfaceRepresentation3D::ReleaseGraphicsResources(vtkWindow *win)
{
  this->Superclass::ReleaseGraphicsResources(win);
  this->BezierSurfaceActor->ReleaseGraphicsResources(win);
  this->ControlPolygonActor->ReleaseGraphicsResources(win);
}

//----------------------------------------------------------------------
int vtkSlicerGridSurfaceRepresentation3D::RenderOverlay(vtkViewport *viewport)
{
  int count = this->Superclass::RenderOverlay(viewport);
  if (this->BezierSurfaceActor->GetVisibility())
  {
    count += this->BezierSurfaceActor->RenderOverlay(viewport);
    count += this->ControlPolygonActor->RenderOverlay(viewport);
  }
  return count;
}

//-----------------------------------------------------------------------------
int vtkSlicerGridSurfaceRepresentation3D::RenderOpaqueGeometry(vtkViewport *viewport)
{
  int count=0;
  count = this->Superclass::RenderOpaqueGeometry(viewport);
  if (this->BezierSurfaceActor->GetVisibility())
  {
    count += this->BezierSurfaceActor->RenderOpaqueGeometry(viewport);
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
  if (this->BezierSurfaceActor->GetVisibility())
  {
    // The internal actor needs to share property keys.
    // This ensures the mapper state is consistent and allows depth peeling to work as expected.
    this->BezierSurfaceActor->SetPropertyKeys(this->GetPropertyKeys());
    count += this->BezierSurfaceActor->RenderTranslucentPolygonalGeometry(viewport);
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
  if (this->BezierSurfaceActor->GetVisibility() && this->BezierSurfaceActor->HasTranslucentPolygonalGeometry())
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
  const std::vector<vtkProp*> actors({this->BezierSurfaceActor, this->ControlPolygonActor});
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

  if (this->BezierSurfaceActor)
  {
    os << indent << "BezierSurface Visibility: " << this->BezierSurfaceActor->GetVisibility() << "\n";
  }
  else
  {
    os << indent << "BezierSurface Visibility: (none)\n";
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
// void vtkSlicerGridSurfaceRepresentation3D::UpdateInteractionPipeline()
// {
//   if (!this->MarkupsNode || this->MarkupsNode->GetNumberOfDefinedControlPoints(true) < 16) //TODO:
//     {
//     this->InteractionPipeline->Actor->SetVisibility(false);
//     return;
//     }
//   // Final visibility handled by superclass in vtkSlicerMarkupsWidgetRepresentation
//   Superclass::UpdateInteractionPipeline();
// }

//-----------------------------------------------------------------------------
void vtkSlicerGridSurfaceRepresentation3D::UpdateBezierSurface(vtkMRMLMarkupsGridSurfaceNode *node)
{
  if (!node)
  {
    return;
  }

  // Change internal grid if grid resolution has changed
  int gridResolution[2] = {0};
  node->GetGridResolution(gridResolution);
  if ( gridResolution[0] != this->BezierSurfaceSource->GetNumberOfControlPointsX()
    || gridResolution[1] != this->BezierSurfaceSource->GetNumberOfControlPointsY() )
  {
    this->InitializeBezierSurfaceControlPoints(gridResolution[0], gridResolution[1]);
    this->BezierSurfaceSource->SetNumberOfControlPoints(gridResolution[0], gridResolution[1]);
    this->BezierSurfaceSource->SetResolution(gridResolution[0] * 3, gridResolution[1] * 3);

    node->ResampleToNewGridResolution();
  }

  if (node->GetNumberOfControlPoints() == gridResolution[0] * gridResolution[1])
  {
    for (int i=0; i<gridResolution[0]; ++i)
    {
      for (int j=0; j<gridResolution[1]; ++j)
      {
        // Transpose point array for Bezier source
        double point[3] = {0.0};
        node->GetNthControlPointPosition(j*gridResolution[0] + i, point);
        this->BezierSurfaceControlPoints->SetPoint(i*gridResolution[1] + j,
          static_cast<float>(point[0]),
          static_cast<float>(point[1]),
          static_cast<float>(point[2]));
      }
    }

    this->BezierSurfaceSource->SetControlPoints(this->BezierSurfaceControlPoints);
  }
}

//-----------------------------------------------------------------------------
void vtkSlicerGridSurfaceRepresentation3D::UpdateControlPolygon(vtkMRMLMarkupsGridSurfaceNode *node)
{
  int gridResolution[2] = {0};
  node->GetGridResolution(gridResolution);
  if (node->GetNumberOfControlPoints() == gridResolution[0] * gridResolution[1])
  {
    // Generate topology;
    vtkSmartPointer<vtkCellArray> planeCells = vtkSmartPointer<vtkCellArray>::New();
    for (int i=0; i<gridResolution[0]-1; ++i)
    {
      for (int j=0; j<gridResolution[1]-1; ++j)
      {
        vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
        polyLine->GetPointIds()->SetNumberOfIds(5);
        polyLine->GetPointIds()->SetId(0, i*gridResolution[1] + j);
        polyLine->GetPointIds()->SetId(1, i*gridResolution[1] + j+1);
        polyLine->GetPointIds()->SetId(2, (i+1)*gridResolution[1] + j+1);
        polyLine->GetPointIds()->SetId(3, (i+1)*gridResolution[1] + j);
        polyLine->GetPointIds()->SetId(4, i*gridResolution[1] + j);
        planeCells->InsertNextCell(polyLine);
      }
    }

    this->ControlPolygonPolyData->SetPoints(this->BezierSurfaceControlPoints);
    this->ControlPolygonPolyData->SetLines(planeCells);
  }
}

//-----------------------------------------------------------------------------
void vtkSlicerGridSurfaceRepresentation3D::InitializeBezierSurfaceControlPoints(int resX, int resY)
{
  // Fill the initial point array of the Bezier surface with the appropriate number of points
  auto planeSource = vtkSmartPointer<vtkPlaneSource>::New();
  planeSource->SetResolution(resX - 1, resY - 1);
  planeSource->Update();

  if (!this->BezierSurfaceControlPoints.GetPointer())
  {
    this->BezierSurfaceControlPoints = vtkSmartPointer<vtkPoints>::New();
  }
  this->BezierSurfaceControlPoints->SetNumberOfPoints(resX * resY);
  this->BezierSurfaceControlPoints->DeepCopy(planeSource->GetOutput()->GetPoints());
}
