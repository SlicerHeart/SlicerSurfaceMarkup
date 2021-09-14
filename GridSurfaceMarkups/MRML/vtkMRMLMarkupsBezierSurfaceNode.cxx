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

#include "vtkMRMLMarkupsBezierSurfaceNode.h"

// MRML includes
#include <vtkMRMLScene.h>
#include <vtkMRMLTransformNode.h>

// VTK includes
#include <vtkCallbackCommand.h>
#include <vtkCommand.h>
#include <vtkMath.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPoints.h>

const int NUMBER_OF_PLANE_CONTROL_POINTS = 3; // 3 points used for initial plane definition, then filled with the rest

//--------------------------------------------------------------------------------
vtkMRMLNodeNewMacro(vtkMRMLMarkupsBezierSurfaceNode);

//--------------------------------------------------------------------------------
vtkMRMLMarkupsBezierSurfaceNode::vtkMRMLMarkupsBezierSurfaceNode()
  :Superclass()
{
  this->RequiredNumberOfControlPoints = NUMBER_OF_PLANE_CONTROL_POINTS;
  this->MaximumNumberOfControlPoints = 0;

  this->CurveInputPoly->GetPoints()->AddObserver(vtkCommand::ModifiedEvent, this->MRMLCallbackCommand);

  this->AddObserver(vtkMRMLMarkupsNode::PointPositionDefinedEvent, this->MRMLCallbackCommand);
}

//----------------------------------------------------------------------------
void vtkMRMLMarkupsBezierSurfaceNode::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf(os,indent);
}

//---------------------------------------------------------------------------
void vtkMRMLMarkupsBezierSurfaceNode::ProcessMRMLEvents(vtkObject* caller, unsigned long event, void* callData)
{
  if (caller == this->CurveInputPoly->GetPoints() || caller == this->GetParentTransformNode())
  {
    this->UpdateSurfaceFromControlPoints();
    //this->UpdateObjectToWorldMatrix();
  }
  else if (caller == this && event == vtkMRMLMarkupsNode::PointPositionDefinedEvent)
  {
    this->UpdateControlPointsFromSurface();
  }
  Superclass::ProcessMRMLEvents(caller, event, callData);
}

//----------------------------------------------------------------------------
void vtkMRMLMarkupsBezierSurfaceNode::SetGridSize(const int gridSize[2])
{
  this->SetGridSize(gridSize[0], gridSize[1]);
}

//----------------------------------------------------------------------------
void vtkMRMLMarkupsBezierSurfaceNode::SetGridSize(int a, int b)
{
  if (this->GridSize[0] == a && this->GridSize[1] == b)
  {
    return;
  }

  MRMLNodeModifyBlocker blocker(this);
  this->GridSize[0] = a;
  this->GridSize[1] = b;
  this->UpdateControlPointsFromSurface();
  this->Modified();
}

//---------------------------------------------------------------------------
void vtkMRMLMarkupsBezierSurfaceNode::UpdateSurfaceFromControlPoints()
{
  if (this->IsUpdatingControlPointsFromSurface || this->IsUpdatingSurfaceFromControlPoints)
  {
    return;
  }

  this->IsUpdatingSurfaceFromControlPoints = true;

  // Block events in this scope
  MRMLNodeModifyBlocker blocker(this);

  //TODO:

  this->IsUpdatingSurfaceFromControlPoints = false;
}

//----------------------------------------------------------------------------
void vtkMRMLMarkupsBezierSurfaceNode::UpdateControlPointsFromSurface()
{
  if (this->IsUpdatingControlPointsFromSurface || this->IsUpdatingSurfaceFromControlPoints)
  {
    return;
  }

  this->IsUpdatingControlPointsFromSurface = true;

  // Block events in this scope
  MRMLNodeModifyBlocker blocker(this);

  int numberOfControlPoints = this->GetNumberOfControlPoints();
  if (numberOfControlPoints == NUMBER_OF_PLANE_CONTROL_POINTS)
  {
    // Auto-fill control points for the surface if the first three points defining
    // the plane from scratch have been defined.

    //
    // Calculate fourth point of plane rectangle from first three control points:
    //   Get vector between second and third point and add it to first point.
    //
    //     0----3
    //     |    |
    //     1----2
    //
    double position_0[3] = {0.0};
    this->GetNthControlPointPositionWorld(0, position_0);
    double position_1[3] = {0.0};
    this->GetNthControlPointPositionWorld(1, position_1);
    double position_2[3] = {0.0};
    this->GetNthControlPointPositionWorld(2, position_2);
    double vector_1_2[3] = {0.0};
    vtkMath::Subtract(position_2, position_1, vector_1_2);
    double position_3[3] = {0.0};
    vtkMath::Add(position_0, vector_1_2, position_3);

    //
    // Calculate grid points based on specified grid size
    //

    // Vectors in both directions for displacement between adjacent control points
    double vector_a[3] = {0.0}; // = vector_0_1 / GridSize[0]
    vtkMath::Subtract(position_1, position_0, vector_a);
    vtkMath::MultiplyScalar(vector_a, 1.0 / (double)(this->GridSize[0] - 1));
    double vector_b[3] = { vector_1_2[0], vector_1_2[1], vector_1_2[2] }; // = vector_1_2 / GridSize[0]
    vtkMath::MultiplyScalar(vector_b, 1.0 / (double)(this->GridSize[1] - 1));

    // Fill in control point list
    vtkNew<vtkPoints> controlPoints;
    for (int b=0; b<this->GridSize[1]; ++b)
    {
      double vector_b_scaled[3] ={vector_b[0], vector_b[1], vector_b[2]};
      vtkMath::MultiplyScalar(vector_b_scaled, b);

      for (int a=0; a<this->GridSize[0]; ++a)
      {
        double vector_a_scaled[3] ={vector_a[0], vector_a[1], vector_a[2]};
        vtkMath::MultiplyScalar(vector_a_scaled, a);

        double position_a_b[3] ={0.0};
        vtkMath::Add(position_0, vector_a_scaled, position_a_b);
        vtkMath::Add(position_a_b, vector_b_scaled, position_a_b);

        controlPoints->InsertNextPoint(position_a_b);
      }
    }

    this->SetControlPointPositionsWorld(controlPoints);
  }
  else if (numberOfControlPoints > NUMBER_OF_PLANE_CONTROL_POINTS
    && numberOfControlPoints != this->GridSize[0] * this->GridSize[1])
  {
    // Re-generate new control point positions from current control points
    // after the grid size has been changed.

    //TODO:
    vtkNew<vtkPoints> pointsWorld;
    //pointsWorld->InsertNextPoint(center_World);
    //this->SetControlPointPositionsWorld(pointsWorld);
    this->MaximumNumberOfControlPoints = this->GridSize[0] * this->GridSize[1];
    this->RequiredNumberOfControlPoints = this->GridSize[0] * this->GridSize[1];
  }

  this->IsUpdatingControlPointsFromSurface = false;
}
