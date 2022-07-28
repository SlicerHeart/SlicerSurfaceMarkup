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

#include "vtkMRMLMarkupsGridSurfaceNode.h"

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

const char* vtkMRMLMarkupsGridSurfaceNode::OutputSurfaceModelNodeReferenceRole = "outputSurfaceModel";
const char* vtkMRMLMarkupsGridSurfaceNode::OutputSurfaceModelNodeReferenceMRMLAttributeName = "outputSurfaceModelNodeRef";

//--------------------------------------------------------------------------------
vtkMRMLNodeNewMacro(vtkMRMLMarkupsGridSurfaceNode);

//--------------------------------------------------------------------------------
vtkMRMLMarkupsGridSurfaceNode::vtkMRMLMarkupsGridSurfaceNode()
  :Superclass()
{
  this->RequiredNumberOfControlPoints = NUMBER_OF_PLANE_CONTROL_POINTS;
  this->MaximumNumberOfControlPoints = -1;

  this->CurveInputPoly->GetPoints()->AddObserver(vtkCommand::ModifiedEvent, this->MRMLCallbackCommand);

  this->AddObserver(vtkMRMLMarkupsNode::PointPositionDefinedEvent, this->MRMLCallbackCommand);

  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkCommand::ModifiedEvent);
  this->AddNodeReferenceRole(
    this->GetOutputSurfaceModelNodeReferenceRole(), this->GetOutputSurfaceModelNodeReferenceMRMLAttributeName(), events.GetPointer());
}

//----------------------------------------------------------------------------
void vtkMRMLMarkupsGridSurfaceNode::WriteXML(ostream& of, int nIndent)
{
  Superclass::WriteXML(of,nIndent);

  vtkMRMLWriteXMLBeginMacro(of);
  vtkMRMLWriteXMLVectorMacro(gridResolution, GridResolution, int, 2);
  vtkMRMLWriteXMLFloatMacro(expansionFactor, ExpansionFactor);
  vtkMRMLWriteXMLEnumMacro(wrapAround, WrapAround);
  vtkMRMLWriteXMLFloatMacro(samplingResolution, SamplingResolution);
  vtkMRMLWriteXMLBooleanMacro(iterativeParameterSpaceCalculation, IterativeParameterSpaceCalculation);
  vtkMRMLWriteXMLBooleanMacro(generateQuadMesh, GenerateQuadMesh);
  vtkMRMLWriteXMLEndMacro();
}

//----------------------------------------------------------------------------
void vtkMRMLMarkupsGridSurfaceNode::ReadXMLAttributes(const char** atts)
{
  MRMLNodeModifyBlocker blocker(this);

  Superclass::ReadXMLAttributes(atts);

  vtkMRMLReadXMLBeginMacro(atts);
  vtkMRMLReadXMLVectorMacro(gridResolution, GridResolution, int, 2);
  vtkMRMLReadXMLFloatMacro(expansionFactor, ExpansionFactor);
  vtkMRMLReadXMLEnumMacro(wrapAround, WrapAround);
  vtkMRMLReadXMLFloatMacro(samplingResolution, SamplingResolution);
  vtkMRMLReadXMLBooleanMacro(iterativeParameterSpaceCalculation, IterativeParameterSpaceCalculation);
  vtkMRMLReadXMLBooleanMacro(generateQuadMesh, GenerateQuadMesh);
  vtkMRMLReadXMLEndMacro();

  // In scenes created by Slicer version version 4.13.0 revision 30287 (built 2021-10-05).
  // The value used to represent unlimited control points has been changed to -1.
  if (this->MaximumNumberOfControlPoints == 0)
    {
    this->MaximumNumberOfControlPoints = -1;
    }
}

//----------------------------------------------------------------------------
void vtkMRMLMarkupsGridSurfaceNode::CopyContent(vtkMRMLNode* anode, bool deepCopy/*=true*/)
{
  MRMLNodeModifyBlocker blocker(this);
  Superclass::CopyContent(anode, deepCopy);

  vtkMRMLCopyBeginMacro(anode);
  vtkMRMLCopyVectorMacro(GridResolution, int, 2);
  vtkMRMLCopyFloatMacro(ExpansionFactor);
  vtkMRMLCopyEnumMacro(WrapAround);
  vtkMRMLCopyFloatMacro(SamplingResolution);
  vtkMRMLCopyBooleanMacro(IterativeParameterSpaceCalculation);
  vtkMRMLCopyBooleanMacro(GenerateQuadMesh);
  vtkMRMLCopyEndMacro();
}

//----------------------------------------------------------------------------
void vtkMRMLMarkupsGridSurfaceNode::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf(os,indent);

  vtkMRMLPrintBeginMacro(os,indent);
  vtkMRMLPrintVectorMacro(GridResolution, int, 2);
  vtkMRMLPrintFloatMacro(ExpansionFactor);
  vtkMRMLPrintEnumMacro(WrapAround);
  vtkMRMLPrintFloatMacro(SamplingResolution);
  vtkMRMLPrintBooleanMacro(IterativeParameterSpaceCalculation);
  vtkMRMLPrintBooleanMacro(GenerateQuadMesh);
  vtkMRMLPrintEndMacro();
}

//---------------------------------------------------------------------------
void vtkMRMLMarkupsGridSurfaceNode::ProcessMRMLEvents(vtkObject* caller, unsigned long event, void* callData)
{
  if (caller == this->CurveInputPoly->GetPoints() || caller == this->GetParentTransformNode())
  {
    this->UpdateGridSurfaceFromControlPoints();
    //this->UpdateObjectToWorldMatrix();
  }
  else if (caller == this && event == vtkMRMLMarkupsNode::PointPositionDefinedEvent)
  {
    this->UpdateControlPointsFromGridSurface();
  }
  Superclass::ProcessMRMLEvents(caller, event, callData);
}

//----------------------------------------------------------------------------
void vtkMRMLMarkupsGridSurfaceNode::SetGridSurfaceType(int gridSurfaceType)
{
  if (this->GridSurfaceType == gridSurfaceType)
  {
    return;
  }

  this->GridSurfaceType = gridSurfaceType;

  this->UpdateGridSurfaceFromControlPoints();
  this->Modified();
}

//----------------------------------------------------------------------------
const char* vtkMRMLMarkupsGridSurfaceNode::GetGridSurfaceTypeAsString(int gridSurfaceType)
{
  switch (gridSurfaceType)
  {
  case vtkMRMLMarkupsGridSurfaceNode::NURBS:
    return "NURBS";
  case vtkMRMLMarkupsGridSurfaceNode::Bezier:
    return "Bezier";
  //case vtkMRMLMarkupsGridSurfaceNode::ThinPlate:
  //  return "ThinPlate";
  default:
    // invalid id
    return "Invalid";
  }
}

//-----------------------------------------------------------
int vtkMRMLMarkupsGridSurfaceNode::GetGridSurfaceTypeFromString(const char* name)
{
  if (name == nullptr)
  {
    // invalid name
    return -1;
  }
  for (int i = 0; i < vtkMRMLMarkupsGridSurfaceNode::GridSurfaceType_Last; i++)
  {
    if (strcmp(name, vtkMRMLMarkupsGridSurfaceNode::GetGridSurfaceTypeAsString(i)) == 0)
    {
      // found a matching name
      return i;
    }
  }
  // unknown name
  return -1;
}

//----------------------------------------------------------------------------
void vtkMRMLMarkupsGridSurfaceNode::SetWrapAround(int wrapAround)
{
  if (this->WrapAround == wrapAround)
  {
    return;
  }

  this->WrapAround = wrapAround;

  this->UpdateGridSurfaceFromControlPoints();
  this->Modified();
}

//----------------------------------------------------------------------------
const char* vtkMRMLMarkupsGridSurfaceNode::GetWrapAroundAsString(int wrapAround)
{
  switch (wrapAround)
  {
  case vtkMRMLMarkupsGridSurfaceNode::NoWrap:
    return "No wrap";
  case vtkMRMLMarkupsGridSurfaceNode::AlongU:
    return "Along u axis";
  case vtkMRMLMarkupsGridSurfaceNode::AlongV:
    return "Along v axis";
  default:
    // invalid id
    return "Invalid";
  }
}

//-----------------------------------------------------------
int vtkMRMLMarkupsGridSurfaceNode::GetWrapAroundFromString(const char* name)
{
  if (name == nullptr)
  {
    // invalid name
    return -1;
  }
  for (int i = 0; i < vtkMRMLMarkupsGridSurfaceNode::WrapAround_Last; i++)
  {
    if (strcmp(name, vtkMRMLMarkupsGridSurfaceNode::GetWrapAroundAsString(i)) == 0)
    {
      // found a matching name
      return i;
    }
  }
  // unknown name
  return -1;
}

//----------------------------------------------------------------------------
void vtkMRMLMarkupsGridSurfaceNode::SetGridResolution(const int gridResolution[2])
{
  this->SetGridResolution(gridResolution[0], gridResolution[1]);
}

//----------------------------------------------------------------------------
void vtkMRMLMarkupsGridSurfaceNode::SetGridResolution(int x, int y)
{
  if (this->GridResolution[0] == x && this->GridResolution[1] == y)
  {
    return;
  }

  this->PreviousGridResolution[0] = this->GridResolution[0];
  this->PreviousGridResolution[1] = this->GridResolution[1];

  MRMLNodeModifyBlocker blocker(this);
  this->GridResolution[0] = x;
  this->GridResolution[1] = y;

  // Modified event triggers update in the representation, which then calls ResampleToNewGridResolution
  this->Modified();
}

//----------------------------------------------------------------------------
const char* vtkMRMLMarkupsGridSurfaceNode::GetOutputSurfaceModelNodeReferenceMRMLAttributeName()
{
  return vtkMRMLMarkupsGridSurfaceNode::OutputSurfaceModelNodeReferenceMRMLAttributeName;
}

//----------------------------------------------------------------------------
const char* vtkMRMLMarkupsGridSurfaceNode::GetOutputSurfaceModelNodeReferenceRole()
{
  return vtkMRMLMarkupsGridSurfaceNode::OutputSurfaceModelNodeReferenceRole;
}

//----------------------------------------------------------------------------
void vtkMRMLMarkupsGridSurfaceNode::SetOutputSurfaceModelNodeID(const char* modelNodeId)
{
  // Disconnect previous output surface model node
  vtkMRMLModelNode* currentOutputSurfaceModelNode = this->GetOutputSurfaceModelNode();
  if (currentOutputSurfaceModelNode)
  {
    currentOutputSurfaceModelNode->SetDisplayVisibility(false);
    currentOutputSurfaceModelNode->SetPolyDataConnection(nullptr);
  }

  this->SetAndObserveNodeReferenceID(this->GetOutputSurfaceModelNodeReferenceRole(), modelNodeId);

  vtkMRMLModelNode* newOutputSurfaceModelNode = this->GetOutputSurfaceModelNode();
  if (newOutputSurfaceModelNode)
  {
    // Make sure new output surface is visible
    newOutputSurfaceModelNode->CreateDefaultDisplayNodes();
    newOutputSurfaceModelNode->SetDisplayVisibility(true);
    // Prevent control point snapping on the output surface by default
    newOutputSurfaceModelNode->SetSelectable(false);
  }
}

//----------------------------------------------------------------------------
const char* vtkMRMLMarkupsGridSurfaceNode::GetOutputSurfaceModelNodeID()
{
  return this->GetNodeReferenceID(this->GetOutputSurfaceModelNodeReferenceRole());
}

//----------------------------------------------------------------------------
vtkMRMLModelNode* vtkMRMLMarkupsGridSurfaceNode::GetOutputSurfaceModelNode()
{
  return vtkMRMLModelNode::SafeDownCast(this->GetNodeReference(this->GetOutputSurfaceModelNodeReferenceRole()));
}

//---------------------------------------------------------------------------
void vtkMRMLMarkupsGridSurfaceNode::UpdateGridSurfaceFromControlPoints()
{
  if (this->IsUpdatingControlPointsFromGridSurface || this->IsUpdatingGridSurfaceFromControlPoints)
  {
    return;
  }

  this->IsUpdatingGridSurfaceFromControlPoints = true;

  // Block events in this scope
  MRMLNodeModifyBlocker blocker(this);

  //TODO:

  this->IsUpdatingGridSurfaceFromControlPoints = false;
}

//----------------------------------------------------------------------------
void vtkMRMLMarkupsGridSurfaceNode::UpdateControlPointsFromGridSurface()
{
  if (this->IsUpdatingControlPointsFromGridSurface || this->IsUpdatingGridSurfaceFromControlPoints)
  {
    return;
  }

  this->IsUpdatingControlPointsFromGridSurface = true;

  // Block events in this scope
  MRMLNodeModifyBlocker blocker(this);

  int numberOfControlPoints = this->GetNumberOfControlPoints();
  if (numberOfControlPoints == NUMBER_OF_PLANE_CONTROL_POINTS)
  {
    // Auto-fill control points for the surface if the first three points defining
    // the plane from scratch have been defined.
    double position_0[3] = {0.0};
    this->GetNthControlPointPositionWorld(0, position_0);
    double position_1[3] = {0.0};
    this->GetNthControlPointPositionWorld(1, position_1);
    double position_2[3] = {0.0};
    this->GetNthControlPointPositionWorld(2, position_2);

    vtkNew<vtkPoints> controlPoints;
    this->FillControlPointGridFromCorners(position_0, position_1, position_2, controlPoints);

    this->SetControlPointPositionsWorld(controlPoints);
  }
  else if (numberOfControlPoints > NUMBER_OF_PLANE_CONTROL_POINTS
    && numberOfControlPoints != this->GridResolution[0] * this->GridResolution[1])
  {
    // Re-generate new control point positions from current control points
    // after the grid size has been changed.

    //TODO:
    vtkNew<vtkPoints> pointsWorld;
    //pointsWorld->InsertNextPoint(center_World);
    //this->SetControlPointPositionsWorld(pointsWorld);
    this->MaximumNumberOfControlPoints = this->GridResolution[0] * this->GridResolution[1];
    this->RequiredNumberOfControlPoints = this->GridResolution[0] * this->GridResolution[1];
  }

  this->IsUpdatingControlPointsFromGridSurface = false;
}

//----------------------------------------------------------------------------
void vtkMRMLMarkupsGridSurfaceNode::FillControlPointGridFromCorners(double position_0[3], double position_1[3], double position_2[3], vtkPoints* controlPoints)
{
  if (!controlPoints)
  {
    return;
  }
  controlPoints->Reset();

  //
  // Calculate fourth point of plane rectangle from first three control points (numbers according to order of placement):
  //   Get vector between second and third point and add it to first point.
  //
  //             vector_a *
  //          GridResolution[0]
  //     0------------------------>1
  //     |                         |
  //     |vector_b *               |
  //     |GridResolution[1]        |
  //     |                         v
  //     3-------------------------2
  //
  double vector_1_2[3] = {0.0};
  vtkMath::Subtract(position_2, position_1, vector_1_2);
  double position_3[3] = {0.0};
  vtkMath::Add(position_0, vector_1_2, position_3);

  //
  // Calculate grid points based on specified grid size
  //

  // Vectors in both directions for displacement between adjacent control points
  double vector_a[3] = {0.0}; // = vector_0_1 / GridResolution[0]
  vtkMath::Subtract(position_1, position_0, vector_a);
  vtkMath::MultiplyScalar(vector_a, 1.0 / (double)(this->GridResolution[0] - 1));
  double vector_b[3] = { vector_1_2[0], vector_1_2[1], vector_1_2[2] }; // = vector_1_2 / GridResolution[1]
  vtkMath::MultiplyScalar(vector_b, 1.0 / (double)(this->GridResolution[1] - 1));

  // Fill in control point list
  for (int b = 0; b < this->GridResolution[1]; ++b)
  {
    double vector_b_scaled[3] = {vector_b[0], vector_b[1], vector_b[2]};
    vtkMath::MultiplyScalar(vector_b_scaled, b);

    for (int a = 0; a < this->GridResolution[0]; ++a)
    {
      double vector_a_scaled[3] = {vector_a[0], vector_a[1], vector_a[2]};
      vtkMath::MultiplyScalar(vector_a_scaled, a);

      double position_a_b[3] = {0.0};
      vtkMath::Add(position_0, vector_a_scaled, position_a_b);
      vtkMath::Add(position_a_b, vector_b_scaled, position_a_b);

      controlPoints->InsertNextPoint(position_a_b);
    }
  }
}

//----------------------------------------------------------------------------
void vtkMRMLMarkupsGridSurfaceNode::ResampleToNewGridResolution()
{
  if (this->PreviousGridResolution[0] == 0 || this->PreviousGridResolution[1] == 0)
  {
    vtkErrorMacro("ResampleToNewGridResolution: Cannot resample to new grid resolution because previous resolution is invalid");
    return;
  }

  vtkNew<vtkPoints> oldPoints;
  this->GetControlPointPositionsWorld(oldPoints);

  // Get corners and re-generate flat control points
  //TODO: This is a very basic way to change grid resolution and a proper resampling method will be needed!

  double position_0[3] = {0.0};
  this->GetNthControlPointPositionWorld(0, position_0);
  double position_1[3] = {0.0};
  this->GetNthControlPointPositionWorld(this->PreviousGridResolution[0]-1, position_1);
  double position_2[3] = {0.0};
  this->GetNthControlPointPositionWorld((this->PreviousGridResolution[0]*this->PreviousGridResolution[1])-1, position_2);

  vtkNew<vtkPoints> controlPoints;
  this->FillControlPointGridFromCorners(position_0, position_1, position_2, controlPoints);

  this->SetControlPointPositionsWorld(controlPoints);

  // Reset previous grid resolution to indicate that resampling is done
  this->PreviousGridResolution[0] = 0;
  this->PreviousGridResolution[1] = 0;
}
