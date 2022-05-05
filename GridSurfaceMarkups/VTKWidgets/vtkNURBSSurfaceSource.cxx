/*==============================================================================

  Distributed under the OSI-approved BSD 3-Clause License.

  Copyright (c) Children's Hospital of Philadelphia. All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.

  * Neither the name of Kitware, Inc. nor the names of Contributors
    may be used to endorse or promote products derived from this
    software without specific prior written permission.

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

  This file was originally developed by Csaba Pinter (Pixel Medical /
  Ebatinca).

==============================================================================*/

#include "vtkNURBSSurfaceSource.h"

// VTK includes
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkExecutive.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkMath.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>

//-------------------------------------------------------------------------------
vtkStandardNewMacro(vtkNURBSSurfaceSource);

//-------------------------------------------------------------------------------
vtkNURBSSurfaceSource::vtkNURBSSurfaceSource()
{
  this->SetNumberOfInputPorts(1);
}

//-------------------------------------------------------------------------------
vtkNURBSSurfaceSource::~vtkNURBSSurfaceSource()
{
}

//-------------------------------------------------------------------------------
void vtkNURBSSurfaceSource::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkPolyDataAlgorithm::PrintSelf(os, indent);

  os << "Input resolution: " << this->InputResolution[0] << ", " << this->InputResolution[1] << "\n";
  os << "Interpolation degrees: " << this->InterpolationDegrees[0] << ", " << this->InterpolationDegrees[1] << "\n";

  os << "Delta: " << this->Delta << "\n";
  os << "Use centripetal: " << this->UseCentripetal << "\n";
  os << "Wrap around: " << vtkMRMLMarkupsGridSurfaceNode::GetWrapAroundAsString(this->WrapAround) << "\n";
}

//----------------------------------------------------------------------------
int vtkNURBSSurfaceSource::FillInputPortInformation(int port, vtkInformation* info)
{
  if (port == 0)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    }
  else
    {
    vtkErrorMacro("Cannot set input info for port " << port);
    return 0;
    }

  return 1;
}

//-------------------------------------------------------------------------------
int vtkNURBSSurfaceSource::RequestData(
  vtkInformation* vtkNotUsed(request), vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkPointSet* inputPointSet = vtkPointSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  if (!inputPointSet)
  {
    return 1;
  }

  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkPolyData* outputPolyData = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  this->ComputeNurbsPolyData(inputPointSet->GetPoints(), outputPolyData);

  return 1;
}

//---------------------------------------------------------------------------
void vtkNURBSSurfaceSource::ComputeNurbsPolyData(vtkPoints* inputPoints, vtkPolyData* outputPolyData) // (points, size_u, size_v, degree_u, degree_v, **kwargs)
{
  if (outputPolyData == nullptr)
  {
    return;
  }

  if (!inputPoints)
  {
    vtkErrorMacro("ComputeNurbsPolyData: No valid input point list is set");
    return;
  }
  if (this->InterpolationDegrees[0] < 2 || this->InterpolationDegrees[1] < 2)
  {
    vtkErrorMacro("ComputeNurbsPolyData: No valid interpolation degree is set - they need to be at least 2");
    return;
  }
  if (this->InputResolution[0] < InterpolationDegrees[0] || this->InputResolution[1] < this->InterpolationDegrees[1])
  {
    vtkErrorMacro("ComputeNurbsPolyData: No valid input resolution is set - they need to be at least 2");
    return;
  }

  // Account for points included for wrapping around
  int interpolatingOverlap[2] = {0};
  this->GetInterpolatingOverlap(interpolatingOverlap);
  int interpolatingGridResolution[2] = {0};
  this->GetInterpolatingGridResolution(interpolatingGridResolution);
  
  // Get parameter arrays in the two directions
  vtkNew<vtkDoubleArray> ukParams;
  vtkNew<vtkDoubleArray> vlParams;
  this->ComputeParamsSurface(inputPoints, ukParams, vlParams);

  // Compute knot vectors
  vtkNew<vtkDoubleArray> uKnots;
  this->ComputeKnotVector(this->InterpolationDegrees[0], interpolatingGridResolution[0], ukParams, uKnots);
  vtkNew<vtkDoubleArray> vKnots;
  this->ComputeKnotVector(this->InterpolationDegrees[1], interpolatingGridResolution[1], vlParams, vKnots);

  //
  // Do global interpolation along the u direction
  //
  vtkNew<vtkPoints> controlPointsR;
  double** matrixA = this->AllocateMatrix(interpolatingGridResolution[0]);
  for (int v = -interpolatingOverlap[1]; v < this->InputResolution[1] + interpolatingOverlap[1]; ++v)
  {
    // Collect column point indices
    vtkNew<vtkPoints> points;
    for (int u = -interpolatingOverlap[0]; u < this->InputResolution[0] + interpolatingOverlap[0]; ++u)
    {
      points->InsertNextPoint(inputPoints->GetPoint(this->GetPointIndexUV(u,v)));
    }

    this->BuildCoeffMatrix(this->InterpolationDegrees[0], uKnots, ukParams, points, matrixA);

    // Insert control points for each point column
    this->LuSolve(matrixA, interpolatingGridResolution[0], points, controlPointsR);
  }
  this->DestructMatrix(matrixA, interpolatingGridResolution[0]);

  //
  // Do global interpolation along the v direction
  //
  vtkNew<vtkPoints> controlPoints;
  matrixA = this->AllocateMatrix(interpolatingGridResolution[1]);
  for (int u = -interpolatingOverlap[0]; u < this->InputResolution[0] + interpolatingOverlap[0]; ++u)
  {
    // Collect row point indices
    vtkNew<vtkPoints> points;
    for (int v = -interpolatingOverlap[1]; v < this->InputResolution[1] + interpolatingOverlap[1]; ++v)
    {
      int controlPointRIndex = (v + interpolatingOverlap[1]) * interpolatingGridResolution[0] + (u + interpolatingOverlap[0]);
      points->InsertNextPoint(controlPointsR->GetPoint(controlPointRIndex));
    }

    this->BuildCoeffMatrix(this->InterpolationDegrees[1], vKnots, vlParams, points, matrixA);

    // Insert control points for each point row
    this->LuSolve(matrixA, interpolatingGridResolution[1], points, controlPoints);
  }
  this->DestructMatrix(matrixA, interpolatingGridResolution[1]);

  //
  // Construct surface
  //
  vtkNew<vtkPoints> evalPoints;
  this->EvaluateSurface(uKnots, vKnots, controlPoints, evalPoints);
  outputPolyData->SetPoints(evalPoints);

  this->TriangulateSurface(outputPolyData);

  // Add point indices as scalar array for debugging purposes
  vtkNew<vtkIntArray> indexArray;
  indexArray->SetName("PointIndices");
  for (int i=0; i<evalPoints->GetNumberOfPoints(); ++i)
  {
    indexArray->InsertNextValue(i);
  }
  outputPolyData->GetPointData()->AddArray(indexArray);
}

//---------------------------------------------------------------------------
void vtkNURBSSurfaceSource::EvaluateSurface(vtkDoubleArray* uKnots, vtkDoubleArray* vKnots, vtkPoints* controlPoints, vtkPoints* outEvalPoints)
{
  if (this->Delta < 0.001)
  {
    vtkErrorMacro("EvaluateSurface: Delta value too small");
    return;
  }

  int dimension = 3; // We only work in three dimensions

  int interpolatingOverlap[2] = {0};
  this->GetInterpolatingOverlap(interpolatingOverlap);
  int interpolatingGridResolution[2] = {0};
  this->GetInterpolatingGridResolution(interpolatingGridResolution);

  // Calculate parameter space to evaluate (expand / shrink, remove wrap around overlap)
  double minLinSpaceU = 0, maxLinSpaceU = 0, minLinSpaceV = 0, maxLinSpaceV = 0;
  this->CalculateEvaluatedParameterSpace(minLinSpaceU, maxLinSpaceU, minLinSpaceV, maxLinSpaceV);
  int sampleSizeU = 0, sampleSizeV = 0;
  this->CalculateSampleSize(sampleSizeU, sampleSizeV);

  vtkNew<vtkDoubleArray> knotsU;
  this->LinSpace(minLinSpaceU, maxLinSpaceU, sampleSizeU, knotsU);
  vtkNew<vtkDoubleArray> knotsV;
  this->LinSpace(minLinSpaceV, maxLinSpaceV, sampleSizeV, knotsV);

  vtkNew<vtkIntArray> uSpans;
  vtkNew<vtkDoubleArray> uBasis;
  this->FindSpans(this->InterpolationDegrees[0], uKnots, interpolatingGridResolution[0], knotsU, uSpans);
  this->BasisFunctions(this->InterpolationDegrees[0], uKnots, uSpans, knotsU, uBasis);

  vtkNew<vtkIntArray> vSpans;
  vtkNew<vtkDoubleArray> vBasis;
  this->FindSpans(this->InterpolationDegrees[1], vKnots, interpolatingGridResolution[1], knotsV, vSpans);
  this->BasisFunctions(this->InterpolationDegrees[1], vKnots, vSpans, knotsV, vBasis);

  // Store first of the two overlapping sections separately in order to unify them afterwards
  vtkNew<vtkPoints> firstOverlappingPoints;

  outEvalPoints->Initialize();
  for (int i=0; i<uSpans->GetNumberOfValues(); ++i)
  {
    int idxU = uSpans->GetValue(i) - this->InterpolationDegrees[0];
    for (int j=0; j<vSpans->GetNumberOfValues(); ++j)
    {
      int idxV = vSpans->GetValue(j) - this->InterpolationDegrees[1];
      double spt[3] = {0.0};
      for (int k=0; k<this->InterpolationDegrees[0]+1; ++k)
      {
        double temp[3] = {0.0};
        for (int l=0; l<this->InterpolationDegrees[1]+1; ++l)
        {
          double* controlPoint = controlPoints->GetPoint(idxV + l + (interpolatingGridResolution[1] * (idxU + k)));
          for (int d=0; d<dimension; ++d)
          {
            temp[d] = temp[d] + vBasis->GetValue(j * (this->InterpolationDegrees[1]+1) + l) * controlPoint[d];
          }
        }
        for (int d=0; d<dimension; ++d)
        {
          spt[d] = spt[d] + uBasis->GetValue(i * (this->InterpolationDegrees[0]+1) + k) * temp[d];
        }
      }
      outEvalPoints->InsertNextPoint(spt);
    } // v direction
  } // u direction
}

//---------------------------------------------------------------------------
void vtkNURBSSurfaceSource::TriangulateSurface(vtkPolyData* outputPolyData)
{
  if (!outputPolyData)
  {
    vtkErrorMacro("TriangulateSurface: Invalid poly data given");
    return;
  }

  int sampleSizeU = 0, sampleSizeV = 0;
  this->CalculateSampleSize(sampleSizeU, sampleSizeV);

  vtkNew<vtkCellArray> cells;
  for (unsigned int i=0; i<sampleSizeU-1; i++)
  {
    for (unsigned int j=0; j<sampleSizeV-1; j++)
    {
      unsigned int base = i*sampleSizeV + j;
      unsigned int a = base;
      unsigned int b = base + 1;
      unsigned int c = base + sampleSizeV + 1;
      unsigned int d = base + sampleSizeV;
      vtkIdType triangle[3] = {0};

      triangle[0] = c;
      triangle[1] = b;
      triangle[2] = a;
      cells->InsertNextCell(3, triangle);

      triangle[0] = d;
      triangle[1] = c;
      triangle[2] = a;
      cells->InsertNextCell(3, triangle);
    }
  }

  // Insert strip of triangles between the meeting wrapped around edges
  if (this->WrapAround == vtkMRMLMarkupsGridSurfaceNode::AlongU)
  {
    for (unsigned int v=0; v<sampleSizeV-1; v++)
    {
      unsigned int base = v;
      unsigned int a = base;
      unsigned int b = base + 1;
      unsigned int c = base + (sampleSizeU-1) * sampleSizeV + 1;
      unsigned int d = base + (sampleSizeU-1) * sampleSizeV;
      vtkIdType triangle[3] = {0};

      triangle[0] = c;
      triangle[1] = b;
      triangle[2] = d;
      cells->InsertNextCell(3, triangle);

      triangle[0] = d;
      triangle[1] = b;
      triangle[2] = a;
      cells->InsertNextCell(3, triangle);
    }
  }
  else if (this->WrapAround == vtkMRMLMarkupsGridSurfaceNode::AlongV)
  {
    for (unsigned int u=0; u<sampleSizeU-1; u++)
    {
      unsigned int base = u*sampleSizeV;
      unsigned int a = base;
      unsigned int b = base + sampleSizeV;
      unsigned int c = base + sampleSizeV * 2 - 1;
      unsigned int d = base + sampleSizeV - 1;
      vtkIdType triangle[3] = {0};

      triangle[0] = c;
      triangle[1] = b;
      triangle[2] = d;
      cells->InsertNextCell(3, triangle);

      triangle[0] = d;
      triangle[1] = b;
      triangle[2] = a;
      cells->InsertNextCell(3, triangle);
    }
  }
  outputPolyData->SetPolys(cells);
}

//---------------------------------------------------------------------------
void vtkNURBSSurfaceSource::ComputeParamsSurface(vtkPoints* inputPoints, vtkDoubleArray* ukParams, vtkDoubleArray* vlParams)
{
  // Account for points included for wrapping around
  int interpolatingOverlap[2] = {0};
  this->GetInterpolatingOverlap(interpolatingOverlap);
  int interpolatingGridResolution[2] = {0};
  this->GetInterpolatingGridResolution(interpolatingGridResolution);

  // Compute parameters for each curve on the v direction. Concatenate into one parameter array
  vtkNew<vtkDoubleArray> ukParamsTemp;
  for (int v = -interpolatingOverlap[1]; v < this->InputResolution[1] + interpolatingOverlap[1]; ++v)
  {
    // Collect column point indices
    vtkNew<vtkIdList> columnIds;
    for (int u = -interpolatingOverlap[0]; u < this->InputResolution[0] + interpolatingOverlap[0]; ++u)
    {
      columnIds->InsertNextId(this->GetPointIndexUV(u,v));
    }

    // Compute parameters for column
    this->ComputeParamsCurve(inputPoints, columnIds, ukParamsTemp);
  }

  // Do averaging on the u direction
  ukParams->Initialize();
  for (int u = 0; u < interpolatingGridResolution[0]; ++u)
  {
    double sumKnots_v = 0.0;
    for (int v = 0; v < interpolatingGridResolution[1]; ++v)
    {
      sumKnots_v += ukParamsTemp->GetValue(u + (interpolatingGridResolution[0] * v));
    }
    ukParams->InsertNextValue(sumKnots_v / interpolatingGridResolution[1]);
  }

  // Compute parameters for each curve on the u direction. Concatenate into one parameter array
  vtkNew<vtkDoubleArray> vlParamsTemp;
  for (int u = -interpolatingOverlap[0]; u < this->InputResolution[0] + interpolatingOverlap[0]; ++u)
  {
    // Collect row point indices
    vtkNew<vtkIdList> rowIds;
    for (int v = -interpolatingOverlap[1]; v < this->InputResolution[1] + interpolatingOverlap[1]; ++v)
    {
      rowIds->InsertNextId(this->GetPointIndexUV(u,v));
    }

    // Compute parameters for row
    this->ComputeParamsCurve(inputPoints, rowIds, vlParamsTemp);
  }

  // Do averaging on the v direction
  vlParams->Initialize();
  for (int v = 0; v < interpolatingGridResolution[1]; ++v)
  {
    double sumKnots_u = 0.0;
    for (int u = 0; u < interpolatingGridResolution[0]; ++u)
    {
      sumKnots_u += vlParamsTemp->GetValue(v + (interpolatingGridResolution[1] * u));
    }
    vlParams->InsertNextValue(sumKnots_u / interpolatingGridResolution[0]);
  }
}

//---------------------------------------------------------------------------
void vtkNURBSSurfaceSource::ComputeParamsCurve(
  vtkPoints* inputPoints, vtkIdList* pointIndexList, vtkDoubleArray* outParametersArray)
{
  if (!pointIndexList || !outParametersArray)
  {
    vtkErrorMacro("ComputeParamsCurve: Invalid arguments");
    return;
  }

  // Calculate chord lengths between the curve points
  vtkNew<vtkDoubleArray> chordLengthsArray;
  double prevPoint[3] = {0.0};
  double currentPoint[3] = {0.0};
  double distance = 0.0;
  double sumAllChordLengths = 0.0;
  for (int i=1; i<pointIndexList->GetNumberOfIds(); ++i)
  {
    inputPoints->GetPoint(pointIndexList->GetId(i-1), prevPoint);
    inputPoints->GetPoint(pointIndexList->GetId(i), currentPoint);
    double distance = sqrt(vtkMath::Distance2BetweenPoints(currentPoint, prevPoint));
    double chordLength = this->UseCentripetal ? sqrt(distance) : distance;
    chordLengthsArray->InsertNextValue(chordLength);
    sumAllChordLengths += chordLength;
  }

  // Insert first parameter (always with value 0) to the array
  outParametersArray->InsertNextValue(0.0);

  // Divide individual chord lengths by the total chord length and insert parameter into output array
  double currentSumChordLengths = 0.0;
  for (int i=0; i<chordLengthsArray->GetNumberOfValues()-1; ++i)
  {
    currentSumChordLengths += chordLengthsArray->GetValue(i);
    outParametersArray->InsertNextValue(currentSumChordLengths / sumAllChordLengths);
  }

  // Insert last parameter (always with value 1) to the array
  outParametersArray->InsertNextValue(1.0);
}

//---------------------------------------------------------------------------
void vtkNURBSSurfaceSource::ComputeKnotVector(
  int degree, int numOfPoints, vtkDoubleArray* inParams, vtkDoubleArray* outKnotVector)
{
  if (degree == 0 || numOfPoints == 0)
  {
    vtkErrorMacro("ComputeKnotVector: Input values should be different than zero");
    return;
  }
  if (!inParams)
  {
    vtkErrorMacro("ComputeKnotVector: Invalid parameters array");
    return;
  }
  if (!outKnotVector)
  {
    vtkErrorMacro("ComputeKnotVector: Invalid output array");
    return;
  }

  outKnotVector->Initialize();

  // Start knot vector
  for (int i=0; i<degree+1; ++i)
  {
    outKnotVector->InsertNextValue(0.0);
  }

  // Use averaging method (Eqn 9.8) to compute internal knots in the knot vector
  for (int i=0; i<numOfPoints-degree-1; ++i)
  {
    double currentSumParams = 0.0;
    for (int j=i+1; j<i+degree+1; ++j)
    {
      currentSumParams += inParams->GetValue(j);
    }

    outKnotVector->InsertNextValue((1.0/degree) * currentSumParams);
  }

  // End knot vector
  for (int i=0; i<degree+1; ++i)
  {
    outKnotVector->InsertNextValue(1.0);
  }
}

//---------------------------------------------------------------------------
void vtkNURBSSurfaceSource::BuildCoeffMatrix(
  int degree, vtkDoubleArray* knotVector, vtkDoubleArray* params, vtkPoints* points, double** outCoeffMatrix)
{
  if (degree == 0)
  {
    vtkErrorMacro("BuildCoeffMatrix: Invalid input degree");
    return;
  }
  if (!knotVector || !params || !points)
  {
    vtkErrorMacro("BuildCoeffMatrix: Invalid input arrays");
    return;
  }
  if (!outCoeffMatrix)
  {
    vtkErrorMacro("BuildCoeffMatrix: Invalid output array");
    return;
  }

  int numPoints = points->GetNumberOfPoints();

  // Initialize coefficient matrix
  for (int i=0; i<numPoints; ++i)
  {
    for (int j=0; j<numPoints; ++j)
    {
      outCoeffMatrix[i][j] = 0.0;
    }
  }

  for (int i=0; i<numPoints; ++i)
  {
    double knot = params->GetValue(i);
    int span = this->FindSpanLinear(degree, knotVector, numPoints, knot);
    vtkNew<vtkDoubleArray> basisFunction;
    this->BasisFunction(degree, knotVector, span, knot, basisFunction);
    for (int n=0; n<basisFunction->GetNumberOfValues(); ++n)
    {
      outCoeffMatrix[i][span-degree+n] = basisFunction->GetValue(n);
    }
  }
}

//
// helpers
//

//---------------------------------------------------------------------------
void vtkNURBSSurfaceSource::BasisFunction(int degree, vtkDoubleArray* knotVector, int span, double knot, vtkDoubleArray* outBasisFunctions)
{
  if (!outBasisFunctions)
  {
    vtkErrorMacro("BasisFunction: Invalid output array");
    return;
  }

  double* left = new double[degree + 1];
  double* right = new double[degree + 1];
  vtkDoubleArray* N = outBasisFunctions;
  N->Initialize();
  N->SetNumberOfValues(degree + 1);
  for (int i=0; i<degree+1; ++i)
  {
    left[i] = 0.0;
    right[i] = 0.0;
    N->SetValue(i, 1.0); // N[0] = 1.0 by definition
  }

  for (int j=1; j<degree+1; ++j)
  {
    left[j] = knot - knotVector->GetValue(span + 1 - j);
    right[j] = knotVector->GetValue(span + j) - knot;
    double saved = 0.0;
    for (int r=0; r<j; ++r)
    {
      double temp = N->GetValue(r) / (right[r + 1] + left[j - r]);
      N->SetValue(r, saved + right[r + 1] * temp);
      saved = left[j - r] * temp;
    }

    N->SetValue(j, saved);
  }

  delete[] left;
  delete[] right;
}

//---------------------------------------------------------------------------
void vtkNURBSSurfaceSource::BasisFunctions(
  int degree, vtkDoubleArray* knotVector, vtkIntArray* spans, vtkDoubleArray* knots, vtkDoubleArray* outBasisFunctions)
{
  if (!knotVector || !spans || !knots)
  {
    vtkErrorMacro("BasisFunctions: Invalid input arrays");
    return;
  }
  if (spans->GetNumberOfValues() != knots->GetNumberOfValues())
  {
    vtkErrorMacro("BasisFunctions: Spans and knots count mismatch");
    return;
  }
  if (!outBasisFunctions)
  {
    vtkErrorMacro("BasisFunctions: Invalid output array");
    return;
  }

  outBasisFunctions->Initialize();
  for (int i=0; i<knots->GetNumberOfValues(); ++i)
  {
    vtkNew<vtkDoubleArray> currentBasisFunctions;
    this->BasisFunction(degree, knotVector, spans->GetValue(i), knots->GetValue(i), currentBasisFunctions);
    for (int j=0; j<currentBasisFunctions->GetNumberOfValues(); ++j)
    {
      outBasisFunctions->InsertNextValue(currentBasisFunctions->GetValue(j));
    }
  }
}

//---------------------------------------------------------------------------
int vtkNURBSSurfaceSource::FindSpanLinear(int degree, vtkDoubleArray* knotVector, int numControlPoints, double knot)
{
  // Knot span index starts from zero
  int span = degree + 1;

  while (span < numControlPoints && knotVector->GetValue(span) <= knot)
  {
    span++;
  }

  return span - 1;
}

//---------------------------------------------------------------------------
void vtkNURBSSurfaceSource::FindSpans(int degree, vtkDoubleArray* knotVector, int numControlPoints, vtkDoubleArray* knots, vtkIntArray* outSpans)
{
  if (!knotVector || !knots)
  {
    vtkErrorMacro("FindSpans: Invalid input arrays");
    return;
  }
  if (!outSpans)
  {
    vtkErrorMacro("FindSpans: Invalid output array");
    return;
  }

  outSpans->Initialize();
  for (int i=0; i<knots->GetNumberOfValues(); ++i)
  {
    outSpans->InsertNextValue(this->FindSpanLinear(degree, knotVector, numControlPoints, knots->GetValue(i)));
  }
}

//
// linalg
//

//---------------------------------------------------------------------------
void vtkNURBSSurfaceSource::LuSolve(double** coeffMatrix, int matrixSize, vtkPoints* points, vtkPoints* outControlPointsR)
{
  if (!coeffMatrix || !points)
  {
    vtkErrorMacro("LuSolve: Invalid input coefficients matrix or points");
    return;
  }
  if (!outControlPointsR)
  {
    vtkErrorMacro("LuSolve: Invalid output point array");
    return;
  }

  int dim = 3; // We only work in three dimensions
  int numX = points->GetNumberOfPoints();

  double** matrixL = this->AllocateMatrix(matrixSize);
  double** matrixU = this->AllocateMatrix(matrixSize);
  double** x = this->AllocateMatrix(numX, dim);

  this->LuDecomposition(coeffMatrix, matrixL, matrixU, matrixSize);

  for (int i=0; i<dim; ++i)
  {
    // Get i'th coordinate of each point into a column vector
    double* b = new double[numX];
    for (int j=0; j<numX; ++j)
    {
      b[j] = points->GetPoint(j)[i];
    }

    double* y = new double[numX];
    this->ForwardSubstitution(matrixL, b, numX, y);

    double* xt = new double[numX];
    this->BackwardSubstitution(matrixU, y, numX, xt);

    for (int j=0; j<numX; ++j)
    {
      x[j][i] = xt[j];
    }

    delete[] b;
    delete[] y;
    delete[] xt;
  }

  // Insert control points from solution matrix
  for (int j=0; j<numX; ++j)
  {
    outControlPointsR->InsertNextPoint(x[j][0], x[j][1], x[j][2]);
  }

  this->DestructMatrix(matrixL, matrixSize);
  this->DestructMatrix(matrixU, matrixSize);
  this->DestructMatrix(x, numX, dim);
}

//---------------------------------------------------------------------------
void vtkNURBSSurfaceSource::LuDecomposition(double** matrixA, double** matrixL, double** matrixU, int size)
{
  if (!matrixA || !matrixL || !matrixU)
  {
    vtkErrorMacro("LuDecomposition: Invalid input matrices");
    return;
  }

  int n = size;
  double** a = matrixA;
  double** l = matrixL;
  double** u = matrixU;

  // From https://titanwolf.org/Network/Articles/Article?AID=af9bbd90-f722-4bec-803c-523bd0ca1d9e
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      if (j < i)
        l[j][i] = 0;
      else
      {
        l[j][i] = a[j][i];
        for (int k = 0; k < i; k++)
        {
          l[j][i] = l[j][i] - l[j][k] * u[k][i];
        }
      }
    }
    for (int j = 0; j < n; j++)
    {
      if (j < i)
        u[i][j] = 0;
      else if (j == i)
        u[i][j] = 1;
      else
      {
        u[i][j] = a[i][j]/l[i][i];
        for (int k = 0; k < i; k++)
        {
          u[i][j] = u[i][j] - ((l[i][k] * u[k][j])/l[i][i]);
        }
      }
    }
  }
}

//---------------------------------------------------------------------------
void vtkNURBSSurfaceSource::ForwardSubstitution(double** matrixL, double* b, int size, double* outY)
{
  if (!matrixL || !b || !outY)
  {
    vtkErrorMacro("ForwardSubstitution: Invalid input/output arguments");
    return;
  }

  outY[0] = b[0] / matrixL[0][0];
  for (int i=1; i<size; ++i)
  {
    outY[i] = 0.0;
  }

  for (int i=1; i<size; ++i)
  {
    double sum = 0.0;
    for (int j=0; j<i; ++j)
    {
      sum += matrixL[i][j] * outY[j];
    }
    outY[i] = b[i] - sum;
    outY[i] /= matrixL[i][i];
  }
}

//---------------------------------------------------------------------------
void vtkNURBSSurfaceSource::BackwardSubstitution(double** matrixU, double* y, int size, double* outX)
{
  if (!matrixU || !y || !outX)
  {
    vtkErrorMacro("BackwardSubstitution: Invalid input/output arguments");
    return;
  }

  outX[size-1] = y[size-1] / matrixU[size-1][size-1];
  for (int i=0; i<size-1; ++i)
  {
    outX[i] = 0.0;
  }

  for (int i=size-2; i>=0; --i)
  {
    double sum = 0.0;
    for (int j=i; j<size; ++j)
    {
      sum += matrixU[i][j] * outX[j];
    }
    outX[i] = y[i] - sum;
    outX[i] /= matrixU[i][i];
  }
}

//---------------------------------------------------------------------------
void vtkNURBSSurfaceSource::LinSpace(double start, double stop, int numOfSamples, vtkDoubleArray* outSpace)
{
  if (!outSpace)
  {
    vtkErrorMacro("LinSpace: Invalid output array");
    return;
  }

  outSpace->Initialize();
  if (fabs(start - stop) <= 1e-8)
  {
    outSpace->InsertNextValue(start);
    return;
  }

  if (numOfSamples > 1)
  {
    int div = numOfSamples - 1;
    double delta = stop - start;

    for (int i=0; i<numOfSamples; ++i)
    {
      outSpace->InsertNextValue(start + (double)i * delta / div);
    }
  }
}

//
// Utility functions
//

//---------------------------------------------------------------------------
unsigned int vtkNURBSSurfaceSource::GetPointIndexUV(int u, int v)
{
  int interpolatingOverlap[2] = {0};
  this->GetInterpolatingOverlap(interpolatingOverlap);

  if ( u < -interpolatingOverlap[0] || v < -interpolatingOverlap[1] ||
       u >= this->InputResolution[0] + interpolatingOverlap[0] || v >= this->InputResolution[1] + interpolatingOverlap[1] )
  {
    vtkErrorMacro("GetPointUV: Index pair (" << u << ", " << v << ") is out of range");
    return -1;
  }

  if (u < 0)
  {
    u += this->InputResolution[0];
  }
  else if (u >= this->InputResolution[0])
  {
    u -= this->InputResolution[0];
  }

  if (v < 0)
  {
    v += this->InputResolution[1];
  }
  else if (v >= this->InputResolution[1])
  {
    v -= this->InputResolution[1];
  }

  return u * this->InputResolution[1] + v;
}

//---------------------------------------------------------------------------
void vtkNURBSSurfaceSource::GetInterpolatingGridResolution(int (&resolutionUV)[2])
{
  int interpolatingOverlap[2] = {0};
  this->GetInterpolatingOverlap(interpolatingOverlap);

  resolutionUV[0] = this->InputResolution[0] + 2 * interpolatingOverlap[0];
  resolutionUV[1] = this->InputResolution[1] + 2 * interpolatingOverlap[1];
}

//---------------------------------------------------------------------------
void vtkNURBSSurfaceSource::GetInterpolatingOverlap(int (&overlapUV)[2])
{
  if (this->WrapAround == vtkMRMLMarkupsGridSurfaceNode::AlongU)
  {
    overlapUV[0] = this->InterpolationDegrees[0] + 1;
    overlapUV[1] = 0;
  }
  else if (this->WrapAround == vtkMRMLMarkupsGridSurfaceNode::AlongV)
  {
    overlapUV[0] = 0;
    overlapUV[1] = this->InterpolationDegrees[1] + 1;
  }
  else // Disabled
  {
    overlapUV[0] = 0;
    overlapUV[1] = 0;
    if (this->WrapAround != vtkMRMLMarkupsGridSurfaceNode::NoWrap) // Valid value for disabled
    {
      vtkErrorMacro("GetInterpolatingOverlap: Invalid WrapAround value " << this->WrapAround);
    }
  }
}

//---------------------------------------------------------------------------
double** vtkNURBSSurfaceSource::AllocateMatrix(int m, int n/*=0*/)
{
  if (n == 0)
  {
    n = m;
  }
  double** matrix = new double*[m];
  for (int i=0; i<m; ++i)
  {
    matrix[i] = new double[n];
  }

  return matrix;
}

//---------------------------------------------------------------------------
void vtkNURBSSurfaceSource::DestructMatrix(double** matrix, int m, int n/*=0*/)
{
  if (n == 0)
  {
    n = m;
  }
  for (int i=0; i<m; ++i)
  {
    delete[] matrix[i];
  }
  delete[] matrix;
}

//---------------------------------------------------------------------------
void vtkNURBSSurfaceSource::CalculateEvaluatedParameterSpace(
  double& minLinSpaceU, double& maxLinSpaceU, double& minLinSpaceV, double& maxLinSpaceV)
{
  // Calculate parameter space to evaluate (expand / shrink, remove wrap around overlap)
  if (this->WrapAround == vtkMRMLMarkupsGridSurfaceNode::NoWrap)
  {
    minLinSpaceU = -this->ExpansionFactor;
    maxLinSpaceU = 1.0 + this->ExpansionFactor;
    minLinSpaceV = -this->ExpansionFactor;
    maxLinSpaceV = 1.0 + this->ExpansionFactor;
  }
  else
  {
    int interpolatingOverlap[2] = {0};
    this->GetInterpolatingOverlap(interpolatingOverlap);
    int interpolatingGridResolution[2] = {0};
    this->GetInterpolatingGridResolution(interpolatingGridResolution);

    if (this->WrapAround == vtkMRMLMarkupsGridSurfaceNode::AlongU)
    {
      minLinSpaceU = 0.0; // Do not expand along the wrapped direction
      maxLinSpaceU = 1.0;
      minLinSpaceV = -this->ExpansionFactor;
      maxLinSpaceV = 1.0 + this->ExpansionFactor;

      // Leave the stitching face from the min side (arbitrary side selection).
      // Also leave two triangles worth of space for overreaching (which may happen when the spacing between the control points is very irregular)
      //TODO: Try to fix the root cause instead, which is currently unknown (why having irregularly spaced points change the linear space that much)
      int samplesPerGridCell = vtkMath::Floor((1.0 / this->Delta) + 0.5);
      double linSpaceUPerTriangle = 2.0 / ((interpolatingGridResolution[0] - 1) * samplesPerGridCell);
      minLinSpaceU += (double)(interpolatingOverlap[0] - 1) / (interpolatingGridResolution[0] - 1) + linSpaceUPerTriangle * 2;
      maxLinSpaceU -= (double)interpolatingOverlap[0] / (interpolatingGridResolution[0] - 1);
    }
    else if (this->WrapAround == vtkMRMLMarkupsGridSurfaceNode::AlongV)
    {
      minLinSpaceU = -this->ExpansionFactor;
      maxLinSpaceU = 1.0 + this->ExpansionFactor;
      minLinSpaceV = 0.0; // Do not expand along the wrapped direction
      maxLinSpaceV = 1.0;

      // Leave the stitching face from the min side (arbitrary side selection).
      // Also leave two triangles worth of space for overreaching (which may happen when the spacing between the control points is very irregular)
      //TODO: Try to fix the root cause instead, which is currently unknown (why having irregularly spaced points change the linear space that much)
      int samplesPerGridCell = vtkMath::Floor((1.0 / this->Delta) + 0.5);
      double linSpaceVPerTriangle = 2.0 / ((interpolatingGridResolution[1] - 1) * samplesPerGridCell);
      minLinSpaceV += (double)(interpolatingOverlap[1] - 1) / (interpolatingGridResolution[1] - 1) + linSpaceVPerTriangle * 2;
      maxLinSpaceV -= (double)interpolatingOverlap[1] / (interpolatingGridResolution[1] - 1);
    }
  }
}

//---------------------------------------------------------------------------
void vtkNURBSSurfaceSource::CalculateSampleSize(int& sampleSizeU, int& sampleSizeV)
{
  int interpolatingGridResolution[2] = {0};
  this->GetInterpolatingGridResolution(interpolatingGridResolution);

  double minLinSpaceU = 0, maxLinSpaceU = 0, minLinSpaceV = 0, maxLinSpaceV = 0;
  this->CalculateEvaluatedParameterSpace(minLinSpaceU, maxLinSpaceU, minLinSpaceV, maxLinSpaceV);
  double linSpaceSizeU = maxLinSpaceU - minLinSpaceU;
  double linSpaceSizeV = maxLinSpaceV - minLinSpaceV;

  int samplesPerGridCell = vtkMath::Floor((1.0 / this->Delta) + 0.5);
  sampleSizeU = samplesPerGridCell * (int)((interpolatingGridResolution[0] - 1) * linSpaceSizeU + 0.5) + 1;
  sampleSizeV = samplesPerGridCell * (int)((interpolatingGridResolution[1] - 1) * linSpaceSizeV + 0.5) + 1;
}
