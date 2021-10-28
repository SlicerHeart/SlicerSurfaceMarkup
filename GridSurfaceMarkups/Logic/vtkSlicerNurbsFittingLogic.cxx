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

#include "vtkSlicerNurbsFittingLogic.h"

// GridSurface Markups MRML includes
#include "vtkMRMLMarkupsGridSurfaceNode.h"

// MRML includes
//#include <vtkMRMLModelNode.h>
#include <vtkMRMLScene.h>

// Markups MRML includes
//#include <vtkMRMLMarkupsDisplayNode.h>

// VTK includes
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkExecutive.h>
#include <vtkIdList.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkMath.h>
#include <vtkObjectFactory.h>

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerNurbsFittingLogic);

//---------------------------------------------------------------------------
vtkSlicerNurbsFittingLogic::vtkSlicerNurbsFittingLogic()
{

}

//---------------------------------------------------------------------------
vtkSlicerNurbsFittingLogic::~vtkSlicerNurbsFittingLogic() = default;

//---------------------------------------------------------------------------
void vtkSlicerNurbsFittingLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
void vtkSlicerNurbsFittingLogic::SetInputPoints(vtkPoints* points)
{
  if (this->InputPoints == points)
  {
    return;
  }

  this->InputPoints = points;
  this->Modified();
}

//-------------------------------------------------------------------------------
int vtkSlicerNurbsFittingLogic::RequestData(
  vtkInformation* vtkNotUsed(request), vtkInformationVector** vtkNotUsed(inputVector), vtkInformationVector* outputVector)
{
  vtkInformation* bezierSurfaceOutputInfo = outputVector->GetInformationObject(0);
  if (bezierSurfaceOutputInfo)
  {
    vtkPolyData* nurbsSurfaceOutput =
      vtkPolyData::SafeDownCast(bezierSurfaceOutputInfo->Get(vtkDataObject::DATA_OBJECT()));
    this->UpdateNurbsPolyData(nurbsSurfaceOutput);
    //nurbsSurfaceOutput->SetPolys(this->Topology);
  }

  return 1;
}

//---------------------------------------------------------------------------
void vtkSlicerNurbsFittingLogic::UpdateNurbsPolyData(vtkPolyData* polyData) // (points, size_u, size_v, degree_u, degree_v, **kwargs)
{
  // # Keyword arguments
  // use_centripetal = kwargs.get('centripetal', False)
  //
  // # Get uk and vl
  // uk, vl = compute_params_surface(points, size_u, size_v, use_centripetal)
  //
  // # Compute knot vectors
  // kv_u = compute_knot_vector(degree_u, size_u, uk)
  // kv_v = compute_knot_vector(degree_v, size_v, vl)
  //
  // # Do global interpolation on the u-direction
  // ctrlpts_r = []
  // for v in range(size_v):
  //     pts = [points[v + (size_v * u)] for u in range(size_u)] #CP_Note: Get vth column of points
  //     matrix_a = _build_coeff_matrix(degree_u, kv_u, uk, pts)
  //     ctrlpts_r += linalg.lu_solve(matrix_a, pts)
  //
  // # Do global interpolation on the v-direction
  // ctrlpts = []
  // for u in range(size_u):
  //     pts = [ctrlpts_r[u + (size_u * v)] for v in range(size_v)]
  //     matrix_a = _build_coeff_matrix(degree_v, kv_v, vl, pts)
  // #CP_Note: linalg.lu_solve
  //     ctrlpts += linalg.lu_solve(matrix_a, pts)
  //
  // # Generate B-spline surface
  // surf = BSpline.Surface()
  // surf.degree_u = degree_u
  // surf.degree_v = degree_v
  // surf.ctrlpts_size_u = size_u
  // surf.ctrlpts_size_v = size_v
  // surf.ctrlpts = ctrlpts
  // surf.knotvector_u = kv_u
  // surf.knotvector_v = kv_v
  //
  // return surf

  if (polyData == nullptr)
  {
    return;
  }

  if (!this->InputPoints)
  {
    vtkErrorMacro("UpdateNurbsPolyData: No valid input point list is set");
    return;
  }
  if (this->InputResolution[0] < 2 || this->InputResolution[1] < 2)
  {
    vtkErrorMacro("UpdateNurbsPolyData: No valid input resolution is set");
    return;
  }

  // Get parameter arrays in the two directions
  vtkNew<vtkDoubleArray> ukParams;
  vtkNew<vtkDoubleArray> vlParams;
  this->ComputeParamsSurface(ukParams, vlParams);

  // Compute knot vectors
  vtkNew<vtkDoubleArray> uKnots;
  this->ComputeKnotVector(this->InterpolationDegrees[0], this->InputResolution[0], ukParams, uKnots);
  vtkNew<vtkDoubleArray> vKnots;
  this->ComputeKnotVector(this->InterpolationDegrees[1], this->InputResolution[1], vlParams, vKnots);

  //
  // Do global interpolation along the u direction
  //
  vtkNew<vtkPoints> controlPointsR;
  double** matrixA = this->AllocateMatrix(this->InputResolution[0]);
  for (int v=0; v<this->InputResolution[1]; ++v)
  {
    // Collect column point indices
    vtkNew<vtkPoints> points;
    for (int u=0; u<this->InputResolution[0]; ++u)
    {
      points->InsertNextPoint(this->InputPoints->GetPoint(this->GetPointIndexUV(u,v)));
    }

    this->BuildCoeffMatrix(this->InterpolationDegrees[0], uKnots, ukParams, points, matrixA);

    // Insert control points for each point column
    this->LuSolve(matrixA, this->InputResolution[0], points, controlPointsR);
  }
  this->DestructMatrix(matrixA, this->InputResolution[0]);

  //
  // Do global interpolation along the v direction
  //
  vtkNew<vtkPoints> controlPoints;
  matrixA = this->AllocateMatrix(this->InputResolution[1]);
  for (int u=0; u<this->InputResolution[0]; ++u)
  {
    // Collect row point indices
    vtkNew<vtkPoints> points;
    for (int v=0; v<this->InputResolution[1]; ++v)
    {
      points->InsertNextPoint(controlPointsR->GetPoint(v * this->InputResolution[0] + u));
    }

    this->BuildCoeffMatrix(this->InterpolationDegrees[1], vKnots, vlParams, points, matrixA);

    // Insert control points for each point row
    this->LuSolve(matrixA, this->InputResolution[1], points, controlPoints);
  }
  this->DestructMatrix(matrixA, this->InputResolution[1]);

  // Fill output
  vtkSmartPointer<vtkPoints> surfacePoints = vtkSmartPointer<vtkPoints>::New();
  //TODO: Fill output

  // Note from abstract.Surface: smaller the delta value, smoother the surface.
}

//---------------------------------------------------------------------------
void vtkSlicerNurbsFittingLogic::EvaluateSurface() // (self, datadict, **kwargs):
{
  // """ Evaluates the surface.
  // 
  // Keyword Arguments:
  //     * ``start``: starting parametric position for evaluation
  //     * ``stop``: ending parametric position for evaluation
  // 
  // :param datadict: data dictionary containing the necessary variables
  // :type datadict: dict
  // :return: evaluated points
  // :rtype: list
  // """
  // # Geometry data from datadict
  // sample_size = datadict['sample_size']
  // degree = datadict['degree']
  // knotvector = datadict['knotvector']
  // ctrlpts = datadict['control_points']
  // size = datadict['size']
  // dimension = datadict['dimension'] + 1 if datadict['rational'] else datadict['dimension']
  // pdimension = datadict['pdimension']
  // precision = datadict['precision']
  // 
  // # Keyword arguments
  // start = kwargs.get('start', [0.0 for _ in range(pdimension)])
  // stop = kwargs.get('stop', [1.0 for _ in range(pdimension)])
  // 
  // # Algorithm A3.5
  // spans = [[] for _ in range(pdimension)]
  // basis = [[] for _ in range(pdimension)]
  // for idx in range(pdimension):
  //     knots = linalg.linspace(start[idx], stop[idx], sample_size[idx], decimals=precision)
  //     spans[idx] = helpers.find_spans(degree[idx], knotvector[idx], size[idx], knots, self._span_func)
  //     basis[idx] = helpers.basis_functions(degree[idx], knotvector[idx], spans[idx], knots)
  // 
  // eval_points = []
  // for i in range(len(spans[0])):
  //     idx_u = spans[0][i] - degree[0]
  //     for j in range(len(spans[1])):
  //         idx_v = spans[1][j] - degree[1]
  //         spt = [0.0 for _ in range(dimension)]
  //         for k in range(0, degree[0] + 1):
  //             temp = [0.0 for _ in range(dimension)]
  //             for l in range(0, degree[1] + 1):
  //                 temp[:] = [tmp + (basis[1][j][l] * cp) for tmp, cp in
  //                            zip(temp, ctrlpts[idx_v + l + (size[1] * (idx_u + k))])]
  //             spt[:] = [pt + (basis[0][i][k] * tmp) for pt, tmp in zip(spt, temp)]
  // 
  //         eval_points.append(spt)
  // 
  // return eval_points
   
  //TODO:
}

//---------------------------------------------------------------------------
void vtkSlicerNurbsFittingLogic::ComputeParamsSurface(vtkDoubleArray* ukParams, vtkDoubleArray* vlParams)
{
  // :param points: data points
  // :type points: list, tuple
  // :param size_u: number of points on the u-direction
  // :type size_u: int
  // :param size_v: number of points on the v-direction
  // :type size_v: int
  // :param centripetal: activates centripetal parametrization method
  // :type centripetal: bool
  // :return: :math:`\\overline{u}_{k}` and :math:`\\overline{u}_{l}` parameter arrays as a tuple
  // :rtype: tuple
  // """
  // # Compute uk
  // uk = [0.0 for _ in range(size_u)]
  //
  // # Compute for each curve on the v-direction
  // uk_temp = []
  // for v in range(size_v):
  // #CP_Note: Get v'th column
  //     pts_u = [points[v + (size_v * u)] for u in range(size_u)]
  // #CP_Note: Append parameter list to one long list
  //     uk_temp += compute_params_curve(pts_u, centripetal)
  //
  // # Do averaging on the u-direction
  // for u in range(size_u):
  //     knots_v = [uk_temp[u + (size_u * v)] for v in range(size_v)]
  //     uk[u] = sum(knots_v) / size_v
  //
  // # Compute vl
  // vl = [0.0 for _ in range(size_v)]
  //
  // # Compute for each curve on the u-direction
  // vl_temp = []
  // #CP_Note: Get u'th row
  // for u in range(size_u):
  //     pts_v = [points[v + (size_v * u)] for v in range(size_v)]
  //     vl_temp += compute_params_curve(pts_v, centripetal)
  //
  // # Do averaging on the v-direction
  // for v in range(size_v):
  //     knots_u = [vl_temp[v + (size_v * u)] for u in range(size_u)]
  //     vl[v] = sum(knots_u) / size_u
  //
  // return uk, vl

  // Compute parameters for each curve on the v direction. Concatenate into one parameter array
  vtkNew<vtkDoubleArray> ukParamsTemp;
  for (int v=0; v<this->InputResolution[1]; ++v)
  {
    // Collect column point indices
    vtkNew<vtkIdList> columnIds;
    for (int u=0; u<this->InputResolution[0]; ++u)
    {
      columnIds->InsertNextId(this->GetPointIndexUV(u,v));
    }

    // Compute parameters for column
    this->ComputeParamsCurve(columnIds, ukParamsTemp);
  }

  // Do averaging on the u direction
  ukParams->Initialize();
  for (int u=0; u<this->InputResolution[0]; ++u)
  {
    double sumKnots_v = 0.0;
    for (int v=0; v<this->InputResolution[1]; ++v)
    {
      sumKnots_v += ukParamsTemp->GetValue(u + (this->InputResolution[0] * v));
    }
    ukParams->InsertNextValue(sumKnots_v / this->InputResolution[1]);
  }

  // Compute parameters for each curve on the u direction. Concatenate into one parameter array
  vtkNew<vtkDoubleArray> vlParamsTemp;
  for (int u=0; u<this->InputResolution[0]; ++u)
  {
    // Collect row point indices
    vtkNew<vtkIdList> rowIds;
    for (int v=0; v<this->InputResolution[1]; ++v)
    {
      rowIds->InsertNextId(this->GetPointIndexUV(u,v));
    }

    // Compute parameters for row
    this->ComputeParamsCurve(rowIds, vlParamsTemp);
  }

  // Do averaging on the v direction
  vlParams->Initialize();
  for (int v=0; v<this->InputResolution[1]; ++v)
  {
    double sumKnots_u = 0.0;
    for (int u=0; u<this->InputResolution[0]; ++u)
    {
      sumKnots_u += vlParamsTemp->GetValue(v + (this->InputResolution[1] * u));
    }
    vlParams->InsertNextValue(sumKnots_u / this->InputResolution[0]);
  }
}

//---------------------------------------------------------------------------
void vtkSlicerNurbsFittingLogic::ComputeParamsCurve(
  vtkIdList* pointIndexList, vtkDoubleArray* outParametersArray)
{
  // :param points: data points
  // :type points: list, tuple
  // :param centripetal: activates centripetal parametrization method
  // :type centripetal: bool
  // :return: parameter array, :math:`\\overline{u}_{k}`
  // :rtype: list
  // """
  // if not isinstance(points, (list, tuple)):
  //     raise TypeError("Data points must be a list or a tuple")
  //
  // # Length of the points array
  // num_points = len(points)
  //
  // # Calculate chord lengths
  // cds = [0.0 for _ in range(num_points + 1)]
  // cds[-1] = 1.0
  // for i in range(1, num_points):
  // #CP_Note: linalg.point_distance
  //     distance = linalg.point_distance(points[i], points[i - 1])
  //     cds[i] = math.sqrt(distance) if centripetal else distance
  //
  // # Find the total chord length
  // d = sum(cds[1:-1])
  //
  // # Divide individual chord lengths by the total chord length
  // uk = [0.0 for _ in range(num_points)]
  // for i in range(num_points):
  //     uk[i] = sum(cds[0:i + 1]) / d
  //
  // return uk
  
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
    this->InputPoints->GetPoint(pointIndexList->GetId(i-1), prevPoint);
    this->InputPoints->GetPoint(pointIndexList->GetId(i), currentPoint);
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
void vtkSlicerNurbsFittingLogic::ComputeKnotVector(
  int degree, int numOfPoints, vtkDoubleArray* inParams, vtkDoubleArray* outKnotVector)
{
  // :param degree: degree
  // :type degree: int
  // :param num_points: number of data points
  // :type num_points: int
  // :param params: list of parameters, :math:`\\overline{u}_{k}`
  // :type params: list, tuple
  // :return: knot vector
  // :rtype: list
  // """
  // # Start knot vector
  // kv = [0.0 for _ in range(degree + 1)]
  //
  // # Use averaging method (Eqn 9.8) to compute internal knots in the knot vector
  // for i in range(num_points - degree - 1):
  //     temp_kv = (1.0 / degree) * sum([params[j] for j in range(i + 1, i + degree + 1)])
  //     kv.append(temp_kv)
  //
  // # End knot vector
  // kv += [1.0 for _ in range(degree + 1)]
  //
  // return kv

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
void vtkSlicerNurbsFittingLogic::BuildCoeffMatrix(
  int degree, vtkDoubleArray* knotVector, vtkDoubleArray* params, vtkPoints* points, double** outCoeffMatrix)
{
  // :param degree: degree
  // :type degree: int
  // :param knotvector: knot vector
  // :type knotvector: list, tuple
  // :param params: list of parameters
  // :type params: list, tuple
  // :param points: data points
  // :type points: list, tuple
  // :return: coefficient matrix
  // :rtype: list
  // """
  // # Number of data points
  // num_points = len(points)
  //
  // # Set up coefficient matrix
  // matrix_a = [[0.0 for _ in range(num_points)] for _ in range(num_points)]
  // for i in range(num_points):
  // #CP_Note: helpers.find_span_linear
  //     span = helpers.find_span_linear(degree, knotvector, num_points, params[i])
  // #CP_Note: helpers.basis_function
  //     matrix_a[i][span-degree:span+1] = helpers.basis_function(degree, knotvector, span, params[i])
  //
  // # Return coefficient matrix
  // return matrix_a

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
void vtkSlicerNurbsFittingLogic::BasisFunction(int degree, vtkDoubleArray* knotVector, int span, double knot, vtkDoubleArray* outBasisFunctions)
{
  // :param degree: degree, :math:`p`
  // :type degree: int
  // :param knot_vector: knot vector, :math:`U`
  // :type knot_vector: list, tuple
  // :param span: knot span, :math:`i`
  // :type span: int
  // :param knot: knot or parameter, :math:`u`
  // :type knot: float
  // :return: basis functions
  // :rtype: list
  // """
  // left = [0.0 for _ in range(degree + 1)]
  // right = [0.0 for _ in range(degree + 1)]
  // N = [1.0 for _ in range(degree + 1)]  # N[0] = 1.0 by definition
  //
  // for j in range(1, degree + 1):
  //     left[j] = knot - knot_vector[span + 1 - j]
  //     right[j] = knot_vector[span + j] - knot
  //     saved = 0.0
  //     for r in range(0, j):
  //         temp = N[r] / (right[r + 1] + left[j - r])
  //         N[r] = saved + right[r + 1] * temp
  //         saved = left[j - r] * temp
  //     N[j] = saved
  //
  // return N

  if (!outBasisFunctions)
  {
    vtkErrorMacro("GenerateKnotVector: Invalid output array");
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
int vtkSlicerNurbsFittingLogic::FindSpanLinear(int degree, vtkDoubleArray* knotVector, int numControlPoints, double knot)
{
  // :param degree: degree, :math:`p`
  // :type degree: int
  // :param knot_vector: knot vector, :math:`U`
  // :type knot_vector: list, tuple
  // :param num_ctrlpts: number of control points, :math:`n + 1`
  // :type num_ctrlpts: int
  // :param knot: knot or parameter, :math:`u`
  // :type knot: float
  // :return: knot span
  // :rtype: int
  // """
  // span = degree + 1  # Knot span index starts from zero
  // while span < num_ctrlpts and knot_vector[span] <= knot:
  //     span += 1
  //
  // return span - 1
 
  // Knot span index starts from zero
  int span = degree + 1;
  
  while (span < numControlPoints && knotVector->GetValue(span) <= knot)
  {
    span++;
  }

  return span - 1;
}

//
// linalg
//

//---------------------------------------------------------------------------
void vtkSlicerNurbsFittingLogic::LuSolve(double** coeffMatrix, int matrixSize, vtkPoints* points, vtkPoints* outControlPointsR)
{
  // :param matrix_a: matrix A
  // :type m_l: list
  // :param b: matrix of M column vectors
  // :type b: list
  // :return: x, the solution matrix
  // :rtype: list
  // """
  // # Variable initialization
  // dim = len(b[0])
  // num_x = len(b)
  // x = [[0.0 for _ in range(dim)] for _ in range(num_x)]
  //
  // # LU decomposition
  // m_l, m_u = lu_decomposition(matrix_a)
  //
  // # Solve the system of linear equations
  // for i in range(dim):
  //     bt = [b1[i] for b1 in b]
  //     y = forward_substitution(m_l, bt)
  //     xt = backward_substitution(m_u, y)
  //     for j in range(num_x):
  //         x[j][i] = xt[j]
  //
  // # Return the solution
  // return x

  if (!coeffMatrix || !points)
  {
    vtkErrorMacro("LuSolve: Invalid input coefficients matrix or points");
    return;
  }

  int dim = 3;
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
void vtkSlicerNurbsFittingLogic::LuDecomposition(double** matrixA, double** matrixL, double** matrixU, int size)
{
  // :param matrix_a: Input matrix (must be a square matrix)
  // :type matrix_a: list, tuple
  // :return: a tuple containing matrices L and U
  // :rtype: tuple
  // """
  // # Check if the 2-dimensional input matrix is a square matrix
  // q = len(matrix_a)
  // for idx, m_a in enumerate(matrix_a):
  //     if len(m_a) != q:
  //         raise ValueError("The input must be a square matrix. " +
  //                          "Row " + str(idx + 1) + " has a size of " + str(len(m_a)) + ".")
  //
  // # Return L and U matrices
  // return _linalg.doolittle(matrix_a)

  // Porting note: L and U matrices differ from the Python implementation but their product is exactly matrix A

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
void vtkSlicerNurbsFittingLogic::ForwardSubstitution(double** matrixL, double* b, int size, double* outY)
{
  // :param matrix_l: L, lower triangular matrix
  // :type matrix_l: list, tuple
  // :param matrix_b: b, column matrix
  // :type matrix_b: list, tuple
  // :return: y, column matrix
  // :rtype: list
  // """
  // q = len(matrix_b)
  // matrix_y = [0.0 for _ in range(q)]
  // matrix_y[0] = float(matrix_b[0]) / float(matrix_l[0][0])
  // for i in range(1, q):
  //     matrix_y[i] = float(matrix_b[i]) - sum([matrix_l[i][j] * matrix_y[j] for j in range(0, i)])
  //     matrix_y[i] /= float(matrix_l[i][i])
  // return matrix_y

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
void vtkSlicerNurbsFittingLogic::BackwardSubstitution(double** matrixU, double* y, int size, double* outX)
{
  // :param matrix_u: U, upper triangular matrix
  // :type matrix_u: list, tuple
  // :param matrix_y: y, column matrix
  // :type matrix_y: list, tuple
  // :return: x, column matrix
  // :rtype: list
  // """
  // q = len(matrix_y)
  // matrix_x = [0.0 for _ in range(q)]
  // matrix_x[q - 1] = float(matrix_y[q - 1]) / float(matrix_u[q - 1][q - 1])
  // for i in range(q - 2, -1, -1):
  //     matrix_x[i] = float(matrix_y[i]) - sum([matrix_u[i][j] * matrix_x[j] for j in range(i, q)])
  //     matrix_x[i] /= float(matrix_u[i][i])
  // return matrix_x

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
unsigned int vtkSlicerNurbsFittingLogic::GetPointIndexUV(unsigned int u, unsigned int v)
{
  if (u >= this->InputResolution[0] || v >= this->InputResolution[1])
  {
    vtkErrorMacro("GetPointUV: Index pair (" << u << ", " << v << ") is out of range");
    return -1;
  }

  return u * this->InputResolution[1] + v;
}

//---------------------------------------------------------------------------
double** vtkSlicerNurbsFittingLogic::AllocateMatrix(int m, int n/*=0*/)
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
void vtkSlicerNurbsFittingLogic::DestructMatrix(double** matrix, int m, int n/*=0*/)
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
