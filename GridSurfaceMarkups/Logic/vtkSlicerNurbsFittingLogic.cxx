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
#include <vtkExecutive.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
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
void vtkSlicerNurbsFittingLogic::UpdateNurbsPolyData(vtkPolyData* polyData)
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
  //     pts = [points[v + (size_v * u)] for u in range(size_u)]
  //     matrix_a = _build_coeff_matrix(degree_u, kv_u, uk, pts)
  //     ctrlpts_r += linalg.lu_solve(matrix_a, pts)
  //
  // # Do global interpolation on the v-direction
  // ctrlpts = []
  // for u in range(size_u):
  //     pts = [ctrlpts_r[u + (size_u * v)] for v in range(size_v)]
  //     matrix_a = _build_coeff_matrix(degree_v, kv_v, vl, pts)
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

  vtkSmartPointer<vtkPoints> surfacePoints = vtkSmartPointer<vtkPoints>::New();

}

//---------------------------------------------------------------------------
void vtkSlicerNurbsFittingLogic::ComputeParamsSurface() // (points, size_u, size_v, degree_u, degree_v, **kwargs)
{
  // """ Computes :math:`\\overline{u}_{k}` and :math:`\\overline{u}_{l}` for surfaces.
  //
  // The data points array has a row size of ``size_v`` and column size of ``size_u`` and it is 1-dimensional. Please
  // refer to The NURBS Book (2nd Edition), pp.366-367 for details on how to compute :math:`\\overline{u}_{k}` and
  // :math:`\\overline{u}_{l}` arrays for global surface interpolation.
  //
  // Please note that this function is not a direct implementation of Algorithm A9.3 which can be found on The NURBS Book
  // (2nd Edition), pp.377-378. However, the output is the same.
  //
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
  //     pts_u = [points[v + (size_v * u)] for u in range(size_u)]
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
}

//---------------------------------------------------------------------------
void vtkSlicerNurbsFittingLogic::ComputeParamsCurve() // (points, centripetal=False)
{
  // """ Computes :math:`\\overline{u}_{k}` for curves.
  //
  // Please refer to the Equations 9.4 and 9.5 for chord length parametrization, and Equation 9.6 for centripetal method
  // on The NURBS Book (2nd Edition), pp.364-365.
  //
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
}

//---------------------------------------------------------------------------
void vtkSlicerNurbsFittingLogic::ComputeKnotVector() // (degree, num_points, params)
{
  // """ Computes a knot vector from the parameter list using averaging method.
  //
  // Please refer to the Equation 9.8 on The NURBS Book (2nd Edition), pp.365 for details.
  //
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
}

//---------------------------------------------------------------------------
void vtkSlicerNurbsFittingLogic::BuildCoeffMatrix() // (degree, knotvector, params, points)
{
  // """ Builds the coefficient matrix for global interpolation.
  //
  // This function only uses data points to build the coefficient matrix. Please refer to The NURBS Book (2nd Edition),
  // pp364-370 for details.
  //
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
  //     span = helpers.find_span_linear(degree, knotvector, num_points, params[i])
  //     matrix_a[i][span-degree:span+1] = helpers.basis_function(degree, knotvector, span, params[i])
  //
  // # Return coefficient matrix
  // return matrix_a
}
