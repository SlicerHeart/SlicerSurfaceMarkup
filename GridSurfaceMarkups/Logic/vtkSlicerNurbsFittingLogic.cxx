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
// linalg.lu_solve
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
// linalg.point_distance
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
// helpers.find_span_linear
  //     span = helpers.find_span_linear(degree, knotvector, num_points, params[i])
// helpers.basis_function
  //     matrix_a[i][span-degree:span+1] = helpers.basis_function(degree, knotvector, span, params[i])
  //
  // # Return coefficient matrix
  // return matrix_a
}

//
// helpers
//

//---------------------------------------------------------------------------
void vtkSlicerNurbsFittingLogic::BasisFunction() // (degree, knot_vector, span, knot)
{
  // """ Computes the non-vanishing basis functions for a single parameter.
  //
  // Implementation of Algorithm A2.2 from The NURBS Book by Piegl & Tiller.
  // Uses recurrence to compute the basis functions, also known as Cox - de
  // Boor recursion formula.
  //
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

}

//---------------------------------------------------------------------------
void vtkSlicerNurbsFittingLogic::FindSpanLinear() // (degree, knot_vector, num_ctrlpts, knot, **kwargs)
{
  // """ Finds the span of a single knot over the knot vector using linear search.
  //
  // Alternative implementation for the Algorithm A2.1 from The NURBS Book by Piegl & Tiller.
  //
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

}

//
// linalg
//

//---------------------------------------------------------------------------
void vtkSlicerNurbsFittingLogic::LuSolve() // (matrix_a, b)
{
  // """ Computes the solution to a system of linear equations.
  //
  // This function solves :math:`Ax = b` using LU decomposition. :math:`A` is a
  // :math:`N \\times N` matrix, :math:`b` is :math:`N \\times M` matrix of
  // :math:`M` column vectors. Each column of :math:`x` is a solution for
  // corresponding column of :math:`b`.
  //
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

}

//---------------------------------------------------------------------------
void vtkSlicerNurbsFittingLogic::LuDecomposition() // (matrix_a)
{
  // """ LU-Factorization method using Doolittle's Method for solution of linear systems.
  //
  // Decomposes the matrix :math:`A` such that :math:`A = LU`.
  //
  // The input matrix is represented by a list or a tuple. The input matrix is **2-dimensional**, i.e. list of lists of
  // integers and/or floats.
  //
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

}

//---------------------------------------------------------------------------
void vtkSlicerNurbsFittingLogic::ForwardSubstitution()
{
  // """ Forward substitution method for the solution of linear systems.
  //
  // Solves the equation :math:`Ly = b` using forward substitution method
  // where :math:`L` is a lower triangular matrix and :math:`b` is a column matrix.
  //
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

}

//---------------------------------------------------------------------------
void vtkSlicerNurbsFittingLogic::BackwardSubstitution()
{
  // """ Backward substitution method for the solution of linear systems.
  //
  // Solves the equation :math:`Ux = y` using backward substitution method
  // where :math:`U` is a upper triangular matrix and :math:`y` is a column matrix.
  //
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

}

//---------------------------------------------------------------------------
void vtkSlicerNurbsFittingLogic::LinSpace()
{
  // """ Returns a list of evenly spaced numbers over a specified interval.
  //
  // Inspired from Numpy's linspace function: https://github.com/numpy/numpy/blob/master/numpy/core/function_base.py
  //
  // :param start: starting value
  // :type start: float
  // :param stop: end value
  // :type stop: float
  // :param num: number of samples to generate
  // :type num: int
  // :param decimals: number of significands
  // :type decimals: int
  // :return: a list of equally spaced numbers
  // :rtype: list
  // """
  // start = float(start)
  // stop = float(stop)
  // if abs(start - stop) <= 10e-8:
  //     return [start]
  // num = int(num)
  // if num > 1:
  //     div = num - 1
  //     delta = stop - start
  //     return [float(("{:." + str(decimals) + "f}").format((start + (float(x) * float(delta) / float(div)))))
  //             for x in range(num)]
  // return [float(("{:." + str(decimals) + "f}").format(start))]


}

void vtkSlicerNurbsFittingLogic::DooLittle() // (matrix_a):
{
  // """ Doolittle's Method for LU-factorization.
  //
  // :param matrix_a: Input matrix (must be a square matrix)
  // :type matrix_a: list, tuple
  // :return: a tuple containing matrices (L,U)
  // :rtype: tuple
  // """
  // # Initialize L and U matrices
  // matrix_u = [[0.0 for _ in range(len(matrix_a))] for _ in range(len(matrix_a))]
  // matrix_l = [[0.0 for _ in range(len(matrix_a))] for _ in range(len(matrix_a))]
  //
  // # Doolittle Method
  // for i in range(0, len(matrix_a)):
  //     for k in range(i, len(matrix_a)):
  //         # Upper triangular (U) matrix
  //         matrix_u[i][k] = float(matrix_a[i][k] - sum([matrix_l[i][j] * matrix_u[j][k] for j in range(0, i)]))
  //         # Lower triangular (L) matrix
  //         if i == k:
  //             matrix_l[i][i] = 1.0
  //         else:
  //             matrix_l[k][i] = float(matrix_a[k][i] - sum([matrix_l[k][j] * matrix_u[j][i] for j in range(0, i)]))
  //             # Handle zero division error
  //             try:
  //                 matrix_l[k][i] /= float(matrix_u[i][i])
  //             except ZeroDivisionError:
  //                 matrix_l[k][i] = 0.0
  //
  // return matrix_l, matrix_u

}

//---------------------------------------------------------------------------
void vtkSlicerNurbsFittingLogic::GenerateKnotVector() // knotvector(degree, num_ctrlpts, **kwargs)
{
  // """ Generates an equally spaced knot vector.
  //
  // It uses the following equality to generate knot vector: :math:`m = n + p + 1`
  //
  // where;
  //
  // * :math:`p`, degree
  // * :math:`n + 1`, number of control points
  // * :math:`m + 1`, number of knots
  //
  // Keyword Arguments:
  //
  //     * ``clamped``: Flag to choose from clamped or unclamped knot vector options. *Default: True*
  //
  // :param degree: degree
  // :type degree: int
  // :param num_ctrlpts: number of control points
  // :type num_ctrlpts: int
  // :return: knot vector
  // :rtype: list
  // """
  // if degree == 0 or num_ctrlpts == 0:
  //     raise ValueError("Input values should be different than zero.")
  //
  // # Get keyword arguments
  // clamped = kwargs.get('clamped', True)
  //
  // # Number of repetitions at the start and end of the array
  // num_repeat = degree
  //
  // # Number of knots in the middle
  // num_segments = num_ctrlpts - (degree + 1)
  //
  // if not clamped:
  //     # No repetitions at the start and end
  //     num_repeat = 0
  //     # Should conform the rule: m = n + p + 1
  //     num_segments = degree + num_ctrlpts - 1
  //
  // # First knots
  // knot_vector = [0.0 for _ in range(0, num_repeat)]
  //
  // # Middle knots
  // knot_vector += linspace(0.0, 1.0, num_segments + 2)
  //
  // # Last knots
  // knot_vector += [1.0 for _ in range(0, num_repeat)]
  //
  // # Return auto-generated knot vector
  // return knot_vector

}

//---------------------------------------------------------------------------
//void vtkSlicerNurbsFittingLogic::PointDistance() // (pt1, pt2)
//{
  // """ Computes distance between two points.
  //
  // :param pt1: point 1
  // :type pt1: list, tuple
  // :param pt2: point 2
  // :type pt2: list, tuple
  // :return: distance between input points
  // :rtype: float
  // """
  // if len(pt1) != len(pt2):
  //     raise ValueError("The input points should have the same dimension")
  //
  // dist_vector = vector_generate(pt1, pt2, normalize=False)
  // distance = vector_magnitude(dist_vector)
  // return distance

//}

