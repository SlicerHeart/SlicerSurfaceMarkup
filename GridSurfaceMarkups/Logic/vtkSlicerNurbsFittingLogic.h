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

#ifndef __vtkSlicerNurbsFittingLogic_h_
#define __vtkSlicerNurbsFittingLogic_h_

#include <vtkPolyDataAlgorithm.h>

#include "vtkSlicerGridSurfaceMarkupsModuleLogicExport.h"

class vtkDoubleArray;
class vtkIdList;

class VTK_SLICER_GRIDSURFACEMARKUPS_MODULE_LOGIC_EXPORT vtkSlicerNurbsFittingLogic : public vtkPolyDataAlgorithm
{
public:
  static vtkSlicerNurbsFittingLogic* New();
  vtkTypeMacro(vtkSlicerNurbsFittingLogic, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  /// TODO:
  void SetInputPoints(vtkPoints* points);
  /// TODO:
  vtkSmartPointer<vtkPoints> GetInputPoints() const { return this->InputPoints.GetPointer(); };

  /// TODO:
  vtkSetVector2Macro(InputResolution, unsigned int);
  /// TODO:
  vtkGetVector2Macro(InputResolution, unsigned int);
  /// TODO:
  vtkSetVector2Macro(InterpolationDegrees, unsigned int);
  /// TODO:
  vtkGetVector2Macro(InterpolationDegrees, unsigned int);
  /// TODO:
  vtkSetMacro(UseCentripetal, bool);
  vtkBooleanMacro(UseCentripetal, bool);
  /// TODO:
  vtkGetMacro(UseCentripetal, bool);

  /// Function computing the NURBS surface according to the pipeline architecture of VTK.
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

//protected: //TODO: FOR DEBUGGING EVERYTHING IS PULIC, UNCOMMENT THIS LINE
  /// Compute NURBS surface poly data from the input points according to input resolution and degrees
  void UpdateNurbsPolyData(vtkPolyData* polyData);

  /// TODO:
  void ComputeParamsSurface(vtkDoubleArray* ukParams, vtkDoubleArray* vlParams);
  /// TODO:
  void ComputeParamsCurve(vtkIdList* indexList, vtkDoubleArray* parametersArray);
  /// TODO:
  void ComputeKnotVector(int degree, int numberOfControlPoints, vtkDoubleArray* params, vtkDoubleArray* outKnotVector);
  /// TODO:
  void BuildCoeffMatrix();

  /// TODO: helpers
  void BasisFunction();
  /// TODO: helpers
  void FindSpanLinear();

  /// TODO: linalg
  void LuSolve();
  /// TODO: linalg
  void LuDecomposition();
  /// TODO: linalg
  void ForwardSubstitution();
  /// TODO: linalg
  void BackwardSubstitution();
  /// TODO: linalg
  //void LinSpace(double start, double stop, int numOfSamples, vtkDoubleArray* inOutArray);
  /// TODO: _linalg
  void DooLittle();
  /// TODO: knotvector
  //void GenerateKnotVector(int degree, int numberOfControlPoints, bool clamped=true);
  /// TODO: linalg
  //void PointDistance();
   
  /// Convenience function to get point index from input point list with the two (u,v) indices
  unsigned int GetPointIndexUV(unsigned int u, unsigned int v);

protected:
  /// Input control points. The number of points is u*v, and the strides contain the rows (0:[0,0], 1:[0,1], ...)
  vtkSmartPointer<vtkPoints> InputPoints;

  /// Number of data points along the u and v directions, respectively
  unsigned int InputResolution[2] = {4,4};
  /// Degree of the output surface for the u and v directions, respectively
  unsigned int InterpolationDegrees[2] = {3,3};

  /// Activate centripetal parametrization method. Default: false
  bool UseCentripetal = false;

protected:
  vtkSlicerNurbsFittingLogic();
  ~vtkSlicerNurbsFittingLogic() override;

private:
  vtkSlicerNurbsFittingLogic(const vtkSlicerNurbsFittingLogic&) = delete;
  void operator=(const vtkSlicerNurbsFittingLogic&) = delete;
};

#endif // __vtkSlicerNurbsFittingLogic_h_