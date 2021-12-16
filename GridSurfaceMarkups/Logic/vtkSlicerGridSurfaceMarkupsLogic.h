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

#ifndef __vtkSlicerGridSurfaceMarkupsLogic_h_
#define __vtkSlicerGridSurfaceMarkupsLogic_h_

#include <vtkSlicerModuleLogic.h>

#include "vtkSlicerGridSurfaceMarkupsModuleLogicExport.h"

class VTK_SLICER_GRIDSURFACEMARKUPS_MODULE_LOGIC_EXPORT vtkSlicerGridSurfaceMarkupsLogic:
  public vtkSlicerModuleLogic
{
public:
  static vtkSlicerGridSurfaceMarkupsLogic* New();
  vtkTypeMacro(vtkSlicerGridSurfaceMarkupsLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent) override;

protected:
  vtkSlicerGridSurfaceMarkupsLogic();
  ~vtkSlicerGridSurfaceMarkupsLogic() override;

  void RegisterNodes() override;

private:
  vtkSlicerGridSurfaceMarkupsLogic(const vtkSlicerGridSurfaceMarkupsLogic&) = delete;
  void operator=(const vtkSlicerGridSurfaceMarkupsLogic&) = delete;
};

#endif // __vtkSlicerGridSurfaceMarkupsLogic_h_
