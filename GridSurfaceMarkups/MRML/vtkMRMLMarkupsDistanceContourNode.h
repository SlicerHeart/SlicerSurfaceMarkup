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

#ifndef __vtkmrmlmarkupsdistancecontournode_h_
#define __vtkmrmlmarkupsdistancecontournode_h_

#include "vtkSlicerGridSurfaceMarkupsModuleMRMLExport.h"

// MRML includes
#include <vtkMRMLMarkupsLineNode.h>
#include <vtkMRMLModelNode.h>

//VTK includes
#include <vtkWeakPointer.h>

//-----------------------------------------------------------------------------
class VTK_SLICER_GRIDSURFACEMARKUPS_MODULE_MRML_EXPORT vtkMRMLMarkupsDistanceContourNode
: public vtkMRMLMarkupsLineNode
{
public:
  static vtkMRMLMarkupsDistanceContourNode* New();
  vtkTypeMacro(vtkMRMLMarkupsDistanceContourNode, vtkMRMLMarkupsLineNode);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  //--------------------------------------------------------------------------------
  // MRMLNode methods
  //--------------------------------------------------------------------------------
  const char* GetIcon() override {return ":/Icons/MarkupsGeneric.png";}
  const char* GetAddIcon() override {return ":/Icons/MarkupsGenericMouseModePlace.png";}
  const char* GetPlaceAddIcon() override {return ":/Icons/MarkupsGenericMouseModePlaceAdd.png";}

  vtkMRMLNode* CreateNodeInstance() override;

  /// Get node XML tag name (like Volume, Model)
  ///
  const char* GetNodeTagName() override {return "MarkupsDistanceContour";}

  /// Get markup name
  const char* GetMarkupType() override {return "DistanceContour";}

  /// Get markup short name
  const char* GetDefaultNodeNamePrefix() override {return "SC";}

  /// \sa vtkMRMLNode::CopyContent
  vtkMRMLCopyContentDefaultMacro(vtkMRMLMarkupsDistanceContourNode);

  vtkMRMLModelNode* GetTarget() const {return this->Target;}
  void SetTarget(vtkMRMLModelNode* target) {this->Target = target; this->Modified();}

protected:
  vtkMRMLMarkupsDistanceContourNode();
  ~vtkMRMLMarkupsDistanceContourNode() override = default;

private:
 vtkWeakPointer<vtkMRMLModelNode> Target;

private:
 vtkMRMLMarkupsDistanceContourNode(const vtkMRMLMarkupsDistanceContourNode&);
 void operator=(const vtkMRMLMarkupsDistanceContourNode&);


};

#endif //__vtkmrmlmarkupsdistancecontournode_h_
