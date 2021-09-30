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

#ifndef __qSlicerMarkupsAdditionalGridSurfaceOptionsWidget_h
#define __qSlicerMarkupsAdditionalGridSurfaceOptionsWidget_h

// Qt includes
#include <QWidget>

// Markups widgets includes
#include "qSlicerMarkupsAdditionalOptionsWidget.h"

// Grid Surface widgets includes
#include "qSlicerGridSurfaceMarkupsModuleWidgetsExport.h"

// CTK includes
#include <ctkPimpl.h>
#include <ctkVTKObject.h>

class vtkMRMLNode;
class vtkMRMLMarkupsGridSurfaceNode;
class qSlicerMarkupsAdditionalGridSurfaceOptionsWidget;
class qSlicerMarkupsAdditionalGridSurfaceOptionsWidgetPrivate;

class Q_SLICER_MODULE_GRIDSURFACEMARKUPS_WIDGETS_EXPORT qSlicerMarkupsAdditionalGridSurfaceOptionsWidget : public qSlicerMarkupsAdditionalOptionsWidget
{
  Q_OBJECT
  QVTK_OBJECT

public:
  typedef qSlicerMarkupsAdditionalOptionsWidget Superclass;
  qSlicerMarkupsAdditionalGridSurfaceOptionsWidget(QWidget* parent=nullptr);
  ~qSlicerMarkupsAdditionalGridSurfaceOptionsWidget() override;

  /// Returns the current MRML Grid Surface node
  vtkMRMLMarkupsGridSurfaceNode* mrmlGridSurfaceNode()const;

  /// Gets the name of the additional options widget type
  const QString getAdditionalOptionsWidgetTypeName() override { return "Grid Surface"; }

  /// Updates the widget based on information from MRML.
  void updateWidgetFromMRML() override;

  /// Checks whether a given node can be handled by the widget
  bool canManageMRMLMarkupsNode(vtkMRMLMarkupsNode *markupsNode) const override;

public slots:

  /// Set the MRML node of interest
  void setMRMLMarkupsNode(vtkMRMLMarkupsNode* node) override;

  /// Sets the vtkMRMLMarkupsNode to operate on.
  void setMRMLMarkupsNode(vtkMRMLNode* node) override;

protected slots:
  /// Internal function to update the widgets based on the Grid Surface node
  void onMRMLNodeModified();
  /// Internal function to update type of Grid Surface
  void onGridSurfaceTypeParameterChanged();
  /// Handle apply grid surface button click
  void onApplyGridResolution();

protected:
  qSlicerMarkupsAdditionalGridSurfaceOptionsWidget(qSlicerMarkupsAdditionalGridSurfaceOptionsWidgetPrivate &d, QWidget* parent=nullptr);
  void setup();

private:
  Q_DECLARE_PRIVATE(qSlicerMarkupsAdditionalGridSurfaceOptionsWidget);
  Q_DISABLE_COPY(qSlicerMarkupsAdditionalGridSurfaceOptionsWidget);
};

#endif
