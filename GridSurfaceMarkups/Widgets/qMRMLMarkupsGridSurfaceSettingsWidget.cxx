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

// qMRML includes
#include "qMRMLMarkupsGridSurfaceSettingsWidget.h"
#include "ui_qMRMLMarkupsGridSurfaceSettingsWidget.h"

// MRML includes
#include <vtkMRMLDisplayNode.h>

// GridSurfaceMarkups includes
#include <vtkMRMLMarkupsGridSurfaceNode.h>

// --------------------------------------------------------------------------
class qMRMLMarkupsGridSurfaceSettingsWidgetPrivate:
  public Ui_qMRMLMarkupsGridSurfaceSettingsWidget
{
  Q_DECLARE_PUBLIC(qMRMLMarkupsGridSurfaceSettingsWidget);

public:
  qMRMLMarkupsGridSurfaceSettingsWidgetPrivate(qMRMLMarkupsGridSurfaceSettingsWidget* object);
  void setupUi(qMRMLMarkupsGridSurfaceSettingsWidget* widget);

protected:
  qMRMLMarkupsGridSurfaceSettingsWidget* const q_ptr;
};

// --------------------------------------------------------------------------
qMRMLMarkupsGridSurfaceSettingsWidgetPrivate::qMRMLMarkupsGridSurfaceSettingsWidgetPrivate(qMRMLMarkupsGridSurfaceSettingsWidget* object)
  : q_ptr(object)
{
}

// --------------------------------------------------------------------------
void qMRMLMarkupsGridSurfaceSettingsWidgetPrivate::setupUi(qMRMLMarkupsGridSurfaceSettingsWidget* widget)
{
  Q_Q(qMRMLMarkupsGridSurfaceSettingsWidget);

  this->Ui_qMRMLMarkupsGridSurfaceSettingsWidget::setupUi(widget);

  this->surfaceTypeComboBox->clear();
  for (int gridSurfaceType = 0; gridSurfaceType < vtkMRMLMarkupsGridSurfaceNode::GridSurfaceType_Last; ++gridSurfaceType)
  {
    this->surfaceTypeComboBox->addItem(vtkMRMLMarkupsGridSurfaceNode::GetGridSurfaceTypeAsString(gridSurfaceType), gridSurfaceType);
  }

  QObject::connect(this->surfaceTypeComboBox, SIGNAL(currentIndexChanged(int)), q, SLOT(onGridSurfaceTypeParameterChanged()));
  QObject::connect(this->applyGridResolutionButton, SIGNAL(clicked()), q, SLOT(onApplyGridResolution()));
  q->setEnabled(q->MarkupsNode != nullptr);
}

// --------------------------------------------------------------------------
// qMRMLMarkupsGridSurfaceSettingsWidget methods

// --------------------------------------------------------------------------
qMRMLMarkupsGridSurfaceSettingsWidget::qMRMLMarkupsGridSurfaceSettingsWidget(QWidget* parent)
  : Superclass(parent)
  , d_ptr(new qMRMLMarkupsGridSurfaceSettingsWidgetPrivate(this))
{
  this->setup();
}

// --------------------------------------------------------------------------
qMRMLMarkupsGridSurfaceSettingsWidget::~qMRMLMarkupsGridSurfaceSettingsWidget() = default;

// --------------------------------------------------------------------------
void qMRMLMarkupsGridSurfaceSettingsWidget::setup()
{
  Q_D(qMRMLMarkupsGridSurfaceSettingsWidget);
  d->setupUi(this);
}

// --------------------------------------------------------------------------
void qMRMLMarkupsGridSurfaceSettingsWidget::setMRMLMarkupsNode(vtkMRMLMarkupsNode* markupsNode)
{
  Q_D(qMRMLMarkupsGridSurfaceSettingsWidget);

  this->MarkupsNode = vtkMRMLMarkupsGridSurfaceNode::SafeDownCast(markupsNode);
  this->setEnabled(this->MarkupsNode!= nullptr);
}

//-----------------------------------------------------------------------------
bool qMRMLMarkupsGridSurfaceSettingsWidget::canManageMRMLMarkupsNode(vtkMRMLMarkupsNode* markupsNode)const
{
  Q_D(const qMRMLMarkupsGridSurfaceSettingsWidget);

  vtkMRMLMarkupsGridSurfaceNode* gridSurfaceNode = vtkMRMLMarkupsGridSurfaceNode::SafeDownCast(markupsNode);
  if (!gridSurfaceNode)
  {
    return false;
  }

  return true;
}

// --------------------------------------------------------------------------
void qMRMLMarkupsGridSurfaceSettingsWidget::updateWidgetFromMRML()
{
  Q_D(qMRMLMarkupsGridSurfaceSettingsWidget);

  vtkMRMLMarkupsGridSurfaceNode* markupsGridSurfaceNode = vtkMRMLMarkupsGridSurfaceNode::SafeDownCast(this->MarkupsNode);
  if (!markupsGridSurfaceNode)
  {
    return;
  }

  bool wasBlocked = d->surfaceTypeComboBox->blockSignals(true);
  d->surfaceTypeComboBox->setCurrentIndex(d->surfaceTypeComboBox->findData(markupsGridSurfaceNode->GetGridSurfaceType()));
  d->surfaceTypeComboBox->blockSignals(wasBlocked);

  int gridResolution[2] = {0};
  markupsGridSurfaceNode->GetGridResolution(gridResolution);
  wasBlocked = d->gridResolutionSpinBox_X->blockSignals(true);
  d->gridResolutionSpinBox_X->setValue(gridResolution[0]);
  d->gridResolutionSpinBox_X->blockSignals(wasBlocked);
  wasBlocked = d->gridResolutionSpinBox_Y->blockSignals(true);
  d->gridResolutionSpinBox_Y->setValue(gridResolution[1]);
  d->gridResolutionSpinBox_Y->blockSignals(wasBlocked);
}

//-----------------------------------------------------------------------------
void qMRMLMarkupsGridSurfaceSettingsWidget::onGridSurfaceTypeParameterChanged()
{
  Q_D(qMRMLMarkupsGridSurfaceSettingsWidget);
  vtkMRMLMarkupsGridSurfaceNode* gridSurfaceNode = vtkMRMLMarkupsGridSurfaceNode::SafeDownCast(this->MarkupsNode);
  if (!gridSurfaceNode)
  {
    return;
  }
  MRMLNodeModifyBlocker blocker(gridSurfaceNode);
  gridSurfaceNode->SetGridSurfaceType(d->surfaceTypeComboBox->currentData().toInt());
}

//-----------------------------------------------------------------------------
void qMRMLMarkupsGridSurfaceSettingsWidget::onApplyGridResolution()
{
  Q_D(qMRMLMarkupsGridSurfaceSettingsWidget);
  vtkMRMLMarkupsGridSurfaceNode* gridSurfaceNode = vtkMRMLMarkupsGridSurfaceNode::SafeDownCast(this->MarkupsNode);
  if (!gridSurfaceNode)
  {
    return;
  }
  MRMLNodeModifyBlocker blocker(gridSurfaceNode);
  gridSurfaceNode->SetGridResolution(d->gridResolutionSpinBox_X->value(), d->gridResolutionSpinBox_Y->value());
}
