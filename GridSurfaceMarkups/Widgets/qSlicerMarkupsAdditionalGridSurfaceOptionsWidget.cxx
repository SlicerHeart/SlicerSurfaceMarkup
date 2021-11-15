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
#include "qSlicerMarkupsAdditionalGridSurfaceOptionsWidget.h"
#include "ui_qSlicerMarkupsAdditionalGridSurfaceOptionsWidget.h"

// MRML includes
#include <vtkMRMLDisplayNode.h>

// GridSurfaceMarkups includes
#include <vtkMRMLMarkupsGridSurfaceNode.h>

// --------------------------------------------------------------------------
class qSlicerMarkupsAdditionalGridSurfaceOptionsWidgetPrivate:
  public qSlicerMarkupsAdditionalOptionsWidgetPrivate,
  public Ui_qSlicerMarkupsAdditionalGridSurfaceOptionsWidget
{
  Q_DECLARE_PUBLIC(qSlicerMarkupsAdditionalGridSurfaceOptionsWidget);

protected:
  qSlicerMarkupsAdditionalGridSurfaceOptionsWidget* const q_ptr;

public:
  qSlicerMarkupsAdditionalGridSurfaceOptionsWidgetPrivate(qSlicerMarkupsAdditionalGridSurfaceOptionsWidget* object);
  void setupUi(qSlicerMarkupsAdditionalGridSurfaceOptionsWidget* widget);

public:
  bool IsProcessingOnMRMLNodeModified;
};

// --------------------------------------------------------------------------
qSlicerMarkupsAdditionalGridSurfaceOptionsWidgetPrivate::qSlicerMarkupsAdditionalGridSurfaceOptionsWidgetPrivate(qSlicerMarkupsAdditionalGridSurfaceOptionsWidget* object)
  : q_ptr(object)
{
  this->MarkupsNode = nullptr;
  this->IsProcessingOnMRMLNodeModified = false;
}

// --------------------------------------------------------------------------
void qSlicerMarkupsAdditionalGridSurfaceOptionsWidgetPrivate::setupUi(qSlicerMarkupsAdditionalGridSurfaceOptionsWidget* widget)
{
  Q_Q(qSlicerMarkupsAdditionalGridSurfaceOptionsWidget);

  this->Ui_qSlicerMarkupsAdditionalGridSurfaceOptionsWidget::setupUi(widget);

  this->gridSurfaceSettingsCollapsibleButton->setVisible(false);
  this->surfaceTypeComboBox->clear();
  for (int gridSurfaceType = 0; gridSurfaceType < vtkMRMLMarkupsGridSurfaceNode::GridSurfaceType_Last; ++gridSurfaceType)
    {
    this->surfaceTypeComboBox->addItem(vtkMRMLMarkupsGridSurfaceNode::GetGridSurfaceTypeAsString(gridSurfaceType), gridSurfaceType);
    }

  QObject::connect(this->surfaceTypeComboBox, SIGNAL(currentIndexChanged(int)),
                   q, SLOT(onGridSurfaceTypeParameterChanged()));
  QObject::connect(this->applyGridResolutionButton, SIGNAL(clicked()),
                   q, SLOT(onApplyGridResolution()));
  q->setEnabled(this->MarkupsNode != nullptr);
}

// --------------------------------------------------------------------------
// qSlicerMarkupsAdditionalGridSurfaceOptionsWidget methods

// --------------------------------------------------------------------------
qSlicerMarkupsAdditionalGridSurfaceOptionsWidget::
qSlicerMarkupsAdditionalGridSurfaceOptionsWidget(QWidget* parent)
  : Superclass(*new qSlicerMarkupsAdditionalGridSurfaceOptionsWidgetPrivate(this), parent)
{
  this->setup();
}

// --------------------------------------------------------------------------
qSlicerMarkupsAdditionalGridSurfaceOptionsWidget::
qSlicerMarkupsAdditionalGridSurfaceOptionsWidget(qSlicerMarkupsAdditionalGridSurfaceOptionsWidgetPrivate &d, QWidget* parent)
  : Superclass(d, parent)
{
  this->setup();
}

// --------------------------------------------------------------------------
qSlicerMarkupsAdditionalGridSurfaceOptionsWidget::~qSlicerMarkupsAdditionalGridSurfaceOptionsWidget() = default;

// --------------------------------------------------------------------------
void qSlicerMarkupsAdditionalGridSurfaceOptionsWidget::setup()
{
  Q_D(qSlicerMarkupsAdditionalGridSurfaceOptionsWidget);
  d->setupUi(this);
}

// --------------------------------------------------------------------------
vtkMRMLMarkupsGridSurfaceNode* qSlicerMarkupsAdditionalGridSurfaceOptionsWidget::mrmlGridSurfaceNode()const
{
  Q_D(const qSlicerMarkupsAdditionalGridSurfaceOptionsWidget);
  return vtkMRMLMarkupsGridSurfaceNode::SafeDownCast(d->MarkupsNode);
}

// --------------------------------------------------------------------------
void qSlicerMarkupsAdditionalGridSurfaceOptionsWidget::setMRMLMarkupsNode(vtkMRMLMarkupsNode* markupsNode)
{
  Q_D(qSlicerMarkupsAdditionalGridSurfaceOptionsWidget);

  Superclass::setMRMLMarkupsNode(markupsNode);

  this->qvtkReconnect(d->MarkupsNode, markupsNode, vtkCommand::ModifiedEvent,
                      this, SLOT(onMRMLNodeModified()));

  this->onMRMLNodeModified();
  this->setEnabled(markupsNode != nullptr);
}

// --------------------------------------------------------------------------
void qSlicerMarkupsAdditionalGridSurfaceOptionsWidget::setMRMLMarkupsNode(vtkMRMLNode* node)
{
  this->setMRMLMarkupsNode(vtkMRMLMarkupsGridSurfaceNode::SafeDownCast(node));
}

// --------------------------------------------------------------------------
void qSlicerMarkupsAdditionalGridSurfaceOptionsWidget::onMRMLNodeModified()
{
  Q_D(qSlicerMarkupsAdditionalGridSurfaceOptionsWidget);

  vtkMRMLMarkupsGridSurfaceNode* gridSurfaceNode = vtkMRMLMarkupsGridSurfaceNode::SafeDownCast(d->MarkupsNode);
  if (!gridSurfaceNode)
  {
    return;
  }

  d->IsProcessingOnMRMLNodeModified = true;

  int gridResolution[2] = {0};
  gridSurfaceNode->GetGridResolution(gridResolution);

  d->gridResolutionSpinBox_X->setValue(gridResolution[0]);
  d->gridResolutionSpinBox_Y->setValue(gridResolution[1]);

  d->IsProcessingOnMRMLNodeModified = false;
}

//-----------------------------------------------------------------------------
bool qSlicerMarkupsAdditionalGridSurfaceOptionsWidget::canManageMRMLMarkupsNode(vtkMRMLMarkupsNode* markupsNode) const
{
  Q_D(const qSlicerMarkupsAdditionalGridSurfaceOptionsWidget);

  vtkMRMLMarkupsGridSurfaceNode* gridSurfaceNode = vtkMRMLMarkupsGridSurfaceNode::SafeDownCast(markupsNode);
  if (!gridSurfaceNode)
  {
    return false;
  }

  return true;
}

// --------------------------------------------------------------------------
void qSlicerMarkupsAdditionalGridSurfaceOptionsWidget::updateWidgetFromMRML()
{
  Q_D(qSlicerMarkupsAdditionalGridSurfaceOptionsWidget);

  if (!this->canManageMRMLMarkupsNode(d->MarkupsNode))
  {
    d->gridSurfaceSettingsCollapsibleButton->setVisible(false);
    return;
  }

  d->gridSurfaceSettingsCollapsibleButton->setVisible(true);
  vtkMRMLMarkupsGridSurfaceNode* markupsGridSurfaceNode = vtkMRMLMarkupsGridSurfaceNode::SafeDownCast(d->MarkupsNode);
  if (markupsGridSurfaceNode)
  {
    bool wasBlocked = d->surfaceTypeComboBox->blockSignals(true);
    d->surfaceTypeComboBox->setCurrentIndex(d->surfaceTypeComboBox->findData(markupsGridSurfaceNode->GetGridSurfaceType()));
    d->surfaceTypeComboBox->blockSignals(wasBlocked);
    this->setMRMLMarkupsNode(markupsGridSurfaceNode);
  }
}

//-----------------------------------------------------------------------------
void qSlicerMarkupsAdditionalGridSurfaceOptionsWidget::onGridSurfaceTypeParameterChanged()
{
  Q_D(qSlicerMarkupsAdditionalGridSurfaceOptionsWidget);
  vtkMRMLMarkupsGridSurfaceNode* gridSurfaceNode = vtkMRMLMarkupsGridSurfaceNode::SafeDownCast(d->MarkupsNode);
  if (!gridSurfaceNode)
  {
    return;
  }
  MRMLNodeModifyBlocker blocker(gridSurfaceNode);
  gridSurfaceNode->SetGridSurfaceType(d->surfaceTypeComboBox->currentData().toInt());
}

//-----------------------------------------------------------------------------
void qSlicerMarkupsAdditionalGridSurfaceOptionsWidget::onApplyGridResolution()
{
  Q_D(qSlicerMarkupsAdditionalGridSurfaceOptionsWidget);
  vtkMRMLMarkupsGridSurfaceNode* gridSurfaceNode = vtkMRMLMarkupsGridSurfaceNode::SafeDownCast(d->MarkupsNode);
  if (!gridSurfaceNode)
  {
    return;
  }
  MRMLNodeModifyBlocker blocker(gridSurfaceNode);
  gridSurfaceNode->SetGridResolution(d->gridResolutionSpinBox_X->value(), d->gridResolutionSpinBox_Y->value());
}
