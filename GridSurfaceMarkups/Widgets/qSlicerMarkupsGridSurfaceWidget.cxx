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

// Markups widgets includes
#include "qSlicerMarkupsAdditionalOptionsWidget_p.h"

// qMRML includes
#include "qSlicerMarkupsGridSurfaceWidget.h"
#include "ui_qSlicerMarkupsGridSurfaceWidget.h"

// MRML includes
#include <vtkMRMLDisplayNode.h>

// GridSurfaceMarkups includes
#include <vtkMRMLMarkupsGridSurfaceNode.h>

// --------------------------------------------------------------------------
class qSlicerMarkupsGridSurfaceWidgetPrivate:
  public qSlicerMarkupsAdditionalOptionsWidgetPrivate,
  public Ui_qSlicerMarkupsGridSurfaceWidget
{
  Q_DECLARE_PUBLIC(qSlicerMarkupsGridSurfaceWidget);

protected:
  qSlicerMarkupsGridSurfaceWidget* const q_ptr;

public:
  qSlicerMarkupsGridSurfaceWidgetPrivate(qSlicerMarkupsGridSurfaceWidget* object);
  void setupUi(qSlicerMarkupsGridSurfaceWidget* widget);

public:
  bool IsProcessingOnMRMLNodeModified;
};

// --------------------------------------------------------------------------
qSlicerMarkupsGridSurfaceWidgetPrivate::qSlicerMarkupsGridSurfaceWidgetPrivate(qSlicerMarkupsGridSurfaceWidget* object)
  : q_ptr(object)
{
  this->MarkupsNode = nullptr;
  this->IsProcessingOnMRMLNodeModified = false;
}

// --------------------------------------------------------------------------
void qSlicerMarkupsGridSurfaceWidgetPrivate::setupUi(qSlicerMarkupsGridSurfaceWidget* widget)
{
  Q_Q(qSlicerMarkupsGridSurfaceWidget);

  this->Ui_qSlicerMarkupsGridSurfaceWidget::setupUi(widget);

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
// qSlicerMarkupsGridSurfaceWidget methods

// --------------------------------------------------------------------------
qSlicerMarkupsGridSurfaceWidget::
qSlicerMarkupsGridSurfaceWidget(QWidget* parent)
  : Superclass(*new qSlicerMarkupsGridSurfaceWidgetPrivate(this), parent)
{
  this->setup();
}

// --------------------------------------------------------------------------
qSlicerMarkupsGridSurfaceWidget::
qSlicerMarkupsGridSurfaceWidget(qSlicerMarkupsGridSurfaceWidgetPrivate &d, QWidget* parent)
  : Superclass(d, parent)
{
  this->setup();
}

// --------------------------------------------------------------------------
qSlicerMarkupsGridSurfaceWidget::~qSlicerMarkupsGridSurfaceWidget() = default;

// --------------------------------------------------------------------------
void qSlicerMarkupsGridSurfaceWidget::setup()
{
  Q_D(qSlicerMarkupsGridSurfaceWidget);
  d->setupUi(this);
}

// --------------------------------------------------------------------------
vtkMRMLMarkupsGridSurfaceNode* qSlicerMarkupsGridSurfaceWidget::mrmlGridSurfaceNode()const
{
  Q_D(const qSlicerMarkupsGridSurfaceWidget);
  return vtkMRMLMarkupsGridSurfaceNode::SafeDownCast(d->MarkupsNode);
}

// --------------------------------------------------------------------------
void qSlicerMarkupsGridSurfaceWidget::setMRMLMarkupsNode(vtkMRMLMarkupsNode* markupsNode)
{
  Q_D(qSlicerMarkupsGridSurfaceWidget);

  Superclass::setMRMLMarkupsNode(markupsNode);

  this->qvtkReconnect(d->MarkupsNode, markupsNode, vtkCommand::ModifiedEvent,
                      this, SLOT(onMRMLNodeModified()));

  this->onMRMLNodeModified();
  this->setEnabled(markupsNode != nullptr);
}

// --------------------------------------------------------------------------
void qSlicerMarkupsGridSurfaceWidget::setMRMLMarkupsNode(vtkMRMLNode* node)
{
  this->setMRMLMarkupsNode(vtkMRMLMarkupsGridSurfaceNode::SafeDownCast(node));
}

// --------------------------------------------------------------------------
void qSlicerMarkupsGridSurfaceWidget::onMRMLNodeModified()
{
  Q_D(qSlicerMarkupsGridSurfaceWidget);

  vtkMRMLMarkupsGridSurfaceNode* gridSurfaceNode = vtkMRMLMarkupsGridSurfaceNode::SafeDownCast(d->MarkupsNode);
  if (!gridSurfaceNode)
  {
    return;
  }

  d->IsProcessingOnMRMLNodeModified = true;

  int gridResolution[2] = {0};
  gridSurfaceNode->GetGridResolution(gridResolution);

  d->gridResolutionSpinBox_A->setValue(gridResolution[0]);
  d->gridResolutionSpinBox_B->setValue(gridResolution[1]);

  d->IsProcessingOnMRMLNodeModified = false;
}

//-----------------------------------------------------------------------------
bool qSlicerMarkupsGridSurfaceWidget::canManageMRMLMarkupsNode(vtkMRMLMarkupsNode* markupsNode) const
{
  Q_D(const qSlicerMarkupsGridSurfaceWidget);

  vtkMRMLMarkupsGridSurfaceNode* gridSurfaceNode = vtkMRMLMarkupsGridSurfaceNode::SafeDownCast(markupsNode);
  if (!gridSurfaceNode)
  {
    return false;
  }

  return true;
}

// --------------------------------------------------------------------------
void qSlicerMarkupsGridSurfaceWidget::updateWidgetFromMRML()
{
  Q_D(qSlicerMarkupsGridSurfaceWidget);

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
void qSlicerMarkupsGridSurfaceWidget::onGridSurfaceTypeParameterChanged()
{
  Q_D(qSlicerMarkupsGridSurfaceWidget);
  vtkMRMLMarkupsGridSurfaceNode* gridSurfaceNode = vtkMRMLMarkupsGridSurfaceNode::SafeDownCast(d->MarkupsNode);
  if (!gridSurfaceNode)
  {
    return;
  }
  MRMLNodeModifyBlocker blocker(gridSurfaceNode);
  gridSurfaceNode->SetGridSurfaceType(d->surfaceTypeComboBox->currentData().toInt());
}

//-----------------------------------------------------------------------------
void qSlicerMarkupsGridSurfaceWidget::applyGridResolutionButton()
{
  Q_D(qSlicerMarkupsGridSurfaceWidget);
  vtkMRMLMarkupsGridSurfaceNode* gridSurfaceNode = vtkMRMLMarkupsGridSurfaceNode::SafeDownCast(d->MarkupsNode);
  if (!gridSurfaceNode)
  {
    return;
  }
  MRMLNodeModifyBlocker blocker(gridSurfaceNode);
  gridSurfaceNode->SetGridResolution(d->gridResolutionSpinBox_A->value(), d->gridResolutionSpinBox_B->value());
}
