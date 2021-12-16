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

#include "qSlicerGridSurfaceMarkupsModule.h"

// Qt includes
#include <QDebug>

// GridSurface Logic includes
#include "vtkSlicerGridSurfaceMarkupsLogic.h"

// GridSurface Widgets
#include "qMRMLMarkupsGridSurfaceSettingsWidget.h"

// Markups Widgets includes
#include "qMRMLMarkupsOptionsWidgetsFactory.h" 

#include <qSlicerModuleManager.h>
#include <qSlicerCoreApplication.h>

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerGridSurfaceMarkupsModulePrivate
{
public:
  qSlicerGridSurfaceMarkupsModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerGridSurfaceMarkupsModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerGridSurfaceMarkupsModulePrivate::qSlicerGridSurfaceMarkupsModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerGridSurfaceMarkupsModule methods

//-----------------------------------------------------------------------------
qSlicerGridSurfaceMarkupsModule::qSlicerGridSurfaceMarkupsModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerGridSurfaceMarkupsModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerGridSurfaceMarkupsModule::~qSlicerGridSurfaceMarkupsModule()
{
}

//-----------------------------------------------------------------------------
bool qSlicerGridSurfaceMarkupsModule::isHidden() const
{
  // The module has no GUI.
  // Widget options will be shown in Markups module.
  return true;
}

//-----------------------------------------------------------------------------
QString qSlicerGridSurfaceMarkupsModule::helpText() const
{
  return "";
}

//-----------------------------------------------------------------------------
QString qSlicerGridSurfaceMarkupsModule::acknowledgementText() const
{
  return "";
}

//-----------------------------------------------------------------------------
QStringList qSlicerGridSurfaceMarkupsModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString("Csaba Pinter (Pixel Medical / Ebatinca)");
  moduleContributors << QString("Andras Lasso (PerkLab, Queen's)");
  moduleContributors << QString("Rafael Palomar (Oslo University Hospital / NTNU)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerGridSurfaceMarkupsModule::icon() const
{
  return QIcon(":/Icons/GridSurfaceMarkups.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerGridSurfaceMarkupsModule::categories() const
{
  return QStringList() << "Informatics";
}

//-----------------------------------------------------------------------------
QStringList qSlicerGridSurfaceMarkupsModule::dependencies() const
{
  return QStringList() << "Markups";
}

//-----------------------------------------------------------------------------
void qSlicerGridSurfaceMarkupsModule::setup()
{
  this->Superclass::setup();

  // Create and configure the options widgets
  auto optionsWidgetFactory = qMRMLMarkupsOptionsWidgetsFactory::instance();
  optionsWidgetFactory->registerOptionsWidget(new qMRMLMarkupsGridSurfaceSettingsWidget()); 
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerGridSurfaceMarkupsModule::createWidgetRepresentation()
{
  return nullptr;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerGridSurfaceMarkupsModule::createLogic()
{
  return vtkSlicerGridSurfaceMarkupsLogic::New();
}
