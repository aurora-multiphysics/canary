#include "CanaryApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
CanaryApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

CanaryApp::CanaryApp(const InputParameters & parameters) : MooseApp(parameters)
{
  CanaryApp::registerAll(_factory, _action_factory, _syntax);
}

CanaryApp::~CanaryApp() {}

void
CanaryApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<CanaryApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"CanaryApp"});
  Registry::registerActionsTo(af, {"CanaryApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
CanaryApp::registerApps()
{
  registerApp(CanaryApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
CanaryApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  CanaryApp::registerAll(f, af, s);
}
extern "C" void
CanaryApp__registerApps()
{
  CanaryApp::registerApps();
}
