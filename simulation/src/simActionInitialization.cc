#include "simActionInitialization.hh"
#include "simPrimaryGeneratorAction.hh"

#include "simEvent.hh"
#include "simRun.hh"
#include "simStepping.hh"

simActionInitialization::simActionInitialization()
    : G4VUserActionInitialization()
{
}

simActionInitialization::~simActionInitialization() {}

void simActionInitialization::BuildForMaster() const
{
    simRunAction *runAction = new simRunAction();
    SetUserAction(runAction);
}

void simActionInitialization::Build() const
{
    simPrimaryGeneratorAction *generator = new simPrimaryGeneratorAction();
    SetUserAction(generator);

    simRunAction *runAction = new simRunAction();
    SetUserAction(runAction);

    simEventAction *eventAction = new simEventAction(runAction);
    SetUserAction(eventAction);

    simSteppingAction *steppingAction = new simSteppingAction(eventAction);
    SetUserAction(steppingAction);
}
