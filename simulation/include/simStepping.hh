#ifndef STEPPING_HH
#define STEPPING_HH

#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Step.hh"
#include "G4UserSteppingAction.hh"

#include "simDetectorConstruction.hh"
#include "simEvent.hh"

class simSteppingAction : public G4UserSteppingAction
{
  public:
    simSteppingAction(simEventAction *eventAction);
    ~simSteppingAction();

    virtual void UserSteppingAction(const G4Step *);

  private:
    simEventAction *fEventAction;
};

#endif
