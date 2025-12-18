#ifndef EVENT_HH
#define EVENT_HH

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4UserEventAction.hh"
#include "simRun.hh"

class simEventAction : public G4UserEventAction
{
  public:
    simEventAction(simRunAction *);
    ~simEventAction();

    virtual void BeginOfEventAction(const G4Event *);
    virtual void EndOfEventAction(const G4Event *);

    void AddEdep(G4double edep)
    {
        fEdep += edep;
    } // transfer energy from stepping class to it//{what the function should
      // do}

  private:
    //  simRunAction *fRunAction;
    G4double fEdep; // deposite of energy
};

#endif
