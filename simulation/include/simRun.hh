#ifndef simRun_HH
#define simRun_HH

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4UserRunAction.hh"

class simRunAction : public G4UserRunAction
{
  public:
    simRunAction();
    ~simRunAction();

    virtual void BeginOfRunAction(const G4Run *);
    virtual void EndOfRunAction(const G4Run *);
};

#endif
