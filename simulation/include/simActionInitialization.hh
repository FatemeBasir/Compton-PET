#ifndef simActionInitialization_h
#define simActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "simPrimaryGeneratorAction.hh"

#include "simEvent.hh"
#include "simRun.hh"
#include "simStepping.hh"

class simActionInitialization : public G4VUserActionInitialization
{
  public:
    simActionInitialization();
    virtual ~simActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;
};

#endif
