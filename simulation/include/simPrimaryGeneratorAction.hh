#ifndef simPrimaryGeneratorAction_h
#define simPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4AnalysisManager.hh"
#include "G4ChargedGeantino.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "globals.hh"

// The primary generator action class with particle gum.
/// It defines an ion (F18), at rest, randomly distribued within a zone

class simPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    simPrimaryGeneratorAction();
    ~simPrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event *);

  private:
    G4ParticleGun *fParticleGun;
};

#endif
