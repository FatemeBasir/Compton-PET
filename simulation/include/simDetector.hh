#ifndef simDetector_HH
#define simDetector_HH

#include "G4AnalysisManager.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VSensitiveDetector.hh"

class simSensitiveDetector : public G4VSensitiveDetector
{
  public:
    simSensitiveDetector(G4String);
    ~simSensitiveDetector();

  private:
    virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);

    G4PhysicsOrderedFreeVector *quEff;
    
    G4double SmearEnergy(G4double energy);
};

#endif
