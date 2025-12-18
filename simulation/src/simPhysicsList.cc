#include "simPhysicsList.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "QGSP_BERT.hh"

simPhysicsList::simPhysicsList()
{
    // SetVerboseLevel(1);

    // Default physics->any particle decay
    RegisterPhysics(new G4DecayPhysics());

    // EM physics
    RegisterPhysics(new G4EmStandardPhysics());

    // Optical photons
    RegisterPhysics(new G4OpticalPhysics());

    // Radioactive decay->decays of R-ion
    RegisterPhysics(new G4RadioactiveDecayPhysics());
    
    // Handling radioactive decays
    RegisterPhysics(new G4DecayPhysics());
    
    // Extra physics processes
    RegisterPhysics(new G4EmExtraPhysics());
};

simPhysicsList::~simPhysicsList(){};
