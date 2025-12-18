#include "simPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"

simPrimaryGeneratorAction::simPrimaryGeneratorAction()
{
    double phi = 2. * M_PI * G4UniformRand();    // Uniformly from 0 to 2π
    double cosTheta = 2. * G4UniformRand() - 1.; // Uniformly from -1 to 1
    double sinTheta = std::sqrt(1. - cosTheta * cosTheta);

    G4double px = sinTheta * std::cos(phi);
    G4double py = sinTheta * std::sin(phi);
    G4double pz = cosTheta;

    fParticleGun = new G4ParticleGun(
        1); // in (): number of particles per event, 1 primary vertex per event

    // Particle type:
    G4ParticleTable *particleTable = G4ParticleTable ::GetParticleTable();
    G4ParticleDefinition *particle = particleTable->FindParticle("e+");

    // Particle position
    G4ThreeVector pos(-19.4395, -2.9528, 0.);

    // Particle direction
    G4ThreeVector mom(px, py, pz);

    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(mom);
    fParticleGun->SetParticleEnergy(0. * keV);
    fParticleGun->SetParticleDefinition(particle);
}

simPrimaryGeneratorAction::~simPrimaryGeneratorAction() { delete fParticleGun; }

void simPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
    G4AnalysisManager *man = G4AnalysisManager::Instance();

    G4ParticleDefinition *particle = fParticleGun->GetParticleDefinition();

    if (particle == G4ChargedGeantino::ChargedGeantino())
    {
        G4int Z = 9; // fluorine
        G4int A = 18;

        G4double charge = 0. * eplus; // atom
        G4double energy = 0. * keV;

        G4ParticleDefinition *ion = G4IonTable::GetIonTable()->GetIon(
            Z, A, energy); // related to PhysicsList
        fParticleGun->SetParticleDefinition(ion);
        fParticleGun->SetParticleCharge(charge);

        // G4ThreeVector pos(G4UniformRand() * 20,G4UniformRand() * 20,-20);
    }

    G4int eventID = anEvent->GetEventID();
    G4int runID = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
    G4int numberOfEvents = G4RunManager::GetRunManager()
                               ->GetCurrentRun()
                               ->GetNumberOfEventToBeProcessed();
    G4long uniqueEventID =
        runID * numberOfEvents + eventID; // Caluculate an uniqe event ID

    G4ThreeVector pos = fParticleGun->GetParticlePosition();

    man->FillNtupleIColumn(3, 0, uniqueEventID); // event number
    man->FillNtupleDColumn(3, 1, pos[0]);        // event number
    man->FillNtupleDColumn(3, 2, pos[1]);        // event number
    man->FillNtupleDColumn(3, 3, pos[2]);        // event number
    man->AddNtupleRow(3);

    // create vertex
    fParticleGun->GeneratePrimaryVertex(anEvent);
}
