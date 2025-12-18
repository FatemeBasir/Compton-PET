#include "simDetector.hh"

#include "G4PhysicsOrderedFreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "simActionInitialization.hh"

simSensitiveDetector::simSensitiveDetector(G4String name)
    : G4VSensitiveDetector(name) // name of the detector
{
    /*quEff = new G4PhysicsOrderedFreeVector();

    std::ifstream datafile;
    datafile.open("eff.dat");

    while (1)
    {
        G4double wlen, queff;

        datafile >> wlen >> queff;

        if (datafile.eof())
            break;

        // G4cout << wlen << " " << queff << std::endl;

        quEff->InsertValues(wlen, queff / 100.); // 100: values between 0 - 1
    }

    datafile.close();

    quEff->GetSpline(); // it has an Erorr*/
}

simSensitiveDetector::~simSensitiveDetector() {}

G4double simSensitiveDetector::SmearEnergy(G4double energy)
{
    // Define your detector resolution (e.g., 2% at 511 keV)
    G4double resolution = 0.05;
    G4double sigma = resolution * energy;
    
    return G4RandGauss::shoot(energy, sigma);
}

G4bool simSensitiveDetector::ProcessHits(G4Step *aStep,
                                         G4TouchableHistory *ROhist)
{
    G4Track *track = aStep->GetTrack();

    // track->SetTrackStatus(fStopAndKill);

    G4StepPoint *preStepPoint =
        aStep->GetPreStepPoint(); // prePoint=when photons enter to the
                                  // sensitive detector
    G4StepPoint *postStepPoint =
        aStep->GetPostStepPoint(); // posPoint=when photons leave from the
                                   // sensitivedetector

    //--------------------------------------position------------------------------------//
    G4ThreeVector posPhoton = preStepPoint->GetPosition();
    G4ThreeVector momPhoton =
        preStepPoint->GetMomentum(); // to calculate wavelength
    G4double momPhotonMag = momPhoton.mag();

    G4double globalTime = preStepPoint->GetGlobalTime();

    G4double wlen = (1.239841939 * eV / momPhotonMag) * 1E+03; // mag->magnitude

    track->SetTrackStatus(fStopAndKill);

    const G4ParticleDefinition *pd = track->GetParticleDefinition();
    G4int pdg = pd->GetPDGEncoding();

    // G4cout << "Photon position: " << posPhoton << G4endl;

    //-----------------------------------copy-number-------------------------------------//
    // reconstruction of photon & info
    const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();
    G4int copyNo = touchable->GetCopyNumber();

    // G4cout << "copy number: " << copyNo << G4endl;

    //-----------------------position of detector in world
    // volume------------------------//
    G4VPhysicalVolume *physVol = touchable->GetVolume();
    G4ThreeVector posDetector = physVol->GetTranslation();
#ifndef G4MULTITHREADED // in order to be fast
    G4cout << "Detector position: " << posDetector << G4endl;
    G4cout << "Photon wavelength: " << posDetector << G4endl;
#endif

    G4int eventID =
        G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    G4int runID = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
    G4int numberOfEvents = G4RunManager::GetRunManager()
                               ->GetCurrentRun()
                               ->GetNumberOfEventToBeProcessed();

    G4long uniqueEventID =
        runID * numberOfEvents + eventID; // Caluculate an uniqe event ID

    G4AnalysisManager *man = G4AnalysisManager::Instance();

    globalTime += uniqueEventID * 100 * ns; // Adding event time

    man->FillNtupleIColumn(0, 0, uniqueEventID); // Event number
    man->FillNtupleDColumn(0, 1, posPhoton[0]);
    man->FillNtupleDColumn(0, 2, posPhoton[1]);
    man->FillNtupleDColumn(0, 3, posPhoton[2]);
    man->FillNtupleDColumn(0, 4, globalTime);
    man->FillNtupleDColumn(0, 5, wlen);
    man->AddNtupleRow(0);

    // if(G4UniformRand() < quEff->Value(wlen))
    {
        man->FillNtupleIColumn(1, 0, uniqueEventID);
        man->FillNtupleIColumn(1, 1, pdg);
        man->FillNtupleDColumn(1, 2, posDetector[0]);
        man->FillNtupleDColumn(1, 3, posDetector[1]);
        man->FillNtupleDColumn(1, 4, posDetector[2]);
        man->FillNtupleDColumn(1, 5, globalTime);
        man->FillNtupleDColumn(1, 6, SmearEnergy(momPhotonMag));
        man->AddNtupleRow(1);
    }

    return true;
};
