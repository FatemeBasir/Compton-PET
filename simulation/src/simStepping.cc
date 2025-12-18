#include "simStepping.hh"

simSteppingAction::simSteppingAction(simEventAction *eventAction)
{
    fEventAction = eventAction; // get access to the object we created
}

simSteppingAction::~simSteppingAction() {}

void simSteppingAction::UserSteppingAction(const G4Step *step)
{
    G4LogicalVolume *volume = step->GetPreStepPoint()
                                  ->GetTouchableHandle()
                                  ->GetVolume()
                                  ->GetLogicalVolume();

    // check scoring volume
    const simDetectorConstruction *detectorConstruction =
        static_cast<const simDetectorConstruction *>(
            G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    G4LogicalVolume *fScoringVolume = detectorConstruction->GetScoringVolume();

    G4int eventID =
        G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    G4int runID = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
    G4int numberOfEvents = G4RunManager::GetRunManager()
                               ->GetCurrentRun()
                               ->GetNumberOfEventToBeProcessed();

    G4long uniqueEventID =
        runID * numberOfEvents + eventID; // Caluculate an uniqe event ID

    // G4cout << volume->GetName() << G4endl;

    G4Track* track = step->GetTrack();
        
    if (track->GetDefinition() == G4Gamma::GammaDefinition())
    {
        // Ensure this is the first step
        if (track->GetCurrentStepNumber() == 1)
        {
            // Get the direction of the photon
            G4ThreeVector photonDirection = track->GetMomentumDirection();

            // Calculate the angle between the photon and the primary particle
            G4double theta = photonDirection.getTheta();
            G4double phi = photonDirection.getPhi();
            
            G4double thetaDeg = theta * 180.0 / CLHEP::pi;
            G4double phiDeg = phi * 180.0 / CLHEP::pi;

            G4AnalysisManager *man = G4AnalysisManager::Instance();

            man->FillNtupleIColumn(4, 0, uniqueEventID);
            man->FillNtupleDColumn(4, 1, thetaDeg);
            man->FillNtupleDColumn(4, 2, phiDeg);
            man->AddNtupleRow(4);
        }
    }

    if (volume != fScoringVolume)
        return;

    G4double edep = step->GetTotalEnergyDeposit();
    fEventAction->AddEdep(edep); // then accumulate
}
