#include "simEvent.hh"
#include "G4RunManager.hh"
#include "simRun.hh"

simEventAction::simEventAction(simRunAction *)
{
    fEdep = 0.; // start value we should define
}

simEventAction::~simEventAction() {}

void simEventAction::BeginOfEventAction(const G4Event *event)
{
    fEdep = 0.; // reset it when event start

    G4cout << "Event number: " << event->GetEventID() << G4endl;
}

void simEventAction::EndOfEventAction(const G4Event *event)
{
    // fRunAction->AddEdep(fEdep);

    // G4cout << "Energy deposition: " << fEdep << G4endl;

    // root file
    G4AnalysisManager *man = G4AnalysisManager::Instance();

    man->FillNtupleDColumn(2, 0, fEdep);

    man->AddNtupleRow(2);
}
