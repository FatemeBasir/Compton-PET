#include "simRun.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

simRunAction::simRunAction()
{
    ////////////////////      root File     //////////////////
    G4AnalysisManager *man = G4AnalysisManager::Instance();

    ////position of detector:  D-column

    man->CreateNtuple("Photons", "Photons");
    man->CreateNtupleIColumn("fEvent");       // I : integer
    man->CreateNtupleDColumn("fX");           // D : double
    man->CreateNtupleDColumn("fY");
    man->CreateNtupleDColumn("fZ");
    man->CreateNtupleDColumn("fTime");
    man->CreateNtupleDColumn("fWlen");        // wavelength of photon
    man->FinishNtuple(0);                     // 0: Ntuple number zero that we created

    man->CreateNtuple("Hits", "Hits");
    man->CreateNtupleIColumn("fEvent");
    man->CreateNtupleIColumn("fPDG");
    man->CreateNtupleDColumn("fX");
    man->CreateNtupleDColumn("fY");
    man->CreateNtupleDColumn("fZ");
    man->CreateNtupleDColumn("fTime");
    man->CreateNtupleDColumn("fEnergy");
    man->FinishNtuple(1);

    man->CreateNtuple("Scoring", "Scoring");
    man->CreateNtupleDColumn("fEdep");
    man->FinishNtuple(2);

    man->CreateNtuple("Header", "Header");
    man->CreateNtupleIColumn("mcEvent");
    man->CreateNtupleDColumn("mcX");
    man->CreateNtupleDColumn("mcY");
    man->CreateNtupleDColumn("mcZ");
    man->FinishNtuple(3);

    man->CreateNtuple("MCTruth", "MCTruth");
    man->CreateNtupleIColumn("fEvent");
    man->CreateNtupleDColumn("fThetaDeg");
    man->CreateNtupleDColumn("fPhiDeg");
    man->FinishNtuple(4);
}

simRunAction::~simRunAction() {}

void simRunAction::BeginOfRunAction(const G4Run *run)
{
    ////////////////////      root File     //////////////////
    G4AnalysisManager *man = G4AnalysisManager::Instance();

    G4int runID = run->GetRunID();

    std::stringstream strRunID; // integer -> string
    strRunID << runID;

    man->OpenFile("output" + strRunID.str() + ".root");
    
    G4cout << "Starting run " << runID << G4endl;
}

void simRunAction::EndOfRunAction(const G4Run * run)
{
    G4AnalysisManager *man = G4AnalysisManager::Instance();

    man->Write(); // before closeing rootFile we should add write cause of
                  // recovery &...
    man->CloseFile();
    
    G4int runID = run->GetRunID();
    
    G4cout << "Finishing run " << runID << G4endl;
}
