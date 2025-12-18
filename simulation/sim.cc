#include <iostream> //c++

#include "G4MTRunManager.hh"
#include "G4RunManager.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4VisManager.hh"

#include "simActionInitialization.hh"
#include "simDetectorConstruction.hh"
#include "simPhysicsList.hh"

int main(int argc, char **argv)
{

    G4UIExecutive *ui = 0;

#ifdef G4MULTITHREADED
    G4MTRunManager *runManager = new G4MTRunManager;
#else
    G4RunManager *runManager = new G4RunManager;
#endif

    // Detector construction
    runManager->SetUserInitialization(new simDetectorConstruction());

    // Physics list
    runManager->SetUserInitialization(new simPhysicsList());

    // User action initialization
    runManager->SetUserInitialization(new simActionInitialization());

    if (argc == 1)
    {
        ui = new G4UIExecutive(argc, argv);
    }

    // Initialize visualization
    G4VisManager *visManager = new G4VisExecutive();
    visManager->Initialize();

    // Get the pointer to the User Interface manager
    G4UImanager *UImanager = G4UImanager::GetUIpointer();

    if (ui)
    {
        UImanager->ApplyCommand("/control/execute vis.mac");
        ui->SessionStart();
    }
    else
    {
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    }

    return 0;
}
