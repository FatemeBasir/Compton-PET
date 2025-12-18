#ifndef simDetectorConstruction_h
#define simDetectorConstruction_h 1

#include "simDetector.hh"

#include "G4Box.hh"
#include "G4GenericMessenger.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

class simDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    simDetectorConstruction();
    virtual ~simDetectorConstruction();

    G4LogicalVolume *GetScoringVolume() const { return fScoringVolume; }

    virtual G4VPhysicalVolume *Construct();
    void ConstructScintillator();
    // void ConstructTOF();

  private:
    // with below commands we can use geometry codes without refering...//
    G4Box *solidWorld, *solidDetector, *solidScintillator, *solidScaterer,
        *solidDetectorS;
    G4LogicalVolume *logicWorld, *logicDetector, *logicScintillator,
        *logicScaterer, *logicDetectorS;
    G4VPhysicalVolume *physWorld, *physDetector, *physScintillator,
        *physScaterer;

    G4OpticalSurface *mirrorSurface;

    virtual void ConstructSDandField();

    G4int nRows, nCols; // related to fMessenger//

    G4LogicalVolume *fScoringVolume;

    G4Material *worldMat, *LYSO, *GAGG, *polyethylene, *hole_mat;
    G4Element *O, *Si, *Lu, *Y, *Gd, *Al, *Ga;

    void DefineMaterials();

    G4double xWorld, yWorld, zWorld;

    G4bool isScintillator; // isTOF;
};

#endif
