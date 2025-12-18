#include "simDetectorConstruction.hh"
#include "simDetector.hh"
// add geometry libraries
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
// add units libraries
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
// add material libraries
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"

#include "G4GenericMessenger.hh"
// add visualize libraries
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

simDetectorConstruction::simDetectorConstruction()
{
    DefineMaterials();

    isScintillator = true;

    xWorld = 2 * m;
    yWorld = 2 * m;
    zWorld = 2 * m;
}

simDetectorConstruction::~simDetectorConstruction() {}

//....oooOO0OOooo.....oooOO0OOooo.......Materials.......oooOO0OOooo.......oooOO0OOooo........//

void simDetectorConstruction::DefineMaterials()
{
    G4NistManager *nist = G4NistManager::Instance();

    worldMat = nist->FindOrBuildMaterial("G4_AIR");

    // Define Absorber detector, scintillator//
    GAGG = new G4Material("GAGG", 6.63 * g / cm3, 4);
    GAGG->AddElement(nist->FindOrBuildElement("Gd"), 3);
    GAGG->AddElement(nist->FindOrBuildElement("Al"), 2);
    GAGG->AddElement(nist->FindOrBuildElement("Ga"), 3);
    GAGG->AddElement(nist->FindOrBuildElement("O"), 12);

    // Define Scaterer detector ,compton detector//
    LYSO = new G4Material("LYSO", 7.1 * g / cm3, 4);
    LYSO->AddElement(nist->FindOrBuildElement("Lu"), 18);
    LYSO->AddElement(nist->FindOrBuildElement("Y"), 2);
    LYSO->AddElement(nist->FindOrBuildElement("Si"), 10);
    LYSO->AddElement(nist->FindOrBuildElement("O"), 50);

    // Define Derenzo phantom//
    hole_mat = nist->FindOrBuildMaterial("G4_AIR");
    polyethylene = nist->FindOrBuildMaterial("G4_POLYETHYLENE");

    // LYSO property//
    std::vector<G4double> LYSO_Energy = {
        2.2545 * eV, 2.2963 * eV, 2.3396 * eV, 2.3846 * eV, 2.4314 * eV,
        2.48 * eV,   2.5306 * eV, 2.5833 * eV, 2.6383 * eV, 2.6957 * eV,
        2.7556 * eV, 2.8182 * eV, 2.8837 * eV, 2.9524 * eV, 3.0244 * eV,
        3.1 * eV,    3.1795 * eV, 3.2632 * eV, 3.3514 * eV, 3.4444 * eV,
        3.5429 * eV};

    std::vector<G4double> LYSO_SCINT = {
        0.0123, 0.0186, 0.0313, 0.0419, 0.0567, 0.0757, 0.1115,
        0.1643, 0.2424, 0.3543, 0.5106, 0.6963, 0.8906, 0.9645,
        0.9856, 0.9898, 0.7766, 0.2171, 0.0376, 0.0039, 0.0039};

    std::vector<G4double> LYSO_RIND = {
        1.8060, 1.8078, 1.8115, 1.8133, 1.8133, 1.8133, 1.8169,
        1.8206, 1.8206, 1.8206, 1.8224, 1.8242, 1.8279, 1.8315,
        1.8352, 1.8370, 1.8370, 1.8352, 1.8352, 1.8370, 1.8388};

    std::vector<G4double> LYSO_fraction = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                           1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                           1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    std::vector<G4double> reflectivity = {
        1.0, 1.0}; // all photons should be reflected

    G4MaterialPropertiesTable *mptLYSO = new G4MaterialPropertiesTable();
    mptLYSO->AddProperty("RINDEX", LYSO_Energy, LYSO_RIND);
    mptLYSO->AddProperty("FASTCOMPONENT", LYSO_Energy, LYSO_fraction,
                         true); // fraction of wavelength, its tell to G4 how
                                // many photons for each wavelenges create
    mptLYSO->AddConstProperty(
        "SCINTILLATIONYIELD",
        30. / MeV); // average scintillation yield is 30K/MeV, :how many photon
                    // per loss energy create , using 'const': we haven't array
    mptLYSO->AddConstProperty("RESOLUTIONSCALE",
                              1.0); // regarding of the photon ratio
    mptLYSO->AddConstProperty("SCINTILLATIONTIMECONSTANT", 45. * ns,
                              true); // scintillation decay time is 45ns
    mptLYSO->AddConstProperty(
        "YIELDRATIO", 1.0,
        true); // YIELDRATIO: related to distribution of photons ,
    LYSO->SetMaterialPropertiesTable(mptLYSO);

    G4MaterialPropertiesTable *mptWorld = new G4MaterialPropertiesTable();

    // mptWorld-> AddProperty("RINDEX",energy,rindexWorld,2);
    // worldMat->SetMaterialPropertiesTable(mptWorld);

    // GAGG property//
    G4double rindexWorld[2] = {1.0, 1.0};
    std::vector<G4double> GAGG_Energy = {662 * keV / 0.07, 662 * keV / 0.04};
    std::vector<G4double> GAGG_RIND = {1.91, 1.91};
    std::vector<G4double> GAGG_fraction = {1.0, 1.0};

    G4MaterialPropertiesTable *mptGAGG = new G4MaterialPropertiesTable();
    mptGAGG->AddProperty("RINDEX", GAGG_Energy, GAGG_RIND);
    mptGAGG->AddProperty("FASTCOMPONENT", GAGG_Energy, GAGG_fraction, true);
    mptGAGG->AddConstProperty("SCINTILLATIONYIELD", 60.0 / keV);
    mptGAGG->AddConstProperty("RESOLUTIONSCALE", 1.0);
    mptGAGG->AddConstProperty("SCINTILLATIONTIMECONSTANT", 50. * ns, true);
    mptGAGG->AddConstProperty("YIELDRATIO", 1.0, true);
    GAGG->SetMaterialPropertiesTable(mptGAGG);

    mirrorSurface = new G4OpticalSurface("mirrorSurface");
    mirrorSurface->SetType(dielectric_metal);
    mirrorSurface->SetFinish(ground);
    mirrorSurface->SetModel(unified);

    G4MaterialPropertiesTable *mptMirror = new G4MaterialPropertiesTable();
    mptMirror->AddProperty("REFLECTIVITY", GAGG_Energy, reflectivity,
                           2); // in G4 REFLECTIVITY: what is the fraction of
                               // the light which is reflected from surface

    mirrorSurface->SetMaterialPropertiesTable(mptMirror);
}

//.....oooOO0OOooo.........oooOO0OOooo.....Scintillator(absorber).....oooOO0OOooo.........oooOO0OOooo.....//

void simDetectorConstruction::ConstructScintillator()
{
    G4bool checkOverlaps = false;

    fScoringVolume = logicScintillator;

    // Compton detector(Scaterer)//

    G4double length = 60. * cm;

    G4double radius = 40 * cm;

    G4double xCryst = 2 * mm;
    G4double yCryst = 2 * mm;
    G4double zCryst = 2 * mm;

    G4double gap = 0.1 * mm;

    G4double xDet = xCryst;
    G4double yDet = yCryst;
    G4double zDet = zCryst;

    solidScaterer = new G4Box("Scaterer", xCryst, yCryst, zCryst);
    logicScaterer = new G4LogicalVolume(solidScaterer, LYSO, "logicalScaterer");

    solidScintillator = new G4Box("solidScintillator", xDet, yDet, zDet);

    logicScintillator =
        new G4LogicalVolume(solidScintillator, GAGG, "logicalScintillator");

    G4LogicalSkinSurface *skin = new G4LogicalSkinSurface(
        "skin", logicWorld,
        mirrorSurface); // applying mirror surface(reflectivity) to Scintillator

    const double pi = 3.14159265358979323846;

    G4int nCrystR = (pi * radius) / (yCryst + gap);
    G4double dAngle = 360. / nCrystR;

    G4int nCrystL = length / (zCryst + gap);

    for (G4int i = 0; i < nCrystL; i++)
    {
        for (G4int j = 0; j < nCrystR; j++)
        {
            G4ThreeVector trans = G4ThreeVector(
                radius, 0., 2 * (-nCrystL / 2 + i) * (zCryst + gap));

            G4Rotate3D rotZ(j * dAngle * deg, G4ThreeVector(0, 0, 1));
            G4Translate3D transScaterer(trans);
            G4Transform3D transformScaterer = (rotZ) * (transScaterer);

            physScaterer = new G4PVPlacement(transformScaterer, logicScaterer,
                                             "physScaterer", logicWorld, false,
                                             i * 16 + j, checkOverlaps);

            G4Translate3D transScint(trans +
                                     G4ThreeVector(xCryst + xDet, 0., 0.));
            G4Transform3D transformScint = (rotZ) * (transScint);

            physScintillator = new G4PVPlacement(
                transformScint, logicScintillator, "physScintillator",
                logicWorld, false, i * 16 + j, checkOverlaps);
        }
    }
}

//.....oooOO0OOooo.........oooOO0OOooo.....WORD.....oooOO0OOooo.........oooOO0OOooo.......//

G4VPhysicalVolume *simDetectorConstruction::Construct()
{
    solidWorld = new G4Box("solidWorld", xWorld, yWorld, zWorld);
    logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
    physWorld = new G4PVPlacement(0,                      // no rotation
                                  G4ThreeVector(0, 0, 0), // at (0,0,0)
                                  logicWorld,  // its logical volume name
                                  "physWorld", // its name
                                  0,           // its mother  volume
                                  false,       // no boolean operation
                                  0,           // copy number
                                  true);       // overlaps checking

    //.....oooOO0OOooo.........oooOO0OOooo.....Derenzo
    // phantom.....oooOO0OOooo.........oooOO0OOooo.......//

    // define a cylynder
    G4Tubs *solidPhantom = new G4Tubs("SPhantom",   // name
                                      0. * mm,      // min raduis
                                      25. * mm,     // max raduis
                                      18.5 * mm,    // half-z
                                      0.0 * deg,    // phi-1
                                      360.0 * deg); // phi-2
    G4LogicalVolume *LogicPhantom =
        new G4LogicalVolume(solidPhantom, polyethylene, "SPhantom");
    G4VPhysicalVolume *PhysPhantom =
        new G4PVPlacement(0, G4ThreeVector(0, 0, 0 * cm), LogicPhantom,
                          "PhysPhantom", logicWorld, false, 3, true);

    // define holes
    G4Tubs *solidHole1 = new G4Tubs("SHole1", 0. * mm, 1.1947 * mm, 18.5 * mm,
                                    0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole1 =
        new G4LogicalVolume(solidHole1, hole_mat, "SHole1");
    G4VPhysicalVolume *PhysHole1 = new G4PVPlacement(
        0, G4ThreeVector(-19.4395 * mm, -2.9528 * mm, 0. * mm), LogicHole1,
        "PhysHole1", LogicPhantom, false, 0, true);

    G4Tubs *solidHole2 = new G4Tubs("SHole2", 0. * mm, 1.5881 * mm, 18.5 * mm,
                                    0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole2 =
        new G4LogicalVolume(solidHole2, hole_mat, "SHole2");
    G4VPhysicalVolume *PhysHole2 = new G4PVPlacement(
        0, G4ThreeVector(-17.8383 * mm, 3.0455 * mm, 0. * mm), LogicHole2,
        "PhysHole2", LogicPhantom, false, 0, true);

    G4Tubs *solidHole3 = new G4Tubs("SHole3", 0. * mm, 1.1947 * mm, 18.5 * mm,
                                    0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole3 =
        new G4LogicalVolume(solidHole3, hole_mat, "SHole3");
    G4VPhysicalVolume *PhysHole3 = new G4PVPlacement(
        0, G4ThreeVector(-17.0353 * mm, -7.1047 * mm, 0. * mm), LogicHole3,
        "PhysHole3", LogicPhantom, false, 0, true);

    G4Tubs *solidHole4 = new G4Tubs("SHole4", 0. * mm, 1.5881 * mm, 18.5 * mm,
                                    0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole4 =
        new G4LogicalVolume(solidHole4, hole_mat, "SHole4");
    G4VPhysicalVolume *PhysHole4 = new G4PVPlacement(
        0, G4ThreeVector(-14.6379 * mm, 8.5824 * mm, 0. * mm), LogicHole4,
        "PhysHole4", LogicPhantom, false, 0, true);

    G4Tubs *solidHole5 = new G4Tubs("SHole5", 0. * mm, 1.1959 * mm, 18.5 * mm,
                                    0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole5 =
        new G4LogicalVolume(solidHole5, hole_mat, "SHole5");
    G4VPhysicalVolume *PhysHole5 = new G4PVPlacement(
        0, G4ThreeVector(-14.6368 * mm, -11.2587 * mm, 0. * mm), LogicHole5,
        "PhysHole5", LogicPhantom, false, 0, true);

    G4Tubs *solidHole6 = new G4Tubs("SHole6", 0. * mm, 1.1898 * mm, 18.5 * mm,
                                    0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole6 =
        new G4LogicalVolume(solidHole6, hole_mat, "SHole6");
    G4VPhysicalVolume *PhysHole6 = new G4PVPlacement(
        0, G4ThreeVector(-14.6446 * mm, -2.9494 * mm, 0. * mm), LogicHole6,
        "PhysHole6", LogicPhantom, false, 0, true);

    G4Tubs *solidHole7 = new G4Tubs("SHole7", 0. * mm, 1.1898 * mm, 18.5 * mm,
                                    0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole7 =
        new G4LogicalVolume(solidHole7, hole_mat, "SHole7");
    G4VPhysicalVolume *PhysHole7 = new G4PVPlacement(
        0, G4ThreeVector(-12.2429 * mm, -15.4138 * mm, 0. * mm), LogicHole7,
        "PhysHole7", LogicPhantom, false, 0, true);

    G4Tubs *solidHole8 = new G4Tubs("SHole8", 0. * mm, 1.1935 * mm, 18.5 * mm,
                                    0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole8 =
        new G4LogicalVolume(solidHole8, hole_mat, "SHole8");
    G4VPhysicalVolume *PhysHole8 = new G4PVPlacement(
        0, G4ThreeVector(-12.2431 * mm, -7.1035 * mm, 0. * mm), LogicHole8,
        "PhysHole8", LogicPhantom, false, 0, true);

    G4Tubs *solidHole9 = new G4Tubs("SHole9", 0. * mm, 1.5937 * mm, 18.5 * mm,
                                    0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole9 =
        new G4LogicalVolume(solidHole9, hole_mat, "SHole9");
    G4VPhysicalVolume *PhysHole9 = new G4PVPlacement(
        0, G4ThreeVector(-11.4423 * mm, 3.0511 * mm, 0. * mm), LogicHole9,
        "PhysHole9", LogicPhantom, false, 0, true);

    G4Tubs *solidHole10 = new G4Tubs("SHole10", 0. * mm, 1.5918 * mm, 18.5 * mm,
                                     0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole10 =
        new G4LogicalVolume(solidHole10, hole_mat, "SHole10");
    G4VPhysicalVolume *PhysHole10 = new G4PVPlacement(
        0, G4ThreeVector(-11.445 * mm, 14.1224 * mm, 0. * mm), LogicHole10,
        "PhysHole10", LogicPhantom, false, 0, true);

    G4Tubs *solidHole11 = new G4Tubs("SHole11", 0. * mm, 1.191 * mm, 18.5 * mm,
                                     0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole11 =
        new G4LogicalVolume(solidHole11, hole_mat, "SHole11");
    G4VPhysicalVolume *PhysHole11 = new G4PVPlacement(
        0, G4ThreeVector(-9.846 * mm, -11.2564 * mm, 0. * mm), LogicHole11,
        "PhysHole11", LogicPhantom, false, 0, true);

    G4Tubs *solidHole12 = new G4Tubs("SHole12", 0. * mm, 1.1898 * mm, 18.5 * mm,
                                     0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole12 =
        new G4LogicalVolume(solidHole12, hole_mat, "SHole12");
    G4VPhysicalVolume *PhysHole12 = new G4PVPlacement(
        0, G4ThreeVector(-9.8453 * mm, -2.9458 * mm, 0. * mm), LogicHole12,
        "PhysHole12", LogicPhantom, false, 0, true);

    G4Tubs *solidHole13 = new G4Tubs("SHole13", 0. * mm, 1.5881 * mm, 18.5 * mm,
                                     0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole13 =
        new G4LogicalVolume(solidHole13, hole_mat, "SHole13");
    G4VPhysicalVolume *PhysHole13 = new G4PVPlacement(
        0, G4ThreeVector(-8.2467 * mm, 8.5824 * mm, 0. * mm), LogicHole13,
        "PhysHole13", LogicPhantom, false, 0, true);

    G4Tubs *solidHole14 = new G4Tubs("SHole14", 0. * mm, 0.79147 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole14 =
        new G4LogicalVolume(solidHole14, hole_mat, "SHole14");
    G4VPhysicalVolume *PhysHole14 = new G4PVPlacement(
        0, G4ThreeVector(-7.8485 * mm, -19.7874 * mm, 0. * mm), LogicHole14,
        "PhysHole14", LogicPhantom, false, 0, true);

    G4Tubs *solidHole15 = new G4Tubs("SHole15", 0. * mm, 1.1898 * mm, 18.5 * mm,
                                     0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole15 =
        new G4LogicalVolume(solidHole15, hole_mat, "SHole15");
    G4VPhysicalVolume *PhysHole15 = new G4PVPlacement(
        0, G4ThreeVector(-7.4464 * mm, -7.1032 * mm, 0. * mm), LogicHole15,
        "PhysHole15", LogicPhantom, false, 0, true);

    G4Tubs *solidHole16 = new G4Tubs("SHole16", 0. * mm, 0.80281 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole16 =
        new G4LogicalVolume(solidHole16, hole_mat, "SHole16");
    G4VPhysicalVolume *PhysHole16 = new G4PVPlacement(
        0, G4ThreeVector(-6.25 * mm, -17.0263 * mm, 0. * mm), LogicHole16,
        "PhysHole16", LogicPhantom, false, 0, true);

    G4Tubs *solidHole17 = new G4Tubs("SHole17", 0. * mm, 1.5872 * mm, 18.5 * mm,
                                     0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole17 =
        new G4LogicalVolume(solidHole17, hole_mat, "SHole17");
    G4VPhysicalVolume *PhysHole17 = new G4PVPlacement(
        0, G4ThreeVector(-5.0481 * mm, 3.0461 * mm, 0. * mm), LogicHole17,
        "PhysHole17", LogicPhantom, false, 0, true);

    G4Tubs *solidHole18 = new G4Tubs("SHole18", 0. * mm, 1.1861 * mm, 18.5 * mm,
                                     0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole18 =
        new G4LogicalVolume(solidHole18, hole_mat, "SHole18");
    G4VPhysicalVolume *PhysHole18 = new G4PVPlacement(
        0, G4ThreeVector(-5.0481 * mm, -2.951 * mm, 0. * mm), LogicHole18,
        "PhysHole18", LogicPhantom, false, 0, true);

    G4Tubs *solidHole19 = new G4Tubs("SHole19", 0. * mm, 1.9888 * mm, 18.5 * mm,
                                     0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole19 =
        new G4LogicalVolume(solidHole19, hole_mat, "SHole19");
    G4VPhysicalVolume *PhysHole19 = new G4PVPlacement(
        0, G4ThreeVector(-3.8512 * mm, 14.9687 * mm, 0. * mm), LogicHole19,
        "PhysHole19", LogicPhantom, false, 0, true);

    G4Tubs *solidHole20 = new G4Tubs("SHole20", 0. * mm, 0.79359 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole20 =
        new G4LogicalVolume(solidHole20, hole_mat, "SHole20");
    G4VPhysicalVolume *PhysHole20 = new G4PVPlacement(
        0, G4ThreeVector(-4.6549 * mm, -19.7888 * mm, 0. * mm), LogicHole20,
        "PhysHole20", LogicPhantom, false, 0, true);

    G4Tubs *solidHole21 = new G4Tubs("SHole21", 0. * mm, 0.79174 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole21 =
        new G4LogicalVolume(solidHole21, hole_mat, "SHole21");
    G4VPhysicalVolume *PhysHole21 = new G4PVPlacement(
        0, G4ThreeVector(-4.6574 * mm, -14.2637 * mm, 0. * mm), LogicHole21,
        "PhysHole21", LogicPhantom, false, 0, true);

    G4Tubs *solidHole22 = new G4Tubs("SHole22", 0. * mm, 0.79729 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole22 =
        new G4LogicalVolume(solidHole22, hole_mat, "SHole22");
    G4VPhysicalVolume *PhysHole22 = new G4PVPlacement(
        0, G4ThreeVector(-3.056 * mm, -17.0272 * mm, 0. * mm), LogicHole22,
        "PhysHole22", LogicPhantom, false, 0, true);

    G4Tubs *solidHole23 = new G4Tubs("SHole23", 0. * mm, 0.79147 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole23 =
        new G4LogicalVolume(solidHole23, hole_mat, "SHole23");
    G4VPhysicalVolume *PhysHole23 = new G4PVPlacement(
        0, G4ThreeVector(-3.0426 * mm, -11.4902 * mm, 0. * mm), LogicHole23,
        "PhysHole23", LogicPhantom, false, 0, true);

    G4Tubs *solidHole24 = new G4Tubs("SHole24", 0. * mm, 0.80098 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole24 =
        new G4LogicalVolume(solidHole24, hole_mat, "SHole24");
    G4VPhysicalVolume *PhysHole24 = new G4PVPlacement(
        0, G4ThreeVector(-1.452 * mm, -19.798 * mm, 0. * mm), LogicHole24,
        "PhysHole24", LogicPhantom, false, 0, true);

    G4Tubs *solidHole25 = new G4Tubs("SHole25", 0. * mm, 0.79359 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole25 =
        new G4LogicalVolume(solidHole25, hole_mat, "SHole25");
    G4VPhysicalVolume *PhysHole25 = new G4PVPlacement(
        0, G4ThreeVector(-1.471 * mm, -14.2528 * mm, 0. * mm), LogicHole25,
        "PhysHole25", LogicPhantom, false, 0, true);

    G4Tubs *solidHole26 = new G4Tubs("SHole26", 0. * mm, 0.79147 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole26 =
        new G4LogicalVolume(solidHole26, hole_mat, "SHole26");
    G4VPhysicalVolume *PhysHole26 = new G4PVPlacement(
        0, G4ThreeVector(-1.4568 * mm, -8.7175 * mm, 0. * mm), LogicHole26,
        "PhysHole26", LogicPhantom, false, 0, true);

    G4Tubs *solidHole27 = new G4Tubs("SHole27", 0. * mm, 1.9858 * mm, 18.5 * mm,
                                     0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole27 =
        new G4LogicalVolume(solidHole27, hole_mat, "SHole27");
    G4VPhysicalVolume *PhysHole27 = new G4PVPlacement(
        0, G4ThreeVector(0.14423 * mm, 8.0448 * mm, 0. * mm), LogicHole27,
        "PhysHole27", LogicPhantom, false, 0, true);

    G4Tubs *solidHole28 = new G4Tubs("SHole28", 0. * mm, 0.70729 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole28 =
        new G4LogicalVolume(solidHole28, hole_mat, "SHole28");
    G4VPhysicalVolume *PhysHole28 = new G4PVPlacement(
        0, G4ThreeVector(0.14423 * mm, -17.0102 * mm, 0. * mm), LogicHole28,
        "PhysHole28", LogicPhantom, false, 0, true);

    G4Tubs *solidHole29 = new G4Tubs("SHole29", 0. * mm, 0.78614 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole29 =
        new G4LogicalVolume(solidHole29, hole_mat, "SHole29");
    G4VPhysicalVolume *PhysHole29 = new G4PVPlacement(
        0, G4ThreeVector(0.14423 * mm, -11.4835 * mm, 0. * mm), LogicHole29,
        "PhysHole29", LogicPhantom, false, 0, true);

    G4Tubs *solidHole30 = new G4Tubs("SHole30", 0. * mm, 0.79729 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole30 =
        new G4LogicalVolume(solidHole30, hole_mat, "SHole30");
    G4VPhysicalVolume *PhysHole30 = new G4PVPlacement(
        0, G4ThreeVector(0.14423 * mm, -5.9366 * mm, 0. * mm), LogicHole30,
        "PhysHole30", LogicPhantom, false, 0, true);

    G4Tubs *solidHole31 = new G4Tubs("SHole31", 0. * mm, 0.80098 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole31 =
        new G4LogicalVolume(solidHole31, hole_mat, "SHole31");
    G4VPhysicalVolume *PhysHole31 = new G4PVPlacement(
        0, G4ThreeVector(1.7405 * mm, -19.798 * mm, 0. * mm), LogicHole31,
        "PhysHole31", LogicPhantom, false, 0, true);

    G4Tubs *solidHole32 = new G4Tubs("SHole32", 0. * mm, 0.79359 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole32 =
        new G4LogicalVolume(solidHole32, hole_mat, "SHole32");
    G4VPhysicalVolume *PhysHole32 = new G4PVPlacement(
        0, G4ThreeVector(1.7456 * mm, -14.2528 * mm, 0. * mm), LogicHole32,
        "PhysHole32", LogicPhantom, false, 0, true);

    G4Tubs *solidHole33 = new G4Tubs("SHole33", 0. * mm, 0.79147 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole33 =
        new G4LogicalVolume(solidHole33, hole_mat, "SHole33");
    G4VPhysicalVolume *PhysHole33 = new G4PVPlacement(
        0, G4ThreeVector(1.7452 * mm, -8.7175 * mm, 0. * mm), LogicHole33,
        "PhysHole33", LogicPhantom, false, 0, true);

    G4Tubs *solidHole34 = new G4Tubs("SHole34", 0. * mm, 1.9888 * mm, 18.5 * mm,
                                     0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole34 =
        new G4LogicalVolume(solidHole34, hole_mat, "SHole34");
    G4VPhysicalVolume *PhysHole34 = new G4PVPlacement(
        0, G4ThreeVector(4.1396 * mm, 14.9687 * mm, 0. * mm), LogicHole34,
        "PhysHole34", LogicPhantom, false, 0, true);

    G4Tubs *solidHole35 = new G4Tubs("SHole35", 0. * mm, 0.79729 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole35 =
        new G4LogicalVolume(solidHole35, hole_mat, "SHole35");
    G4VPhysicalVolume *PhysHole35 = new G4PVPlacement(
        0, G4ThreeVector(3.3445 * mm, -17.0272 * mm, 0. * mm), LogicHole35,
        "PhysHole35", LogicPhantom, false, 0, true);

    G4Tubs *solidHole36 = new G4Tubs("SHole36", 0. * mm, 0.79147 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole36 =
        new G4LogicalVolume(solidHole36, hole_mat, "SHole36");
    G4VPhysicalVolume *PhysHole36 = new G4PVPlacement(
        0, G4ThreeVector(3.3311 * mm, -11.4902 * mm, 0. * mm), LogicHole36,
        "PhysHole36", LogicPhantom, false, 0, true);

    G4Tubs *solidHole37 = new G4Tubs("SHole37", 0. * mm, 0.79359 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole37 =
        new G4LogicalVolume(solidHole37, hole_mat, "SHole37");
    G4VPhysicalVolume *PhysHole37 = new G4PVPlacement(
        0, G4ThreeVector(4.9434 * mm, -19.7888 * mm, 0. * mm), LogicHole37,
        "PhysHole37", LogicPhantom, false, 0, true);

    G4Tubs *solidHole38 = new G4Tubs("SHole38", 0. * mm, 0.79147 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole38 =
        new G4LogicalVolume(solidHole38, hole_mat, "SHole38");
    G4VPhysicalVolume *PhysHole38 = new G4PVPlacement(
        0, G4ThreeVector(4.9458 * mm, -14.2637 * mm, 0. * mm), LogicHole38,
        "PhysHole38", LogicPhantom, false, 0, true);

    G4Tubs *solidHole39 = new G4Tubs("SHole39", 0. * mm, 2.3833 * mm, 18.5 * mm,
                                     0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole39 =
        new G4LogicalVolume(solidHole39, hole_mat, "SHole39");
    G4VPhysicalVolume *PhysHole39 = new G4PVPlacement(
        0, G4ThreeVector(7.0673 * mm, 4.0426 * mm, 0. * mm), LogicHole39,
        "PhysHole39", LogicPhantom, false, 0, true);

    G4Tubs *solidHole40 = new G4Tubs("SHole40", 0. * mm, 0.58929 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole40 =
        new G4LogicalVolume(solidHole40, hole_mat, "SHole40");
    G4VPhysicalVolume *PhysHole40 = new G4PVPlacement(
        0, G4ThreeVector(5.3365 * mm, -2.9417 * mm, 0. * mm), LogicHole40,
        "PhysHole40", LogicPhantom, false, 0, true);

    G4Tubs *solidHole41 = new G4Tubs("SHole41", 0. * mm, 0.80281 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole41 =
        new G4LogicalVolume(solidHole41, hole_mat, "SHole41");
    G4VPhysicalVolume *PhysHole41 = new G4PVPlacement(
        0, G4ThreeVector(6.5385 * mm, -17.0263 * mm, 0. * mm), LogicHole41,
        "PhysHole41", LogicPhantom, false, 0, true);

    G4Tubs *solidHole42 = new G4Tubs("SHole42", 0. * mm, 0.5992 * mm, 18.5 * mm,
                                     0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole42 =
        new G4LogicalVolume(solidHole42, hole_mat, "SHole42");
    G4VPhysicalVolume *PhysHole42 = new G4PVPlacement(
        0, G4ThreeVector(6.5306 * mm, -5.0292 * mm, 0. * mm), LogicHole42,
        "PhysHole42", LogicPhantom, false, 0, true);

    G4Tubs *solidHole43 = new G4Tubs("SHole43", 0. * mm, 0.5992 * mm, 18.5 * mm,
                                     0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole43 =
        new G4LogicalVolume(solidHole43, hole_mat, "SHole43");
    G4VPhysicalVolume *PhysHole43 = new G4PVPlacement(
        0, G4ThreeVector(7.7333 * mm, -7.0996 * mm, 0. * mm), LogicHole43,
        "PhysHole43", LogicPhantom, false, 0, true);

    G4Tubs *solidHole44 = new G4Tubs("SHole44", 0. * mm, 0.59427 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole44 =
        new G4LogicalVolume(solidHole44, hole_mat, "SHole44");
    G4VPhysicalVolume *PhysHole44 = new G4PVPlacement(
        0, G4ThreeVector(7.7348 * mm, -2.9407 * mm, 0. * mm), LogicHole44,
        "PhysHole44", LogicPhantom, false, 0, true);

    G4Tubs *solidHole45 = new G4Tubs("SHole45", 0. * mm, 0.79174 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole45 =
        new G4LogicalVolume(solidHole45, hole_mat, "SHole45");
    G4VPhysicalVolume *PhysHole45 = new G4PVPlacement(
        0, G4ThreeVector(8.137 * mm, -19.7874 * mm, 0. * mm), LogicHole45,
        "PhysHole45", LogicPhantom, false, 0, true);

    G4Tubs *solidHole46 = new G4Tubs("SHole46", 0. * mm, 0.60165 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole46 =
        new G4LogicalVolume(solidHole46, hole_mat, "SHole46");
    G4VPhysicalVolume *PhysHole46 = new G4PVPlacement(
        0, G4ThreeVector(8.9306 * mm, -9.1784 * mm, 0. * mm), LogicHole46,
        "PhysHole46", LogicPhantom, false, 0, true);

    G4Tubs *solidHole47 = new G4Tubs("SHole47", 0. * mm, 0.5992 * mm, 18.5 * mm,
                                     0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole47 =
        new G4LogicalVolume(solidHole47, hole_mat, "SHole47");
    G4VPhysicalVolume *PhysHole47 = new G4PVPlacement(
        0, G4ThreeVector(8.9258 * mm, -5.0331 * mm, 0. * mm), LogicHole47,
        "PhysHole47", LogicPhantom, false, 0, true);

    G4Tubs *solidHole48 = new G4Tubs("SHole48", 0. * mm, 2.3845 * mm, 18.5 * mm,
                                     0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole48 =
        new G4LogicalVolume(solidHole48, hole_mat, "SHole48");
    G4VPhysicalVolume *PhysHole48 = new G4PVPlacement(
        0, G4ThreeVector(11.8615 * mm, 12.3535 * mm, 0. * mm), LogicHole48,
        "PhysHole48", LogicPhantom, false, 0, true);

    G4Tubs *solidHole49 = new G4Tubs("SHole49", 0. * mm, 0.60165 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole49 =
        new G4LogicalVolume(solidHole49, hole_mat, "SHole49");
    G4VPhysicalVolume *PhysHole49 = new G4PVPlacement(
        0, G4ThreeVector(10.1329 * mm, -11.2633 * mm, 0. * mm), LogicHole49,
        "PhysHole49", LogicPhantom, false, 0, true);

    G4Tubs *solidHole50 = new G4Tubs("SHole50", 0. * mm, 0.60165 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole50 =
        new G4LogicalVolume(solidHole50, hole_mat, "SHole50");
    G4VPhysicalVolume *PhysHole50 = new G4PVPlacement(
        0, G4ThreeVector(10.1329 * mm, -7.1021 * mm, 0. * mm), LogicHole50,
        "PhysHole50", LogicPhantom, false, 0, true);

    G4Tubs *solidHole51 = new G4Tubs("SHole51", 0. * mm, 0.59179 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole51 =
        new G4LogicalVolume(solidHole51, hole_mat, "SHole51");
    G4VPhysicalVolume *PhysHole51 = new G4PVPlacement(
        0, G4ThreeVector(10.1341 * mm, -2.9428 * mm, 0. * mm), LogicHole51,
        "PhysHole51", LogicPhantom, false, 0, true);

    G4Tubs *solidHole52 = new G4Tubs("SHole52", 0. * mm, 0.59674 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole52 =
        new G4LogicalVolume(solidHole52, hole_mat, "SHole52");
    G4VPhysicalVolume *PhysHole52 = new G4PVPlacement(
        0, G4ThreeVector(11.3247 * mm, -13.3328 * mm, 0. * mm), LogicHole52,
        "PhysHole52", LogicPhantom, false, 0, true);

    G4Tubs *solidHole53 = new G4Tubs("SHole53", 0. * mm, 0.60165 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole53 =
        new G4LogicalVolume(solidHole53, hole_mat, "SHole53");
    G4VPhysicalVolume *PhysHole53 = new G4PVPlacement(
        0, G4ThreeVector(11.3263 * mm, -9.1783 * mm, 0. * mm), LogicHole53,
        "PhysHole53", LogicPhantom, false, 0, true);

    G4Tubs *solidHole54 = new G4Tubs("SHole54", 0. * mm, 0.59674 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole54 =
        new G4LogicalVolume(solidHole54, hole_mat, "SHole54");
    G4VPhysicalVolume *PhysHole54 = new G4PVPlacement(
        0, G4ThreeVector(11.3247 * mm, -5.0326 * mm, 0. * mm), LogicHole54,
        "PhysHole54", LogicPhantom, false, 0, true);

    G4Tubs *solidHole55 = new G4Tubs("SHole55", 0. * mm, 0.59179 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole55 =
        new G4LogicalVolume(solidHole55, hole_mat, "SHole55");
    G4VPhysicalVolume *PhysHole55 = new G4PVPlacement(
        0, G4ThreeVector(12.538 * mm, -15.4226 * mm, 0. * mm), LogicHole55,
        "PhysHole55", LogicPhantom, false, 0, true);

    G4Tubs *solidHole56 = new G4Tubs("SHole56", 0. * mm, 0.5992 * mm, 18.5 * mm,
                                     0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole56 =
        new G4LogicalVolume(solidHole56, hole_mat, "SHole56");
    G4VPhysicalVolume *PhysHole56 = new G4PVPlacement(
        0, G4ThreeVector(12.5331 * mm, -11.2666 * mm, 0. * mm), LogicHole56,
        "PhysHole56", LogicPhantom, false, 0, true);

    G4Tubs *solidHole57 = new G4Tubs("SHole57", 0. * mm, 0.5992 * mm, 18.5 * mm,
                                     0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole57 =
        new G4LogicalVolume(solidHole57, hole_mat, "SHole57");
    G4VPhysicalVolume *PhysHole57 = new G4PVPlacement(
        0, G4ThreeVector(12.5331 * mm, -7.0988 * mm, 0. * mm), LogicHole57,
        "PhysHole57", LogicPhantom, false, 0, true);

    G4Tubs *solidHole58 = new G4Tubs("SHole58", 0. * mm, 0.59179 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole58 =
        new G4LogicalVolume(solidHole58, hole_mat, "SHole58");
    G4VPhysicalVolume *PhysHole58 = new G4PVPlacement(
        0, G4ThreeVector(12.538 * mm, -2.9428 * mm, 0. * mm), LogicHole58,
        "PhysHole58", LogicPhantom, false, 0, true);

    G4Tubs *solidHole59 = new G4Tubs("SHole50", 0. * mm, 0.60165 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole59 =
        new G4LogicalVolume(solidHole59, hole_mat, "SHole59");
    G4VPhysicalVolume *PhysHole59 = new G4PVPlacement(
        0, G4ThreeVector(13.7219 * mm, -17.4836 * mm, 0. * mm), LogicHole59,
        "PhysHole59", LogicPhantom, false, 0, true);

    G4Tubs *solidHole60 = new G4Tubs("SHole60", 0. * mm, 0.59427 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole60 =
        new G4LogicalVolume(solidHole60, hole_mat, "SHole60");
    G4VPhysicalVolume *PhysHole60 = new G4PVPlacement(
        0, G4ThreeVector(13.7236 * mm, -13.3325 * mm, 0. * mm), LogicHole60,
        "PhysHole60", LogicPhantom, false, 0, true);

    G4Tubs *solidHole61 = new G4Tubs("SHole61", 0. * mm, 0.59427 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole61 =
        new G4LogicalVolume(solidHole61, hole_mat, "SHole61");
    G4VPhysicalVolume *PhysHole61 = new G4PVPlacement(
        0, G4ThreeVector(13.7179 * mm, -9.1827 * mm, 0. * mm), LogicHole61,
        "PhysHole61", LogicPhantom, false, 0, true);

    G4Tubs *solidHole62 = new G4Tubs("SHole62", 0. * mm, 0.59427 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole62 =
        new G4LogicalVolume(solidHole62, hole_mat, "SHole62");
    G4VPhysicalVolume *PhysHole62 = new G4PVPlacement(
        0, G4ThreeVector(13.7236 * mm, -5.0329 * mm, 0. * mm), LogicHole62,
        "PhysHole62", LogicPhantom, false, 0, true);

    G4Tubs *solidHole63 = new G4Tubs("SHole63", 0. * mm, 2.3833 * mm, 18.5 * mm,
                                     0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole63 =
        new G4LogicalVolume(solidHole63, hole_mat, "SHole63");
    G4VPhysicalVolume *PhysHole63 = new G4PVPlacement(
        0, G4ThreeVector(16.6608 * mm, 4.0447 * mm, 0. * mm), LogicHole63,
        "PhysHole63", LogicPhantom, false, 0, true);

    G4Tubs *solidHole64 = new G4Tubs("SHole64", 0. * mm, 0.59674 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole64 =
        new G4LogicalVolume(solidHole64, hole_mat, "SHole64");
    G4VPhysicalVolume *PhysHole64 = new G4PVPlacement(
        0, G4ThreeVector(14.9348 * mm, -15.4156 * mm, 0. * mm), LogicHole64,
        "PhysHole64", LogicPhantom, false, 0, true);

    G4Tubs *solidHole65 = new G4Tubs("SHole65", 0. * mm, 0.5992 * mm, 18.5 * mm,
                                     0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole65 =
        new G4LogicalVolume(solidHole65, hole_mat, "SHole65");
    G4VPhysicalVolume *PhysHole65 = new G4PVPlacement(
        0, G4ThreeVector(14.9369 * mm, -11.2666 * mm, 0. * mm), LogicHole65,
        "PhysHole65", LogicPhantom, false, 0, true);

    G4Tubs *solidHole66 = new G4Tubs("SHole66", 0. * mm, 0.5992 * mm, 18.5 * mm,
                                     0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole66 =
        new G4LogicalVolume(solidHole66, hole_mat, "SHole66");
    G4VPhysicalVolume *PhysHole66 = new G4PVPlacement(
        0, G4ThreeVector(14.9369 * mm, -7.0988 * mm, 0. * mm), LogicHole66,
        "PhysHole66", LogicPhantom, false, 0, true);

    G4Tubs *solidHole67 = new G4Tubs("SHole67", 0. * mm, 0.58929 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole67 =
        new G4LogicalVolume(solidHole67, hole_mat, "SHole67");
    G4VPhysicalVolume *PhysHole67 = new G4PVPlacement(
        0, G4ThreeVector(14.9373 * mm, -2.9408 * mm, 0. * mm), LogicHole67,
        "PhysHole67", LogicPhantom, false, 0, true);

    G4Tubs *solidHole68 = new G4Tubs("SHole68", 0. * mm, 0.58929 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole68 =
        new G4LogicalVolume(solidHole68, hole_mat, "SHole68");
    G4VPhysicalVolume *PhysHole68 = new G4PVPlacement(
        0, G4ThreeVector(16.1204 * mm, -13.3255 * mm, 0. * mm), LogicHole68,
        "PhysHole68", LogicPhantom, false, 0, true);

    G4Tubs *solidHole69 = new G4Tubs("SHole69", 0. * mm, 0.59179 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole69 =
        new G4LogicalVolume(solidHole69, hole_mat, "SHole69");
    G4VPhysicalVolume *PhysHole69 = new G4PVPlacement(
        0, G4ThreeVector(16.1183 * mm, -9.1791 * mm, 0. * mm), LogicHole69,
        "PhysHole69", LogicPhantom, false, 0, true);

    G4Tubs *solidHole70 = new G4Tubs("SHole70", 0. * mm, 0.59674 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole70 =
        new G4LogicalVolume(solidHole70, hole_mat, "SHole70");
    G4VPhysicalVolume *PhysHole70 = new G4PVPlacement(
        0, G4ThreeVector(16.1229 * mm, -5.031 * mm, 0. * mm), LogicHole70,
        "PhysHole70", LogicPhantom, false, 0, true);

    G4Tubs *solidHole71 = new G4Tubs("SHole71", 0. * mm, 0.59674 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole71 =
        new G4LogicalVolume(solidHole71, hole_mat, "SHole71");
    G4VPhysicalVolume *PhysHole71 = new G4PVPlacement(
        0, G4ThreeVector(17.3291 * mm, -11.2595 * mm, 0. * mm), LogicHole71,
        "PhysHole71", LogicPhantom, false, 0, true);

    G4Tubs *solidHole72 = new G4Tubs("SHole72", 0. * mm, 0.60165 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole72 =
        new G4LogicalVolume(solidHole72, hole_mat, "SHole72");
    G4VPhysicalVolume *PhysHole72 = new G4PVPlacement(
        0, G4ThreeVector(17.3358 * mm, -7.099 * mm, 0. * mm), LogicHole72,
        "PhysHole72", LogicPhantom, false, 0, true);

    G4Tubs *solidHole73 = new G4Tubs("SHole73", 0. * mm, 0.59427 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole73 =
        new G4LogicalVolume(solidHole73, hole_mat, "SHole73");
    G4VPhysicalVolume *PhysHole73 = new G4PVPlacement(
        0, G4ThreeVector(17.3341 * mm, -2.9479 * mm, 0. * mm), LogicHole73,
        "PhysHole73", LogicPhantom, false, 0, true);

    G4Tubs *solidHole74 = new G4Tubs("SHole74", 0. * mm, 0.58929 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole74 =
        new G4LogicalVolume(solidHole74, hole_mat, "SHole74");
    G4VPhysicalVolume *PhysHole74 = new G4PVPlacement(
        0, G4ThreeVector(18.5186 * mm, -9.1827 * mm, 0. * mm), LogicHole74,
        "PhysHole74", LogicPhantom, false, 0, true);

    G4Tubs *solidHole75 = new G4Tubs("SHole75", 0. * mm, 0.59179 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole75 =
        new G4LogicalVolume(solidHole75, hole_mat, "SHole75");
    G4VPhysicalVolume *PhysHole75 = new G4PVPlacement(
        0, G4ThreeVector(18.5197 * mm, -5.038 * mm, 0. * mm), LogicHole75,
        "PhysHole75", LogicPhantom, false, 0, true);

    G4Tubs *solidHole76 = new G4Tubs("SHole76", 0. * mm, 0.59179 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole76 =
        new G4LogicalVolume(solidHole76, hole_mat, "SHole76");
    G4VPhysicalVolume *PhysHole76 = new G4PVPlacement(
        0, G4ThreeVector(19.7277 * mm, -7.0992 * mm, 0. * mm), LogicHole76,
        "PhysHole76", LogicPhantom, false, 0, true);

    G4Tubs *solidHole77 = new G4Tubs("SHole77", 0. * mm, 0.59674 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole77 =
        new G4LogicalVolume(solidHole77, hole_mat, "SHole77");
    G4VPhysicalVolume *PhysHole77 = new G4PVPlacement(
        0, G4ThreeVector(19.733 * mm, -2.9482 * mm, 0. * mm), LogicHole77,
        "PhysHole77", LogicPhantom, false, 0, true);

    G4Tubs *solidHole78 = new G4Tubs("SHole78", 0. * mm, 0.58929 * mm,
                                     18.5 * mm, 0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole78 =
        new G4LogicalVolume(solidHole78, hole_mat, "SHole78");
    G4VPhysicalVolume *PhysHole78 = new G4PVPlacement(
        0, G4ThreeVector(20.9216 * mm, -5.0334 * mm, 0. * mm), LogicHole78,
        "PhysHole78", LogicPhantom, false, 0, true);

    G4Tubs *solidHole79 = new G4Tubs("SHole79", 0. * mm, 0.5992 * mm, 18.5 * mm,
                                     0.0 * deg, 360.0 * deg);
    G4LogicalVolume *LogicHole79 =
        new G4LogicalVolume(solidHole79, hole_mat, "SHole79");
    G4VPhysicalVolume *PhysHole79 = new G4PVPlacement(
        0, G4ThreeVector(22.1319 * mm, -2.9477 * mm, 0. * mm), LogicHole79,
        "PhysHole79", LogicPhantom, false, 0, true);

    if (isScintillator)
        ConstructScintillator();

    G4VisAttributes *cyan = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
    G4VisAttributes *magenta = new G4VisAttributes(G4Colour(1.0, 0, 1.0));
    G4VisAttributes *gray = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    G4VisAttributes *green = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
    cyan->SetVisibility(true);
    magenta->SetVisibility(true);
    // cyan->SetForceSolid(true);
    // magenta->SetForceSolid(true);
    gray->SetVisibility(true);
    green->SetVisibility(true);
    // gray->SetForceSolid(true);
    // green->SetForceSolid(true);

    logicScaterer->SetVisAttributes(cyan);
    logicScintillator->SetVisAttributes(magenta);
    LogicPhantom->SetVisAttributes(gray);

    return physWorld;
}

//.....oooOO0OOooo.........oooOO0OOooo.....SensitiveDetector:
// absorbor.....oooOO0OOooo.........oooOO0OOooo.......

void simDetectorConstruction::ConstructSDandField()
{
    simSensitiveDetector *sensDet1 =
        new simSensitiveDetector("SensitiveDetector");

    if (logicScintillator !=
        NULL) // detector define to any set up we can use it
        logicScintillator->SetSensitiveDetector(sensDet1);

    //.....oooOO0OOooo.........oooOO0OOooo.....SensitiveDetector:
    // Scaterer.....oooOO0OOooo.........oooOO0OOooo.......

    simSensitiveDetector *sensDet2 =
        new simSensitiveDetector("SensitiveDetector");

    if (logicScaterer != NULL)
        logicScaterer->SetSensitiveDetector(sensDet2);
}
