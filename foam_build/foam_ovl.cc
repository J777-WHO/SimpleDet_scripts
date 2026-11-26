// Read a GDML with Geant4 and run G4PVPlacement::CheckOverlaps on every world daughter whose name matches argv[2]
// (default "PassiveTrackerFoam"). Prints overlaps (warning only, no abort).
#include "G4GDMLParser.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
int main(int argc, char** argv) {
  G4GDMLParser p; p.SetOverlapCheck(false); p.Read(argv[1], false);
  G4String key = argc > 2 ? argv[2] : "PassiveTrackerFoam";
  int npts = argc > 3 ? atoi(argv[3]) : 10000;
  G4LogicalVolume* w = p.GetWorldVolume()->GetLogicalVolume();
  int n = 0, bad = 0;
  for (size_t i = 0; i < w->GetNoDaughters(); i++) {
    auto* pv = dynamic_cast<G4PVPlacement*>(w->GetDaughter(i));
    if (!pv || pv->GetName().find(key) == std::string::npos) continue;
    n++;
    bool o = pv->CheckOverlaps(npts, 0.01 * mm, true, 3);
    if (o) { bad++; std::cout << "### OVERLAP in " << pv->GetName() << std::endl; }
  }
  std::cout << "checked " << n << " volumes matching '" << key << "', with overlaps: " << bad << std::endl;
}
