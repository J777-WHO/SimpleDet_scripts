// check_shape_overlap.C: strip shape bounding boxes (local, before node rotation) + ROOT overlap check
void check_shape_overlap(const char* fn){
  TGeoManager::LockDefaultUnits(false); TGeoManager::SetDefaultUnits(TGeoManager::kG4Units); TGeoManager::LockDefaultUnits(true);
  TGeoManager* g = TGeoManager::Import(fn); gErrorIgnoreLevel = kWarning;
  TGeoNode* w = g->GetTopNode(); int shown = 0;
  for (int i = 0; i < w->GetNdaughters() && shown < 16; i++) {
    TGeoNode* nd = w->GetDaughter(i); TString nm = nd->GetName();
    if (!nm.BeginsWith("ActiveSiLiDetCylStripVolId20000000") && !nm.BeginsWith("ActiveSiLiDetCylStripVolId20100000")) continue;
    auto* bb = (TGeoBBox*)nd->GetVolume()->GetShape(); const double* o = bb->GetOrigin();
    printf("%-36.36s LV %-18.18s bbox local centre (%6.2f,%6.2f) half (%5.2f,%5.2f) | vol %.1f mm3\n", nm.Data(), nd->GetVolume()->GetName(), o[0], o[1], bb->GetDX(), bb->GetDY(), nd->GetVolume()->GetShape()->Capacity());
    shown++;
  }
  g->CheckOverlaps(0.01);   // 0.01 = 10 um in G4 units (mm)
  TObjArray* ov = g->GetListOfOverlaps(); int nstrip = 0;
  for (int i = 0; i < ov->GetEntriesFast(); i++) { TString s = ((TNamed*)ov->At(i))->GetTitle(); if (s.Contains("Strip")) { nstrip++; if (nstrip <= 5) printf("OVL %s\n", s.Data()); } }
  printf("overlaps total %d, involving strips %d\n", ov->GetEntriesFast(), nstrip);
}
