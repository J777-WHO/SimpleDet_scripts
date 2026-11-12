// singleEvent.C
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

using namespace std;

void singleEvent(const char* inputFile, Long64_t eventId, const char* outputFile = "") {
    
    TFile *f = TFile::Open(inputFile);
    if (!f || f->IsZombie()) {
        cerr << "Error: Cannot open file " << inputFile << endl;
        return;
    }
    
    TTree *tree = (TTree*)f->Get("TreeMc");
    if (!tree) {
        cerr << "Error: Cannot find TreeMc" << endl;
        f->Close();
        return;
    }
    
    CEventMc *event = 0;
    tree->SetBranchAddress("Mc", &event);
    
    Long64_t nEntries = tree->GetEntries();
    if (eventId >= nEntries) {
        cerr << "Error: EventId " << eventId << " out of range (0-" << nEntries-1 << ")" << endl;
        f->Close();
        return;
    }
    
    tree->GetEntry(eventId);
    
    ostream* out = &cout;
    ofstream outFile;
    if (outputFile && outputFile[0] != '\0') {
        outFile.open(outputFile);
        if (outFile.is_open()) out = &outFile;
    }
    
    out->setf(ios::fixed);
    out->precision(8);
    
    vector<CTrackMc*> tracks = event->GetTracks();
    CTrackMc* primary = nullptr;
    for (auto track : tracks) {
        if (track && track->GetParentId() == 0) {
            primary = track;
            break;
        }
    }
    
    *out << "========================================" << endl;
    *out << "Event ID: " << eventId << endl;
    *out << "========================================" << endl;
    
    int stopVol = event->GetPrimaryStoppingVolume();
    TVector3 stopPos = event->GetPrimaryStoppingPosition();
    
    *out << "PrimaryPDG: " << (primary ? primary->GetPdg() : 0) << endl;
    *out << "PrimaryTrackID: " << (primary ? primary->GetTrackId() : 0) << endl;
    *out << "PrimaryParentID: " << (primary ? primary->GetParentId() : 0) << endl;
    *out << "PrimaryProcessType: " << (primary ? (int)primary->GetProcessType() : 0) << endl;
    *out << "PrimaryBeta: " << event->GetPrimaryBeta() << endl;
    *out << "PrimaryBetaGen: " << event->GetPrimaryBetaGenerated() << endl;
    *out << "PrimaryKineticEnergy: " << event->GetPrimaryKineticEnergy() << " MeV" << endl;
    *out << "PrimaryKineticEnergyGenerated: " << event->GetPrimaryKineticEnergyGenerated() << " MeV" << endl;
    *out << "StoppingVolume: " << stopVol << endl;
    *out << "StoppingPosition(mm): (" << stopPos.x() << "," << stopPos.y() << "," << stopPos.z() << ")" << endl;
    
    // ====================== 表格1: 初级粒子（全部动能输出） ======================
    if (primary) {
        *out << "\n========== Table 1: Primary Particle (All Kinetic Energy) ==========" << endl;
        
        vector<double> preKinE = primary->GetPreKineticEnergy();
        vector<double> kinE = primary->GetKineticEnergy();
        vector<unsigned int> vols = primary->GetVolumeId();
        vector<TVector3> positions = primary->GetPosition();
        vector<TVector3> momenta = primary->GetMomentumDirection();
        double vertexE = primary->GetVertexKineticEnergy();
        
        *out << "Vertex (诞生) Kinetic Energy: " << vertexE << " MeV" << endl;
        
        *out << "\n+----+----------+-------------------+-------------------+-------------------+-----------+-----------+-----------+-----------+-----------+" << endl;
        *out << "|Step| VolumeId | VertexE           | PreKinE           | PostKinE          |    Px     |    Py     |    Pz     |    X      |    Y      |    Z      |" << endl;
        *out << "+----+----------+-------------------+-------------------+-------------------+-----------+-----------+-----------+-----------+-----------+-----------+" << endl;
        
        for (size_t step = 0; step < vols.size(); step++) {
            *out << "|" << setw(3) << step << " |" << setw(8) << vols[step] << "|";
            
            *out << setw(18) << vertexE << "|";
            if (step < preKinE.size()) *out << setw(18) << preKinE[step] << "|";
            else *out << setw(18) << "N/A |";
            if (step < kinE.size()) *out << setw(18) << kinE[step] << "|";
            else *out << setw(18) << "N/A |";
            
            if (step < momenta.size()) {
                *out << setw(11) << momenta[step].x() << "|"
                     << setw(11) << momenta[step].y() << "|"
                     << setw(11) << momenta[step].z() << "|";
            } else {
                *out << setw(11) << "N/A |" << setw(11) << "N/A |" << setw(11) << "N/A |";
            }
            
            if (step < positions.size()) {
                *out << setw(11) << positions[step].x() << "|"
                     << setw(11) << positions[step].y() << "|"
                     << setw(11) << positions[step].z() << "|";
            } else {
                *out << setw(11) << "N/A |" << setw(11) << "N/A |" << setw(11) << "N/A |";
            }
            *out << endl;
        }
        *out << "+----+----------+-------------------+-------------------+-------------------+-----------+-----------+-----------+-----------+-----------+-----------+" << endl;
    }
    
    // ====================== 表格2: All Tracks Summary ======================
    *out << "\n========== Table 2: All Tracks Summary ==========" << endl;
    *out << "+-----+--------+--------+-------------+-------------------+-------------------+-------------------+-------------------+----------+" << endl;
    *out << "|Index|  PDG   |ParentId| ProcessType | VertexE(MeV)      |    VertexVolId    | FirstPreE(MeV)    | LastPostE(MeV)    | LastVol  |" << endl;
    *out << "+-----+--------+--------+-------------+-------------------+-------------------+-------------------+-------------------+----------+" << endl;
    
    for (size_t i = 0; i < tracks.size(); i++) {
        CTrackMc* track = tracks[i];
        if (!track) continue;
        
        double vertexE = track->GetVertexKineticEnergy();
        unsigned int vertexVolId = track->GetVertexVolumeId();
        vector<double> preE = track->GetPreKineticEnergy();
        vector<double> postE = track->GetKineticEnergy();
        vector<unsigned int> vols = track->GetVolumeId();
        
        *out << "|" << setw(4) << i << " |" << setw(6) << track->GetPdg() << " |"
             << setw(6) << track->GetParentId() << " |" << setw(11) << (int)track->GetProcessType() << " |";
        
        *out << setw(18) << vertexE << "|";
        *out << setw(18) << vertexVolId << "|"; 
        if (!preE.empty()) *out << setw(18) << preE.front() << "|";
        else *out << setw(18) << "N/A |";
        if (!postE.empty()) *out << setw(18) << postE.back() << "|";
        else *out << setw(18) << "N/A |";
        
        if (!vols.empty()) *out << setw(8) << vols.back() << " |";
        else *out << setw(8) << "N/A |";
        *out << endl;
    }
    *out << "+-----+--------+--------+-------------+-------------------+-------------------+-------------------+----------+" << endl;
    
    // ====================== 表格3: 湮灭光子（逐步骤输出全部三种动能） ======================
    *out << "\n========== Table 3: Annihilation Photons (ProcType=40 & PDG=22) ==========" << endl;
    
    int photonCount = 0;
    for (auto track : tracks) {
        if (!track) continue;
        if (track->GetParentId() == 0) continue;
        if ((int)track->GetProcessType() != 40) continue;
        if (track->GetPdg() != 22) continue;
        
        photonCount++;
        vector<double> preE = track->GetPreKineticEnergy();
        vector<double> postE = track->GetKineticEnergy();
        vector<unsigned int> vols = track->GetVolumeId();
        double vertexE = track->GetVertexKineticEnergy();
        
        *out << "\n--- Photon " << photonCount << " | ParentID=" << track->GetParentId() << " ---" << endl;
        *out << "+----+----------+-------------------+-------------------+-------------------+" << endl;
        *out << "|Step| VolumeId | VertexE           | PreKinE           | PostKinE          |" << endl;
        *out << "+----+----------+-------------------+-------------------+-------------------+" << endl;
        
        for (size_t s = 0; s < vols.size(); s++) {
            *out << "|" << setw(3) << s << " |" << setw(8) << vols[s] << "|";
            *out << setw(18) << vertexE << "|";
            if (s < preE.size()) *out << setw(18) << preE[s] << "|";
            else *out << setw(18) << "N/A |";
            if (s < postE.size()) *out << setw(18) << postE[s] << "|";
            else *out << setw(18) << "N/A |";
            *out << endl;
        }
        *out << "+----+----------+-------------------+-------------------+-------------------+" << endl;
    }
    
    if (photonCount == 0) *out << "  No annihilation photons found." << endl;
    
    *out << "\n========================================" << endl;
    *out << "End of Event" << endl;
    *out << "========================================" << endl;
    
    if (outFile.is_open()) outFile.close();
    f->Close();
}