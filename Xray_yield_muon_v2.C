// Xray_yield_muon_v2.C 批量遍历文件夹，多root文件累加汇总
// Function: Batch‑process multiple GAPS simulation ROOT files.
// Count muonic X‑ray photons produced by stopped negative muons in silicon,
// accumulate histograms and perform yield normalization, write aggregated output ROOT file.
// Input directory: contains multiple ecoMug_FTFP_BERT_HP_*.root simulation output files.
// Output: aggregated ROOT file including photon spectrum, yield spectra,
//         TTree with per‑qualified‑event information, and metadata records.
// Selection criteria: Negative muon (PDG = 13) stops inside silicon volume (volumeId starts with '2'),
//                     kinetic energy close to zero at stopping point.
// Photon selection: Only secondary photons directly produced by the stopped muon,
//                   filtered by matching photon parentId to muon trackId.
// Two normalization schemes:
//   hYieldDensity         : Normalized by total number of actual qualified stopped‑muon events
//                           Unit: photons / event / keV
//   hYieldDensityFixedDenom : Normalized by predefined fixed total event denominator
//                             Unit: photons / FIXED_DENOMINATOR / keV
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <TH1D.h>
#include <TAxis.h>
#include <TNamed.h>
#include <dirent.h>
#include <algorithm>
#include <cstdlib>

#include <gaps/CEventMc.hh>
#include <gaps/GDataEvent.hh>
#include <gaps/GGeometryTools.hh>
#include <gaps/GReco.hh>
using namespace std;

// 【修复】全局定义事件标识结构体，必须放在所有函数之前
struct EventIdentifier {
    Long64_t entry;
    Long64_t eventNumber;
    int trackId;
};

Long64_t extractRunId(const string& filename)
{
    size_t pos = filename.find("ecoMug_FTFP_BERT_HP_");
    if(pos == string::npos) return -1;
    pos += string("ecoMug_FTFP_BERT_HP_").size();
    size_t endpos = pos;
    while(endpos<filename.size() && isdigit(filename[endpos])){
        endpos++;
    }
    string numStr = filename.substr(pos, endpos-pos);
    return atoll(numStr.c_str());
}

vector<string> getRootFiles(const string& directory)
{
    vector<string> files;
    DIR* dir = opendir(directory.c_str());
    if (!dir) return files;

    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        string filename = entry->d_name;
        if (filename.find("ecoMug_FTFP_BERT_HP_") != string::npos &&
            filename.find(".root") != string::npos) {
            files.push_back(directory + "/" + filename);
        }
    }
    closedir(dir);

    sort(files.begin(), files.end(), [](const string& a, const string& b) {
        return extractRunId(a) < extractRunId(b);
    });

    return files;
}

// 处理单个输入文件，返回本文件停止事件列表
vector<EventIdentifier> processSingleFile_Pass1(const string& inFile, Long64_t &nTotalEntry_file, Long64_t &nHaveStopMu_file)
{
    vector<EventIdentifier> qualifiedEvents;
    nTotalEntry_file = 0;
    nHaveStopMu_file = 0;

    TFile *f = TFile::Open(inFile.c_str());
    if (!f || f->IsZombie()) {
        cerr<<"[WARN] Cannot open "<<inFile<<endl;
        if(f) f->Close();
        return qualifiedEvents;
    }
    TTree *tree = (TTree*)f->Get("TreeMc");
    if(!tree){
        cerr<<"[WARN] No TreeMc in "<<inFile<<endl;
        f->Close();
        return qualifiedEvents;
    }

    CEventMc *event = nullptr;
    tree->SetBranchAddress("Mc", &event);
    Long64_t nEntries = tree->GetEntries();
    nTotalEntry_file = nEntries;
    cout<<"\n==== Pass1 file:"<<inFile<<" entries="<<nEntries<<endl;

    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        CTrackMc* stoppingMuon = nullptr;
        std::vector<CTrackMc*> tracks = event->GetTracks();
        Long64_t currentEventNum = event->GetEventNumber();

        for (auto track : tracks) {
            if (track && track->GetParentId() == 0) {
                break;
            }
        }
        // select stopping muon  mu+ pdg= -13；mu‑ pdg=13
        for (auto track : tracks) {
            if (track && track->GetPdg() == 13 && track->GetParentId() == 0) {
                std::vector<unsigned int> vols = track->GetVolumeId();
                std::vector<double> kinE = track->GetKineticEnergy();
                if (vols.empty() || kinE.empty()) continue;

                double lastVol = vols.back();
                double lastE = kinE.back();
                std::string lastVolStr = std::to_string((int)lastVol);
                bool stopInSi = (lastVolStr[0] == '2');
                bool lowEnergy = (lastE < 1e-7);

                if (stopInSi && lowEnergy) {
                    stoppingMuon = track;
                    break;
                }
            }
        }

        if (stoppingMuon) {
            nHaveStopMu_file++;
            EventIdentifier id;
            id.entry = i;
            id.eventNumber = currentEventNum;
            id.trackId = stoppingMuon->GetTrackId();
            qualifiedEvents.push_back(id);
        }
        if ((i + 1) % 100000 == 0) {
            cout << "Progress: " << (i+1) << "/" << nEntries << " | Qualified in this file: " << nHaveStopMu_file << endl;
        }
    }
    f->Close();
    cout<<"[Pass1 done] this file qualified muon events = "<<nHaveStopMu_file<<endl;
    return qualifiedEvents;
}

void photon_save_with_yield_batch(const char* inputDir = "./",
                                  const char* outputRootFile = "xray_batch_sum.root")
{
    vector<string> rootFiles = getRootFiles(inputDir);
    if(rootFiles.empty()){
        cerr<<"No matched root files in directory "<<inputDir<<endl;
        return;
    }
    cout<<"Found "<<rootFiles.size()<<" root files:"<<endl;
    for(auto &f:rootFiles) cout<<"  "<<f<<" runId="<<extractRunId(f)<<endl;

    // 全局累积直方图，所有文件累加
    const int MAX_ENERGY = 120;
    const int NBINS = 120;
    const double binWidth = (double)MAX_ENERGY / NBINS;

    TH1D *hPhotonCounts = new TH1D("hPhotonCounts", ";Energy (keV);Counts", NBINS, 0, MAX_ENERGY);  // 次级光子直接计数
    TH1D *hYieldDensity = new TH1D("hYieldDensity", ";Energy (keV);Yield (photons/event/keV)", NBINS, 0, MAX_ENERGY);   // 次级光子对停止muon事件数归一
    TH1D *hYieldDensityFixedDenom = new TH1D("hYieldDensityFixedDenom", ";Energy (keV);Yield (photons/1.96M events/keV)", NBINS, 0, MAX_ENERGY);    // 次级光子对固定总事例归一
    TH1D *hPhotonsPerEvent = new TH1D("hPhotonsPerEvent", ";N_{photons};Number of Events", 30, -0.5, 29.5);     // 每个停止muon事件里面产生光子数

    TTree *eventInfoTree = new TTree("eventInfoTree", "Information for qualified events (batch sum)");
    Long64_t entry;
    Long64_t eventNumber;
    int muonTrackId;
    int primaryPdg;
    int nPhotons;
    double totalPhotonEnergy;
    double avgPhotonEnergy;
    Long64_t runId;

    eventInfoTree->Branch("runId",&runId,"runId/L");
    eventInfoTree->Branch("entry", &entry, "entry/L");
    eventInfoTree->Branch("eventNumber", &eventNumber, "eventNumber/L");
    eventInfoTree->Branch("muonTrackId", &muonTrackId, "muonTrackId/I");
    eventInfoTree->Branch("primaryPdg", &primaryPdg, "primaryPdg/I");
    eventInfoTree->Branch("nPhotons", &nPhotons, "nPhotons/I");
    eventInfoTree->Branch("totalPhotonEnergy_keV", &totalPhotonEnergy, "totalPhotonEnergy_keV/D");
    eventInfoTree->Branch("avgPhotonEnergy_keV", &avgPhotonEnergy, "avgPhotonEnergy_keV/D");

    Long64_t totalQualifiedAll = 0;     // 文件中的所有停止的事件数
    Long64_t totalEntriesAll = 0;       // 文件中的所有事件数
    const double FIXED_DENOMINATOR = 382 * 10000.0;

    // 循环每个root文件
    for(const string &fname : rootFiles)
    {
        Long64_t run = extractRunId(fname);
        cout<<"\n>>>>>>>>>> Process file: "<<fname<<" runId="<<run<<" <<<<<<<<<<"<<endl;

        Long64_t nEntry_file=0, nQual_file=0;
        vector<EventIdentifier> qualList = processSingleFile_Pass1(fname, nEntry_file, nQual_file);
        totalEntriesAll += nEntry_file;
        totalQualifiedAll += nQual_file;
        if(qualList.empty()){
            cout<<"No qualified events for this file, skip pass2."<<endl;
            continue;
        }

        //==== Pass2 当前文件重新打开，填充直方图与Tree
        TFile *f = TFile::Open(fname.c_str());
        if (!f || f->IsZombie()) {
            cerr<<"Cannot reopen "<<fname<<" skip pass2"<<endl;
            if(f)f->Close();
            continue;
        }
        TTree *tree = (TTree*)f->Get("TreeMc");
        CEventMc *event = nullptr;
        tree->SetBranchAddress("Mc", &event);

        Long64_t processed = 0;
        for(auto &id : qualList)
        {
            tree->GetEntry(id.entry);
            runId = run;
            entry = id.entry;
            eventNumber = id.eventNumber;
            muonTrackId = id.trackId;

            std::vector<CTrackMc*> tracks = event->GetTracks();
            CTrackMc* primaryParticle = nullptr;
            for (auto track : tracks) {
                if (track && track->GetParentId() == 0) {
                    primaryParticle = track;
                    break;
                }
            }
            primaryPdg = primaryParticle ? primaryParticle->GetPdg() : 0;

            nPhotons = 0;
            totalPhotonEnergy = 0.0;
            for (auto track : tracks) {
                if (!track) continue;
                if (track->GetPdg() != 22) continue;
                if (track->GetParentId() != muonTrackId) continue;

                double energy_keV = track->GetVertexKineticEnergy() * 1000.0;
                hPhotonCounts->Fill(energy_keV);
                hYieldDensity->Fill(energy_keV);
                hYieldDensityFixedDenom->Fill(energy_keV);
                nPhotons++;
                totalPhotonEnergy += energy_keV;
            }
            avgPhotonEnergy = (nPhotons>0) ? totalPhotonEnergy / nPhotons : 0.0;
            hPhotonsPerEvent->Fill(nPhotons);
            eventInfoTree->Fill();

            processed++;
            if (processed % 10000 == 0) {
                cout << "File processed: " << processed << "/" << qualList.size()
                     << " | Photons so far(total): " << hPhotonCounts->GetEntries() << endl;
            }
        }
        f->Close();
        cout<<"File done, this file qualified:"<<qualList.size()<<endl;
    }

    cout<<"\n======== All files finished ========="<<endl;
    cout<<"Total entries over all files: "<<totalEntriesAll<<endl;
    cout<<"Total qualified muon events(sum): "<<totalQualifiedAll<<endl;
    if(totalQualifiedAll <=0){
        cerr<<"Zero qualified events in all files, exit."<<endl;
        delete hPhotonCounts;
        delete hYieldDensity;
        delete hYieldDensityFixedDenom;
        delete hPhotonsPerEvent;
        delete eventInfoTree;
        return;
    }

    //===== 全局归一化
    double scaleFactorQualified = 1.0 / (totalQualifiedAll * binWidth);
    hYieldDensity->Scale(scaleFactorQualified);

    double scaleFactorFixed = 1.0 / (FIXED_DENOMINATOR * binWidth);
    hYieldDensityFixedDenom->Scale(scaleFactorFixed);

    //===== 写出汇总ROOT
    TFile *outRoot = TFile::Open(outputRootFile, "RECREATE");
    hPhotonCounts->Write();
    hYieldDensity->Write();
    hYieldDensityFixedDenom->Write();
    hPhotonsPerEvent->Write();
    eventInfoTree->Write();

    TNamed *nEventsInfo = new TNamed("QualifiedMuonEvents", Form("%lld", totalQualifiedAll));
    nEventsInfo->Write();
    TNamed *fixedDenomInfo = new TNamed("FixedDenominator", Form("%.0f", FIXED_DENOMINATOR));
    fixedDenomInfo->Write();

    TH1D *hPrimaryPdg = new TH1D("hPrimaryPdg", "Primary Particle PDG Distribution;PDG Code;Counts", 50, -50, 50);
    eventInfoTree->Draw("primaryPdg>>hPrimaryPdg", "", "goff");
    hPrimaryPdg->Write();

    outRoot->Close();

    cout<<"\nAll done! Output: "<<outputRootFile<<endl;
    cout<<"  hYieldDensity normalized by total qualified = "<<totalQualifiedAll<<endl;
    cout<<"  hYieldDensityFixedDenom normalized by FIXED_DENOMINATOR = "<<FIXED_DENOMINATOR<<endl;
}

int main(int argc, char** argv)
{
    if(argc!=3){
        fprintf(stderr,"Usage: Xray_yield_muon_v2 input_dir output.root\n");
        fprintf(stderr,"Example: ./Xray_yield_muon_v2 /data/xxx/simu_0817_NY/ ../out_08/out.root\n");
        return 1;
    }
    photon_save_with_yield_batch(argv[1], argv[2]);
    return 0;
}