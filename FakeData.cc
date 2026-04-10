#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TVector3.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"
#include "TEfficiency.h"
#include "TMVA/Reader.h"

#include <vector>

#include "Helpers.cc"
#include "FakeDataHelper.h"

void FakeData(int sample = 1) {
    // Set defaults
    gStyle->SetOptStat(0); // get rid of stats box
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    TH1::AddDirectory(false);
    gStyle->SetPalette(kRainBow);

    TGraph* gProton = new TGraph();
    TGraph* gPion   = new TGraph();
    TGraph* gMuonTG = new TGraph();

    // Initialize points
    initializeProtonPoints(gProton);
    initializePionPoints(gPion);
    initializeMuonNoBraggPoints(gMuonTG);

    int FontStyle = 132;
    double TextSize = 0.06;
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis/figs/FakeData/Sample" + TString::Itoa(sample, 10) + "/";
    std::cout << "Generating Sample: " << sample << std::endl;

    // Load file with data products
    TString RootFilePath = "/exp/lariat/app/users/epelaez/files/anatree_60a/histo.root";
    std::unique_ptr<TFile> File(TFile::Open(RootFilePath));
    TDirectory* Directory = (TDirectory*)File->Get("anatree");

    // Fake data configuration
    FakeDataFD::Scenario FDScenario;
    if (sample == 1) {
        FDScenario = FakeDataFD::Scenario::FDa;
    } else if (sample == 2) {
        FDScenario = FakeDataFD::Scenario::FDb;
    } else if (sample == 3) {
        FDScenario = FakeDataFD::Scenario::FDc;
    } else if (sample == 4) {
        FDScenario = FakeDataFD::Scenario::FDd;
    }

    ///////////////////
    // Load branches //
    ///////////////////

    // Load tree and branches
    TTree* tree = (TTree*) Directory->Get<TTree>("anatree");

    // Run information
    int run, subrun, event; bool isData = false;
    tree->SetBranchAddress("run", &run); 
    tree->SetBranchAddress("subrun", &subrun); 
    tree->SetBranchAddress("event", &event);

    // Track information
    const int kMaxTrack = 100; // 10000
    int   ntracks_reco;                       tree->SetBranchAddress("ntracks_reco",    &ntracks_reco);
    float trkvtxx[kMaxTrack];                 tree->SetBranchAddress("trkvtxx",          &trkvtxx);
    float trkvtxy[kMaxTrack];                 tree->SetBranchAddress("trkvtxy",          &trkvtxy);
    float trkvtxz[kMaxTrack];                 tree->SetBranchAddress("trkvtxz",          &trkvtxz);
    float trkendx[kMaxTrack];                 tree->SetBranchAddress("trkendx",          &trkendx);
    float trkendy[kMaxTrack];                 tree->SetBranchAddress("trkendy",          &trkendy);
    float trkendz[kMaxTrack];                 tree->SetBranchAddress("trkendz",          &trkendz);
    // float trkstartdcosx[kMaxTrack];           tree->SetBranchAddress("trkstartdcosx",    &trkstartdcosx);
    // float trkstartdcosy[kMaxTrack];           tree->SetBranchAddress("trkstartdcosy",    &trkstartdcosy);
    // float trkstartdcosz[kMaxTrack];           tree->SetBranchAddress("trkstartdcosz",    &trkstartdcosz);
    // float trkenddcosx[kMaxTrack];             tree->SetBranchAddress("trkenddcosx",      &trkenddcosx);
    // float trkenddcosy[kMaxTrack];             tree->SetBranchAddress("trkenddcosy",      &trkenddcosy);
    // float trkenddcosz[kMaxTrack];             tree->SetBranchAddress("trkenddcosz",      &trkenddcosz);
    // float trklength[kMaxTrack];               tree->SetBranchAddress("trklength",        &trklength);
    int   trkWCtoTPCMatch[kMaxTrack];         tree->SetBranchAddress("trkWCtoTPCMatch",  &trkWCtoTPCMatch);

    // Wire-chamber track information
    const int kMaxWCTracks = 1; // 1000
    int   nwctrks;                              tree->SetBranchAddress("nwctrks",           &nwctrks);
    // float wctrk_XFaceCoor[kMaxWCTracks];        tree->SetBranchAddress("wctrk_XFaceCoor",   &wctrk_XFaceCoor);
    // float wctrk_YFaceCoor[kMaxWCTracks];        tree->SetBranchAddress("wctrk_YFaceCoor",   &wctrk_YFaceCoor);
    float wctrk_momentum[kMaxWCTracks];         tree->SetBranchAddress("wctrk_momentum",    &wctrk_momentum);
    // float wctrk_Px[kMaxWCTracks];               tree->SetBranchAddress("wctrk_Px",          &wctrk_Px);
    // float wctrk_Py[kMaxWCTracks];               tree->SetBranchAddress("wctrk_Py",          &wctrk_Py);
    // float wctrk_Pz[kMaxWCTracks];               tree->SetBranchAddress("wctrk_Pz",          &wctrk_Pz);
    float wctrk_theta[kMaxWCTracks];            tree->SetBranchAddress("wctrk_theta",        &wctrk_theta);
    float wctrk_phi[kMaxWCTracks];              tree->SetBranchAddress("wctrk_phi",          &wctrk_phi);
    // float wctrk_residual[kMaxWCTracks];         tree->SetBranchAddress("wctrk_residual",     &wctrk_residual);
    // int   wctrk_wcmissed[kMaxWCTracks];         tree->SetBranchAddress("wctrk_wcmissed",     &wctrk_wcmissed);
    int   wctrk_picky[kMaxWCTracks];            tree->SetBranchAddress("wctrk_picky",        &wctrk_picky);
    // float wctrk_XDist[kMaxWCTracks];            tree->SetBranchAddress("wctrk_XDist",        &wctrk_XDist);
    // float wctrk_YDist[kMaxWCTracks];            tree->SetBranchAddress("wctrk_YDist",        &wctrk_YDist);
    // float wctrk_ZDist[kMaxWCTracks];            tree->SetBranchAddress("wctrk_ZDist",        &wctrk_ZDist);
    // float wctrk_YKink[kMaxWCTracks];            tree->SetBranchAddress("wctrk_YKink",        &wctrk_YKink);
    // int   wctrk_WC1XMult[kMaxWCTracks];         tree->SetBranchAddress("wctrk_WC1XMult",     &wctrk_WC1XMult);
    // int   wctrk_WC1YMult[kMaxWCTracks];         tree->SetBranchAddress("wctrk_WC1YMult",     &wctrk_WC1YMult);
    // int   wctrk_WC2XMult[kMaxWCTracks];         tree->SetBranchAddress("wctrk_WC2XMult",     &wctrk_WC2XMult);
    // int   wctrk_WC2YMult[kMaxWCTracks];         tree->SetBranchAddress("wctrk_WC2YMult",     &wctrk_WC2YMult);
    // int   wctrk_WC3XMult[kMaxWCTracks];         tree->SetBranchAddress("wctrk_WC3XMult",     &wctrk_WC3XMult);
    // int   wctrk_WC3YMult[kMaxWCTracks];         tree->SetBranchAddress("wctrk_WC3YMult",     &wctrk_WC3YMult);
    // int   wctrk_WC4XMult[kMaxWCTracks];         tree->SetBranchAddress("wctrk_WC4XMult",     &wctrk_WC4XMult);
    // int   wctrk_WC4YMult[kMaxWCTracks];         tree->SetBranchAddress("wctrk_WC4YMult",     &wctrk_WC4YMult);
    // float XWireHist[kMaxWCTracks][1000];         tree->SetBranchAddress("XWireHist",          &XWireHist);
    // float YWireHist[kMaxWCTracks][1000];         tree->SetBranchAddress("YWireHist",          &YWireHist);
    // float XAxisHist[kMaxWCTracks][1000];         tree->SetBranchAddress("XAxisHist",          &XAxisHist);
    // float YAxisHist[kMaxWCTracks][1000];         tree->SetBranchAddress("YAxisHist",          &YAxisHist);
    // float WC1xPos[kMaxWCTracks];                 tree->SetBranchAddress("WC1xPos",            &WC1xPos);
    // float WC1yPos[kMaxWCTracks];                 tree->SetBranchAddress("WC1yPos",            &WC1yPos);
    // float WC1zPos[kMaxWCTracks];                 tree->SetBranchAddress("WC1zPos",            &WC1zPos);
    // float WC2xPos[kMaxWCTracks];                 tree->SetBranchAddress("WC2xPos",            &WC2xPos);
    // float WC2yPos[kMaxWCTracks];                 tree->SetBranchAddress("WC2yPos",            &WC2yPos);
    // float WC2zPos[kMaxWCTracks];                 tree->SetBranchAddress("WC2zPos",            &WC2zPos);
    // float WC3xPos[kMaxWCTracks];                 tree->SetBranchAddress("WC3xPos",            &WC3xPos);
    // float WC3yPos[kMaxWCTracks];                 tree->SetBranchAddress("WC3yPos",            &WC3yPos);
    // float WC3zPos[kMaxWCTracks];                 tree->SetBranchAddress("WC3zPos",            &WC3zPos);
    float WC4xPos[kMaxWCTracks];                 tree->SetBranchAddress("WC4xPos",            &WC4xPos);
    // float WC4yPos[kMaxWCTracks];                 tree->SetBranchAddress("WC4yPos",            &WC4yPos);
    // float WC4zPos[kMaxWCTracks];                 tree->SetBranchAddress("WC4zPos",            &WC4zPos);

    // Calorimetry information
    const int kMaxTrackHits = 1000; // 1000
    int   ntrkcalopts[kMaxTrack][2];                    tree->SetBranchAddress("ntrkcalopts", &ntrkcalopts);
    // float trkpida[kMaxTrack][2];                        tree->SetBranchAddress("trkpida",     &trkpida);
    // float trkke[kMaxTrack][2];                          tree->SetBranchAddress("trkke",       &trkke);
    float trkdedx[kMaxTrack][2][kMaxTrackHits];         tree->SetBranchAddress("trkdedx",     &trkdedx);
    // float trkdqdx[kMaxTrack][2][kMaxTrackHits];         tree->SetBranchAddress("trkdqdx",     &trkdqdx);
    float trkrr[kMaxTrack][2][kMaxTrackHits];           tree->SetBranchAddress("trkrr",       &trkrr);
    float trkpitch[kMaxTrack][2][kMaxTrackHits];        tree->SetBranchAddress("trkpitch",    &trkpitch);
    float trkxyz[kMaxTrack][2][kMaxTrackHits][3];       tree->SetBranchAddress("trkxyz",      &trkxyz);

    // Trajectory information for tracks
    const int kMaxTrajHits = 1000; // 1000
    int   nTrajPoint[kMaxTrack];                  tree->SetBranchAddress("nTrajPoint", &nTrajPoint);
    // float pHat0_X[kMaxTrack][kMaxTrajHits];       tree->SetBranchAddress("pHat0_X",    &pHat0_X);
    // float pHat0_Y[kMaxTrack][kMaxTrajHits];       tree->SetBranchAddress("pHat0_Y",    &pHat0_Y);
    // float pHat0_Z[kMaxTrack][kMaxTrajHits];       tree->SetBranchAddress("pHat0_Z",    &pHat0_Z);
    float trjPt_X[kMaxTrack][kMaxTrajHits];       tree->SetBranchAddress("trjPt_X",    &trjPt_X);
    float trjPt_Y[kMaxTrack][kMaxTrajHits];       tree->SetBranchAddress("trjPt_Y",    &trjPt_Y);
    float trjPt_Z[kMaxTrack][kMaxTrajHits];       tree->SetBranchAddress("trjPt_Z",    &trjPt_Z);

    // Geant4 information for truth tracks
    const int kMaxPrimaries      = 4000; // 20000
    const int kMaxPrimaryPart    = 60;   // 50
    const int kMaxTruePrimaryPts = 1000; // 5000
    int   no_primaries;                                   tree->SetBranchAddress("no_primaries",        &no_primaries);
    int   geant_list_size;                                tree->SetBranchAddress("geant_list_size",      &geant_list_size);
    int   pdg[kMaxPrimaries];                             tree->SetBranchAddress("pdg",                  &pdg);
    float Mass[kMaxPrimaries];                            tree->SetBranchAddress("Mass",                 &Mass);
    // float StartPointx[kMaxPrimaries];                     tree->SetBranchAddress("StartPointx",          &StartPointx);
    // float StartPointy[kMaxPrimaries];                     tree->SetBranchAddress("StartPointy",          &StartPointy);
    float StartPointz[kMaxPrimaries];                     tree->SetBranchAddress("StartPointz",          &StartPointz);
    float Eng[kMaxPrimaries];                             tree->SetBranchAddress("Eng",                  &Eng);
    float Px[kMaxPrimaries];                              tree->SetBranchAddress("Px",                   &Px);
    float Py[kMaxPrimaries];                              tree->SetBranchAddress("Py",                   &Py);
    float Pz[kMaxPrimaries];                              tree->SetBranchAddress("Pz",                   &Pz);
    float EndPointx[kMaxPrimaries];                       tree->SetBranchAddress("EndPointx",            &EndPointx);
    float EndPointy[kMaxPrimaries];                       tree->SetBranchAddress("EndPointy",            &EndPointy);
    float EndPointz[kMaxPrimaries];                       tree->SetBranchAddress("EndPointz",            &EndPointz);
    float EndEng[kMaxPrimaries];                          tree->SetBranchAddress("EndEng",               &EndEng);
    float EndPx[kMaxPrimaries];                           tree->SetBranchAddress("EndPx",                &EndPx);
    float EndPy[kMaxPrimaries];                           tree->SetBranchAddress("EndPy",                &EndPy);
    float EndPz[kMaxPrimaries];                           tree->SetBranchAddress("EndPz",                &EndPz);
    // float StartT[kMaxPrimaries];                          tree->SetBranchAddress("StartT",               &StartT);
    // float EndT[kMaxPrimaries];                            tree->SetBranchAddress("EndT",                 &EndT);
    // float PathLenInTpcAV[kMaxPrimaries];                  tree->SetBranchAddress("PathLenInTpcAV",       &PathLenInTpcAV);
    // bool  StartInTpcAV[kMaxPrimaries];                    tree->SetBranchAddress("StartInTpcAV",         &StartInTpcAV);
    // bool  EndInTpcAV[kMaxPrimaries];                      tree->SetBranchAddress("EndInTpcAV",           &EndInTpcAV);
    int   Process[kMaxPrimaries];                         tree->SetBranchAddress("Process",              &Process);
    int   NumberDaughters[kMaxPrimaries];                 tree->SetBranchAddress("NumberDaughters",      &NumberDaughters);
    int   TrackId[kMaxPrimaries];                         tree->SetBranchAddress("TrackId",              &TrackId);
    int   Mother[kMaxPrimaries];                          tree->SetBranchAddress("Mother",               &Mother);
    int   process_primary[kMaxPrimaries];                 tree->SetBranchAddress("process_primary",      &process_primary);
    // std::vector<std::string> G4Process;                   tree->SetBranchAddress("G4Process",            &G4Process);
    // std::vector<std::string> G4FinalProcess;              tree->SetBranchAddress("G4FinalProcess",       &G4FinalProcess);
    std::vector<int>*         InteractionPoint = nullptr;            tree->SetBranchAddress("InteractionPoint",     &InteractionPoint);
    std::vector<int>*         InteractionPointType = nullptr;        tree->SetBranchAddress("InteractionPointType", &InteractionPointType);
    int   NTrTrajPts[kMaxPrimaryPart];                    tree->SetBranchAddress("NTrTrajPts",           &NTrTrajPts);
    float MidPosX[kMaxPrimaryPart][kMaxTruePrimaryPts];   tree->SetBranchAddress("MidPosX",              &MidPosX);
    float MidPosY[kMaxPrimaryPart][kMaxTruePrimaryPts];   tree->SetBranchAddress("MidPosY",              &MidPosY);
    float MidPosZ[kMaxPrimaryPart][kMaxTruePrimaryPts];   tree->SetBranchAddress("MidPosZ",              &MidPosZ);
    float MidPx[kMaxPrimaryPart][kMaxTruePrimaryPts];     tree->SetBranchAddress("MidPx",                &MidPx);
    float MidPy[kMaxPrimaryPart][kMaxTruePrimaryPts];     tree->SetBranchAddress("MidPy",                &MidPy);
    float MidPz[kMaxPrimaryPart][kMaxTruePrimaryPts];     tree->SetBranchAddress("MidPz",                &MidPz);

    // Information about wire plane hits
    const int kMaxHits = 1000;
    int    nhits;                              tree->SetBranchAddress("nhits",              &nhits);
    int    hit_plane[kMaxHits];                tree->SetBranchAddress("hit_plane",           hit_plane);
    // int    hit_wire[kMaxHits];                 tree->SetBranchAddress("hit_wire",            hit_wire);
    int    hit_channel[kMaxHits];              tree->SetBranchAddress("hit_channel",         hit_channel);
    int    hit_trkid[kMaxHits];               tree->SetBranchAddress("hit_trkid",           hit_trkid);
    // float  hit_peakT[kMaxHits];               tree->SetBranchAddress("hit_peakT",           hit_peakT);
    // float  hit_charge[kMaxHits];              tree->SetBranchAddress("hit_charge",           hit_charge);
    // float  hit_electrons[kMaxHits];           tree->SetBranchAddress("hit_electrons",        hit_electrons);
    // float  hit_ph[kMaxHits];                  tree->SetBranchAddress("hit_ph",               hit_ph);
    // float  hit_rms[kMaxHits];                 tree->SetBranchAddress("hit_rms",              hit_rms);
    // float  hit_tstart[kMaxHits];              tree->SetBranchAddress("hit_tstart",           hit_tstart);
    // float  hit_tend[kMaxHits];                tree->SetBranchAddress("hit_tend",             hit_tend);
    float  hit_driftT[kMaxHits];              tree->SetBranchAddress("hit_driftT",           hit_driftT);
    // float  hit_dQds[kMaxHits];               tree->SetBranchAddress("hit_dQds",             hit_dQds);
    // float  hit_dEds[kMaxHits];               tree->SetBranchAddress("hit_dEds",             hit_dEds);
    // float  hit_ds[kMaxHits];                 tree->SetBranchAddress("hit_ds",               hit_ds);
    // float  hit_resrange[kMaxHits];            tree->SetBranchAddress("hit_resrange",         hit_resrange);
    float  hit_x[kMaxHits];                  tree->SetBranchAddress("hit_x",                hit_x);
    float  hit_y[kMaxHits];                  tree->SetBranchAddress("hit_y",                hit_y);
    float  hit_z[kMaxHits];                  tree->SetBranchAddress("hit_z",                hit_z);
    // int    hit_g4id[kMaxHits];               tree->SetBranchAddress("hit_g4id",             hit_g4id);
    // float  hit_g4frac[kMaxHits];             tree->SetBranchAddress("hit_g4frac",           hit_g4frac);
    // float  hit_g4nelec[kMaxHits];            tree->SetBranchAddress("hit_g4nelec",          hit_g4nelec);
    // float  hit_g4energy[kMaxHits];           tree->SetBranchAddress("hit_g4energy",         hit_g4energy);

    // Simchannel information
    const int kMaxIDE = 5000;
    int maxTrackIDE;          tree->SetBranchAddress("maxTrackIDE", &maxTrackIDE);
    float IDEEnergy[kMaxIDE]; tree->SetBranchAddress("IDEEnergy",   &IDEEnergy);
    float IDEPos[kMaxIDE][3]; tree->SetBranchAddress("IDEPos",      &IDEPos);

    ///////////////////////
    // Create histograms //
    ///////////////////////

    // Histograms with classified events
    TH1D* hPionAbs0p   = new TH1D("hPionAbs0p", "hPionAbs0p;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbsNp   = new TH1D("hPionAbsNp", "hPionAbsNp;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionScatter = new TH1D("hPionScatter", "hPionScatter;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    std::vector<TH1*> RecoSignals = {hPionAbs0p, hPionAbsNp, hPionScatter};

    // Reconstructed incident energy
    TH1D* hPionIncidentKE = new TH1D("hPionIncidentKE", "hPionIncidentKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    
    // True interacting events
    TH1D* hTruePionAbs0p   = new TH1D("hTruePionAbs0p", "hTruePionAbs0p;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTruePionAbsNp   = new TH1D("hTruePionAbsNp", "hTruePionAbsNp;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTruePionScatter = new TH1D("hTruePionScatter", "hTruePionScatter;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    std::vector<TH1*> TotalEventsHistos = {hTruePionAbs0p, hTruePionAbsNp, hTruePionScatter};

    // True incident energy
    TH1D* hTrueIncidentKE = new TH1D("hTrueIncidentKE", "hTrueIncidentKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    //////////////////////
    // Loop over events //
    //////////////////////

    Int_t NumEntries = (Int_t) tree->GetEntries();
    std::cout << "Num entries: " << NumEntries << std::endl;

    bool verbose = false;

    for (Int_t i = 0; i < NumEntries; ++i) {
        // For some reason crashes
        if (i == 1402) continue;

        // Get tree entry and reset variables
        std::cout << std::endl;
        std::cout << "=================================" << std::endl;
        std::cout << "Getting tree entry: " << i << std::endl;
        tree->GetEntry(i);
        std::cout << "Got tree entry" << std::endl;
        std::cout << "Reseting variables" << std::endl;
        EventVariables ev;
        std::cout << "Variables reset" << std::endl;
        std::cout << "=================================" << std::endl;

        //////////////////////////////////////
        // Construct variables I care about //
        //////////////////////////////////////

        // Load track information
        if (verbose) std::cout << std::endl;
        if (verbose) std::cout << "=================================" << std::endl;
        if (verbose) std::cout << "Event: " << event << std::endl;
        if (verbose) std::cout << "Geant4 list size: " << geant_list_size << std::endl;
        if (verbose) std::cout << "Number of primaries: " << no_primaries << std::endl;
        
        // First, we just want to grab WC to TPC match
        int primaryTrackIdx = -1;
        for (size_t trk_idx = 0; trk_idx < ntracks_reco; ++trk_idx) {
            if (trkWCtoTPCMatch[trk_idx]) {
                // Grab position information
                ev.WC2TPCPrimaryBeginX = trkvtxx[trk_idx];
                ev.WC2TPCPrimaryBeginY = trkvtxy[trk_idx];
                ev.WC2TPCPrimaryBeginZ = trkvtxz[trk_idx];
                ev.WC2TPCPrimaryEndX   = trkendx[trk_idx];
                ev.WC2TPCPrimaryEndY   = trkendy[trk_idx];
                ev.WC2TPCPrimaryEndZ   = trkendz[trk_idx];

                // Grab calorimetry information (in collection plane)
                ev.wcMatchResR.assign(trkrr[trk_idx][1], trkrr[trk_idx][1] + ntrkcalopts[trk_idx][1]);
                ev.wcMatchDEDX.assign(trkdedx[trk_idx][1], trkdedx[trk_idx][1] + ntrkcalopts[trk_idx][1]);

                for (size_t dep_idx = 0; dep_idx < ntrkcalopts[trk_idx][1]; ++dep_idx) {
                    ev.wcMatchEDep.push_back(trkdedx[trk_idx][1][dep_idx] * trkpitch[trk_idx][1][dep_idx]);
                    ev.wcMatchXPos.push_back(trkxyz[trk_idx][1][dep_idx][0]);
                    ev.wcMatchYPos.push_back(trkxyz[trk_idx][1][dep_idx][1]);
                    ev.wcMatchZPos.push_back(trkxyz[trk_idx][1][dep_idx][2]);
                }

                // Get location information
                ev.WC2TPCLocationsX.assign(trjPt_X[trk_idx], trjPt_X[trk_idx] + nTrajPoint[trk_idx]);
                ev.WC2TPCLocationsY.assign(trjPt_Y[trk_idx], trjPt_Y[trk_idx] + nTrajPoint[trk_idx]);
                ev.WC2TPCLocationsZ.assign(trjPt_Z[trk_idx], trjPt_Z[trk_idx] + nTrajPoint[trk_idx]);

                // Set flag and index
                primaryTrackIdx = trk_idx;
                ev.WC2TPCMatch     = true;
                ev.WC2TPCsize++;
            }
        }

        if (verbose) std::cout << "Found WC2TPC match: " << ev.WC2TPCMatch  << std::endl;

        // Copy vertex and end for all tracks
        ev.recoEndX.assign(trkendx,  trkendx  + ntracks_reco);
        ev.recoEndY.assign(trkendy,  trkendy  + ntracks_reco);
        ev.recoEndZ.assign(trkendz,  trkendz  + ntracks_reco);
        ev.recoBeginX.assign(trkvtxx, trkvtxx + ntracks_reco);
        ev.recoBeginY.assign(trkvtxy, trkvtxy + ntracks_reco);
        ev.recoBeginZ.assign(trkvtxz, trkvtxz + ntracks_reco);

        // Now, we want to loop through all tracks
        for (size_t trk_idx = 0; trk_idx < ntracks_reco; ++trk_idx) {
            // Grab calorimetry information
            ev.recoResR.push_back(std::vector<double>(trkrr[trk_idx][1], trkrr[trk_idx][1] + ntrkcalopts[trk_idx][1]));
            ev.recoDEDX.push_back(std::vector<double>(trkdedx[trk_idx][1], trkdedx[trk_idx][1] + ntrkcalopts[trk_idx][1]));

            // Check reversed
            double startDistance = distance(trkvtxx[trk_idx], ev.WC2TPCPrimaryEndX, trkvtxy[trk_idx], ev.WC2TPCPrimaryEndY, trkvtxz[trk_idx], ev.WC2TPCPrimaryEndZ);
            double endDistance   = distance(trkendx[trk_idx], ev.WC2TPCPrimaryEndX, trkendy[trk_idx], ev.WC2TPCPrimaryEndY, trkendz[trk_idx], ev.WC2TPCPrimaryEndZ);

            if (startDistance > endDistance && !trkWCtoTPCMatch[trk_idx]) {
                ev.isTrackInverted.push_back(true);

                std::swap(ev.recoEndX[trk_idx], ev.recoBeginX[trk_idx]);
                std::swap(ev.recoEndY[trk_idx], ev.recoBeginY[trk_idx]);
                std::swap(ev.recoEndZ[trk_idx], ev.recoBeginZ[trk_idx]);

                std::reverse(ev.recoResR[trk_idx].begin(), ev.recoResR[trk_idx].end());
                std::reverse(ev.recoDEDX[trk_idx].begin(), ev.recoDEDX[trk_idx].end());
            } else {
                ev.isTrackInverted.push_back(false);
            }
        }

        // Load wire-chamber track information
        ev.WCTrackMomentum = wctrk_momentum[0];
        ev.WCTheta         = wctrk_theta[0];
        ev.WCPhi           = wctrk_phi[0];
        ev.WC4PrimaryX     = WC4xPos[0];

        // In case we overwrite for elastic scattering
        double primaryGeantEndX, primaryGeantEndY, primaryGeantEndZ; 

        // Get information about Geant4 primary particles
        int iPrimary = 0;
        TVector3 vtxPion, outPion;
        for (size_t true_idx = 0; true_idx < geant_list_size; ++true_idx) {
            // Get truth primary
            if (process_primary[true_idx] && StartPointz[true_idx] == -100.) {
                ev.truthPrimaryPDG        = pdg[true_idx];
                ev.truthPrimaryVertexX    = EndPointx[true_idx];
                ev.truthPrimaryVertexY    = EndPointy[true_idx];
                ev.truthPrimaryVertexZ    = EndPointz[true_idx];
                ev.truthPrimaryVertexKE   = EndEng[true_idx] - Mass[true_idx];
                vtxPion                   = TVector3(EndPx[true_idx], EndPy[true_idx], EndPz[true_idx]);

                primaryGeantEndX = EndPointx[true_idx]; primaryGeantEndY = EndPointy[true_idx]; primaryGeantEndZ = EndPointz[true_idx];

                if (verbose) std::cout << "Interacting pion mom.: " << EndPx[true_idx] << "   " << EndPy[true_idx] << "   " << EndPz[true_idx] << std::endl;
                if (verbose) std::cout << "Primary id: " << TrackId[true_idx] << std::endl;
                if (verbose) std::cout << "Primary pdg: " << pdg[true_idx] << std::endl;
                if (verbose) std::cout << "Number traj points " << NTrTrajPts[true_idx] << std::endl;
                if (verbose) std::cout << "Number daughters " << NumberDaughters[true_idx] << std::endl;
                int found = 0;

                // Get daughters
                for (size_t inner_idx = 0; inner_idx < geant_list_size; ++inner_idx) {
                    // Find daughters of primary
                    if (Mother[inner_idx] == TrackId[true_idx]) {
                        found++;
                        ev.truthPrimaryDaughtersPDG.push_back(pdg[inner_idx]);
                        ev.truthPrimaryDaughtersProcess.push_back(ProcessToString(Process[inner_idx]));
                        ev.truthPrimaryDaughtersKE.push_back(Eng[inner_idx] - Mass[inner_idx]);

                        // Get daughter -211
                        if (pdg[inner_idx] == -211) {
                            ev.truthScatteredPionKE = Eng[inner_idx] - Mass[inner_idx];
                            outPion              = TVector3(Px[inner_idx], Py[inner_idx], Pz[inner_idx]);
                        }
                    }
                }
                if (verbose) std::cout << "Found primary daughters: " << found << std::endl;
                
                // Get initial momentum for energy loss
                ev.trajectoryInitialMomentumX = 1000 * MidPx[iPrimary][0];

                // Get trajectory information
                ev.truthPrimaryLocationX.assign(MidPosX[iPrimary], MidPosX[iPrimary] + NTrTrajPts[true_idx]);
                ev.truthPrimaryLocationY.assign(MidPosY[iPrimary], MidPosY[iPrimary] + NTrTrajPts[true_idx]);
                ev.truthPrimaryLocationZ.assign(MidPosZ[iPrimary], MidPosZ[iPrimary] + NTrTrajPts[true_idx]);

                // Get interaction in trajectory
                float lastTPCX = -999, lastTPCY = -999, lastTPCZ = -999;
                ev.interactionInTrajectory = false;
                if (verbose) std::cout << "Looking at trajectory" << std::endl;
                if (InteractionPoint->at(0) != NTrTrajPts[true_idx] - 1) {
                    for (int i_point = 0; i_point < InteractionPoint->size(); ++i_point) {
                        if (i_point >= InteractionPointType->size()) continue;

                        if (verbose) std::cout << "   Point: " << InteractionPoint->at(i_point) << "   Interaction: " << InteractionPointType->at(i_point) << std::endl;

                        // "hadElastic" is type 3
                        if (InteractionPointType->at(i_point) == 3) {
                            // Get kinetic energy at this point
                            float p = std::sqrt(
                                MidPx[iPrimary][InteractionPoint->at(i_point)]*MidPx[iPrimary][InteractionPoint->at(i_point)] + 
                                MidPy[iPrimary][InteractionPoint->at(i_point)]*MidPy[iPrimary][InteractionPoint->at(i_point)] + 
                                MidPz[iPrimary][InteractionPoint->at(i_point)]*MidPz[iPrimary][InteractionPoint->at(i_point)]
                            );
                            float KE = std::sqrt(p*p + Mass[true_idx]*Mass[true_idx]) - Mass[true_idx];

                            if (!ev.interactionInTrajectory) {
                                // First interaction in trajectory
                                ev.interactionInTrajectory    = true;
                                ev.trajectoryInteractionLabel = "hadElastic";
                                ev.trajectoryInteractionKE    = KE;

                                lastTPCX = MidPosX[iPrimary][InteractionPoint->at(i_point)];
                                lastTPCY = MidPosY[iPrimary][InteractionPoint->at(i_point)];
                                lastTPCZ = MidPosZ[iPrimary][InteractionPoint->at(i_point)];

                                ev.truthPrimaryVertexX = lastTPCX;
                                ev.truthPrimaryVertexY = lastTPCY;
                                ev.truthPrimaryVertexZ = lastTPCZ;

                                TVector3 momBefore, momAfter;
                                if (InteractionPoint->at(i_point) - 1 >= 0) {
                                    momBefore = TVector3(
                                        MidPx[iPrimary][InteractionPoint->at(i_point)-1], 
                                        MidPy[iPrimary][InteractionPoint->at(i_point)-1], 
                                        MidPz[iPrimary][InteractionPoint->at(i_point)-1]
                                    );
                                } else {
                                    momBefore = TVector3(
                                        MidPx[iPrimary][InteractionPoint->at(i_point)], 
                                        MidPy[iPrimary][InteractionPoint->at(i_point)],
                                        MidPz[iPrimary][InteractionPoint->at(i_point)]
                                    );
                                }
                                momAfter                      = TVector3(MidPx[iPrimary][InteractionPoint->at(i_point)], MidPy[iPrimary][InteractionPoint->at(i_point)], MidPz[iPrimary][InteractionPoint->at(i_point)]);
                                ev.trajectoryInteractionAngle = momBefore.Angle(momAfter);
                            } else if (ev.interactionInTrajectory) {
                                // Subsequent interactions in trajectory
                                ev.secondaryInteractionTypes.push_back(12);
                                ev.secondaryInteractionXPosition.push_back(MidPosX[iPrimary][InteractionPoint->at(i_point)]);
                                ev.secondaryInteractionYPosition.push_back(MidPosY[iPrimary][InteractionPoint->at(i_point)]);
                                ev.secondaryInteractionZPosition.push_back(MidPosZ[iPrimary][InteractionPoint->at(i_point)]);
                                ev.secondaryInteractionInteractingKE.push_back(KE);

                                ev.secondaryInteractionDaughtersPDG.push_back({});
                                ev.secondaryInteractionDaughtersKE.push_back({});

                                TVector3 momBefore, momAfter;
                                if (InteractionPoint->at(i_point) - 1 >= 0) {
                                    momBefore = TVector3(
                                        MidPx[iPrimary][InteractionPoint->at(i_point)-1], 
                                        MidPy[iPrimary][InteractionPoint->at(i_point)-1], 
                                        MidPz[iPrimary][InteractionPoint->at(i_point)-1]
                                    );
                                } else {
                                    momBefore = TVector3(
                                        MidPx[iPrimary][InteractionPoint->at(i_point)], 
                                        MidPy[iPrimary][InteractionPoint->at(i_point)],
                                        MidPz[iPrimary][InteractionPoint->at(i_point)]
                                    );
                                }
                                momAfter                   = TVector3(MidPx[iPrimary][InteractionPoint->at(i_point)], MidPy[iPrimary][InteractionPoint->at(i_point)], MidPz[iPrimary][InteractionPoint->at(i_point)]);
                                ev.secondaryInteractionAngle.push_back(momBefore.Angle(momAfter));
                            }
                        }
                    }
                    if (verbose) std::cout << "Finished trajectory" << std::endl;
                }

                // If no interaction in trajectory, last point found by looping backwards
                if (!ev.interactionInTrajectory) {
                    for (int i = ev.truthPrimaryLocationX.size() - 1; i >= 0; i--) {
                        if (isWithinActiveVolume(ev.truthPrimaryLocationX.at(i), ev.truthPrimaryLocationY.at(i), ev.truthPrimaryLocationZ.at(i))) {
                            lastTPCX = ev.truthPrimaryLocationX.at(i);
                            lastTPCY = ev.truthPrimaryLocationY.at(i);
                            lastTPCZ = ev.truthPrimaryLocationZ.at(i);
                            break;
                        }
                    }
                }

                // Find first point in TPC
                float firstTPCX  = -999, firstTPCY  = -999, firstTPCZ  = -999;
                float firstTPCPx = -999, firstTPCPy = -999, firstTPCPz = -999;
                for (size_t i = 0; i < ev.truthPrimaryLocationX.size() - 1; i++) {
                    if (isWithinActiveVolume(ev.truthPrimaryLocationX.at(i), ev.truthPrimaryLocationY.at(i), ev.truthPrimaryLocationZ.at(i))) {
                        firstTPCX = ev.truthPrimaryLocationX.at(i);
                        firstTPCY = ev.truthPrimaryLocationY.at(i);
                        firstTPCZ = ev.truthPrimaryLocationZ.at(i);

                        firstTPCPx = MidPx[iPrimary][i];
                        firstTPCPy = MidPy[iPrimary][i];
                        firstTPCPz = MidPz[iPrimary][i];
                        break;
                    }
                }

                // Get true incident KE
                ev.validTrueIncidentKE = true;
                if (ev.truthPrimaryPDG != -211) ev.validTrueIncidentKE = false;
                if (firstTPCX == -999)          ev.validTrueIncidentKE = false;
                if (firstTPCX == lastTPCX)      ev.validTrueIncidentKE = false;

                double totalLength = distance(firstTPCX, lastTPCX, firstTPCY, lastTPCY, firstTPCZ, lastTPCZ);
                
                if (verbose) std::cout << "First TPC point: " << firstTPCX << "  " << firstTPCY << "  " << firstTPCZ << std::endl;
                if (verbose) std::cout << "Last TPC point:  " << lastTPCX << "  " << lastTPCY << "  " << lastTPCZ << std::endl;
                if (verbose) std::cout << "Total length:    " << totalLength << std::endl;
                if (verbose) std::cout << "Is true incident valid: " << ev.validTrueIncidentKE << std::endl;

                // If length more than track pitch, get contributions
                if (ev.validTrueIncidentKE && totalLength >= TRACK_PITCH) {
                    std::map<double, TVector3> orderedUniformTrjPts;
                    auto positionVector0 = TVector3(firstTPCX, firstTPCY, firstTPCZ);
                    auto positionVector1 = TVector3(lastTPCX, lastTPCY, lastTPCZ);
                    orderedUniformTrjPts[positionVector0.Z()] = positionVector0;
                    orderedUniformTrjPts[positionVector1.Z()] = positionVector1;

                    int numberPts = (int) (totalLength / TRACK_PITCH);
                    for (int iPoint = 1; iPoint <= numberPts; ++iPoint) {
                        auto newPoint = positionVector0 + iPoint * (TRACK_PITCH / totalLength) * (positionVector1 - positionVector0);
                        orderedUniformTrjPts[newPoint.Z()] = newPoint;
                    }

                    // If distance between last point and second to last is less than 0.235, eliminate second to last
                    auto lastPt         = (orderedUniformTrjPts.rbegin())->second;
                    auto secondtoLastPt = (std::next(orderedUniformTrjPts.rbegin()))->second;
                    double lastDist     = distance(lastPt.X(), secondtoLastPt.X(), lastPt.Y(), secondtoLastPt.Y(), lastPt.Z(), secondtoLastPt.Z());
                    if (lastDist < (TRACK_PITCH / 2)) orderedUniformTrjPts.erase((std::next(orderedUniformTrjPts.rbegin()))->first);

                    double trueKineticEnergy = 1000 * (
                        TMath::Sqrt(
                            firstTPCPx * firstTPCPx + 
                            firstTPCPy * firstTPCPy +
                            firstTPCPz * firstTPCPz +
                            Mass[true_idx] * Mass[true_idx]
                        ) - Mass[true_idx]
                    );
                    for (auto it = std::next(orderedUniformTrjPts.begin()), old_it = orderedUniformTrjPts.begin(); it != orderedUniformTrjPts.end(); it++, old_it++) {
                        auto oldPos        = old_it->second;
                        auto currentPos    = it->second;
                        double uniformDist = (currentPos - oldPos).Mag();

                        double currentDepEnergy = 0.;
                        for (int i = 0; i < maxTrackIDE; i++) {
                            if (IDEPos[i][2] < oldPos.Z())     continue;
                            if (IDEPos[i][2] > currentPos.Z()) continue;
                            currentDepEnergy += IDEEnergy[i];
                        }

                        // Skip tiny depositions
                        if (currentDepEnergy / uniformDist < 0.1) continue;

                        // Calculate current KE
                        trueKineticEnergy -= currentDepEnergy;

                        if (isWithinReducedVolume(currentPos.X(), currentPos.Y(), currentPos.Z())) {
                            ev.trueIncidentKEContributions.push_back(trueKineticEnergy);
                        }

                        // std::cout << "Step " << std::distance(orderedUniformTrjPts.begin(), it)
                        // << ": oldPos=("     << oldPos.X()     << ", " << oldPos.Y()     << ", " << oldPos.Z()     << ")"
                        // << "  currentPos=(" << currentPos.X() << ", " << currentPos.Y() << ", " << currentPos.Z() << ")"
                        // << "  dist="        << uniformDist
                        // << "  KE="          << trueKineticEnergy
                        // << "  depEnergy="   << currentDepEnergy
                        // << std::endl;
                    }
                } else if (ev.validTrueIncidentKE) {
                    // If length less than track pitch, empty contributions
                    ev.trueIncidentKEContributions.push_back({});
                }

                // Only one primary starting at -100 (particle gun)
                break;
            } else if (process_primary[true_idx]) iPrimary++;
        }

        if (!(ev.truthPrimaryPDG == -211 || ev.truthPrimaryPDG == 11 || ev.truthPrimaryPDG == 13)) {
            std::cout << "Weird primary PDG: " << ev.truthPrimaryPDG << ", skipping event" << std::endl;
            continue;
        }

        // Get scattering angle
        ev.truthScatteringAngle = vtxPion.Angle(outPion);

        // Get information about wire hits
        std::map<int, std::vector<int>>    trackHitMap;
        std::map<int, std::vector<double>> trackHitXMap;
        std::map<int, std::vector<double>> trackHitYMap;
        std::map<int, std::vector<double>> trackHitZMap;

        for (size_t i_hit = 0; i_hit < nhits; ++i_hit) {
            ev.fHitPlane.push_back(hit_plane[i_hit]);
            ev.fHitT.push_back(hit_driftT[i_hit]);
            ev.fHitX.push_back(hit_driftT[i_hit] * DRIFT_VELOCITY);
            ev.fHitW.push_back(hit_channel[i_hit]);

            // If it is -9, no match to a track
            if (hit_trkid[i_hit] != -9) {
                ev.hitRecoAsTrackKey.push_back(i_hit);
                if (hit_trkid[i_hit] == primaryTrackIdx) ev.hitWC2TPCKey.push_back(i_hit);

                trackHitMap[hit_trkid[i_hit]].push_back(i_hit);
                trackHitXMap[hit_trkid[i_hit]].push_back(hit_x[i_hit]);
                trackHitYMap[hit_trkid[i_hit]].push_back(hit_y[i_hit]);
                trackHitZMap[hit_trkid[i_hit]].push_back(hit_z[i_hit]);
            }
        }

        // Fill out vectors
        for (auto& [trkid, hits] : trackHitMap) {
            if (trkid >= (int) ev.recoTrackHitIndices.size()) {
                ev.recoTrackHitIndices.resize(trkid + 1);
                ev.recoTrackHitX.resize(trkid + 1);
                ev.recoTrackHitY.resize(trkid + 1);
                ev.recoTrackHitZ.resize(trkid + 1);
            }
            ev.recoTrackHitIndices[trkid] = hits;
            ev.recoTrackHitX[trkid]       = trackHitXMap[trkid];
            ev.recoTrackHitY[trkid]       = trackHitYMap[trkid];
            ev.recoTrackHitZ[trkid]       = trackHitZMap[trkid];
        }

        // Sanity check
        // for (size_t i_trk = 0; i_trk < recoTrackHitIndices->size(); ++i_trk) {
        //     std::cout << "Track " << i_trk << ":" << std::endl;
        //     for (int idx : (*recoTrackHitIndices)[i_trk]) {
        //         std::cout << "  hit index=" << idx << "  hit_trkid=" << hit_trkid[idx] << std::endl;
        //     }
        // }

        if (verbose) std::cout << "Hits associated to WC2TPC match: " << ev.hitWC2TPCKey.size() << std::endl;

        // Get truth background type
        ev.backgroundType = fillSignalInformation(
            ev.truthPrimaryPDG,
            ev.truthPrimaryVertexX,
            ev.truthPrimaryVertexY,
            ev.truthPrimaryVertexZ,
            ev.interactionInTrajectory,
            ev.trajectoryInteractionLabel,
            ev.truthPrimaryDaughtersPDG,
            ev.truthPrimaryDaughtersProcess,
            ev.truthPrimaryDaughtersKE,
            true,
            ev.numVisibleProtons
        );

        if (verbose) std::cout << "Background type: " << ev.backgroundType << std::endl;

        // For elastic scattering, record what happens at end of track
        if (ev.backgroundType == 12) {
            ev.secondaryInteractionTypes.push_back(
                fillSignalInformation(
                    ev.truthPrimaryPDG,
                    primaryGeantEndX,
                    primaryGeantEndY,
                    primaryGeantEndZ,
                    false,
                    "",
                    ev.truthPrimaryDaughtersPDG,
                    ev.truthPrimaryDaughtersProcess,
                    ev.truthPrimaryDaughtersKE,
                    false,
                    ev.numVisibleProtons
                )
            );
            ev.secondaryInteractionInteractingKE.push_back(ev.truthPrimaryVertexKE);
            ev.secondaryInteractionAngle.push_back(ev.truthScatteringAngle);
            ev.secondaryInteractionXPosition.push_back(primaryGeantEndX);
            ev.secondaryInteractionYPosition.push_back(primaryGeantEndY);
            ev.secondaryInteractionZPosition.push_back(primaryGeantEndZ);
            ev.secondaryInteractionDaughtersPDG.push_back(ev.truthPrimaryDaughtersPDG);
            ev.secondaryInteractionDaughtersKE.push_back(ev.truthPrimaryDaughtersKE);
        }

        // Grab incident KE for secondary interactions along primary track
        for (int i_seg = 0; i_seg < ev.secondaryInteractionXPosition.size(); i_seg++) {
            TVector3 segStart;
            if (i_seg == 0) {
                segStart = TVector3(
                    ev.truthPrimaryVertexX,
                    ev.truthPrimaryVertexY,
                    ev.truthPrimaryVertexZ
                );
            } else {
                segStart = TVector3(
                    ev.secondaryInteractionXPosition.at(i_seg - 1),
                    ev.secondaryInteractionYPosition.at(i_seg - 1),
                    ev.secondaryInteractionZPosition.at(i_seg - 1)
                );
            }
            
            auto segEnd   = TVector3(
                ev.secondaryInteractionXPosition.at(i_seg),
                ev.secondaryInteractionYPosition.at(i_seg),
                ev.secondaryInteractionZPosition.at(i_seg)
            );

            if ((segEnd - segStart).Mag() < TRACK_PITCH) {
                ev.secondaryIncidentKEContributions.push_back({});
                continue;
            }

            std::map<double, TVector3> orderedUniformTrjPtsSecondary;
            orderedUniformTrjPtsSecondary[segStart.Z()] = segStart;
            orderedUniformTrjPtsSecondary[segEnd.Z()]   = segEnd;

            double segLength = (segEnd - segStart).Mag();
            int numberPts    = (int)(segLength / TRACK_PITCH);
            for (int iPoint = 1; iPoint <= numberPts; ++iPoint) {
                auto newPoint = segStart + iPoint * (TRACK_PITCH / segLength) * (segEnd - segStart);
                orderedUniformTrjPtsSecondary[newPoint.Z()] = newPoint;
            }

            // If distance between last point and second to last is less than half pitch, eliminate second to last
            auto lastPt         = (orderedUniformTrjPtsSecondary.rbegin())->second;
            auto secondtoLastPt = (std::next(orderedUniformTrjPtsSecondary.rbegin()))->second;
            double lastDist     = distance(lastPt.X(), secondtoLastPt.X(), lastPt.Y(), secondtoLastPt.Y(), lastPt.Z(), secondtoLastPt.Z());
            if (lastDist < (TRACK_PITCH / 2)) orderedUniformTrjPtsSecondary.erase((std::next(orderedUniformTrjPtsSecondary.rbegin()))->first);

            std::vector<double> segIncidentKE = {};
            double segKineticEnergy = ev.secondaryInteractionInteractingKE.at(i_seg) * 1000;
            for (auto it = std::next(orderedUniformTrjPtsSecondary.begin()), old_it = orderedUniformTrjPtsSecondary.begin(); it != orderedUniformTrjPtsSecondary.end(); it++, old_it++) {
                auto oldPos     = old_it->second;
                auto currentPos = it->second;
                double uniformDist = (currentPos - oldPos).Mag();

                double currentDepEnergy = 0.;
                for (int i = 0; i < maxTrackIDE; i++) {
                    if (IDEPos[i][2] < oldPos.Z())     continue;
                    if (IDEPos[i][2] > currentPos.Z()) continue;
                    currentDepEnergy += IDEEnergy[i];
                }

                if (currentDepEnergy / uniformDist < 0.1) continue;

                segKineticEnergy -= currentDepEnergy;

                if (isWithinReducedVolume(currentPos.X(), currentPos.Y(), currentPos.Z())) {
                    segIncidentKE.push_back(segKineticEnergy);
                }
            }
            ev.secondaryIncidentKEContributions.push_back(segIncidentKE);
        }

        // Make script go faster
        // if (i > USE_NUM_EVENTS) break;

        std::unordered_set<int> hitsInTracks(ev.hitRecoAsTrackKey.begin(), ev.hitRecoAsTrackKey.end());

        // Extend reco cylinder
        removeRepeatedPoints(&ev.WC2TPCLocationsX, &ev.WC2TPCLocationsY, &ev.WC2TPCLocationsZ);
        std::vector<double>* wcX = new std::vector<double>(ev.WC2TPCLocationsX);
        std::vector<double>* wcY = new std::vector<double>(ev.WC2TPCLocationsY);
        std::vector<double>* wcZ = new std::vector<double>(ev.WC2TPCLocationsZ);

        // Get direction to end cylinder
        int numPoints = wcX->size();
        int numTail   = std::min(10, numPoints - 1);
        std::vector<std::vector<double>> points;
        for (int j = numPoints - numTail - 1; j < numPoints; ++j) {
            points.push_back({
                wcX->at(j),
                wcY->at(j),
                wcZ->at(j)
            });
        }

        std::vector<double> avgDir(3, 0);
        if (numTail > 0) {
            avgDir = getAverageDir(points);

            // Extrapolate track to end
            double scale = (maxZ - points.back()[2]) / avgDir[2];
            wcX->push_back(points.back()[0] + scale * avgDir[0]);
            wcY->push_back(points.back()[1] + scale * avgDir[1]);
            wcZ->push_back(points.back()[2] + scale * avgDir[2]);
        }

        // Extend truth cylinder
        std::vector<double>* truthCylinderLocationX = new std::vector<double>(ev.truthPrimaryLocationX);
        std::vector<double>* truthCylinderLocationY = new std::vector<double>(ev.truthPrimaryLocationY);
        std::vector<double>* truthCylinderLocationZ = new std::vector<double>(ev.truthPrimaryLocationZ);

        if (truthCylinderLocationX->size() == 0) {
            std::cout << "Primary trajectory has no size" << std::endl;
        }

        // Get direction to end cylinder
        int truthNumPoints = truthCylinderLocationX->size();
        int truthNumTail   = std::min(10, truthNumPoints - 1);
        std::vector<std::vector<double>> truth_points;
        for (int j = truthNumPoints - truthNumTail; j < truthNumPoints; ++j) {
            truth_points.push_back({
                truthCylinderLocationX->at(j),
                truthCylinderLocationY->at(j),
                truthCylinderLocationZ->at(j)
            });
        }

        //////////////////////////////////
        // Look at unreconstructed hits //
        //////////////////////////////////

        // Units for coordinates:
        //     W: (channel) * (channel width) = cm
        //     X: (drift velocity) * (time) = (cm / us) * us = cm

        // First, find hits near the primary track that were not
        // reconstructed into any of the already existing tracks
        std::vector<int> candidateHits;

        int numTotalHitsNearPrimary            = 0;
        int numUnRecoHitsNearPrimaryInduction  = 0;
        int numUnRecoHitsNearPrimaryCollection = 0;
        for (size_t iHit = 0; iHit < nhits; ++iHit) {
            double hitX     = ev.fHitX.at(iHit);
            double hitW     = ev.fHitW.at(iHit);
            int    hitPlane = ev.fHitPlane.at(iHit);

            // Check if hit is near vertex of the primary
            if (isHitNearPrimary(
                &ev.hitWC2TPCKey,
                &ev.fHitX,
                &ev.fHitW,
                &ev.fHitPlane,
                hitX,
                hitW,
                hitPlane,
                DISTANCE_TO_PRIMARY_THRESHOLD,
                true
            )) {
                numTotalHitsNearPrimary++;

                // Skip hits already in tracks
                if (hitsInTracks.count(iHit) > 0) continue;

                if (hitPlane == 0) numUnRecoHitsNearPrimaryInduction++;
                else if (hitPlane == 1) numUnRecoHitsNearPrimaryCollection++;

                candidateHits.push_back(iHit); 
            }
        }

        // Now cluster using those hits as starting points
        std::unordered_set<int> usedHits;
        std::vector<HitCluster> hitClusters;
        int nCandidateHits = candidateHits.size();

        // Hits in the same cluster must be separated by at most some number of wires
        for (int iHit = 0; iHit < nCandidateHits; ++iHit) {
            // The candidate hits are only STARTING points, as this could make up 
            // really long tracks in the induction plane that are no longer near 
            // the ending point of the primary track
            
            // First, check if we have already used this hit
            int   thisHitKey       = candidateHits.at(iHit);
            int   thisHitPlane     = ev.fHitPlane.at(thisHitKey);
            float thisHitW         = ev.fHitW.at(thisHitKey);
            float thisHitX         = ev.fHitX.at(thisHitKey);
            float thisHitCharge    = -1;
            float thisHitChargeCol = -1;
            
            if (usedHits.count(thisHitKey)) continue;           
            
            std::vector<int> clusterKeys;
            std::vector<float> clusterX;
            std::vector<float> clusterW;
            std::vector<float> clusterCharge;
            std::vector<float> clusterChargeCol;

            clusterKeys.push_back(thisHitKey);
            clusterX.push_back(thisHitX);
            clusterW.push_back(thisHitW);
            clusterCharge.push_back(thisHitCharge);
            clusterChargeCol.push_back(thisHitChargeCol);
            
            for (int iAllHit = 0; iAllHit < nhits; ++iAllHit) {
                // Skip already used hits, and those reconstructed in tracks
                if (usedHits.count(iAllHit) || hitsInTracks.count(iAllHit)) continue;

                // Clusters have to be in same plane
                if (ev.fHitPlane.at(iAllHit) != thisHitPlane) continue;

                float internalHitW  = ev.fHitW.at(iAllHit);
                float internalHitX  = ev.fHitX.at(iAllHit);
                float dW            = std::abs(internalHitW - thisHitW);
                float dX            = std::abs(internalHitX - thisHitX);
                float distance      = std::sqrt(std::pow(dW, 2) + std::pow(dX, 2));

                int nClusterSoFar = clusterW.size();
                for (int iCluster = 0; iCluster < nClusterSoFar; ++iCluster) {
                    float tempdW       = std::abs(internalHitW - clusterW.at(iCluster));
                    float tempdX       = std::abs(internalHitX - clusterX.at(iCluster));
                    float tempDistance = std::sqrt(std::pow(tempdW, 2) + std::pow(tempdX, 2));

                    if (tempDistance < distance) distance = tempDistance;
                    if (distance < MAX_IN_CLUSTER_SEPARATION) break;
                }

                if (distance < MAX_IN_CLUSTER_SEPARATION) {
                    usedHits.insert(iAllHit);
                    clusterKeys.push_back(iAllHit);
                    clusterX.push_back(ev.fHitX.at(iAllHit));
                    clusterW.push_back(ev.fHitW.at(iAllHit));
                    clusterCharge.push_back(-1);
                    clusterChargeCol.push_back(-1);
                }
            }

            if (clusterKeys.size() > MINIMUM_HITS_FOR_CLUSTER) {
                HitCluster thisCluster;

                double clusterSize = 0;
                for (int j = 0; j < clusterW.size() - 1; ++j) {
                    for (int k = j + 1; k < clusterW.size(); ++k) {
                        double wSeparation = clusterW[j] - clusterW[k];
                        double xSeparation = clusterX[j] - clusterX[k];

                        double thisDiameter = std::sqrt(
                            std::pow(wSeparation, 2) + 
                            std::pow(xSeparation, 2)
                        );

                        if (thisDiameter > clusterSize) clusterSize = thisDiameter;
                    }
                }

                thisCluster.plane        = thisHitPlane;
                thisCluster.hitKeys      = clusterKeys;
                thisCluster.hitX         = clusterX;
                thisCluster.hitW         = clusterW;
                thisCluster.hitCharge    = clusterCharge;
                thisCluster.hitChargeCol = clusterChargeCol;
                thisCluster.clusterSize  = clusterSize;
                
                usedHits.insert(thisHitKey);
                hitClusters.push_back(thisCluster);
            }
        }


        ////////////////////////////
        // Save truth information //
        ////////////////////////////

        // Scattering only if degree > THRESHOLD and energy > THRESHOLD
        double scatteringAngle  = -9999;
        double scatteringEnergy = -9999;

        // Modify scatterings
        if (ev.backgroundType == 12 || ev.backgroundType == 6) {
            if (ev.backgroundType == 12) {
                scatteringAngle         = ev.trajectoryInteractionAngle;
                scatteringEnergy        = ev.trajectoryInteractionKE;
                ev.truthPrimaryVertexKE = ev.trajectoryInteractionKE; // in case we do not modify anything
            }
            else if (ev.backgroundType == 6) {
                scatteringAngle  = ev.truthScatteringAngle;
                scatteringEnergy = ev.truthScatteredPionKE; 
            }

            if (verbose) std::cout << "Checking if elastic scattering is above threshold, with angle " << scatteringAngle << " and KE " << scatteringEnergy << std::endl;

            // If outgoing pion below threshold, absorption
            if (scatteringEnergy < PION_SCATTERING_ENERGY_THRESHOLD) {
                if (ev.backgroundType == 12) ev.backgroundType = 0;
                else if (ev.backgroundType == 6) {
                    int tempNumVisibleProtons = 0;
                    for (int i = 0; i < ev.truthPrimaryDaughtersPDG.size(); ++i) {
                        if (ev.truthPrimaryDaughtersPDG.at(i) == 2212) {
                            if (
                                ev.truthPrimaryDaughtersKE.at(i) > PROTON_ENERGY_LOWER_BOUND &&
                                ev.truthPrimaryDaughtersKE.at(i) < PROTON_ENERGY_UPPER_BOUND
                            ) tempNumVisibleProtons++;
                        }
                    }
                    if (tempNumVisibleProtons == 0) ev.backgroundType = 0;
                    else ev.backgroundType = 1;
                    ev.numVisibleProtons = tempNumVisibleProtons;
                }
            }
            // If pion above threshold but angle not large enough, go to next interaction and check there
            else if (scatteringAngle < SCATTERING_ANGLE_THRESHOLD) {
                // Use secondary interaction for 
                for (int iInteraction = 0; iInteraction < ev.secondaryInteractionTypes.size(); ++iInteraction) {
                    int currentInteraction = ev.secondaryInteractionTypes.at(iInteraction);
                    scatteringAngle        = ev.secondaryInteractionAngle.at(iInteraction);
                    scatteringEnergy       = ev.secondaryInteractionInteractingKE.at(iInteraction);

                    // Get scattering energy for inelastic scattering from outgoing pion kinematics
                    if (currentInteraction == 6) {
                        for (int i = 0; i < ev.secondaryInteractionDaughtersPDG.at(iInteraction).size(); ++i) {
                            if (ev.secondaryInteractionDaughtersPDG.at(iInteraction)[i] == -211) {
                                scatteringEnergy = ev.secondaryInteractionDaughtersKE.at(iInteraction)[i];
                                break;
                            }
                        }
                    }

                    // Add incident slices true contributions
                    for (int iContribution = 0; iContribution < ev.secondaryIncidentKEContributions.at(iInteraction).size(); ++iContribution) {
                        ev.trueIncidentKEContributions.push_back(ev.secondaryIncidentKEContributions.at(iInteraction)[iContribution]);
                    }

                    // If scattering but outgoing below threshold, absorption
                    if (
                        (currentInteraction == 6 || currentInteraction == 12) &&
                        scatteringEnergy < PION_SCATTERING_ENERGY_THRESHOLD
                    ) {
                        if (currentInteraction == 12) ev.backgroundType = 0;
                        else if (currentInteraction == 6) {
                            int tempNumVisibleProtons = 0;
                            for (int i = 0; i < ev.secondaryInteractionDaughtersPDG.at(iInteraction).size(); ++i) {
                                if (ev.secondaryInteractionDaughtersPDG.at(iInteraction)[i] == 2212) {
                                    if (
                                        ev.secondaryInteractionDaughtersKE.at(iInteraction)[i] > PROTON_ENERGY_LOWER_BOUND &&
                                        ev.secondaryInteractionDaughtersKE.at(iInteraction)[i] < PROTON_ENERGY_UPPER_BOUND
                                    ) tempNumVisibleProtons++;
                                }
                            }
                            if (tempNumVisibleProtons == 0) ev.backgroundType = 0;
                            else ev.backgroundType = 1;
                            ev.numVisibleProtons = tempNumVisibleProtons;
                        }
                        ev.truthPrimaryVertexKE = ev.secondaryInteractionInteractingKE.at(iInteraction);
                        ev.truthPrimaryVertexX  = ev.secondaryInteractionXPosition.at(iInteraction);
                        ev.truthPrimaryVertexY  = ev.secondaryInteractionYPosition.at(iInteraction);
                        ev.truthPrimaryVertexZ  = ev.secondaryInteractionZPosition.at(iInteraction);
                        break;
                    }
                    // If scattering energy above threshold but angle not large enough, keep going
                    else if (
                        (currentInteraction == 6 || currentInteraction == 12) &&
                        scatteringAngle < SCATTERING_ANGLE_THRESHOLD
                    ) {
                        // We want to keep going
                        continue;
                    } else {
                        // We found non-scattering interaction or scattering interaction with angle and energy above our threshold value
                        ev.backgroundType       = currentInteraction;
                        ev.truthPrimaryVertexKE = ev.secondaryInteractionInteractingKE.at(iInteraction);
                        ev.truthPrimaryVertexX  = ev.secondaryInteractionXPosition.at(iInteraction);
                        ev.truthPrimaryVertexY  = ev.secondaryInteractionYPosition.at(iInteraction);
                        ev.truthPrimaryVertexZ  = ev.secondaryInteractionZPosition.at(iInteraction);
                        break;
                    }
                }
            }
        }

        ///////////////////////////////
        // Get weight for this event //
        ///////////////////////////////

        double fd_weight= FakeDataFD::GetFDWeight(FDScenario, ev.backgroundType, ev.truthPrimaryVertexKE * 1000);
        // if (fd_weight != 1) {
        //     std::cout << "Background type: " << ev.backgroundType << std::endl;
        //     std::cout << "Truth primary vertex KE: " << ev.truthPrimaryVertexKE * 1000 << " MeV" << std::endl;
        //     std::cout << "Fake data weight applied: " << fd_weight << std::endl;
        //     std::cout << std::endl;
        // }

        ///////////////////
        // True spectrum //
        ///////////////////

        if (ev.trueIncidentKEContributions.size() == 0) ev.truthPrimaryVertexKE = 0;
        else ev.truthPrimaryVertexKE = ev.trueIncidentKEContributions.back() / 1000.0;

        // Get true energy bin
        int TrueEnergyBin = getBin(ev.truthPrimaryVertexKE * 1000, ARRAY_KE_BINS);

        // Add true incident KE
        if (ev.validTrueIncidentKE) {
            for (double x : ev.trueIncidentKEContributions) {
                hTrueIncidentKE->Fill(x);
            }
        }

        // True interaction histograms
        if (ev.backgroundType == 0) {
            hTruePionAbs0p->Fill(ev.truthPrimaryVertexKE * 1000, fd_weight);
        } else if (ev.backgroundType == 1) {
            hTruePionAbsNp->Fill(ev.truthPrimaryVertexKE * 1000, fd_weight);
        } else if (ev.backgroundType == 6 || ev.backgroundType == 12) {
            hTruePionScatter->Fill(ev.truthPrimaryVertexKE * 1000, fd_weight);
        }

        //////////////////////////////
        // Classification algorithm //
        //////////////////////////////

        // Reject stuff with no data products
        if (!ev.WC2TPCMatch || ev.WC2TPCsize != 1) continue;

        ////////////////////////////////
        // Cylinder and TG track cuts //
        ////////////////////////////////

        int numSmallTracksInCylinder = 0;
        int numTracksInCylinder      = 0;
        int numTGTracks              = 0;
        for (int trk_idx = 0; trk_idx < ev.recoBeginX.size(); ++trk_idx) {
            if (trkWCtoTPCMatch[trk_idx]) continue;

            // Check if tracks is through-going
            if (
                !isWithinReducedVolume(ev.recoBeginX.at(trk_idx), ev.recoBeginY.at(trk_idx), ev.recoBeginZ.at(trk_idx)) &&
                !isWithinReducedVolume(ev.recoEndX.at(trk_idx), ev.recoEndY.at(trk_idx), ev.recoEndZ.at(trk_idx))
            ) numTGTracks++;

            double trackLength = sqrt(
                pow(ev.recoEndX.at(trk_idx) - ev.recoBeginX.at(trk_idx), 2) +
                pow(ev.recoEndY.at(trk_idx) - ev.recoBeginY.at(trk_idx), 2) +
                pow(ev.recoEndZ.at(trk_idx) - ev.recoBeginZ.at(trk_idx), 2)
            );

            bool startInCylinder = IsPointInsideTrackCylinder(
                wcX, wcY, wcZ,
                ev.recoBeginX.at(trk_idx), ev.recoBeginY.at(trk_idx), ev.recoBeginZ.at(trk_idx),
                CYLINDER_RADIUS
            );
            bool endInCylinder = IsPointInsideTrackCylinder(
                wcX, wcY, wcZ,
                ev.recoEndX.at(trk_idx), ev.recoEndY.at(trk_idx), ev.recoEndZ.at(trk_idx),
                CYLINDER_RADIUS
            );

            if (startInCylinder && endInCylinder) numTracksInCylinder++;

            if (
                startInCylinder && endInCylinder &&
                trackLength < CYLINDER_SMALL_TRACK
            ) {
                numSmallTracksInCylinder++;
            }
        }

        // Reject events with too many TG tracks
        if (numTGTracks > MAX_NUM_TG_TRACKS) continue;

        // Reject events with too many small tracks inside the cylinder
        if (numSmallTracksInCylinder > ALLOWED_CYLINDER_SMALL_TRACKS) continue;

        ///////////////////////
        // Primary track PID //
        ///////////////////////

        int totalCaloPoints = ev.wcMatchDEDX.size();
        int nRemoveOutliers = 2;
        int nRemoveEnds     = 3;
        int minPoints       = 5;

        // Get chi^2 fits, primary tracks are already checked for reversal in first module
        double pionChi2   = computeReducedChi2(  gPion, ev.wcMatchResR, ev.wcMatchDEDX, false, totalCaloPoints, nRemoveOutliers, nRemoveEnds);
        double MIPChi2    = computeReducedChi2(gMuonTG, ev.wcMatchResR, ev.wcMatchDEDX, false, totalCaloPoints, nRemoveOutliers, nRemoveEnds);
        double protonChi2 = computeReducedChi2(gProton, ev.wcMatchResR, ev.wcMatchDEDX, false, totalCaloPoints, nRemoveOutliers, nRemoveEnds);

        double minStitchedChi2 = std::numeric_limits<double>::max();
        int bestBreakPoint = -1;
        if (totalCaloPoints >= 4 * nRemoveEnds + 2 * nRemoveOutliers + 2 * minPoints) {
            for (int caloBreakPoint = 2 * nRemoveEnds + nRemoveOutliers + minPoints; caloBreakPoint < totalCaloPoints - (2 * nRemoveEnds + nRemoveOutliers + minPoints); ++caloBreakPoint) {
                std::vector<double> leftResR(ev.wcMatchResR.begin(), ev.wcMatchResR.begin() + caloBreakPoint);
                std::vector<double> leftDEDX(ev.wcMatchDEDX.begin(), ev.wcMatchDEDX.begin() + caloBreakPoint);

                std::vector<double> rightResR(ev.wcMatchResR.begin() + caloBreakPoint, ev.wcMatchResR.end());
                std::vector<double> rightDEDX(ev.wcMatchDEDX.begin() + caloBreakPoint, ev.wcMatchDEDX.end());

                // Shift right-hand side to fix r.r.
                for (int i = 0; i < rightResR.size(); ++i) {
                    rightResR[i] = rightResR[i] - rightResR[0];
                }

                double chi2LHS = computeReducedChi2(gProton, leftResR, leftDEDX, false, leftResR.size(), nRemoveOutliers, nRemoveEnds);
                double chi2RHS = computeReducedChi2(gMuonTG, rightResR, rightDEDX, false, rightResR.size(), nRemoveOutliers, nRemoveEnds);

                double totalChi2 = (chi2LHS * leftResR.size() + chi2RHS * rightResR.size()) / totalCaloPoints;
                
                if (totalChi2 < minStitchedChi2) {
                    minStitchedChi2 = totalChi2;
                    bestBreakPoint  = caloBreakPoint;
                }
            }
        }

        // Get smallest chi^2 value for primary track
        double minChi2 = std::min({minStitchedChi2, pionChi2, MIPChi2, protonChi2});

        // If primary track stitched, get break point, otherwise break point is end of track
        double breakPointX = ev.WC2TPCPrimaryEndX; 
        double breakPointY = ev.WC2TPCPrimaryEndY; 
        double breakPointZ = ev.WC2TPCPrimaryEndZ;
        if (minChi2 == minStitchedChi2) {
            breakPointX = ev.wcMatchXPos.at(bestBreakPoint);
            breakPointY = ev.wcMatchYPos.at(bestBreakPoint);
            breakPointZ = ev.wcMatchZPos.at(bestBreakPoint);
        }

        // At this point, we want to fill the incident kinetic energy histograms
        double WCKE             = TMath::Sqrt(ev.WCTrackMomentum * ev.WCTrackMomentum + PionMass * PionMass) - PionMass;
        double calculatedEnLoss = energyLossCalculation(ev.WC4PrimaryX, ev.trajectoryInitialMomentumX, isData);
        const double initialKE  = WCKE - calculatedEnLoss;

        double energyDeposited = 0;
        for (size_t iDep = 0; iDep < ev.wcMatchDEDX.size(); ++iDep) {
            // If we are past detected breaking point, exit loop
            if (ev.wcMatchZPos.at(iDep) > breakPointZ) break;

            // If larger than threshold, continue
            if (ev.wcMatchDEDX.at(iDep) > HIT_DEDX_THRESHOLD) continue;

            // Else, add to energy deposited so far
            energyDeposited += ev.wcMatchEDep.at(iDep);
            
            // Add to incident KE if inside reduced volume
            if (isWithinReducedVolume(ev.wcMatchXPos.at(iDep), ev.wcMatchYPos.at(iDep), ev.wcMatchZPos.at(iDep))) {
                hPionIncidentKE->Fill(initialKE - energyDeposited);
            }
        }
        double energyAtVertex = initialKE - energyDeposited;

        // Reject if vertex outside reduced volume
        if (!isWithinReducedVolume(breakPointX, breakPointY, breakPointZ)) continue;

        // Reject based on primary PID 
        if (minChi2 == pionChi2 || minChi2 == protonChi2) continue;

        /////////////////////////
        // Secondary track PID //
        /////////////////////////

        int secondaryTaggedPion   = 0;
        int secondaryTaggedProton = 0;
        int secondaryTaggedOther  = 0;

        int otherTaggedPion   = 0;
        int otherTaggedProton = 0;

        for (size_t trk_idx = 0; trk_idx < ev.recoBeginX.size(); ++trk_idx) {
            if (trkWCtoTPCMatch[trk_idx]) continue;

            // Have to re-check track ordering for stitched case
            double distanceFromStart = distance(
                ev.recoBeginX.at(trk_idx), breakPointX, 
                ev.recoBeginY.at(trk_idx), breakPointY,
                ev.recoBeginZ.at(trk_idx), breakPointZ
            );
            double distanceFromEnd = distance(
                ev.recoEndX.at(trk_idx), breakPointX, 
                ev.recoEndY.at(trk_idx), breakPointY,
                ev.recoEndZ.at(trk_idx), breakPointZ
            );
            double thisTrackLength = sqrt(
                pow(ev.recoBeginX.at(trk_idx) - ev.recoEndX.at(trk_idx), 2) +
                pow(ev.recoBeginY.at(trk_idx) - ev.recoEndY.at(trk_idx), 2) + 
                pow(ev.recoBeginZ.at(trk_idx) - ev.recoEndZ.at(trk_idx), 2)
            );

            if ((distanceFromStart < VERTEX_RADIUS || distanceFromEnd < VERTEX_RADIUS)) {
                std::vector<double> secondaryDEDX = ev.recoDEDX.at(trk_idx);
                std::vector<double> secondaryResR = ev.recoResR.at(trk_idx);

                bool secondaryReversed  = false;
                bool originallyReversed = ev.isTrackInverted.at(trk_idx);
                if (distanceFromEnd < distanceFromStart) secondaryReversed = true;
                // If it was originally reversed, we do not reverse again
                if (originallyReversed && secondaryReversed) secondaryReversed = false; 

                double pionChi2          = computeReducedChi2(gPion, secondaryResR, secondaryDEDX, secondaryReversed, secondaryDEDX.size(), nRemoveOutliers, nRemoveEnds);
                double protonChi2        = computeReducedChi2(gProton, secondaryResR, secondaryDEDX, secondaryReversed, secondaryDEDX.size(), nRemoveOutliers, nRemoveEnds);
                double secondaryMeanDEDX = meanDEDX(secondaryDEDX, secondaryReversed, MEAN_DEDX_NUM_TRAJ_POINTS);

                // First, try classifying track using chi^2 fits
                if ((pionChi2 < PION_CHI2_PION_VALUE) && (protonChi2 > PROTON_CHI2_PION_VALUE)) {
                    // Tagged as pion
                    secondaryTaggedPion++;
                } else if ((pionChi2 > PION_CHI2_PROTON_VALUE) && (protonChi2 < PROTON_CHI2_PROTON_VALUE)) {
                    // Tagged as proton
                    secondaryTaggedProton++;
                } else {
                    // Not tagged with chi^2, use mean dE/dx
                    secondaryTaggedOther++;
                    if (secondaryMeanDEDX <= MEAN_DEDX_THRESHOLD) {
                        otherTaggedPion++;
                    } else {
                        otherTaggedProton++;
                    }
                }
            }
        }    
        
        // For particles where we stitched, we also need to analyze the second part of the primary track
        // If the new secondary track looks like a pion, the stitching is non-sensical, and we should revert
        bool newSecondaryPion = false;
        if (minChi2 == minStitchedChi2) {
            std::vector<double> newSecondaryResR(ev.wcMatchResR.begin(), ev.wcMatchResR.begin() + bestBreakPoint);
            std::vector<double> newSecondaryDEDX(ev.wcMatchDEDX.begin(), ev.wcMatchDEDX.begin() + bestBreakPoint);

            double newPionChi2   = computeReducedChi2(gPion, newSecondaryResR, newSecondaryDEDX, false, newSecondaryDEDX.size(), nRemoveOutliers, nRemoveEnds);
            double newProtonChi2 = computeReducedChi2(gProton, newSecondaryResR, newSecondaryDEDX, false, newSecondaryDEDX.size(), nRemoveOutliers, nRemoveEnds);
            double newMeanDEDX   = meanDEDX(newSecondaryDEDX, false, MEAN_DEDX_NUM_TRAJ_POINTS);

            if ((newPionChi2 < PION_CHI2_PION_VALUE) && (newProtonChi2 > PROTON_CHI2_PION_VALUE)) {
                // Tagged as pion
                newSecondaryPion = true;
                secondaryTaggedPion++;
            } else if ((newPionChi2 > PION_CHI2_PROTON_VALUE) && (newProtonChi2 < PROTON_CHI2_PROTON_VALUE)) {
                // Tagged as proton
                secondaryTaggedProton++;
            } else {
                // Not tagged with chi^2, use mean dE/dx
                secondaryTaggedOther++;
                if (newMeanDEDX <= MEAN_DEDX_THRESHOLD) {
                    newSecondaryPion = true;
                    otherTaggedPion++;
                } else {
                    otherTaggedProton++;
                }
            }
        }

        int totalTaggedPions   = secondaryTaggedPion + otherTaggedPion;
        int totalTaggedProtons = secondaryTaggedProton + otherTaggedProton;

        // Fill scattering
        if (totalTaggedPions > 0) {
            if (totalTaggedPions > 1 || newSecondaryPion) continue;

            hPionScatter->Fill(energyAtVertex, fd_weight);
            continue;
        }

        // Fill abs Np
        if (totalTaggedProtons > 0) {
            hPionAbsNp->Fill(energyAtVertex, fd_weight);
            continue;
        }

        ////////////////////////////////////////
        // Cluster non-reconstructed hits cut //
        ////////////////////////////////////////

        // Get data for cut
        int numLargeClustersInduction  = 0;
        int numLargeClustersCollection = 0;

        int numClustersInduction  = 0;
        int numClustersCollection = 0;

        double largestClusterSizeInduction  = 0;
        double largestClusterSizeCollection = 0;

        for (int i = 0; i < hitClusters.size(); ++i) {
            double clusterSize = hitClusters[i].clusterSize;

            if (hitClusters[i].plane == 0) {
                if (clusterSize > largestClusterSizeInduction) largestClusterSizeInduction = clusterSize;
                numClustersInduction++;
            }
            else if (hitClusters[i].plane == 1) {
                if (clusterSize > largestClusterSizeCollection) largestClusterSizeCollection = clusterSize;
                numClustersCollection++;
            }

            if (clusterSize > LARGE_CLUSTER_THRESHOLD) {
                if (hitClusters[i].plane == 0) numLargeClustersInduction++;
                else if (hitClusters[i].plane == 1) numLargeClustersCollection++;
            }
        }

        // Fill abs 0p
        if (numClustersInduction < MAX_NUM_CLUSTERS_INDUCTION) {
            hPionAbs0p->Fill(energyAtVertex, fd_weight);
            continue;
        }
    }

    ////////////////////////
    // True cross-section //
    ////////////////////////

    TH1D* hSignal = new TH1D("hSignal", "hSignal;;", NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index   = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            double xsec =  XSEC_UNITS * (TotalEventsHistos[iSignal]->GetBinContent(iBin + 1) / hTrueIncidentKE->GetBinContent(iBin + 1));
            hSignal->SetBinContent(index + 1,  xsec);
        }
    }

    ////////////////////////
    // Reco cross-section //
    ////////////////////////

    TH1D* hMeasure = new TH1D("hMeasure", "hMeasure;;", NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            double xsec = XSEC_UNITS * (RecoSignals[iSignal]->GetBinContent(iBin + 1) / hPionIncidentKE->GetBinContent(iBin + 1));
            hMeasure->SetBinContent(index + 1, xsec);
        }
    }

    //////////////////
    // Save to file //
    //////////////////

    TFile* saveFile = new TFile(("/exp/lariat/app/users/epelaez/histos/fake_data/FD" + TString::Itoa(sample, 10) + ".root").Data(), "RECREATE");
    saveFile->cd();
    hSignal->Write("", TObject::kOverwrite);
    hMeasure->Write("", TObject::kOverwrite);
    hPionIncidentKE->Write("", TObject::kOverwrite);
    saveFile->Close();

    ///////////////////////
    // Create histograms //
    ///////////////////////

    std::vector<int> Colors = {
        kBlack,
        kBlue,
        kRed,
        kGreen,
        kOrange+1,
        kMagenta,
        kCyan+1,
        kViolet+1,
        kAzure+1,
        kPink+6
    };

    std::vector<std::vector<TH1*>> PlotGroups = {
        {hMeasure, hSignal}
    };

    std::vector<std::vector<TString>> PlotLabelGroups = {
        {"Reco", "True"}
    };

    std::vector<TString> PlotTitles = {
        "Sample"
    };

    std::vector<TString> XLabels = {
        "Bin"
    };

    std::vector<TString> YLabels = {
        "Cross-Section [barn]"
    };

    std::vector<bool> PlotStacked = {
        false
    };

    std::vector<std::vector<bool>> PlotsAsPoints = {
        {true, false}
    };

    printOneDPlots(
        SaveDir, FontStyle, TextSize,
        PlotGroups,
        Colors,
        PlotLabelGroups,
        PlotTitles,
        XLabels,
        YLabels,
        PlotStacked,
        PlotsAsPoints
    );
}