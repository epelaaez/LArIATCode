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

void BeamlineElCov() {
    // This module generates all relevant covariance matrices coming from beamline content variations

    // Set defaults
    gStyle->SetOptStat(0); // get rid of stats box
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    TH1::AddDirectory(false);
    gStyle->SetPalette(kRainBow);

    TGraph* gProton = new TGraph();
    TGraph* gPion   = new TGraph();
    TGraph* gMuonTG = new TGraph();

    // Number of universes for this systematic
    int NUM_UNIVERSES = 2; // +/- mu unc

    // Initialize points
    initializeProtonPoints(gProton);
    initializePionPoints(gPion);
    initializeMuonNoBraggPoints(gMuonTG);

    int FontStyle = 132;
    double TextSize = 0.06;
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis/figs/Systematics/";

    // Load files
    TChain* Chain = new TChain("anatree/anatree");
    Chain->Add("/exp/lariat/app/users/epelaez/files/anatree_60a/chunks/*.root");
    std::cout << "Files:   " << Chain->GetListOfFiles()->GetEntries() << std::endl;

    // Load weights
    TH1D* hWeightsFrontFace = dynamic_cast<TH1D*>(fWeights->Get(WEIGHTS_NAME));
    hWeightsFrontFace->SetDirectory(nullptr);
    fWeights->Close();

    ///////////////////
    // Load branches //
    ///////////////////

    // Run information
    int run, subrun, event; bool isData = false;
    Chain->SetBranchAddress("run", &run); 
    Chain->SetBranchAddress("subrun", &subrun); 
    Chain->SetBranchAddress("event", &event);

    // Track information
    int   ntracks_reco;                       Chain->SetBranchAddress("ntracks_reco",    &ntracks_reco);
    static float trkvtxx[kMaxTrack];                 Chain->SetBranchAddress("trkvtxx",          &trkvtxx);
    static float trkvtxy[kMaxTrack];                 Chain->SetBranchAddress("trkvtxy",          &trkvtxy);
    static float trkvtxz[kMaxTrack];                 Chain->SetBranchAddress("trkvtxz",          &trkvtxz);
    static float trkendx[kMaxTrack];                 Chain->SetBranchAddress("trkendx",          &trkendx);
    static float trkendy[kMaxTrack];                 Chain->SetBranchAddress("trkendy",          &trkendy);
    static float trkendz[kMaxTrack];                 Chain->SetBranchAddress("trkendz",          &trkendz);
    static int   trkWCtoTPCMatch[kMaxTrack];         Chain->SetBranchAddress("trkWCtoTPCMatch",  &trkWCtoTPCMatch);

    // Wire-chamber track information
    int   nwctrks;                              Chain->SetBranchAddress("nwctrks",           &nwctrks);
    static float wctrk_momentum[kMaxWCTracks];         Chain->SetBranchAddress("wctrk_momentum",    &wctrk_momentum);
    static float wctrk_theta[kMaxWCTracks];            Chain->SetBranchAddress("wctrk_theta",        &wctrk_theta);
    static float wctrk_phi[kMaxWCTracks];              Chain->SetBranchAddress("wctrk_phi",          &wctrk_phi);
    static int   wctrk_picky[kMaxWCTracks];            Chain->SetBranchAddress("wctrk_picky",        &wctrk_picky);
    static float WC4xPos[kMaxWCTracks];                 Chain->SetBranchAddress("WC4xPos",            &WC4xPos);

    // Calorimetry information
    static int   ntrkcalopts[kMaxTrack][2];                    Chain->SetBranchAddress("ntrkcalopts", &ntrkcalopts);
    static float trkdedx[kMaxTrack][2][kMaxTrackHits];         Chain->SetBranchAddress("trkdedx",     &trkdedx);
    static float trkrr[kMaxTrack][2][kMaxTrackHits];           Chain->SetBranchAddress("trkrr",       &trkrr);
    static float trkpitch[kMaxTrack][2][kMaxTrackHits];        Chain->SetBranchAddress("trkpitch",    &trkpitch);
    static float trkxyz[kMaxTrack][2][kMaxTrackHits][3];       Chain->SetBranchAddress("trkxyz",      &trkxyz);

    // Trajectory information for tracks
    int   nTrajPoint[kMaxTrack];                  Chain->SetBranchAddress("nTrajPoint", &nTrajPoint);
    static float trjPt_X[kMaxTrack][kMaxTrajHits];       Chain->SetBranchAddress("trjPt_X",    &trjPt_X);
    static float trjPt_Y[kMaxTrack][kMaxTrajHits];       Chain->SetBranchAddress("trjPt_Y",    &trjPt_Y);
    static float trjPt_Z[kMaxTrack][kMaxTrajHits];       Chain->SetBranchAddress("trjPt_Z",    &trjPt_Z);

    // Geant4 information for truth tracks
    int   no_primaries;                                   Chain->SetBranchAddress("no_primaries",        &no_primaries);
    int   geant_list_size;                                Chain->SetBranchAddress("geant_list_size",      &geant_list_size);
    static int   pdg[kMaxPrimaries];                             Chain->SetBranchAddress("pdg",                  &pdg);
    static float Mass[kMaxPrimaries];                            Chain->SetBranchAddress("Mass",                 &Mass);
    static float StartPointz[kMaxPrimaries];                     Chain->SetBranchAddress("StartPointz",          &StartPointz);
    static float Eng[kMaxPrimaries];                             Chain->SetBranchAddress("Eng",                  &Eng);
    static float Px[kMaxPrimaries];                              Chain->SetBranchAddress("Px",                   &Px);
    static float Py[kMaxPrimaries];                              Chain->SetBranchAddress("Py",                   &Py);
    static float Pz[kMaxPrimaries];                              Chain->SetBranchAddress("Pz",                   &Pz);
    static float EndPointx[kMaxPrimaries];                       Chain->SetBranchAddress("EndPointx",            &EndPointx);
    static float EndPointy[kMaxPrimaries];                       Chain->SetBranchAddress("EndPointy",            &EndPointy);
    static float EndPointz[kMaxPrimaries];                       Chain->SetBranchAddress("EndPointz",            &EndPointz);
    static float EndEng[kMaxPrimaries];                          Chain->SetBranchAddress("EndEng",               &EndEng);
    static float EndPx[kMaxPrimaries];                           Chain->SetBranchAddress("EndPx",                &EndPx);
    static float EndPy[kMaxPrimaries];                           Chain->SetBranchAddress("EndPy",                &EndPy);
    static float EndPz[kMaxPrimaries];                           Chain->SetBranchAddress("EndPz",                &EndPz);
    static int   Process[kMaxPrimaries];                         Chain->SetBranchAddress("Process",              &Process);
    static int   NumberDaughters[kMaxPrimaries];                 Chain->SetBranchAddress("NumberDaughters",      &NumberDaughters);
    static int   TrackId[kMaxPrimaries];                         Chain->SetBranchAddress("TrackId",              &TrackId);
    static int   Mother[kMaxPrimaries];                          Chain->SetBranchAddress("Mother",               &Mother);
    static int   process_primary[kMaxPrimaries];                 Chain->SetBranchAddress("process_primary",      &process_primary);
    std::vector<int>* InteractionPoint = nullptr;                Chain->SetBranchAddress("InteractionPoint",     &InteractionPoint);
    std::vector<int>* InteractionPointType = nullptr;            Chain->SetBranchAddress("InteractionPointType", &InteractionPointType);
    static int   NTrTrajPts[kMaxPrimaryPart];                    Chain->SetBranchAddress("NTrTrajPts",           &NTrTrajPts);
    static float MidPosX[kMaxPrimaryPart][kMaxTruePrimaryPts];   Chain->SetBranchAddress("MidPosX",              &MidPosX);
    static float MidPosY[kMaxPrimaryPart][kMaxTruePrimaryPts];   Chain->SetBranchAddress("MidPosY",              &MidPosY);
    static float MidPosZ[kMaxPrimaryPart][kMaxTruePrimaryPts];   Chain->SetBranchAddress("MidPosZ",              &MidPosZ);
    static float MidPx[kMaxPrimaryPart][kMaxTruePrimaryPts];     Chain->SetBranchAddress("MidPx",                &MidPx);
    static float MidPy[kMaxPrimaryPart][kMaxTruePrimaryPts];     Chain->SetBranchAddress("MidPy",                &MidPy);
    static float MidPz[kMaxPrimaryPart][kMaxTruePrimaryPts];     Chain->SetBranchAddress("MidPz",                &MidPz);

    // Information about wire plane hits
    int    nhits;                              Chain->SetBranchAddress("nhits",              &nhits);
    static int    hit_plane[kMaxHits];                Chain->SetBranchAddress("hit_plane",           hit_plane);
    static int    hit_channel[kMaxHits];              Chain->SetBranchAddress("hit_channel",         hit_channel);
    static int    hit_trkid[kMaxHits];               Chain->SetBranchAddress("hit_trkid",           hit_trkid);
    static float  hit_driftT[kMaxHits];              Chain->SetBranchAddress("hit_driftT",           hit_driftT);
    static float  hit_x[kMaxHits];                  Chain->SetBranchAddress("hit_x",                hit_x);
    static float  hit_y[kMaxHits];                  Chain->SetBranchAddress("hit_y",                hit_y);
    static float  hit_z[kMaxHits];                  Chain->SetBranchAddress("hit_z",                hit_z);

    // Simchannel information
    int maxTrackIDE;          Chain->SetBranchAddress("maxTrackIDE", &maxTrackIDE);
    static float IDEEnergy[kMaxIDE]; Chain->SetBranchAddress("IDEEnergy",   &IDEEnergy);
    static float IDEPos[kMaxIDE][3]; Chain->SetBranchAddress("IDEPos",      &IDEPos);

    ///////////////////////
    // Create histograms //
    ///////////////////////

    /////////////////
    // Incident KE //
    /////////////////

    // Selected incident pion
    std::vector<TH1D*> PionIncidentKEUnivs = MakeUniverseHists("hIncidentKE", "IncidentKE;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hPionIncidentKENom = new TH1D("hPionIncidentKENom", "hPionIncidentKENom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Selected incident pion that is true pion
    std::vector<TH1D*> PionIncidentKETrueUnivs = MakeUniverseHists("hIncidentKETrue", "IncidentKETrue;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hPionIncidentKETrueNom = new TH1D("hPionIncidentKETrueNom", "hPionIncidentKETrueNom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // True incident KE
    std::vector<TH1D*> TrueIncidentKEUnivs = MakeUniverseHists("hTrueIncidentKE", "TrueIncidentKE;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hTrueIncidentKENom = new TH1D("hTrueIncidentKENom", "hTrueIncidentKENom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Psi corrections
    std::vector<TH1D*> PsiIncidentUnivs = MakeUniverseHists("PsiIncidentUnivs", "PsiIncidentUnivs;Kinetic Energy [MeV];Factor", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hPsiIncidentNom = new TH1D("hPsiIncidentNom", "hPsiIncidentNom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // C corrections
    std::vector<TH1D*> CIncidentUnivs = MakeUniverseHists("CIncidentUnivs", "CIncidentUnivs;Kinetic Energy [MeV];Factor", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hCIncidentNom = new TH1D("hCIncidentNom", "hCIncidentNom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Corrected selected incident pion
    std::vector<TH1D*> PionIncidentCorrectedKEUnivs = MakeUniverseHists("hPionIncidentCorrectedKE", "PionIncidentCorrectedKE;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hPionIncidentCorrectedKENom = new TH1D("hPionIncidentCorrectedKENom", "hPionIncidentCorrectedKENom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    /////////////////
    // Selected KE // 
    /////////////////

    // Selected pion abs 0p for universes
    std::vector<TH1D*> PionAbs0pKEUnivs = MakeUniverseHists("hPionAbs0pKE", "PionAbs0pKE;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hPionAbs0pKENom = new TH1D("hPionAbs0pKENom", "hPionAbs0pKENom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Pion abs 0p estimated backgrounds
    std::vector<TH1D*> PionAbs0pKEBkgUnivs = MakeUniverseHists("hPionAbs0pKEBkg", "hPionAbs0pKEBkg;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hPionAbs0pKEBkgNom = new TH1D("hPionAbs0pKEBkgNom", "hPionAbs0pKEBkgNom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Selected pion abs Np for universes
    std::vector<TH1D*> PionAbsNpKEUnivs = MakeUniverseHists("hPionAbsNpKE", "PionAbsNpKE;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hPionAbsNpKENom = new TH1D("hPionAbsNpKENom", "hPionAbsNpKENom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Pion abs Np estimated backgrounds
    std::vector<TH1D*> PionAbsNpKEBkgUnivs = MakeUniverseHists("hPionAbsNpKEBkg", "hPionAbsNpKEBkg;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hPionAbsNpKEBkgNom = new TH1D("hPionAbsNpKEBkgNom", "hPionAbsNpKEBkgNom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Selected pion scattering for universes
    std::vector<TH1D*> PionScatterKEUnivs = MakeUniverseHists("hPionScatterKE", "PionScatterKE;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hPionScatterKENom = new TH1D("hPionScatterKENom", "hPionScatterKENom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Pion scattering estimated backgrounds
    std::vector<TH1D*> PionScatterKEBkgUnivs = MakeUniverseHists("hPionScatterKEBkg", "hPionScatterKEBkg;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hPionScatterKEBkgNom = new TH1D("hPionScatterKEBkgNom", "hPionScatterKEBkgNom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    std::vector<TH1D*> RecoSignals = {
        hPionAbs0pKENom, hPionAbsNpKENom, hPionScatterKENom
    };

    std::vector<TH1D*> RecoBkgSignals = {
        hPionAbs0pKEBkgNom, hPionAbsNpKEBkgNom, hPionScatterKEBkgNom
    };

    /////////////
    // True KE //
    /////////////

    // True abs 0p KE
    std::vector<TH1D*> TruePionAbs0pKEUnivs = MakeUniverseHists("hTruePionAbs0pKE", "TruePionAbs0pKE;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hTruePionAbs0pKENom = new TH1D("hTruePionAbs0pKENom", "hTruePionAbs0pKENom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // True abs Np KE
    std::vector<TH1D*> TruePionAbsNpKEUnivs = MakeUniverseHists("hTruePionAbsNpKE", "TruePionAbsNpKE;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hTruePionAbsNpKENom = new TH1D("hTruePionAbsNpKENom", "hTruePionAbsNpKENom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // True abs scattering KE
    std::vector<TH1D*> TruePionScatterKEUnivs = MakeUniverseHists("hTruePionScatterKE", "TruePionScatterKE;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hTruePionScatterKENom = new TH1D("hTruePionScatterKENom", "hTruePionScatterKENom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    std::vector<TH1D*> TotalEventsHistosNom = {
        hTruePionAbs0pKENom, hTruePionAbsNpKENom, hTruePionScatterKENom
    };

    /////////////////////////
    // For response matrix //
    /////////////////////////

    // For all universes
    std::vector<std::vector<std::vector<TH1D*>>> TrueAbs0pAsByBinUnivs = MakeUniverseHistBlock(
        "hTrueAbs0p", "True Abs0p (U=%d): As %d in true bin %d",
        NUM_UNIVERSES, NUM_BINS_KE, NUM_SIGNAL_TYPES,
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    std::vector<std::vector<std::vector<TH1D*>>> TrueAbsNpAsByBinUnivs = MakeUniverseHistBlock(
        "hTrueAbsNp", "True AbsNp (U=%d): As %d in true bin %d",
        NUM_UNIVERSES, NUM_BINS_KE, NUM_SIGNAL_TYPES,
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    std::vector<std::vector<std::vector<TH1D*>>> TrueScatterAsByBinUnivs = MakeUniverseHistBlock(
        "hTrueScatter", "True Scatter (U=%d): As %d in true bin %d",
        NUM_UNIVERSES, NUM_BINS_KE, NUM_SIGNAL_TYPES,
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );

    // Nominal
    std::vector<std::vector<TH1D*>> TrueAbs0pAsByBinNom;
    for (int iEnergyBin = 0; iEnergyBin < NUM_BINS_KE; ++iEnergyBin) {
        std::vector<TH1D*> TempVec;
        for (int iInt = 0; iInt < NUM_SIGNAL_TYPES; ++iInt) {
            TH1D* hTempHist = new TH1D(Form("hTrueAbs0p_%d_Bin_As_%d", iEnergyBin, iInt), Form("True Abs 0p KE As %d in bin %d", iInt, iEnergyBin), NUM_BINS_KE, ARRAY_KE_BINS.data());
            TempVec.push_back(hTempHist);
        }
        TrueAbs0pAsByBinNom.push_back(TempVec);
    }
    std::vector<std::vector<TH1D*>> TrueAbsNpAsByBinNom;
    for (int iEnergyBin = 0; iEnergyBin < NUM_BINS_KE; ++iEnergyBin) {
        std::vector<TH1D*> TempVec;
        for (int iInt = 0; iInt < NUM_SIGNAL_TYPES; ++iInt) {
            TH1D* hTempHist = new TH1D(Form("hTrueAbsNp_%d_Bin_As_%d", iEnergyBin, iInt), Form("True Abs Np KE As %d in bin %d", iInt, iEnergyBin), NUM_BINS_KE, ARRAY_KE_BINS.data());
            TempVec.push_back(hTempHist);
        }
        TrueAbsNpAsByBinNom.push_back(TempVec);
    }
    std::vector<std::vector<TH1D*>> TrueScatterAsByBinNom;
    for (int iEnergyBin = 0; iEnergyBin < NUM_BINS_KE; ++iEnergyBin) {
        std::vector<TH1D*> TempVec;
        for (int iInt = 0; iInt < NUM_SIGNAL_TYPES; ++iInt) {
            TH1D* hTempHist = new TH1D(Form("hTrueScatter_%d_Bin_As_%d", iEnergyBin, iInt), Form("True Scatter KE As %d in bin %d", iInt, iEnergyBin), NUM_BINS_KE, ARRAY_KE_BINS.data());
            TempVec.push_back(hTempHist);
        }
        TrueScatterAsByBinNom.push_back(TempVec);
    }

    std::vector<std::vector<std::vector<TH1D*>>> TrueRecoAsByBinNom = {
        TrueAbs0pAsByBinNom, TrueAbsNpAsByBinNom, TrueScatterAsByBinNom
    };

    //////////////////////
    // Loop over events //
    //////////////////////

    bool verbose = false;

    Long64_t i = 0;
    while (true) {
        if (SKIP_INDICES.count(i)) { i++; continue; }
        if (Chain->GetEntry(i++) <= 0) break;

        // Go faster
        // if (i > USE_NUM_EVENTS) break;

        // Reset variables
        if (verbose) std::cout << std::endl;
        if (verbose) std::cout << "=================================" << std::endl;
        if (verbose) std::cout << "Got tree entry: " << i << std::endl;
        EventVariables ev;
        if (verbose) std::cout << "Variables reset" << std::endl;
        if (verbose) std::cout << "=================================" << std::endl;

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
        for (size_t trk_idx = 0; trk_idx < std::min(ntracks_reco, kMaxTrack); ++trk_idx) {
            if (trkWCtoTPCMatch[trk_idx]) {
                // Grab position information
                ev.WC2TPCPrimaryBeginX = trkvtxx[trk_idx];
                ev.WC2TPCPrimaryBeginY = trkvtxy[trk_idx];
                ev.WC2TPCPrimaryBeginZ = trkvtxz[trk_idx];
                ev.WC2TPCPrimaryEndX   = trkendx[trk_idx];
                ev.WC2TPCPrimaryEndY   = trkendy[trk_idx];
                ev.WC2TPCPrimaryEndZ   = trkendz[trk_idx];

                // Grab calorimetry information (in collection plane)
                int npts_dedx = std::min(ntrkcalopts[trk_idx][1], kMaxTrackHits);
                ev.wcMatchResR.assign(trkrr[trk_idx][1], trkrr[trk_idx][1] + npts_dedx);
                ev.wcMatchDEDX.assign(trkdedx[trk_idx][1], trkdedx[trk_idx][1] + npts_dedx);

                for (size_t dep_idx = 0; dep_idx < std::min(ntrkcalopts[trk_idx][1], kMaxTrackHits); ++dep_idx) {
                    ev.wcMatchEDep.push_back(trkdedx[trk_idx][1][dep_idx] * trkpitch[trk_idx][1][dep_idx]);
                    ev.wcMatchXPos.push_back(trkxyz[trk_idx][1][dep_idx][0]);
                    ev.wcMatchYPos.push_back(trkxyz[trk_idx][1][dep_idx][1]);
                    ev.wcMatchZPos.push_back(trkxyz[trk_idx][1][dep_idx][2]);
                }

                // Get location information
                int npts_wc2tpc = std::min(nTrajPoint[trk_idx], kMaxTrajHits);
                ev.WC2TPCLocationsX.assign(trjPt_X[trk_idx], trjPt_X[trk_idx] + npts_wc2tpc);
                ev.WC2TPCLocationsY.assign(trjPt_Y[trk_idx], trjPt_Y[trk_idx] + npts_wc2tpc);
                ev.WC2TPCLocationsZ.assign(trjPt_Z[trk_idx], trjPt_Z[trk_idx] + npts_wc2tpc);

                // Set flag and index
                primaryTrackIdx = trk_idx;
                ev.WC2TPCMatch     = true;
                ev.WC2TPCsize++;
            }
        }

        if (verbose) std::cout << "Found WC2TPC match: " << ev.WC2TPCMatch  << std::endl;

        // Copy vertex and end for all tracks
        int npts_trk = std::min(ntracks_reco, kMaxTrack);
        ev.recoEndX.assign(trkendx,  trkendx  + npts_trk);
        ev.recoEndY.assign(trkendy,  trkendy  + npts_trk);
        ev.recoEndZ.assign(trkendz,  trkendz  + npts_trk);
        ev.recoBeginX.assign(trkvtxx, trkvtxx + npts_trk);
        ev.recoBeginY.assign(trkvtxy, trkvtxy + npts_trk);
        ev.recoBeginZ.assign(trkvtxz, trkvtxz + npts_trk);

        // Now, we want to loop through all tracks
        for (size_t trk_idx = 0; trk_idx < std::min(ntracks_reco, kMaxTrack); ++trk_idx) {
            // Grab calorimetry information
            ev.recoResR.push_back(std::vector<double>(trkrr[trk_idx][1], trkrr[trk_idx][1] + std::min(ntrkcalopts[trk_idx][1], kMaxTrackHits)));
            ev.recoDEDX.push_back(std::vector<double>(trkdedx[trk_idx][1], trkdedx[trk_idx][1] + std::min(ntrkcalopts[trk_idx][1], kMaxTrackHits)));

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
                            outPion                 = TVector3(Px[inner_idx], Py[inner_idx], Pz[inner_idx]);
                        }
                    }
                }
                if (verbose) std::cout << "Found primary daughters: " << found << std::endl;
                
                // Get initial momentum for energy loss
                ev.trajectoryInitialMomentumX = 1000 * MidPx[iPrimary][0];

                // Get trajectory information
                int npts_truth = std::min(NTrTrajPts[true_idx], kMaxTruePrimaryPts);
                ev.truthPrimaryLocationX.assign(MidPosX[iPrimary], MidPosX[iPrimary] + npts_truth);
                ev.truthPrimaryLocationY.assign(MidPosY[iPrimary], MidPosY[iPrimary] + npts_truth);
                ev.truthPrimaryLocationZ.assign(MidPosZ[iPrimary], MidPosZ[iPrimary] + npts_truth);

                ev.truthPrimaryVertexX = MidPosX[iPrimary][NTrTrajPts[true_idx] - 1];
                ev.truthPrimaryVertexY = MidPosY[iPrimary][NTrTrajPts[true_idx] - 1];
                ev.truthPrimaryVertexZ = MidPosZ[iPrimary][NTrTrajPts[true_idx] - 1];

                // Get interaction in trajectory
                float lastTPCX = -999, lastTPCY = -999, lastTPCZ = -999;
                ev.interactionInTrajectory = false;
                if (verbose) std::cout << "Looking at trajectory" << std::endl;
                if (InteractionPoint->at(0) != NTrTrajPts[true_idx] - 1) {
                    for (int i_point = 0; i_point < InteractionPoint->size(); ++i_point) {
                        if (i_point >= InteractionPointType->size()) continue;
                        if (!isWithinReducedVolume(
                            MidPosX[iPrimary][InteractionPoint->at(i_point)], 
                            MidPosY[iPrimary][InteractionPoint->at(i_point)], 
                            MidPosZ[iPrimary][InteractionPoint->at(i_point)]
                        )) continue;

                        if (verbose) std::cout << "   Point: " << InteractionPoint->at(i_point) << "   Interaction: " << InteractionPointType->at(i_point) << std::endl;

                        // "hadElastic" is type 3
                        if (InteractionPointType->at(i_point) == 3) {
                            // Get kinetic energy at this point
                            float KE = std::sqrt(
                                MidPx[iPrimary][InteractionPoint->at(i_point)] * MidPx[iPrimary][InteractionPoint->at(i_point)] + 
                                MidPy[iPrimary][InteractionPoint->at(i_point)] * MidPy[iPrimary][InteractionPoint->at(i_point)] + 
                                MidPz[iPrimary][InteractionPoint->at(i_point)] * MidPz[iPrimary][InteractionPoint->at(i_point)] +
                                Mass[true_idx] * Mass[true_idx]
                            ) - Mass[true_idx];

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
                                momAfter = TVector3(
                                    MidPx[iPrimary][InteractionPoint->at(i_point)], 
                                    MidPy[iPrimary][InteractionPoint->at(i_point)], 
                                    MidPz[iPrimary][InteractionPoint->at(i_point)]
                                );
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
                                momAfter = TVector3(
                                    MidPx[iPrimary][InteractionPoint->at(i_point)], 
                                    MidPy[iPrimary][InteractionPoint->at(i_point)], 
                                    MidPz[iPrimary][InteractionPoint->at(i_point)]
                                );
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
                if (firstTPCX == -999 || firstTPCY == -999 || firstTPCZ == -999) ev.validTrueIncidentKE = false;
                if (
                    firstTPCX == lastTPCX ||
                    firstTPCY == lastTPCY ||
                    firstTPCZ == lastTPCZ
                ) ev.validTrueIncidentKE = false;

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
                        for (int i = 0; i < std::min(maxTrackIDE, kMaxIDE); i++) {
                            if (IDEPos[i][2] < oldPos.Z())     continue;
                            if (IDEPos[i][2] > currentPos.Z()) continue;
                            currentDepEnergy += IDEEnergy[i];
                        }

                        // Skip tiny depositions
                        // if (currentDepEnergy / uniformDist < 0.1) continue;

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

        for (size_t i_hit = 0; i_hit < std::min(nhits, kMaxHits); ++i_hit) {
            ev.fHitPlane.push_back(hit_plane[i_hit]);
            ev.fHitT.push_back(SAMPLING_RATE * hit_driftT[i_hit]);
            ev.fHitX.push_back(SAMPLING_RATE * hit_driftT[i_hit] * DRIFT_VELOCITY);
            if (hit_channel[i_hit] < 240) {
                ev.fHitW.push_back(hit_channel[i_hit] * 0.4);
            } else {
                ev.fHitW.push_back((hit_channel[i_hit] - 240) * 0.4);
            }

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

            ev.truthPrimaryVertexKE = ev.trajectoryInteractionKE;
        }

        // Grab incident KE for secondary interactions along primary track
        for (int i_seg = 0; i_seg < ev.secondaryInteractionXPosition.size(); i_seg++) {
            TVector3 segStart;
            if (i_seg == 0) {
                segStart = TVector3(
                    primaryGeantEndX,
                    primaryGeantEndY,
                    primaryGeantEndZ
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
            double segKineticEnergy;
            if (i_seg == 0) {
                segKineticEnergy = ev.truthPrimaryVertexKE * 1000;
            } else {
                segKineticEnergy = ev.secondaryInteractionInteractingKE.at(i_seg-1) * 1000;
            }
            for (auto it = std::next(orderedUniformTrjPtsSecondary.begin()), old_it = orderedUniformTrjPtsSecondary.begin(); it != orderedUniformTrjPtsSecondary.end(); it++, old_it++) {
                auto oldPos     = old_it->second;
                auto currentPos = it->second;
                double uniformDist = (currentPos - oldPos).Mag();

                double currentDepEnergy = 0.;
                for (int i = 0; i < std::min(maxTrackIDE, kMaxIDE); i++) {
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

        // Setup hits in tracks
        std::unordered_set<int> hitsInTracks(ev.hitRecoAsTrackKey.begin(), ev.hitRecoAsTrackKey.end());

        // Extend reco cylinder
        removeRepeatedPoints(&ev.WC2TPCLocationsX, &ev.WC2TPCLocationsY, &ev.WC2TPCLocationsZ);
        std::vector<double> wcX(ev.WC2TPCLocationsX);
        std::vector<double> wcY(ev.WC2TPCLocationsY);
        std::vector<double> wcZ(ev.WC2TPCLocationsZ);

        // Get direction to end cylinder
        int numPoints = wcX.size();
        int numTail   = std::min(10, numPoints - 1);
        std::vector<std::vector<double>> points;
        for (int j = numPoints - numTail - 1; j < numPoints; ++j) {
            points.push_back({
                wcX.at(j),
                wcY.at(j),
                wcZ.at(j)
            });
        }

        std::vector<double> avgDir(3, 0);
        if (numTail > 0) {
            avgDir = getAverageDir(points);

            // Extrapolate track to end
            double scale = (maxZ - points.back()[2]) / avgDir[2];
            wcX.push_back(points.back()[0] + scale * avgDir[0]);
            wcY.push_back(points.back()[1] + scale * avgDir[1]);
            wcZ.push_back(points.back()[2] + scale * avgDir[2]);
        }

        //////////////////////////////////
        // Look at unreconstructed hits //
        //////////////////////////////////

        std::vector<int> candidateHits;
        for (size_t iHit = 0; iHit < std::min(nhits, kMaxHits); ++iHit) {
            double hitX     = ev.fHitX.at(iHit);
            double hitW     = ev.fHitW.at(iHit);
            int    hitPlane = ev.fHitPlane.at(iHit);

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
                // Skip hits already in tracks
                if (hitsInTracks.count(iHit) > 0) continue;
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
            
            for (int iAllHit = 0; iAllHit < std::min(nhits, kMaxHits); ++iAllHit) {
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

        ///////////////////////////////
        // Correct signal definition //
        ///////////////////////////////

        // Scattering only if degree > THRESHOLD and energy > THRESHOLD
        double scatteringAngle  = -9999;
        double scatteringEnergy = -9999;
        
        // Modify scatterings
        if (ev.backgroundType == 12 || ev.backgroundType == 6) {
            if (ev.backgroundType == 12) {
                scatteringAngle  = ev.trajectoryInteractionAngle;
                scatteringEnergy = ev.trajectoryInteractionKE;
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

        ///////////////////////////////////
        // Get front face KE and reweigh //
        ///////////////////////////////////

        // At this point, we want to fill the incident kinetic energy histograms
        double WCKE             = TMath::Sqrt(ev.WCTrackMomentum * ev.WCTrackMomentum + PionMass * PionMass) - PionMass;
        double calculatedEnLoss = energyLossCalculation(ev.WC4PrimaryX, ev.trajectoryInitialMomentumX, isData);
        const double initialKE  = WCKE - calculatedEnLoss;

        ev.weight *= GetKEWeight(hWeightsFrontFace, initialKE);

        ///////////////////////////
        // Get truth information //
        ///////////////////////////

        if (ev.trueIncidentKEContributions.size() == 0) ev.truthPrimaryVertexKE = 0;
        else ev.truthPrimaryVertexKE = ev.trueIncidentKEContributions.back() / 1000.0;

        int TrueEnergyBin = getBin(ev.truthPrimaryVertexKE * 1000, ARRAY_KE_BINS);

        if (ev.validTrueIncidentKE) {
            for (double x : ev.trueIncidentKEContributions) {
                // Add to nominal
                hTrueIncidentKENom->Fill(x, ev.weight);

                // Add to universes. For generator systematics, we do not vary the 
                // incident flux, just the interacting flux
                for (int iUniv = 0; iUniv < NUM_UNIVERSES; ++iUniv) {
                    TrueIncidentKEUnivs[iUniv]->Fill(x, ev.weight);
                }
            }
        }

        if (ev.backgroundType == 0) {
            hTruePionAbs0pKENom->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
            for (int iUniv = 0; iUniv < NUM_UNIVERSES; ++iUniv) TruePionAbs0pKEUnivs[iUniv]->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
        } else if (ev.backgroundType == 1) {
            hTruePionAbsNpKENom->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
            for (int iUniv = 0; iUniv < NUM_UNIVERSES; ++iUniv) TruePionAbsNpKEUnivs[iUniv]->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
        } else if (ev.backgroundType == 6 || ev.backgroundType == 12) {
            hTruePionScatterKENom->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
            for (int iUniv = 0; iUniv < NUM_UNIVERSES; ++iUniv) TruePionScatterKEUnivs[iUniv]->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
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
                &wcX, &wcY, &wcZ,
                ev.recoBeginX.at(trk_idx), ev.recoBeginY.at(trk_idx), ev.recoBeginZ.at(trk_idx),
                CYLINDER_RADIUS
            );
            bool endInCylinder = IsPointInsideTrackCylinder(
                &wcX, &wcY, &wcZ,
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

        // Check these are in order
        if (ev.wcMatchZPos.size() > 1 && ev.wcMatchZPos.front() > ev.wcMatchZPos.back()) {
            std::reverse(ev.wcMatchZPos.begin(), ev.wcMatchZPos.end());
            std::reverse(ev.wcMatchDEDX.begin(), ev.wcMatchDEDX.end());
            std::reverse(ev.wcMatchEDep.begin(), ev.wcMatchEDep.end());
            std::reverse(ev.wcMatchXPos.begin(), ev.wcMatchXPos.end());
            std::reverse(ev.wcMatchYPos.begin(), ev.wcMatchYPos.end());
        }

        double energyDeposited  = 0;
        for (size_t iDep = 0; iDep < ev.wcMatchDEDX.size(); ++iDep) {
            // If we are past detected breaking point, exit loop
            if (ev.wcMatchZPos.at(iDep) > breakPointZ) break;

            // If larger than threshold, continue
            if (ev.wcMatchDEDX.at(iDep) > HIT_DEDX_THRESHOLD) continue;

            // Else, add to energy deposited so far
            energyDeposited += ev.wcMatchEDep.at(iDep);

            // Add to incident KE if inside reduced volume
            if (isWithinReducedVolume(ev.wcMatchXPos.at(iDep), ev.wcMatchYPos.at(iDep), ev.wcMatchZPos.at(iDep))) {
                hPionIncidentKENom->Fill(initialKE - energyDeposited, ev.weight);
                PionIncidentKEUnivs[0]->Fill(initialKE - energyDeposited, ev.weight);
                PionIncidentKEUnivs[1]->Fill(initialKE - energyDeposited, ev.weight);

                if (ev.truthPrimaryPDG == 11) {
                    PionIncidentKEUnivs[0]->Fill(initialKE - energyDeposited, -ELECTRON_COMP_UNC * ev.weight);
                    PionIncidentKEUnivs[1]->Fill(initialKE - energyDeposited, +ELECTRON_COMP_UNC * ev.weight);
                }

                if (ev.truthPrimaryPDG == -211) {
                    hPionIncidentKETrueNom->Fill(initialKE - energyDeposited, ev.weight);
                    PionIncidentKETrueUnivs[0]->Fill(initialKE - energyDeposited, ev.weight);
                    PionIncidentKETrueUnivs[1]->Fill(initialKE - energyDeposited, ev.weight);
                }
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
                std::pow(ev.recoBeginX.at(trk_idx) - ev.recoEndX.at(trk_idx), 2) +
                std::pow(ev.recoBeginY.at(trk_idx) - ev.recoEndY.at(trk_idx), 2) + 
                std::pow(ev.recoBeginZ.at(trk_idx) - ev.recoEndZ.at(trk_idx), 2)
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

        if (totalTaggedPions > 0) {
            if (totalTaggedPions > 1 || newSecondaryPion) continue;

            // Add weights to scatter
            hPionScatterKENom->Fill(energyAtVertex, ev.weight);

            PionScatterKEUnivs[0]->Fill(energyAtVertex, ev.weight);
            PionScatterKEUnivs[1]->Fill(energyAtVertex, ev.weight);
            if (ev.truthPrimaryPDG == 11) {
                PionScatterKEUnivs[0]->Fill(energyAtVertex, -ELECTRON_COMP_UNC * ev.weight); 
                PionScatterKEUnivs[1]->Fill(energyAtVertex, +ELECTRON_COMP_UNC * ev.weight);
            }

            if (ev.backgroundType == 0) {
                if (TrueEnergyBin != -1) TrueAbs0pAsByBinNom.at(TrueEnergyBin).at(2)->Fill(energyAtVertex, ev.weight);
            } else if (ev.backgroundType == 1) {
                if (TrueEnergyBin != -1) TrueAbsNpAsByBinNom.at(TrueEnergyBin).at(2)->Fill(energyAtVertex, ev.weight);
            } else if (ev.backgroundType == 6 || ev.backgroundType == 12) {
                if (TrueEnergyBin != -1) TrueScatterAsByBinNom.at(TrueEnergyBin).at(2)->Fill(energyAtVertex, ev.weight);
            } else {
                hPionScatterKEBkgNom->Fill(energyAtVertex, ev.weight);

                PionScatterKEBkgUnivs[0]->Fill(energyAtVertex, ev.weight);
                PionScatterKEBkgUnivs[1]->Fill(energyAtVertex, ev.weight);
                if (ev.truthPrimaryPDG == 11) {
                    PionScatterKEBkgUnivs[0]->Fill(energyAtVertex, -ELECTRON_COMP_UNC * ev.weight);
                    PionScatterKEBkgUnivs[1]->Fill(energyAtVertex, +ELECTRON_COMP_UNC * ev.weight);
                }
            }

            continue;
        }

        if (totalTaggedProtons > 0) {
            // Add weights to abs Np
            hPionAbsNpKENom->Fill(energyAtVertex, ev.weight);

            PionAbsNpKEUnivs[0]->Fill(energyAtVertex, ev.weight);
            PionAbsNpKEUnivs[1]->Fill(energyAtVertex, ev.weight);
            if (ev.truthPrimaryPDG == 11) {
                PionAbsNpKEUnivs[0]->Fill(energyAtVertex, -ELECTRON_COMP_UNC * ev.weight); 
                PionAbsNpKEUnivs[1]->Fill(energyAtVertex, +ELECTRON_COMP_UNC * ev.weight);
            }

            if (ev.backgroundType == 0) {
                if (TrueEnergyBin != -1) TrueAbs0pAsByBinNom.at(TrueEnergyBin).at(1)->Fill(energyAtVertex, ev.weight);
            } else if (ev.backgroundType == 1) {
                if (TrueEnergyBin != -1) TrueAbsNpAsByBinNom.at(TrueEnergyBin).at(1)->Fill(energyAtVertex, ev.weight);
            } else if (ev.backgroundType == 6 || ev.backgroundType == 12) {
                if (TrueEnergyBin != -1) TrueScatterAsByBinNom.at(TrueEnergyBin).at(1)->Fill(energyAtVertex, ev.weight);
            } else {
                hPionAbsNpKEBkgNom->Fill(energyAtVertex, ev.weight);

                PionAbsNpKEBkgUnivs[0]->Fill(energyAtVertex, ev.weight);
                PionAbsNpKEBkgUnivs[1]->Fill(energyAtVertex, ev.weight);
                if (ev.truthPrimaryPDG == 11) {
                    PionAbsNpKEBkgUnivs[0]->Fill(energyAtVertex, -ELECTRON_COMP_UNC * ev.weight);
                    PionAbsNpKEBkgUnivs[1]->Fill(energyAtVertex, +ELECTRON_COMP_UNC * ev.weight);
                }
            }
            
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

        if (
            // numClustersInduction < MAX_NUM_CLUSTERS_INDUCTION
            numLargeClustersInduction < MAX_NUM_LARGE_CLUSTERS_INDUCTION
        ) {
            // Add weights to abs 0p
            hPionAbs0pKENom->Fill(energyAtVertex, ev.weight);

            PionAbs0pKEUnivs[0]->Fill(energyAtVertex, ev.weight);
            PionAbs0pKEUnivs[1]->Fill(energyAtVertex, ev.weight);
            if (ev.truthPrimaryPDG == 11) {
                PionAbs0pKEUnivs[0]->Fill(energyAtVertex, -ELECTRON_COMP_UNC * ev.weight); 
                PionAbs0pKEUnivs[1]->Fill(energyAtVertex, +ELECTRON_COMP_UNC * ev.weight);
            }

            if (ev.backgroundType == 0) {
                if (TrueEnergyBin != -1) TrueAbs0pAsByBinNom.at(TrueEnergyBin).at(0)->Fill(energyAtVertex, ev.weight);
            } else if (ev.backgroundType == 1) {
                if (TrueEnergyBin != -1) TrueAbsNpAsByBinNom.at(TrueEnergyBin).at(0)->Fill(energyAtVertex, ev.weight);
            } else if (ev.backgroundType == 6 || ev.backgroundType == 12) {
                if (TrueEnergyBin != -1) TrueScatterAsByBinNom.at(TrueEnergyBin).at(0)->Fill(energyAtVertex, ev.weight);
            } else {
                hPionAbs0pKEBkgNom->Fill(energyAtVertex, ev.weight);

                PionAbs0pKEBkgUnivs[0]->Fill(energyAtVertex, ev.weight);
                PionAbs0pKEBkgUnivs[1]->Fill(energyAtVertex, ev.weight);
                if (ev.truthPrimaryPDG == 11) {
                    PionAbs0pKEBkgUnivs[0]->Fill(energyAtVertex, -ELECTRON_COMP_UNC * ev.weight);
                    PionAbs0pKEBkgUnivs[1]->Fill(energyAtVertex, +ELECTRON_COMP_UNC * ev.weight);
                }
            }

            continue;
        }
    }

    //////////////////////////////
    // Cross-section extraction //
    //////////////////////////////

    ///////////////////
    // Incident flux //
    ///////////////////
    
    // Get nominal incident flux and corrections
    for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
        double psi = hPionIncidentKETrueNom->GetBinContent(iBin + 1) / hTrueIncidentKENom->GetBinContent(iBin + 1);
        double   c = hPionIncidentKETrueNom->GetBinContent(iBin + 1) / hPionIncidentKENom->GetBinContent(iBin + 1);

        hPsiIncidentNom->SetBinContent(iBin + 1, psi); hCIncidentNom->SetBinContent(iBin + 1, c);
        hPionIncidentCorrectedKENom->SetBinContent(iBin + 1, (c / psi) * hPionIncidentKENom->GetBinContent(iBin + 1));
    }

    // Get incident flux and corrections for all universes
    for (int iUniv = 0; iUniv < NUM_UNIVERSES; ++iUniv) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            double psi = PionIncidentKETrueUnivs[iUniv]->GetBinContent(iBin + 1) / TrueIncidentKEUnivs[iUniv]->GetBinContent(iBin + 1);
            double   c = PionIncidentKETrueUnivs[iUniv]->GetBinContent(iBin + 1) / PionIncidentKEUnivs[iUniv]->GetBinContent(iBin + 1);

            PsiIncidentUnivs[iUniv]->SetBinContent(iBin + 1, psi); CIncidentUnivs[iUniv]->SetBinContent(iBin + 1, c);
            PionIncidentCorrectedKEUnivs[iUniv]->SetBinContent(iBin + 1, (c / psi) * PionIncidentKEUnivs[iUniv]->GetBinContent(iBin + 1));
        }
    }

    ///////////////////////
    // Response matrices //
    ///////////////////////

    // Get nominal response matrix
    TH2D* hResponseMatrixNominal = new TH2D(
        "hResponseMatrix", "hResponseMatrix;Reco (j, #beta);True (i, #alpha)",
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE,
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    );
    GetResponseMatrix(
        NUM_SIGNAL_TYPES, NUM_BINS_KE,
        TotalEventsHistosNom,
        TrueRecoAsByBinNom,
        hPionIncidentKENom,
        hTrueIncidentKENom,
        hResponseMatrixNominal
    );

    // Compute response matrix for all universes
    std::vector<TH2D*> ResponseMatricesUnivs;
    for (int iUniv = 0; iUniv < NUM_UNIVERSES; ++iUniv) {
        TH2D* hResponseMatrixUniv = new TH2D(
            "hResponseMatrix", "hResponseMatrix;Reco (j, #beta);True (i, #alpha)",
            NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE,
            NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
        );
        std::vector<TH1D*> TotalEventsUniv = {hTruePionAbs0pKENom, hTruePionAbsNpKENom, hTruePionScatterKENom};
        std::vector<std::vector<std::vector<TH1D*>>> TrueRecoAsByBinUniv = {TrueAbs0pAsByBinNom, TrueAbsNpAsByBinNom, TrueScatterAsByBinNom};
        TH1D* RecoIncidentFluxUniv = PionIncidentKEUnivs[iUniv];
        TH1D* TrueIncidentFluxUniv = hTrueIncidentKENom;
        
        GetResponseMatrix(
            NUM_SIGNAL_TYPES, NUM_BINS_KE,
            TotalEventsUniv,
            TrueRecoAsByBinUniv,
            RecoIncidentFluxUniv,
            TrueIncidentFluxUniv,
            hResponseMatrixUniv
        );

        ResponseMatricesUnivs.push_back(hResponseMatrixUniv);
    }

    //////////////////////////
    // Signal cross-section //
    //////////////////////////

    // Get nominal signal vector
    TH1D* hSignalNominal = new TH1D("hSignalNominal", "hSignalNominal;;", NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            double xsec =  XSEC_UNITS * (TotalEventsHistosNom[iSignal]->GetBinContent(iBin + 1) / hTrueIncidentKENom->GetBinContent(iBin + 1));
            hSignalNominal->SetBinContent(index + 1,  xsec); TotalEventsHistosNom[iSignal]->SetBinContent(iBin + 1,  xsec);
        }
    }

    //////////////////////////////
    // Background cross-section //
    //////////////////////////////

    // Get nominal estimated background vector
    TH1D* hBackgroundNominal = new TH1D("hBackgroundNominal", "hBackgroundNominal;;", NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            double xsec = XSEC_UNITS * (RecoBkgSignals[iSignal]->GetBinContent(iBin + 1) / hPionIncidentKENom->GetBinContent(iBin + 1));
            hBackgroundNominal->SetBinContent(index + 1, xsec); RecoBkgSignals[iSignal]->SetBinContent(iBin + 1, xsec);
        }
    }

    // Compose background vectors for all universes
    std::vector<TH1D*> BackgroundVectorUnivs;
    BackgroundVectorUnivs.reserve(NUM_UNIVERSES);
    for (int iUniv = 0; iUniv < NUM_UNIVERSES; ++iUniv) {
        TH1D* hBackground = new TH1D(Form("hBackground_Univ%d", iUniv), Form("Background Vector for Universe %d", iUniv), NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE);
        for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
            for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
                int index = flattenIndex(iSignal, iBin, NUM_BINS_KE);

                double content = 0;
                if (iSignal == 0) {
                    content = XSEC_UNITS * (PionAbs0pKEBkgUnivs[iUniv]->GetBinContent(iBin + 1) / PionIncidentKEUnivs[iUniv]->GetBinContent(iBin + 1));
                    PionAbs0pKEBkgUnivs[iUniv]->SetBinContent(iBin + 1, content);
                } else if (iSignal == 1) {
                    content = XSEC_UNITS * (PionAbsNpKEBkgUnivs[iUniv]->GetBinContent(iBin + 1) / PionIncidentKEUnivs[iUniv]->GetBinContent(iBin + 1));
                    PionAbsNpKEBkgUnivs[iUniv]->SetBinContent(iBin + 1, content);
                } else if (iSignal == 2) {
                    content = XSEC_UNITS * (PionScatterKEBkgUnivs[iUniv]->GetBinContent(iBin + 1) / PionIncidentKEUnivs[iUniv]->GetBinContent(iBin + 1));
                    PionScatterKEBkgUnivs[iUniv]->SetBinContent(iBin + 1, content);
                }

                hBackground->SetBinContent(index + 1, content);
                hBackground->SetBinError(index + 1, 0.0);
            }
        }
        BackgroundVectorUnivs.push_back(hBackground);
    }

    ////////////////////////////
    // Measured cross-section //
    ////////////////////////////

    // Get nominal measure vector
    TH1D* hMeasureNominal = new TH1D("hMeasureVectorNominal", "hMeasureVectorNominal;;", NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            double xsec = XSEC_UNITS * (RecoSignals[iSignal]->GetBinContent(iBin + 1) / hPionIncidentKENom->GetBinContent(iBin + 1));
            hMeasureNominal->SetBinContent(index + 1, xsec); RecoSignals[iSignal]->SetBinContent(iBin + 1, xsec);
        }
    }

    // Get nominal true cross-section in vector form
    TVectorD SignalNom(NUM_SIGNAL_TYPES * NUM_BINS_KE); H2V(hSignalNominal, SignalNom);

    // Compose measured vectors for all universes
    std::vector<TH1D*> MeasureVectorUnivs;
    MeasureVectorUnivs.reserve(NUM_UNIVERSES);
    for (int iUniv = 0; iUniv < NUM_UNIVERSES; ++iUniv) {
        // One TH1D holding all signals flattened
        TH1D* hMeasure = new TH1D(Form("hMeasure_Univ%d", iUniv), Form("Measured Vector for Universe %d", iUniv), NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE);

        // Compute M = R * S + B with S fixed
        TMatrixD R(NUM_SIGNAL_TYPES * NUM_BINS_KE, NUM_SIGNAL_TYPES * NUM_BINS_KE); H2M(static_cast<const TH2D*>(ResponseMatricesUnivs[iUniv]), R, kTRUE);
        TVectorD B(NUM_SIGNAL_TYPES * NUM_BINS_KE); H2V(BackgroundVectorUnivs[iUniv], B);

        TVectorD M = (R * SignalNom) + B;

        // Fill it with flattened (signal, bin) content
        for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
            for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
                int index     = flattenIndex(iSignal, iBin, NUM_BINS_KE);

                // double content = 0;
                if (iSignal == 0) {
                    // content = XSEC_UNITS * (PionAbs0pKEUnivs[iUniv]->GetBinContent(iBin + 1) / PionIncidentCorrectedKEUnivs[iUniv]->GetBinContent(iBin + 1));
                    PionAbs0pKEUnivs[iUniv]->SetBinContent(iBin + 1, M(index));
                    // std::cout << "measured: " << content << " vs. vector: " << M(index) << std::endl;
                } else if (iSignal == 1) {
                    // content = XSEC_UNITS * (PionAbsNpKEUnivs[iUniv]->GetBinContent(iBin + 1) / PionIncidentCorrectedKEUnivs[iUniv]->GetBinContent(iBin + 1));
                    PionAbsNpKEUnivs[iUniv]->SetBinContent(iBin + 1, M(index));
                    // std::cout << "measured: " << content << " vs. vector: " << M(index) << std::endl;
                } else if (iSignal == 2) {
                    // content = XSEC_UNITS * (PionScatterKEUnivs[iUniv]->GetBinContent(iBin + 1) / PionIncidentCorrectedKEUnivs[iUniv]->GetBinContent(iBin + 1));
                    PionScatterKEUnivs[iUniv]->SetBinContent(iBin + 1, M(index));
                    // std::cout << "measured: " << content << " vs. vector: " << M(index) << std::endl;
                }
                
                hMeasure->SetBinContent(index + 1, M(index));
                hMeasure->SetBinError(index + 1, 0.0);
            }
        }
        MeasureVectorUnivs.push_back(hMeasure);
    }

    // Use TMatrices/Vectors to compute response equation
    TMatrixD ResponseNom(NUM_SIGNAL_TYPES * NUM_BINS_KE, NUM_SIGNAL_TYPES * NUM_BINS_KE); H2M(static_cast<const TH2D*>(hResponseMatrixNominal), ResponseNom, kTRUE);
    TVectorD BackgroundNom(NUM_SIGNAL_TYPES * NUM_BINS_KE); H2V(hBackgroundNominal, BackgroundNom);

    // Get measured vector by "folding" true vector
    TVectorD MeasureNom = (ResponseNom * SignalNom) + BackgroundNom;

    // Sanity check, for CV, these two should match
    std::cout << "Cross-section: " << std::endl;
    TH1D* hCrossSectionNominal = new TH1D("hCrossSectionNominal", "hCrossSectionNominal;;", NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            std::cout << "  Reco: " << index << " " << hMeasureNominal->GetBinContent(index + 1) << std::endl;
            std::cout << "  Fold: " << index << " " << MeasureNom(index) << std::endl;
            std::cout << std::endl;
        }
    }

    /////////////////////
    // Get covariances //
    /////////////////////

    // Meausured vector
    TH2D* hMeasureCovMatrix = new TH2D(
        "hMeasureCovMatrix", "hMeasureCovMatrix;Reco (j, #beta);True (i, #alpha)",
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE,
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    );
    TH2D* hMeasureFracCovMatrix = new TH2D(
        "hMeasureFracCovMatrix", "hMeasureFracCovMatrix;Reco (j, #beta);True (i, #alpha)",
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE,
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    );
    TH2D* hMeasureCorrMatrix = new TH2D(
        "hMeasureCorrMatrix", "hMeasureCorrMatrix;Reco (j, #beta);True (i, #alpha)",
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE,
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    );
    GetCovMatrix(hMeasureNominal, MeasureVectorUnivs, hMeasureCovMatrix);
    GetFracCovAndCorrMatrix(hMeasureNominal, hMeasureCovMatrix, hMeasureFracCovMatrix, hMeasureCorrMatrix);
    DrawHistosWithErrorBands(
        hMeasureNominal,
        MeasureVectorUnivs,
        SaveDir,
        "BeamlineEl",
        "Measure",
        TextSize,
        FontStyle,
        "Measured Cross Section",
        "Bin number",
        "Cross Section [barn]"
    );

    // Abs 0p
    TH2D* hAbs0pCovMatrix = new TH2D(
        "hAbs0pCovMatrix", "hAbs0pCovMatrix;Kinetic Energy [MeV];Kinetic Energy [MeV]",
        NUM_BINS_KE, ARRAY_KE_BINS.data(),
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    TH2D* hAbs0pFracCovMatrix = new TH2D(
        "hAbs0pFracCovMatrix", "hAbs0pFracCovMatrix;Kinetic Energy [MeV];Kinetic Energy [MeV]",
        NUM_BINS_KE, ARRAY_KE_BINS.data(),
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    TH2D* hAbs0pCorrMatrix = new TH2D(
        "hAbs0pCorrMatrix", "hAbs0pCorrMatrix;Kinetic Energy [MeV];Kinetic Energy [MeV]",
        NUM_BINS_KE, ARRAY_KE_BINS.data(),
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    GetCovMatrix(hPionAbs0pKENom, PionAbs0pKEUnivs, hAbs0pCovMatrix);
    GetFracCovAndCorrMatrix(hPionAbs0pKENom, hAbs0pCovMatrix, hAbs0pFracCovMatrix, hAbs0pCorrMatrix);
    DrawHistosWithErrorBands(
        hPionAbs0pKENom, 
        PionAbs0pKEUnivs, 
        SaveDir, 
        "BeamlineEl", 
        "Abs0p",
        TextSize,
        FontStyle,
        "Abs 0p Cross Section",
        "Kinetic Energy [MeV]",
        "Cross Section [barn]"
    );

    // Abs Np
    TH2D* hAbsNpCovMatrix = new TH2D(
        "hAbsNpCovMatrix", "hAbsNpCovMatrix;Kinetic Energy [MeV];Kinetic Energy [MeV]",
        NUM_BINS_KE, ARRAY_KE_BINS.data(),
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    TH2D* hAbsNpFracCovMatrix = new TH2D(
        "hAbsNpFracCovMatrix", "hAbsNpFracCovMatrix;Kinetic Energy [MeV];Kinetic Energy [MeV]",
        NUM_BINS_KE, ARRAY_KE_BINS.data(),
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    TH2D* hAbsNpCorrMatrix = new TH2D(
        "hAbsNpCorrMatrix", "hAbsNpCorrMatrix;Kinetic Energy [MeV];Kinetic Energy [MeV]",
        NUM_BINS_KE, ARRAY_KE_BINS.data(),
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    GetCovMatrix(hPionAbsNpKENom, PionAbsNpKEUnivs, hAbsNpCovMatrix);
    GetFracCovAndCorrMatrix(hPionAbsNpKENom, hAbsNpCovMatrix, hAbsNpFracCovMatrix, hAbsNpCorrMatrix);
    DrawHistosWithErrorBands(
        hPionAbsNpKENom, 
        PionAbsNpKEUnivs, 
        SaveDir,
        "BeamlineEl", 
        "AbsNp",
        TextSize,
        FontStyle,
        "Abs Np Cross Section",
        "Kinetic Energy [MeV]",
        "Cross Section [barn]"
    );

    // Scatter
    TH2D* hScatterCovMatrix = new TH2D(
        "hScatterCovMatrix", "hScatterCovMatrix;Kinetic Energy [MeV];Kinetic Energy [MeV]",
        NUM_BINS_KE, ARRAY_KE_BINS.data(),
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    TH2D* hScatterFracCovMatrix = new TH2D(
        "hScatterFracCovMatrix", "hScatterFracCovMatrix;Kinetic Energy [MeV];Kinetic Energy [MeV]",
        NUM_BINS_KE, ARRAY_KE_BINS.data(),
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    TH2D* hScatterCorrMatrix = new TH2D(
        "hScatterCorrMatrix", "hScatterCorrMatrix;Kinetic Energy [MeV];Kinetic Energy [MeV]",
        NUM_BINS_KE, ARRAY_KE_BINS.data(),
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    GetCovMatrix(hPionScatterKENom, PionScatterKEUnivs, hScatterCovMatrix);
    GetFracCovAndCorrMatrix(hPionScatterKENom, hScatterCovMatrix, hScatterFracCovMatrix, hScatterCorrMatrix);
    DrawHistosWithErrorBands(
        hPionScatterKENom, 
        PionScatterKEUnivs, 
        SaveDir, 
        "BeamlineEl",
        "Scatter",
        TextSize,
        FontStyle,
        "Scatter Cross Section",
        "Kinetic Energy [MeV]",
        "Cross Section [barn]"
    );

    ///////////////////
    // Draw 2D Plots //
    ///////////////////

    std::vector<TH2*> TwoDPlots = {
        hMeasureCovMatrix,
        hMeasureFracCovMatrix, 
        hMeasureCorrMatrix,
        hAbs0pCovMatrix,
        hAbs0pFracCovMatrix,
        hAbs0pCorrMatrix,
        hAbsNpCovMatrix,
        hAbsNpFracCovMatrix,
        hAbsNpCorrMatrix,
        hScatterCovMatrix,
        hScatterFracCovMatrix,
        hScatterCorrMatrix
    };

    std::vector<TString> TwoDTitles = {
        "BeamlineEl/MeasureCovariance",
        "BeamlineEl/MeasureFracCovariance",
        "BeamlineEl/MeasureCorrelation",
        "BeamlineEl/Abs0pCovariance",
        "BeamlineEl/Abs0pFracCovariance",
        "BeamlineEl/Abs0pCorrelation",
        "BeamlineEl/AbsNpCovariance",
        "BeamlineEl/AbsNpFracCovariance",
        "BeamlineEl/AbsNpCorrelation",
        "BeamlineEl/ScatterCovariance",
        "BeamlineEl/ScatterFracCovariance",
        "BeamlineEl/ScatterCorrelation"
    };

    std::vector<std::pair<double,double>> TwoDRanges = {
        {0, 0},
        {0, 0},
        {-1, 1},
        {0, 0},
        {0, 0},
        {-1, 1},
        {0, 0},
        {0, 0},
        {-1, 1},
        {0, 0},
        {0, 0},
        {-1, 1}
    };

    std::vector<bool> TwoDDisplayNumbers = {
        false,
        false,
        false,
        false,
        false,
        false,
        false,
        false,
        false,
        false,
        false,
        false,
        false,
        false,
        false
    };

    printTwoDPlots(SaveDir, TwoDPlots, TwoDTitles, TwoDRanges, TwoDDisplayNumbers);

    ////////////////////////////
    // Save covariance matrix //
    ////////////////////////////

    TFile* saveFile = new TFile("/exp/lariat/app/users/epelaez/histos/beamlineel/Measure.root", "RECREATE");
    saveFile->cd();
    hMeasureCovMatrix->Write("", TObject::kOverwrite);
    saveFile->Close();
}