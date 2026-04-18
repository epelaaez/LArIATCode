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

void RecoClassify3Cat() {
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
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis_abs_scatt/figs/Classify3Cat/";

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

    // Primary track index
    static int primarytrkkey; Chain->SetBranchAddress("primarytrkkey", &primarytrkkey);

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

    /////////////////////////////////
    // Files for event information //
    /////////////////////////////////

    std::ofstream outFileEvents("files/Classify3Cat/GoodEvents.txt");

    TFile* nominalFile = new TFile("/exp/lariat/app/users/epelaez/analysis_abs_scatt/histos/nominal/Reco.root", "RECREATE");

    ///////////////////////
    // Create histograms //
    ///////////////////////

    TH1D* hTotalEvents = new TH1D("hTotalEvents", "hTotalEvents", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);

    // Histograms with classified events
    TH1D* hPionAbs     = new TH1D("hPionAbs", "hPionAbs;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D* hPionScatter = new TH1D("hPionScatter", "hPionScatter;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);

    // Events composition after each cut
    TH1D* hDataProdsAndWC2TPC = new TH1D("hDataProdsAndWC2TPC", "hDataProdsAndWC2TPC;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D* hNotManyTGTracks    = new TH1D("hNotManyTGTracks", "hNotManyTGTracks;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D* hNotAnElectron      = new TH1D("hNotAnElectron", "hNotAnElectron;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D* hPrimaryInRedVol    = new TH1D("hPrimaryInRedVol", "hPrimaryInRedVol;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D* hPrimaryPID         = new TH1D("hPrimaryPID", "hPrimaryPID;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D* hNotScatter         = new TH1D("hNotScatter", "hNotScatter;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D* hNotPionAbsNp       = new TH1D("hNotPionAbsNp", "hNotPionAbsNp;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D* hNotPionAbs0p       = new TH1D("hNotPionAbs0p", "hNotPionAbs0p;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);

    // Incident kinetic energy
    TH1D* hIncidentKE         = new TH1D("hIncidentKE", "hIncidentKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hIncidentKEPion     = new TH1D("hIncidentKEPion", "hIncidentKEPion;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hIncidentKEElectron = new TH1D("hIncidentKEElectron", "hIncidentKEElectron;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hIncidentKEMuon     = new TH1D("hIncidentKEMuon", "hIncidentKEMuon;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueIncidentKE     = new TH1D("hTrueIncidentKE", "hTrueIncidentKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Incident kinetic energy (finer binning)
    TH1D* hIncidentKEFine         = new TH1D("hIncidentKEFine", "hIncidentKEFine;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());
    TH1D* hIncidentKEPionFine     = new TH1D("hIncidentKEPionFine", "hIncidentKEPionFine;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());
    TH1D* hIncidentKEElectronFine = new TH1D("hIncidentKEElectronFine", "hIncidentKEElectronFine;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());
    TH1D* hIncidentKEMuonFine     = new TH1D("hIncidentKEMuonFine", "hIncidentKEMuonFine;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());

    // Front-face kinetic energy
    TH1D* hFrontFaceKE         = new TH1D("hFrontFaceKE", "hFrontFaceKE;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());
    TH1D* hFrontFaceKEPion     = new TH1D("hFrontFaceKEPion", "hFrontFaceKEPion;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());
    TH1D* hFrontFaceKEElectron = new TH1D("hFrontFaceKEElectron", "hFrontFaceKEElectron;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());
    TH1D* hFrontFaceKEMuon     = new TH1D("hFrontFaceKEMuon", "hFrontFaceKEMuon;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());

    TH1D* hFrontFacePionKE         = new TH1D("hFrontFacePionKE", "hFrontFacePionKE;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());
    TH1D* hFrontFacePionKEPion     = new TH1D("hFrontFacePionKEPion", "hFrontFacePionKEPion;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());
    TH1D* hFrontFacePionKEElectron = new TH1D("hFrontFacePionKEElectron", "hFrontFacePionKEElectron;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());
    TH1D* hFrontFacePionKEMuon     = new TH1D("hFrontFacePionKEMuon", "hFrontFacePionKEMuon;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());

    TH1D* hFrontFaceKENoWeight         = new TH1D("hFrontFaceKENoWeight", "hFrontFaceKENoWeight;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());
    TH1D* hFrontFaceKEPionNoWeight     = new TH1D("hFrontFaceKEPionNoWeight", "hFrontFaceKEPionNoWeight;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());
    TH1D* hFrontFaceKEElectronNoWeight = new TH1D("hFrontFaceKEElectronNoWeight", "hFrontFaceKEElectronNoWeight;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());
    TH1D* hFrontFaceKEMuonNoWeight     = new TH1D("hFrontFaceKEMuonNoWeight", "hFrontFaceKEMuonNoWeight;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());

    TH1D* hFrontFaceKENoWeightPre         = new TH1D("hFrontFaceKENoWeightPre", "hFrontFaceKENoWeightPre;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());
    TH1D* hFrontFaceKEPionNoWeightPre     = new TH1D("hFrontFaceKEPionNoWeightPre", "hFrontFaceKEPionNoWeightPre;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());
    TH1D* hFrontFaceKEElectronNoWeightPre = new TH1D("hFrontFaceKEElectronNoWeightPre", "hFrontFaceKEElectronNoWeightPre;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());
    TH1D* hFrontFaceKEMuonNoWeightPre     = new TH1D("hFrontFaceKEMuonNoWeightPre", "hFrontFaceKEMuonNoWeightPre;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());

    // Wire-chamber momentum
    TH1D* hWCKE         = new TH1D("hWCKE", "hWCKE;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());  
    TH1D* hWCKEPion     = new TH1D("hWCKEPion", "hWCKEPion;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());  
    TH1D* hWCKEMuon     = new TH1D("hWCKEMuon", "hWCKEMuon;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());  
    TH1D* hWCKEElectron = new TH1D("hWCKEElectron", "hWCKEElectron;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());  

    // True interacting flux
    TH1D* hTrueAllKE   = new TH1D("hTrueAllKE", "hTrueAllKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueOtherKE = new TH1D("hTrueOtherKE", "hTrueOtherKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Interacting pion absorption energy
    TH1D* hPionAbsKE         = new TH1D("hPionAbsKE", "hPionAbsKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbsKETrue     = new TH1D("hPionAbsKETrue", "hPionAbsKETrue;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbsKEScatter  = new TH1D("hPionAbsKEScatter", "hPionAbsKEScatter;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbsKEChExch   = new TH1D("hPionAbsKEChExch", "hPionAbsKEChExch;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbsKEMuon     = new TH1D("hPionAbsKEMuon", "hPionAbsKEMuon;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbsKEElectron = new TH1D("hPionAbsKEElectron", "hPionAbsKEElectron;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbsKEOther    = new TH1D("hPionAbsKEOther", "hPionAbsKEOther;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    TH1D* hPionAbsKEMuonTG    = new TH1D("hPionAbsKEMuonTG", "hPionAbsKEMuonTG;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbsKEMuonDecay = new TH1D("hPionAbsKEMuonDecay", "hPionAbsKEMuonDecay;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbsKEMuonCAR   = new TH1D("hPionAbsKEMuonCAR", "hPionAbsKEMuonCAR;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Estimated backgrounds for pion absorption
    std::vector<TH1*> PionAbsBkg = {
        hPionAbsKEChExch, hPionAbsKEMuon, hPionAbsKEElectron, hPionAbsKEOther
    };

    // Interacting pion scattering energy
    TH1D* hPionScatterKE         = new TH1D("hPionScatterKE", "hPionScatterKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionScatterKETrue     = new TH1D("hPionScatterKETrue", "hPionScatterKETrue;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionScatterKEAbs      = new TH1D("hPionScatterKEAbs", "hPionScatterKEAbs;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionScatterKEChExch   = new TH1D("hPionScatterKEChExch", "hPionScatterKEChExch;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionScatterKEMuon     = new TH1D("hPionScatterKEMuon", "hPionScatterKEMuon;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionScatterKEElectron = new TH1D("hPionScatterKEElectron", "hPionScatterKEElectron;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionScatterKEOther    = new TH1D("hPionScatterKEOther", "hPionScatterKEOther;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    TH1D* hPionScatterKEMuonTG    = new TH1D("hPionScatterKEMuonTG", "hPionScatterKEMuonTG;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionScatterKEMuonDecay = new TH1D("hPionScatterKEMuonDecay", "hPionScatterKEMuonDecay;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionScatterKEMuonCAR   = new TH1D("hPionScatterKEMuonCAR", "hPionScatterKEMuonCAR;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Estimated backgrounds for pion scattering
    std::vector<TH1*> PionScatterBkg = {
        hPionScatterKEChExch, hPionScatterKEMuon, hPionScatterKEElectron, hPionScatterKEOther
    };

    // All our reconstructed backgrounds
    std::vector<std::vector<TH1*>> RecoSignalBackgrounds = {
        PionAbsBkg, PionScatterBkg
    };

    // All our reconstructed signals
    std::vector<TH1*> RecoSignals = {
        hPionAbsKE, hPionScatterKE
    };

    // True absorption (abs0p + absNp merged)
    TH1D* hTrueAbsKE          = new TH1D("hTrueAbsKE", "hTrueAbsKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbsKEAsAbs     = new TH1D("hTrueAbsKEAsAbs", "hTrueAbsKEAsAbs;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbsKEAsScatter = new TH1D("hTrueAbsKEAsScatter", "hTrueAbsKEAsScatter;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbsKERejected  = new TH1D("hTrueAbsKERejected", "hTrueAbsKERejected;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // For each TRUE energy bin, what are absorption events reconstructed as?
    std::vector<std::vector<TH1D*>> TrueAbsAsByBin;
    for (int iEnergyBin = 0; iEnergyBin < NUM_BINS_KE; ++iEnergyBin) {
        std::vector<TH1D*> TempVec;
        for (int iInt = 0; iInt < NUM_SIGNAL_TYPES; ++iInt) {
            TH1D* hTempHist = new TH1D(Form("hTrueAbs_%d_Bin_As_%d", iEnergyBin, iInt), Form("True Abs KE As %d in bin %d", iInt, iEnergyBin), NUM_BINS_KE, ARRAY_KE_BINS.data());
            TempVec.push_back(hTempHist);
        }
        TrueAbsAsByBin.push_back(TempVec);
    }

    TH1D* hTrueAbsKERejDataProds = new TH1D("hTrueAbsKERejDataProds", "hTrueAbsKERejDataProds;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbsKERejElectron  = new TH1D("hTrueAbsKERejElectron", "hTrueAbsKERejElectron;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbsKERejRedVol    = new TH1D("hTrueAbsKERejRedVol", "hTrueAbsKERejRedVol;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbsKERejPID       = new TH1D("hTrueAbsKERejPID", "hTrueAbsKERejPID;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbsKERejManyPions = new TH1D("hTrueAbsKERejManyPions", "hTrueAbsKERejManyPions;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbsKERejClusters  = new TH1D("hTrueAbsKERejClusters", "hTrueAbsKERejClusters;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // True scatter
    TH1D* hTrueScatterKE          = new TH1D("hTrueScatterKE", "hTrueScatterKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueScatterKEAsAbs     = new TH1D("hTrueScatterKEAsAbs", "hTrueScatterKEAsAbs;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueScatterKEAsScatter = new TH1D("hTrueScatterKEAsScatter", "hTrueScatterKEAsScatter;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueScatterKERejected  = new TH1D("hTrueScatterKERejected", "hTrueScatterKERejected;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // For each TRUE energy bin, what are scatter events reconstructed as?
    std::vector<std::vector<TH1D*>> TrueScatterAsByBin;
    for (int iEnergyBin = 0; iEnergyBin < NUM_BINS_KE; ++iEnergyBin) {
        std::vector<TH1D*> TempVec;
        for (int iInt = 0; iInt < NUM_SIGNAL_TYPES; ++iInt) {
            TH1D* hTempHist = new TH1D(Form("hTrueScatter_%d_Bin_As_%d", iEnergyBin, iInt), Form("True Scatter KE As %d in bin %d", iInt, iEnergyBin), NUM_BINS_KE, ARRAY_KE_BINS.data());
            TempVec.push_back(hTempHist);
        }
        TrueScatterAsByBin.push_back(TempVec);
    }

    TH1D* hTrueScatterKERejDataProds = new TH1D("hTrueScatterKERejDataProds", "hTrueScatterKERejDataProds;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueScatterKERejElectron  = new TH1D("hTrueScatterKERejElectron", "hTrueScatterKERejElectron;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueScatterKERejRedVol    = new TH1D("hTrueScatterKERejRedVol", "hTrueScatterKERejRedVol;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueScatterKERejPID       = new TH1D("hTrueScatterKERejPID", "hTrueScatterKERejPID;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueScatterKERejManyPions = new TH1D("hTrueScatterKERejManyPions", "hTrueScatterKERejManyPions;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueScatterKERejClusters  = new TH1D("hTrueScatterKERejClusters", "hTrueScatterKERejClusters;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // True charge exchange
    TH1D* hTrueChExchKE          = new TH1D("hTrueChExchKE", "hTrueChExchKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueChExchKEAsAbs     = new TH1D("hTrueChExchKEAsAbs", "hTrueChExchKEAsAbs;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueChExchKEAsScatter = new TH1D("hTrueChExchKEAsScatter", "hTrueChExchKEAsScatter;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueChExchKERejected  = new TH1D("hTrueChExchKERejected", "hTrueChExchKERejected;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    //////////////////////////
    // Through-going tracks //
    //////////////////////////

    TH1D* hNumTGTracksAbs      = new TH1D("hNumTGTracksAbs", "hNumTGTracksAbs", 10, 0, 10);
    TH1D* hNumTGTracksMuon     = new TH1D("hNumTGTracksMuon", "hNumTGTracksMuon", 10, 0, 10);
    TH1D* hNumTGTracksElectron = new TH1D("hNumTGTracksElectron", "hNumTGTracksElectron", 10, 0, 10);
    TH1D* hNumTGTracksScatter  = new TH1D("hNumTGTracksScatter", "hNumTGTracksScatter", 10, 0, 10);
    TH1D* hNumTGTracksChExch   = new TH1D("hNumTGTracksChExch", "hNumTGTracksChExch", 10, 0, 10);
    TH1D* hNumTGTracksOther    = new TH1D("hNumTGTracksOther", "hNumTGTracksOther", 10, 0, 10);

    ////////////////////////////////
    // Study unreconstructed hits //
    ////////////////////////////////

    TH1D* hNegativeTimePrimaryHits = new TH1D("hNegativeTimePrimaryHits", "hNegativeTimePrimaryHits", 30, -20, 10);
    TH1D* hTimePrimaryHits         = new TH1D("hTimePrimaryHits", "hTimePrimaryHits", 60, -20, 300);

    TH1D* hUnRecoHitsInductionAbs      = new TH1D("hUnRecoHitsInductionAbs", "hUnRecoHitsInductionAbs", 25, 0, 50);
    TH1D* hUnRecoHitsInductionMuon     = new TH1D("hUnRecoHitsInductionMuon", "hUnRecoHitsInductionMuon", 25, 0, 50);
    TH1D* hUnRecoHitsInductionElectron = new TH1D("hUnRecoHitsInductionElectron", "hUnRecoHitsInductionElectron", 25, 0, 50);
    TH1D* hUnRecoHitsInductionScatter  = new TH1D("hUnRecoHitsInductionScatter", "hUnRecoHitsInductionScatter", 25, 0, 50);
    TH1D* hUnRecoHitsInductionChExch   = new TH1D("hUnRecoHitsInductionChExch", "hUnRecoHitsInductionChExch", 25, 0, 50);
    TH1D* hUnRecoHitsInductionOther    = new TH1D("hUnRecoHitsInductionOther", "hUnRecoHitsInductionOther", 25, 0, 50);

    TH1D* hUnRecoHitsCollectionAbs      = new TH1D("hUnRecoHitsCollectionAbs", "hUnRecoHitsCollectionAbs", 25, 0, 50);
    TH1D* hUnRecoHitsCollectionMuon     = new TH1D("hUnRecoHitsCollectionMuon", "hUnRecoHitsCollectionMuon", 25, 0, 50);
    TH1D* hUnRecoHitsCollectionElectron = new TH1D("hUnRecoHitsCollectionElectron", "hUnRecoHitsCollectionElectron", 25, 0, 50);
    TH1D* hUnRecoHitsCollectionScatter  = new TH1D("hUnRecoHitsCollectionScatter", "hUnRecoHitsCollectionScatter", 25, 0, 50);
    TH1D* hUnRecoHitsCollectionChExch   = new TH1D("hUnRecoHitsCollectionChExch", "hUnRecoHitsCollectionChExch", 25, 0, 50);
    TH1D* hUnRecoHitsCollectionOther    = new TH1D("hUnRecoHitsCollectionOther", "hUnRecoHitsCollectionOther", 25, 0, 50);

    TH1D* hHitClusterInductionSizesAbs      = new TH1D("hHitClusterInductionSizesAbs", "hHitClusterInductionSizesAbs;;", 20, 0, 20);
    TH1D* hHitClusterInductionSizesMuon     = new TH1D("hHitClusterInductionSizesMuon", "hHitClusterInductionSizesMuon;;", 20, 0, 20);
    TH1D* hHitClusterInductionSizesElectron = new TH1D("hHitClusterInductionSizesElectron", "hHitClusterInductionSizesElectron;;", 20, 0, 20);
    TH1D* hHitClusterInductionSizesScatter  = new TH1D("hHitClusterInductionSizesScatter", "hHitClusterInductionSizesScatter;;", 20, 0, 20);
    TH1D* hHitClusterInductionSizesChExch   = new TH1D("hHitClusterInductionSizesChExch", "hHitClusterInductionSizesChExch;;", 20, 0, 20);
    TH1D* hHitClusterInductionSizesOther    = new TH1D("hHitClusterInductionSizesOther", "hHitClusterInductionSizesOther;;", 20, 0, 20);

    TH1D* hHitClusterCollectionSizesAbs      = new TH1D("hHitClusterCollectionSizesAbs", "hHitClusterCollectionSizesAbs;;", 20, 0, 20);
    TH1D* hHitClusterCollectionSizesMuon     = new TH1D("hHitClusterCollectionSizesMuon", "hHitClusterCollectionSizesMuon;;", 20, 0, 20);
    TH1D* hHitClusterCollectionSizesElectron = new TH1D("hHitClusterCollectionSizesElectron", "hHitClusterCollectionSizesElectron;;", 20, 0, 20);
    TH1D* hHitClusterCollectionSizesScatter  = new TH1D("hHitClusterCollectionSizesScatter", "hHitClusterCollectionSizesScatter;;", 20, 0, 20);
    TH1D* hHitClusterCollectionSizesChExch   = new TH1D("hHitClusterCollectionSizesChExch", "hHitClusterCollectionSizesChExch;;", 20, 0, 20);
    TH1D* hHitClusterCollectionSizesOther    = new TH1D("hHitClusterCollectionSizesOther", "hHitClusterCollectionSizesOther;;", 20, 0, 20);

    TH1D* hLargestHitClusterInductionAbs      = new TH1D("hLargestHitClusterInductionAbs", "hLargestHitClusterInductionAbs;;", 40, 0, 20);
    TH1D* hLargestHitClusterInductionMuon     = new TH1D("hLargestHitClusterInductionMuon", "hLargestHitClusterInductionMuon;;", 40, 0, 20);
    TH1D* hLargestHitClusterInductionElectron = new TH1D("hLargestHitClusterInductionElectron", "hLargestHitClusterInductionElectron;;", 40, 0, 20);
    TH1D* hLargestHitClusterInductionScatter  = new TH1D("hLargestHitClusterInductionScatter", "hLargestHitClusterInductionScatter;;", 40, 0, 20);
    TH1D* hLargestHitClusterInductionChExch   = new TH1D("hLargestHitClusterInductionChExch", "hLargestHitClusterInductionChExch;;", 40, 0, 20);
    TH1D* hLargestHitClusterInductionOther    = new TH1D("hLargestHitClusterInductionOther", "hLargestHitClusterInductionOther;;", 40, 0, 20);

    TH1D* hLargestHitClusterCollectionAbs      = new TH1D("hLargestHitClusterCollectionAbs", "hLargestHitClusterCollectionAbs;;", 40, 0, 20);
    TH1D* hLargestHitClusterCollectionMuon     = new TH1D("hLargestHitClusterCollectionMuon", "hLargestHitClusterCollectionMuon;;", 40, 0, 20);
    TH1D* hLargestHitClusterCollectionElectron = new TH1D("hLargestHitClusterCollectionElectron", "hLargestHitClusterCollectionElectron;;", 40, 0, 20);
    TH1D* hLargestHitClusterCollectionScatter  = new TH1D("hLargestHitClusterCollectionScatter", "hLargestHitClusterCollectionScatter;;", 40, 0, 20);
    TH1D* hLargestHitClusterCollectionChExch   = new TH1D("hLargestHitClusterCollectionChExch", "hLargestHitClusterCollectionChExch;;", 40, 0, 20);
    TH1D* hLargestHitClusterCollectionOther    = new TH1D("hLargestHitClusterCollectionOther", "hLargestHitClusterCollectionOther;;", 40, 0, 20);

    TH1D* hNumClustersInductionAbs      = new TH1D("hNumClustersInductionAbs", "hNumClustersInductionAbs;;", 10, 0, 10);
    TH1D* hNumClustersInductionMuon     = new TH1D("hNumClustersInductionMuon", "hNumClustersInductionMuon;;", 10, 0, 10);
    TH1D* hNumClustersInductionElectron = new TH1D("hNumClustersInductionElectron", "hNumClustersInductionElectron;;", 10, 0, 10);
    TH1D* hNumClustersInductionScatter  = new TH1D("hNumClustersInductionScatter", "hNumClustersInductionScatter;;", 10, 0, 10);
    TH1D* hNumClustersInductionChExch   = new TH1D("hNumClustersInductionChExch", "hNumClustersInductionChExch;;", 10, 0, 10);
    TH1D* hNumClustersInductionOther    = new TH1D("hNumClustersInductionOther", "hNumClustersInductionOther;;", 10, 0, 10);

    TH1D* hNumClustersCollectionAbs      = new TH1D("hNumClustersCollectionAbs", "hNumClustersCollectionAbs;;", 10, 0, 10);
    TH1D* hNumClustersCollectionMuon     = new TH1D("hNumClustersCollectionMuon", "hNumClustersCollectionMuon;;", 10, 0, 10);
    TH1D* hNumClustersCollectionElectron = new TH1D("hNumClustersCollectionElectron", "hNumClustersCollectionElectron;;", 10, 0, 10);
    TH1D* hNumClustersCollectionScatter  = new TH1D("hNumClustersCollectionScatter", "hNumClustersCollectionScatter;;", 10, 0, 10);
    TH1D* hNumClustersCollectionChExch   = new TH1D("hNumClustersCollectionChExch", "hNumClustersCollectionChExch;;", 10, 0, 10);
    TH1D* hNumClustersCollectionOther    = new TH1D("hNumClustersCollectionOther", "hNumClustersCollectionOther;;", 10, 0, 10);

    TH1D* hLargeHitClusterInductionAbs      = new TH1D("hLargeHitClusterInductionAbs", "hLargeHitClusterInductionAbs;;", 5, 0, 5);
    TH1D* hLargeHitClusterInductionMuon     = new TH1D("hLargeHitClusterInductionMuon", "hLargeHitClusterInductionMuon;;", 5, 0, 5);
    TH1D* hLargeHitClusterInductionElectron = new TH1D("hLargeHitClusterInductionElectron", "hLargeHitClusterInductionElectron;;", 5, 0, 5);
    TH1D* hLargeHitClusterInductionScatter  = new TH1D("hLargeHitClusterInductionScatter", "hLargeHitClusterInductionScatter;;", 5, 0, 5);
    TH1D* hLargeHitClusterInductionChExch   = new TH1D("hLargeHitClusterInductionChExch", "hLargeHitClusterInductionChExch;;", 5, 0, 5);
    TH1D* hLargeHitClusterInductionOther    = new TH1D("hLargeHitClusterInductionOther", "hLargeHitClusterInductionOther;;", 5, 0, 5);

    TH1D* hLargeHitClusterCollectionAbs      = new TH1D("hLargeHitClusterCollectionAbs", "hLargeHitClusterCollectionAbs;;", 5, 0, 5);
    TH1D* hLargeHitClusterCollectionMuon     = new TH1D("hLargeHitClusterCollectionMuon", "hLargeHitClusterCollectionMuon;;", 5, 0, 5);
    TH1D* hLargeHitClusterCollectionElectron = new TH1D("hLargeHitClusterCollectionElectron", "hLargeHitClusterCollectionElectron;;", 5, 0, 5);
    TH1D* hLargeHitClusterCollectionScatter  = new TH1D("hLargeHitClusterCollectionScatter", "hLargeHitClusterCollectionScatter;;", 5, 0, 5);
    TH1D* hLargeHitClusterCollectionChExch   = new TH1D("hLargeHitClusterCollectionChExch", "hLargeHitClusterCollectionChExch;;", 5, 0, 5);
    TH1D* hLargeHitClusterCollectionOther    = new TH1D("hLargeHitClusterCollectionOther", "hLargeHitClusterCollectionOther;;", 5, 0, 5);

    //////////////
    // Matrices //
    //////////////

    std::vector<TH1D*> TotalEventsHistos = {
        hTrueAbsKE, hTrueScatterKE
    };

    std::vector<std::vector<std::vector<TH1D*>>> TrueRecoAsByBin = {
        TrueAbsAsByBin, TrueScatterAsByBin
    };

    ////////////////
    // Muon types //
    ////////////////

    TH1D* hTrueMuonTypes = new TH1D("hTrueMuonTypes", "hTrueMuonTypes", 4, 0, 4);

    //////////////////////
    // Loop over events //
    //////////////////////

    int scatteringsModified = 0;

    int numEventsPrimaryReco  = 0;
    int numEventsNoCollection = 0;
    int numEventsNoInduction  = 0;
    int numEventsNoEither     = 0;

    int numEventsPrimaryHitNegativeTime = 0;

    // Keep track of event counts for stat estimations
    int eventCount0TG = 0;
    int eventCount1TG = 0;
    int eventCount2TG = 0;

    bool verbose = false;

    Long64_t i = 0;
    while (true) {
        if (SKIP_INDICES.count(i)) { i++; continue; }
        if (Chain->GetEntry(i++) <= 0) break;

        // Make script go faster
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
            if (trkWCtoTPCMatch[trk_idx] == 1) {
                // Grab position information
                ev.WC2TPCPrimaryBeginX = trkvtxx[trk_idx];
                ev.WC2TPCPrimaryBeginY = trkvtxy[trk_idx];
                ev.WC2TPCPrimaryBeginZ = trkvtxz[trk_idx];
                ev.WC2TPCPrimaryEndX   = trkendx[trk_idx];
                ev.WC2TPCPrimaryEndY   = trkendy[trk_idx];
                ev.WC2TPCPrimaryEndZ   = trkendz[trk_idx];

                if (verbose) std::cout << "Begin: " << ev.WC2TPCPrimaryBeginX << " " << ev.WC2TPCPrimaryBeginY << " " << ev.WC2TPCPrimaryBeginZ << std::endl;
                if (verbose) std::cout << "End:   " << ev.WC2TPCPrimaryEndX << " " << ev.WC2TPCPrimaryEndY << " " << ev.WC2TPCPrimaryEndZ << std::endl;

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

                if (verbose) std::cout << nTrajPoint[trk_idx] << "  " << kMaxTrajHits << std::endl;
                if (verbose) std::cout << "(loc.) Begin: " << ev.WC2TPCLocationsX[0] << " " << ev.WC2TPCLocationsY[0] << " " << ev.WC2TPCLocationsZ[0] << std::endl;
                if (verbose) std::cout << "(loc.) End:   " << ev.WC2TPCLocationsX[npts_wc2tpc - 1] << " " << ev.WC2TPCLocationsY[npts_wc2tpc - 1] << " " << ev.WC2TPCLocationsZ[npts_wc2tpc - 1] << std::endl;

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
        float lastTPCX = -999, lastTPCY = -999, lastTPCZ = -999;
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
                }
                if (verbose) std::cout << "Finished trajectory" << std::endl;

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
                        // << "  within red vol? " << isWithinReducedVolume(currentPos.X(), currentPos.Y(), currentPos.Z())
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

        if (verbose) std::cout << "Filling out hit information" << std::endl;
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
            if (trkid < 0 || trkid >= ntracks_reco) {
                std::cerr << "WARNING: unexpected trkid=" << trkid << " for entry " << i << ", skipping\n";
                continue;
            }
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

        // Get unordered set for hits in tracks
        std::unordered_set<int> hitsInTracks(ev.hitRecoAsTrackKey.begin(), ev.hitRecoAsTrackKey.end());

        // Sanity check
        removeRepeatedPoints(&ev.WC2TPCLocationsX, &ev.WC2TPCLocationsY, &ev.WC2TPCLocationsZ);

        // Extend reco cylinder
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

        // Units for coordinates:
        //     W: (channel) * (channel width) = cm
        //     X: (drift velocity) * (time) = (cm / us) * us = cm

        // First, find hits near the primary track that were not
        // reconstructed into any of the already existing tracks
        std::vector<int> candidateHits;

        int numTotalHitsNearPrimary            = 0;
        int numUnRecoHitsNearPrimaryInduction  = 0;
        int numUnRecoHitsNearPrimaryCollection = 0;
        for (size_t iHit = 0; iHit < std::min(nhits, kMaxHits); ++iHit) {
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

        // Separate primary hits into collection and induction
        std::vector<int> hitWC2TPCKeyInduction;
        std::vector<int> hitWC2TPCKeyCollection;

        bool sawPrimaryInduction  = false; 
        bool sawPrimaryCollection = false;

        bool foundHitNegativeTime = false;
        for (size_t i = 0; i < ev.hitWC2TPCKey.size(); ++i) {
            // Only care about hits inside the reduced volume
            auto [i_hit, j_hit] = find_unique_position(&ev.recoTrackHitIndices, ev.hitWC2TPCKey.at(i));

            if (ev.fHitPlane.at(ev.hitWC2TPCKey.at(i)) == 0) sawPrimaryInduction = true;
            else if (ev.fHitPlane.at(ev.hitWC2TPCKey.at(i)) == 1) sawPrimaryCollection = true;

            if (!isWithinReducedVolume(
                ev.recoTrackHitX.at(i_hit)[j_hit],
                ev.recoTrackHitY.at(i_hit)[j_hit],
                ev.recoTrackHitZ.at(i_hit)[j_hit]
            )) continue;

            if (ev.fHitPlane.at(ev.hitWC2TPCKey.at(i)) == 0) hitWC2TPCKeyInduction.push_back(ev.hitWC2TPCKey.at(i));
            else if (ev.fHitPlane.at(ev.hitWC2TPCKey.at(i)) == 1) hitWC2TPCKeyCollection.push_back(ev.hitWC2TPCKey.at(i));

            // Look at hits with negative time
            if (ev.fHitT.at(ev.hitWC2TPCKey.at(i)) < 0) {
                foundHitNegativeTime = true;
                hNegativeTimePrimaryHits->Fill(ev.fHitT.at(ev.hitWC2TPCKey.at(i)), ev.weight);
            }
            hTimePrimaryHits->Fill(ev.fHitT.at(ev.hitWC2TPCKey.at(i)), ev.weight);
        }
        if (foundHitNegativeTime) numEventsPrimaryHitNegativeTime++;

        if (!sawPrimaryCollection) numEventsNoCollection++;
        if (!sawPrimaryInduction) numEventsNoInduction++;
        if (!sawPrimaryCollection && !sawPrimaryInduction) numEventsNoEither++;

        ///////////////////////
        // Get front face KE //
        ///////////////////////

        double WCKE             = TMath::Sqrt(ev.WCTrackMomentum * ev.WCTrackMomentum + PionMass * PionMass) - PionMass;
        double calculatedEnLoss = energyLossCalculation(); 
        if (isData) {
            double tanThetaCosPhi = TMath::Tan(ev.WCTheta) * TMath::Cos(ev.WCPhi);
            double tanThetaSinPhi = TMath::Tan(ev.WCTheta) * TMath::Sin(ev.WCPhi);
            double den            = TMath::Sqrt(1 + tanThetaCosPhi * tanThetaCosPhi);
            double onTheFlyPz     = ev.WCTrackMomentum / den;
            double onTheFlyPx     = onTheFlyPz * tanThetaSinPhi;
            calculatedEnLoss      = energyLossCalculation(ev.WC4PrimaryX, onTheFlyPx, isData);
        } else { calculatedEnLoss = energyLossCalculation(ev.WC4PrimaryX, ev.trajectoryInitialMomentumX, isData); }
        const double initialKE = WCKE - calculatedEnLoss;

        hFrontFaceKENoWeightPre->Fill(initialKE);
        if (ev.truthPrimaryPDG == -211) {
            hFrontFaceKEPionNoWeightPre->Fill(initialKE);
        } else if (ev.truthPrimaryPDG == 13) {
            hFrontFaceKEMuonNoWeightPre->Fill(initialKE);
        } else if (ev.truthPrimaryPDG == 11) {
            hFrontFaceKEElectronNoWeightPre->Fill(initialKE);
        }

        // Reweigh!
        ev.weight *= GetKEWeight(hWeightsFrontFace, initialKE);

        ////////////////////////////
        // Save truth information //
        ////////////////////////////

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
                        if (
                            ev.truthPrimaryDaughtersPDG.at(i) == 2212 &&
                            ev.truthPrimaryDaughtersProcess.at(i) == "pi-Inelastic" &&
                            ev.truthPrimaryDaughtersKE.at(i) >= PROTON_ENERGY_LOWER_BOUND &&
                            ev.truthPrimaryDaughtersKE.at(i) <= PROTON_ENERGY_UPPER_BOUND
                        ) {
                            tempNumVisibleProtons++;
                        }
                    }
                    if (tempNumVisibleProtons == 0) ev.backgroundType = 0;
                    else ev.backgroundType = 1;
                    ev.numVisibleProtons = tempNumVisibleProtons;
                }
            }
            // If pion above threshold but angle not large enough, go to next interaction and check there
            else if (scatteringAngle < SCATTERING_ANGLE_THRESHOLD) {
                scatteringsModified++;
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
                                if (
                                    ev.secondaryInteractionDaughtersPDG.at(iInteraction)[i] == 2212 &&
                                    ev.secondaryInteractionDaughtersKE.at(iInteraction)[i] >= PROTON_ENERGY_LOWER_BOUND &&
                                    ev.secondaryInteractionDaughtersKE.at(iInteraction)[i] <= PROTON_ENERGY_UPPER_BOUND
                                ) {
                                    tempNumVisibleProtons++;
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

        // Further classify muon backgrounds
        int muonType = -1; // 0: through-going, 1: decaying, 2: capture at rest
        if (ev.backgroundType == 2) {
            if (!isWithinReducedVolume(ev.truthPrimaryVertexX, ev.truthPrimaryVertexY, ev.truthPrimaryVertexZ)) {
                // through-going muon
                muonType = 0;
            } else {
                bool sawElectron = false;
                for (int i = 0; i < ev.truthPrimaryDaughtersPDG.size(); ++i) {
                    if (ev.truthPrimaryDaughtersProcess.at(i) == "muMinusCaptureAtRest") {
                        // muon capture at rest
                        muonType = 2;
                        break;
                    }
                    if (ev.truthPrimaryDaughtersPDG.at(i) == 11) {
                        sawElectron = true;
                    }
                }

                if (muonType == -1 && sawElectron) {
                    // muon decay
                    muonType = 1;
                }
            }
            
            if (muonType != 0 && muonType != 1 && muonType != 2) {
                // I don't know what this event is
                muonType = 3;
            }

            hTrueMuonTypes->Fill(muonType, ev.weight);
        }

        ///////////////////////////
        // Fill true interacting //
        ///////////////////////////

        if (verbose) std::cout << "Original interacting KE: " << ev.truthPrimaryVertexKE << std::endl;
        if (ev.trueIncidentKEContributions.size() == 0) ev.truthPrimaryVertexKE = 0;
        else ev.truthPrimaryVertexKE = ev.trueIncidentKEContributions.back() / 1000.0;
        if (verbose) for (const auto& ke : ev.trueIncidentKEContributions) {std::cout << ke << " ";}
        if (verbose) if (ev.trueIncidentKEContributions.size() > 0) std::cout << std::endl;
        if (verbose) std::cout << "Final interacting KE: " << ev.truthPrimaryVertexKE << std::endl;

        // Get true energy bin
        int TrueEnergyBin = getBin(ev.truthPrimaryVertexKE * 1000, ARRAY_KE_BINS);

        // Add true incident KE
        if (ev.validTrueIncidentKE) {
            for (double x : ev.trueIncidentKEContributions) hTrueIncidentKE->Fill(x, ev.weight);
        }

        // Fill in true interaction energy histograms
        if (ev.backgroundType == 0 || ev.backgroundType == 1) {
            hTrueAbsKE->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
            hTrueAllKE->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
        } else if (ev.backgroundType == 6 || ev.backgroundType == 12) {
            hTrueScatterKE->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
            hTrueAllKE->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
        } else if (ev.backgroundType == 7) {
            hTrueChExchKE->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
            hTrueAllKE->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
        } else if (
            ev.backgroundType == 8 ||
            ev.backgroundType == 9 ||
            ev.backgroundType == 10 ||
            ev.backgroundType == 11
        ) {
            hTrueOtherKE->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
            hTrueAllKE->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
        }
        hTotalEvents->Fill(ev.backgroundType, ev.weight);

        //////////////////////
        // WC2TPC match cut //
        //////////////////////

        // If no track matched to wire-chamber, skip
        if (!ev.WC2TPCMatch || ev.WC2TPCsize != 1) {
            if (ev.backgroundType == 0 || ev.backgroundType == 1) {
                hTrueAbsKERejected->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
                hTrueAbsKERejDataProds->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
            } else if (ev.backgroundType == 6 || ev.backgroundType == 12) {
                hTrueScatterKERejected->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
                hTrueScatterKERejDataProds->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
            } else if (ev.backgroundType == 7) {
                hTrueChExchKERejected->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
            }
            continue;
        }
        hDataProdsAndWC2TPC->Fill(ev.backgroundType, ev.weight);

        //////////////////////////////////
        // Check WC/front-face momentum //
        //////////////////////////////////

        hWCKE->Fill(WCKE, ev.weight);
        if (ev.truthPrimaryPDG == -211) {
            hWCKEPion->Fill(WCKE, ev.weight);
        } else if (ev.truthPrimaryPDG == 13) {
            hWCKEMuon->Fill(WCKE, ev.weight);
        } else if (ev.truthPrimaryPDG == 11) {
            hWCKEElectron->Fill(WCKE, ev.weight);
        }

        hFrontFaceKE->Fill(initialKE, ev.weight);
        hFrontFaceKENoWeight->Fill(initialKE);
        if (ev.truthPrimaryPDG == -211) {
            hFrontFaceKEPion->Fill(initialKE, ev.weight);
            hFrontFaceKEPionNoWeight->Fill(initialKE);
        } else if (ev.truthPrimaryPDG == 13) {
            hFrontFaceKEMuon->Fill(initialKE, ev.weight);
            hFrontFaceKEMuonNoWeight->Fill(initialKE);
        } else if (ev.truthPrimaryPDG == 11) {
            hFrontFaceKEElectron->Fill(initialKE, ev.weight);
            hFrontFaceKEElectronNoWeight->Fill(initialKE);
        }

        //////////////////////////////////////
        // Back to classification algorithm //
        //////////////////////////////////////

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
        if (verbose) std::cout << "Chi2 results:" << std::endl;
        if (verbose) std::cout << "    Pion: " << pionChi2 << std::endl;
        if (verbose) std::cout << "    MIP: " << MIPChi2 << std::endl;
        if (verbose) std::cout << "    Proton: " << protonChi2 << std::endl;
        if (verbose) std::cout << "    Stitch: " << minStitchedChi2 << std::endl;

        // If primary track stitched, get break point, otherwise break point is end of track
        double breakPointX = ev.WC2TPCPrimaryEndX; 
        double breakPointY = ev.WC2TPCPrimaryEndY; 
        double breakPointZ = ev.WC2TPCPrimaryEndZ;
        if (minChi2 == minStitchedChi2) {
            breakPointX = ev.wcMatchXPos.at(bestBreakPoint);
            breakPointY = ev.wcMatchYPos.at(bestBreakPoint);
            breakPointZ = ev.wcMatchZPos.at(bestBreakPoint);
        }

        ////////////////////////////////
        // Cylinder and TG track cuts //
        ////////////////////////////////

        int numTGTracks              = 0;
        int numSmallTracksInCylinder = 0;
        int numTracksInCylinder      = 0;

        for (int trk_idx = 0; trk_idx < ev.recoBeginX.size(); ++trk_idx) {
            if (trkWCtoTPCMatch[trk_idx]) continue;

            // Check if track is through-going
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
        
        // Save data about TG tracks
        if (ev.backgroundType == 0 || ev.backgroundType == 1) {
            hNumTGTracksAbs->Fill(numTGTracks, ev.weight);
        } else if (ev.backgroundType == 2) {
            hNumTGTracksMuon->Fill(numTGTracks, ev.weight);
        } else if (ev.backgroundType == 3) {
            hNumTGTracksElectron->Fill(numTGTracks, ev.weight);
        } else if (ev.backgroundType == 6 || ev.backgroundType == 12) {
            hNumTGTracksScatter->Fill(numTGTracks, ev.weight);
        } else if (ev.backgroundType == 7) {
            hNumTGTracksChExch->Fill(numTGTracks, ev.weight);
        } else {
            hNumTGTracksOther->Fill(numTGTracks, ev.weight);
        }

        // Grab data about number of events with each cutoff
        if (numTGTracks <= 0) eventCount0TG++;
        if (numTGTracks <= 1) eventCount1TG++;
        if (numTGTracks <= 2) eventCount2TG++;

        // Perform TG track cut
        if (numTGTracks > MAX_NUM_TG_TRACKS) continue;
        hNotManyTGTracks->Fill(ev.backgroundType, ev.weight);

        if (numSmallTracksInCylinder > ALLOWED_CYLINDER_SMALL_TRACKS) continue;
        hNotAnElectron->Fill(ev.backgroundType, ev.weight);

        hFrontFacePionKE->Fill(initialKE, ev.weight);
        if (ev.truthPrimaryPDG == -211) {
            hFrontFacePionKEPion->Fill(initialKE, ev.weight);
        } else if (ev.truthPrimaryPDG == 13) {
            hFrontFacePionKEMuon->Fill(initialKE, ev.weight);
        } else if (ev.truthPrimaryPDG == 11) {
            hFrontFacePionKEElectron->Fill(initialKE, ev.weight);
        }

        //////////////////////
        // Incident KE fill //
        //////////////////////
        
        // Check these are in order
        if (ev.wcMatchZPos.size() > 1 && ev.wcMatchZPos.front() > ev.wcMatchZPos.back()) {
            std::reverse(ev.wcMatchZPos.begin(), ev.wcMatchZPos.end());
            std::reverse(ev.wcMatchDEDX.begin(), ev.wcMatchDEDX.end());
            std::reverse(ev.wcMatchEDep.begin(), ev.wcMatchEDep.end());
            std::reverse(ev.wcMatchXPos.begin(), ev.wcMatchXPos.end());
            std::reverse(ev.wcMatchYPos.begin(), ev.wcMatchYPos.end());
        }

        double energyDeposited = 0.0;
        for (size_t iDep = 0; iDep < ev.wcMatchDEDX.size(); ++iDep) {
            // If we are past detected breaking point, exit loop
            if (ev.wcMatchZPos.at(iDep) > breakPointZ) break;

            // If larger than threshold, continue
            if (ev.wcMatchDEDX.at(iDep) > HIT_DEDX_THRESHOLD) continue;

            // Else, add to energy deposited so far
            energyDeposited += ev.wcMatchEDep.at(iDep);

            // Add to incident KE if inside reduced volume
            if (isWithinReducedVolume(ev.wcMatchXPos.at(iDep), ev.wcMatchYPos.at(iDep), ev.wcMatchZPos.at(iDep))) {
                hIncidentKE->Fill(initialKE - energyDeposited, ev.weight);
                hIncidentKEFine->Fill(initialKE - energyDeposited, ev.weight);

                // Background breakdown
                if (ev.truthPrimaryPDG == -211) {
                    hIncidentKEPion->Fill(initialKE - energyDeposited, ev.weight);
                    hIncidentKEPionFine->Fill(initialKE - energyDeposited, ev.weight);
                } else if (ev.truthPrimaryPDG == 13) {
                    hIncidentKEMuon->Fill(initialKE - energyDeposited, ev.weight);
                    hIncidentKEMuonFine->Fill(initialKE - energyDeposited, ev.weight);
                } else if (ev.truthPrimaryPDG == 11) {
                    hIncidentKEElectron->Fill(initialKE - energyDeposited, ev.weight);
                    hIncidentKEElectronFine->Fill(initialKE - energyDeposited, ev.weight);
                }
            }
        }
        double energyAtVertex = initialKE - energyDeposited;

        ////////////////////////
        // Reduced volume cut //
        ////////////////////////

        if (!isWithinReducedVolume(breakPointX, breakPointY, breakPointZ)) {
            if (ev.backgroundType == 0 || ev.backgroundType == 1) {
                hTrueAbsKERejected->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
                hTrueAbsKERejRedVol->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
            } else if (ev.backgroundType == 6 || ev.backgroundType == 12) {
                hTrueScatterKERejected->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
                hTrueScatterKERejRedVol->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
            } else if (ev.backgroundType == 7) {
                hTrueChExchKERejected->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
            }
            continue;
        }
        hPrimaryInRedVol->Fill(ev.backgroundType, ev.weight);

        /////////////////////
        // Primary PID cut //
        /////////////////////

        if (minChi2 == pionChi2 || minChi2 == protonChi2) {
            if (ev.backgroundType == 0 || ev.backgroundType == 1) {
                hTrueAbsKERejected->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
                hTrueAbsKERejPID->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
            } else if (ev.backgroundType == 6 || ev.backgroundType == 12) {
                hTrueScatterKERejected->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
                hTrueScatterKERejPID->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
            } else if (ev.backgroundType == 7) {
                hTrueChExchKERejected->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
            }
            continue;
        }
        hPrimaryPID->Fill(ev.backgroundType, ev.weight);

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

        if (totalTaggedPions > 0) {
            if (totalTaggedPions > 1 || newSecondaryPion) {
                // reject events with > 1 tagged pion
                if (ev.backgroundType == 0 || ev.backgroundType == 1) {
                    hTrueAbsKERejected->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
                    hTrueAbsKERejManyPions->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
                } else if (ev.backgroundType == 6 || ev.backgroundType == 12) {
                    hTrueScatterKERejected->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
                    hTrueScatterKERejManyPions->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
                } else if (ev.backgroundType == 7) {
                    hTrueChExchKERejected->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
                }
                continue;
            }

            // Select as scatter
            hPionScatter->Fill(ev.backgroundType, ev.weight);

            hPionScatterKE->Fill(energyAtVertex, ev.weight);
            if (ev.backgroundType == 0 || ev.backgroundType == 1) {
                hPionScatterKEAbs->Fill(energyAtVertex, ev.weight);
            } else if (ev.backgroundType == 6 || ev.backgroundType == 12) {
                hPionScatterKETrue->Fill(energyAtVertex, ev.weight);
                outFileEvents << "True scatter as scatter: " << run << " " << subrun << " " << event << std::endl;
            } else if (ev.backgroundType == 2) {
                hPionScatterKEMuon->Fill(energyAtVertex, ev.weight);

                if (muonType == 0) {
                    hPionScatterKEMuonTG->Fill(energyAtVertex, ev.weight);
                } else if (muonType == 1) {
                    hPionScatterKEMuonDecay->Fill(energyAtVertex, ev.weight);
                } else if (muonType == 2) {
                    hPionScatterKEMuonCAR->Fill(energyAtVertex, ev.weight);
                }
            } else if (ev.backgroundType == 3) {
                hPionScatterKEElectron->Fill(energyAtVertex, ev.weight);
            } else if (ev.backgroundType == 7) {
                hPionScatterKEChExch->Fill(energyAtVertex, ev.weight);
            } else {
                hPionScatterKEOther->Fill(energyAtVertex, ev.weight);
            }

            if (ev.backgroundType == 0 || ev.backgroundType == 1) {
                hTrueAbsKEAsScatter->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
                if (TrueEnergyBin != -1) TrueAbsAsByBin.at(TrueEnergyBin).at(1)->Fill(energyAtVertex, ev.weight);
            } else if (ev.backgroundType == 6 || ev.backgroundType == 12) {
                hTrueScatterKEAsScatter->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
                if (TrueEnergyBin != -1) TrueScatterAsByBin.at(TrueEnergyBin).at(1)->Fill(energyAtVertex, ev.weight);
            } else if (ev.backgroundType == 7) {
                hTrueChExchKEAsScatter->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
            }

            // TODO: maybe check not rejecting more than 1 pion

            continue;
        }
        hNotScatter->Fill(ev.backgroundType, ev.weight);

        if (totalTaggedProtons > 0) {
            // Select as absorption (Np)
            hPionAbs->Fill(ev.backgroundType, ev.weight);

            hPionAbsKE->Fill(energyAtVertex, ev.weight);
            if (ev.backgroundType == 0 || ev.backgroundType == 1) {
                hPionAbsKETrue->Fill(energyAtVertex, ev.weight);
                outFileEvents << "True abs np as abs: " << run << " " << subrun << " " << event << std::endl;
            } else if (ev.backgroundType == 6 || ev.backgroundType == 12) {
                hPionAbsKEScatter->Fill(energyAtVertex, ev.weight);
            } else if (ev.backgroundType == 2) {
                hPionAbsKEMuon->Fill(energyAtVertex, ev.weight);
            } else if (ev.backgroundType == 3) {
                hPionAbsKEElectron->Fill(energyAtVertex, ev.weight);
            } else if (ev.backgroundType == 7) {
                hPionAbsKEChExch->Fill(energyAtVertex, ev.weight);
            } else {
                hPionAbsKEOther->Fill(energyAtVertex, ev.weight);
            }

            if (ev.backgroundType == 0 || ev.backgroundType == 1) {
                hTrueAbsKEAsAbs->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
                if (TrueEnergyBin != -1) TrueAbsAsByBin.at(TrueEnergyBin).at(0)->Fill(energyAtVertex, ev.weight);
            } else if (ev.backgroundType == 6 || ev.backgroundType == 12) {
                hTrueScatterKEAsAbs->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
                if (TrueEnergyBin != -1) TrueScatterAsByBin.at(TrueEnergyBin).at(0)->Fill(energyAtVertex, ev.weight);
            } else if (ev.backgroundType == 7) {
                hTrueChExchKEAsAbs->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
            }

            continue;
        }
        hNotPionAbsNp->Fill(ev.backgroundType, ev.weight);

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

            if (ev.backgroundType == 0 || ev.backgroundType == 1) {
                if (hitClusters[i].plane == 0) hHitClusterInductionSizesAbs->Fill(clusterSize, ev.weight);
                else if (hitClusters[i].plane == 1) hHitClusterCollectionSizesAbs->Fill(clusterSize, ev.weight);
            } else if (ev.backgroundType == 2) {
                if (hitClusters[i].plane == 0) hHitClusterInductionSizesMuon->Fill(clusterSize, ev.weight);
                else if (hitClusters[i].plane == 1) hHitClusterCollectionSizesMuon->Fill(clusterSize, ev.weight);
            } else if (ev.backgroundType == 3) {
                if (hitClusters[i].plane == 0) hHitClusterInductionSizesElectron->Fill(clusterSize, ev.weight);
                else if (hitClusters[i].plane == 1) hHitClusterCollectionSizesElectron->Fill(clusterSize, ev.weight);
            } else if (ev.backgroundType == 6 || ev.backgroundType == 12) {
                if (hitClusters[i].plane == 0) hHitClusterInductionSizesScatter->Fill(clusterSize, ev.weight);
                else if (hitClusters[i].plane == 1) hHitClusterCollectionSizesScatter->Fill(clusterSize, ev.weight);
            } else if (ev.backgroundType == 7) {
                if (hitClusters[i].plane == 0) hHitClusterInductionSizesChExch->Fill(clusterSize, ev.weight);
                else if (hitClusters[i].plane == 1) hHitClusterCollectionSizesChExch->Fill(clusterSize, ev.weight);
            } else {
                if (hitClusters[i].plane == 0) hHitClusterInductionSizesOther->Fill(clusterSize, ev.weight);
                else if (hitClusters[i].plane == 1) hHitClusterCollectionSizesOther->Fill(clusterSize, ev.weight);
            }
        }

        if (ev.backgroundType == 0 || ev.backgroundType == 1) {
            hLargeHitClusterInductionAbs->Fill(numLargeClustersInduction, ev.weight);
            hLargestHitClusterInductionAbs->Fill(largestClusterSizeInduction, ev.weight);
            hNumClustersInductionAbs->Fill(numClustersInduction, ev.weight);

            hLargeHitClusterCollectionAbs->Fill(numLargeClustersCollection, ev.weight);
            hLargestHitClusterCollectionAbs->Fill(largestClusterSizeCollection, ev.weight);
            hNumClustersCollectionAbs->Fill(numClustersCollection, ev.weight);
        } else if (ev.backgroundType == 7) {
            hLargeHitClusterInductionChExch->Fill(numLargeClustersInduction, ev.weight);
            hLargestHitClusterInductionChExch->Fill(largestClusterSizeInduction, ev.weight);
            hNumClustersInductionChExch->Fill(numClustersInduction, ev.weight);

            hLargeHitClusterCollectionChExch->Fill(numLargeClustersCollection, ev.weight);
            hLargestHitClusterCollectionChExch->Fill(largestClusterSizeCollection, ev.weight);
            hNumClustersCollectionChExch->Fill(numClustersCollection, ev.weight);
        } else if (ev.backgroundType == 6 || ev.backgroundType == 12) {
            hLargeHitClusterInductionScatter->Fill(numLargeClustersInduction, ev.weight);
            hLargestHitClusterInductionScatter->Fill(largestClusterSizeInduction, ev.weight);
            hNumClustersInductionScatter->Fill(numClustersInduction, ev.weight);

            hLargeHitClusterCollectionScatter->Fill(numLargeClustersCollection, ev.weight);
            hLargestHitClusterCollectionScatter->Fill(largestClusterSizeCollection, ev.weight);
            hNumClustersCollectionScatter->Fill(numClustersCollection, ev.weight);
        } else if (ev.backgroundType == 2) {
            hLargeHitClusterInductionMuon->Fill(numLargeClustersInduction, ev.weight);
            hLargestHitClusterInductionMuon->Fill(largestClusterSizeInduction, ev.weight);
            hNumClustersInductionMuon->Fill(numClustersInduction, ev.weight);

            hLargeHitClusterCollectionMuon->Fill(numLargeClustersCollection, ev.weight);
            hLargestHitClusterCollectionMuon->Fill(largestClusterSizeCollection, ev.weight);
            hNumClustersCollectionMuon->Fill(numClustersCollection, ev.weight);
        } else if (ev.backgroundType == 3) {
            hLargeHitClusterInductionElectron->Fill(numLargeClustersInduction, ev.weight);
            hLargestHitClusterInductionElectron->Fill(largestClusterSizeInduction, ev.weight);
            hNumClustersInductionElectron->Fill(numClustersInduction, ev.weight);

            hLargeHitClusterCollectionElectron->Fill(numLargeClustersCollection, ev.weight);
            hLargestHitClusterCollectionElectron->Fill(largestClusterSizeCollection, ev.weight);
            hNumClustersCollectionElectron->Fill(numClustersCollection, ev.weight);
        } else {
            hLargeHitClusterInductionOther->Fill(numLargeClustersInduction, ev.weight);
            hLargestHitClusterInductionOther->Fill(largestClusterSizeInduction, ev.weight);
            hNumClustersInductionOther->Fill(numClustersInduction, ev.weight);

            hLargeHitClusterCollectionOther->Fill(numLargeClustersCollection, ev.weight);
            hLargestHitClusterCollectionOther->Fill(largestClusterSizeCollection, ev.weight);
            hNumClustersCollectionOther->Fill(numClustersCollection, ev.weight);
        }

        if (
            // numClustersInduction < MAX_NUM_CLUSTERS_INDUCTION
            // numLargeClustersCollection < MAX_NUM_LARGE_CLUSTERS_COLLECTION
            numLargeClustersInduction < MAX_NUM_LARGE_CLUSTERS_INDUCTION
        ) {
            hPionAbs->Fill(ev.backgroundType, ev.weight);

            hPionAbsKE->Fill(energyAtVertex, ev.weight);
            if (ev.backgroundType == 0 || ev.backgroundType == 1) {
                hPionAbsKETrue->Fill(energyAtVertex, ev.weight);
                outFileEvents << "True abs 0p as abs: " << run << " " << subrun << " " << event << std::endl;
            } else if (ev.backgroundType == 7) {
                hPionAbsKEChExch->Fill(energyAtVertex, ev.weight);
            } else if (ev.backgroundType == 6 || ev.backgroundType == 12) {
                hPionAbsKEScatter->Fill(energyAtVertex, ev.weight);
            } else if (ev.backgroundType == 2) {
                hPionAbsKEMuon->Fill(energyAtVertex, ev.weight);
                if (muonType == 0) {
                    hPionAbsKEMuonTG->Fill(energyAtVertex, ev.weight);
                } else if (muonType == 1) {
                    hPionAbsKEMuonDecay->Fill(energyAtVertex, ev.weight);
                } else if (muonType == 2) {
                    hPionAbsKEMuonCAR->Fill(energyAtVertex, ev.weight);
                }
            } else if (ev.backgroundType == 3) {
                hPionAbsKEElectron->Fill(energyAtVertex, ev.weight);
            } else {
                hPionAbsKEOther->Fill(energyAtVertex, ev.weight);
            }

            if (ev.backgroundType == 0 || ev.backgroundType == 1) {
                hTrueAbsKEAsAbs->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
                if (TrueEnergyBin != -1) TrueAbsAsByBin.at(TrueEnergyBin).at(0)->Fill(energyAtVertex, ev.weight);
            } else if (ev.backgroundType == 6 || ev.backgroundType == 12) {
                hTrueScatterKEAsAbs->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
                if (TrueEnergyBin != -1) TrueScatterAsByBin.at(TrueEnergyBin).at(0)->Fill(energyAtVertex, ev.weight);
            } else if (ev.backgroundType == 7) {
                hTrueChExchKEAsAbs->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
            }

            continue;
        }
        hNotPionAbs0p->Fill(ev.backgroundType, ev.weight);

        if (ev.backgroundType == 0 || ev.backgroundType == 1) {
            hUnRecoHitsInductionAbs->Fill(numUnRecoHitsNearPrimaryInduction, ev.weight);
            hUnRecoHitsCollectionAbs->Fill(numUnRecoHitsNearPrimaryCollection, ev.weight);
        } else if (ev.backgroundType == 2) {
            hUnRecoHitsInductionMuon->Fill(numUnRecoHitsNearPrimaryInduction, ev.weight);
            hUnRecoHitsCollectionMuon->Fill(numUnRecoHitsNearPrimaryCollection, ev.weight);
        } else if (ev.backgroundType == 3) {
            hUnRecoHitsInductionElectron->Fill(numUnRecoHitsNearPrimaryInduction, ev.weight);
            hUnRecoHitsCollectionElectron->Fill(numUnRecoHitsNearPrimaryCollection, ev.weight);
        } else if (ev.backgroundType == 6 || ev.backgroundType == 12) {
            hUnRecoHitsInductionScatter->Fill(numUnRecoHitsNearPrimaryInduction, ev.weight);
            hUnRecoHitsCollectionScatter->Fill(numUnRecoHitsNearPrimaryCollection, ev.weight);
        } else if (ev.backgroundType == 7) {
            hUnRecoHitsInductionChExch->Fill(numUnRecoHitsNearPrimaryInduction, ev.weight);
            hUnRecoHitsCollectionChExch->Fill(numUnRecoHitsNearPrimaryCollection, ev.weight);
        } else {
            hUnRecoHitsInductionOther->Fill(numUnRecoHitsNearPrimaryInduction, ev.weight);
            hUnRecoHitsCollectionOther->Fill(numUnRecoHitsNearPrimaryCollection, ev.weight);
        }

        // Anything left here is rejected
        if (ev.backgroundType == 0 || ev.backgroundType == 1) {
            hTrueAbsKERejected->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
            hTrueAbsKERejClusters->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
        } else if (ev.backgroundType == 6 || ev.backgroundType == 12) {
            hTrueScatterKERejected->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
            hTrueScatterKERejClusters->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
        } else if (ev.backgroundType == 7) {
            hTrueChExchKERejected->Fill(ev.truthPrimaryVertexKE * 1000, ev.weight);
        }
    }

    std::cout << std::endl;
    std::cout << "Events with reconstructed primary track: " << numEventsPrimaryReco << std::endl;
    std::cout << "Events with no primary hits in induction: " << numEventsNoInduction << std::endl;
    std::cout << "Events with no primary hits in collection: " << numEventsNoCollection << std::endl;
    std::cout << "Events with no primary hits in either: " << numEventsNoEither << std::endl;

    std::cout << std::endl;
    std::cout << "Events with primary reconstructed with negative time hits: " << numEventsPrimaryHitNegativeTime << std::endl;

    std::cout << std::endl;
    std::cout << "Original sample composition: " << std::endl;
    printBackgroundInfo(hTotalEvents, std::cout);

    std::cout << std::endl;
    std::cout << "Original muon event composition: " << std::endl;
    std::cout << "  Total muon events: " << hTrueMuonTypes->Integral() << std::endl;
    std::cout << "  Through-going:     " << hTrueMuonTypes->GetBinContent(1) << std::endl;
    std::cout << "  Decay:             " << hTrueMuonTypes->GetBinContent(2) << std::endl;
    std::cout << "  Capture at rest:   " << hTrueMuonTypes->GetBinContent(3) << std::endl;
    std::cout << "  Other:             " << hTrueMuonTypes->GetBinContent(4) << std::endl;

    std::cout << std::endl;
    std::cout << "Number of scattering events modified: " << scatteringsModified << std::endl;

    std::cout << std::endl;
    std::cout << "Number of events with at most 0 TG tracks: " << eventCount0TG << std::endl;
    std::cout << "Number of events with at most 1 TG tracks: " << eventCount1TG << std::endl;
    std::cout << "Number of events with at most 2 TG tracks: " << eventCount2TG << std::endl;

    //////////////////////////////////
    // Breakdown of rejected events //
    //////////////////////////////////

    std::cout << std::endl;
    std::cout << "=======================================================" << std::endl;
    std::cout << "Rejection breakdown" << std::endl;
    std::cout << "=======================================================" << std::endl;

    std::cout << std::endl;
    std::cout << "Events with sufficiency data products and WC to TPC match: " << std::endl;
    printBackgroundInfo(hDataProdsAndWC2TPC, std::cout);

    std::cout << std::endl;
    std::cout << "Events with number of TG tracks below threshold: " << std::endl;
    printBackgroundInfo(hNotManyTGTracks, std::cout);

    std::cout << std::endl;
    std::cout << "Events not classified as electron (with cylinder small tracks): " << std::endl;
    printBackgroundInfo(hNotAnElectron, std::cout);

    std::cout << std::endl;
    std::cout << "Events with vertex identified inside reduced volume: " << std::endl;
    printBackgroundInfo(hPrimaryInRedVol, std::cout);

    std::cout << std::endl;
    std::cout << "Events with primary matched to a negative pion: " << std::endl;
    printBackgroundInfo(hPrimaryPID, std::cout);

    std::cout << std::endl;
    std::cout << "Events not classified as scatter (no secondary pion): " << std::endl;
    printBackgroundInfo(hNotScatter, std::cout);

    std::cout << std::endl;
    std::cout << "Events not classified as pion Np absorption (no secondary particles): " << std::endl;
    printBackgroundInfo(hNotPionAbsNp, std::cout);

    std::cout << std::endl;
    std::cout << "Events not classified as pion 0p absorption (too many non-reconstructed hit clusters): " << std::endl;
    printBackgroundInfo(hNotPionAbs0p, std::cout);

    ///////////////////////////////////////////
    // Information about events at each step //
    ///////////////////////////////////////////

    std::cout << std::endl;
    std::cout << "=======================================================" << std::endl;
    std::cout << "Final classification breakdown" << std::endl;
    std::cout << "=======================================================" << std::endl;

    std::cout << std::endl;
    
    std::cout << "Pion abs total reco " <<  hPionAbs->Integral() << " with composition: " << std::endl;
    printBackgroundInfo(hPionAbs, std::cout);
    std::cout << "Purity: " << (hPionAbs->GetBinContent(1) + hPionAbs->GetBinContent(2)) / hPionAbs->Integral() << std::endl;
    std::cout << "Efficiency: " << (hPionAbs->GetBinContent(1) + hPionAbs->GetBinContent(2)) / (hTotalEvents->GetBinContent(1) + hTotalEvents->GetBinContent(2)) << std::endl;

    std::cout << std::endl;
    std::cout << "Pion scattering total reco " << hPionScatter->Integral() << " with composition: " << std::endl;
    printBackgroundInfo(hPionScatter, std::cout);
    std::cout << "Purity: " << (hPionScatter->GetBinContent(7) + hPionScatter->GetBinContent(13)) / hPionScatter->Integral() << std::endl;
    std::cout << "Efficiency: " << (hPionScatter->GetBinContent(7) + hPionScatter->GetBinContent(13)) / (hTotalEvents->GetBinContent(7) + hTotalEvents->GetBinContent(13)) << std::endl;


    ///////////////////////////////////
    // Corrections for incident flux //
    ///////////////////////////////////

    // Set errors to original incident KE histogram
    for (int iBin = 1; iBin <= hIncidentKE->GetNbinsX(); ++iBin) {
        double entries = hIncidentKE->GetBinContent(iBin);
        hIncidentKE->SetBinError(iBin, std::sqrt(entries)); // stat error
    }

    // Compute corrections for use when computing cross-section
    TH1D* hPsiInc = (TH1D*) hIncidentKEPion->Clone("hPsiInc");
    hPsiInc->Divide(hTrueIncidentKE);

    TH1D* hCInc = (TH1D*) hIncidentKEPion->Clone("hCInc");
    hCInc->Divide(hIncidentKE);

    // Get corrected incident KE
    TH1D* hIncidentKECorrected = (TH1D*) hIncidentKE->Clone("hIncidentKECorrected");

    // Set errors to corrected incident KE histogram
    for (int iBin = 1; iBin <= hIncidentKECorrected->GetNbinsX(); ++iBin) {
        double entries = hIncidentKECorrected->GetBinContent(iBin);
        hIncidentKECorrected->SetBinError(iBin, std::sqrt(entries));
    }

    // Apply corrections
    hIncidentKECorrected->Divide(hPsiInc);
    hIncidentKECorrected->Multiply(hCInc);

    //////////////////////////////////////////////
    // Perform unfolding for interacting slices //
    //////////////////////////////////////////////

    // Response matrix with components: 
    //    R[(i, \alpha), (j, \beta)] = P(reco signal \beta in reco energy bin j | true signal \alpha in true energy bin i)
    TH2D* hResponseMatrix = new TH2D(
        "hResponseMatrix", "Response;Reco (j, #beta);True (i, #alpha)",
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE,
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    );

    // Construct response matrix
    GetResponseMatrix(
        NUM_SIGNAL_TYPES, NUM_BINS_KE,
        TotalEventsHistos,
        TrueRecoAsByBin,
        hIncidentKE,
        hTrueIncidentKE,
        hResponseMatrix
    );

    // Construct large background vector
    TVectorD Background(NUM_SIGNAL_TYPES * NUM_BINS_KE);
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            double reco = RecoSignals[iSignal]->GetBinContent(iBin + 1);
            double inc  = hIncidentKE->GetBinContent(iBin + 1);
            double bkg  = 0;
            for (int iBkg = 0; iBkg < RecoSignalBackgrounds[iSignal].size(); ++iBkg) {
                bkg += RecoSignalBackgrounds[iSignal][iBkg]->GetBinContent(iBin + 1);
            }
            int index = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            Background(index) = XSEC_UNITS * (bkg / inc);
        }
    }

    // Construct large measured vector
    TVectorD Measure(NUM_SIGNAL_TYPES * NUM_BINS_KE);
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index      = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            Measure(index) = XSEC_UNITS * (RecoSignals[iSignal]->GetBinContent(iBin + 1) / hIncidentKE->GetBinContent(iBin + 1));
        }
    }

    // Construct large signal vector
    TVectorD Signal(NUM_SIGNAL_TYPES * NUM_BINS_KE);
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index     = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            Signal(index) = XSEC_UNITS * (TotalEventsHistos[iSignal]->GetBinContent(iBin + 1) / hTrueIncidentKE->GetBinContent(iBin + 1));
        }
    }

    // Save nominal measure and signal
    TH1D* hSignalVectorNominal = new TH1D(
        "hSignalVectorNominal", "hSignalVectorNominal;;", 
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    ); V2H(Signal, hSignalVectorNominal);
    TH1D* hMeasureVectorNominal = new TH1D(
        "hMeasureVectorNominal", "hMeasureVectorNominal;;", 
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    ); V2H(Measure, hMeasureVectorNominal);
    TH1D* hBackgroundVectorNominal = new TH1D(
        "hBackgroundVectorNominal", "hBackgroundVectorNominal;;", 
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    ); V2H(Background, hBackgroundVectorNominal);

    // Subtract background from measured
    TVectorD MeasureMinusBackground = Measure - Background;

    // Construct statistical covariance matrix
    TMatrixD StatCovariance(NUM_SIGNAL_TYPES * NUM_BINS_KE, NUM_SIGNAL_TYPES * NUM_BINS_KE); StatCovariance.Zero();
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index   = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            double Ninc = hIncidentKE->GetBinContent(iBin + 1);
            double N    = Measure(index) * Ninc / XSEC_UNITS;

            // Scale to data to estimate stat systematics
            Ninc = Ninc * MC_SCALING;
            N    = N * MC_SCALING;

            double numSigma = (Ninc>0.0 && N>0.0 ? std::sqrt(N*(1.0 - N/Ninc)) : 0.0);
            double denSigma = (Ninc>0.0 ? std::sqrt(Ninc) : 0.0);

            double xs       = (N / Ninc) * XSEC_UNITS;
            double relVarXS = (N>0.0 && Ninc>0.0 ? std::pow(numSigma/N + denSigma/Ninc, 2) : 0.0);

            StatCovariance(index, index) = xs * xs * relVarXS;
        }
    }

    // Add all covariances (in this script, only care about MC stat)
    TMatrixD Covariance(NUM_SIGNAL_TYPES * NUM_BINS_KE, NUM_SIGNAL_TYPES * NUM_BINS_KE); Covariance.Zero();
    Covariance += StatCovariance;

    // Convert response histogram to matrix
    TMatrixD Response(NUM_SIGNAL_TYPES * NUM_BINS_KE, NUM_SIGNAL_TYPES * NUM_BINS_KE); 
    H2M(static_cast<const TH2D*>(hResponseMatrix), Response, kTRUE);

    // Objects to store stuff
    TMatrixD AddSmear(NUM_SIGNAL_TYPES * NUM_BINS_KE, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    TMatrixD AddSmearInverse(NUM_SIGNAL_TYPES * NUM_BINS_KE, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    TVectorD WF(NUM_SIGNAL_TYPES * NUM_BINS_KE);
    TMatrixD UnfoldCov(NUM_SIGNAL_TYPES * NUM_BINS_KE, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    TMatrixD CovRotation(NUM_SIGNAL_TYPES * NUM_BINS_KE, NUM_SIGNAL_TYPES * NUM_BINS_KE);

    // Invert response matrix with Wiener-SVD w/o regularization
    TVectorD UnfoldedReco = WienerSVD(
        Response,
        Signal,
        MeasureMinusBackground,
        Covariance,
        0, // 0: unit, 2: second derivative
        0, // smoothing parameter
        AddSmear,
        WF,
        UnfoldCov,
        CovRotation,
        AddSmearInverse
    );
    TH1D* hUnfReco = new TH1D("hUnfReco", "Unfolded Reco;;", NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    V2H(UnfoldedReco, hUnfReco);

    // Sanity check: CovRotation = RInv?
    TDecompSVD svd(Response); svd.SetTol(1e-8);
    TMatrixD Rinv = svd.Invert();
    if (!EqualApprox(CovRotation, Rinv, 1e-8, 1e-8)) {
        std::cerr << "Warning: CovRotation matrix does not match the inverse of the response matrix within tolerance." << std::endl;
    }

    // Organize results
    std::vector<TH1*> UnfoldedRecoHistos;
    for (int iBin = 0; iBin < NUM_SIGNAL_TYPES; ++iBin) {
        TH1D* unfoldedHist = new TH1D(Form("hUnfoldedRecoVec_%d", iBin), Form("Unfolded Reco Vec %d", iBin), NUM_BINS_KE, ARRAY_KE_BINS.data());
        UnfoldedRecoHistos.push_back(unfoldedHist);
    }

    for (int iFlatIndex = 0; iFlatIndex < NUM_SIGNAL_TYPES * NUM_BINS_KE; ++iFlatIndex) {
        auto [signalBin, energyBin] = unflattenIndex(iFlatIndex, NUM_BINS_KE);

        // Get error for unfolded reco histogram from unfolded covariance
        double err = std::sqrt(UnfoldCov[iFlatIndex][iFlatIndex]);
        UnfoldedRecoHistos[signalBin]->SetBinContent(energyBin + 1, UnfoldedReco(iFlatIndex));
        UnfoldedRecoHistos[signalBin]->SetBinError(energyBin + 1, err);
    }

    // Copy covariance and inverse response matrices into histos
    TH2D* hCovariance = new TH2D(
        "Covariance", "Covariance Matrix;(i, #alpha);(j, #beta)",
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE, 
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE 
    ); M2H(Covariance, hCovariance);
    TH2D* hUnfCovariance = new TH2D(
        "Unfolded Covariance", "Unfolded Covariance Matrix;(i, #alpha);(j, #beta)",
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE, 
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE 
    ); M2H(UnfoldCov, hUnfCovariance);
    TH2D* hResponseInvMatrix = new TH2D("Response Inverse", "Response Inverse;Reco (j, #beta);True (i, #alpha)",
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE,
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    ); M2H(CovRotation, hResponseInvMatrix);

    TH2D* hStatCovariance = new TH2D(
        "hStatCovariance", "Statistical Covariance Matrix;(i, #alpha);(j, #beta)",
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE, 
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE 
    ); M2H(StatCovariance, hStatCovariance);

    // Save stuff to nominal file
    nominalFile->cd();
    hMeasureVectorNominal->Write("", TObject::kOverwrite);
    hBackgroundVectorNominal->Write("", TObject::kOverwrite);
    hSignalVectorNominal->Write("", TObject::kOverwrite);
    hResponseMatrix->Write("", TObject::kOverwrite);
    hStatCovariance->Write("", TObject::kOverwrite);
    nominalFile->Close();

    ///////////////////////////////////////////////
    // Get cross-section using corrrected fluxes //
    ///////////////////////////////////////////////

    TH1D* hPionAbsCrossSection     = new TH1D("hPionAbsCrossSection", "hPionAbsCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTruePionAbsCrossSection = new TH1D("hTruePionAbsCrossSection", "hTruePionAbsCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    TH1D* hPionScatterCrossSection     = new TH1D("hPionScatterCrossSection", "hPionScatterCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTruePionScatterCrossSection = new TH1D("hTruePionScatterCrossSection", "hTruePionScatterCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    TH1D* hTruePionChExchCrossSection = new TH1D("hTruePionChExchCrossSection", "hTruePionChExchCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTruePionOtherCrossSection  = new TH1D("hTruePionOtherCrossSection", "hTruePionOtherCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueTotalCrossSection      = new TH1D("hTrueTotalCrossSection", "hTrueTotalCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    // TODO: some sort of check with the total?

    TH1D* hTruePionAbsNoThinCrossSection     = new TH1D("hTruePionAbsNoThinCrossSection", "hTruePionAbsNoThinCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTruePionScatterNoThinCrossSection = new TH1D("hTruePionScatterNoThinCrossSection", "hTruePionScatterNoThinCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    std::vector<TH1D*> UnfoldedCrossSections = {
        hPionAbsCrossSection,
        hPionScatterCrossSection
    };

    std::vector<TH1D*> TrueCrossSections = {
        hTruePionAbsCrossSection,
        hTruePionScatterCrossSection
    };

    std::vector<TH1D*> TrueNoThinCrossSections = {
        hTruePionAbsNoThinCrossSection,
        hTruePionScatterNoThinCrossSection
    };

    for (int i = 0; i < UnfoldedCrossSections.size(); ++i) {
        TH1D* unfXSec  = UnfoldedCrossSections[i];
        TH1D* trueXSec = TrueCrossSections[i];
        TH1D* trueNoThinXSec = TrueNoThinCrossSections[i];

        for (int iBin = 1; iBin <= NUM_BINS_KE; ++iBin) {
            double xsecErr     = UnfoldedRecoHistos[i]->GetBinError(iBin);
            double xsecContent = UnfoldedRecoHistos[i]->GetBinContent(iBin);

            unfXSec->SetBinContent(iBin, xsecContent);
            unfXSec->SetBinError(iBin, xsecErr);

            // True cross-section, no error bars
            trueXSec->SetBinContent(iBin, XSEC_UNITS * (TotalEventsHistos[i]->GetBinContent(iBin) / hTrueIncidentKE->GetBinContent(iBin)));
            trueNoThinXSec->SetBinContent(iBin, XSEC_UNITS * (-std::log(1 - (TotalEventsHistos[i]->GetBinContent(iBin) / hTrueIncidentKE->GetBinContent(iBin)))));
        }

        // Make contents per 100 MeV
        reweightOneDHisto(unfXSec, 100.);
        reweightOneDHisto(trueXSec, 100.);
        reweightOneDHisto(trueNoThinXSec, 100.);
    }

    // True "other" cross-section
    for (int iBin = 1; iBin <= NUM_BINS_KE; ++iBin) {
        hTruePionOtherCrossSection->SetBinContent(iBin, XSEC_UNITS * (hTrueOtherKE->GetBinContent(iBin) / hTrueIncidentKE->GetBinContent(iBin)));
        hTruePionChExchCrossSection->SetBinContent(iBin, XSEC_UNITS * (hTrueChExchKE->GetBinContent(iBin) / hTrueIncidentKE->GetBinContent(iBin)));
    }
    reweightOneDHisto(hTruePionOtherCrossSection, 100.); reweightOneDHisto(hTruePionChExchCrossSection, 100.);

    ///////////////////////////
    // Make plots per 50 MeV //
    ///////////////////////////

    // TODO: create function that scales histogram bins per 50 MeV,
    //       this is already done above for cross-sections, but have to do 
    //       it for everything

    ///////////////////////////////////
    // Save histograms to text files //
    ///////////////////////////////////

    histoToText(hTrueAllKE, "/exp/lariat/app/users/epelaez/analysis_abs_scatt/files/Flux/InteractingFlux.txt");
    histoToText(hTrueIncidentKE, "/exp/lariat/app/users/epelaez/analysis_abs_scatt/files/Flux/IncidentFlux.txt");

    //////////////////
    // Create plots //
    //////////////////

    std::vector<std::vector<TH1*>> PlotGroups = {
        // Incident KE
        {hIncidentKE, hIncidentKEPion, hIncidentKEMuon, hIncidentKEElectron},
        {hIncidentKEFine, hIncidentKEPionFine, hIncidentKEMuonFine, hIncidentKEElectronFine},
        {hIncidentKECorrected, hTrueIncidentKE},
        {hFrontFaceKE, hFrontFaceKEPion, hFrontFaceKEMuon, hFrontFaceKEElectron},
        {hFrontFaceKENoWeight, hFrontFaceKEPionNoWeight, hFrontFaceKEMuonNoWeight, hFrontFaceKEElectronNoWeight},
        {hFrontFaceKENoWeightPre, hFrontFaceKEPionNoWeightPre, hFrontFaceKEMuonNoWeightPre, hFrontFaceKEElectronNoWeightPre},
        {hFrontFacePionKE, hFrontFacePionKEPion, hFrontFacePionKEElectron, hFrontFacePionKEMuon},
        {hWCKE, hWCKEPion, hWCKEMuon, hWCKEElectron},

        // Interacting KE
        {hPionAbsKE, hPionAbsKETrue, hPionAbsKEScatter, hPionAbsKEChExch, hPionAbsKEMuon, hPionAbsKEElectron, hPionAbsKEOther},
        {hPionAbsKE, hPionAbsKETrue, hPionAbsKEScatter, hPionAbsKEChExch, hPionAbsKEMuonTG, hPionAbsKEMuonDecay, hPionAbsKEMuonCAR, hPionAbsKEElectron, hPionAbsKEOther},
        {hPionScatterKE, hPionScatterKETrue, hPionScatterKEAbs, hPionScatterKEChExch, hPionScatterKEMuon, hPionScatterKEElectron, hPionScatterKEOther},
        {hPionScatterKE, hPionScatterKETrue, hPionScatterKEAbs, hPionScatterKEChExch, hPionScatterKEMuonTG, hPionScatterKEMuonDecay, hPionScatterKEMuonCAR, hPionScatterKEElectron, hPionScatterKEOther},

        // True interacting KE
        {hTrueAllKE},
        {hTrueAbsKE, hTrueScatterKE, hTrueChExchKE, hTrueOtherKE},

        // True events classified breakdown
        {hTrueAbsKEAsAbs, hTrueAbsKEAsScatter, hTrueAbsKERejected},
        {hTrueScatterKEAsAbs, hTrueScatterKEAsScatter, hTrueScatterKERejected},
        {hTrueChExchKEAsAbs, hTrueChExchKEAsScatter, hTrueChExchKERejected},

        // Rejected events
        {hTrueAbsKERejDataProds, hTrueAbsKERejElectron, hTrueAbsKERejRedVol, hTrueAbsKERejPID, hTrueAbsKERejManyPions, hTrueAbsKERejClusters},
        {hTrueScatterKERejDataProds, hTrueScatterKERejElectron, hTrueScatterKERejRedVol, hTrueScatterKERejPID, hTrueScatterKERejManyPions, hTrueScatterKERejClusters},

        // Cross-sections (unfolded)
        {hTruePionAbsCrossSection, hPionAbsCrossSection},
        {hTruePionScatterCrossSection, hPionScatterCrossSection},

        // Total true-cross section
        {hTruePionAbsCrossSection, hTruePionScatterCrossSection, hTruePionChExchCrossSection, hTruePionOtherCrossSection},

        // Total true-cross section with no thin-slice approximation
        {hTruePionAbsCrossSection, hTruePionAbsNoThinCrossSection},
        {hTruePionScatterCrossSection, hTruePionScatterNoThinCrossSection},

        // Unreconstructed hits
        {hHitClusterInductionSizesAbs, hHitClusterInductionSizesMuon, hHitClusterInductionSizesElectron, hHitClusterInductionSizesScatter, hHitClusterInductionSizesChExch, hHitClusterInductionSizesOther},
        {hLargeHitClusterInductionAbs, hLargeHitClusterInductionMuon, hLargeHitClusterInductionElectron, hLargeHitClusterInductionScatter, hLargeHitClusterInductionChExch, hLargeHitClusterInductionOther},
        {hLargestHitClusterInductionAbs, hLargestHitClusterInductionMuon, hLargestHitClusterInductionElectron, hLargestHitClusterInductionScatter, hLargestHitClusterInductionChExch, hLargestHitClusterInductionOther},
        {hNumClustersInductionAbs, hNumClustersInductionMuon, hNumClustersInductionElectron, hNumClustersInductionScatter, hNumClustersInductionChExch, hNumClustersInductionOther},
        {hUnRecoHitsInductionAbs, hUnRecoHitsInductionMuon, hUnRecoHitsInductionElectron, hUnRecoHitsInductionScatter, hUnRecoHitsInductionChExch, hUnRecoHitsInductionOther},

        {hHitClusterCollectionSizesAbs, hHitClusterCollectionSizesMuon, hHitClusterCollectionSizesElectron, hHitClusterCollectionSizesScatter, hHitClusterCollectionSizesChExch, hHitClusterCollectionSizesOther},
        {hLargeHitClusterCollectionAbs, hLargeHitClusterCollectionMuon, hLargeHitClusterCollectionElectron, hLargeHitClusterCollectionScatter, hLargeHitClusterCollectionChExch, hLargeHitClusterCollectionOther},
        {hLargestHitClusterCollectionAbs, hLargestHitClusterCollectionMuon, hLargestHitClusterCollectionElectron, hLargestHitClusterCollectionScatter, hLargestHitClusterCollectionChExch, hLargestHitClusterCollectionOther},
        {hNumClustersCollectionAbs, hNumClustersCollectionMuon, hNumClustersCollectionElectron, hNumClustersCollectionScatter, hNumClustersCollectionChExch, hNumClustersCollectionOther},
        {hUnRecoHitsCollectionAbs, hUnRecoHitsCollectionMuon, hUnRecoHitsCollectionElectron, hUnRecoHitsCollectionScatter, hUnRecoHitsCollectionChExch, hUnRecoHitsCollectionOther},

        {hNegativeTimePrimaryHits},
        {hTimePrimaryHits},

        // Through-going tracks
        {hNumTGTracksAbs, hNumTGTracksMuon, hNumTGTracksElectron, hNumTGTracksScatter, hNumTGTracksChExch, hNumTGTracksOther},

    };

    std::vector<std::vector<TString>> PlotLabelGroups = {
        // Incident KE
        {"All", "Pions", "Muons", "Electrons"},
        {"All", "Pions", "Muons", "Electrons"},
        {"Corrected", "True"},
        {"All", "Pions", "Muons", "Electrons"},
        {"All", "Pions", "Muons", "Electrons"},
        {"All", "Pions", "Muons", "Electrons"},
        {"All", "Pions", "Muons", "Electrons"},
        {"All", "Pions", "Muons", "Electrons"},

        // Interacting KE
        {"All", "True", "Scatter", "Ch. exch.", "Muon", "Electron", "Other"},
        {"All", "True", "Scatter", "Ch. exch.", "Muon TG", "Muon decay", "Muon CAR", "Electron", "Other"},
        {"All", "True", "Abs", "Ch. exch.", "Muon", "Electron", "Other"},
        {"All", "True", "Abs", "Ch. exch.", "Muon TG", "Muon decay", "Muon CAR", "Electron", "Other"},

        // True interacting KE
        {"All"},
        {"Abs", "Scatter", "Ch. exch.", "Other"},

        // True events classified breakdown
        {"Abs", "Scatter", "Rejected"},
        {"Abs", "Scatter", "Rejected"},
        {"Abs", "Scatter", "Rejected"},

        // Rejected events
        {"Data-prods", "Shower-like", "Red. vol.", "PID reject", "> 1 pion", "Hit clusters"},
        {"Data-prods", "Shower-like", "Red. vol.", "PID reject", "> 1 pion", "Hit clusters"},

        // Cross-sections (unfolded)
        {"True", "Unf."},
        {"True", "Unf."},

        // Total true-cross section
        {"Abs", "Scatter", "Ch. exch.", "Other"},

        // Total true-cross section with no thin-slice approximation
        {"Thin-slice", "Log."},
        {"Thin-slice", "Log."},

        // Unreconstructed hits
        {"Abs", "Muon", "Electron", "Scatter", "Ch. exch.", "Other"},
        {"Abs", "Muon", "Electron", "Scatter", "Ch. exch.", "Other"},
        {"Abs", "Muon", "Electron", "Scatter", "Ch. exch.", "Other"},
        {"Abs", "Muon", "Electron", "Scatter", "Ch. exch.", "Other"},
        {"Abs", "Muon", "Electron", "Scatter", "Ch. exch.", "Other"},
        {"Abs", "Muon", "Electron", "Scatter", "Ch. exch.", "Other"},
        {"Abs", "Muon", "Electron", "Scatter", "Ch. exch.", "Other"},
        {"Abs", "Muon", "Electron", "Scatter", "Ch. exch.", "Other"},
        {"Abs", "Muon", "Electron", "Scatter", "Ch. exch.", "Other"},
        {"Abs", "Muon", "Electron", "Scatter", "Ch. exch.", "Other"},
        {"Negative time hits"},
        {"All time hits"},

        // Through-going tracks
        {"Abs", "Muon", "Electron", "Scatter", "Ch. exch.", "Other"},

    };

    std::vector<TString> PlotTitles = {
        // Incident KE
        "Incident/IncidentKE",
        "Incident/IncidentKEFine",
        "Incident/IncidentKECorrected",
        "Incident/FrontFaceKE",
        "Incident/FrontFaceKENoWeight",
        "Incident/FrontFaceKENoWeightPre",
        "Incident/FrontFacePionCutKE",
        "Incident/WireChamberKE",

        // Interacting KE
        "RecoInteracting/AbsInteractingKE",
        "RecoInteracting/AbsInteractingKEDetailed",
        "RecoInteracting/ScatterInteractingKE",
        "RecoInteracting/ScatterInteractingKEDetailed",

        // True interacting KE
        "TrueInteracting/AllInteracting",
        "TrueInteracting/AllInteractingBreakdown",

        // True events classified breakdown
        "TrueInteracting/TrueAbsBreakdown",
        "TrueInteracting/TrueScatterBreakdown",
        "TrueInteracting/TrueChExchBreakdown",

        // Rejected events
        "Rejected/TrueAbsRej",
        "Rejected/TrueScatterRej",

        // True cross-section (unfolded)
        "CrossSection/PionAbsCrossSection",
        "CrossSection/PionScatterCrossSection",

        // Total true-cross section
        "CrossSection/TotalTrueCrossSection",

        // Total true-cross section with no thin-slice approximation
        "CrossSection/AbsLogComp",
        "CrossSection/ScatterLogComp",

        // Unreconstructed hits
        "Hits/ClusterSizesInduction",
        "Hits/NumLargeClustersInduction",
        "Hits/LargestClusterInduction",
        "Hits/NumClustersInduction",
        "Hits/UnReconstructedInduction",
        "Hits/ClusterSizesCollection",
        "Hits/NumLargeClustersCollection",
        "Hits/LargestClusterCollection",
        "Hits/NumClustersCollection",
        "Hits/UnReconstructedCollection",
        "Hits/PrimaryNegativeTimeHits",
        "Hits/PrimaryTimeHits",

        // Through-going tracks
        "TGTracks/NumTGTracks",

    };

    std::vector<TString> XLabels = {
        // Incident KE
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",

        // Interacting KE
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",

        // True interacting KE
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",

        // True events classified breakdown
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",

        // Rejected events
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",

        // Cross-sections (unfolded)
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",

        // Total true-cross section
        "Kinetic energy [MeV]",

        // Total true-cross section with no thin-slice approximation
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",

        // Unreconstructed hits
        "Cluster size [cm]",
        "Number of large clusters",
        "Cluster size [cm]",
        "# of clusters",
        "Unreconstructed hits",
        "Cluster size [cm]",
        "Number of large clusters",
        "Cluster size [cm]",
        "# of clusters",
        "Unreconstructed hits",
        "Hit time [us]",
        "Hit time [us]",

        // Through-going tracks
        "# of TG tracks",

    };

    std::vector<TString> YLabels = {
        // Incident KE
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",

        // Interacting KE
        "Counts",
        "Counts",
        "Counts",
        "Counts",

        // True interacting KE
        "Counts",
        "Counts",

        // True events classified breakdown
        "Counts",
        "Counts",
        "Counts",

        // Rejected events
        "Counts",
        "Counts",

        // Cross-sections (unfolded)
        "Cross section [barn] per 100 MeV",
        "Cross section [barn] per 100 MeV",

        // Total true-cross section
        "Cross section [barn] per 100 MeV",

        // Total true-cross section with no thin-slice approximation
        "Cross section [barn] per 100 MeV",
        "Cross section [barn] per 100 MeV",

        // Unreconstructed hits
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",

        // Through-going tracks
        "Counts",

    };

    std::vector<bool> PlotStacked = {
        // Incident KE
        true,
        true,
        false,
        true,
        true,
        true,
        true,
        true,

        // Interacting KE
        true,
        true,
        true,
        true,

        // True interacting KE
        true,
        true,

        // True events classified breakdown
        true,
        true,
        true,

        // Rejected events
        true,
        true,

        // Cross-sections (unfolded)
        false,
        false,

        // Total true-cross section
        true,

        // Total true-cross section with no thin-slice approximation
        false,
        false,

        // Unreconstructed hits
        true,
        true,
        true,
        true,
        true,
        true,
        true,
        true,
        true,
        true,
        true,
        true,

        // Through-going tracks
        true,

    };

    std::vector<std::vector<bool>> PlotsAsPoints = {
        // Incident KE
        {true, false, false, false},
        {true, false, false, false},
        {true, false},
        {true, false, false, false},
        {true, false, false, false},
        {true, false, false, false},
        {true, false, false, false},
        {true, false, false, false},

        // Interacting KE
        {true, false, false, false, false, false, false},
        {true, false, false, false, false, false, false, false, false},
        {true, false, false, false, false, false, false},
        {true, false, false, false, false, false, false, false, false},

        // True interacting KE
        {false},
        {false, false, false, false},

        // True events classified breakdown
        {false, false, false},
        {false, false, false},
        {false, false, false},

        // Rejected events
        {false, false, false, false, false, false},
        {false, false, false, false, false, false},

        // Cross-sections (unfolded)
        {false, true},
        {false, true},

        // Total true-cross section
        {false, false, false, false},

        // Total true-cross section with no thin-slice approximation
        {false, false},
        {false, false},

        // Unreconstructed hits
        {false, false, false, false, false, false},
        {false, false, false, false, false, false},
        {false, false, false, false, false, false},
        {false, false, false, false, false, false},
        {false, false, false, false, false, false},
        {false, false, false, false, false, false},
        {false, false, false, false, false, false},
        {false, false, false, false, false, false},
        {false, false, false, false, false, false},
        {false, false, false, false, false, false},
        {false},
        {false},

        // Through-going tracks
        {false, false, false, false, false, false},

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

    ////////////////////////////////////////////////
    // Get universes for statistical fluctuations //
    ////////////////////////////////////////////////

    TH1D* hMeasure1 = new TH1D(Form("hMeasure_Univ%d", 1), Form("Measured Vector for Universe %d", 1), NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    TH1D* hMeasure2 = new TH1D(Form("hMeasure_Univ%d", 2), Form("Measured Vector for Universe %d", 2), NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    for (int i = 1; i <= NUM_SIGNAL_TYPES * NUM_BINS_KE; ++i) {
        double xsec  = hMeasureVectorNominal->GetBinContent(i);
        double sigma = std::sqrt(hCovariance->GetBinContent(i, i));
        hMeasure1->SetBinContent(i, xsec + sigma);
        hMeasure2->SetBinContent(i, xsec - sigma);
    }
    std::vector<TH1D*> MeasureVectorUnivs = {hMeasure1, hMeasure2};

    DrawHistosWithErrorBands(
        hMeasureVectorNominal,
        MeasureVectorUnivs,
        SaveDir,
        "Covariance",
        "Measure",
        TextSize,
        FontStyle,
        "Measured Cross Section",
        "Bin number",
        "Cross Section [barn]"
    );

    ///////////////////////////
    // Two-dimensional plots //
    ///////////////////////////

    // Get extra matrices
    TH2D* hFracCovMatrix = new TH2D(
        "Fractional Covariance", "Fractional Covariance;Reco (j, #beta);True (i, #alpha)",
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE,
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    );
    TH2D* hCorrMatrix = new TH2D(
        "Correlation", "Correlation;Reco (j, #beta);True (i, #alpha)",
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE,
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    );
    GetFracCovAndCorrMatrix(hMeasureVectorNominal, hCovariance, hFracCovMatrix, hCorrMatrix);

    TH2D* hUnfFracCovMatrix = new TH2D(
        "Unfolded Fractional Covariance", "Unfolded Fractional Covariance;Reco (j, #beta);True (i, #alpha)",
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE,
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    );
    TH2D* hUnfCorrMatrix = new TH2D(
        "Unfolded Correlation", "Unfolded Correlation;Reco (j, #beta);True (i, #alpha)",
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE,
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    );
    GetFracCovAndCorrMatrix(hUnfReco, hUnfCovariance, hUnfFracCovMatrix, hUnfCorrMatrix);

    std::vector<TH2*> TwoDPlots = {
        hResponseMatrix,
        hResponseInvMatrix,
        hCovariance,
        hFracCovMatrix,
        hCorrMatrix,
        hUnfCovariance,
        hUnfFracCovMatrix,
        hUnfCorrMatrix
    };

    std::vector<TString> TwoDTitles = {
        // Response
        "Response/ResponseMatrix",
        "Response/ResponseInverseMatrix",

        // Covariances
        "Covariance/Covariance",
        "Covariance/FracCovariance",
        "Covariance/Correlation",
        "Covariance/UnfoldedCovariance",
        "Covariance/UnfoldedFracCovariance",
        "Covariance/UnfoldedCorrelation"
    };

    std::vector<std::pair<double,double>> TwoDRanges = {
        // Response
        {0, 0},
        {0, 0},

        // Covariances
        {0, 0},
        {0, 0},
        {-1, -1},
        {0, 0},
        {0, 0},
        {-1, -1}
    };

    std::vector<bool> TwoDDisplayNumbers = {
        // Response
        false,
        false,

        // Covariances
        false,
        false,
        false,
        false,
        false,
        false
    };

    printTwoDPlots(SaveDir, TwoDPlots, TwoDTitles, TwoDRanges, TwoDDisplayNumbers);

    ////////////////////////////////////
    // Save everything from later use //
    ////////////////////////////////////

    TString outPath = "/exp/lariat/app/users/epelaez/analysis_abs_scatt/histos/nominal/RecoClassify3Cat_AllHists.root";
    TFile outAll(outPath, "RECREATE");
    outAll.cd();

    std::unordered_set<std::string> written;

    // Helper lambda: write once by name
    auto writeOnce = [&](TObject* obj) {
        if (!obj) return;
        const std::string name = obj->GetName();
        if (name.empty()) return;
        if (written.insert(name).second) {
            obj->Write(name.c_str(), TObject::kOverwrite);
        }
    };

    // 1D groups
    for (auto& group : PlotGroups) {
        for (auto* h : group) writeOnce(h);
    }

    // 2D plots
    for (auto* h2 : TwoDPlots) writeOnce(h2);

    // Other important hist collections
    for (auto* h : UnfoldedRecoHistos) writeOnce(h);

    writeOnce(hUnfReco);
    writeOnce(hResponseMatrix);
    writeOnce(hResponseInvMatrix);
    writeOnce(hCovariance);
    writeOnce(hUnfCovariance);
    writeOnce(hFracCovMatrix);
    writeOnce(hCorrMatrix);
    writeOnce(hUnfFracCovMatrix);
    writeOnce(hUnfCorrMatrix);

    outAll.Write();
    outAll.Close();

    std::cout << "\nWrote " << written.size() << " objects to " << outPath << std::endl;
}