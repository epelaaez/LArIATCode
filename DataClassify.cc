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

#include <vector>

#include "Helpers.cc"

void DataClassify() {
    // Set defaults
    gStyle->SetOptStat(0); // get rid of stats box
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
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
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis/figs/DataClassify/";

    // Load files
    TChain* Chain = new TChain("anatree/anatree");
    Chain->Add("/exp/lariat/app/users/epelaez/files/anatree_60a_data/chunks/*.root");
    std::cout << "Files:   " << Chain->GetListOfFiles()->GetEntries() << std::endl;

    // Load nominal MC histos
    TString NominalHistsPath = "/exp/lariat/app/users/epelaez/histos/nominal/RecoClassify3Cat_AllHists.root";
    std::unique_ptr<TFile> fNom(TFile::Open(NominalHistsPath, "READ"));

    // Load nominal MC histos (abs + scatt)
    TString NominalAbsScattHistsPath = "/exp/lariat/app/users/epelaez/analysis_abs_scatt/histos/nominal/RecoClassify3Cat_AllHists.root";
    std::unique_ptr<TFile> fNomAbsScatt(TFile::Open(NominalAbsScattHistsPath, "READ"));
    
    ///////////////////////
    // Counters for cuts //
    ///////////////////////

    int TotalEvents = 0;
    int EventsPassingProj = 0;
    int EventsPassingAperture = 0;
    int EventsPassingTOFMass = 0;
    int EventsWithWCMatch = 0;
    int EventsWithPickyWC = 0;
    int EventsPassingTG = 0;
    int EventsPassingSmall = 0;
    int EventsInRedVol = 0;
    int EventsPassingPrimaryPID = 0;
    int EventsSelectedAsScatter = 0;
    int EventsSelectedAsAbsNp = 0;
    int EventsNoSecondary = 0;
    int EventsSelectedAsAbs0p = 0;

    //////////////////////////
    // Get histos for plots //
    //////////////////////////

    TH1D* hMCIncidentKEPion     = dynamic_cast<TH1D*>(fNom->Get("hIncidentKEPion"));
    TH1D* hMCIncidentKEElectron = dynamic_cast<TH1D*>(fNom->Get("hIncidentKEElectron"));
    TH1D* hMCIncidentKEMuon     = dynamic_cast<TH1D*>(fNom->Get("hIncidentKEMuon"));

    TH1D* hMCIncidentKEPionFine     = dynamic_cast<TH1D*>(fNom->Get("hIncidentKEPionFine"));
    TH1D* hMCIncidentKEElectronFine = dynamic_cast<TH1D*>(fNom->Get("hIncidentKEElectronFine"));
    TH1D* hMCIncidentKEMuonFine     = dynamic_cast<TH1D*>(fNom->Get("hIncidentKEMuonFine"));

    TH1D* hMCSmallTrksInCylinderPions     = dynamic_cast<TH1D*>(fNom->Get("hSmallTrksInCylinderPions"));
    TH1D* hMCSmallTrksInCylinderElectrons = dynamic_cast<TH1D*>(fNom->Get("hSmallTrksInCylinderElectrons"));
    TH1D* hMCSmallTrksInCylinderMuons     = dynamic_cast<TH1D*>(fNom->Get("hSmallTrksInCylinderMuons"));

    TH1D* hMCFrontFaceKEPion     = dynamic_cast<TH1D*>(fNom->Get("hFrontFaceKEPion"));
    TH1D* hMCFrontFaceKEElectron = dynamic_cast<TH1D*>(fNom->Get("hFrontFaceKEElectron"));
    TH1D* hMCFrontFaceKEMuon     = dynamic_cast<TH1D*>(fNom->Get("hFrontFaceKEMuon"));

    TH1D* hMCWCKEPion     = dynamic_cast<TH1D*>(fNom->Get("hWCKEPion"));
    TH1D* hMCWCKEMuon     = dynamic_cast<TH1D*>(fNom->Get("hWCKEMuon"));
    TH1D* hMCWCKEElectron = dynamic_cast<TH1D*>(fNom->Get("hWCKEElectron"));

    // Abs 0p interacting KE
    TH1D* hMCPionAbs0pKETrue     = dynamic_cast<TH1D*>(fNom->Get("hPionAbs0pKETrue"));
    TH1D* hMCPionAbs0pKEAbsNp    = dynamic_cast<TH1D*>(fNom->Get("hPionAbs0pKEAbsNp"));
    TH1D* hMCPionAbs0pKEScatter  = dynamic_cast<TH1D*>(fNom->Get("hPionAbs0pKEScatter"));
    TH1D* hMCPionAbs0pKEChExch   = dynamic_cast<TH1D*>(fNom->Get("hPionAbs0pKEChExch"));
    TH1D* hMCPionAbs0pKEMuon     = dynamic_cast<TH1D*>(fNom->Get("hPionAbs0pKEMuon"));
    TH1D* hMCPionAbs0pKEElectron = dynamic_cast<TH1D*>(fNom->Get("hPionAbs0pKEElectron"));
    TH1D* hMCPionAbs0pKEOther    = dynamic_cast<TH1D*>(fNom->Get("hPionAbs0pKEOther"));

    // Abs Np interacting KE
    TH1D* hMCPionAbsNpKETrue     = dynamic_cast<TH1D*>(fNom->Get("hPionAbsNpKETrue"));
    TH1D* hMCPionAbsNpKEAbs0p    = dynamic_cast<TH1D*>(fNom->Get("hPionAbsNpKEAbs0p"));
    TH1D* hMCPionAbsNpKEScatter  = dynamic_cast<TH1D*>(fNom->Get("hPionAbsNpKEScatter"));
    TH1D* hMCPionAbsNpKEChExch   = dynamic_cast<TH1D*>(fNom->Get("hPionAbsNpKEChExch"));
    TH1D* hMCPionAbsNpKEMuon     = dynamic_cast<TH1D*>(fNom->Get("hPionAbsNpKEMuon"));
    TH1D* hMCPionAbsNpKEElectron = dynamic_cast<TH1D*>(fNom->Get("hPionAbsNpKEElectron"));
    TH1D* hMCPionAbsNpKEOther    = dynamic_cast<TH1D*>(fNom->Get("hPionAbsNpKEOther"));

    // Scatter interacting KE
    TH1D* hMCPionScatterKETrue     = dynamic_cast<TH1D*>(fNom->Get("hPionScatterKETrue"));
    TH1D* hMCPionScatterKEAbs0p    = dynamic_cast<TH1D*>(fNom->Get("hPionScatterKEAbs0p"));
    TH1D* hMCPionScatterKEAbsNp    = dynamic_cast<TH1D*>(fNom->Get("hPionScatterKEAbsNp"));
    TH1D* hMCPionScatterKEChExch   = dynamic_cast<TH1D*>(fNom->Get("hPionScatterKEChExch"));
    TH1D* hMCPionScatterKEMuon     = dynamic_cast<TH1D*>(fNom->Get("hPionScatterKEMuon"));
    TH1D* hMCPionScatterKEElectron = dynamic_cast<TH1D*>(fNom->Get("hPionScatterKEElectron"));
    TH1D* hMCPionScatterKEOther    = dynamic_cast<TH1D*>(fNom->Get("hPionScatterKEOther"));

    // Abs interacting KE (abs + scatt)
    TH1D* hMCPionAbsKETrue     = dynamic_cast<TH1D*>(fNomAbsScatt->Get("hPionAbsKETrue"));
    TH1D* hMCPionAbsKEScatter  = dynamic_cast<TH1D*>(fNomAbsScatt->Get("hPionAbsKEScatter"));
    TH1D* hMCPionAbsKEChExch   = dynamic_cast<TH1D*>(fNomAbsScatt->Get("hPionAbsKEChExch"));
    TH1D* hMCPionAbsKEMuon     = dynamic_cast<TH1D*>(fNomAbsScatt->Get("hPionAbsKEMuon"));
    TH1D* hMCPionAbsKEElectron = dynamic_cast<TH1D*>(fNomAbsScatt->Get("hPionAbsKEElectron"));
    TH1D* hMCPionAbsKEOther    = dynamic_cast<TH1D*>(fNomAbsScatt->Get("hPionAbsKEOther"));

    // Abs interacting KE
    TH1D* hMCPionScatterKETrue2     = dynamic_cast<TH1D*>(fNomAbsScatt->Get("hPionScatterKETrue"));
    TH1D* hMCPionScatterKEAbs2      = dynamic_cast<TH1D*>(fNomAbsScatt->Get("hPionScatterKEAbs"));
    TH1D* hMCPionScatterKEChExch2   = dynamic_cast<TH1D*>(fNomAbsScatt->Get("hPionScatterKEChExch"));
    TH1D* hMCPionScatterKEMuon2     = dynamic_cast<TH1D*>(fNomAbsScatt->Get("hPionScatterKEMuon"));
    TH1D* hMCPionScatterKEElectron2 = dynamic_cast<TH1D*>(fNomAbsScatt->Get("hPionScatterKEElectron"));
    TH1D* hMCPionScatterKEOther2    = dynamic_cast<TH1D*>(fNomAbsScatt->Get("hPionScatterKEOther"));

    ///////////////////
    // Load branches //
    ///////////////////

    int run, subrun, event; bool isData = true;
    Chain->SetBranchAddress("run", &run);
    Chain->SetBranchAddress("subrun", &subrun);
    Chain->SetBranchAddress("event", &event);

    // Track information
    int   ntracks_reco;                       Chain->SetBranchAddress("ntracks_reco",    &ntracks_reco);
    static float trkvtxx[kMaxTrackData];                 Chain->SetBranchAddress("trkvtxx",          &trkvtxx);
    static float trkvtxy[kMaxTrackData];                 Chain->SetBranchAddress("trkvtxy",          &trkvtxy);
    static float trkvtxz[kMaxTrackData];                 Chain->SetBranchAddress("trkvtxz",          &trkvtxz);
    static float trkendx[kMaxTrackData];                 Chain->SetBranchAddress("trkendx",          &trkendx);
    static float trkendy[kMaxTrackData];                 Chain->SetBranchAddress("trkendy",          &trkendy);
    static float trkendz[kMaxTrackData];                 Chain->SetBranchAddress("trkendz",          &trkendz);
    static int   trkWCtoTPCMatch[kMaxTrackData];         Chain->SetBranchAddress("trkWCtoTPCMatch",  &trkWCtoTPCMatch);

    // Wire-chamber track information
    float beamline_mass;                        Chain->SetBranchAddress("beamline_mass",     &beamline_mass);
    int   nwctrks;                              Chain->SetBranchAddress("nwctrks",           &nwctrks);
    static float wctrk_momentum[kMaxWCTracksData];         Chain->SetBranchAddress("wctrk_momentum",    &wctrk_momentum);
    static float wctrk_theta[kMaxWCTracksData];            Chain->SetBranchAddress("wctrk_theta",        &wctrk_theta);
    static float wctrk_phi[kMaxWCTracksData];              Chain->SetBranchAddress("wctrk_phi",          &wctrk_phi);
    static int   wctrk_picky[kMaxWCTracksData];            Chain->SetBranchAddress("wctrk_picky",        &wctrk_picky);
    static float WC1xPos[kMaxWCTracksData];                 Chain->SetBranchAddress("WC1xPos",            &WC1xPos);
    static float WC1yPos[kMaxWCTracksData];                 Chain->SetBranchAddress("WC1yPos",            &WC1yPos);
    static float WC1zPos[kMaxWCTracksData];                 Chain->SetBranchAddress("WC1zPos",            &WC1zPos);
    static float WC2xPos[kMaxWCTracksData];                 Chain->SetBranchAddress("WC2xPos",            &WC2xPos);
    static float WC2yPos[kMaxWCTracksData];                 Chain->SetBranchAddress("WC2yPos",            &WC2yPos);
    static float WC2zPos[kMaxWCTracksData];                 Chain->SetBranchAddress("WC2zPos",            &WC2zPos);
    static float WC3xPos[kMaxWCTracksData];                 Chain->SetBranchAddress("WC3xPos",            &WC3xPos);
    static float WC3yPos[kMaxWCTracksData];                 Chain->SetBranchAddress("WC3yPos",            &WC3yPos);
    static float WC3zPos[kMaxWCTracksData];                 Chain->SetBranchAddress("WC3zPos",            &WC3zPos);
    static float WC4xPos[kMaxWCTracksData];                 Chain->SetBranchAddress("WC4xPos",            &WC4xPos);
    static float WC4yPos[kMaxWCTracksData];                 Chain->SetBranchAddress("WC4yPos",            &WC4yPos);
    static float WC4zPos[kMaxWCTracksData];                 Chain->SetBranchAddress("WC4zPos",            &WC4zPos);

    // Calorimetry information
    static int   ntrkcalopts[kMaxTrackData][2];                    Chain->SetBranchAddress("ntrkcalopts", &ntrkcalopts);
    static float trkdedx[kMaxTrackData][2][kMaxTrackHitsData];         Chain->SetBranchAddress("trkdedx",     &trkdedx);
    static float trkrr[kMaxTrackData][2][kMaxTrackHitsData];           Chain->SetBranchAddress("trkrr",       &trkrr);
    static float trkpitch[kMaxTrackData][2][kMaxTrackHitsData];        Chain->SetBranchAddress("trkpitch",    &trkpitch);
    static float trkxyz[kMaxTrackData][2][kMaxTrackHitsData][3];       Chain->SetBranchAddress("trkxyz",      &trkxyz);

    // Trajectory information for tracks
    static int   nTrajPoint[kMaxTrackData];                  Chain->SetBranchAddress("nTrajPoint", &nTrajPoint);
    static float trjPt_X[kMaxTrackData][kMaxTrajHitsData];       Chain->SetBranchAddress("trjPt_X",    &trjPt_X);
    static float trjPt_Y[kMaxTrackData][kMaxTrajHitsData];       Chain->SetBranchAddress("trjPt_Y",    &trjPt_Y);
    static float trjPt_Z[kMaxTrackData][kMaxTrajHitsData];       Chain->SetBranchAddress("trjPt_Z",    &trjPt_Z);

    // Information about wire plane hits
    int    nhits;                              Chain->SetBranchAddress("nhits",              &nhits);
    static int    hit_plane[kMaxHitsData];         Chain->SetBranchAddress("hit_plane",           hit_plane);
    static int    hit_channel[kMaxHitsData];       Chain->SetBranchAddress("hit_channel",         hit_channel);
    static int    hit_trkid[kMaxHitsData];         Chain->SetBranchAddress("hit_trkid",           hit_trkid);
    static float  hit_driftT[kMaxHitsData];        Chain->SetBranchAddress("hit_driftT",           hit_driftT);
    static float  hit_x[kMaxHitsData];             Chain->SetBranchAddress("hit_x",                hit_x);
    static float  hit_y[kMaxHitsData];             Chain->SetBranchAddress("hit_y",                hit_y);
    static float  hit_z[kMaxHitsData];             Chain->SetBranchAddress("hit_z",                hit_z);

    // TOF object
    int ntof;                     Chain->SetBranchAddress("ntof",          &ntof);
    static float tof[kMaxTOFData];           Chain->SetBranchAddress("tof",           tof);
    static float tof_timestamp[kMaxTOFData]; Chain->SetBranchAddress("tof_timestamp", tof_timestamp);

    //////////////////////
    // Create histogram //
    //////////////////////

    // Small tracks in cylinder
    TH1D* hSmallTracksInCylinder = new TH1D("hSmallTracksInCylinder ", "hSmallTracksInCylinder ;;", 10, 0, 10);

    // Selected interaction type
    TH1D* hPionAbs0pKE   = new TH1D("hPionAbs0pKE", "hPionAbs0pKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbsNpKE   = new TH1D("hPionAbsNpKE", "hPionAbsNpKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionScatterKE = new TH1D("hPionScatterKE", "hPionScatterKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    std::vector<TH1*> RecoSignals = {
        hPionAbs0pKE, hPionAbsNpKE, hPionScatterKE
    };

    // Also want to keep track of total abs
    TH1D* hPionAbsKE = new TH1D("hPionAbsKE", "hPionAbsKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    std::vector<TH1*> RecoSignalsAbsScatt = {
        hPionAbsKE, hPionScatterKE
    };

    // Incident kinetic energy
    TH1D* hIncidentKE     = new TH1D("hIncidentKE", "hIncidentKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hIncidentKEFine = new TH1D("hIncidentKEFine", "hIncidentKEFine;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());

    // Kinetic energy at the front face of the TPC
    TH1D* hFrontFaceKE = new TH1D("hFrontFaceKE", "hFrontFaceKE;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());

    // Wire chamber kinetic energy (before energy loss correction)
    TH1D* hWCKE = new TH1D("hWCKE", "hWCKE;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());

    // Wire-chamber quality
    TH1D* hWCHits          = new TH1D("hWCHits", "hWCHits", 10, 0, 10);
    TH1D* hRadDistMidPlane = new TH1D("hRadDistMidPlane", "hRadDistMidPlane", 60, 0.0, 20.);
    TH1D* hRadDistWC4      = new TH1D("hRadDistWC4", "RadDistWC4", 90, 0.0, 30.);

    //////////////////////
    // Loop over events //
    //////////////////////

    bool verbose = false;

    Long64_t i = 0;
    while (true) {
        if (SKIP_INDICES_DATA.count(i)) { i++; continue; }
        if (Chain->GetEntry(i++) <= 0) break;

        // Reset variables
        if (verbose) std::cout << std::endl;
        if (verbose) std::cout << "=================================" << std::endl;
        if (verbose) std::cout << "Got tree entry: " << i << std::endl;
        EventVariablesData ev;
        if (verbose) std::cout << "Variables reset" << std::endl;
        if (verbose) std::cout << "=================================" << std::endl;

        // Make script go faster
        // if (i > USE_NUM_EVENTS) break;

        TotalEvents++;

        ////////////////////////////////
        // Load variables of interest //
        ////////////////////////////////

        // Load track information
        if (verbose) std::cout << std::endl;
        if (verbose) std::cout << "=================================" << std::endl;
        if (verbose) std::cout << "Event: " << event << std::endl;
        
        // First, we just want to grab WC to TPC match
        int primaryTrackIdx = -1;
        for (size_t trk_idx = 0; trk_idx < std::min(ntracks_reco, kMaxTrackData); ++trk_idx) {
            if (trkWCtoTPCMatch[trk_idx]) {
                // Grab position information
                ev.WC2TPCPrimaryBeginX = trkvtxx[trk_idx];
                ev.WC2TPCPrimaryBeginY = trkvtxy[trk_idx];
                ev.WC2TPCPrimaryBeginZ = trkvtxz[trk_idx];
                ev.WC2TPCPrimaryEndX   = trkendx[trk_idx];
                ev.WC2TPCPrimaryEndY   = trkendy[trk_idx];
                ev.WC2TPCPrimaryEndZ   = trkendz[trk_idx];

                // Grab calorimetry information (in collection plane)
                int npts_dedx = std::min(ntrkcalopts[trk_idx][1], kMaxTrackHitsData);
                ev.wcMatchResR.assign(trkrr[trk_idx][1], trkrr[trk_idx][1] + npts_dedx);
                ev.wcMatchDEDX.assign(trkdedx[trk_idx][1], trkdedx[trk_idx][1] + npts_dedx);

                for (size_t dep_idx = 0; dep_idx < std::min(ntrkcalopts[trk_idx][1], kMaxTrackHitsData); ++dep_idx) {
                    ev.wcMatchEDep.push_back(trkdedx[trk_idx][1][dep_idx] * trkpitch[trk_idx][1][dep_idx]);
                    ev.wcMatchXPos.push_back(trkxyz[trk_idx][1][dep_idx][0]);
                    ev.wcMatchYPos.push_back(trkxyz[trk_idx][1][dep_idx][1]);
                    ev.wcMatchZPos.push_back(trkxyz[trk_idx][1][dep_idx][2]);
                }

                // Get location information
                int npts_wc2tpc = std::min(nTrajPoint[trk_idx], kMaxTrajHitsData);
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
        if (verbose) std::cout << "Number of matches: " << ev.WC2TPCsize  << std::endl;

        // Set beamline information
        ev.TOFMass   = std::abs(beamline_mass);
        if (ntof == 1) ev.tofObject = tof[0];

        if (verbose) std::cout << "Beamline mass: " << ev.TOFMass  << std::endl;
        if (verbose) std::cout << "TOF: " << ev.tofObject  << std::endl;

        // Copy vertex and end for all tracks
        int npts_trk = std::min(ntracks_reco, kMaxTrackData);
        ev.recoEndX.assign(trkendx,  trkendx  + npts_trk);
        ev.recoEndY.assign(trkendy,  trkendy  + npts_trk);
        ev.recoEndZ.assign(trkendz,  trkendz  + npts_trk);
        ev.recoBeginX.assign(trkvtxx, trkvtxx + npts_trk);
        ev.recoBeginY.assign(trkvtxy, trkvtxy + npts_trk);
        ev.recoBeginZ.assign(trkvtxz, trkvtxz + npts_trk);

        // Now, we want to loop through all tracks
        for (size_t trk_idx = 0; trk_idx < std::min(ntracks_reco, kMaxTrackData); ++trk_idx) {
            // Grab calorimetry information
            ev.recoResR.push_back(std::vector<double>(trkrr[trk_idx][1], trkrr[trk_idx][1] + std::min(ntrkcalopts[trk_idx][1], kMaxTrackHitsData)));
            ev.recoDEDX.push_back(std::vector<double>(trkdedx[trk_idx][1], trkdedx[trk_idx][1] + std::min(ntrkcalopts[trk_idx][1], kMaxTrackHitsData)));

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
        if (verbose) std::cout << "Number of WC tracks : " << nwctrks << std::endl;

        if (nwctrks == 1) {
            ev.wcTrackPicky = wctrk_picky[0];

            ev.WCTrackMomentum = wctrk_momentum[0];
            ev.WCTheta         = wctrk_theta[0];
            ev.WCPhi           = wctrk_phi[0];
            ev.WC4PrimaryX     = WC4xPos[0];

            ev.wcHit0 = {WC1xPos[0], WC1yPos[0], WC1zPos[0]};
            ev.wcHit1 = {WC2xPos[0], WC2yPos[0], WC2zPos[0]};
            ev.wcHit2 = {WC3xPos[0], WC3yPos[0], WC3zPos[0]};
            ev.wcHit3 = {WC4xPos[0], WC4yPos[0], WC4zPos[0]};
        }
        if (verbose) std::cout << "Is WC track picky? " << ev.wcTrackPicky  << std::endl;

        // Get information about wire hits
        std::map<int, std::vector<int>>    trackHitMap;
        std::map<int, std::vector<double>> trackHitXMap;
        std::map<int, std::vector<double>> trackHitYMap;
        std::map<int, std::vector<double>> trackHitZMap;

        for (size_t i_hit = 0; i_hit < std::min(nhits, kMaxHitsData); ++i_hit) {
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

        /////////////////////////////////////////
        // Wire-chamber and other initial cuts //
        /////////////////////////////////////////

        // Project downstream to midplane @ -437.97 (without angular corrections)
        std::vector<double> midUp = projToZ(ev.wcHit0, ev.wcHit1, -437.97);
        // Use this point and WC3 to project up to WC4
        std::vector<double> projDown = projToZ(midUp, ev.wcHit2, -95.0);
        // Requires some corrections because magnets are not the same
        projDown[0] -= tan(1.32 * TMath::Pi() / 180.0) * (-95.0 - -437.97);

        // Compare x and y coordinate in projection and real hit for WC4
        double radDistWC4 = TMath::Sqrt(pow(projDown[0] - ev.wcHit3.at(0), 2.) + pow(projDown[1] - ev.wcHit3.at(1), 2.));

        // Project upstream to midplane @ -437.97
        std::vector<double> midDown = projToZ(ev.wcHit2, ev.wcHit3, -437.97);
        midDown[0] -= tan(1.32 * TMath::Pi() / 180.0) * (-339.57 - -437.97);
        double midPlaneDist = TMath::Sqrt(pow(midUp[0] - midDown[0] + 0.75, 2) + pow(midUp[1] - midDown[1], 2));
        // double midPlaneDist = TMath::Sqrt(pow(midUp[0] - midDown[0], 2) + pow(midUp[1] - midDown[1], 2));

        hRadDistMidPlane->Fill(midPlaneDist);
        hRadDistWC4->Fill(radDistWC4);

        // Cuts
        if (radDistWC4 > 8.0) continue;
        if (midPlaneDist > 3.0) continue;
        EventsPassingProj++;

        // Check projected tracks go through all apertures
        bool Magnet1ApertureCheck = CheckUpstreamMagnetAperture(ev.wcHit0, ev.wcHit1);
        bool Magnet2ApertureCheck = CheckDownstreamMagnetAperture(ev.wcHit2, ev.wcHit3);
        bool DSColApertureCheck   = CheckDownstreamCollimatorAperture(ev.wcHit2, ev.wcHit3);

        if (!Magnet1ApertureCheck || !Magnet2ApertureCheck || !DSColApertureCheck) continue;
        EventsPassingAperture++;

        // Candidate mass cut, keep pions, muons and electrons
        if (std::abs(ev.TOFMass) > PI_MU_EL_MASS_CUTOFF) continue;
        EventsPassingTOFMass++;

        // If no or multiple tracks matched to wire-chamber, skip
        if (!(ev.WC2TPCMatch && ev.WC2TPCsize == 1)) continue;
        EventsWithWCMatch++;

        // If not picky track, skip
        if (!ev.wcTrackPicky) continue;
        EventsWithPickyWC++;

        // Check WC and front-face momentum
        double WCKE             = TMath::Sqrt(ev.WCTrackMomentum * ev.WCTrackMomentum + PionMass * PionMass) - PionMass;
        double calculatedEnLoss = energyLossCalculation();
        if (isData) {
            double tanThetaCosPhi = TMath::Tan(ev.WCTheta) * TMath::Cos(ev.WCPhi);
            double tanThetaSinPhi = TMath::Tan(ev.WCTheta) * TMath::Sin(ev.WCPhi);
            double den            = TMath::Sqrt(1 + tanThetaCosPhi * tanThetaCosPhi);
            double onTheFlyPz     = ev.WCTrackMomentum / den;
            double onTheFlyPx     = onTheFlyPz * tanThetaSinPhi;
            calculatedEnLoss      = energyLossCalculation(ev.WC4PrimaryX, onTheFlyPx, isData);
        }
        const double initialKE = WCKE - calculatedEnLoss;
        hFrontFaceKE->Fill(initialKE);
        hWCKE->Fill(WCKE);

        /////////////////////////////////////////
        // All this are valid events, analyze! //
        /////////////////////////////////////////

        // Sanity check
        removeRepeatedPoints(&ev.WC2TPCLocationsX, &ev.WC2TPCLocationsY, &ev.WC2TPCLocationsZ);

        // Copy WC2TPCLocations
        std::vector<double> wcX(ev.WC2TPCLocationsX);
        std::vector<double> wcY(ev.WC2TPCLocationsY);
        std::vector<double> wcZ(ev.WC2TPCLocationsZ);

        // Get direction to end cylinder
        int numPoints = wcX.size();
        int numTail   = std::min(10, numPoints - 1);
        std::vector<std::vector<double>> points;
        for (int j = numPoints - numTail; j < numPoints; ++j) {
            points.push_back({
                wcX.at(j),
                wcY.at(j),
                wcZ.at(j)
            });
        }
        if (numTail > 0) {
            std::vector<double> avgDir = getAverageDir(points);

            // Extrapolate track to end
            double scale = (maxZ - points.back()[2]) / avgDir[2];
            wcX.push_back(points.back()[0] + scale * avgDir[0]);
            wcY.push_back(points.back()[1] + scale * avgDir[1]);
            wcZ.push_back(points.back()[2] + scale * avgDir[2]);
        }

        // First, number of non-primary TG tracks and electron cut
        int numTGTracks = 0; int numSmallTracksInCylinder = 0;
        for (size_t trk_idx = 0; trk_idx < ev.recoBeginX.size(); ++trk_idx) {
            if (trkWCtoTPCMatch[trk_idx]) continue;

            // Get track length
            double trackLength = sqrt(
                pow(ev.recoEndX.at(trk_idx) - ev.recoBeginX.at(trk_idx), 2) +
                pow(ev.recoEndY.at(trk_idx) - ev.recoBeginY.at(trk_idx), 2) +
                pow(ev.recoEndZ.at(trk_idx) - ev.recoBeginZ.at(trk_idx), 2)
            );

            // Is track contained in 10 cm cylinder?
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
            if (startInCylinder && endInCylinder && (trackLength < CYLINDER_SMALL_TRACK)) numSmallTracksInCylinder++;

            if (
                !isWithinReducedVolume(ev.recoBeginX.at(trk_idx), ev.recoBeginY.at(trk_idx), ev.recoBeginZ.at(trk_idx)) &&
                !isWithinReducedVolume(ev.recoEndX.at(trk_idx), ev.recoEndY.at(trk_idx), ev.recoEndZ.at(trk_idx))
            ) numTGTracks++;
        }

        // Apply cuts
        if (numTGTracks > MAX_NUM_TG_TRACKS) continue;
        EventsPassingTG++;

        hSmallTracksInCylinder->Fill(numSmallTracksInCylinder);
        if (numSmallTracksInCylinder > ALLOWED_CYLINDER_SMALL_TRACKS) continue;
        EventsPassingSmall++;

        ///////////////////////
        // Primary track PID //
        ///////////////////////

        int totalCaloPoints = ev.wcMatchDEDX.size();
        int nRemoveOutliers = 2;
        int nRemoveEnds     = 3;
        int minPoints       = 5;

        // Get chi^2 fits, primary tracks are already checked for reversal in first module
        double pionChi2   = computeReducedChi2(gPion, ev.wcMatchResR,  ev.wcMatchDEDX, false, totalCaloPoints, nRemoveOutliers, nRemoveEnds);
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

        //////////////////////
        // Incident KE fill //
        //////////////////////

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
                hIncidentKE->Fill(initialKE - energyDeposited);
                hIncidentKEFine->Fill(initialKE - energyDeposited);
            }
        }
        double energyAtVertex = initialKE - energyDeposited;

        ////////////////////////
        // Reduced volume cut //
        ////////////////////////

        if (!isWithinReducedVolume(breakPointX, breakPointY, breakPointZ)) continue;
        EventsInRedVol++;

        /////////////////////
        // Primary PID cut //
        /////////////////////

        if (minChi2 == pionChi2 || minChi2 == protonChi2) continue;
        EventsPassingPrimaryPID++;

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

        //////////////////////////////////////////
        // Selection based off secondary tracks //
        //////////////////////////////////////////

        if (totalTaggedPions > 0) {
            // Reject events with > 1 pion
            if (totalTaggedPions > 1 || newSecondaryPion) continue;

            // Select as scatter
            EventsSelectedAsScatter++;
            hPionScatterKE->Fill(energyAtVertex);
            continue;
        }

        if (totalTaggedProtons > 0) {
            // Select as Np absorption
            EventsSelectedAsAbsNp++;
            hPionAbsNpKE->Fill(energyAtVertex);
            hPionAbsKE->Fill(energyAtVertex);
            continue;
        }

        EventsNoSecondary++;

        ////////////////////////////////////
        // Cluster non-reconstructed hits //
        ////////////////////////////////////

        // Get unordered set for hits in tracks
        std::unordered_set<int> hitsInTracks(ev.hitRecoAsTrackKey.begin(), ev.hitRecoAsTrackKey.end());

        // First, we construct the clusters 
        std::vector<int> candidateHits;
        for (size_t iHit = 0; iHit < std::min(nhits, kMaxHitsData); ++iHit) {
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
            
            for (int iAllHit = 0; iAllHit < std::min(nhits, kMaxHitsData); ++iAllHit) {
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

        // Now, get data for cut
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

        // Perform cut
        if (numClustersInduction < MAX_NUM_CLUSTERS_INDUCTION) {
            EventsSelectedAsAbs0p++;
            hPionAbs0pKE->Fill(energyAtVertex);
            hPionAbsKE->Fill(energyAtVertex);
            continue;
        }

        // We are done, all events that reach here are
        // rejected as background
    }

    // Diagnostics
    std::cout << std::endl;
    std::cout << "DIAGNOSTICS" << std::endl;
    std::cout << "=================================" << std::endl;

    std::cout << std::endl;
    std::cout << hMCSmallTrksInCylinderMuons->Integral() + hMCSmallTrksInCylinderPions->Integral() + hMCSmallTrksInCylinderElectrons->Integral() << " events in MC" << std::endl;
    std::cout << hSmallTracksInCylinder->Integral() << " events in data" << std::endl;
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout << "All events in data: " << TotalEvents << std::endl;
    std::cout << "    Events passing wire-chamber projection cuts: " << EventsPassingProj << "(" << ((double) EventsPassingProj / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events passing aperture checks:              " << EventsPassingAperture << "(" << ((double) EventsPassingAperture / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events passing TOF mass check:               " << EventsPassingTOFMass << "(" << ((double) EventsPassingTOFMass / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events with wire-chamber to TPC match:       " << EventsWithWCMatch << "(" << ((double) EventsWithWCMatch / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events with picky WC track:                  " << EventsWithPickyWC << "(" << ((double) EventsWithPickyWC / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events under TG tracks threshold:            " << EventsPassingTG << "(" << ((double) EventsPassingTG / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events under small tracks threshold:         " << EventsPassingSmall << "(" << ((double) EventsPassingSmall / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events with vertex in reduced volume:        " << EventsInRedVol << "(" << ((double) EventsInRedVol / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events passing primary PID:                  " << EventsPassingPrimaryPID << "(" << ((double) EventsPassingPrimaryPID / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events selected as scattering:               " << EventsSelectedAsScatter << "(" << ((double) EventsSelectedAsScatter / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events selected as absorption Np:            " << EventsSelectedAsAbsNp << "(" << ((double) EventsSelectedAsAbsNp / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events with no secondary tracks:             " << EventsNoSecondary << "(" << ((double) EventsNoSecondary / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events selected as absorption 0p:            " << EventsSelectedAsAbs0p << "(" << ((double) EventsSelectedAsAbs0p / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << std::endl;

    std::cout << "=================================" << std::endl;
    std::cout << std::endl;

    double SCALING_FACTOR     = hSmallTracksInCylinder->Integral() / (hMCSmallTrksInCylinderMuons->Integral() + hMCSmallTrksInCylinderPions->Integral() + hMCSmallTrksInCylinderElectrons->Integral());
    double SCALING_FACTOR_PRE = hFrontFaceKE->Integral() / (hMCFrontFaceKEPion->Integral() + hMCFrontFaceKEMuon->Integral() + hMCFrontFaceKEElectron->Integral());

    std::vector<TH1*> scaleByNormal = {
        hMCIncidentKEPion, hMCIncidentKEElectron, hMCIncidentKEMuon,
        hMCSmallTrksInCylinderPions, hMCSmallTrksInCylinderElectrons, hMCSmallTrksInCylinderMuons,
        hMCIncidentKEPionFine, hMCIncidentKEElectronFine, hMCIncidentKEMuonFine,
        hMCPionAbs0pKETrue, hMCPionAbs0pKEAbsNp, hMCPionAbs0pKEScatter, hMCPionAbs0pKEChExch, hMCPionAbs0pKEMuon, hMCPionAbs0pKEElectron, hMCPionAbs0pKEOther,
        hMCPionAbsNpKETrue, hMCPionAbsNpKEAbs0p, hMCPionAbsNpKEScatter, hMCPionAbsNpKEChExch, hMCPionAbsNpKEMuon, hMCPionAbsNpKEElectron, hMCPionAbsNpKEOther,
        hMCPionScatterKETrue, hMCPionScatterKEAbs0p, hMCPionScatterKEAbsNp, hMCPionScatterKEChExch, hMCPionScatterKEMuon, hMCPionScatterKEElectron, hMCPionScatterKEOther,
        hMCPionAbsKETrue, hMCPionAbsKEScatter, hMCPionAbsKEChExch, hMCPionAbsKEMuon, hMCPionAbsKEElectron, hMCPionAbsKEOther,
        hMCPionScatterKETrue2, hMCPionScatterKEAbs2, hMCPionScatterKEChExch2, hMCPionScatterKEMuon2, hMCPionScatterKEElectron2, hMCPionScatterKEOther2
    };

    std::vector<TH1*> scaleByPre = {
        hMCFrontFaceKEPion, hMCFrontFaceKEElectron, hMCFrontFaceKEMuon,
        hMCWCKEPion, hMCWCKEMuon, hMCWCKEElectron
    };

    auto scaleAll = [](const std::vector<TH1*>& hists, double factor) {
        for (auto* h : hists) {
            h->Scale(factor);
        }
    };

    scaleAll(scaleByNormal, SCALING_FACTOR);
    scaleAll(scaleByPre, SCALING_FACTOR_PRE);

    ///////////////////////////////////
    // Save histograms for unfolding //
    ///////////////////////////////////

    // Construct measured vector
    TVectorD Measure(NUM_SIGNAL_TYPES * NUM_BINS_KE);
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index      = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            Measure(index) = XSEC_UNITS * (RecoSignals[iSignal]->GetBinContent(iBin + 1) / hIncidentKE->GetBinContent(iBin + 1));
        }
    }
    TH1D* hMeasureVector = new TH1D(
        "hMeasureVector", "hMeasureVector;;", 
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    ); V2H(Measure, hMeasureVector);

    // Construct stat covariance matrix
    TMatrixD StatCovariance(NUM_SIGNAL_TYPES * NUM_BINS_KE, NUM_SIGNAL_TYPES * NUM_BINS_KE); StatCovariance.Zero();
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index   = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            double Ninc = hIncidentKE->GetBinContent(iBin + 1);
            double N    = Measure(index) * Ninc / XSEC_UNITS;

            double numSigma = (Ninc>0.0 && N>0.0 ? std::sqrt(N*(1.0 - N/Ninc)) : 0.0);
            double denSigma = (Ninc>0.0 ? std::sqrt(Ninc) : 0.0);

            double xs       = (N / Ninc) * XSEC_UNITS;
            double relVarXS = (N>0.0 && Ninc>0.0 ? std::pow(numSigma/N + denSigma/Ninc, 2) : 0.0);

            StatCovariance(index, index) = xs * xs * relVarXS;
        }
    }
    TH2D* hStatCovariance = new TH2D(
        "hStatCovariance", "Statistical Covariance Matrix;(i, #alpha);(j, #beta)",
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE, 
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE 
    ); M2H(StatCovariance, hStatCovariance);

    TFile* dataFile = new TFile("/exp/lariat/app/users/epelaez/histos/data/Reco.root", "RECREATE");

    // Construct measure and stat cov for abs + scatt measurement
    TVectorD MeasureAbsScatt(2 * NUM_BINS_KE);
    for (int iSignal = 0; iSignal < 2; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            MeasureAbsScatt(index) = XSEC_UNITS * (RecoSignalsAbsScatt[iSignal]->GetBinContent(iBin + 1) / hIncidentKE->GetBinContent(iBin + 1));
        }
    }
    TH1D* hMeasureVectorAbsScatt = new TH1D(
        "hMeasureVectorAbsScatt", "hMeasureVectorAbsScatt;;", 
        2 * NUM_BINS_KE, 0, 2 * NUM_BINS_KE
    ); V2H(MeasureAbsScatt, hMeasureVectorAbsScatt);

    TMatrixD StatCovarianceAbsScatt(2 * NUM_BINS_KE, 2 * NUM_BINS_KE); StatCovarianceAbsScatt.Zero();
    for (int iSignal = 0; iSignal < 2; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index   = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            double Ninc = hIncidentKE->GetBinContent(iBin + 1);
            double N    = MeasureAbsScatt(index) * Ninc / XSEC_UNITS;

            double numSigma = (Ninc>0.0 && N>0.0 ? std::sqrt(N*(1.0 - N/Ninc)) : 0.0);
            double denSigma = (Ninc>0.0 ? std::sqrt(Ninc) : 0.0);

            double xs       = (N / Ninc) * XSEC_UNITS;
            double relVarXS = (N>0.0 && Ninc>0.0 ? std::pow(numSigma/N + denSigma/Ninc, 2) : 0.0);

            StatCovarianceAbsScatt(index, index) = xs * xs * relVarXS;
        }
    }
    TH2D* hStatCovarianceAbsScatt = new TH2D(
        "hStatCovarianceAbsScatt", "Statistical Covariance Matrix;(i, #alpha);(j, #beta)",
        2 * NUM_BINS_KE, 0, 2 * NUM_BINS_KE, 
        2 * NUM_BINS_KE, 0, 2 * NUM_BINS_KE 
    ); M2H(StatCovarianceAbsScatt, hStatCovarianceAbsScatt);

    dataFile->cd();
    hMeasureVector->Write("", TObject::kOverwrite);
    hStatCovariance->Write("", TObject::kOverwrite);

    hMeasureVectorAbsScatt->Write("", TObject::kOverwrite);
    hStatCovarianceAbsScatt->Write("", TObject::kOverwrite);
    dataFile->Close();

    ////////////////////////////
    // Get raw cross-sections //
    ////////////////////////////

    //////////////////
    // Create plots //
    //////////////////

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

    std::vector<TH1*> PlotDataGroups = {
        // Cylinder
        hSmallTracksInCylinder,

        // Incident KE
        hIncidentKE,
        hIncidentKEFine,
        hFrontFaceKE,
        hWCKE,

        // Interacting KE
        hPionAbs0pKE,
        hPionAbsNpKE,
        hPionScatterKE,

        // Interacting KE (abs + scatt)
        hPionAbsKE,
        hPionScatterKE
    };

    std::vector<std::vector<TH1*>> PlotMCGroups = {
        // Cylinder
        {hMCSmallTrksInCylinderPions, hMCSmallTrksInCylinderMuons, hMCSmallTrksInCylinderElectrons},

        // Incident KE
        {hMCIncidentKEPion, hMCIncidentKEMuon, hMCIncidentKEElectron},
        {hMCIncidentKEPionFine, hMCIncidentKEMuonFine, hMCIncidentKEElectronFine},
        {hMCFrontFaceKEPion, hMCFrontFaceKEMuon, hMCFrontFaceKEElectron},
        {hMCWCKEPion, hMCWCKEMuon, hMCWCKEElectron},

        // Interacting KE
        {hMCPionAbs0pKETrue, hMCPionAbs0pKEAbsNp, hMCPionAbs0pKEScatter, hMCPionAbs0pKEChExch, hMCPionAbs0pKEMuon, hMCPionAbs0pKEElectron, hMCPionAbs0pKEOther},
        {hMCPionAbsNpKETrue, hMCPionAbsNpKEAbs0p, hMCPionAbsNpKEScatter, hMCPionAbsNpKEChExch, hMCPionAbsNpKEMuon, hMCPionAbsNpKEElectron, hMCPionAbsNpKEOther},
        {hMCPionScatterKETrue, hMCPionScatterKEAbs0p, hMCPionScatterKEAbsNp, hMCPionScatterKEChExch, hMCPionScatterKEMuon, hMCPionScatterKEElectron, hMCPionScatterKEOther},

        // Interacting KE (abs + scatt)
        {hMCPionAbsKETrue, hMCPionAbsKEScatter, hMCPionAbsKEChExch, hMCPionAbsKEMuon, hMCPionAbsKEElectron, hMCPionAbsKEOther},
        {hMCPionScatterKETrue2, hMCPionScatterKEAbs2, hMCPionScatterKEChExch2, hMCPionScatterKEMuon2, hMCPionScatterKEElectron2, hMCPionScatterKEOther2}
    };

    std::vector<std::vector<TString>> PlotMCLabelGroups = {
        // Cylinder
        {"Pion", "Muon", "Electron"},
        
        // Incident KE
        {"Pion", "Muon", "Electron"},
        {"Pion", "Muon", "Electron"},
        {"Pion", "Muon", "Electron"},
        {"Pion", "Muon", "Electron"},

        // Interacting KE
        {"True", "Abs Np", "Scatter", "Ch. Exch.", "Muon", "Electron", "Other"},
        {"True", "Abs 0p", "Scatter", "Ch. Exch.", "Muon", "Electron", "Other"},
        {"True", "Abs 0p", "Abs Np", "Ch. Exch.", "Muon", "Electron", "Other"},

        // Interacting KE (abs + scatt)
        {"True", "Scatter", "Ch. Exch.", "Muon", "Electron", "Other"},
        {"True", "Abs", "Ch. Exch.", "Muon", "Electron", "Other"},
    };

    std::vector<TString> PlotName = {
        // Cylinder 
        "Cylinder/SmallTracks",

        // Incident KE
        "Incident/IncidentKE",
        "Incident/IncidentKEFine",
        "Incident/FrontFaceKE",
        "Incident/WireChamberKE",

        // Interacting KE
        "Interacting/PionAbs0pKE",
        "Interacting/PionAbsNpKE",
        "Interacting/PionScatterKE",

        // Interacting KE (abs + scatt)
        "InteractingAbsScatt/PionAbsKE",
        "InteractingAbsScatt/PionScatterKE"
    };

    std::vector<TString> PlotTitle = {
        // Cylinder
        "Small Tracks in Cylinder",

        // Incident KE
        "Incident Kinetic Energy",
        "Incident Kinetic Energy (Fine Binning)",
        "Front-Face Kinetic Energy",
        "Wire-Chamber Kinetic Energy",

        // Interacting KE
        "Abs 0p Kinetic Energy",
        "Abs Np Kinetic Energy",
        "Scatter Kinetic Energy",

        // Interacting KE (abs + scatt)
        "Abs Kinetic Energy",
        "Scatter Kinetic Energy"
    };

    std::vector<TString> XLabels = {
        // Cylinder
        "# of small tracks in cylinder",

        // Incident KE
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",

        // Interacting KE
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",

        // Interacting KE (abs + scatt)
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]"
    };

    std::vector<TString> YLabels = {
        // Cylinder
        "Counts",

        // Incident KE
        "Counts",
        "Counts",
        "Counts",
        "Counts",

        // Interacting KE
        "Counts",
        "Counts",
        "Counts",

        // Interacting KE (abs + scatt)
        "Counts",
        "Counts"
    };

    for (size_t i = 0; i < PlotDataGroups.size(); ++i) {
        PrintDataVsMCContribPlot(
            SaveDir,
            PlotName[i],
            PlotDataGroups[i],
            PlotMCGroups[i],
            PlotMCLabelGroups[i],
            Colors,
            PlotTitle[i],
            XLabels[i],
            YLabels[i],
            FontStyle, TextSize,
            true,
            nullptr,
            true
        );
    }

    ////////////////////////////////////
    // Save everything from later use //
    ////////////////////////////////////

    TString outPath = "/exp/lariat/app/users/epelaez/histos/data/All.root";
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

    // 1D data plots
    for (auto* h : PlotDataGroups) writeOnce(h);

    outAll.Write();
    outAll.Close();

    std::cout << "\nWrote " << written.size() << " objects to " << outPath << std::endl;
}