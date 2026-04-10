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

struct EventVariables {
    // Picky track
    int wcTrackPicky;

    // WC match calorimetry
    std::vector<double> wcMatchResR;
    std::vector<double> wcMatchDEDX;
    std::vector<double> wcMatchEDep;
    std::vector<double> wcMatchXPos;
    std::vector<double> wcMatchYPos;
    std::vector<double> wcMatchZPos;

    // WC2TPC location
    std::vector<double> WC2TPCLocationsX;
    std::vector<double> WC2TPCLocationsY;
    std::vector<double> WC2TPCLocationsZ;

    // Reco track endpoints
    std::vector<double> recoEndX;
    std::vector<double> recoEndY;
    std::vector<double> recoEndZ;
    std::vector<double> recoBeginX;
    std::vector<double> recoBeginY;
    std::vector<double> recoBeginZ;
    std::vector<bool>   isTrackInverted;

    // Reco calorimetry
    std::vector<std::vector<double>> recoResR;
    std::vector<std::vector<double>> recoDEDX;

    // WC track scalars
    int    WC2TPCsize          = 0;
    bool   WC2TPCMatch         = false;
    double WCTrackMomentum     = -999;
    double WCTheta             = -999;
    double WCPhi               = -999;
    double WC4PrimaryX         = -999;
    double WC2TPCPrimaryBeginX = -999;
    double WC2TPCPrimaryBeginY = -999;
    double WC2TPCPrimaryBeginZ = -999;
    double WC2TPCPrimaryEndX   = -999;
    double WC2TPCPrimaryEndY   = -999;
    double WC2TPCPrimaryEndZ   = -999;

    // WC quality data
    std::vector<double> wcHit0;
    std::vector<double> wcHit1;
    std::vector<double> wcHit2;
    std::vector<double> wcHit3;

    // Hit information
    std::vector<int>   fHitPlane;
    std::vector<int>   hitRecoAsTrackKey;
    std::vector<int>   hitWC2TPCKey;
    std::vector<float> fHitT;
    std::vector<float> fHitX;
    std::vector<float> fHitW;

    // Reco track hits
    std::vector<std::vector<int>>    recoTrackHitIndices;
    std::vector<std::vector<double>> recoTrackHitX;
    std::vector<std::vector<double>> recoTrackHitY;
    std::vector<std::vector<double>> recoTrackHitZ;

    // TOF information
    double TOFMass   = -999;
    double tofObject = -999;
};

void RecoDataAnalysis() {
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
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis/figs/RecoDataAnalysis/";

    // Load file with NN data products
    TString RootFilePath = "/exp/lariat/app/users/epelaez/files/anatree_60a_data/histo.root";
    std::unique_ptr<TFile> File(TFile::Open(RootFilePath));
    TDirectory* Directory = (TDirectory*)File->Get("anatree");

    // Load file with MC comparison histograms
    TString MCCompFilePath = "/exp/lariat/app/users/epelaez/files/DataMCComparisons.root";
    std::unique_ptr<TFile> MCFile(TFile::Open(MCCompFilePath));

    ///////////////////
    // Load branches //
    ///////////////////

    ///////////////////
    // Load branches //
    ///////////////////

    // Load tree and branches
    TTree* tree = (TTree*) Directory->Get<TTree>("anatree");

    int run, subrun, event; bool isData = true;
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
    const int kMaxWCTracks = 100; // 1000
    float beamline_mass;                        tree->SetBranchAddress("beamline_mass",     &beamline_mass);
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
    float WC1xPos[kMaxWCTracks];                 tree->SetBranchAddress("WC1xPos",            &WC1xPos);
    float WC1yPos[kMaxWCTracks];                 tree->SetBranchAddress("WC1yPos",            &WC1yPos);
    float WC1zPos[kMaxWCTracks];                 tree->SetBranchAddress("WC1zPos",            &WC1zPos);
    float WC2xPos[kMaxWCTracks];                 tree->SetBranchAddress("WC2xPos",            &WC2xPos);
    float WC2yPos[kMaxWCTracks];                 tree->SetBranchAddress("WC2yPos",            &WC2yPos);
    float WC2zPos[kMaxWCTracks];                 tree->SetBranchAddress("WC2zPos",            &WC2zPos);
    float WC3xPos[kMaxWCTracks];                 tree->SetBranchAddress("WC3xPos",            &WC3xPos);
    float WC3yPos[kMaxWCTracks];                 tree->SetBranchAddress("WC3yPos",            &WC3yPos);
    float WC3zPos[kMaxWCTracks];                 tree->SetBranchAddress("WC3zPos",            &WC3zPos);
    float WC4xPos[kMaxWCTracks];                 tree->SetBranchAddress("WC4xPos",            &WC4xPos);
    float WC4yPos[kMaxWCTracks];                 tree->SetBranchAddress("WC4yPos",            &WC4yPos);
    float WC4zPos[kMaxWCTracks];                 tree->SetBranchAddress("WC4zPos",            &WC4zPos);

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

    // TOF object
    const int kMaxTOF = 10;
    int ntof;                     tree->SetBranchAddress("ntof",          &ntof);
    float tof[kMaxTOF];           tree->SetBranchAddress("tof",           tof);
    float tof_timestamp[kMaxTOF]; tree->SetBranchAddress("tof_timestamp", tof_timestamp);

    /////////////////////
    // Load histograms //
    /////////////////////

    TH1D* hMCNumWC2TPCMatch             = (TH1D*) MCFile->Get("hMCNumWC2TPCMatch");
    TH1D* hMCTGTrackLengths             = (TH1D*) MCFile->Get("hMCTGTrackLengths");
    TH1D* hMCTracksNearVertex           = (TH1D*) MCFile->Get("hMCTracksNearVertex");
    TH1D* hMCTrackLengthsNearVertex     = (TH1D*) MCFile->Get("hMCTrackLengthsNearVertex");
    TH1D* hMCNumTGTracks                = (TH1D*) MCFile->Get("hMCNumTGTracks");
    TH2D* hMCSmallVsTGTracks            = (TH2D*) MCFile->Get("hMCSmallVsTGTracks");
    TH2D* hMCTGNumSmallTracksVsThresh   = (TH2D*) MCFile->Get("hMCTGNumSmallTracksVsThresh");
    TH2D* hMCPrimaryTrackPosition       = (TH2D*) MCFile->Get("hMCPrimaryTrackPosition");

    TH1D* hMCNumTracksInCylinder    = (TH1D*) MCFile->Get("hMCNumTracksInCylinder");
    TH1D* hMCNumTracksInCylinder0TG = (TH1D*) MCFile->Get("hMCNumTracksInCylinder0TG");
    TH1D* hMCNumTracksInCylinder1TG = (TH1D*) MCFile->Get("hMCNumTracksInCylinder1TG");
    TH1D* hMCNumTracksInCylinder2TG = (TH1D*) MCFile->Get("hMCNumTracksInCylinder2TG");
    TH1D* hMCNumTracksInCylinderNTG = (TH1D*) MCFile->Get("hMCNumTracksInCylinderNTG");

    TH1D* hMCNumSmallTracksInCylinder    = (TH1D*) MCFile->Get("hMCNumSmallTracksInCylinder");
    TH1D* hMCNumSmallTracksInCylinder0TG = (TH1D*) MCFile->Get("hMCNumSmallTracksInCylinder0TG");
    TH1D* hMCNumSmallTracksInCylinder1TG = (TH1D*) MCFile->Get("hMCNumSmallTracksInCylinder1TG");
    TH1D* hMCNumSmallTracksInCylinder2TG = (TH1D*) MCFile->Get("hMCNumSmallTracksInCylinder2TG");
    TH1D* hMCNumSmallTracksInCylinderNTG = (TH1D*) MCFile->Get("hMCNumSmallTracksInCylinderNTG");

    TH1D* hMCTGSmallTracks    = (TH1D*) MCFile->Get("hMCTGSmallTracks");
    TH1D* hMCTGSmallTracks0TG = (TH1D*) MCFile->Get("hMCTGSmallTracks0TG");
    TH1D* hMCTGSmallTracks1TG = (TH1D*) MCFile->Get("hMCTGSmallTracks1TG");
    TH1D* hMCTGSmallTracks2TG = (TH1D*) MCFile->Get("hMCTGSmallTracks2TG");
    TH1D* hMCTGSmallTracksNTG = (TH1D*) MCFile->Get("hMCTGSmallTracksNTG");

    TH1D* hMCTGUnreconstructedHitsInduction    = (TH1D*) MCFile->Get("hMCTGUnreconstructedHitsInduction");
    TH1D* hMCTGUnreconstructedHitsInduction0TG = (TH1D*) MCFile->Get("hMCTGUnreconstructedHitsInduction0TG");
    TH1D* hMCTGUnreconstructedHitsInduction1TG = (TH1D*) MCFile->Get("hMCTGUnreconstructedHitsInduction1TG");
    TH1D* hMCTGUnreconstructedHitsInduction2TG = (TH1D*) MCFile->Get("hMCTGUnreconstructedHitsInduction2TG");
    TH1D* hMCTGUnreconstructedHitsInductionNTG = (TH1D*) MCFile->Get("hMCTGUnreconstructedHitsInductionNTG");

    TH1D* hMCTGUnreconstructedHitsCollection    = (TH1D*) MCFile->Get("hMCTGUnreconstructedHitsCollection");
    TH1D* hMCTGUnreconstructedHitsCollection0TG = (TH1D*) MCFile->Get("hMCTGUnreconstructedHitsCollection0TG");
    TH1D* hMCTGUnreconstructedHitsCollection1TG = (TH1D*) MCFile->Get("hMCTGUnreconstructedHitsCollection1TG");
    TH1D* hMCTGUnreconstructedHitsCollection2TG = (TH1D*) MCFile->Get("hMCTGUnreconstructedHitsCollection2TG");
    TH1D* hMCTGUnreconstructedHitsCollectionNTG = (TH1D*) MCFile->Get("hMCTGUnreconstructedHitsCollectionNTG");

    TH1D* hMCTGNumClustersInduction    = (TH1D*) MCFile->Get("hMCTGNumClustersInduction");
    TH1D* hMCTGNumClustersInduction0TG = (TH1D*) MCFile->Get("hMCTGNumClustersInduction0TG");
    TH1D* hMCTGNumClustersInduction1TG = (TH1D*) MCFile->Get("hMCTGNumClustersInduction1TG");
    TH1D* hMCTGNumClustersInduction2TG = (TH1D*) MCFile->Get("hMCTGNumClustersInduction2TG");
    TH1D* hMCTGNumClustersInductionNTG = (TH1D*) MCFile->Get("hMCTGNumClustersInductionNTG");

    TH1D* hMCTGClusterSizesInduction    = (TH1D*) MCFile->Get("hMCTGClusterSizesInduction");
    TH1D* hMCTGClusterSizesInduction0TG = (TH1D*) MCFile->Get("hMCTGClusterSizesInduction0TG");
    TH1D* hMCTGClusterSizesInduction1TG = (TH1D*) MCFile->Get("hMCTGClusterSizesInduction1TG");
    TH1D* hMCTGClusterSizesInduction2TG = (TH1D*) MCFile->Get("hMCTGClusterSizesInduction2TG");
    TH1D* hMCTGClusterSizesInductionNTG = (TH1D*) MCFile->Get("hMCTGClusterSizesInductionNTG");

    TH1D* hMCTGNumClustersCollection    = (TH1D*) MCFile->Get("hMCTGNumClustersCollection");
    TH1D* hMCTGNumClustersCollection0TG = (TH1D*) MCFile->Get("hMCTGNumClustersCollection0TG");
    TH1D* hMCTGNumClustersCollection1TG = (TH1D*) MCFile->Get("hMCTGNumClustersCollection1TG");
    TH1D* hMCTGNumClustersCollection2TG = (TH1D*) MCFile->Get("hMCTGNumClustersCollection2TG");
    TH1D* hMCTGNumClustersCollectionNTG = (TH1D*) MCFile->Get("hMCTGNumClustersCollectionNTG");

    TH1D* hMCTGClusterSizesCollection    = (TH1D*) MCFile->Get("hMCTGClusterSizesCollection");
    TH1D* hMCTGClusterSizesCollection0TG = (TH1D*) MCFile->Get("hMCTGClusterSizesCollection0TG");
    TH1D* hMCTGClusterSizesCollection1TG = (TH1D*) MCFile->Get("hMCTGClusterSizesCollection1TG");
    TH1D* hMCTGClusterSizesCollection2TG = (TH1D*) MCFile->Get("hMCTGClusterSizesCollection2TG");
    TH1D* hMCTGClusterSizesCollectionNTG = (TH1D*) MCFile->Get("hMCTGClusterSizesCollectionNTG");

    TH1D* hMCNumCandidateProtons = (TH1D*) MCFile->Get("hMCNumCandidateProtons");
    TH1D* hMCLengthCandidateProtons = (TH1D*) MCFile->Get("hMCLengthCandidateProtons");

    ///////////////////////
    // Create histograms //
    ///////////////////////

    // TOF
    TH1D* hTOFMass = new TH1D("hTOFMass", "TOF Mass Distribution", 50, 0, 1200);
    TH1D* hTOF     = new TH1D("hTOF", "TOF Distribution", 50, 0, 80);
    
    // WC
    TH1D* hNumWC2TPCMatch = new TH1D("hNumWC2TPCMatch", "NumWC2TPCMatch", 10, 0, 10);
    TH1D* hNumWCHits      = new TH1D("hNumWCHits", "hNumWCHits", 10, 0, 10);

    TH1D* hTGTrackLengths      = new TH1D("hTGTrackLengths", "TGTrackLengths", 25, 0, 50);

    TH1D* hNumTracksInCylinder    = new TH1D("hNumTracksInCylinder", "NumTracksInCylinder", 10, 0, 10);
    TH1D* hNumTracksInCylinder0TG = new TH1D("hNumTracksInCylinder0TG", "NumTracksInCylinder0TG", 10, 0, 10);
    TH1D* hNumTracksInCylinder1TG = new TH1D("hNumTracksInCylinder1TG", "NumTracksInCylinder1TG", 10, 0, 10);
    TH1D* hNumTracksInCylinder2TG = new TH1D("hNumTracksInCylinder2TG", "NumTracksInCylinder2TG", 10, 0, 10);
    TH1D* hNumTracksInCylinderNTG = new TH1D("hNumTracksInCylinderNTG", "NumTracksInCylinderNTG", 10, 0, 10);
    
    TH1D* hTGSmallTracks    = new TH1D("hTGSmallTracks", "TGSmallTracks", 10, 0, 10);
    TH1D* hTGSmallTracks0TG = new TH1D("hTGSmallTracks0TG", "TGSmallTracks0TG", 10, 0, 10);
    TH1D* hTGSmallTracks1TG = new TH1D("hTGSmallTracks1TG", "TGSmallTracks1TG", 10, 0, 10);
    TH1D* hTGSmallTracks2TG = new TH1D("hTGSmallTracks2TG", "TGSmallTracks2TG", 10, 0, 10);
    TH1D* hTGSmallTracksNTG = new TH1D("hTGSmallTracksNTG", "TGSmallTracksNTG", 10, 0, 10);

    TH1D* hNumSmallTracksInCylinder    = new TH1D("hNumSmallTracksInCylinder", "hNumSmallTracksInCylinder", 10, 0, 10);
    TH1D* hNumSmallTracksInCylinder0TG = new TH1D("hNumSmallTracksInCylinder0TG", "NumSmallTracksInCylinder0TG", 10, 0, 10);
    TH1D* hNumSmallTracksInCylinder1TG = new TH1D("hNumSmallTracksInCylinder1TG", "NumSmallTracksInCylinder1TG", 10, 0, 10);
    TH1D* hNumSmallTracksInCylinder2TG = new TH1D("hNumSmallTracksInCylinder2TG", "NumSmallTracksInCylinder2TG", 10, 0, 10);
    TH1D* hNumSmallTracksInCylinderNTG = new TH1D("hNumSmallTracksInCylinderNTG", "NumSmallTracksInCylinderNTG", 10, 0, 10);

    TH1D* hTGUnreconstructedHitsInduction    = new TH1D("hTGUnreconstructedHitsInduction", "hTGUnreconstructedHitsInduction;;", 30, 0, 30);
    TH1D* hTGUnreconstructedHitsInduction0TG = new TH1D("hTGUnreconstructedHitsInduction0TG", "UnreconstructedInductionHits0TG;;", 30, 0, 30);
    TH1D* hTGUnreconstructedHitsInduction1TG = new TH1D("hTGUnreconstructedHitsInduction1TG", "UnreconstructedInductionHits1TG;;", 30, 0, 30);
    TH1D* hTGUnreconstructedHitsInduction2TG = new TH1D("hTGUnreconstructedHitsInduction2TG", "UnreconstructedInductionHits2TG;;", 30, 0, 30);
    TH1D* hTGUnreconstructedHitsInductionNTG = new TH1D("hTGUnreconstructedHitsInductionNTG", "TGUnreconstructedHitsInductionNTG;;", 30, 0, 30);

    TH1D* hTGUnreconstructedHitsCollection    = new TH1D("hTGUnreconstructedHitsCollection", "hTGUnreconstructedHitsCollection;;", 30, 0, 30);
    TH1D* hTGUnreconstructedHitsCollection0TG = new TH1D("hTGUnreconstructedHitsCollection0TG", "UnreconstructedCollectionHits0TG;;", 30, 0, 30);
    TH1D* hTGUnreconstructedHitsCollection1TG = new TH1D("hTGUnreconstructedHitsCollection1TG", "UnreconstructedCollectionHits1TG;;", 30, 0, 30);
    TH1D* hTGUnreconstructedHitsCollection2TG = new TH1D("hTGUnreconstructedHitsCollection2TG", "UnreconstructedCollectionHits2TG;;", 30, 0, 30);
    TH1D* hTGUnreconstructedHitsCollectionNTG = new TH1D("hTGUnreconstructedHitsCollectionNTG", "TGUnreconstructedHitsCollectionNTG;;", 30, 0, 30);

    TH1D* hTGNumClustersInduction    = new TH1D("hTGNumClustersInduction", "hTGNumClustersInduction;;", 10, 0, 10);
    TH1D* hTGNumClustersInduction0TG = new TH1D("hTGNumClustersInduction0TG", "NumClustersInduction0TG;;", 10, 0, 10);
    TH1D* hTGNumClustersInduction1TG = new TH1D("hTGNumClustersInduction1TG", "NumClustersInduction1TG;;", 10, 0, 10);
    TH1D* hTGNumClustersInduction2TG = new TH1D("hTGNumClustersInduction2TG", "NumClustersInduction2TG;;", 10, 0, 10);
    TH1D* hTGNumClustersInductionNTG = new TH1D("hTGNumClustersInductionNTG", "TGNumClustersInductionNTG;;", 10, 0, 10);

    TH1D* hTGNumClustersCollection    = new TH1D("hTGNumClustersCollection", "hTGNumClustersCollection;;", 10, 0, 10);
    TH1D* hTGNumClustersCollection0TG = new TH1D("hTGNumClustersCollection0TG", "NumClustersCollection0TG;;", 10, 0, 10);
    TH1D* hTGNumClustersCollection1TG = new TH1D("hTGNumClustersCollection1TG", "NumClustersCollection1TG;;", 10, 0, 10);
    TH1D* hTGNumClustersCollection2TG = new TH1D("hTGNumClustersCollection2TG", "NumClustersCollection2TG;;", 10, 0, 10);
    TH1D* hTGNumClustersCollectionNTG = new TH1D("hTGNumClustersCollectionNTG", "TGNumClustersCollectionNTG;;", 10, 0, 10);

    TH1D* hTGClusterSizesInduction    = new TH1D("hTGClusterSizesInduction", "hTGClusterSizesInduction;;", 15, 0, 30);
    TH1D* hTGClusterSizesInduction0TG = new TH1D("hTGClusterSizesInduction0TG", "ClusterSizesInduction0TG;;", 15, 0, 30);
    TH1D* hTGClusterSizesInduction1TG = new TH1D("hTGClusterSizesInduction1TG", "ClusterSizesInduction1TG;;", 15, 0, 30);
    TH1D* hTGClusterSizesInduction2TG = new TH1D("hTGClusterSizesInduction2TG", "ClusterSizesInduction2TG;;", 15, 0, 30);
    TH1D* hTGClusterSizesInductionNTG = new TH1D("hTGClusterSizesInductionNTG", "TGClusterSizesInductionNTG;;", 15, 0, 30);

    TH1D* hTGClusterSizesCollection    = new TH1D("hTGClusterSizesCollection", "hTGClusterSizesCollection;;", 15, 0, 30);
    TH1D* hTGClusterSizesCollection0TG = new TH1D("hTGClusterSizesCollection0TG", "ClusterSizesCollection0TG;;", 15, 0, 30);
    TH1D* hTGClusterSizesCollection1TG = new TH1D("hTGClusterSizesCollection1TG", "ClusterSizesCollection1TG;;", 15, 0, 30);
    TH1D* hTGClusterSizesCollection2TG = new TH1D("hTGClusterSizesCollection2TG", "ClusterSizesCollection2TG;;", 15, 0, 30);
    TH1D* hTGClusterSizesCollectionNTG = new TH1D("hTGClusterSizesCollectionNTG", "TGClusterSizesCollectionNTG;;", 15, 0, 30);

    // With primary interacting
    TH1D* hTracksNearVertex       = new TH1D("hTracksNearVertex", "TracksNearVertex", 10, 0, 10);
    TH1D* hTrackLengthsNearVertex = new TH1D("hTrackLengthsNearVertex", "TrackLengthsNearVertex", 50, 0, 100);

    // TG and small tracks
    TH1D* hNumTGTracks              = new TH1D("hNumTGTracks", "NumTGTracks", 10, 0, 10);
    TH2D* hSmallVsTGTracks          = new TH2D("hSmallVsTGTracks", "SmallVsTGTracks;Small Tracks;TG Tracks", 15, 0, 15, 15, 0, 15);
    TH2D* hTGNumSmallTracksVsThresh = new TH2D("hTGNumSmallTracksVsThresh", "TGNumSmallTracksVsThresh;Small Track Length Threshold (cm);Num Small Tracks", 10, 0, 40, 15, 0, 15);
    TH2D* hTOFVsTOFMass             = new TH2D("hTOFVsTOFMass", "TOFVsTOFMass;TOF [ns];TOF Mass [MeV/c^2]", 35, 10, 80, 50, 0, 1200);

    // Kinematics
    TH1D* hNumCandidateProtons    = new TH1D("hNumCandidateProtons", "NumCandidateProtons;;", 10, 0, 10);
    TH1D* hLengthCandidateProtons = new TH1D("hLengthCandidateProtons", "LengthCandidateProtons;;", 40, 0, 80);

    ///////////////////////////////////
    // Distribution of primary track //
    ///////////////////////////////////

    TH2D* hPrimaryTrackPosition = new TH2D(
        "hPrimaryTrackPosition", "PrimaryTrackPosition;x-position;y-position",
        20, minX, maxX,
        20, minY, maxY
    );

    ///////////////////////////////////////
    // Distribution of background tracks //
    ///////////////////////////////////////

    TH2D* hBackgroundTracksPosition = new TH2D(
        "hBackgroundTracksPosition", "BackgroundTracksPosition;x-position;y-position",
        20, minX, maxX,
        20, minY, maxY
    );

    TH2D* hBackgroundTracksDirection = new TH2D(
        "hBackgroundTracksDirection", "BackgroundTracksDirection;#phi (x-y) [rad];#theta (forward) [rad]",
        20, -TMath::Pi(), TMath::Pi(),
        20, 0, TMath::Pi() / 2
    );

    /////////////////////////////
    // Distribution of WC hits //
    /////////////////////////////

    TH2D* hProjVsRealWC4 = new TH2D(
        "hProjVsRealWC4", "ProjVsRealWC4;p_{x} - WC_{x};p_{y} - WC_{y}",
        100, -40., 40.,
        100, -15., 15.
    );
    TH1D* hRadDistWC4 = new TH1D("hRadDistWC4", "RadDistWC4", 90, 0.0, 30.);

    TH2D* hMidPlaneCoinc = new TH2D(
        "hMidPlaneCoinc", "MidPlaneCoinc;up_{x} - down_{x};up_{y} - down_{y}",
        100, -20., 20.,
        100, -5., 5.
    );
    TH1D* hRadDistMidPlane = new TH1D("hRadDistMidPlane", "RadDistMidPlane", 60, 0.0, 20.);

    ///////////////////////////////////////////////////
    // Distribution of hits matched to primary track //
    ///////////////////////////////////////////////////

    TH1D* hTimePrimaryHits = new TH1D("hTimePrimaryHits", "hTimePrimaryHits", 60, -20, 300);

    //////////////////////
    // Loop over events //
    //////////////////////

    int candidateInteractingEvents = 0;

    Int_t NumEntries = (Int_t) tree->GetEntries();
    std::cout << "Num entries: " << NumEntries << std::endl;

    bool verbose = true;

    // Keep track of event counts for stat estimations
    int eventCount0TG = 0;
    int eventCount1TG = 0;
    int eventCount2TG = 0;
    int eventCountNTG = 0;

    int numValidEvents = 0;
    for (Int_t i = 0; i < NumEntries; ++i) {
        // For some reason, these events crash
        if (
            i == 199835 ||
            i == 311635
        ) continue;

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

        // Make script go faster
        if (i > USE_NUM_EVENTS) break;

        ////////////////////////////////
        // Load variables of interest //
        ////////////////////////////////

        // Load track information
        if (verbose) std::cout << std::endl;
        if (verbose) std::cout << "=================================" << std::endl;
        if (verbose) std::cout << "Event: " << event << std::endl;
        
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
        if (verbose) std::cout << "Number of matches: " << ev.WC2TPCsize  << std::endl;

        // Set beamline information
        ev.TOFMass   = std::abs(beamline_mass);
        if (ntof == 1) ev.tofObject = tof[0];

        if (verbose) std::cout << "Beamline mass: " << ev.TOFMass  << std::endl;
        if (verbose) std::cout << "TOF: " << ev.tofObject  << std::endl;

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

        // Fill for all events
        hTOFMass->Fill(std::abs(ev.TOFMass));
        hNumWC2TPCMatch->Fill(ev.WC2TPCsize);
        hTOF->Fill(ev.tofObject);
        hTOFVsTOFMass->Fill(ev.tofObject, std::abs(ev.TOFMass));

        /////////////////////////////////////////////////
        // Wire-chamber quality and other initial cuts //
        /////////////////////////////////////////////////

        // Candidate mass cut, keep pions, muons and electrons
        if (std::abs(ev.TOFMass) > PI_MU_EL_MASS_CUTOFF) continue;

        // If no track matched to wire-chamber, skip
        if (!(ev.WC2TPCMatch && ev.WC2TPCsize == 1)) continue;

        // If not picky track, skip
        if (!ev.wcTrackPicky) continue;

        // Project downstream to midplane @ -437.97 (without angular corrections)
        std::vector<double> midUp = projToZ(ev.wcHit0, ev.wcHit1, -437.97);
        // Use this point and WC3 to project up to WC4
        std::vector<double> projDown = projToZ(midUp, ev.wcHit2, -95.0);
        // Requires some corrections because magnets are not the same
        projDown[0] -= tan(1.32 * TMath::Pi() / 180.0) * (-95.0 - -437.97);

        // Compare x and y coordinate in projection and real hit for WC4
        hProjVsRealWC4->Fill(projDown[0] - ev.wcHit3.at(0), projDown[1] - ev.wcHit3.at(1));
        double radDistWC4 = TMath::Sqrt(pow(projDown[0] - ev.wcHit3.at(0), 2.) + pow(projDown[1] - ev.wcHit3.at(1), 2.));
        hRadDistWC4->Fill(radDistWC4);

        // Project upstream to midplane @ -437.97
        std::vector<double> midDown = projToZ(ev.wcHit2, ev.wcHit3, -437.97);
        midDown[0] -= tan(1.32 * TMath::Pi() / 180.0) * (-339.57 - -437.97);

        hMidPlaneCoinc->Fill(midUp[0] - midDown[0], midUp[1] - midDown[1]);
        double midPlaneDist = TMath::Sqrt(pow(midUp[0] - midDown[0] + 0.75, 2) + pow(midUp[1] - midDown[1], 2));
        // double midPlaneDist = TMath::Sqrt(pow(midUp[0] - midDown[0], 2) + pow(midUp[1] - midDown[1], 2));
        hRadDistMidPlane->Fill(midPlaneDist);

        if (radDistWC4 > MAX_RAD_DIST_WC4) continue;
        if (midPlaneDist > MAX_MID_PLANE_DIST) continue;

        // Check projected tracks go through all apertures
        bool Magnet1ApertureCheck = CheckUpstreamMagnetAperture(ev.wcHit0, ev.wcHit1);
        bool Magnet2ApertureCheck = CheckDownstreamMagnetAperture(ev.wcHit2, ev.wcHit3);
        bool DSColApertureCheck   = CheckDownstreamCollimatorAperture(ev.wcHit2, ev.wcHit3);

        if (!Magnet1ApertureCheck || !Magnet2ApertureCheck || !DSColApertureCheck) continue;

        //////////////////////////
        // Analyze valid events //
        //////////////////////////

        numValidEvents++;
        hPrimaryTrackPosition->Fill(ev.WC2TPCPrimaryBeginX, ev.WC2TPCPrimaryBeginY);

        // Check if WC2TPC is through-going
        bool isPrimaryTG = !isWithinReducedVolume(ev.WC2TPCPrimaryEndX, ev.WC2TPCPrimaryEndY, ev.WC2TPCPrimaryEndZ);

        // Loop over reconstructed tracks
        int numSmallTracks = 0;
        int numTracksNearVertex = 0;
        int smallTracksTPCStart = 0;
        int numTGTracks = 0;
        int numTracksInCylinder = 0;
        int numSmallTracksInCylinder = 0;

        for (size_t trk_idx = 0; trk_idx < ev.recoBeginX.size(); ++trk_idx) {
            if (trkWCtoTPCMatch[trk_idx]) continue;

            double distanceFromStart = distance(
                ev.recoBeginX.at(trk_idx), ev.WC2TPCPrimaryEndX, 
                ev.recoBeginY.at(trk_idx), ev.WC2TPCPrimaryEndY,
                ev.recoBeginZ.at(trk_idx), ev.WC2TPCPrimaryEndZ
            );
            double distanceFromEnd = distance(
                ev.recoEndX.at(trk_idx), ev.WC2TPCPrimaryEndX, 
                ev.recoEndY.at(trk_idx), ev.WC2TPCPrimaryEndY,
                ev.recoEndZ.at(trk_idx), ev.WC2TPCPrimaryEndZ
            );

            double trackLength = sqrt(
                pow(ev.recoEndX.at(trk_idx) - ev.recoBeginX.at(trk_idx), 2) +
                pow(ev.recoEndY.at(trk_idx) - ev.recoBeginY.at(trk_idx), 2) +
                pow(ev.recoEndZ.at(trk_idx) - ev.recoBeginZ.at(trk_idx), 2)
            );

            // Copy WC2TPCLocations
            std::vector<double>* wcX = new std::vector<double>(ev.WC2TPCLocationsX);
            std::vector<double>* wcY = new std::vector<double>(ev.WC2TPCLocationsY);
            std::vector<double>* wcZ = new std::vector<double>(ev.WC2TPCLocationsZ);

            // Get direction to end cylinder
            int numPoints = wcX->size();
            int numTail   = std::min(10, numPoints - 1);
            std::vector<std::vector<double>> points;
            for (int j = numPoints - numTail; j < numPoints; ++j) {
                points.push_back({
                    wcX->at(j),
                    wcY->at(j),
                    wcZ->at(j)
                });
            }
            if (numTail > 0) {
                std::vector<double> avgDir = getAverageDir(points);

                // Extrapolate track to end
                double scale = (maxZ - points.back()[2]) / avgDir[2];
                wcX->push_back(points.back()[0] + scale * avgDir[0]);
                wcY->push_back(points.back()[1] + scale * avgDir[1]);
                wcZ->push_back(points.back()[2] + scale * avgDir[2]);
            }

            // Is track contained in 10 cm cylinder?
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
            if (startInCylinder && endInCylinder) {
                numTracksInCylinder++;
                if (trackLength < CYLINDER_SMALL_TRACK) numSmallTracksInCylinder++;
            }

            if (isPrimaryTG) hTGTrackLengths->Fill(trackLength);

            if (
                ev.recoEndZ.at(trk_idx) < 30.0 && 
                ev.recoBeginZ.at(trk_idx) < 30.0
            ) {
                if (trackLength < SMALL_TRACK_LENGTH_CHEX) smallTracksTPCStart++;
            }

            if (
                !isWithinReducedVolume(ev.recoBeginX.at(trk_idx), ev.recoBeginY.at(trk_idx), ev.recoBeginZ.at(trk_idx)) &&
                !isWithinReducedVolume(ev.recoEndX.at(trk_idx), ev.recoEndY.at(trk_idx), ev.recoEndZ.at(trk_idx))
            ) {
                numTGTracks++;

                if (isPrimaryTG) {
                    // Track information about background tracks
                    bool isBeginTrue = ev.recoBeginZ.at(trk_idx) < ev.recoEndZ.at(trk_idx);

                    const double xs = isBeginTrue ? ev.recoBeginX.at(trk_idx) : ev.recoEndX.at(trk_idx);
                    const double ys = isBeginTrue ? ev.recoBeginY.at(trk_idx) : ev.recoEndY.at(trk_idx);
                    const double zs = isBeginTrue ? ev.recoBeginZ.at(trk_idx) : ev.recoEndZ.at(trk_idx);

                    const double xe = isBeginTrue ? ev.recoEndX.at(trk_idx)   : ev.recoBeginX.at(trk_idx);
                    const double ye = isBeginTrue ? ev.recoEndY.at(trk_idx)   : ev.recoBeginY.at(trk_idx);
                    const double ze = isBeginTrue ? ev.recoEndZ.at(trk_idx)   : ev.recoBeginZ.at(trk_idx);

                    hBackgroundTracksPosition->Fill(xs, ys);

                    auto [phi, theta] = azimuth_polar_from_points(
                        xs, ys, zs,
                        xe, ye, ze
                    );
                    hBackgroundTracksDirection->Fill(phi, theta);
                }
            }

            if (
                !isPrimaryTG &&
                (distanceFromStart < VERTEX_RADIUS || distanceFromEnd < VERTEX_RADIUS)
            ) {
                numTracksNearVertex++;
                hTrackLengthsNearVertex->Fill(trackLength);
            }
            if (trackLength < SMALL_TRACK_LENGTH_CHEX) numSmallTracks++;
        }

        // Look at relevant quantities with different number of TG tracks
        if (numTGTracks == 0 && isPrimaryTG) {
            hNumTracksInCylinder0TG->Fill(numTracksInCylinder);
            hNumSmallTracksInCylinder0TG->Fill(numSmallTracksInCylinder);
            hTGSmallTracks0TG->Fill(numSmallTracks);
        } 
        if (numTGTracks <= 1 && isPrimaryTG) {
            hNumTracksInCylinder1TG->Fill(numTracksInCylinder);
            hNumSmallTracksInCylinder1TG->Fill(numSmallTracksInCylinder);
            hTGSmallTracks1TG->Fill(numSmallTracks);
        }
        if (numTGTracks <= 2 && isPrimaryTG) {
            hNumTracksInCylinder2TG->Fill(numTracksInCylinder);
            hNumSmallTracksInCylinder2TG->Fill(numSmallTracksInCylinder);
            hTGSmallTracks2TG->Fill(numSmallTracks);
        }
        if (numTGTracks <= N_TG_TRACKS && isPrimaryTG) {
            hNumTracksInCylinderNTG->Fill(numTracksInCylinder);
            hNumSmallTracksInCylinderNTG->Fill(numSmallTracksInCylinder);
            hTGSmallTracksNTG->Fill(numSmallTracks);
        }

        hNumTGTracks->Fill(numTGTracks);

        if (!isPrimaryTG) hTracksNearVertex->Fill(numTracksNearVertex);

        // Add to histogram of small tracks if primary is throughgoing
        if (isPrimaryTG) {
            hSmallVsTGTracks->Fill(numSmallTracks, numTGTracks);

            hTGSmallTracks->Fill(numSmallTracks);
            hNumTracksInCylinder->Fill(numTracksInCylinder);
            hNumSmallTracksInCylinder->Fill(numSmallTracksInCylinder);

            // Scan over small track length thresholds and fill 2D histogram
            for (int threshBin = 1; threshBin <= hTGNumSmallTracksVsThresh->GetNbinsX(); ++threshBin) {
                double threshold = hTGNumSmallTracksVsThresh->GetXaxis()->GetBinCenter(threshBin);
                int nSmallTracks = 0;
                for (size_t trk_idx = 0; trk_idx < ev.recoBeginX.size(); ++trk_idx) {
                    if (trkWCtoTPCMatch[trk_idx]) continue;
                    double trackLength = sqrt(
                        pow(ev.recoEndX.at(trk_idx) - ev.recoBeginX.at(trk_idx), 2) +
                        pow(ev.recoEndY.at(trk_idx) - ev.recoBeginY.at(trk_idx), 2) +
                        pow(ev.recoEndZ.at(trk_idx) - ev.recoBeginZ.at(trk_idx), 2)
                    );
                    if (trackLength < threshold) nSmallTracks++;
                }
                hTGNumSmallTracksVsThresh->Fill(threshold, nSmallTracks);
            }
        }

        //////////////////////////
        // Unreconstructed hits //
        //////////////////////////

        // Look at unreconstructed hits
        std::unordered_set<int> hitsInTracks(ev.hitRecoAsTrackKey.begin(), ev.hitRecoAsTrackKey.end());

        // Reconstruct hit clusters
        std::vector<int> candidateHits;
        for (size_t iHit = 0; iHit < nhits; ++iHit) {
            // Skip hits already in tracks
            if (hitsInTracks.count(iHit) > 0) continue;

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
            )) candidateHits.push_back(iHit);
        }

        // Now cluster using those hits as starting points
        std::unordered_set<int> usedHits;
        std::vector<HitCluster> hitClusters;
        int nCandidateHits = candidateHits.size();

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

        // Separate primary hits into collection and induction
        std::vector<int> hitWC2TPCKeyInduction;
        std::vector<int> hitWC2TPCKeyCollection;

        for (size_t i = 0; i < ev.hitWC2TPCKey.size(); ++i) {
            // Only care about hits inside the reduced volume
            auto [i_hit, j_hit] = find_unique_position(&ev.recoTrackHitIndices, ev.hitWC2TPCKey.at(i));
            if (!isWithinReducedVolume(
                ev.recoTrackHitX.at(i_hit)[j_hit],
                ev.recoTrackHitY.at(i_hit)[j_hit],
                ev.recoTrackHitZ.at(i_hit)[j_hit]
            )) continue;

            if (ev.fHitPlane.at(ev.hitWC2TPCKey.at(i)) == 0) hitWC2TPCKeyInduction.push_back(ev.hitWC2TPCKey.at(i));
            else if (ev.fHitPlane.at(ev.hitWC2TPCKey.at(i)) == 1) hitWC2TPCKeyCollection.push_back(ev.hitWC2TPCKey.at(i));

            hTimePrimaryHits->Fill(ev.fHitT.at(ev.hitWC2TPCKey.at(i)));
        }

        // First, get a random point in the primary track
        std::vector<int> randomInduction; std::vector<int> randomCollection;

        if (hitWC2TPCKeyInduction.size() != 0) {
            randomInduction.push_back(0);

            for (int i = 1; i < (int) hitWC2TPCKeyInduction.size(); ++i) {
                double dX = std::abs(ev.fHitX.at(hitWC2TPCKeyInduction[randomInduction[randomInduction.size() - 1]]) - ev.fHitX.at(hitWC2TPCKeyInduction[i]));
                double dW = std::abs(ev.fHitW.at(hitWC2TPCKeyInduction[randomInduction[randomInduction.size() - 1]]) - ev.fHitW.at(hitWC2TPCKeyInduction[i]));

                if (dX > MAX_IN_CLUSTER_X_SEPARATION && dW > MAX_IN_CLUSTER_W_SEPARATION) randomInduction.push_back(i);
            }
        } 
        if (hitWC2TPCKeyCollection.size() != 0) {
            randomCollection.push_back(0);

            for (int i = 1; i < (int) hitWC2TPCKeyCollection.size(); ++i) {
                double dX = std::abs(ev.fHitX.at(hitWC2TPCKeyCollection[randomCollection[randomCollection.size() - 1]]) - ev.fHitX.at(hitWC2TPCKeyCollection[i]));
                double dW = std::abs(ev.fHitW.at(hitWC2TPCKeyCollection[randomCollection[randomCollection.size() - 1]]) - ev.fHitW.at(hitWC2TPCKeyCollection[i]));

                if (dX > MAX_IN_CLUSTER_X_SEPARATION && dW > MAX_IN_CLUSTER_W_SEPARATION) randomCollection.push_back(i);
            }
        }

        int numUnrecoHitsInduction  = 0;
        int numUnrecoHitsCollection = 0;

        std::vector<int> candidateRandomHitsInduction;
        std::vector<int> candidateRandomHitsCollection;

        for (size_t iHit = 0; iHit < nhits; ++iHit) {
            // Skip hits already in tracks
            if (hitsInTracks.count(iHit) > 0) continue;

            double hitX     = ev.fHitX.at(iHit);
            double hitW     = ev.fHitW.at(iHit);
            int    hitPlane = ev.fHitPlane.at(iHit);

            if (hitPlane == 0 && hitWC2TPCKeyInduction.size() > 0) {
                for (int idx : randomInduction) {
                    double dW = (hitW - ev.fHitW.at(hitWC2TPCKeyInduction[idx]));
                    double dX = (hitX - ev.fHitX.at(hitWC2TPCKeyInduction[idx]));
                    double d  = std::sqrt(std::pow(dW, 2) + std::pow(dX, 2));
                    if (d < DISTANCE_TO_PRIMARY_THRESHOLD) {
                        numUnrecoHitsInduction++;
                        candidateRandomHitsInduction.push_back(iHit);
                    }
                }
            } else if (hitPlane == 1 && hitWC2TPCKeyCollection.size() > 0) {
                for (int idx : randomCollection) {
                    double dW = (hitW - ev.fHitW.at(hitWC2TPCKeyCollection[idx]));
                    double dX = (hitX - ev.fHitX.at(hitWC2TPCKeyCollection[idx]));
                    double d  = std::sqrt(std::pow(dW, 2) + std::pow(dX, 2));
                    if (d < DISTANCE_TO_PRIMARY_THRESHOLD) {
                        numUnrecoHitsCollection++;
                        candidateRandomHitsCollection.push_back(iHit);
                    }
                }
            }
        }

        // Loop through clusters and see which would be close to random point
        int numClustersInduction = 0; int numClustersCollection = 0;

        for (size_t iCluster = 0; iCluster < hitClusters.size(); ++iCluster) {
            HitCluster thisCluster = hitClusters[iCluster];

            for (size_t iHit = 0; iHit < thisCluster.hitKeys.size(); ++iHit) {
                if (std::find(candidateRandomHitsInduction.begin(), candidateRandomHitsInduction.end(), thisCluster.hitKeys[iHit]) != candidateRandomHitsInduction.end()) {
                    numClustersInduction++;
                    if (isPrimaryTG) hTGClusterSizesInduction->Fill(thisCluster.clusterSize);
                    if (numTGTracks == 0 && isPrimaryTG) hTGClusterSizesInduction0TG->Fill(thisCluster.clusterSize);
                    if (numTGTracks <= 1 && isPrimaryTG) hTGClusterSizesInduction1TG->Fill(thisCluster.clusterSize);
                    if (numTGTracks <= 2 && isPrimaryTG) hTGClusterSizesInduction2TG->Fill(thisCluster.clusterSize);
                    if (numTGTracks <= N_TG_TRACKS && isPrimaryTG) hTGClusterSizesInductionNTG->Fill(thisCluster.clusterSize);
                    break;
                } else if (std::find(candidateRandomHitsCollection.begin(), candidateRandomHitsCollection.end(), thisCluster.hitKeys[iHit]) != candidateRandomHitsCollection.end()) {
                    numClustersCollection++;
                    if (isPrimaryTG) hTGClusterSizesCollection->Fill(thisCluster.clusterSize);
                    if (numTGTracks == 0 && isPrimaryTG) hTGClusterSizesCollection0TG->Fill(thisCluster.clusterSize);
                    if (numTGTracks <= 1 && isPrimaryTG) hTGClusterSizesCollection1TG->Fill(thisCluster.clusterSize);
                    if (numTGTracks <= 2 && isPrimaryTG) hTGClusterSizesCollection2TG->Fill(thisCluster.clusterSize);
                    if (numTGTracks <= N_TG_TRACKS && isPrimaryTG) hTGClusterSizesCollectionNTG->Fill(thisCluster.clusterSize);
                    break;
                }
            }
        }

        if (isPrimaryTG) {
            hTGNumClustersInduction->Fill(numClustersInduction);
            hTGNumClustersCollection->Fill(numClustersCollection);
            hTGUnreconstructedHitsInduction->Fill(numUnrecoHitsInduction);
            hTGUnreconstructedHitsCollection->Fill(numUnrecoHitsCollection);
        }
        if (numTGTracks == 0 && isPrimaryTG) {
            hTGNumClustersInduction0TG->Fill(numClustersInduction);
            hTGNumClustersCollection0TG->Fill(numClustersCollection);
            hTGUnreconstructedHitsInduction0TG->Fill(numUnrecoHitsInduction);
            hTGUnreconstructedHitsCollection0TG->Fill(numUnrecoHitsCollection);
        }
        if (numTGTracks <= 1 && isPrimaryTG) {
            hTGNumClustersInduction1TG->Fill(numClustersInduction);
            hTGNumClustersCollection1TG->Fill(numClustersCollection);
            hTGUnreconstructedHitsInduction1TG->Fill(numUnrecoHitsInduction);
            hTGUnreconstructedHitsCollection1TG->Fill(numUnrecoHitsCollection);
        }
        if (numTGTracks <= 2 && isPrimaryTG) {
            hTGNumClustersInduction2TG->Fill(numClustersInduction);
            hTGNumClustersCollection2TG->Fill(numClustersCollection);
            hTGUnreconstructedHitsInduction2TG->Fill(numUnrecoHitsInduction);
            hTGUnreconstructedHitsCollection2TG->Fill(numUnrecoHitsCollection);
        }
        if (numTGTracks <= N_TG_TRACKS && isPrimaryTG) {
            hTGNumClustersInductionNTG->Fill(numClustersInduction);
            hTGNumClustersCollectionNTG->Fill(numClustersCollection);
            hTGUnreconstructedHitsInductionNTG->Fill(numUnrecoHitsInduction);
            hTGUnreconstructedHitsCollectionNTG->Fill(numUnrecoHitsCollection);
        }

        //////////////////
        // TG track cut //
        //////////////////

        // Grab data about number of events with each cutoff
        if (numTGTracks <= 0) eventCount0TG++;
        if (numTGTracks <= 1) eventCount1TG++;
        if (numTGTracks <= 2) eventCount2TG++;
        if (numTGTracks <= 2) eventCountNTG++;

        if (numTGTracks > MAX_NUM_TG_TRACKS) continue;
        
        //////////////////
        // Cylinder cut //
        //////////////////

        if (numSmallTracksInCylinder > ALLOWED_CYLINDER_SMALL_TRACKS) continue;
        
        ////////////////////////
        // Reduced volume cut //
        ////////////////////////

        candidateInteractingEvents++;

        if (isPrimaryTG) continue;

        ///////////////////////////////////
        // Secondary particle kinematics //
        ///////////////////////////////////

        int numCandidateProtons = 0;

        for (size_t trk_idx = 0; trk_idx < ev.recoBeginX.size(); ++trk_idx) {
            if (trkWCtoTPCMatch[trk_idx]) continue;

            double distanceFromStart = distance(
                ev.recoBeginX.at(trk_idx), ev.WC2TPCPrimaryEndX, 
                ev.recoBeginY.at(trk_idx), ev.WC2TPCPrimaryEndY,
                ev.recoBeginZ.at(trk_idx), ev.WC2TPCPrimaryEndZ
            );
            double distanceFromEnd = distance(
                ev.recoEndX.at(trk_idx), ev.WC2TPCPrimaryEndX, 
                ev.recoEndY.at(trk_idx), ev.WC2TPCPrimaryEndY,
                ev.recoEndZ.at(trk_idx), ev.WC2TPCPrimaryEndZ
            );

            double trackLength = sqrt(
                pow(ev.recoEndX.at(trk_idx) - ev.recoBeginX.at(trk_idx), 2) +
                pow(ev.recoEndY.at(trk_idx) - ev.recoBeginY.at(trk_idx), 2) +
                pow(ev.recoEndZ.at(trk_idx) - ev.recoBeginZ.at(trk_idx), 2)
            );

            if (distanceFromStart < VERTEX_RADIUS || distanceFromEnd < VERTEX_RADIUS) {
                numCandidateProtons++;
                hLengthCandidateProtons->Fill(trackLength);
            }
        }
        hNumCandidateProtons->Fill(numCandidateProtons);
    }

    std::cout << std::endl;
    std::cout << "Number of events with at most 0 TG tracks: " << eventCount0TG << std::endl;
    std::cout << "Number of events with at most 1 TG tracks: " << eventCount1TG << std::endl;
    std::cout << "Number of events with at most 2 TG tracks: " << eventCount2TG << std::endl;
    std::cout << "Number of events with at most " << N_TG_TRACKS << " TG tracks: " << eventCountNTG << std::endl;
    std::cout << std::endl;

    std::cout << "Number of events before RV cut: " << candidateInteractingEvents << std::endl;
    std::cout << std::endl;
    
    double numMCEvents   = hMCNumTGTracks->Integral(0, hMCNumTGTracks->GetNbinsX() + 1);
    double numDataEvents = hNumTGTracks->Integral(0, hNumTGTracks->GetNbinsX() + 1);
    double scaling       = numDataEvents / numMCEvents;
    std::cout << "Num valid MC events: " << numMCEvents << std::endl;
    std::cout << "Num valid data events: " << numDataEvents << std::endl;
    std::cout << "Scaling: " << scaling << std::endl;
    std::cout << std::endl;

    double numMCTGEvents   = hMCTGSmallTracks->Integral(0, hMCTGSmallTracks->GetNbinsX() + 1);
    double numDataTGEvents = hTGSmallTracks->Integral(0, hTGSmallTracks->GetNbinsX() + 1);
    double scalingTG       = numDataTGEvents / numMCTGEvents;
    std::cout << "Num valid MC TG events: " << numMCTGEvents << std::endl;
    std::cout << "Num valid data TG events: " << numDataTGEvents << std::endl;
    std::cout << "Scaling TG: " << scalingTG << std::endl;
    std::cout << std::endl;

    double numMCNotTGEvents = hMCTracksNearVertex->Integral(0, hMCTracksNearVertex->GetNbinsX() + 1);
    double numDataNoTGvents = hTracksNearVertex->Integral(0, hTracksNearVertex->GetNbinsX() + 1);
    double scalingNoTG      = numDataNoTGvents / numMCNotTGEvents;
    std::cout << "Num valid MC not TG events: " << numMCNotTGEvents << std::endl;
    std::cout << "Num valid data not TG events: " << numDataNoTGvents << std::endl;
    std::cout << "Scaling not TG: " << scalingNoTG << std::endl;
    std::cout << std::endl;

    std::cout << "Num valid data TG events with no other TG: " << hNumTracksInCylinder0TG->Integral() << std::endl;
    std::cout << "Num valid data TG events with at most one more TG: " << hNumTracksInCylinder1TG->Integral() << std::endl;
    std::cout << "Num valid data TG events with at most two more TG: " << hNumTracksInCylinder2TG->Integral() << std::endl;
    std::cout << std::endl;

    double scaling0TG = hNumTracksInCylinder0TG->Integral() / hMCNumTracksInCylinder0TG->Integral();
    double scaling1TG = hNumTracksInCylinder1TG->Integral() / hMCNumTracksInCylinder1TG->Integral();
    double scaling2TG = hNumTracksInCylinder2TG->Integral() / hMCNumTracksInCylinder2TG->Integral();
    double scalingNTG = hNumTracksInCylinderNTG->Integral() / hMCNumTracksInCylinderNTG->Integral();
    double scalingKin = hNumCandidateProtons->Integral() / hMCNumCandidateProtons->Integral();

    std::cout << "Scaling 0 TG: " << scaling0TG << std::endl;
    std::cout << "Scaling 1 TG: " << scaling1TG << std::endl;
    std::cout << "Scaling 2 TG: " << scaling2TG << std::endl;
    std::cout << "Scaling N TG: " << scalingNTG << std::endl;
    std::cout << "Scaling Kin: " << scalingKin << std::endl;
    std::cout << std::endl;

    // Scale MC histograms to data histograms using event counts
    hMCNumWC2TPCMatch->Scale(hNumWC2TPCMatch->Integral() / hMCNumWC2TPCMatch->Integral());

    // all events
    hMCNumTGTracks->Scale(scaling);

    // primary not TG
    hMCTracksNearVertex->Scale(scalingNoTG);
    hMCTrackLengthsNearVertex->Scale(scalingNoTG);

    // primary TG
    hMCTGTrackLengths->Scale(scalingTG);

    hMCTGSmallTracks->Scale(scalingTG);
    hMCTGSmallTracks0TG->Scale(scaling0TG);
    hMCTGSmallTracks1TG->Scale(scaling1TG);
    hMCTGSmallTracks2TG->Scale(scaling2TG);
    hMCTGSmallTracksNTG->Scale(scalingNTG);

    hMCNumTracksInCylinder->Scale(scaling);
    hMCNumTracksInCylinder0TG->Scale(scaling0TG);
    hMCNumTracksInCylinder1TG->Scale(scaling1TG);
    hMCNumTracksInCylinder2TG->Scale(scaling2TG);
    hMCNumTracksInCylinderNTG->Scale(scalingNTG);

    hMCNumSmallTracksInCylinder->Scale(scaling);
    hMCNumSmallTracksInCylinder0TG->Scale(scaling0TG);
    hMCNumSmallTracksInCylinder1TG->Scale(scaling1TG);
    hMCNumSmallTracksInCylinder2TG->Scale(scaling2TG);
    hMCNumSmallTracksInCylinderNTG->Scale(scalingNTG);

    hMCTGUnreconstructedHitsInduction->Scale(scaling);
    hMCTGUnreconstructedHitsInduction0TG->Scale(scaling0TG);
    hMCTGUnreconstructedHitsInduction1TG->Scale(scaling1TG);
    hMCTGUnreconstructedHitsInduction2TG->Scale(scaling2TG);
    hMCTGUnreconstructedHitsInductionNTG->Scale(scalingNTG);

    hMCTGUnreconstructedHitsCollection->Scale(scaling);
    hMCTGUnreconstructedHitsCollection0TG->Scale(scaling0TG);
    hMCTGUnreconstructedHitsCollection1TG->Scale(scaling1TG);
    hMCTGUnreconstructedHitsCollection2TG->Scale(scaling2TG);
    hMCTGUnreconstructedHitsCollectionNTG->Scale(scalingNTG);

    hMCTGClusterSizesInduction->Scale(scaling);
    hMCTGClusterSizesInduction0TG->Scale(scaling0TG);
    hMCTGClusterSizesInduction1TG->Scale(scaling1TG);
    hMCTGClusterSizesInduction2TG->Scale(scaling2TG);
    hMCTGClusterSizesInductionNTG->Scale(scalingNTG);

    hMCTGClusterSizesCollection->Scale(scaling);
    hMCTGClusterSizesCollection0TG->Scale(scaling0TG);
    hMCTGClusterSizesCollection1TG->Scale(scaling1TG);
    hMCTGClusterSizesCollection2TG->Scale(scaling2TG);
    hMCTGClusterSizesCollectionNTG->Scale(scalingNTG);

    hMCTGNumClustersCollection->Scale(scaling);
    hMCTGNumClustersCollection0TG->Scale(scaling0TG);
    hMCTGNumClustersCollection1TG->Scale(scaling1TG);
    hMCTGNumClustersCollection2TG->Scale(scaling2TG);
    hMCTGNumClustersCollectionNTG->Scale(scalingNTG);

    hMCTGNumClustersInduction->Scale(scaling);
    hMCTGNumClustersInduction0TG->Scale(scaling0TG);
    hMCTGNumClustersInduction1TG->Scale(scaling1TG);
    hMCTGNumClustersInduction2TG->Scale(scaling2TG);
    hMCTGNumClustersInductionNTG->Scale(scalingNTG);

    // Secondary particle kinematics
    hMCNumCandidateProtons->Scale(scalingKin);
    hMCLengthCandidateProtons->Scale(scalingKin);

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

    std::vector<std::vector<TH1*>> PlotGroups = {
        // TOF
        {hTOFMass},
        {hTOF},

        // WC quality
        {hNumWC2TPCMatch, hMCNumWC2TPCMatch},
        {hNumWCHits},
        {hRadDistWC4},
        {hRadDistMidPlane},

        // Cylinder
        {hNumTracksInCylinder, hMCNumTracksInCylinder},
        {hNumTracksInCylinder0TG, hMCNumTracksInCylinder0TG},
        {hNumTracksInCylinder1TG, hMCNumTracksInCylinder1TG},
        {hNumTracksInCylinder2TG, hMCNumTracksInCylinder2TG},
        {hNumTracksInCylinderNTG, hMCNumTracksInCylinderNTG},

        {hNumSmallTracksInCylinder, hMCNumSmallTracksInCylinder},
        {hNumSmallTracksInCylinder0TG, hMCNumSmallTracksInCylinder0TG},
        {hNumSmallTracksInCylinder1TG, hMCNumSmallTracksInCylinder1TG},
        {hNumSmallTracksInCylinder2TG, hMCNumSmallTracksInCylinder2TG},
        {hNumSmallTracksInCylinderNTG, hMCNumSmallTracksInCylinderNTG},

        // Hits in primary
        {hTimePrimaryHits},

        // Comparisons with MC
        {hTracksNearVertex, hMCTracksNearVertex},
        {hTrackLengthsNearVertex, hMCTrackLengthsNearVertex},
        {hNumTGTracks, hMCNumTGTracks},
        {hTGTrackLengths, hMCTGTrackLengths},

        // Comparing TG threshold
        {hTGSmallTracks, hMCTGSmallTracks},
        {hTGSmallTracks0TG, hMCTGSmallTracks0TG},
        {hTGSmallTracks1TG, hMCTGSmallTracks1TG},
        {hTGSmallTracks2TG, hMCTGSmallTracks2TG},
        {hTGSmallTracksNTG, hMCTGSmallTracksNTG},

        // Unreconstructed hits
        {hTGUnreconstructedHitsInduction, hMCTGUnreconstructedHitsInduction},
        {hTGUnreconstructedHitsInduction0TG, hMCTGUnreconstructedHitsInduction0TG},
        {hTGUnreconstructedHitsInduction1TG, hMCTGUnreconstructedHitsInduction1TG},
        {hTGUnreconstructedHitsInduction2TG, hMCTGUnreconstructedHitsInduction2TG},
        {hTGUnreconstructedHitsInductionNTG, hMCTGUnreconstructedHitsInductionNTG},

        {hTGUnreconstructedHitsCollection, hMCTGUnreconstructedHitsCollection},
        {hTGUnreconstructedHitsCollection0TG, hMCTGUnreconstructedHitsCollection0TG},
        {hTGUnreconstructedHitsCollection1TG, hMCTGUnreconstructedHitsCollection1TG},
        {hTGUnreconstructedHitsCollection2TG, hMCTGUnreconstructedHitsCollection2TG},
        {hTGUnreconstructedHitsCollectionNTG, hMCTGUnreconstructedHitsCollectionNTG},

        {hTGClusterSizesCollection, hMCTGClusterSizesCollection},
        {hTGClusterSizesCollection0TG, hMCTGClusterSizesCollection0TG},
        {hTGClusterSizesCollection1TG, hMCTGClusterSizesCollection1TG},
        {hTGClusterSizesCollection2TG, hMCTGClusterSizesCollection2TG},
        {hTGClusterSizesCollectionNTG, hMCTGClusterSizesCollectionNTG},

        {hTGClusterSizesInduction, hMCTGClusterSizesInduction},
        {hTGClusterSizesInduction0TG, hMCTGClusterSizesInduction0TG},
        {hTGClusterSizesInduction1TG, hMCTGClusterSizesInduction1TG},
        {hTGClusterSizesInduction2TG, hMCTGClusterSizesInduction2TG},
        {hTGClusterSizesInductionNTG, hMCTGClusterSizesInductionNTG},

        {hTGNumClustersCollection, hMCTGNumClustersCollection},
        {hTGNumClustersCollection0TG, hMCTGNumClustersCollection0TG},
        {hTGNumClustersCollection1TG, hMCTGNumClustersCollection1TG},
        {hTGNumClustersCollection2TG, hMCTGNumClustersCollection2TG},
        {hTGNumClustersCollectionNTG, hMCTGNumClustersCollectionNTG},

        {hTGNumClustersInduction, hMCTGNumClustersInduction},
        {hTGNumClustersInduction0TG, hMCTGNumClustersInduction0TG},
        {hTGNumClustersInduction1TG, hMCTGNumClustersInduction1TG},
        {hTGNumClustersInduction2TG, hMCTGNumClustersInduction2TG},
        {hTGNumClustersInductionNTG, hMCTGNumClustersInductionNTG},

        // Secondary particle kinematics
        {hNumCandidateProtons, hMCNumCandidateProtons},
        {hLengthCandidateProtons, hMCLengthCandidateProtons}
    };

    std::vector<std::vector<TString>> PlotLabelGroups = {
        // TOF
        {"Data"},
        {"Data"},

        // WC quality
        {"Data", "MC (scaled)"},
        {"Data"},
        {"Data"},
        {"Data"},

        // Cylinder
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},

        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},

        // Hits in primary
        {"Data"},

        // Comparisons with MC
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data before", "Data after", "MC before (scaled)", "MC after (scaled)"},
        {"Data", "MC (scaled)"},

        // Comparing TG threshold
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},

        // Unreconstructed hits
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},

        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},

        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},

        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},

        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},

        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},

        // Secondary particle kinematics
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"}
    };

    std::vector<TString> PlotTitles = {
        // TOF
        "TOF/TOFMass",
        "TOF/TOF",

        // WC quality
        "WCQuality/NumMatches",
        "WCQuality/NumWCHits",
        "WCQuality/RadDistWC4Proj",
        "WCQuality/RadDistMidPlane",

        // Cylinder
        "Cylinder/NumTracksInCylinder",
        "Cylinder/0TG/NumTracksInCylinder0TG",
        "Cylinder/1TG/NumTracksInCylinder1TG",
        "Cylinder/2TG/NumTracksInCylinder2TG",
        "Cylinder/NTG/NumTracksInCylinderNTG",

        "Cylinder/NumSmallTracksInCylinder",
        "Cylinder/0TG/NumSmallTracksInCylinder0TG",
        "Cylinder/1TG/NumSmallTracksInCylinder1TG",
        "Cylinder/2TG/NumSmallTracksInCylinder2TG",
        "Cylinder/NTG/NumSmallTracksInCylinderNTG",

        // Hits in primary
        "Primary/HitTime",

        // Comparisons with MC
        "NearVertex/TracksNearVertex",
        "NearVertex/TrackLengthsNearVertex",
        "TGTracks/NumTGTracks",
        "Primary/TGTrackLengths",

        // Small tracks when primary is TG
        "Primary/TGSmallTracks",
        "Primary/0TG/TGSmallTracks0TG",
        "Primary/1TG/TGSmallTracks1TG",
        "Primary/2TG/TGSmallTracks2TG",
        "Primary/NTG/TGSmallTracksNTG",

        // Unreconstructed hits
        "UnrecoHits/UnreconstructedInductionHits",
        "UnrecoHits/0TG/UnreconstructedInductionHits0TG",
        "UnrecoHits/1TG/UnreconstructedInductionHits1TG",
        "UnrecoHits/2TG/UnreconstructedInductionHits2TG",
        "UnrecoHits/NTG/UnreconstructedInductionHitsNTG",

        "UnrecoHits/UnreconstructedCollectionHits",
        "UnrecoHits/0TG/UnreconstructedCollectionHits0TG",
        "UnrecoHits/1TG/UnreconstructedCollectionHits1TG",
        "UnrecoHits/2TG/UnreconstructedCollectionHits2TG",
        "UnrecoHits/NTG/UnreconstructedCollectionHitsNTG",

        "UnrecoHits/ClusterSizesCollection",
        "UnrecoHits/0TG/ClusterSizesCollection0TG",
        "UnrecoHits/1TG/ClusterSizesCollection1TG",
        "UnrecoHits/2TG/ClusterSizesCollection2TG",
        "UnrecoHits/NTG/ClusterSizesCollectionNTG",

        "UnrecoHits/ClusterSizesInduction",
        "UnrecoHits/0TG/ClusterSizesInduction0TG",
        "UnrecoHits/1TG/ClusterSizesInduction1TG",
        "UnrecoHits/2TG/ClusterSizesInduction2TG",
        "UnrecoHits/NTG/ClusterSizesInductionNTG",

        "UnrecoHits/NumClustersCollection",
        "UnrecoHits/0TG/NumClustersCollection0TG",
        "UnrecoHits/1TG/NumClustersCollection1TG",
        "UnrecoHits/2TG/NumClustersCollection2TG",
        "UnrecoHits/NTG/NumClustersCollectionNTG",

        "UnrecoHits/NumClustersInduction",
        "UnrecoHits/0TG/NumClustersInduction0TG",
        "UnrecoHits/1TG/NumClustersInduction1TG",
        "UnrecoHits/2TG/NumClustersInduction2TG",
        "UnrecoHits/NTG/NumClustersInductionNTG",

        // Secondary particle kinematics
        "Kinematics/NumCandidateProtons",
        "Kinematics/LengthsCandidateProtons"
    };

    std::vector<TString> XLabels = {
        // TOF
        "Mass [MeV/c^2]",
        "Time of flight [ns]",

        // WC quality
        "# of WC to TPC matches",
        "# of WC hits",
        "Proj. to WC hit distance [cm]",
        "Upstream to downstream proj. distance [cm]",

        // Cylinder
        "# of tracks",
        "# of tracks",
        "# of tracks",
        "# of tracks",
        "# of tracks",

        "# of small tracks",
        "# of small tracks",
        "# of small tracks",
        "# of small tracks",
        "# of small tracks",

        // Hits in primary
        "Hit time [us]",

        // Comparisons with MC
        "# of tracks near vertex",
        "Track length [cm]",
        "# of throughgoing tracks",
        "Track length [cm]",

        // Small tracks when primary is TG
        "# of small tracks",
        "# of small tracks",
        "# of small tracks",
        "# of small tracks",
        "# of small tracks",

        // Unreconstructed hits
        "# of unreconstructed hits",
        "# of unreconstructed hits",
        "# of unreconstructed hits",
        "# of unreconstructed hits",
        "# of unreconstructed hits",

        "# of unreconstructed hits",
        "# of unreconstructed hits",
        "# of unreconstructed hits",
        "# of unreconstructed hits",
        "# of unreconstructed hits",

        "Cluster size [cm]",
        "Cluster size [cm]",
        "Cluster size [cm]",
        "Cluster size [cm]",
        "Cluster size [cm]",

        "Cluster size [cm]",
        "Cluster size [cm]",
        "Cluster size [cm]",
        "Cluster size [cm]",
        "Cluster size [cm]",

        "Number of clusters",
        "Number of clusters",
        "Number of clusters",
        "Number of clusters",
        "Number of clusters",

        "Number of clusters",
        "Number of clusters",
        "Number of clusters",
        "Number of clusters",
        "Number of clusters",

        // Secondary particle kinematics
        "# of candidate proton",
        "Candidate proton length [cm]"
    };

    std::vector<TString> YLabels = {
        // TOF
        "Counts",
        "Counts",

        // WC quality
        "Counts",
        "Counts",
        "Counts",
        "Counts",

        // Cylinder
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

        // Hits in primary
        "Counts",

        // Comparisons with MC
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",

        // Small tracks when primary is TG
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",

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
        "Counts",

        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",

        // Secondary particle kinematics
        "Counts",
        "Counts"
    };

    std::vector<bool> PlotStacked = {
        // TOF
        false,
        false,

        // WC quality
        false,
        false,
        false,
        false,

        // Cylinder
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

        // Hits in primary
        false,

        // Comparisons with MC
        false,
        false,
        false,
        false,

        // Small tracks when primary is TG
        false,
        false,
        false,
        false,
        false,

        // Unreconstructed hits
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
        false,
        false,

        // Secondary particle kinematics
        false,
        false
    };

    std::vector<std::vector<bool>> PlotsAsPoints = {
        // Data plots
        {true},
        {true},

        // WC quality
        {true, false},
        {true},
        {true},
        {true},

        // Cylinder
        {true, false},
        {true, false},
        {true, false},
        {true, false},
        {true, false},

        {true, false},
        {true, false},
        {true, false},
        {true, false},
        {true, false},

        // Hits in primary
        {true},

        // Comparisons with MC
        {true, false},
        {true, false},
        {true, false},
        {true, false},

        // Small tracks when primary is TG
        {true, false},
        {true, false},
        {true, false},
        {true, false},
        {true, false},

        // Unreconstructed hits
        {true, false},
        {true, false},
        {true, false},
        {true, false},
        {true, false},

        {true, false},
        {true, false},
        {true, false},
        {true, false},
        {true, false},

        {true, false},
        {true, false},
        {true, false},
        {true, false},
        {true, false},

        {true, false},
        {true, false},
        {true, false},
        {true, false},
        {true, false},

        {true, false},
        {true, false},
        {true, false},
        {true, false},
        {true, false},

        {true, false},
        {true, false},
        {true, false},
        {true, false},
        {true, false},

        // Secondary particle kinematics
        {true, false},
        {true, false}
    };

    /////////////////////////////////
    // Fractional difference plots //
    /////////////////////////////////
    
    std::vector<std::pair<TH1*,TH1*>> PlotGroupsFracDiff = {
        // Primary
        {hTGSmallTracks, hMCTGSmallTracks},
        {hTGSmallTracks0TG, hMCTGSmallTracks0TG},
        {hTGSmallTracks1TG, hMCTGSmallTracks1TG},
        {hTGSmallTracks2TG, hMCTGSmallTracks2TG},
        {hTGSmallTracksNTG, hMCTGSmallTracksNTG},
        {hTGTrackLengths, hMCTGTrackLengths},

        // Cylinder
        {hNumTracksInCylinder, hMCNumTracksInCylinder},
        {hNumTracksInCylinder0TG, hMCNumTracksInCylinder0TG},
        {hNumTracksInCylinder1TG, hMCNumTracksInCylinder1TG},
        {hNumTracksInCylinder2TG, hMCNumTracksInCylinder2TG},
        {hNumTracksInCylinderNTG, hMCNumTracksInCylinderNTG},

        {hNumSmallTracksInCylinder, hMCNumSmallTracksInCylinder},
        {hNumSmallTracksInCylinder0TG, hMCNumSmallTracksInCylinder0TG},
        {hNumSmallTracksInCylinder1TG, hMCNumSmallTracksInCylinder1TG},
        {hNumSmallTracksInCylinder2TG, hMCNumSmallTracksInCylinder2TG},
        {hNumSmallTracksInCylinderNTG, hMCNumSmallTracksInCylinderNTG},

        // Near vertex
        {hTracksNearVertex, hMCTracksNearVertex},
        {hTrackLengthsNearVertex, hMCTrackLengthsNearVertex},

        // TG tracks
        {hNumTGTracks, hMCNumTGTracks},

        // Unreconstruced hits
        {hTGUnreconstructedHitsInduction, hMCTGUnreconstructedHitsInduction},
        {hTGUnreconstructedHitsInduction0TG, hMCTGUnreconstructedHitsInduction0TG},
        {hTGUnreconstructedHitsInduction1TG, hMCTGUnreconstructedHitsInduction1TG},
        {hTGUnreconstructedHitsInduction2TG, hMCTGUnreconstructedHitsInduction2TG},
        {hTGUnreconstructedHitsInductionNTG, hMCTGUnreconstructedHitsInductionNTG},

        {hTGUnreconstructedHitsCollection, hMCTGUnreconstructedHitsCollection},
        {hTGUnreconstructedHitsCollection0TG, hMCTGUnreconstructedHitsCollection0TG},
        {hTGUnreconstructedHitsCollection1TG, hMCTGUnreconstructedHitsCollection1TG},
        {hTGUnreconstructedHitsCollection2TG, hMCTGUnreconstructedHitsCollection2TG},
        {hTGUnreconstructedHitsCollectionNTG, hMCTGUnreconstructedHitsCollectionNTG},

        {hTGClusterSizesCollection, hMCTGClusterSizesCollection},
        {hTGClusterSizesCollection0TG, hMCTGClusterSizesCollection0TG},
        {hTGClusterSizesCollection1TG, hMCTGClusterSizesCollection1TG},
        {hTGClusterSizesCollection2TG, hMCTGClusterSizesCollection2TG},
        {hTGClusterSizesCollectionNTG, hMCTGClusterSizesCollectionNTG},

        {hTGClusterSizesInduction, hMCTGClusterSizesInduction},
        {hTGClusterSizesInduction0TG, hMCTGClusterSizesInduction0TG},
        {hTGClusterSizesInduction1TG, hMCTGClusterSizesInduction1TG},
        {hTGClusterSizesInduction2TG, hMCTGClusterSizesInduction2TG},
        {hTGClusterSizesInductionNTG, hMCTGClusterSizesInductionNTG},

        {hTGNumClustersCollection, hMCTGNumClustersCollection},
        {hTGNumClustersCollection0TG, hMCTGNumClustersCollection0TG},
        {hTGNumClustersCollection1TG, hMCTGNumClustersCollection1TG},
        {hTGNumClustersCollection2TG, hMCTGNumClustersCollection2TG},
        {hTGNumClustersCollectionNTG, hMCTGNumClustersCollectionNTG},

        {hTGNumClustersInduction, hMCTGNumClustersInduction},
        {hTGNumClustersInduction0TG, hMCTGNumClustersInduction0TG},
        {hTGNumClustersInduction1TG, hMCTGNumClustersInduction1TG},
        {hTGNumClustersInduction2TG, hMCTGNumClustersInduction2TG},
        {hTGNumClustersInductionNTG, hMCTGNumClustersInductionNTG},

        // Secondary particle kinematics
        {hNumCandidateProtons, hMCNumCandidateProtons},
        {hLengthCandidateProtons, hMCLengthCandidateProtons}
    };

    std::vector<std::string> FracDiffDirectory = {
        "Primary",
        "Primary/0TG",
        "Primary/1TG",
        "Primary/2TG",
        "Primary/NTG",
        "Primary",

        "Cylinder",
        "Cylinder/0TG",
        "Cylinder/1TG",
        "Cylinder/2TG",
        "Cylinder/NTG",

        "Cylinder",
        "Cylinder/0TG",
        "Cylinder/1TG",
        "Cylinder/2TG",
        "Cylinder/NTG",

        "NearVertex",
        "NearVertex",

        "TGTracks",

        "UnrecoHits",
        "UnrecoHits/0TG",
        "UnrecoHits/1TG",
        "UnrecoHits/2TG",
        "UnrecoHits/NTG",

        "UnrecoHits",
        "UnrecoHits/0TG",
        "UnrecoHits/1TG",
        "UnrecoHits/2TG",
        "UnrecoHits/NTG",

        "UnrecoHits",
        "UnrecoHits/0TG",
        "UnrecoHits/1TG",
        "UnrecoHits/2TG",
        "UnrecoHits/NTG",

        "UnrecoHits",
        "UnrecoHits/0TG",
        "UnrecoHits/1TG",
        "UnrecoHits/2TG",
        "UnrecoHits/NTG",

        "UnrecoHits",
        "UnrecoHits/0TG",
        "UnrecoHits/1TG",
        "UnrecoHits/2TG",
        "UnrecoHits/NTG",

        "UnrecoHits",
        "UnrecoHits/0TG",
        "UnrecoHits/1TG",
        "UnrecoHits/2TG",
        "UnrecoHits/NTG",

        "Kinematics",
        "Kinematics"
    };

    for (int i = 0; i < (int)PlotGroupsFracDiff.size(); ++i) {
        TH1* hData = PlotGroupsFracDiff[i].first;
        TH1* hMC   = PlotGroupsFracDiff[i].second;
        if (!hData || !hMC) continue;

        TH1D* hFrac = (TH1D*)hData->Clone(
            (std::string(hData->GetName()) + "_FracVsMC").c_str()
        );
        hFrac->Reset("ICE");
        hFrac->Sumw2(kTRUE);
        hFrac->GetYaxis()->SetTitle("(D - MC) / MC");

        const int nb = hData->GetNbinsX();
        for (int b = 1; b <= nb; ++b) {
            double D  = hData->GetBinContent(b);
            double S  = hMC->GetBinContent(b);
            double sD = hData->GetBinError(b);
            double sS = hMC->GetBinError(b);

            if (S != 0.0) {
                double f = (D - S) / S;
                // df/dD = 1/S ; df/dS = -D/S^2
                double dfdD = 1.0 / S;
                double dfdS = -D / (S*S);
                double var = dfdD*dfdD*sD*sD + dfdS*dfdS*sS*sS; // rho=0
                hFrac->SetBinContent(b, f);
                hFrac->SetBinError(b, (var > 0 ? std::sqrt(var) : 0.0));
            } else {
                hFrac->SetBinContent(b, 0.0);
                hFrac->SetBinError(b,   0.0);
            }
        }

        PlotGroups.push_back({hFrac});
        PlotLabelGroups.push_back({"Frac. diff."});
        PlotTitles.push_back(FracDiffDirectory[i] + "/" + std::string(hData->GetTitle()) + "FracDiff");
        XLabels.push_back(hData->GetXaxis()->GetTitle());
        YLabels.push_back("(D - MC) / MC");
        PlotStacked.push_back(false);
        PlotsAsPoints.push_back({true});
    }

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

    ///////////////////////////
    // Two-dimensional plots //
    ///////////////////////////

    hMCSmallVsTGTracks->Scale(scalingTG);
    hMCTGNumSmallTracksVsThresh->Scale(scalingTG);

    std::vector<TH2*> TwoDPlots = {
        hSmallVsTGTracks,
        hMCSmallVsTGTracks,
        hTGNumSmallTracksVsThresh,
        hMCTGNumSmallTracksVsThresh,
        hTOFVsTOFMass,
        hPrimaryTrackPosition,
        hMCPrimaryTrackPosition,
        hBackgroundTracksPosition,
        hBackgroundTracksDirection,
        hProjVsRealWC4,
        hMidPlaneCoinc
    };

    std::vector<TString> TwoDTitles = {
        "TGTracks/SmallVsTGTracks",
        "TGTracks/MCSmallVsTGTracks",
        "TGTracks/TGNumSmallTracksVsThresh",
        "TGTracks/MCTGNumSmallTracksVsThresh",
        "TOF/TOFVsTOFMass",
        "Primary/TrackIncidentPosition",
        "Primary/MCTrackIncidentPosition",
        "BkgTracks/IncidentPosition",
        "BkgTracks/IncidentDirection",
        "WCQuality/ProjVsRealWC4",
        "WCQuality/MidPlaneCoinc"
    };
    std::vector<std::pair<double,double>> TwoDRanges = {
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0}
    };

    std::vector<bool> TwoDDisplayNumbers = {
        true,
        true,
        true,
        true,
        false,
        false,
        false,
        false,
        false,
        false
    };

    printTwoDPlots(SaveDir, TwoDPlots, TwoDTitles, TwoDRanges, TwoDDisplayNumbers);
}