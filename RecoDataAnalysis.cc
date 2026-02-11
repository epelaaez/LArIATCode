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
    TString RootFilePath = "/exp/lariat/app/users/epelaez/files/DataAll_histo.root";
    std::unique_ptr<TFile> File(TFile::Open(RootFilePath));
    TDirectory* Directory = (TDirectory*)File->Get("RecoNNDataEval");

    // Load file with MC comparison histograms
    TString MCCompFilePath = "/exp/lariat/app/users/epelaez/files/DataMCComparisons.root";
    std::unique_ptr<TFile> MCFile(TFile::Open(MCCompFilePath));

    ///////////////////
    // Load branches //
    ///////////////////

    // Load tree and branches
    TTree* tree = (TTree*) Directory->Get<TTree>("RecoNNDataEvalTree");

    int run, subrun, event; bool isData;
    tree->SetBranchAddress("run", &run); 
    tree->SetBranchAddress("subrun", &subrun); 
    tree->SetBranchAddress("event", &event);
    tree->SetBranchAddress("isData", &isData);

    // WC match information
    int WC2TPCtrkID, WC2TPCsize;
    double WCTrackMomentum, WCTheta, WCPhi, WC4PrimaryX;
    double WC2TPCPrimaryBeginX, WC2TPCPrimaryBeginY, WC2TPCPrimaryBeginZ;
    double WC2TPCPrimaryEndX, WC2TPCPrimaryEndY, WC2TPCPrimaryEndZ;
    std::vector<double>* wcMatchResR = nullptr;
    std::vector<double>* wcMatchDEDX = nullptr;
    std::vector<double>* wcMatchEDep = nullptr;
    std::vector<double>* wcMatchXPos = nullptr;
    std::vector<double>* wcMatchYPos = nullptr;
    std::vector<double>* wcMatchZPos = nullptr;
    tree->SetBranchAddress("WC2TPCtrkID", &WC2TPCtrkID);
    tree->SetBranchAddress("WC2TPCsize", &WC2TPCsize);
    tree->SetBranchAddress("WCTrackMomentum", &WCTrackMomentum);
    tree->SetBranchAddress("WCTheta", &WCTheta);
    tree->SetBranchAddress("WCPhi", &WCPhi);
    tree->SetBranchAddress("WC4PrimaryX", &WC4PrimaryX);
    tree->SetBranchAddress("WC2TPCPrimaryBeginX", &WC2TPCPrimaryBeginX);
    tree->SetBranchAddress("WC2TPCPrimaryBeginY", &WC2TPCPrimaryBeginY);
    tree->SetBranchAddress("WC2TPCPrimaryBeginZ", &WC2TPCPrimaryBeginZ);
    tree->SetBranchAddress("WC2TPCPrimaryEndX", &WC2TPCPrimaryEndX);
    tree->SetBranchAddress("WC2TPCPrimaryEndY", &WC2TPCPrimaryEndY);
    tree->SetBranchAddress("WC2TPCPrimaryEndZ", &WC2TPCPrimaryEndZ);
    tree->SetBranchAddress("wcMatchResR", &wcMatchResR);
    tree->SetBranchAddress("wcMatchDEDX", &wcMatchDEDX);
    tree->SetBranchAddress("wcMatchEDep", &wcMatchEDep);
    tree->SetBranchAddress("wcMatchXPos", &wcMatchXPos);
    tree->SetBranchAddress("wcMatchYPos", &wcMatchYPos);
    tree->SetBranchAddress("wcMatchZPos", &wcMatchZPos);

    // WC match location information
    std::vector<double>* WC2TPCLocationsX = nullptr;
    std::vector<double>* WC2TPCLocationsY = nullptr;
    std::vector<double>* WC2TPCLocationsZ = nullptr;
    tree->SetBranchAddress("WC2TPCLocationsX", &WC2TPCLocationsX);
    tree->SetBranchAddress("WC2TPCLocationsY", &WC2TPCLocationsY);
    tree->SetBranchAddress("WC2TPCLocationsZ", &WC2TPCLocationsZ);

    // WC quality data
    int wcNumHits;
    std::vector<double>* wcHit0 = nullptr;
    std::vector<double>* wcHit1 = nullptr;
    std::vector<double>* wcHit2 = nullptr;
    std::vector<double>* wcHit3 = nullptr;
    tree->SetBranchAddress("wcNumHits", &wcNumHits);
    tree->SetBranchAddress("wcHit0", &wcHit0);
    tree->SetBranchAddress("wcHit1", &wcHit1);
    tree->SetBranchAddress("wcHit2", &wcHit2);
    tree->SetBranchAddress("wcHit3", &wcHit3);

    // Reco information
    std::vector<double>* recoMeanDEDX      = nullptr;
    std::vector<double>* recoEndX          = nullptr;
    std::vector<double>* recoEndY          = nullptr;
    std::vector<double>* recoEndZ          = nullptr;
    std::vector<double>* recoBeginX        = nullptr;
    std::vector<double>* recoBeginY        = nullptr;
    std::vector<double>* recoBeginZ        = nullptr;
    std::vector<int>*    recoTrkID         = nullptr;
    std::vector<bool>*   isTrackNearVertex = nullptr;
    std::vector<bool>*   isTrackInverted   = nullptr;
    tree->SetBranchAddress("recoMeanDEDX", &recoMeanDEDX);
    tree->SetBranchAddress("recoEndX", &recoEndX);
    tree->SetBranchAddress("recoEndY", &recoEndY);
    tree->SetBranchAddress("recoEndZ", &recoEndZ);
    tree->SetBranchAddress("recoBeginX", &recoBeginX);
    tree->SetBranchAddress("recoBeginY", &recoBeginY);
    tree->SetBranchAddress("recoBeginZ", &recoBeginZ);
    tree->SetBranchAddress("recoTrkID", &recoTrkID);
    tree->SetBranchAddress("isTrackNearVertex", &isTrackNearVertex);
    tree->SetBranchAddress("isTrackInverted", &isTrackInverted);

    // Calorimetry information
    std::vector<std::vector<double>>* recoResR = nullptr;
    std::vector<std::vector<double>>* recoDEDX = nullptr;
    tree->SetBranchAddress("recoResR", &recoResR);
    tree->SetBranchAddress("recoDEDX", &recoDEDX);

    // Individual hit information
    double primaryEndPointHitW, primaryEndPointHitX;
    std::vector<int>*   fHitKey = nullptr;
    std::vector<int>*   fHitPlane = nullptr;
    std::vector<int>*   hitRecoAsTrackKey = nullptr;
    std::vector<float>* fHitT = nullptr;
    std::vector<float>* fHitX = nullptr;
    std::vector<float>* fHitW = nullptr;
    std::vector<float>* fHitCharge = nullptr;
    std::vector<float>* fHitChargeCol = nullptr;
    std::vector<int>*   hitWC2TPCKey = nullptr;
    std::vector<int>*   hitThroughTrack = nullptr;
    tree->SetBranchAddress("primaryEndPointHitW", &primaryEndPointHitW);
    tree->SetBranchAddress("primaryEndPointHitX", &primaryEndPointHitX);
    tree->SetBranchAddress("fHitKey", &fHitKey);
    tree->SetBranchAddress("fHitPlane", &fHitPlane);
    tree->SetBranchAddress("hitRecoAsTrackKey", &hitRecoAsTrackKey);
    tree->SetBranchAddress("fHitT", &fHitT);
    tree->SetBranchAddress("fHitX", &fHitX);
    tree->SetBranchAddress("fHitW", &fHitW);
    tree->SetBranchAddress("fHitCharge", &fHitCharge);
    tree->SetBranchAddress("fHitChargeCol", &fHitChargeCol);
    tree->SetBranchAddress("hitWC2TPCKey", &hitWC2TPCKey);
    tree->SetBranchAddress("hitThroughTrack", &hitThroughTrack);

    // Information about hits in tracks
    std::vector<std::vector<int>>* recoTrackHitIndices  = nullptr;
    std::vector<std::vector<double>>*     recoTrackHitX = nullptr; 
    std::vector<std::vector<double>>*     recoTrackHitY = nullptr; 
    std::vector<std::vector<double>>*     recoTrackHitZ = nullptr; 
    tree->SetBranchAddress("recoTrackHitIndices", &recoTrackHitIndices);
    tree->SetBranchAddress("recoTrackHitX", &recoTrackHitX);
    tree->SetBranchAddress("recoTrackHitY", &recoTrackHitY);
    tree->SetBranchAddress("recoTrackHitZ", &recoTrackHitZ);

    // TOF information
    double TOFMass, tofObject;
    tree->SetBranchAddress("TOFMass", &TOFMass);
    tree->SetBranchAddress("tofObject", &tofObject);

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

    TH1D* hMCNumTracksInCylinder0TG = (TH1D*) MCFile->Get("hMCNumTracksInCylinder0TG");
    TH1D* hMCNumTracksInCylinder1TG = (TH1D*) MCFile->Get("hMCNumTracksInCylinder1TG");
    TH1D* hMCNumTracksInCylinder2TG = (TH1D*) MCFile->Get("hMCNumTracksInCylinder2TG");

    TH1D* hMCNumSmallTracksInCylinder0TG = (TH1D*) MCFile->Get("hMCNumSmallTracksInCylinder0TG");
    TH1D* hMCNumSmallTracksInCylinder1TG = (TH1D*) MCFile->Get("hMCNumSmallTracksInCylinder1TG");
    TH1D* hMCNumSmallTracksInCylinder2TG = (TH1D*) MCFile->Get("hMCNumSmallTracksInCylinder2TG");

    TH1D* hMCTGSmallTracks    = (TH1D*) MCFile->Get("hMCTGSmallTracks");
    TH1D* hMCTGSmallTracks0TG = (TH1D*) MCFile->Get("hMCTGSmallTracks0TG");
    TH1D* hMCTGSmallTracks1TG = (TH1D*) MCFile->Get("hMCTGSmallTracks1TG");
    TH1D* hMCTGSmallTracks2TG = (TH1D*) MCFile->Get("hMCTGSmallTracks2TG");

    TH1D* hMCTGUnreconstructedHitsInduction0TG = (TH1D*) MCFile->Get("hMCTGUnreconstructedHitsInduction0TG");
    TH1D* hMCTGUnreconstructedHitsInduction1TG = (TH1D*) MCFile->Get("hMCTGUnreconstructedHitsInduction1TG");
    TH1D* hMCTGUnreconstructedHitsInduction2TG = (TH1D*) MCFile->Get("hMCTGUnreconstructedHitsInduction2TG");

    TH1D* hMCTGUnreconstructedHitsCollection0TG = (TH1D*) MCFile->Get("hMCTGUnreconstructedHitsCollection0TG");
    TH1D* hMCTGUnreconstructedHitsCollection1TG = (TH1D*) MCFile->Get("hMCTGUnreconstructedHitsCollection1TG");
    TH1D* hMCTGUnreconstructedHitsCollection2TG = (TH1D*) MCFile->Get("hMCTGUnreconstructedHitsCollection2TG");

    TH1D* hMCTGNumClustersInduction0TG = (TH1D*) MCFile->Get("hMCTGNumClustersInduction0TG");
    TH1D* hMCTGNumClustersInduction1TG = (TH1D*) MCFile->Get("hMCTGNumClustersInduction1TG");
    TH1D* hMCTGNumClustersInduction2TG = (TH1D*) MCFile->Get("hMCTGNumClustersInduction2TG");

    TH1D* hMCTGClusterSizesInduction0TG = (TH1D*) MCFile->Get("hMCTGClusterSizesInduction0TG");
    TH1D* hMCTGClusterSizesInduction1TG = (TH1D*) MCFile->Get("hMCTGClusterSizesInduction1TG");
    TH1D* hMCTGClusterSizesInduction2TG = (TH1D*) MCFile->Get("hMCTGClusterSizesInduction2TG");

    TH1D* hMCTGNumClustersCollection0TG = (TH1D*) MCFile->Get("hMCTGNumClustersCollection0TG");
    TH1D* hMCTGNumClustersCollection1TG = (TH1D*) MCFile->Get("hMCTGNumClustersCollection1TG");
    TH1D* hMCTGNumClustersCollection2TG = (TH1D*) MCFile->Get("hMCTGNumClustersCollection2TG");

    TH1D* hMCTGClusterSizesCollection0TG = (TH1D*) MCFile->Get("hMCTGClusterSizesCollection0TG");
    TH1D* hMCTGClusterSizesCollection1TG = (TH1D*) MCFile->Get("hMCTGClusterSizesCollection1TG");
    TH1D* hMCTGClusterSizesCollection2TG = (TH1D*) MCFile->Get("hMCTGClusterSizesCollection2TG");

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

    // With primary TG
    TH1D* hNumTracksInCylinder = new TH1D("hNumTracksInCylinder", "NumTracksInCylinder", 10, 0, 10);
    TH1D* hTGTrackLengths      = new TH1D("hTGTrackLengths", "TGTrackLengths", 25, 0, 50);
    TH1D* hTGSmallTracks       = new TH1D("hTGSmallTracks", "TGSmallTracks", 10, 0, 10);

    TH1D* hNumTracksInCylinder0TG      = new TH1D("hNumTracksInCylinder0TG", "NumTracksInCylinder0TG", 10, 0, 10);
    TH1D* hNumSmallTracksInCylinder0TG = new TH1D("hNumSmallTracksInCylinder0TG", "NumSmallTracksInCylinder0TG", 10, 0, 10);
    TH1D* hTGSmallTracks0TG            = new TH1D("hTGSmallTracks0TG", "TGSmallTracks0TG", 10, 0, 10);

    TH1D* hNumTracksInCylinder1TG      = new TH1D("hNumTracksInCylinder1TG", "NumTracksInCylinder1TG", 10, 0, 10);
    TH1D* hNumSmallTracksInCylinder1TG = new TH1D("hNumSmallTracksInCylinder1TG", "NumSmallTracksInCylinder1TG", 10, 0, 10);
    TH1D* hTGSmallTracks1TG            = new TH1D("hTGSmallTracks1TG", "TGSmallTracks1TG", 10, 0, 10);

    TH1D* hNumTracksInCylinder2TG      = new TH1D("hNumTracksInCylinder2TG", "NumTracksInCylinder2TG", 10, 0, 10);
    TH1D* hNumSmallTracksInCylinder2TG = new TH1D("hNumSmallTracksInCylinder2TG", "NumSmallTracksInCylinder2TG", 10, 0, 10);
    TH1D* hTGSmallTracks2TG            = new TH1D("hTGSmallTracks2TG", "TGSmallTracks2TG", 10, 0, 10);

    TH1D* hTGUnreconstructedHitsInduction0TG = new TH1D("hTGUnreconstructedHitsInduction0TG", "UnreconstructedInductionHits0TG;;", 30, 0, 30);
    TH1D* hTGUnreconstructedHitsInduction1TG = new TH1D("hTGUnreconstructedHitsInduction1TG", "UnreconstructedInductionHits1TG;;", 30, 0, 30);
    TH1D* hTGUnreconstructedHitsInduction2TG = new TH1D("hTGUnreconstructedHitsInduction2TG", "UnreconstructedInductionHits2TG;;", 30, 0, 30);

    TH1D* hTGUnreconstructedHitsCollection0TG = new TH1D("hTGUnreconstructedHitsCollection0TG", "UnreconstructedCollectionHits0TG;;", 30, 0, 30);
    TH1D* hTGUnreconstructedHitsCollection1TG = new TH1D("hTGUnreconstructedHitsCollection1TG", "UnreconstructedCollectionHits1TG;;", 30, 0, 30);
    TH1D* hTGUnreconstructedHitsCollection2TG = new TH1D("hTGUnreconstructedHitsCollection2TG", "UnreconstructedCollectionHits2TG;;", 30, 0, 30);

    TH1D* hTGNumClustersInduction0TG = new TH1D("hTGNumClustersInduction0TG", "NumClustersInduction0TG;;", 10, 0, 10);
    TH1D* hTGNumClustersInduction1TG = new TH1D("hTGNumClustersInduction1TG", "NumClustersInduction1TG;;", 10, 0, 10);
    TH1D* hTGNumClustersInduction2TG = new TH1D("hTGNumClustersInduction2TG", "NumClustersInduction2TG;;", 10, 0, 10);

    TH1D* hTGNumClustersCollection0TG = new TH1D("hTGNumClustersCollection0TG", "NumClustersCollection0TG;;", 10, 0, 10);
    TH1D* hTGNumClustersCollection1TG = new TH1D("hTGNumClustersCollection1TG", "NumClustersCollection1TG;;", 10, 0, 10);
    TH1D* hTGNumClustersCollection2TG = new TH1D("hTGNumClustersCollection2TG", "NumClustersCollection2TG;;", 10, 0, 10);

    TH1D* hTGClusterSizesInduction0TG = new TH1D("hTGClusterSizesInduction0TG", "ClusterSizesInduction0TG;;", 15, 0, 30);
    TH1D* hTGClusterSizesInduction1TG = new TH1D("hTGClusterSizesInduction1TG", "ClusterSizesInduction1TG;;", 15, 0, 30);
    TH1D* hTGClusterSizesInduction2TG = new TH1D("hTGClusterSizesInduction2TG", "ClusterSizesInduction2TG;;", 15, 0, 30);

    TH1D* hTGClusterSizesCollection0TG = new TH1D("hTGClusterSizesCollection0TG", "ClusterSizesCollection0TG;;", 15, 0, 30);
    TH1D* hTGClusterSizesCollection1TG = new TH1D("hTGClusterSizesCollection1TG", "ClusterSizesCollection1TG;;", 15, 0, 30);
    TH1D* hTGClusterSizesCollection2TG = new TH1D("hTGClusterSizesCollection2TG", "ClusterSizesCollection2TG;;", 15, 0, 30);

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
    TH1D* hRadDistWC4 = new TH1D("hRadDistWC4", "RadDistWC4", 80, 0.0, 40.);

    TH2D* hMidPlaneCoinc = new TH2D(
        "hMidPlaneCoinc", "MidPlaneCoinc;up_{x} - down_{x};up_{y} - down_{y}",
        100, -20., 20.,
        100, -5., 5.
    );
    TH1D* hRadDistMidPlane = new TH1D("hRadDistMidPlane", "RadDistMidPlane", 60, 0.0, 30.);

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

    // Keep track of event counts for stat estimations
    int eventCount0TG = 0;
    int eventCount1TG = 0;
    int eventCount2TG = 0;

    int numValidEvents = 0;
    for (Int_t i = 0; i < NumEntries; ++i) {
        tree->GetEntry(i);

        // Make script go faster
        if (i > 100000) break;

        // Fill for all events
        hTOFMass->Fill(std::abs(TOFMass));
        hNumWC2TPCMatch->Fill(WC2TPCsize);
        hNumWCHits->Fill(wcNumHits);
        hTOF->Fill(tofObject);
        hTOFVsTOFMass->Fill(tofObject, std::abs(TOFMass));

        /////////////////////////////////////////
        // Wire-chamber and other initial cuts //
        /////////////////////////////////////////

        if (wcNumHits != 8) continue;

        // Project downstream to midplane @ -437.97 (without angular corrections)
        std::vector<double> midUp = projToZ(*wcHit0, *wcHit1, -437.97);
        // Use this point and WC3 to project up to WC4
        std::vector<double> projDown = projToZ(midUp, *wcHit2, -95.0);
        // Requires some corrections because magnets are not the same
        projDown[0] -= tan(1.32 * TMath::Pi() / 180.0) * (-95.0 - -437.97);

        // Compare x and y coordinate in projection and real hit for WC4
        hProjVsRealWC4->Fill(projDown[0] - wcHit3->at(0), projDown[1] - wcHit3->at(1));
        double radDistWC4 = TMath::Sqrt(pow(projDown[0] - wcHit3->at(0), 2.) + pow(projDown[1] - wcHit3->at(1), 2.));
        hRadDistWC4->Fill(radDistWC4);

        // if (radDistWC4 > 8.0) continue; // TODO: figure out number

        // Project upstream to midplane @ -437.97
        std::vector<double> midDown = projToZ(*wcHit3, *wcHit2, -437.97);
        midDown[0] -= tan(1.32 * TMath::Pi() / 180.0) * (-339.57 - -437.97);

        hMidPlaneCoinc->Fill(midUp[0] - midDown[0], midUp[1] - midDown[1]);
        double midPlaneDist = TMath::Sqrt(pow(midUp[0] - midDown[0], 2) + pow(midUp[1] - midDown[1], 2));
        hRadDistMidPlane->Fill(midPlaneDist);

        // if (midPlaneDist > 3.0) continue;

        // Check projected tracks go through all apertures
        bool Magnet1ApertureCheck = CheckUpstreamMagnetAperture(*wcHit0, *wcHit1);
        bool Magnet2ApertureCheck = CheckDownstreamMagnetAperture(*wcHit3, *wcHit2);
        bool DSColApertureCheck   = CheckDownstreamCollimatorAperture(*wcHit3, *wcHit2);

        if (!Magnet1ApertureCheck || !Magnet2ApertureCheck || !DSColApertureCheck) continue;

        // Candidate mass cut, keep pions, muons and electrons
        if (std::abs(TOFMass) > PI_MU_EL_MASS_CUTOFF) continue;

        // If no track matched to wire-chamber, skip
        if (WC2TPCtrkID == -99999) continue;

        //////////////////////////
        // Analyze valid events //
        //////////////////////////

        numValidEvents++;
        hPrimaryTrackPosition->Fill(WC2TPCPrimaryBeginX, WC2TPCPrimaryBeginY);

        // Check if WC2TPC is through-going
        bool isPrimaryTG = !isWithinReducedVolume(WC2TPCPrimaryEndX, WC2TPCPrimaryEndY, WC2TPCPrimaryEndZ);

        // Loop over reconstructed tracks
        int numSmallTracks = 0; 
        int numTracksNearVertex = 0; 
        int smallTracksTPCStart = 0; 
        int numTGTracks = 0;
        int numTracksInCylinder = 0;
        int numSmallTracksInCylinder = 0;

        for (size_t trk_idx = 0; trk_idx < recoBeginX->size(); ++trk_idx) {
            if (recoTrkID->at(trk_idx) == WC2TPCtrkID) continue;

            double distanceFromStart = distance(
                recoBeginX->at(trk_idx), WC2TPCPrimaryEndX, 
                recoBeginY->at(trk_idx), WC2TPCPrimaryEndY,
                recoBeginZ->at(trk_idx), WC2TPCPrimaryEndZ
            );
            double distanceFromEnd = distance(
                recoEndX->at(trk_idx), WC2TPCPrimaryEndX, 
                recoEndY->at(trk_idx), WC2TPCPrimaryEndY,
                recoEndZ->at(trk_idx), WC2TPCPrimaryEndZ
            );

            double trackLength = sqrt(
                pow(recoEndX->at(trk_idx) - recoBeginX->at(trk_idx), 2) +
                pow(recoEndY->at(trk_idx) - recoBeginY->at(trk_idx), 2) +
                pow(recoEndZ->at(trk_idx) - recoBeginZ->at(trk_idx), 2)
            );

            // Copy WC2TPCLocations
            std::vector<double>* wcX = new std::vector<double>(*WC2TPCLocationsX);
            std::vector<double>* wcY = new std::vector<double>(*WC2TPCLocationsY);
            std::vector<double>* wcZ = new std::vector<double>(*WC2TPCLocationsZ);

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
                recoBeginX->at(trk_idx), recoBeginY->at(trk_idx), recoBeginZ->at(trk_idx),
                CYLINDER_RADIUS
            );
            bool endInCylinder = IsPointInsideTrackCylinder(
                wcX, wcY, wcZ,
                recoEndX->at(trk_idx), recoEndY->at(trk_idx), recoEndZ->at(trk_idx),
                CYLINDER_RADIUS
            );
            if (startInCylinder && endInCylinder) {
                numTracksInCylinder++;
                if (trackLength < CYLINDER_SMALL_TRACK) numSmallTracksInCylinder++;
            }

            if (isPrimaryTG) hTGTrackLengths->Fill(trackLength);

            if (
                recoEndZ->at(trk_idx) < 30.0 && 
                recoBeginZ->at(trk_idx) < 30.0
            ) {
                if (trackLength < SMALL_TRACK_LENGTH_CHEX) smallTracksTPCStart++;
            }

            if (
                !isWithinReducedVolume(recoBeginX->at(trk_idx), recoBeginY->at(trk_idx), recoBeginZ->at(trk_idx)) &&
                !isWithinReducedVolume(recoEndX->at(trk_idx), recoEndY->at(trk_idx), recoEndZ->at(trk_idx))
            ) {
                numTGTracks++;

                if (isPrimaryTG) {
                    // Track information about background tracks
                    bool isBeginTrue = recoBeginZ->at(trk_idx) < recoEndZ->at(trk_idx);

                    const double xs = isBeginTrue ? recoBeginX->at(trk_idx) : recoEndX->at(trk_idx);
                    const double ys = isBeginTrue ? recoBeginY->at(trk_idx) : recoEndY->at(trk_idx);
                    const double zs = isBeginTrue ? recoBeginZ->at(trk_idx) : recoEndZ->at(trk_idx);

                    const double xe = isBeginTrue ? recoEndX->at(trk_idx)   : recoBeginX->at(trk_idx);
                    const double ye = isBeginTrue ? recoEndY->at(trk_idx)   : recoBeginY->at(trk_idx);
                    const double ze = isBeginTrue ? recoEndZ->at(trk_idx)   : recoBeginZ->at(trk_idx);

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

        hNumTGTracks->Fill(numTGTracks);

        if (!isPrimaryTG) hTracksNearVertex->Fill(numTracksNearVertex);

        // Add to histogram of small tracks if primary is throughgoing
        if (isPrimaryTG) {
            hTGSmallTracks->Fill(numSmallTracks);
            hSmallVsTGTracks->Fill(numSmallTracks, numTGTracks);
            hNumTracksInCylinder->Fill(numTracksInCylinder);

            // Scan over small track length thresholds and fill 2D histogram
            for (int threshBin = 1; threshBin <= hTGNumSmallTracksVsThresh->GetNbinsX(); ++threshBin) {
                double threshold = hTGNumSmallTracksVsThresh->GetXaxis()->GetBinCenter(threshBin);
                int nSmallTracks = 0;
                for (size_t trk_idx = 0; trk_idx < recoBeginX->size(); ++trk_idx) {
                    if (recoTrkID->at(trk_idx) == WC2TPCtrkID) continue;
                    double trackLength = sqrt(
                        pow(recoEndX->at(trk_idx) - recoBeginX->at(trk_idx), 2) +
                        pow(recoEndY->at(trk_idx) - recoBeginY->at(trk_idx), 2) +
                        pow(recoEndZ->at(trk_idx) - recoBeginZ->at(trk_idx), 2)
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
        std::unordered_set<int> hitsInTracks(hitRecoAsTrackKey->begin(), hitRecoAsTrackKey->end());

        // Reconstruct hit clusters
        std::vector<int> candidateHits;
        for (size_t iHit = 0; iHit < fHitKey->size(); ++iHit) {
            // Skip hits already in tracks
            if (hitsInTracks.count(iHit) > 0) continue;

            double hitX     = fHitX->at(iHit);
            double hitW     = fHitW->at(iHit);
            int    hitPlane = fHitPlane->at(iHit);

            if (isHitNearPrimary(
                hitWC2TPCKey,
                fHitX,
                fHitW,
                fHitPlane,
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
            int   thisHitPlane     = fHitPlane->at(thisHitKey);
            float thisHitW         = fHitW->at(thisHitKey);
            float thisHitX         = fHitX->at(thisHitKey);
            float thisHitCharge    = fHitCharge->at(thisHitKey);
            float thisHitChargeCol = fHitChargeCol->at(thisHitKey);
            
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
            
            for (int iAllHit = 0; iAllHit < fHitKey->size(); ++iAllHit) {
                // Skip already used hits, and those reconstructed in tracks
                if (usedHits.count(iAllHit) || hitsInTracks.count(iAllHit)) continue;

                // Clusters have to be in same plane
                if (fHitPlane->at(iAllHit) != thisHitPlane) continue;

                float internalHitW  = fHitW->at(iAllHit);
                float internalHitX  = fHitX->at(iAllHit);
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
                    clusterX.push_back(fHitX->at(iAllHit));
                    clusterW.push_back(fHitW->at(iAllHit));
                    clusterCharge.push_back(fHitCharge->at(iAllHit));
                    clusterChargeCol.push_back(fHitChargeCol->at(iAllHit));
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

        for (size_t i = 0; i < hitWC2TPCKey->size(); ++i) {
            // Only care about hits inside the reduced volume
            auto [i_hit, j_hit] = find_unique_position(recoTrackHitIndices, hitWC2TPCKey->at(i));
            if (!isWithinReducedVolume(
                recoTrackHitX->at(i_hit)[j_hit],
                recoTrackHitY->at(i_hit)[j_hit],
                recoTrackHitZ->at(i_hit)[j_hit]
            )) continue;

            if (fHitPlane->at(hitWC2TPCKey->at(i)) == 0) hitWC2TPCKeyInduction.push_back(hitWC2TPCKey->at(i));
            else if (fHitPlane->at(hitWC2TPCKey->at(i)) == 1) hitWC2TPCKeyCollection.push_back(hitWC2TPCKey->at(i));

            hTimePrimaryHits->Fill(fHitT->at(hitWC2TPCKey->at(i)));
        }

        // First, get a random point in the primary track
        std::vector<int> randomInduction; std::vector<int> randomCollection;

        if (hitWC2TPCKeyInduction.size() != 0) {
            randomInduction.push_back(0);

            for (int i = 1; i < (int) hitWC2TPCKeyInduction.size(); ++i) {
                double dX = std::abs(fHitX->at(hitWC2TPCKeyInduction[randomInduction[randomInduction.size() - 1]]) - fHitX->at(hitWC2TPCKeyInduction[i]));
                double dW = std::abs(fHitW->at(hitWC2TPCKeyInduction[randomInduction[randomInduction.size() - 1]]) - fHitW->at(hitWC2TPCKeyInduction[i]));

                if (dX > MAX_IN_CLUSTER_X_SEPARATION && dW > MAX_IN_CLUSTER_W_SEPARATION) randomInduction.push_back(i);
            }
        } 
        if (hitWC2TPCKeyCollection.size() != 0) {
            randomCollection.push_back(0);

            for (int i = 1; i < (int) hitWC2TPCKeyCollection.size(); ++i) {
                double dX = std::abs(fHitX->at(hitWC2TPCKeyCollection[randomCollection[randomCollection.size() - 1]]) - fHitX->at(hitWC2TPCKeyCollection[i]));
                double dW = std::abs(fHitW->at(hitWC2TPCKeyCollection[randomCollection[randomCollection.size() - 1]]) - fHitW->at(hitWC2TPCKeyCollection[i]));

                if (dX > MAX_IN_CLUSTER_X_SEPARATION && dW > MAX_IN_CLUSTER_W_SEPARATION) randomCollection.push_back(i);
            }
        }

        int numUnrecoHitsInduction  = 0;
        int numUnrecoHitsCollection = 0;

        std::vector<int> candidateRandomHitsInduction;
        std::vector<int> candidateRandomHitsCollection;

        for (size_t iHit = 0; iHit < fHitKey->size(); ++iHit) {
            // Skip hits already in tracks
            if (hitsInTracks.count(iHit) > 0) continue;

            double hitX     = fHitX->at(iHit);
            double hitW     = fHitW->at(iHit);
            int    hitPlane = fHitPlane->at(iHit);

            if (hitPlane == 0 && hitWC2TPCKeyInduction.size() > 0) {
                for (int idx : randomInduction) {
                    double dW = (hitW - fHitW->at(hitWC2TPCKeyInduction[idx]));
                    double dX = (hitX - fHitX->at(hitWC2TPCKeyInduction[idx]));
                    double d  = std::sqrt(std::pow(dW, 2) + std::pow(dX, 2));
                    if (d < DISTANCE_TO_PRIMARY_THRESHOLD) {
                        numUnrecoHitsInduction++;
                        candidateRandomHitsInduction.push_back(iHit);
                    }
                }
            } else if (hitPlane == 1 && hitWC2TPCKeyCollection.size() > 0) {
                for (int idx : randomCollection) {
                    double dW = (hitW - fHitW->at(hitWC2TPCKeyCollection[idx]));
                    double dX = (hitX - fHitX->at(hitWC2TPCKeyCollection[idx]));
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
                    if (numTGTracks == 0 && isPrimaryTG) hTGClusterSizesInduction0TG->Fill(thisCluster.clusterSize);
                    else if (numTGTracks == 1 && isPrimaryTG) hTGClusterSizesInduction1TG->Fill(thisCluster.clusterSize);
                    else if (numTGTracks == 2 && isPrimaryTG) hTGClusterSizesInduction2TG->Fill(thisCluster.clusterSize);
                    break;
                } else if (std::find(candidateRandomHitsCollection.begin(), candidateRandomHitsCollection.end(), thisCluster.hitKeys[iHit]) != candidateRandomHitsCollection.end()) {
                    numClustersCollection++;
                    if (numTGTracks == 0 && isPrimaryTG) hTGClusterSizesCollection0TG->Fill(thisCluster.clusterSize);
                    else if (numTGTracks == 1 && isPrimaryTG) hTGClusterSizesCollection1TG->Fill(thisCluster.clusterSize);
                    else if (numTGTracks == 2 && isPrimaryTG) hTGClusterSizesCollection2TG->Fill(thisCluster.clusterSize);
                    break;
                }
            }
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

        //////////////////
        // TG track cut //
        //////////////////

        // Grab data about number of events with each cutoff
        if (numTGTracks <= 0) eventCount0TG++;
        if (numTGTracks <= 1) eventCount1TG++;
        if (numTGTracks <= 2) eventCount2TG++;

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

        for (size_t trk_idx = 0; trk_idx < recoBeginX->size(); ++trk_idx) {
            if (recoTrkID->at(trk_idx) == WC2TPCtrkID) continue;

            double distanceFromStart = distance(
                recoBeginX->at(trk_idx), WC2TPCPrimaryEndX, 
                recoBeginY->at(trk_idx), WC2TPCPrimaryEndY,
                recoBeginZ->at(trk_idx), WC2TPCPrimaryEndZ
            );
            double distanceFromEnd = distance(
                recoEndX->at(trk_idx), WC2TPCPrimaryEndX, 
                recoEndY->at(trk_idx), WC2TPCPrimaryEndY,
                recoEndZ->at(trk_idx), WC2TPCPrimaryEndZ
            );

            double trackLength = sqrt(
                pow(recoEndX->at(trk_idx) - recoBeginX->at(trk_idx), 2) +
                pow(recoEndY->at(trk_idx) - recoBeginY->at(trk_idx), 2) +
                pow(recoEndZ->at(trk_idx) - recoBeginZ->at(trk_idx), 2)
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

    double scaling0TG = hNumTracksInCylinder0TG->Integral() / hMCNumTracksInCylinder0TG->Integral();
    double scaling1TG = hNumTracksInCylinder1TG->Integral() / hMCNumTracksInCylinder1TG->Integral();
    double scaling2TG = hNumTracksInCylinder2TG->Integral() / hMCNumTracksInCylinder2TG->Integral();
    double scalingKin = hNumCandidateProtons->Integral() / hMCNumCandidateProtons->Integral();

    // Scale MC histograms to data histograms using event counts
    hMCNumWC2TPCMatch->Scale(hNumWC2TPCMatch->Integral() / hMCNumWC2TPCMatch->Integral());

    // all events
    hMCNumTGTracks->Scale(scaling);

    // primary not TG
    hMCTracksNearVertex->Scale(scalingNoTG);
    hMCTrackLengthsNearVertex->Scale(scalingNoTG);

    // primary TG
    hMCTGSmallTracks->Scale(scalingTG);
    hMCTGTrackLengths->Scale(scalingTG);

    hMCTGSmallTracks0TG->Scale(scaling0TG);
    hMCTGSmallTracks1TG->Scale(scaling1TG);
    hMCTGSmallTracks2TG->Scale(scaling2TG);

    hMCNumTracksInCylinder0TG->Scale(scaling0TG);
    hMCNumTracksInCylinder1TG->Scale(scaling1TG);
    hMCNumTracksInCylinder2TG->Scale(scaling2TG);

    hMCNumSmallTracksInCylinder0TG->Scale(scaling0TG);
    hMCNumSmallTracksInCylinder1TG->Scale(scaling1TG);
    hMCNumSmallTracksInCylinder2TG->Scale(scaling2TG);

    hMCTGUnreconstructedHitsInduction0TG->Scale(scaling0TG);
    hMCTGUnreconstructedHitsInduction1TG->Scale(scaling1TG);
    hMCTGUnreconstructedHitsInduction2TG->Scale(scaling2TG);

    hMCTGUnreconstructedHitsCollection0TG->Scale(scaling0TG);
    hMCTGUnreconstructedHitsCollection1TG->Scale(scaling1TG);
    hMCTGUnreconstructedHitsCollection2TG->Scale(scaling2TG);

    hMCTGClusterSizesInduction0TG->Scale(scaling0TG);
    hMCTGClusterSizesInduction1TG->Scale(scaling1TG);
    hMCTGClusterSizesInduction2TG->Scale(scaling2TG);

    hMCTGClusterSizesCollection0TG->Scale(scaling0TG);
    hMCTGClusterSizesCollection1TG->Scale(scaling1TG);
    hMCTGClusterSizesCollection2TG->Scale(scaling2TG);

    hMCTGNumClustersCollection0TG->Scale(scaling0TG);
    hMCTGNumClustersCollection1TG->Scale(scaling1TG);
    hMCTGNumClustersCollection2TG->Scale(scaling2TG);

    hMCTGNumClustersInduction0TG->Scale(scaling0TG);
    hMCTGNumClustersInduction1TG->Scale(scaling1TG);
    hMCTGNumClustersInduction2TG->Scale(scaling2TG);

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
        {hNumTracksInCylinder0TG, hMCNumTracksInCylinder0TG},
        {hNumTracksInCylinder1TG, hMCNumTracksInCylinder1TG},
        {hNumTracksInCylinder2TG, hMCNumTracksInCylinder2TG},

        {hNumSmallTracksInCylinder0TG, hMCNumSmallTracksInCylinder0TG},
        {hNumSmallTracksInCylinder1TG, hMCNumSmallTracksInCylinder1TG},
        {hNumSmallTracksInCylinder2TG, hMCNumSmallTracksInCylinder2TG},

        // Hits in primary
        {hTimePrimaryHits},

        // Comparisons with MC
        {hTracksNearVertex, hMCTracksNearVertex},
        {hTrackLengthsNearVertex, hMCTrackLengthsNearVertex},
        {hNumTGTracks, hMCNumTGTracks},
        {hTGTrackLengths, hMCTGTrackLengths},

        // Comparing TG threshold
        {hTGSmallTracks0TG, hMCTGSmallTracks0TG},
        {hTGSmallTracks1TG, hMCTGSmallTracks1TG},
        {hTGSmallTracks2TG, hMCTGSmallTracks2TG},

        // Unreconstructed hits
        {hTGUnreconstructedHitsInduction0TG, hMCTGUnreconstructedHitsInduction0TG},
        {hTGUnreconstructedHitsInduction1TG, hMCTGUnreconstructedHitsInduction1TG},
        {hTGUnreconstructedHitsInduction2TG, hMCTGUnreconstructedHitsInduction2TG},

        {hTGUnreconstructedHitsCollection0TG, hMCTGUnreconstructedHitsCollection0TG},
        {hTGUnreconstructedHitsCollection1TG, hMCTGUnreconstructedHitsCollection1TG},
        {hTGUnreconstructedHitsCollection2TG, hMCTGUnreconstructedHitsCollection2TG},

        {hTGClusterSizesCollection0TG, hMCTGClusterSizesCollection0TG},
        {hTGClusterSizesCollection1TG, hMCTGClusterSizesCollection1TG},
        {hTGClusterSizesCollection2TG, hMCTGClusterSizesCollection2TG},

        {hTGClusterSizesInduction0TG, hMCTGClusterSizesInduction0TG},
        {hTGClusterSizesInduction1TG, hMCTGClusterSizesInduction1TG},
        {hTGClusterSizesInduction2TG, hMCTGClusterSizesInduction2TG},

        {hTGNumClustersCollection0TG, hMCTGNumClustersCollection0TG},
        {hTGNumClustersCollection1TG, hMCTGNumClustersCollection1TG},
        {hTGNumClustersCollection2TG, hMCTGNumClustersCollection2TG},

        {hTGNumClustersInduction0TG, hMCTGNumClustersInduction0TG},
        {hTGNumClustersInduction1TG, hMCTGNumClustersInduction1TG},
        {hTGNumClustersInduction2TG, hMCTGNumClustersInduction2TG},

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
        "Cylinder/NumTracksInCylinder0TG",
        "Cylinder/NumTracksInCylinder1TG",
        "Cylinder/NumTracksInCylinder2TG",

        "Cylinder/NumSmallTracksInCylinder0TG",
        "Cylinder/NumSmallTracksInCylinder1TG",
        "Cylinder/NumSmallTracksInCylinder2TG",

        // Hits in primary
        "Primary/HitTime",

        // Comparisons with MC
        "NearVertex/TracksNearVertex",
        "NearVertex/TrackLengthsNearVertex",
        "TGTracks/NumTGTracks",
        "Primary/TGTrackLengths",

        // Small tracks when primary is TG
        "Primary/TGSmallTracks0TG",
        "Primary/TGSmallTracks1TG",
        "Primary/TGSmallTracks2TG",

        // Unreconstructed hits
        "UnrecoHits/UnreconstructedInductionHits0TG",
        "UnrecoHits/UnreconstructedInductionHits1TG",
        "UnrecoHits/UnreconstructedInductionHits2TG",

        "UnrecoHits/UnreconstructedCollectionHits0TG",
        "UnrecoHits/UnreconstructedCollectionHits1TG",
        "UnrecoHits/UnreconstructedCollectionHits2TG",

        "UnrecoHits/ClusterSizesCollection0TG",
        "UnrecoHits/ClusterSizesCollection1TG",
        "UnrecoHits/ClusterSizesCollection2TG",

        "UnrecoHits/ClusterSizesInduction0TG",
        "UnrecoHits/ClusterSizesInduction1TG",
        "UnrecoHits/ClusterSizesInduction2TG",

        "UnrecoHits/NumClustersCollection0TG",
        "UnrecoHits/NumClustersCollection1TG",
        "UnrecoHits/NumClustersCollection2TG",

        "UnrecoHits/NumClustersInduction0TG",
        "UnrecoHits/NumClustersInduction1TG",
        "UnrecoHits/NumClustersInduction2TG",

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

        // Unreconstructed hits
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

        // Secondary particle kinematics
        {true, false},
        {true, false}
    };

    /////////////////////////////////
    // Fractional difference plots //
    /////////////////////////////////
    
    std::vector<std::pair<TH1*,TH1*>> PlotGroupsFracDiff = {
        // Primary
        {hTGSmallTracks0TG, hMCTGSmallTracks0TG},
        {hTGSmallTracks1TG, hMCTGSmallTracks1TG},
        {hTGSmallTracks2TG, hMCTGSmallTracks2TG},
        {hTGTrackLengths, hMCTGTrackLengths},

        // Cylinder
        {hNumTracksInCylinder0TG, hMCNumTracksInCylinder0TG},
        {hNumTracksInCylinder1TG, hMCNumTracksInCylinder1TG},
        {hNumTracksInCylinder2TG, hMCNumTracksInCylinder2TG},

        {hNumSmallTracksInCylinder0TG, hMCNumSmallTracksInCylinder0TG},
        {hNumSmallTracksInCylinder1TG, hMCNumSmallTracksInCylinder1TG},
        {hNumSmallTracksInCylinder2TG, hMCNumSmallTracksInCylinder2TG},

        // Near vertex
        {hTracksNearVertex, hMCTracksNearVertex},
        {hTrackLengthsNearVertex, hMCTrackLengthsNearVertex},

        // TG tracks
        {hNumTGTracks, hMCNumTGTracks},

        // Unreconstruced hits
        {hTGUnreconstructedHitsInduction0TG, hMCTGUnreconstructedHitsInduction0TG},
        {hTGUnreconstructedHitsInduction1TG, hMCTGUnreconstructedHitsInduction1TG},
        {hTGUnreconstructedHitsInduction2TG, hMCTGUnreconstructedHitsInduction2TG},

        {hTGUnreconstructedHitsCollection0TG, hMCTGUnreconstructedHitsCollection0TG},
        {hTGUnreconstructedHitsCollection1TG, hMCTGUnreconstructedHitsCollection1TG},
        {hTGUnreconstructedHitsCollection2TG, hMCTGUnreconstructedHitsCollection2TG},

        {hTGClusterSizesCollection0TG, hMCTGClusterSizesCollection0TG},
        {hTGClusterSizesCollection1TG, hMCTGClusterSizesCollection1TG},
        {hTGClusterSizesCollection2TG, hMCTGClusterSizesCollection2TG},

        {hTGClusterSizesInduction0TG, hMCTGClusterSizesInduction0TG},
        {hTGClusterSizesInduction1TG, hMCTGClusterSizesInduction1TG},
        {hTGClusterSizesInduction2TG, hMCTGClusterSizesInduction2TG},

        {hTGNumClustersCollection0TG, hMCTGNumClustersCollection0TG},
        {hTGNumClustersCollection1TG, hMCTGNumClustersCollection1TG},
        {hTGNumClustersCollection2TG, hMCTGNumClustersCollection2TG},

        {hTGNumClustersInduction0TG, hMCTGNumClustersInduction0TG},
        {hTGNumClustersInduction1TG, hMCTGNumClustersInduction1TG},
        {hTGNumClustersInduction2TG, hMCTGNumClustersInduction2TG},

        // Secondary particle kinematics
        {hNumCandidateProtons, hMCNumCandidateProtons},
        {hLengthCandidateProtons, hMCLengthCandidateProtons}
    };

    std::vector<std::string> FracDiffDirectory = {
        "Primary",
        "Primary",
        "Primary",
        "Primary",

        "Cylinder",
        "Cylinder",
        "Cylinder",

        "Cylinder",
        "Cylinder",
        "Cylinder",

        "NearVertex",
        "NearVertex",

        "TGTracks",

        "UnrecoHits",
        "UnrecoHits",
        "UnrecoHits",

        "UnrecoHits",
        "UnrecoHits",
        "UnrecoHits",

        "UnrecoHits",
        "UnrecoHits",
        "UnrecoHits",

        "UnrecoHits",
        "UnrecoHits",
        "UnrecoHits",

        "UnrecoHits",
        "UnrecoHits",
        "UnrecoHits",

        "UnrecoHits",
        "UnrecoHits",
        "UnrecoHits",

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