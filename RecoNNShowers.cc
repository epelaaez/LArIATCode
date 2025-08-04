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

//////////////////////
// Global constants //
//////////////////////

std::map<int, std::string> backgroundTypes = {
    {-1, "Not flagged as background"},
    {0, "Abs 0p"},
    {1, "Abs Np"},
    {2, "Primary muon"},
    {3, "Primary electron"},
    {4, "Other primary"},
    {5, "Outside reduced volume"},
    {6, "Inelastic scattering"},
    {7, "Charge exchange"},
    {8, "Double charge exchange"},
    {9, "Capture at rest"},
    {10, "Decay"},
    {11, "Other"},
    {12, "Elastic scattering"}
};

double SMALL_TRACK_LENGTH = 10.;
double LARGE_TRACK_LENGTH = 30.;
double VERTEX_RADIUS      = 5.;
double FAR_TRACK_DISTANCE = 5.;
double SHOWER_PROB_CUT    = 0.6;

double PION_CHI2_PION_VALUE     = 1.125;
double PION_CHI2_PROTON_VALUE   = 1.125;
double PROTON_CHI2_PION_VALUE   = 1.375;
double PROTON_CHI2_PROTON_VALUE = 0.125;

double MEAN_DEDX_THRESHOLD       = 4.0;
int    MEAN_DEDX_NUM_TRAJ_POINTS = 5;

//////////////////////
// Helper functions //
//////////////////////

void printOneDPlots(
    const TString& dir, 
    int fontStyle, 
    double textSize,
    std::vector<std::vector<TH1*>>& groups,
    std::vector<int>& colors, 
    std::vector<std::vector<TString>>& labels, 
    std::vector<TString>& titles, 
    std::vector<TString>& xlabels, 
    std::vector<TString>& ylabels,
    std::vector<bool>& stack
);

void printTwoDPlots(
    const TString& dir, 
    const std::vector<TH2*>& plots, 
    const std::vector<TString>& titles
);

double computeReducedChi2(const TGraph* theory, std::vector<double> xData, std::vector<double> yData,  bool dataReversed, int nPoints, int nOutliersToDiscard = 0, int nTrim = 0);
double meanDEDX(std::vector<double>& trackDEDX, bool isTrackReversed, int pointsToUse);
double distance(double x1, double x2, double y1, double y2, double z1, double z2);
void initializeProtonPoints(TGraph* gProton);
void initializePionPoints(TGraph* gPion);
void initializeMuonNoBraggPoints(TGraph* gMuonTG);

///////////////////
// Main function //
///////////////////

void RecoNNShowers() {
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
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis/figs/RecoNNShowers/";

    // Load NN root file
    TString RootFilePath = "/exp/lariat/app/users/epelaez/files/RecoNNAllEval_histo.root";
    std::unique_ptr<TFile> File(TFile::Open(RootFilePath));
    TDirectory* Directory = (TDirectory*)File->Get("RecoNNAllEval");

    ///////////////////
    // Load branches //
    ///////////////////

    // Load tree and branches
    TTree* tree = (TTree*) Directory->Get<TTree>("RecoNNAllEvalTree");

    // Signal information
    bool isPionAbsorptionSignal; int backgroundType, numVisibleProtons; 
    tree->SetBranchAddress("isPionAbsorptionSignal", &isPionAbsorptionSignal);
    tree->SetBranchAddress("backgroundType", &backgroundType);
    tree->SetBranchAddress("numVisibleProtons", &numVisibleProtons);

    // Shower probability information
    double trackProb, showerProb;
    bool obtainedProbabilities;
    tree->SetBranchAddress("trackProb", &trackProb);
    tree->SetBranchAddress("showerProb", &showerProb);
    tree->SetBranchAddress("obtainedProbabilities", &obtainedProbabilities);

    // WC match information
    int WC2TPCtrkID;
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

    // WC match truth information
    int                       wcMatchPDG;
    std::string*              wcMatchProcess          = new std::string();
    std::vector<int>*         wcMatchDaughtersPDG     = nullptr;
    std::vector<std::string>* wcMatchDaughtersProcess = nullptr;
    tree->SetBranchAddress("wcMatchDaughtersPDG", &wcMatchDaughtersPDG);
    tree->SetBranchAddress("wcMatchDaughtersProcess", &wcMatchDaughtersProcess);
    tree->SetBranchAddress("wcMatchProcess", &wcMatchProcess);
    tree->SetBranchAddress("wcMatchPDG", &wcMatchPDG);

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

    // Reco truth-match information
    std::vector<int>*         matchedIdentity = nullptr;
    std::vector<int>*         matchedTrkID    = nullptr;
    std::vector<double>*      matchedKEnergy  = nullptr;
    std::vector<std::string>* matchedProcess  = nullptr;
    tree->SetBranchAddress("matchedIdentity", &matchedIdentity);
    tree->SetBranchAddress("matchedTrkID", &matchedTrkID);
    tree->SetBranchAddress("matchedKEnergy", &matchedKEnergy);
    tree->SetBranchAddress("matchedProcess", &matchedProcess);

    // Calorimetry information
    std::vector<std::vector<double>>* recoResR = nullptr;
    std::vector<std::vector<double>>* recoDEDX = nullptr;
    tree->SetBranchAddress("recoResR", &recoResR);
    tree->SetBranchAddress("recoDEDX", &recoDEDX);

    // Truth-level information
    int truthPrimaryPDG;
    double truthPrimaryVertexX, truthPrimaryVertexY, truthPrimaryVertexZ;
    double truthPrimaryIncidentKE, truthPrimaryVertexKE;
    std::vector<int>*         truthPrimaryDaughtersPDG     = nullptr;
    std::vector<std::string>* truthPrimaryDaughtersProcess = nullptr;
    std::vector<double>*      truthPrimaryDaughtersKE      = nullptr;
    tree->SetBranchAddress("truthPrimaryPDG", &truthPrimaryPDG);
    tree->SetBranchAddress("truthPrimaryIncidentKE", &truthPrimaryIncidentKE);
    tree->SetBranchAddress("truthPrimaryVertexKE", &truthPrimaryVertexKE);
    tree->SetBranchAddress("truthPrimaryVertexX", &truthPrimaryVertexX);
    tree->SetBranchAddress("truthPrimaryVertexY", &truthPrimaryVertexY);
    tree->SetBranchAddress("truthPrimaryVertexZ", &truthPrimaryVertexZ);
    tree->SetBranchAddress("truthPrimaryDaughtersPDG", &truthPrimaryDaughtersPDG);
    tree->SetBranchAddress("truthPrimaryDaughtersProcess", &truthPrimaryDaughtersProcess);
    tree->SetBranchAddress("truthPrimaryDaughtersKE", &truthPrimaryDaughtersKE);

    // Truth-level interaction information
    bool         interactionInTrajectory;
    std::string* trajectoryInteractionLabel = new std::string();
    double       trajectoryInteractionAngle;
    double       trajectoryInteractionX, trajectoryInteractionY, trajectoryInteractionZ;
    double       trajectoryInteractionKE;
    double       trajectoryInitialMomentumX;
    tree->SetBranchAddress("interactionInTrajectory", &interactionInTrajectory);
    tree->SetBranchAddress("trajectoryInteractionLabel", &trajectoryInteractionLabel);
    tree->SetBranchAddress("trajectoryInteractionAngle", &trajectoryInteractionAngle);
    tree->SetBranchAddress("trajectoryInteractionX", &trajectoryInteractionX);
    tree->SetBranchAddress("trajectoryInteractionY", &trajectoryInteractionY);
    tree->SetBranchAddress("trajectoryInteractionZ", &trajectoryInteractionZ);
    tree->SetBranchAddress("trajectoryInteractionKE", &trajectoryInteractionKE);
    tree->SetBranchAddress("trajectoryInitialMomentumX", &trajectoryInitialMomentumX);

    // Truth-level secondary pion interaction information
    std::vector<int>*         truthSecondaryPionDaughtersPDG     = nullptr;
    std::vector<std::string>* truthSecondaryPionDaughtersProcess = nullptr;
    std::vector<double>*      truthSecondaryPionDaughtersKE      = nullptr;
    double truthScatteringAngle, truthSecondaryVertexX, truthSecondaryVertexY, truthSecondaryVertexZ, truthScatteredPionLength, truthScatteredPionKE;
    tree->SetBranchAddress("truthSecondaryPionDaughtersPDG", &truthSecondaryPionDaughtersPDG);
    tree->SetBranchAddress("truthSecondaryPionDaughtersProcess", &truthSecondaryPionDaughtersProcess);
    tree->SetBranchAddress("truthSecondaryPionDaughtersKE", &truthSecondaryPionDaughtersKE);
    tree->SetBranchAddress("truthScatteringAngle", &truthScatteringAngle);
    tree->SetBranchAddress("truthScatteredPionLength", &truthScatteredPionLength);
    tree->SetBranchAddress("truthScatteredPionKE", &truthScatteredPionKE);
    tree->SetBranchAddress("truthSecondaryVertexX", &truthSecondaryVertexX);
    tree->SetBranchAddress("truthSecondaryVertexY", &truthSecondaryVertexY);
    tree->SetBranchAddress("truthSecondaryVertexZ", &truthSecondaryVertexZ);

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

    ///////////////////////
    // Create histograms //
    ///////////////////////

    // Electron shower probability
    TH1D* hElectronShowerProb = new TH1D("hElectronShowerProb", "hElectronShowerProb;;", 20, 0., 1.);
    TH1D* hPionShowerProb     = new TH1D("hPionShowerProb", "hPionShowerProb;;", 20, 0., 1.);
    TH1D* hMuonShowerProb     = new TH1D("hMuonShowerProb", "hMuonShowerProb;;", 20, 0., 1.);

    TH1D* hPionAbs0pShowerProb          = new TH1D("hPionAbs0pShowerProb", "hPionAbs0pShowerProb;;", 20, 0., 1.);
    TH1D* hPionAbsNpShowerProb          = new TH1D("hPionAbsNpShowerProb", "hPionAbsNpShowerProb;;", 20, 0., 1.);
    TH1D* hPion0pScatterShowerProb      = new TH1D("hPion0pScatterShowerProb", "hPion0pScatterShowerProb;;", 20, 0., 1.);
    TH1D* hPionNpScatterShowerProb      = new TH1D("hPionNpScatterShowerProb", "hPionNpScatterShowerProb;;", 20, 0., 1.);
    TH1D* hPionChargeExchangeShowerProb = new TH1D("hPionChargeExchangeShowerProb", "hPionChargeExchangeShowerProb;;", 20, 0., 1.);
    TH1D* hPionOtherShowerProb          = new TH1D("hPionOtherShowerProb", "hPionOtherShowerProb;;", 20, 0., 1.);

    // Tracks near vertex
    TH1D* hChargeExchangeTracksNVertex = new TH1D("hChargeExchangeTracksNVertex", "hChargeExchangeTracksNVertex;;", 5, 0, 5);
    TH1D* hPionAbs0pTracksNVertex     = new TH1D("hPionAbs0pTracksNVertex", "hPionAbs0pTracksNVertex;;", 5, 0, 5);
    TH1D* hPionAbsNpTracksNVertex     = new TH1D("hPionAbsNpTracksNVertex", "hPionAbsNpTracksNVertex;;", 5, 0, 5);
    TH1D* hPion0pScatterTracksNVertex = new TH1D("hPion0pScatterTracksNVertex", "hPion0pScatterTracksNVertex;;", 5, 0, 5);
    TH1D* hPionNpScatterTracksNVertex = new TH1D("hPionNpScatterTracksNVertex", "hPionNpScatterTracksNVertex;;", 5, 0, 5);
    TH1D* hOtherPionTracksNVertex      = new TH1D("hOtherPionTracksNVertex", "hOtherPionTracksNVertex;;", 5, 0, 5);
    TH1D* hElectronTracksNVertex      = new TH1D("hElectronTracksNVertex", "hElectronTracksNVertex;;", 5, 0, 5);
    TH1D* hMuonTracksNVertex          = new TH1D("hMuonTracksNVertex", "hMuonTracksNVertex;;", 5, 0, 5);

    // Tracks away from vertex
    TH1D* hChargeExchangeTracksFarVertex = new TH1D("hChargeExchangeTracksFarVertex", "hChargeExchangeTracksFarVertex;;", 5, 0, 5);
    TH1D* hPionAbs0pTracksFarVertex     = new TH1D("hPionAbs0pTracksFarVertex", "hPionAbs0pTracksFarVertex;;", 5, 0, 5);
    TH1D* hPionAbsNpTracksFarVertex     = new TH1D("hPionAbsNpTracksFarVertex", "hPionAbsNpTracksFarVertex;;", 5, 0, 5);
    TH1D* hPion0pScatterTracksFarVertex = new TH1D("hPion0pScatterTracksFarVertex", "hPion0pScatterTracksFarVertex;;", 5, 0, 5);
    TH1D* hPionNpScatterTracksFarVertex = new TH1D("hPionNpScatterTracksFarVertex", "hPionNpScatterTracksFarVertex;;", 5, 0, 5);
    TH1D* hOtherPionTracksFarVertex      = new TH1D("hOtherPionTracksFarVertex", "hOtherPionTracksFarVertex;;", 5, 0, 5);
    TH1D* hElectronTracksFarVertex      = new TH1D("hElectronTracksFarVertex", "hElectronTracksFarVertex;;", 5, 0, 5);
    TH1D* hMuonTracksFarVertex          = new TH1D("hMuonTracksFarVertex", "hMuonTracksFarVertex;;", 5, 0, 5);

    // Track length analysis
    TH1D* hChargeExchangeTrackLengths = new TH1D("hChargeExchangeTrackLengths", "hChargeExchangeTrackLengths;;", 40, 0., 90.);
    TH1D* hOtherPionTrackLengths      = new TH1D("hOtherPionTrackLengths", "hOtherPionTrackLengths;;", 40, 0., 90.);
    TH1D* hElectronTrackLengths       = new TH1D("hElectronTrackLengths", "hElectronTrackLengths;;", 40, 0., 90.);
    TH1D* hMuonTrackLengths           = new TH1D("hMuonTrackLengths", "hMuonTrackLengths;;", 40, 0., 90.);

    TH1D* hChargeExchangeSmallTracksCount = new TH1D("hChargeExchangeSmallTracksCount", "hChargeExchangeSmallTracksCount;;", 20, 0, 20);
    TH1D* hOtherPionSmallTracksCount      = new TH1D("hOtherPionSmallTracksCount", "hOtherPionSmallTracksCount;;", 20, 0, 20);
    TH1D* hElectronSmallTracksCount       = new TH1D("hElectronSmallTracksCount", "hElectronSmallTracksCount;;", 20, 0, 20);
    TH1D* hMuonSmallTracksCount           = new TH1D("hMuonSmallTracksCount", "hMuonSmallTracksCount;;", 20, 0, 20);

    TH1D* hChargeExchangeLargeTracksCount = new TH1D("hChargeExchangeLargeTracksCount", "hChargeExchangeLargeTracksCount;;", 20, 0, 20);
    TH1D* hOtherPionLargeTracksCount      = new TH1D("hOtherPionLargeTracksCount", "hOtherPionLargeTracksCount;;", 20, 0, 20);
    TH1D* hElectronLargeTracksCount       = new TH1D("hElectronLargeTracksCount", "hElectronLargeTracksCount;;", 20, 0, 20);
    TH1D* hMuonLargeTracksCount           = new TH1D("hMuonLargeTracksCount", "hMuonLargeTracksCount;;", 20, 0, 20);

    double TRACK_THRESHOLD_MIN  = 1;
    double TRACK_THRESHOLD_MAX  = 30;
    double TRACK_THRESHOLD_STEP = 1;

    int NUMBER_TRACKS_MIN  = 1;
    int NUMBER_TRACKS_MAX  = 10;
    int NUMBER_TRACKS_STEP = 1;

    TH2D* hChargeExchangeMisIDs = new TH2D(
        "hChargeExchangeMisIDs", 
        ";Small Track Threshold; Number of small tracks",
        (TRACK_THRESHOLD_MAX - TRACK_THRESHOLD_MIN) / TRACK_THRESHOLD_STEP, TRACK_THRESHOLD_MIN, TRACK_THRESHOLD_MAX,
        (NUMBER_TRACKS_MAX - NUMBER_TRACKS_MIN) / NUMBER_TRACKS_STEP, NUMBER_TRACKS_MIN, NUMBER_TRACKS_MAX
    );

    TH2D* hChargeExchangeReco = new TH2D(
        "hChargeExchangeReco", 
        ";Small Track Threshold; Number of small tracks",
        (TRACK_THRESHOLD_MAX - TRACK_THRESHOLD_MIN) / TRACK_THRESHOLD_STEP, TRACK_THRESHOLD_MIN, TRACK_THRESHOLD_MAX,
        (NUMBER_TRACKS_MAX - NUMBER_TRACKS_MIN) / NUMBER_TRACKS_STEP, NUMBER_TRACKS_MIN, NUMBER_TRACKS_MAX
    );

    TH2D* hChargeExchangeRecoTrue = new TH2D(
        "hChargeExchangeRecoTrue",
        ";Small Track Threshold; Number of small tracks",
        (TRACK_THRESHOLD_MAX - TRACK_THRESHOLD_MIN) / TRACK_THRESHOLD_STEP, TRACK_THRESHOLD_MIN, TRACK_THRESHOLD_MAX,
        (NUMBER_TRACKS_MAX - NUMBER_TRACKS_MIN) / NUMBER_TRACKS_STEP, NUMBER_TRACKS_MIN, NUMBER_TRACKS_MAX
    );

    //////////////////////
    // Loop over events //
    //////////////////////

    Int_t NumEntries = (Int_t) tree->GetEntries();
    std::cout << "Num entries: " << NumEntries << std::endl;

    for (Int_t i = 0; i < NumEntries; ++i) {
        tree->GetEntry(i);

        // Categorize scatterings into 0p and Np scatterings
        int scatteringType = -1; // -1: not scattering, 0: 0p scattering, 1: Np scattering
        if (backgroundType == 12 || (backgroundType == 6 && numVisibleProtons == 0)) {
            scatteringType = 0;
        }
        else if (backgroundType == 6 && numVisibleProtons > 0) {
            scatteringType = 1;
        }

        // If no track matched to wire-chamber, skip
        if (WC2TPCtrkID == -99999) continue;

        // If no shower probability obtained, skip
        if (!obtainedProbabilities) continue;

        if (backgroundType == 2) {
            hMuonShowerProb->Fill(showerProb);
        } else if (backgroundType == 3) {
            hElectronShowerProb->Fill(showerProb);
        } else {
            hPionShowerProb->Fill(showerProb);
        }

        if (backgroundType == 7) {
            hPionChargeExchangeShowerProb->Fill(showerProb);
        } else if (backgroundType == 0) {
            hPionAbs0pShowerProb->Fill(showerProb);
        } else if (backgroundType == 1) {
            hPionAbsNpShowerProb->Fill(showerProb);
        } else if (scatteringType == 0) {
            hPion0pScatterShowerProb->Fill(showerProb);
        } else if (scatteringType == 1) {
            hPionNpScatterShowerProb->Fill(showerProb);
        } else if (backgroundType != 2 && backgroundType != 3) {
            hPionOtherShowerProb->Fill(showerProb);
        }

        // Apply shower probability cut
        if (showerProb >= SHOWER_PROB_CUT) continue;

        ///////////////////////
        // Primary track PID //
        ///////////////////////

        int totalCaloPoints = wcMatchDEDX->size();
        int nRemoveOutliers = 2;
        int nRemoveEnds     = 3;
        int minPoints       = 5;

        // Get chi^2 fits, primary tracks are already checked for reversal in first module
        double pionChi2   = computeReducedChi2(gPion, *wcMatchResR,  *wcMatchDEDX, false, totalCaloPoints, nRemoveOutliers, nRemoveEnds);
        double MIPChi2    = computeReducedChi2(gMuonTG, *wcMatchResR, *wcMatchDEDX, false, totalCaloPoints, nRemoveOutliers, nRemoveEnds);
        double protonChi2 = computeReducedChi2(gProton, *wcMatchResR, *wcMatchDEDX, false, totalCaloPoints, nRemoveOutliers, nRemoveEnds);

        double minStitchedChi2 = std::numeric_limits<double>::max();
        int bestBreakPoint = -1;
        if (totalCaloPoints >= 4 * nRemoveEnds + 2 * nRemoveOutliers + 2 * minPoints) {
            for (int caloBreakPoint = 2 * nRemoveEnds + nRemoveOutliers + minPoints; caloBreakPoint < totalCaloPoints - (2 * nRemoveEnds + nRemoveOutliers + minPoints); ++caloBreakPoint) {
                std::vector<double> leftResR(wcMatchResR->begin(), wcMatchResR->begin() + caloBreakPoint);
                std::vector<double> leftDEDX(wcMatchDEDX->begin(), wcMatchDEDX->begin() + caloBreakPoint);

                std::vector<double> rightResR(wcMatchResR->begin() + caloBreakPoint, wcMatchResR->end());
                std::vector<double> rightDEDX(wcMatchDEDX->begin() + caloBreakPoint, wcMatchDEDX->end());

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
        double breakPointX = WC2TPCPrimaryEndX; 
        double breakPointY = WC2TPCPrimaryEndY; 
        double breakPointZ = WC2TPCPrimaryEndZ;
        if (minChi2 == minStitchedChi2) {
            breakPointX = wcMatchXPos->at(bestBreakPoint); 
            breakPointY = wcMatchYPos->at(bestBreakPoint); 
            breakPointZ = wcMatchZPos->at(bestBreakPoint);    
        }

        /////////////////////////
        // Secondary track PID //
        /////////////////////////

        // Get reco tracks near vertex
        int secondaryTaggedPion   = 0;
        int secondaryTaggedProton = 0;
        int secondaryTaggedOther  = 0;

        int otherTaggedPion   = 0;
        int otherTaggedProton = 0;

        int numTracksNearVertex = 0;
        for (size_t trk_idx = 0; trk_idx < recoBeginX->size(); ++trk_idx) {
            // Have to re-check track ordering for stitched case
            double distanceFromStart = distance(
                recoBeginX->at(trk_idx), breakPointX, 
                recoBeginY->at(trk_idx), breakPointY,
                recoBeginZ->at(trk_idx), breakPointZ
            );
            double distanceFromEnd = distance(
                recoEndX->at(trk_idx), breakPointX, 
                recoEndY->at(trk_idx), breakPointY,
                recoEndZ->at(trk_idx), breakPointZ
            );

            if ((distanceFromStart < VERTEX_RADIUS || distanceFromEnd < VERTEX_RADIUS)) {
                numTracksNearVertex++;

                std::vector<double> secondaryDEDX = recoDEDX->at(trk_idx);
                std::vector<double> secondaryResR = recoResR->at(trk_idx);

                bool secondaryReversed  = false;
                bool originallyReversed = isTrackInverted->at(trk_idx);
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

        if (minChi2 == pionChi2) {
            continue; 
        } else if (minChi2 == protonChi2) {
            continue; // reject events matching to proton primary
        }

        int totalTaggedPions   = secondaryTaggedPion + otherTaggedPion;
        int totalTaggedProtons = secondaryTaggedProton + otherTaggedProton;

        if (backgroundType == 7) {
            hChargeExchangeTracksNVertex->Fill(numTracksNearVertex);
        } else if (backgroundType == 0) {
            hPionAbs0pTracksNVertex->Fill(numTracksNearVertex);
        } else if (backgroundType == 1) {
            hPionAbsNpTracksNVertex->Fill(numTracksNearVertex);
        } else if (scatteringType == 0) {
            hPion0pScatterTracksNVertex->Fill(numTracksNearVertex);
        } else if (scatteringType == 1) {
            hPionNpScatterTracksNVertex->Fill(numTracksNearVertex);
        } else if (backgroundType == 3) {
            hElectronTracksNVertex->Fill(numTracksNearVertex);
        } else if (backgroundType == 2) {
            hMuonTracksNVertex->Fill(numTracksNearVertex);
        } else if (backgroundType != 2 && backgroundType != 3) {
            hOtherPionTracksNVertex->Fill(numTracksNearVertex);
        }

        // Not charge exchange if we found a proton or pion near the vertex
        if (totalTaggedPions > 0 || totalTaggedProtons > 0) continue;

        ///////////////////////////
        // Track length analysis //
        ///////////////////////////

        int numSmallTracks = 0; int numLargeTracks = 0; int numFarTracks = 0;
        for (int iTrk = 0; iTrk < (int) recoTrkID->size(); ++iTrk) {
            double distanceFromStart = distance(
                recoBeginX->at(iTrk), breakPointX, 
                recoBeginY->at(iTrk), breakPointY,
                recoBeginZ->at(iTrk), breakPointZ
            );
            double distanceFromEnd = distance(
                recoEndX->at(iTrk), breakPointX, 
                recoEndY->at(iTrk), breakPointY,
                recoEndZ->at(iTrk), breakPointZ
            );

            if ((distanceFromStart < FAR_TRACK_DISTANCE || distanceFromEnd < FAR_TRACK_DISTANCE)) continue;
            numFarTracks++;

            double thisTrackLength = sqrt(
                pow(recoBeginX->at(iTrk) - recoEndX->at(iTrk), 2) +
                pow(recoBeginY->at(iTrk) - recoEndY->at(iTrk), 2) + 
                pow(recoBeginZ->at(iTrk) - recoEndZ->at(iTrk), 2)
            );
            if (thisTrackLength < SMALL_TRACK_LENGTH) numSmallTracks++;
            if (thisTrackLength > LARGE_TRACK_LENGTH) numLargeTracks++;

            if (backgroundType == 7) {
                hChargeExchangeTrackLengths->Fill(thisTrackLength);
            } else if (backgroundType != 2 && backgroundType != 3) {
                hOtherPionTrackLengths->Fill(thisTrackLength);
            } else if (backgroundType == 3) {
                hElectronTrackLengths->Fill(thisTrackLength);
            } else if (backgroundType == 2) {
                hMuonTrackLengths->Fill(thisTrackLength);
            }
        }

        if (backgroundType == 7) {
            hChargeExchangeSmallTracksCount->Fill(numSmallTracks);
            hChargeExchangeLargeTracksCount->Fill(numLargeTracks);
        } else if (backgroundType != 2 && backgroundType != 3) {
            hOtherPionSmallTracksCount->Fill(numSmallTracks);
            hOtherPionLargeTracksCount->Fill(numLargeTracks);
        } else if (backgroundType == 3) {
            hElectronSmallTracksCount->Fill(numSmallTracks);
            hElectronLargeTracksCount->Fill(numLargeTracks);
        } else if (backgroundType == 2) {
            hMuonSmallTracksCount->Fill(numSmallTracks);
            hMuonLargeTracksCount->Fill(numLargeTracks);
        }

        if (backgroundType == 7) {
            hChargeExchangeTracksFarVertex->Fill(numFarTracks);
        } else if (backgroundType == 0) {
            hPionAbs0pTracksFarVertex->Fill(numFarTracks);
        } else if (backgroundType == 1) {
            hPionAbsNpTracksFarVertex->Fill(numFarTracks);
        } else if (scatteringType == 0) {
            hPion0pScatterTracksFarVertex->Fill(numFarTracks);
        } else if (scatteringType == 1) {
            hPionNpScatterTracksFarVertex->Fill(numFarTracks);
        } else if (backgroundType == 3) {
            hElectronTracksFarVertex->Fill(numFarTracks);
        } else if (backgroundType == 2) {
            hMuonTracksFarVertex->Fill(numFarTracks);
        } else if (backgroundType != 2 && backgroundType != 3) {
            hOtherPionTracksFarVertex->Fill(numFarTracks);
        } 

        // Fill misID histogram
        for (double threshold = TRACK_THRESHOLD_MIN; threshold <= TRACK_THRESHOLD_MAX; threshold += TRACK_THRESHOLD_STEP) {
            for (int iTrk = 0; iTrk < (int) recoTrkID->size(); ++iTrk) {
                double thisTrackLength = sqrt(
                    pow(recoBeginX->at(iTrk) - recoEndX->at(iTrk), 2) +
                    pow(recoBeginY->at(iTrk) - recoEndY->at(iTrk), 2) + 
                    pow(recoBeginZ->at(iTrk) - recoEndZ->at(iTrk), 2)
                );
                if (thisTrackLength < threshold) numSmallTracks++;
            }

            for (int numTrksBound = NUMBER_TRACKS_MIN; numTrksBound <= NUMBER_TRACKS_MAX; numTrksBound += NUMBER_TRACKS_STEP) {
                // Look for misidentifications
                if (((numSmallTracks > numTrksBound) && (backgroundType != 7)) || ((numSmallTracks <= numTrksBound) && (backgroundType == 7))) {
                    hChargeExchangeMisIDs->Fill(threshold, numTrksBound, 1. / NumEntries);
                }

                if (numSmallTracks > numTrksBound) {
                    hChargeExchangeReco->Fill(threshold, numTrksBound, 1);
                    if (backgroundType == 7) {
                        hChargeExchangeRecoTrue->Fill(threshold, numTrksBound, 1);
                    }
                }
            }
        }
    }
    
    TH2D* hChargeExchangePurity = (TH2D*) hChargeExchangeRecoTrue->Clone("hChargeExchangePurity");
    hChargeExchangePurity->SetTitle(";Small Track Threshold; Number of small tracks;Purity");
    hChargeExchangePurity->Divide(hChargeExchangeRecoTrue, hChargeExchangeReco, 1.0, 1.0);

    /////////////////////////////////
    // Print one dimensional plots //
    /////////////////////////////////

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
        // Shower probability
        {hElectronShowerProb, hPionShowerProb, hMuonShowerProb},
        {hPionChargeExchangeShowerProb, hPionAbs0pShowerProb, hPionAbsNpShowerProb, hPion0pScatterShowerProb, hPionNpScatterShowerProb, hPionOtherShowerProb},

        // Tracks near vertex
        {hChargeExchangeTracksNVertex, hPionAbs0pTracksNVertex, hPionAbsNpTracksNVertex, hPion0pScatterTracksNVertex, hPionNpScatterTracksNVertex, hOtherPionTracksNVertex, hElectronTracksNVertex, hMuonTracksNVertex},

        // Tracks away from vertex
        {hChargeExchangeTracksFarVertex, hPionAbs0pTracksFarVertex, hPionAbsNpTracksFarVertex, hPion0pScatterTracksFarVertex, hPionNpScatterTracksFarVertex, hOtherPionTracksFarVertex, hElectronTracksFarVertex, hMuonTracksFarVertex},

        // Track lengths
        {hChargeExchangeTrackLengths, hOtherPionTrackLengths, hElectronTrackLengths, hMuonTrackLengths},
        {hChargeExchangeSmallTracksCount, hOtherPionSmallTracksCount, hElectronSmallTracksCount, hMuonSmallTracksCount},
        {hChargeExchangeLargeTracksCount, hOtherPionLargeTracksCount, hElectronLargeTracksCount, hMuonLargeTracksCount}
    };

    std::vector<std::vector<TString>> PlotLabelGroups = {
        // Shower probability
        {"Electron", "Pion", "Muon"},
        {"Ch. exch.", "Abs. 0p", "Abs. Np", "0p scatter", "Np scatter", "Other"},

        // Tracks near vertex
        {"Ch. exch.", "Abs. 0p", "Abs. Np", "0p scatter", "Np scatter", "Other pion", "Electron", "Muon"},

        // Tracks away from vertex
        {"Ch. exch.", "Abs. 0p", "Abs. Np", "0p scatter", "Np scatter", "Other pion", "Electron", "Muon"},

        // Track lengths
        {"Ch. exch.", "Other pion", "Electron", "Muon"},
        {"Ch. exch.", "Other pion", "Electron", "Muon"},
        {"Ch. exch.", "Other pion", "Electron", "Muon"}
    };

    std::vector<TString> PlotTitles = {
        // Shower probability
        "PrimaryTrackShowerProb",
        "PionInteractionShowerProb",

        // Tracks near vertex
        "TracksNearVertex",

        // Tracks away from vertex
        "TracksAwayFromVertex",

        // Track lengths
        "TrackLengths",
        "SmallTrackCounts",
        "LargeTrackCounts"
    };

    std::vector<TString> XLabels = {
        // Shower probability
        "Shower probability",
        "Shower probability",

        // Tracks near vertex
        "Number of tracks near vertex",

        // Tracks away from vertex
        "Number of tracks away from vertex",

        // Track lengths
        "Track length (cm)",
        "Number of small tracks",
        "Number of large tracks"
    };

    std::vector<TString> YLabels = {
        // Shower probability
        "Number of events",
        "Number of events",

        // Tracks near vertex
        "Number of events",

        // Tracks away from vertex
        "Number of events",

        // Track lengths
        "Number of tracks",
        "Number of events",
        "Number of events"
    };

    std::vector<bool> PlotStacked = {
        // Shower probability
        true,
        true,

        // Tracks near vertex
        true,

        // Tracks away from vertex
        true,

        // Track lengths
        true,
        true,
        true
    };

    printOneDPlots(
        SaveDir, FontStyle, TextSize,
        PlotGroups,
        Colors,
        PlotLabelGroups,
        PlotTitles,
        XLabels,
        YLabels,
        PlotStacked
    );

    /////////////////////////////////
    // Print two dimensional plots //
    /////////////////////////////////

    std::vector<TH2*> TwoDPlots = {
        hChargeExchangeMisIDs,
        hChargeExchangeReco,
        hChargeExchangeRecoTrue,
        hChargeExchangePurity
    };

    std::vector<TString> TwoDTitles = {
        "ChargeExchangeMisIDs",
        "ChargeExchangeReco",
        "ChargeExchangeRecoTrue",
        "ChargeExchangePurity"
    };

    printTwoDPlots(SaveDir, TwoDPlots, TwoDTitles);
}

void printOneDPlots(
    const TString& dir, 
    int fontStyle, 
    double textSize,
    std::vector<std::vector<TH1*>>& groups,
    std::vector<int>& colors, 
    std::vector<std::vector<TString>>& labels, 
    std::vector<TString>& titles, 
    std::vector<TString>& xlabels, 
    std::vector<TString>& ylabels,
    std::vector<bool>& stack
) {
    int numPlots = groups.size();
    for (int iPlot = 0; iPlot < numPlots; ++iPlot) {
        // Set up canvas
        TCanvas* PlotCanvas = new TCanvas("Canvas", "Canvas", 205, 34, 1300, 768);

        TPad* mainPad = new TPad("mainPad","",0.0, 0.0, 0.80, 1.0);
        mainPad->SetRightMargin(0.05);
        mainPad->SetLeftMargin(0.15);
        mainPad->SetBottomMargin(0.15);
        mainPad->Draw();
        mainPad->cd();

        TLegend* leg = new TLegend(0.0, 0.0, 1.0, 1.0);
        leg->SetTextSize(textSize * 2.5);
        leg->SetTextFont(fontStyle);

        TPad* legendPad = new TPad("legendPad", "", 0.80, 0.4, 1.0, 0.8);
        legendPad->SetLeftMargin(0.05);
        legendPad->SetRightMargin(0.05);
        legendPad->SetBottomMargin(0.15);
        legendPad->SetFillStyle(0);
        legendPad->SetBorderMode(0);

        // Get histograms and labels
        std::vector<TH1*> Plots = groups.at(iPlot);
        std::vector<TString> Labels = labels.at(iPlot);	

        // If stacked
        if (stack.at(iPlot)) {
            THStack stack("stack", titles.at(iPlot));

            // Style and add each histogram to the stack
            for (int iSubPlot = 0; iSubPlot < (int) Plots.size(); ++iSubPlot) {
                TH1* h = Plots[iSubPlot];
                leg->AddEntry(h, Labels[iSubPlot], "f");
                h->SetLineWidth(2);
                h->SetLineColor(colors.at(iSubPlot));
                h->SetFillColor(colors.at(iSubPlot));
                h->SetFillColorAlpha(colors.at(iSubPlot), 0.2);
                h->SetFillStyle(3001);
                stack.Add(h, "H");
            }

            // Style the stack
            stack.SetTitle(titles.at(iPlot));

            // Draw stack
            stack.Draw("hist");

            stack.GetXaxis()->SetTitleFont(fontStyle);
            stack.GetXaxis()->SetLabelFont(fontStyle);
            stack.GetXaxis()->SetNdivisions(8);
            stack.GetXaxis()->SetLabelSize(textSize);
            stack.GetXaxis()->SetTitleSize(textSize);
            stack.GetXaxis()->SetTitle(xlabels.at(iPlot));
            stack.GetXaxis()->SetTitleOffset(1.1);
            stack.GetXaxis()->CenterTitle();

            stack.GetYaxis()->SetTitleFont(fontStyle);
            stack.GetYaxis()->SetLabelFont(fontStyle);
            stack.GetYaxis()->SetNdivisions(6);
            stack.GetYaxis()->SetLabelSize(textSize);
            stack.GetYaxis()->SetTitleSize(textSize);
            stack.GetYaxis()->SetTitle(ylabels.at(iPlot));
            stack.GetYaxis()->SetTitleOffset(1.1);
            stack.GetYaxis()->CenterTitle();

            // Determine proper max from stack object
            double stackMax = stack.GetMaximum();
            double YAxisRange = 1.15 * stackMax;
            stack.SetMaximum(YAxisRange);
            stack.SetMinimum(0.0);

            gPad->Update();
            TAxis *x = Plots[0]->GetXaxis();
            x->SetMaxDigits(3);
            gPad->Modified();
            gPad->Update();

            PlotCanvas->cd();
            legendPad->Draw();
            legendPad->cd();
            leg->Draw();

            PlotCanvas->SaveAs(dir + titles.at(iPlot) + ".png");
        } 
        // If not stacked
        else {
            // Style the first plot
            Plots[0]->SetTitle(titles.at(iPlot));

            Plots[0]->GetXaxis()->SetTitleFont(fontStyle);
            Plots[0]->GetXaxis()->SetLabelFont(fontStyle);
            Plots[0]->GetXaxis()->SetNdivisions(8);
            Plots[0]->GetXaxis()->SetLabelSize(textSize);
            Plots[0]->GetXaxis()->SetTitleSize(textSize);
            Plots[0]->GetXaxis()->SetTitle(xlabels.at(iPlot));
            Plots[0]->GetXaxis()->SetTitleOffset(1.1);
            Plots[0]->GetXaxis()->CenterTitle();

            Plots[0]->GetYaxis()->SetTitleFont(fontStyle);
            Plots[0]->GetYaxis()->SetLabelFont(fontStyle);
            Plots[0]->GetYaxis()->SetNdivisions(6);
            Plots[0]->GetYaxis()->SetLabelSize(textSize);
            Plots[0]->GetYaxis()->SetTitleSize(textSize);
            Plots[0]->GetYaxis()->SetTitle(ylabels.at(iPlot));
            Plots[0]->GetYaxis()->SetTitleOffset(1.1);
            Plots[0]->GetYaxis()->CenterTitle();

            double imax;
            for (int iSubPlot = 0; iSubPlot < (int) Plots.size(); ++iSubPlot) {
                leg->AddEntry(Plots[iSubPlot], Labels[iSubPlot], "f");
                Plots[iSubPlot]->SetLineWidth(2);
                Plots[iSubPlot]->SetLineColor(colors.at(iSubPlot));
                Plots[iSubPlot]->Draw("hist same");

                imax = TMath::Max(Plots[iSubPlot]->GetMaximum(), Plots[0]->GetMaximum());

                double YAxisRange = 1.15 * imax;
                Plots[iSubPlot]->GetYaxis()->SetRangeUser(0., YAxisRange);
                Plots[0]->GetYaxis()->SetRangeUser(0., YAxisRange);	
            }

            gPad->Update();
            TAxis *x = Plots[0]->GetXaxis();
            x->SetMaxDigits(3);
            gPad->Modified();
            gPad->Update();

            PlotCanvas->cd();
            legendPad->Draw();
            legendPad->cd();
            leg->Draw();

            PlotCanvas->SaveAs(dir + titles.at(iPlot) + ".png");
        }
        delete PlotCanvas;
    }
}

void printTwoDPlots(
    const TString& dir, 
    const std::vector<TH2*>& plots, 
    const std::vector<TString>& titles
) {
    int nPlots = plots.size();

    TCanvas* c1 = new TCanvas("c1", "EnergyLossPlots", 800, 600);
    for (int iPlot = 0; iPlot < nPlots; ++iPlot) {
        TH1* hPlot = plots.at(iPlot);
        hPlot->SetMinimum(0);
        hPlot->SetMaximum(hPlot->GetMaximum());
        hPlot->Draw("COLZ");
        c1->SaveAs(dir + titles.at(iPlot) + ".png");
    }
}

double distance(double x1, double x2, double y1, double y2, double z1, double z2) {
    return sqrt(
        pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)
    );
}

double computeReducedChi2(const TGraph* theory, std::vector<double> xData, std::vector<double> yData, bool dataReversed, int nPoints, int nOutliersToDiscard = 0, int nTrim = 0) {
    std::vector<double> chi2Contributions;
    if (dataReversed) {
        std::reverse(xData.begin(), xData.end());
        std::reverse(yData.begin(), yData.end());
    }

    int startIdx = 0;
    int endIdx   = nPoints;
    if (2 * nTrim < nPoints) {
        startIdx = nTrim;
        endIdx   = nPoints - nTrim;
    }

    for (int i = startIdx; i < endIdx; ++i) {
        double theoryY = theory->Eval(xData[i]); // interpolate the theory at xData[i]
        double deltaY  = yData[i] - theoryY;
        double chi2    = (deltaY * deltaY) / std::abs(theoryY);
        chi2Contributions.push_back(chi2);
    }

    std::sort(chi2Contributions.begin(), chi2Contributions.end());
    int totalContrib = chi2Contributions.size();
    int nUsed        = std::max(0, totalContrib - nOutliersToDiscard);
    double chi2Sum   = 0.0;
    for (int i = 0; i < nUsed; ++i) {
        chi2Sum += chi2Contributions[i];
    }

    // Return reduced chi²
    return nUsed > 0 ? chi2Sum / nUsed : 0.0;
}

double meanDEDX(std::vector<double>& trackDEDX, bool isTrackReversed, int pointsToUse) {
    double dEdx = 0.;
    if (isTrackReversed) std::reverse(trackDEDX.begin(), trackDEDX.end());

    unsigned int bound = pointsToUse;
    if (pointsToUse > trackDEDX.size()) bound = trackDEDX.size();
    for (unsigned int i = 0; i < bound; ++i) dEdx += trackDEDX.at(i);
    dEdx /= bound;
    return dEdx;
}

void initializeProtonPoints(TGraph* gProton) {
    double protonData[107][2] = {
        {31.95, 4.14}, {31.65, 4.16}, {31.35, 4.17}, {31.05, 4.18}, {30.75, 4.20},
        {30.45, 4.21}, {30.15, 4.23}, {29.85, 4.25}, {29.55, 4.26}, {29.25, 4.28},
        {28.95, 4.29}, {28.65, 4.31}, {28.35, 4.33}, {28.05, 4.34}, {27.75, 4.36},
        {27.45, 4.38}, {27.15, 4.40}, {26.85, 4.42}, {26.55, 4.43}, {26.25, 4.45},
        {25.95, 4.47}, {25.65, 4.49}, {25.35, 4.51}, {25.05, 4.53}, {24.75, 4.55},
        {24.45, 4.57}, {24.15, 4.60}, {23.85, 4.62}, {23.55, 4.64}, {23.25, 4.66},
        {22.95, 4.69}, {22.65, 4.71}, {22.35, 4.73}, {22.05, 4.76}, {21.75, 4.78},
        {21.45, 4.81}, {21.15, 4.83}, {20.85, 4.86}, {20.55, 4.89}, {20.25, 4.92},
        {19.95, 4.94}, {19.65, 4.97}, {19.35, 5.00}, {19.05, 5.03}, {18.75, 5.07},
        {18.45, 5.10}, {18.15, 5.13}, {17.85, 5.16}, {17.55, 5.20}, {17.25, 5.23},
        {16.95, 5.27}, {16.65, 5.31}, {16.35, 5.35}, {16.05, 5.39}, {15.75, 5.43},
        {15.45, 5.47}, {15.15, 5.51}, {14.85, 5.56}, {14.55, 5.60}, {14.25, 5.65},
        {13.95, 5.70}, {13.65, 5.75}, {13.35, 5.80}, {13.05, 5.85}, {12.75, 5.91},
        {12.45, 5.97}, {12.15, 6.03}, {11.85, 6.09}, {11.55, 6.15}, {11.25, 6.22},
        {10.95, 6.29}, {10.65, 6.36}, {10.35, 6.44}, {10.05, 6.52}, {9.75, 6.60},
        {9.45, 6.68}, {9.15, 6.77}, {8.85, 6.87}, {8.55, 6.97}, {8.25, 7.08},
        {7.95, 7.19}, {7.65, 7.30}, {7.35, 7.43}, {7.05, 7.56}, {6.75, 7.70},
        {6.45, 7.85}, {6.15, 8.02}, {5.85, 8.19}, {5.55, 8.38}, {5.25, 8.58},
        {4.95, 8.81}, {4.65, 9.05}, {4.35, 9.32}, {4.05, 9.61}, {3.75, 9.94},
        {3.45, 10.32}, {3.15, 10.74}, {2.85, 11.23}, {2.55, 11.80}, {2.25, 12.48},
        {1.95, 13.31}, {1.65, 14.35}, {1.35, 15.71}, {1.05, 17.59}, {0.75, 20.44},
        {0.45, 25.48}, {0.15, 38.12}
    };

    for (int i = 0; i < 107; ++i) {
        gProton->SetPoint(i, protonData[i][0], protonData[i][1]);
    }
}

void initializePionPoints(TGraph* gPion) {
    double pionData[107][2] = {
        {31.95, 2.4}, {31.65, 2.4}, {31.35, 2.4}, {31.05, 2.4}, {30.75, 2.4},
        {30.45, 2.4}, {30.15, 2.4}, {29.85, 2.4}, {29.55, 2.4}, {29.25, 2.4},
        {28.95, 2.4}, {28.65, 2.4}, {28.35, 2.4}, {28.05, 2.4}, {27.75, 2.5},
        {27.45, 2.5}, {27.15, 2.5}, {26.85, 2.5}, {26.55, 2.5}, {26.25, 2.5},
        {25.95, 2.5}, {25.65, 2.5}, {25.35, 2.5}, {25.05, 2.5}, {24.75, 2.5},
        {24.45, 2.5}, {24.15, 2.5}, {23.85, 2.5}, {23.55, 2.5}, {23.25, 2.6},
        {22.95, 2.6}, {22.65, 2.6}, {22.35, 2.6}, {22.05, 2.6}, {21.75, 2.6},
        {21.45, 2.6}, {21.15, 2.6}, {20.85, 2.6}, {20.55, 2.6}, {20.25, 2.6},
        {19.95, 2.6}, {19.65, 2.7}, {19.35, 2.7}, {19.05, 2.7}, {18.75, 2.7},
        {18.45, 2.7}, {18.15, 2.7}, {17.85, 2.7}, {17.55, 2.7}, {17.25, 2.8},
        {16.95, 2.8}, {16.65, 2.8}, {16.35, 2.8}, {16.05, 2.8}, {15.75, 2.8},
        {15.45, 2.8}, {15.15, 2.9}, {14.85, 2.9}, {14.55, 2.9}, {14.25, 2.9},
        {13.95, 2.9}, {13.65, 2.9}, {13.35, 3.0}, {13.05, 3.0}, {12.75, 3.0},
        {12.45, 3.0}, {12.15, 3.0}, {11.85, 3.1}, {11.55, 3.1}, {11.25, 3.1},
        {10.95, 3.1}, {10.65, 3.2}, {10.35, 3.2}, {10.05, 3.2}, {9.75, 3.3},
        {9.45, 3.3}, {9.15, 3.3}, {8.85, 3.4}, {8.55, 3.4}, {8.25, 3.4},
        {7.95, 3.5}, {7.65, 3.5}, {7.35, 3.6}, {7.05, 3.6}, {6.75, 3.7},
        {6.45, 3.7}, {6.15, 3.8}, {5.85, 3.9}, {5.55, 3.9}, {5.25, 4.0},
        {4.95, 4.1}, {4.65, 4.2}, {4.35, 4.3}, {4.05, 4.4}, {3.75, 4.6},
        {3.45, 4.7}, {3.15, 4.9}, {2.85, 5.1}, {2.55, 5.3}, {2.25, 5.6},
        {1.95, 5.9}, {1.65, 6.4}, {1.35, 6.9}, {1.05, 7.7}, {0.75, 8.9},
        {0.45, 11.0}, {0.15, 16.5}
    };

    for (int i = 0; i < 107; ++i) {
        gPion->SetPoint(i, pionData[i][0], pionData[i][1]);
    }
}

void initializeMuonNoBraggPoints(TGraph* gMuonTG) {
    double muonTGData[107][2] = {
        {31.95, 2.3}, {31.65, 2.3}, {31.35, 2.3}, {31.05, 2.3}, {30.75, 2.3},
        {30.45, 2.3}, {30.15, 2.3}, {29.85, 2.3}, {29.55, 2.3}, {29.25, 2.3},
        {28.95, 2.3}, {28.65, 2.3}, {28.35, 2.3}, {28.05, 2.3}, {27.75, 2.3},
        {27.45, 2.3}, {27.15, 2.3}, {26.85, 2.3}, {26.55, 2.3}, {26.25, 2.3},
        {25.95, 2.3}, {25.65, 2.3}, {25.35, 2.3}, {25.05, 2.3}, {24.75, 2.3},
        {24.45, 2.3}, {24.15, 2.3}, {23.85, 2.3}, {23.55, 2.3}, {23.25, 2.3},
        {22.95, 2.3}, {22.65, 2.3}, {22.35, 2.3}, {22.05, 2.3}, {21.75, 2.3},
        {21.45, 2.3}, {21.15, 2.3}, {20.85, 2.3}, {20.55, 2.3}, {20.25, 2.3},
        {19.95, 2.3}, {19.65, 2.3}, {19.35, 2.3}, {19.05, 2.3}, {18.75, 2.3},
        {18.45, 2.3}, {18.15, 2.3}, {17.85, 2.3}, {17.55, 2.3}, {17.25, 2.3},
        {16.95, 2.3}, {16.65, 2.3}, {16.35, 2.3}, {16.05, 2.3}, {15.75, 2.3},
        {15.45, 2.3}, {15.15, 2.3}, {14.85, 2.3}, {14.55, 2.3}, {14.25, 2.3},
        {13.95, 2.3}, {13.65, 2.3}, {13.35, 2.3}, {13.05, 2.3}, {12.75, 2.3},
        {12.45, 2.3}, {12.15, 2.3}, {11.85, 2.3}, {11.55, 2.3}, {11.25, 2.3},
        {10.95, 2.3}, {10.65, 2.3}, {10.35, 2.3}, {10.05, 2.3}, {9.75, 2.3},
        {9.45, 2.3}, {9.15, 2.3}, {8.85, 2.3}, {8.55, 2.3}, {8.25, 2.3},
        {7.95, 2.3}, {7.65, 2.3}, {7.35, 2.3}, {7.05, 2.3}, {6.75, 2.3},
        {6.45, 2.3}, {6.15, 2.3}, {5.85, 2.3}, {5.55, 2.3}, {5.25, 2.3},
        {4.95, 2.3}, {4.65, 2.3}, {4.35, 2.3}, {4.05, 2.3}, {3.75, 2.3},
        {3.45, 2.3}, {3.15, 2.3}, {2.85, 2.3}, {2.55, 2.3}, {2.25, 2.3},
        {1.95, 2.3}, {1.65, 2.3}, {1.35, 2.3}, {1.05, 2.3}, {0.75, 2.3},
        {0.45, 2.3}, {0.15, 2.3}
    };

    for (int i = 0; i < 107; ++i) {
        gMuonTG->SetPoint(i, muonTGData[i][0], muonTGData[i][1]);
    }
}