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

void RecoClassifyAllSimplified() {
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
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis/figs/ClassifyAllSimplified/";

    // Load file with NN data products
    TString RootFilePath = "/exp/lariat/app/users/epelaez/files/RecoNNAllEval_histo.root";
    std::unique_ptr<TFile> File(TFile::Open(RootFilePath));
    TDirectory* Directory = (TDirectory*)File->Get("RecoNNAllEval");

    ///////////////////
    // Load branches //
    ///////////////////

    // Load tree and branches
    TTree* tree = (TTree*) Directory->Get<TTree>("RecoNNAllEvalTree");
    TH1D* hTotalEvents = (TH1D*) Directory->Get<TH1D>("hTotalEvents");

    int run, subrun, event; bool isData;
    tree->SetBranchAddress("run", &run); 
    tree->SetBranchAddress("subrun", &subrun); 
    tree->SetBranchAddress("event", &event);
    tree->SetBranchAddress("isData", &isData);

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

    double showerNoBoxProb;
    bool obtainedNoBoxProbabilities;
    tree->SetBranchAddress("showerNoBoxProb", &showerNoBoxProb);
    tree->SetBranchAddress("obtainedNoBoxProbabilities", &obtainedNoBoxProbabilities);

    double showerOutsideBoxProb;
    bool obtainedOutsideBoxProbabilities;
    tree->SetBranchAddress("showerOutsideBoxProb", &showerOutsideBoxProb);
    tree->SetBranchAddress("obtainedOutsideBoxProbabilities", &obtainedOutsideBoxProbabilities);

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

    // Truth level incident KE information
    bool validTrueIncidentKE;
    std::vector<double>* trueIncidentKEContributions = nullptr;
    tree->SetBranchAddress("validTrueIncidentKE", &validTrueIncidentKE);
    tree->SetBranchAddress("trueIncidentKEContributions", &trueIncidentKEContributions);

    /////////////////////////////////
    // Files for event information //
    /////////////////////////////////

    std::ofstream outFileChExchBkg("files/ClassifyAll/ChExchBackground.txt");
    TFile* comparisonsFile = new TFile("/exp/lariat/app/users/epelaez/files/DataMCComparisons.root", "UPDATE");

    ///////////////////////
    // Create histograms //
    ///////////////////////

    // Histograms with classified events
    TH1D* hPionAbs0p     = new TH1D("hPassShowerProb", "hPassShowerProb;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D* hPionAbsNp     = new TH1D("hPionAbsNp", "hPionAbsNp;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D* hPionScatter   = new TH1D("hPionScatter", "hPionScatter;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D* hPionChExch    = new TH1D("hPionChExch", "hPionChExch;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);

    // Events composition after each cut
    TH1D* hDataProdsAndWC2TPC = new TH1D("hDataProdsAndWC2TPC", "hDataProdsAndWC2TPC;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D* hNotAnElectron      = new TH1D("hNotAnElectron", "hNotAnElectron;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D* hPrimaryInRedVol    = new TH1D("hPrimaryInRedVol", "hPrimaryInRedVol;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D* hPrimaryPID         = new TH1D("hPrimaryPID", "hPrimaryPID;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D* hNotScatter         = new TH1D("hNotScatter", "hNotScatter;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D* hNotPionAbsNp       = new TH1D("hNotPionAbsNp", "hNotPionAbsNp;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D* hNotChExch          = new TH1D("hNotChExch", "hNotChExch;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D* hNotPionAbs0p       = new TH1D("hNotPionAbs0p", "hNotPionAbs0p;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);

    // Incident kinetic energy
    TH1D* hIncidentKE         = new TH1D("hIncidentKE", "hIncidentKE;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hIncidentKEPion     = new TH1D("hIncidentKEPion", "hIncidentKEPion;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hIncidentKEElectron = new TH1D("hIncidentKEElectron", "hIncidentKEElectron;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hIncidentKEMuon     = new TH1D("hIncidentKEMuon", "hIncidentKEMuon;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueIncidentKE     = new TH1D("hTrueIncidentKE", "hTrueIncidentKE;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    // Interacting pion abs 0p energy
    TH1D* hPionAbs0pKE         = new TH1D("hPionAbs0pKE", "hPionAbs0pKE;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbs0pKETrue     = new TH1D("hPionAbs0PKETrue", "hPionAbs0PKETrue;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbs0pKEAbsNp    = new TH1D("hPionAbs0pKEAbsNp", "hPionAbs0pKEAbsNp;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbs0pKEScatter  = new TH1D("hPionAbs0pKEScatter", "hPionAbs0pKEScatter;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbs0pKEChExch   = new TH1D("hPionAbs0pKEChExch", "hPionAbs0pKEChExch;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbs0pKEMuon     = new TH1D("hPionAbs0pKEMuon", "hPionAbs0pKEMuon;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbs0pKEElectron = new TH1D("hPionAbs0pKEElectron", "hPionAbs0pKEElectron;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbs0pKEOther    = new TH1D("hPionAbs0pKEOther", "hPionAbs0pKEOther;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    // Estimated backgrounds for pion abs 0p
    std::vector<TH1*> PionAbs0pBkg = {
        hPionAbs0pKEMuon, hPionAbs0pKEElectron, hPionAbs0pKEOther
    };

    // Interacting pion abs Np energy
    TH1D* hPionAbsNpKE         = new TH1D("hPionAbsNpKE", "hPionAbsNpKE;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbsNpKETrue     = new TH1D("hPionAbsNpKETrue", "hPionAbsNpKETrue;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbsNpKEAbs0p    = new TH1D("hPionAbsNpKEAbs0p", "hPionAbsNpKEAbs0p;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbsNpKEScatter  = new TH1D("hPionAbsNpKEScatter", "hPionAbsNpKEScatter;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbsNpKEChExch   = new TH1D("hPionAbsNpKEChExch", "hPionAbsNpKEChExch;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbsNpKEMuon     = new TH1D("hPionAbsNpKEMuon", "hPionAbsNpKEMuon;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbsNpKEElectron = new TH1D("hPionAbsNpKEElectron", "hPionAbsNpKEElectron;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbsNpKEOther    = new TH1D("hPionAbsNpKEOther", "hPionAbsNpKEOther;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    // Estimated backgrounds for pion abs Np
    std::vector<TH1*> PionAbsNpBkg = {
        hPionAbsNpKEMuon, hPionAbsNpKEElectron, hPionAbsNpKEOther
    };

    // Interacting pion scattering energy
    TH1D* hPionScatterKE         = new TH1D("hPionScatterKE", "hPionScatterKE;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionScatterKETrue     = new TH1D("hPionScatterKETrue", "hPionScatterKETrue;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionScatterKEAbs0p    = new TH1D("hPionScatterKEAbs0p", "hPionScatterKEAbs0p;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionScatterKEAbsNp    = new TH1D("hPionScatterKEAbsNp", "hPionScatterKEAbsNp;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionScatterKEChExch   = new TH1D("hPionScatterKEChExch", "hPionScatterKEChExch;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionScatterKEMuon     = new TH1D("hPionScatterKEMuon", "hPionScatterKEMuon;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionScatterKEElectron = new TH1D("hPionScatterKEElectron", "hPionScatterKEElectron;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionScatterKEOther    = new TH1D("hPionScatterKEOther", "hPionScatterKEOther;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    // Estimated backgrounds for pion scattering
    std::vector<TH1*> PionScatterBkg = {
        hPionScatterKEMuon, hPionScatterKEElectron, hPionScatterKEOther
    };

    // Interacting pion charge exchange energy
    TH1D* hPionChExchKE         = new TH1D("hPionChExchKE", "hPionChExchKE;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionChExchKETrue     = new TH1D("hPionChExchKETrue", "hPionChExchKETrue;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionChExchKEAbs0p    = new TH1D("hPionChExchKEAbs0p", "hPionChExchKEAbs0p;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionChExchKEAbsNp    = new TH1D("hPionChExchKEAbsNp", "hPionChExchKEAbsNp;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionChExchKEScatter  = new TH1D("hPionChExchKEScatter", "hPionChExchKEScatter;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionChExchKEMuon     = new TH1D("hPionChExchKEMuon", "hPionChExchKEMuon;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionChExchKEElectron = new TH1D("hPionChExchKEElectron", "hPionChExchKEElectron;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionChExchKEOther    = new TH1D("hPionChExchKEOther", "hPionChExchKEOther;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    // Estimated backgrounds for pion charge exchange
    std::vector<TH1*> PionChExchBkg = {
        hPionChExchKEMuon, hPionChExchKEElectron, hPionChExchKEOther
    };

    // All our reconstructed backgrounds 
    std::vector<std::vector<TH1*>> RecoSignalBackgrounds = {
        PionAbs0pBkg, PionAbsNpBkg, PionScatterBkg, PionChExchBkg
    };

    // All our reconstructed signals
    std::vector<TH1*> RecoSignals = {
        hPionAbs0pKE, hPionAbsNpKE, hPionScatterKE, hPionChExchKE
    };

    // True abs 0p
    TH1D* hTrueAbs0pKE            = new TH1D("hTrueAbs0pKE", "hTrueAbs0pKE;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbs0pKEAsAbs0p     = new TH1D("hTrueAbs0pKEAsAbs0p", "hTrueAbs0pKEAsAbs0p;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbs0pKEAsAbsNp     = new TH1D("hTrueAbs0pKEAsAbsNp", "hTrueAbs0pKEAsAbsNp;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbs0pKEAsScatter   = new TH1D("hTrueAbs0pKEAsScatter", "hTrueAbs0pKEAsScatter;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbs0pKEAsChExch    = new TH1D("hTrueAbs0pKEAsChExch", "hTrueAbs0pKEAsChExch;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbs0pKERejected    = new TH1D("hTrueAbs0pKERejected", "hTrueAbs0pKERejected;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    // For each TRUE energy bin, what are abs 0p events reconstructed as?
    std::vector<std::vector<TH1*>> TrueAbs0pAsByBin;
    for (int iEnergyBin = 0; iEnergyBin < NUM_BINS_KE; ++iEnergyBin) {
        std::vector<TH1*> TempVec;
        for (int iInt = 0; iInt < NUM_SIGNAL_TYPES; ++iInt) {
            TH1* hTempHist = new TH1D(Form("hTrueAbs0p_%d_Bin_As_%d", iEnergyBin, iInt), Form("True Abs 0p KE As %d in bin %d", iInt, iEnergyBin), NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
            TempVec.push_back(hTempHist);
        }
        TrueAbs0pAsByBin.push_back(TempVec);
    }

    TH1D* hTrueAbs0pKERejDataProds = new TH1D("hTrueAbs0pKERejDataProds", "hTrueAbs0pKERejDataProds;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbs0pKERejElectron  = new TH1D("hTrueAbs0pKERejElectron", "hTrueAbs0pKERejElectron;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbs0pKERejRedVol    = new TH1D("hTrueAbs0pKERejRedVol", "hTrueAbs0pKERejRedVol;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbs0pKERejPID       = new TH1D("hTrueAbs0pKERejPID", "hTrueAbs0pKERejPID;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbs0pKERejManyPions = new TH1D("hTrueAbs0pKERejManyPions", "hTrueAbs0pKERejManyPions;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbs0pKERejClusters  = new TH1D("hTrueAbs0pKERejClusters", "hTrueAbs0pKERejClusters;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    // True abs Np
    TH1D* hTrueAbsNpKE          = new TH1D("hTrueAbsNpKE", "hTrueAbsNpKE;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbsNpKEAsAbs0p   = new TH1D("hTrueAbsNpKEAsAbs0p", "hTrueAbsNpKEAsAbs0p;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbsNpKEAsAbsNp   = new TH1D("hTrueAbsNpKEAsAbsNp", "hTrueAbsNpKEAsAbsNp;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbsNpKEAsScatter = new TH1D("hTrueAbsNpKEAsScatter", "hTrueAbsNpKEAsScatter;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbsNpKEAsChExch  = new TH1D("hTrueAbsNpKEAsChExch", "hTrueAbsNpKEAsChExch;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbsNpKERejected  = new TH1D("hTrueAbsNpKERejected", "hTrueAbsNpKERejected;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    // For each TRUE energy bin, what are abs Np events reconstructed as?
    std::vector<std::vector<TH1*>> TrueAbsNpAsByBin;
    for (int iEnergyBin = 0; iEnergyBin < NUM_BINS_KE; ++iEnergyBin) {
        std::vector<TH1*> TempVec;
        for (int iInt = 0; iInt < NUM_SIGNAL_TYPES; ++iInt) {
            TH1* hTempHist = new TH1D(Form("hTrueAbsNp_%d_Bin_As_%d", iEnergyBin, iInt), Form("True Abs Np KE As %d in bin %d", iInt, iEnergyBin), NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
            TempVec.push_back(hTempHist);
        }
        TrueAbsNpAsByBin.push_back(TempVec);
    }

    TH1D* hTrueAbsNpKERejDataProds = new TH1D("hTrueAbsNpKERejDataProds", "hTrueAbsNpKERejDataProds;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbsNpKERejElectron  = new TH1D("hTrueAbsNpKERejElectron", "hTrueAbsNpKERejElectron;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbsNpKERejRedVol    = new TH1D("hTrueAbsNpKERejRedVol", "hTrueAbsNpKERejRedVol;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbsNpKERejPID       = new TH1D("hTrueAbsNpKERejPID", "hTrueAbsNpKERejPID;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbsNpKERejManyPions = new TH1D("hTrueAbsNpKERejManyPions", "hTrueAbsNpKERejManyPions;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbsNpKERejClusters  = new TH1D("hTrueAbsNpKERejClusters", "hTrueAbsNpKERejClusters;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    // True scatter
    TH1D* hTrueScatterKE          = new TH1D("hTrueScatterKE", "hTrueScatterKE;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueScatterKEAsAbs0p   = new TH1D("hTrueScatterKEAsAbs0p", "hTrueScatterKEAsAbs0p;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueScatterKEAsAbsNp   = new TH1D("hTrueScatterKEAsAbsNp", "hTrueScatterKEAsAbsNp;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueScatterKEAsScatter = new TH1D("hTrueScatterKEAsScatter", "hTrueScatterKEAsScatter;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueScatterKEAsChExch  = new TH1D("hTrueScatterKEAsChExch", "hTrueScatterKEAsChExch;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueScatterKERejected  = new TH1D("hTrueScatterKERejected", "hTrueScatterKERejected;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    // For each TRUE energy bin, what are scatter events reconstructed as?
    std::vector<std::vector<TH1*>> TrueScatterAsByBin;
    for (int iEnergyBin = 0; iEnergyBin < NUM_BINS_KE; ++iEnergyBin) {
        std::vector<TH1*> TempVec;
        for (int iInt = 0; iInt < NUM_SIGNAL_TYPES; ++iInt) {
            TH1* hTempHist = new TH1D(Form("hTrueAbsScatter_%d_Bin_As_%d", iEnergyBin, iInt), Form("True Abs Scatter KE As %d in bin %d", iInt, iEnergyBin), NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
            TempVec.push_back(hTempHist);
        }
        TrueScatterAsByBin.push_back(TempVec);
    }

    TH1D* hTrueScatterKERejDataProds = new TH1D("hTrueScatterKERejDataProds", "hTrueScatterKERejDataProds;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueScatterKERejElectron  = new TH1D("hTrueScatterKERejElectron", "hTrueScatterKERejElectron;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueScatterKERejRedVol    = new TH1D("hTrueScatterKERejRedVol", "hTrueScatterKERejRedVol;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueScatterKERejPID       = new TH1D("hTrueScatterKERejPID", "hTrueScatterKERejPID;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueScatterKERejManyPions = new TH1D("hTrueScatterKERejManyPions", "hTrueScatterKERejManyPions;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueScatterKERejClusters  = new TH1D("hTrueScatterKERejClusters", "hTrueScatterKERejClusters;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    // True charge exchange
    TH1D* hTrueChExchKE          = new TH1D("hTrueChExchKE", "hTrueChExchKE;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueChExchKEAsAbs0p   = new TH1D("hTrueChExchKEAsAbs0p", "hTrueChExchKEAsAbs0p;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueChExchKEAsAbsNp   = new TH1D("hTrueChExchKEAsAbsNp", "hTrueChExchKEAsAbsNp;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueChExchKEAsScatter = new TH1D("hTrueChExchKEAsScatter", "hTrueChExchKEAsScatter;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueChExchKEAsChExch  = new TH1D("hTrueChExchKEAsChExch", "hTrueChExchKEAsChExch;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueChExchKERejected  = new TH1D("hTrueChExchKERejected", "hTrueChExchKERejected;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    // For each TRUE energy bin, what are charge exchange events reconstructed as?
    std::vector<std::vector<TH1*>> TrueChExchAsByBin;
    for (int iEnergyBin = 0; iEnergyBin < NUM_BINS_KE; ++iEnergyBin) {
        std::vector<TH1*> TempVec;
        for (int iInt = 0; iInt < NUM_SIGNAL_TYPES; ++iInt) {
            TH1* hTempHist = new TH1D(Form("hTrueChExch_%d_Bin_As_%d", iEnergyBin, iInt), Form("True Ch Exch KE As %d in bin %d", iInt, iEnergyBin), NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
            TempVec.push_back(hTempHist);
        }
        TrueChExchAsByBin.push_back(TempVec);
    }

    TH1D* hTrueChExchKERejDataProds = new TH1D("hTrueChExchKERejDataProds", "hTrueChExchKERejDataProds;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueChExchKERejElectron  = new TH1D("hTrueChExchKERejElectron", "hTrueChExchKERejElectron;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueChExchKERejRedVol    = new TH1D("hTrueChExchKERejRedVol", "hTrueChExchKERejRedVol;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueChExchKERejPID       = new TH1D("hTrueChExchKERejPID", "hTrueChExchKERejPID;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueChExchKERejManyPions = new TH1D("hTrueChExchKERejManyPions", "hTrueChExchKERejManyPions;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueChExchKERejClusters  = new TH1D("hTrueChExchKERejClusters", "hTrueChExchKERejClusters;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    //////////////
    // Matrices //
    //////////////

    NUM_SIGNAL_TYPES = 4; // override, for now

    std::vector<TH1*> TotalEventsHistos = {
        hTrueAbs0pKE, hTrueAbsNpKE, hTrueScatterKE, hTrueChExchKE
    };

    std::vector<std::vector<std::vector<TH1*>>> TrueRecoAsByBin = {
        TrueAbs0pAsByBin, TrueAbsNpAsByBin, TrueScatterAsByBin, TrueChExchAsByBin
    };

    ///////////////////////////////////
    // Data-MC comparison histograms //
    ///////////////////////////////////

    TH1D* hMCTGSmallTracks          = new TH1D("hMCTGSmallTracks", "hMCTGSmallTracks;;", 10, 0, 10);
    TH1D* hMCTracksNearVertex       = new TH1D("hMCTracksNearVertex", "hMCTracksNearVertex;;", 10, 0, 10);
    TH1D* hMCTrackLengthsNearVertex = new TH1D("hMCTrackLengthsNearVertex", "hMCTrackLengthsNearVertex;;", 50, 0, 100);

    //////////////////////
    // Loop over events //
    //////////////////////

    Int_t NumEntries = (Int_t) tree->GetEntries();
    std::cout << "Num entries: " << NumEntries << std::endl;

    for (Int_t i = 0; i < NumEntries; ++i) {
        tree->GetEntry(i);

        // Make script go faster
        // if (i > 10000) break;

        ////////////////////////////////////////
        // Histograms for data-MC comparisons //
        ////////////////////////////////////////

        // For data-MC comparisons
        if (
            WC2TPCtrkID != -99999 &&
            obtainedProbabilities &&
            showerProb < SHOWER_PROB_CUT
        ) {
            int smallTracksComp = 0; int tracksNearVertexComp = 0;
            for (size_t trk_idx = 0; trk_idx < recoBeginX->size(); ++trk_idx) {
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
                    tracksNearVertexComp++;
                    hMCTrackLengthsNearVertex->Fill(trackLength);
                }
                if (trackLength < SMALL_TRACK_LENGTH_CHEX) smallTracksComp++;
            }

            hMCTracksNearVertex->Fill(tracksNearVertexComp);
            if (!isWithinReducedVolume(WC2TPCPrimaryEndX, WC2TPCPrimaryEndY, WC2TPCPrimaryEndZ)) {
                hMCTGSmallTracks->Fill(smallTracksComp);
            }
        }

        // Back to classification algorithm //
        /////////////

        // Set correct primary vertex KE for elastic scattering
        if (backgroundType == 12) {
            truthPrimaryVertexKE = trajectoryInteractionKE;
        }

        // Get true energy bin
        int TrueEnergyBin = getBin(truthPrimaryVertexKE * 1000, LOWER_BOUND_KE, UPPER_BOUND_KE, NUM_BINS_KE);

        // Add true incident KE
        if (validTrueIncidentKE) {
            for (double x : *trueIncidentKEContributions) {
                hTrueIncidentKE->Fill(x);
            }
        }

        // In this script, we replace background types 6 and 12 with background types 13 and 14,
        // because we are interested in contamination across the bins that we are measuring
        if (backgroundType == 12 || (backgroundType == 6 && numVisibleProtons == 0)) {
            backgroundType = 13;
        } else if (backgroundType == 6 && numVisibleProtons > 0) {
            backgroundType = 14;
        }

        // Fill in true interaction energy histograms
        if (backgroundType == 0) {
            hTrueAbs0pKE->Fill(truthPrimaryVertexKE * 1000);
        } else if (backgroundType == 1) {
            hTrueAbsNpKE->Fill(truthPrimaryVertexKE * 1000);
        } else if (backgroundType == 13 || backgroundType == 14) {
            hTrueScatterKE->Fill(truthPrimaryVertexKE * 1000);
        } else if (backgroundType == 7) {
            hTrueChExchKE->Fill(truthPrimaryVertexKE * 1000);
        }

        // If no track matched to wire-chamber, skip
        if (WC2TPCtrkID == -99999) {
            if (backgroundType == 0) {
                hTrueAbs0pKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueAbs0pKERejDataProds->Fill(truthPrimaryVertexKE * 1000);
            } else if (backgroundType == 1) {
                hTrueAbsNpKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueAbsNpKERejDataProds->Fill(truthPrimaryVertexKE * 1000);
            } else if (backgroundType == 13 || backgroundType == 14) {
                hTrueScatterKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueScatterKERejDataProds->Fill(truthPrimaryVertexKE * 1000);
            } else if (backgroundType == 7) {
                hTrueChExchKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueChExchKERejDataProds->Fill(truthPrimaryVertexKE * 1000);
            }
            continue;
        }
        hDataProdsAndWC2TPC->Fill(backgroundType);

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

        // At this point, we want to fill the incident kinetic energy histograms
        double WCKE             = TMath::Sqrt(WCTrackMomentum * WCTrackMomentum + PionMass * PionMass) - PionMass;
        double calculatedEnLoss = energyLossCalculation(); 
        if (isData) {
            double tanThetaCosPhi = TMath::Tan(WCTheta) * TMath::Cos(WCPhi);
            double tanThetaSinPhi = TMath::Tan(WCTheta) * TMath::Sin(WCPhi);
            double den            = TMath::Sqrt(1 + tanThetaCosPhi * tanThetaCosPhi);
            double onTheFlyPz     = WCTrackMomentum / den;
            double onTheFlyPx     = onTheFlyPz * tanThetaSinPhi;
            calculatedEnLoss      = energyLossCalculation(WC4PrimaryX, onTheFlyPx, isData);
        } else { calculatedEnLoss = energyLossCalculation(WC4PrimaryX, trajectoryInitialMomentumX, isData); }
        const double initialKE = WCKE - calculatedEnLoss;

        double energyDeposited = 0;
        for (size_t iDep = 0; iDep < wcMatchDEDX->size(); ++iDep) {
            // If we are past detected breaking point, exit loop
            if (wcMatchZPos->at(iDep) > breakPointZ) break;

            // If larger than threshold, continue
            if (wcMatchDEDX->at(iDep) > HIT_DEDX_THRESHOLD) continue;

            // Else, add to energy deposited so far
            energyDeposited += wcMatchEDep->at(iDep);
            
            // Add to incident KE if inside reduced volume
            if (isWithinReducedVolume(wcMatchXPos->at(iDep), wcMatchYPos->at(iDep), wcMatchZPos->at(iDep))) {
                hIncidentKE->Fill(initialKE - energyDeposited);

                // Background breakdown
                if (truthPrimaryPDG == -211) {
                    hIncidentKEPion->Fill(initialKE - energyDeposited);
                } else if (truthPrimaryPDG == 13) {
                    hIncidentKEMuon->Fill(initialKE - energyDeposited);
                } else if (truthPrimaryPDG == 11) {
                    hIncidentKEElectron->Fill(initialKE - energyDeposited);
                }
            }
        }
        double energyAtVertex = initialKE - energyDeposited;

        // If did not obtain neural network electron shower probabilities, skip 
        if (!obtainedProbabilities) showerProb = 1.;
        if (showerProb >= SHOWER_PROB_CUT){
            if (backgroundType == 0) {
                hTrueAbs0pKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueAbs0pKERejElectron->Fill(truthPrimaryVertexKE * 1000);
            } else if (backgroundType == 1) {
                hTrueAbsNpKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueAbsNpKERejElectron->Fill(truthPrimaryVertexKE * 1000);
            } else if (backgroundType == 13 || backgroundType == 14) {
                hTrueScatterKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueScatterKERejElectron->Fill(truthPrimaryVertexKE * 1000);
            } else if (backgroundType == 7) {
                hTrueChExchKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueChExchKERejElectron->Fill(truthPrimaryVertexKE * 1000);
            }
            continue;
        }
        hNotAnElectron->Fill(backgroundType);

        if (!isWithinReducedVolume(breakPointX, breakPointY, breakPointZ)) {
            if (backgroundType == 0) {
                hTrueAbs0pKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueAbs0pKERejRedVol->Fill(truthPrimaryVertexKE * 1000);
            } else if (backgroundType == 1) {
                hTrueAbsNpKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueAbsNpKERejRedVol->Fill(truthPrimaryVertexKE * 1000);
            } else if (backgroundType == 13 || backgroundType == 14) {
                hTrueScatterKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueScatterKERejRedVol->Fill(truthPrimaryVertexKE * 1000);
            } else if (backgroundType == 7) {
                hTrueChExchKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueChExchKERejRedVol->Fill(truthPrimaryVertexKE * 1000);
            }
            continue;
        }
        hPrimaryInRedVol->Fill(backgroundType);

        if (minChi2 == pionChi2 || minChi2 == protonChi2) {
            if (backgroundType == 0) {
                hTrueAbs0pKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueAbs0pKERejPID->Fill(truthPrimaryVertexKE * 1000);
            } else if (backgroundType == 1) {
                hTrueAbsNpKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueAbsNpKERejPID->Fill(truthPrimaryVertexKE * 1000);
            } else if (backgroundType == 13 || backgroundType == 14) {
                hTrueScatterKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueScatterKERejPID->Fill(truthPrimaryVertexKE * 1000);
            } else if (backgroundType == 7) {
                hTrueChExchKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueChExchKERejPID->Fill(truthPrimaryVertexKE * 1000);
            }
            continue;
        }
        hPrimaryPID->Fill(backgroundType);

        /////////////////////////
        // Secondary track PID //
        /////////////////////////

        int secondaryTaggedPion   = 0;
        int secondaryTaggedProton = 0;
        int secondaryTaggedOther  = 0;

        int otherTaggedPion   = 0;
        int otherTaggedProton = 0;

        int numSmallTracks = 0;
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

            
            double thisTrackLength = sqrt(
                pow(recoBeginX->at(trk_idx) - recoEndX->at(trk_idx), 2) +
                pow(recoBeginY->at(trk_idx) - recoEndY->at(trk_idx), 2) + 
                pow(recoBeginZ->at(trk_idx) - recoEndZ->at(trk_idx), 2)
            );
            if (thisTrackLength < SMALL_TRACK_LENGTH_CHEX) numSmallTracks++;

            if ((distanceFromStart < VERTEX_RADIUS || distanceFromEnd < VERTEX_RADIUS)) {
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

        // For particles where we stitched, we also need to analyze the second part of the primary track
        if (minChi2 == minStitchedChi2) {
            std::vector<double> newSecondaryResR(wcMatchResR->begin(), wcMatchResR->begin() + bestBreakPoint);
            std::vector<double> newSecondaryDEDX(wcMatchDEDX->begin(), wcMatchDEDX->begin() + bestBreakPoint);

            double newPionChi2   = computeReducedChi2(gPion, newSecondaryResR, newSecondaryDEDX, false, newSecondaryDEDX.size(), nRemoveOutliers, nRemoveEnds);
            double newProtonChi2 = computeReducedChi2(gProton, newSecondaryResR, newSecondaryDEDX, false, newSecondaryDEDX.size(), nRemoveOutliers, nRemoveEnds);
            double newMeanDEDX   = meanDEDX(newSecondaryDEDX, false, MEAN_DEDX_NUM_TRAJ_POINTS);

            if ((newPionChi2 < PION_CHI2_PION_VALUE) && (newProtonChi2 > PROTON_CHI2_PION_VALUE)) {
                // Tagged as pion
                secondaryTaggedPion++;
            } else if ((newPionChi2 > PION_CHI2_PROTON_VALUE) && (newProtonChi2 < PROTON_CHI2_PROTON_VALUE)) {
                // Tagged as proton
                secondaryTaggedProton++;
            } else {
                // Not tagged with chi^2, use mean dE/dx
                secondaryTaggedOther++;
                if (newMeanDEDX <= MEAN_DEDX_THRESHOLD) {
                    otherTaggedPion++;
                } else {
                    otherTaggedProton++;
                }
            }
        }

        int totalTaggedPions   = secondaryTaggedPion + otherTaggedPion;
        int totalTaggedProtons = secondaryTaggedProton + otherTaggedProton;

        if (totalTaggedPions > 0) {
            if (totalTaggedPions > 1) {
                // reject events with > 1 tagged pion
                if (backgroundType == 0) {
                    hTrueAbs0pKERejected->Fill(truthPrimaryVertexKE * 1000);
                    hTrueAbs0pKERejManyPions->Fill(truthPrimaryVertexKE * 1000);
                } else if (backgroundType == 1) {
                    hTrueAbsNpKERejected->Fill(truthPrimaryVertexKE * 1000);
                    hTrueAbsNpKERejManyPions->Fill(truthPrimaryVertexKE * 1000);
                } else if (backgroundType == 13 || backgroundType == 14) {
                    hTrueScatterKERejected->Fill(truthPrimaryVertexKE * 1000);
                    hTrueScatterKERejManyPions->Fill(truthPrimaryVertexKE * 1000);
                } else if (backgroundType == 7) {
                    hTrueChExchKERejected->Fill(truthPrimaryVertexKE * 1000);
                    hTrueChExchKERejManyPions->Fill(truthPrimaryVertexKE * 1000);
                }
                continue;
            }

            // Select as scatter
            hPionScatter->Fill(backgroundType);

            hPionScatterKE->Fill(energyAtVertex);
            if (backgroundType == 0) {
                hPionScatterKEAbs0p->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                hPionScatterKEAbsNp->Fill(energyAtVertex);
            } else if (backgroundType == 7) {
                hPionScatterKEChExch->Fill(energyAtVertex);
            } else if (backgroundType == 13 || backgroundType == 14) {
                hPionScatterKETrue->Fill(energyAtVertex);
            } else if (backgroundType == 2) {
                hPionScatterKEMuon->Fill(energyAtVertex);
            } else if (backgroundType == 3) {
                hPionScatterKEElectron->Fill(energyAtVertex);
            } else {
                hPionScatterKEOther->Fill(energyAtVertex);
            }

            if (backgroundType == 0) {
                hTrueAbs0pKEAsScatter->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueAbs0pAsByBin.at(TrueEnergyBin).at(2)->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                hTrueAbsNpKEAsScatter->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueAbsNpAsByBin.at(TrueEnergyBin).at(2)->Fill(energyAtVertex);
            } else if (backgroundType == 13 || backgroundType == 14) {
                hTrueScatterKEAsScatter->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueScatterAsByBin.at(TrueEnergyBin).at(2)->Fill(energyAtVertex);
            } else if (backgroundType == 7) {
                hTrueChExchKEAsScatter->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueChExchAsByBin.at(TrueEnergyBin).at(2)->Fill(energyAtVertex);
            }

            // TODO: maybe check not rejecting more than 1 pion

            continue;
        }
        hNotScatter->Fill(backgroundType);

        if (totalTaggedProtons > 0) {
            // Select as Np absorption
            hPionAbsNp->Fill(backgroundType);

            hPionAbsNpKE->Fill(energyAtVertex);
            if (backgroundType == 0) {
                hPionAbsNpKEAbs0p->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                hPionAbsNpKETrue->Fill(energyAtVertex);
            } else if (backgroundType == 7) {
                hPionAbsNpKEChExch->Fill(energyAtVertex);
            } else if (backgroundType == 13 || backgroundType == 14) {
                hPionAbsNpKEScatter->Fill(energyAtVertex);
            } else if (backgroundType == 2) {
                hPionAbsNpKEMuon->Fill(energyAtVertex);
            } else if (backgroundType == 3) {
                hPionAbsNpKEElectron->Fill(energyAtVertex);
            } else {
                hPionAbsNpKEOther->Fill(energyAtVertex);
            }

            if (backgroundType == 0) {
                hTrueAbs0pKEAsAbsNp->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueAbs0pAsByBin.at(TrueEnergyBin).at(1)->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                hTrueAbsNpKEAsAbsNp->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueAbsNpAsByBin.at(TrueEnergyBin).at(1)->Fill(energyAtVertex);
            } else if (backgroundType == 13 || backgroundType == 14) {
                hTrueScatterKEAsAbsNp->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueScatterAsByBin.at(TrueEnergyBin).at(1)->Fill(energyAtVertex);
            } else if (backgroundType == 7) {
                hTrueChExchKEAsAbsNp->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueChExchAsByBin.at(TrueEnergyBin).at(1)->Fill(energyAtVertex);
            }

            continue;
        }
        hNotPionAbsNp->Fill(backgroundType);

        /////////////////////
        // Small track cut //
        /////////////////////

        if (numSmallTracks > SMALL_TRACK_CUT_CHEX) {
            // Tag as charge exchange
            hPionChExch->Fill(backgroundType);

            hPionChExchKE->Fill(energyAtVertex);
            if (backgroundType == 0) {
                hPionChExchKEAbs0p->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                hPionChExchKEAbsNp->Fill(energyAtVertex);
            } else if (backgroundType == 7) {
                hPionChExchKETrue->Fill(energyAtVertex);
            } else if (backgroundType == 13 || backgroundType == 14) {
                hPionChExchKEScatter->Fill(energyAtVertex);
            } else if (backgroundType == 2) {
                hPionChExchKEMuon->Fill(energyAtVertex);
            } else if (backgroundType == 3) {
                hPionChExchKEElectron->Fill(energyAtVertex);
            } else {
                hPionChExchKEOther->Fill(energyAtVertex);
            }

            if (backgroundType == 0) {
                hTrueAbs0pKEAsChExch->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueAbs0pAsByBin.at(TrueEnergyBin).at(3)->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                hTrueAbsNpKEAsChExch->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueAbsNpAsByBin.at(TrueEnergyBin).at(3)->Fill(energyAtVertex);
            } else if (backgroundType == 13 || backgroundType == 14) {
                hTrueScatterKEAsChExch->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueScatterAsByBin.at(TrueEnergyBin).at(3)->Fill(energyAtVertex);
            } else if (backgroundType == 7) {
                hTrueChExchKEAsChExch->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueChExchAsByBin.at(TrueEnergyBin).at(3)->Fill(energyAtVertex);
            }

            outFileChExchBkg << "Run: " << run << " subrun: " << subrun << " event: " << event << std::endl;
            outFileChExchBkg << " Classified as: " << backgroundTypes[backgroundType] << " with energy at vertex: " << truthPrimaryVertexKE * 1000 << std::endl;
            outFileChExchBkg << std::endl;

            continue;
        }
        hNotChExch->Fill(backgroundType);

        /////////////////////////////////////////
        // Selection on non-reconstructed hits //
        /////////////////////////////////////////

        // First, find hits in the induction plane near the primary track that
        // were not reconstructed into any of the already existing tracks
        const float xThreshold = 2.0;
        const float wThreshold = 5.0 * HIT_WIRE_SEPARATION;
        std::unordered_set<int> hitsInTracks(hitRecoAsTrackKey->begin(), hitRecoAsTrackKey->end());
        std::vector<int> candidateInductionHits;

        int nTotalHits = fHitKey->size();
        for (size_t iHit = 0; iHit < nTotalHits; ++iHit) {
            // Skip hits already in tracks and collection plane
            // if ((hitsInTracks.count(iHit) > 0) || (fHitPlane->at(iHit) == 1)) continue;
            if (hitsInTracks.count(iHit) > 0) continue;

            float hitX = fHitX->at(iHit);
            float hitW = fHitW->at(iHit);

            if (isHitNearPrimary(
                hitWC2TPCKey, 
                fHitX, 
                fHitW, 
                hitX, 
                hitW,
                xThreshold, 
                wThreshold
            )) { candidateInductionHits.push_back(iHit); }
        }

        // Now cluster using those hits as starting points
        std::unordered_set<int> usedHits;
        std::vector<HitCluster> hitClusters;
        int nCandidateHits    = candidateInductionHits.size();

        // Hits in the same cluster must be separated by at most some number of wires
        float maxHitClusterSeparation = HIT_WIRE_SEPARATION * 3.0; 
        float maxHitXSeparation       = 1.0;
        
        for (int iHit = 0; iHit < nCandidateHits; ++iHit) {
            // The candidate hits are only STARTING points, as this could make up 
            // really long tracks in the induction plane that are no longer near 
            // the ending point of the primary track
            
            // First, check if we have already used this hit
            int   thisHitKey       = candidateInductionHits.at(iHit);
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
            
            for (int iAllHit = 0; iAllHit < nTotalHits; ++iAllHit) {
                // Skip already used hits, those reconstructed in tracks, and those in collection plane
                // if (usedHits.count(iAllHit) || hitsInTracks.count(iAllHit) || (fHitPlane->at(iHit) != 0)) continue;
                if (usedHits.count(iAllHit) || hitsInTracks.count(iAllHit)) continue;
                float internalHitW  = fHitW->at(iAllHit);
                float internalHitX  = fHitX->at(iAllHit);
                float dW            = std::abs(internalHitW - thisHitW);
                float dX            = std::abs(internalHitX - thisHitX);

                int   nClusterSoFar = clusterW.size();
                for (int iCluster = 0; iCluster < nClusterSoFar; ++iCluster) {
                    float tempdW = std::abs(internalHitW - clusterW.at(iCluster));
                    if (tempdW < dW) dW = tempdW;
                    float tempdX = std::abs(internalHitX - clusterX.at(iCluster));
                    if (tempdX < dX) dX = tempdX;
                }

                if ((dW < maxHitClusterSeparation) && (dX < maxHitXSeparation)) {
                    usedHits.insert(iAllHit);
                    clusterKeys.push_back(iAllHit);
                    clusterX.push_back(fHitX->at(iAllHit));
                    clusterW.push_back(fHitW->at(iAllHit));
                    clusterCharge.push_back(fHitCharge->at(iAllHit));
                    clusterChargeCol.push_back(fHitChargeCol->at(iAllHit));
                }
            }

            HitCluster thisCluster;
            thisCluster.hitKeys      = clusterKeys;
            thisCluster.hitX         = clusterX;
            thisCluster.hitW         = clusterW;
            thisCluster.hitCharge    = clusterCharge;
            thisCluster.hitChargeCol = clusterChargeCol;
            
            usedHits.insert(thisHitKey);
            hitClusters.push_back(thisCluster);
        }

        // Get data for cut
        int   numLargeClusters = 0;
        for (int i = 0; i < hitClusters.size(); ++i) {
            float maxWire = *std::max_element(hitClusters[i].hitW.begin(), hitClusters[i].hitW.end());
            float minWire = *std::min_element(hitClusters[i].hitW.begin(), hitClusters[i].hitW.end());
            if (std::abs(maxWire - minWire) > LARGE_CLUSTER_THRESHOLD) numLargeClusters++;
        }

        if (numLargeClusters < NUM_CLUSTERS_THRESHOLD) {
            hPionAbs0p->Fill(backgroundType);

            hPionAbs0pKE->Fill(energyAtVertex);
            if (backgroundType == 0) {
                hPionAbs0pKETrue->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                hPionAbs0pKEAbsNp->Fill(energyAtVertex);
            } else if (backgroundType == 7) {
                hPionAbs0pKEChExch->Fill(energyAtVertex);
            } else if (backgroundType == 13 || backgroundType == 14) {
                hPionAbs0pKEScatter->Fill(energyAtVertex);
            } else if (backgroundType == 2) {
                hPionAbs0pKEMuon->Fill(energyAtVertex);
            } else if (backgroundType == 3) {
                hPionAbs0pKEElectron->Fill(energyAtVertex);
            } else {
                hPionAbs0pKEOther->Fill(energyAtVertex);
            }

            if (backgroundType == 0) {
                hTrueAbs0pKEAsAbs0p->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueAbs0pAsByBin.at(TrueEnergyBin).at(0)->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                hTrueAbsNpKEAsAbs0p->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueAbsNpAsByBin.at(TrueEnergyBin).at(0)->Fill(energyAtVertex);
            } else if (backgroundType == 13 || backgroundType == 14) {
                hTrueScatterKEAsAbs0p->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueScatterAsByBin.at(TrueEnergyBin).at(0)->Fill(energyAtVertex);
            } else if (backgroundType == 7) {
                hTrueChExchKEAsAbs0p->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueChExchAsByBin.at(TrueEnergyBin).at(0)->Fill(energyAtVertex);
            }

            continue;
        }
        hNotPionAbs0p->Fill(backgroundType);

        // Anything left here is rejected
        if (backgroundType == 0) {
            hTrueAbs0pKERejected->Fill(truthPrimaryVertexKE * 1000);
            hTrueAbs0pKERejClusters->Fill(truthPrimaryVertexKE * 1000);
        } else if (backgroundType == 1) {
            hTrueAbsNpKERejected->Fill(truthPrimaryVertexKE * 1000);
            hTrueAbsNpKERejClusters->Fill(truthPrimaryVertexKE * 1000);
        } else if (backgroundType == 13 || backgroundType == 14) {
            hTrueScatterKERejected->Fill(truthPrimaryVertexKE * 1000);
            hTrueScatterKERejClusters->Fill(truthPrimaryVertexKE * 1000);
        } else if (backgroundType == 7) {
            hTrueChExchKERejected->Fill(truthPrimaryVertexKE * 1000);
            hTrueChExchKERejClusters->Fill(truthPrimaryVertexKE * 1000);
        }
    }

    std::cout << std::endl;
    std::cout << "Original sample composition: " << std::endl;
    printBackgroundInfo(hTotalEvents, std::cout);

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
    std::cout << "Events not classified as electron: " << std::endl;
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
    std::cout << "Events not classified as charge exchange (not enough small tracks): " << std::endl;
    printBackgroundInfo(hNotChExch, std::cout);

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
    std::cout << "Pion abs 0p total reco " <<  hPionAbs0p->Integral() << " with composition: " << std::endl;
    printBackgroundInfo(hPionAbs0p, std::cout);
    std::cout << "Purity: " << hPionAbs0p->GetBinContent(1) / hPionAbs0p->Integral() << std::endl;
    std::cout << "Efficiency: " << hPionAbs0p->GetBinContent(1) / hTotalEvents->GetBinContent(1) << std::endl;

    std::cout << std::endl;
    std::cout << "Pion abs Np total reco " <<  hPionAbsNp->Integral() << " with composition: " << std::endl;
    printBackgroundInfo(hPionAbsNp, std::cout);
    std::cout << "Purity: " << hPionAbsNp->GetBinContent(2) / hPionAbsNp->Integral() << std::endl;
    std::cout << "Efficiency: " << hPionAbsNp->GetBinContent(2) / hTotalEvents->GetBinContent(2) << std::endl;

    std::cout << std::endl;
    std::cout << "Pion scattering total reco " << hPionScatter->Integral() << " with composition: " << std::endl;
    printBackgroundInfo(hPionScatter, std::cout);
    std::cout << "Purity: " << (hPionScatter->GetBinContent(14) + hPionScatter->GetBinContent(15)) / hPionScatter->Integral() << std::endl;
    std::cout << "Efficiency: " << (hPionScatter->GetBinContent(14) + hPionScatter->GetBinContent(15)) / (hTotalEvents->GetBinContent(14) + hTotalEvents->GetBinContent(15)) << std::endl;

    std::cout << std::endl;
    std::cout << "Pion charge exchange total reco " << hPionChExch->Integral() << " with composition: " << std::endl;
    printBackgroundInfo(hPionChExch, std::cout);
    std::cout << "Purity: " << hPionChExch->GetBinContent(8) / hPionChExch->Integral() << std::endl;
    std::cout << "Efficiency: " << hPionChExch->GetBinContent(8) / hTotalEvents->GetBinContent(8) << std::endl;

    ////////////////////////////////
    // Save comparison histograms //
    ////////////////////////////////

    hMCTGSmallTracks->SetDirectory(comparisonsFile);
    hMCTGSmallTracks->Write();

    hMCTracksNearVertex->SetDirectory(comparisonsFile);
    hMCTracksNearVertex->Write();

    hMCTrackLengthsNearVertex->SetDirectory(comparisonsFile);
    hMCTrackLengthsNearVertex->Write();

    //////////////////////////////////////////////
    // Perform unfolding for interacting slices //
    //////////////////////////////////////////////

    // Response matrix with components: 
    //    R[(i, \alpha), (j, \beta)] = P(reco signal \beta in reco energy bin j | true signal \alpha in true energy bin i)
    TH2D* hResponseMatrix = new TH2D(
        "hResponseMatrix", "hResponseMatrix;Reco (j, #beta);True (i, #alpha)",
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE,
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    );

    // Construct response matrix
    for (int iOuterSignalBin = 0; iOuterSignalBin < NUM_SIGNAL_TYPES; ++iOuterSignalBin) {
        for (int iOuterEnergyBin = 0; iOuterEnergyBin < NUM_BINS_KE; ++iOuterEnergyBin) {
            // Get index for row (i, \alpha) -> k
            int row = flattenIndex(iOuterSignalBin, iOuterEnergyBin, NUM_BINS_KE);

            // Get total number of true events of true signal type \alpha with true energy bin i
            double denom = TotalEventsHistos.at(iOuterSignalBin)->GetBinContent(iOuterEnergyBin + 1);

            for (int iInnerSignalBin = 0; iInnerSignalBin < NUM_SIGNAL_TYPES; ++iInnerSignalBin) {
                for (int iInnerEnergyBin = 0; iInnerEnergyBin < NUM_BINS_KE; ++iInnerEnergyBin) {
                    // Get index for column (j, \beta) -> l
                    int  column = flattenIndex(iInnerSignalBin, iInnerEnergyBin, NUM_BINS_KE);
                    double prob = 0;
                    
                    // Get total number of true (i, \alpha) events that were reconstructed as (j, \beta)
                    if (denom > 0) prob = TrueRecoAsByBin.at(iOuterSignalBin).at(iOuterEnergyBin).at(iInnerSignalBin)->GetBinContent(iInnerEnergyBin + 1) / denom;
                    hResponseMatrix->SetBinContent(column + 1, row + 1, prob);
                }
            }
        }
    }

    // Before unfolding, we have to remove the estimated backgrounds using ratios
    std::vector<TH1D*> hSubtractEstimatedBkgs(NUM_SIGNAL_TYPES);
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        hSubtractEstimatedBkgs[iSignal] = (TH1D*) RecoSignals[iSignal]->Clone(Form("hSubtractEstimatedBkgs_%d", iSignal));
        hSubtractEstimatedBkgs[iSignal]->SetTitle(Form("Ratio for non-bkgs for signal %d", iSignal));
        for (int iBin = 1; iBin <= NUM_BINS_KE; ++iBin) {
            double total = RecoSignals[iSignal]->GetBinContent(iBin);

            double bkg = 0;
            for (int iBkg = 0; iBkg < RecoSignalBackgrounds[iSignal].size(); ++iBkg) {
                bkg += RecoSignalBackgrounds[iSignal][iBkg]->GetBinContent(iBin);
            }

            double ratio = (total > 0) ? (total - bkg) / total : 0;
            hSubtractEstimatedBkgs[iSignal]->SetBinContent(iBin, ratio);
        }

        // Apply ratio to reco signal
        for (int iBin = 1; iBin <= NUM_BINS_KE; ++iBin) {
            double reco  = RecoSignals[iSignal]->GetBinContent(iBin);
            double ratio = hSubtractEstimatedBkgs[iSignal]->GetBinContent(iBin);
            RecoSignals[iSignal]->SetBinContent(iBin, reco * ratio);
        }
    }

    // Construct large measured vector
    TVectorD Measure(NUM_SIGNAL_TYPES * NUM_BINS_KE);
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index      = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            Measure(index) = RecoSignals[iSignal]->GetBinContent(iBin + 1);
        }
    }

    // Construct large signal vector
    TVectorD Signal(NUM_SIGNAL_TYPES * NUM_BINS_KE);
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index     = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            Signal(index) = TotalEventsHistos[iSignal]->GetBinContent(iBin + 1);
        }
    }

    // Construct covariance matrix (only statistical variance for now)
    TMatrixD Covariance(NUM_SIGNAL_TYPES * NUM_BINS_KE, NUM_SIGNAL_TYPES * NUM_BINS_KE); Covariance.Zero();
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            double N    = RecoSignals[iSignal]->GetBinContent(iBin + 1);
            double Ninc = hIncidentKE->GetBinContent(iBin + 1);

            if (N > 0) {
                Covariance(index, index) = N * (1 - (N / Ninc)); // variance = N_int(1-N_int/N_inc)
            } else {
                Covariance(index, index) = 1.0; // small floor to avoid zero variance
            }
        }
    }

    // Convert histo to matrix
    TMatrixD Response(NUM_SIGNAL_TYPES * NUM_BINS_KE, NUM_SIGNAL_TYPES * NUM_BINS_KE); H2M(static_cast<const TH2D*>(hResponseMatrix), Response, kTRUE);

    // Objects to store stuff
    TMatrixD AddSmear(NUM_SIGNAL_TYPES * NUM_BINS_KE, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    TVectorD WF(NUM_SIGNAL_TYPES * NUM_BINS_KE);
    TMatrixD UnfoldCov(NUM_SIGNAL_TYPES * NUM_BINS_KE, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    TMatrixD CovRotation(NUM_SIGNAL_TYPES * NUM_BINS_KE, NUM_SIGNAL_TYPES * NUM_BINS_KE);

    TVectorD UnfoldedReco = WienerSVD(
        Response,
        Signal,
        Measure,
        Covariance,
        2,
        0.5,
        AddSmear,
        WF,
        UnfoldCov,
        CovRotation
    );
    TVectorD SmearedTrue = AddSmear * Signal;

    // Organize results
    std::vector<TH1*> UnfoldedRecoHistos;
    std::vector<TH1*> SmearedTrueHistos;
    for (int iBin = 0; iBin < NUM_SIGNAL_TYPES; ++iBin) {
        TH1D* unfoldedHist = new TH1D(Form("hUnfoldedRecoVec_%d", iBin), Form("Unfolded Reco Vec %d", iBin), NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
        TH1D* smearedHist = new TH1D(Form("hSmearedTrueVec_%d", iBin), Form("Smeared True Vec %d", iBin), NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
        UnfoldedRecoHistos.push_back(unfoldedHist);
        SmearedTrueHistos.push_back(smearedHist);
    }

    for (int iFlatIndex = 0; iFlatIndex < NUM_SIGNAL_TYPES * NUM_BINS_KE; ++iFlatIndex) {
        auto [signalBin, energyBin] = unflattenIndex(iFlatIndex, NUM_BINS_KE);

        // Get error for unfolded reco histogram from unfolded covariance
        double err = std::sqrt(UnfoldCov[iFlatIndex][iFlatIndex]);
        UnfoldedRecoHistos[signalBin]->SetBinContent(energyBin + 1, UnfoldedReco(iFlatIndex));
        UnfoldedRecoHistos[signalBin]->SetBinError(energyBin + 1, err);

        SmearedTrueHistos[signalBin]->SetBinContent(energyBin + 1, SmearedTrue(iFlatIndex));
    }

    // Copy covariance and smearing matrices into histos
    TH2D* hCovariance = new TH2D(
        "hCovariance", "Covariance Matrix;(i, #alpha);(j, #beta)",
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE, 
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE 
    );
    TH2D* hSmearing = new TH2D("hSmearing", "Smearing Matrix;True (i, #alpha); Reco (j, #beta)", 
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE, 
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    );
    M2H(Covariance, hCovariance); M2H(AddSmear, hSmearing);

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

    // Only for plotting corrected incident KE
    TH1D* hIncidentKECorrected = (TH1D*) hIncidentKE->Clone("hIncidentKECorrected");
    hIncidentKECorrected->Divide(hPsiInc);
    hIncidentKECorrected->Multiply(hCInc);

    ///////////////////////////////////////////////
    // Get cross-section using corrrected fluxes //
    ///////////////////////////////////////////////

    TH1D* hPionAbs0pCrossSection     = new TH1D("hPionAbs0pCrossSection", "hPionAbs0pCrossSection;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTruePionAbs0pCrossSection = new TH1D("hTruePionAbs0pCrossSection", "hTruePionAbs0pCrossSection;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    TH1D* hPionAbsNpCrossSection     = new TH1D("hPionAbsNpCrossSection", "hPionAbsNpCrossSection;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTruePionAbsNpCrossSection = new TH1D("hTruePionAbsNpCrossSection", "hTruePionAbsNpCrossSection;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    TH1D* hPionScatterCrossSection     = new TH1D("hPionScatterCrossSection", "hPionScatterCrossSection;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTruePionScatterCrossSection = new TH1D("hTruePionScatterCrossSection", "hTruePionScatterCrossSection;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    TH1D* hPionChExchCrossSection     = new TH1D("hPionChExchCrossSection", "hPionChExchCrossSection;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTruePionChExchCrossSection = new TH1D("hTruePionChExchCrossSection", "hTruePionChExchCrossSection;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    std::vector<TH1D*> UnfoldedCrossSections = {
        hPionAbs0pCrossSection,
        hPionAbsNpCrossSection,
        hPionScatterCrossSection,
        hPionChExchCrossSection
    };

    std::vector<TH1D*> TrueCrossSections = {
        hTruePionAbs0pCrossSection,
        hTruePionAbsNpCrossSection,
        hTruePionScatterCrossSection,
        hTruePionChExchCrossSection
    };

    for (int i = 0; i < UnfoldedCrossSections.size(); ++i) {
        TH1D* unfXSec  = UnfoldedCrossSections[i];
        TH1D* trueXSec = TrueCrossSections[i];

        for (int iBin = 1; iBin <= NUM_BINS_KE; ++iBin) {
            double incidentErr     = hIncidentKE->GetBinError(iBin);
            double incidentContent = hIncidentKE->GetBinContent(iBin);

            double interactingErr     = UnfoldedRecoHistos[i]->GetBinError(iBin);
            double interactingContent = UnfoldedRecoHistos[i]->GetBinContent(iBin);

            // we compute error as δσ / σ = δN_int / N_int + δN_inc / N_inc (conservative upper bound)
            double crossSection    = interactingContent / incidentContent;
            double crossSectionErr = crossSection * ((interactingErr / interactingContent) + (incidentErr / incidentContent));
            unfXSec->SetBinContent(iBin, crossSection);
            unfXSec->SetBinError(iBin, crossSectionErr);

            // True cross-section, no error bars
            trueXSec->SetBinContent(iBin, SmearedTrueHistos[i]->GetBinContent(iBin) / hTrueIncidentKE->GetBinContent(iBin));
        }

        // Incident KE corrections
        unfXSec->Multiply(hPsiInc);
        unfXSec->Divide(hCInc);

        // Scale units
        unfXSec->Scale(1.0 / (number_density * slab_width * 1e-28));
        trueXSec->Scale(1.0 / (number_density * slab_width * 1e-28));
    }

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
        // Incident KE
        {hIncidentKEPion, hIncidentKEMuon, hIncidentKEElectron},
        {hIncidentKECorrected, hTrueIncidentKE},

        // Interacting KE
        {hPionAbs0pKETrue, hPionAbs0pKEAbsNp, hPionAbs0pKEScatter, hPionAbs0pKEChExch, hPionAbs0pKEMuon, hPionAbs0pKEElectron, hPionAbs0pKEOther},
        {hPionAbsNpKETrue, hPionAbsNpKEAbs0p, hPionAbsNpKEScatter, hPionAbsNpKEChExch, hPionAbsNpKEMuon, hPionAbsNpKEElectron, hPionAbsNpKEOther},
        {hPionScatterKETrue, hPionScatterKEAbs0p, hPionScatterKEAbsNp, hPionScatterKEChExch, hPionScatterKEMuon, hPionScatterKEElectron, hPionScatterKEOther},
        {hPionChExchKETrue, hPionChExchKEAbs0p, hPionChExchKEAbsNp, hPionChExchKEScatter, hPionChExchKEMuon, hPionChExchKEElectron, hPionChExchKEOther},

        // True events classified breakdown
        {hTrueAbs0pKEAsAbs0p, hTrueAbs0pKEAsAbsNp, hTrueAbs0pKEAsScatter, hTrueAbs0pKEAsChExch, hTrueAbs0pKERejected},
        {hTrueAbsNpKEAsAbs0p, hTrueAbsNpKEAsAbsNp, hTrueAbsNpKEAsScatter, hTrueAbsNpKEAsChExch, hTrueAbsNpKERejected},
        {hTrueScatterKEAsAbs0p, hTrueScatterKEAsAbsNp, hTrueScatterKEAsScatter, hTrueScatterKEAsChExch, hTrueScatterKERejected},
        {hTrueChExchKEAsAbs0p, hTrueChExchKEAsAbsNp, hTrueChExchKEAsScatter, hTrueChExchKEAsChExch, hTrueChExchKERejected},

        // Rejected events
        {hTrueAbs0pKERejDataProds, hTrueAbs0pKERejElectron, hTrueAbs0pKERejRedVol, hTrueAbs0pKERejPID, hTrueAbs0pKERejManyPions, hTrueAbs0pKERejClusters},
        {hTrueAbsNpKERejDataProds, hTrueAbsNpKERejElectron, hTrueAbsNpKERejRedVol, hTrueAbsNpKERejPID, hTrueAbsNpKERejManyPions, hTrueAbsNpKERejClusters},
        {hTrueScatterKERejDataProds, hTrueScatterKERejElectron, hTrueScatterKERejRedVol, hTrueScatterKERejPID, hTrueScatterKERejManyPions, hTrueScatterKERejClusters},
        {hTrueChExchKERejDataProds, hTrueChExchKERejElectron, hTrueChExchKERejRedVol, hTrueChExchKERejPID, hTrueChExchKERejManyPions, hTrueChExchKERejClusters},

        // Cross-sections
        {hTruePionAbs0pCrossSection, hPionAbs0pCrossSection},
        {hTruePionAbsNpCrossSection, hPionAbsNpCrossSection},
        {hTruePionScatterCrossSection, hPionScatterCrossSection},
        {hTruePionChExchCrossSection, hPionChExchCrossSection}
    };

    std::vector<std::vector<TString>> PlotLabelGroups = {
        // Incident KE
        {"Pions", "Muons", "Electrons"},
        {"Corrected", "True"},

        // Interacting KE
        {"True", "Abs Np", "Scatter", "Ch. exch.", "Muon", "Electron", "Other"},
        {"True", "Abs Np", "Scatter", "Ch. exch.", "Muon", "Electron", "Other"},
        {"True", "Abs 0p", "Abs Np", "Ch. exch.", "Muon", "Electron", "Other"},
        {"True", "Abs 0p", "Abs Np", "Scatter", "Muon", "Electron", "Other"},

        // True events classified breakdown
        {"Abs 0p", "Abs Np", "Scatter", "Ch. exch.", "Rejected"},
        {"Abs 0p", "Abs Np", "Scatter", "Ch. exch.", "Rejected"},
        {"Abs 0p", "Abs Np", "Scatter", "Ch. exch.", "Rejected"},
        {"Abs 0p", "Abs Np", "Scatter", "Ch. exch.", "Rejected"},

        // Rejected events
        {"Data-prods", "Shower-like", "Red. vol.", "PID reject", "> 1 pion", "Hit clusters"},
        {"Data-prods", "Shower-like", "Red. vol.", "PID reject", "> 1 pion", "Hit clusters"},
        {"Data-prods", "Shower-like", "Red. vol.", "PID reject", "> 1 pion", "Hit clusters"},
        {"Data-prods", "Shower-like", "Red. vol.", "PID reject", "> 1 pion", "Hit clusters"},

        // Cross-sections
        {"True", "Unfolded"},
        {"True", "Unfolded"},
        {"True", "Unfolded"},
        {"True", "Unfolded"}
    };

    std::vector<TString> PlotTitles = {
        // Incident KE
        "Incident/IncidentKE",
        "Incident/IncidentKECorrected",

        // Interacting KE
        "RecoInteracting/Abs0pInteractingKE",
        "RecoInteracting/AbsNpInteractingKE",
        "RecoInteracting/ScatterInteractingKE",
        "RecoInteracting/ChExchInteractingKE",

        // True events classified breakdown
        "TrueInteracting/TrueAbs0pBreakdown",
        "TrueInteracting/TrueAbsNpBreakdown",
        "TrueInteracting/TrueScatterBreakdown",
        "TrueInteracting/TrueChExchBreakdown",

        // Rejected events
        "Rejected/TrueAbs0pRej",
        "Rejected/TrueAbsNpRej",
        "Rejected/TrueScatterRej",
        "Rejected/TrueChExchRej",

        // Cross-sections
        "CrossSection/PionAbs0pCrossSection",
        "CrossSection/PionAbsNpCrossSection",
        "CrossSection/PionScatterCrossSection",
        "CrossSection/PionChExchCrossSection"
    };

    std::vector<TString> XLabels = {
        // Incident KE
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",

        // Interacting KE
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",

        // True events classified breakdown
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",

        // Rejected events
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",

        // Cross-sections
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]"
    };

    std::vector<TString> YLabels = {
        // Incident KE
        "Counts",
        "Counts",

        // Interacting KE
        "Counts",
        "Counts",
        "Counts",
        "Counts",

        // True events classified breakdown
        "Counts",
        "Counts",
        "Counts",
        "Counts",

        // Rejected events
        "Counts",
        "Counts",
        "Counts",
        "Counts",

        // Cross-sections
        "Cross section [barn]",
        "Cross section [barn]",
        "Cross section [barn]",
        "Cross section [barn]"
    };

    std::vector<bool> PlotStacked = {
        // Incident KE
        true,
        false,

        // Interacting KE
        true,
        true,
        true,
        true,

        // True events classified breakdown
        true,
        true,
        true,
        true,

        // Rejected events
        true,
        true,
        true,
        true,

        // Cross-sections
        false,
        false,
        false,
        false
    };

    std::vector<std::vector<bool>> PlotsAsPoints = {
        // Incident KE
        {false, false, false},
        {true, false},

        // Interacting KE
        {false, false, false, false, false, false, false},
        {false, false, false, false, false, false, false},
        {false, false, false, false, false, false, false},
        {false, false, false, false, false, false, false},

        // True events classified breakdown
        {false, false, false, false, false},
        {false, false, false, false, false},
        {false, false, false, false, false},
        {false, false, false, false, false},

        // Rejected events
        {false, false, false, false, false, false},
        {false, false, false, false, false, false},
        {false, false, false, false, false, false},
        {false, false, false, false, false, false},

        // Cross-sections
        {false, true},
        {false, true},
        {false, true},
        {false, true}
    };

    // Add each unfolded histogram as a single plot group for plotting
    std::vector<TString> UnfHistTitles = {
        "UnfoldedAbs0p",
        "UnfoldedAbsNp",
        "UnfoldedScatter",
        "UnfoldedChExch"
    };

    for (size_t i = 0; i < UnfoldedRecoHistos.size(); ++i) {
        PlotGroups.push_back({SmearedTrueHistos[i], UnfoldedRecoHistos[i]});
        PlotsAsPoints.push_back({false, true});
        PlotLabelGroups.push_back({"Smeared true", "Unfolded"});
        PlotTitles.push_back("Unfolded/" + UnfHistTitles[i]);
        XLabels.push_back("Kinetic energy [MeV]");
        YLabels.push_back("Counts");
        PlotStacked.push_back(false);
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

    std::vector<TH2*> TwoDPlots = {
        hResponseMatrix,
        hCovariance, 
        hSmearing
    };

    std::vector<TString> TwoDTitles = {
        "Response/ResponseMatrix",
        "Covariance/UnfoldedCovariance",
        "Smearing/AddSmearing"
    };

    std::vector<std::pair<double,double>> TwoDRanges = {
        {0, 0},
        {0, 0},
        {0, 0}
    };

    std::vector<bool> TwoDDisplayNumbers = {
        false,
        false,
        false
    };

    printTwoDPlots(SaveDir, TwoDPlots, TwoDTitles, TwoDRanges, TwoDDisplayNumbers);
}