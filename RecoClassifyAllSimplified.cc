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

    /////////////////////////////////
    // Files for event information //
    /////////////////////////////////

    std::ofstream outFileChExchBkg("files/ClassifyAll/ChExchBackground.txt");

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

    // Interacting pion abs 0p energy
    TH1D* hPionAbs0pKE         = new TH1D("hPionAbs0pKE", "hPionAbs0pKE;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbs0pKETrue     = new TH1D("hPionAbs0PKETrue", "hPionAbs0PKETrue;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbs0pKEAbsNp    = new TH1D("hPionAbs0pKEAbsNp", "hPionAbs0pKEAbsNp;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbs0pKEScatter  = new TH1D("hPionAbs0pKEScatter", "hPionAbs0pKEScatter;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbs0pKEChExch   = new TH1D("hPionAbs0pKEChExch", "hPionAbs0pKEChExch;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbs0pKEMuon     = new TH1D("hPionAbs0pKEMuon", "hPionAbs0pKEMuon;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbs0pKEElectron = new TH1D("hPionAbs0pKEElectron", "hPionAbs0pKEElectron;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbs0pKEOther    = new TH1D("hPionAbs0pKEOther", "hPionAbs0pKEOther;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    // Interacting pion abs Np energy
    TH1D* hPionAbsNpKE         = new TH1D("hPionAbsNpKE", "hPionAbsNpKE;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbsNpKETrue     = new TH1D("hPionAbsNpKETrue", "hPionAbsNpKETrue;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbsNpKEAbs0p    = new TH1D("hPionAbsNpKEAbs0p", "hPionAbsNpKEAbs0p;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbsNpKEScatter  = new TH1D("hPionAbsNpKEScatter", "hPionAbsNpKEScatter;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbsNpKEChExch   = new TH1D("hPionAbsNpKEChExch", "hPionAbsNpKEChExch;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbsNpKEMuon     = new TH1D("hPionAbsNpKEMuon", "hPionAbsNpKEMuon;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbsNpKEElectron = new TH1D("hPionAbsNpKEElectron", "hPionAbsNpKEElectron;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionAbsNpKEOther    = new TH1D("hPionAbsNpKEOther", "hPionAbsNpKEOther;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    // Interacting pion scattering energy
    TH1D* hPionScatterKE         = new TH1D("hPionScatterKE", "hPionScatterKE;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionScatterKETrue     = new TH1D("hPionScatterKETrue", "hPionScatterKETrue;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionScatterKEAbs0p    = new TH1D("hPionScatterKEAbs0p", "hPionScatterKEAbs0p;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionScatterKEAbsNp    = new TH1D("hPionScatterKEAbsNp", "hPionScatterKEAbsNp;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionScatterKEChExch   = new TH1D("hPionScatterKEChExch", "hPionScatterKEChExch;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionScatterKEMuon     = new TH1D("hPionScatterKEMuon", "hPionScatterKEMuon;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionScatterKEElectron = new TH1D("hPionScatterKEElectron", "hPionScatterKEElectron;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionScatterKEOther    = new TH1D("hPionScatterKEOther", "hPionScatterKEOther;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    // Interacting pion charge exchange energy
    TH1D* hPionChExchKE         = new TH1D("hPionChExchKE", "hPionChExchKE;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionChExchKETrue     = new TH1D("hPionChExchKETrue", "hPionChExchKETrue;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionChExchKEAbs0p    = new TH1D("hPionChExchKEAbs0p", "hPionChExchKEAbs0p;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionChExchKEAbsNp    = new TH1D("hPionChExchKEAbsNp", "hPionChExchKEAbsNp;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionChExchKEScatter  = new TH1D("hPionChExchKEScatter", "hPionChExchKEScatter;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionChExchKEMuon     = new TH1D("hPionChExchKEMuon", "hPionChExchKEMuon;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionChExchKEElectron = new TH1D("hPionChExchKEElectron", "hPionChExchKEElectron;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hPionChExchKEOther    = new TH1D("hPionChExchKEOther", "hPionChExchKEOther;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

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

    std::vector<TH1*> TrueAbs0pAs = {
        hTrueAbs0pKEAsAbs0p, hTrueAbs0pKEAsAbsNp, hTrueAbs0pKEAsScatter, hTrueAbs0pKEAsChExch, hTrueAbs0pKERejected
    };

    // True abs Np
    TH1D* hTrueAbsNpKE          = new TH1D("hTrueAbsNpKE", "hTrueAbsNpKE;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbsNpKEAsAbs0p   = new TH1D("hTrueAbsNpKEAsAbs0p", "hTrueAbsNpKEAsAbs0p;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbsNpKEAsAbsNp   = new TH1D("hTrueAbsNpKEAsAbsNp", "hTrueAbsNpKEAsAbsNp;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbsNpKEAsScatter = new TH1D("hTrueAbsNpKEAsScatter", "hTrueAbsNpKEAsScatter;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbsNpKEAsChExch  = new TH1D("hTrueAbsNpKEAsChExch", "hTrueAbsNpKEAsChExch;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueAbsNpKERejected  = new TH1D("hTrueAbsNpKERejected", "hTrueAbsNpKERejected;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    std::vector<TH1*> TrueAbsNpAs = {
        hTrueAbsNpKEAsAbs0p, hTrueAbsNpKEAsAbsNp, hTrueAbsNpKEAsScatter, hTrueAbsNpKEAsChExch, hTrueAbsNpKERejected
    };

    // True scatter 0p
    TH1D* hTrueScatterKE          = new TH1D("hTrueScatterKE", "hTrueScatterKE;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueScatterKEAsAbs0p   = new TH1D("hTrueScatterKEAsAbs0p", "hTrueScatterKEAsAbs0p;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueScatterKEAsAbsNp   = new TH1D("hTrueScatterKEAsAbsNp", "hTrueScatterKEAsAbsNp;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueScatterKEAsScatter = new TH1D("hTrueScatterKEAsScatter", "hTrueScatterKEAsScatter;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueScatterKEAsChExch  = new TH1D("hTrueScatterKEAsChExch", "hTrueScatterKEAsChExch;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueScatterKERejected  = new TH1D("hTrueScatterKERejected", "hTrueScatterKERejected;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    std::vector<TH1*> TrueScatterAs = {
        hTrueScatterKEAsAbs0p, hTrueScatterKEAsAbsNp, hTrueScatterKEAsScatter, hTrueScatterKEAsChExch, hTrueScatterKERejected
    };

    // True charge exchange
    TH1D* hTrueChExchKE          = new TH1D("hTrueChExchKE", "hTrueChExchKE;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueChExchKEAsAbs0p   = new TH1D("hTrueChExchKEAsAbs0p", "hTrueChExchKEAsAbs0p;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueChExchKEAsAbsNp   = new TH1D("hTrueChExchKEAsAbsNp", "hTrueChExchKEAsAbsNp;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueChExchKEAsScatter = new TH1D("hTrueChExchKEAsScatter", "hTrueChExchKEAsScatter;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueChExchKEAsChExch  = new TH1D("hTrueChExchKEAsChExch", "hTrueChExchKEAsChExch;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
    TH1D* hTrueChExchKERejected  = new TH1D("hTrueChExchKERejected", "hTrueChExchKERejected;;", NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);

    std::vector<TH1*> TrueChExchAs = {
        hTrueChExchKEAsAbs0p, hTrueChExchKEAsAbsNp, hTrueChExchKEAsScatter, hTrueChExchKEAsChExch, hTrueChExchKERejected
    };

    //////////////
    // Matrices //
    //////////////

    NUM_SIGNAL_TYPES = 4; // override, for now

    std::vector<TH2*> ProbabilityMatrices;
    for (int i = 1; i <= NUM_BINS_KE; ++i) {
        double low_bin  = LOWER_BOUND_KE + (i - 1) * (UPPER_BOUND_KE - LOWER_BOUND_KE) / NUM_BINS_KE;
        double high_bin = LOWER_BOUND_KE + i * (UPPER_BOUND_KE - LOWER_BOUND_KE) / NUM_BINS_KE;
        TH2D* hPMatrix = new TH2D(
            Form("hPMatrix_Bin_%d_%d", (int)low_bin, (int)high_bin), Form("hPMatrix_Bin_%d_%d;Reco interaction;True interaction", (int)low_bin, (int)high_bin), 
            NUM_SIGNAL_TYPES, 0, NUM_SIGNAL_TYPES, 
            NUM_SIGNAL_TYPES, 0, NUM_SIGNAL_TYPES
        );
        ProbabilityMatrices.push_back(hPMatrix);
    }

    std::vector<TH1*> TotalEventsHistos = {
        hTrueAbs0pKE, hTrueAbsNpKE, hTrueScatterKE, hTrueChExchKE
    };

    std::vector<std::vector<TH1*>> TrueRecoAs = {
        TrueAbs0pAs, TrueAbsNpAs, TrueScatterAs, TrueChExchAs
    };

    //////////////////////
    // Loop over events //
    //////////////////////

    Int_t NumEntries = (Int_t) tree->GetEntries();
    std::cout << "Num entries: " << NumEntries << std::endl;

    for (Int_t i = 0; i < NumEntries; ++i) {
        tree->GetEntry(i);

        // In this script, we replace background types 6 and 12 with background types 13 and 14,
        // because we are interested in contamination across the bins that we are measuring
        if (backgroundType == 12 || (backgroundType == 6 && numVisibleProtons == 0)) {
            backgroundType = 13;
        } else if (backgroundType == 6 && numVisibleProtons > 0) {
            backgroundType = 14;
        }

        // If no track matched to wire-chamber, skip
        if (WC2TPCtrkID == -99999) continue;
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

        // Fill in true interaction energy histograms
        if (backgroundType == 0) {
            hTrueAbs0pKE->Fill(energyAtVertex);
        } else if (backgroundType == 1) {
            hTrueAbsNpKE->Fill(energyAtVertex);
        } else if (backgroundType == 13 || backgroundType == 14) {
            hTrueScatterKE->Fill(energyAtVertex);
        } else if (backgroundType == 7) {
            hTrueChExchKE->Fill(energyAtVertex);
        }

        // If did not obtain neural network electron shower probabilities, skip 
        if (!obtainedProbabilities) {
            if (backgroundType == 0) {
                hTrueAbs0pKERejected->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                hTrueAbsNpKERejected->Fill(energyAtVertex);
            } else if (backgroundType == 13 || backgroundType == 14) {
                hTrueScatterKERejected->Fill(energyAtVertex);
            } else if (backgroundType == 7) {
                hTrueChExchKERejected->Fill(energyAtVertex);
            }
            continue;
        }
        if (showerProb >= SHOWER_PROB_CUT){
            if (backgroundType == 0) {
                hTrueAbs0pKERejected->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                hTrueAbsNpKERejected->Fill(energyAtVertex);
            } else if (backgroundType == 13 || backgroundType == 14) {
                hTrueScatterKERejected->Fill(energyAtVertex);
            } else if (backgroundType == 7) {
                hTrueChExchKERejected->Fill(energyAtVertex);
            }
            continue;
        }
        hNotAnElectron->Fill(backgroundType);

        if (!isWithinReducedVolume(breakPointX, breakPointY, breakPointZ)) {
            if (backgroundType == 0) {
                hTrueAbs0pKERejected->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                hTrueAbsNpKERejected->Fill(energyAtVertex);
            } else if (backgroundType == 13 || backgroundType == 14) {
                hTrueScatterKERejected->Fill(energyAtVertex);
            } else if (backgroundType == 7) {
                hTrueChExchKERejected->Fill(energyAtVertex);
            }
            continue;
        }
        hPrimaryInRedVol->Fill(backgroundType);

        if (minChi2 == pionChi2) {
            if (backgroundType == 0) {
                hTrueAbs0pKERejected->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                hTrueAbsNpKERejected->Fill(energyAtVertex);
            } else if (backgroundType == 13 || backgroundType == 14) {
                hTrueScatterKERejected->Fill(energyAtVertex);
            } else if (backgroundType == 7) {
                hTrueChExchKERejected->Fill(energyAtVertex);
            }
            continue;
        } else if (minChi2 == protonChi2) {
            if (backgroundType == 0) {
                hTrueAbs0pKERejected->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                hTrueAbsNpKERejected->Fill(energyAtVertex);
            } else if (backgroundType == 13 || backgroundType == 14) {
                hTrueScatterKERejected->Fill(energyAtVertex);
            } else if (backgroundType == 7) {
                hTrueChExchKERejected->Fill(energyAtVertex);
            }
            continue; // reject events matching to proton primary
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
                    hTrueAbs0pKERejected->Fill(energyAtVertex);
                } else if (backgroundType == 1) {
                    hTrueAbsNpKERejected->Fill(energyAtVertex);
                } else if (backgroundType == 13 || backgroundType == 14) {
                    hTrueScatterKERejected->Fill(energyAtVertex);
                } else if (backgroundType == 7) {
                    hTrueChExchKERejected->Fill(energyAtVertex);
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
                hTrueAbs0pKEAsScatter->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                hTrueAbsNpKEAsScatter->Fill(energyAtVertex);
            } else if (backgroundType == 13 || backgroundType == 14) {
                hTrueScatterKEAsScatter->Fill(energyAtVertex);
            } else if (backgroundType == 7) {
                hTrueChExchKEAsScatter->Fill(energyAtVertex);
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
                hTrueAbs0pKEAsAbsNp->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                hTrueAbsNpKEAsAbsNp->Fill(energyAtVertex);
            } else if (backgroundType == 13 || backgroundType == 14) {
                hTrueScatterKEAsAbsNp->Fill(energyAtVertex);
            } else if (backgroundType == 7) {
                hTrueChExchKEAsAbsNp->Fill(energyAtVertex);
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
                hTrueAbs0pKEAsChExch->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                hTrueAbsNpKEAsChExch->Fill(energyAtVertex);
            } else if (backgroundType == 13 || backgroundType == 14) {
                hTrueScatterKEAsChExch->Fill(energyAtVertex);
            } else if (backgroundType == 7) {
                hTrueChExchKEAsChExch->Fill(energyAtVertex);
            }

            outFileChExchBkg << "Run: " << run << " subrun: " << subrun << " event: " << event << std::endl;
            outFileChExchBkg << " Classified as: " << backgroundTypes[backgroundType] << " with energy at vertex: " << energyAtVertex << std::endl;
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
                hTrueAbs0pKEAsAbs0p->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                hTrueAbsNpKEAsAbs0p->Fill(energyAtVertex);
            } else if (backgroundType == 13 || backgroundType == 14) {
                hTrueScatterKEAsAbs0p->Fill(energyAtVertex);
            } else if (backgroundType == 7) {
                hTrueChExchKEAsAbs0p->Fill(energyAtVertex);
            }

            continue;
        }
        hNotPionAbs0p->Fill(backgroundType);

        // Anything left here is rejected
        if (backgroundType == 0) {
            hTrueAbs0pKERejected->Fill(energyAtVertex);
        } else if (backgroundType == 1) {
            hTrueAbsNpKERejected->Fill(energyAtVertex);
        } else if (backgroundType == 13 || backgroundType == 14) {
            hTrueScatterKERejected->Fill(energyAtVertex);
        } else if (backgroundType == 7) {
            hTrueChExchKERejected->Fill(energyAtVertex);
        }
    }

    ///////////////////////////////////////////
    // Information about events at each step //
    ///////////////////////////////////////////

    std::cout << std::endl;
    std::cout << "Original sample composition: " << std::endl;
    printBackgroundInfo(hTotalEvents, std::cout);

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
    // Compute probability matrix //
    ////////////////////////////////

    for (int iBin = 1; iBin <= NUM_BINS_KE; ++iBin) {
        TH2* currentMatrix = ProbabilityMatrices.at(iBin - 1);
        for (int column = 0; column < NUM_SIGNAL_TYPES; ++column) {
            double denom = TotalEventsHistos.at(column)->GetBinContent(iBin);
            for (int row = 0; row < NUM_SIGNAL_TYPES; ++row) {
                double num = TrueRecoAs.at(column).at(row)->GetBinContent(iBin);

                if (denom == 0) currentMatrix->SetBinContent(row + 1, column + 1, 0);
                else currentMatrix->SetBinContent(row + 1, column + 1, num / denom);
            }
        }
    }

    // Convert TH2 into TMatrix
    std::vector<TMatrixD> PInvMatrices;
    std::vector<TH2*>     PInvMatricesHistos;
    for (int iBin = 1; iBin <= NUM_BINS_KE; ++iBin) {
        TMatrixD mat(NUM_SIGNAL_TYPES, NUM_SIGNAL_TYPES);
        for (int i = 1; i <= NUM_SIGNAL_TYPES; ++i) {
            for (int j = 1; j <= NUM_SIGNAL_TYPES; ++j) {
            mat(i-1, j-1) = ProbabilityMatrices.at(iBin - 1)->GetBinContent(i, j);
            }
        }

        TMatrixD inv = mat.Invert();
        PInvMatrices.push_back(inv);

        // Save the inverse matrix as a histogram for plots
        double low_bin  = LOWER_BOUND_KE + (iBin - 1) * (UPPER_BOUND_KE - LOWER_BOUND_KE) / NUM_BINS_KE;
        double high_bin = LOWER_BOUND_KE + iBin * (UPPER_BOUND_KE - LOWER_BOUND_KE) / NUM_BINS_KE;
        TH2D* hPInv = new TH2D(
            Form("hPInvMatrix_Bin_%d_%d", (int)low_bin, (int)high_bin),
            Form("hPInvMatrix_Bin_%d_%d;Reco interaction;True interaction", (int)low_bin, (int)high_bin),
            NUM_SIGNAL_TYPES, 0, NUM_SIGNAL_TYPES,
            NUM_SIGNAL_TYPES, 0, NUM_SIGNAL_TYPES
        );
        for (int i = 0; i < NUM_SIGNAL_TYPES; ++i) {
            for (int j = 0; j < NUM_SIGNAL_TYPES; ++j) {
                hPInv->SetBinContent(i + 1, j + 1, inv(i, j));
            }
        }
        PInvMatricesHistos.push_back(hPInv);
    }

    // Get vectors from reconstructed signals
    std::vector<TVectorD> SignalUnfVectors;
    for (int iBin = 1; iBin <= NUM_BINS_KE; ++iBin) {
        TVectorD vec(NUM_SIGNAL_TYPES);
        for (int i = 0; i < NUM_SIGNAL_TYPES; ++i) {
            vec(i) = RecoSignals[i]->GetBinContent(iBin);
        }
        
        TVectorD unfoldedVec = PInvMatrices.at(iBin - 1) * vec;
        SignalUnfVectors.push_back(unfoldedVec);
    }

    // Organize unfolded results
    std::vector<TH1*> UnfHistograms;
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        TH1* hist = new TH1D(Form("hUnfSignal_%d", iSignal), Form("Unfolded Signal %d", iSignal), NUM_BINS_KE, LOWER_BOUND_KE, UPPER_BOUND_KE);
        for (int iBin = 1; iBin <= NUM_BINS_KE; ++iBin) {
            hist->SetBinContent(iBin, SignalUnfVectors.at(iBin - 1)(iSignal));
        }
        UnfHistograms.push_back(hist);
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

        // Interacting KE
        {hPionAbs0pKETrue, hPionAbs0pKEAbsNp, hPionAbs0pKEScatter, hPionAbs0pKEChExch, hPionAbs0pKEMuon, hPionAbs0pKEElectron, hPionAbs0pKEOther},
        {hPionAbsNpKETrue, hPionAbsNpKEAbs0p, hPionAbsNpKEScatter, hPionAbsNpKEChExch, hPionAbsNpKEMuon, hPionAbsNpKEElectron, hPionAbsNpKEOther},
        {hPionScatterKETrue, hPionScatterKEAbs0p, hPionScatterKEAbsNp, hPionScatterKEChExch, hPionScatterKEMuon, hPionScatterKEElectron, hPionScatterKEOther},
        {hPionChExchKETrue, hPionChExchKEAbs0p, hPionChExchKEAbsNp, hPionChExchKEScatter, hPionChExchKEMuon, hPionChExchKEElectron, hPionChExchKEOther},

        // True events classified breakdown
        {hTrueAbs0pKEAsAbs0p, hTrueAbs0pKEAsAbsNp, hTrueAbs0pKEAsScatter, hTrueAbs0pKEAsChExch, hTrueAbs0pKERejected},
        {hTrueAbsNpKEAsAbs0p, hTrueAbsNpKEAsAbsNp, hTrueAbsNpKEAsScatter, hTrueAbsNpKEAsChExch, hTrueAbsNpKERejected},
        {hTrueScatterKEAsAbs0p, hTrueScatterKEAsAbsNp, hTrueScatterKEAsScatter, hTrueScatterKEAsChExch, hTrueScatterKERejected},
        {hTrueChExchKEAsAbs0p, hTrueChExchKEAsAbsNp, hTrueChExchKEAsScatter, hTrueChExchKEAsChExch, hTrueChExchKERejected}
    };

    std::vector<std::vector<TString>> PlotLabelGroups = {
        // Incident KE
        {"Pions", "Muons", "Electrons"},

        // Interacting KE
        {"True", "Abs Np", "Scatter", "Ch. exch.", "Muon", "Electron", "Other"},
        {"True", "Abs Np", "Scatter", "Ch. exch.", "Muon", "Electron", "Other"},
        {"True", "Abs 0p", "Abs Np", "Ch. exch.", "Muon", "Electron", "Other"},
        {"True", "Abs 0p", "Abs Np", "Scatter", "Muon", "Electron", "Other"},

        // True events classified breakdown
        {"Abs 0p", "Abs Np", "Scatter", "Ch. exch.", "Rejected"},
        {"Abs 0p", "Abs Np", "Scatter", "Ch. exch.", "Rejected"},
        {"Abs 0p", "Abs Np", "Scatter", "Ch. exch.", "Rejected"},
        {"Abs 0p", "Abs Np", "Scatter", "Ch. exch.", "Rejected"}
    };

    std::vector<TString> PlotTitles = {
        // Incident KE
        "Incident/IncidentKE",

        // Interacting KE
        "RecoInteracting/Abs0pInteractingKE",
        "RecoInteracting/AbsNpInteractingKE",
        "RecoInteracting/ScatterInteractingKE",
        "RecoInteracting/ChExchInteractingKE",

        // True events classified breakdown
        "TrueInteracting/TrueAbs0pBreakdown",
        "TrueInteracting/TrueAbsNpBreakdown",
        "TrueInteracting/TrueScatterBreakdown",
        "TrueInteracting/TrueChExchBreakdown"
    };

    std::vector<TString> XLabels = {
        // Incident KE
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
        "Kinetic energy [MeV]"
    };

    std::vector<TString> YLabels = {
        // Incident KE
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
        "Counts"
    };

    std::vector<bool> PlotStacked = {
        // Incident KE
        true,

        // Interacting KE
        true,
        true,
        true,
        true,

        // True events classified breakdown
        true,
        true,
        true,
        true
    };

        // Add each unfolded histogram as a single plot group for plotting
        std::vector<TString> UnfHistTitles = {
            "UnfoldedAbs0p",
            "UnfoldedAbsNp",
            "UnfoldedScatter",
            "UnfoldedChExch"
        };

        for (size_t i = 0; i < UnfHistograms.size(); ++i) {
            PlotGroups.push_back({UnfHistograms[i], TotalEventsHistos[i]});
            PlotLabelGroups.push_back({"Unfolded", "True signal"});
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
        PlotStacked
    );

    ///////////////////////////
    // Two-dimensional plots //
    ///////////////////////////

    std::vector<TH2*> TwoDPlots;
    std::vector<TString> TwoDTitles;
    std::vector<std::pair<double,double>> TwoDRanges;
    std::vector<bool> TwoDDisplayNumbers;

    for (int i = 0; i < NUM_BINS_KE; ++i) {
        TwoDPlots.push_back(ProbabilityMatrices.at(i));
        TwoDTitles.push_back("ProbMatrices/ProbabilityMatrix_Bin" + std::to_string(i+1));
        TwoDRanges.push_back({0, 1});
        TwoDDisplayNumbers.push_back(true);

        TwoDPlots.push_back(PInvMatricesHistos.at(i));
        TwoDTitles.push_back("ProbInvMatrices/PInvMatrix_Bin" + std::to_string(i+1));
        TwoDRanges.push_back({0, 1});
        TwoDDisplayNumbers.push_back(true);
    }

    printTwoDPlots(SaveDir, TwoDPlots, TwoDTitles, TwoDRanges, TwoDDisplayNumbers);
}