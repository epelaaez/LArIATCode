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
    TString RootFilePath = "/exp/lariat/app/users/epelaez/files/RecoNNAllEval_histo2.root";
    std::unique_ptr<TFile> File(TFile::Open(RootFilePath));
    TDirectory* Directory = (TDirectory*)File->Get("RecoNNAllEval");

    ///////////////////
    // Load branches //
    ///////////////////

    // Load tree and branches
    TTree* tree = (TTree*) Directory->Get<TTree>("RecoNNAllEvalTree");

    int run, subrun, event; bool isData;
    tree->SetBranchAddress("run", &run); 
    tree->SetBranchAddress("subrun", &subrun); 
    tree->SetBranchAddress("event", &event);
    tree->SetBranchAddress("isData", &isData);

    // Signal information
    int backgroundType, numVisibleProtons; 
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
    int truthPrimaryPDG, truthPrimaryID;
    double truthPrimaryVertexX, truthPrimaryVertexY, truthPrimaryVertexZ;
    double truthPrimaryIncidentKE, truthPrimaryVertexKE;
    std::vector<int>*         truthPrimaryDaughtersPDG     = nullptr;
    std::vector<std::string>* truthPrimaryDaughtersProcess = nullptr;
    std::vector<double>*      truthPrimaryDaughtersKE      = nullptr;
    tree->SetBranchAddress("truthPrimaryPDG", &truthPrimaryPDG);
    tree->SetBranchAddress("truthPrimaryID", &truthPrimaryID);
    tree->SetBranchAddress("truthPrimaryIncidentKE", &truthPrimaryIncidentKE);
    tree->SetBranchAddress("truthPrimaryVertexKE", &truthPrimaryVertexKE);
    tree->SetBranchAddress("truthPrimaryVertexX", &truthPrimaryVertexX);
    tree->SetBranchAddress("truthPrimaryVertexY", &truthPrimaryVertexY);
    tree->SetBranchAddress("truthPrimaryVertexZ", &truthPrimaryVertexZ);
    tree->SetBranchAddress("truthPrimaryDaughtersPDG", &truthPrimaryDaughtersPDG);
    tree->SetBranchAddress("truthPrimaryDaughtersProcess", &truthPrimaryDaughtersProcess);
    tree->SetBranchAddress("truthPrimaryDaughtersKE", &truthPrimaryDaughtersKE);

    // True particle location information
    std::vector<double>* truthPrimaryLocationX = nullptr;
    std::vector<double>* truthPrimaryLocationY = nullptr;
    std::vector<double>* truthPrimaryLocationZ = nullptr;
    tree->SetBranchAddress("truthPrimaryLocationX", &truthPrimaryLocationX);
    tree->SetBranchAddress("truthPrimaryLocationY", &truthPrimaryLocationY);
    tree->SetBranchAddress("truthPrimaryLocationZ", &truthPrimaryLocationZ);

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

    // Truth level charge exchange information
    std::vector<int>* chExchShowerPDGs = nullptr;
    std::vector<int>* chExchShowerIDs = nullptr;
    std::vector<std::vector<double>>* chExchShowerStart = nullptr;
    std::vector<std::vector<double>>* chExchShowerEnd = nullptr;
    std::vector<int>* chExchShowerNeutralPionDaughtersID = nullptr;
    tree->SetBranchAddress("chExchShowerPDGs", &chExchShowerPDGs);
    tree->SetBranchAddress("chExchShowerIDs", &chExchShowerIDs);
    tree->SetBranchAddress("chExchShowerStart", &chExchShowerStart);
    tree->SetBranchAddress("chExchShowerEnd", &chExchShowerEnd);
    tree->SetBranchAddress("chExchShowerNeutralPionDaughtersID", &chExchShowerNeutralPionDaughtersID);

    // Primaries information
    std::vector<double>* primariesStartX = nullptr;
    std::vector<double>* primariesStartY = nullptr;
    std::vector<double>* primariesStartZ = nullptr;
    std::vector<double>* primariesEndX = nullptr;
    std::vector<double>* primariesEndY = nullptr;
    std::vector<double>* primariesEndZ = nullptr;
    std::vector<int>*    primariesID = nullptr;
    tree->SetBranchAddress("primariesStartX", &primariesStartX);
    tree->SetBranchAddress("primariesStartY", &primariesStartY);
    tree->SetBranchAddress("primariesStartZ", &primariesStartZ);
    tree->SetBranchAddress("primariesEndX", &primariesEndX);
    tree->SetBranchAddress("primariesEndY", &primariesEndY);
    tree->SetBranchAddress("primariesEndZ", &primariesEndZ);
    tree->SetBranchAddress("primariesID", &primariesID);

    // Information about post-scattering interactions
    std::vector<int>*                 secondaryInteractionTypes = nullptr;
    std::vector<int>*                 secondaryInteractionTrkID = nullptr;
    std::vector<double>*              secondaryInteractionInteractingKE = nullptr;
    std::vector<double>*              secondaryInteractionAngle = nullptr;
    std::vector<double>*              secondaryInteractionXPosition = nullptr;
    std::vector<double>*              secondaryInteractionYPosition = nullptr;
    std::vector<double>*              secondaryInteractionZPosition = nullptr;
    std::vector<std::vector<double>>* secondaryIncidentKEContributions = nullptr;
    tree->SetBranchAddress("secondaryInteractionTypes", &secondaryInteractionTypes);
    tree->SetBranchAddress("secondaryInteractionTrkID", &secondaryInteractionTrkID);
    tree->SetBranchAddress("secondaryInteractionInteractingKE", &secondaryInteractionInteractingKE);
    tree->SetBranchAddress("secondaryInteractionAngle", &secondaryInteractionAngle);
    tree->SetBranchAddress("secondaryInteractionXPosition", &secondaryInteractionXPosition);
    tree->SetBranchAddress("secondaryInteractionYPosition", &secondaryInteractionYPosition);
    tree->SetBranchAddress("secondaryInteractionZPosition", &secondaryInteractionZPosition);
    tree->SetBranchAddress("secondaryIncidentKEContributions", &secondaryIncidentKEContributions);

    /////////////////////////////////
    // Files for event information //
    /////////////////////////////////

    std::ofstream outFileChExchBkg("files/ClassifyAll/ChExchBackground.txt");
    std::ofstream outFileAbs0pBkg("files/ClassifyAll/Abs0pBackground.txt");
    TFile* comparisonsFile = new TFile("/exp/lariat/app/users/epelaez/files/DataMCComparisons.root", "UPDATE");

    ///////////////////////
    // Create histograms //
    ///////////////////////

    TH1D* hTotalEvents = new TH1D("hTotalEvents", "hTotalEvents", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);

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
    TH1D* hIncidentKE         = new TH1D("hIncidentKE", "hIncidentKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hIncidentKEPion     = new TH1D("hIncidentKEPion", "hIncidentKEPion;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hIncidentKEElectron = new TH1D("hIncidentKEElectron", "hIncidentKEElectron;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hIncidentKEMuon     = new TH1D("hIncidentKEMuon", "hIncidentKEMuon;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueIncidentKE     = new TH1D("hTrueIncidentKE", "hTrueIncidentKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Interacting pion abs 0p energy
    TH1D* hPionAbs0pKE         = new TH1D("hPionAbs0pKE", "hPionAbs0pKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbs0pKETrue     = new TH1D("hPionAbs0PKETrue", "hPionAbs0PKETrue;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbs0pKEAbsNp    = new TH1D("hPionAbs0pKEAbsNp", "hPionAbs0pKEAbsNp;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbs0pKEScatter  = new TH1D("hPionAbs0pKEScatter", "hPionAbs0pKEScatter;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbs0pKEChExch   = new TH1D("hPionAbs0pKEChExch", "hPionAbs0pKEChExch;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbs0pKEMuon     = new TH1D("hPionAbs0pKEMuon", "hPionAbs0pKEMuon;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbs0pKEElectron = new TH1D("hPionAbs0pKEElectron", "hPionAbs0pKEElectron;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbs0pKEOther    = new TH1D("hPionAbs0pKEOther", "hPionAbs0pKEOther;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    TH1D* hPionAbs0pKEMuonTG    = new TH1D("hPionAbs0pKEMuonTG", "hPionAbs0pKEMuonTG;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbs0pKEMuonDecay = new TH1D("hPionAbs0pKEMuonDecay", "hPionAbs0pKEMuonDecay;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbs0pKEMuonCAR   = new TH1D("hPionAbs0pKEMuonCAR", "hPionAbs0pKEMuonCAR;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Estimated backgrounds for pion abs 0p
    std::vector<TH1*> PionAbs0pBkg = {
        hPionAbs0pKEMuon, hPionAbs0pKEElectron, hPionAbs0pKEOther
    };

    // Interacting pion abs Np energy
    TH1D* hPionAbsNpKE         = new TH1D("hPionAbsNpKE", "hPionAbsNpKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbsNpKETrue     = new TH1D("hPionAbsNpKETrue", "hPionAbsNpKETrue;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbsNpKEAbs0p    = new TH1D("hPionAbsNpKEAbs0p", "hPionAbsNpKEAbs0p;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbsNpKEScatter  = new TH1D("hPionAbsNpKEScatter", "hPionAbsNpKEScatter;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbsNpKEChExch   = new TH1D("hPionAbsNpKEChExch", "hPionAbsNpKEChExch;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbsNpKEMuon     = new TH1D("hPionAbsNpKEMuon", "hPionAbsNpKEMuon;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbsNpKEElectron = new TH1D("hPionAbsNpKEElectron", "hPionAbsNpKEElectron;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbsNpKEOther    = new TH1D("hPionAbsNpKEOther", "hPionAbsNpKEOther;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Estimated backgrounds for pion abs Np
    std::vector<TH1*> PionAbsNpBkg = {
        hPionAbsNpKEMuon, hPionAbsNpKEElectron, hPionAbsNpKEOther
    };

    // Interacting pion scattering energy
    TH1D* hPionScatterKE         = new TH1D("hPionScatterKE", "hPionScatterKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionScatterKETrue     = new TH1D("hPionScatterKETrue", "hPionScatterKETrue;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionScatterKEAbs0p    = new TH1D("hPionScatterKEAbs0p", "hPionScatterKEAbs0p;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionScatterKEAbsNp    = new TH1D("hPionScatterKEAbsNp", "hPionScatterKEAbsNp;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionScatterKEChExch   = new TH1D("hPionScatterKEChExch", "hPionScatterKEChExch;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionScatterKEMuon     = new TH1D("hPionScatterKEMuon", "hPionScatterKEMuon;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionScatterKEElectron = new TH1D("hPionScatterKEElectron", "hPionScatterKEElectron;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionScatterKEOther    = new TH1D("hPionScatterKEOther", "hPionScatterKEOther;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    TH1D* hPionScatterKEMuonTG    = new TH1D("hPionScatterKEMuonTG", "hPionScatterKEMuonTG;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionScatterKEMuonDecay = new TH1D("hPionScatterKEMuonDecay", "hPionScatterKEMuonDecay;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionScatterKEMuonCAR   = new TH1D("hPionScatterKEMuonCAR", "hPionScatterKEMuonCAR;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Estimated backgrounds for pion scattering
    std::vector<TH1*> PionScatterBkg = {
        hPionScatterKEMuon, hPionScatterKEElectron, hPionScatterKEOther
    };

    // Interacting pion charge exchange energy
    TH1D* hPionChExchKE         = new TH1D("hPionChExchKE", "hPionChExchKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionChExchKETrue     = new TH1D("hPionChExchKETrue", "hPionChExchKETrue;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionChExchKEAbs0p    = new TH1D("hPionChExchKEAbs0p", "hPionChExchKEAbs0p;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionChExchKEAbsNp    = new TH1D("hPionChExchKEAbsNp", "hPionChExchKEAbsNp;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionChExchKEScatter  = new TH1D("hPionChExchKEScatter", "hPionChExchKEScatter;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionChExchKEMuon     = new TH1D("hPionChExchKEMuon", "hPionChExchKEMuon;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionChExchKEElectron = new TH1D("hPionChExchKEElectron", "hPionChExchKEElectron;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionChExchKEOther    = new TH1D("hPionChExchKEOther", "hPionChExchKEOther;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

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
    TH1D* hTrueAbs0pKE            = new TH1D("hTrueAbs0pKE", "hTrueAbs0pKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbs0pKEAsAbs0p     = new TH1D("hTrueAbs0pKEAsAbs0p", "hTrueAbs0pKEAsAbs0p;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbs0pKEAsAbsNp     = new TH1D("hTrueAbs0pKEAsAbsNp", "hTrueAbs0pKEAsAbsNp;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbs0pKEAsScatter   = new TH1D("hTrueAbs0pKEAsScatter", "hTrueAbs0pKEAsScatter;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbs0pKEAsChExch    = new TH1D("hTrueAbs0pKEAsChExch", "hTrueAbs0pKEAsChExch;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbs0pKERejected    = new TH1D("hTrueAbs0pKERejected", "hTrueAbs0pKERejected;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // For each TRUE energy bin, what are abs 0p events reconstructed as?
    std::vector<std::vector<TH1*>> TrueAbs0pAsByBin;
    for (int iEnergyBin = 0; iEnergyBin < NUM_BINS_KE; ++iEnergyBin) {
        std::vector<TH1*> TempVec;
        for (int iInt = 0; iInt < NUM_SIGNAL_TYPES; ++iInt) {
            TH1* hTempHist = new TH1D(Form("hTrueAbs0p_%d_Bin_As_%d", iEnergyBin, iInt), Form("True Abs 0p KE As %d in bin %d", iInt, iEnergyBin), NUM_BINS_KE, ARRAY_KE_BINS.data());
            TempVec.push_back(hTempHist);
        }
        TrueAbs0pAsByBin.push_back(TempVec);
    }

    TH1D* hTrueAbs0pKERejDataProds = new TH1D("hTrueAbs0pKERejDataProds", "hTrueAbs0pKERejDataProds;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbs0pKERejElectron  = new TH1D("hTrueAbs0pKERejElectron", "hTrueAbs0pKERejElectron;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbs0pKERejRedVol    = new TH1D("hTrueAbs0pKERejRedVol", "hTrueAbs0pKERejRedVol;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbs0pKERejPID       = new TH1D("hTrueAbs0pKERejPID", "hTrueAbs0pKERejPID;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbs0pKERejManyPions = new TH1D("hTrueAbs0pKERejManyPions", "hTrueAbs0pKERejManyPions;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbs0pKERejClusters  = new TH1D("hTrueAbs0pKERejClusters", "hTrueAbs0pKERejClusters;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // True abs Np
    TH1D* hTrueAbsNpKE          = new TH1D("hTrueAbsNpKE", "hTrueAbsNpKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbsNpKEAsAbs0p   = new TH1D("hTrueAbsNpKEAsAbs0p", "hTrueAbsNpKEAsAbs0p;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbsNpKEAsAbsNp   = new TH1D("hTrueAbsNpKEAsAbsNp", "hTrueAbsNpKEAsAbsNp;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbsNpKEAsScatter = new TH1D("hTrueAbsNpKEAsScatter", "hTrueAbsNpKEAsScatter;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbsNpKEAsChExch  = new TH1D("hTrueAbsNpKEAsChExch", "hTrueAbsNpKEAsChExch;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbsNpKERejected  = new TH1D("hTrueAbsNpKERejected", "hTrueAbsNpKERejected;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // For each TRUE energy bin, what are abs Np events reconstructed as?
    std::vector<std::vector<TH1*>> TrueAbsNpAsByBin;
    for (int iEnergyBin = 0; iEnergyBin < NUM_BINS_KE; ++iEnergyBin) {
        std::vector<TH1*> TempVec;
        for (int iInt = 0; iInt < NUM_SIGNAL_TYPES; ++iInt) {
            TH1* hTempHist = new TH1D(Form("hTrueAbsNp_%d_Bin_As_%d", iEnergyBin, iInt), Form("True Abs Np KE As %d in bin %d", iInt, iEnergyBin), NUM_BINS_KE, ARRAY_KE_BINS.data());
            TempVec.push_back(hTempHist);
        }
        TrueAbsNpAsByBin.push_back(TempVec);
    }

    TH1D* hTrueAbsNpKERejDataProds = new TH1D("hTrueAbsNpKERejDataProds", "hTrueAbsNpKERejDataProds;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbsNpKERejElectron  = new TH1D("hTrueAbsNpKERejElectron", "hTrueAbsNpKERejElectron;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbsNpKERejRedVol    = new TH1D("hTrueAbsNpKERejRedVol", "hTrueAbsNpKERejRedVol;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbsNpKERejPID       = new TH1D("hTrueAbsNpKERejPID", "hTrueAbsNpKERejPID;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbsNpKERejManyPions = new TH1D("hTrueAbsNpKERejManyPions", "hTrueAbsNpKERejManyPions;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbsNpKERejClusters  = new TH1D("hTrueAbsNpKERejClusters", "hTrueAbsNpKERejClusters;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // True scatter
    TH1D* hTrueScatterKE          = new TH1D("hTrueScatterKE", "hTrueScatterKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueScatterKEAsAbs0p   = new TH1D("hTrueScatterKEAsAbs0p", "hTrueScatterKEAsAbs0p;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueScatterKEAsAbsNp   = new TH1D("hTrueScatterKEAsAbsNp", "hTrueScatterKEAsAbsNp;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueScatterKEAsScatter = new TH1D("hTrueScatterKEAsScatter", "hTrueScatterKEAsScatter;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueScatterKEAsChExch  = new TH1D("hTrueScatterKEAsChExch", "hTrueScatterKEAsChExch;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueScatterKERejected  = new TH1D("hTrueScatterKERejected", "hTrueScatterKERejected;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // For each TRUE energy bin, what are scatter events reconstructed as?
    std::vector<std::vector<TH1*>> TrueScatterAsByBin;
    for (int iEnergyBin = 0; iEnergyBin < NUM_BINS_KE; ++iEnergyBin) {
        std::vector<TH1*> TempVec;
        for (int iInt = 0; iInt < NUM_SIGNAL_TYPES; ++iInt) {
            TH1* hTempHist = new TH1D(Form("hTrueAbsScatter_%d_Bin_As_%d", iEnergyBin, iInt), Form("True Abs Scatter KE As %d in bin %d", iInt, iEnergyBin), NUM_BINS_KE, ARRAY_KE_BINS.data());
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
    TH1D* hTrueChExchKEAsAbs0p   = new TH1D("hTrueChExchKEAsAbs0p", "hTrueChExchKEAsAbs0p;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueChExchKEAsAbsNp   = new TH1D("hTrueChExchKEAsAbsNp", "hTrueChExchKEAsAbsNp;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueChExchKEAsScatter = new TH1D("hTrueChExchKEAsScatter", "hTrueChExchKEAsScatter;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueChExchKEAsChExch  = new TH1D("hTrueChExchKEAsChExch", "hTrueChExchKEAsChExch;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueChExchKERejected  = new TH1D("hTrueChExchKERejected", "hTrueChExchKERejected;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // For each TRUE energy bin, what are charge exchange events reconstructed as?
    std::vector<std::vector<TH1*>> TrueChExchAsByBin;
    for (int iEnergyBin = 0; iEnergyBin < NUM_BINS_KE; ++iEnergyBin) {
        std::vector<TH1*> TempVec;
        for (int iInt = 0; iInt < NUM_SIGNAL_TYPES; ++iInt) {
            TH1* hTempHist = new TH1D(Form("hTrueChExch_%d_Bin_As_%d", iEnergyBin, iInt), Form("True Ch Exch KE As %d in bin %d", iInt, iEnergyBin), NUM_BINS_KE, ARRAY_KE_BINS.data());
            TempVec.push_back(hTempHist);
        }
        TrueChExchAsByBin.push_back(TempVec);
    }

    TH1D* hTrueChExchKERejDataProds = new TH1D("hTrueChExchKERejDataProds", "hTrueChExchKERejDataProds;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueChExchKERejElectron  = new TH1D("hTrueChExchKERejElectron", "hTrueChExchKERejElectron;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueChExchKERejRedVol    = new TH1D("hTrueChExchKERejRedVol", "hTrueChExchKERejRedVol;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueChExchKERejPID       = new TH1D("hTrueChExchKERejPID", "hTrueChExchKERejPID;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueChExchKERejManyPions = new TH1D("hTrueChExchKERejManyPions", "hTrueChExchKERejManyPions;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueChExchKERejClusters  = new TH1D("hTrueChExchKERejClusters", "hTrueChExchKERejClusters;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    ////////////////
    // BDT scores //
    ////////////////

    int PROB_DIVISIONS = 30;

    TH1D* hBestAbs0pAbs0p    = new TH1D("hBestAbs0pAbs0p", "hBestAbs0pAbs0p;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestAbs0pAbsNp    = new TH1D("hBestAbs0pAbsNp", "hBestAbs0pAbsNp;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestAbs0pMuon     = new TH1D("hBestAbs0pMuon", "hBestAbs0pMuon;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestAbs0pElectron = new TH1D("hBestAbs0pElectron", "hBestAbs0pElectron;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestAbs0pScatter  = new TH1D("hBestAbs0pScatter", "hBestAbs0pScatter;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestAbs0pChExch   = new TH1D("hBestAbs0pChExch", "hBestAbs0pChExch;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestAbs0pOther    = new TH1D("hBestAbs0pOther", "hBestAbs0pOther;;", PROB_DIVISIONS, 0, 1);

    TH1D* hBestChExchAbs0p    = new TH1D("hBestChExchAbs0p", "hBestChExchAbs0p;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestChExchAbsNp    = new TH1D("hBestChExchAbsNp", "hBestChExchAbsNp;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestChExchMuon     = new TH1D("hBestChExchMuon", "hBestChExchMuon;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestChExchElectron = new TH1D("hBestChExchElectron", "hBestChExchElectron;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestChExchScatter  = new TH1D("hBestChExchScatter", "hBestChExchScatter;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestChExchChExch   = new TH1D("hBestChExchChExch", "hBestChExchChExch;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestChExchOther    = new TH1D("hBestChExchOther", "hBestChExchOther;;", PROB_DIVISIONS, 0, 1);

    TH1D* hBestScatterAbs0p    = new TH1D("hBestScatterAbs0p", "hBestScatterAbs0p;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestScatterAbsNp    = new TH1D("hBestScatterAbsNp", "hBestScatterAbsNp;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestScatterMuon     = new TH1D("hBestScatterMuon", "hBestScatterMuon;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestScatterElectron = new TH1D("hBestScatterElectron", "hBestScatterElectron;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestScatterScatter  = new TH1D("hBestScatterScatter", "hBestScatterScatter;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestScatterChExch   = new TH1D("hBestScatterChExch", "hBestScatterChExch;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestScatterOther    = new TH1D("hBestScatterOther", "hBestScatterOther;;", PROB_DIVISIONS, 0, 1);

    TH1D* hBestElectronAbs0p    = new TH1D("hBestElectronAbs0p", "hBestElectronAbs0p;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestElectronAbsNp    = new TH1D("hBestElectronAbsNp", "hBestElectronAbsNp;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestElectronMuon     = new TH1D("hBestElectronMuon", "hBestElectronMuon;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestElectronElectron = new TH1D("hBestElectronElectron", "hBestElectronElectron;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestElectronScatter  = new TH1D("hBestElectronScatter", "hBestElectronScatter;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestElectronChExch   = new TH1D("hBestElectronChExch", "hBestElectronChExch;;", PROB_DIVISIONS, 0, 1);
    TH1D* hBestElectronOther    = new TH1D("hBestElectronOther", "hBestElectronOther;;", PROB_DIVISIONS, 0, 1);

    ////////////////////////////////
    // Study unreconstructed hits //
    ////////////////////////////////

    TH1D* hHitClusterSizesAbs0p    = new TH1D("hHitClusterSizesAbs0p", "hHitClusterSizesAbs0p;;", 25, 0, 20);
    TH1D* hHitClusterSizesAbsNp    = new TH1D("hHitClusterSizesAbsNp", "hHitClusterSizesAbsNp;;", 25, 0, 20);
    TH1D* hHitClusterSizesMuon     = new TH1D("hHitClusterSizesMuon", "hHitClusterSizesMuon;;", 25, 0, 20);
    TH1D* hHitClusterSizesElectron = new TH1D("hHitClusterSizesElectron", "hHitClusterSizesElectron;;", 25, 0, 20);
    TH1D* hHitClusterSizesScatter  = new TH1D("hHitClusterSizesScatter", "hHitClusterSizesScatter;;", 25, 0, 20);
    TH1D* hHitClusterSizesChExch   = new TH1D("hHitClusterSizesChExch", "hHitClusterSizesChExch;;", 25, 0, 20);
    TH1D* hHitClusterSizesOther    = new TH1D("hHitClusterSizesOther", "hHitClusterSizesOther;;", 25, 0, 20);

    TH1D* hLargeHitClusterAbs0p    = new TH1D("hLargeHitClusterAbs0p", "hLargeHitClusterAbs0p;;", 5, 0, 5);
    TH1D* hLargeHitClusterAbsNp    = new TH1D("hLargeHitClusterAbsNp", "hLargeHitClusterAbsNp;;", 5, 0, 5);
    TH1D* hLargeHitClusterMuon     = new TH1D("hLargeHitClusterMuon", "hLargeHitClusterMuon;;", 5, 0, 5);
    TH1D* hLargeHitClusterElectron = new TH1D("hLargeHitClusterElectron", "hLargeHitClusterElectron;;", 5, 0, 5);
    TH1D* hLargeHitClusterScatter  = new TH1D("hLargeHitClusterScatter", "hLargeHitClusterScatter;;", 5, 0, 5);
    TH1D* hLargeHitClusterChExch   = new TH1D("hLargeHitClusterChExch", "hLargeHitClusterChExch;;", 5, 0, 5);
    TH1D* hLargeHitClusterOther    = new TH1D("hLargeHitClusterOther", "hLargeHitClusterOther;;", 5, 0, 5);

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

    TH1D* hMCNumWC2TPCMatch      = new TH1D("hMCNumWC2TPCMatch", "hMCNumWC2TPCMatch", 10, 0, 10);
    TH1D* hMCNumTracksInCylinder = new TH1D("hMCNumTracksInCylinder", "hMCNumTracksInCylinder", 10, 0, 10);

    TH1D* hMCTGTrackLengths         = new TH1D("hMCTGTrackLengths", "hMCTGTrackLengths;;", 25, 0, 50);
    TH1D* hMCTGSmallTracks          = new TH1D("hMCTGSmallTracks", "hMCTGSmallTracks;;", 10, 0, 10);
    TH1D* hMCTracksNearVertex       = new TH1D("hMCTracksNearVertex", "hMCTracksNearVertex;;", 10, 0, 10);
    TH1D* hMCTrackLengthsNearVertex = new TH1D("hMCTrackLengthsNearVertex", "hMCTrackLengthsNearVertex;;", 50, 0, 100);
    TH1D* hMCNumTGTracks            = new TH1D("hMCNumTGTracks", "hMCNumTGTracks;;", 10, 0, 10);
    TH1D* hMCShowerProb             = new TH1D("hMCShowerProb", "hMCShowerProb;;", 20, 0, 1.);

    TH1D* hMCBeforeShowerCutSmallTracks = new TH1D("hMCBeforeShowerCutSmallTracks", "hMCBeforeShowerCutSmallTracks;;", 10, 0, 10);
    TH1D* hMCAfterShowerCutSmallTracks  = new TH1D("hMCAfterShowerCutSmallTracks", "hMCAfterShowerCutSmallTracks;;", 10, 0, 10);

    TH2D* hMCSmallVsTGTracks          = new TH2D("hMCSmallVsTGTracks", "MCSmallVsTGTracks;Small Tracks;TG Tracks", 15, 0, 15, 15, 0, 15);
    TH2D* hMCTGNumSmallTracksVsThresh = new TH2D("hMCTGNumSmallTracksVsThresh", "MCTGNumSmallTracksVsThresh;Small Track Length Threshold (cm);Num Small Tracks", 10, 0, 40, 15, 0, 15);

    ///////////////////////////////////
    // Distribution of primary track //
    ///////////////////////////////////

    TH2D* hMCPrimaryTrackPosition = new TH2D(
        "hMCPrimaryTrackPosition", "MCPrimaryTrackPosition;x-position;y-position",
        20, minX, maxX,
        20, minY, maxY
    );

    ////////////////
    // Muon types //
    ////////////////

    TH1D* hTrueMuonTypes = new TH1D("hTrueMuonTypes", "hTrueMuonTypes", 4, 0, 4);

    ////////////////////
    // Load BDT model //
    ////////////////////

    TMVA::Reader BDTReader;
    
    Float_t bdt_WC2TPCTrackLength;
    Float_t bdt_WC2TPCBeginX;
    Float_t bdt_WC2TPCBeginY;
    Float_t bdt_WC2TPCBeginZ;
    Float_t bdt_WC2TPCEndX;
    Float_t bdt_WC2TPCEndY;
    Float_t bdt_WC2TPCEndZ;
    Float_t bdt_WC2TPCdEdx;
    Float_t bdt_numRecoTrksInCylinder;

    BDTReader.AddVariable("WC2TPCTrackLength", &bdt_WC2TPCTrackLength);
    BDTReader.AddVariable("WC2TPCBeginX", &bdt_WC2TPCBeginX);
    BDTReader.AddVariable("WC2TPCBeginY", &bdt_WC2TPCBeginY);
    BDTReader.AddVariable("WC2TPCBeginZ", &bdt_WC2TPCBeginZ);
    BDTReader.AddVariable("WC2TPCEndX", &bdt_WC2TPCEndX);
    BDTReader.AddVariable("WC2TPCEndY", &bdt_WC2TPCEndY);
    BDTReader.AddVariable("WC2TPCEndZ", &bdt_WC2TPCEndZ);
    BDTReader.AddVariable("WC2TPCdEdx", &bdt_WC2TPCdEdx);

    Float_t bdt_recoTrkBeginX[BDT_NUM_RECO_TRKS];
    Float_t bdt_recoTrkBeginY[BDT_NUM_RECO_TRKS];
    Float_t bdt_recoTrkBeginZ[BDT_NUM_RECO_TRKS];
    Float_t bdt_recoTrkEndX[BDT_NUM_RECO_TRKS];
    Float_t bdt_recoTrkEndY[BDT_NUM_RECO_TRKS];
    Float_t bdt_recoTrkEndZ[BDT_NUM_RECO_TRKS];
    Float_t bdt_recoTrkLen[BDT_NUM_RECO_TRKS];
    Float_t bdt_recoTrkdEdx[BDT_NUM_RECO_TRKS];

    for (int j = 0; j < BDT_NUM_RECO_TRKS; ++j) {
        BDTReader.AddVariable(Form("recoTrkBeginX_%d", j), &bdt_recoTrkBeginX[j]);
        BDTReader.AddVariable(Form("recoTrkBeginY_%d", j), &bdt_recoTrkBeginY[j]);
        BDTReader.AddVariable(Form("recoTrkBeginZ_%d", j), &bdt_recoTrkBeginZ[j]);
        BDTReader.AddVariable(Form("recoTrkEndX_%d", j), &bdt_recoTrkEndX[j]);
        BDTReader.AddVariable(Form("recoTrkEndY_%d", j), &bdt_recoTrkEndY[j]);
        BDTReader.AddVariable(Form("recoTrkEndZ_%d", j), &bdt_recoTrkEndZ[j]);
        BDTReader.AddVariable(Form("recoTrkLen_%d", j), &bdt_recoTrkLen[j]);
        BDTReader.AddVariable(Form("recoTrkdEdx_%d", j), &bdt_recoTrkdEdx[j]);
    }
    
    BDTReader.AddVariable("numRecoTrksInCylinder", &bdt_numRecoTrksInCylinder);

    for (int i = 0; i < 4; ++i) {
        BDTReader.BookMVA("BDT" + std::to_string(i), "/exp/lariat/app/users/epelaez/analysis/python/model/chexch_abs_model_class_" + std::to_string(i) + ".xml");
    }

    //////////////////////
    // Loop over events //
    //////////////////////

    Int_t NumEntries = (Int_t) tree->GetEntries();
    std::cout << "Num entries: " << NumEntries << std::endl;

    int chExchAll               = 0;
    int chExchContainedRedVol   = 0;
    int chExchContainedCylinder = 0;

    int chExchContainedRedVolReco   = 0;
    int chExchContainedCylinderReco = 0;

    int scatteringsModified = 0;

    for (Int_t i = 0; i < NumEntries; ++i) {
        tree->GetEntry(i);

        // Events to debug BDT
        // if (!(
        //     event == 149610 || 
        //     event == 149626 ||
        //     event == 149639 ||
        //     event == 149650 ||
        //     event == 149655 ||
        //     event == 149666 ||
        //     event == 149694 ||
        //     event == 149696 ||
        //     event == 149702 ||
        //     event == 149706
        // )) continue;

        // Use first 50k events, last 50k events were used for BDT training
        if (i > 50000) break;

        // Make script go faster
        // if (i > 10000) break;

        // Sanity check
        removeRepeatedPoints(WC2TPCLocationsX, WC2TPCLocationsY, WC2TPCLocationsZ);

        // Extend reco cylinder
        std::vector<double>* wcX = new std::vector<double>(*WC2TPCLocationsX);
        std::vector<double>* wcY = new std::vector<double>(*WC2TPCLocationsY);
        std::vector<double>* wcZ = new std::vector<double>(*WC2TPCLocationsZ);

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

        int validPrimaryIdx = -1;
        for (size_t i = 0; i < primariesID->size(); ++i) {
            if (primariesID->at(i) == truthPrimaryID) {
                validPrimaryIdx = i;
                break;
            }
        }

        // Extend truth cylinder
        std::vector<double>* truthCylinderLocationX = new std::vector<double>(*truthPrimaryLocationX);
        std::vector<double>* truthCylinderLocationY = new std::vector<double>(*truthPrimaryLocationY);
        std::vector<double>* truthCylinderLocationZ = new std::vector<double>(*truthPrimaryLocationZ);

        if (truthCylinderLocationX->size() == 0) {
            truthCylinderLocationX->push_back(primariesStartX->at(validPrimaryIdx));
            truthCylinderLocationY->push_back(primariesStartY->at(validPrimaryIdx));
            truthCylinderLocationZ->push_back(primariesStartZ->at(validPrimaryIdx));

            truthCylinderLocationX->push_back(primariesEndX->at(validPrimaryIdx));
            truthCylinderLocationY->push_back(primariesEndY->at(validPrimaryIdx));
            truthCylinderLocationZ->push_back(primariesEndZ->at(validPrimaryIdx));
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

        ////////////////////////////////////////
        // Histograms for data-MC comparisons //
        ////////////////////////////////////////

        // For data-MC comparisons
        hMCNumWC2TPCMatch->Fill(WC2TPCsize);
        if (WC2TPCtrkID != -99999) {
            hMCPrimaryTrackPosition->Fill(WC2TPCPrimaryBeginX, WC2TPCPrimaryBeginY);

            bool isPrimaryTG = !isWithinReducedVolume(WC2TPCPrimaryEndX, WC2TPCPrimaryEndY, WC2TPCPrimaryEndZ);

            if (obtainedProbabilities) hMCShowerProb->Fill(showerProb);

            int smallTracksComp = 0;
            int tracksNearVertexComp = 0;
            int numTGTracksComp = 0;
            int smallTracksTPCStart = 0;
            int numTracksInCylinder = 0;

            for (size_t trk_idx = 0; trk_idx < recoBeginX->size(); ++trk_idx) {
                // Skip WC2TPC match itself
                if (recoTrkID->at(trk_idx) == WC2TPCtrkID) continue;

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
                if (startInCylinder && endInCylinder) numTracksInCylinder++;

                if (
                    !isWithinReducedVolume(recoBeginX->at(trk_idx), recoBeginY->at(trk_idx), recoBeginZ->at(trk_idx)) &&
                    !isWithinReducedVolume(recoEndX->at(trk_idx), recoEndY->at(trk_idx), recoEndZ->at(trk_idx))
                ) {
                    numTGTracksComp++;
                }

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
                if (isPrimaryTG) {
                    hMCTGTrackLengths->Fill(trackLength);
                }

                // Start of TPC
                if (
                    recoEndZ->at(trk_idx) < 30.0 && 
                    recoBeginZ->at(trk_idx) < 30.0
                ) {
                    if (trackLength < SMALL_TRACK_LENGTH_CHEX) smallTracksTPCStart++;
                }

                if (
                    !isPrimaryTG &&
                    (distanceFromStart < VERTEX_RADIUS || distanceFromEnd < VERTEX_RADIUS)
                ) {
                    tracksNearVertexComp++;
                    hMCTrackLengthsNearVertex->Fill(trackLength);
                }
                if (trackLength < SMALL_TRACK_LENGTH_CHEX) smallTracksComp++;
            }

            hMCBeforeShowerCutSmallTracks->Fill(smallTracksTPCStart);
            hMCNumTGTracks->Fill(numTGTracksComp);

            if (obtainedProbabilities && showerProb < SHOWER_PROB_CUT) hMCAfterShowerCutSmallTracks->Fill(smallTracksTPCStart);
            if (!isPrimaryTG) hMCTracksNearVertex->Fill(tracksNearVertexComp);

            if (isPrimaryTG) {
                hMCTGSmallTracks->Fill(smallTracksComp);
                hMCSmallVsTGTracks->Fill(smallTracksComp, numTGTracksComp);
                hMCNumTracksInCylinder->Fill(numTracksInCylinder);

                // Scan over small track length thresholds and fill 2D histogram
                for (int threshBin = 1; threshBin <= hMCTGNumSmallTracksVsThresh->GetNbinsX(); ++threshBin) {
                    double threshold = hMCTGNumSmallTracksVsThresh->GetXaxis()->GetBinCenter(threshBin);
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
                    hMCTGNumSmallTracksVsThresh->Fill(threshold, nSmallTracks);
                }
            }
        }

        ////////////////////////////
        // Save truth information //
        ////////////////////////////

        // Scattering only if degree > THRESHOLD
        double scatteringAngle = -9999;
        if (backgroundType == 12 || backgroundType == 6) {
            if (backgroundType == 12) {
                scatteringAngle = trajectoryInteractionAngle;
            } else if (backgroundType == 6) {
                scatteringAngle = truthScatteringAngle;
            }

            // std::cout << "Initial scattering angle: " << scatteringAngle * (180 / TMath::Pi()) << " for interaction " << backgroundType << std::endl;

            if (scatteringAngle < SCATTERING_ANGLE_THRESHOLD) {
                scatteringsModified++;

                // Use secondary interaction
                for (int iInteraction = 0; iInteraction < secondaryInteractionTypes->size(); ++iInteraction) {
                    int currentInteraction = secondaryInteractionTypes->at(iInteraction);
                    scatteringAngle        = secondaryInteractionAngle->at(iInteraction);

                    for (int iContribution = 0; iContribution < secondaryIncidentKEContributions->at(iInteraction).size(); ++iContribution) {
                        trueIncidentKEContributions->push_back(secondaryIncidentKEContributions->at(iInteraction)[iContribution]);
                    }

                    // std::cout << "  Current interaction: " << currentInteraction << ", scattering angle: " << scatteringAngle * (180 / TMath::Pi()) << std::endl;

                    if (
                        (currentInteraction == 6 || currentInteraction == 12) &&
                        scatteringAngle < SCATTERING_ANGLE_THRESHOLD
                    ) {
                        // We want to keep going
                        continue;
                    } else {
                        // We found non-scattering interaction or scattering interaction 
                        // with angle above our threshold value
                        backgroundType       = currentInteraction;
                        truthPrimaryVertexKE = secondaryInteractionInteractingKE->at(iInteraction);

                        truthPrimaryVertexX = secondaryInteractionXPosition->at(iInteraction);
                        truthPrimaryVertexY = secondaryInteractionYPosition->at(iInteraction);
                        truthPrimaryVertexZ = secondaryInteractionZPosition->at(iInteraction);
                        break;
                    }
                }
                // std::cout << "new interaction : " << backgroundType << std::endl;
                // std::cout << std::endl;
            } else {
                // If initial interaction has valid angle, just correct interaction KE for 
                // elastic scattering interactions
                if (backgroundType == 12) truthPrimaryVertexKE = trajectoryInteractionKE;
            }
        }

        // More information about charge exchange events
        bool chExchInsideRedVol   = false;
        bool chExchInsideCylinder = false;
        if (backgroundType == 7) {
            // If the photons end outside the TPC, there is no way for us to detect the charge exchange
            // If the photons end outside the cylinder, it is also pretty hard for us to detect it reliably
            chExchAll++;

            bool firstPhotonInsideRedVol  = false;
            bool secondPhotonInsideRedVol = false;

            bool firstPhotonInsideCylinder  = false;
            bool secondPhotonInsideCylinder = false;

            int  photonCount = 0;

            for (int i = 0; i < chExchShowerIDs->size(); ++i) {
                if (std::find(chExchShowerNeutralPionDaughtersID->begin(), chExchShowerNeutralPionDaughtersID->end(), chExchShowerIDs->at(i)) != chExchShowerNeutralPionDaughtersID->end()) {
                    // Found direct daughter of neutral pion
                    if (chExchShowerPDGs->at(i) == 22) {
                        // Found photon from neutral pion decay
                        photonCount++;

                        bool endInCylinder = IsPointInsideTrackCylinder(
                            truthCylinderLocationX, truthCylinderLocationY, truthCylinderLocationZ,
                            chExchShowerEnd->at(i)[0], chExchShowerEnd->at(i)[1], chExchShowerEnd->at(i)[2],
                            CYLINDER_RADIUS
                        );
                        bool endInRedVol = isWithinReducedVolume(chExchShowerEnd->at(i)[0], chExchShowerEnd->at(i)[1], chExchShowerEnd->at(i)[2]);

                        if (photonCount == 1) {
                            firstPhotonInsideCylinder = endInCylinder;
                            firstPhotonInsideRedVol   = endInRedVol;
                        } else if (photonCount == 2) {
                            secondPhotonInsideCylinder = endInCylinder;
                            secondPhotonInsideRedVol   = endInRedVol;
                        }
                    }
                }
            }

            // We can only hope to reconstruct it if there is a WC2TPC match
            if (WC2TPCtrkID != -99999) {
                if (firstPhotonInsideRedVol || secondPhotonInsideRedVol) {
                    chExchContainedRedVol++;
                    chExchInsideRedVol = true;
                }

                if (firstPhotonInsideCylinder || secondPhotonInsideCylinder) {
                    chExchContainedCylinder++;
                    chExchInsideCylinder = true;
                }
            }
        }

        // Further classify muon backgrounds
        int muonType = -1; // 0: through-going, 1: decaying, 2: capture at rest
        if (backgroundType == 2) {
            if (!isWithinReducedVolume(truthPrimaryVertexX, truthPrimaryVertexY, truthPrimaryVertexZ)) {
                // through-going muon
                muonType = 0;
            } else {
                bool sawElectron = false;
                for (int i = 0; i < truthPrimaryDaughtersPDG->size(); ++i) {
                    if (truthPrimaryDaughtersProcess->at(i) == "muMinusCaptureAtRest") {
                        // muon capture at rest
                        muonType = 2;
                        break;
                    }
                    if (truthPrimaryDaughtersPDG->at(i) == 11) {
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

            hTrueMuonTypes->Fill(muonType);
        }

        // Get true energy bin
        int TrueEnergyBin = getBin(truthPrimaryVertexKE * 1000, ARRAY_KE_BINS);

        // Add true incident KE
        if (validTrueIncidentKE) {
            for (double x : *trueIncidentKEContributions) {
                hTrueIncidentKE->Fill(x);
            }
        }

        // In this script, we replace background types 6 and 12 with background types 13 and 14,
        // because we are interested in contamination across the bins that we are measuring
        // if (backgroundType == 12 || (backgroundType == 6 && numVisibleProtons == 0)) {
        //     backgroundType = 13;
        // } else if (backgroundType == 6 && numVisibleProtons > 0) {
        //     backgroundType = 14;
        // }

        // Fill in true interaction energy histograms
        if (backgroundType == 0) {
            hTrueAbs0pKE->Fill(truthPrimaryVertexKE * 1000);
        } else if (backgroundType == 1) {
            hTrueAbsNpKE->Fill(truthPrimaryVertexKE * 1000);
        } else if (backgroundType == 6 || backgroundType == 12) {
            hTrueScatterKE->Fill(truthPrimaryVertexKE * 1000);
        } else if (backgroundType == 7) {
            hTrueChExchKE->Fill(truthPrimaryVertexKE * 1000);
        }

        hTotalEvents->Fill(backgroundType);

        // If no track matched to wire-chamber, skip
        if (WC2TPCtrkID == -99999) {
            if (backgroundType == 0) {
                hTrueAbs0pKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueAbs0pKERejDataProds->Fill(truthPrimaryVertexKE * 1000);
            } else if (backgroundType == 1) {
                hTrueAbsNpKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueAbsNpKERejDataProds->Fill(truthPrimaryVertexKE * 1000);
            } else if (backgroundType == 6 || backgroundType == 12) {
                hTrueScatterKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueScatterKERejDataProds->Fill(truthPrimaryVertexKE * 1000);
            } else if (backgroundType == 7) {
                hTrueChExchKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueChExchKERejDataProds->Fill(truthPrimaryVertexKE * 1000);
            }
            continue;
        }
        hDataProdsAndWC2TPC->Fill(backgroundType);

        //////////////////////////////////////
        // Back to classification algorithm //
        //////////////////////////////////////

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

        //////////////////
        // Cylinder cut //
        //////////////////

        int numSmallTracksInCylinder = 0;
        int numTracksInCylinder      = 0;
        for (int trk_idx = 0; trk_idx < recoTrkID->size(); ++trk_idx) {
            if (recoTrkID->at(trk_idx) == WC2TPCtrkID) continue;

            double trackLength = sqrt(
                pow(recoEndX->at(trk_idx) - recoBeginX->at(trk_idx), 2) +
                pow(recoEndY->at(trk_idx) - recoBeginY->at(trk_idx), 2) +
                pow(recoEndZ->at(trk_idx) - recoBeginZ->at(trk_idx), 2)
            );

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

            if (startInCylinder && endInCylinder) numTracksInCylinder++;

            if (
                startInCylinder && endInCylinder &&
                trackLength < CYLINDER_SMALL_TRACK
            ) {
                numSmallTracksInCylinder++;
            }
        }

        // Perform cut

        // Cut on tracks
        // if (numTracksInCylinder > ALLOWED_CYLINDER_TRACKS) {
        // Cut on small tracks
        if (numSmallTracksInCylinder > ALLOWED_CYLINDER_SMALL_TRACKS) {
            if (backgroundType == 0) {
                hTrueAbs0pKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueAbs0pKERejElectron->Fill(truthPrimaryVertexKE * 1000);
            } else if (backgroundType == 1) {
                hTrueAbsNpKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueAbsNpKERejElectron->Fill(truthPrimaryVertexKE * 1000);
            } else if (backgroundType == 6 || backgroundType == 12) {
                hTrueScatterKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueScatterKERejElectron->Fill(truthPrimaryVertexKE * 1000);
            } else if (backgroundType == 7) {
                hTrueChExchKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueChExchKERejElectron->Fill(truthPrimaryVertexKE * 1000);
            }
            continue;
        }
        hNotAnElectron->Fill(backgroundType);

        ////////////////////////
        // Neural network cut //
        ////////////////////////

        // If did not obtain neural network electron shower probabilities, skip 
        // if (!obtainedProbabilities) showerProb = 1.;
        // if (showerProb >= SHOWER_PROB_CUT){
        //     if (backgroundType == 0) {
        //         hTrueAbs0pKERejected->Fill(truthPrimaryVertexKE * 1000);
        //         hTrueAbs0pKERejElectron->Fill(truthPrimaryVertexKE * 1000);
        //     } else if (backgroundType == 1) {
        //         hTrueAbsNpKERejected->Fill(truthPrimaryVertexKE * 1000);
        //         hTrueAbsNpKERejElectron->Fill(truthPrimaryVertexKE * 1000);
        //     } else if (backgroundType == 6 || backgroundType == 12) {
        //         hTrueScatterKERejected->Fill(truthPrimaryVertexKE * 1000);
        //         hTrueScatterKERejElectron->Fill(truthPrimaryVertexKE * 1000);
        //     } else if (backgroundType == 7) {
        //         hTrueChExchKERejected->Fill(truthPrimaryVertexKE * 1000);
        //         hTrueChExchKERejElectron->Fill(truthPrimaryVertexKE * 1000);
        //     }
        //     continue;
        // }
        // hNotAnElectron->Fill(backgroundType);

        ///////////////////////
        // Setup BDT for use //
        ///////////////////////

        getBDTVariables(
            WC2TPCtrkID, 
            WC2TPCLocationsX, WC2TPCLocationsY, WC2TPCLocationsZ,
            wcX, wcY, wcZ,
            recoBeginX, recoBeginY, recoBeginZ,
            recoEndX,   recoEndY,   recoEndZ,
            recoDEDX, recoTrkID, BDT_NUM_RECO_TRKS,
            &bdt_WC2TPCTrackLength,
            &bdt_WC2TPCBeginX, &bdt_WC2TPCBeginY, &bdt_WC2TPCBeginZ,
            &bdt_WC2TPCEndX,   &bdt_WC2TPCEndY,   &bdt_WC2TPCEndZ,
            &bdt_WC2TPCdEdx,
            bdt_recoTrkBeginX, bdt_recoTrkBeginY, bdt_recoTrkBeginZ,
            bdt_recoTrkEndX,   bdt_recoTrkEndY,   bdt_recoTrkEndZ,
            bdt_recoTrkLen,    bdt_recoTrkdEdx,
            &bdt_numRecoTrksInCylinder
        );

        // Compute BDT scores
        double s0 = BDTReader.EvaluateMVA("BDT0");
        double s1 = BDTReader.EvaluateMVA("BDT1");
        double s2 = BDTReader.EvaluateMVA("BDT2");
        double s3 = BDTReader.EvaluateMVA("BDT3");

        // Make sure raw output is in (-1, 1)
        const double eps = 1e-12; 
        s0 = std::max(-1.0 + eps, std::min(1.0 - eps, s0));
        s1 = std::max(-1.0 + eps, std::min(1.0 - eps, s1));
        s2 = std::max(-1.0 + eps, std::min(1.0 - eps, s2));
        s3 = std::max(-1.0 + eps, std::min(1.0 - eps, s3));
        
        // Map scores in (-1,1) -> margins in (-inf, +inf) via atanh
        double m0 = 0.5 * std::log((1.0 + s0) / (1.0 - s0)); // atanh(s0)
        double m1 = 0.5 * std::log((1.0 + s1) / (1.0 - s1)); // atanh(s1)
        double m2 = 0.5 * std::log((1.0 + s2) / (1.0 - s2)); // atanh(s2)
        double m3 = 0.5 * std::log((1.0 + s3) / (1.0 - s3)); // atanh(s3)

        // Convert margins -> class probabilities via softmax
        double mmax = std::max({m0, m1, m2, m3});
        double e0 = std::exp(m0 - mmax);
        double e1 = std::exp(m1 - mmax);
        double e2 = std::exp(m2 - mmax);
        double e3 = std::exp(m3 - mmax);
        double Z  = e0 + e1 + e2 + e3;
        if (Z == 0) Z = 1e-12;

        double p0 = e0 / Z;
        double p1 = e1 / Z;
        double p2 = e2 / Z;
        double p3 = e3 / Z;

        // Test few events
        if (
            event == 149610 || 
            event == 149626 ||
            event == 149639 ||
            event == 149650 ||
            event == 149655 ||
            event == 149666 ||
            event == 149694 ||
            event == 149696 ||
            event == 149702 ||
            event == 149706
        ) {
            // std::cout << "Event " << event << " BDT scores: " << p0 << ", " << p1 << ", " << p2 << ", " << p3 << std::endl;
            // std::cout << "           " << "Raw scores: " << s0 << ", " << s1 << ", " << s2 << ", " << s3 << std::endl;

            // std::cout << "BDT input variables for event " << event << ":" << std::endl;
            // std::cout << "  WC2TPCTrackLength: " << bdt_WC2TPCTrackLength << std::endl;
            // std::cout << "  WC2TPCBeginX: " << bdt_WC2TPCBeginX << std::endl;
            // std::cout << "  WC2TPCBeginY: " << bdt_WC2TPCBeginY << std::endl;
            // std::cout << "  WC2TPCBeginZ: " << bdt_WC2TPCBeginZ << std::endl;
            // std::cout << "  WC2TPCEndX: " << bdt_WC2TPCEndX << std::endl;
            // std::cout << "  WC2TPCEndY: " << bdt_WC2TPCEndY << std::endl;
            // std::cout << "  WC2TPCEndZ: " << bdt_WC2TPCEndZ << std::endl;
            // std::cout << "  WC2TPCdEdx: " << bdt_WC2TPCdEdx << std::endl;
            // for (int j = 0; j < BDT_NUM_RECO_TRKS; ++j) {
            //     std::cout << "  recoTrkBeginX[" << j << "]: " << bdt_recoTrkBeginX[j] << std::endl;
            //     std::cout << "  recoTrkBeginY[" << j << "]: " << bdt_recoTrkBeginY[j] << std::endl;
            //     std::cout << "  recoTrkBeginZ[" << j << "]: " << bdt_recoTrkBeginZ[j] << std::endl;
            //     std::cout << "  recoTrkEndX[" << j << "]: " << bdt_recoTrkEndX[j] << std::endl;
            //     std::cout << "  recoTrkEndY[" << j << "]: " << bdt_recoTrkEndY[j] << std::endl;
            //     std::cout << "  recoTrkEndZ[" << j << "]: " << bdt_recoTrkEndZ[j] << std::endl;
            //     std::cout << "  recoTrkLen[" << j << "]: " << bdt_recoTrkLen[j] << std::endl;
            //     std::cout << "  recoTrkdEdx[" << j << "]: " << bdt_recoTrkdEdx[j] << std::endl;
            // }
            // std::cout << "  numRecoTrksInCylinder: " << bdt_numRecoTrksInCylinder << std::endl;
        }

        ////////////////////////
        // Reduced volume cut //
        ////////////////////////

        if (!isWithinReducedVolume(breakPointX, breakPointY, breakPointZ)) {
            if (backgroundType == 0) {
                hTrueAbs0pKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueAbs0pKERejRedVol->Fill(truthPrimaryVertexKE * 1000);
            } else if (backgroundType == 1) {
                hTrueAbsNpKERejected->Fill(truthPrimaryVertexKE * 1000);
                hTrueAbsNpKERejRedVol->Fill(truthPrimaryVertexKE * 1000);
            } else if (backgroundType == 6 || backgroundType == 12) {
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
            } else if (backgroundType == 6 || backgroundType == 12) {
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
            if (recoTrkID->at(trk_idx) == WC2TPCtrkID) continue;

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
        // If the new secondary track looks like a pion, the stitching is non-sensical, and we should revert
        bool newSecondaryPion = false;
        if (minChi2 == minStitchedChi2) {
            std::vector<double> newSecondaryResR(wcMatchResR->begin(), wcMatchResR->begin() + bestBreakPoint);
            std::vector<double> newSecondaryDEDX(wcMatchDEDX->begin(), wcMatchDEDX->begin() + bestBreakPoint);

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
                if (backgroundType == 0) {
                    hTrueAbs0pKERejected->Fill(truthPrimaryVertexKE * 1000);
                    hTrueAbs0pKERejManyPions->Fill(truthPrimaryVertexKE * 1000);
                } else if (backgroundType == 1) {
                    hTrueAbsNpKERejected->Fill(truthPrimaryVertexKE * 1000);
                    hTrueAbsNpKERejManyPions->Fill(truthPrimaryVertexKE * 1000);
                } else if (backgroundType == 6 || backgroundType == 12) {
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
            } else if (backgroundType == 6 || backgroundType == 12) {
                hPionScatterKETrue->Fill(energyAtVertex);
            } else if (backgroundType == 2) {
                hPionScatterKEMuon->Fill(energyAtVertex);

                if (muonType == 0) {
                    hPionScatterKEMuonTG->Fill(energyAtVertex);
                } else if (muonType == 1) {
                    hPionScatterKEMuonDecay->Fill(energyAtVertex);
                } else if (muonType == 2) {
                    hPionScatterKEMuonCAR->Fill(energyAtVertex);
                }
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
            } else if (backgroundType == 6 || backgroundType == 12) {
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
            } else if (backgroundType == 6 || backgroundType == 12) {
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
            } else if (backgroundType == 6 || backgroundType == 12) {
                hTrueScatterKEAsAbsNp->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueScatterAsByBin.at(TrueEnergyBin).at(1)->Fill(energyAtVertex);
            } else if (backgroundType == 7) {
                hTrueChExchKEAsAbsNp->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueChExchAsByBin.at(TrueEnergyBin).at(1)->Fill(energyAtVertex);
            }

            continue;
        }
        hNotPionAbsNp->Fill(backgroundType);

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
        int numLargeClusters = 0;
        for (int i = 0; i < hitClusters.size(); ++i) {
            double maxWire = *std::max_element(hitClusters[i].hitW.begin(), hitClusters[i].hitW.end());
            double minWire = *std::min_element(hitClusters[i].hitW.begin(), hitClusters[i].hitW.end());
            if (std::abs(maxWire - minWire) > LARGE_CLUSTER_THRESHOLD) numLargeClusters++;

            double clusterSize = std::abs(maxWire - minWire);
            if (backgroundType == 0) {
                hHitClusterSizesAbs0p->Fill(clusterSize);
            } else if (backgroundType == 1) {
                hHitClusterSizesAbsNp->Fill(clusterSize);
            } else if (backgroundType == 2) {
                hHitClusterSizesMuon->Fill(clusterSize);
            } else if (backgroundType == 3) {
                hHitClusterSizesElectron->Fill(clusterSize);
            } else if (backgroundType == 6 || backgroundType == 12) {
                hHitClusterSizesScatter->Fill(clusterSize);
            } else if (backgroundType == 7) {
                hHitClusterSizesChExch->Fill(clusterSize);
            } else {
                hHitClusterSizesOther->Fill(clusterSize);
            }
        }

        if (backgroundType == 0) {
            hLargeHitClusterAbs0p->Fill(numLargeClusters);
        } else if (backgroundType == 1) {
            hLargeHitClusterAbsNp->Fill(numLargeClusters);
        } else if (backgroundType == 7) {
            hLargeHitClusterChExch->Fill(numLargeClusters);
        } else if (backgroundType == 6 || backgroundType == 12) {
            hLargeHitClusterScatter->Fill(numLargeClusters);
        } else if (backgroundType == 2) {
            hLargeHitClusterMuon->Fill(numLargeClusters);
        } else if (backgroundType == 3) {
            hLargeHitClusterElectron->Fill(numLargeClusters);
        } else {
            hLargeHitClusterOther->Fill(numLargeClusters);
        }

        // if (numLargeClusters < NUM_CLUSTERS_THRESHOLD) {
        //     hPionAbs0p->Fill(backgroundType);

        //     hPionAbs0pKE->Fill(energyAtVertex);
        //     if (backgroundType == 0) {
        //         hPionAbs0pKETrue->Fill(energyAtVertex);
        //     } else if (backgroundType == 1) {
        //         hPionAbs0pKEAbsNp->Fill(energyAtVertex);
        //     } else if (backgroundType == 7) {
        //         hPionAbs0pKEChExch->Fill(energyAtVertex);
        //     } else if (backgroundType == 6 || backgroundType == 12) {
        //         hPionAbs0pKEScatter->Fill(energyAtVertex);
        //     } else if (backgroundType == 2) {
        //         hPionAbs0pKEMuon->Fill(energyAtVertex);
        //     } else if (backgroundType == 3) {
        //         hPionAbs0pKEElectron->Fill(energyAtVertex);
        //     } else {
        //         hPionAbs0pKEOther->Fill(energyAtVertex);
        //     }

        //     if (backgroundType == 0) {
        //         hTrueAbs0pKEAsAbs0p->Fill(truthPrimaryVertexKE * 1000);
        //         if (TrueEnergyBin != -1) TrueAbs0pAsByBin.at(TrueEnergyBin).at(0)->Fill(energyAtVertex);
        //     } else if (backgroundType == 1) {
        //         hTrueAbsNpKEAsAbs0p->Fill(truthPrimaryVertexKE * 1000);
        //         if (TrueEnergyBin != -1) TrueAbsNpAsByBin.at(TrueEnergyBin).at(0)->Fill(energyAtVertex);
        //     } else if (backgroundType == 6 || backgroundType == 12) {
        //         hTrueScatterKEAsAbs0p->Fill(truthPrimaryVertexKE * 1000);
        //         if (TrueEnergyBin != -1) TrueScatterAsByBin.at(TrueEnergyBin).at(0)->Fill(energyAtVertex);
        //     } else if (backgroundType == 7) {
        //         hTrueChExchKEAsAbs0p->Fill(truthPrimaryVertexKE * 1000);
        //         if (TrueEnergyBin != -1) TrueChExchAsByBin.at(TrueEnergyBin).at(0)->Fill(energyAtVertex);
        //     }

        //     continue;
        // }
        // hNotPionAbs0p->Fill(backgroundType);

        /////////////////////////////////////////////////
        // Discriminate remaining from charge exchange //
        /////////////////////////////////////////////////

        double max_prob = std::max({p0, p1, p2, p3});

        if (p0 == max_prob) {
            // classify as pion abs 0p
            if (p0 < 0.5) continue;

            hPionAbs0p->Fill(backgroundType);

            hPionAbs0pKE->Fill(energyAtVertex);
            if (backgroundType == 0) {
                hPionAbs0pKETrue->Fill(energyAtVertex);
                hBestAbs0pAbs0p->Fill(p0);
            } else if (backgroundType == 1) {
                hPionAbs0pKEAbsNp->Fill(energyAtVertex);
                hBestAbs0pAbsNp->Fill(p0);
            } else if (backgroundType == 7) {
                hPionAbs0pKEChExch->Fill(energyAtVertex);
                hBestAbs0pChExch->Fill(p0);
            } else if (backgroundType == 6 || backgroundType == 12) {
                hPionAbs0pKEScatter->Fill(energyAtVertex);
                hBestAbs0pScatter->Fill(p0);
            } else if (backgroundType == 2) {
                hPionAbs0pKEMuon->Fill(energyAtVertex);
                hBestAbs0pMuon->Fill(p0);

                if (muonType == 0) {
                    hPionAbs0pKEMuonTG->Fill(energyAtVertex);
                } else if (muonType == 1) {
                    hPionAbs0pKEMuonDecay->Fill(energyAtVertex);
                } else if (muonType == 2) {
                    hPionAbs0pKEMuonCAR->Fill(energyAtVertex);
                }
            } else if (backgroundType == 3) {
                hPionAbs0pKEElectron->Fill(energyAtVertex);
                hBestAbs0pElectron->Fill(p0);
            } else {
                hPionAbs0pKEOther->Fill(energyAtVertex);
                hBestAbs0pOther->Fill(p0);
            }

            if (backgroundType == 0) {
                hTrueAbs0pKEAsAbs0p->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueAbs0pAsByBin.at(TrueEnergyBin).at(0)->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                hTrueAbsNpKEAsAbs0p->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueAbsNpAsByBin.at(TrueEnergyBin).at(0)->Fill(energyAtVertex);
            } else if (backgroundType == 6 || backgroundType == 12) {
                hTrueScatterKEAsAbs0p->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueScatterAsByBin.at(TrueEnergyBin).at(0)->Fill(energyAtVertex);
            } else if (backgroundType == 7) {
                hTrueChExchKEAsAbs0p->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueChExchAsByBin.at(TrueEnergyBin).at(0)->Fill(energyAtVertex);
            }

            outFileAbs0pBkg << "Run: " << run << " subrun: " << subrun << " event: " << event << std::endl;
            outFileAbs0pBkg << " Classified as: " << backgroundTypes[backgroundType] << " with energy at vertex: " << truthPrimaryVertexKE * 1000 << std::endl;
            outFileAbs0pBkg << " Scattering angle (degrees): " << scatteringAngle * (180 / TMath::Pi()) << std::endl;
            outFileAbs0pBkg << " Tagged pions: " << totalTaggedPions << std::endl;
            outFileAbs0pBkg << " Tagged protons: " << totalTaggedProtons << std::endl;
            outFileAbs0pBkg << " Reco trks in cylinder: " << bdt_numRecoTrksInCylinder << std::endl;
            outFileAbs0pBkg << std::endl;

            continue;
        } else if (p1 == max_prob) {
            // classify as electron shower
            if (backgroundType == 0) {
                hBestElectronAbs0p->Fill(p1);
            } else if (backgroundType == 1) {
                hBestElectronAbsNp->Fill(p1);
            } else if (backgroundType == 7) {
                hBestElectronChExch->Fill(p1);
            } else if (backgroundType == 6 || backgroundType == 12) {
                hBestElectronScatter->Fill(p1);
            } else if (backgroundType == 2) {
                hBestElectronMuon->Fill(p1);
            } else if (backgroundType == 3) {
                hBestElectronElectron->Fill(p1);
            } else {
                hBestElectronOther->Fill(p1);
            }

            continue;
        } else if (p2 == max_prob) {
            // classify as charge exchange
            hPionChExch->Fill(backgroundType);

            hPionChExchKE->Fill(energyAtVertex);
            if (backgroundType == 0) {
                hPionChExchKEAbs0p->Fill(energyAtVertex);
                hBestChExchAbs0p->Fill(p2);
            } else if (backgroundType == 1) {
                hPionChExchKEAbsNp->Fill(energyAtVertex);
                hBestChExchAbsNp->Fill(p2);
            } else if (backgroundType == 7) {
                hPionChExchKETrue->Fill(energyAtVertex);
                hBestChExchChExch->Fill(p2);
            } else if (backgroundType == 6 || backgroundType == 12) {
                hPionChExchKEScatter->Fill(energyAtVertex);
                hBestChExchScatter->Fill(p2);
            } else if (backgroundType == 2) {
                hPionChExchKEMuon->Fill(energyAtVertex);
                hBestChExchMuon->Fill(p2);
            } else if (backgroundType == 3) {
                hPionChExchKEElectron->Fill(energyAtVertex);
                hBestChExchElectron->Fill(p2);
            } else {
                hPionChExchKEOther->Fill(energyAtVertex);
                hBestChExchOther->Fill(p2);
            }

            if (backgroundType == 0) {
                hTrueAbs0pKEAsChExch->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueAbs0pAsByBin.at(TrueEnergyBin).at(3)->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                hTrueAbsNpKEAsChExch->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueAbsNpAsByBin.at(TrueEnergyBin).at(3)->Fill(energyAtVertex);
            } else if (backgroundType == 6 || backgroundType == 12) {
                hTrueScatterKEAsChExch->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueScatterAsByBin.at(TrueEnergyBin).at(3)->Fill(energyAtVertex);
            } else if (backgroundType == 7) {
                hTrueChExchKEAsChExch->Fill(truthPrimaryVertexKE * 1000);
                if (TrueEnergyBin != -1) TrueChExchAsByBin.at(TrueEnergyBin).at(3)->Fill(energyAtVertex);
            }

            if (chExchInsideCylinder) {
                chExchContainedCylinderReco++;
            } else if (chExchInsideRedVol) {
                chExchContainedRedVolReco++;
            }

            continue;
        } else if (p3 == max_prob) {
            // classify as remaining scattering
            if (backgroundType == 0) {
                hBestScatterAbs0p->Fill(p3);
            } else if (backgroundType == 1) {
                hBestScatterAbsNp->Fill(p3);
            } else if (backgroundType == 7) {
                hBestScatterChExch->Fill(p3);
            } else if (backgroundType == 6 || backgroundType == 12) {
                hBestScatterScatter->Fill(p3);
            } else if (backgroundType == 2) {
                hBestScatterMuon->Fill(p3);
            } else if (backgroundType == 3) {
                hBestScatterElectron->Fill(p3);
            } else {
                hBestScatterOther->Fill(p3);
            }

            continue;
        }

        /////////////////////
        // Small track cut //
        /////////////////////

        // if (numSmallTracks > SMALL_TRACK_CUT_CHEX) {
        //     // Tag as charge exchange
        //     hPionChExch->Fill(backgroundType);

        //     hPionChExchKE->Fill(energyAtVertex);
        //     if (backgroundType == 0) {
        //         hPionChExchKEAbs0p->Fill(energyAtVertex);
        //     } else if (backgroundType == 1) {
        //         hPionChExchKEAbsNp->Fill(energyAtVertex);
        //     } else if (backgroundType == 7) {
        //         hPionChExchKETrue->Fill(energyAtVertex);
        //     } else if (backgroundType == 6 || backgroundType == 12) {
        //         hPionChExchKEScatter->Fill(energyAtVertex);
        //     } else if (backgroundType == 2) {
        //         hPionChExchKEMuon->Fill(energyAtVertex);
        //     } else if (backgroundType == 3) {
        //         hPionChExchKEElectron->Fill(energyAtVertex);
        //     } else {
        //         hPionChExchKEOther->Fill(energyAtVertex);
        //     }

        //     if (backgroundType == 0) {
        //         hTrueAbs0pKEAsChExch->Fill(truthPrimaryVertexKE * 1000);
        //         if (TrueEnergyBin != -1) TrueAbs0pAsByBin.at(TrueEnergyBin).at(3)->Fill(energyAtVertex);
        //     } else if (backgroundType == 1) {
        //         hTrueAbsNpKEAsChExch->Fill(truthPrimaryVertexKE * 1000);
        //         if (TrueEnergyBin != -1) TrueAbsNpAsByBin.at(TrueEnergyBin).at(3)->Fill(energyAtVertex);
        //     } else if (backgroundType == 6 || backgroundType == 12) {
        //         hTrueScatterKEAsChExch->Fill(truthPrimaryVertexKE * 1000);
        //         if (TrueEnergyBin != -1) TrueScatterAsByBin.at(TrueEnergyBin).at(3)->Fill(energyAtVertex);
        //     } else if (backgroundType == 7) {
        //         hTrueChExchKEAsChExch->Fill(truthPrimaryVertexKE * 1000);
        //         if (TrueEnergyBin != -1) TrueChExchAsByBin.at(TrueEnergyBin).at(3)->Fill(energyAtVertex);
        //     }

        //     outFileChExchBkg << "Run: " << run << " subrun: " << subrun << " event: " << event << std::endl;
        //     outFileChExchBkg << " Classified as: " << backgroundTypes[backgroundType] << " with energy at vertex: " << truthPrimaryVertexKE * 1000 << std::endl;
        //     outFileChExchBkg << std::endl;

        //     continue;
        // }
        // hNotChExch->Fill(backgroundType);

        // Anything left here is rejected
        if (backgroundType == 0) {
            hTrueAbs0pKERejected->Fill(truthPrimaryVertexKE * 1000);
            hTrueAbs0pKERejClusters->Fill(truthPrimaryVertexKE * 1000);
        } else if (backgroundType == 1) {
            hTrueAbsNpKERejected->Fill(truthPrimaryVertexKE * 1000);
            hTrueAbsNpKERejClusters->Fill(truthPrimaryVertexKE * 1000);
        } else if (backgroundType == 6 || backgroundType == 12) {
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

    std::cout << std::endl;
    std::cout << "Original muon event composition: " << std::endl;
    std::cout << "  Total muon events: " << hTrueMuonTypes->Integral() << std::endl;
    std::cout << "  Through-going:     " << hTrueMuonTypes->GetBinContent(1) << std::endl;
    std::cout << "  Decay:             " << hTrueMuonTypes->GetBinContent(2) << std::endl;
    std::cout << "  Capture at rest:   " << hTrueMuonTypes->GetBinContent(3) << std::endl;
    std::cout << "  Other:             " << hTrueMuonTypes->GetBinContent(4) << std::endl;

    std::cout << std::endl;
    std::cout << "Total charge exchange events: " << chExchAll << std::endl;
    std::cout << "  At least one photon in reduced volume: " << chExchContainedRedVol << ", " << ((double) chExchContainedRedVol / (double) chExchAll) * 100 << "%" << std::endl;
    std::cout << "  At least one photon contained in reduced volume correctly tagged: " << chExchContainedRedVolReco << ", " << ((double) chExchContainedRedVolReco / (double) chExchContainedRedVol) * 100 << "%" << std::endl;
    std::cout << std::endl;
    std::cout << "  At least one photon in cylinder: " << chExchContainedCylinder << ", " << ((double) chExchContainedCylinder / (double) chExchAll) * 100 << "%" << std::endl;
    std::cout << "  At least one photon contained in cylinder correctly tagged: " << chExchContainedCylinderReco << ", " << ((double) chExchContainedCylinderReco / (double) chExchContainedCylinder) * 100 << "%" << std::endl;

    std::cout << std::endl;
    std::cout << "Number of scattering events modified: " << scatteringsModified << std::endl;

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
    std::cout << "Events not classified as pion 0p absorption (too many non-reconstructed hit clusters): " << std::endl;
    printBackgroundInfo(hNotPionAbs0p, std::cout);

    // std::cout << std::endl;
    // std::cout << "Events not classified as charge exchange (not enough small tracks): " << std::endl;
    // printBackgroundInfo(hNotChExch, std::cout);

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
    std::cout << "Purity: " << (hPionScatter->GetBinContent(7) + hPionScatter->GetBinContent(13)) / hPionScatter->Integral() << std::endl;
    std::cout << "Efficiency: " << (hPionScatter->GetBinContent(7) + hPionScatter->GetBinContent(13)) / (hTotalEvents->GetBinContent(7) + hTotalEvents->GetBinContent(13)) << std::endl;

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

    hMCNumTGTracks->SetDirectory(comparisonsFile);
    hMCNumTGTracks->Write();

    hMCSmallVsTGTracks->SetDirectory(comparisonsFile);
    hMCSmallVsTGTracks->Write();

    hMCShowerProb->SetDirectory(comparisonsFile);
    hMCShowerProb->Write();

    hMCBeforeShowerCutSmallTracks->SetDirectory(comparisonsFile);
    hMCBeforeShowerCutSmallTracks->Write();

    hMCAfterShowerCutSmallTracks->SetDirectory(comparisonsFile);
    hMCAfterShowerCutSmallTracks->Write();

    hMCTGTrackLengths->SetDirectory(comparisonsFile);
    hMCTGTrackLengths->Write();

    hMCTGNumSmallTracksVsThresh->SetDirectory(comparisonsFile);
    hMCTGNumSmallTracksVsThresh->Write();

    hMCNumWC2TPCMatch->SetDirectory(comparisonsFile);
    hMCNumWC2TPCMatch->Write();

    hMCNumTracksInCylinder->SetDirectory(comparisonsFile);
    hMCNumTracksInCylinder->Write();

    hMCPrimaryTrackPosition->SetDirectory(comparisonsFile);
    hMCPrimaryTrackPosition->Write();

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
    TMatrixD AddSmearInverse(NUM_SIGNAL_TYPES * NUM_BINS_KE, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    TVectorD WF(NUM_SIGNAL_TYPES * NUM_BINS_KE);
    TMatrixD UnfoldCov(NUM_SIGNAL_TYPES * NUM_BINS_KE, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    TMatrixD CovRotation(NUM_SIGNAL_TYPES * NUM_BINS_KE, NUM_SIGNAL_TYPES * NUM_BINS_KE);

    TVectorD UnfoldedReco = WienerSVD(
        Response,
        Signal,
        Measure,
        Covariance,
        0, // 0: unit, 2: second derivative
        0, // smoothing parameter
        AddSmear,
        WF,
        UnfoldCov,
        CovRotation,
        AddSmearInverse
    );
    TVectorD SmearedTrue  = AddSmear * Signal;
    TVectorD UnSmearedUnf = AddSmearInverse * UnfoldedReco;

    // Organize results
    std::vector<TH1*> UnfoldedRecoHistos;
    std::vector<TH1*> UnSmearedUnfRecoHistos;
    std::vector<TH1*> SmearedTrueHistos;
    for (int iBin = 0; iBin < NUM_SIGNAL_TYPES; ++iBin) {
        TH1D* unfoldedHist          = new TH1D(Form("hUnfoldedRecoVec_%d", iBin), Form("Unfolded Reco Vec %d", iBin), NUM_BINS_KE, ARRAY_KE_BINS.data());
        TH1D* unSmearedUnfoldedHist = new TH1D(Form("hUnSmearedUnfoldedRecoVec_%d", iBin), Form("Unsmeared Unfolded Reco Vec %d", iBin), NUM_BINS_KE, ARRAY_KE_BINS.data());
        TH1D* smearedHist           = new TH1D(Form("hSmearedTrueVec_%d", iBin), Form("Smeared True Vec %d", iBin), NUM_BINS_KE, ARRAY_KE_BINS.data());
        UnfoldedRecoHistos.push_back(unfoldedHist);
        UnSmearedUnfRecoHistos.push_back(unSmearedUnfoldedHist);
        SmearedTrueHistos.push_back(smearedHist);
    }

    for (int iFlatIndex = 0; iFlatIndex < NUM_SIGNAL_TYPES * NUM_BINS_KE; ++iFlatIndex) {
        auto [signalBin, energyBin] = unflattenIndex(iFlatIndex, NUM_BINS_KE);

        // Get error for unfolded reco histogram from unfolded covariance
        double err = std::sqrt(UnfoldCov[iFlatIndex][iFlatIndex]);
        UnfoldedRecoHistos[signalBin]->SetBinContent(energyBin + 1, UnfoldedReco(iFlatIndex));
        UnfoldedRecoHistos[signalBin]->SetBinError(energyBin + 1, err);

        // TODO: get error for un-smeared histos
        UnSmearedUnfRecoHistos[signalBin]->SetBinContent(energyBin + 1, UnSmearedUnf(iFlatIndex));

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
    TH2D* hSmearingInv = new TH2D("hSmearingInverse", "Smearing Inverse Matrix;True (i, #alpha); Reco (j, #beta)", 
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE, 
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    );
    M2H(Covariance, hCovariance); M2H(AddSmear, hSmearing); M2H(AddSmearInverse, hSmearingInv);

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

    // Set errors to corrected incident KE histogram (TODO: errors on corrections)
    for (int iBin = 1; iBin <= hIncidentKECorrected->GetNbinsX(); ++iBin) {
        double entries = hIncidentKECorrected->GetBinContent(iBin);
        hIncidentKECorrected->SetBinError(iBin, std::sqrt(entries));
    }

    ///////////////////////////////////////////////
    // Get cross-section using corrrected fluxes //
    ///////////////////////////////////////////////

    TH1D* hPionAbs0pCrossSection            = new TH1D("hPionAbs0pCrossSection", "hPionAbs0pCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hUnSmearedPionAbs0pCrossSection   = new TH1D("hUnSmearedPionAbs0pCrossSection", "hUnSmearedPionAbs0pCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hSmearedTruePionAbs0pCrossSection = new TH1D("hSmearedTruePionAbs0pCrossSection", "hSmearedTruePionAbs0pCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTruePionAbs0pCrossSection        = new TH1D("hTruePionAbs0pCrossSection", "hTruePionAbs0pCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    TH1D* hPionAbsNpCrossSection            = new TH1D("hPionAbsNpCrossSection", "hPionAbsNpCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hUnSmearedPionAbsNpCrossSection   = new TH1D("hUnSmearedPionAbsNpCrossSection", "hUnSmearedPionAbsNpCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hSmearedTruePionAbsNpCrossSection = new TH1D("hSmearedTruePionAbsNpCrossSection", "hSmearedTruePionAbsNpCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTruePionAbsNpCrossSection        = new TH1D("hTruePionAbsNpCrossSection", "hTruePionAbsNpCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    TH1D* hPionScatterCrossSection            = new TH1D("hPionScatterCrossSection", "hPionScatterCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hUnSmearedPionScatterCrossSection   = new TH1D("hUnSmearedPionScatterCrossSection", "hUnSmearedPionScatterCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hSmearedTruePionScatterCrossSection = new TH1D("hSmearedTruePionScatterCrossSection", "hSmearedTruePionScatterCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTruePionScatterCrossSection        = new TH1D("hTruePionScatterCrossSection", "hTruePionScatterCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    TH1D* hPionChExchCrossSection            = new TH1D("hPionChExchCrossSection", "hPionChExchCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hUnSmearedPionChExchCrossSection   = new TH1D("hUnSmearedPionChExchCrossSection", "hUnSmearedPionChExchCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hSmearedTruePionChExchCrossSection = new TH1D("hSmearedTruePionChExchCrossSection", "hSmearedTruePionChExchCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTruePionChExchCrossSection        = new TH1D("hTruePionChExchCrossSection", "hTruePionChExchCrossSection;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    std::vector<TH1D*> UnfoldedCrossSections = {
        hPionAbs0pCrossSection,
        hPionAbsNpCrossSection,
        hPionScatterCrossSection,
        hPionChExchCrossSection
    };

    std::vector<TH1D*> UnSmearedUnfoldedCrossSections = {
        hUnSmearedPionAbs0pCrossSection,
        hUnSmearedPionAbsNpCrossSection,
        hUnSmearedPionScatterCrossSection,
        hUnSmearedPionChExchCrossSection
    };

    std::vector<TH1D*> SmearedTrueCrossSections = {
        hSmearedTruePionAbs0pCrossSection,
        hSmearedTruePionAbsNpCrossSection,
        hSmearedTruePionScatterCrossSection,
        hSmearedTruePionChExchCrossSection
    };

    std::vector<TH1D*> TrueCrossSections = {
        hTruePionAbs0pCrossSection,
        hTruePionAbsNpCrossSection,
        hTruePionScatterCrossSection,
        hTruePionChExchCrossSection
    };

    for (int i = 0; i < UnfoldedCrossSections.size(); ++i) {
        TH1D* unfXSec          = UnfoldedCrossSections[i];
        TH1D* unSmearedUnfXSec = UnSmearedUnfoldedCrossSections[i];
        TH1D* smearedTrueXSec  = SmearedTrueCrossSections[i];
        TH1D* trueXSec         = TrueCrossSections[i];

        for (int iBin = 1; iBin <= NUM_BINS_KE; ++iBin) {
            double incidentErr     = hIncidentKE->GetBinError(iBin);
            double incidentContent = hIncidentKE->GetBinContent(iBin);

            double interactingErr     = UnfoldedRecoHistos[i]->GetBinError(iBin);
            double interactingContent = UnfoldedRecoHistos[i]->GetBinContent(iBin);

            double unSmearedInteractingContent = UnSmearedUnfRecoHistos[i]->GetBinContent(iBin);

            // we compute error as δσ / σ = δN_int / N_int + δN_inc / N_inc (conservative upper bound)
            double crossSection    = interactingContent / incidentContent;
            double crossSectionErr = crossSection * ((interactingErr / interactingContent) + (incidentErr / incidentContent));
            unfXSec->SetBinContent(iBin, crossSection);
            unfXSec->SetBinError(iBin, crossSectionErr);

            // TODO: error
            double unsmearedCrossSection = unSmearedInteractingContent / incidentContent;
            unSmearedUnfXSec->SetBinContent(iBin, unsmearedCrossSection);

            // True cross-sections, no error bars
            smearedTrueXSec->SetBinContent(iBin, SmearedTrueHistos[i]->GetBinContent(iBin) / hTrueIncidentKE->GetBinContent(iBin));
            trueXSec->SetBinContent(iBin, TotalEventsHistos[i]->GetBinContent(iBin) / hTrueIncidentKE->GetBinContent(iBin));
        }

        // Incident KE corrections
        unfXSec->Multiply(hPsiInc);
        unfXSec->Divide(hCInc);

        unSmearedUnfXSec->Multiply(hPsiInc);
        unSmearedUnfXSec->Divide(hCInc);

        // Scale units
        unfXSec->Scale(1.0 / (number_density * slab_width * 1e-28));
        unSmearedUnfXSec->Scale(1.0 / (number_density * slab_width * 1e-28));
        smearedTrueXSec->Scale(1.0 / (number_density * slab_width * 1e-28));
        trueXSec->Scale(1.0 / (number_density * slab_width * 1e-28));

        // Make contents per 50 MeV
        reweightOneDHisto(unfXSec, 50.);
        reweightOneDHisto(unSmearedUnfXSec, 50.);
        reweightOneDHisto(smearedTrueXSec, 50.);
        reweightOneDHisto(trueXSec, 50.);
    }

    ///////////////////////////
    // Make plots per 50 MeV //
    ///////////////////////////

    // TODO: create function that scales histogram bins per 50 MeV

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
        {hPionAbs0pKETrue, hPionAbs0pKEAbsNp, hPionAbs0pKEScatter, hPionAbs0pKEChExch, hPionAbs0pKEMuonTG, hPionAbs0pKEMuonDecay, hPionAbs0pKEMuonCAR, hPionAbs0pKEElectron, hPionAbs0pKEOther},
        {hPionAbsNpKETrue, hPionAbsNpKEAbs0p, hPionAbsNpKEScatter, hPionAbsNpKEChExch, hPionAbsNpKEMuon, hPionAbsNpKEElectron, hPionAbsNpKEOther},
        {hPionScatterKETrue, hPionScatterKEAbs0p, hPionScatterKEAbsNp, hPionScatterKEChExch, hPionScatterKEMuon, hPionScatterKEElectron, hPionScatterKEOther},
        {hPionScatterKETrue, hPionScatterKEAbs0p, hPionScatterKEAbsNp, hPionScatterKEChExch, hPionScatterKEMuonTG, hPionScatterKEMuonDecay, hPionScatterKEMuonCAR, hPionScatterKEElectron, hPionScatterKEOther},
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

        // True cross-section comparison (reg vs true space)
        {hTruePionAbs0pCrossSection, hSmearedTruePionAbs0pCrossSection},
        {hTruePionAbsNpCrossSection, hSmearedTruePionAbsNpCrossSection},
        {hTruePionScatterCrossSection, hSmearedTruePionScatterCrossSection},
        {hTruePionChExchCrossSection, hSmearedTruePionChExchCrossSection},

        // Cross-sections (reg space)
        {hSmearedTruePionAbs0pCrossSection, hPionAbs0pCrossSection},
        {hSmearedTruePionAbsNpCrossSection, hPionAbsNpCrossSection},
        {hSmearedTruePionScatterCrossSection, hPionScatterCrossSection},
        {hSmearedTruePionChExchCrossSection, hPionChExchCrossSection},

        // Cross-sections (true space)
        {hTruePionAbs0pCrossSection, hUnSmearedPionAbs0pCrossSection},
        {hTruePionAbsNpCrossSection, hUnSmearedPionAbsNpCrossSection},
        {hTruePionScatterCrossSection, hUnSmearedPionScatterCrossSection},
        {hTruePionChExchCrossSection, hUnSmearedPionChExchCrossSection},

        // BDT probabilities
        {hBestAbs0pAbs0p, hBestAbs0pAbsNp, hBestAbs0pMuon, hBestAbs0pElectron, hBestAbs0pScatter, hBestAbs0pChExch, hBestAbs0pOther},
        {hBestChExchAbs0p, hBestChExchAbsNp, hBestChExchMuon, hBestChExchElectron, hBestChExchScatter, hBestChExchChExch, hBestChExchOther},
        {hBestScatterAbs0p, hBestScatterAbsNp, hBestScatterMuon, hBestScatterElectron, hBestScatterScatter, hBestScatterChExch, hBestScatterOther},
        {hBestElectronAbs0p, hBestElectronAbsNp, hBestElectronMuon, hBestElectronElectron, hBestElectronScatter, hBestElectronChExch, hBestElectronOther},

        // Unreconstructed hits
        {hHitClusterSizesAbs0p, hHitClusterSizesAbsNp, hHitClusterSizesMuon, hHitClusterSizesElectron, hHitClusterSizesScatter, hHitClusterSizesChExch, hHitClusterSizesOther},
        {hLargeHitClusterAbs0p, hLargeHitClusterAbsNp, hLargeHitClusterMuon, hLargeHitClusterElectron, hLargeHitClusterScatter, hLargeHitClusterChExch, hLargeHitClusterOther}
    };

    std::vector<std::vector<TString>> PlotLabelGroups = {
        // Incident KE
        {"Pions", "Muons", "Electrons"},
        {"Corrected", "True"},

        // Interacting KE
        {"True", "Abs Np", "Scatter", "Ch. exch.", "Muon", "Electron", "Other"},
        {"True", "Abs Np", "Scatter", "Ch. exch.", "Muon TG", "Muon decay", "Muon CAR", "Electron", "Other"},
        {"True", "Abs 0p", "Scatter", "Ch. exch.", "Muon", "Electron", "Other"},
        {"True", "Abs 0p", "Abs Np", "Ch. exch.", "Muon", "Electron", "Other"},
        {"True", "Abs 0p", "Abs Np", "Ch. exch.", "Muon TG", "Muon decay", "Muon CAR", "Electron", "Other"},
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

        // True cross-section comparison (reg vs true space)
        {"True (t)", "True (r)"},
        {"True (t)", "True (r)"},
        {"True (t)", "True (r)"},
        {"True (t)", "True (r)"},

        // Cross-sections (reg space)
        {"True (r)", "Unf. (r)"},
        {"True (r)", "Unf. (r)"},
        {"True (r)", "Unf. (r)"},
        {"True (r)", "Unf. (r)"},

        // Cross-sections (true space)
        {"True (t)", "Unf. (t)"},
        {"True (t)", "Unf. (t)"},
        {"True (t)", "Unf. (t)"},
        {"True (t)", "Unf. (t)"},

        // BDT probabilities
        {"Abs 0p", "Abs Np", "Muon", "Electron", "Scatter", "Ch. exch.", "Other"},
        {"Abs 0p", "Abs Np", "Muon", "Electron", "Scatter", "Ch. exch.", "Other"},
        {"Abs 0p", "Abs Np", "Muon", "Electron", "Scatter", "Ch. exch.", "Other"},
        {"Abs 0p", "Abs Np", "Muon", "Electron", "Scatter", "Ch. exch.", "Other"},

        // Unreconstructed hits
        {"Abs 0p", "Abs Np", "Muon", "Electron", "Scatter", "Ch. exch.", "Other"},
        {"Abs 0p", "Abs Np", "Muon", "Electron", "Scatter", "Ch. exch.", "Other"}
    };

    std::vector<TString> PlotTitles = {
        // Incident KE
        "Incident/IncidentKE",
        "Incident/IncidentKECorrected",

        // Interacting KE
        "RecoInteracting/Abs0pInteractingKE",
        "RecoInteracting/Abs0pInteractingKEDetailed",
        "RecoInteracting/AbsNpInteractingKE",
        "RecoInteracting/ScatterInteractingKE",
        "RecoInteracting/ScatterInteractingKEDetailed",
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

        // True cross-section comparison (reg vs true space)
        "CrossSection/ComparePionAbs0pCrossSection",
        "CrossSection/ComparePionAbsNpCrossSection",
        "CrossSection/ComparePionScatterCrossSection",
        "CrossSection/ComparePionChExchCrossSection",

        // Cross-sections (reg space)
        "CrossSection/PionAbs0pCrossSection",
        "CrossSection/PionAbsNpCrossSection",
        "CrossSection/PionScatterCrossSection",
        "CrossSection/PionChExchCrossSection",

        // Cross-sections (true space)
        "CrossSection/UnSmearedPionAbs0pCrossSection",
        "CrossSection/UnSmearedPionAbsNpCrossSection",
        "CrossSection/UnSmearedPionScatterCrossSection",
        "CrossSection/UnSmearedPionChExchCrossSection",

        // BDT probabilities
        "BDTScores/Abs0pTrue",
        "BDTScores/ChExchTrue",
        "BDTScores/ScatterTrue",
        "BDTScores/ElectronTrue",

        // Unreconstructed hits
        "Hits/ClusterSizes",
        "Hits/NumLargeClusters"
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

        // True cross-section comparison (reg vs true space)
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",

        // Cross-sections (reg space)
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",

        // Cross-sections (true space)
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",

        // BDT probabilities
        "Probability",
        "Probability",
        "Probability",
        "Probability",

        // Unreconstructed hits
        "Cluster size [cm]",
        "Number of large clusters"
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

        // True cross-section comparison (reg vs true space)
        "Cross section [barn] per 50 MeV",
        "Cross section [barn] per 50 MeV",
        "Cross section [barn] per 50 MeV",
        "Cross section [barn] per 50 MeV",

        // Cross-sections (reg space)
        "Cross section [barn] per 50 MeV",
        "Cross section [barn] per 50 MeV",
        "Cross section [barn] per 50 MeV",
        "Cross section [barn] per 50 MeV",

        // Cross-sections (true space)
        "Cross section [barn] per 50 MeV",
        "Cross section [barn] per 50 MeV",
        "Cross section [barn] per 50 MeV",
        "Cross section [barn] per 50 MeV",

        // BDT probabilities
        "Counts",
        "Counts",
        "Counts",
        "Counts",

        // Unreconstructed hits
        "Counts",
        "Counts"
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

        // True cross-section comparison (reg vs true space)
        false,
        false,
        false,
        false,

        // Cross-sections (reg space)
        false,
        false,
        false,
        false,

        // Cross-sections (true space)
        false,
        false,
        false,
        false,

        // BDT probabilities
        true,
        true,
        true,
        true,

        // Unreconstructed hits
        true,
        true
    };

    std::vector<std::vector<bool>> PlotsAsPoints = {
        // Incident KE
        {false, false, false},
        {true, false},

        // Interacting KE
        {false, false, false, false, false, false, false},
        {false, false, false, false, false, false, false, false, false},
        {false, false, false, false, false, false, false},
        {false, false, false, false, false, false, false},
        {false, false, false, false, false, false, false, false, false},
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

        // True cross-section comparison (reg vs true space)
        {false, false},
        {false, false},
        {false, false},
        {false, false},

        // Cross-sections (reg space)
        {false, true},
        {false, true},
        {false, true},
        {false, true},

        // Cross-sections (true space)
        {false, true},
        {false, true},
        {false, true},
        {false, true},

        // BDT probabilities
        {false, false, false, false, false, false, false},
        {false, false, false, false, false, false, false},
        {false, false, false, false, false, false, false},
        {false, false, false, false, false, false, false},

        // Unreconstructed hits
        {false, false, false, false, false, false, false},
        {false, false, false, false, false, false, false}
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
        PlotLabelGroups.push_back({"True (r)", "Unf. (r)"});
        PlotTitles.push_back("Unfolded/" + UnfHistTitles[i]);
        XLabels.push_back("Kinetic energy [MeV]");
        YLabels.push_back("Counts");
        PlotStacked.push_back(false);

        PlotGroups.push_back({TotalEventsHistos[i], UnSmearedUnfRecoHistos[i]});
        PlotsAsPoints.push_back({false, true});
        PlotLabelGroups.push_back({"True (t)", "Unf. (t)"});
        PlotTitles.push_back("Unfolded/UnSmeared" + UnfHistTitles[i] );
        XLabels.push_back("Kinetic energy [MeV]");
        YLabels.push_back("Counts");
        PlotStacked.push_back(false);

        PlotGroups.push_back({TotalEventsHistos[i], SmearedTrueHistos[i]});
        PlotsAsPoints.push_back({false, false});
        PlotLabelGroups.push_back({"True (t)", "True (r)"});
        TString modifiedTitle = UnfHistTitles[i];
        modifiedTitle.ReplaceAll("Unfolded", "");
        PlotTitles.push_back("Unfolded/TrueComparison" + modifiedTitle);
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
        hSmearing,
        hSmearingInv
    };

    std::vector<TString> TwoDTitles = {
        "Response/ResponseMatrix",
        "Covariance/UnfoldedCovariance",
        "Smearing/AddSmearing",
        "Smearing/AddSmearingInverse"
    };

    std::vector<std::pair<double,double>> TwoDRanges = {
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0}
    };

    std::vector<bool> TwoDDisplayNumbers = {
        false,
        false,
        false,
        false
    };

    printTwoDPlots(SaveDir, TwoDPlots, TwoDTitles, TwoDRanges, TwoDDisplayNumbers);
}