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

void RecoAllAnalysis() {
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
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis/figs/RecoAllAnalysis/";

    // Load root file
    TString RootFilePath = "/exp/lariat/app/users/epelaez/files/RecoNNAllEval_histo2.root"; 
    std::unique_ptr<TFile> File(TFile::Open(RootFilePath));
    TDirectory* Directory = (TDirectory*)File->Get("RecoNNAllEval");

    // Load root file with true cross sections
    TString TrueXSRootFilePath = "/exp/lariat/app/users/epelaez/files/TrueXSPionAbs_histo.root";
    std::unique_ptr<TFile> TrueXSFile(TFile::Open(TrueXSRootFilePath));
    TDirectory* TrueXSDirectory = (TDirectory*) TrueXSFile->Get("TrueXSPionAbs");

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

    // Cut information
    bool passesPionInRedVolume, passesNoOutgoingPion, passesSmallTracksCut, passesMeanCurvatureCut;
    tree->SetBranchAddress("passesPionInRedVolume", &passesPionInRedVolume);
    tree->SetBranchAddress("passesNoOutgoingPion", &passesNoOutgoingPion);
    tree->SetBranchAddress("passesSmallTracksCut", &passesSmallTracksCut);
    tree->SetBranchAddress("passesMeanCurvatureCut", &passesMeanCurvatureCut);

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

    // Truth-level charge exchange information
    std::vector<int>* chExchShowerIDs = nullptr;
    std::vector<std::string>* chExchShowerProcesses = nullptr;
    std::vector<int>* chExchShowerPDGs = nullptr;
    std::vector<double>* chExchShowerLengths = nullptr;
    std::vector<std::vector<double>>* chExchShowerStart = nullptr;
    std::vector<std::vector<double>>* chExchShowerEnd = nullptr;
    std::vector<int>* chExchShowerNeutralPionDaughtersID = nullptr;
    std::vector<std::vector<double>>* chExchShowerNeutralPionDaughtersMom = nullptr;
    tree->SetBranchAddress("chExchShowerIDs", &chExchShowerIDs);
    tree->SetBranchAddress("chExchShowerProcesses", &chExchShowerProcesses);
    tree->SetBranchAddress("chExchShowerPDGs", &chExchShowerPDGs);
    tree->SetBranchAddress("chExchShowerLengths", &chExchShowerLengths);
    tree->SetBranchAddress("chExchShowerStart", &chExchShowerStart);
    tree->SetBranchAddress("chExchShowerEnd", &chExchShowerEnd);
    tree->SetBranchAddress("chExchShowerNeutralPionDaughtersID", &chExchShowerNeutralPionDaughtersID);
    tree->SetBranchAddress("chExchShowerNeutralPionDaughtersMom", &chExchShowerNeutralPionDaughtersMom);

    // Truth-level beamline electron shower information
    std::vector<int>* electronShowerIDs = nullptr;
    std::vector<std::string>* electronShowerProcesses = nullptr;
    std::vector<int>* electronShowerPDGs = nullptr;
    std::vector<double>* electronShowerLengths = nullptr;
    std::vector<std::vector<double>>* electronShowerStart = nullptr;
    std::vector<std::vector<double>>* electronShowerEnd = nullptr;
    tree->SetBranchAddress("electronShowerIDs", &electronShowerIDs);
    tree->SetBranchAddress("electronShowerProcesses", &electronShowerProcesses);
    tree->SetBranchAddress("electronShowerPDGs", &electronShowerPDGs);
    tree->SetBranchAddress("electronShowerLengths", &electronShowerLengths);
    tree->SetBranchAddress("electronShowerStart", &electronShowerStart);
    tree->SetBranchAddress("electronShowerEnd", &electronShowerEnd);

    // Primaries information
    int                  numPrimaries;
    int                  numValidPrimaries;
    std::vector<double>* primariesStartX = nullptr;
    std::vector<double>* primariesStartY = nullptr;
    std::vector<double>* primariesStartZ = nullptr;
    std::vector<double>* primariesEndX = nullptr;
    std::vector<double>* primariesEndY = nullptr;
    std::vector<double>* primariesEndZ = nullptr;
    std::vector<int>*    primariesPDG = nullptr;
    std::vector<int>*    primariesID = nullptr;
    tree->SetBranchAddress("numPrimaries", &numPrimaries);
    tree->SetBranchAddress("numValidPrimaries", &numValidPrimaries);
    tree->SetBranchAddress("primariesStartX", &primariesStartX);
    tree->SetBranchAddress("primariesStartY", &primariesStartY);
    tree->SetBranchAddress("primariesStartZ", &primariesStartZ);
    tree->SetBranchAddress("primariesEndX", &primariesEndX);
    tree->SetBranchAddress("primariesEndY", &primariesEndY);
    tree->SetBranchAddress("primariesEndZ", &primariesEndZ);
    tree->SetBranchAddress("primariesPDG", &primariesPDG);
    tree->SetBranchAddress("primariesID", &primariesID);

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

    // Information about post-scattering interactions
    std::vector<int>*                 secondaryInteractionTypes = nullptr;
    std::vector<int>*                 secondaryInteractionTrkID = nullptr;
    std::vector<double>*              secondaryInteractionInteractingKE = nullptr;
    std::vector<double>*              secondaryInteractionAngle = nullptr;
    std::vector<double>*              secondaryInteractionXPosition = nullptr;
    std::vector<double>*              secondaryInteractionYPosition = nullptr;
    std::vector<double>*              secondaryInteractionZPosition = nullptr;
    std::vector<std::vector<double>>* secondaryIncidentKEContributions = nullptr;
    std::vector<std::vector<int>>* secondaryInteractionDaughtersPDG = nullptr;
    std::vector<std::vector<double>>* secondaryInteractionDaughtersKE = nullptr;
    // std::vector<std::vector<std::string>>* secondaryInteractionDaughtersProcess = nullptr;
    tree->SetBranchAddress("secondaryInteractionTypes", &secondaryInteractionTypes);
    tree->SetBranchAddress("secondaryInteractionTrkID", &secondaryInteractionTrkID);
    tree->SetBranchAddress("secondaryInteractionInteractingKE", &secondaryInteractionInteractingKE);
    tree->SetBranchAddress("secondaryInteractionAngle", &secondaryInteractionAngle);
    tree->SetBranchAddress("secondaryInteractionXPosition", &secondaryInteractionXPosition);
    tree->SetBranchAddress("secondaryInteractionYPosition", &secondaryInteractionYPosition);
    tree->SetBranchAddress("secondaryInteractionZPosition", &secondaryInteractionZPosition);
    tree->SetBranchAddress("secondaryIncidentKEContributions", &secondaryIncidentKEContributions);
    tree->SetBranchAddress("secondaryInteractionDaughtersPDG", &secondaryInteractionDaughtersPDG);
    tree->SetBranchAddress("secondaryInteractionDaughtersKE", &secondaryInteractionDaughtersKE);
    // tree->SetBranchAddress("secondaryInteractionDaughtersProcess", &secondaryInteractionDaughtersProcess);

    ///////////////////////
    // Create histograms //
    ///////////////////////

    TH1D* hTotalEvents = new TH1D("hTotalEvents", "hTotalEvents", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);

    TH1D *hRecoAbsorption   = new TH1D("hRecoAbsorption", "hRecoAbsorption;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D *hRecoAbsorption0p = new TH1D("hRecoAbsorption0p", "hRecoAbsorption0p;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D *hRecoAbsorptionNp = new TH1D("hRecoAbsorptionNp", "hRecoAbsorptionNp;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);

    //////////////////////////////////////////////
    // Histograms for primary stitching studies //
    //////////////////////////////////////////////

    TH1D *hStitchedDistanceFromVertex         = new TH1D("hStitchedDistanceFromVertex", "StitchedDistanceFromVertex;;", 20, 0, 70);
    TH1D *hStitchedOriginalDistanceFromVertex = new TH1D("hStitchedOriginalDistanceFromVertex", "hStitchedOriginalDistanceFromVertex;;", 20, 0, 70);

    TH1D *hStitchAsPionAndProton = new TH1D("hStitchAsPionAndProton", "hStitchAsPionAndProton;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D *hStitchAsPionBraggPeak = new TH1D("hStitchAsPionBraggPeak", "hStitchAsPionBraggPeak;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D *hStitchAsMIP  = new TH1D("hStitchAsMIP", "hStitchAsMIP;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D *hStitchAsProton        = new TH1D("hStitchAsProton", "hStitchAsProton;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D *hStitchFracBreakPoints = new TH1D("hStitchFracBreakPoints", "hStitchFracBreakPoints;;", 50, 0, 1);

    TH1D *hInelScatteringBackgroundOriginalDistanceToPrimaryVertex   = new TH1D("hInelScatteringBackgroundOriginalDistanceToPrimaryVertex", "hInelScatteringBackgroundOriginalDistanceToPrimaryVertex;;", 20, 0, 20);
    TH1D *hInelScatteringBackgroundOriginalDistanceToSecondaryVertex = new TH1D("hInelScatteringBackgroundOriginalDistanceToSecondaryVertex", "hInelScatteringBackgroundOriginalDistanceToSecondaryVertex;;", 20, 0, 20);
    TH1D *hInelScatteringBackgroundDistanceToPrimaryVertex           = new TH1D("hInelScatteringBackgroundDistanceToPrimaryVertex", "hInelScatteringBackgroundDistanceToPrimaryVertex;;", 20, 0, 20);
    TH1D *hInelScatteringBackgroundDistanceToSecondaryVertex         = new TH1D("hInelScatteringBackgroundDistanceToSecondaryVertex", "hInelScatteringBackgroundDistanceToSecondaryVertex;;", 20, 0, 20);

    /////////////////////////////////////////////////////
    // Histograms for secondary track PID optimization //
    /////////////////////////////////////////////////////

    TH1D *hSecondaryPionChi2Protons = new TH1D("hSecondaryPionChi2Protons", "hSecondaryPionChi2Protons;;", 20, 0, 10);
    TH1D *hSecondaryPionChi2Pions   = new TH1D("hSecondaryPionChi2Pions", "hSecondaryPionChi2Pions;;", 20, 0, 10);
    TH1D *hSecondaryPionChi2Others  = new TH1D("hSecondaryPionChi2Others", "hSecondaryPionChi2Others;;", 20, 0, 10);

    TH1D *hSecondaryProtonChi2Protons = new TH1D("hSecondaryProtonChi2Protons", "hSecondaryProtonChi2Protons;;", 20, 0, 10);
    TH1D *hSecondaryProtonChi2Pions   = new TH1D("hSecondaryProtonChi2Pions", "hSecondaryProtonChi2Pions;;", 20, 0, 10);
    TH1D *hSecondaryProtonChi2Others  = new TH1D("hSecondaryProtonChi2Others", "hSecondaryProtonChi2Others;;", 20, 0, 10);

    TH1D *hSecondaryMeanDEDXProtons = new TH1D("hSecondaryMeanDEDXProtons", "hSecondaryMeanDEDXProtons;;", 10, 0, 10);
    TH1D *hSecondaryMeanDEDXPions   = new TH1D("hSecondaryMeanDEDXPions", "hSecondaryMeanDEDXPions;;", 10, 0, 10);
    TH1D *hSecondaryMeanDEDXOthers  = new TH1D("hSecondaryMeanDEDXOthers", "hSecondaryMeanDEDXOthers;;", 10, 0, 10);

    //////////////////////////////////////////////////
    // Reconstruction power for secondary particles //
    //////////////////////////////////////////////////

    TH1D* hAllSecondaryPions  = new TH1D("hAllSecondaryPions", "hAllSecondaryPions;;", 20, 0, 0.5);
    TH1D* hRecoSecondaryPions = new TH1D("hRecoSecondaryPions", "hRecoSecondaryPions;;", 20, 0, 0.5);

    TH1D* hAllSecondaryProtons  = new TH1D("hAllSecondaryProtons", "hAllSecondaryProtons;;", 20, 0, 0.4);
    TH1D* hRecoSecondaryProtons = new TH1D("hRecoSecondaryProtons", "hRecoSecondaryProtons", 20, 0, 0.4); 

    TH2D* hScatteringVertexKEVsOutKE = new TH2D(
        "hScatteringVertexKEVsOutKE", 
        "hScatteringVertexKEVsOutKE;Vertex KE [GeV];Outgoing KE [GeV]", 
        20, 0, 0.6,
        20, 0, 0.6
    );

    /////////////////////////////////////////////////////////////////////////////
    // Variables and histograms for local linearity derivative cut optimizaton //
    /////////////////////////////////////////////////////////////////////////////

    double INIT_LIN_DERIV_THRESHOLD = 0.06e-3;
    double LIN_DERIV_STEP_SIZE      = 0.02e-3;
    int    LIN_DERIV_NUM_STEPS      = 10;

    TH1D *hLocalLinearityDerivativeReco     = new TH1D("hLocalLinearityDerivativeReco", "hLocalLinearityDerivativeReco;;", LIN_DERIV_NUM_STEPS, INIT_LIN_DERIV_THRESHOLD, INIT_LIN_DERIV_THRESHOLD + LIN_DERIV_NUM_STEPS * LIN_DERIV_STEP_SIZE);
    TH1D *hLocalLinearityDerivativeRecoTrue = new TH1D("hLocalLinearityDerivativeRecoTrue", "hLocalLinearityDerivativeRecoTrue;;", LIN_DERIV_NUM_STEPS, INIT_LIN_DERIV_THRESHOLD, INIT_LIN_DERIV_THRESHOLD + LIN_DERIV_NUM_STEPS * LIN_DERIV_STEP_SIZE);
    TH1D *hLocalLinearityDerivativeCutPur   = new TH1D("hLocalLinearityDerivativeCutPur", "hLocalLinearityDerivativeCutPur;;", LIN_DERIV_NUM_STEPS, INIT_LIN_DERIV_THRESHOLD, INIT_LIN_DERIV_THRESHOLD + LIN_DERIV_NUM_STEPS * LIN_DERIV_STEP_SIZE);
    TH1D *hLocalLinearityDerivativeCutEff   = new TH1D("hLocalLinearityDerivativeCutEff", "hLocalLinearityDerivativeCutEff;;", LIN_DERIV_NUM_STEPS, INIT_LIN_DERIV_THRESHOLD, INIT_LIN_DERIV_THRESHOLD + LIN_DERIV_NUM_STEPS * LIN_DERIV_STEP_SIZE);

    double INIT_LIN_STDDEV_THRESHOLD = 0.01e-3;
    double LIN_STDDEV_STEP_SIZE      = 0.005e-3;
    int    LIN_STDDEV_NUM_STEPS      = 20;

    TH1D *hLocalLinearityStdDevReco     = new TH1D("hLocalLinearityStdDevReco", "hLocalLinearityStdDevReco;;", LIN_STDDEV_NUM_STEPS, INIT_LIN_STDDEV_THRESHOLD, INIT_LIN_STDDEV_THRESHOLD + LIN_STDDEV_NUM_STEPS * LIN_STDDEV_STEP_SIZE);
    TH1D *hLocalLinearityStdDevRecoTrue = new TH1D("hLocalLinearityStdDevRecoTrue", "hLocalLinearityStdDevRecoTrue;;", LIN_STDDEV_NUM_STEPS, INIT_LIN_STDDEV_THRESHOLD, INIT_LIN_STDDEV_THRESHOLD + LIN_STDDEV_NUM_STEPS * LIN_STDDEV_STEP_SIZE);
    TH1D *hLocalLinearityStdDevCutPur   = new TH1D("hLocalLinearityStdDevCutPur", "hLocalLinearityStdDevCutPur;;", LIN_STDDEV_NUM_STEPS, INIT_LIN_STDDEV_THRESHOLD, INIT_LIN_STDDEV_THRESHOLD + LIN_STDDEV_NUM_STEPS * LIN_STDDEV_STEP_SIZE);
    TH1D *hLocalLinearityStdDevCutEff   = new TH1D("hLocalLinearityStdDevCutEff", "hLocalLinearityStdDevCutEff;;", LIN_STDDEV_NUM_STEPS, INIT_LIN_STDDEV_THRESHOLD, INIT_LIN_STDDEV_THRESHOLD + LIN_STDDEV_NUM_STEPS * LIN_STDDEV_STEP_SIZE);

    double INIT_LIN_MIN_THRESHOLD = 0.9995;
    double LIN_MIN_STEP_SIZE      = -0.0002;
    int    LIN_MIN_NUM_STEPS      = 20;

    TH1D *hLocalLinearityMinReco     = new TH1D("hLocalLinearityMinReco", "hLocalLinearityMinReco;;", LIN_MIN_NUM_STEPS, 1 - INIT_LIN_MIN_THRESHOLD, 1 - (INIT_LIN_MIN_THRESHOLD + LIN_MIN_NUM_STEPS * LIN_MIN_STEP_SIZE));
    TH1D *hLocalLinearityMinRecoTrue = new TH1D("hLocalLinearityMinRecoTrue", "hLocalLinearityMinRecoTrue;;", LIN_MIN_NUM_STEPS, 1 - INIT_LIN_MIN_THRESHOLD, 1 - (INIT_LIN_MIN_THRESHOLD + LIN_MIN_NUM_STEPS * LIN_MIN_STEP_SIZE));
    TH1D *hLocalLinearityMinCutPur   = new TH1D("hLocalLinearityMinCutPur", "hLocalLinearityMinCutPur;;", LIN_MIN_NUM_STEPS, 1 - INIT_LIN_MIN_THRESHOLD, 1 - (INIT_LIN_MIN_THRESHOLD + LIN_MIN_NUM_STEPS * LIN_MIN_STEP_SIZE));
    TH1D *hLocalLinearityMinCutEff   = new TH1D("hLocalLinearityMinCutEff", "hLocalLinearityMinCutEff;;", LIN_MIN_NUM_STEPS, 1 - INIT_LIN_MIN_THRESHOLD, 1 - (INIT_LIN_MIN_THRESHOLD + LIN_MIN_NUM_STEPS * LIN_MIN_STEP_SIZE));

    TH1D *hMinimumLinearity0p           = new TH1D("hMinimumLinearity0p", "hMinimumLinearity0p;;", 20, 0.998, 1);
    TH1D *hMinimumLinearity0pBackground = new TH1D("hMinimumLinearity0pBackground", "hMinimumLinearity0pBackground;;", 20, 0.998, 1);
    TH1D *hMinimumLinearityNp           = new TH1D("hMinimumLinearityNp", "hMinimumLinearityNp;;", 20, 0.995, 1);
    TH1D *hMinimumLinearityNpBackground = new TH1D("hMinimumLinearityNpBackground", "hMinimumLinearityNpBackground;;", 20, 0.995, 1);

    TH1D *hStdDevLinearity0p           = new TH1D("hStdDevLinearity0p", "hStdDevLinearity0p;;", 20, 0, 0.0001);
    TH1D *hStdDevLinearity0pBackground = new TH1D("hStdDevLinearity0pBackground", "hStdDevLinearity0pBackground;;", 20, 0, 0.0001);
    TH1D *hStdDevLinearityNp           = new TH1D("hStdDevLinearityNp", "hStdDevLinearityNp;;", 20, 0, 0.0002);
    TH1D *hStdDevLinearityNpBackground = new TH1D("hStdDevLinearityNpBackground", "hStdDevLinearityNpBackground;;", 20, 0, 0.0002);

    TH1D *hMaxLinearityD0p           = new TH1D("hMaxLinearityD0p", "hMaxLinearityD0p;;", 20, 0, 0.0004);
    // TH1D *hMaxLinearityD0pScattering = new TH1D("hMaxLinearityD0pScattering", "hMaxLinearityD0pScattering;;", 20, 0, 0.001);
    TH1D *hMaxLinearityD0pBackground = new TH1D("hMaxLinearityD0pBackground", "hMaxLinearityD0pBackground;;", 20, 0, 0.0004);

    TH1D *hMaxLinearityDD0p           = new TH1D("hMaxLinearityDD0p", "hMaxLinearityDD0p;;", 20, 0, 0.0002);
    TH1D *hMaxLinearityDD0pBackground = new TH1D("hMaxLinearityDD0pBackground", "hMaxLinearityDD0pBackground;;", 20, 0, 0.0002);

    /////////////////////////////////////////////////////
    // Variables and histograms for chi^2 optimization //
    /////////////////////////////////////////////////////

    int protonChiNumSteps = 20;
    double protonChiStart = 0.;
    double protonChiEnd   = 5.;
    double protonChiStep  = (protonChiEnd - protonChiStart) / ((double) protonChiNumSteps);

    int pionChiNumSteps = 20;
    double pionChiStart = 0.;
    double pionChiEnd   = 5.;
    double pionChiStep  = (pionChiEnd - pionChiStart) / ((double) pionChiNumSteps);

    TH2D* hProtonChi2TruePositives = new TH2D(
        "hProtonChi2TruePositives", 
        "hProtonChi2TruePositives;Proton chi squared;Pion chi squared", 
        protonChiNumSteps, protonChiStart, protonChiEnd, 
        pionChiNumSteps, pionChiStart, pionChiEnd
    );
    TH2D* hProtonChi2FalsePositives = new TH2D(
        "hProtonChi2FalsePositives", 
        "hProtonChi2FalsePositives;Proton chi squared;Pion chi squared", 
        protonChiNumSteps, protonChiStart, protonChiEnd, 
        pionChiNumSteps, pionChiStart, pionChiEnd
    );
    TH2D* hProtonChi2TrueNegatives = new TH2D(
        "hProtonChi2TrueNegatives", 
        "hProtonChi2TrueNegatives;Proton chi squared;Pion chi squared", 
        protonChiNumSteps, protonChiStart, protonChiEnd, 
        pionChiNumSteps, pionChiStart, pionChiEnd
    );
    TH2D* hProtonChi2FalseNegatives = new TH2D(
        "hProtonChi2FalseNegatives", 
        "hProtonChi2FalseNegatives;Proton chi squared;Pion chi squared", 
        protonChiNumSteps, protonChiStart, protonChiEnd, 
        pionChiNumSteps, pionChiStart, pionChiEnd
    );

    TH2D* hPionChi2TruePositives = new TH2D(
        "hPionChi2TruePositives", 
        "hPionChi2TruePositives;Proton chi squared;Pion chi squared", 
        protonChiNumSteps, protonChiStart, protonChiEnd, 
        pionChiNumSteps, pionChiStart, pionChiEnd
    );
    TH2D* hPionChi2FalsePositives = new TH2D(
        "hPionChi2FalsePositives", 
        "hPionChi2FalsePositives;Proton chi squared;Pion chi squared", 
        protonChiNumSteps, protonChiStart, protonChiEnd, 
        pionChiNumSteps, pionChiStart, pionChiEnd
    );
    TH2D* hPionChi2TrueNegatives = new TH2D(
        "hPionChi2TrueNegatives", 
        "hPionChi2TrueNegatives;Proton chi squared;Pion chi squared", 
        protonChiNumSteps, protonChiStart, protonChiEnd, 
        pionChiNumSteps, pionChiStart, pionChiEnd
    );
    TH2D* hPionChi2FalseNegatives = new TH2D(
        "hPionChi2FalseNegatives", 
        "hPionChi2FalseNegatives;Proton chi squared;Pion chi squared", 
        protonChiNumSteps, protonChiStart, protonChiEnd, 
        pionChiNumSteps, pionChiStart, pionChiEnd
    );

    ///////////////////////////////////////////////////////////
    // Variables and histograms for hit cluster optimization //
    ///////////////////////////////////////////////////////////

    double STARTING_CLUSTER_SIZE  = 2. * HIT_WIRE_SEPARATION;
    double CLUSTER_SIZE_STEP      = 1. * HIT_WIRE_SEPARATION;
    int    CLUSTER_SIZE_NUM_STEPS = 5;
    std::vector<int> NUM_CLUSTERS;

    int STARTING_NUM_CLUSTERS  = 1;
    int NUM_CLUSTERS_STEP      = 1;
    int NUM_CLUSTERS_NUM_STEPS = 5;

    TH2D* hHitClusterCut0pReco = new TH2D(
        "hHitClusterCut0pReco", 
        "hHitClusterCut0pReco;Large cluster size;Number of large clusters", 
        CLUSTER_SIZE_NUM_STEPS, STARTING_CLUSTER_SIZE, STARTING_CLUSTER_SIZE + (CLUSTER_SIZE_STEP * CLUSTER_SIZE_NUM_STEPS), 
        NUM_CLUSTERS_NUM_STEPS, STARTING_NUM_CLUSTERS, STARTING_NUM_CLUSTERS + (NUM_CLUSTERS_STEP * NUM_CLUSTERS_NUM_STEPS)
    );
    TH2D* hHitClusterCut0pRecoTrue = new TH2D(
        "hHitClusterCut0pRecoTrue", 
        "hHitClusterCut0pRecoTrue;Large cluster size;Number of large clusters", 
        CLUSTER_SIZE_NUM_STEPS, STARTING_CLUSTER_SIZE, STARTING_CLUSTER_SIZE + (CLUSTER_SIZE_STEP * CLUSTER_SIZE_NUM_STEPS), 
        NUM_CLUSTERS_NUM_STEPS, STARTING_NUM_CLUSTERS, STARTING_NUM_CLUSTERS + (NUM_CLUSTERS_STEP * NUM_CLUSTERS_NUM_STEPS)
    );

    TH1D *hHitClusters0p           = new TH1D("hHitClusters0p", "hHitClusters0p;;", 10, 0, 10);
    TH1D *hHitClusters0pBackground = new TH1D("hHitClusters0pBackground", "hHitClusters0pBackground;;", 10, 0, 10);
    TH1D *hHitClustersNp           = new TH1D("hHitClustersNp", "hHitClustersNp;;", 10, 0, 10);
    TH1D *hHitClustersNpBackground = new TH1D("hHitClustersNpBackground", "hHitClustersNpBackground;;", 10, 0, 10);

    TH1D *hHitClustersSize0p           = new TH1D("hHitClustersSize0p", "hHitClustersSize0p;;", 10, 2, 8);
    TH1D *hHitClustersSize0pBackground = new TH1D("hHitClustersSize0pBackground", "hHitClustersSize0pBackground;;", 10, 2, 8);
    TH1D *hHitClustersSizeNp           = new TH1D("hHitClustersSizeNp", "hHitClustersSizeNp;;", 10, 2, 8);
    TH1D *hHitClustersSizeNpBackground = new TH1D("hHitClustersSizeNpBackground", "hHitClustersSizeNpBackground;;", 10, 2, 8);

    TH1D *hHitLargeClusters0p           = new TH1D("hHitLargeClusters0p", "hHitClustersSize0p;;", 10, 0, 10);
    TH1D *hHitLargeClusters0pBackground = new TH1D("hHitLargeClusters0pBackground", "hHitClustersSize0pBackground;;", 10, 0, 10);
    TH1D *hHitLargeClustersNp           = new TH1D("hHitLargeClustersNp", "hHitClustersSizeNp;;", 10, 0, 10);
    TH1D *hHitLargeClustersNpBackground = new TH1D("hHitLargeClustersNpBackground", "hHitClustersSizeNpBackground;;", 10, 0, 10);

    //////////////////////////////////////////
    // Data for inelastic scattering events //
    //////////////////////////////////////////

    TH2D *hTotalBackgroundInelScatteringLengthVSAngle = new TH2D(
        "hTotalBackgroundInelScatteringLengthVSAngle",
        "hTotalBackgroundInelScatteringLengthVSAngle;Length (cm);Angle (rad)",
        10, 0, 40, 
        10, 0, TMath::Pi()
    );

    TH2D *hTotalBackgroundInelScatteringLengthVSVertexKE = new TH2D(
        "hTotalBackgroundInelScatteringLengthVSVertexKE",
        "hTotalBackgroundInelScatteringLengthVSVertexKE;Length (cm); Vertex energy (GeV/c)",
        10, 0, 40,
        8, 0, 0.4
    );

    TH2D *hTotalBackgroundInelScatteringVertexKEVSAngle = new TH2D(
        "hTotalBackgroundInelScatteringVertexKEVSAngle",
        "hTotalBackgroundInelScatteringVertexKEVSAngle;Vertex KE (GeV/c);Angle (rad)",
        8, 0, 0.4,
        10, 0, TMath::Pi()
    );

    TH2D *h0pBackgroundInelScatteringVertexKEVSAngle = new TH2D(
        "h0pBackgroundInelScatteringVertexKEVSAngle",
        "h0pBackgroundInelScatteringVertexKEVSAngle;Vertex KE (GeV/c);Angle (rad)",
        8, 0, 0.4,
        10, 0, TMath::Pi()
    );

    TH2D *hNpBackgroundInelScatteringVertexKEVSAngle = new TH2D(
        "hNpBackgroundInelScatteringVertexKEVSAngle",
        "hNpBackgroundInelScatteringVertexKEVSAngle;Vertex KE (GeV/c);Angle (rad)",
        8, 0, 0.4,
        10, 0, TMath::Pi()
    );

    std::vector<double> InelasticScatteringIncidentKE;
    std::vector<double> InelasticScatteringVertexKE;
    std::vector<double> InelasticScatteringOutgoingKE;

    std::vector<double> Background0pInelasticScatteringIncidentKE;
    std::vector<double> Background0pInelasticScatteringVertexKE;
    std::vector<double> Background0pInelasticScatteringOutgoingKE;

    std::vector<double> BackgroundNpInelasticScatteringIncidentKE;
    std::vector<double> BackgroundNpInelasticScatteringVertexKE;
    std::vector<double> BackgroundNpInelasticScatteringOutgoingKE;

    TH1D *hAllScatteringTotal         = new TH1D("hAllScatteringReconstructionTotal", "hAllScatteringReconstructionTotal;;", 50, 0, 0.4);
    TH1D *hAllScatteringReconstructed = new TH1D("hAllScatteringReconstructionEfficiency", "hAllScatteringReconstructionEfficiency;;", 50, 0, 0.4);

    TH1D *hInelasticScatteringReconstructed = new TH1D("hInelasticScatteringReconstructionEfficiency", "hInelasticScatteringReconstructionEfficiency;;", 20, 0, 0.4);
    TH1D *hInelasticScatteringTotal         = new TH1D("hInelasticScatteringTotal", "hInelasticScatteringTotal;;", 20, 0, 0.4);
    TH1D *hInelasticScatteringVertexKE      = new TH1D("hInelasticScatteringVertexKE", "hInelasticScatteringVertexKE;;", 20, 0, 0.5);

    TH1D *hInelasticScatteringReconstructed0pBkg = new TH1D("hInelasticScatteringReconstructionEfficiency0pBkg", "hInelasticScatteringReconstructionEfficiency0pBkg;;", 20, 0, 0.4);
    TH1D *hInelasticScatteringTotal0pBkg         = new TH1D("hInelasticScatteringTotal0pBkg", "hInelasticScatteringTotal0pBkg;;", 20, 0, 0.4);

    TH1D *hInelasticScatteringReconstructedNpBkg = new TH1D("hInelasticScatteringReconstructionEfficiencyNpBkg", "hInelasticScatteringReconstructionEfficiencyNpBkg;;", 20, 0, 0.4);
    TH1D *hInelasticScatteringTotalNpBkg         = new TH1D("hInelasticScatteringTotalNpBkg", "hInelasticScatteringTotalNpBkg;;", 20, 0, 0.4);

    TH1D *hTotalBackgroundInelScatteringAngle = new TH1D("hTotalBackgroundInelScatteringAngle", "hTotalBackgroundInelScatteringAngle;;", 20, 0, TMath::Pi());
    TH1D *h0pBackgroundInelScatteringAngle    = new TH1D("h0pBackgroundInelScatteringAngle", "h0pBackgroundInelScatteringAngle;;", 20, 0, TMath::Pi());
    TH1D *hNpBackgroundInelScatteringAngle    = new TH1D("hNpBackgroundInelScatteringAngle", "hNpBackgroundInelScatteringAngle;;", 20, 0, TMath::Pi());

    TH1D *hTotalBackgroundInelScatteringLength = new TH1D("hTotalBackgroundInelScatteringLength", "hTotalBackgroundInelScatteringLength;;", 15, 0, 40);
    TH1D *h0pBackgroundInelScatteringLength    = new TH1D("h0pBackgroundInelScatteringLength", "h0pBackgroundInelScatteringLength;;", 15, 0, 40);
    TH1D *hNpBackgroundInelScatteringLength    = new TH1D("hNpBackgroundInelScatteringLength", "hNpBackgroundInelScatteringLength;;", 15, 0, 40);

    TH1D *hTotalBackgroundInelScatteringVertexKE = new TH1D("hTotalBackgroundInelScatteringVertexKE", "hTotalBackgroundInelScatteringVertexKE;;", 20, 0, 0.5);
    TH1D *h0pBackgroundInelScatteringVertexKE    = new TH1D("h0pBackgroundInelScatteringVertexKE", "h0pBackgroundInelScatteringVertexKE;;", 20, 0, 0.5);
    TH1D *hNpBackgroundInelScatteringVertexKE    = new TH1D("hNpBackgroundInelScatteringVertexKE", "hNpBackgroundInelScatteringVertexKE;;", 20, 0, 0.5);

    TH1D *h0pInelasticBackgroundSecondaryInteraction = new TH1D("h0pInelasticBackgroundSecondaryInteraction", "h0pInelasticBackgroundSecondaryInteraction;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D *hNpInelasticBackgroundSecondaryInteraction = new TH1D("hNpInelasticBackgroundSecondaryInteraction", "hNpInelasticBackgroundSecondaryInteraction;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);

    TH1D *hInelasticScatteringAngle     = new TH1D("hInelasticScatteringAngle", "hInelasticScatteringAngle;;", 20, 0, 180);
    TH1D *hInelasticScatteringRecoAngle = new TH1D("hInelasticScatteringRecoAngle", "hInelasthInelasticScatteringRecoAngleicScatteringAngle;;", 20, 0, 180);

    TH1D *hInelasticScatteringEnergyDiff     = new TH1D("hInelasticScatteringEnergyDiff", "hInelasticScatteringEnergyDiff;;", 20, 0, 0.4);
    TH1D *hInelasticScatteringRecoEnergyDiff = new TH1D("hInelasticScatteringRecoEnergyDiff", "hInelasticScatteringRecoEnergyDiff;;", 20, 0, 0.4);

    TH1D *hTotalBackgroundInelasticScatteringEnergyDiff = new TH1D("hTotalBackgroundInelasticScatteringEnergyDiff", "hTotalBackgroundInelasticScatteringEnergyDiff;;", 20, 0, 0.4);
    TH1D *h0pBackgroundInelasticScatteringEnergyDiff    = new TH1D("h0pBackgroundInelasticScatteringEnergyDiff", "h0pBackgroundInelasticScatteringEnergyDiff;;", 20, 0, 0.4);
    TH1D *hNpBackgroundInelasticScatteringEnergyDiff    = new TH1D("hNpBackgroundInelasticScatteringEnergyDiff", "hNpBackgroundInelasticScatteringEnergyDiff;;", 20, 0, 0.4);

    ////////////////////////////////////////
    // Data for elastic scattering events //
    ////////////////////////////////////////

    TH2D *hTotalBackgroundElasticScatteringVertexKEVSAngle = new TH2D(
        "hTotalBackgroundElasticScatteringVertexKEVSAngle",
        "hTotalBackgroundElasticScatteringVertexKEVSAngle;Vertex KE (GeV/c);Angle (rad)",
        8, 0, 0.4,
        10, 0, TMath::Pi()
    );

    TH2D *h0pBackgroundElasticScatteringVertexKEVSAngle = new TH2D(
        "h0pBackgroundElasticScatteringVertexKEVSAngle",
        "h0pBackgroundElasticScatteringVertexKEVSAngle;Vertex KE (GeV/c);Angle (rad)",
        8, 0, 0.4,
        10, 0, TMath::Pi()
    );

    TH2D *hNpBackgroundElasticScatteringVertexKEVSAngle = new TH2D(
        "hNpBackgroundElasticScatteringVertexKEVSAngle",
        "hNpBackgroundElasticScatteringVertexKEVSAngle;Vertex KE (GeV/c);Angle (rad)",
        8, 0, 0.4,
        10, 0, TMath::Pi()
    );

    TH1D *hElasticScatteringRecoAngle = new TH1D("hElasticScatteringRecoAngle", "hElasticScatteringRecoAngle;;", 20, 0, 180);
    TH1D *hElasticScatteringAngle     = new TH1D("hElasticScatteringAngle", "hElasticScatteringAngle;;", 20, 0, 180);

    TH1D *hElasticScatteringVertexKE     = new TH1D("hElasticScatteringVertexKE", "hElasticScatteringVertexKE;;", 20, 0, 0.5);
    TH1D *hElasticScatteringRecoVertexKE = new TH1D("hElasticScatteringRecoVertexKE", "hElasticScatteringRecoVertexKE;;", 20, 0, 0.5);

    TH1D *hTotalBackgroundElasticScatteringAngle = new TH1D("hTotalBackgroundElasticScatteringAngle", "hTotalBackgroundElasticScatteringAngle;;", 20, 0, TMath::Pi());
    TH1D *h0pBackgroundElasticScatteringAngle    = new TH1D("h0pBackgroundElasticScatteringAngle", "h0pBackgroundElasticScatteringAngle;;", 20, 0, TMath::Pi());
    TH1D *hNpBackgroundElasticScatteringAngle    = new TH1D("hNpBackgroundElasticScatteringAngle", "hNpBackgroundElasticScatteringAngle;;", 20, 0, TMath::Pi());

    TH1D *hTotalBackgroundElasticScatteringVertexKE = new TH1D("hTotalBackgroundElasticScatteringVertexKE", "hTotalBackgroundElasticScatteringVertexKE;;", 20, 0, 0.5);
    TH1D *h0pBackgroundElasticScatteringVertexKE    = new TH1D("h0pBackgroundElasticScatteringVertexKE", "h0pBackgroundElasticScatteringVertexKE;;", 20, 0, 0.5);
    TH1D *hNpBackgroundElasticScatteringVertexKE    = new TH1D("hNpBackgroundElasticScatteringVertexKE", "hNpBackgroundElasticScatteringVertexKE;;", 20, 0, 0.5);

    ///////////////////////////////////
    // Histograms for ALL scattering //
    ///////////////////////////////////

    TH2D *hTotalBackgroundAllScatteringVertexKEVSAngle = new TH2D(
        "hTotalBackgroundAllScatteringVertexKEVSAngle",
        "hTotalBackgroundAllScatteringVertexKEVSAngle;Vertex KE (GeV/c);Angle (rad)",
        8, 0, 0.4,
        10, 0, TMath::Pi()
    );

    TH2D *h0pBackgroundAllScatteringVertexKEVSAngle = new TH2D(
        "h0pBackgroundAllScatteringVertexKEVSAngle",
        "h0pBackgroundAllScatteringVertexKEVSAngle;Vertex KE (GeV/c);Angle (rad)",
        8, 0, 0.4,
        10, 0, TMath::Pi()
    );

    TH2D *hNpBackgroundAllScatteringVertexKEVSAngle = new TH2D(
        "hNpBackgroundAllScatteringVertexKEVSAngle",
        "hNpBackgroundAllScatteringVertexKEVSAngle;Vertex KE (GeV/c);Angle (rad)",
        8, 0, 0.4,
        10, 0, TMath::Pi()
    );

    TH1D *hAllScatteringAngle      = new TH1D("hAllScatteringAngle", "hAllScatteringAngle;;", 20, 0, 180);
    TH1D *hAllScatteringRecoAngle  = new TH1D("hAllScatteringRecoAngle", "hAllScatteringRecoAngle;;", 20, 0, 180);

    TH1D *hAllScatteringVertexKE     = new TH1D("hAllScatteringVertexKE", "hAllScatteringVertexKE;;", 20, 0, 0.5);
    TH1D *hAllScatteringRecoVertexKE = new TH1D("hAllScatteringRecoVertexKE", "hAllScatteringRecoVertexKE;;", 20, 0, 0.5);

    TH1D *hTotalBackgroundAllScatteringAngle = new TH1D("hTotalBackgroundAllScatteringAngle", "hTotalBackgroundAllScatteringAngle;;", 20, 0, TMath::Pi());
    TH1D *h0pBackgroundAllScatteringAngle    = new TH1D("h0pBackgroundAllScatteringAngle", "h0pBackgroundAllScatteringAngle;;", 20, 0, TMath::Pi());
    TH1D *hNpBackgroundAllScatteringAngle    = new TH1D("hNpBackgroundAllScatteringAngle", "hNpBackgroundAllScatteringAngle;;", 20, 0, TMath::Pi());

    TH1D *hTotalBackgroundAllScatteringVertexKE = new TH1D("hTotalBackgroundAllScatteringVertexKE", "hTotalBackgroundAllScatteringVertexKE;;", 20, 0, 0.5);
    TH1D *h0pBackgroundAllScatteringVertexKE    = new TH1D("h0pBackgroundAllScatteringVertexKE", "h0pBackgroundAllScatteringVertexKE;;", 20, 0, 0.5);
    TH1D *hNpBackgroundAllScatteringVertexKE    = new TH1D("hNpBackgroundAllScatteringVertexKE", "hNpBackgroundAllScatteringVertexKE;;", 20, 0, 0.5);

    /////////////////////////////////////
    // Data for charge exchange events //
    /////////////////////////////////////

    TH1D* hChExchShowerRecoTrkLengths  = new TH1D("hChExchShowerRecoTrkLengths", "Charge Exchange Shower Reconstructed Track Lengths;Length (cm);Entries", 20, 0, 40);
    TH1D* hChExchShowerTruthTrkLengths = new TH1D("hChExchShowerTruthTrkLengths", "Charge Exchange Shower Truth Track Lengths;Length (cm);Entries", 20, 0, 40);

    TH2D* hChExchShowerDaughtersTrueIncident = new TH2D(
        "hChExchShowerDaughtersTrueIncident",
        "hChExchShowerDaughtersTrueIncident;true daughter_{x} - WC_{x};true daughter_{y} - WC_{y}",
        50, -15, 15,
        50, -15, 15
    );

    TH2D* hChExchShowerDaughtersRecoIncident = new TH2D(
        "hChExchShowerDaughtersRecoIncident",
        "hChExchShowerDaughtersRecoIncident;reco daughter_{x} - WC_{x};reco daughter_{y} - WC_{y}",
        30, -15, 15,
        30, -15, 15
    );

    TH2D* hChExchShowerDaughtersExitAngles = new TH2D(
        "hChExchShowerDaughtersExitAngles",
        "hChExchShowerDaughtersExitAngles;Azimuth (x-y plane) [rad];Polar (away from z)",
        20, -TMath::Pi(), TMath::Pi(),
        20, 0, TMath::Pi()
    );
    TH1D* hChExchShowerDaughtersForwardMom = new TH1D("hChExchShowerDaughtersForwardMom", "hChExchShowerDaughtersForwardMom", 20, 0, 0.5);

    // Sliced cone 
    TH1D* hChExchShowerTrksConeContained       = new TH1D("hChExchShowerTrksConeContained", "hChExchShowerTrksConeContained", 20, 0, 20);
    TH1D* hChExchShowerTrksConeUnContained     = new TH1D("hChExchShowerTrksConeUnContained", "hChExchShowerTrksConeUnContained", 20, 0, 20);
    TH1D* hChExchShowerRecoTrksConeContained   = new TH1D("hChExchShowerRecoTrksConeContained", "hChExchShowerRecoTrksConeContained", 20, 0, 20);
    TH1D* hChExchShowerRecoTrksConeUnContained = new TH1D("hChExchShowerRecoTrksConeUnContained", "hChExchShowerRecoTrksConeUnContained", 20, 0, 20);

    // Cylinder
    TH1D* hChExchShowerTrksCylinderContained       = new TH1D("hChExchShowerTrksCylinderContained", "hChExchShowerTrksCylinderContained", 20, 0, 20);
    TH1D* hChExchShowerTrksCylinderUnContained     = new TH1D("hChExchShowerTrksCylinderUnContained", "hChExchShowerTrksCylinderUnContained", 20, 0, 20);
    TH1D* hChExchShowerRecoTrksCylinderContained   = new TH1D("hChExchShowerRecoTrksCylinderContained", "hChExchShowerRecoTrksCylinderContained", 20, 0, 20);
    TH1D* hChExchShowerRecoTrksCylinderUnContained = new TH1D("hChExchShowerRecoTrksCylinderUnContained", "hChExchShowerRecoTrksCylinderUnContained", 20, 0, 20);

    //////////////////////////////
    // Data for sliced cone cut //
    //////////////////////////////

    TH1D* hSlicedConeTracksPiAbs0p = new TH1D("hSlicedConeTracksPiAbs0p", "hSlicedConeTracksPiAbs0p", 20, 0, 20);
    TH1D* hSlicedConeTracksPiAbsNp = new TH1D("hSlicedConeTracksPiAbsNp", "hSlicedConeTracksPiAbsNp", 20, 0, 20);
    TH1D* hSlicedConeTracksPiScat  = new TH1D("hSlicedConeTracksPiScat", "hSlicedConeTracksPiScat", 20, 0, 20);
    TH1D* hSlicedConeTracksPiChEx  = new TH1D("hSlicedConeTracksPiChEx", "hSlicedConeTracksPiChEx", 20, 0, 20);
    TH1D* hSlicedConeTracksPiOther = new TH1D("hSlicedConeTracksPiOther", "hSlicedConeTracksPiOther", 20, 0, 20);
    TH1D* hSlicedConeTracksMu      = new TH1D("hSlicedConeTracksMu", "hSlicedConeTracksMu", 20, 0, 20);
    TH1D* hSlicedConeTracksE       = new TH1D("hSlicedConeTracksE", "hSlicedConeTracksE", 20, 0, 20);

    TH1D* hSlicedConeNumTrksPiAbs0p = new TH1D("hSlicedConeNumTrksPiAbs0p", "hSlicedConeNumTrksPiAbs0p", 10, 0, 10);
    TH1D* hSlicedConeNumTrksPiAbsNp = new TH1D("hSlicedConeNumTrksPiAbsNp", "hSlicedConeNumTrksPiAbsNp", 10, 0, 10);
    TH1D* hSlicedConeNumTrksPiScat  = new TH1D("hSlicedConeNumTrksPiScat", "hSlicedConeNumTrksPiScat", 10, 0, 10);
    TH1D* hSlicedConeNumTrksPiChEx  = new TH1D("hSlicedConeNumTrksPiChEx", "hSlicedConeNumTrksPiChEx", 10, 0, 10);
    TH1D* hSlicedConeNumTrksPiOther = new TH1D("hSlicedConeNumTrksPiOther", "hSlicedConeNumTrksPiOther", 10, 0, 10);
    TH1D* hSlicedConeNumTrksMu      = new TH1D("hSlicedConeNumTrksMu", "hSlicedConeNumTrksMu", 10, 0, 10);
    TH1D* hSlicedConeNumTrksE       = new TH1D("hSlicedConeNumTrksE", "hSlicedConeNumTrksE", 10, 0, 10);

    TH1D* hSlicedConeNumSmallTrksPiAbs0p = new TH1D("hSlicedConeNumSmallTrksPiAbs0p", "hSlicedConeNumSmallTrksPiAbs0p", 10, 0, 10);
    TH1D* hSlicedConeNumSmallTrksPiAbsNp = new TH1D("hSlicedConeNumSmallTrksPiAbsNp", "hSlicedConeNumSmallTrksPiAbsNp", 10, 0, 10);
    TH1D* hSlicedConeNumSmallTrksPiScat  = new TH1D("hSlicedConeNumSmallTrksPiScat", "hSlicedConeNumSmallTrksPiScat", 10, 0, 10);
    TH1D* hSlicedConeNumSmallTrksPiChEx  = new TH1D("hSlicedConeNumSmallTrksPiChEx", "hSlicedConeNumSmallTrksPiChEx", 10, 0, 10);
    TH1D* hSlicedConeNumSmallTrksPiOther = new TH1D("hSlicedConeNumSmallTrksPiOther", "hSlicedConeNumSmallTrksPiOther", 10, 0, 10);
    TH1D* hSlicedConeNumSmallTrksMu      = new TH1D("hSlicedConeNumSmallTrksMu", "hSlicedConeNumSmallTrksMu", 10, 0, 10);
    TH1D* hSlicedConeNumSmallTrksE       = new TH1D("hSlicedConeNumSmallTrksE", "hSlicedConeNumSmallTrksE", 10, 0, 10);

    ///////////////////////////
    // Data for cylinder cut //
    ///////////////////////////

    int numSurvivingPions     = 0;
    int numSurvivingMuons     = 0;
    int numSurvivingElectrons = 0;

    TH1D* hNumTracksInCylinder          = new TH1D("hNumTracksInCylinder", "hNumTracksInCylinder", 10, 0, 10);
    TH1D* hNumTracksInCylinderPions     = new TH1D("hNumTracksInCylinderPions", "hNumTracksInCylinderPions", 10, 0, 10);
    TH1D* hNumTracksInCylinderMuons     = new TH1D("hNumTracksInCylinderMuons", "hNumTracksInCylinderMuons", 10, 0, 10);
    TH1D* hNumTracksInCylinderElectrons = new TH1D("hNumTracksInCylinderElectrons", "hNumTracksInCylinderElectrons", 10, 0, 10);

    TH1D* hTrkLengthsInCylinder          = new TH1D("hTrkLengthsInCylinder", "hTrkLengthsInCylinder", 50, 0, 10);
    TH1D* hTrkLengthsInCylinderPions     = new TH1D("hTrkLengthsInCylinderPions", "hTrkLengthsInCylinderPions", 50, 0, 10);
    TH1D* hTrkLengthsInCylinderMuons     = new TH1D("hTrkLengthsInCylinderMuons", "hTrkLengthsInCylinderMuons", 50, 0, 10);
    TH1D* hTrkLengthsInCylinderElectrons = new TH1D("hTrkLengthsInCylinderElectrons", "hTrkLengthsInCylinderElectrons", 50, 0, 10);

    TH1D* hTrkLengthsInCylinderPionAbs0p      = new TH1D("hTrkLengthsInCylinderPionAbs0p", "hTrkLengthsInCylinderPionAbs0p", 50, 0, 10);
    TH1D* hTrkLengthsInCylinderPionAbsNp      = new TH1D("hTrkLengthsInCylinderPionAbsNp", "hTrkLengthsInCylinderPionAbsNp", 50, 0, 10);
    TH1D* hTrkLengthsInCylinderPionScattering = new TH1D("hTrkLengthsInCylinderPionScattering", "hTrkLengthsInCylinderPionScattering", 50, 0, 10);
    TH1D* hTrkLengthsInCylinderPionChExch     = new TH1D("hTrkLengthsInCylinderPionChExch", "hTrkLengthsInCylinderPionChExch", 50, 0, 10);
    TH1D* hTrkLengthsInCylinderPionOther      = new TH1D("hTrkLengthsInCylinderPionOther", "hTrkLengthsInCylinderPionOther", 50, 0, 10);

    TH1D* hSmallTrksInCylinder          = new TH1D("hSmallTrksInCylinder", "hSmallTrksInCylinder", 10, 0, 10);
    TH1D* hSmallTrksInCylinderPions     = new TH1D("hSmallTrksInCylinderPions", "hSmallTrksInCylinderPions", 10, 0, 10);
    TH1D* hSmallTrksInCylinderMuons     = new TH1D("hSmallTrksInCylinderMuons", "hSmallTrksInCylinderMuons", 10, 0, 10);
    TH1D* hSmallTrksInCylinderElectrons = new TH1D("hSmallTrksInCylinderElectrons", "hSmallTrksInCylinderElectrons", 10, 0, 10);

    TH1D* hSmallTrksInCylinderPionAbs0p      = new TH1D("hSmallTrksInCylinderPionAbs0p", "hSmallTrksInCylinderPionAbs0p", 10, 0, 10);
    TH1D* hSmallTrksInCylinderPionAbsNp      = new TH1D("hSmallTrksInCylinderPionAbsNp", "hSmallTrksInCylinderPionAbsNp", 10, 0, 10);
    TH1D* hSmallTrksInCylinderPionScattering = new TH1D("hSmallTrksInCylinderPionScattering", "hSmallTrksInCylinderPionScattering", 10, 0, 10);
    TH1D* hSmallTrksInCylinderPionChExch     = new TH1D("hSmallTrksInCylinderPionChExch", "hSmallTrksInCylinderPionChExch", 10, 0, 10);
    TH1D* hSmallTrksInCylinderPionOther      = new TH1D("hSmallTrksInCylinderPionOther", "hSmallTrksInCylinderPionOther", 10, 0, 10);

    //////////////////////
    // Beamline showers //
    //////////////////////

    TH1D* hElectronShowerTrueTrkLengths    = new TH1D("hElectronShowerTrueTrkLengths", "Electron Shower Reconstructed Track Lengths;Length (cm);Entries", 20, 0, 20);
    TH1D* hElectronShowerRecoTrkLengths    = new TH1D("hElectronShowerRecoTrkLengths", "Electron Shower Reconstructed Track Lengths;Length (cm);Entries", 20, 0, 20);
    TH1D* hElectronShowerTrkContained      = new TH1D("hElectronShowerTrkContained", "Electron Shower Reconstructed Track Contained;Length (cm);Entries", 20, 0, 20);
    TH1D* hElectronShowerTrkUnContained    = new TH1D("hElectronShowerTrkUnContained", "Electron Shower Reconstructed Track UnContained;Length (cm);Entries", 20, 0, 20);
    TH1D* hElectronShowerTrkContainedRatio = new TH1D("hElectronShowerTrkContainedRatio", "Electron Shower Reconstructed Track Contained Ratio;Length (cm);Entries", 20, 0, 20);

    TH2D* hElectronShowerDaughtersTrueIncident = new TH2D(
        "hElectronShowerDaughtersTrueIncident",
        "hElectronShowerDaughtersTrueIncident;true daughter_{x} - WC_{x} [cm];true daughter_{y} - WC_{y} [cm]",
        50, -15, 15,
        50, -15, 15
    );

    TH2D* hElectronShowerDaughtersRecoIncident = new TH2D(
        "hElectronShowerDaughtersRecoIncident",
        "hElectronShowerDaughtersRecoIncident;reco daughter_{x} - WC_{x} [cm];reco daughter_{y} - WC_{y} [cm]",
        50, -15, 15,
        50, -15, 15
    );

    ////////////////////////
    // Multiple primaries //
    ////////////////////////

    TH1D* hNumPrimaries      = new TH1D("hNumPrimaries", "Number of Primaries;N_{primaries};Entries", 10, 0, 10);
    TH1D* hNumValidPrimaries = new TH1D("hNumValidPrimaries", "Number of Valid Primaries;N_{valid primaries};Entries", 10, 0, 10);

    //////////////////////////////
    // Cross-section histograms //
    //////////////////////////////

    // True cross sections
    TH1D* hTrueIncidentKE  = (TH1D*) TrueXSDirectory->Get("hIncidentKE");
    TH1D* hTruePionAbsKE   = (TH1D*) TrueXSDirectory->Get("hInteractingKEPionAbs");
    TH1D* hTruePionAbs0pKE = (TH1D*) TrueXSDirectory->Get("hInteractingKEPionAbs0p");
    TH1D* hTruePionAbsNpKE = (TH1D*) TrueXSDirectory->Get("hInteractingKEPionAbsNp");

    TH1D* hTruePionScatteringKE   = (TH1D*) TrueXSDirectory->Get("hInteractingKEScattering");
    TH1D* hTruePion0pScatteringKE = (TH1D*) TrueXSDirectory->Get("hInteractingKE0pScattering");
    TH1D* hTruePionNpScatteringKE = (TH1D*) TrueXSDirectory->Get("hInteractingKENpScattering");

    TH1D* hPionAbsTrueCrossSection   = (TH1D*) TrueXSDirectory->Get("hCrossSectionPionAbs");
    TH1D* hPion0pAbsTrueCrossSection = (TH1D*) TrueXSDirectory->Get("hCrossSectionPionAbs0p");
    TH1D* hPionNpAbsTrueCrossSection = (TH1D*) TrueXSDirectory->Get("hCrossSectionPionAbsNp");

    // Reco cross sections
    TH1D *hIncidentKE  = new TH1D("hRecoIncidentKE", "Incident KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hPionAbsKE   = new TH1D("hPionAbsKE", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *h0pPionAbsKE = new TH1D("h0pPionAbsKE", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hNpPionAbsKE = new TH1D("hNpPionAbsKE", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    
    TH1D *hIncidentKECorrect   = new TH1D("hRecoIncidentKECorrect", "Incident KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hIncidentKEMuons     = new TH1D("hIncidentKEMuons", "Incident KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hIncidentKEElectrons = new TH1D("hIncidentKEElectrons", "Incident KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Background breakdown for all absorption
    TH1D *hPionAbsKECorrect      = new TH1D("hPionAbsKECorrect", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hPionAbsKE0pScattering = new TH1D("hPionAbsKE0pScattering", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hPionAbsKENpScattering = new TH1D("hPionAbsKENpScattering", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hPionAbsKEChEx         = new TH1D("hPionAbsKEChEx", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hPionAbsKEMuons        = new TH1D("hPionAbsKEMuons", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hPionAbsKEElectrons    = new TH1D("hPionAbsKEElectrons", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hPionAbsKEOther        = new TH1D("hPionAbsKEOther", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Background breakdown for 0p absorption
    TH1D *h0pPionAbsKECorrect      = new TH1D("h0pPionAbsKECorrect", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *h0pPionAbsKENpAbs        = new TH1D("h0pPionAbsKENpAbs", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *h0pPionAbsKE0pScattering = new TH1D("h0pPionAbsKE0pScattering", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *h0pPionAbsKENpScattering = new TH1D("h0pPionAbsKENpScattering", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *h0pPionAbsKEChEx         = new TH1D("h0pPionAbsKEChEx", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *h0pPionAbsKEMuons        = new TH1D("h0pPionAbsKEMuons", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *h0pPionAbsKEElectrons    = new TH1D("h0pPionAbsKEElectrons", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *h0pPionAbsKEOther        = new TH1D("h0pPionAbsKEOther", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Background breakdown for Np absorption
    TH1D *hNpPionAbsKECorrect      = new TH1D("hNpPionAbsKECorrect", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hNpPionAbsKE0pAbs        = new TH1D("hNpPionAbsKE0pAbs", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hNpPionAbsKE0pScattering = new TH1D("hNpPionAbsKE0pScattering", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hNpPionAbsKENpScattering = new TH1D("hNpPionAbsKENpScattering", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hNpPionAbsKEChEx         = new TH1D("hNpPionAbsKEChEx", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hNpPionAbsKEMuons        = new TH1D("hNpPionAbsKEMuons", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hNpPionAbsKEElectrons    = new TH1D("hNpPionAbsKEElectrons", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hNpPionAbsKEOther        = new TH1D("hNpPionAbsKEOther", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());

    TH1D *hPionAbsRecoCrossSection   = new TH1D("hPionAbsRecoCrossSection", "Cross section [barn]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hPion0pAbsRecoCrossSection = new TH1D("hPion0pAbsRecoCrossSection", "Cross section [barn]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hPionNpAbsRecoCrossSection = new TH1D("hPionNpAbsRecoCrossSection", "Cross section [barn]", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Histograms for corrections
    TH1D *hIncidentKEOnlyPions  = new TH1D("hIncidentKEOnlyPions", "Incident KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());

    TH1D *hPionAbsKEOnlyAbs     = new TH1D("hPionAbsKEOnlyAbs", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *h0pPionAbsKEOnly0pAbs = new TH1D("h0pPionAbsKEOnly0pAbs", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hNpPionAbsKEOnlyNpAbs = new TH1D("hNpPionAbsKEOnlyNpAbs", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());

    ///////////////////////////
    // Scattering histograms //
    ///////////////////////////

    TH1D *hPionScatteringKE   = new TH1D("hPionScatteringKE", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hPionScatteringKE0p = new TH1D("hPionScatteringKE0p", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hPionScatteringKENp = new TH1D("hPionScatteringKENp", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Background breakdown for all scattering
    TH1D *hPionScatteringKECorrect   = new TH1D("hPionScatteringKECorrect", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hPionScatteringKE0pAbs     = new TH1D("hPionScatteringKE0pAbs", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hPionScatteringKENpAbs     = new TH1D("hPionScatteringKENpAbs", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hPionScatteringKEChEx      = new TH1D("hPionScatteringKEChEx", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hPionScatteringKEMuons     = new TH1D("hPionScatteringKEMuons", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hPionScatteringKEElectrons = new TH1D("hPionScatteringKEElectrons", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hPionScatteringKEOther     = new TH1D("hPionScatteringKEOther", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Background breakdown for 0p scattering
    TH1D *h0pPionScatteringKECorrect   = new TH1D("h0pPionScatteringKECorrect", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *h0pPionScatteringKENpScatter = new TH1D("h0pPionScatteringKENpScatter", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *h0pPionScatteringKE0pAbs     = new TH1D("h0pPionScatteringKE0pAbs", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *h0pPionScatteringKENpAbs     = new TH1D("h0pPionScatteringKENpAbs", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *h0pPionScatteringKEChEx      = new TH1D("h0pPionScatteringKEChEx", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *h0pPionScatteringKEMuons     = new TH1D("h0pPionScatteringKEMuons", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *h0pPionScatteringKEElectrons = new TH1D("h0pPionScatteringKEElectrons", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *h0pPionScatteringKEOther     = new TH1D("h0pPionScatteringKEOther", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    
    // Background breakdown for Np scattering
    TH1D *hNpPionScatteringKECorrect   = new TH1D("hNpPionScatteringKECorrect", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hNpPionScatteringKE0pScatter = new TH1D("hNpPionScatteringKE0pScatter", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hNpPionScatteringKE0pAbs     = new TH1D("hNpPionScatteringKE0pAbs", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hNpPionScatteringKENpAbs     = new TH1D("hNpPionScatteringKENpAbs", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hNpPionScatteringKEChEx      = new TH1D("hNpPionScatteringKEChEx", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hNpPionScatteringKEMuons     = new TH1D("hNpPionScatteringKEMuons", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hNpPionScatteringKEElectrons = new TH1D("hNpPionScatteringKEElectrons", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D *hNpPionScatteringKEOther     = new TH1D("hNpPionScatteringKEOther", "Interacting KE [MeV]", NUM_BINS_KE, ARRAY_KE_BINS.data());

    /////////////////////////////////
    // Files for event information //
    /////////////////////////////////

    std::ofstream outFile0pTrue("files/RecoAllAnalysis/0pRecoTrue.txt");
    std::ofstream outFileNpTrue("files/RecoAllAnalysis/NpRecoTrue.txt");
    std::ofstream outFile0pBackground("files/RecoAllAnalysis/0pBackground.txt");
    std::ofstream outFileNpBackground("files/RecoAllAnalysis/NpBackground.txt");
    std::ofstream outStitchedFile("files/RecoAllAnalysis/WCMatchStitching.txt");
    std::ofstream outWCAll("files/RecoAllAnalysis/WCMatchAllChi2.txt");

    /////////////////////////////////////
    // Plot energy deposition profiles //
    /////////////////////////////////////

    gProton->SetLineColor(kRed+1);
    gProton->SetLineWidth(2);
    gProton->SetTitle("Proton");

    gPion->SetLineColor(kRed+3);
    gPion->SetLineWidth(2);
    gPion->SetTitle("Pion");

    gMuonTG->SetLineColor(kOrange+7);
    gMuonTG->SetLineWidth(2);
    gMuonTG->SetTitle("MIP");

    // Set axis labels using one of the graphs
    gProton->SetTitle("dE/dx vs Residual Range;Residual Range [cm];dE/dx [MeV/cm]");

    // Create a canvas and draw the graphs
    TCanvas* c1 = new TCanvas("c1", "dE/dx vs Residual Range", 800, 600);
    gProton->Draw("AL");
    gPion->Draw("L SAME");
    gMuonTG->Draw("L SAME");

    // Add a legend
    TLegend* legend = new TLegend(0.65, 0.70, 0.88, 0.88);
    legend->AddEntry(gProton, "Proton", "l");
    legend->AddEntry(gPion,   "Pion",   "l");
    legend->AddEntry(gMuonTG, "MIP",   "l");
    legend->Draw();

    // Optional: set grid
    c1->SetGrid();
    c1->SaveAs(SaveDir +  "dEdxProfiles/AllProfiles.png");
    delete c1;

    // Absorption Np vertex reco
    int betterVertexNp = 0;
    int worseVertexNp  = 0;

    //////////////////////
    // Loop over events //
    //////////////////////

    Int_t NumEntries = (Int_t) tree->GetEntries();
    std::cout << "Num entries: " << NumEntries << std::endl;

    for (Int_t i = 0; i < NumEntries; ++i) {
        tree->GetEntry(i);

        // Make it go faster
        // if (i > 10000) break;

        std::vector<int>             thisEventPrimaryDaughterPDG = *truthPrimaryDaughtersPDG;
        std::vector<double>           thisEventPrimaryDaughterKE = *truthPrimaryDaughtersKE;
        // std::vector<std::string> thisEventPrimaryDaughterProcess = *truthPrimaryDaughtersProcess;

        // Sanity check
        removeRepeatedPoints(WC2TPCLocationsX, WC2TPCLocationsY, WC2TPCLocationsZ);

        // Scattering only if degree > THRESHOLD
        double scatteringAngle = -9999;
        if (backgroundType == 12 || backgroundType == 6) {
            if (backgroundType == 12) {
                scatteringAngle = trajectoryInteractionAngle;
            } else if (backgroundType == 6) {
                scatteringAngle = truthScatteringAngle;
            }

            if (scatteringAngle < SCATTERING_ANGLE_THRESHOLD) {
                // Use secondary interaction
                for (int iInteraction = 0; iInteraction < secondaryInteractionTypes->size(); ++iInteraction) {
                    int currentInteraction = secondaryInteractionTypes->at(iInteraction);
                    scatteringAngle        = secondaryInteractionAngle->at(iInteraction);

                    // for (int iContribution = 0; iContribution < secondaryIncidentKEContributions->at(iInteraction).size(); ++iContribution) {
                    //     trueIncidentKEContributions->push_back(secondaryIncidentKEContributions->at(iInteraction)[iContribution]);
                    // }

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

                        thisEventPrimaryDaughterPDG     = secondaryInteractionDaughtersPDG->at(iInteraction);
                        thisEventPrimaryDaughterKE      = secondaryInteractionDaughtersKE->at(iInteraction);
                        // thisEventPrimaryDaughterProcess = secondaryInteractionDaughtersProcess->at(iInteraction);

                        // Re-do variables (hacky)
                        if (currentInteraction == 12) trajectoryInteractionAngle = scatteringAngle;
                        else if (currentInteraction == 6) {
                            truthScatteringAngle = scatteringAngle;
                            for (int i = 0; i < thisEventPrimaryDaughterPDG.size(); ++i) {
                                if (thisEventPrimaryDaughterPDG[i] == -211) {
                                    truthScatteredPionKE = secondaryInteractionDaughtersKE->at(iInteraction)[i];
                                    break;
                                }
                            }
                        }
                        break;
                    }
                }
            } else {
                // If initial interaction has valid angle, just correct interaction KE for 
                // elastic scattering interactions
                if (backgroundType == 12) truthPrimaryVertexKE = trajectoryInteractionKE;
            }
        }
        hTotalEvents->Fill(backgroundType);

        // Get data about truth-level secondary particles: pions and protons
        for (int i = 0; i < thisEventPrimaryDaughterPDG.size(); ++i) {
            if (thisEventPrimaryDaughterPDG[i] == -211) {
                hAllSecondaryPions->Fill(thisEventPrimaryDaughterKE[i]);
            } else if (thisEventPrimaryDaughterPDG[i] == 2212) {
                hAllSecondaryProtons->Fill(thisEventPrimaryDaughterKE[i]);
            }
        }

        // Categorize scatterings into 0p and Np scatterings
        int scatteringType = -1; // -1: not scattering, 0: 0p scattering, 1: Np scattering
        if (backgroundType == 12 || (backgroundType == 6 && numVisibleProtons == 0)) {
            scatteringType = 0;
        }
        else if (backgroundType == 6 && numVisibleProtons > 0) {
            scatteringType = 1;
        }

        // Multiple primaries
        hNumPrimaries->Fill(numPrimaries);
        hNumValidPrimaries->Fill(numValidPrimaries);

        int validPrimaryIdx = -1;
        for (size_t i = 0; i < primariesID->size(); ++i) {
            if (primariesID->at(i) == truthPrimaryID) {
                validPrimaryIdx = i;
                break;
            }
        }

        /////////////////////////
        // Construct cylinders //
        /////////////////////////

        // Get truth cylinder for cuts
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

        if (truthNumTail > 0) {
            std::vector<double> avgDir = getAverageDir(truth_points);

            // Extrapolate track to end
            double scale = (maxZ - truth_points.back()[2]) / avgDir[2];
            truthCylinderLocationX->push_back(truth_points.back()[0] + scale * avgDir[0]);
            truthCylinderLocationY->push_back(truth_points.back()[1] + scale * avgDir[1]);
            truthCylinderLocationZ->push_back(truth_points.back()[2] + scale * avgDir[2]);
        }

        // Get reco cylinder for cuts 

        // Copy WC2TPCLocations
        std::vector<double>* wcX = new std::vector<double>(*WC2TPCLocationsX);
        std::vector<double>* wcY = new std::vector<double>(*WC2TPCLocationsY);
        std::vector<double>* wcZ = new std::vector<double>(*WC2TPCLocationsZ);

        // Get direction to end cylinder
        int wcNumPoints = wcX->size();
        int wcNumTail   = std::min(10, wcNumPoints - 1);
        std::vector<std::vector<double>> wc_points;
        for (int j = wcNumPoints - wcNumTail; j < wcNumPoints; ++j) {
            wc_points.push_back({
                wcX->at(j),
                wcY->at(j),
                wcZ->at(j)
            });
        }
        if (wcNumTail > 0) {
            std::vector<double> avgDir = getAverageDir(wc_points);

            // Extrapolate track to end
            double scale = (maxZ - wc_points.back()[2]) / avgDir[2];
            wcX->push_back(wc_points.back()[0] + scale * avgDir[0]);
            wcY->push_back(wc_points.back()[1] + scale * avgDir[1]);
            wcZ->push_back(wc_points.back()[2] + scale * avgDir[2]);
        }

        //////////////////////////////
        // Study truth interactions //
        //////////////////////////////

        // Study elastic and inelastic scattering events
        if (backgroundType == 12) {
            // Elastic scattering event
            hElasticScatteringAngle->Fill(trajectoryInteractionAngle * (180. / TMath::Pi()));
            hElasticScatteringVertexKE->Fill(trajectoryInteractionKE);
            hAllScatteringTotal->Fill(trajectoryInteractionKE);

            hAllScatteringAngle->Fill(trajectoryInteractionAngle * (180. / TMath::Pi()));
            hAllScatteringVertexKE->Fill(trajectoryInteractionKE);

            for (int iRecoTrk = 0; iRecoTrk < matchedIdentity->size(); ++iRecoTrk) {
                if (recoTrkID->at(iRecoTrk) == WC2TPCtrkID) continue;

                if (
                    (matchedIdentity->at(iRecoTrk) == -211) && // pion
                    (distance(
                        recoBeginX->at(iRecoTrk),
                        trajectoryInteractionX,
                        recoBeginY->at(iRecoTrk),
                        trajectoryInteractionY,
                        recoBeginZ->at(iRecoTrk),
                        trajectoryInteractionZ
                    ) < VERTEX_RADIUS) // within 5 cm of scattering vertex
                ) {
                    hAllScatteringReconstructed->Fill(trajectoryInteractionKE);
                    hElasticScatteringRecoAngle->Fill(trajectoryInteractionAngle * (180. / TMath::Pi()));
                    hElasticScatteringRecoVertexKE->Fill(trajectoryInteractionKE);
                    hAllScatteringRecoAngle->Fill(trajectoryInteractionAngle * (180. / TMath::Pi()));
                    hAllScatteringRecoVertexKE->Fill(trajectoryInteractionKE);
                    break;
                }
            }
        } else if (backgroundType == 6) {
            // Inelastic scattering event
            hInelasticScatteringTotal->Fill(truthScatteredPionKE);
            hAllScatteringTotal->Fill(truthScatteredPionKE);
            hInelasticScatteringVertexKE->Fill(truthPrimaryVertexKE);
            hInelasticScatteringAngle->Fill(truthScatteringAngle * (180. / TMath::Pi()));
            hInelasticScatteringEnergyDiff->Fill(truthPrimaryVertexKE - truthScatteredPionKE);

            hAllScatteringAngle->Fill(truthScatteringAngle * (180. / TMath::Pi()));
            hAllScatteringVertexKE->Fill(truthPrimaryVertexKE);

            InelasticScatteringIncidentKE.push_back(truthPrimaryIncidentKE);
            InelasticScatteringVertexKE.push_back(truthPrimaryVertexKE);
            InelasticScatteringOutgoingKE.push_back(truthScatteredPionKE);

            hScatteringVertexKEVsOutKE->Fill(truthPrimaryVertexKE, truthScatteredPionKE);

            for (int iRecoTrk = 0; iRecoTrk < matchedIdentity->size(); ++iRecoTrk) {
                if (recoTrkID->at(iRecoTrk) == WC2TPCtrkID) continue;
                
                if (
                    (matchedIdentity->at(iRecoTrk) == -211) &&          // pion
                    (matchedProcess->at(iRecoTrk) == "pi-Inelastic") && // inelastic process
                    (distance(
                        recoBeginX->at(iRecoTrk),
                        truthPrimaryVertexX,
                        recoBeginY->at(iRecoTrk),
                        truthPrimaryVertexY,
                        recoBeginZ->at(iRecoTrk),
                        truthPrimaryVertexZ
                    ) < VERTEX_RADIUS)
                ) {
                    hInelasticScatteringReconstructed->Fill(truthScatteredPionKE);
                    hAllScatteringReconstructed->Fill(truthScatteredPionKE);
                    hInelasticScatteringRecoAngle->Fill(truthScatteringAngle * (180. / TMath::Pi()));
                    hInelasticScatteringRecoEnergyDiff->Fill(truthPrimaryVertexKE - truthScatteredPionKE);

                    hAllScatteringRecoAngle->Fill(truthScatteringAngle * (180. / TMath::Pi()));
                    hAllScatteringRecoVertexKE->Fill(truthScatteredPionKE);
                    break;
                }
            }
        }

        // Study charge exchange events
        if (backgroundType == 7) {
            // Look at truth tracks that make up the shower
            for (int iTruthTrk = 0; iTruthTrk < chExchShowerIDs->size(); ++iTruthTrk) {
                hChExchShowerDaughtersTrueIncident->Fill(
                    chExchShowerStart->at(iTruthTrk)[0] - primariesStartX->at(validPrimaryIdx), 
                    chExchShowerStart->at(iTruthTrk)[1] - primariesStartY->at(validPrimaryIdx)
                );

                // Check if track is contained inside sliced cone
                double coneDirX = primariesEndX->at(validPrimaryIdx) - primariesStartX->at(validPrimaryIdx);
                double coneDirY = primariesEndY->at(validPrimaryIdx) - primariesStartY->at(validPrimaryIdx);
                double coneDirZ = primariesEndZ->at(validPrimaryIdx) - primariesStartZ->at(validPrimaryIdx);

                bool startsInCone = IsPointInsideSlicedCone(
                    chExchShowerStart->at(iTruthTrk)[0], chExchShowerStart->at(iTruthTrk)[1], chExchShowerStart->at(iTruthTrk)[2],
                    primariesEndX->at(validPrimaryIdx), primariesEndY->at(validPrimaryIdx), primariesEndZ->at(validPrimaryIdx),
                    coneDirX, coneDirY, coneDirZ,
                    SLICED_CONE_HEIGHT, SLICED_CONE_MIN_RADIUS, SLICED_CONE_MAX_RADIUS
                );

                bool endsInCone = IsPointInsideSlicedCone(
                    chExchShowerEnd->at(iTruthTrk)[0], chExchShowerEnd->at(iTruthTrk)[1], chExchShowerEnd->at(iTruthTrk)[2],
                    primariesEndX->at(validPrimaryIdx), primariesEndY->at(validPrimaryIdx), primariesEndZ->at(validPrimaryIdx),
                    coneDirX, coneDirY, coneDirZ,
                    SLICED_CONE_HEIGHT, SLICED_CONE_MIN_RADIUS, SLICED_CONE_MAX_RADIUS
                );

                // Check if track is contained inside cylinder
                bool startInCylinder = IsPointInsideTrackCylinder(
                    truthCylinderLocationX, truthCylinderLocationY, truthCylinderLocationZ,
                    chExchShowerStart->at(iTruthTrk)[0], chExchShowerStart->at(iTruthTrk)[1], chExchShowerStart->at(iTruthTrk)[2],
                    CYLINDER_RADIUS
                );
                bool endInCylinder = IsPointInsideTrackCylinder(
                    truthCylinderLocationX, truthCylinderLocationY, truthCylinderLocationZ,
                    chExchShowerEnd->at(iTruthTrk)[0], chExchShowerEnd->at(iTruthTrk)[1], chExchShowerEnd->at(iTruthTrk)[2],
                    CYLINDER_RADIUS
                );

                // We do not care about photons because they won't give us anything to reconstruct themselves
                if (chExchShowerPDGs->at(iTruthTrk) != 22) {
                    hChExchShowerTruthTrkLengths->Fill(chExchShowerLengths->at(iTruthTrk));

                    if (startsInCone && endsInCone) {
                        hChExchShowerTrksConeContained->Fill(chExchShowerLengths->at(iTruthTrk));
                    } else {
                        hChExchShowerTrksConeUnContained->Fill(chExchShowerLengths->at(iTruthTrk));
                    }

                    if (startInCylinder && endInCylinder) {
                        hChExchShowerTrksCylinderContained->Fill(chExchShowerLengths->at(iTruthTrk));
                    } else {
                        hChExchShowerTrksCylinderUnContained->Fill(chExchShowerLengths->at(iTruthTrk));
                    }
                }
            }

            for (int iNeutralPionDaughter = 0; iNeutralPionDaughter < chExchShowerNeutralPionDaughtersID->size(); ++iNeutralPionDaughter) {
                double px = chExchShowerNeutralPionDaughtersMom->at(iNeutralPionDaughter)[0];
                double py = chExchShowerNeutralPionDaughtersMom->at(iNeutralPionDaughter)[1];
                double pz = chExchShowerNeutralPionDaughtersMom->at(iNeutralPionDaughter)[2];

                double azimuth = TMath::ATan2(py, px);
                double polar   = TMath::ATan2(TMath::Sqrt(px*px + py*py), pz);

                hChExchShowerDaughtersExitAngles->Fill(azimuth, polar);
                hChExchShowerDaughtersForwardMom->Fill(pz);
            }

            // Look at reco tracks that match to tracks in the shower
            if (WC2TPCtrkID != -99999) {
                for (int iRecoTrk = 0; iRecoTrk < matchedTrkID->size(); ++iRecoTrk) {
                    if (recoTrkID->at(iRecoTrk) == WC2TPCtrkID) continue;

                    if (std::find(chExchShowerIDs->begin(), chExchShowerIDs->end(), matchedTrkID->at(iRecoTrk)) != chExchShowerIDs->end()) {
                        double trk_length = TMath::Sqrt(
                            TMath::Power(recoEndX->at(iRecoTrk) - recoBeginX->at(iRecoTrk), 2) +
                            TMath::Power(recoEndY->at(iRecoTrk) - recoBeginY->at(iRecoTrk), 2) +
                            TMath::Power(recoEndZ->at(iRecoTrk) - recoBeginZ->at(iRecoTrk), 2)
                        );
                        hChExchShowerRecoTrkLengths->Fill(trk_length);

                        hChExchShowerDaughtersRecoIncident->Fill(
                            recoBeginX->at(iRecoTrk) - WC2TPCPrimaryBeginX,
                            recoBeginY->at(iRecoTrk) - WC2TPCPrimaryBeginY
                        );
                    }
                }
            }
        }

        // Study electron showers
        if (backgroundType == 3) {
            // Start at iTruthTrk = 1 to avoid counting electron itself
            for (int iTruthTrk = 1; iTruthTrk < electronShowerIDs->size(); ++iTruthTrk) {
                hElectronShowerDaughtersTrueIncident->Fill(
                    electronShowerStart->at(iTruthTrk)[0] - primariesStartX->at(validPrimaryIdx), 
                    electronShowerStart->at(iTruthTrk)[1] - primariesStartY->at(validPrimaryIdx)
                );

                bool startInCylinder = IsPointInsideTrackCylinder(
                    truthCylinderLocationX, truthCylinderLocationY, truthCylinderLocationZ,
                    electronShowerStart->at(iTruthTrk)[0], electronShowerStart->at(iTruthTrk)[1], electronShowerStart->at(iTruthTrk)[2],
                    CYLINDER_RADIUS
                );
                bool endInCylinder = IsPointInsideTrackCylinder(
                    truthCylinderLocationX, truthCylinderLocationY, truthCylinderLocationZ,
                    electronShowerEnd->at(iTruthTrk)[0], electronShowerEnd->at(iTruthTrk)[1], electronShowerEnd->at(iTruthTrk)[2],
                    CYLINDER_RADIUS
                );

                if (electronShowerPDGs->at(iTruthTrk) != 22) {
                    hElectronShowerTrueTrkLengths->Fill(electronShowerLengths->at(iTruthTrk));
                    if (startInCylinder && endInCylinder) {
                        hElectronShowerTrkContained->Fill(electronShowerLengths->at(iTruthTrk));
                    } else {
                        hElectronShowerTrkUnContained->Fill(electronShowerLengths->at(iTruthTrk));
                    }
                }
            }

            // Look at reco tracks that match to tracks in the shower
            if (WC2TPCtrkID != -99999) {
                for (int iRecoTrk = 0; iRecoTrk < matchedTrkID->size(); ++iRecoTrk) {
                    if (recoTrkID->at(iRecoTrk) == WC2TPCtrkID) continue;

                    if (std::find(electronShowerIDs->begin(), electronShowerIDs->end(), matchedTrkID->at(iRecoTrk)) != electronShowerIDs->end()) {
                        double trk_length = TMath::Sqrt(
                            TMath::Power(recoEndX->at(iRecoTrk) - recoBeginX->at(iRecoTrk), 2) +
                            TMath::Power(recoEndY->at(iRecoTrk) - recoBeginY->at(iRecoTrk), 2) +
                            TMath::Power(recoEndZ->at(iRecoTrk) - recoBeginZ->at(iRecoTrk), 2)
                        );
                        hElectronShowerRecoTrkLengths->Fill(trk_length);

                        hElectronShowerDaughtersRecoIncident->Fill(
                            recoBeginX->at(iRecoTrk) - WC2TPCPrimaryBeginX,
                            recoBeginY->at(iRecoTrk) - WC2TPCPrimaryBeginY
                        );
                    }
                }
            }
        }

        // Look at energy profile for events we are interested in
        if (event == 200920) {
            int caloPoints = wcMatchResR->size();

            TH2D *hDEDXProfile = new TH2D(
                "hDEDXProfile", 
                "hDEDXProfile;Residual range (cm);dE/dx (MeV/cm)",
                caloPoints, 0, 0, 
                caloPoints, 0, 0
            );

            if (event == 200920) {
                std::cout << "Event " << event << " has breakpoint at " << wcMatchResR->at(73) << " rr (cm)" << std::endl; 
            }

            for (int iCalo = 0; iCalo < caloPoints; iCalo++) {
                hDEDXProfile->Fill(wcMatchResR->at(iCalo), wcMatchDEDX->at(iCalo));
            }

            TCanvas* c1 = new TCanvas("c1", "DEDXProfiles", 800, 600);
            hDEDXProfile->SetMinimum(0);
            hDEDXProfile->SetMaximum(1);
            hDEDXProfile->Draw("COLZ");
            Overlay_dEdx_RR_Reference_PP(gProton, gPion, gMuonTG, true, wcMatchResR->at(73));
            c1->SaveAs(SaveDir +  "dEdxProfiles/" + event + ".png");
            
            delete c1;
            delete hDEDXProfile;
        }

        /////////////////////
        // Start selection //
        /////////////////////

        // If no track matched to wire-chamber, skip
        if (WC2TPCtrkID == -99999) continue;

        // Study all reconstructed tracks
        for (int iRecoTrk = 0; iRecoTrk < matchedIdentity->size(); ++iRecoTrk) {
            if (recoTrkID->at(iRecoTrk) == WC2TPCtrkID) continue;

            if (distance(
                    recoBeginX->at(iRecoTrk),
                    truthPrimaryVertexX,
                    recoBeginY->at(iRecoTrk),
                    truthPrimaryVertexY,
                    recoBeginZ->at(iRecoTrk),
                    truthPrimaryVertexZ
                ) < VERTEX_RADIUS) {
                    if (matchedIdentity->at(iRecoTrk) == -211) hRecoSecondaryPions->Fill(matchedKEnergy->at(iRecoTrk));
                    if (matchedIdentity->at(iRecoTrk) == 2212) hRecoSecondaryProtons->Fill(matchedKEnergy->at(iRecoTrk));
                }
        }

        // Study tracks inside 10 cm cylinder around primary track
        int numTracksInCylinder = 0; int numSmallTracksInCylinder = 0;
        if (WC2TPCtrkID != -99999) {
            for (int trk_idx = 0; trk_idx < matchedTrkID->size(); ++trk_idx) {
                if (recoTrkID->at(trk_idx) == WC2TPCtrkID) continue;

                double trackLength = sqrt(
                    pow(recoEndX->at(trk_idx) - recoBeginX->at(trk_idx), 2) +
                    pow(recoEndY->at(trk_idx) - recoBeginY->at(trk_idx), 2) +
                    pow(recoEndZ->at(trk_idx) - recoBeginZ->at(trk_idx), 2)
                );

                // Check if ends are inside cylinder
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

                // Check cylinder
                if (startInCylinder && endInCylinder) {
                    numTracksInCylinder++;

                    hTrkLengthsInCylinder->Fill(trackLength);
                    if (backgroundType == 2) {
                        hTrkLengthsInCylinderMuons->Fill(trackLength);
                    } else if (backgroundType == 3) {
                        hTrkLengthsInCylinderElectrons->Fill(trackLength);
                    } else {
                        hTrkLengthsInCylinderPions->Fill(trackLength);

                        if (backgroundType == 0) {
                            hTrkLengthsInCylinderPionAbs0p->Fill(trackLength);
                        } else if (backgroundType == 1) {
                            hTrkLengthsInCylinderPionAbsNp->Fill(trackLength);
                        } else if (scatteringType == 0 || scatteringType == 1) {
                            hTrkLengthsInCylinderPionScattering->Fill(trackLength);
                        } else if (backgroundType == 7) {
                            hTrkLengthsInCylinderPionChExch->Fill(trackLength);
                            hChExchShowerRecoTrksCylinderContained->Fill(trackLength);
                        } else {
                            hTrkLengthsInCylinderPionOther->Fill(trackLength);
                        }
                    }

                    if (trackLength < CYLINDER_SMALL_TRACK) numSmallTracksInCylinder++;
                } else {
                    if (backgroundType == 7) {
                        hChExchShowerRecoTrksCylinderUnContained->Fill(trackLength);
                    }
                }
            }
            hNumTracksInCylinder->Fill(numTracksInCylinder);
            hSmallTrksInCylinder->Fill(numSmallTracksInCylinder);
            if (backgroundType == 2) {
                hNumTracksInCylinderMuons->Fill(numTracksInCylinder);
                hSmallTrksInCylinderMuons->Fill(numSmallTracksInCylinder);
            } else if (backgroundType == 3) {
                hNumTracksInCylinderElectrons->Fill(numTracksInCylinder);
                hSmallTrksInCylinderElectrons->Fill(numSmallTracksInCylinder);
            } else {
                hNumTracksInCylinderPions->Fill(numTracksInCylinder);
                hSmallTrksInCylinderPions->Fill(numSmallTracksInCylinder);

                if (backgroundType == 0) {
                    hSmallTrksInCylinderPionAbs0p->Fill(numSmallTracksInCylinder);
                } else if (backgroundType == 1) {
                    hSmallTrksInCylinderPionAbsNp->Fill(numSmallTracksInCylinder);
                } else if (scatteringType == 0 || scatteringType == 1) {
                    hSmallTrksInCylinderPionScattering->Fill(numSmallTracksInCylinder);
                } else if (backgroundType == 7) {
                    hSmallTrksInCylinderPionChExch->Fill(numSmallTracksInCylinder);
                } else {
                    hSmallTrksInCylinderPionOther->Fill(numSmallTracksInCylinder);
                }
            }
            
            // Characterize cut performance
            if (numTracksInCylinder <= ALLOWED_CYLINDER_SMALL_TRACKS) {
                if (backgroundType == 2) {
                    numSurvivingMuons++;
                } else if (backgroundType == 3) {
                    numSurvivingElectrons++;
                } else {
                    numSurvivingPions++;
                }
            }
        }

        // Get linearity profile for wire-chamber match
        std::vector<double> WC2TPCLinearity = calcLinearityProfile(*WC2TPCLocationsX, *WC2TPCLocationsY, *WC2TPCLocationsZ, WINDOW_SIZE_LINEARITY);

        if ((event == 19977) || (event == 6197) || (event == 5899) || (event == 5966) || (event == 200213) || (event == 200494) || (event == 54275)) {
            TCanvas* c1 = new TCanvas("c1", "DEDXProfiles", 900, 600);
            c1->SetLeftMargin(0.15);
            c1->SetBottomMargin(0.15);

            int N = WC2TPCLinearity.size();
            std::vector<double> x(N);
            for (int i = 0; i < N; ++i) x[i] = i;

            TGraph *gLinearity = new TGraph(N, x.data(), WC2TPCLinearity.data());
            
            gLinearity->GetYaxis()->SetRangeUser(0.995, 1.001);
            gLinearity->SetTitle("Linearity;Position index;Local linearity");
            gLinearity->SetMarkerStyle(20);
            gLinearity->SetMarkerColor(kBlack);
            gLinearity->Draw("AP");
            c1->SaveAs(SaveDir + "linearityProfiles/" + event + ".png");

            delete c1;
            delete gLinearity;
        }

        // Apply beamline electron shower cut
        // We can do this before incident KE is added since we do not care about electrons
        // if (!passesSmallTracksCut)  continue;
        if (numSmallTracksInCylinder > ALLOWED_CYLINDER_SMALL_TRACKS) continue;

        // Perform chi^2 stitching for primary track
        // For this analysis, we get the following chi^2 values:
        //   - Pion with Bragg peak chi^2
        //   - Proton with Bragg peak chi^2
        //   - MIP chi^2
        //   - Scanning fit with rhs to MIP and lhs to proton

        int totalCaloPoints = wcMatchDEDX->size();
        int nRemoveOutliers = 2;
        int nRemoveEnds     = 3;
        int minPoints       = 5;

        outWCAll << "Event number: " << event << std::endl; 
        outWCAll << "  Calorimetry data points: " << totalCaloPoints << std::endl;

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

        // Get fractional break point
        double fracBreakPoint = (double) bestBreakPoint / totalCaloPoints;
        
        outWCAll << "  Fractional break point: " << (double) bestBreakPoint / totalCaloPoints << std::endl;
        outWCAll << std::endl;
        outWCAll << "  Pion with Bragg peak chi^2: " << pionChi2 << std::endl; 
        outWCAll << "  Proton with Bragg peak chi^2: " << protonChi2 << std::endl; 
        outWCAll << "  MIP chi^2: " << MIPChi2 << std::endl; 
        outWCAll << "  Best break point: " << bestBreakPoint << " with stitched chi^2: " << minStitchedChi2 << std::endl;
        outWCAll << std::endl;
        outWCAll << "  Primary particle PDG: " << truthPrimaryPDG << std::endl;
        outWCAll << "  Event classified as: " << backgroundTypes[backgroundType] << std::endl;
        int numDaughters = truthPrimaryDaughtersPDG->size();
        outWCAll << "  Daughters: " << std::endl;
        for (int i = 0; i < numDaughters; i++) {
            if (truthPrimaryDaughtersPDG->at(i) == 11) continue;
            outWCAll << "    " << truthPrimaryDaughtersPDG->at(i) << " with process: " << truthPrimaryDaughtersProcess->at(i) << std::endl;
        }
        outWCAll << std::endl;

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

        // Before rejecting anything, we want to add energy depositions to incident
        // KE histogram for cross-section computation using the found vertex
        double WCKE             = TMath::Sqrt(WCTrackMomentum * WCTrackMomentum + PionMass * PionMass) - PionMass;
        double calculatedEnLoss = energyLossCalculation(); 
        if (isData) {
            double tanThetaCosPhi = TMath::Tan(WCTheta) * TMath::Cos(WCPhi);
            double tanThetaSinPhi = TMath::Tan(WCTheta) * TMath::Sin(WCPhi);
            double den            = TMath::Sqrt(1 + tanThetaCosPhi * tanThetaCosPhi);
            double onTheFlyPz     = WCTrackMomentum / den;
            double onTheFlyPx     = onTheFlyPz * tanThetaSinPhi;

            calculatedEnLoss = energyLossCalculation(WC4PrimaryX, onTheFlyPx, isData);
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
                if (truthPrimaryPDG == -211) hIncidentKEOnlyPions->Fill(initialKE - energyDeposited);

                // Background breakdown
                if (truthPrimaryPDG == -211) {
                    hIncidentKECorrect->Fill(initialKE - energyDeposited);
                } else if (truthPrimaryPDG == 13) {
                    hIncidentKEMuons->Fill(initialKE - energyDeposited);
                } else if (truthPrimaryPDG == 11) {
                    hIncidentKEElectrons->Fill(initialKE - energyDeposited);
                }
            }
        }
        double energyAtVertex = initialKE - energyDeposited;

        // Check if vertex is inside reduced volume
        if (!isWithinReducedVolume(breakPointX, breakPointY, breakPointZ)) continue;

        double distanceFromVertex         = distance(breakPointX, truthPrimaryVertexX, breakPointY, truthPrimaryVertexY, breakPointZ, truthPrimaryVertexZ);
        double originalDistanceFromVertex = distance(WC2TPCPrimaryEndX, truthPrimaryVertexX, WC2TPCPrimaryEndY, truthPrimaryVertexY, WC2TPCPrimaryEndZ, truthPrimaryVertexZ);
        if (minChi2 == minStitchedChi2) {
            // Track looks the most like a pion + proton, background
            outStitchedFile << "Event number: " << event << std::endl; 
            outStitchedFile << "  Pion with Bragg peak chi^2: " << pionChi2 << std::endl; 
            outStitchedFile << "  Proton with Bragg peak chi^2: " << protonChi2 << std::endl; 
            outStitchedFile << "  MIP chi^2: " << MIPChi2 << std::endl; 
            outStitchedFile << "  Best break point: " << bestBreakPoint << " with stitched chi^2: " << minStitchedChi2 << std::endl;
            outStitchedFile << std::endl;
            outStitchedFile << "  Original point distance from vertex: " << originalDistanceFromVertex << std::endl;
            outStitchedFile << "  Break point distance from vertex: " << distanceFromVertex << std::endl;
            outStitchedFile << "  Is new better than original: " << std::boolalpha << (distanceFromVertex < originalDistanceFromVertex) << std::endl;
            outStitchedFile << std::endl;
            outStitchedFile << "  Primary particle PDG: " << truthPrimaryPDG << std::endl;
            outStitchedFile << "  Event classified as: " << backgroundTypes[backgroundType] << std::endl;
            outStitchedFile << "  Daughters: " << std::endl;
            for (int i = 0; i < numDaughters; i++) {
                if (truthPrimaryDaughtersPDG->at(i) == 11) continue;
                outStitchedFile << "    " << truthPrimaryDaughtersPDG->at(i) << " with process: " << truthPrimaryDaughtersProcess->at(i) << std::endl;
            }
            outStitchedFile << std::endl;

            // Gauge performance for Np events
            if (backgroundType == 1) {
                hStitchedDistanceFromVertex->Fill(distanceFromVertex);
                hStitchedOriginalDistanceFromVertex->Fill(originalDistanceFromVertex);

                if (distanceFromVertex < originalDistanceFromVertex) { betterVertexNp++; }
                else { worseVertexNp++; }
            }
            hStitchAsPionAndProton->Fill(backgroundType);
            hStitchFracBreakPoints->Fill(fracBreakPoint);
            
        } else if (minChi2 == protonChi2) {
            // Track looks the most like a proton, background
            hStitchAsProton->Fill(backgroundType);
        } else if (minChi2 == pionChi2) {
            // Track looks like pion with Bragg peak, background
            hStitchAsPionBraggPeak->Fill(backgroundType);
        } else {
            // Track looks like mip, could be signal
            hStitchAsMIP->Fill(backgroundType);
        }

        ////////////////////////
        // Continue selection //
        ////////////////////////

        // If the particle is categorized as a pion w/ bragg peak:
        //   - Kill
        // If the particle is categorized as proton:
        //   - Kill
        // If the particle is categorized as MIP:
        //   - Chi^2 selection on secondary tracks
        //   - If no secondary tracks, accept as 0p event
        // If the particle is categorized as stitched:
        //   - Chi^2 selection on secondary tracks from new vertex

        // Get reco tracks near vertex
        int secondaryTaggedPion   = 0;
        int secondaryTaggedProton = 0;
        int secondaryTaggedOther  = 0;

        int otherTaggedPion   = 0;
        int otherTaggedProton = 0;

        int numTracksNearVertex = 0;

        int numTracksInSlicedCone      = 0;
        int numSmallTracksInSlicedCone = 0;

        for (size_t trk_idx = 0; trk_idx < recoBeginX->size(); trk_idx++) {
            if (recoTrkID->at(trk_idx) == WC2TPCtrkID) continue;

            int secondaryMatchedPDG = matchedIdentity->at(trk_idx);

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

            // Track is near vertex
            if ((distanceFromStart < VERTEX_RADIUS) || (distanceFromEnd < VERTEX_RADIUS)) {
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

                // Using truth-matched, fill corresponding histogram
                if (secondaryMatchedPDG == -211) {
                    hSecondaryPionChi2Pions->Fill(pionChi2);
                    hSecondaryProtonChi2Pions->Fill(protonChi2);
                } else if ((secondaryMatchedPDG == 2212) && (matchedKEnergy->at(trk_idx) >= PROTON_ENERGY_LOWER_BOUND) && (matchedKEnergy->at(trk_idx) <= PROTON_ENERGY_UPPER_BOUND)) {
                    hSecondaryPionChi2Protons->Fill(pionChi2);
                    hSecondaryProtonChi2Protons->Fill(protonChi2);
                } else {
                    hSecondaryPionChi2Others->Fill(pionChi2);
                    hSecondaryProtonChi2Others->Fill(protonChi2);
                }

                // Fill chi2 plots for optimization
                for (int iPionChiStep = 0; iPionChiStep < pionChiNumSteps; iPionChiStep++) {
                    double currentPionChiValue = pionChiStart + (iPionChiStep * pionChiStep);
                    for (int iProtonChiStep = 0; iProtonChiStep < protonChiNumSteps; iProtonChiStep++) {
                        double currentProtonChiValue = protonChiStart + (iProtonChiStep * protonChiStep);
                        
                        if ((pionChi2 < currentPionChiValue) && (protonChi2 > currentProtonChiValue)) {
                            // Tagged as pion
                            if (secondaryMatchedPDG == -211) {
                                hPionChi2TruePositives->Fill(currentProtonChiValue, currentPionChiValue);
                            } else {
                                hPionChi2FalsePositives->Fill(currentProtonChiValue, currentPionChiValue);
                            }
                        } else if ((pionChi2 > currentPionChiValue) && (protonChi2 < currentProtonChiValue)) {
                            // Tagged as proton
                            if (secondaryMatchedPDG == 2212) {
                                hProtonChi2TruePositives->Fill(currentProtonChiValue, currentPionChiValue);
                            } else {
                                hProtonChi2FalsePositives->Fill(currentProtonChiValue, currentPionChiValue);
                            }
                        }

                        if (!((pionChi2 < currentPionChiValue) && (protonChi2 > currentProtonChiValue))) {
                            if (secondaryMatchedPDG == -211) {
                                hPionChi2FalseNegatives->Fill(currentProtonChiValue, currentPionChiValue);
                            } else {
                                hPionChi2TrueNegatives->Fill(currentProtonChiValue, currentPionChiValue);
                            }
                        }

                        if (!((pionChi2 > currentPionChiValue) && (protonChi2 < currentProtonChiValue))) {
                            if (secondaryMatchedPDG == 2212) {
                                hProtonChi2FalseNegatives->Fill(currentProtonChiValue, currentPionChiValue);
                            } else {
                                hProtonChi2TrueNegatives->Fill(currentProtonChiValue, currentPionChiValue);
                            }
                        }
                    }
                }

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
                    if (secondaryMatchedPDG == -211) {
                        hSecondaryMeanDEDXPions->Fill(secondaryMeanDEDX);
                    } else if ((secondaryMatchedPDG == 2212) && (matchedKEnergy->at(trk_idx) >= PROTON_ENERGY_LOWER_BOUND) && (matchedKEnergy->at(trk_idx) <= PROTON_ENERGY_UPPER_BOUND)) {
                        hSecondaryMeanDEDXProtons->Fill(secondaryMeanDEDX);
                    } else {
                        hSecondaryMeanDEDXOthers->Fill(secondaryMeanDEDX);
                    }

                    if (secondaryMeanDEDX <= MEAN_DEDX_THRESHOLD) {
                        otherTaggedPion++;
                    } else {
                        otherTaggedProton++;
                    }
                }
            } 

            // Apart from checking if track is near vertex, check if it is in ROI for ch exch
            bool startsInCone = IsPointInsideSlicedCone(
                recoBeginX->at(trk_idx), recoBeginY->at(trk_idx), recoBeginZ->at(trk_idx),
                breakPointX, breakPointY, breakPointZ,
                breakPointX - WC2TPCPrimaryBeginX, breakPointY - WC2TPCPrimaryBeginY, breakPointZ - WC2TPCPrimaryBeginZ,
                SLICED_CONE_HEIGHT, SLICED_CONE_MIN_RADIUS, SLICED_CONE_MAX_RADIUS
            );

            bool endsInCone = IsPointInsideSlicedCone(
                recoEndX->at(trk_idx), recoEndY->at(trk_idx), recoEndZ->at(trk_idx),
                breakPointX, breakPointY, breakPointZ,
                breakPointX - WC2TPCPrimaryBeginX, breakPointY - WC2TPCPrimaryBeginY, breakPointZ - WC2TPCPrimaryBeginZ,
                SLICED_CONE_HEIGHT, SLICED_CONE_MIN_RADIUS, SLICED_CONE_MAX_RADIUS
            );

            double trackLength = sqrt(
                pow(recoEndX->at(trk_idx) - recoBeginX->at(trk_idx), 2) +
                pow(recoEndY->at(trk_idx) - recoBeginY->at(trk_idx), 2) +
                pow(recoEndZ->at(trk_idx) - recoBeginZ->at(trk_idx), 2)
            );

            if (startsInCone && endsInCone) {
                numTracksInSlicedCone++;
                if (trackLength < 5) numSmallTracksInSlicedCone++;
                if (backgroundType == 7) hChExchShowerRecoTrksConeContained->Fill(trackLength);

                if (backgroundType == 0) {
                    hSlicedConeTracksPiAbs0p->Fill(trackLength);
                } else if (backgroundType == 1) {
                    hSlicedConeTracksPiAbsNp->Fill(trackLength);
                } else if (scatteringType == 0 || scatteringType == 1) {
                    hSlicedConeTracksPiScat->Fill(trackLength);
                } else if (backgroundType == 7) {
                    hSlicedConeTracksPiChEx->Fill(trackLength);
                } else if (backgroundType == 2) {
                    hSlicedConeTracksMu->Fill(trackLength);
                } else if (backgroundType == 3) {
                    hSlicedConeTracksE->Fill(trackLength);
                } else {
                    hSlicedConeTracksPiOther->Fill(trackLength);
                }
            } else {
                if (backgroundType == 7) hChExchShowerRecoTrksConeUnContained->Fill(trackLength);
            }
        }

        // Fill data for tracks in sliced cone
        if (backgroundType == 0) {
            hSlicedConeNumTrksPiAbs0p->Fill(numTracksInSlicedCone);
            hSlicedConeNumSmallTrksPiAbs0p->Fill(numSmallTracksInSlicedCone);
        } else if (backgroundType == 1) {
            hSlicedConeNumTrksPiAbsNp->Fill(numTracksInSlicedCone);
            hSlicedConeNumSmallTrksPiAbsNp->Fill(numSmallTracksInSlicedCone);
        } else if (scatteringType == 0 || scatteringType == 1) {
            hSlicedConeNumTrksPiScat->Fill(numTracksInSlicedCone);
            hSlicedConeNumSmallTrksPiScat->Fill(numSmallTracksInSlicedCone);
        } else if (backgroundType == 7) {
            hSlicedConeNumTrksPiChEx->Fill(numTracksInSlicedCone);
            hSlicedConeNumSmallTrksPiChEx->Fill(numSmallTracksInSlicedCone);
        } else if (backgroundType == 2) {
            hSlicedConeNumTrksMu->Fill(numTracksInSlicedCone);
            hSlicedConeNumSmallTrksMu->Fill(numSmallTracksInSlicedCone);
        } else if (backgroundType == 3) {
            hSlicedConeNumTrksE->Fill(numTracksInSlicedCone);
            hSlicedConeNumSmallTrksE->Fill(numSmallTracksInSlicedCone);
        } else {
            hSlicedConeNumTrksPiOther->Fill(numTracksInSlicedCone);
            hSlicedConeNumSmallTrksPiOther->Fill(numSmallTracksInSlicedCone);
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

        // If primary is MIP:
        //   - Usual selection
        // If primary is pion with Bragg peak:
        //   - Reject
        // If primary is proton:
        //   - Reject
        // If primary is stitched:
        //   - Usual selection

        if (minChi2 == pionChi2) {
            continue; 
        } else if (minChi2 == protonChi2) {
            continue; // reject events matching to proton primary
        }

        int totalTaggedPions   = secondaryTaggedPion + otherTaggedPion;
        int totalTaggedProtons = secondaryTaggedProton + otherTaggedProton;

        if (totalTaggedPions > 0) {
            if (totalTaggedPions > 1) continue; // reject events with more than one tagged pion
            
            // If any secondary tagged as pion, event is not pion absorption signal
            // However, we are still interested in tagging scattering events
            hPionScatteringKE->Fill(energyAtVertex);

            // Background breakdown for pion scattering
            if (scatteringType != -1) hPionScatteringKECorrect->Fill(energyAtVertex);
            else if (backgroundType == 0) hPionScatteringKE0pAbs->Fill(energyAtVertex);
            else if (backgroundType == 1) hPionScatteringKENpAbs->Fill(energyAtVertex);
            else if (backgroundType == 7) hPionScatteringKEChEx->Fill(energyAtVertex);
            else if (backgroundType == 2) hPionScatteringKEMuons->Fill(energyAtVertex);
            else if (backgroundType == 3) hPionScatteringKEElectrons->Fill(energyAtVertex);
            else hPionScatteringKEOther->Fill(energyAtVertex);

            if (totalTaggedProtons == 0) {
                hPionScatteringKE0p->Fill(energyAtVertex);

                // Background breakdown for pion 0p scattering
                if (scatteringType == 0) h0pPionScatteringKECorrect->Fill(energyAtVertex);
                else if (scatteringType == 1) h0pPionScatteringKENpScatter->Fill(energyAtVertex);
                else if (backgroundType == 0) h0pPionScatteringKE0pAbs->Fill(energyAtVertex);
                else if (backgroundType == 1) h0pPionScatteringKENpAbs->Fill(energyAtVertex);
                else if (backgroundType == 7) h0pPionScatteringKEChEx->Fill(energyAtVertex);
                else if (backgroundType == 2) h0pPionScatteringKEMuons->Fill(energyAtVertex);
                else if (backgroundType == 3) h0pPionScatteringKEElectrons->Fill(energyAtVertex);
                else h0pPionScatteringKEOther->Fill(energyAtVertex);
            } else if (totalTaggedProtons > 0) {
                hPionScatteringKENp->Fill(energyAtVertex);

                // Background breakdown for pion Np scattering
                if (scatteringType == 1) hNpPionScatteringKECorrect->Fill(energyAtVertex);
                else if (scatteringType == 0) hNpPionScatteringKE0pScatter->Fill(energyAtVertex);
                else if (backgroundType == 0) hNpPionScatteringKE0pAbs->Fill(energyAtVertex);
                else if (backgroundType == 1) hNpPionScatteringKENpAbs->Fill(energyAtVertex);
                else if (backgroundType == 7) hNpPionScatteringKEChEx->Fill(energyAtVertex);
                else if (backgroundType == 2) hNpPionScatteringKEMuons->Fill(energyAtVertex);
                else if (backgroundType == 3) hNpPionScatteringKEElectrons->Fill(energyAtVertex);
                else hNpPionScatteringKEOther->Fill(energyAtVertex);
            }
            continue;
        }

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

            float hitX     = fHitX->at(iHit);
            float hitW     = fHitW->at(iHit);
            int   hitPlane = fHitPlane->at(iHit);

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

        // std::cout << "Found clusters: " << hitClusters.size() << std::endl;
        // std::cout << std::endl;

        // Optimize cut; note, this cut is only for events so far categorized as 0p
        if (totalTaggedProtons == 0) {
            double currentClusterSize = STARTING_CLUSTER_SIZE;
            for (int clusterStep = 0; clusterStep <= CLUSTER_SIZE_NUM_STEPS; ++clusterStep) {                
                int largeClusterCounter = 0;
                for (int i = 0; i < hitClusters.size(); ++i) {
                    double maxWire = *std::max_element(hitClusters[i].hitW.begin(), hitClusters[i].hitW.end());
                    double minWire = *std::min_element(hitClusters[i].hitW.begin(), hitClusters[i].hitW.end());
                    if (std::abs(maxWire - minWire) > currentClusterSize) largeClusterCounter++;
                }

                int currentNumClusters = STARTING_NUM_CLUSTERS;
                for (int numClusterStep = 0; numClusterStep <= NUM_CLUSTERS_NUM_STEPS; ++numClusterStep) {                    
                    if (largeClusterCounter < currentNumClusters) {
                        // ID event as 0p
                        hHitClusterCut0pReco->Fill(currentClusterSize, currentNumClusters);
                        if (backgroundType == 0) hHitClusterCut0pRecoTrue->Fill(currentClusterSize, currentNumClusters);
                    }
                    currentNumClusters += NUM_CLUSTERS_STEP;
                }
                currentClusterSize += CLUSTER_SIZE_STEP;
            }
        }

        // Get data for cut
        float totalClusterSize = 0.;
        int   numLargeClusters = 0;
        for (int i = 0; i < hitClusters.size(); ++i) {
            float maxWire = *std::max_element(hitClusters[i].hitW.begin(), hitClusters[i].hitW.end());
            float minWire = *std::min_element(hitClusters[i].hitW.begin(), hitClusters[i].hitW.end());
            if (std::abs(maxWire - minWire) > LARGE_CLUSTER_THRESHOLD) numLargeClusters++;
            totalClusterSize += maxWire - minWire;
        }
        float averageClusterSize = totalClusterSize / ((float) hitClusters.size());

        // Fill histograms
        if (totalTaggedProtons == 0) {
            if (backgroundType == 0) {
                hHitClusters0p->Fill(hitClusters.size());
                hHitLargeClusters0p->Fill(numLargeClusters);
                hHitClustersSize0p->Fill(averageClusterSize);
            } else {
                hHitClusters0pBackground->Fill(hitClusters.size());
                hHitLargeClusters0pBackground->Fill(numLargeClusters);
                hHitClustersSize0pBackground->Fill(averageClusterSize);
            }
        } else if (totalTaggedProtons > 0) {
            if (backgroundType == 1) {
                hHitClustersNp->Fill(hitClusters.size());
                hHitLargeClustersNp->Fill(numLargeClusters);
                hHitClustersSizeNp->Fill(averageClusterSize);
            } else {
                hHitClustersNpBackground->Fill(hitClusters.size());
                hHitLargeClustersNpBackground->Fill(numLargeClusters);
                hHitClustersSizeNpBackground->Fill(averageClusterSize);
            }
        }

        // If tagged as 0p and too many large clusters, reject
        if ((totalTaggedProtons == 0) && (numLargeClusters >= NUM_CLUSTERS_THRESHOLD)) continue;

        ////////////////////////////////////////////////////
        // Study local spatial linearity of primary track //
        ////////////////////////////////////////////////////

        // If we did not detect a kink with the energy profile, we look
        // for a spatial kink using the local linearity of the primary track
        double maxLocalLinearityD, maxLocalLinearityDD, minLocalLinearity, stdevLocalLinearity;

        if (minChi2 == MIPChi2) {
            minLocalLinearity          = *std::min_element(WC2TPCLinearity.begin(), WC2TPCLinearity.end());
            double sumLocalLinearity   = std::accumulate(WC2TPCLinearity.begin(), WC2TPCLinearity.end(), 0.0);
            double meanLocalLinearity  = sumLocalLinearity / WC2TPCLinearity.size();
            double sqsumLocalLinearity = std::inner_product(WC2TPCLinearity.begin(), WC2TPCLinearity.end(), WC2TPCLinearity.begin(), 0.0);
            stdevLocalLinearity        = std::sqrt(sqsumLocalLinearity / WC2TPCLinearity.size() - meanLocalLinearity * meanLocalLinearity);

            std::vector<double> localLinearityD;
            for (int i = 1; i < WC2TPCLinearity.size(); ++i) {
                double diff = WC2TPCLinearity.at(i) - WC2TPCLinearity.at(i - 1);
                localLinearityD.push_back(diff);
            }

            std::vector<double> localLinearityDAbs;
            for (int i = 0; i < localLinearityD.size(); ++i) localLinearityDAbs.push_back(std::abs(localLinearityD[i]));
            maxLocalLinearityD = *std::max_element(localLinearityDAbs.begin(), localLinearityDAbs.end());

            std::vector<double> localLinearityDD;
            for (int i = 1; i < localLinearityD.size(); ++i) {
                double diff = localLinearityD.at(i) - localLinearityD.at(i - 1);
                localLinearityDD.push_back(std::abs(diff));
            }
            maxLocalLinearityDD = *std::max_element(localLinearityDD.begin(), localLinearityDD.end());

            if (totalTaggedProtons == 0) {
                if (backgroundType == 0) {
                    hMinimumLinearity0p->Fill(minLocalLinearity);
                    hStdDevLinearity0p->Fill(stdevLocalLinearity);
                    hMaxLinearityD0p->Fill(maxLocalLinearityD);
                    hMaxLinearityDD0p->Fill(maxLocalLinearityDD);
                } else {
                    hMinimumLinearity0pBackground->Fill(minLocalLinearity);
                    hStdDevLinearity0pBackground->Fill(stdevLocalLinearity);
                    hMaxLinearityD0pBackground->Fill(maxLocalLinearityD);
                    hMaxLinearityDD0pBackground->Fill(maxLocalLinearityDD);

                    // if (backgroundType == 6){ hMaxLinearityD0pScattering->Fill(maxLocalLinearityD); }
                    // else { hMaxLinearityD0pBackground->Fill(maxLocalLinearityD); }
                }
            } else if (totalTaggedProtons > 0) {
                if (backgroundType == 1) {
                    hMinimumLinearityNp->Fill(minLocalLinearity);
                    hStdDevLinearityNp->Fill(stdevLocalLinearity);
                } else {
                    hMinimumLinearityNpBackground->Fill(minLocalLinearity);
                    hStdDevLinearityNpBackground->Fill(stdevLocalLinearity);
                }
            }

            // Study efficiency and purity curves
            double currentLinDerivThreshold = INIT_LIN_DERIV_THRESHOLD;
            for (int iLinDeriv = 0; iLinDeriv < LIN_DERIV_NUM_STEPS; ++iLinDeriv) {
                if ((totalTaggedProtons == 0) && (maxLocalLinearityD < currentLinDerivThreshold)) {
                    hLocalLinearityDerivativeReco->Fill(currentLinDerivThreshold + 0.01e-3);
                    if (backgroundType == 0) hLocalLinearityDerivativeRecoTrue->Fill(currentLinDerivThreshold + 0.01e-3);
                }
                currentLinDerivThreshold += LIN_DERIV_STEP_SIZE;
            }

            double currentLinStdDevThreshold = INIT_LIN_STDDEV_THRESHOLD;
            for (int iLinStdDev = 0; iLinStdDev < LIN_STDDEV_NUM_STEPS; ++iLinStdDev) {
                if ((totalTaggedProtons == 0) && (stdevLocalLinearity < currentLinStdDevThreshold)) {
                    hLocalLinearityStdDevReco->Fill(currentLinStdDevThreshold + 0.001e-3);
                    if (backgroundType == 0) hLocalLinearityStdDevRecoTrue->Fill(currentLinStdDevThreshold + 0.001e-3);
                }
                currentLinStdDevThreshold += LIN_STDDEV_STEP_SIZE;
            }

            double currentLinMinThreshold = INIT_LIN_MIN_THRESHOLD;
            for (int iLinMin = 0; iLinMin < LIN_MIN_NUM_STEPS; ++iLinMin) {
                if ((totalTaggedProtons == 0) && (minLocalLinearity > currentLinMinThreshold)) {
                    hLocalLinearityMinReco->Fill((1 - currentLinMinThreshold) + 0.0001e-3);
                    if (backgroundType == 0) hLocalLinearityMinRecoTrue->Fill((1 - currentLinMinThreshold) + 0.0001e-3);
                }
                currentLinMinThreshold += LIN_MIN_STEP_SIZE;
            }

            // Actually perform cut 
            // if ((totalTaggedProtons == 0) && (maxLocalLinearityD >= LINEARITY_DERIVATIVE_THRESHOLD)) continue;
            // if ((totalTaggedProtons == 0) && (stdevLocalLinearity >= STD_DEV_LINEARITY_THRESHOLD)) continue;
        }

        ///////////////////////////
        // Fill final histograms //
        ///////////////////////////

        hRecoAbsorption->Fill(backgroundType);
        if (totalTaggedProtons == 0) { 
            hRecoAbsorption0p->Fill(backgroundType);
        } else if (totalTaggedProtons > 0) {
            hRecoAbsorptionNp->Fill(backgroundType);
        }

        hPionAbsKE->Fill(energyAtVertex);
        if (totalTaggedProtons == 0) {
            h0pPionAbsKE->Fill(energyAtVertex);
        } else if (totalTaggedProtons > 0) {
            hNpPionAbsKE->Fill(energyAtVertex);
        }

        // We want to see if we filled the interacting bin correctly using truth info
        if (backgroundType == 0 || backgroundType == 1) hPionAbsKECorrect->Fill(energyAtVertex);
        else if (scatteringType == 0) hPionAbsKE0pScattering->Fill(energyAtVertex);
        else if (scatteringType == 1) hPionAbsKENpScattering->Fill(energyAtVertex);
        else if (backgroundType == 7) hPionAbsKEChEx->Fill(energyAtVertex);
        else if (backgroundType == 2) hPionAbsKEMuons->Fill(energyAtVertex);
        else if (backgroundType == 3) hPionAbsKEElectrons->Fill(energyAtVertex);
        else hPionAbsKEOther->Fill(energyAtVertex);

        // Same thing for 0p absorption
        if (totalTaggedProtons == 0) {
            if (backgroundType == 0) h0pPionAbsKECorrect->Fill(energyAtVertex);
            else if (backgroundType == 1) h0pPionAbsKENpAbs->Fill(energyAtVertex);
            else if (scatteringType == 0) h0pPionAbsKE0pScattering->Fill(energyAtVertex);
            else if (scatteringType == 1) h0pPionAbsKENpScattering->Fill(energyAtVertex);
            else if (backgroundType == 7) h0pPionAbsKEChEx->Fill(energyAtVertex);
            else if (backgroundType == 2) h0pPionAbsKEMuons->Fill(energyAtVertex);
            else if (backgroundType == 3) h0pPionAbsKEElectrons->Fill(energyAtVertex);
            else h0pPionAbsKEOther->Fill(energyAtVertex);
        } 
        // And for Np absorption
        else if (totalTaggedProtons > 0) {
            if (backgroundType == 1) hNpPionAbsKECorrect->Fill(energyAtVertex);
            else if (backgroundType == 0) hNpPionAbsKE0pAbs->Fill(energyAtVertex);
            else if (scatteringType == 0) hNpPionAbsKE0pScattering->Fill(energyAtVertex);
            else if (scatteringType == 1) hNpPionAbsKENpScattering->Fill(energyAtVertex);
            else if (backgroundType == 7) hNpPionAbsKEChEx->Fill(energyAtVertex);
            else if (backgroundType == 2) hNpPionAbsKEMuons->Fill(energyAtVertex);
            else if (backgroundType == 3) hNpPionAbsKEElectrons->Fill(energyAtVertex);
            else hNpPionAbsKEOther->Fill(energyAtVertex);
        }

        if (backgroundType == 0 || backgroundType == 1) {
            hPionAbsKEOnlyAbs->Fill(energyAtVertex);
            if (backgroundType == 0 && totalTaggedProtons == 0) {
                h0pPionAbsKEOnly0pAbs->Fill(energyAtVertex);
            } else if (backgroundType == 1 && totalTaggedProtons > 0) {
                hNpPionAbsKEOnlyNpAbs->Fill(energyAtVertex);
            }
        }

        /////////////////////////
        // Analyze backgrounds //
        /////////////////////////

        // At this point, any event still alive was tagged as pion absorption
        
        // First, just print all the data for the background events
        EventInfo thisEventInfo;
        thisEventInfo.run    = run;
        thisEventInfo.subrun = subrun;
        thisEventInfo.event  = event;

        thisEventInfo.vertexX = truthPrimaryVertexX;
        thisEventInfo.vertexY = truthPrimaryVertexY;
        thisEventInfo.vertexZ = truthPrimaryVertexZ;

        thisEventInfo.wcMatchPDG              = wcMatchPDG;
        thisEventInfo.wcMatchProcess          = *wcMatchProcess;
        thisEventInfo.wcMatchDaughtersPDG     = *wcMatchDaughtersPDG;
        thisEventInfo.wcMatchDaughtersProcess = *wcMatchDaughtersProcess;

        thisEventInfo.truthPrimaryPDG              = truthPrimaryPDG;
        thisEventInfo.truthPrimaryDaughtersPDG     = *truthPrimaryDaughtersPDG;
        thisEventInfo.truthPrimaryDaughtersProcess = *truthPrimaryDaughtersProcess;

        thisEventInfo.visibleProtons   = numVisibleProtons;
        thisEventInfo.backgroundNum    = backgroundType;
        thisEventInfo.recoProtonCount  = totalTaggedProtons;
        thisEventInfo.tracksNearVertex = numTracksNearVertex;

        thisEventInfo.truthSecondaryDaughtersPDG     = *truthSecondaryPionDaughtersPDG;
        thisEventInfo.truthSecondaryDaughtersProcess = *truthSecondaryPionDaughtersProcess;
        thisEventInfo.truthSecondaryDaughtersKE      = *truthSecondaryPionDaughtersKE;

        int secondaryInteractionTag           = isSecondaryInteractionAbsorption(*truthSecondaryPionDaughtersPDG, *truthSecondaryPionDaughtersProcess, *truthSecondaryPionDaughtersKE);
        thisEventInfo.secondaryInteractionTag = secondaryInteractionTag;

        thisEventInfo.numHitClusters     = hitClusters.size();
        thisEventInfo.minLocalLinearity  = minLocalLinearity;
        thisEventInfo.maxLocalLinearityD = maxLocalLinearityD;

        // Get info about scattering events that survive
        if (backgroundType == 6) {
            hTotalBackgroundInelScatteringAngle->Fill(truthScatteringAngle);
            hTotalBackgroundInelScatteringLength->Fill(truthScatteredPionLength);
            hTotalBackgroundInelScatteringVertexKE->Fill(truthPrimaryVertexKE);

            hTotalBackgroundInelScatteringLengthVSAngle->Fill(truthScatteredPionLength, truthScatteringAngle);
            hTotalBackgroundInelScatteringLengthVSVertexKE->Fill(truthScatteredPionLength, truthPrimaryVertexKE);
            hTotalBackgroundInelScatteringVertexKEVSAngle->Fill(truthPrimaryVertexKE, truthScatteringAngle);

            hTotalBackgroundInelasticScatteringEnergyDiff->Fill(truthPrimaryVertexKE - truthScatteredPionKE);

            // For all scattering
            hTotalBackgroundAllScatteringAngle->Fill(truthScatteringAngle);
            hTotalBackgroundAllScatteringVertexKE->Fill(truthPrimaryVertexKE);
            hTotalBackgroundAllScatteringVertexKEVSAngle->Fill(truthPrimaryVertexKE, truthScatteringAngle);

            // Compute distance of reco'ed vertex to secondary vertex
            double distanceFromSecondaryVertex         = distance(breakPointX, truthSecondaryVertexX, breakPointY, truthSecondaryVertexY, breakPointZ, truthSecondaryVertexZ);
            double originalDistanceFromSecondaryVertex = distance(WC2TPCPrimaryEndX, truthSecondaryVertexX, WC2TPCPrimaryEndY, truthSecondaryVertexY, WC2TPCPrimaryEndZ, truthSecondaryVertexZ);

            hInelScatteringBackgroundOriginalDistanceToPrimaryVertex->Fill(originalDistanceFromVertex);
            hInelScatteringBackgroundOriginalDistanceToSecondaryVertex->Fill(originalDistanceFromSecondaryVertex);
            hInelScatteringBackgroundDistanceToPrimaryVertex->Fill(distanceFromVertex);
            hInelScatteringBackgroundDistanceToSecondaryVertex->Fill(distanceFromSecondaryVertex);

            if (totalTaggedProtons == 0) {
                Background0pInelasticScatteringIncidentKE.push_back(truthPrimaryIncidentKE);
                Background0pInelasticScatteringVertexKE.push_back(truthPrimaryVertexKE);
                Background0pInelasticScatteringOutgoingKE.push_back(truthScatteredPionKE);
                
                h0pBackgroundInelScatteringVertexKEVSAngle->Fill(truthPrimaryVertexKE, truthScatteringAngle);
                h0pBackgroundInelasticScatteringEnergyDiff->Fill(truthPrimaryVertexKE - truthScatteredPionKE);

                hInelasticScatteringTotal0pBkg->Fill(truthScatteredPionKE);

                // If reco'ed
                for (int iRecoTrk = 0; iRecoTrk < matchedIdentity->size(); ++iRecoTrk) {
                    if (recoTrkID->at(iRecoTrk) == WC2TPCtrkID) continue;

                    if ((matchedIdentity->at(iRecoTrk) == -211) && (matchedTrkID->at(iRecoTrk) != WC2TPCtrkID) && (matchedProcess->at(iRecoTrk) == "pi-Inelastic")) {
                        hInelasticScatteringReconstructed0pBkg->Fill(truthScatteredPionKE);
                        break;
                    }
                }

                // For all scattering
                h0pBackgroundAllScatteringAngle->Fill(truthScatteringAngle);
                h0pBackgroundAllScatteringVertexKE->Fill(truthPrimaryVertexKE);
                h0pBackgroundAllScatteringVertexKEVSAngle->Fill(truthPrimaryVertexKE, truthScatteringAngle);
            } else if (totalTaggedProtons > 0) {
                BackgroundNpInelasticScatteringIncidentKE.push_back(truthPrimaryIncidentKE);
                BackgroundNpInelasticScatteringVertexKE.push_back(truthPrimaryVertexKE);
                BackgroundNpInelasticScatteringOutgoingKE.push_back(truthScatteredPionKE);

                hNpBackgroundInelScatteringVertexKEVSAngle->Fill(truthPrimaryVertexKE, truthScatteringAngle);
                hNpBackgroundInelasticScatteringEnergyDiff->Fill(truthPrimaryVertexKE - truthScatteredPionKE);

                hInelasticScatteringTotalNpBkg->Fill(truthScatteredPionKE);

                // If reco'ed
                for (int iRecoTrk = 0; iRecoTrk < matchedIdentity->size(); ++iRecoTrk) {
                    if (recoTrkID->at(iRecoTrk) == WC2TPCtrkID) continue;

                    if ((matchedIdentity->at(iRecoTrk) == -211) && (matchedTrkID->at(iRecoTrk) != WC2TPCtrkID) && (matchedProcess->at(iRecoTrk) == "pi-Inelastic")) {
                        hInelasticScatteringReconstructedNpBkg->Fill(truthScatteredPionKE);
                        break;
                    }
                }

                // For all scattering
                hNpBackgroundAllScatteringAngle->Fill(truthScatteringAngle);
                hNpBackgroundAllScatteringVertexKE->Fill(truthPrimaryVertexKE);
                hNpBackgroundAllScatteringVertexKEVSAngle->Fill(truthPrimaryVertexKE, truthScatteringAngle);
            }
        } else if (backgroundType == 12) {
            hTotalBackgroundElasticScatteringAngle->Fill(trajectoryInteractionAngle);
            hTotalBackgroundElasticScatteringVertexKE->Fill(trajectoryInteractionKE);
            hTotalBackgroundElasticScatteringVertexKEVSAngle->Fill(trajectoryInteractionKE, trajectoryInteractionAngle);
            
            // For all scattering
            hTotalBackgroundAllScatteringAngle->Fill(trajectoryInteractionAngle);
            hTotalBackgroundAllScatteringVertexKE->Fill(trajectoryInteractionKE);
            hTotalBackgroundAllScatteringVertexKEVSAngle->Fill(trajectoryInteractionKE, trajectoryInteractionAngle);

            if (totalTaggedProtons == 0) {
                h0pBackgroundElasticScatteringAngle->Fill(trajectoryInteractionAngle);
                h0pBackgroundElasticScatteringVertexKE->Fill(trajectoryInteractionKE);
                h0pBackgroundElasticScatteringVertexKEVSAngle->Fill(trajectoryInteractionKE, trajectoryInteractionAngle);

                // For all scattering
                h0pBackgroundAllScatteringAngle->Fill(trajectoryInteractionAngle);
                h0pBackgroundAllScatteringVertexKE->Fill(trajectoryInteractionKE);
                h0pBackgroundAllScatteringVertexKEVSAngle->Fill(trajectoryInteractionKE, trajectoryInteractionAngle);
            } else if (totalTaggedProtons > 0) {
                hNpBackgroundElasticScatteringAngle->Fill(trajectoryInteractionAngle);
                hNpBackgroundElasticScatteringVertexKE->Fill(trajectoryInteractionKE);
                hNpBackgroundElasticScatteringVertexKEVSAngle->Fill(trajectoryInteractionKE, trajectoryInteractionAngle);

                // For all scattering
                hNpBackgroundAllScatteringAngle->Fill(trajectoryInteractionAngle);
                hNpBackgroundAllScatteringVertexKE->Fill(trajectoryInteractionKE);
                hNpBackgroundAllScatteringVertexKEVSAngle->Fill(trajectoryInteractionKE, trajectoryInteractionAngle);
            }
        }

        if (totalTaggedProtons == 0) {
            if (backgroundType == 0) {
                printEventInfo(thisEventInfo, outFile0pTrue);
            } else if (backgroundType == 6) {
                printEventInfo(thisEventInfo, outFile0pBackground);
                h0pInelasticBackgroundSecondaryInteraction->Fill(secondaryInteractionTag);
                h0pBackgroundInelScatteringAngle->Fill(truthScatteringAngle);
                h0pBackgroundInelScatteringLength->Fill(truthScatteredPionLength);
                h0pBackgroundInelScatteringVertexKE->Fill(truthPrimaryVertexKE);
            } else { printEventInfo(thisEventInfo, outFile0pBackground); }
        } else if (totalTaggedProtons > 0) {
            if (backgroundType == 1) {
                printEventInfo(thisEventInfo, outFileNpTrue);
            } else if (backgroundType == 6) {
                printEventInfo(thisEventInfo, outFileNpBackground);
                hNpInelasticBackgroundSecondaryInteraction->Fill(secondaryInteractionTag);
                hNpBackgroundInelScatteringAngle->Fill(truthScatteringAngle);
                hNpBackgroundInelScatteringLength->Fill(truthScatteredPionLength);
                hNpBackgroundInelScatteringVertexKE->Fill(truthPrimaryVertexKE);
            } else { printEventInfo(thisEventInfo, outFileNpBackground); }
        }
    }

    std::cout << std::endl;
    std::cout << "Total events: " << hTotalEvents->Integral() << std::endl;
    printBackgroundInfo(hTotalEvents, std::cout);
    std::cout << std::endl;
    std::cout << "Total absorption reco'ed events: " << hRecoAbsorption->Integral() << std::endl;
    std::cout << "  Abs overall purity: " << hRecoAbsorption->Integral(1,2) / hRecoAbsorption->Integral() << std::endl;
    std::cout << "  Abs overall efficiency: " << hRecoAbsorption->Integral(1,2) / hTotalEvents->Integral(1,2) << std::endl;
    std::cout << std::endl;
    printBackgroundInfo(hRecoAbsorption, std::cout);
    std::cout << std::endl;
    std::cout << "Total 0p absorption reco'ed events: " << hRecoAbsorption0p->Integral() << std::endl;
    std::cout << "  Abs 0p purity: " << hRecoAbsorption0p->Integral(1,1) / hRecoAbsorption0p->Integral() << std::endl;
    std::cout << "  Abs 0p efficiency: " << hRecoAbsorption0p->Integral(1,1) / hTotalEvents->Integral(1,1) << std::endl;
    std::cout << std::endl;
    printBackgroundInfo(hRecoAbsorption0p, std::cout);
    std::cout << std::endl;
    std::cout << "0p inelastic scattering background secondary interactions: " << std::endl;
    printBackgroundInfo(h0pInelasticBackgroundSecondaryInteraction, std::cout);
    std::cout << std::endl;
    std::cout << "Total Np absorption reco'ed events: " << hRecoAbsorptionNp->Integral() << std::endl;
    std::cout << "  Abs Np purity: " << hRecoAbsorptionNp->Integral(2,2) / hRecoAbsorptionNp->Integral() << std::endl;
    std::cout << "  Abs Np efficiency: " << hRecoAbsorptionNp->Integral(2,2) / hTotalEvents->Integral(2,2) << std::endl;
    std::cout << std::endl;
    printBackgroundInfo(hRecoAbsorptionNp, std::cout);
    std::cout << std::endl;
    std::cout << "Np inelastic scattering background secondary interactions: " << std::endl;
    printBackgroundInfo(hNpInelasticBackgroundSecondaryInteraction, std::cout);
    std::cout << std::endl;

    ///////////////////////////
    // Compute FOM for chi^2 //
    ///////////////////////////

    TH2D* hProtonChi2FOM = new TH2D(
        "hProtonChi2FOM", 
        "hProtonChi2FOM;Proton chi squared;Pion chi squared", 
        protonChiNumSteps, protonChiStart, protonChiEnd, 
        pionChiNumSteps, pionChiStart, pionChiEnd
    );
    TH2D* hPionChi2FOM = new TH2D(
        "hPionChi2FOM", 
        "hPionChi2FOM;Proton chi squared;Pion chi squared", 
        protonChiNumSteps, protonChiStart, protonChiEnd, 
        pionChiNumSteps, pionChiStart, pionChiEnd
    );

    for (int iPionChiStep = 0; iPionChiStep < pionChiNumSteps; iPionChiStep++) {
        double currentPionChiValue = pionChiStart + (iPionChiStep * pionChiStep);
        for (int iProtonChiStep = 0; iProtonChiStep < protonChiNumSteps; iProtonChiStep++) {
            double currentProtonChiValue = protonChiStart + (iProtonChiStep * protonChiStep);

            int binX = hPionChi2TruePositives->GetXaxis()->FindFixBin(currentProtonChiValue);
            int binY = hPionChi2TruePositives->GetYaxis()->FindFixBin(currentPionChiValue);

            double pion_tp     = hPionChi2TruePositives->GetBinContent(binX, binY);
            double pion_tn     = hPionChi2TrueNegatives->GetBinContent(binX, binY);
            double pion_fp     = hPionChi2FalsePositives->GetBinContent(binX, binY);
            double pion_fn     = hPionChi2FalseNegatives->GetBinContent(binX, binY);
            double pion_eff    = pion_tp / (pion_tp + pion_tn);
            double pion_purity = pion_tp / (pion_tp + pion_fp);

            // double pion_fom = pion_tp / std::sqrt(pion_tp + pion_fp);
            double pion_fom = (2 * pion_tp) / (2 * pion_tp + pion_fp + pion_fn);
            // double pion_fom = ((pion_tp * pion_tn) - (pion_fp * pion_fn)) / std::sqrt((pion_tp + pion_fp) * (pion_tp + pion_fn) * (pion_tn * pion_fp) * (pion_tn * pion_fn));
            hPionChi2FOM->SetBinContent(binX, binY, pion_fom);

            double proton_tp     = hProtonChi2TruePositives->GetBinContent(binX, binY);
            double proton_tn     = hProtonChi2TrueNegatives->GetBinContent(binX, binY);
            double proton_fp     = hProtonChi2FalsePositives->GetBinContent(binX, binY);
            double proton_fn     = hProtonChi2FalseNegatives->GetBinContent(binX, binY);
            double proton_eff    = proton_tp / (proton_tp + proton_tn);
            double proton_purity = proton_tp / (proton_tp + proton_fp);

            // double proton_fom = proton_tp / std::sqrt(proton_tp + proton_fp);
            double proton_fom = (2 * proton_tp) / (2 * proton_tp + proton_fp + proton_fn);
            // double proton_fom = ((proton_tp * proton_tn) - (proton_fp * proton_fn)) / std::sqrt((proton_tp + proton_fp) * (proton_tp + proton_fn) * (proton_tn * proton_fp) * (proton_tn * proton_fn));
            hProtonChi2FOM->SetBinContent(binX, binY, proton_fom);
        }
    }

    Int_t binxPion, binyPion, binzPion;
    Int_t globalBinPion   = hPionChi2FOM->GetMaximumBin(binxPion, binyPion, binzPion);

    double maxContentPion = hPionChi2FOM->GetBinContent(globalBinPion);
    double xAtMaxPion     = hPionChi2FOM->GetXaxis()->GetBinCenter(binxPion);
    double yAtMaxPion     = hPionChi2FOM->GetYaxis()->GetBinCenter(binyPion);

    std::cout << std::endl;
    std::cout << "Pion χ² max‑FOM = " << maxContentPion << " at (χ²_p, χ²_π) = (" << xAtMaxPion << ", " << yAtMaxPion << ")" << std::endl;

    Int_t binxProt, binyProt, binzProt;
    Int_t globalBinProt   = hProtonChi2FOM->GetMaximumBin(binxProt, binyProt, binzProt);

    double maxContentProt = hProtonChi2FOM->GetBinContent(globalBinProt);
    double xAtMaxProt     = hProtonChi2FOM->GetXaxis()->GetBinCenter(binxProt);
    double yAtMaxProt     = hProtonChi2FOM->GetYaxis()->GetBinCenter(binyProt);

    std::cout << "Proton χ² max‑FOM = " << maxContentProt << " at (χ²_p, χ²_π) = (" << xAtMaxProt << ", " << yAtMaxProt << ")" << std::endl;
    std::cout << std::endl;

    ///////////////////////////////////////
    // Compute purity for hit clustering //
    ///////////////////////////////////////

    TH2D* hHitClusterCut0pPurity = (TH2D*) hHitClusterCut0pRecoTrue->Clone("hHitClusterCut0pPurity");
    hHitClusterCut0pPurity->SetTitle("0p purity;Large cluster size;Number of large clusters");
    hHitClusterCut0pPurity->Divide(hHitClusterCut0pReco);

    Int_t hitClusterBinX, hitClusterBinY, hitClusterBinZ;
    Int_t globalHitClusterBin = hHitClusterCut0pPurity->GetMaximumBin(hitClusterBinX, hitClusterBinY, hitClusterBinZ);

    double maxContentHitCluster = hHitClusterCut0pPurity->GetBinContent(globalHitClusterBin);
    double xAtMaxHitCluster     = hHitClusterCut0pPurity->GetXaxis()->GetBinLowEdge(hitClusterBinX);
    double yAtMaxHitCluster     = hHitClusterCut0pPurity->GetYaxis()->GetBinLowEdge(hitClusterBinY);

    std::cout << "Hit cluster cut max purity = " << maxContentHitCluster << " at (cluster size, num clusters) = (" << xAtMaxHitCluster << ", " << yAtMaxHitCluster << ")" << std::endl;
    std::cout << std::endl;

    //////////////////////////////////////////////////////////////////
    // Compute purity and efficiency for local linearity derivative //
    //////////////////////////////////////////////////////////////////

    hLocalLinearityDerivativeCutPur = (TH1D*) hLocalLinearityDerivativeRecoTrue->Clone("hLocalLinearityDerivativeCutPur");
    hLocalLinearityDerivativeCutPur->Divide(hLocalLinearityDerivativeReco);

    hLocalLinearityDerivativeCutEff = (TH1D*) hLocalLinearityDerivativeRecoTrue->Clone("hLocalLinearityDerivativeCutPur");
    hLocalLinearityDerivativeCutEff->Scale(1.0 / hTotalEvents->Integral(2,2));

    hLocalLinearityStdDevCutPur = (TH1D*) hLocalLinearityStdDevRecoTrue->Clone("hLocalLinearityStdDevCutPur");
    hLocalLinearityStdDevCutPur->Divide(hLocalLinearityStdDevReco);

    hLocalLinearityStdDevCutEff = (TH1D*) hLocalLinearityStdDevRecoTrue->Clone("hLocalLinearityStdDevCutEff");
    hLocalLinearityStdDevCutEff->Scale(1.0 / hTotalEvents->Integral(2,2));

    hLocalLinearityMinCutPur = (TH1D*) hLocalLinearityMinRecoTrue->Clone("hLocalLinearityMinCutPur");
    hLocalLinearityMinCutPur->Divide(hLocalLinearityMinReco);

    hLocalLinearityMinCutEff = (TH1D*) hLocalLinearityMinRecoTrue->Clone("hLocalLinearityMinCutEff");
    hLocalLinearityMinCutEff->Scale(1.0 / hTotalEvents->Integral(2,2));

    ////////////////////////////////////////////
    // Compute reco efficiency for scattering //
    ////////////////////////////////////////////

    TEfficiency* hAllScatteringReconstructionEfficiency = nullptr;
    hAllScatteringReconstructionEfficiency = new TEfficiency(*hAllScatteringReconstructed, *hAllScatteringTotal);

    TEfficiency* hInelasticScatteringReconstructionEfficiency = nullptr;
    hInelasticScatteringReconstructionEfficiency = new TEfficiency(*hInelasticScatteringReconstructed, *hInelasticScatteringTotal);

    TEfficiency* hInelasticScatteringReconstructionEfficiency0pBkg = nullptr;
    hInelasticScatteringReconstructionEfficiency0pBkg = new TEfficiency(*hInelasticScatteringReconstructed0pBkg, *hInelasticScatteringTotal0pBkg);

    TEfficiency* hInelasticScatteringReconstructionEfficiencyNpBkg = nullptr;
    hInelasticScatteringReconstructionEfficiencyNpBkg = new TEfficiency(*hInelasticScatteringReconstructedNpBkg, *hInelasticScatteringTotalNpBkg);

    TEfficiency* hInelasticScatteringReconstructionEfficiencyAngle = nullptr;
    hInelasticScatteringReconstructionEfficiencyAngle = new TEfficiency(*hInelasticScatteringRecoAngle, *hInelasticScatteringAngle);

    TEfficiency* hElasticScatteringReconstructionEfficiencyAngle = nullptr;
    hElasticScatteringReconstructionEfficiencyAngle = new TEfficiency(*hElasticScatteringRecoAngle, *hElasticScatteringAngle);

    TEfficiency* hAllScatteringReconstructionEfficiencyAngle = nullptr;
    hAllScatteringReconstructionEfficiencyAngle = new TEfficiency(*hAllScatteringRecoAngle, *hAllScatteringAngle);

    TEfficiency* hInelasticScatteringReconstructionEnergyDiff = nullptr;
    hInelasticScatteringReconstructionEnergyDiff = new TEfficiency(*hInelasticScatteringRecoEnergyDiff, *hInelasticScatteringEnergyDiff);

    std::vector<TEfficiency*> EfficiencyPlots = {
        hAllScatteringReconstructionEfficiency,
        hInelasticScatteringReconstructionEfficiency,
        hInelasticScatteringReconstructionEfficiency0pBkg,
        hInelasticScatteringReconstructionEfficiencyNpBkg,
        hInelasticScatteringReconstructionEfficiencyAngle,
        hElasticScatteringReconstructionEfficiencyAngle,
        hAllScatteringReconstructionEfficiencyAngle,
        hInelasticScatteringReconstructionEnergyDiff
    };

    std::vector<TString> EfficiencyPlotsTitles = {
        "Scattering/RecoEfficiencyAllScatteredPions",
        "Scattering/RecoEfficiencyInelScatteredPions",
        "Scattering/RecoEfficiencyInelScatteredPions0pBkg",
        "Scattering/RecoEfficiencyInelScatteredPionsNpBkg",
        "Scattering/RecoEfficiencyInelScatteredAngle",
        "Scattering/RecoEfficiencyElasticScatteredAngle",
        "Scattering/RecoEfficiencyAllScatteredAngle",
        "Scattering/RecoEfficiencyInelScatteredEnergyDiff"
    };

    std::vector<TString> EfficiencyPlotsXLabels = {
        "Outgoing pion KE (GeV/c)",
        "Outgoing pion KE (GeV/c)",
        "Outgoing pion KE (GeV/c)",
        "Outgoing pion KE (GeV/c)",
        "Scattering angle (deg)",
        "Scattering angle (deg)",
        "Scattering angle (deg)",
        "Energy difference (GeV/c)"
    };

    printEfficiencyPlots(
        SaveDir, FontStyle, TextSize,
        EfficiencyPlots,
        EfficiencyPlotsTitles,
        EfficiencyPlotsXLabels
    );

    ////////////////////////////
    // Compute cross-sections //
    ////////////////////////////

    float rho            = 1396; // kg / m^3
    float molar_mass     = 39.95; // g / mol
    float g_per_kg       = 1000; 
    float avogadro       = 6.022e+23; // number/mol
    float number_density = rho * g_per_kg / molar_mass * avogadro;
    float slab_width     = 0.0047; // in m

    for (int iBin = 1; iBin <= hIncidentKE->GetNbinsX(); ++iBin) {
        if (hIncidentKE->GetBinContent(iBin) == 0) continue;

        float pionAbsCrossSection   = ((hPionAbsKE->GetBinContent(iBin) / hIncidentKE->GetBinContent(iBin)) * (1 / number_density) * (1 / slab_width)) * (1 / 1e-28);
        float pion0pAbsCrossSection = ((h0pPionAbsKE->GetBinContent(iBin) / hIncidentKE->GetBinContent(iBin)) * (1 / number_density) * (1 / slab_width)) * (1 / 1e-28);
        float pionNpAbsCrossSection = ((hNpPionAbsKE->GetBinContent(iBin) / hIncidentKE->GetBinContent(iBin)) * (1 / number_density) * (1 / slab_width)) * (1 / 1e-28);

        hPionAbsRecoCrossSection->SetBinContent(iBin, pionAbsCrossSection);
        hPion0pAbsRecoCrossSection->SetBinContent(iBin, pion0pAbsCrossSection);
        hPionNpAbsRecoCrossSection->SetBinContent(iBin, pionNpAbsCrossSection);

        float denomError = pow(hIncidentKE->GetBinContent(iBin),0.5);
        float denom      = hIncidentKE->GetBinContent(iBin);
        if (denom == 0) { continue; }
        float termDenom = denomError / denom;

        float numErrorAbs = pow(hPionAbsKE->GetBinContent(iBin), 0.5);
        float numAbs      = hPionAbsKE->GetBinContent(iBin);
        if(numAbs == 0){ continue; }
        float termAbs       = numErrorAbs / numAbs;
        float totalAbsError = (pionAbsCrossSection) * (std::pow(((termAbs * termAbs) + (termDenom * termDenom)), 0.5)) * (1 / number_density) * (1 / slab_width) * (1e26);
        hPionAbsRecoCrossSection->SetBinError(iBin, totalAbsError);

        float numErrorAbs0p = pow(h0pPionAbsKE->GetBinContent(iBin), 0.5);
        float numAbs0p      = h0pPionAbsKE->GetBinContent(iBin);
        if(numAbs0p == 0){ continue; }
        float termAbs0p       = numErrorAbs0p / numAbs0p;
        float totalAbs0pError = (pion0pAbsCrossSection) * (std::pow(((termAbs0p * termAbs0p) + (termDenom * termDenom)), 0.5)) * (1 / number_density) * (1 / slab_width) * (1e26);
        hPion0pAbsRecoCrossSection->SetBinError(iBin, totalAbs0pError);

        float numErrorAbsNp = pow(hNpPionAbsKE->GetBinContent(iBin), 0.5);
        float numAbsNp      = hNpPionAbsKE->GetBinContent(iBin);
        if(numAbsNp == 0){ continue; }
        float termAbsNp       = numErrorAbsNp / numAbsNp;
        float totalAbsNpError = (pionNpAbsCrossSection) * (std::pow(((termAbsNp * termAbsNp) + (termDenom * termDenom)), 0.5)) * (1 / number_density) * (1 / slab_width) * (1e26);
        hPionNpAbsRecoCrossSection->SetBinError(iBin, totalAbsNpError);
    }

    // Compute corrections
    TH1D* hPsiInc = (TH1D*) hIncidentKEOnlyPions->Clone("hPsiInc");
    hPsiInc->Divide(hTrueIncidentKE);

    TH1D* hPsiIntPionAbs = (TH1D*) hPionAbsKEOnlyAbs->Clone("hPsiIntPionAbs");
    hPsiIntPionAbs->Divide(hTruePionAbsKE);

    TH1D* hPsiIntPionAbs0p = (TH1D*) h0pPionAbsKEOnly0pAbs->Clone("hPsiIntPionAbs0p");
    hPsiIntPionAbs0p->Divide(hTruePionAbs0pKE);

    TH1D* hPsiIntPionAbsNp = (TH1D*) hNpPionAbsKEOnlyNpAbs->Clone("hPsiIntPionAbsNp");
    hPsiIntPionAbsNp->Divide(hTruePionAbsNpKE);

    TH1D* hCInc = (TH1D*) hIncidentKEOnlyPions->Clone("hCInc");
    hCInc->Divide(hIncidentKE);

    TH1D* hCIntPionAbs = (TH1D*) hPionAbsKEOnlyAbs->Clone("hCIntPionAbs");
    hCIntPionAbs->Divide(hPionAbsKE);

    TH1D* hCIntPionAbs0p = (TH1D*) h0pPionAbsKEOnly0pAbs->Clone("hCIntPionAbs0p");
    hCIntPionAbs0p->Divide(h0pPionAbsKE);

    TH1D* hCIntPionAbsNp = (TH1D*) hNpPionAbsKEOnlyNpAbs->Clone("hCIntPionAbsNp");
    hCIntPionAbsNp->Divide(hNpPionAbsKE);

    // Apply corrections
    hPionAbsRecoCrossSection->Multiply(hPsiInc);
    hPionAbsRecoCrossSection->Divide(hPsiIntPionAbs);
    hPionAbsRecoCrossSection->Multiply(hCIntPionAbs);
    hPionAbsRecoCrossSection->Divide(hCInc);

    hPion0pAbsRecoCrossSection->Multiply(hPsiInc);
    hPion0pAbsRecoCrossSection->Divide(hPsiIntPionAbs0p);
    hPion0pAbsRecoCrossSection->Multiply(hCIntPionAbs0p);
    hPion0pAbsRecoCrossSection->Divide(hCInc);

    hPionNpAbsRecoCrossSection->Multiply(hPsiInc);
    hPionNpAbsRecoCrossSection->Divide(hPsiIntPionAbsNp);
    hPionNpAbsRecoCrossSection->Multiply(hCIntPionAbsNp);
    hPionNpAbsRecoCrossSection->Divide(hCInc);

    std::vector<TH1D*> XSecHistos = {
        hPionAbsRecoCrossSection, 
        hPion0pAbsRecoCrossSection,
        hPionNpAbsRecoCrossSection
    };

    TCanvas* XSecCanvas = new TCanvas("Canvas", "Canvas", 205, 34, 1024, 768);

    ////////////////////////////
    // Draw them individually //
    ////////////////////////////

    // gStyle->SetOptStat(0);
    for (size_t i = 0; i < XSecHistos.size(); ++i) {
        if (!XSecHistos[i]) continue;

        double maxY = XSecHistos[i]->GetMaximum();
        for (int bin = 1; bin <= XSecHistos[i]->GetNbinsX(); ++bin) {
            double y = XSecHistos[i]->GetBinContent(bin) + XSecHistos[i]->GetBinError(bin);
            if (y > maxY) maxY = y;
        }
        XSecHistos[i]->GetYaxis()->SetRangeUser(0, 1.1 * maxY);

        XSecHistos[i]->SetLineColor(kBlack);
        XSecHistos[i]->SetMarkerColor(kBlack);
        XSecHistos[i]->SetMarkerStyle(20);
        XSecHistos[i]->SetMarkerSize(1.2);

        XSecHistos[i]->SetFillColor(kBlack);
        XSecHistos[i]->SetFillStyle(3001);        // e.g. 3001 = hatched; 1001 = solid
        XSecHistos[i]->SetFillColorAlpha(kBlack, 0.2);

        XSecHistos[i]->SetTitle(XSecHistos[i]->GetName());
        XSecHistos[i]->GetXaxis()->SetTitle("Kinetic Energy [MeV]");
        XSecHistos[i]->GetYaxis()->SetTitle("Cross Section [barn]");
        // XSecHistos[i]->Draw("HISTO E1");
        XSecHistos[i]->Draw("E1");

        XSecCanvas->SetLeftMargin(0.13);
        XSecCanvas->SetBottomMargin(0.13);

        TString histName = XSecHistos[i]->GetName();
        histName.Remove(0, 1); // Remove the first character
        XSecCanvas->SaveAs(SaveDir + "CrossSection/" + histName + ".png");
    }

    ///////////////////////
    // Draw them stacked //
    ///////////////////////

    // Create a stacked histogram for pion absorption 0p and Np
    THStack* hStackPionAbs = new THStack("hStackPionAbs", "Pion Absorption Reco Cross Section;Kinetic Energy [MeV];Cross Section [barn]");

    // Set colors for stack
    hPion0pAbsRecoCrossSection->SetFillColor(kAzure+1);
    hPion0pAbsRecoCrossSection->SetLineColor(kAzure+1);

    hPionNpAbsRecoCrossSection->SetFillColor(kOrange+7);
    hPionNpAbsRecoCrossSection->SetLineColor(kOrange+7);

    // Add to stack
    hStackPionAbs->Add(hPion0pAbsRecoCrossSection, "H");
    hStackPionAbs->Add(hPionNpAbsRecoCrossSection, "H");

    hStackPionAbs->Draw("hist");
    hStackPionAbs->SetMaximum(1.1 * std::max(hPionAbsRecoCrossSection->GetMaximum(), hStackPionAbs->GetMaximum()));
    hStackPionAbs->GetXaxis()->SetTitle("Kinetic Energy [MeV]");
    hStackPionAbs->GetYaxis()->SetTitle("Cross Section [barn]");

    // Draw total pion absorption cross-section on top
    hPionAbsRecoCrossSection->SetLineColor(kBlack);
    hPionAbsRecoCrossSection->SetLineWidth(2);
    hPionAbsRecoCrossSection->SetMarkerSize(1.2);
    hPionAbsRecoCrossSection->SetMarkerStyle(20);
    hPionAbsRecoCrossSection->SetMarkerColor(kBlack);
    hPionAbsRecoCrossSection->Draw("E1 SAME");

    // Add legend
    TLegend* leg = new TLegend(0.60, 0.65, 0.85, 0.85);
    leg->SetTextFont(FontStyle);
    leg->SetTextSize(TextSize * 0.9);
    leg->AddEntry(hPion0pAbsRecoCrossSection, "Reco 0p", "f");
    leg->AddEntry(hPionNpAbsRecoCrossSection, "Reco Np", "f");
    leg->AddEntry(hPionAbsRecoCrossSection, "Reco Total", "lep");

    leg->Draw();
    XSecCanvas->SaveAs(SaveDir + "CrossSection/PionAbsorptionStacked.png");
    delete leg;

    ///////////////////////////////////
    // Draw them overlaid with truth //
    ///////////////////////////////////

    hPionAbsTrueCrossSection->SetMaximum(1.15 * std::max({
        hPionAbsTrueCrossSection->GetMaximum(),
        hPionAbsRecoCrossSection->GetMaximum()
    }));

    hPionAbsTrueCrossSection->GetXaxis()->SetTitle("Kinetic Energy [MeV]");
    hPionAbsTrueCrossSection->GetYaxis()->SetTitle("Cross Section [barn]");
    hPionAbsTrueCrossSection->SetTitle("Pion Absorption Cross-Section");

    // Draw true cross-sections
    hPionAbsTrueCrossSection->SetFillColorAlpha(kGray+2, 0.1);
    hPionAbsTrueCrossSection->SetMarkerStyle(20);
    hPionAbsTrueCrossSection->SetMarkerColor(kGray+2);
    hPionAbsTrueCrossSection->SetLineColor(kGray+2);
    hPionAbsTrueCrossSection->SetMarkerSize(1.2);
    hPionAbsTrueCrossSection->Draw("HIST E1");

    hPion0pAbsTrueCrossSection->SetFillColorAlpha(kAzure+4, 0.1);
    hPion0pAbsTrueCrossSection->SetMarkerStyle(20);
    hPion0pAbsTrueCrossSection->SetMarkerColor(kAzure+4);
    hPion0pAbsTrueCrossSection->SetLineColor(kAzure+4);
    hPion0pAbsTrueCrossSection->SetMarkerSize(1.2);
    hPion0pAbsTrueCrossSection->Draw("HIST E1 SAME");

    hPionNpAbsTrueCrossSection->SetFillColorAlpha(kOrange+9, 0.2);
    hPionNpAbsTrueCrossSection->SetMarkerStyle(20);
    hPionNpAbsTrueCrossSection->SetMarkerColor(kOrange+9);
    hPionNpAbsTrueCrossSection->SetLineColor(kOrange+9);
    hPionNpAbsTrueCrossSection->SetMarkerSize(1.2);
    hPionNpAbsTrueCrossSection->Draw("HIST E1 SAME");

    // Draw reco cross sections
    hPionAbsRecoCrossSection->SetMarkerColor(kBlack);
    hPionAbsRecoCrossSection->SetLineColor(kBlack);
    hPionAbsRecoCrossSection->SetLineWidth(2);
    hPionAbsRecoCrossSection->Draw("E1 SAME");

    hPion0pAbsRecoCrossSection->SetMarkerColor(kAzure+1);
    hPion0pAbsRecoCrossSection->SetLineColor(kAzure+1);
    hPion0pAbsRecoCrossSection->SetLineWidth(2);
    hPion0pAbsRecoCrossSection->Draw("E1 SAME");

    hPionNpAbsRecoCrossSection->SetMarkerColor(kOrange+7);
    hPionNpAbsRecoCrossSection->SetLineColor(kOrange+7);
    hPionNpAbsRecoCrossSection->SetLineWidth(2);
    hPionNpAbsRecoCrossSection->Draw("E1 SAME");

    // Add legend
    TLegend* leg2 = new TLegend(0.60, 0.55, 0.85, 0.85);
    leg2->SetTextFont(FontStyle);
    leg2->SetTextSize(TextSize * 0.9);
    leg2->AddEntry(hPionAbsTrueCrossSection, "True Total", "f");
    leg2->AddEntry(hPionAbsRecoCrossSection, "Reco Total", "lep");
    leg2->AddEntry(hPion0pAbsTrueCrossSection, "True 0p", "f");
    leg2->AddEntry(hPion0pAbsRecoCrossSection, "Reco 0p", "lep");
    leg2->AddEntry(hPionNpAbsTrueCrossSection, "True Np", "f");
    leg2->AddEntry(hPionNpAbsRecoCrossSection, "Reco Np", "lep");

    leg2->Draw();
    XSecCanvas->SaveAs(SaveDir + "CrossSection/PionAbsorptionOverlay.png");
    delete leg2;

    delete XSecCanvas;

    // gStyle->SetOptStat(1);

    //////////////////////////////////////////
    // Get bin-by-bin efficiency and purity //
    //////////////////////////////////////////

    auto AbsEffPur = getBinEfficiencyAndPurity(hTruePionAbsKE, hPionAbsKE, hPionAbsKECorrect);
    TH1* hAbsEff = AbsEffPur.first;
    TH1* hAbsPur = AbsEffPur.second;

    auto Abs0pEffPur = getBinEfficiencyAndPurity(hTruePionAbs0pKE, h0pPionAbsKE, h0pPionAbsKECorrect);
    TH1* hAbs0pEff = Abs0pEffPur.first;
    TH1* hAbs0pPur = Abs0pEffPur.second;

    auto AbsNpEffPur = getBinEfficiencyAndPurity(hTruePionAbsNpKE, hNpPionAbsKE, hNpPionAbsKECorrect);
    TH1* hAbsNpEff = AbsNpEffPur.first;
    TH1* hAbsNpPur = AbsNpEffPur.second;

    auto ScatterEffPur = getBinEfficiencyAndPurity(hTruePionScatteringKE, hPionScatteringKE, hPionScatteringKECorrect);
    TH1* hScatterEff = ScatterEffPur.first;
    TH1* hScatterPur = ScatterEffPur.second;

    auto Scatter0pEffPur = getBinEfficiencyAndPurity(hTruePion0pScatteringKE, hPionScatteringKE0p, h0pPionScatteringKECorrect);
    TH1* hScatter0pEff = Scatter0pEffPur.first;
    TH1* hScatter0pPur = Scatter0pEffPur.second;

    auto ScatterNpEffPur = getBinEfficiencyAndPurity(hTruePionNpScatteringKE, hPionScatteringKENp, hNpPionScatteringKECorrect);
    TH1* hScatterNpEff = ScatterNpEffPur.first;
    TH1* hScatterNpPur = ScatterNpEffPur.second;

    ////////////////////////////////////////
    // Print information for cylinder cut //
    ////////////////////////////////////////

    std::cout << std::endl;
    std::cout << "Results from cylinder cut" << std::endl;
    std::cout << "  Number of pions passing: " << numSurvivingPions << std::endl;
    std::cout << "  Number of muons passing: " << numSurvivingMuons << std::endl;
    std::cout << "  Number of electrons passing: " << numSurvivingElectrons << std::endl;
    std::cout << std::endl;


    for (int i = 1; i <= hElectronShowerTrkContained->GetNbinsX(); ++i) {
        double contained   = hElectronShowerTrkContained->GetBinContent(i);
        double uncontained = hElectronShowerTrkUnContained->GetBinContent(i);
        double total = contained + uncontained;
        double ratio = (total > 0) ? (contained / total) * 100.0 : 0.0;
        hElectronShowerTrkContainedRatio->SetBinContent(i, ratio);
        
        // Set error using binomial error propagation
        double err = (total > 0) ? 100.0 * std::sqrt(contained * (total - contained)) / (total * total) : 0.0;
        hElectronShowerTrkContainedRatio->SetBinError(i, err);
    }

    /////////////////////////////////////
    // Print information for vertex Np //
    /////////////////////////////////////

    std::cout << std::endl;
    std::cout << "Percentage of abs Np events with improved vertex: " << (((float) betterVertexNp) / ((float) betterVertexNp + (float) worseVertexNp) * 100.0) << "%" << std::endl;
    std::cout << std::endl;

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
        // Stitching
        {hStitchedDistanceFromVertex, hStitchedOriginalDistanceFromVertex},
        {hStitchAsPionAndProton, hStitchAsPionBraggPeak, hStitchAsMIP, hStitchAsProton},
        {hStitchFracBreakPoints},

        // Secondary PID
        {hSecondaryPionChi2Pions, hSecondaryPionChi2Protons, hSecondaryPionChi2Others},
        {hSecondaryProtonChi2Pions, hSecondaryProtonChi2Protons, hSecondaryProtonChi2Others},
        {hSecondaryMeanDEDXPions, hSecondaryMeanDEDXProtons, hSecondaryMeanDEDXOthers},

        // Scattering
        {hTotalBackgroundAllScatteringAngle, h0pBackgroundAllScatteringAngle, hNpBackgroundAllScatteringAngle},
        {hTotalBackgroundElasticScatteringAngle, h0pBackgroundElasticScatteringAngle, hNpBackgroundElasticScatteringAngle},
        {hTotalBackgroundInelScatteringAngle, h0pBackgroundInelScatteringAngle, hNpBackgroundInelScatteringAngle},
        {hTotalBackgroundInelScatteringLength, h0pBackgroundInelScatteringLength, hNpBackgroundInelScatteringLength},
        {hTotalBackgroundAllScatteringVertexKE, h0pBackgroundAllScatteringVertexKE, hNpBackgroundAllScatteringVertexKE},
        {hTotalBackgroundElasticScatteringVertexKE, h0pBackgroundElasticScatteringVertexKE, hNpBackgroundElasticScatteringVertexKE},
        {hTotalBackgroundInelScatteringVertexKE, h0pBackgroundInelScatteringVertexKE, hNpBackgroundInelScatteringVertexKE},
        {hTotalBackgroundInelasticScatteringEnergyDiff, h0pBackgroundInelasticScatteringEnergyDiff, hNpBackgroundInelasticScatteringEnergyDiff},
        {hInelScatteringBackgroundOriginalDistanceToPrimaryVertex, hInelScatteringBackgroundOriginalDistanceToSecondaryVertex, hInelScatteringBackgroundDistanceToPrimaryVertex, hInelScatteringBackgroundDistanceToSecondaryVertex},

        // Hit clustering
        {hHitClusters0p, hHitClusters0pBackground, hHitClustersNp, hHitClustersNpBackground},
        {hHitClustersSize0p, hHitClustersSize0pBackground, hHitClustersSizeNp, hHitClustersSizeNpBackground},
        {hHitLargeClusters0p, hHitLargeClusters0pBackground, hHitLargeClustersNp, hHitLargeClustersNpBackground},

        // Scattering reconstruction
        {hInelasticScatteringTotal, hInelasticScatteringReconstructed},
        {hInelasticScatteringTotal0pBkg, hInelasticScatteringReconstructed0pBkg},
        {hInelasticScatteringTotalNpBkg, hInelasticScatteringReconstructedNpBkg},
        {hInelasticScatteringEnergyDiff, hInelasticScatteringRecoEnergyDiff},
        {hAllScatteringTotal, hAllScatteringReconstructed},

        // Local linearity
        {hMinimumLinearity0p, hMinimumLinearity0pBackground},
        {hStdDevLinearity0p, hStdDevLinearity0pBackground},
        {hMinimumLinearityNp, hMinimumLinearityNpBackground},
        {hStdDevLinearityNp, hStdDevLinearityNpBackground},
        {hMaxLinearityD0p, hMaxLinearityD0pBackground},
        {hMaxLinearityDD0p, hMaxLinearityDD0pBackground},
        {hLocalLinearityDerivativeCutPur, hLocalLinearityDerivativeCutEff},
        {hLocalLinearityStdDevCutPur, hLocalLinearityStdDevCutEff},
        {hLocalLinearityMinCutPur, hLocalLinearityMinCutEff},

        // All scattering 
        {hAllScatteringAngle, hElasticScatteringAngle, hInelasticScatteringAngle},
        {hAllScatteringVertexKE, hElasticScatteringVertexKE, hInelasticScatteringVertexKE},

        // Incident histograms
        {hIncidentKE},
        {hIncidentKECorrect, hIncidentKEMuons, hIncidentKEElectrons},

        // Pion absorption interacting slices
        {hPionAbsKE, h0pPionAbsKE, hNpPionAbsKE},
        {hPionAbsKECorrect, hPionAbsKE0pScattering, hPionAbsKENpScattering, hPionAbsKEChEx, hPionAbsKEMuons, hPionAbsKEElectrons, hPionAbsKEOther},
        {h0pPionAbsKECorrect, h0pPionAbsKENpAbs, h0pPionAbsKE0pScattering, h0pPionAbsKENpScattering, h0pPionAbsKEChEx, h0pPionAbsKEMuons, h0pPionAbsKEElectrons, h0pPionAbsKEOther},
        {hNpPionAbsKECorrect, hNpPionAbsKE0pAbs, hNpPionAbsKE0pScattering, hNpPionAbsKENpScattering, hNpPionAbsKEChEx, hNpPionAbsKEMuons, hNpPionAbsKEElectrons, hNpPionAbsKEOther},

        // Scattering interacting slices
        {hPionScatteringKE, hPionScatteringKE0p, hPionScatteringKENp},
        {hPionScatteringKECorrect, hPionScatteringKE0pAbs, hPionScatteringKENpAbs, hPionScatteringKEChEx, hPionScatteringKEMuons, hPionScatteringKEElectrons, hPionScatteringKEOther},
        {h0pPionScatteringKECorrect, h0pPionScatteringKENpScatter, h0pPionScatteringKE0pAbs, h0pPionScatteringKENpAbs, h0pPionScatteringKEChEx, h0pPionScatteringKEMuons, h0pPionScatteringKEElectrons, h0pPionScatteringKEOther},
        {hNpPionScatteringKECorrect, hNpPionScatteringKE0pScatter, hNpPionScatteringKE0pAbs, hNpPionScatteringKENpAbs, hNpPionScatteringKEChEx, hNpPionScatteringKEMuons, hNpPionScatteringKEElectrons, hNpPionScatteringKEOther},

        // Bin-by-bin efficiency/purity
        {hAbsEff, hAbsPur},
        {hAbs0pEff, hAbs0pPur},
        {hAbsNpEff, hAbsNpPur},
        {hScatterEff, hScatterPur},
        {hScatter0pEff, hScatter0pPur},
        {hScatterNpEff, hScatterNpPur},

        // Charge exchange events
        {hChExchShowerRecoTrkLengths, hChExchShowerTruthTrkLengths},
        {hChExchShowerDaughtersForwardMom},
        {hChExchShowerTrksConeContained, hChExchShowerTrksConeUnContained},
        {hChExchShowerRecoTrksConeContained, hChExchShowerRecoTrksConeUnContained},
        {hChExchShowerTrksCylinderContained, hChExchShowerTrksCylinderUnContained},
        {hChExchShowerRecoTrksCylinderContained, hChExchShowerRecoTrksCylinderUnContained},

        // Cylinder cuts
        {hTrkLengthsInCylinder},
        {hTrkLengthsInCylinderPions, hTrkLengthsInCylinderMuons, hTrkLengthsInCylinderElectrons},
        {hTrkLengthsInCylinderPionAbs0p, hTrkLengthsInCylinderPionAbsNp, hTrkLengthsInCylinderPionScattering, hTrkLengthsInCylinderPionChExch, hTrkLengthsInCylinderPionOther, hTrkLengthsInCylinderMuons, hTrkLengthsInCylinderElectrons},
        {hNumTracksInCylinder},
        {hNumTracksInCylinderPions, hNumTracksInCylinderMuons, hNumTracksInCylinderElectrons},
        {hSmallTrksInCylinder},
        {hSmallTrksInCylinderPions, hSmallTrksInCylinderMuons, hSmallTrksInCylinderElectrons},
        {hSmallTrksInCylinderPionAbs0p, hSmallTrksInCylinderPionAbsNp, hSmallTrksInCylinderPionScattering, hSmallTrksInCylinderPionChExch, hSmallTrksInCylinderPionOther, hSmallTrksInCylinderMuons, hSmallTrksInCylinderElectrons},

        // Sliced cone cut
        {hSlicedConeTracksPiAbs0p, hSlicedConeTracksPiAbsNp, hSlicedConeTracksPiScat, hSlicedConeTracksPiChEx, hSlicedConeTracksPiOther, hSlicedConeTracksMu, hSlicedConeTracksE},
        {hSlicedConeNumTrksPiAbs0p, hSlicedConeNumTrksPiAbsNp, hSlicedConeNumTrksPiScat, hSlicedConeNumTrksPiChEx, hSlicedConeNumTrksPiOther, hSlicedConeNumTrksMu, hSlicedConeNumTrksE},
        {hSlicedConeNumSmallTrksPiAbs0p, hSlicedConeNumSmallTrksPiAbsNp, hSlicedConeNumSmallTrksPiScat, hSlicedConeNumSmallTrksPiChEx, hSlicedConeNumSmallTrksPiOther, hSlicedConeNumSmallTrksMu, hSlicedConeNumSmallTrksE},

        // Beamline showers
        {hElectronShowerRecoTrkLengths, hElectronShowerTrueTrkLengths},
        {hElectronShowerTrkContained, hElectronShowerTrkUnContained},
        {hElectronShowerTrkContainedRatio},

        // Multiple primaries
        {hNumPrimaries, hNumValidPrimaries},

        // Reconstructed secondary particles
        {hRecoSecondaryPions, hAllSecondaryPions},
        {hRecoSecondaryProtons, hAllSecondaryProtons}
    };

    std::vector<std::vector<TString>> PlotLabelGroups = {
        // Stitching
        {"Detected", "Original"},
        {"Stitched", "Bragg pion", "MIP", "Proton"},
        {"Stitched tracks"},

        // Secondary PID
        {"Pions", "Protons", "Others"},
        {"Pions", "Protons", "Others"},
        {"Pions", "Protons", "Others"},

        // Scattering
        {"All bkg", "0p bkg", "Np bkg"},
        {"All bkg", "0p bkg", "Np bkg"},
        {"All bkg", "0p bkg", "Np bkg"},
        {"All bkg", "0p bkg", "Np bkg"},
        {"All bkg", "0p bkg", "Np bkg"},
        {"All bkg", "0p bkg", "Np bkg"},
        {"All bkg", "0p bkg", "Np bkg"},
        {"All bkg", "0p bkg", "Np bkg"},
        {"Original to primary", "Original to secondary", "Detected to primary", "Detected to secondary"},

        // Hit clustering
        {"Reco 0p true", "Reco 0p bkg", "Reco Np true", "Reco Np bkg"},
        {"Reco 0p true", "Reco 0p bkg", "Reco Np true", "Reco Np bkg"},
        {"Reco 0p true", "Reco 0p bkg", "Reco Np true", "Reco Np bkg"},

        // Scattering reconstruction
        {"All", "Reconstructed"},
        {"All", "Reconstructed"},
        {"All", "Reconstructed"},
        {"All", "Reconstructed"},
        {"All", "Reconstructed"},

        // Local linearity
        {"Reco 0p true", "Reco 0p bkg"},
        {"Reco 0p true", "Reco 0p bkg"},
        {"Reco 0p true", "Reco 0p bkg"},
        {"Reco Np true", "Reco Np bkg"},
        {"Reco 0p true", "Reco 0p bkg"},
        {"Reco 0p true", "Reco 0p bkg"},
        {"Purity", "Efficiency"},
        {"Purity", "Efficiency"},
        {"Purity", "Efficiency"},

        // All scattering
        {"All", "Elastic", "Inelastic"},
        {"All", "Elastic", "Inelastic"},

        // Incident histograms
        {"Incident"},
        {"Pions", "Muons", "Electrons"},

        // Pion absorption interacting slices
        {"All abs", "0p abs", "Np abs"},
        {"Abs", "0p scatter", "Np scatter", "Ch. exch.", "Muon", "Electron", "Other"},
        {"0p abs", "Np abs", "0p scatter", "Np scatter", "Ch. exch.", "Muon", "Electron", "Other"},
        {"Np abs", "0p abs", "0p scatter", "Np scatter", "Ch. exch.", "Muon", "Electron", "Other"},

        // Scattering interacting slices        
        {"All scatter", "0p Scatter", "Np scatter"},
        {"Scatter", "0p abs", "Np abs", "Ch. exch.", "Muon", "Electron", "Other"},
        {"0p scatter", "Np scatter", "0p abs", "Np abs", "Ch. exch.", "Muon", "Electron", "Other"},
        {"Np scatter", "0p scatter", "0p abs", "Np abs", "Ch. exch.", "Muon", "Electron", "Other"},

        // Bin-by-bin efficiency/purity
        {"Efficiency", "Purity"},
        {"Efficiency", "Purity"},
        {"Efficiency", "Purity"},
        {"Efficiency", "Purity"},
        {"Efficiency", "Purity"},
        {"Efficiency", "Purity"},

        // Charge exchange events
        {"Reco", "Truth"},
        {"#pi^{0} daughters"},
        {"Contained", "Uncontained"},
        {"Contained", "Uncontained"},
        {"Contained", "Uncontained"},
        {"Contained", "Uncontained"},

        // Cylinder cuts
        {"All"},
        {"Pions", "Muons", "Electrons"},
        {"Abs 0p", "Abs Np", "Scatter", "Ch. exch.", "Other", "Muon", "Electron"},
        {"All"},
        {"Pions", "Muons", "Electrons"},
        {"All"},
        {"Pions", "Muons", "Electrons"},
        {"Abs 0p", "Abs Np", "Scatter", "Ch. exch.", "Other", "Muon", "Electron"},

        // Sliced cone cut
        {"Abs 0p", "Abs Np", "Scatter", "Ch. exch.", "Other", "Muon", "Electron"},
        {"Abs 0p", "Abs Np", "Scatter", "Ch. exch.", "Other", "Muon", "Electron"},
        {"Abs 0p", "Abs Np", "Scatter", "Ch. exch.", "Other", "Muon", "Electron"},

        // Beamline showers
        {"Reco", "Truth"},
        {"Contained", "Uncontained"},
        {"\% contained"},

        // Multiple primaries
        {"Total", "Valid"},

        // Reconstructed secondary particles
        {"Reco", "All"},
        {"Reco", "All"}
    };

    std::vector<TString> PlotTitles = {
        // Stitching
        "PrimaryStitching/StitchedDistanceFromVertexNp",
        "PrimaryStitching/StitchedBackgroundTypes",
        "PrimaryStitching/StitchedFracBreakPoints",

        // Secondary PID
        "SecondaryPID/PionChi2SecondaryTracks",
        "SecondaryPID/ProtonChi2SecondaryTracks",
        "SecondaryPID/MeanDEDXSecondaryTracks",

        // Scattering
        "Scattering/AllScatteringBackgroundAngle",
        "Scattering/ElasticScatteringBackgroundAngle",
        "Scattering/InelasticScatteringBackgroundAngle",
        "Scattering/InelasticScatteringBackgroundLength",
        "Scattering/AllScatteringBackgroundVertexKE",
        "Scattering/ElasticScatteringBackgroundVertexKE",
        "Scattering/InelasticScatteringBackgroundVertexKE",
        "Scattering/InelasticScatteringBackgroundEnergyDiff",
        "Scattering/InelasticScatteringBackgroundDistanceFromVertex",

        // Hit clustering
        "HitClustering/NumberOfHitClusters",
        "HitClustering/AverageHitClusterSize",
        "HitClustering/NumberOfLargeHitClusters",

        // Scattering reconstruction
        "Scattering/RecoInelScatteredPions",
        "Scattering/RecoInelScatteredPions0pBkg",
        "Scattering/RecoInelScatteredPionsNpBkg",
        "Scattering/RecoInelScatteredPionsEnergyDiff",
        "Scattering/RecoAllScatteredPions",

        // Local linearity
        "LocalLinearity/MinLocalLinearity0pReco",
        "LocalLinearity/StdDevLinearity0pReco",
        "LocalLinearity/MinLocalLinearityNpReco",
        "LocalLinearity/StdDevLinearityNpReco",
        "LocalLinearity/MaxLocalLinearityDeriv0pReco",
        "LocalLinearity/MaxLocalLinearitySecondDeriv0pReco",
        "LocalLinearity/MaxLocalLinearityDerivCut",
        "LocalLinearity/MaxLocalLinearityStdDevCut",
        "LocalLinearity/MinLocalLinearityCut",

        // All scattering
        "Scattering/AllScatteringAngle",
        "Scattering/AllScatteringVertexKE",

        // Incident histograms
        "CrossSection/IncidentKE",
        "CrossSection/IncidentKEBreakdown",

        // Pion absorption interacting slices
        "CrossSection/AbsInteractingKE",
        "CrossSection/AbsInteractingKEBreakdown",
        "CrossSection/0pAbsInteractingKEBreakdown",
        "CrossSection/NpAbsInteractingKEBreakdown",

        // Scattering interacting slices
        "CrossSection/ScatteringInteractingKE",
        "CrossSection/ScatteringInteractingKEBreakdown",
        "CrossSection/0pScatteringInteractingKEBreakdown",
        "CrossSection/NpScatteringInteractingKEBreakdown",

        // Bin-by-bin efficiency/purity
        "EffPur/PionAbsorption",
        "EffPur/Pion0pAbsorption",
        "EffPur/PionNpAbsorption",
        "EffPur/PionScattering",
        "EffPur/Pion0pScattering",
        "EffPur/PionNpScattering",

        // Charge exchange events
        "ChExch/ShowerTrksLengths",
        "ChExch/NeutralPionDaughtersMom",
        "ChExch/ShowerTrueTrksConeContainment",
        "ChExch/ShowerRecoTrksConeContainment",
        "ChExch/ShowerTrueTrksCylinderContainment",
        "ChExch/ShowerRecoTrksCylinderContainment",

        // Cylinder cuts
        "Cylinder/TrackLengths",
        "Cylinder/TrackLengthsBreakdown",
        "Cylinder/TrackLengthsPiBreakdown",
        "Cylinder/NumTracks",
        "Cylinder/NumTracksBreakdown",
        "Cylinder/SmallTracks",
        "Cylinder/SmallTracksBreakdown",
        "Cylinder/SmallTracksPiBreakdown",

        // Sliced cone cut
        "SlicedCone/SlicedConeTracksPiBreakdown",
        "SlicedCone/SlicedConeNumTracksPiBreakdown",
        "SlicedCone/SlicedConeNumSmallTracksPiBreakdown",

        // Beamline showers
        "BeamlineShower/ElectronShowerTrkLengths",
        "BeamlineShower/ElectronShowerCylinderContainedTrks",
        "BeamlineShower/ElectronShowerCylinderContainedRatio",

        // Multiple primaries
        "Primaries/NumPrimaries",

        // Reconstructed secondary particles
        "SecondaryReco/RecoPions",
        "SecondaryReco/RecoProtons"
    };

    std::vector<TString> XLabels = {
        // Stitching
        "Distance from true vertex (cm)",
        "Background type",
        "Fractional break point",

        // Secondary PID
        "Pion reduced #chi^{2}",
        "Proton reduced #chi^{2}",
        "Mean dE/dx (MeV/cm)",

        // Scattering
        "Scattering angle (rad)",
        "Scattering angle (rad)",
        "Scattering angle (rad)",
        "Scattered pion length (cm)",
        "Vertex energy (GeV/c)",
        "Vertex energy (GeV/c)",
        "Vertex energy (GeV/c)",
        "Energy difference (GeV/c)",
        "Distance from true vertex (cm)",

        // Hit clustering
        "Number of induction plane hit clusters",
        "Induction plane hit cluster average size",
        "Number of induction plane large hit clusters",

        // Scattering reconstruction
        "Outgoing pion energy (GeV/c)",
        "Outgoing pion energy (GeV/c)",
        "Outgoing pion energy (GeV/c)",
        "Energy diff. (GeV/c)",
        "Outgoing pion energy (GeV/c)",

        // Local linearity
        "Minimum local linearity",
        "Standard dev. local linearity",
        "Minimum local linearity",
        "Standard dev. local linearity",
        "Maximum local linearity derivative",
        "Maximum local linearity second derivative",
        "Maximum local linearity derivative allowed",
        "Maximum local linearity std. dev. allowed",
        "1 - minimum local linearity allowed",

        // All scattering
        "Scattering angle (deg)",
        "Scattering vertex KE (GeV/c)",

        // Incident histograms
        "Kinetic Energy [MeV]",
        "Kinetic Energy [MeV]",

        // Pion absorption interacting slices
        "Kinetic Energy [MeV]",
        "Kinetic Energy [MeV]",
        "Kinetic Energy [MeV]",
        "Kinetic Energy [MeV]",

        // Scattering interacting slices
        "Kinetic Energy [MeV]",
        "Kinetic Energy [MeV]",
        "Kinetic Energy [MeV]",
        "Kinetic Energy [MeV]",

        // Bin-by-bin efficiency/purity
        "Kinetic Energy [MeV]",
        "Kinetic Energy [MeV]",
        "Kinetic Energy [MeV]",
        "Kinetic Energy [MeV]",
        "Kinetic Energy [MeV]",
        "Kinetic Energy [MeV]",

        // Charge exchange events
        "Track length (cm)",
        "Forward momentum [GeV/c]",
        "Track length (cm)",
        "Track length (cm)",
        "Track length (cm)",
        "Track length (cm)",

        // Cylinder cuts
        "Track length (cm)",
        "Track length (cm)",
        "Track length (cm)",
        "Number of tracks",
        "Number of tracks",
        "Number of small tracks",
        "Number of small tracks",
        "Number of small tracks",

        // Sliced cone cut
        "Track length (cm)",
        "# of tracks",
        "# of small tracks",

        // Beamline showers
        "Track length (cm)",
        "Track length (cm)",
        "Track length (cm)",

        // Multiple primaries
        "# of primaries",

        // Reconstructed secondary particles
        "Initial KE",
        "Initial KE"
    };

    std::vector<TString> YLabels = {
        // Stitching
        "Number of tracks",
        "Number of tracks",
        "Number of tracks",

        // Secondary PID
        "Number of tracks",
        "Number of tracks",
        "Number of tracks",

        // Scattering
        "Number of events",
        "Number of events",
        "Number of events",
        "Number of events",
        "Number of events",
        "Number of events",
        "Number of events",
        "Number of events",
        "Number of events",

        // Hit clustering
        "Number of events",
        "Number of events",
        "Number of events",

        // Scattering reconstruction
        "Number of tracks",
        "Number of tracks",
        "Number of tracks",
        "Number of tracks",
        "Number of tracks",

        // Local linearity
        "Number of tracks",
        "Number of tracks",
        "Number of tracks",
        "Number of tracks",
        "Number of tracks",
        "Number of tracks",
        "Purity/Efficiency",
        "Purity/Efficiency",
        "Purity/Efficiency",

        // All scattering
        "Number of events",
        "Number of events",

        // Incident histograms
        "Counts",
        "Counts",

        // Pion absorption interacting slices
        "Counts",
        "Counts",
        "Counts",
        "Counts",

        // Scattering interacting slices
        "Counts",
        "Counts",
        "Counts",
        "Counts",

        // Bin-by-bin efficiency/purity
        "Eff/Pur",
        "Eff/Pur",
        "Eff/Pur",
        "Eff/Pur",
        "Eff/Pur",
        "Eff/Pur",

        // Charge exchange events
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",

        // Cylinder cuts
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",

        // Sliced cone cut
        "Counts",
        "Counts",
        "Counts",

        // Beamline showers
        "Counts",
        "Counts",
        "Counts",

        // Multiple primaries
        "Counts",

        // Reconstructed secondary particles
        "Counts",
        "Counts"
    };

    std::vector<bool> PlotStacked = {
        // Stitching
        false,
        true,
        false,

        // Secondary PID
        false,
        false,
        false,

        // Scattering
        false,
        false,
        false,
        false,
        false,
        false,
        false,
        false,
        false,

        // Hit clustering
        false,
        false,
        false,

        // Scattering reconstruction
        false,
        false,
        false,
        false,
        false,

        // Local linearity
        false,
        false,
        false,
        false,
        false,
        false,
        false,
        false,
        false,

        // All scattering
        false,
        false,

        // Incident histograms
        false,
        true,

        // Pion absorption interacting slices
        false,
        true,
        true,
        true,

        // Scattering interacting slices
        false,
        true,
        true,
        true,

        // Bin-by-bin efficiency/purity
        false,
        false,
        false,
        false,
        false,
        false,

        // Charge exchange events
        false,
        false,
        true,
        true,
        true,
        true,

        // Cylinder cuts
        false, 
        true,
        true,
        false,
        true,
        false,
        true,
        true,

        // Sliced cone cut
        true,
        true,
        true,

        // Beamline showers
        false,
        true,
        false,

        // Multiple primaries
        false,

        // Reconstructed secondary particles
        false,
        false
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

    ///////////////////////////
    // Two-dimensional plots //
    ///////////////////////////

    std::vector<TH2*> TwoDPlots = {
        // Chi^2 optimization
        hProtonChi2TruePositives, 
        hProtonChi2FalsePositives,
        hProtonChi2TrueNegatives,
        hProtonChi2FalseNegatives,
        hPionChi2TruePositives, 
        hPionChi2FalsePositives,
        hPionChi2TrueNegatives,
        hPionChi2FalseNegatives,
        hProtonChi2FOM,
        hPionChi2FOM,

        // Scattering
        hTotalBackgroundInelScatteringLengthVSAngle,
        hTotalBackgroundInelScatteringLengthVSVertexKE,
        hTotalBackgroundInelScatteringVertexKEVSAngle,
        h0pBackgroundInelScatteringVertexKEVSAngle,
        hNpBackgroundInelScatteringVertexKEVSAngle,
        hTotalBackgroundElasticScatteringVertexKEVSAngle,
        h0pBackgroundElasticScatteringVertexKEVSAngle,
        hNpBackgroundElasticScatteringVertexKEVSAngle,
        hTotalBackgroundAllScatteringVertexKEVSAngle,
        h0pBackgroundAllScatteringVertexKEVSAngle,
        hNpBackgroundAllScatteringVertexKEVSAngle,
        hScatteringVertexKEVsOutKE,

        // Hit clustering
        hHitClusterCut0pReco,
        hHitClusterCut0pRecoTrue,
        hHitClusterCut0pPurity,

        // Beamline shower
        hElectronShowerDaughtersTrueIncident,
        hElectronShowerDaughtersRecoIncident,

        // Charge exchange
        hChExchShowerDaughtersTrueIncident,
        hChExchShowerDaughtersRecoIncident,
        hChExchShowerDaughtersExitAngles
    };

    std::vector<TString> TwoDTitles = {
        // Chi^2 optimization
        "SecondaryPID/2DProtonChi2TruePositives",
        "SecondaryPID/2DProtonChi2FalsePositives",
        "SecondaryPID/2DProtonChi2TrueNegatives",
        "SecondaryPID/2DProtonChi2FalseNegatives",
        "SecondaryPID/2DPionChi2TruePositives",
        "SecondaryPID/2DPionChi2FalsePositives",
        "SecondaryPID/2DPionChi2TrueNegatives",
        "SecondaryPID/2DPionChi2FalseNegatives",
        "SecondaryPID/2DProtonChi2FOM",
        "SecondaryPID/2DPionChi2FOM",

        // Scattering
        "Scattering/2DTotalBackgroundInelScatteringLengthVSAngle",
        "Scattering/2DTotalBackgroundInelScatteringLengthVSVertexKE",
        "Scattering/2DTotalBackgroundInelScatteringVertexKEVSAngle",
        "Scattering/2D0pBackgroundInelScatteringVertexKEVSAngle",
        "Scattering/2DNpBackgroundInelScatteringVertexKEVSAngle",
        "Scattering/2DTotalBackgroundElasticScatteringVertexKEVSAngle",
        "Scattering/2D0pBackgroundElasticScatteringVertexKEVSAngle",
        "Scattering/2DNpBackgroundElasticScatteringVertexKEVSAngle",
        "Scattering/2DTotalBackgroundAllScatteringVertexKEVSAngle",
        "Scattering/2D0pBackgroundAllScatteringVertexKEVSAngle",
        "Scattering/2DNpBackgroundAllScatteringVertexKEVSAngle",
        "Scattering/2DInelasticVertexVSOutKE",

        // Hit clustering
        "HitClustering/2DHitClusterCut0pReco",
        "HitClustering/2DHitClusterCut0pRecoTrue",
        "HitClustering/2DHitClusterCut0pPurity",

        // Beamline shower
        "BeamlineShower/2DElectronShowerDaughtersTrueIncident",
        "BeamlineShower/2DElectronShowerDaughtersRecoIncident",

        // Charge exchange
        "ChExch/2DChExchShowerDaughtersTrueIncident",
        "ChExch/2DChExchShowerDaughtersRecoIncident",
        "ChExch/2DChExchShowerNeutralPiDaughtersAngles",
    };

    printTwoDPlots(SaveDir, TwoDPlots, TwoDTitles);

    //////////////////////////////////////////////////////
    // Create scatter plots for inelastic scattering KE //
    //////////////////////////////////////////////////////

    TCanvas* c = new TCanvas("Canvas", "Canvas", 205, 34, 1024, 768);

    int nDataPoints               = InelasticScatteringIncidentKE.size();
    TGraph *gIncidentVSOutgoingKE = new TGraph(nDataPoints, InelasticScatteringIncidentKE.data(), InelasticScatteringOutgoingKE.data());
    TGraph *gVertexVSOutgoingKE   = new TGraph(nDataPoints, InelasticScatteringVertexKE.data(), InelasticScatteringOutgoingKE.data());

    gIncidentVSOutgoingKE->SetTitle("Outgoing KE vs Incident KE;Incident KE [MeV];Outgoing KE [MeV]");
    gIncidentVSOutgoingKE->SetMarkerStyle(20);
    gIncidentVSOutgoingKE->SetMarkerColor(kBlack);
    gIncidentVSOutgoingKE->Draw("AP");
    c->SaveAs(SaveDir + "Scattering/InelasticScatteringIncidentVSOutgoingKE.png");

    gVertexVSOutgoingKE->SetTitle("Outgoing KE vs Vertex KE;Vertex KE [MeV];Outgoing KE [MeV]");
    gVertexVSOutgoingKE->SetMarkerStyle(20);
    gVertexVSOutgoingKE->SetMarkerColor(kBlack);
    gVertexVSOutgoingKE->Draw("AP");
    c->SaveAs(SaveDir + "Scattering/InelasticScatteringVertexVSOutgoingKE.png");

    int n0pDataPoints               = Background0pInelasticScatteringIncidentKE.size();
    TGraph *g0pIncidentVSOutgoingKE = new TGraph(n0pDataPoints, Background0pInelasticScatteringIncidentKE.data(), Background0pInelasticScatteringOutgoingKE.data());
    TGraph *g0pVertexVSOutgoingKE   = new TGraph(n0pDataPoints, Background0pInelasticScatteringVertexKE.data(), Background0pInelasticScatteringOutgoingKE.data());

    g0pIncidentVSOutgoingKE->SetTitle("Outgoing KE vs Incident KE;Incident KE [MeV];Outgoing KE [MeV]");
    g0pIncidentVSOutgoingKE->SetMarkerStyle(20);
    g0pIncidentVSOutgoingKE->SetMarkerColor(kBlack);
    g0pIncidentVSOutgoingKE->Draw("AP");
    c->SaveAs(SaveDir + "Scattering/Background0pInelasticScatteringIncidentVSOutgoingKE.png");

    g0pVertexVSOutgoingKE->SetTitle("Outgoing KE vs Vertex KE;Vertex KE [MeV];Outgoing KE [MeV]");
    g0pVertexVSOutgoingKE->SetMarkerStyle(20);
    g0pVertexVSOutgoingKE->SetMarkerColor(kBlack);
    g0pVertexVSOutgoingKE->Draw("AP");
    c->SaveAs(SaveDir + "Scattering/Background0pInelasticScatteringVertexVSOutgoingKE.png");

    int nNpDataPoints               = BackgroundNpInelasticScatteringIncidentKE.size();
    TGraph *gNpIncidentVSOutgoingKE = new TGraph(nNpDataPoints, BackgroundNpInelasticScatteringIncidentKE.data(), BackgroundNpInelasticScatteringOutgoingKE.data());
    TGraph *gNpVertexVSOutgoingKE   = new TGraph(nNpDataPoints, BackgroundNpInelasticScatteringVertexKE.data(), BackgroundNpInelasticScatteringOutgoingKE.data());

    gNpIncidentVSOutgoingKE->SetTitle("Outgoing KE vs Incident KE;Incident KE [MeV];Outgoing KE [MeV]");
    gNpIncidentVSOutgoingKE->SetMarkerStyle(20);
    gNpIncidentVSOutgoingKE->SetMarkerColor(kBlack);
    gNpIncidentVSOutgoingKE->Draw("AP");
    c->SaveAs(SaveDir + "Scattering/BackgroundNpInelasticScatteringIncidentVSOutgoingKE.png");

    gNpVertexVSOutgoingKE->SetTitle("Outgoing KE vs Vertex KE;Vertex KE [MeV];Outgoing KE [MeV]");
    gNpVertexVSOutgoingKE->SetMarkerStyle(20);
    gNpVertexVSOutgoingKE->SetMarkerColor(kBlack);
    gNpVertexVSOutgoingKE->Draw("AP");
    c->SaveAs(SaveDir + "Scattering/BackgroundNpInelasticScatteringVertexVSOutgoingKE.png");
}