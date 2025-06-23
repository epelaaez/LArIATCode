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
    {11, "Other"}
};

struct EventInfo {
    int run; int subrun; int event;
    float vertexX; float vertexY; float vertexZ;

    int                      wcMatchPDG;
    std::string              wcMatchProcess;
    std::vector<int>         wcMatchDaughtersPDG;
    std::vector<std::string> wcMatchDaughtersProcess;

    int                      truthPrimaryPDG;
    std::vector<int>         truthPrimaryDaughtersPDG;
    std::vector<std::string> truthPrimaryDaughtersProcess;

    std::vector<int>         truthSecondaryDaughtersPDG;
    std::vector<std::string> truthSecondaryDaughtersProcess;
    std::vector<double>      truthSecondaryDaughtersKE;

    int visibleProtons;
    int recoProtonCount;
    int backgroundNum;
    int tracksNearVertex;
    int secondaryInteractionTag;

    int    numHitClusters;
    double minLocalLinearity;
    double maxLocalLinearityD;
};

struct HitCluster {
    std::vector<int>   hitKeys;
    std::vector<float> hitX;
    std::vector<float> hitW;
    std::vector<float> hitCharge;
    std::vector<float> hitChargeCol;
};

// Reduced volume for interactions
const double RminX =  5.0;
const double RmaxX = 42.0;
const double RminY =-15.0; 
const double RmaxY = 15.0;
const double RminZ =  8.0;
const double RmaxZ = 82.0;

// Detector dimensions
const double minX =  0.0;
const double maxX = 47.0;
const double minY =-20.0; 
const double maxY = 20.0; 
const double minZ =  3.0;
const double maxZ = 87.0;

// Background types
int NUM_BACKGROUND_TYPES = 11;
//    0:  0p pion absorption
//    1:  Np pion absorption
//    2:  primary muon event
//    3:  primary electron event
//    4:  other primary event
//    5:  primary pion outside reduced volume
//    6:  pion inelastic scattering
//    7:  charge exchange
//    8:  double charge exchange
//    9:  capture at rest
//    10: decay
//    11: other

// Values for chi^2 secondary fits
double PION_CHI2_PION_VALUE     = 1.375;
double PION_CHI2_PROTON_VALUE   = 2.625;
double PROTON_CHI2_PION_VALUE   = 1.375;
double PROTON_CHI2_PROTON_VALUE = 0.125;

// Mean dE/dx threshold
double MEAN_DEDX_THRESHOLD = 4.0;

// Vertex radius
double VERTEX_RADIUS = 5.0;

// Number of points to use in mean dE/dx calculation
int MEAN_DEDX_NUM_TRAJ_POINTS = 5;

// Proton energy bounds
double PROTON_ENERGY_LOWER_BOUND = 0.075;
double PROTON_ENERGY_UPPER_BOUND = 1.0;

// Hit wire separation
double HIT_WIRE_SEPARATION = 0.4; // converts wire # to cm

// Window size for local linearity
int WINDOW_SIZE_LINEARITY = 8;

// Threshold(s) used to cut based on primary track local linearity
double LINEARITY_DERIVATIVE_THRESHOLD = 0.06e-3;
double STD_DEV_LINEARITY_THRESHOLD    = 0.045E-3;

//////////////////////
// Helper functions //
//////////////////////

void initializeProtonPoints(TGraph* gProton);
void initializePionPoints(TGraph* gPion);
void initializeMuonNoBraggPoints(TGraph* gMuonTG);
void printEventInfo(EventInfo event, std::ostream& os);
bool isWithinReducedVolume(double x, double y, double z);
double computeReducedChi2(const TGraph* theory, std::vector<double> xData, std::vector<double> yData,  bool dataReversed, int nPoints, int nOutliersToDiscard = 0, int nTrim = 0);
double distance(double x1, double x2, double y1, double y2, double z1, double z2);
double meanDEDX(std::vector<double> trackDEDX, bool isTrackReversed, int pointsToUse);
void printBackgroundInfo(TH1D* background_histo, std::ostream& os);
int isSecondaryInteractionAbsorption(std::vector<int> daughtersPDG, std::vector<string> daughtersProcess, std::vector<double> daughtersKE);
void printOneDPlots(
    TString dir, int fontStyle, double textSize,
    std::vector<std::vector<TH1*>> groups,
    std::vector<int> colors, 
    std::vector<std::vector<TString>> labels, 
    std::vector<TString> titles, 
    std::vector<TString> xlabels, 
    std::vector<TString> ylabels
);
void printTwoDPlots(TString dir, std::vector<TH2*> plots, std::vector<TString> titles);
bool isHitNearPrimary(std::vector<int>* primaryKey, std::vector<float>* hitX, std::vector<float>* hitW, float thisHitX, float thisHitW, float xThreshold, float wThreshold);
void Overlay_dEdx_RR_Reference_PP(TGraph* gProton, TGraph* gPion, bool addLegend = true, TVirtualPad* pad = gPad);
void printEfficiencyPlots(TString dir, int fontStyle, double textSize, std::vector<TEfficiency*> efficiencies, std::vector<TString> titles, std::vector<TString> xlabels);
std::vector<double> calcLinearityProfile(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, int nb);

///////////////////
// Main function //
///////////////////

void RecoAllAnalysis() {
    // Set defaults
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
    TString RootFilePath = "/exp/lariat/app/users/epelaez/files/RecoAllEval_histo.root"; 
    std::unique_ptr<TFile> File(TFile::Open(RootFilePath));
    TDirectory* Directory = (TDirectory*)File->Get("RecoAllEval");

    ///////////////////
    // Load branches //
    ///////////////////

    // Load tree and branches
    TTree* tree = (TTree*) Directory->Get<TTree>("RecoEvalTree");

    int run, subrun, event;
    tree->SetBranchAddress("run", &run); 
    tree->SetBranchAddress("subrun", &subrun); 
    tree->SetBranchAddress("event", &event);

    // Signal information
    bool isPionAbsorptionSignal; int backgroundType, numVisibleProtons; 
    tree->SetBranchAddress("isPionAbsorptionSignal", &isPionAbsorptionSignal);
    tree->SetBranchAddress("backgroundType", &backgroundType);
    tree->SetBranchAddress("numVisibleProtons", &numVisibleProtons);

    // Cut information
    bool passesPionInRedVolume, passesNoOutgoingPion, passesSmallTracksCut, passesMeanCurvatureCut;
    tree->SetBranchAddress("passesPionInRedVolume", &passesPionInRedVolume);
    tree->SetBranchAddress("passesNoOutgoingPion", &passesNoOutgoingPion);
    tree->SetBranchAddress("passesSmallTracksCut", &passesSmallTracksCut);
    tree->SetBranchAddress("passesMeanCurvatureCut", &passesMeanCurvatureCut);

    // WC match information
    int WC2TPCtrkID;
    double WC2TPCPrimaryBeginX, WC2TPCPrimaryBeginY, WC2TPCPrimaryBeginZ;
    double WC2TPCPrimaryEndX, WC2TPCPrimaryEndY, WC2TPCPrimaryEndZ;
    std::vector<double>* wcMatchResR = nullptr;
    std::vector<double>* wcMatchDEDX = nullptr;
    std::vector<double>* wcMatchXPos = nullptr;
    std::vector<double>* wcMatchYPos = nullptr;
    std::vector<double>* wcMatchZPos = nullptr;
    tree->SetBranchAddress("WC2TPCtrkID", &WC2TPCtrkID);
    tree->SetBranchAddress("WC2TPCPrimaryBeginX", &WC2TPCPrimaryBeginX);
    tree->SetBranchAddress("WC2TPCPrimaryBeginY", &WC2TPCPrimaryBeginY);
    tree->SetBranchAddress("WC2TPCPrimaryBeginZ", &WC2TPCPrimaryBeginZ);
    tree->SetBranchAddress("WC2TPCPrimaryEndX", &WC2TPCPrimaryEndX);
    tree->SetBranchAddress("WC2TPCPrimaryEndY", &WC2TPCPrimaryEndY);
    tree->SetBranchAddress("WC2TPCPrimaryEndZ", &WC2TPCPrimaryEndZ);
    tree->SetBranchAddress("wcMatchResR", &wcMatchResR);
    tree->SetBranchAddress("wcMatchDEDX", &wcMatchDEDX);
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

    // Retrieve histograms
    TH1D* hTotalEvents = (TH1D*) Directory->Get<TH1D>("hTotalEvents");

    ///////////////////////
    // Create histograms //
    ///////////////////////

    TH1D *hRecoAbsorption   = new TH1D("hRecoAbsorption", "hRecoAbsorption;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D *hRecoAbsorption0p = new TH1D("hRecoAbsorption0p", "hRecoAbsorption0p;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D *hRecoAbsorptionNp = new TH1D("hRecoAbsorptionNp", "hRecoAbsorptionNp;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);

    TH1D *hStitchedDistanceFromVertex         = new TH1D("hStitchedDistanceFromVertex", "StitchedDistanceFromVertex;;", 20, 0, 70);
    TH1D *hStitchedOriginalDistanceFromVertex = new TH1D("hStitchedOriginalDistanceFromVertex", "hStitchedOriginalDistanceFromVertex;;", 20, 0, 70);

    TH1D *hStitchAsPionAndProton = new TH1D("hStitchAsPionAndProton", "hStitchAsPionAndProton;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D *hStitchAsPionBraggPeak = new TH1D("hStitchAsPionBraggPeak", "hStitchAsPionBraggPeak;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D *hStitchAsMIP  = new TH1D("hStitchAsMIP", "hStitchAsMIP;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D *hStitchAsProton        = new TH1D("hStitchAsProton", "hStitchAsProton;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D *hStitchFracBreakPoints = new TH1D("hStitchFracBreakPoints", "hStitchFracBreakPoints;;", 50, 0, 1);

    TH1D *hSecondaryPionChi2Protons = new TH1D("hSecondaryPionChi2Protons", "hSecondaryPionChi2Protons;;", 20, 0, 10);
    TH1D *hSecondaryPionChi2Pions   = new TH1D("hSecondaryPionChi2Pions", "hSecondaryPionChi2Pions;;", 20, 0, 10);
    TH1D *hSecondaryPionChi2Others  = new TH1D("hSecondaryPionChi2Others", "hSecondaryPionChi2Others;;", 20, 0, 10);

    TH1D *hSecondaryProtonChi2Protons = new TH1D("hSecondaryProtonChi2Protons", "hSecondaryProtonChi2Protons;;", 20, 0, 10);
    TH1D *hSecondaryProtonChi2Pions   = new TH1D("hSecondaryProtonChi2Pions", "hSecondaryProtonChi2Pions;;", 20, 0, 10);
    TH1D *hSecondaryProtonChi2Others  = new TH1D("hSecondaryProtonChi2Others", "hSecondaryProtonChi2Others;;", 20, 0, 10);

    TH1D *hSecondaryMeanDEDXProtons = new TH1D("hSecondaryMeanDEDXProtons", "hSecondaryMeanDEDXProtons;;", 10, 0, 10);
    TH1D *hSecondaryMeanDEDXPions   = new TH1D("hSecondaryMeanDEDXPions", "hSecondaryMeanDEDXPions;;", 10, 0, 10);
    TH1D *hSecondaryMeanDEDXOthers  = new TH1D("hSecondaryMeanDEDXOthers", "hSecondaryMeanDEDXOthers;;", 10, 0, 10);

    TH1D *hTotalBackgroundScatteringAngle = new TH1D("hTotalBackgroundScatteringAngle", "hTotalBackgroundScatteringAngle;;", 20, 0, TMath::Pi());
    TH1D *h0pBackgroundScatteringAngle    = new TH1D("h0pBackgroundScatteringAngle", "h0pBackgroundScatteringAngle;;", 20, 0, TMath::Pi());
    TH1D *hNpBackgroundScatteringAngle    = new TH1D("hNpBackgroundScatteringAngle", "hNpBackgroundScatteringAngle;;", 20, 0, TMath::Pi());

    TH1D *hTotalBackgroundScatteringLength = new TH1D("hTotalBackgroundScatteringLength", "hTotalBackgroundScatteringLength;;", 15, 0, 40);
    TH1D *h0pBackgroundScatteringLength    = new TH1D("h0pBackgroundScatteringLength", "h0pBackgroundScatteringLength;;", 15, 0, 40);
    TH1D *hNpBackgroundScatteringLength    = new TH1D("hNpBackgroundScatteringLength", "hNpBackgroundScatteringLength;;", 15, 0, 40);

    TH1D *hTotalBackgroundScatteringKE = new TH1D("hTotalBackgroundScatteringKE", "hTotalBackgroundScatteringKE;;", 20, 0, 0.5);
    TH1D *h0pBackgroundScatteringKE    = new TH1D("h0pBackgroundScatteringKE", "h0pBackgroundScatteringKE;;", 20, 0, 0.5);
    TH1D *hNpBackgroundScatteringKE    = new TH1D("hNpBackgroundScatteringKE", "hNpBackgroundScatteringKE;;", 20, 0, 0.5);

    TH1D *hScatteringBackgroundOriginalDistanceToPrimaryVertex = new TH1D("hScatteringBackgroundOriginalDistanceToPrimaryVertex", "hScatteringBackgroundOriginalDistanceToPrimaryVertex;;", 20, 0, 20);
    TH1D *hScatteringBackgroundOriginalDistanceToSecondaryVertex = new TH1D("hScatteringBackgroundOriginalDistanceToSecondaryVertex", "hScatteringBackgroundOriginalDistanceToSecondaryVertex;;", 20, 0, 20);
    TH1D *hScatteringBackgroundDistanceToPrimaryVertex = new TH1D("hScatteringBackgroundDistanceToPrimaryVertex", "hScatteringBackgroundDistanceToPrimaryVertex;;", 20, 0, 20);
    TH1D *hScatteringBackgroundDistanceToSecondaryVertex = new TH1D("hScatteringBackgroundDistanceToSecondaryVertex", "hScatteringBackgroundDistanceToSecondaryVertex;;", 20, 0, 20);

    TH1D *h0pInelasticBackgroundSecondaryInteraction = new TH1D("h0pInelasticBackgroundSecondaryInteraction", "h0pInelasticBackgroundSecondaryInteraction;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D *hNpInelasticBackgroundSecondaryInteraction = new TH1D("hNpInelasticBackgroundSecondaryInteraction", "hNpInelasticBackgroundSecondaryInteraction;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);

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

    TH2D *hTotalBackgroundScatteringLengthVSAngle = new TH2D(
        "hTotalBackgroundScatteringLengthVSAngle",
        "hTotalBackgroundScatteringLengthVSAngle;Length (cm);Angle (rad)",
        10, 0, 40, 
        10, 0, TMath::Pi()
    );

    TH2D *hTotalBackgroundScatteringLengthVSKEnergy = new TH2D(
        "hTotalBackgroundScatteringLengthVSKEnergy",
        "hTotalBackgroundScatteringLengthVSKEnergy;Length (cm);Energy (GeV/c)",
        10, 0, 40,
        8, 0, 0.4
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
    int NUM_CLUSTERS_NUM_STEPS = 3;

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

    int    NUM_CLUSTERS_THRESHOLD  = 1;
    double LARGE_CLUSTER_THRESHOLD = HIT_WIRE_SEPARATION * 3;

    //////////////////////////////////////////
    // Data for inelastic scattering events //
    //////////////////////////////////////////

    std::vector<double> InelasticScatteringIncidentKE;
    std::vector<double> InelasticScatteringVertexKE;
    std::vector<double> InelasticScatteringOutgoingKE;

    std::vector<double> Background0pInelasticScatteringIncidentKE;
    std::vector<double> Background0pInelasticScatteringVertexKE;
    std::vector<double> Background0pInelasticScatteringOutgoingKE;

    std::vector<double> BackgroundNpInelasticScatteringIncidentKE;
    std::vector<double> BackgroundNpInelasticScatteringVertexKE;
    std::vector<double> BackgroundNpInelasticScatteringOutgoingKE;

    TH1D *hInelasticScatteringReconstructed = new TH1D("hInelasticScatteringReconstructionEfficiency", "hInelasticScatteringReconstructionEfficiency;;", 20, 0, 0.4);
    TH1D *hInelasticScatteringTotal         = new TH1D("hInelasticScatteringTotal", "hInelasticScatteringTotal;;", 20, 0, 0.4);

    TH1D *hInelasticScatteringReconstructed0pBkg = new TH1D("hInelasticScatteringReconstructionEfficiency0pBkg", "hInelasticScatteringReconstructionEfficiency0pBkg;;", 20, 0, 0.4);
    TH1D *hInelasticScatteringTotal0pBkg         = new TH1D("hInelasticScatteringTotal0pBkg", "hInelasticScatteringTotal0pBkg;;", 20, 0, 0.4);

    TH1D *hInelasticScatteringReconstructedNpBkg = new TH1D("hInelasticScatteringReconstructionEfficiencyNpBkg", "hInelasticScatteringReconstructionEfficiencyNpBkg;;", 20, 0, 0.4);
    TH1D *hInelasticScatteringTotalNpBkg         = new TH1D("hInelasticScatteringTotalNpBkg", "hInelasticScatteringTotalNpBkg;;", 20, 0, 0.4);

    /////////////////////////////////
    // Files for event information //
    /////////////////////////////////

    std::ofstream outFile0pTrue("files/RecoAllAnalysis/0pRecoTrue.txt");
    std::ofstream outFileNpTrue("files/RecoAllAnalysis/NpRecoTrue.txt");
    std::ofstream outFile0pBackground("files/RecoAllAnalysis/0pBackground.txt");
    std::ofstream outFileNpBackground("files/RecoAllAnalysis/NpBackground.txt");
    std::ofstream outStitchedFile("files/RecoAllAnalysis/WCMatchStitching.txt");
    std::ofstream outWCAll("files/RecoAllAnalysis/WCMatchAllChi2.txt");

    //////////////////////
    // Loop over events //
    //////////////////////

    Int_t NumEntries = (Int_t) tree->GetEntries();
    std::cout << "Num entries: " << NumEntries << std::endl;

    for (Int_t i = 0; i < NumEntries; ++i) {
        tree->GetEntry(i);

        // Label background type as 0 for 0p signal and 1 for Np signal
        if (isPionAbsorptionSignal) {
            if (numVisibleProtons == 0) backgroundType = 0;
            if (numVisibleProtons > 0)  backgroundType = 1;
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
            Overlay_dEdx_RR_Reference_PP(gProton, gPion);
            c1->SaveAs(SaveDir +  "dEdxProfiles/" + event + ".png");
            
            delete c1;
            delete hDEDXProfile;
        }

        // Study ALL inelastic scattering events
        if (backgroundType == 6) {
            InelasticScatteringIncidentKE.push_back(truthPrimaryIncidentKE);
            InelasticScatteringVertexKE.push_back(truthPrimaryVertexKE);
            InelasticScatteringOutgoingKE.push_back(truthScatteredPionKE);

            hInelasticScatteringTotal->Fill(truthScatteredPionKE);
            for (int iRecoTrk = 0; iRecoTrk < matchedIdentity->size(); ++iRecoTrk) {
                if (
                    (matchedIdentity->at(iRecoTrk) == -211) && 
                    (matchedTrkID->at(iRecoTrk) != WC2TPCtrkID) && 
                    (matchedProcess->at(iRecoTrk) == "pi-Inelastic")
                ) {
                    hInelasticScatteringReconstructed->Fill(truthScatteredPionKE);
                    break;
                }
            }
        }

        // If no track matched to wire-chamber, skip
        if (WC2TPCtrkID == -99999) continue;

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

        // Apply small tracks cut
        if (!passesSmallTracksCut)  continue;

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
            } else {
                // TODO: if track is not near vertex
                //       First, check at least one of the two endpoints is inside reduced volume
                //       Then, we want to check that this track is not a pion
                //       If pion, reject event
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
            // If any secondary tagged as pion, not signal
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
            if ((hitsInTracks.count(iHit) > 0) || (fHitPlane->at(iHit) != 0)) continue;

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

        // std::cout << "Candidate hits: " << nCandidateHits << std::endl;
        // std::cout << "Vertex hit X: " << primaryEndPointHitX << "  W: " << primaryEndPointHitW << std::endl;
        
        for (int iHit = 0; iHit < nCandidateHits; ++iHit) {
            // std::cout << "  X: " << fHitX->at(candidateInductionHits.at(iHit));
            // std::cout << "  W: " << fHitW->at(candidateInductionHits.at(iHit)) << std::endl;

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
                if (usedHits.count(iAllHit) || hitsInTracks.count(iAllHit) || (fHitPlane->at(iHit) != 0)) continue;
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
        float totalClusterSize        = 0.;
        int   numLargeClusters        = 0;
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

        if (backgroundType == 6) {
            hTotalBackgroundScatteringAngle->Fill(truthScatteringAngle);
            hTotalBackgroundScatteringLength->Fill(truthScatteredPionLength);
            hTotalBackgroundScatteringKE->Fill(truthScatteredPionKE);

            hTotalBackgroundScatteringLengthVSAngle->Fill(truthScatteredPionLength, truthScatteringAngle);
            hTotalBackgroundScatteringLengthVSKEnergy->Fill(truthScatteredPionLength, truthScatteredPionKE);

            // Compute distance of reco'ed vertex to secondary vertex
            double distanceFromSecondaryVertex         = distance(breakPointX, truthSecondaryVertexX, breakPointY, truthSecondaryVertexY, breakPointZ, truthSecondaryVertexZ);
            double originalDistanceFromSecondaryVertex = distance(WC2TPCPrimaryEndX, truthSecondaryVertexX, WC2TPCPrimaryEndY, truthSecondaryVertexY, WC2TPCPrimaryEndZ, truthSecondaryVertexZ);

            hScatteringBackgroundOriginalDistanceToPrimaryVertex->Fill(originalDistanceFromVertex);
            hScatteringBackgroundOriginalDistanceToSecondaryVertex->Fill(originalDistanceFromSecondaryVertex);
            hScatteringBackgroundDistanceToPrimaryVertex->Fill(distanceFromVertex);
            hScatteringBackgroundDistanceToSecondaryVertex->Fill(distanceFromSecondaryVertex);

            if (totalTaggedProtons == 0) {
                Background0pInelasticScatteringIncidentKE.push_back(truthPrimaryIncidentKE);
                Background0pInelasticScatteringVertexKE.push_back(truthPrimaryVertexKE);
                Background0pInelasticScatteringOutgoingKE.push_back(truthScatteredPionKE);

                hInelasticScatteringTotal0pBkg->Fill(truthScatteredPionKE);
                for (int iRecoTrk = 0; iRecoTrk < matchedIdentity->size(); ++iRecoTrk) {
                    if ((matchedIdentity->at(iRecoTrk) == -211) && (matchedTrkID->at(iRecoTrk) != WC2TPCtrkID) && (matchedProcess->at(iRecoTrk) == "pi-Inelastic")) {
                        hInelasticScatteringReconstructed0pBkg->Fill(truthScatteredPionKE);
                        break;
                    }
                }

            } else if (totalTaggedProtons > 0) {
                BackgroundNpInelasticScatteringIncidentKE.push_back(truthPrimaryIncidentKE);
                BackgroundNpInelasticScatteringVertexKE.push_back(truthPrimaryVertexKE);
                BackgroundNpInelasticScatteringOutgoingKE.push_back(truthScatteredPionKE);

                hInelasticScatteringTotalNpBkg->Fill(truthScatteredPionKE);
                for (int iRecoTrk = 0; iRecoTrk < matchedIdentity->size(); ++iRecoTrk) {
                    if ((matchedIdentity->at(iRecoTrk) == -211) && (matchedTrkID->at(iRecoTrk) != WC2TPCtrkID) && (matchedProcess->at(iRecoTrk) == "pi-Inelastic")) {
                        hInelasticScatteringReconstructedNpBkg->Fill(truthScatteredPionKE);
                        break;
                    }
                }
            }
        }

        if (totalTaggedProtons == 0) {
            if (backgroundType == 0) {
                printEventInfo(thisEventInfo, outFile0pTrue);
            } else if (backgroundType == 6) {
                printEventInfo(thisEventInfo, outFile0pBackground);
                h0pInelasticBackgroundSecondaryInteraction->Fill(secondaryInteractionTag);
                h0pBackgroundScatteringAngle->Fill(truthScatteringAngle);
                h0pBackgroundScatteringLength->Fill(truthScatteredPionLength);
                h0pBackgroundScatteringKE->Fill(truthScatteredPionKE);
            } else { printEventInfo(thisEventInfo, outFile0pBackground); }
        } else if (totalTaggedProtons > 0) {
            if (backgroundType == 1) {
                printEventInfo(thisEventInfo, outFileNpTrue);
            } else if (backgroundType == 6) {
                printEventInfo(thisEventInfo, outFileNpBackground);
                hNpInelasticBackgroundSecondaryInteraction->Fill(secondaryInteractionTag);
                hNpBackgroundScatteringAngle->Fill(truthScatteringAngle);
                hNpBackgroundScatteringLength->Fill(truthScatteredPionLength);
                hNpBackgroundScatteringKE->Fill(truthScatteredPionKE);
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

    ////////////////////////////////////////////
    // Compute reco efficiency for scattering //
    ////////////////////////////////////////////

    TEfficiency* hInelasticScatteringReconstructionEfficiency = nullptr;
    hInelasticScatteringReconstructionEfficiency = new TEfficiency(*hInelasticScatteringReconstructed, *hInelasticScatteringTotal);

    TEfficiency* hInelasticScatteringReconstructionEfficiency0pBkg = nullptr;
    hInelasticScatteringReconstructionEfficiency0pBkg = new TEfficiency(*hInelasticScatteringReconstructed0pBkg, *hInelasticScatteringTotal0pBkg);

    TEfficiency* hInelasticScatteringReconstructionEfficiencyNpBkg = nullptr;
    hInelasticScatteringReconstructionEfficiencyNpBkg = new TEfficiency(*hInelasticScatteringReconstructedNpBkg, *hInelasticScatteringTotalNpBkg);

    std::vector<TEfficiency*> EfficiencyPlots = {
        hInelasticScatteringReconstructionEfficiency,
        hInelasticScatteringReconstructionEfficiency0pBkg,
        hInelasticScatteringReconstructionEfficiencyNpBkg
    };

    std::vector<TString> EfficiencyPlotsTitles = {
        "RecoEfficiencyScatteredPions",
        "RecoEfficiencyScatteredPions0pBkg",
        "RecoEfficiencyScatteredPionsNpBkg"
    };

    std::vector<TString> EfficiencyPlotsXLabels = {
        "Outgoing pion KE (GeV/c)",
        "Outgoing pion KE (GeV/c)",
        "Outgoing pion KE (GeV/c)"
    };

    printEfficiencyPlots(
        SaveDir, FontStyle, TextSize,
        EfficiencyPlots,
        EfficiencyPlotsTitles,
        EfficiencyPlotsXLabels
    );

    //////////////////
    // Create plots //
    //////////////////

    std::vector<int> Colors = {
        kBlack, kBlue, kRed, kGreen
    };

    std::vector<std::vector<TH1*>> PlotGroups = {
        {hStitchedDistanceFromVertex, hStitchedOriginalDistanceFromVertex},
        {hStitchAsPionAndProton, hStitchAsPionBraggPeak, hStitchAsMIP, hStitchAsProton},
        {hStitchFracBreakPoints},
        {hSecondaryPionChi2Pions, hSecondaryPionChi2Protons, hSecondaryPionChi2Others},
        {hSecondaryProtonChi2Pions, hSecondaryProtonChi2Protons, hSecondaryProtonChi2Others},
        {hSecondaryMeanDEDXPions, hSecondaryMeanDEDXProtons, hSecondaryMeanDEDXOthers},
        {hTotalBackgroundScatteringAngle, h0pBackgroundScatteringAngle, hNpBackgroundScatteringAngle},
        {hTotalBackgroundScatteringLength, h0pBackgroundScatteringLength, hNpBackgroundScatteringLength},
        {hTotalBackgroundScatteringKE, h0pBackgroundScatteringKE, hNpBackgroundScatteringKE},
        {hScatteringBackgroundOriginalDistanceToPrimaryVertex, hScatteringBackgroundOriginalDistanceToSecondaryVertex, hScatteringBackgroundDistanceToPrimaryVertex, hScatteringBackgroundDistanceToSecondaryVertex},
        {hHitClusters0p, hHitClusters0pBackground, hHitClustersNp, hHitClustersNpBackground},
        {hHitClustersSize0p, hHitClustersSize0pBackground, hHitClustersSizeNp, hHitClustersSizeNpBackground},
        {hHitLargeClusters0p, hHitLargeClusters0pBackground, hHitLargeClustersNp, hHitLargeClustersNpBackground},
        {hInelasticScatteringTotal, hInelasticScatteringReconstructed},
        {hInelasticScatteringTotal0pBkg, hInelasticScatteringReconstructed0pBkg},
        {hInelasticScatteringTotalNpBkg, hInelasticScatteringReconstructedNpBkg},
        {hMinimumLinearity0p, hMinimumLinearity0pBackground},
        {hStdDevLinearity0p, hStdDevLinearity0pBackground},
        {hMinimumLinearityNp, hMinimumLinearityNpBackground},
        {hStdDevLinearityNp, hStdDevLinearityNpBackground},
        {hMaxLinearityD0p, hMaxLinearityD0pBackground},
        {hMaxLinearityDD0p, hMaxLinearityDD0pBackground},
        {hLocalLinearityDerivativeCutPur, hLocalLinearityDerivativeCutEff}
    };

    std::vector<std::vector<TString>> PlotLabelGroups = {
        {"Detected vertex distance", "Original distance"},
        {"Stitched MIP-proton", "Bragg pion", "MIP", "Proton"},
        {"Stitched tracks"},
        {"Pions", "Protons", "Others"},
        {"Pions", "Protons", "Others"},
        {"Pions", "Protons", "Others"},
        {"All background", "0p background", "Np background"},
        {"All background", "0p background", "Np background"},
        {"All background", "0p background", "Np background"},
        {"Original to primary", "Original to secondary", "Detected to primary", "Detected to secondary"},
        {"Reco 0p true", "Reco 0p background", "Reco Np true", "Reco Np background"},
        {"Reco 0p true", "Reco 0p background", "Reco Np true", "Reco Np background"},
        {"Reco 0p true", "Reco 0p background", "Reco Np true", "Reco Np background"},
        {"All", "Reconstructed"},
        {"All", "Reconstructed"},
        {"All", "Reconstructed"},
        {"Reco 0p true", "Reco 0p background"},
        {"Reco 0p true", "Reco 0p background"},
        {"Reco 0p true", "Reco 0p background"},
        {"Reco Np true", "Reco Np background"},
        {"Reco 0p true", "Reco 0p background"},
        {"Reco 0p true", "Reco 0p background"},
        {"Purity", "Efficiency"}
    };

    std::vector<TString> PlotTitles = {
        "StitchedDistanceFromVertex",
        "StitchedBackgroundTypes",
        "StitchedFracBreakPoints",
        "PionChi2SecondaryTracks",
        "ProtonChi2SecondaryTracks",
        "MeanDEDXSecondaryTracks",
        "InelasticScatteringBackgroundAngle",
        "InelasticScatteringBackgroundLength",
        "InelasticScatteringBackgroundKE",
        "InelasticScatteringBackgroundDistanceFromVertex",
        "NumberOfHitClusters",
        "AverageHitClusterSize",
        "NumberOfLargeHitClusters",
        "RecoScatteredPions",
        "RecoScatteredPions0pBkg",
        "RecoScatteredPionsNpBkg",
        "MinLocalLinearity0pReco",
        "StdDevLinearity0pReco",
        "MinLocalLinearityNpReco",
        "StdDevLinearityNpReco",
        "MaxLocalLinearityDeriv0pReco",
        "MaxLocalLinearitySecondDeriv0pReco",
        "MaxLocalLinearityDerivCut"
    };

    std::vector<TString> XLabels = {
        "Distance from true vertex (cm)",
        "Background type",
        "Fractional break point",
        "Pion reduced #chi^{2}",
        "Proton reduced #chi^{2}",
        "Mean dE/dx (MeV/cm)",
        "Scattering angle (rad)",
        "Scattered pion length (cm)",
        "Kinetic energy (MeV/c)",
        "Distance from true vertex (cm)",
        "Number of induction plane hit clusters",
        "Induction plane hit cluster average size",
        "Number of induction plane large hit clusters",
        "Outgoing pion energy (GeV/c)",
        "Outgoing pion energy (GeV/c)",
        "Outgoing pion energy (GeV/c)",
        "Minimum local linearity",
        "Standard dev. local linearity",
        "Minimum local linearity",
        "Standard dev. local linearity",
        "Maximum local linearity derivative",
        "Maximum local linearity second derivative",
        "Maximum local linearity derivative allowed"
    };

    std::vector<TString> YLabels = {
        "Number of tracks",
        "Number of tracks",
        "Number of tracks",
        "Number of tracks",
        "Number of tracks",
        "Number of tracks",
        "Number of events",
        "Number of tracks",
        "Number of tracks",
        "Number of tracks",
        "Number of events",
        "Number of events",
        "Number of events",
        "Number of events",
        "Number of events",
        "Number of events",
        "Number of tracks",
        "Number of tracks",
        "Number of tracks",
        "Number of tracks",
        "Number of tracks",
        "Number of tracks",
        "Purity/Efficiency"
    };

    printOneDPlots(
        SaveDir, FontStyle, TextSize,
        PlotGroups,
        Colors,
        PlotLabelGroups,
        PlotTitles,
        XLabels,
        YLabels
    );

    std::vector<TH2*> TwoDPlots = {
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
        hTotalBackgroundScatteringLengthVSAngle,
        hTotalBackgroundScatteringLengthVSKEnergy,
        hHitClusterCut0pReco,
        hHitClusterCut0pRecoTrue,
        hHitClusterCut0pPurity
    };

    std::vector<TString> TwoDTitles = {
        "ProtonChi2TruePositives",
        "ProtonChi2FalsePositives",
        "ProtonChi2TrueNegatives",
        "ProtonChi2FalseNegatives",
        "PionChi2TruePositives",
        "PionChi2FalsePositives",
        "PionChi2TrueNegatives",
        "PionChi2FalseNegatives",
        "ProtonChi2FOM",
        "PionChi2FOM",
        "TotalBackgroundScatteringLengthVSAngle",
        "TotalBackgroundScatteringLengthVSKEnergy",
        "HitClusterCut0pReco",
        "HitClusterCut0pRecoTrue",
        "HitClusterCut0pPurity"
    };

    printTwoDPlots(SaveDir, TwoDPlots, TwoDTitles);

    //////////////////////////
    // Create stacked plots //
    //////////////////////////

    THStack* hBackgroundTypesStack = new THStack("hBackgroundTypesStack", "BackgroundTypesStack");

    double alpha = 0.2;
    hStitchAsPionAndProton->SetTitle("Stitched MIP-proton");
    hStitchAsPionBraggPeak->SetTitle("Bragg pion");
    hStitchAsMIP->SetTitle("MIP");
    hStitchAsProton->SetTitle("Bragg proton");

    hStitchAsPionAndProton->SetFillColorAlpha(Colors[0], alpha);
    hStitchAsPionBraggPeak->SetFillColorAlpha(Colors[1], alpha);
    hStitchAsMIP->SetFillColorAlpha(Colors[2], alpha);
    hStitchAsProton->SetFillColorAlpha(Colors[3], alpha);

    hStitchAsPionAndProton->SetLineColor(Colors[0]);
    hStitchAsPionBraggPeak->SetLineColor(Colors[1]);
    hStitchAsMIP->SetLineColor(Colors[2]);
    hStitchAsProton->SetLineColor(Colors[3]);

    hBackgroundTypesStack->Add(hStitchAsMIP);
    hBackgroundTypesStack->Add(hStitchAsPionAndProton);
    hBackgroundTypesStack->Add(hStitchAsPionBraggPeak);
    hBackgroundTypesStack->Add(hStitchAsProton);

    TCanvas* c = new TCanvas("Canvas", "Canvas", 205, 34, 1024, 768);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);
    hBackgroundTypesStack->Draw("HIST");
    c->BuildLegend(0.65, 0.7, 0.95, 0.85);

    hBackgroundTypesStack->GetXaxis()->SetTitleFont(FontStyle);
    hBackgroundTypesStack->GetXaxis()->SetLabelFont(FontStyle);
    hBackgroundTypesStack->GetXaxis()->SetNdivisions(8);
    hBackgroundTypesStack->GetXaxis()->SetLabelSize(TextSize);
    hBackgroundTypesStack->GetXaxis()->SetTitleSize(TextSize);
    // hBackgroundTypesStack->GetXaxis()->SetTitle("Background type");
    hBackgroundTypesStack->GetXaxis()->SetTitleOffset(1.1);
    hBackgroundTypesStack->GetXaxis()->CenterTitle();

    hBackgroundTypesStack->GetYaxis()->SetTitleFont(FontStyle);
    hBackgroundTypesStack->GetYaxis()->SetLabelFont(FontStyle);
    hBackgroundTypesStack->GetYaxis()->SetNdivisions(6);
    hBackgroundTypesStack->GetYaxis()->SetLabelSize(TextSize);
    hBackgroundTypesStack->GetYaxis()->SetTitleSize(TextSize);
    hBackgroundTypesStack->GetYaxis()->SetTitle("Number of tracks");
    hBackgroundTypesStack->GetYaxis()->SetTitleOffset(1.1);
    hBackgroundTypesStack->GetYaxis()->CenterTitle();

    for (const auto& entry : backgroundTypes) {
        int bin = entry.first;
        const std::string& label = entry.second;
        if ((bin+1 == 0) || (bin+1 == 12)) continue; // extra bins, hacky
        hBackgroundTypesStack->GetXaxis()->SetBinLabel(bin + 1, label.c_str());
    }

    c->Update();
    c->SaveAs(SaveDir + "StitchedBackgroundTypesStacked.png");

    //////////////////////////////////////////////////////
    // Create scatter plots for inelastic scattering KE //
    //////////////////////////////////////////////////////

    int nDataPoints               = InelasticScatteringIncidentKE.size();
    TGraph *gIncidentVSOutgoingKE = new TGraph(nDataPoints, InelasticScatteringIncidentKE.data(), InelasticScatteringOutgoingKE.data());
    TGraph *gVertexVSOutgoingKE   = new TGraph(nDataPoints, InelasticScatteringVertexKE.data(), InelasticScatteringOutgoingKE.data());

    gIncidentVSOutgoingKE->SetTitle("Outgoing KE vs Incident KE;Incident KE [MeV];Outgoing KE [MeV]");
    gIncidentVSOutgoingKE->SetMarkerStyle(20);
    gIncidentVSOutgoingKE->SetMarkerColor(kBlack);
    gIncidentVSOutgoingKE->Draw("AP");
    c->SaveAs(SaveDir + "InelasticScatteringIncidentVSOutgoingKE.png");

    gVertexVSOutgoingKE->SetTitle("Outgoing KE vs Vertex KE;Vertex KE [MeV];Outgoing KE [MeV]");
    gVertexVSOutgoingKE->SetMarkerStyle(20);
    gVertexVSOutgoingKE->SetMarkerColor(kBlack);
    gVertexVSOutgoingKE->Draw("AP");
    c->SaveAs(SaveDir + "InelasticScatteringVertexVSOutgoingKE.png");

    int n0pDataPoints               = Background0pInelasticScatteringIncidentKE.size();
    TGraph *g0pIncidentVSOutgoingKE = new TGraph(n0pDataPoints, Background0pInelasticScatteringIncidentKE.data(), Background0pInelasticScatteringOutgoingKE.data());
    TGraph *g0pVertexVSOutgoingKE   = new TGraph(n0pDataPoints, Background0pInelasticScatteringVertexKE.data(), Background0pInelasticScatteringOutgoingKE.data());

    g0pIncidentVSOutgoingKE->SetTitle("Outgoing KE vs Incident KE;Incident KE [MeV];Outgoing KE [MeV]");
    g0pIncidentVSOutgoingKE->SetMarkerStyle(20);
    g0pIncidentVSOutgoingKE->SetMarkerColor(kBlack);
    g0pIncidentVSOutgoingKE->Draw("AP");
    c->SaveAs(SaveDir + "Background0pInelasticScatteringIncidentVSOutgoingKE.png");

    g0pVertexVSOutgoingKE->SetTitle("Outgoing KE vs Vertex KE;Vertex KE [MeV];Outgoing KE [MeV]");
    g0pVertexVSOutgoingKE->SetMarkerStyle(20);
    g0pVertexVSOutgoingKE->SetMarkerColor(kBlack);
    g0pVertexVSOutgoingKE->Draw("AP");
    c->SaveAs(SaveDir + "Background0pInelasticScatteringVertexVSOutgoingKE.png");

    int nNpDataPoints               = BackgroundNpInelasticScatteringIncidentKE.size();
    TGraph *gNpIncidentVSOutgoingKE = new TGraph(nNpDataPoints, BackgroundNpInelasticScatteringIncidentKE.data(), BackgroundNpInelasticScatteringOutgoingKE.data());
    TGraph *gNpVertexVSOutgoingKE   = new TGraph(nNpDataPoints, BackgroundNpInelasticScatteringVertexKE.data(), BackgroundNpInelasticScatteringOutgoingKE.data());

    gNpIncidentVSOutgoingKE->SetTitle("Outgoing KE vs Incident KE;Incident KE [MeV];Outgoing KE [MeV]");
    gNpIncidentVSOutgoingKE->SetMarkerStyle(20);
    gNpIncidentVSOutgoingKE->SetMarkerColor(kBlack);
    gNpIncidentVSOutgoingKE->Draw("AP");
    c->SaveAs(SaveDir + "BackgroundNpInelasticScatteringIncidentVSOutgoingKE.png");

    gNpVertexVSOutgoingKE->SetTitle("Outgoing KE vs Vertex KE;Vertex KE [MeV];Outgoing KE [MeV]");
    gNpVertexVSOutgoingKE->SetMarkerStyle(20);
    gNpVertexVSOutgoingKE->SetMarkerColor(kBlack);
    gNpVertexVSOutgoingKE->Draw("AP");
    c->SaveAs(SaveDir + "BackgroundNpInelasticScatteringVertexVSOutgoingKE.png");
}

bool isHitNearPrimary(std::vector<int>* primaryKey, std::vector<float>* hitX, std::vector<float>* hitW, float thisHitX, float thisHitW, float xThreshold, float wThreshold) {
    int nPrimaryHits = primaryKey->size();
    for (int iHit = 0; iHit < nPrimaryHits; ++iHit) {
        float dX = std::abs(thisHitX - hitX->at(primaryKey->at(iHit)));
        float dW = std::abs(thisHitW - hitW->at(primaryKey->at(iHit)));
        if ((dX < xThreshold) && (dW < wThreshold)) return true;
    }
    return false;
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

double meanDEDX(std::vector<double> trackDEDX, bool isTrackReversed, int pointsToUse) {
    double dEdx = 0.;
    if (isTrackReversed) std::reverse(trackDEDX.begin(), trackDEDX.end());

    unsigned int bound = pointsToUse;
    if (pointsToUse > trackDEDX.size()) bound = trackDEDX.size();
    for (unsigned int i = 0; i < bound; ++i) dEdx += trackDEDX.at(i);
    dEdx /= bound;
    return dEdx;
}

int isSecondaryInteractionAbsorption(std::vector<int> daughtersPDG, std::vector<string> daughtersProcess, std::vector<double> daughtersKE) {
    int numDaughters = daughtersPDG.size();
    if (numDaughters == 0) return 11;

    int tempNumProtons = 0;
    for (int i = 0; i < numDaughters; ++i) {
        if ((daughtersPDG[i] == 11) && (daughtersProcess[i] == "hIoni")) continue;
        if (daughtersPDG[i] == -211) return 6;
        if (daughtersPDG[i] == 111) return 7;
        if (daughtersPDG[i] == 211) return 8;
        if (daughtersProcess[i] == "Decay") return 10;
        if (daughtersProcess[i] == "hBertiniCaptureAtRest") return 9;

        if (daughtersProcess[i] == "pi-Inelastic") {
            if ((daughtersPDG[i] == 13) || (daughtersPDG[i] == -13)) return 11;
            else if ((daughtersPDG[i] == 321) || (daughtersPDG[i] == -321) || (daughtersPDG[i] == 311)) return 11;
            else if (daughtersPDG[i] == 2212) {
                if ((daughtersKE[i] >= PROTON_ENERGY_LOWER_BOUND) && (daughtersKE[i] <= PROTON_ENERGY_UPPER_BOUND)) {
                    tempNumProtons++;
                }
            }
        }
    }

    if (tempNumProtons == 0) return 0;
    return 1;
}

std::vector<double> calcLinearityProfile(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, int nb) {
    size_t N = vx.size();
    std::vector<double> linearity(N, 0.0);

    if (N < 3 || vx.size() != vy.size() || vx.size() != vz.size()) return linearity;

    for (size_t index = 0; index < N; ++index) {
        int k1 = std::max(int(0), int(index - nb));
        int k2 = std::min(int(N - 1), int(index + nb));
        int M  = k2 - k1 + 1;

        double mx = 0., my = 0., mz = 0.;
        for (int i = k1; i <= k2; ++i) { mx += vx[i]; my += vy[i]; mz += vz[i]; }
        mx /= M; my /= M; mz /= M;

        double Cxx = 0., Cxy = 0., Cxz = 0.;
        double Cyy = 0., Cyz = 0., Czz = 0.;

        for (int i = k1; i <= k2; ++i) {
            double dx = vx[i] - mx;
            double dy = vy[i] - my;
            double dz = vz[i] - mz;

            Cxx += dx * dx; Cxy += dx * dy; Cxz += dx * dz;
            Cyy += dy * dy; Cyz += dy * dz;
            Czz += dz * dz;
        }

        Cxx /= M; Cxy /= M; Cxz /= M;
        Cyy /= M; Cyz /= M;
        Czz /= M;

        TMatrixDSym cov(3);
        cov(0,0) = Cxx; cov(0,1) = Cxy; cov(0,2) = Cxz;
        cov(1,0) = Cxy; cov(1,1) = Cyy; cov(1,2) = Cyz;
        cov(2,0) = Cxz; cov(2,1) = Cyz; cov(2,2) = Czz;

        TMatrixDSymEigen eigenSolver(cov);
        TVectorD evals = eigenSolver.GetEigenValues();
        std::vector<double> eigs = { evals[0], evals[1], evals[2] };

        double lambda_sum = eigs[0] + eigs[1] + eigs[2];
        if (lambda_sum == 0.0f) {
            linearity[index] = 1.0;
        } else {
            linearity[index] = eigs[0] / lambda_sum;
        }
    }
    return linearity;
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

void Overlay_dEdx_RR_Reference_PP(
    TGraph* gProton,
    TGraph* gPion,
    bool addLegend = true, 
    TVirtualPad* pad = gPad
) {
    const int N = 107;

    /* ---------- Proton ---------- */
    gProton->SetLineColor(kRed+1);
    gProton->SetLineWidth(3);
    gProton->SetName("Proton");

    /* ---------- Pion ---------- */
    gPion->SetLineColor(kRed+3);
    gPion->SetLineWidth(3);
    gPion->SetName("Pion");

    if (pad) pad->cd();
    gProton->Draw("L same");
    gPion  ->Draw("L same");
    
    if (addLegend) {
        static TLegend *leg = nullptr;
        if (!leg) {
            leg = new TLegend(0.50,0.80,0.75,0.88);
            leg->AddEntry(gProton, "Proton", "l");
            leg->AddEntry(gPion,   "Pion",   "l");
            leg->SetFillStyle(1001);
            leg->SetBorderSize(0);
        }
        leg->Draw();
    }
}


double distance(double x1, double x2, double y1, double y2, double z1, double z2) {
    return sqrt(
        pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)
    );
}

bool isWithinReducedVolume(double x, double y, double z) {
    return (
        (x > RminX) && (x < RmaxX) && 
        (y > RminY) && (y < RmaxY) && 
        (z > RminZ) && (z < RmaxZ)
    );
}

void printEfficiencyPlots(
    TString dir, int fontStyle, double textSize,
    std::vector<TEfficiency*> efficiencies,
    std::vector<TString> titles,
    std::vector<TString> xlabels
) {
    for (size_t iPlot = 0; iPlot < efficiencies.size(); ++iPlot) {
        TEfficiency* eff = efficiencies.at(iPlot);

        TCanvas* PlotCanvas = new TCanvas("Canvas","Canvas",205,34,1024,768);
        PlotCanvas->cd();
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetBottomMargin(0.15);

        eff->SetTitle(titles.at(iPlot) + ";" + xlabels.at(iPlot) + ";Efficiency");
        eff->SetLineWidth(2);
        eff->SetMarkerStyle(20);
        eff->SetMarkerSize(1.0);
        eff->Draw("AP");

        PlotCanvas->Update();
        TGraph* gr = eff->GetPaintedGraph();
        if (gr) {
            gr->GetXaxis()->SetTitleFont(fontStyle);
            gr->GetXaxis()->SetLabelFont(fontStyle);
            gr->GetXaxis()->SetNdivisions(8);
            gr->GetXaxis()->SetLabelSize(textSize);
            gr->GetXaxis()->SetTitleSize(textSize);
            gr->GetXaxis()->SetTitle(xlabels[iPlot]);
            gr->GetYaxis()->SetTitleOffset(1.1);
            gr->GetYaxis()->CenterTitle();

            gr->GetYaxis()->SetTitleFont(fontStyle);
            gr->GetYaxis()->SetLabelFont(fontStyle);
            gr->GetYaxis()->SetNdivisions(8);
            gr->GetYaxis()->SetLabelSize(textSize);
            gr->GetYaxis()->SetTitleSize(textSize);
            gr->GetYaxis()->SetTitle("Efficiency");
            gr->GetYaxis()->SetTitleOffset(1.1);
            gr->GetYaxis()->CenterTitle();
        }
        PlotCanvas->Update();

        // Save plot and delete canvas
        PlotCanvas->SaveAs(dir + titles.at(iPlot) + ".png");
        delete PlotCanvas;
    }
}

void printOneDPlots(
    TString dir, int fontStyle, double textSize, 
    std::vector<std::vector<TH1*>> groups,
    std::vector<int> colors,
    std::vector<std::vector<TString>> labels, 
    std::vector<TString> titles, 
    std::vector<TString> xlabels,
    std::vector<TString> ylabels
) {
    int numPlots = groups.size();
    for (int iPlot = 0; iPlot < numPlots; ++iPlot) {
        // Set up canvas
        TCanvas* PlotCanvas = new TCanvas("Canvas","Canvas",205,34,1024,768);
        PlotCanvas->cd();
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetBottomMargin(0.15);

        TLegend* leg = new TLegend(0.65,0.65,0.85,0.80);
        leg->SetBorderSize(0);
        leg->SetTextSize(textSize * 0.8);
        leg->SetTextFont(fontStyle);

        // Get histograms and labels
        std::vector<TH1*> Plots = groups.at(iPlot);
        std::vector<TString> Labels = labels.at(iPlot);

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

        for (int iSubPlot = 0; iSubPlot < Plots.size(); ++iSubPlot) {
            leg->AddEntry(Plots[iSubPlot], Labels[iSubPlot], "l");
            Plots[iSubPlot]->SetLineWidth(2);
            Plots[iSubPlot]->SetLineColor(colors.at(iSubPlot));
            Plots[iSubPlot]->Draw("hist same");

            double imax = TMath::Max(Plots[iSubPlot]->GetMaximum(), Plots[0]->GetMaximum());

            double YAxisRange = 1.15*imax;
            Plots[iSubPlot]->GetYaxis()->SetRangeUser(0., YAxisRange);
            Plots[0]->GetYaxis()->SetRangeUser(0., YAxisRange);	
        }

        leg->Draw();
        PlotCanvas->SaveAs(dir + titles.at(iPlot) + ".png");
        delete PlotCanvas;
    }
}

void printTwoDPlots(TString dir, std::vector<TH2*> plots, std::vector<TString> titles) {
    int nPlots = plots.size();

    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas("c1", "EnergyLossPlots", 800, 600);
    for (int iPlot = 0; iPlot < nPlots; ++iPlot) {
        TH1* hPlot = plots.at(iPlot);
        hPlot->SetMinimum(0);
        hPlot->SetMaximum(hPlot->GetMaximum());
        hPlot->Draw("COLZ");
        c1->SaveAs(dir + titles.at(iPlot) + ".png");
    }
}

void printBackgroundInfo(TH1D* background_histo, std::ostream& os) {
    int pionAbs0p            = background_histo->GetBinContent(1);
    int pionAbsNp            = background_histo->GetBinContent(2);
    int primaryMuon          = background_histo->GetBinContent(3);
    int primaryElectron      = background_histo->GetBinContent(4);
    int otherPrimary         = background_histo->GetBinContent(5);
    int pionOutRedVol        = background_histo->GetBinContent(6);
    int pionInelScatter      = background_histo->GetBinContent(7);
    int chargeExchange       = background_histo->GetBinContent(8);
    int doubleChargeExchange = background_histo->GetBinContent(9);
    int captureAtRest        = background_histo->GetBinContent(10);
    int pionDecay            = background_histo->GetBinContent(11);
    int other                = background_histo->GetBinContent(12);

    os << "  Abs 0p: " << pionAbs0p << std::endl;;
    os << "  Abs Np: " << pionAbsNp << std::endl;;
    os << "  Primary muon: " << primaryMuon << std::endl;;
    os << "  Primary electron: " << primaryElectron << std::endl;;
    os << "  Other primary: " << otherPrimary << std::endl;;
    os << "  Outside reduced volume: " << pionOutRedVol << std::endl;;
    os << "  Inelastic scattering: " << pionInelScatter << std::endl;;
    os << "  Charge exchange: " << chargeExchange << std::endl;;
    os << "  Double charge exchange: " << doubleChargeExchange << std::endl;;
    os << "  Capture at rest: " << captureAtRest << std::endl;;
    os << "  Decay: " << pionDecay << std::endl;;
    os << "  Other: " << other << std::endl;
}

void printEventInfo(EventInfo event, std::ostream& os) {
    os << "Run: " << event.run << " subrun: " << event.subrun << " event: " << event.event << std::endl;
    os << "  Background type: " << backgroundTypes[event.backgroundNum] << std::endl;
    os << "  WC match PDG: " << event.wcMatchPDG << std::endl;
    os << "  WC match process: " << event.wcMatchProcess << std::endl;
    os << "  WC match daughters: " << std::endl;
    for (int n = 0; n < event.wcMatchDaughtersPDG.size(); ++n) {
        if (event.wcMatchDaughtersPDG[n] == 11 && event.wcMatchDaughtersProcess[n] == "hIoni") continue;
        os << "      PDG: " << event.wcMatchDaughtersPDG[n];
        os << " Process: " << event.wcMatchDaughtersProcess[n];
        os << std::endl;
    }
    os << std::endl;
    os << "  Primary PDG: " << event.truthPrimaryPDG << std::endl;
    os << "  Primary vertex: " << event.vertexX << ", " << event.vertexY << ", " << event.vertexZ << std::endl; 
    os << "  Primary daughters: " << std::endl;
    for (int n = 0; n < event.truthPrimaryDaughtersPDG.size(); ++n) {
        if (event.truthPrimaryDaughtersPDG[n] == 11 && event.truthPrimaryDaughtersProcess[n] == "hIoni") continue;
        os << "      PDG: " << event.truthPrimaryDaughtersPDG[n];
        os << " Process: " << event.truthPrimaryDaughtersProcess[n];
        os << std::endl;
    }
    os << "  Secondary pion daughters: " << std::endl;
    for (int n = 0; n < event.truthSecondaryDaughtersPDG.size(); ++n) {
        if (event.truthSecondaryDaughtersPDG[n] == 11 && event.truthSecondaryDaughtersProcess[n] == "hIoni") continue;
        os << "      PDG: " << event.truthSecondaryDaughtersPDG[n];
        os << " Process: " << event.truthSecondaryDaughtersProcess[n];
        os << std::endl;
    }
    os << "  Secondary interaction tagged as: " << backgroundTypes[event.secondaryInteractionTag] << std::endl;
    os << std::endl;
    os << "  Tracks reco'ed near vertex: " << event.tracksNearVertex << std::endl;
    os << "  Tracks reco'ed as protons: " << event.recoProtonCount << std::endl;
    os << "  True visible protons: " << event.visibleProtons << std::endl;
    os << std::endl;
    os << "  Number of hit clusters found near primary track: " << event.numHitClusters << std::endl;
    os << std::endl;
    os << "  Minimum primary track linearity: " << event.minLocalLinearity << std::endl;
    os << "  Maximum primary track linearity derivative: " << event.maxLocalLinearityD << " above threshold: " << (event.maxLocalLinearityD >= LINEARITY_DERIVATIVE_THRESHOLD) << std::endl;
    os << std::endl;
}