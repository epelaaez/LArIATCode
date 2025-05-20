#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TVector3.h>
#include <TEfficiency.h>

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
    std::vector<int>*     matchedIdentity = nullptr;
    std::vector<double> * matchedKEnergy  = nullptr;
    tree->SetBranchAddress("matchedIdentity", &matchedIdentity);
    tree->SetBranchAddress("matchedKEnergy", &matchedKEnergy);

    // Calorimetry information
    std::vector<std::vector<double>>* recoResR = nullptr;
    std::vector<std::vector<double>>* recoDEDX = nullptr;
    tree->SetBranchAddress("recoResR", &recoResR);
    tree->SetBranchAddress("recoDEDX", &recoDEDX);

    // Truth-level information
    int truthPrimaryPDG;
    double truthPrimaryVertexX, truthPrimaryVertexY, truthPrimaryVertexZ;
    std::vector<int>*         truthPrimaryDaughtersPDG     = nullptr;
    std::vector<std::string>* truthPrimaryDaughtersProcess = nullptr;
    std::vector<double>*      truthPrimaryDaughtersKE      = nullptr;
    tree->SetBranchAddress("truthPrimaryPDG", &truthPrimaryPDG);
    tree->SetBranchAddress("truthPrimaryVertexX", &truthPrimaryVertexX);
    tree->SetBranchAddress("truthPrimaryVertexY", &truthPrimaryVertexY);
    tree->SetBranchAddress("truthPrimaryVertexZ", &truthPrimaryVertexZ);
    tree->SetBranchAddress("truthPrimaryDaughtersPDG", &truthPrimaryDaughtersPDG);
    tree->SetBranchAddress("truthPrimaryDaughtersProcess", &truthPrimaryDaughtersProcess);
    tree->SetBranchAddress("truthPrimaryDaughtersKE", &truthPrimaryDaughtersKE);

    std::vector<int>*         truthSecondaryPionDaughtersPDG     = nullptr;
    std::vector<std::string>* truthSecondaryPionDaughtersProcess = nullptr;
    std::vector<double>*      truthSecondaryPionDaughtersKE      = nullptr;
    tree->SetBranchAddress("truthSecondaryPionDaughtersPDG", &truthSecondaryPionDaughtersPDG);
    tree->SetBranchAddress("truthSecondaryPionDaughtersProcess", &truthSecondaryPionDaughtersProcess);
    tree->SetBranchAddress("truthSecondaryPionDaughtersKE", &truthSecondaryPionDaughtersKE);

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
    TH1D *hStitchAsThroughgoing  = new TH1D("hStitchAsThroughgoing", "hStitchAsThroughgoing;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
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

    TH1D *h0pInelasticBackgroundSecondaryInteraction = new TH1D("h0pInelasticBackgroundSecondaryInteraction", "h0pInelasticBackgroundSecondaryInteraction;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);
    TH1D *hNpInelasticBackgroundSecondaryInteraction = new TH1D("hNpInelasticBackgroundSecondaryInteraction", "hNpInelasticBackgroundSecondaryInteraction;;", NUM_BACKGROUND_TYPES, 0, NUM_BACKGROUND_TYPES);

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

    /////////////////////////////////
    // Files for event information //
    /////////////////////////////////

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

        // If no track matched to wire-chamber, skip
        if (WC2TPCtrkID == -99999) continue;

        // Apply small tracks and reduced volume cuts
        if (!passesSmallTracksCut)  continue;

        // Perform chi^2 stitching for primary track
        // For this analysis, we get the following chi^2 values:
        //   - Pion with Bragg peak chi^2
        //   - Through-going chi^2
        //   - Proton chi^2
        //   - Scanning fit with rhs to through-going and lhs to proton

        int totalCaloPoints = wcMatchDEDX->size();
        int nRemoveOutliers = 2;
        int nRemoveEnds     = 3;
        int minPoints       = 5;

        outWCAll << "Event number: " << event << std::endl; 
        outWCAll << "  Calorimetry data points: " << totalCaloPoints << std::endl;

        // Get chi^2 fits, primary tracks are already checked for reversal in first module
        double pionChi2         = computeReducedChi2(gPion, *wcMatchResR,  *wcMatchDEDX, false, totalCaloPoints, nRemoveOutliers, nRemoveEnds);
        double throughGoingChi2 = computeReducedChi2(gMuonTG, *wcMatchResR, *wcMatchDEDX, false, totalCaloPoints, nRemoveOutliers, nRemoveEnds);
        double protonChi2       = computeReducedChi2(gProton, *wcMatchResR, *wcMatchDEDX, false, totalCaloPoints, nRemoveOutliers, nRemoveEnds);

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
        outWCAll << "  Through-going pion chi^2: " << throughGoingChi2 << std::endl; 
        outWCAll << "  Proton with Bragg peak chi^2: " << protonChi2 << std::endl; 
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
        double minChi2 = std::min({minStitchedChi2, pionChi2, throughGoingChi2, protonChi2});

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

        if (minChi2 == minStitchedChi2) {
            double distanceFromVertex         = distance(breakPointX, truthPrimaryVertexX, breakPointY, truthPrimaryVertexY, breakPointZ, truthPrimaryVertexZ);
            double originalDistanceFromVertex = distance(WC2TPCPrimaryEndX, truthPrimaryVertexX, WC2TPCPrimaryEndY, truthPrimaryVertexY, WC2TPCPrimaryEndZ, truthPrimaryVertexZ);

            // Track looks the most like a pion + proton, background
            outStitchedFile << "Event number: " << event << std::endl; 
            outStitchedFile << "  Pion with Bragg peak chi^2: " << pionChi2 << std::endl; 
            outStitchedFile << "  Through-going pion chi^2: " << throughGoingChi2 << std::endl; 
            outStitchedFile << "  Proton with Bragg peak chi^2: " << protonChi2 << std::endl; 
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

            if ((backgroundType == 1) || (backgroundType == 7)) {
                hStitchedDistanceFromVertex->Fill(distanceFromVertex);
                hStitchedOriginalDistanceFromVertex->Fill(originalDistanceFromVertex);
            }
            // hStitchedDistanceFromVertex->Fill(distanceFromVertex);
            // hStitchedOriginalDistanceFromVertex->Fill(originalDistanceFromVertex);
            hStitchAsPionAndProton->Fill(backgroundType);
            hStitchFracBreakPoints->Fill(fracBreakPoint);
            
        } else if (minChi2 == protonChi2) {
            // Track looks the most like a proton, background
            hStitchAsProton->Fill(backgroundType);
        } else if (minChi2 == pionChi2) {
            // Track looks like pion with Bragg peak, background
            hStitchAsPionBraggPeak->Fill(backgroundType);
            continue;
        } else {
            // Track looks like through-going pion, could be signal
            hStitchAsThroughgoing->Fill(backgroundType);
        }

        ////////////////////////
        // Continue selection //
        ////////////////////////

        // If the particle is categorized as throughgoing:
        //   - Chi^2 selection on secondary tracks
        //   - If no secondary tracks, accept as 0p event
        // If the particle is categorized as a pion w/ bragg peak:
        //   - Chi^2 selection on secondary tracks
        //   - If no secondary tracks, kill the event 
        // If the particle is categorized as proton:
        //   - Test if kill/keep
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

        // If primary is through-going:
        //   - Usual selection
        // If primary is pion with Bragg peak:
        //   - If no secondaries, reject
        //   - If secondaries, usual selection
        // If primary is proton:
        //   - Reject
        // If primary is stitched:
        //   - Usual selection

        if (minChi2 == pionChi2) {
            if (numTracksNearVertex == 0) continue;
        } else if (minChi2 == protonChi2) {
            continue;
        }

        int totalTaggedPions   = secondaryTaggedPion + otherTaggedPion;
        int totalTaggedProtons = secondaryTaggedProton + otherTaggedProton;

        if (totalTaggedPions > 0) {
            // If any secondary tagged as pion, not signal
            continue;
        } else {
            hRecoAbsorption->Fill(backgroundType);
            if (totalTaggedProtons == 0) { 
                hRecoAbsorption0p->Fill(backgroundType);
            } else if (totalTaggedProtons > 0) {
                hRecoAbsorptionNp->Fill(backgroundType);
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

        if ((totalTaggedProtons == 0) && (backgroundType != 0)) {
            printEventInfo(thisEventInfo, outFile0pBackground);
            if (backgroundType == 6) h0pInelasticBackgroundSecondaryInteraction->Fill(secondaryInteractionTag);
        } else if ((totalTaggedProtons > 0) && (backgroundType != 1)) {
            printEventInfo(thisEventInfo, outFileNpBackground);
            if (backgroundType == 6) hNpInelasticBackgroundSecondaryInteraction->Fill(secondaryInteractionTag);
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

    //////////////////
    // Create plots //
    //////////////////

    std::vector<int> Colors = {
        kBlack, kBlue, kRed, kGreen
    };

    std::vector<std::vector<TH1*>> PlotGroups = {
        {hStitchedDistanceFromVertex, hStitchedOriginalDistanceFromVertex},
        {hStitchAsPionAndProton, hStitchAsPionBraggPeak, hStitchAsThroughgoing, hStitchAsProton},
        {hStitchFracBreakPoints},
        {hSecondaryPionChi2Pions, hSecondaryPionChi2Protons, hSecondaryPionChi2Others},
        {hSecondaryProtonChi2Pions, hSecondaryProtonChi2Protons, hSecondaryProtonChi2Others},
        {hSecondaryMeanDEDXPions, hSecondaryMeanDEDXProtons, hSecondaryMeanDEDXOthers}
    };

    std::vector<std::vector<TString>> PlotLabelGroups = {
        {"Detected vertex distance", "Original distance"},
        {"Stitched pion-proton", "Bragg pion", "Through-going", "Proton"},
        {"Stitched tracks"},
        {"Pions", "Protons", "Others"},
        {"Pions", "Protons", "Others"},
        {"Pions", "Protons", "Others"}
    };

    std::vector<TString> PlotTitles = {
        "StitchedDistanceFromVertex",
        "StitchedBackgroundTypes",
        "StitchedFracBreakPoints",
        "PionChi2SecondaryTracks",
        "ProtonChi2SecondaryTracks",
        "MeanDEDXSecondaryTracks"
    };

    std::vector<TString> XLabels = {
        "Distance from true vertex (cm)",
        "Background type",
        "Fractional break point",
        "Pion reduced #chi^{2}",
        "Proton reduced #chi^{2}",
        "Mean dE/dx (MeV/cm)"

    };

    std::vector<TString> YLabels = {
        "Number of tracks",
        "Number of tracks",
        "Number of tracks",
        "Number of tracks",
        "Number of tracks",
        "Number of tracks"
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
        hPionChi2FOM
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
        "PionChi2FOM"
    };

    printTwoDPlots(SaveDir, TwoDPlots, TwoDTitles);

    //////////////////////////
    // Create stacked plots //
    //////////////////////////

    THStack* hBackgroundTypesStack = new THStack("hBackgroundTypesStack", "BackgroundTypesStack");

    double alpha = 0.2;
    hStitchAsPionAndProton->SetTitle("Stitched pion-proton");
    hStitchAsPionBraggPeak->SetTitle("Bragg pion");
    hStitchAsThroughgoing->SetTitle("Through-going");
    hStitchAsProton->SetTitle("Proton");

    hStitchAsPionAndProton->SetFillColorAlpha(Colors[0], alpha);
    hStitchAsPionBraggPeak->SetFillColorAlpha(Colors[1], alpha);
    hStitchAsThroughgoing->SetFillColorAlpha(Colors[2], alpha);
    hStitchAsProton->SetFillColorAlpha(Colors[3], alpha);

    hStitchAsPionAndProton->SetLineColor(Colors[0]);
    hStitchAsPionBraggPeak->SetLineColor(Colors[1]);
    hStitchAsThroughgoing->SetLineColor(Colors[2]);
    hStitchAsProton->SetLineColor(Colors[3]);

    hBackgroundTypesStack->Add(hStitchAsThroughgoing);
    hBackgroundTypesStack->Add(hStitchAsPionAndProton);
    hBackgroundTypesStack->Add(hStitchAsPionBraggPeak);
    hBackgroundTypesStack->Add(hStitchAsProton);

    TCanvas* c = new TCanvas("Canvas", "Canvas", 205, 34, 1024, 768);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);
    hBackgroundTypesStack->Draw("HIST");
    // c->BuildLegend();

    hBackgroundTypesStack->GetXaxis()->SetTitleFont(FontStyle);
    hBackgroundTypesStack->GetXaxis()->SetLabelFont(FontStyle);
    hBackgroundTypesStack->GetXaxis()->SetNdivisions(8);
    hBackgroundTypesStack->GetXaxis()->SetLabelSize(TextSize);
    hBackgroundTypesStack->GetXaxis()->SetTitleSize(TextSize);
    hBackgroundTypesStack->GetXaxis()->SetTitle("Background type");
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

    c->Update();
    c->SaveAs(SaveDir + "StitchedBackgroundTypesStacked.png");
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
}