#pragma once

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include "TGraph.h"

//////////////////////
// Global constants //
//////////////////////

std::map<int, std::string> backgroundTypes = {
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
    {12, "Elastic scattering"},
    {13, "Scattering 0p"},
    {14, "Scattering Np"}
};

static const std::unordered_set<int> SKIP_INDICES = {
    1402,
    15345,
    46200,
    152596,
    255737,
    315186,
    370495,
};

static const std::unordered_set<int> SKIP_INDICES_DATA = {
    311635
};

std::vector<int> Colors = {
    TColor::GetColor(64,  64,  64),   // dark charcoal
    TColor::GetColor(0,   114, 178),  // deep blue
    TColor::GetColor(213, 94,  0),    // burnt orange
    TColor::GetColor(0,   158, 115),  // teal green
    TColor::GetColor(200, 185, 0),    // bright yellow
    TColor::GetColor(86,  180, 233),  // sky blue
    TColor::GetColor(204, 121, 167),  // mauve pink
    TColor::GetColor(0,   158, 115),  // teal green
    TColor::GetColor(86,  180, 233),  // sky blue
    TColor::GetColor(204, 121, 167),  // mauve pink
};

// Loading products
const int kMaxTrack          = 10000; // 10000
const int kMaxWCTracks       = 1;     // 1000
const int kMaxTrackHits      = 1000;  // 1000
const int kMaxTrajHits       = 1000;  // 1000
const int kMaxPrimaries      = 20000; // 20000
const int kMaxPrimaryPart    = 50;    // 50
const int kMaxTruePrimaryPts = 5000;  // 5000
const int kMaxHits           = 20000; // 20000
const int kMaxIDE            = 5000;  // 5000

// Loading products (data)
const int kMaxTrackData     = 10000; // 10000
const int kMaxWCTracksData  = 1000;  // 1000
const int kMaxTrackHitsData = 1000;  // 1000
const int kMaxTrajHitsData  = 1000;  // 1000
const int kMaxHitsData      = 20000; // 20000
const int kMaxTOFData       = 100;   // 100

// Background types
int NUM_BACKGROUND_TYPES = backgroundTypes.size();

// Signal types
int NUM_SIGNAL_TYPES = 2; // abs, scattering

// Run through this many events
int USE_NUM_EVENTS = 100000;

// Test threshold for TG tracks
int N_TG_TRACKS = 3;

// Number of events at each data quality cut
double TG0_MC_EVENTS = 118999;
double TG1_MC_EVENTS = 156432;
double TG2_MC_EVENTS = 176123;

double TG0_DATA_EVENTS = 4742;
double TG1_DATA_EVENTS = 9041;
double TG2_DATA_EVENTS = 12867;

// For scaling: events before RV cut (beam-line mass, WC2TPC match, TG tracks)
double TOTAL_DATA_EVENTS      = 346349;
double RESTRICTED_DATA_EVENTS = 100000;
double TOTAL_MC_EVENTS        = 493581;
double NUM_DATA_EVENTS        = TG1_DATA_EVENTS * (TOTAL_DATA_EVENTS / RESTRICTED_DATA_EVENTS);
double NUM_MC_EVENTS          = TG1_MC_EVENTS * (TOTAL_DATA_EVENTS / TOTAL_MC_EVENTS);
double MC_SCALING             = NUM_DATA_EVENTS / NUM_MC_EVENTS;

// File with weights
extern TFile* fWeights;
const TString WEIGHTS_NAME = "hWeightsFrontFacePre";

// Sampling rate
double SAMPLING_RATE = 0.128;

// Drift velocity
double DRIFT_VELOCITY = 0.14732; // cm / \mu s

// Track pitch
const double TRACK_PITCH = 0.47;

// Wire-chamber quality
double MAX_RAD_DIST_WC4   = 8.0; // original 8.0
double MAX_MID_PLANE_DIST = 3.0; // original 3.0

struct EventInfo {
    int run; int subrun; int event; bool isData;
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
    int                plane;
    double             clusterSize;
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
const double RminZ = 30.0;
const double RmaxZ = 82.0;

// Detector dimensions
const double minX =  0.0;
const double maxX = 47.0;
const double minY =-20.0;
const double maxY = 20.0;
const double minZ =  3.0;
const double maxZ = 87.0;

// Uncertainty on beamline muon and electron composition
const double MUON_COMP_UNC     = 0.26;
const double ELECTRON_COMP_UNC = 1.00;

// Classify pions 
double PION_CHI2_PION_VALUE     = 0.750; // 1.125
double PROTON_CHI2_PION_VALUE   = 1.375; // 1.375

// Classify protons
double PION_CHI2_PROTON_VALUE   = 1.375; // 1.125
double PROTON_CHI2_PROTON_VALUE = 0.750; // 0.125

// Mean dE/dx threshold
double MEAN_DEDX_THRESHOLD = 4.0;

// Vertex radius
double VERTEX_RADIUS = 10.0;

// Number of points to use in mean dE/dx calculation
int MEAN_DEDX_NUM_TRAJ_POINTS = 5;

// Proton energy bounds
double PROTON_ENERGY_LOWER_BOUND = 0.075;
double PROTON_ENERGY_UPPER_BOUND = 1.0;

// Hit wire separation
double HIT_WIRE_SEPARATION = 0.4; // converts wire # to cm

// Window size for local linearity
int WINDOW_SIZE_LINEARITY = 8;

// Shower probability cut
double SHOWER_PROB_CUT = 0.5;

// Track length cuts charge exchange
double SMALL_TRACK_LENGTH_CHEX = 28.;
double LARGE_TRACK_LENGTH_CHEX = 30.;
double FAR_TRACK_DISTANCE_CHEX = 10.;
int    SMALL_TRACK_CUT_CHEX    = 3;

// Hit cluster cut
double DISTANCE_TO_PRIMARY_THRESHOLD = 10.;

double MAX_IN_CLUSTER_W_SEPARATION = 1.5;
double MAX_IN_CLUSTER_X_SEPARATION = 1.5;
double MAX_IN_CLUSTER_SEPARATION   = 2.5;

int MINIMUM_HITS_FOR_CLUSTER = 5;

int    NUM_CLUSTERS_THRESHOLD  = 2;
double LARGE_CLUSTER_THRESHOLD = 2.0;

// Cluster cut thresholds
int MAX_NUM_CLUSTERS_INDUCTION        = 1;
int MAX_NUM_CLUSTERS_COLLECTION       = 1;
int MAX_NUM_LARGE_CLUSTERS_INDUCTION  = 1;
int MAX_NUM_LARGE_CLUSTERS_COLLECTION = 1;
double MAX_LARGEST_CLUSTER_INDUCTION  = 2.;
double MAX_LARGEST_CLUSTER_COLLECTION = 2.;

// Threshold(s) used to cut based on primary track local linearity
double LINEARITY_DERIVATIVE_THRESHOLD = 0.06e-3;
double STD_DEV_LINEARITY_THRESHOLD    = 0.045E-3;

// Masses
const double PionMass    = 139.57018;    // in MeV
const double ProtonMass  = 938.27208816; // in MeV
const double NeutronMass = 939.5654133;  // in MeV

// Threshold for individual dE/dx hits
double HIT_DEDX_THRESHOLD = 40.;

// Threshold for scattering angle
double SCATTERING_ANGLE_THRESHOLD = 5.0 * (TMath::Pi() / 180);

// Threshold for pion scattering energy 
double PION_SCATTERING_ENERGY_THRESHOLD = 0.05; // in GeV / c

// Bins for energy bins
static const std::vector<double> ARRAY_KE_BINS = {0., 100., 200., 300., 500.};
static const int                   NUM_BINS_KE = static_cast<int>(ARRAY_KE_BINS.size()) - 1;

// Bins for energy bins (finer)
static const std::vector<double> ARRAY_KE_FINE_BINS = {
      0., 10.,   20.,  30.,  40., 50.,
     60., 70.,   80.,  90., 100.,
    110., 120., 130., 140., 150.,
    160., 170., 180., 190., 200.,
    210., 220., 230., 240., 250.,
    260., 270., 280., 290., 300.,
    310., 320., 330., 340., 350.,
    360., 370., 380., 390., 400.,
    410., 420., 430., 440., 450.,
    460., 470., 480., 490., 500.,
    510., 520., 530., 540., 550.,
    560., 570., 580., 590., 600.
};
static const int NUM_BINS_KE_FINE = static_cast<int>(ARRAY_KE_FINE_BINS.size()) - 1;

// TOF mass cut
double PI_MU_EL_MASS_CUTOFF = 350.;

// Cylinder radius for shower identification
double CYLINDER_RADIUS               = 10.;
double CYLINDER_SMALL_TRACK          = 5.;
int    ALLOWED_CYLINDER_SMALL_TRACKS = 1;
int    ALLOWED_CYLINDER_TRACKS       = 1;

// Sliced cone for ch exch events parameters
double SLICED_CONE_MIN_RADIUS           = 10.0;
double SLICED_CONE_MAX_RADIUS           = 30.0;
double SLICED_CONE_HEIGHT               = 30.0;
double SLICED_CONE_SMALL_TRACK          = 5.0;
int    SLICED_CONE_ALLOWED_SMALL_TRACKS = 0;

// Number of allowed TG tracks
int MAX_NUM_TG_TRACKS = 2;

// Parameters for BDT model
int BDT_NUM_RECO_TRKS = 5;

// Cross section calculation
double rho            = 1396; // kg / m^3
double molar_mass     = 39.95; // g / mol
double g_per_kg       = 1000; 
double avogadro       = 6.022e+23; // number / mol
double number_density = rho * g_per_kg / molar_mass * avogadro;
double slab_width     = 0.0047; // in m

double XSEC_UNITS = 1 / (number_density * slab_width * 1e-28);

// For each "collimator", the bounds of the face of both aperatures [xlow_frontface, xhigh_frontface, xlow_backface, xhigh_backface], 
// similarly for y. In cm, in TPC coordinates. Taken from survey.
double xboundMagnet1[4] = {45.74, 75.52, 35.09, 64.87};
double yboundMagnet1[4] = {-13.12, 13.59, -13.16, 13.55};

double xboundMagnet2[4] = {34.88, 65.12, 27.90, 58.15};
double yboundMagnet2[4] = {-13.10, 13.61, -13.08, 13.63};

double xboundDSCol[4] = {30.33, 45.42, 26.65, 41.73};
double yboundDSCol[4] = {-15.70, 14.91, -15.53, 15.08};

// Z Position of the center of the aperatures of each collimator, found by taking the average of the z bounds of the aperature. [zcent_US, zcent_DS]
double zcentMagnet1[2] = { (-501.95-494.98)/2, (-449.49-456.46)/2};
double zcentMagnet2[2] = { (-432.04-427.50)/2, (-381.27-385.81)/2};
double zcentDSCol[2]   = { (-296.67-297.36)/2, (-205.94-206.63)/2};

/////////////////////
// Event variables //
/////////////////////

struct EventVariables {
    // Weight
    double weight = 1.0;

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

    // Truth primary daughters
    std::vector<int>         truthPrimaryDaughtersPDG;
    std::vector<std::string> truthPrimaryDaughtersProcess;
    std::vector<double>      truthPrimaryDaughtersKE;

    // Truth primary trajectory
    std::vector<double> truthPrimaryLocationX;
    std::vector<double> truthPrimaryLocationY;
    std::vector<double> truthPrimaryLocationZ;

    // Truth primary scalars
    int    truthPrimaryPDG        = -999;
    double truthPrimaryVertexX    = -999;
    double truthPrimaryVertexY    = -999;
    double truthPrimaryVertexZ    = -999;
    double truthPrimaryVertexKE   = -999;

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

    // Hit information
    std::vector<int>   fHitPlane;
    std::vector<int>   hitRecoAsTrackKey;
    std::vector<int>   hitWC2TPCKey;
    std::vector<float> fHitT;
    std::vector<float> fHitX;
    std::vector<float> fHitW;

    // Trajectory interaction
    bool        interactionInTrajectory    = false;
    std::string trajectoryInteractionLabel = "";
    double      trajectoryInteractionAngle = -999;
    double      trajectoryInteractionKE    = -999;
    double      trajectoryInitialMomentumX = -999;

    // Reco track hits
    std::vector<std::vector<int>>    recoTrackHitIndices;
    std::vector<std::vector<double>> recoTrackHitX;
    std::vector<std::vector<double>> recoTrackHitY;
    std::vector<std::vector<double>> recoTrackHitZ;

    // Scattering truth
    double truthScatteringAngle = -999;
    double truthScatteredPionKE = -999;

    // Secondary interactions
    std::vector<int>                 secondaryInteractionTypes;
    std::vector<double>              secondaryInteractionInteractingKE;
    std::vector<double>              secondaryInteractionAngle;
    std::vector<double>              secondaryInteractionXPosition;
    std::vector<double>              secondaryInteractionYPosition;
    std::vector<double>              secondaryInteractionZPosition;
    std::vector<std::vector<int>>    secondaryInteractionDaughtersPDG;
    std::vector<std::vector<double>> secondaryInteractionDaughtersKE;
    std::vector<std::vector<double>> secondaryIncidentKEContributions;

    // Incident KE
    bool                validTrueIncidentKE      = true;
    std::vector<double> trueIncidentKEContributions;

    // Classification
    int backgroundType    = -1;
    int numVisibleProtons = 0;
};

struct EventVariablesData {
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

//////////////////////
// Helper functions //
//////////////////////

void initializeProtonPoints(TGraph* gProton);
void initializePionPoints(TGraph* gPion);
void initializeMuonNoBraggPoints(TGraph* gMuonTG);
void printEventInfo(EventInfo event, std::ostream& os);
void printOneDPlots(const TString& dir, int fontStyle, double textSize, std::vector<std::vector<TH1*>>& groups, std::vector<int>& colors, std::vector<std::vector<TString>>& labels, std::vector<TString>& titles, std::vector<TString>& xlabels, std::vector<TString>& ylabels, std::vector<bool>& stack, std::vector<std::vector<bool>> asPoints = {}, std::vector<std::vector<bool>> addText = {});
void printTwoDPlots(const TString& dir, const std::vector<TH2*>& plots, const std::vector<TString>& titles, const std::vector<std::pair<double,double>>& zRanges = {}, const std::vector<bool>& displayNumbers = {});
void Overlay_dEdx_RR_Reference_PP(TGraph* gProton, TGraph* gPion, TGraph* gMIP, bool truncate = false, double MIPStart = 0., bool addLegend = true, TVirtualPad* pad = gPad);
void printEfficiencyPlots(TString dir, int fontStyle, double textSize, std::vector<TEfficiency*> efficiencies, std::vector<TString> titles, std::vector<TString> xlabels);
void printBackgroundInfo(TH1D* background_histo, std::ostream& os);
void reweightOneDHisto(TH1D* histo, double weight);

bool isWithinReducedVolume(double x, double y, double z);
bool isWithinActiveVolume(double x, double y, double z);
bool isHitNearPrimary(
    std::vector<int>* primaryKey, 
    std::vector<float>* hitX, 
    std::vector<float>* hitW, 
    std::vector<int>* hitPlane,
    float thisHitX, 
    float thisHitW, 
    int thisHitPlane,
    float threshold,
    bool onlyVertex = false
);

double computeReducedChi2(const TGraph* theory, std::vector<double> xData, std::vector<double> yData,  bool dataReversed, int nPoints, int nOutliersToDiscard = 0, int nTrim = 0);
double distance(double x1, double x2, double y1, double y2, double z1, double z2);
double meanDEDX(std::vector<double> trackDEDX, bool isTrackReversed, int pointsToUse);
double energyLossCalculation();
double energyLossCalculation(double x, double px, bool isData);
double getClusterElongation(HitCluster cluster);
double getClusterWidth(HitCluster cluster);

int isSecondaryInteractionAbsorption(std::vector<int> daughtersPDG, std::vector<string> daughtersProcess, std::vector<double> daughtersKE);
int flattenIndex (int beta, int j, int S);
int getBin(double data, const std::vector<double>& edges);

std::pair<int,int> unflattenIndex(int f, int S);

std::vector<double> calcLinearityProfile(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, int nb);
std::pair<TH1*, TH1*> getBinEfficiencyAndPurity(TH1* hTrue, TH1* hReco, TH1* hRecoTrue);

void H2M(const TH2D* histo, TMatrixD& mat, bool rowcolumn);
void H2V(const TH1D* histo, TVectorD& vec);
void M2H(const TMatrixD mat, TH2D* histo);
void V2H(const TVectorD vec, TH1D* histo);

void removeRepeatedPoints(std::vector<double>* x, std::vector<double>* y, std::vector<double>* z);
void getBDTVariables(
    int WC2TPCtrkID,
    const std::vector<double>* wcLocationX,
    const std::vector<double>* wcLocationY,
    const std::vector<double>* wcLocationZ,
    const std::vector<double>* wcX,
    const std::vector<double>* wcY,
    const std::vector<double>* wcZ,
    const std::vector<double>* recoBeginX,
    const std::vector<double>* recoBeginY,
    const std::vector<double>* recoBeginZ,
    const std::vector<double>* recoEndX,
    const std::vector<double>* recoEndY,
    const std::vector<double>* recoEndZ,
    const std::vector<std::vector<double>>* recoDEDX,
    const std::vector<int>*    recoTrkID,
    int maxRecoTrks,

    Float_t& bdt_WC2TPCTrackLength,
    Float_t& bdt_WC2TPCBeginX,
    Float_t& bdt_WC2TPCBeginY,
    Float_t& bdt_WC2TPCBeginZ,
    Float_t& bdt_WC2TPCEndX,
    Float_t& bdt_WC2TPCEndY,
    Float_t& bdt_WC2TPCEndZ,
    Float_t& bdt_WC2TPCdEdx,

    Float_t* bdt_recoTrkBeginX,
    Float_t* bdt_recoTrkBeginY,
    Float_t* bdt_recoTrkBeginZ,
    Float_t* bdt_recoTrkEndX,
    Float_t* bdt_recoTrkEndY,
    Float_t* bdt_recoTrkEndZ,
    Float_t* bdt_recoTrkLen,
    Float_t* bdt_recoTrkdEdx,
    Float_t& bdt_numRecoTrksInCylinder
);

void histoToText(const TH1D* hist, const std::string& filename);

bool IsPointInsideTrackCylinder(
    const std::vector<double>* primaryX,
    const std::vector<double>* primaryY,
    const std::vector<double>* primaryZ,
    double x, double y, double z,
    double radius,
    double* minDist2_out = nullptr,
    std::size_t* segIndex_out = nullptr
);

bool IsPointInsideSlicedCone(
    double px, double py, double pz,
    double ax, double ay, double az,
    double vx, double vy, double vz,
    double height,
    double r0,
    double r1
);

std::tuple<double, double> azimuth_polar_from_points(
    double xs, double ys, double zs,
    double xe, double ye, double ze
);

std::vector<double> projToZ(const std::vector<double>& hit0, const std::vector<double>& hit1, double zpos);

std::vector<double> getAverageDir(const std::vector<double>& points);

std::pair<size_t, size_t> find_unique_position(const std::vector<std::vector<int>>* v, int n);

bool EqualApprox(const TMatrixD& A, const TMatrixD& B, double rtol = 1e-12, double atol = 1e-15);

///////////////
// WC checks //
///////////////

bool CheckUpstreamMagnetAperture(const std::vector<double>& hit1, const std::vector<double>& hit2);
bool CheckDownstreamMagnetAperture(const std::vector<double>& hit1, const std::vector<double>& hit2);
bool CheckDownstreamCollimatorAperture(const std::vector<double>& hit1, const std::vector<double>& hit2);

//////////////////////////
// Wiener SVD unfolding //
// Author: Hanyu WEI    //
//////////////////////////

TMatrixD Matrix_C(Int_t n, Int_t type);

// Wiener filter unfolding:
//   - first five parameters are inputs as names read
//   - AddSmear is additional smearing matrix after unfolding
//   - WF is the elements of Wiener Filter
//   - Covariance matrix of the unfolded spectrum
TVectorD WienerSVD(TMatrixD Response, TVectorD Signal, TVectorD Measure, TMatrixD Covariance, Int_t C_type, Float_t Norm_type, TMatrixD& AddSmear, TVectorD& WF, TMatrixD& UnfoldCov, TMatrixD& CovRotation, TMatrixD& AddSmearInverse);

/////////////////
// Systematics //
/////////////////

std::vector<TH1D*> MakeUniverseHists(const std::string& baseName, const std::string& title, int nBins, const double* binEdges, int nUniverses);
std::vector<std::vector<std::vector<TH1D*>>> MakeUniverseHistBlock(const std::string& baseName, const std::string& titleFmt, int nUniverses, int nOuter, int nInner, int nbins, const double* binEdges);

void GetResponseMatrix(int SizeOuter, int SizeInner, const std::vector<TH1D*>& TotalEventsHistos, const std::vector<std::vector<std::vector<TH1D*>>>& TrueRecoAsByBin, const TH1D* RecoIncidentFlux, const TH1D* TrueIncidentFlux, TH2D* ReponseMatrix);
void GetCovMatrix(TH1* RecoHisto, std::vector<TH1D*> UnivRecoHisto, TH2* CovMatrix);
void GetFracCovAndCorrMatrix(TH1* RecoHisto, TH2* CovMatrix, TH2* FracCovMatrix, TH2* CorrMatrix);
double FindQuantile(double frac, std::vector<double>& xs_in);
void DrawHistosWithErrorBands(TH1D* RecoHisto, std::vector<TH1D*> UnivRecoHisto, const TString& dir, const TString& SystName, const TString& PlotName, double textSize, int fontStyle, const TString& title, const TString& xlabel, const TString& ylabel);
void DrawErrorBand(TH1* nom, TGraphAsymmErrors* band, int bandCol, double alpha);
void PrintTruthUnfoldedStatShapePlot(const TString& SaveDir, const TString& Name, TH1* hTrue, TH1* hUnfolded, TH1* hStat, TH1* hShape, const TString& PlotTitle, const TString& XTitle, const TString& YTitle, int FontStyle = 132, double TextSize = 0.06);
void PrintDataUnfoldedStatShapePlot(const TString& SaveDir, const TString& Name, TH1* hTrue, TH1* hUnfolded, TH1* hStat, TH1* hShape, const TString& PlotTitle, const TString& XTitle, const TString& YTitle, int FontStyle = 132, double TextSize = 0.06);
void PrintFDPlot(const TString& SaveDir, const TString& Name, TH1* hTrue, TH1* hNomTrue, TH1* hUnfolded, TH1* hSysts, const TString& PlotTitle, const TString& XTitle, const TString& YTitle, std::pair<double, double> chi, std::pair<int, int> ndof, std::pair<double, double> pval, std::pair<double, double> sigma, int FontStyle = 132, double TextSize = 0.06);
void CalcChiSquared(TH1D* h_model, TH1D* h_data, TH2D* cov, double &chi, int &ndof, double &pval, double &sigma);

void PrintDataVsMCContribPlot(const TString& SaveDir, const TString& Name, TH1* hData, const std::vector<TH1*>& mcHists, const std::vector<TString>& mcLabels, const std::vector<int>& mcColors, const TString& PlotTitle, const TString& XTitle, const TString& YTitle, int FontStyle, double TextSize, bool UsePoissonDataErrors = true, TH1* hMCAbsUnc = nullptr, bool DrawRatio = true);
void PrintDataVsTwoMCPlot(const TString& SaveDir, const TString& Name, TH1* hData, TH1* hMC1, TH1* hMC2, const TString& label1, const TString& label2, int color1, int color2, const TString& PlotTitle, const TString& XTitle, const TString& YTitle, int FontStyle, double TextSize);

int fillSignalInformation(int pdg, double vx, double vy, double vz, bool interactionInTrajectory, std::string trajectoryInteractionLabel, std::vector<int> daughtersPDG, std::vector<std::string> daughtersProcess,  std::vector<double> daughtersKE, bool saveNumProtons, int& numVisibleProtons);
int getBackgroundInteractionType(int pdg, double vx, double vy, double vz, bool interactionInTrajectory, std::string trajectoryInteractionLabel, std::vector<int> daughtersPDG, std::vector<std::string> daughtersProcess, std::vector<double> daughtersKE);
std::string ProcessToString(int proc);

void printEventVariables(const EventVariables& t, std::ofstream& out);

// Reweighing
void ApplyWeights(TH1* h, const TH1* weights);
void ComputeReweightingWeights(const TH1* source, const TH1* target, TH1* weights);
double GetKEWeight(TH1D* hWeights, double ke);