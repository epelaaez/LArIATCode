#pragma once

#include <vector>
#include <string>
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

// Background types
int NUM_BACKGROUND_TYPES = backgroundTypes.size();

// Signal types
int NUM_SIGNAL_TYPES = 5; // for now: 0p abs, Np abs, 0p scatter, Np scatter, charge exchange

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

// Values for chi^2 secondary fits
double PION_CHI2_PION_VALUE     = 1.125;
double PION_CHI2_PROTON_VALUE   = 1.125;
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

// Shower probability cut
double SHOWER_PROB_CUT = 0.5;

// Track length cuts charge exchange
double SMALL_TRACK_LENGTH_CHEX = 28.;
double LARGE_TRACK_LENGTH_CHEX = 30.;
double FAR_TRACK_DISTANCE_CHEX = 10.;
int    SMALL_TRACK_CUT_CHEX    = 3;

// Hit cluster cut
int    NUM_CLUSTERS_THRESHOLD  = 2;
double LARGE_CLUSTER_THRESHOLD = HIT_WIRE_SEPARATION * 3;

// Threshold(s) used to cut based on primary track local linearity
double LINEARITY_DERIVATIVE_THRESHOLD = 0.06e-3;
double STD_DEV_LINEARITY_THRESHOLD    = 0.045E-3;

// Masses
const double PionMass    = 139.57018;    // in MeV
const double ProtonMass  = 938.27208816; // in MeV
const double NeutronMass = 939.5654133;  // in MeV

// Threshold for individual dE/dx hits
double HIT_DEDX_THRESHOLD = 40.;

// Bins for energy bins
double LOWER_BOUND_KE = 0.;
double UPPER_BOUND_KE = 500.;
int    NUM_BINS_KE    = 10;

// TOF mass cut
double PI_MU_EL_MASS_CUTOFF = 350.;

// Cross section calculation
double rho            = 1396; // kg / m^3
double molar_mass     = 39.95; // g / mol
double g_per_kg       = 1000; 
double avogadro       = 6.022e+23; // number / mol
double number_density = rho * g_per_kg / molar_mass * avogadro;
double slab_width     = 0.0047; // in m

//////////////////////
// Helper functions //
//////////////////////

void initializeProtonPoints(TGraph* gProton);
void initializePionPoints(TGraph* gPion);
void initializeMuonNoBraggPoints(TGraph* gMuonTG);
void printEventInfo(EventInfo event, std::ostream& os);
void printOneDPlots(const TString& dir, int fontStyle, double textSize, std::vector<std::vector<TH1*>>& groups, std::vector<int>& colors, std::vector<std::vector<TString>>& labels, std::vector<TString>& titles, std::vector<TString>& xlabels, std::vector<TString>& ylabels, std::vector<bool>& stack, std::vector<std::vector<bool>> asPoints = {});
void printTwoDPlots(const TString& dir, const std::vector<TH2*>& plots, const std::vector<TString>& titles, const std::vector<std::pair<double,double>>& zRanges = {}, const std::vector<bool>& displayNumbers = {});
void Overlay_dEdx_RR_Reference_PP(TGraph* gProton, TGraph* gPion, TGraph* gMIP, bool truncate = false, double MIPStart = 0., bool addLegend = true, TVirtualPad* pad = gPad);
void printEfficiencyPlots(TString dir, int fontStyle, double textSize, std::vector<TEfficiency*> efficiencies, std::vector<TString> titles, std::vector<TString> xlabels);
void printBackgroundInfo(TH1D* background_histo, std::ostream& os);

bool isWithinReducedVolume(double x, double y, double z);
bool isHitNearPrimary(std::vector<int>* primaryKey, std::vector<float>* hitX, std::vector<float>* hitW, float thisHitX, float thisHitW, float xThreshold, float wThreshold);

double computeReducedChi2(const TGraph* theory, std::vector<double> xData, std::vector<double> yData,  bool dataReversed, int nPoints, int nOutliersToDiscard = 0, int nTrim = 0);
double distance(double x1, double x2, double y1, double y2, double z1, double z2);
double meanDEDX(std::vector<double> trackDEDX, bool isTrackReversed, int pointsToUse);
double energyLossCalculation();
double energyLossCalculation(double x, double px, bool isData);
double getClusterElongation(HitCluster cluster);
double getClusterWidth(HitCluster cluster);

int isSecondaryInteractionAbsorption(std::vector<int> daughtersPDG, std::vector<string> daughtersProcess, std::vector<double> daughtersKE);
int getCorrespondingBin(double value, int num_bins, double low, double high);
int flattenIndex (int beta, int j, int S);
std::pair<int,int> unflattenIndex(int f, int S);

std::vector<double> calcLinearityProfile(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, int nb);
std::pair<TH1*, TH1*> getBinEfficiencyAndPurity(TH1* hTrue, TH1* hReco, TH1* hRecoTrue);

void H2M(const TH2D* histo, TMatrixD& mat, bool rowcolumn);
void H2V(const TH1D* histo, TVectorD& vec);
void M2H(const TMatrixD mat, TH2D* histo);
void V2H(const TVectorD vec, TH1D* histo);

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
TVectorD WienerSVD(TMatrixD Response, TVectorD Signal, TVectorD Measure, TMatrixD Covariance, Int_t C_type, Float_t Norm_type, TMatrixD& AddSmear, TVectorD& WF, TMatrixD& UnfoldCov, TMatrixD& CovRotation);