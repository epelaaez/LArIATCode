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

void BeamlineElCov() {
    // This module generates all relevant covariance matrices coming from beamline content variations

    // Set defaults
    gStyle->SetOptStat(0); // get rid of stats box
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    TH1::AddDirectory(false);
    gStyle->SetPalette(kRainBow);

    TGraph* gProton = new TGraph();
    TGraph* gPion   = new TGraph();
    TGraph* gMuonTG = new TGraph();

    // Number of universes for this systematic
    int NUM_UNIVERSES = 2; // +/- mu unc

    // Initialize points
    initializeProtonPoints(gProton);
    initializePionPoints(gPion);
    initializeMuonNoBraggPoints(gMuonTG);

    int FontStyle = 132;
    double TextSize = 0.06;
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis/figs/Systematics/";

    // Load file with data products
    TString RootFilePath = "/exp/lariat/app/users/epelaez/files/RecoAll_histo.root"; // RV at z = 30
    std::unique_ptr<TFile> File(TFile::Open(RootFilePath));
    TDirectory* Directory = (TDirectory*)File->Get("RecoNNAllEval");

    ///////////////////
    // Load branches //
    ///////////////////

    // TODO: delete unnecessary data products to make this faster

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
    std::vector<std::vector<int>>*    secondaryInteractionDaughtersPDG = nullptr;
    std::vector<std::vector<double>>* secondaryInteractionDaughtersKE = nullptr;
    std::vector<std::vector<double>>* secondaryIncidentKEContributions = nullptr;
    tree->SetBranchAddress("secondaryInteractionTypes", &secondaryInteractionTypes);
    tree->SetBranchAddress("secondaryInteractionTrkID", &secondaryInteractionTrkID);
    tree->SetBranchAddress("secondaryInteractionInteractingKE", &secondaryInteractionInteractingKE);
    tree->SetBranchAddress("secondaryInteractionAngle", &secondaryInteractionAngle);
    tree->SetBranchAddress("secondaryInteractionXPosition", &secondaryInteractionXPosition);
    tree->SetBranchAddress("secondaryInteractionYPosition", &secondaryInteractionYPosition);
    tree->SetBranchAddress("secondaryInteractionZPosition", &secondaryInteractionZPosition);
    tree->SetBranchAddress("secondaryInteractionDaughtersPDG", &secondaryInteractionDaughtersPDG);
    tree->SetBranchAddress("secondaryInteractionDaughtersKE", &secondaryInteractionDaughtersKE);
    tree->SetBranchAddress("secondaryIncidentKEContributions", &secondaryIncidentKEContributions);

    ///////////////////////
    // Create histograms //
    ///////////////////////

    /////////////////
    // Incident KE //
    /////////////////

    // Selected incident pion
    std::vector<TH1D*> PionIncidentKEUnivs = MakeUniverseHists("hIncidentKE", "IncidentKE;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hPionIncidentKENom = new TH1D("hPionIncidentKENom", "hPionIncidentKENom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Selected incident pion that is true pion
    std::vector<TH1D*> PionIncidentKETrueUnivs = MakeUniverseHists("hIncidentKETrue", "IncidentKETrue;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hPionIncidentKETrueNom = new TH1D("hPionIncidentKETrueNom", "hPionIncidentKETrueNom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // True incident KE
    std::vector<TH1D*> TrueIncidentKEUnivs = MakeUniverseHists("hTrueIncidentKE", "TrueIncidentKE;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hTrueIncidentKENom = new TH1D("hTrueIncidentKENom", "hTrueIncidentKENom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Psi corrections
    std::vector<TH1D*> PsiIncidentUnivs = MakeUniverseHists("PsiIncidentUnivs", "PsiIncidentUnivs;Kinetic Energy [MeV];Factor", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hPsiIncidentNom = new TH1D("hPsiIncidentNom", "hPsiIncidentNom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // C corrections
    std::vector<TH1D*> CIncidentUnivs = MakeUniverseHists("CIncidentUnivs", "CIncidentUnivs;Kinetic Energy [MeV];Factor", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hCIncidentNom = new TH1D("hCIncidentNom", "hCIncidentNom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Corrected selected incident pion
    std::vector<TH1D*> PionIncidentCorrectedKEUnivs = MakeUniverseHists("hPionIncidentCorrectedKE", "PionIncidentCorrectedKE;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hPionIncidentCorrectedKENom = new TH1D("hPionIncidentCorrectedKENom", "hPionIncidentCorrectedKENom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    /////////////////
    // Selected KE // 
    /////////////////

    // Selected pion abs 0p for universes
    std::vector<TH1D*> PionAbs0pKEUnivs = MakeUniverseHists("hPionAbs0pKE", "PionAbs0pKE;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hPionAbs0pKENom = new TH1D("hPionAbs0pKENom", "hPionAbs0pKENom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Pion abs 0p estimated backgrounds
    std::vector<TH1D*> PionAbs0pKEBkgUnivs = MakeUniverseHists("hPionAbs0pKEBkg", "hPionAbs0pKEBkg;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hPionAbs0pKEBkgNom = new TH1D("hPionAbs0pKEBkgNom", "hPionAbs0pKEBkgNom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Selected pion abs Np for universes
    std::vector<TH1D*> PionAbsNpKEUnivs = MakeUniverseHists("hPionAbsNpKE", "PionAbsNpKE;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hPionAbsNpKENom = new TH1D("hPionAbsNpKENom", "hPionAbsNpKENom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Pion abs Np estimated backgrounds
    std::vector<TH1D*> PionAbsNpKEBkgUnivs = MakeUniverseHists("hPionAbsNpKEBkg", "hPionAbsNpKEBkg;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hPionAbsNpKEBkgNom = new TH1D("hPionAbsNpKEBkgNom", "hPionAbsNpKEBkgNom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Selected pion scattering for universes
    std::vector<TH1D*> PionScatterKEUnivs = MakeUniverseHists("hPionScatterKE", "PionScatterKE;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hPionScatterKENom = new TH1D("hPionScatterKENom", "hPionScatterKENom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // Pion scattering estimated backgrounds
    std::vector<TH1D*> PionScatterKEBkgUnivs = MakeUniverseHists("hPionScatterKEBkg", "hPionScatterKEBkg;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hPionScatterKEBkgNom = new TH1D("hPionScatterKEBkgNom", "hPionScatterKEBkgNom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    std::vector<TH1D*> RecoSignals = {
        hPionAbs0pKENom, hPionAbsNpKENom, hPionScatterKENom
    };

    std::vector<TH1D*> RecoBkgSignals = {
        hPionAbs0pKEBkgNom, hPionAbsNpKEBkgNom, hPionScatterKEBkgNom
    };

    /////////////
    // True KE //
    /////////////

    // True abs 0p KE
    std::vector<TH1D*> TruePionAbs0pKEUnivs = MakeUniverseHists("hTruePionAbs0pKE", "TruePionAbs0pKE;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hTruePionAbs0pKENom = new TH1D("hTruePionAbs0pKENom", "hTruePionAbs0pKENom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // True abs Np KE
    std::vector<TH1D*> TruePionAbsNpKEUnivs = MakeUniverseHists("hTruePionAbsNpKE", "TruePionAbsNpKE;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hTruePionAbsNpKENom = new TH1D("hTruePionAbsNpKENom", "hTruePionAbsNpKENom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    // True abs scattering KE
    std::vector<TH1D*> TruePionScatterKEUnivs = MakeUniverseHists("hTruePionScatterKE", "TruePionScatterKE;Kinetic Energy [MeV];Counts", NUM_BINS_KE, ARRAY_KE_BINS.data(), NUM_UNIVERSES);
    TH1D* hTruePionScatterKENom = new TH1D("hTruePionScatterKENom", "hTruePionScatterKENom;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    std::vector<TH1D*> TotalEventsHistosNom = {
        hTruePionAbs0pKENom, hTruePionAbsNpKENom, hTruePionScatterKENom
    };

    /////////////////////////
    // For response matrix //
    /////////////////////////

    // For all universes
    std::vector<std::vector<std::vector<TH1D*>>> TrueAbs0pAsByBinUnivs = MakeUniverseHistBlock(
        "hTrueAbs0p", "True Abs0p (U=%d): As %d in true bin %d",
        NUM_UNIVERSES, NUM_BINS_KE, NUM_SIGNAL_TYPES,
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    std::vector<std::vector<std::vector<TH1D*>>> TrueAbsNpAsByBinUnivs = MakeUniverseHistBlock(
        "hTrueAbsNp", "True AbsNp (U=%d): As %d in true bin %d",
        NUM_UNIVERSES, NUM_BINS_KE, NUM_SIGNAL_TYPES,
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    std::vector<std::vector<std::vector<TH1D*>>> TrueScatterAsByBinUnivs = MakeUniverseHistBlock(
        "hTrueScatter", "True Scatter (U=%d): As %d in true bin %d",
        NUM_UNIVERSES, NUM_BINS_KE, NUM_SIGNAL_TYPES,
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );

    // Nominal
    std::vector<std::vector<TH1D*>> TrueAbs0pAsByBinNom;
    for (int iEnergyBin = 0; iEnergyBin < NUM_BINS_KE; ++iEnergyBin) {
        std::vector<TH1D*> TempVec;
        for (int iInt = 0; iInt < NUM_SIGNAL_TYPES; ++iInt) {
            TH1D* hTempHist = new TH1D(Form("hTrueAbs0p_%d_Bin_As_%d", iEnergyBin, iInt), Form("True Abs 0p KE As %d in bin %d", iInt, iEnergyBin), NUM_BINS_KE, ARRAY_KE_BINS.data());
            TempVec.push_back(hTempHist);
        }
        TrueAbs0pAsByBinNom.push_back(TempVec);
    }
    std::vector<std::vector<TH1D*>> TrueAbsNpAsByBinNom;
    for (int iEnergyBin = 0; iEnergyBin < NUM_BINS_KE; ++iEnergyBin) {
        std::vector<TH1D*> TempVec;
        for (int iInt = 0; iInt < NUM_SIGNAL_TYPES; ++iInt) {
            TH1D* hTempHist = new TH1D(Form("hTrueAbsNp_%d_Bin_As_%d", iEnergyBin, iInt), Form("True Abs Np KE As %d in bin %d", iInt, iEnergyBin), NUM_BINS_KE, ARRAY_KE_BINS.data());
            TempVec.push_back(hTempHist);
        }
        TrueAbsNpAsByBinNom.push_back(TempVec);
    }
    std::vector<std::vector<TH1D*>> TrueScatterAsByBinNom;
    for (int iEnergyBin = 0; iEnergyBin < NUM_BINS_KE; ++iEnergyBin) {
        std::vector<TH1D*> TempVec;
        for (int iInt = 0; iInt < NUM_SIGNAL_TYPES; ++iInt) {
            TH1D* hTempHist = new TH1D(Form("hTrueScatter_%d_Bin_As_%d", iEnergyBin, iInt), Form("True Scatter KE As %d in bin %d", iInt, iEnergyBin), NUM_BINS_KE, ARRAY_KE_BINS.data());
            TempVec.push_back(hTempHist);
        }
        TrueScatterAsByBinNom.push_back(TempVec);
    }

    std::vector<std::vector<std::vector<TH1D*>>> TrueRecoAsByBinNom = {
        TrueAbs0pAsByBinNom, TrueAbsNpAsByBinNom, TrueScatterAsByBinNom
    };

    //////////////////////
    // Loop over events //
    //////////////////////

    Int_t NumEntries = (Int_t) tree->GetEntries();
    std::cout << "Num entries: " << NumEntries << std::endl;

    for (Int_t i = 0; i < NumEntries; ++i) {
        tree->GetEntry(i);

        // Go faster
        // if (i > USE_NUM_EVENTS) break;

        // Setup hits in tracks
        std::unordered_set<int> hitsInTracks(hitRecoAsTrackKey->begin(), hitRecoAsTrackKey->end());

        // Extend reco cylinder
        removeRepeatedPoints(WC2TPCLocationsX, WC2TPCLocationsY, WC2TPCLocationsZ);
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

        //////////////////////////////////
        // Look at unreconstructed hits //
        //////////////////////////////////

        std::vector<int> candidateHits;
        for (size_t iHit = 0; iHit < fHitKey->size(); ++iHit) {
            // if (fHitPlane->at(iHit) == 1) continue;

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
            )) {
                // Skip hits already in tracks
                if (hitsInTracks.count(iHit) > 0) continue;
                candidateHits.push_back(iHit); 
            }
        }

        // Now cluster using those hits as starting points
        std::unordered_set<int> usedHits;
        std::vector<HitCluster> hitClusters;
        int nCandidateHits = candidateHits.size();

        // Hits in the same cluster must be separated by at most some number of wires
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

        ///////////////////////////////
        // Correct signal definition //
        ///////////////////////////////

        // Scattering only if degree > THRESHOLD and energy > THRESHOLD
        double scatteringAngle = -9999;
        if (backgroundType == 12 || backgroundType == 6) {
            if (backgroundType == 12) scatteringAngle = trajectoryInteractionAngle;
            else if (backgroundType == 6) scatteringAngle = truthScatteringAngle;

            if (scatteringAngle < SCATTERING_ANGLE_THRESHOLD) {
                // Use secondary interaction
                for (int iInteraction = 0; iInteraction < secondaryInteractionTypes->size(); ++iInteraction) {
                    int currentInteraction = secondaryInteractionTypes->at(iInteraction);
                    scatteringAngle        = secondaryInteractionAngle->at(iInteraction);

                    for (int iContribution = 0; iContribution < secondaryIncidentKEContributions->at(iInteraction).size(); ++iContribution) {
                        trueIncidentKEContributions->push_back(secondaryIncidentKEContributions->at(iInteraction)[iContribution]);
                    }

                    if (
                        (currentInteraction == 6 || currentInteraction == 12) &&
                        scatteringAngle < SCATTERING_ANGLE_THRESHOLD
                    ) {
                        // We want to keep going
                        continue;
                    } else {
                        // We found non-scattering interaction or scattering interaction 
                        // with angle above our threshold value
                        backgroundType = currentInteraction;

                        if (backgroundType == 6 || backgroundType == 12) {
                            // Look at outgoing particles
                            int secondaryVisibleProtons = 0; double scatteringEnergy = 0;
                            for (int i = 0; i < secondaryInteractionDaughtersPDG->at(iInteraction).size(); ++i) {
                                if (secondaryInteractionDaughtersPDG->at(iInteraction)[i] == -211) {
                                    scatteringEnergy = secondaryInteractionDaughtersKE->at(iInteraction)[i];
                                } else if (secondaryInteractionDaughtersPDG->at(iInteraction)[i] == 2212) {
                                    if (
                                        secondaryInteractionDaughtersKE->at(iInteraction)[i] > PROTON_ENERGY_LOWER_BOUND &&
                                        secondaryInteractionDaughtersKE->at(iInteraction)[i] < PROTON_ENERGY_UPPER_BOUND
                                    ) secondaryVisibleProtons++;
                                }
                            }

                            // If scattering energy is below threshold, absorption
                            if (scatteringEnergy < PION_SCATTERING_ENERGY_THRESHOLD) {
                                if (secondaryVisibleProtons == 0) backgroundType = 0;
                                else backgroundType = 1;
                            }
                        }
                        truthPrimaryVertexKE = secondaryInteractionInteractingKE->at(iInteraction);
                        truthPrimaryVertexX  = secondaryInteractionXPosition->at(iInteraction);
                        truthPrimaryVertexY  = secondaryInteractionYPosition->at(iInteraction);
                        truthPrimaryVertexZ  = secondaryInteractionZPosition->at(iInteraction);
                        break;
                    }
                }
            } else {
                // If initial interaction has valid angle, first correct for vertex energy 
                if (backgroundType == 12) truthPrimaryVertexKE = trajectoryInteractionKE;

                // Check if outgoing pion is above threshold
                if (backgroundType == 6 && truthScatteredPionKE < PION_SCATTERING_ENERGY_THRESHOLD) {
                    int secondaryVisibleProtons = 0;
                    for (int i = 0; i < truthPrimaryDaughtersPDG->size(); ++i) {
                        if (truthPrimaryDaughtersPDG->at(i) == 2212) {
                            if (
                                truthPrimaryDaughtersKE->at(i) > PROTON_ENERGY_LOWER_BOUND &&
                                truthPrimaryDaughtersKE->at(i) < PROTON_ENERGY_UPPER_BOUND
                            ) secondaryVisibleProtons++;
                        }
                    }
                    if (secondaryVisibleProtons == 0) backgroundType = 0;
                    else backgroundType = 1;
                } else if (backgroundType == 12 && trajectoryInteractionKE < PION_SCATTERING_ENERGY_THRESHOLD) {
                    // If outgoing pion from elastic scattering has energy lower than threshold, we
                    // label it as 0p absorption since there are no other products in an elastic scattering
                    backgroundType = 0;
                }

            }
        }

        ///////////////////////////
        // Get truth information //
        ///////////////////////////

        int TrueEnergyBin = getBin(truthPrimaryVertexKE * 1000, ARRAY_KE_BINS);

        if (validTrueIncidentKE) {
            for (double x : *trueIncidentKEContributions) {
                // Add to nominal
                hTrueIncidentKENom->Fill(x);

                // Add to universes. For generator systematics, we do not vary the 
                // incident flux, just the interacting flux
                for (int iUniv = 0; iUniv < NUM_UNIVERSES; ++iUniv) {
                    TrueIncidentKEUnivs[iUniv]->Fill(x);
                }
            }
        }

        if (backgroundType == 0) {
            hTruePionAbs0pKENom->Fill(truthPrimaryVertexKE * 1000);
            for (int iUniv = 0; iUniv < NUM_UNIVERSES; ++iUniv) TruePionAbs0pKEUnivs[iUniv]->Fill(truthPrimaryVertexKE * 1000, 1);
        } else if (backgroundType == 1) {
            hTruePionAbsNpKENom->Fill(truthPrimaryVertexKE * 1000);
            for (int iUniv = 0; iUniv < NUM_UNIVERSES; ++iUniv) TruePionAbsNpKEUnivs[iUniv]->Fill(truthPrimaryVertexKE * 1000, 1);
        } else if (backgroundType == 6 || backgroundType == 12) {
            hTruePionScatterKENom->Fill(truthPrimaryVertexKE * 1000);
            for (int iUniv = 0; iUniv < NUM_UNIVERSES; ++iUniv) TruePionScatterKEUnivs[iUniv]->Fill(truthPrimaryVertexKE * 1000, 1);
        }

        //////////////////////////////
        // Classification algorithm //
        //////////////////////////////

        // Reject stuff with no data products
        if (WC2TPCtrkID == -99999) continue;

        ////////////////////////////////
        // Cylinder and TG track cuts //
        ////////////////////////////////

        int numSmallTracksInCylinder = 0;
        int numTracksInCylinder      = 0;
        int numTGTracks              = 0;
        for (int trk_idx = 0; trk_idx < recoTrkID->size(); ++trk_idx) {
            if (recoTrkID->at(trk_idx) == WC2TPCtrkID) continue;

            // Check if tracks is through-going
            if (
                !isWithinReducedVolume(recoBeginX->at(trk_idx), recoBeginY->at(trk_idx), recoBeginZ->at(trk_idx)) &&
                !isWithinReducedVolume(recoEndX->at(trk_idx), recoEndY->at(trk_idx), recoEndZ->at(trk_idx))
            ) numTGTracks++;

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

        // Reject events with too many TG tracks
        if (numTGTracks > MAX_NUM_TG_TRACKS) continue;

        // Reject events with too many small tracks inside the cylinder
        if (numSmallTracksInCylinder > ALLOWED_CYLINDER_SMALL_TRACKS) continue;

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
        double calculatedEnLoss = energyLossCalculation(WC4PrimaryX, trajectoryInitialMomentumX, isData);
        const double initialKE  = WCKE - calculatedEnLoss;
        double energyDeposited  = 0;

        for (size_t iDep = 0; iDep < wcMatchDEDX->size(); ++iDep) {
            // If we are past detected breaking point, exit loop
            if (wcMatchZPos->at(iDep) > breakPointZ) break;

            // If larger than threshold, continue
            if (wcMatchDEDX->at(iDep) > HIT_DEDX_THRESHOLD) continue;

            // Else, add to energy deposited so far
            energyDeposited += wcMatchEDep->at(iDep);

            // Add to incident KE if inside reduced volume
            if (isWithinReducedVolume(wcMatchXPos->at(iDep), wcMatchYPos->at(iDep), wcMatchZPos->at(iDep))) {
                hPionIncidentKENom->Fill(initialKE - energyDeposited);
                PionIncidentKEUnivs[0]->Fill(initialKE - energyDeposited);
                PionIncidentKEUnivs[1]->Fill(initialKE - energyDeposited);

                if (truthPrimaryPDG == 11) {
                    PionIncidentKEUnivs[0]->Fill(initialKE - energyDeposited, -ELECTRON_COMP_UNC);
                    PionIncidentKEUnivs[1]->Fill(initialKE - energyDeposited, +ELECTRON_COMP_UNC);
                }

                if (truthPrimaryPDG == -211) {
                    hPionIncidentKETrueNom->Fill(initialKE - energyDeposited);
                    PionIncidentKETrueUnivs[0]->Fill(initialKE - energyDeposited);
                    PionIncidentKETrueUnivs[1]->Fill(initialKE - energyDeposited);
                }
            }
        }

        double energyAtVertex = initialKE - energyDeposited;

        // Reject if vertex outside reduced volume
        if (!isWithinReducedVolume(breakPointX, breakPointY, breakPointZ)) continue;

        // Reject based on primary PID 
        if (minChi2 == pionChi2 || minChi2 == protonChi2) continue;

        /////////////////////////
        // Secondary track PID //
        /////////////////////////

        int secondaryTaggedPion   = 0;
        int secondaryTaggedProton = 0;
        int secondaryTaggedOther  = 0;

        int otherTaggedPion   = 0;
        int otherTaggedProton = 0;

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
                std::pow(recoBeginX->at(trk_idx) - recoEndX->at(trk_idx), 2) +
                std::pow(recoBeginY->at(trk_idx) - recoEndY->at(trk_idx), 2) + 
                std::pow(recoBeginZ->at(trk_idx) - recoEndZ->at(trk_idx), 2)
            );

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
            if (totalTaggedPions > 1 || newSecondaryPion) continue;

            // Add weights to scatter
            hPionScatterKENom->Fill(energyAtVertex);

            PionScatterKEUnivs[0]->Fill(energyAtVertex);
            PionScatterKEUnivs[1]->Fill(energyAtVertex);
            if (truthPrimaryPDG == 11) {
                PionScatterKEUnivs[0]->Fill(energyAtVertex, -ELECTRON_COMP_UNC); 
                PionScatterKEUnivs[1]->Fill(energyAtVertex, +ELECTRON_COMP_UNC);
            }

            if (backgroundType == 0) {
                if (TrueEnergyBin != -1) TrueAbs0pAsByBinNom.at(TrueEnergyBin).at(2)->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                if (TrueEnergyBin != -1) TrueAbsNpAsByBinNom.at(TrueEnergyBin).at(2)->Fill(energyAtVertex);
            } else if (backgroundType == 6 || backgroundType == 12) {
                if (TrueEnergyBin != -1) TrueScatterAsByBinNom.at(TrueEnergyBin).at(2)->Fill(energyAtVertex);
            } else {
                hPionScatterKEBkgNom->Fill(energyAtVertex);

                PionScatterKEBkgUnivs[0]->Fill(energyAtVertex);
                PionScatterKEBkgUnivs[1]->Fill(energyAtVertex);
                if (truthPrimaryPDG == 11) {
                    PionScatterKEBkgUnivs[0]->Fill(energyAtVertex, -ELECTRON_COMP_UNC);
                    PionScatterKEBkgUnivs[1]->Fill(energyAtVertex, +ELECTRON_COMP_UNC);
                }
            }

            continue;
        }

        if (totalTaggedProtons > 0) {
            // Add weights to abs Np
            hPionAbsNpKENom->Fill(energyAtVertex);

            PionAbsNpKEUnivs[0]->Fill(energyAtVertex);
            PionAbsNpKEUnivs[1]->Fill(energyAtVertex);
            if (truthPrimaryPDG == 11) {
                PionAbsNpKEUnivs[0]->Fill(energyAtVertex, -ELECTRON_COMP_UNC); 
                PionAbsNpKEUnivs[1]->Fill(energyAtVertex, +ELECTRON_COMP_UNC);
            }

            if (backgroundType == 0) {
                if (TrueEnergyBin != -1) TrueAbs0pAsByBinNom.at(TrueEnergyBin).at(1)->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                if (TrueEnergyBin != -1) TrueAbsNpAsByBinNom.at(TrueEnergyBin).at(1)->Fill(energyAtVertex);
            } else if (backgroundType == 6 || backgroundType == 12) {
                if (TrueEnergyBin != -1) TrueScatterAsByBinNom.at(TrueEnergyBin).at(1)->Fill(energyAtVertex);
            } else {
                hPionAbsNpKEBkgNom->Fill(energyAtVertex);

                PionAbsNpKEBkgUnivs[0]->Fill(energyAtVertex);
                PionAbsNpKEBkgUnivs[1]->Fill(energyAtVertex);
                if (truthPrimaryPDG == 11) {
                    PionAbsNpKEBkgUnivs[0]->Fill(energyAtVertex, -ELECTRON_COMP_UNC);
                    PionAbsNpKEBkgUnivs[1]->Fill(energyAtVertex, +ELECTRON_COMP_UNC);
                }
            }
            
            continue;
        }

        ////////////////////////////////////////
        // Cluster non-reconstructed hits cut //
        ////////////////////////////////////////

        // Get data for cut
        int numLargeClustersInduction  = 0;
        int numLargeClustersCollection = 0;

        int numClustersInduction  = 0;
        int numClustersCollection = 0;

        double largestClusterSizeInduction  = 0;
        double largestClusterSizeCollection = 0;

        for (int i = 0; i < hitClusters.size(); ++i) {
            double clusterSize = hitClusters[i].clusterSize;

            if (hitClusters[i].plane == 0) {
                if (clusterSize > largestClusterSizeInduction) largestClusterSizeInduction = clusterSize;
                numClustersInduction++;
            }
            else if (hitClusters[i].plane == 1) {
                if (clusterSize > largestClusterSizeCollection) largestClusterSizeCollection = clusterSize;
                numClustersCollection++;
            }

            if (clusterSize > LARGE_CLUSTER_THRESHOLD) {
                if (hitClusters[i].plane == 0) numLargeClustersInduction++;
                else if (hitClusters[i].plane == 1) numLargeClustersCollection++;
            }
        }

        if (numClustersInduction < MAX_NUM_CLUSTERS_INDUCTION) {
            // Add weights to abs 0p
            hPionAbs0pKENom->Fill(energyAtVertex);

            PionAbs0pKEUnivs[0]->Fill(energyAtVertex);
            PionAbs0pKEUnivs[1]->Fill(energyAtVertex);
            if (truthPrimaryPDG == 11) {
                PionAbs0pKEUnivs[0]->Fill(energyAtVertex, -ELECTRON_COMP_UNC); 
                PionAbs0pKEUnivs[1]->Fill(energyAtVertex, +ELECTRON_COMP_UNC);
            }

            if (backgroundType == 0) {
                if (TrueEnergyBin != -1) TrueAbs0pAsByBinNom.at(TrueEnergyBin).at(0)->Fill(energyAtVertex);
            } else if (backgroundType == 1) {
                if (TrueEnergyBin != -1) TrueAbsNpAsByBinNom.at(TrueEnergyBin).at(0)->Fill(energyAtVertex);
            } else if (backgroundType == 6 || backgroundType == 12) {
                if (TrueEnergyBin != -1) TrueScatterAsByBinNom.at(TrueEnergyBin).at(0)->Fill(energyAtVertex);
            } else {
                hPionAbs0pKEBkgNom->Fill(energyAtVertex);

                PionAbs0pKEBkgUnivs[0]->Fill(energyAtVertex);
                PionAbs0pKEBkgUnivs[1]->Fill(energyAtVertex);
                if (truthPrimaryPDG == 11) {
                    PionAbs0pKEBkgUnivs[0]->Fill(energyAtVertex, -ELECTRON_COMP_UNC);
                    PionAbs0pKEBkgUnivs[1]->Fill(energyAtVertex, +ELECTRON_COMP_UNC);
                }
            }

            continue;
        }
    }

    //////////////////////////////
    // Cross-section extraction //
    //////////////////////////////

    ///////////////////
    // Incident flux //
    ///////////////////
    
    // Get nominal incident flux and corrections
    for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
        double psi = hPionIncidentKETrueNom->GetBinContent(iBin + 1) / hTrueIncidentKENom->GetBinContent(iBin + 1);
        double   c = hPionIncidentKETrueNom->GetBinContent(iBin + 1) / hPionIncidentKENom->GetBinContent(iBin + 1);

        hPsiIncidentNom->SetBinContent(iBin + 1, psi); hCIncidentNom->SetBinContent(iBin + 1, c);
        hPionIncidentCorrectedKENom->SetBinContent(iBin + 1, (c / psi) * hPionIncidentKENom->GetBinContent(iBin + 1));
    }

    // Get incident flux and corrections for all universes
    for (int iUniv = 0; iUniv < NUM_UNIVERSES; ++iUniv) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            double psi = PionIncidentKETrueUnivs[iUniv]->GetBinContent(iBin + 1) / TrueIncidentKEUnivs[iUniv]->GetBinContent(iBin + 1);
            double   c = PionIncidentKETrueUnivs[iUniv]->GetBinContent(iBin + 1) / PionIncidentKEUnivs[iUniv]->GetBinContent(iBin + 1);

            PsiIncidentUnivs[iUniv]->SetBinContent(iBin + 1, psi); CIncidentUnivs[iUniv]->SetBinContent(iBin + 1, c);
            PionIncidentCorrectedKEUnivs[iUniv]->SetBinContent(iBin + 1, (c / psi) * PionIncidentKEUnivs[iUniv]->GetBinContent(iBin + 1));
        }
    }

    ///////////////////////
    // Response matrices //
    ///////////////////////

    // Get nominal response matrix
    TH2D* hResponseMatrixNominal = new TH2D(
        "hResponseMatrix", "hResponseMatrix;Reco (j, #beta);True (i, #alpha)",
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE,
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    );
    GetResponseMatrix(
        NUM_SIGNAL_TYPES, NUM_BINS_KE,
        TotalEventsHistosNom,
        TrueRecoAsByBinNom,
        hPionIncidentKENom,
        hTrueIncidentKENom,
        hResponseMatrixNominal
    );

    // Compute response matrix for all universes
    std::vector<TH2D*> ResponseMatricesUnivs;
    for (int iUniv = 0; iUniv < NUM_UNIVERSES; ++iUniv) {
        TH2D* hResponseMatrixUniv = new TH2D(
            "hResponseMatrix", "hResponseMatrix;Reco (j, #beta);True (i, #alpha)",
            NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE,
            NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
        );
        std::vector<TH1D*> TotalEventsUniv = {hTruePionAbs0pKENom, hTruePionAbsNpKENom, hTruePionScatterKENom};
        std::vector<std::vector<std::vector<TH1D*>>> TrueRecoAsByBinUniv = {TrueAbs0pAsByBinNom, TrueAbsNpAsByBinNom, TrueScatterAsByBinNom};
        TH1D* RecoIncidentFluxUniv = PionIncidentKEUnivs[iUniv];
        TH1D* TrueIncidentFluxUniv = hTrueIncidentKENom;
        
        GetResponseMatrix(
            NUM_SIGNAL_TYPES, NUM_BINS_KE,
            TotalEventsUniv,
            TrueRecoAsByBinUniv,
            RecoIncidentFluxUniv,
            TrueIncidentFluxUniv,
            hResponseMatrixUniv
        );

        ResponseMatricesUnivs.push_back(hResponseMatrixUniv);
    }

    //////////////////////////
    // Signal cross-section //
    //////////////////////////

    // Get nominal signal vector
    TH1D* hSignalNominal = new TH1D("hSignalNominal", "hSignalNominal;;", NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            double xsec =  XSEC_UNITS * (TotalEventsHistosNom[iSignal]->GetBinContent(iBin + 1) / hTrueIncidentKENom->GetBinContent(iBin + 1));
            hSignalNominal->SetBinContent(index + 1,  xsec); TotalEventsHistosNom[iSignal]->SetBinContent(iBin + 1,  xsec);
        }
    }

    //////////////////////////////
    // Background cross-section //
    //////////////////////////////

    // Get nominal estimated background vector
    TH1D* hBackgroundNominal = new TH1D("hBackgroundNominal", "hBackgroundNominal;;", NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            double xsec = XSEC_UNITS * (RecoBkgSignals[iSignal]->GetBinContent(iBin + 1) / hPionIncidentKENom->GetBinContent(iBin + 1));
            hBackgroundNominal->SetBinContent(index + 1, xsec); RecoBkgSignals[iSignal]->SetBinContent(iBin + 1, xsec);
        }
    }

    // Compose background vectors for all universes
    std::vector<TH1D*> BackgroundVectorUnivs;
    BackgroundVectorUnivs.reserve(NUM_UNIVERSES);
    for (int iUniv = 0; iUniv < NUM_UNIVERSES; ++iUniv) {
        TH1D* hBackground = new TH1D(Form("hBackground_Univ%d", iUniv), Form("Background Vector for Universe %d", iUniv), NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE);
        for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
            for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
                int index = flattenIndex(iSignal, iBin, NUM_BINS_KE);

                double content = 0;
                if (iSignal == 0) {
                    content = XSEC_UNITS * (PionAbs0pKEBkgUnivs[iUniv]->GetBinContent(iBin + 1) / PionIncidentKEUnivs[iUniv]->GetBinContent(iBin + 1));
                    PionAbs0pKEBkgUnivs[iUniv]->SetBinContent(iBin + 1, content);
                } else if (iSignal == 1) {
                    content = XSEC_UNITS * (PionAbsNpKEBkgUnivs[iUniv]->GetBinContent(iBin + 1) / PionIncidentKEUnivs[iUniv]->GetBinContent(iBin + 1));
                    PionAbsNpKEBkgUnivs[iUniv]->SetBinContent(iBin + 1, content);
                } else if (iSignal == 2) {
                    content = XSEC_UNITS * (PionScatterKEBkgUnivs[iUniv]->GetBinContent(iBin + 1) / PionIncidentKEUnivs[iUniv]->GetBinContent(iBin + 1));
                    PionScatterKEBkgUnivs[iUniv]->SetBinContent(iBin + 1, content);
                }

                hBackground->SetBinContent(index + 1, content);
                hBackground->SetBinError(index + 1, 0.0);
            }
        }
        BackgroundVectorUnivs.push_back(hBackground);
    }

    ////////////////////////////
    // Measured cross-section //
    ////////////////////////////

    // Get nominal measure vector
    TH1D* hMeasureNominal = new TH1D("hMeasureVectorNominal", "hMeasureVectorNominal;;", NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            double xsec = XSEC_UNITS * (RecoSignals[iSignal]->GetBinContent(iBin + 1) / hPionIncidentKENom->GetBinContent(iBin + 1));
            hMeasureNominal->SetBinContent(index + 1, xsec); RecoSignals[iSignal]->SetBinContent(iBin + 1, xsec);
        }
    }

    // Get nominal true cross-section in vector form
    TVectorD SignalNom(NUM_SIGNAL_TYPES * NUM_BINS_KE); H2V(hSignalNominal, SignalNom);

    // Compose measured vectors for all universes
    std::vector<TH1D*> MeasureVectorUnivs;
    MeasureVectorUnivs.reserve(NUM_UNIVERSES);
    for (int iUniv = 0; iUniv < NUM_UNIVERSES; ++iUniv) {
        // One TH1D holding all signals flattened
        TH1D* hMeasure = new TH1D(Form("hMeasure_Univ%d", iUniv), Form("Measured Vector for Universe %d", iUniv), NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE);

        // Compute M = R * S + B with S fixed
        TMatrixD R(NUM_SIGNAL_TYPES * NUM_BINS_KE, NUM_SIGNAL_TYPES * NUM_BINS_KE); H2M(static_cast<const TH2D*>(ResponseMatricesUnivs[iUniv]), R, kTRUE);
        TVectorD B(NUM_SIGNAL_TYPES * NUM_BINS_KE); H2V(BackgroundVectorUnivs[iUniv], B);

        TVectorD M = (R * SignalNom) + B;

        // Fill it with flattened (signal, bin) content
        for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
            for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
                int index     = flattenIndex(iSignal, iBin, NUM_BINS_KE);

                // double content = 0;
                if (iSignal == 0) {
                    // content = XSEC_UNITS * (PionAbs0pKEUnivs[iUniv]->GetBinContent(iBin + 1) / PionIncidentCorrectedKEUnivs[iUniv]->GetBinContent(iBin + 1));
                    PionAbs0pKEUnivs[iUniv]->SetBinContent(iBin + 1, M(index));
                    // std::cout << "measured: " << content << " vs. vector: " << M(index) << std::endl;
                } else if (iSignal == 1) {
                    // content = XSEC_UNITS * (PionAbsNpKEUnivs[iUniv]->GetBinContent(iBin + 1) / PionIncidentCorrectedKEUnivs[iUniv]->GetBinContent(iBin + 1));
                    PionAbsNpKEUnivs[iUniv]->SetBinContent(iBin + 1, M(index));
                    // std::cout << "measured: " << content << " vs. vector: " << M(index) << std::endl;
                } else if (iSignal == 2) {
                    // content = XSEC_UNITS * (PionScatterKEUnivs[iUniv]->GetBinContent(iBin + 1) / PionIncidentCorrectedKEUnivs[iUniv]->GetBinContent(iBin + 1));
                    PionScatterKEUnivs[iUniv]->SetBinContent(iBin + 1, M(index));
                    // std::cout << "measured: " << content << " vs. vector: " << M(index) << std::endl;
                }
                
                hMeasure->SetBinContent(index + 1, M(index));
                hMeasure->SetBinError(index + 1, 0.0);
            }
        }
        MeasureVectorUnivs.push_back(hMeasure);
    }

    // Use TMatrices/Vectors to compute response equation
    TMatrixD ResponseNom(NUM_SIGNAL_TYPES * NUM_BINS_KE, NUM_SIGNAL_TYPES * NUM_BINS_KE); H2M(static_cast<const TH2D*>(hResponseMatrixNominal), ResponseNom, kTRUE);
    TVectorD BackgroundNom(NUM_SIGNAL_TYPES * NUM_BINS_KE); H2V(hBackgroundNominal, BackgroundNom);

    // Get measured vector by "folding" true vector
    TVectorD MeasureNom = (ResponseNom * SignalNom) + BackgroundNom;

    // Sanity check, for CV, these two should match
    std::cout << "Cross-section: " << std::endl;
    TH1D* hCrossSectionNominal = new TH1D("hCrossSectionNominal", "hCrossSectionNominal;;", NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            std::cout << "  Reco: " << index << " " << hMeasureNominal->GetBinContent(index + 1) << std::endl;
            std::cout << "  Fold: " << index << " " << MeasureNom(index) << std::endl;
            std::cout << std::endl;
        }
    }

    /////////////////////
    // Get covariances //
    /////////////////////

    // Meausured vector
    TH2D* hMeasureCovMatrix = new TH2D(
        "hMeasureCovMatrix", "hMeasureCovMatrix;Reco (j, #beta);True (i, #alpha)",
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE,
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    );
    TH2D* hMeasureFracCovMatrix = new TH2D(
        "hMeasureFracCovMatrix", "hMeasureFracCovMatrix;Reco (j, #beta);True (i, #alpha)",
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE,
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    );
    TH2D* hMeasureCorrMatrix = new TH2D(
        "hMeasureCorrMatrix", "hMeasureCorrMatrix;Reco (j, #beta);True (i, #alpha)",
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE,
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    );
    GetCovMatrix(hMeasureNominal, MeasureVectorUnivs, hMeasureCovMatrix);
    GetFracCovAndCorrMatrix(hMeasureNominal, hMeasureCovMatrix, hMeasureFracCovMatrix, hMeasureCorrMatrix);
    DrawHistosWithErrorBands(
        hMeasureNominal,
        MeasureVectorUnivs,
        SaveDir,
        "BeamlineEl",
        "Measure",
        TextSize,
        FontStyle,
        "Measured Cross Section",
        "Bin number",
        "Cross Section [barn]"
    );

    // Abs 0p
    TH2D* hAbs0pCovMatrix = new TH2D(
        "hAbs0pCovMatrix", "hAbs0pCovMatrix;Kinetic Energy [MeV];Kinetic Energy [MeV]",
        NUM_BINS_KE, ARRAY_KE_BINS.data(),
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    TH2D* hAbs0pFracCovMatrix = new TH2D(
        "hAbs0pFracCovMatrix", "hAbs0pFracCovMatrix;Kinetic Energy [MeV];Kinetic Energy [MeV]",
        NUM_BINS_KE, ARRAY_KE_BINS.data(),
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    TH2D* hAbs0pCorrMatrix = new TH2D(
        "hAbs0pCorrMatrix", "hAbs0pCorrMatrix;Kinetic Energy [MeV];Kinetic Energy [MeV]",
        NUM_BINS_KE, ARRAY_KE_BINS.data(),
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    GetCovMatrix(hPionAbs0pKENom, PionAbs0pKEUnivs, hAbs0pCovMatrix);
    GetFracCovAndCorrMatrix(hPionAbs0pKENom, hAbs0pCovMatrix, hAbs0pFracCovMatrix, hAbs0pCorrMatrix);
    DrawHistosWithErrorBands(
        hPionAbs0pKENom, 
        PionAbs0pKEUnivs, 
        SaveDir, 
        "BeamlineEl", 
        "Abs0p",
        TextSize,
        FontStyle,
        "Abs 0p Cross Section",
        "Kinetic Energy [MeV]",
        "Cross Section [barn]"
    );

    // Abs Np
    TH2D* hAbsNpCovMatrix = new TH2D(
        "hAbsNpCovMatrix", "hAbsNpCovMatrix;Kinetic Energy [MeV];Kinetic Energy [MeV]",
        NUM_BINS_KE, ARRAY_KE_BINS.data(),
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    TH2D* hAbsNpFracCovMatrix = new TH2D(
        "hAbsNpFracCovMatrix", "hAbsNpFracCovMatrix;Kinetic Energy [MeV];Kinetic Energy [MeV]",
        NUM_BINS_KE, ARRAY_KE_BINS.data(),
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    TH2D* hAbsNpCorrMatrix = new TH2D(
        "hAbsNpCorrMatrix", "hAbsNpCorrMatrix;Kinetic Energy [MeV];Kinetic Energy [MeV]",
        NUM_BINS_KE, ARRAY_KE_BINS.data(),
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    GetCovMatrix(hPionAbsNpKENom, PionAbsNpKEUnivs, hAbsNpCovMatrix);
    GetFracCovAndCorrMatrix(hPionAbsNpKENom, hAbsNpCovMatrix, hAbsNpFracCovMatrix, hAbsNpCorrMatrix);
    DrawHistosWithErrorBands(
        hPionAbsNpKENom, 
        PionAbsNpKEUnivs, 
        SaveDir,
        "BeamlineEl", 
        "AbsNp",
        TextSize,
        FontStyle,
        "Abs Np Cross Section",
        "Kinetic Energy [MeV]",
        "Cross Section [barn]"
    );

    // Scatter
    TH2D* hScatterCovMatrix = new TH2D(
        "hScatterCovMatrix", "hScatterCovMatrix;Kinetic Energy [MeV];Kinetic Energy [MeV]",
        NUM_BINS_KE, ARRAY_KE_BINS.data(),
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    TH2D* hScatterFracCovMatrix = new TH2D(
        "hScatterFracCovMatrix", "hScatterFracCovMatrix;Kinetic Energy [MeV];Kinetic Energy [MeV]",
        NUM_BINS_KE, ARRAY_KE_BINS.data(),
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    TH2D* hScatterCorrMatrix = new TH2D(
        "hScatterCorrMatrix", "hScatterCorrMatrix;Kinetic Energy [MeV];Kinetic Energy [MeV]",
        NUM_BINS_KE, ARRAY_KE_BINS.data(),
        NUM_BINS_KE, ARRAY_KE_BINS.data()
    );
    GetCovMatrix(hPionScatterKENom, PionScatterKEUnivs, hScatterCovMatrix);
    GetFracCovAndCorrMatrix(hPionScatterKENom, hScatterCovMatrix, hScatterFracCovMatrix, hScatterCorrMatrix);
    DrawHistosWithErrorBands(
        hPionScatterKENom, 
        PionScatterKEUnivs, 
        SaveDir, 
        "BeamlineEl",
        "Scatter",
        TextSize,
        FontStyle,
        "Scatter Cross Section",
        "Kinetic Energy [MeV]",
        "Cross Section [barn]"
    );

    ///////////////////
    // Draw 2D Plots //
    ///////////////////

    std::vector<TH2*> TwoDPlots = {
        hMeasureCovMatrix,
        hMeasureFracCovMatrix, 
        hMeasureCorrMatrix,
        hAbs0pCovMatrix,
        hAbs0pFracCovMatrix,
        hAbs0pCorrMatrix,
        hAbsNpCovMatrix,
        hAbsNpFracCovMatrix,
        hAbsNpCorrMatrix,
        hScatterCovMatrix,
        hScatterFracCovMatrix,
        hScatterCorrMatrix
    };

    std::vector<TString> TwoDTitles = {
        "BeamlineEl/MeasureCovariance",
        "BeamlineEl/MeasureFracCovariance",
        "BeamlineEl/MeasureCorrelation",
        "BeamlineEl/Abs0pCovariance",
        "BeamlineEl/Abs0pFracCovariance",
        "BeamlineEl/Abs0pCorrelation",
        "BeamlineEl/AbsNpCovariance",
        "BeamlineEl/AbsNpFracCovariance",
        "BeamlineEl/AbsNpCorrelation",
        "BeamlineEl/ScatterCovariance",
        "BeamlineEl/ScatterFracCovariance",
        "BeamlineEl/ScatterCorrelation"
    };

    std::vector<std::pair<double,double>> TwoDRanges = {
        {0, 0},
        {0, 0},
        {-1, 1},
        {0, 0},
        {0, 0},
        {-1, 1},
        {0, 0},
        {0, 0},
        {-1, 1},
        {0, 0},
        {0, 0},
        {-1, 1}
    };

    std::vector<bool> TwoDDisplayNumbers = {
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
        false
    };

    printTwoDPlots(SaveDir, TwoDPlots, TwoDTitles, TwoDRanges, TwoDDisplayNumbers);

    ////////////////////////////
    // Save covariance matrix //
    ////////////////////////////

    TFile* saveFile = new TFile("/exp/lariat/app/users/epelaez/histos/beamlineel/Measure.root", "RECREATE");
    saveFile->cd();
    hMeasureCovMatrix->Write("", TObject::kOverwrite);
    saveFile->Close();
}