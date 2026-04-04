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

void DataClassify() {
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
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis/figs/DataClassify/";

    // Load file with NN data products
    TString RootFilePath = "/exp/lariat/app/users/epelaez/files/DataNeg60_histo.root";
    std::unique_ptr<TFile> File(TFile::Open(RootFilePath));
    TDirectory* Directory = (TDirectory*)File->Get("RecoNNDataEval");

    // Load nominal MC histos
    TString NominalHistsPath = "/exp/lariat/app/users/epelaez/histos/nominal/RecoClassify3Cat_AllHists.root";
    std::unique_ptr<TFile> fNom(TFile::Open(NominalHistsPath, "READ"));
    
    ///////////////////////
    // Counters for cuts //
    ///////////////////////

    int TotalEvents = 0;
    int EventsPassingProj = 0;
    int EventsPassingAperture = 0;
    int EventsPassingTOFMass = 0;
    int EventsWithWCMatch = 0;
    int EventsPassingTG = 0;
    int EventsPassingSmall = 0;
    int EventsInRedVol = 0;
    int EventsPassingPrimaryPID = 0;
    int EventsSelectedAsScatter = 0;
    int EventsSelectedAsAbsNp = 0;
    int EventsNoSecondary = 0;
    int EventsSelectedAsAbs0p = 0;

    //////////////////////////
    // Get histos for plots //
    //////////////////////////

    TH1D* hMCIncidentKEPion     = dynamic_cast<TH1D*>(fNom->Get("hIncidentKEPion"));
    TH1D* hMCIncidentKEElectron = dynamic_cast<TH1D*>(fNom->Get("hIncidentKEElectron"));
    TH1D* hMCIncidentKEMuon     = dynamic_cast<TH1D*>(fNom->Get("hIncidentKEMuon"));

    TH1D* hMCIncidentKEPionFine     = dynamic_cast<TH1D*>(fNom->Get("hIncidentKEPionFine"));
    TH1D* hMCIncidentKEElectronFine = dynamic_cast<TH1D*>(fNom->Get("hIncidentKEElectronFine"));
    TH1D* hMCIncidentKEMuonFine     = dynamic_cast<TH1D*>(fNom->Get("hIncidentKEMuonFine"));

    TH1D* hMCSmallTrksInCylinderPions     = dynamic_cast<TH1D*>(fNom->Get("hSmallTrksInCylinderPions"));
    TH1D* hMCSmallTrksInCylinderElectrons = dynamic_cast<TH1D*>(fNom->Get("hSmallTrksInCylinderElectrons"));
    TH1D* hMCSmallTrksInCylinderMuons     = dynamic_cast<TH1D*>(fNom->Get("hSmallTrksInCylinderMuons"));

    TH1D* hMCFrontFaceKEPion     = dynamic_cast<TH1D*>(fNom->Get("hFrontFaceKEPion"));
    TH1D* hMCFrontFaceKEElectron = dynamic_cast<TH1D*>(fNom->Get("hFrontFaceKEElectron"));
    TH1D* hMCFrontFaceKEMuon     = dynamic_cast<TH1D*>(fNom->Get("hFrontFaceKEMuon"));

    TH1D* hMCWCKEPion     = dynamic_cast<TH1D*>(fNom->Get("hWCKEPion"));
    TH1D* hMCWCKEMuon     = dynamic_cast<TH1D*>(fNom->Get("hWCKEMuon"));
    TH1D* hMCWCKEElectron = dynamic_cast<TH1D*>(fNom->Get("hWCKEElectron"));

    TH1D* hMCWCHitsPion     = dynamic_cast<TH1D*>(fNom->Get("hWCHitsPion"));
    TH1D* hMCWCHitsMuon     = dynamic_cast<TH1D*>(fNom->Get("hWCHitsMuon"));
    TH1D* hMCWCHitsElectron = dynamic_cast<TH1D*>(fNom->Get("hWCHitsElectron"));

    TH1D* hMCRadDistWC4Pion     = dynamic_cast<TH1D*>(fNom->Get("hRadDistWC4Pion"));
    TH1D* hMCRadDistWC4Muon     = dynamic_cast<TH1D*>(fNom->Get("hRadDistWC4Muon"));
    TH1D* hMCRadDistWC4Electron = dynamic_cast<TH1D*>(fNom->Get("hRadDistWC4Electron"));

    TH1D* hMCRadDistMidPlanePion     = dynamic_cast<TH1D*>(fNom->Get("hRadDistMidPlanePion"));
    TH1D* hMCRadDistMidPlaneMuon     = dynamic_cast<TH1D*>(fNom->Get("hRadDistMidPlaneMuon"));
    TH1D* hMCRadDistMidPlaneElectron = dynamic_cast<TH1D*>(fNom->Get("hRadDistMidPlaneElectron"));

    // Abs 0p interacting KE
    TH1D* hMCPionAbs0pKETrue     = dynamic_cast<TH1D*>(fNom->Get("hPionAbs0pKETrue"));
    TH1D* hMCPionAbs0pKEAbsNp    = dynamic_cast<TH1D*>(fNom->Get("hPionAbs0pKEAbsNp"));
    TH1D* hMCPionAbs0pKEScatter  = dynamic_cast<TH1D*>(fNom->Get("hPionAbs0pKEScatter"));
    TH1D* hMCPionAbs0pKEChExch   = dynamic_cast<TH1D*>(fNom->Get("hPionAbs0pKEChExch"));
    TH1D* hMCPionAbs0pKEMuon     = dynamic_cast<TH1D*>(fNom->Get("hPionAbs0pKEMuon"));
    TH1D* hMCPionAbs0pKEElectron = dynamic_cast<TH1D*>(fNom->Get("hPionAbs0pKEElectron"));
    TH1D* hMCPionAbs0pKEOther    = dynamic_cast<TH1D*>(fNom->Get("hPionAbs0pKEOther"));

    // Abs Np interacting KE
    TH1D* hMCPionAbsNpKETrue     = dynamic_cast<TH1D*>(fNom->Get("hPionAbsNpKETrue"));
    TH1D* hMCPionAbsNpKEAbs0p    = dynamic_cast<TH1D*>(fNom->Get("hPionAbsNpKEAbs0p"));
    TH1D* hMCPionAbsNpKEScatter  = dynamic_cast<TH1D*>(fNom->Get("hPionAbsNpKEScatter"));
    TH1D* hMCPionAbsNpKEChExch   = dynamic_cast<TH1D*>(fNom->Get("hPionAbsNpKEChExch"));
    TH1D* hMCPionAbsNpKEMuon     = dynamic_cast<TH1D*>(fNom->Get("hPionAbsNpKEMuon"));
    TH1D* hMCPionAbsNpKEElectron = dynamic_cast<TH1D*>(fNom->Get("hPionAbsNpKEElectron"));
    TH1D* hMCPionAbsNpKEOther    = dynamic_cast<TH1D*>(fNom->Get("hPionAbsNpKEOther"));

    // Scatter interacting KE
    TH1D* hMCPionScatterKETrue     = dynamic_cast<TH1D*>(fNom->Get("hPionScatterKETrue"));
    TH1D* hMCPionScatterKEAbs0p    = dynamic_cast<TH1D*>(fNom->Get("hPionScatterKEAbs0p"));
    TH1D* hMCPionScatterKEAbsNp    = dynamic_cast<TH1D*>(fNom->Get("hPionScatterKEAbsNp"));
    TH1D* hMCPionScatterKEChExch   = dynamic_cast<TH1D*>(fNom->Get("hPionScatterKEChExch"));
    TH1D* hMCPionScatterKEMuon     = dynamic_cast<TH1D*>(fNom->Get("hPionScatterKEMuon"));
    TH1D* hMCPionScatterKEElectron = dynamic_cast<TH1D*>(fNom->Get("hPionScatterKEElectron"));
    TH1D* hMCPionScatterKEOther    = dynamic_cast<TH1D*>(fNom->Get("hPionScatterKEOther"));

    ///////////////////
    // Load branches //
    ///////////////////

    // Load tree and branches
    TTree* tree = (TTree*) Directory->Get<TTree>("RecoNNDataEvalTree");

    int run, subrun, event; bool isData;
    tree->SetBranchAddress("run", &run);
    tree->SetBranchAddress("subrun", &subrun);
    tree->SetBranchAddress("event", &event);
    tree->SetBranchAddress("isData", &isData);

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

    // WC quality data
    int wcNumHits;
    std::vector<double>* wcHit0 = nullptr;
    std::vector<double>* wcHit1 = nullptr;
    std::vector<double>* wcHit2 = nullptr;
    std::vector<double>* wcHit3 = nullptr;
    tree->SetBranchAddress("wcNumHits", &wcNumHits);
    tree->SetBranchAddress("wcHit0", &wcHit0);
    tree->SetBranchAddress("wcHit1", &wcHit1);
    tree->SetBranchAddress("wcHit2", &wcHit2);
    tree->SetBranchAddress("wcHit3", &wcHit3);

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

    // Calorimetry information
    std::vector<std::vector<double>>* recoResR = nullptr;
    std::vector<std::vector<double>>* recoDEDX = nullptr;
    tree->SetBranchAddress("recoResR", &recoResR);
    tree->SetBranchAddress("recoDEDX", &recoDEDX);

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

    // Information about hits in tracks
    std::vector<std::vector<int>>* recoTrackHitIndices  = nullptr;
    std::vector<std::vector<double>>*     recoTrackHitX = nullptr; 
    std::vector<std::vector<double>>*     recoTrackHitY = nullptr; 
    std::vector<std::vector<double>>*     recoTrackHitZ = nullptr; 
    tree->SetBranchAddress("recoTrackHitIndices", &recoTrackHitIndices);
    tree->SetBranchAddress("recoTrackHitX", &recoTrackHitX);
    tree->SetBranchAddress("recoTrackHitY", &recoTrackHitY);
    tree->SetBranchAddress("recoTrackHitZ", &recoTrackHitZ);

    // TOF information
    double TOFMass, tofObject;
    tree->SetBranchAddress("TOFMass", &TOFMass);
    tree->SetBranchAddress("tofObject", &tofObject);

    //////////////////////
    // Create histogram //
    //////////////////////

    // Small tracks in cylinder
    TH1D* hSmallTracksInCylinder = new TH1D("hSmallTracksInCylinder ", "hSmallTracksInCylinder ;;", 10, 0, 10);

    // Selected interaction type
    TH1D* hPionAbs0pKE   = new TH1D("hPionAbs0pKE", "hPionAbs0pKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbsNpKE   = new TH1D("hPionAbsNpKE", "hPionAbsNpKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionScatterKE = new TH1D("hPionScatterKE", "hPionScatterKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    std::vector<TH1*> RecoSignals = {
        hPionAbs0pKE, hPionAbsNpKE, hPionScatterKE
    };

    // Also want to keep track of total abs
    TH1D* hPionAbsKE = new TH1D("hPionAbsKE", "hPionAbsKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    std::vector<TH1*> RecoSignalsAbsScatt = {
        hPionAbsKE, hPionScatterKE
    };

    // Incident kinetic energy
    TH1D* hIncidentKE     = new TH1D("hIncidentKE", "hIncidentKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hIncidentKEFine = new TH1D("hIncidentKEFine", "hIncidentKEFine;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());

    // Kinetic energy at the front face of the TPC
    TH1D* hFrontFaceKE = new TH1D("hFrontFaceKE", "hFrontFaceKE;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());

    // Wire chamber kinetic energy (before energy loss correction)
    TH1D* hWCKE = new TH1D("hWCKE", "hWCKE;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());

    // Wire-chamber quality
    TH1D* hWCHits          = new TH1D("hWCHits", "hWCHits", 10, 0, 10);
    TH1D* hRadDistMidPlane = new TH1D("hRadDistMidPlane", "hRadDistMidPlane", 60, 0.0, 20.);
    TH1D* hRadDistWC4      = new TH1D("hRadDistWC4", "RadDistWC4", 90, 0.0, 30.);

    //////////////////////
    // Loop over events //
    //////////////////////

    Int_t NumEntries = (Int_t) tree->GetEntries();
    std::cout << "Num entries: " << NumEntries << std::endl;

    for (Int_t i = 0; i < NumEntries; ++i) {
        tree->GetEntry(i);

        // Make script go faster
        // if (i > USE_NUM_EVENTS) break;

        TotalEvents++;

        /////////////////////////////////////////
        // Wire-chamber and other initial cuts //
        /////////////////////////////////////////
        
        hWCHits->Fill(wcNumHits);
        // if (wcNumHits != 8) continue;

        // std::cout << "==== Entry " << i << " ====\n";
        // std::cout << "wcNumHits = " << wcNumHits << "\n";

        // auto printVec = [](const std::string& name, const std::vector<double>* v) {
        //     std::cout << name << " = ";
        //     if (!v) {
        //         std::cout << "nullptr\n";
        //         return;
        //     }

        //     std::cout << "[ ";
        //     for (size_t j = 0; j < v->size(); ++j) {
        //         std::cout << v->at(j);
        //         if (j + 1 < v->size()) std::cout << ", ";
        //     }
        //     std::cout << " ]\n";
        // };

        // printVec("wcHit0", wcHit0);
        // printVec("wcHit1", wcHit1);
        // printVec("wcHit2", wcHit2);
        // printVec("wcHit3", wcHit3);

        // Project downstream to midplane @ -437.97 (without angular corrections)
        std::vector<double> midUp = projToZ(*wcHit0, *wcHit1, -437.97);
        // Use this point and WC3 to project up to WC4
        std::vector<double> projDown = projToZ(midUp, *wcHit2, -95.0);
        // Requires some corrections because magnets are not the same
        projDown[0] -= tan(1.32 * TMath::Pi() / 180.0) * (-95.0 - -437.97);

        // Compare x and y coordinate in projection and real hit for WC4
        double radDistWC4 = TMath::Sqrt(pow(projDown[0] - wcHit3->at(0), 2.) + pow(projDown[1] - wcHit3->at(1), 2.));

        // Project upstream to midplane @ -437.97
        std::vector<double> midDown = projToZ(*wcHit2, *wcHit3, -437.97);
        midDown[0] -= tan(1.32 * TMath::Pi() / 180.0) * (-339.57 - -437.97);
        double midPlaneDist = TMath::Sqrt(pow(midUp[0] - midDown[0] + 0.75, 2) + pow(midUp[1] - midDown[1], 2));
        // double midPlaneDist = TMath::Sqrt(pow(midUp[0] - midDown[0], 2) + pow(midUp[1] - midDown[1], 2));

        hRadDistMidPlane->Fill(midPlaneDist);
        hRadDistWC4->Fill(radDistWC4);

        // Cuts
        if (radDistWC4 > 8.0) continue;
        if (midPlaneDist > 3.0) continue;
        EventsPassingProj++;

        // Check projected tracks go through all apertures
        bool Magnet1ApertureCheck = CheckUpstreamMagnetAperture(*wcHit0, *wcHit1);
        bool Magnet2ApertureCheck = CheckDownstreamMagnetAperture(*wcHit2, *wcHit3);
        bool DSColApertureCheck   = CheckDownstreamCollimatorAperture(*wcHit2, *wcHit3);

        if (!Magnet1ApertureCheck || !Magnet2ApertureCheck || !DSColApertureCheck) continue;
        EventsPassingAperture++;

        // Candidate mass cut, keep pions, muons and electrons
        if (std::abs(TOFMass) > PI_MU_EL_MASS_CUTOFF) continue;
        EventsPassingTOFMass++;

        // If no track matched to wire-chamber, skip
        if (WC2TPCtrkID == -99999) continue;
        EventsWithWCMatch++;

        // Check WC and front-face momentum
        double WCKE             = TMath::Sqrt(WCTrackMomentum * WCTrackMomentum + PionMass * PionMass) - PionMass;
        double calculatedEnLoss = energyLossCalculation();
        if (isData) {
            double tanThetaCosPhi = TMath::Tan(WCTheta) * TMath::Cos(WCPhi);
            double tanThetaSinPhi = TMath::Tan(WCTheta) * TMath::Sin(WCPhi);
            double den            = TMath::Sqrt(1 + tanThetaCosPhi * tanThetaCosPhi);
            double onTheFlyPz     = WCTrackMomentum / den;
            double onTheFlyPx     = onTheFlyPz * tanThetaSinPhi;
            calculatedEnLoss      = energyLossCalculation(WC4PrimaryX, onTheFlyPx, isData);
        }
        const double initialKE = WCKE - calculatedEnLoss;
        hFrontFaceKE->Fill(initialKE);
        hWCKE->Fill(WCKE);

        /////////////////////////////////////////
        // All this are valid events, analyze! //
        /////////////////////////////////////////

        // Sanity check
        removeRepeatedPoints(WC2TPCLocationsX, WC2TPCLocationsY, WC2TPCLocationsZ);

        // Copy WC2TPCLocations
        std::vector<double>* wcX = new std::vector<double>(*WC2TPCLocationsX);
        std::vector<double>* wcY = new std::vector<double>(*WC2TPCLocationsY);
        std::vector<double>* wcZ = new std::vector<double>(*WC2TPCLocationsZ);

        // Get direction to end cylinder
        int numPoints = wcX->size();
        int numTail   = std::min(10, numPoints - 1);
        std::vector<std::vector<double>> points;
        for (int j = numPoints - numTail; j < numPoints; ++j) {
            points.push_back({
                wcX->at(j),
                wcY->at(j),
                wcZ->at(j)
            });
        }
        if (numTail > 0) {
            std::vector<double> avgDir = getAverageDir(points);

            // Extrapolate track to end
            double scale = (maxZ - points.back()[2]) / avgDir[2];
            wcX->push_back(points.back()[0] + scale * avgDir[0]);
            wcY->push_back(points.back()[1] + scale * avgDir[1]);
            wcZ->push_back(points.back()[2] + scale * avgDir[2]);
        }

        // First, number of non-primary TG tracks and electron cut
        int numTGTracks = 0; int numSmallTracksInCylinder = 0;
        for (size_t trk_idx = 0; trk_idx < recoBeginX->size(); ++trk_idx) {
            if (recoTrkID->at(trk_idx) == WC2TPCtrkID) continue;

            // Get track length
            double trackLength = sqrt(
                pow(recoEndX->at(trk_idx) - recoBeginX->at(trk_idx), 2) +
                pow(recoEndY->at(trk_idx) - recoBeginY->at(trk_idx), 2) +
                pow(recoEndZ->at(trk_idx) - recoBeginZ->at(trk_idx), 2)
            );

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
            if (startInCylinder && endInCylinder && (trackLength < CYLINDER_SMALL_TRACK)) numSmallTracksInCylinder++;

            if (
                !isWithinReducedVolume(recoBeginX->at(trk_idx), recoBeginY->at(trk_idx), recoBeginZ->at(trk_idx)) &&
                !isWithinReducedVolume(recoEndX->at(trk_idx), recoEndY->at(trk_idx), recoEndZ->at(trk_idx))
            ) numTGTracks++;
        }

        // Apply cuts
        if (numTGTracks > MAX_NUM_TG_TRACKS) continue;
        EventsPassingTG++;

        hSmallTracksInCylinder->Fill(numSmallTracksInCylinder);
        if (numSmallTracksInCylinder > ALLOWED_CYLINDER_SMALL_TRACKS) continue;
        EventsPassingSmall++;

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

        //////////////////////
        // Incident KE fill //
        //////////////////////

        double energyDeposited = 0.0;
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
                hIncidentKEFine->Fill(initialKE - energyDeposited);
            }
        }
        double energyAtVertex = initialKE - energyDeposited;

        ////////////////////////
        // Reduced volume cut //
        ////////////////////////

        if (!isWithinReducedVolume(breakPointX, breakPointY, breakPointZ)) continue;
        EventsInRedVol++;

        /////////////////////
        // Primary PID cut //
        /////////////////////

        if (minChi2 == pionChi2 || minChi2 == protonChi2) continue;
        EventsPassingPrimaryPID++;

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
                pow(recoBeginX->at(trk_idx) - recoEndX->at(trk_idx), 2) +
                pow(recoBeginY->at(trk_idx) - recoEndY->at(trk_idx), 2) + 
                pow(recoBeginZ->at(trk_idx) - recoEndZ->at(trk_idx), 2)
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

        //////////////////////////////////////////
        // Selection based off secondary tracks //
        //////////////////////////////////////////

        if (totalTaggedPions > 0) {
            // Reject events with > 1 pion
            if (totalTaggedPions > 1 || newSecondaryPion) continue;

            // Select as scatter
            EventsSelectedAsScatter++;
            hPionScatterKE->Fill(energyAtVertex);
            continue;
        }

        if (totalTaggedProtons > 0) {
            // Select as Np absorption
            EventsSelectedAsAbsNp++;
            hPionAbsNpKE->Fill(energyAtVertex);
            hPionAbsKE->Fill(energyAtVertex);
            continue;
        }

        EventsNoSecondary++;

        ////////////////////////////////////
        // Cluster non-reconstructed hits //
        ////////////////////////////////////

        // Get unordered set for hits in tracks
        std::unordered_set<int> hitsInTracks(hitRecoAsTrackKey->begin(), hitRecoAsTrackKey->end());

        // First, we construct the clusters 
        std::vector<int> candidateHits;
        for (size_t iHit = 0; iHit < fHitKey->size(); ++iHit) {
            double hitX     = fHitX->at(iHit);
            double hitW     = fHitW->at(iHit);
            int    hitPlane = fHitPlane->at(iHit);

            // Check if hit is near vertex of the primary
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

        // Now, get data for cut
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

        // Perform cut
        if (numClustersInduction < MAX_NUM_CLUSTERS_INDUCTION) {
            EventsSelectedAsAbs0p++;
            hPionAbs0pKE->Fill(energyAtVertex);
            hPionAbsKE->Fill(energyAtVertex);
            continue;
        }

        // We are done, all events that reach here are
        // rejected as background
    }

    // Diagnostics
    std::cout << std::endl;
    std::cout << "DIAGNOSTICS" << std::endl;
    std::cout << "==================================" << std::endl;
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout << hMCSmallTrksInCylinderMuons->Integral() + hMCSmallTrksInCylinderPions->Integral() + hMCSmallTrksInCylinderElectrons->Integral() << " events in MC" << std::endl;
    std::cout << hSmallTracksInCylinder->Integral() << " events in data" << std::endl;
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout << "All events in data: " << TotalEvents << std::endl;
    std::cout << "    Events passing wire-chamber projection cuts: " << EventsPassingProj << "(" << ((double) EventsPassingProj / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events passing aperture checks:              " << EventsPassingAperture << "(" << ((double) EventsPassingAperture / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events passing TOF mass check:               " << EventsPassingTOFMass << "(" << ((double) EventsPassingTOFMass / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events with wire-chamber to TPC match:       " << EventsWithWCMatch << "(" << ((double) EventsWithWCMatch / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events under TG tracks threshold:            " << EventsPassingTG << "(" << ((double) EventsPassingTG / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events under small tracks threshold:         " << EventsPassingSmall << "(" << ((double) EventsPassingSmall / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events with vertex in reduced volume:        " << EventsInRedVol << "(" << ((double) EventsInRedVol / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events passing primary PID:                  " << EventsPassingPrimaryPID << "(" << ((double) EventsPassingPrimaryPID / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events selected as scattering:               " << EventsSelectedAsScatter << "(" << ((double) EventsSelectedAsScatter / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events selected as absorption Np:            " << EventsSelectedAsAbsNp << "(" << ((double) EventsSelectedAsAbsNp / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events with no secondary tracks:             " << EventsNoSecondary << "(" << ((double) EventsNoSecondary / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << "    Events selected as absorption 0p:            " << EventsSelectedAsAbs0p << "(" << ((double) EventsSelectedAsAbs0p / TotalEvents) * 100. << "%)" << std::endl;
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout << "==================================" << std::endl;
    std::cout << std::endl;

    double SCALING_FACTOR     = hSmallTracksInCylinder->Integral() / (hMCSmallTrksInCylinderMuons->Integral() + hMCSmallTrksInCylinderPions->Integral() + hMCSmallTrksInCylinderElectrons->Integral());
    double SCALING_FACTOR_PRE = hFrontFaceKE->Integral() / (hMCFrontFaceKEPion->Integral() + hMCFrontFaceKEMuon->Integral() + hMCFrontFaceKEElectron->Integral());
    double SCALING_FACTOR_WC  = hWCHits->Integral() / (hMCWCHitsPion->Integral() + hMCWCHitsMuon->Integral() + hMCWCHitsElectron->Integral());

    std::vector<TH1*> scaleByNormal = {
        hMCIncidentKEPion, hMCIncidentKEElectron, hMCIncidentKEMuon,
        hMCSmallTrksInCylinderPions, hMCSmallTrksInCylinderElectrons, hMCSmallTrksInCylinderMuons,
        hMCIncidentKEPionFine, hMCIncidentKEElectronFine, hMCIncidentKEMuonFine,
        hMCPionAbs0pKETrue, hMCPionAbs0pKEAbsNp, hMCPionAbs0pKEScatter, hMCPionAbs0pKEChExch, hMCPionAbs0pKEMuon, hMCPionAbs0pKEElectron, hMCPionAbs0pKEOther,
        hMCPionAbsNpKETrue, hMCPionAbsNpKEAbs0p, hMCPionAbsNpKEScatter, hMCPionAbsNpKEChExch, hMCPionAbsNpKEMuon, hMCPionAbsNpKEElectron, hMCPionAbsNpKEOther,
        hMCPionScatterKETrue, hMCPionScatterKEAbs0p, hMCPionScatterKEAbsNp, hMCPionScatterKEChExch, hMCPionScatterKEMuon, hMCPionScatterKEElectron, hMCPionScatterKEOther
    };

    std::vector<TH1*> scaleByPre = {
        hMCFrontFaceKEPion, hMCFrontFaceKEElectron, hMCFrontFaceKEMuon,
        hMCWCKEPion, hMCWCKEMuon, hMCWCKEElectron
    };

    std::vector<TH1*> scaleByWC = {
        hMCWCHitsPion, hMCWCHitsMuon, hMCWCHitsElectron,
        hMCRadDistMidPlanePion, hMCRadDistMidPlaneMuon, hMCRadDistMidPlaneElectron,
        hMCRadDistWC4Pion, hMCRadDistWC4Muon, hMCRadDistWC4Electron
    };

    auto scaleAll = [](const std::vector<TH1*>& hists, double factor) {
        for (auto* h : hists) {
            h->Scale(factor);
        }
    };

    scaleAll(scaleByNormal, SCALING_FACTOR);
    scaleAll(scaleByPre, SCALING_FACTOR_PRE);
    scaleAll(scaleByWC, SCALING_FACTOR_WC);

    ///////////////////////////////////
    // Save histograms for unfolding //
    ///////////////////////////////////

    // Construct measured vector
    TVectorD Measure(NUM_SIGNAL_TYPES * NUM_BINS_KE);
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index      = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            Measure(index) = XSEC_UNITS * (RecoSignals[iSignal]->GetBinContent(iBin + 1) / hIncidentKE->GetBinContent(iBin + 1));
        }
    }
    TH1D* hMeasureVector = new TH1D(
        "hMeasureVector", "hMeasureVector;;", 
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE
    ); V2H(Measure, hMeasureVector);

    // Construct stat covariance matrix
    TMatrixD StatCovariance(NUM_SIGNAL_TYPES * NUM_BINS_KE, NUM_SIGNAL_TYPES * NUM_BINS_KE); StatCovariance.Zero();
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index   = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            double Ninc = hIncidentKE->GetBinContent(iBin + 1);
            double N    = Measure(index) * Ninc / XSEC_UNITS;

            double numSigma = (Ninc>0.0 && N>0.0 ? std::sqrt(N*(1.0 - N/Ninc)) : 0.0);
            double denSigma = (Ninc>0.0 ? std::sqrt(Ninc) : 0.0);

            double xs       = (N / Ninc) * XSEC_UNITS;
            double relVarXS = (N>0.0 && Ninc>0.0 ? std::pow(numSigma/N + denSigma/Ninc, 2) : 0.0);

            StatCovariance(index, index) = xs * xs * relVarXS;
        }
    }
    TH2D* hStatCovariance = new TH2D(
        "hStatCovariance", "Statistical Covariance Matrix;(i, #alpha);(j, #beta)",
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE, 
        NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE 
    ); M2H(StatCovariance, hStatCovariance);

    TFile* dataFile = new TFile("/exp/lariat/app/users/epelaez/histos/data/Reco.root", "RECREATE");

    // Construct measure and stat cov for abs + scatt measurement
    TVectorD MeasureAbsScatt(2 * NUM_BINS_KE);
    for (int iSignal = 0; iSignal < 2; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            MeasureAbsScatt(index) = XSEC_UNITS * (RecoSignalsAbsScatt[iSignal]->GetBinContent(iBin + 1) / hIncidentKE->GetBinContent(iBin + 1));
        }
    }
    TH1D* hMeasureVectorAbsScatt = new TH1D(
        "hMeasureVectorAbsScatt", "hMeasureVectorAbsScatt;;", 
        2 * NUM_BINS_KE, 0, 2 * NUM_BINS_KE
    ); V2H(MeasureAbsScatt, hMeasureVectorAbsScatt);

    TMatrixD StatCovarianceAbsScatt(2 * NUM_BINS_KE, 2 * NUM_BINS_KE); StatCovarianceAbsScatt.Zero();
    for (int iSignal = 0; iSignal < 2; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index   = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            double Ninc = hIncidentKE->GetBinContent(iBin + 1);
            double N    = MeasureAbsScatt(index) * Ninc / XSEC_UNITS;

            double numSigma = (Ninc>0.0 && N>0.0 ? std::sqrt(N*(1.0 - N/Ninc)) : 0.0);
            double denSigma = (Ninc>0.0 ? std::sqrt(Ninc) : 0.0);

            double xs       = (N / Ninc) * XSEC_UNITS;
            double relVarXS = (N>0.0 && Ninc>0.0 ? std::pow(numSigma/N + denSigma/Ninc, 2) : 0.0);

            StatCovarianceAbsScatt(index, index) = xs * xs * relVarXS;
        }
    }
    TH2D* hStatCovarianceAbsScatt = new TH2D(
        "hStatCovarianceAbsScatt", "Statistical Covariance Matrix;(i, #alpha);(j, #beta)",
        2 * NUM_BINS_KE, 0, 2 * NUM_BINS_KE, 
        2 * NUM_BINS_KE, 0, 2 * NUM_BINS_KE 
    ); M2H(StatCovarianceAbsScatt, hStatCovarianceAbsScatt);

    dataFile->cd();
    hMeasureVector->Write("", TObject::kOverwrite);
    hStatCovariance->Write("", TObject::kOverwrite);

    hMeasureVectorAbsScatt->Write("", TObject::kOverwrite);
    hStatCovarianceAbsScatt->Write("", TObject::kOverwrite);
    dataFile->Close();

    ////////////////////////////
    // Get raw cross-sections //
    ////////////////////////////

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

    std::vector<TH1*> PlotDataGroups = {
        // WC quality
        hWCHits,
        hRadDistWC4,
        hRadDistMidPlane,

        // Cylinder
        hSmallTracksInCylinder,

        // Incident KE
        hIncidentKE,
        hIncidentKEFine,
        hFrontFaceKE,
        hWCKE,

        // Interacting KE
        hPionAbs0pKE,
        hPionAbsNpKE,
        hPionScatterKE
    };

    std::vector<std::vector<TH1*>> PlotMCGroups = {
        // WC quality
        {hMCWCHitsPion, hMCWCHitsMuon, hMCWCHitsElectron},
        {hMCRadDistWC4Pion, hMCRadDistWC4Muon, hMCRadDistWC4Electron},
        {hMCRadDistMidPlanePion, hMCRadDistMidPlaneMuon, hMCRadDistMidPlaneElectron},
        
        // Cylinder
        {hMCSmallTrksInCylinderPions, hMCSmallTrksInCylinderMuons, hMCSmallTrksInCylinderElectrons},

        // Incident KE
        {hMCIncidentKEPion, hMCIncidentKEMuon, hMCIncidentKEElectron},
        {hMCIncidentKEPionFine, hMCIncidentKEMuonFine, hMCIncidentKEElectronFine},
        {hMCFrontFaceKEPion, hMCFrontFaceKEMuon, hMCFrontFaceKEElectron},
        {hMCWCKEPion, hMCWCKEMuon, hMCWCKEElectron},

        // Interacting KE
        {hMCPionAbs0pKETrue, hMCPionAbs0pKEAbsNp, hMCPionAbs0pKEScatter, hMCPionAbs0pKEChExch, hMCPionAbs0pKEMuon, hMCPionAbs0pKEElectron, hMCPionAbs0pKEOther},
        {hMCPionAbsNpKETrue, hMCPionAbsNpKEAbs0p, hMCPionAbsNpKEScatter, hMCPionAbsNpKEChExch, hMCPionAbsNpKEMuon, hMCPionAbsNpKEElectron, hMCPionAbsNpKEOther},
        {hMCPionScatterKETrue, hMCPionScatterKEAbs0p, hMCPionScatterKEAbsNp, hMCPionScatterKEChExch, hMCPionScatterKEMuon, hMCPionScatterKEElectron, hMCPionScatterKEOther}
    };

    std::vector<std::vector<TString>> PlotMCLabelGroups = {
        // WC quality
        {"Pion", "Muon", "Electron"},
        {"Pion", "Muon", "Electron"},
        {"Pion", "Muon", "Electron"},
        
        // Cylinder
        {"Pion", "Muon", "Electron"},
        
        // Incident KE
        {"Pion", "Muon", "Electron"},
        {"Pion", "Muon", "Electron"},
        {"Pion", "Muon", "Electron"},
        {"Pion", "Muon", "Electron"},

        // Interacting KE
        {"True", "Abs Np", "Scatter", "Ch. Exch.", "Muon", "Electron", "Other"},
        {"True", "Abs 0p", "Scatter", "Ch. Exch.", "Muon", "Electron", "Other"},
        {"True", "Abs 0p", "Abs Np", "Ch. Exch.", "Muon", "Electron", "Other"}
    };

    std::vector<TString> PlotName = {
        // WC quality
        "WCQuality/NumHits",
        "WCQuality/RadDistWC4",
        "WCQuality/RadDistMidPlane",

        // Cylinder 
        "Cylinder/SmallTracks",

        // Incident KE
        "Incident/IncidentKE",
        "Incident/IncidentKEFine",
        "Incident/FrontFaceKE",
        "Incident/WireChamberKE",

        // Interacting KE
        "Interacting/PionAbs0pKE",
        "Interacting/PionAbsNpKE",
        "Interacting/PionScatterKE"
    };

    std::vector<TString> PlotTitle = {
        // WC Quality
        "Number of Wire-Chamber Hits",
        "Proj. to WC hit",
        "Up to downstream proj. distance",

        // Cylinder
        "Small Tracks in Cylinder",

        // Incident KE
        "Incident Kinetic Energy",
        "Incident Kinetic Energy (Fine Binning)",
        "Front-Face Kinetic Energy",
        "Wire-Chamber Kinetic Energy",

        // Interacting KE
        "Abs 0p Kinetic Energy",
        "Abs Np Kinetic Energy",
        "Scatter Kinetic Energy"
    };

    std::vector<TString> XLabels = {
        // WC Quality
        "# of hits",
        "Proj. to WC hit distance [cm]",
        "Upstream to downstream proj. distance [cm]",

        // Cylinder
        "# of small tracks in cylinder",

        // Incident KE
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",

        // Interacting KE
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]",
        "Kinetic energy [MeV]"
    };

    std::vector<TString> YLabels = {
        // WC Quality
        "Counts",
        "Counts",
        "Counts",

        // Cylinder
        "Counts",

        // Incident KE
        "Counts",
        "Counts",
        "Counts",
        "Counts",

        // Interacting KE
        "Counts",
        "Counts",
        "Counts"
    };

    for (size_t i = 0; i < PlotDataGroups.size(); ++i) {
        PrintDataVsMCContribPlot(
            SaveDir,
            PlotName[i],
            PlotDataGroups[i],
            PlotMCGroups[i],
            PlotMCLabelGroups[i],
            Colors,
            PlotTitle[i],
            XLabels[i],
            YLabels[i],
            FontStyle, TextSize,
            true,
            nullptr,
            true
        );
    }
}