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
#include "FakeDataHelper.h"

void FakeData(int sample = 1) {
    // Set defaults
    gStyle->SetOptStat(0); // get rid of stats box
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    TH1::AddDirectory(false);
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
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis/figs/FakeData/Sample" + TString::Itoa(sample, 10) + "/";
    std::cout << "Generating Sample: " << sample << std::endl;

    // Load file with data products
    TString RootFilePath = "/exp/lariat/app/users/epelaez/files/RecoAll_histo.root"; // RV at z = 30
    std::unique_ptr<TFile> File(TFile::Open(RootFilePath));
    TDirectory* Directory = (TDirectory*)File->Get("RecoNNAllEval");

    // Fake data configuration
    FakeDataFD::Scenario FDScenario;
    if (sample == 1) {
        FDScenario = FakeDataFD::Scenario::FDa;
    } else if (sample == 2) {
        FDScenario = FakeDataFD::Scenario::FDb;
    } else if (sample == 3) {
        FDScenario = FakeDataFD::Scenario::FDc;
    } else if (sample == 4) {
        FDScenario = FakeDataFD::Scenario::FDd;
    }

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

    // Hits reconstructed into full tracks
    std::vector<std::vector<int>>* recoTrackHitIndices  = nullptr;
    std::vector<std::vector<double>>*     recoTrackHitX = nullptr; 
    std::vector<std::vector<double>>*     recoTrackHitY = nullptr; 
    std::vector<std::vector<double>>*     recoTrackHitZ = nullptr; 
    tree->SetBranchAddress("recoTrackHitIndices", &recoTrackHitIndices);
    tree->SetBranchAddress("recoTrackHitX", &recoTrackHitX);
    tree->SetBranchAddress("recoTrackHitY", &recoTrackHitY);
    tree->SetBranchAddress("recoTrackHitZ", &recoTrackHitZ);

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

    // Histograms with classified events
    TH1D* hPionAbs0p   = new TH1D("hPionAbs0p", "hPionAbs0p;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionAbsNp   = new TH1D("hPionAbsNp", "hPionAbsNp;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hPionScatter = new TH1D("hPionScatter", "hPionScatter;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    std::vector<TH1*> RecoSignals = {hPionAbs0p, hPionAbsNp, hPionScatter};

    // Reconstructed incident energy
    TH1D* hPionIncidentKE = new TH1D("hPionIncidentKE", "hPionIncidentKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    
    // True interacting events
    TH1D* hTruePionAbs0p   = new TH1D("hTruePionAbs0p", "hTruePionAbs0p;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTruePionAbsNp   = new TH1D("hTruePionAbsNp", "hTruePionAbsNp;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTruePionScatter = new TH1D("hTruePionScatter", "hTruePionScatter;;", NUM_BINS_KE, ARRAY_KE_BINS.data());
    std::vector<TH1*> TotalEventsHistos = {hTruePionAbs0p, hTruePionAbsNp, hTruePionScatter};

    // True incident energy
    TH1D* hTrueIncidentKE = new TH1D("hTrueIncidentKE", "hTrueIncidentKE;;", NUM_BINS_KE, ARRAY_KE_BINS.data());

    //////////////////////
    // Loop over events //
    //////////////////////

    Int_t NumEntries = (Int_t) tree->GetEntries();
    std::cout << "Num entries: " << NumEntries << std::endl;

    for (Int_t i = 0; i < NumEntries; ++i) {
        tree->GetEntry(i);

        // Make script go faster
        // if (i > USE_NUM_EVENTS) break;

        // Get unordered set for hits in tracks
        std::unordered_set<int> hitsInTracks(hitRecoAsTrackKey->begin(), hitRecoAsTrackKey->end());

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

        //////////////////////////////////
        // Look at unreconstructed hits //
        //////////////////////////////////

        // Units for coordinates:
        //     W: (channel) * (channel width) = cm
        //     X: (drift velocity) * (time) = (cm / us) * us = cm

        // First, find hits near the primary track that were not
        // reconstructed into any of the already existing tracks
        std::vector<int> candidateHits;

        int numTotalHitsNearPrimary            = 0;
        int numUnRecoHitsNearPrimaryInduction  = 0;
        int numUnRecoHitsNearPrimaryCollection = 0;
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
                numTotalHitsNearPrimary++;

                // Skip hits already in tracks
                if (hitsInTracks.count(iHit) > 0) continue;

                if (hitPlane == 0) numUnRecoHitsNearPrimaryInduction++;
                else if (hitPlane == 1) numUnRecoHitsNearPrimaryCollection++;

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
                // if (usedHits.count(iAllHit) || hitsInTracks.count(iAllHit) || (fHitPlane->at(iAllHit) == 1)) continue;
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


        ////////////////////////////
        // Save truth information //
        ////////////////////////////

        // Scattering only if degree > THRESHOLD and energy > THRESHOLD
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
                if (backgroundType == 12) truthPrimaryVertexKE = trajectoryInteractionKE;
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
                    backgroundType = 0;
                }

            }
        }

        ///////////////////////////////
        // Get weight for this event //
        ///////////////////////////////

        double fd_weight= FakeDataFD::GetFDWeight(FDScenario, backgroundType, truthPrimaryVertexKE * 1000);
        // if (fd_weight != 1) {
        //     std::cout << "Background type: " << backgroundType << std::endl;
        //     std::cout << "Truth primary vertex KE: " << truthPrimaryVertexKE * 1000 << " MeV" << std::endl;
        //     std::cout << "Fake data weight applied: " << fd_weight << std::endl;
        //     std::cout << std::endl;
        // }

        ///////////////////
        // True spectrum //
        ///////////////////

        // Get true energy bin
        int TrueEnergyBin = getBin(truthPrimaryVertexKE * 1000, ARRAY_KE_BINS);

        // Add true incident KE
        if (validTrueIncidentKE) {
            for (double x : *trueIncidentKEContributions) {
                hTrueIncidentKE->Fill(x);
            }
        }

        // True interaction histograms
        if (backgroundType == 0) {
            hTruePionAbs0p->Fill(truthPrimaryVertexKE * 1000, fd_weight);
        } else if (backgroundType == 1) {
            hTruePionAbsNp->Fill(truthPrimaryVertexKE * 1000, fd_weight);
        } else if (backgroundType == 6 || backgroundType == 12) {
            hTruePionScatter->Fill(truthPrimaryVertexKE * 1000, fd_weight);
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
                hPionIncidentKE->Fill(initialKE - energyDeposited);
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

        // Fill scattering
        if (totalTaggedPions > 0) {
            if (totalTaggedPions > 1 || newSecondaryPion) continue;

            hPionScatter->Fill(energyAtVertex, fd_weight);
            continue;
        }

        // Fill abs Np
        if (totalTaggedProtons > 0) {
            hPionAbsNp->Fill(energyAtVertex, fd_weight);
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

        // Fill abs 0p
        if (numClustersInduction < MAX_NUM_CLUSTERS_INDUCTION) {
            hPionAbs0p->Fill(energyAtVertex, fd_weight);
            continue;
        }
    }

    ////////////////////////
    // True cross-section //
    ////////////////////////

    TH1D* hSignal = new TH1D("hSignal", "hSignal;;", NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index   = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            double xsec =  XSEC_UNITS * (TotalEventsHistos[iSignal]->GetBinContent(iBin + 1) / hTrueIncidentKE->GetBinContent(iBin + 1));
            hSignal->SetBinContent(index + 1,  xsec);
        }
    }

    ////////////////////////
    // Reco cross-section //
    ////////////////////////

    TH1D* hMeasure = new TH1D("hMeasure", "hMeasure;;", NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE);
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        for (int iBin = 0; iBin < NUM_BINS_KE; ++iBin) {
            int index = flattenIndex(iSignal, iBin, NUM_BINS_KE);
            double xsec = XSEC_UNITS * (RecoSignals[iSignal]->GetBinContent(iBin + 1) / hPionIncidentKE->GetBinContent(iBin + 1));
            hMeasure->SetBinContent(index + 1, xsec);
        }
    }

    //////////////////
    // Save to file //
    //////////////////

    TFile* saveFile = new TFile(("/exp/lariat/app/users/epelaez/histos/fake_data/FD" + TString::Itoa(sample, 10) + ".root").Data(), "RECREATE");
    saveFile->cd();
    hSignal->Write("", TObject::kOverwrite);
    hMeasure->Write("", TObject::kOverwrite);
    hPionIncidentKE->Write("", TObject::kOverwrite);
    saveFile->Close();

    ///////////////////////
    // Create histograms //
    ///////////////////////

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
        {hMeasure, hSignal}
    };

    std::vector<std::vector<TString>> PlotLabelGroups = {
        {"Reco", "True"}
    };

    std::vector<TString> PlotTitles = {
        "Sample"
    };

    std::vector<TString> XLabels = {
        "Bin"
    };

    std::vector<TString> YLabels = {
        "Cross-Section [barn]"
    };

    std::vector<bool> PlotStacked = {
        false
    };

    std::vector<std::vector<bool>> PlotsAsPoints = {
        {true, false}
    };

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
}