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

void RecoDataAnalysis() {
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
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis/figs/RecoDataAnalysis/";

    // Load file with NN data products
    TString RootFilePath = "/exp/lariat/app/users/epelaez/files/RecoNNDataEval_histo.root";
    std::unique_ptr<TFile> File(TFile::Open(RootFilePath));
    TDirectory* Directory = (TDirectory*)File->Get("RecoNNDataEval");

    // Load file with MC comparison histograms
    TString MCCompFilePath = "/exp/lariat/app/users/epelaez/files/DataMCComparisons.root";
    std::unique_ptr<TFile> MCFile(TFile::Open(MCCompFilePath));

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

    // TOF information
    double TOFMass, tofObject;
    tree->SetBranchAddress("TOFMass", &TOFMass);
    tree->SetBranchAddress("tofObject", &tofObject);

    /////////////////////
    // Load histograms //
    /////////////////////

    TH1D* hMCNumWC2TPCMatch             = (TH1D*) MCFile->Get("hMCNumWC2TPCMatch");
    TH1D* hMCNumTracksInCylinder        = (TH1D*) MCFile->Get("hMCNumTracksInCylinder");
    TH1D* hMCTGTrackLengths             = (TH1D*) MCFile->Get("hMCTGTrackLengths");
    TH1D* hMCTGSmallTracks              = (TH1D*) MCFile->Get("hMCTGSmallTracks");
    TH1D* hMCTracksNearVertex           = (TH1D*) MCFile->Get("hMCTracksNearVertex");
    TH1D* hMCTrackLengthsNearVertex     = (TH1D*) MCFile->Get("hMCTrackLengthsNearVertex");
    TH1D* hMCNumTGTracks                = (TH1D*) MCFile->Get("hMCNumTGTracks");
    TH1D* hMCShowerProb                 = (TH1D*) MCFile->Get("hMCShowerProb");
    TH1D* hMCBeforeShowerCutSmallTracks = (TH1D*) MCFile->Get("hMCBeforeShowerCutSmallTracks");
    TH1D* hMCAfterShowerCutSmallTracks  = (TH1D*) MCFile->Get("hMCAfterShowerCutSmallTracks");
    TH2D* hMCSmallVsTGTracks            = (TH2D*) MCFile->Get("hMCSmallVsTGTracks");
    TH2D* hMCTGNumSmallTracksVsThresh   = (TH2D*) MCFile->Get("hMCTGNumSmallTracksVsThresh");
    TH2D* hMCPrimaryTrackPosition       = (TH2D*) MCFile->Get("hMCPrimaryTrackPosition");

    ///////////////////////
    // Create histograms //
    ///////////////////////

    // TOF
    TH1D* hTOFMass = new TH1D("hTOFMass", "TOF Mass Distribution", 50, 0, 1200);
    TH1D* hTOF     = new TH1D("hTOF", "TOF Distribution", 50, 0, 80);
    
    // WC
    TH1D* hNumWC2TPCMatch = new TH1D("hNumWC2TPCMatch", "NumWC2TPCMatch", 10, 0, 10);
    TH1D* hNumWCHits      = new TH1D("hNumWCHits", "hNumWCHits", 10, 0, 10);

    // Shower cut
    TH1D* hShowerProb                 = new TH1D("hShowerProb", "ShowerProb", 20, 0, 1.);
    TH1D* hBeforeShowerCutSmallTracks = new TH1D("hBeforeShowerCutSmallTracks", "BeforeShowerCutSmallTracks", 10, 0, 10);
    TH1D* hAfterShowerCutSmallTracks  = new TH1D("hAfterShowerCutSmallTracks", "AfterShowerCutSmallTracks", 10, 0, 10);

    // With primary TG
    TH1D* hNumTracksInCylinder = new TH1D("hNumTracksInCylinder", "NumTracksInCylinder", 10, 0, 10);
    TH1D* hTGTrackLengths      = new TH1D("hTGTrackLengths", "TGTrackLengths", 25, 0, 50);
    TH1D* hTGSmallTracks       = new TH1D("hTGSmallTracks", "TGSmallTracks", 10, 0, 10);

    TH1D* hNumTracksInCylinder0TG = new TH1D("hNumTracksInCylinder0TG", "NumTracksInCylinder0TG", 10, 0, 10);
    TH1D* hTGSmallTracks0TG       = new TH1D("hTGSmallTracks0TG", "TGSmallTracks0TG", 10, 0, 10);

    TH1D* hNumTracksInCylinder1TG = new TH1D("hNumTracksInCylinder1TG", "NumTracksInCylinder1TG", 10, 0, 10);
    TH1D* hTGSmallTracks1TG       = new TH1D("hTGSmallTracks1TG", "TGSmallTracks1TG", 10, 0, 10);

    TH1D* hNumTracksInCylinder2TG = new TH1D("hNumTracksInCylinder2TG", "NumTracksInCylinder2TG", 10, 0, 10);
    TH1D* hTGSmallTracks2TG       = new TH1D("hTGSmallTracks2TG", "TGSmallTracks2TG", 10, 0, 10);

    // With primary interacting
    TH1D* hTracksNearVertex       = new TH1D("hTracksNearVertex", "TracksNearVertex", 10, 0, 10);
    TH1D* hTrackLengthsNearVertex = new TH1D("hTrackLengthsNearVertex", "TrackLengthsNearVertex", 50, 0, 100);

    // TG and small tracks
    TH1D* hNumTGTracks              = new TH1D("hNumTGTracks", "NumTGTracks", 10, 0, 10);
    TH2D* hSmallVsTGTracks          = new TH2D("hSmallVsTGTracks", "SmallVsTGTracks;Small Tracks;TG Tracks", 15, 0, 15, 15, 0, 15);
    TH2D* hTGNumSmallTracksVsThresh = new TH2D("hTGNumSmallTracksVsThresh", "TGNumSmallTracksVsThresh;Small Track Length Threshold (cm);Num Small Tracks", 10, 0, 40, 15, 0, 15);
    TH2D* hTOFVsTOFMass             = new TH2D("hTOFVsTOFMass", "TOFVsTOFMass;TOF [ns];TOF Mass [MeV/c^2]", 35, 10, 80, 50, 0, 1200);

    ///////////////////////////////////
    // Distribution of primary track //
    ///////////////////////////////////

    TH2D* hPrimaryTrackPosition = new TH2D(
        "hPrimaryTrackPosition", "PrimaryTrackPosition;x-position;y-position",
        20, minX, maxX,
        20, minY, maxY
    );

    ///////////////////////////////////////
    // Distribution of background tracks //
    ///////////////////////////////////////

    TH2D* hBackgroundTracksPosition = new TH2D(
        "hBackgroundTracksPosition", "BackgroundTracksPosition;x-position;y-position",
        20, minX, maxX,
        20, minY, maxY
    );

    TH2D* hBackgroundTracksDirection = new TH2D(
        "hBackgroundTracksDirection", "BackgroundTracksDirection;#phi (x-y) [rad];#theta (forward) [rad]",
        20, -TMath::Pi(), TMath::Pi(),
        20, 0, TMath::Pi() / 2
    );

    /////////////////////////////
    // Distribution of WC hits //
    /////////////////////////////

    TH2D* hProjVsRealWC4 = new TH2D(
        "hProjVsRealWC4", "ProjVsRealWC4;p_{x} - WC_{x};p_{y} - WC_{y}",
        100, -40., 40.,
        100, -15., 15.
    );
    TH1D* hRadDistWC4 = new TH1D("hRadDistWC4", "RadDistWC4", 80, 0.0, 40.);

    TH2D* hMidPlaneCoinc = new TH2D(
        "hMidPlaneCoinc", "MidPlaneCoinc;up_{x} - down_{x};up_{y} - down_{y}",
        100, -20., 20.,
        100, -5., 5.
    );
    TH1D* hRadDistMidPlane = new TH1D("hRadDistMidPlane", "RadDistMidPlane", 60, 0.0, 30.);

    //////////////////////
    // Loop over events //
    //////////////////////

    Int_t NumEntries = (Int_t) tree->GetEntries();
    std::cout << "Num entries: " << NumEntries << std::endl;

    int numValidEvents = 0;
    for (Int_t i = 0; i < NumEntries; ++i) {
        tree->GetEntry(i);

        // Fill for all events
        hTOFMass->Fill(std::abs(TOFMass));
        hNumWC2TPCMatch->Fill(WC2TPCsize);
        hNumWCHits->Fill(wcNumHits);
        hTOF->Fill(tofObject);
        hTOFVsTOFMass->Fill(tofObject, std::abs(TOFMass));

        /////////////////////////////////////////
        // Wire-chamber and other initial cuts //
        /////////////////////////////////////////

        if (wcNumHits != 8) continue;

        // Project downstream to midplane @ -437.97 (without angular corrections)
        std::vector<double> midUp = projToZ(*wcHit0, *wcHit1, -437.97);
        // Use this point and WC3 to project up to WC4
        std::vector<double> projDown = projToZ(midUp, *wcHit2, -95.0);
        // Requires some corrections because magnets are not the same
        projDown[0] -= tan(1.32 * TMath::Pi() / 180.0) * (-95.0 - -437.97);

        // Compare x and y coordinate in projection and real hit for WC4
        hProjVsRealWC4->Fill(projDown[0] - wcHit3->at(0), projDown[1] - wcHit3->at(1));
        double radDistWC4 = TMath::Sqrt(pow(projDown[0] - wcHit3->at(0), 2.) + pow(projDown[1] - wcHit3->at(1), 2.));
        hRadDistWC4->Fill(radDistWC4);

        // if (radDistWC4 > 8.0) continue; // TODO: figure out number

        // Project upstream to midplane @ -437.97
        std::vector<double> midDown = projToZ(*wcHit3, *wcHit2, -437.97);
        midDown[0] -= tan(1.32 * TMath::Pi() / 180.0) * (-339.57 - -437.97);

        hMidPlaneCoinc->Fill(midUp[0] - midDown[0], midUp[1] - midDown[1]);
        double midPlaneDist = TMath::Sqrt(pow(midUp[0] - midDown[0], 2) + pow(midUp[1] - midDown[1], 2));
        hRadDistMidPlane->Fill(midPlaneDist);

        // if (midPlaneDist > 3.0) continue;

        // Check projected tracks go through all apertures
        bool Magnet1ApertureCheck = CheckUpstreamMagnetAperture(*wcHit0, *wcHit1);
        bool Magnet2ApertureCheck = CheckDownstreamMagnetAperture(*wcHit3, *wcHit2);
        bool DSColApertureCheck   = CheckDownstreamCollimatorAperture(*wcHit3, *wcHit2);

        if (!Magnet1ApertureCheck || !Magnet2ApertureCheck || !DSColApertureCheck) continue;

        // Candidate mass cut, keep pions, muons and electrons
        if (std::abs(TOFMass) > PI_MU_EL_MASS_CUTOFF) continue;

        // If no track matched to wire-chamber, skip
        if (WC2TPCtrkID == -99999) continue;

        //////////////////////////
        // Analyze valid events //
        //////////////////////////

        numValidEvents++;
        hPrimaryTrackPosition->Fill(WC2TPCPrimaryBeginX, WC2TPCPrimaryBeginY);

        // Check if WC2TPC is through-going
        bool isPrimaryTG = !isWithinReducedVolume(WC2TPCPrimaryEndX, WC2TPCPrimaryEndY, WC2TPCPrimaryEndZ);

        // Tracks in first 14 cm of upstream TPC
        int countEarlyTracks = 0;
        for (size_t trk_idx = 0; trk_idx < recoTrkID->size(); ++trk_idx) {
            if (recoTrkID->at(trk_idx) == WC2TPCtrkID) continue;

            if (
                recoBeginZ->at(trk_idx) < 14.0 ||
                recoEndZ->at(trk_idx) < 14.0
            ) {
                countEarlyTracks++;
            }
        }
        // if (countEarlyTracks > 4) continue;

        // If high shower probability, reject
        // if (showerProb >= SHOWER_PROB_CUT) continue;
        if (obtainedProbabilities) hShowerProb->Fill(showerProb);

        // Loop over reconstructed tracks
        int numSmallTracks = 0; 
        int numTracksNearVertex = 0; 
        int smallTracksTPCStart = 0; 
        int numTGTracks = 0;
        int numTracksInCylinder = 0;

        for (size_t trk_idx = 0; trk_idx < recoBeginX->size(); ++trk_idx) {
            if (recoTrkID->at(trk_idx) == WC2TPCtrkID) continue;

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
            if (isPrimaryTG) hTGTrackLengths->Fill(trackLength);

            if (
                recoEndZ->at(trk_idx) < 30.0 && 
                recoBeginZ->at(trk_idx) < 30.0
            ) {
                if (trackLength < SMALL_TRACK_LENGTH_CHEX) smallTracksTPCStart++;
            }

            if (
                !isWithinReducedVolume(recoBeginX->at(trk_idx), recoBeginY->at(trk_idx), recoBeginZ->at(trk_idx)) &&
                !isWithinReducedVolume(recoEndX->at(trk_idx), recoEndY->at(trk_idx), recoEndZ->at(trk_idx))
            ) {
                numTGTracks++;

                if (isPrimaryTG) {
                    // Track information about background tracks
                    bool isBeginTrue = recoBeginZ->at(trk_idx) < recoEndZ->at(trk_idx);

                    const double xs = isBeginTrue ? recoBeginX->at(trk_idx) : recoEndX->at(trk_idx);
                    const double ys = isBeginTrue ? recoBeginY->at(trk_idx) : recoEndY->at(trk_idx);
                    const double zs = isBeginTrue ? recoBeginZ->at(trk_idx) : recoEndZ->at(trk_idx);

                    const double xe = isBeginTrue ? recoEndX->at(trk_idx)   : recoBeginX->at(trk_idx);
                    const double ye = isBeginTrue ? recoEndY->at(trk_idx)   : recoBeginY->at(trk_idx);
                    const double ze = isBeginTrue ? recoEndZ->at(trk_idx)   : recoBeginZ->at(trk_idx);

                    hBackgroundTracksPosition->Fill(xs, ys);

                    auto [phi, theta] = azimuth_polar_from_points(
                        xs, ys, zs,
                        xe, ye, ze
                    );
                    hBackgroundTracksDirection->Fill(phi, theta);
                }
            }

            if (
                !isPrimaryTG &&
                (distanceFromStart < VERTEX_RADIUS || distanceFromEnd < VERTEX_RADIUS)
            ) {
                numTracksNearVertex++;
                hTrackLengthsNearVertex->Fill(trackLength);
            }
            if (trackLength < SMALL_TRACK_LENGTH_CHEX) numSmallTracks++;
        }

        // Look at relevant quantities with different number of TG tracks
        if (numTGTracks == 0 && isPrimaryTG) {
            hNumTracksInCylinder0TG->Fill(numTracksInCylinder);
            hTGSmallTracks0TG->Fill(numSmallTracks);
        } 
        if (numTGTracks <= 1 && isPrimaryTG) {
            hNumTracksInCylinder1TG->Fill(numTracksInCylinder);
            hTGSmallTracks1TG->Fill(numSmallTracks);
        }
        if (numTGTracks <= 2 && isPrimaryTG) {
            hNumTracksInCylinder2TG->Fill(numTracksInCylinder);
            hTGSmallTracks2TG->Fill(numSmallTracks);
        }

        hBeforeShowerCutSmallTracks->Fill(smallTracksTPCStart);
        hNumTGTracks->Fill(numTGTracks);

        if (!isPrimaryTG) hTracksNearVertex->Fill(numTracksNearVertex);
        if (obtainedProbabilities && showerProb < SHOWER_PROB_CUT) hAfterShowerCutSmallTracks->Fill(smallTracksTPCStart);

        // Add to histogram of small tracks if primary is throughgoing
        if (isPrimaryTG) {
            hTGSmallTracks->Fill(numSmallTracks);
            hSmallVsTGTracks->Fill(numSmallTracks, numTGTracks);
            hNumTracksInCylinder->Fill(numTracksInCylinder);

            // Scan over small track length thresholds and fill 2D histogram
            for (int threshBin = 1; threshBin <= hTGNumSmallTracksVsThresh->GetNbinsX(); ++threshBin) {
                double threshold = hTGNumSmallTracksVsThresh->GetXaxis()->GetBinCenter(threshBin);
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
                hTGNumSmallTracksVsThresh->Fill(threshold, nSmallTracks);
            }
        }
    }
    
    double numMCEvents   = hMCNumTGTracks->Integral(0, hMCNumTGTracks->GetNbinsX() + 1);
    double numDataEvents = hNumTGTracks->Integral(0, hNumTGTracks->GetNbinsX() + 1);
    double scaling       = numDataEvents / numMCEvents;
    std::cout << "Num valid MC events: " << numMCEvents << std::endl;
    std::cout << "Num valid data events: " << numDataEvents << std::endl;
    std::cout << "Scaling: " << scaling << std::endl;
    std::cout << std::endl;

    double numMCTGEvents   = hMCTGSmallTracks->Integral(0, hMCTGSmallTracks->GetNbinsX() + 1);
    double numDataTGEvents = hTGSmallTracks->Integral(0, hTGSmallTracks->GetNbinsX() + 1);
    double scalingTG       = numDataTGEvents / numMCTGEvents;
    std::cout << "Num valid MC TG events: " << numMCTGEvents << std::endl;
    std::cout << "Num valid data TG events: " << numDataTGEvents << std::endl;
    std::cout << "Scaling TG: " << scalingTG << std::endl;
    std::cout << std::endl;

    double numMCNotTGEvents = hMCTracksNearVertex->Integral(0, hMCTracksNearVertex->GetNbinsX() + 1);
    double numDataNoTGvents = hTracksNearVertex->Integral(0, hTracksNearVertex->GetNbinsX() + 1);
    double scalingNoTG      = numDataNoTGvents / numMCNotTGEvents;
    std::cout << "Num valid MC not TG events: " << numMCNotTGEvents << std::endl;
    std::cout << "Num valid data not TG events: " << numDataNoTGvents << std::endl;
    std::cout << "Scaling not TG: " << scalingNoTG << std::endl;
    std::cout << std::endl;

    std::cout << "Num valid data TG events with no other TG: " << hNumTracksInCylinder0TG->Integral() << std::endl;
    std::cout << "Num valid data TG events with at most one more TG: " << hNumTracksInCylinder1TG->Integral() << std::endl;
    std::cout << "Num valid data TG events with at most two more TG: " << hNumTracksInCylinder2TG->Integral() << std::endl;

    // Scale MC histograms to data histograms using event counts
    hMCNumTGTracks->Scale(scaling);
    hMCBeforeShowerCutSmallTracks->Scale(scaling);
    hMCAfterShowerCutSmallTracks->Scale(hAfterShowerCutSmallTracks->Integral() / hMCAfterShowerCutSmallTracks->Integral()); // not sure

    hMCTGSmallTracks->Scale(scalingTG);
    hMCTGTrackLengths->Scale(scalingTG);
    hMCNumTracksInCylinder->Scale(scalingTG);
    
    hMCTracksNearVertex->Scale(scalingNoTG);
    hMCTrackLengthsNearVertex->Scale(scalingNoTG);

    // Have to account only for events with valid NN probability
    hMCShowerProb->Scale(hShowerProb->Integral() / hMCShowerProb->Integral());
    hMCNumWC2TPCMatch->Scale(hNumWC2TPCMatch->Integral() / hMCNumWC2TPCMatch->Integral());

    // Have to account for cutting off more events
    TH1D* hMCTGSmallTracks0TG = (TH1D*) hMCTGSmallTracks->Clone("hMCTGSmallTracks0TG");
    TH1D* hMCTGSmallTracks1TG = (TH1D*) hMCTGSmallTracks->Clone("hMCTGSmallTracks1TG");
    TH1D* hMCTGSmallTracks2TG = (TH1D*) hMCTGSmallTracks->Clone("hMCTGSmallTracks2TG");

    hMCTGSmallTracks0TG->Scale(hTGSmallTracks0TG->Integral() / hMCTGSmallTracks0TG->Integral());
    hMCTGSmallTracks1TG->Scale(hTGSmallTracks1TG->Integral() / hMCTGSmallTracks1TG->Integral());
    hMCTGSmallTracks2TG->Scale(hTGSmallTracks2TG->Integral() / hMCTGSmallTracks2TG->Integral());

    TH1D* hMCNumTracksInCylinder0TG = (TH1D*) hMCNumTracksInCylinder->Clone("hMCNumTracksInCylinder0TG");
    TH1D* hMCNumTracksInCylinder1TG = (TH1D*) hMCNumTracksInCylinder->Clone("hMCNumTracksInCylinder1TG");
    TH1D* hMCNumTracksInCylinder2TG = (TH1D*) hMCNumTracksInCylinder->Clone("hMCNumTracksInCylinder2TG");

    hMCNumTracksInCylinder0TG->Scale(hNumTracksInCylinder0TG->Integral() / hMCNumTracksInCylinder0TG->Integral());
    hMCNumTracksInCylinder1TG->Scale(hNumTracksInCylinder1TG->Integral() / hMCNumTracksInCylinder1TG->Integral());
    hMCNumTracksInCylinder2TG->Scale(hNumTracksInCylinder2TG->Integral() / hMCNumTracksInCylinder2TG->Integral());

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
        // TOF
        {hTOFMass},
        {hTOF},

        // WC quality
        {hNumWC2TPCMatch, hMCNumWC2TPCMatch},
        {hNumWCHits},
        {hRadDistWC4},
        {hRadDistMidPlane},

        // Cylinder
        {hNumTracksInCylinder, hMCNumTracksInCylinder},

        // Comparisons with MC
        {hTGSmallTracks, hMCTGSmallTracks},
        {hTracksNearVertex, hMCTracksNearVertex},
        {hTrackLengthsNearVertex, hMCTrackLengthsNearVertex},
        {hNumTGTracks, hMCNumTGTracks},
        {hShowerProb, hMCShowerProb},
        {hBeforeShowerCutSmallTracks, hAfterShowerCutSmallTracks, hMCBeforeShowerCutSmallTracks, hMCAfterShowerCutSmallTracks},
        {hTGTrackLengths, hMCTGTrackLengths},

        // Comparing TG threshold
        {hTGSmallTracks0TG, hMCTGSmallTracks0TG},
        {hTGSmallTracks1TG, hMCTGSmallTracks1TG},
        {hTGSmallTracks2TG, hMCTGSmallTracks2TG},
        {hNumTracksInCylinder0TG, hMCNumTracksInCylinder0TG},
        {hNumTracksInCylinder1TG, hMCNumTracksInCylinder1TG},
        {hNumTracksInCylinder2TG, hMCNumTracksInCylinder2TG}
    };

    std::vector<std::vector<TString>> PlotLabelGroups = {
        // TOF
        {"Data"},
        {"Data"},

        // WC quality
        {"Data", "MC (scaled)"},
        {"Data"},
        {"Data"},
        {"Data"},

        // Cylinder
        {"Data", "MC (scaled)"},

        // Comparisons with MC
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data before", "Data after", "MC before (scaled)", "MC after (scaled)"},
        {"Data", "MC (scaled)"},

        // Comparing TG threshold
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"}
    };

    std::vector<TString> PlotTitles = {
        // TOF
        "TOF/TOFMass",
        "TOF/TOF",

        // WC quality
        "WCQuality/NumMatches",
        "WCQuality/NumWCHits",
        "WCQuality/RadDistWC4Proj",
        "WCQuality/RadDistMidPlane",

        // Cylinder
        "Cylinder/NumTracksInCylinder",

        // Comparisons with MC
        "TGPrimary/TGSmallTracks",
        "NearVertex/TracksNearVertex",
        "NearVertex/TrackLengthsNearVertex",
        "TGTracks/NumTGTracks",
        "ShowerCut/ShowerProb",
        "ShowerCut/SmallTracksShowerCut",
        "TGPrimary/TGTrackLengths",

        // Comparing TG threshold
        "TGThreshold/TGSmallTracks0TG",
        "TGThreshold/TGSmallTracks1TG",
        "TGThreshold/TGSmallTracks2TG",
        "TGThreshold/NumTracksInCylinder0TG",
        "TGThreshold/NumTracksInCylinder1TG",
        "TGThreshold/NumTracksInCylinder2TG"
    };

    std::vector<TString> XLabels = {
        // TOF
        "Mass [MeV/c^2]",
        "Time of flight [ns]",

        // WC quality
        "# of WC to TPC matches",
        "# of WC hits",
        "Proj. to WC hit distance [cm]",
        "Upstream to downstream proj. distance [cm]",

        // Cylinder
        "# of tracks",

        // Comparisons with MC
        "# of small tracks",
        "# of tracks near vertex",
        "Track length [cm]",
        "# of throughgoing tracks",
        "Shower probability",
        "# of small tracks (first 30 cm)",
        "Track length [cm]",

        // Comparing TG threshold
        "# of small tracks",
        "# of small tracks",
        "# of small tracks",
        "# of tracks",
        "# of tracks",
        "# of tracks",
    };

    std::vector<TString> YLabels = {
        // TOF
        "Counts",
        "Counts",

        // WC quality
        "Counts",
        "Counts",
        "Counts",
        "Counts",

        // Cylinder
        "Counts",

        // Comparisons with MC
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",

        // Comparing TG threshold
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts",
        "Counts"
    };

    std::vector<bool> PlotStacked = {
        // TOF
        false,
        false,

        // WC quality
        false,
        false,
        false,
        false,

        // Cylinder
        false,

        // Comparisons with MC
        false,
        false,
        false,
        false,
        false,
        false,
        false,

        // Comparing TG threshold
        false,
        false,
        false,
        false,
        false,
        false
    };

    std::vector<std::vector<bool>> PlotsAsPoints = {
        // Data plots
        {true},
        {true},

        // WC quality
        {true, false},
        {true},
        {true},
        {true},

        // Cylinder
        {true, false},

        // Comparisons with MC
        {true, false},
        {true, false},
        {true, false},
        {true, false},
        {true, false},
        {true, true, false, false},
        {true, false},

        // Comparing TG threshold
        {true, false},
        {true, false},
        {true, false},
        {true, false},
        {true, false},
        {true, false}
    };

    /////////////////////////////////
    // Fractional difference plots //
    /////////////////////////////////
    
    std::vector<std::pair<TH1*,TH1*>> PlotGroupsFracDiff = {
        {hTGSmallTracks, hMCTGSmallTracks},
        {hTracksNearVertex, hMCTracksNearVertex},
        {hTrackLengthsNearVertex, hMCTrackLengthsNearVertex},
        {hNumTGTracks, hMCNumTGTracks},
        {hBeforeShowerCutSmallTracks, hMCBeforeShowerCutSmallTracks},
        {hAfterShowerCutSmallTracks, hMCAfterShowerCutSmallTracks},
        {hTGTrackLengths, hMCTGTrackLengths},
    };

    std::vector<std::string> FracDiffDirectory = {
        "TGPrimary",
        "NearVertex",
        "NearVertex",
        "TGTracks",
        "ShowerCut",
        "ShowerCut",
        "TGPrimary"
    };

    for (int i = 0; i < (int)PlotGroupsFracDiff.size(); ++i) {
        TH1* hData = PlotGroupsFracDiff[i].first;
        TH1* hMC   = PlotGroupsFracDiff[i].second;
        if (!hData || !hMC) continue;

        TH1D* hFrac = (TH1D*)hData->Clone(
            (std::string(hData->GetName()) + "_FracVsMC").c_str()
        );
        hFrac->Reset("ICE");
        hFrac->Sumw2(kTRUE);
        hFrac->GetYaxis()->SetTitle("(D - MC) / MC");

        const int nb = hData->GetNbinsX();
        for (int b = 1; b <= nb; ++b) {
            double D  = hData->GetBinContent(b);
            double S  = hMC->GetBinContent(b);
            double sD = hData->GetBinError(b);
            double sS = hMC->GetBinError(b);

            if (S != 0.0) {
                double f = (D - S) / S;
                // df/dD = 1/S ; df/dS = -D/S^2
                double dfdD = 1.0 / S;
                double dfdS = -D / (S*S);
                double var = dfdD*dfdD*sD*sD + dfdS*dfdS*sS*sS; // rho=0
                hFrac->SetBinContent(b, f);
                hFrac->SetBinError(b, (var > 0 ? std::sqrt(var) : 0.0));
            } else {
                hFrac->SetBinContent(b, 0.0);
                hFrac->SetBinError(b,   0.0);
            }
        }

        PlotGroups.push_back({hFrac});
        PlotLabelGroups.push_back({"Frac. diff."});
        PlotTitles.push_back(FracDiffDirectory[i] + "/" + std::string(hData->GetTitle()) + "FracDiff");
        XLabels.push_back(hData->GetXaxis()->GetTitle());
        YLabels.push_back("(D - MC) / MC");
        PlotStacked.push_back(false);
        PlotsAsPoints.push_back({true});
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

    hMCSmallVsTGTracks->Scale(scalingTG);
    hMCTGNumSmallTracksVsThresh->Scale(scalingTG);

    std::vector<TH2*> TwoDPlots = {
        hSmallVsTGTracks,
        hMCSmallVsTGTracks,
        hTGNumSmallTracksVsThresh,
        hMCTGNumSmallTracksVsThresh,
        hTOFVsTOFMass,
        hPrimaryTrackPosition,
        hMCPrimaryTrackPosition,
        hBackgroundTracksPosition,
        hBackgroundTracksDirection,
        hProjVsRealWC4,
        hMidPlaneCoinc
    };

    std::vector<TString> TwoDTitles = {
        "TGTracks/SmallVsTGTracks",
        "TGTracks/MCSmallVsTGTracks",
        "TGTracks/TGNumSmallTracksVsThresh",
        "TGTracks/MCTGNumSmallTracksVsThresh",
        "TOF/TOFVsTOFMass",
        "PrimaryTrack/IncidentPosition",
        "PrimaryTrack/MCIncidentPosition",
        "BkgTracks/IncidentPosition",
        "BkgTracks/IncidentDirection",
        "WCQuality/ProjVsRealWC4",
        "WCQuality/MidPlaneCoinc"
    };
    std::vector<std::pair<double,double>> TwoDRanges = {
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0}
    };

    std::vector<bool> TwoDDisplayNumbers = {
        true,
        true,
        true,
        true,
        false,
        false,
        false,
        false,
        false,
        false
    };

    printTwoDPlots(SaveDir, TwoDPlots, TwoDTitles, TwoDRanges, TwoDDisplayNumbers);
}