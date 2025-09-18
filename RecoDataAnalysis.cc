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

    TH1D* hMCNumWC2TPCMatch         = (TH1D*) MCFile->Get("hMCNumWC2TPCMatch");
    TH1D* hMCNumTracksInCylinder    = (TH1D*) MCFile->Get("hMCNumTracksInCylinder");
    TH1D* hMCTGTrackLengths         = (TH1D*) MCFile->Get("hMCTGTrackLengths");
    TH1D* hMCTGSmallTracks          = (TH1D*) MCFile->Get("hMCTGSmallTracks");
    TH1D* hMCTracksNearVertex       = (TH1D*) MCFile->Get("hMCTracksNearVertex");
    TH1D* hMCTrackLengthsNearVertex = (TH1D*) MCFile->Get("hMCTrackLengthsNearVertex");
    TH1D* hMCNumTGTracks            = (TH1D*) MCFile->Get("hMCNumTGTracks");
    TH1D* hMCShowerProb             = (TH1D*) MCFile->Get("hMCShowerProb");

    TH1D* hMCBeforeShowerCutSmallTracks = (TH1D*) MCFile->Get("hMCBeforeShowerCutSmallTracks");
    TH1D* hMCAfterShowerCutSmallTracks  = (TH1D*) MCFile->Get("hMCAfterShowerCutSmallTracks");

    TH2D* hMCSmallVsTGTracks          = (TH2D*) MCFile->Get("hMCSmallVsTGTracks");
    TH2D* hMCTGNumSmallTracksVsThresh = (TH2D*) MCFile->Get("hMCTGNumSmallTracksVsThresh");

    TH2D* hMCPrimaryTrackPosition = (TH2D*) MCFile->Get("hMCPrimaryTrackPosition");

    ///////////////////////
    // Create histograms //
    ///////////////////////

    TH1D* hTOFMass = new TH1D("hTOFMass", "TOF Mass Distribution", 50, 0, 1200);
    TH1D* hTOF     = new TH1D("hTOF", "TOF Distribution", 50, 0, 80);
    
    TH1D* hNumWC2TPCMatch      = new TH1D("hNumWC2TPCMatch", "NumWC2TPCMatch", 10, 0, 10);
    TH1D* hNumTracksInCylinder = new TH1D("hNumTracksInCylinder", "NumTracksInCylinder", 10, 0, 10);

    TH1D* hTGTrackLengths         = new TH1D("hTGTrackLengths", "TGTrackLengths", 25, 0, 50);
    TH1D* hTGSmallTracks          = new TH1D("hTGSmallTracks", "TGSmallTracks", 10, 0, 10);
    TH1D* hTracksNearVertex       = new TH1D("hTracksNearVertex", "TracksNearVertex", 10, 0, 10);
    TH1D* hTrackLengthsNearVertex = new TH1D("hTrackLengthsNearVertex", "TrackLengthsNearVertex", 50, 0, 100);
    TH1D* hNumTGTracks            = new TH1D("hNumTGTracks", "NumTGTracks", 10, 0, 10);
    TH1D* hShowerProb             = new TH1D("hShowerProb", "ShowerProb", 20, 0, 1.);

    TH1D* hBeforeShowerCutSmallTracks = new TH1D("hBeforeShowerCutSmallTracks", "BeforeShowerCutSmallTracks", 10, 0, 10);
    TH1D* hAfterShowerCutSmallTracks  = new TH1D("hAfterShowerCutSmallTracks", "AfterShowerCutSmallTracks", 10, 0, 10);

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
        "hBackgroundTracksDirection", "BackgroundTracksDirection;Azimuth;Polar",
        20, -TMath::Pi(), -TMath::Pi(),
        20, 0, TMath::Pi()
    );

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
        hTOF->Fill(tofObject);
        hTOFVsTOFMass->Fill(tofObject, std::abs(TOFMass));

        if (std::abs(TOFMass) > PI_MU_EL_MASS_CUTOFF) {
            // Not a candidate pion, muon, or electron
            continue;
        }

        // If no track matched to wire-chamber, skip
        if (WC2TPCtrkID == -99999) continue;
        hPrimaryTrackPosition->Fill(WC2TPCPrimaryBeginX, WC2TPCPrimaryBeginY);
        numValidEvents++;

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
                // std::cout << "Event " << event << ", trk " << trk_idx << ": beginZ = " << recoBeginZ->at(trk_idx) << ", endZ = " << recoEndZ->at(trk_idx) << std::endl;
                // std::cout << "     begin x = " << recoBeginX->at(trk_idx) << ", y = " << recoBeginY->at(trk_idx) << ", z = " << recoBeginZ->at(trk_idx) << std::endl;
                // std::cout << "     end   x = " << recoEndX->at(trk_idx) << ", y = " << recoEndY->at(trk_idx) << ", z = " << recoEndZ->at(trk_idx) << std::endl;

                countEarlyTracks++;
            }
        }
        // std::cout << "Event " << event << " has " << countEarlyTracks << " tracks in first 14 cm of TPC" << std::endl;
        // std::cout << std::endl;
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

            // Is track contained in 10 cm cylinder?
            bool startInCylinder = IsPointInsideTrackCylinder(
                WC2TPCLocationsX, WC2TPCLocationsY, WC2TPCLocationsZ,
                recoBeginX->at(trk_idx), recoBeginY->at(trk_idx), recoBeginZ->at(trk_idx),
                CYLINDER_RADIUS
            );
            bool endInCylinder = IsPointInsideTrackCylinder(
                WC2TPCLocationsX, WC2TPCLocationsY, WC2TPCLocationsZ,
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
                    double incidentX = isBeginTrue ? recoBeginX->at(trk_idx) : recoEndX->at(trk_idx);
                    double incidentY = isBeginTrue ? recoBeginY->at(trk_idx) : recoEndY->at(trk_idx);

                    hBackgroundTracksPosition->Fill(incidentX, incidentY);

                    auto [phi, theta] = azimuth_polar_from_points(
                        recoBeginX->at(trk_idx), recoBeginY->at(trk_idx), recoBeginZ->at(trk_idx),
                        recoEndX->at(trk_idx), recoEndY->at(trk_idx), recoEndZ->at(trk_idx)
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

        hBeforeShowerCutSmallTracks->Fill(smallTracksTPCStart);
        hNumTracksInCylinder->Fill(numTracksInCylinder);

        if (!isPrimaryTG) hTracksNearVertex->Fill(numTracksNearVertex);
        if (obtainedProbabilities && showerProb < SHOWER_PROB_CUT) hAfterShowerCutSmallTracks->Fill(smallTracksTPCStart);

        // Add to histogram of small tracks if primary is throughgoing
        if (isPrimaryTG) {
            hTGSmallTracks->Fill(numSmallTracks);
            hNumTGTracks->Fill(numTGTracks);
            hSmallVsTGTracks->Fill(numSmallTracks, numTGTracks);

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

    // Scale MC histograms to data histograms using event counts
    hMCNumTGTracks->Scale(scaling);
    hMCBeforeShowerCutSmallTracks->Scale(scaling);
    hMCAfterShowerCutSmallTracks->Scale(scaling); // not sure

    hMCTGSmallTracks->Scale(scalingTG);
    hMCTGTrackLengths->Scale(scalingTG);
    hMCNumTracksInCylinder->Scale(scalingTG);
    
    hMCTracksNearVertex->Scale(scalingNoTG);
    hMCTrackLengthsNearVertex->Scale(scalingNoTG);

    // Have to account only for events with valid NN probability
    hMCShowerProb->Scale(hShowerProb->Integral() / hMCShowerProb->Integral());

    hMCNumWC2TPCMatch->Scale(hNumWC2TPCMatch->Integral() / hMCNumWC2TPCMatch->Integral());

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

        // WC2TPC
        {hNumWC2TPCMatch, hMCNumWC2TPCMatch},

        // Cylinder
        {hNumTracksInCylinder, hMCNumTracksInCylinder},

        // Comparisons with MC
        {hTGSmallTracks, hMCTGSmallTracks},
        {hTracksNearVertex, hMCTracksNearVertex},
        {hTrackLengthsNearVertex, hMCTrackLengthsNearVertex},
        {hNumTGTracks, hMCNumTGTracks},
        {hShowerProb, hMCShowerProb},
        {hBeforeShowerCutSmallTracks, hAfterShowerCutSmallTracks, hMCBeforeShowerCutSmallTracks, hMCAfterShowerCutSmallTracks},
        {hTGTrackLengths, hMCTGTrackLengths}
    };

    std::vector<std::vector<TString>> PlotLabelGroups = {
        // TOF
        {"Data"},
        {"Data"},

        // WC2TPC
        {"Data", "MC (scaled)"},

        // Cylinder
        {"Data", "MC (scaled)"},

        // Comparisons with MC
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data", "MC (scaled)"},
        {"Data before", "Data after", "MC before (scaled)", "MC after (scaled)"},
        {"Data", "MC (scaled)"}
    };

    std::vector<TString> PlotTitles = {
        // TOF
        "TOF/TOFMass",
        "TOF/TOF",

        // WC2TPC
        "WC2TPC/NumMatches",

        // Cylinder
        "Cylinder/NumTracksInCylinder",

        // Comparisons with MC
        "TGPrimary/TGSmallTracks",
        "NearVertex/TracksNearVertex",
        "NearVertex/TrackLengthsNearVertex",
        "TGTracks/NumTGTracks",
        "ShowerCut/ShowerProb",
        "ShowerCut/SmallTracksShowerCut",
        "TGPrimary/TGTrackLengths"
    };

    std::vector<TString> XLabels = {
        // TOF
        "Mass [MeV/c^2]",
        "Time of flight [ns]",

        // WC2TPC
        "# of WC to TPC matches",

        // Cylinder
        "# of tracks",

        // Comparisons with MC
        "# of small tracks",
        "# of tracks near vertex",
        "Track length [cm]",
        "# of throughgoing tracks",
        "Shower probability",
        "# of small tracks (first 30 cm)",
        "Track length [cm]"
    };

    std::vector<TString> YLabels = {
        // TOF
        "Counts",
        "Counts",

        // WC2TPC
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
        "Counts"
    };

    std::vector<bool> PlotStacked = {
        // TOF
        false,
        false,

        // WC2TPC
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
        false
    };

    std::vector<std::vector<bool>> PlotsAsPoints = {
        // Data plots
        {true},
        {true},

        // WC2TPC
        {true, false},

        // Cylinder
        {true, false},

        // Comparisons with MC
        {true, false},
        {true, false},
        {true, false},
        {true, false},
        {true, false},
        {true, true, false, false},
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
        {hTGTrackLengths, hMCTGTrackLengths}
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
        hBackgroundTracksDirection
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
    };
    std::vector<std::pair<double,double>> TwoDRanges = {
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
        false
    };

    printTwoDPlots(SaveDir, TwoDPlots, TwoDTitles, TwoDRanges, TwoDDisplayNumbers);
}