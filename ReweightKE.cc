#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TLegend.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLine.h"

#include <vector>

#include "Helpers.cc"

void ReweightKE() {
    gStyle->SetOptStat(0);
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int    fontStyle = 132;
    double textSize  = 0.06;
    TString saveDir  = "/exp/lariat/app/users/epelaez/analysis/figs/ReweightKE/";

    TGraph* gProton = new TGraph();
    TGraph* gPion   = new TGraph();
    TGraph* gMuonTG = new TGraph();
    initializeProtonPoints(gProton);
    initializePionPoints(gPion);
    initializeMuonNoBraggPoints(gMuonTG);

    // ── Load files ────────────────────────────────────────────────────────────
    TString dataPath = "/exp/lariat/app/users/epelaez/files/anatree_60a_data/histo.root";
    std::unique_ptr<TFile> fData(TFile::Open(dataPath));
    TDirectory* dir = (TDirectory*) fData->Get("anatree");

    TString nomPath = "/exp/lariat/app/users/epelaez/histos/nominal/RecoClassify3Cat_AllHists.root";
    std::unique_ptr<TFile> fNom(TFile::Open(nomPath, "READ"));

    // ── Load MC histograms ────────────────────────────────────────────────────
    TH1D* hMCFrontFaceKEPion     = dynamic_cast<TH1D*>(fNom->Get("hFrontFaceKEPion"));
    TH1D* hMCFrontFaceKEMuon     = dynamic_cast<TH1D*>(fNom->Get("hFrontFaceKEMuon"));
    TH1D* hMCFrontFaceKEElectron = dynamic_cast<TH1D*>(fNom->Get("hFrontFaceKEElectron"));

    TH1D* hMCWCKEPion     = dynamic_cast<TH1D*>(fNom->Get("hWCKEPion"));
    TH1D* hMCWCKEMuon     = dynamic_cast<TH1D*>(fNom->Get("hWCKEMuon"));
    TH1D* hMCWCKEElectron = dynamic_cast<TH1D*>(fNom->Get("hWCKEElectron"));

    // ── Load TTree and branches ───────────────────────────────────────────────
    TTree* tree = (TTree*) dir->Get<TTree>("anatree");

    int run, subrun, event;
    tree->SetBranchAddress("run",    &run);
    tree->SetBranchAddress("subrun", &subrun);
    tree->SetBranchAddress("event",  &event);

    const int kMaxTrack    = 100;
    const int kMaxWCTracks = 100;

    int   ntracks_reco;                tree->SetBranchAddress("ntracks_reco",   &ntracks_reco);
    int   trkWCtoTPCMatch[kMaxTrack];  tree->SetBranchAddress("trkWCtoTPCMatch", &trkWCtoTPCMatch);

    float beamline_mass;               tree->SetBranchAddress("beamline_mass",   &beamline_mass);
    int   nwctrks;                     tree->SetBranchAddress("nwctrks",         &nwctrks);
    float wctrk_momentum[kMaxWCTracks]; tree->SetBranchAddress("wctrk_momentum", &wctrk_momentum);
    float wctrk_theta[kMaxWCTracks];    tree->SetBranchAddress("wctrk_theta",    &wctrk_theta);
    float wctrk_phi[kMaxWCTracks];      tree->SetBranchAddress("wctrk_phi",      &wctrk_phi);
    int   wctrk_picky[kMaxWCTracks];    tree->SetBranchAddress("wctrk_picky",    &wctrk_picky);
    float WC1xPos[kMaxWCTracks];        tree->SetBranchAddress("WC1xPos",        &WC1xPos);
    float WC1yPos[kMaxWCTracks];        tree->SetBranchAddress("WC1yPos",        &WC1yPos);
    float WC1zPos[kMaxWCTracks];        tree->SetBranchAddress("WC1zPos",        &WC1zPos);
    float WC2xPos[kMaxWCTracks];        tree->SetBranchAddress("WC2xPos",        &WC2xPos);
    float WC2yPos[kMaxWCTracks];        tree->SetBranchAddress("WC2yPos",        &WC2yPos);
    float WC2zPos[kMaxWCTracks];        tree->SetBranchAddress("WC2zPos",        &WC2zPos);
    float WC3xPos[kMaxWCTracks];        tree->SetBranchAddress("WC3xPos",        &WC3xPos);
    float WC3yPos[kMaxWCTracks];        tree->SetBranchAddress("WC3yPos",        &WC3yPos);
    float WC3zPos[kMaxWCTracks];        tree->SetBranchAddress("WC3zPos",        &WC3zPos);
    float WC4xPos[kMaxWCTracks];        tree->SetBranchAddress("WC4xPos",        &WC4xPos);
    float WC4yPos[kMaxWCTracks];        tree->SetBranchAddress("WC4yPos",        &WC4yPos);
    float WC4zPos[kMaxWCTracks];        tree->SetBranchAddress("WC4zPos",        &WC4zPos);

    const int kMaxTOF = 10;
    int   ntof;           tree->SetBranchAddress("ntof", &ntof);
    float tof[kMaxTOF];   tree->SetBranchAddress("tof",  tof);

    TH1D* hFrontFaceKE = new TH1D("hFrontFaceKE", "hFrontFaceKE;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());
    TH1D* hWCKE        = new TH1D("hWCKE",        "hWCKE;;",        NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());

    //////////////////////
    // Loop over events //
    //////////////////////

    Int_t numEntries = (Int_t) tree->GetEntries();
    std::cout << "Entries: " << numEntries << std::endl;

    bool verbose = false;

    for (Int_t i = 0; i < numEntries; ++i) {
        // For some reason, these events crash
        if (
            i == 199835 ||
            i == 311635
        ) continue;
        
        // Get tree entry and reset variables
        std::cout << std::endl;
        std::cout << "=================================" << std::endl;
        std::cout << "Getting tree entry: " << i << std::endl;
        tree->GetEntry(i);
        std::cout << "Got tree entry" << std::endl;
        std::cout << "Reseting variables" << std::endl;
        EventVariablesData ev;
        std::cout << "Variables reset" << std::endl;
        std::cout << "=================================" << std::endl;

        // Make script go faster
        // if (i > USE_NUM_EVENTS) break;

        int wcMatchCount = 0;
        for (int trk = 0; trk < ntracks_reco; ++trk) {
            if (trkWCtoTPCMatch[trk]) wcMatchCount++;
        }

        ev.TOFMass = std::abs(beamline_mass);

        if (nwctrks == 1) {
            ev.wcTrackPicky    = wctrk_picky[0];
            ev.WCTrackMomentum = wctrk_momentum[0];
            ev.WCTheta         = wctrk_theta[0];
            ev.WCPhi           = wctrk_phi[0];
            ev.WC4PrimaryX     = WC4xPos[0];

            ev.wcHit0 = {WC1xPos[0], WC1yPos[0], WC1zPos[0]};
            ev.wcHit1 = {WC2xPos[0], WC2yPos[0], WC2zPos[0]};
            ev.wcHit2 = {WC3xPos[0], WC3yPos[0], WC3zPos[0]};
            ev.wcHit3 = {WC4xPos[0], WC4yPos[0], WC4zPos[0]};
        }

        std::vector<double> midUp    = projToZ(ev.wcHit0, ev.wcHit1, -437.97);
        std::vector<double> projDown = projToZ(midUp, ev.wcHit2, -95.0);
        projDown[0] -= tan(1.32 * TMath::Pi() / 180.0) * (-95.0 - -437.97);
        double radDistWC4 = TMath::Sqrt(
            pow(projDown[0] - ev.wcHit3.at(0), 2.) +
            pow(projDown[1] - ev.wcHit3.at(1), 2.)
        );

        std::vector<double> midDown = projToZ(ev.wcHit2, ev.wcHit3, -437.97);
        midDown[0] -= tan(1.32 * TMath::Pi() / 180.0) * (-339.57 - -437.97);
        double midPlaneDist = TMath::Sqrt(
            pow(midUp[0] - midDown[0] + 0.75, 2) +
            pow(midUp[1] - midDown[1], 2)
        );

        if (radDistWC4   > 8.0) continue;
        if (midPlaneDist > 3.0) continue;

        if (!CheckUpstreamMagnetAperture      (ev.wcHit0, ev.wcHit1)) continue;
        if (!CheckDownstreamMagnetAperture    (ev.wcHit2, ev.wcHit3)) continue;
        if (!CheckDownstreamCollimatorAperture(ev.wcHit2, ev.wcHit3)) continue;

        if (std::abs(ev.TOFMass) > PI_MU_EL_MASS_CUTOFF) continue;
        if (wcMatchCount != 1)   continue;
        if (!ev.wcTrackPicky)    continue;

        double WCKE = TMath::Sqrt(
            ev.WCTrackMomentum * ev.WCTrackMomentum + PionMass * PionMass
        ) - PionMass;

        double tanThetaCosPhi = TMath::Tan(ev.WCTheta) * TMath::Cos(ev.WCPhi);
        double tanThetaSinPhi = TMath::Tan(ev.WCTheta) * TMath::Sin(ev.WCPhi);
        double den            = TMath::Sqrt(1 + tanThetaCosPhi * tanThetaCosPhi);
        double onTheFlyPz     = ev.WCTrackMomentum / den;
        double onTheFlyPx     = onTheFlyPz * tanThetaSinPhi;
        double calculatedEnLoss = energyLossCalculation(ev.WC4PrimaryX, onTheFlyPx, true);

        double initialKE = WCKE - calculatedEnLoss;

        hFrontFaceKE->Fill(initialKE);
        hWCKE->Fill(WCKE);
    }

    double scaleFrontFace = hFrontFaceKE->Integral() /
        (hMCFrontFaceKEPion->Integral() + hMCFrontFaceKEMuon->Integral() + hMCFrontFaceKEElectron->Integral());

    double scaleWCKE = hWCKE->Integral() /
        (hMCWCKEPion->Integral() + hMCWCKEMuon->Integral() + hMCWCKEElectron->Integral());

    hMCFrontFaceKEPion    ->Scale(scaleFrontFace);
    hMCFrontFaceKEMuon    ->Scale(scaleFrontFace);
    hMCFrontFaceKEElectron->Scale(scaleFrontFace);
    hMCWCKEPion           ->Scale(scaleWCKE);
    hMCWCKEMuon           ->Scale(scaleWCKE);
    hMCWCKEElectron       ->Scale(scaleWCKE);

    // ── Build total MC for reweighting ────────────────────────────────────────
    TH1D* hMCFrontFaceTotal = (TH1D*) hMCFrontFaceKEPion->Clone("hMCFrontFaceTotal");
    hMCFrontFaceTotal->Add(hMCFrontFaceKEMuon);
    hMCFrontFaceTotal->Add(hMCFrontFaceKEElectron);
    hMCFrontFaceTotal->SetDirectory(nullptr);

    TH1D* hMCWCKETotal = (TH1D*) hMCWCKEPion->Clone("hMCWCKETotal");
    hMCWCKETotal->Add(hMCWCKEMuon);
    hMCWCKETotal->Add(hMCWCKEElectron);
    hMCWCKETotal->SetDirectory(nullptr);

    // ── Compute weights ───────────────────────────────────────────────────────
    TH1D* hWeightsFrontFace = (TH1D*) hFrontFaceKE->Clone("hWeightsFrontFace");
    TH1D* hWeightsWCKE      = (TH1D*) hWCKE->Clone("hWeightsWCKE");
    hWeightsFrontFace->SetDirectory(nullptr);
    hWeightsWCKE->SetDirectory(nullptr);

    ComputeReweightingWeights(hMCFrontFaceTotal, hFrontFaceKE, hWeightsFrontFace);
    ComputeReweightingWeights(hMCWCKETotal,      hWCKE,        hWeightsWCKE);

    // ── Build reweighted MC species (clones with weights applied) ─────────────
    TH1D* hMCFrontFaceKEPionRew     = (TH1D*) hMCFrontFaceKEPion    ->Clone("hMCFrontFaceKEPionRew");
    TH1D* hMCFrontFaceKEMuonRew     = (TH1D*) hMCFrontFaceKEMuon    ->Clone("hMCFrontFaceKEMuonRew");
    TH1D* hMCFrontFaceKEElectronRew = (TH1D*) hMCFrontFaceKEElectron->Clone("hMCFrontFaceKEElectronRew");
    TH1D* hMCWCKEPionRew            = (TH1D*) hMCWCKEPion            ->Clone("hMCWCKEPionRew");
    TH1D* hMCWCKEMuonRew            = (TH1D*) hMCWCKEMuon            ->Clone("hMCWCKEMuonRew");
    TH1D* hMCWCKEElectronRew        = (TH1D*) hMCWCKEElectron        ->Clone("hMCWCKEElectronRew");

    for (auto* h : {hMCFrontFaceKEPionRew, hMCFrontFaceKEMuonRew, hMCFrontFaceKEElectronRew})
        h->SetDirectory(nullptr);
    for (auto* h : {hMCWCKEPionRew, hMCWCKEMuonRew, hMCWCKEElectronRew})
        h->SetDirectory(nullptr);

    ApplyWeights(hMCFrontFaceKEPionRew,     hWeightsFrontFace);
    ApplyWeights(hMCFrontFaceKEMuonRew,     hWeightsFrontFace);
    ApplyWeights(hMCFrontFaceKEElectronRew, hWeightsFrontFace);
    ApplyWeights(hMCWCKEPionRew,            hWeightsWCKE);
    ApplyWeights(hMCWCKEMuonRew,            hWeightsWCKE);
    ApplyWeights(hMCWCKEElectronRew,        hWeightsWCKE);

    // ── Reweighted totals (for the comparison plot) ───────────────────────────
    TH1D* hMCFrontFaceRewTotal = (TH1D*) hMCFrontFaceKEPionRew->Clone("hMCFrontFaceRewTotal");
    hMCFrontFaceRewTotal->Add(hMCFrontFaceKEMuonRew);
    hMCFrontFaceRewTotal->Add(hMCFrontFaceKEElectronRew);
    hMCFrontFaceRewTotal->SetDirectory(nullptr);

    TH1D* hMCWCKERewTotal = (TH1D*) hMCWCKEPionRew->Clone("hMCWCKERewTotal");
    hMCWCKERewTotal->Add(hMCWCKEMuonRew);
    hMCWCKERewTotal->Add(hMCWCKEElectronRew);
    hMCWCKERewTotal->SetDirectory(nullptr);

    // ── Colour palette ────────────────────────────────────────────────────────
    // Species colours (same for original and reweighted breakdowns)
    std::vector<int> speciesColors = {kAzure - 4, kOrange - 3, kGreen - 6};
    std::vector<TString> speciesLabels = {"Pion", "Muon", "Electron"};

    // Comparison colours: original (blue) vs reweighted (red)
    std::vector<int> compColors  = {kAzure + 1, kRed - 4};
    std::vector<TString> compLabels = {"Original MC", "Reweighted MC"};

    // ─────────────────────────────────────────────────────────────────────────
    // FRONT-FACE KE
    // Plot 1: Data vs. original total + reweighted total (comparison)
    // Plot 2: Data vs. original species breakdown
    // Plot 3: Data vs. reweighted species breakdown
    // ─────────────────────────────────────────────────────────────────────────

    // Plot 1 — comparison of totals
    PrintDataVsTwoMCPlot(
        saveDir, "FrontFaceKEComparison",
        hFrontFaceKE,
        hMCFrontFaceTotal, hMCFrontFaceRewTotal,
        "Original MC", "Reweighted MC",
        kAzure + 1, kRed - 4,
        "Front-Face KE: Original vs. Reweighted MC",
        "Kinetic energy [MeV]", "Counts",
        fontStyle, textSize
    );

    // Plot 2 — original species breakdown
    PrintDataVsMCContribPlot(
        saveDir, 
        "FrontFaceKEBreakdown",
        hFrontFaceKE,
        {hMCFrontFaceKEPion, hMCFrontFaceKEMuon, hMCFrontFaceKEElectron},
        speciesLabels,
        speciesColors,
        "Front-Face KE: Original MC Breakdown",
        "Kinetic energy [MeV]", "Counts",
        fontStyle, textSize,
        /*UsePoissonErrors=*/ true,
        /*hMCAbsUnc=*/        nullptr,
        /*DrawRatio=*/        true
    );

    // Plot 3 — reweighted species breakdown
    PrintDataVsMCContribPlot(
        saveDir, "FrontFaceKEReweighBreakdown",
        hFrontFaceKE,
        {hMCFrontFaceKEPionRew, hMCFrontFaceKEMuonRew, hMCFrontFaceKEElectronRew},
        speciesLabels,
        speciesColors,
        "Front-Face KE: Reweighted MC Breakdown",
        "Kinetic energy [MeV]", "Counts",
        fontStyle, textSize,
        /*UsePoissonErrors=*/ true,
        /*hMCAbsUnc=*/        nullptr,
        /*DrawRatio=*/        true
    );

    // ─────────────────────────────────────────────────────────────────────────
    // WIRE-CHAMBER KE — same three plots
    // ─────────────────────────────────────────────────────────────────────────

    // Plot 1 — comparison of totals
    PrintDataVsTwoMCPlot(
        saveDir, "WCKEComparison",
        hWCKE,
        hMCWCKETotal, hMCWCKERewTotal,
        "Original MC", "Reweighted MC",
        kAzure + 1, kRed - 4,
        "Wire-Chamber KE: Original vs. Reweighted MC",
        "Kinetic energy [MeV]", "Counts",
        fontStyle, textSize
    );

    // Plot 2 — original species breakdown
    PrintDataVsMCContribPlot(
        saveDir, "WCKEBreakdown",
        hWCKE,
        {hMCWCKEPion, hMCWCKEMuon, hMCWCKEElectron},
        speciesLabels,
        speciesColors,
        "Wire-Chamber KE: Original MC Breakdown",
        "Kinetic energy [MeV]", "Counts",
        fontStyle, textSize,
        /*UsePoissonErrors=*/ true,
        /*hMCAbsUnc=*/        nullptr,
        /*DrawRatio=*/        true
    );

    // Plot 3 — reweighted species breakdown
    PrintDataVsMCContribPlot(
        saveDir, "WCKEReweighBreakdown",
        hWCKE,
        {hMCWCKEPionRew, hMCWCKEMuonRew, hMCWCKEElectronRew},
        speciesLabels,
        speciesColors,
        "Wire-Chamber KE: Reweighted MC Breakdown",
        "Kinetic energy [MeV]", "Counts",
        fontStyle, textSize,
        /*UsePoissonErrors=*/ true,
        /*hMCAbsUnc=*/        nullptr,
        /*DrawRatio=*/        true
    );

    // Save weights to file
    TString outPath = "/exp/lariat/app/users/epelaez/histos/data/BeamlineWeights.root";
    TFile outFile(outPath, "RECREATE");
    outFile.cd();
    hWeightsFrontFace->Write("hWeightsFrontFace", TObject::kOverwrite);
    hWeightsWCKE     ->Write("hWeightsWCKE",      TObject::kOverwrite);
    outFile.Close();

    std::cout << "\nWeights saved to " << outPath << std::endl;
}