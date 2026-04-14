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

    // Load files
    TChain* Chain = new TChain("anatree/anatree");
    Chain->Add("/exp/lariat/app/users/epelaez/files/anatree_60a_data/chunks/*.root");
    std::cout << "Files:   " << Chain->GetListOfFiles()->GetEntries() << std::endl;

    // Load branches
    int run, subrun, event;
    Chain->SetBranchAddress("run",    &run);
    Chain->SetBranchAddress("subrun", &subrun);
    Chain->SetBranchAddress("event",  &event);

    int   ntracks_reco;                Chain->SetBranchAddress("ntracks_reco",   &ntracks_reco);
    int   trkWCtoTPCMatch[kMaxTrackData];  Chain->SetBranchAddress("trkWCtoTPCMatch", &trkWCtoTPCMatch);

    float beamline_mass;               Chain->SetBranchAddress("beamline_mass",   &beamline_mass);
    int   nwctrks;                     Chain->SetBranchAddress("nwctrks",         &nwctrks);
    float wctrk_momentum[kMaxWCTracksData]; Chain->SetBranchAddress("wctrk_momentum", &wctrk_momentum);
    float wctrk_theta[kMaxWCTracksData];    Chain->SetBranchAddress("wctrk_theta",    &wctrk_theta);
    float wctrk_phi[kMaxWCTracksData];      Chain->SetBranchAddress("wctrk_phi",      &wctrk_phi);
    int   wctrk_picky[kMaxWCTracksData];    Chain->SetBranchAddress("wctrk_picky",    &wctrk_picky);
    float WC1xPos[kMaxWCTracksData];        Chain->SetBranchAddress("WC1xPos",        &WC1xPos);
    float WC1yPos[kMaxWCTracksData];        Chain->SetBranchAddress("WC1yPos",        &WC1yPos);
    float WC1zPos[kMaxWCTracksData];        Chain->SetBranchAddress("WC1zPos",        &WC1zPos);
    float WC2xPos[kMaxWCTracksData];        Chain->SetBranchAddress("WC2xPos",        &WC2xPos);
    float WC2yPos[kMaxWCTracksData];        Chain->SetBranchAddress("WC2yPos",        &WC2yPos);
    float WC2zPos[kMaxWCTracksData];        Chain->SetBranchAddress("WC2zPos",        &WC2zPos);
    float WC3xPos[kMaxWCTracksData];        Chain->SetBranchAddress("WC3xPos",        &WC3xPos);
    float WC3yPos[kMaxWCTracksData];        Chain->SetBranchAddress("WC3yPos",        &WC3yPos);
    float WC3zPos[kMaxWCTracksData];        Chain->SetBranchAddress("WC3zPos",        &WC3zPos);
    float WC4xPos[kMaxWCTracksData];        Chain->SetBranchAddress("WC4xPos",        &WC4xPos);
    float WC4yPos[kMaxWCTracksData];        Chain->SetBranchAddress("WC4yPos",        &WC4yPos);
    float WC4zPos[kMaxWCTracksData];        Chain->SetBranchAddress("WC4zPos",        &WC4zPos);

    int   ntof;           Chain->SetBranchAddress("ntof", &ntof);
    float tof[kMaxTOFData];   Chain->SetBranchAddress("tof",  tof);

    TH1D* hFrontFaceKE         = new TH1D("hFrontFaceKE", "hFrontFaceKE;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());
    TH1D* hFrontFaceKEPreMatch = new TH1D("hFrontFaceKEPreMatch", "hFrontFaceKEPreMatch;;", NUM_BINS_KE_FINE, ARRAY_KE_FINE_BINS.data());

    //////////////////////
    // Loop over events //
    //////////////////////

    bool verbose = false;

    Long64_t i = 0;
    while (true) {
        if (SKIP_INDICES_DATA.count(i)) { i++; continue; }
        if (Chain->GetEntry(i++) <= 0) break;
        
        // Reset variables
        if (verbose) std::cout << std::endl;
        if (verbose) std::cout << "=================================" << std::endl;
        if (verbose) std::cout << "Got tree entry: " << i << std::endl;
        EventVariablesData ev;
        if (verbose) std::cout << "Variables reset" << std::endl;
        if (verbose) std::cout << "=================================" << std::endl;

        // Make script go faster
        // if (i > USE_NUM_EVENTS) break;

        int wcMatchCount = 0;
        for (int trk = 0; trk < std::min(ntracks_reco, kMaxTrackData); ++trk) {
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
        if (!ev.wcTrackPicky) continue;

        double WCKE = TMath::Sqrt(
            ev.WCTrackMomentum * ev.WCTrackMomentum + PionMass * PionMass
        ) - PionMass;

        double tanThetaCosPhi   = TMath::Tan(ev.WCTheta) * TMath::Cos(ev.WCPhi);
        double tanThetaSinPhi   = TMath::Tan(ev.WCTheta) * TMath::Sin(ev.WCPhi);
        double den              = TMath::Sqrt(1 + tanThetaCosPhi * tanThetaCosPhi);
        double onTheFlyPz       = ev.WCTrackMomentum / den;
        double onTheFlyPx       = onTheFlyPz * tanThetaSinPhi;
        double calculatedEnLoss = energyLossCalculation(ev.WC4PrimaryX, onTheFlyPx, true);
        double frontFaceKE      = WCKE - calculatedEnLoss;
        hFrontFaceKEPreMatch->Fill(frontFaceKE);

        if (wcMatchCount != 1)   continue;

        hFrontFaceKE->Fill(frontFaceKE);
    }

    TString nomPath = "/exp/lariat/app/users/epelaez/histos/nominal/RecoClassify3Cat_AllHists.root";
    std::unique_ptr<TFile> fNom(TFile::Open(nomPath, "READ"));

    // Load MC histograms
    TH1D* hMCFrontFaceKEPion     = dynamic_cast<TH1D*>(fNom->Get("hFrontFaceKEPionNoWeight"));
    TH1D* hMCFrontFaceKEMuon     = dynamic_cast<TH1D*>(fNom->Get("hFrontFaceKEMuonNoWeight"));
    TH1D* hMCFrontFaceKEElectron = dynamic_cast<TH1D*>(fNom->Get("hFrontFaceKEElectronNoWeight"));

    TH1D* hMCFrontFaceKEPionPre     = dynamic_cast<TH1D*>(fNom->Get("hFrontFaceKEPionNoWeightPre"));
    TH1D* hMCFrontFaceKEMuonPre     = dynamic_cast<TH1D*>(fNom->Get("hFrontFaceKEMuonNoWeightPre"));
    TH1D* hMCFrontFaceKEElectronPre = dynamic_cast<TH1D*>(fNom->Get("hFrontFaceKEElectronNoWeightPre"));

    double scaleFrontFace    = hFrontFaceKE->Integral() / (hMCFrontFaceKEPion->Integral() + hMCFrontFaceKEMuon->Integral() + hMCFrontFaceKEElectron->Integral());
    double scaleFrontFacePre = hFrontFaceKEPreMatch->Integral() / (hMCFrontFaceKEPionPre->Integral() + hMCFrontFaceKEMuonPre->Integral() + hMCFrontFaceKEElectronPre->Integral());

    hMCFrontFaceKEPion    ->Scale(scaleFrontFace);
    hMCFrontFaceKEMuon    ->Scale(scaleFrontFace);
    hMCFrontFaceKEElectron->Scale(scaleFrontFace);

    hMCFrontFaceKEPionPre    ->Scale(scaleFrontFacePre);
    hMCFrontFaceKEMuonPre    ->Scale(scaleFrontFacePre);
    hMCFrontFaceKEElectronPre->Scale(scaleFrontFacePre);

    // Build total MC for reweighting
    TH1D* hMCFrontFaceTotal = (TH1D*) hMCFrontFaceKEPion->Clone("hMCFrontFaceTotal");
    hMCFrontFaceTotal->Add(hMCFrontFaceKEMuon);
    hMCFrontFaceTotal->Add(hMCFrontFaceKEElectron);
    hMCFrontFaceTotal->SetDirectory(nullptr);

    TH1D* hMCFrontFaceTotalPre = (TH1D*) hMCFrontFaceKEPionPre->Clone("hMCFrontFaceTotalPre");
    hMCFrontFaceTotalPre->Add(hMCFrontFaceKEMuonPre);
    hMCFrontFaceTotalPre->Add(hMCFrontFaceKEElectronPre);
    hMCFrontFaceTotalPre->SetDirectory(nullptr);

    // Compute weights
    TH1D* hWeightsFrontFace    = (TH1D*) hFrontFaceKE->Clone("hWeightsFrontFace");
    TH1D* hWeightsFrontFacePre = (TH1D*) hFrontFaceKEPreMatch->Clone("hWeightsFrontFacePre");
    hWeightsFrontFace->SetDirectory(nullptr);
    hWeightsFrontFacePre->SetDirectory(nullptr);

    ComputeReweightingWeights(hMCFrontFaceTotal, hFrontFaceKE, hWeightsFrontFace);
    ComputeReweightingWeights(hMCFrontFaceTotalPre, hFrontFaceKEPreMatch, hWeightsFrontFacePre);

    // Build reweighted MC species (clones with weights applied)
    TH1D* hMCFrontFaceKEPionRew     = (TH1D*) hMCFrontFaceKEPion    ->Clone("hMCFrontFaceKEPionRew");
    TH1D* hMCFrontFaceKEMuonRew     = (TH1D*) hMCFrontFaceKEMuon    ->Clone("hMCFrontFaceKEMuonRew");
    TH1D* hMCFrontFaceKEElectronRew = (TH1D*) hMCFrontFaceKEElectron->Clone("hMCFrontFaceKEElectronRew");

    TH1D* hMCFrontFaceKEPionRewPre     = (TH1D*) hMCFrontFaceKEPionPre    ->Clone("hMCFrontFaceKEPionRewPre");
    TH1D* hMCFrontFaceKEMuonRewPre     = (TH1D*) hMCFrontFaceKEMuonPre    ->Clone("hMCFrontFaceKEMuonRewPre");
    TH1D* hMCFrontFaceKEElectronRewPre = (TH1D*) hMCFrontFaceKEElectronPre->Clone("hMCFrontFaceKEElectronRewPre");

    for (auto* h : {hMCFrontFaceKEPionRew, hMCFrontFaceKEMuonRew, hMCFrontFaceKEElectronRew}) h->SetDirectory(nullptr);
    for (auto* h : {hMCFrontFaceKEPionRewPre, hMCFrontFaceKEMuonRewPre, hMCFrontFaceKEElectronRewPre}) h->SetDirectory(nullptr);

    ApplyWeights(hMCFrontFaceKEPionRew,     hWeightsFrontFace);
    ApplyWeights(hMCFrontFaceKEMuonRew,     hWeightsFrontFace);
    ApplyWeights(hMCFrontFaceKEElectronRew, hWeightsFrontFace);

    ApplyWeights(hMCFrontFaceKEPionRewPre,     hWeightsFrontFacePre);
    ApplyWeights(hMCFrontFaceKEMuonRewPre,     hWeightsFrontFacePre);
    ApplyWeights(hMCFrontFaceKEElectronRewPre, hWeightsFrontFacePre);

    // Reweighted totals (for the comparison plot)
    TH1D* hMCFrontFaceRewTotal = (TH1D*) hMCFrontFaceKEPionRew->Clone("hMCFrontFaceRewTotal");
    hMCFrontFaceRewTotal->Add(hMCFrontFaceKEMuonRew);
    hMCFrontFaceRewTotal->Add(hMCFrontFaceKEElectronRew);
    hMCFrontFaceRewTotal->SetDirectory(nullptr);

    TH1D* hMCFrontFaceRewTotalPre = (TH1D*) hMCFrontFaceKEPionRewPre->Clone("hMCFrontFaceRewTotalPre");
    hMCFrontFaceRewTotalPre->Add(hMCFrontFaceKEMuonRewPre);
    hMCFrontFaceRewTotalPre->Add(hMCFrontFaceKEElectronRewPre);
    hMCFrontFaceRewTotalPre->SetDirectory(nullptr);

    // Species colours (same for original and reweighted breakdowns)
    std::vector<int> speciesColors = {kAzure - 4, kOrange - 3, kGreen - 6};
    std::vector<TString> speciesLabels = {"Pion", "Muon", "Electron"};

    // Comparison colours: original (blue) vs reweighted (red)
    std::vector<int> compColors  = {kAzure + 1, kRed - 4};
    std::vector<TString> compLabels = {"Original MC", "Reweighted MC"};

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

    // Pre match
    PrintDataVsTwoMCPlot(
        saveDir, "FrontFaceKEComparisonPreMatch",
        hFrontFaceKEPreMatch,
        hMCFrontFaceTotalPre, hMCFrontFaceRewTotalPre,
        "Original MC", "Reweighted MC",
        kAzure + 1, kRed - 4,
        "Front-Face KE: Original vs. Reweighted MC",
        "Kinetic energy [MeV]", "Counts",
        fontStyle, textSize
    );
    PrintDataVsMCContribPlot(
        saveDir, 
        "FrontFaceKEBreakdownPreMatch",
        hFrontFaceKEPreMatch,
        {hMCFrontFaceKEPionPre, hMCFrontFaceKEMuonPre, hMCFrontFaceKEElectronPre},
        speciesLabels,
        speciesColors,
        "Front-Face KE: Original MC Breakdown",
        "Kinetic energy [MeV]", "Counts",
        fontStyle, textSize,
        /*UsePoissonErrors=*/ true,
        /*hMCAbsUnc=*/        nullptr,
        /*DrawRatio=*/        true
    );
    PrintDataVsMCContribPlot(
        saveDir, "FrontFaceKEReweighBreakdownPreMatch",
        hFrontFaceKEPreMatch,
        {hMCFrontFaceKEPionRewPre, hMCFrontFaceKEMuonRewPre, hMCFrontFaceKEElectronRewPre},
        speciesLabels,
        speciesColors,
        "Front-Face KE: Reweighted MC Breakdown",
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
    hWeightsFrontFacePre->Write("hWeightsFrontFacePre", TObject::kOverwrite);

    TH1D* hWeightsFrontFaceDummy = (TH1D*) hWeightsFrontFace->Clone("hWeightsFrontFaceDummy");
    for (int b = 1; b <= hWeightsFrontFaceDummy->GetNbinsX(); ++b) {
        hWeightsFrontFaceDummy->SetBinContent(b, 1.0);
    }
    hWeightsFrontFaceDummy->Write("hWeightsFrontFaceDummy", TObject::kOverwrite);
    outFile.Close();

    std::cout << "\nWeights saved to " << outPath << std::endl;
}