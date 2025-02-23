#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TEfficiency.h>

#include <vector>

void ShowerAnalysis() {
    // Set defaults
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    // gStyle->SetOptStat(0);

    int FontStyle = 132;
    double TextSize = 0.06;
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis/figs/ShowerAnalysis/";

    // Load root file
    TString RootFilePath = "/exp/lariat/app/users/epelaez/files/ShowerRecoEval_histo.root"; 
    std::unique_ptr<TFile> File(TFile::Open(RootFilePath));
    TDirectory* Directory = (TDirectory*)File->Get("ShowerRecoEval");

    // Load tree and branches
    TTree* tree = (TTree*) Directory->Get<TTree>("ShowerRecoEvalTree");

    int run, subrun, event;
    tree->SetBranchAddress("run", &run); 
    tree->SetBranchAddress("subrun", &subrun); 
    tree->SetBranchAddress("event", &event);

    int primaryParticleTrackID, numEmmitedElectrons, numEmmitedPhotons;
    tree->SetBranchAddress("primaryParticleTrackID", &primaryParticleTrackID);
    tree->SetBranchAddress("numEmmitedElectrons", &numEmmitedElectrons);
    tree->SetBranchAddress("numEmmitedPhotons", &numEmmitedPhotons);

    std::vector<double>* truthElectronsLength = nullptr;
    tree->SetBranchAddress("truthElectronsLength", &truthElectronsLength);

    std::vector<double>* recoLength = nullptr;
    std::vector<int>* matchedIdentity = nullptr;
    std::vector<double>* matchedLength = nullptr;
    tree->SetBranchAddress("recoLength", &recoLength);
    tree->SetBranchAddress("matchedIdentity", &matchedIdentity);
    tree->SetBranchAddress("matchedLength", &matchedLength);

    // Declare histograms
    TH1D *hShowerRecoLengths = new TH1D("hShowerRecoLengths", "hShowerRecoLengths;;", 20, -1,-1);
    TH1D *hNoShowerRecoLengths = new TH1D("hNoShowerRecoLengths", "hNoShowerRecoLengths;;", 20, -1,-1);

    TH1D* hShowerSmallTracksCount = new TH1D("hShowerSmallTracksCount", "hShowerSmallTracksCount", 20, 0, 20);
    TH1D* hNoShowerSmallTracksCount = new TH1D("hNoShowerSmallTracksCount", "hNoShowerSmallTracksCount", 20, 0, 20);

    double TRACK_THRESHOLD = 15;

    // Loop over events
    Int_t NumEntries = (Int_t) tree->GetEntries();
    for (Int_t i = 0; i < NumEntries; ++i) {
        tree->GetEntry(i);

        bool isShower = false;
        if ((numEmmitedElectrons + numEmmitedPhotons) > 1) isShower = true;

        int numTrks = recoLength->size();
        int smallTrackCount = 0;
        for (int j = 0; j < numTrks; ++j) {
            if (isShower) hShowerRecoLengths->Fill(recoLength->at(j));
            if (!isShower) hNoShowerRecoLengths->Fill(recoLength->at(j));
            if (recoLength->at(j) < TRACK_THRESHOLD) smallTrackCount++;
        }

        if (isShower) hShowerSmallTracksCount->Fill(smallTrackCount);
        if (!isShower) hNoShowerSmallTracksCount->Fill(smallTrackCount);
    }

    // Setup for drawing plots
    std::vector<int> Colors = {
        kBlack, kBlue, kRed, kGreen
    };

    std::vector<std::vector<TH1*>> PlotGroups = {
        {hShowerRecoLengths, hNoShowerRecoLengths},
        {hShowerSmallTracksCount, hNoShowerSmallTracksCount}
    };

    std::vector<std::vector<TString>> PlotLabelGroups = {
        {"Shower", "No shower"}, 
        {"Shower", "No shower"}
    };

    std::vector<TString> PlotTitles = {
        "TrackLenghts",
        "SmallTrackCounts"
    };
    
    std::vector<TString> XLabels = {
        "Track length (cm)",
        "# of small tracks",
    };

    int numPlots = PlotGroups.size();
    for (int iPlot = 0; iPlot < numPlots; ++iPlot) {
        // Set up canvas
        TCanvas* PlotCanvas = new TCanvas("Canvas","Canvas",205,34,1024,768);
        PlotCanvas->cd();
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetBottomMargin(0.15);

        TLegend* leg = new TLegend(0.65,0.65,0.85,0.75);
        leg->SetBorderSize(0);
        leg->SetTextSize(TextSize * 0.8);
        leg->SetTextFont(FontStyle);

        // Get histograms and labels
        std::vector<TH1*> Plots = PlotGroups.at(iPlot);
        std::vector<TString> Labels = PlotLabelGroups.at(iPlot);

        Plots[0]->SetTitle(PlotTitles.at(iPlot));

        Plots[0]->GetXaxis()->SetTitleFont(FontStyle);
        Plots[0]->GetXaxis()->SetLabelFont(FontStyle);
        Plots[0]->GetXaxis()->SetNdivisions(8);
        Plots[0]->GetXaxis()->SetLabelSize(TextSize);
        Plots[0]->GetXaxis()->SetTitleSize(TextSize);
        Plots[0]->GetXaxis()->SetTitle(XLabels.at(iPlot));
        Plots[0]->GetXaxis()->SetTitleOffset(1.1);
        Plots[0]->GetXaxis()->CenterTitle();

        Plots[0]->GetYaxis()->SetTitleFont(FontStyle);
        Plots[0]->GetYaxis()->SetLabelFont(FontStyle);
        Plots[0]->GetYaxis()->SetNdivisions(6);
        Plots[0]->GetYaxis()->SetLabelSize(TextSize);
        Plots[0]->GetYaxis()->SetTitleSize(TextSize);
        Plots[0]->GetYaxis()->SetTitle("Event counts");
        Plots[0]->GetYaxis()->SetTitleOffset(1.1);
        Plots[0]->GetYaxis()->CenterTitle();	

        for (int iSubPlot = 0; iSubPlot < Plots.size(); ++iSubPlot) {
            double iMax = TMath::Max(Plots[iSubPlot]->GetMaximum(), Plots[0]->GetMaximum());
            Plots[iSubPlot]->GetYaxis()->SetRangeUser(0., 1.1 * iMax);
            Plots[0]->GetYaxis()->SetRangeUser(0., 1.1 * iMax);

            leg->AddEntry(Plots[iSubPlot], Labels[iSubPlot], "l");
            Plots[iSubPlot]->SetLineWidth(2);
            Plots[iSubPlot]->SetLineColor(Colors.at(iSubPlot));
            Plots[iSubPlot]->Draw("hist same");
        }

        leg->Draw();
        PlotCanvas->SaveAs(SaveDir + PlotTitles.at(iPlot) + ".png");
        delete PlotCanvas;
    }
}