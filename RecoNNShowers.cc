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

//////////////////////
// Global constants //
//////////////////////

std::map<int, std::string> backgroundTypes = {
    {-1, "Not flagged as background"},
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
    {12, "Elastic scattering"}
};

double SMALL_TRACK_LENGTH = 35;

//////////////////////
// Helper functions //
//////////////////////

void printOneDPlots(
    const TString& dir, 
    int fontStyle, 
    double textSize,
    std::vector<std::vector<TH1*>>& groups,
    std::vector<int>& colors, 
    std::vector<std::vector<TString>>& labels, 
    std::vector<TString>& titles, 
    std::vector<TString>& xlabels, 
    std::vector<TString>& ylabels,
    std::vector<bool>& stack
);

///////////////////
// Main function //
///////////////////

void RecoNNShowers() {
    // Set defaults
    gStyle->SetOptStat(0); // get rid of stats box
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    gStyle->SetPalette(kRainBow);

    int FontStyle = 132;
    double TextSize = 0.06;
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis/figs/RecoNNShowers/";

    // Load NN root file
    TString RootFilePath = "/exp/lariat/app/users/epelaez/files/NNBeamlineShower_hist.root";
    std::unique_ptr<TFile> File(TFile::Open(RootFilePath));
    TDirectory* Directory = (TDirectory*)File->Get("NNBeamlineShower");

    ///////////////////
    // Load branches //
    ///////////////////

    // Load tree and branches
    TTree* tree = (TTree*) Directory->Get<TTree>("NNShowerTree");

    // Signal information
    bool isPionAbsorptionSignal; int backgroundType; 
    tree->SetBranchAddress("isPionAbsorptionSignal", &isPionAbsorptionSignal);
    tree->SetBranchAddress("backgroundType", &backgroundType);

    // Shower probability information
    double trackProb, showerProb, showerMultProb;
    bool obtainedProbabilities;
    tree->SetBranchAddress("trackProb", &trackProb);
    tree->SetBranchAddress("showerProb", &showerProb);
    tree->SetBranchAddress("showerMultProb", &showerMultProb);
    tree->SetBranchAddress("obtainedProbabilities", &obtainedProbabilities);

    ///////////////////////
    // Create histograms //
    ///////////////////////

    TH1D* hElectronShowerProb = new TH1D("hElectronShowerProb", "hElectronShowerProb;;", 20, 0., 1.);
    TH1D* hPionShowerProb     = new TH1D("hPionShowerProb", "hPionShowerProb;;", 20, 0., 1.);
    TH1D* hMuonShowerProb     = new TH1D("hMuonShowerProb", "hMuonShowerProb;;", 20, 0., 1.);

    TH1D* hPionChargeExchangeShowerProb = new TH1D("hPionChargeExchangeShowerProb", "hPionChargeExchangeShowerProb;;", 20, 0., 1.);
    TH1D* hPionOtherShowerProb          = new TH1D("hPionOtherShowerProb", "hPionOtherShowerProb;;", 20, 0., 1.);

    //////////////////////
    // Loop over events //
    //////////////////////

    Int_t NumEntries = (Int_t) tree->GetEntries();
    std::cout << "Num entries: " << NumEntries << std::endl;

    for (Int_t i = 0; i < NumEntries; ++i) {
        tree->GetEntry(i);

        if (!obtainedProbabilities) continue;

        if (backgroundType == 2) {
            hMuonShowerProb->Fill(showerProb);
        } else if (backgroundType == 3) {
            hElectronShowerProb->Fill(showerProb);
        } else {
            hPionShowerProb->Fill(showerProb);
        }

        if (backgroundType == 7) {
            hPionChargeExchangeShowerProb->Fill(showerProb);
        } else if (backgroundType != 2 && backgroundType != 3) {
            hPionOtherShowerProb->Fill(showerProb);
        }
    }

    /////////////////
    // Print plots //
    /////////////////

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
        {hElectronShowerProb, hPionShowerProb, hMuonShowerProb},
        {hPionChargeExchangeShowerProb, hPionOtherShowerProb}
    };

    std::vector<std::vector<TString>> PlotLabelGroups = {
        {"Electron", "Pion", "Muon"},
        {"Ch. exch.", "Other"}
    };

    std::vector<TString> PlotTitles = {
        "PrimaryTrackShowerProb",
        "PionInteractionShowerProb"
    };

    std::vector<TString> XLabels = {
        "Shower probability",
        "Shower probability"
    };

    std::vector<TString> YLabels = {
        "Number of events",
        "Number of events"
    };

    std::vector<bool> PlotStacked = {
        false,
        false
    };

    printOneDPlots(
        SaveDir, FontStyle, TextSize,
        PlotGroups,
        Colors,
        PlotLabelGroups,
        PlotTitles,
        XLabels,
        YLabels,
        PlotStacked
    );
}

void printOneDPlots(
    const TString& dir, 
    int fontStyle, 
    double textSize,
    std::vector<std::vector<TH1*>>& groups,
    std::vector<int>& colors, 
    std::vector<std::vector<TString>>& labels, 
    std::vector<TString>& titles, 
    std::vector<TString>& xlabels, 
    std::vector<TString>& ylabels,
    std::vector<bool>& stack
) {
    int numPlots = groups.size();
    for (int iPlot = 0; iPlot < numPlots; ++iPlot) {
        // Set up canvas
        TCanvas* PlotCanvas = new TCanvas("Canvas", "Canvas", 205, 34, 1300, 768);

        TPad* mainPad = new TPad("mainPad","",0.0, 0.0, 0.80, 1.0);
        mainPad->SetRightMargin(0.05);
        mainPad->SetLeftMargin(0.15);
        mainPad->SetBottomMargin(0.15);
        mainPad->Draw();
        mainPad->cd();

        TLegend* leg = new TLegend(0.0, 0.0, 1.0, 1.0);
        leg->SetTextSize(textSize * 2.5);
        leg->SetTextFont(fontStyle);

        TPad* legendPad = new TPad("legendPad", "", 0.80, 0.4, 1.0, 0.8);
        legendPad->SetLeftMargin(0.05);
        legendPad->SetRightMargin(0.05);
        legendPad->SetBottomMargin(0.15);
        legendPad->SetFillStyle(0);
        legendPad->SetBorderMode(0);

        // Get histograms and labels
        std::vector<TH1*> Plots = groups.at(iPlot);
        std::vector<TString> Labels = labels.at(iPlot);	

        // If stacked
        if (stack.at(iPlot)) {
            THStack stack("stack", titles.at(iPlot));

            // Style and add each histogram to the stack
            for (int iSubPlot = 0; iSubPlot < (int) Plots.size(); ++iSubPlot) {
                TH1* h = Plots[iSubPlot];
                leg->AddEntry(h, Labels[iSubPlot], "f");
                h->SetLineWidth(2);
                h->SetLineColor(colors.at(iSubPlot));
                h->SetFillColor(colors.at(iSubPlot));
                h->SetFillColorAlpha(colors.at(iSubPlot), 0.2);
                h->SetFillStyle(3001);
                stack.Add(h, "H");
            }

            // Style the stack
            stack.SetTitle(titles.at(iPlot));

            // Draw stack
            stack.Draw("hist");

            stack.GetXaxis()->SetTitleFont(fontStyle);
            stack.GetXaxis()->SetLabelFont(fontStyle);
            stack.GetXaxis()->SetNdivisions(8);
            stack.GetXaxis()->SetLabelSize(textSize);
            stack.GetXaxis()->SetTitleSize(textSize);
            stack.GetXaxis()->SetTitle(xlabels.at(iPlot));
            stack.GetXaxis()->SetTitleOffset(1.1);
            stack.GetXaxis()->CenterTitle();

            stack.GetYaxis()->SetTitleFont(fontStyle);
            stack.GetYaxis()->SetLabelFont(fontStyle);
            stack.GetYaxis()->SetNdivisions(6);
            stack.GetYaxis()->SetLabelSize(textSize);
            stack.GetYaxis()->SetTitleSize(textSize);
            stack.GetYaxis()->SetTitle(ylabels.at(iPlot));
            stack.GetYaxis()->SetTitleOffset(1.1);
            stack.GetYaxis()->CenterTitle();

            // Determine proper max from stack object
            double stackMax = stack.GetMaximum();
            double YAxisRange = 1.15 * stackMax;
            stack.SetMaximum(YAxisRange);
            stack.SetMinimum(0.0);

            gPad->Update();
            TAxis *x = Plots[0]->GetXaxis();
            x->SetMaxDigits(3);
            gPad->Modified();
            gPad->Update();

            PlotCanvas->cd();
            legendPad->Draw();
            legendPad->cd();
            leg->Draw();

            PlotCanvas->SaveAs(dir + titles.at(iPlot) + ".png");
        } 
        // If not stacked
        else {
            // Style the first plot
            Plots[0]->SetTitle(titles.at(iPlot));

            Plots[0]->GetXaxis()->SetTitleFont(fontStyle);
            Plots[0]->GetXaxis()->SetLabelFont(fontStyle);
            Plots[0]->GetXaxis()->SetNdivisions(8);
            Plots[0]->GetXaxis()->SetLabelSize(textSize);
            Plots[0]->GetXaxis()->SetTitleSize(textSize);
            Plots[0]->GetXaxis()->SetTitle(xlabels.at(iPlot));
            Plots[0]->GetXaxis()->SetTitleOffset(1.1);
            Plots[0]->GetXaxis()->CenterTitle();

            Plots[0]->GetYaxis()->SetTitleFont(fontStyle);
            Plots[0]->GetYaxis()->SetLabelFont(fontStyle);
            Plots[0]->GetYaxis()->SetNdivisions(6);
            Plots[0]->GetYaxis()->SetLabelSize(textSize);
            Plots[0]->GetYaxis()->SetTitleSize(textSize);
            Plots[0]->GetYaxis()->SetTitle(ylabels.at(iPlot));
            Plots[0]->GetYaxis()->SetTitleOffset(1.1);
            Plots[0]->GetYaxis()->CenterTitle();

            double imax;
            for (int iSubPlot = 0; iSubPlot < (int) Plots.size(); ++iSubPlot) {
                leg->AddEntry(Plots[iSubPlot], Labels[iSubPlot], "f");
                Plots[iSubPlot]->SetLineWidth(2);
                Plots[iSubPlot]->SetLineColor(colors.at(iSubPlot));
                Plots[iSubPlot]->Draw("hist same");

                imax = TMath::Max(Plots[iSubPlot]->GetMaximum(), Plots[0]->GetMaximum());

                double YAxisRange = 1.15 * imax;
                Plots[iSubPlot]->GetYaxis()->SetRangeUser(0., YAxisRange);
                Plots[0]->GetYaxis()->SetRangeUser(0., YAxisRange);	
            }

            gPad->Update();
            TAxis *x = Plots[0]->GetXaxis();
            x->SetMaxDigits(3);
            gPad->Modified();
            gPad->Update();

            PlotCanvas->cd();
            legendPad->Draw();
            legendPad->cd();
            leg->Draw();

            PlotCanvas->SaveAs(dir + titles.at(iPlot) + ".png");
        }
        delete PlotCanvas;
    }
}