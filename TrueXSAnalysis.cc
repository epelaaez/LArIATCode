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

void TrueXSAnalysis() {
    // Set defaults
    gStyle->SetOptStat(0);
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    gStyle->SetPalette(kRainBow);

    int FontStyle = 132;
    double TextSize = 0.05;
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis/figs/TrueXSAnalysis/";

    // Load root file
    TString RootFilePath = "/exp/lariat/app/users/epelaez/files/TrueXSPionAbs_histo.root";
    std::unique_ptr<TFile> File(TFile::Open(RootFilePath));
    TDirectory* Directory = (TDirectory*) File->Get("TrueXSPionAbs");

    /////////////////////////
    // Retrieve histograms //
    /////////////////////////

    // Kinetic energies
    TH1D* hIncidentKE  = (TH1D*) Directory->Get("hIncidentKE");
    TH1D* hPionAbsKE   = (TH1D*) Directory->Get("hInteractingKEPionAbs");
    TH1D* hPionAbs0pKE = (TH1D*) Directory->Get("hInteractingKEPionAbs0p");
    TH1D* hPionAbsNpKE = (TH1D*) Directory->Get("hInteractingKEPionAbsNp");

    // All hadronic pion cross section
    TH1D* hCrossSection          = (TH1D*) Directory->Get("hCrossSection");

    // Pion absorption cross sections
    TH1D* hCrossSectionPionAbs   = (TH1D*) Directory->Get("hCrossSectionPionAbs");
    TH1D* hCrossSectionPionAbs0p = (TH1D*) Directory->Get("hCrossSectionPionAbs0p");
    TH1D* hCrossSectionPionAbsNp = (TH1D*) Directory->Get("hCrossSectionPionAbsNp");

    // Background cross sections
    TH1D* hCrossSectionPionInelastic       = (TH1D*) Directory->Get("hCrossSectionPionInelastic");
    TH1D* hCrossSectionPionElastic         = (TH1D*) Directory->Get("hCrossSectionPionElastic");
    TH1D* hCrossSectionChargeExchange      = (TH1D*) Directory->Get("hCrossSectionChargeExchange");
    TH1D* hCrossSectionDoubleChargeExchange = (TH1D*) Directory->Get("hCrossSectionDoubleChargeExchange");
    TH1D* hCrossSectionCaptureAtRest       = (TH1D*) Directory->Get("hCrossSectionCaptureAtRest");
    TH1D* hCrossSectionDecay               = (TH1D*) Directory->Get("hCrossSectionDecay");
    TH1D* hCrossSectionOther               = (TH1D*) Directory->Get("hCrossSectionOther");

    std::vector<TH1D*> Histograms = {
        hCrossSection, 
        hCrossSectionPionAbs, 
        hCrossSectionPionAbs0p, 
        hCrossSectionPionAbsNp,
        hCrossSectionPionInelastic,
        hCrossSectionPionElastic,
        hCrossSectionChargeExchange,
        hCrossSectionDoubleChargeExchange,
        hCrossSectionCaptureAtRest,
        hCrossSectionDecay,
        hCrossSectionOther,
        hIncidentKE,
        hPionAbsKE,
        hPionAbs0pKE,
        hPionAbsNpKE
    };

    ///////////////////////////
    // Adjust axes for plots //
    ///////////////////////////

    for (size_t i = 0; i < Histograms.size(); ++i) {
        if (!Histograms[i]) continue;

        // Rebin and modify range
        Histograms[i]->Rebin(1); // NOT SURE THIS WOULD BE RIGHT, BETTER TO ADJUST BINNING IN MODULE
        Histograms[i]->GetXaxis()->SetRangeUser(0, 600);
        double maxY = Histograms[i]->GetMaximum();
        for (int bin = 1; bin <= Histograms[i]->GetNbinsX(); ++bin) {
            double y = Histograms[i]->GetBinContent(bin) + Histograms[i]->GetBinError(bin);
            if (y > maxY) maxY = y;
        }
        Histograms[i]->GetYaxis()->SetRangeUser(0, 1.1 * maxY);
    }

    //////////////////////////////////////////
    // Plot cross-section histograms nicely //
    //////////////////////////////////////////
    
    TCanvas* c = new TCanvas("Canvas", "Canvas", 205, 34, 1024, 768);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);

    for (size_t i = 0; i < Histograms.size(); ++i) {
        if (!Histograms[i]) continue;

        Histograms[i]->SetLineColor(kBlack);
        Histograms[i]->SetMarkerColor(kBlack);
        Histograms[i]->SetMarkerStyle(20);
        Histograms[i]->SetMarkerSize(1.2);

        Histograms[i]->SetFillColor(kBlack);
        Histograms[i]->SetFillStyle(3001);        // e.g. 3001 = hatched; 1001 = solid
        Histograms[i]->SetFillColorAlpha(kBlack, 0.2);

        Histograms[i]->SetTitle(Histograms[i]->GetName());
        Histograms[i]->GetXaxis()->SetTitle("Kinetic Energy [MeV]");
        Histograms[i]->GetYaxis()->SetTitle("Cross Section [barn]");
        // Histograms[i]->Draw("E1");   // only points
        Histograms[i]->Draw("E1 H"); // points with histogram bars

        c->SetLeftMargin(0.13);
        c->SetBottomMargin(0.13);

        TString histName = Histograms[i]->GetName();
        histName.Remove(0, 1); // Remove the first character
        c->SaveAs(SaveDir + histName + ".png");
    }

    ///////////////////////////////////////
    // Pion absorption stacked histogram //
    ///////////////////////////////////////

    // Create a stacked histogram for pion absorption 0p and Np
    THStack* hStackPionAbs = new THStack("hStackPionAbs", "Pion Absorption Cross Section;Kinetic Energy [MeV];Cross Section [barn]");

    // Set colors for stack
    hCrossSectionPionAbs0p->SetFillColor(kAzure+1);
    hCrossSectionPionAbs0p->SetLineColor(kAzure+1);

    hCrossSectionPionAbsNp->SetFillColor(kOrange+7);
    hCrossSectionPionAbsNp->SetLineColor(kOrange+7);

    // Add to stack
    hStackPionAbs->Add(hCrossSectionPionAbs0p, "H");
    hStackPionAbs->Add(hCrossSectionPionAbsNp, "H");

    hStackPionAbs->Draw("hist");
    hStackPionAbs->SetMaximum(1.1 * std::max(hCrossSectionPionAbs->GetMaximum(), hStackPionAbs->GetMaximum()));
    hStackPionAbs->GetXaxis()->SetRangeUser(0, 600);
    hStackPionAbs->GetXaxis()->SetTitle("Kinetic Energy [MeV]");
    hStackPionAbs->GetYaxis()->SetTitle("Cross Section [barn]");

    // Draw total pion absorption cross-section on top
    hCrossSectionPionAbs->SetLineColor(kBlack);
    hCrossSectionPionAbs->SetLineWidth(2);
    hCrossSectionPionAbs->SetMarkerSize(1.2);
    hCrossSectionPionAbs->SetMarkerStyle(20);
    hCrossSectionPionAbs->SetMarkerColor(kBlack);
    hCrossSectionPionAbs->Draw("E1 SAME");

    // Add legend
    TLegend* leg = new TLegend(0.65, 0.65, 0.85, 0.85);
    leg->SetTextFont(FontStyle);
    leg->SetTextSize(TextSize * 0.9);
    leg->AddEntry(hCrossSectionPionAbs0p, "Abs. 0p", "f");
    leg->AddEntry(hCrossSectionPionAbsNp, "Abs. Np", "f");
    leg->AddEntry(hCrossSectionPionAbs, "Total abs.", "lep");
    leg->Draw();

    c->SaveAs(SaveDir + "PionAbsorptionStacked.png");

    /////////////////////////////////////////////////
    // All hadronic interactions stacked histogram //
    /////////////////////////////////////////////////

    THStack* hStackAll = new THStack("hStackAll", "All Interactions Cross Section;Kinetic Energy [MeV];Cross Section [barn]");

    // Set colors for each process
    hCrossSectionPionAbs0p->SetFillColor(kAzure+1);
    hCrossSectionPionAbs0p->SetLineColor(kAzure+1);
    hStackAll->Add(hCrossSectionPionAbs0p, "H");

    hCrossSectionPionAbsNp->SetFillColor(kOrange+7);
    hCrossSectionPionAbsNp->SetLineColor(kOrange+7);
    hStackAll->Add(hCrossSectionPionAbsNp, "H");

    hCrossSectionPionInelastic->SetFillColor(kGreen+2);
    hCrossSectionPionInelastic->SetLineColor(kGreen+2);
    hStackAll->Add(hCrossSectionPionInelastic, "H");

    hCrossSectionPionElastic->SetFillColor(kBlue+2);
    hCrossSectionPionElastic->SetLineColor(kBlue+2);
    hStackAll->Add(hCrossSectionPionElastic, "H");

    hCrossSectionChargeExchange->SetFillColor(kMagenta+1);
    hCrossSectionChargeExchange->SetLineColor(kMagenta+1);
    hStackAll->Add(hCrossSectionChargeExchange, "H");

    hCrossSectionDoubleChargeExchange->SetFillColor(kCyan+2);
    hCrossSectionDoubleChargeExchange->SetLineColor(kCyan+2);
    hStackAll->Add(hCrossSectionDoubleChargeExchange, "H");

    hCrossSectionCaptureAtRest->SetFillColor(kRed+1);
    hCrossSectionCaptureAtRest->SetLineColor(kRed+1);
    hStackAll->Add(hCrossSectionCaptureAtRest, "H");

    hCrossSectionDecay->SetFillColor(kYellow+1);
    hCrossSectionDecay->SetLineColor(kYellow+1);
    hStackAll->Add(hCrossSectionDecay, "H");

    hCrossSectionOther->SetFillColor(kGray+2);
    hCrossSectionOther->SetLineColor(kGray+2);
    hStackAll->Add(hCrossSectionOther, "H");
    
    // Draw the stack
    hStackAll->Draw("hist");
    hStackAll->SetMaximum(1.1 * std::max(hCrossSection->GetMaximum(), hStackAll->GetMaximum()));
    hStackAll->GetXaxis()->SetRangeUser(0, 600);
    hStackAll->GetXaxis()->SetTitle("Kinetic Energy [MeV]");
    hStackAll->GetYaxis()->SetTitle("Cross Section [barn]");

    // Draw total cross-section on top
    hCrossSection->SetLineColor(kBlack);
    hCrossSection->SetLineWidth(2);
    hCrossSection->SetMarkerSize(1.2);
    hCrossSection->SetMarkerStyle(20);
    hCrossSection->SetMarkerColor(kBlack);
    hCrossSection->Draw("E1 SAME");

    // Add legend
    TLegend* legAll = new TLegend(0.65, 0.40, 0.85, 0.85);
    legAll->SetTextFont(FontStyle);
    legAll->SetTextSize(TextSize * 0.5);
    legAll->AddEntry(hCrossSectionPionAbs0p, "Abs. 0p", "f");
    legAll->AddEntry(hCrossSectionPionAbsNp, "Abs. Np", "f");
    legAll->AddEntry(hCrossSectionPionInelastic, "Inelastic", "f");
    legAll->AddEntry(hCrossSectionPionElastic, "Elastic", "f");
    legAll->AddEntry(hCrossSectionChargeExchange, "Charge exc.", "f");
    legAll->AddEntry(hCrossSectionDoubleChargeExchange, "Double charge exc.", "f");
    legAll->AddEntry(hCrossSectionCaptureAtRest, "Cap. at rest", "f");
    legAll->AddEntry(hCrossSectionDecay, "Decay", "f");
    legAll->AddEntry(hCrossSectionOther, "Other", "f");
    legAll->AddEntry(hCrossSection, "Total", "lep");
    legAll->Draw();

    c->SaveAs(SaveDir + "AllInteractionsStacked.png");

    /////////////////////////////////////////////////////
    // Total cross section compared to pion absorption //
    /////////////////////////////////////////////////////

    THStack* hStackTotalAbs = new THStack("hStackTotalAbs", "Total vs. Absorption Cross Section;Kinetic Energy [MeV];Cross Section [barn]");

    hStackTotalAbs->Add(hCrossSectionPionAbs0p, "H");
    hStackTotalAbs->Add(hCrossSectionPionAbsNp, "H");

    // Draw the stack
    hStackTotalAbs->Draw("hist");
    hStackTotalAbs->SetMaximum(1.1 * std::max(hCrossSection->GetMaximum(), hStackTotalAbs->GetMaximum()));
    hStackTotalAbs->GetXaxis()->SetRangeUser(0, 600);
    hStackTotalAbs->GetXaxis()->SetTitle("Kinetic Energy [MeV]");
    hStackTotalAbs->GetYaxis()->SetTitle("Cross Section [barn]");

    // Draw total cross-section on top
    hCrossSection->SetLineColor(kBlack);
    hCrossSection->SetLineWidth(2);
    hCrossSection->SetMarkerSize(1.2);
    hCrossSection->SetMarkerStyle(20);
    hCrossSection->SetMarkerColor(kBlack);
    hCrossSection->Draw("E1 SAME");

    // Add legend
    TLegend* legTotalAbs = new TLegend(0.65, 0.65, 0.85, 0.85);
    legTotalAbs->SetTextFont(FontStyle);
    legTotalAbs->SetTextSize(TextSize * 0.9);
    legTotalAbs->AddEntry(hCrossSectionPionAbs0p, "Abs. 0p", "f");
    legTotalAbs->AddEntry(hCrossSectionPionAbsNp, "Abs. Np", "f");
    legTotalAbs->AddEntry(hCrossSection, "Total", "lep");
    legTotalAbs->Draw();

    c->SaveAs(SaveDir + "TotalComparedAbsStacked.png");
}