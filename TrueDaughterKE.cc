#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"

#include <vector>
#include <string>
#include <iostream>

// True KE of proton and pi+- daughters of the primary pi- at generator (Geant4) level, using the full MC sample
// with no selection, and compared to the daughters in events reaching the secondary PID stage (from SecondaryStudies_MC.root).
// Used to check whether the dip near 200 MeV seen in the selected sample is present in the underlying true spectrum.
void TrueDaughterKE(
    TString SelPath = "/exp/lariat/app/users/epelaez/analysis_matt/files/SecondaryStudies_MC.root",
    TString OutPath = "/exp/lariat/app/users/epelaez/analysis_matt/files/TrueDaughterKE.root",
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis_matt/figs/SecondaryStudies/"
) {
    gStyle->SetOptStat(0);
    TH1::AddDirectory(false);

    int FontStyle = 132;
    double TextSize = 0.05;

    ///////////////////////////////////////
    // Fill true daughter KE (all events) //
    ///////////////////////////////////////

    TChain* Chain = new TChain("anatree/anatree");
    Chain->Add("/exp/lariat/app/users/epelaez/files/anatree_60a/chunks/*.root");
    Chain->SetBranchStatus("*", 0);
    for (const char* b : {"geant_list_size", "pdg", "Mass", "Eng", "TrackId", "Mother", "process_primary", "EndPointz"}) Chain->SetBranchStatus(b, 1);

    const int kMax = 20000;
    int geant_list_size;
    static int pdg[kMax], TrackId[kMax], Mother[kMax], process_primary[kMax];
    static float Mass[kMax], Eng[kMax], EndPointz[kMax];
    Chain->SetBranchAddress("geant_list_size", &geant_list_size);
    Chain->SetBranchAddress("pdg", pdg);
    Chain->SetBranchAddress("TrackId", TrackId);
    Chain->SetBranchAddress("Mother", Mother);
    Chain->SetBranchAddress("process_primary", process_primary);
    Chain->SetBranchAddress("Mass", Mass);
    Chain->SetBranchAddress("Eng", Eng);
    Chain->SetBranchAddress("EndPointz", EndPointz);

    // Fine binning for the standalone spectra, and the same binning as hEffProtonAll / hEffPionAll for the comparison
    TH1D* hProtonFine = new TH1D("hTrueKEProtonFine", ";True proton KE [MeV];Daughters", 80, 0, 400);
    TH1D* hPionFine   = new TH1D("hTrueKEPionFine", ";True #pi^{#pm} KE [MeV];Daughters", 80, 0, 400);
    TH1D* hProton     = new TH1D("hTrueKEProton", ";True proton KE [MeV];", 60, 0, 300);
    TH1D* hPion       = new TH1D("hTrueKEPion", ";True #pi^{#pm} KE [MeV];", 30, 0, 300);

    Long64_t nEntries = Chain->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        Chain->GetEntry(i);
        if (i % 100000 == 0) std::cout << "Processed " << i << " events" << std::endl;

        int primary = -1;
        for (int j = 0; j < geant_list_size; ++j) {
            if (process_primary[j] == 1 && pdg[j] == -211) { primary = j; break; }
        }
        // Primary must end inside the TPC
        if (primary < 0 || EndPointz[primary] < 0 || EndPointz[primary] > 90) continue;

        for (int j = 0; j < geant_list_size; ++j) {
            if (Mother[j] != TrackId[primary]) continue;
            double ke = 1000 * (Eng[j] - Mass[j]);
            if (pdg[j] == 2212) {
                hProtonFine->Fill(ke); hProton->Fill(ke);
            } else if (std::abs(pdg[j]) == 211) {
                hPionFine->Fill(ke); hPion->Fill(ke);
            }
        }
    }

    TFile out(OutPath, "RECREATE");
    for (TH1D* h : {hProtonFine, hPionFine, hProton, hPion}) h->Write();
    out.Close();

    ///////////
    // Plots //
    ///////////

    const int colBlue   = TColor::GetColor(0, 114, 178);
    const int colOrange = TColor::GetColor(213, 94, 0);

    auto styleAxes = [&](TH1* h, double scale = 1.0) {
        for (TAxis* ax : {h->GetXaxis(), h->GetYaxis()}) {
            ax->SetTitleFont(FontStyle); ax->SetLabelFont(FontStyle);
            ax->SetTitleSize(TextSize * scale); ax->SetLabelSize(TextSize * 0.85 * scale);
            ax->CenterTitle();
        }
    };

    // Full-sample spectra on their own
    std::vector<std::pair<TH1D*, TString>> fine = {{hProtonFine, "Proton"}, {hPionFine, "Pion"}};
    for (auto& [h, name] : fine) {
        TCanvas* c = new TCanvas("cFine" + name, "", 1300, 900);
        c->SetLeftMargin(0.14); c->SetRightMargin(0.05); c->SetBottomMargin(0.14); c->SetGridy();
        h->SetTitle((name == "Proton" ? "Proton" : "#pi^{#pm}") + TString(" daughters of the primary #pi^{-}, all MC events (no selection)"));
        styleAxes(h);
        h->GetYaxis()->SetTitleOffset(1.3);
        h->SetLineColor(kBlack); h->SetMarkerColor(kBlack); h->SetMarkerStyle(20); h->SetMarkerSize(0.9);
        h->SetMinimum(0);
        h->Draw("E1P");
        c->SaveAs(SaveDir + "TrueKEAllEvents" + name + ".png");
        delete c;
    }

    // Selected sample (events at secondary PID, weighted) vs all events, both area normalized, with ratio
    std::unique_ptr<TFile> fSel(TFile::Open(SelPath, "READ"));
    std::vector<std::tuple<TH1D*, TString, TString>> comps = {{hProton, "hEffProtonAll", "Proton"}, {hPion, "hEffPionAll", "Pion"}};
    for (auto& [hAll, selName, name] : comps) {
        TH1D* hSel = (TH1D*) fSel->Get(selName);
        if (!hSel) continue;
        hSel = (TH1D*) hSel->Clone("hSel" + name);
        TH1D* hAllN = (TH1D*) hAll->Clone("hAllNorm" + name);
        hAllN->Scale(1. / hAllN->Integral());
        hSel->Scale(1. / hSel->Integral());

        TCanvas* c = new TCanvas("cComp" + name, "", 1300, 1100);
        TPad* pTop = new TPad("pTop" + name, "", 0, 0.30, 1, 1);
        TPad* pBot = new TPad("pBot" + name, "", 0, 0, 1, 0.30);
        pTop->SetLeftMargin(0.14); pTop->SetRightMargin(0.05); pTop->SetBottomMargin(0.02); pTop->SetGridy();
        pBot->SetLeftMargin(0.14); pBot->SetRightMargin(0.05); pBot->SetBottomMargin(0.35); pBot->SetTopMargin(0.05); pBot->SetGridy();
        pTop->Draw(); pBot->Draw();

        pTop->cd();
        hAllN->SetTitle((name == "Proton" ? "Proton" : "#pi^{#pm}") + TString(" daughters of the primary #pi^{-} (area normalized)"));
        hAllN->GetYaxis()->SetTitle("Fraction of daughters");
        styleAxes(hAllN);
        hAllN->GetXaxis()->SetLabelSize(0); hAllN->GetXaxis()->SetTitleSize(0);
        hAllN->SetLineColor(colBlue); hAllN->SetLineWidth(3);
        hAllN->SetMinimum(0);
        hAllN->SetMaximum(1.15 * std::max(hAllN->GetMaximum(), hSel->GetMaximum()));
        hAllN->Draw("HIST E");
        hSel->SetLineColor(colOrange); hSel->SetMarkerColor(colOrange); hSel->SetMarkerStyle(20);
        hSel->Draw("E1P SAME");
        TLegend* leg = new TLegend(0.45, 0.70, 0.93, 0.88);
        leg->SetTextFont(FontStyle); leg->SetTextSize(TextSize * 0.8); leg->SetBorderSize(0);
        leg->AddEntry(hAllN, "All MC events (no selection)", "l");
        leg->AddEntry(hSel, "Events at secondary PID stage (weighted)", "lep");
        leg->Draw();

        pBot->cd();
        TH1D* hRatio = (TH1D*) hSel->Clone("hRatio" + name);
        hRatio->Divide(hAllN);
        hRatio->SetTitle("");
        styleAxes(hRatio, 0.7 / 0.3);
        hRatio->GetYaxis()->SetTitle("Selected / all");
        hRatio->GetYaxis()->SetNdivisions(505);
        hRatio->GetYaxis()->SetTitleOffset(0.55);
        hRatio->GetXaxis()->SetTitleOffset(1.0);
        hRatio->GetYaxis()->SetRangeUser(0, 2);
        hRatio->Draw("E1P");

        c->SaveAs(SaveDir + "TrueKESelectedVsAll" + name + ".png");
        delete c;
    }
}
