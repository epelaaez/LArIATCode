#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TPRegexp.h"

#include <vector>

#include "Helpers.cc"

// Draws several ratio curves (efficiencies / fractions) on one plot with a vertical line at a threshold
void DrawFractionCurves(
    const TString& SaveDir,
    const TString& Name,
    const std::vector<TH1D*>& numerators,
    TH1D* denominator,
    const std::vector<TString>& labels,
    const std::vector<int>& colors,
    double threshold,
    const TString& thresholdLabel,
    const TString& PlotTitle,
    const TString& XTitle,
    const TString& YTitle,
    int FontStyle,
    double TextSize,
    int Rebin = 1
) {
    TCanvas* c = new TCanvas("c_" + Name, "", 1300, 900);
    c->SetLeftMargin(0.13); c->SetRightMargin(0.05);
    c->SetBottomMargin(0.14); c->SetTopMargin(0.08);
    c->SetGridy();

    TLegend* leg = new TLegend(0.15, 0.73, 0.93, 0.90);
    leg->SetNColumns(2);
    leg->SetTextFont(FontStyle);
    leg->SetTextSize(TextSize * 0.75);
    leg->SetBorderSize(0);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);

    TH1D* den = (TH1D*) denominator->Clone(Name + "_den");
    den->Rebin(Rebin);

    std::vector<TH1D*> ratios;
    for (size_t i = 0; i < numerators.size(); ++i) {
        TH1D* num = (TH1D*) numerators[i]->Clone(Name + Form("_num%zu", i));
        num->Rebin(Rebin);
        TH1D* r = (TH1D*) num->Clone(Name + Form("_ratio%zu", i));
        r->Divide(num, den, 1.0, 1.0, "B");
        r->SetLineColor(colors[i]);
        r->SetMarkerColor(colors[i]);
        r->SetMarkerStyle(20);
        r->SetMarkerSize(1.0);
        r->SetLineWidth(2);
        ratios.push_back(r);
    }

    TH1D* frame = ratios[0];
    frame->SetTitle(PlotTitle);
    frame->GetYaxis()->SetRangeUser(0.0, 1.45);
    frame->GetXaxis()->SetTitle(XTitle);
    frame->GetYaxis()->SetTitle(YTitle);
    for (TAxis* ax : {frame->GetXaxis(), frame->GetYaxis()}) {
        ax->SetTitleFont(FontStyle); ax->SetLabelFont(FontStyle);
        ax->SetTitleSize(TextSize);  ax->SetLabelSize(TextSize);
        ax->CenterTitle();
    }
    frame->GetXaxis()->SetTitleOffset(1.0);
    frame->GetYaxis()->SetTitleOffset(1.0);

    frame->Draw("E1P");
    for (size_t i = 1; i < ratios.size(); ++i) ratios[i]->Draw("E1P SAME");
    for (size_t i = 0; i < ratios.size(); ++i) leg->AddEntry(ratios[i], labels[i], "lep");

    if (threshold > frame->GetXaxis()->GetXmin()) {
        TLine* line = new TLine(threshold, 0.0, threshold, 1.2);
        line->SetLineStyle(2); line->SetLineWidth(2); line->SetLineColor(kGray + 2);
        line->Draw();
        leg->AddEntry(line, thresholdLabel, "l");
    }
    leg->Draw();

    c->SaveAs(SaveDir + Name + ".png");
    delete c;
}

void SecondaryStudiesPlots(
    TString MCPath   = "/exp/lariat/app/users/epelaez/analysis_matt/files/SecondaryStudies_MC.root",
    TString DataPath = "/exp/lariat/app/users/epelaez/analysis_matt/files/SecondaryStudies_Data.root",
    TString SaveDir  = "/exp/lariat/app/users/epelaez/analysis_matt/figs/SecondaryStudies/"
) {
    gStyle->SetOptStat(0);
    TH1::AddDirectory(false);

    int FontStyle = 132;
    double TextSize = 0.05;

    std::unique_ptr<TFile> fMC(TFile::Open(MCPath, "READ"));
    std::unique_ptr<TFile> fData(TFile::Open(DataPath, "READ"));
    auto getMC   = [&](const TString& name) { return (TH1D*) fMC->Get(name); };
    auto getData = [&](const TString& name) { return (TH1D*) fData->Get(name); };

    // Data / MC comparisons are drawn with coarser bins than they are filled with (data statistics are small)
    auto rebinned = [&](TH1* h, int n) {
        if (!h) return h;
        TH1* r = (TH1*) h->Clone(TString(h->GetName()) + Form("_rebin%d", n));
        r->Rebin(n);
        return r;
    };

    const int colBlue   = TColor::GetColor(0,   114, 178);
    const int colOrange = TColor::GetColor(213, 94,  0);
    const int colGreen  = TColor::GetColor(0,   158, 115);
    const int colGray   = TColor::GetColor(64,  64,  64);

    ////////////////////////////////
    // Fine-binned efficiencies   //
    ////////////////////////////////

    DrawFractionCurves(
        SaveDir, "ProtonEfficiency",
        {getMC("hEffProtonReco"), getMC("hEffProtonTagChi2"), getMC("hEffProtonTagAny")},
        getMC("hEffProtonAll"),
        {"Reconstructed near vertex", "Tagged proton (#chi^{2})", "Tagged proton (#chi^{2} or mean dE/dx)"},
        {colGray, colBlue, colOrange},
        1000 * PROTON_ENERGY_LOWER_BOUND, Form("Signal threshold (%.0f MeV)", 1000 * PROTON_ENERGY_LOWER_BOUND),
        "Secondary proton efficiency (MC)", "True proton KE [MeV]", "Efficiency",
        FontStyle, TextSize, 3
    );

    DrawFractionCurves(
        SaveDir, "PionEfficiency",
        {getMC("hEffPionReco"), getMC("hEffPionTagChi2"), getMC("hEffPionTagAny")},
        getMC("hEffPionAll"),
        {"Reconstructed near vertex", "Tagged pion (#chi^{2})", "Tagged pion (#chi^{2} or mean dE/dx)"},
        {colGray, colBlue, colOrange},
        1000 * PION_SCATTERING_ENERGY_THRESHOLD, Form("Signal threshold (%.0f MeV)", 1000 * PION_SCATTERING_ENERGY_THRESHOLD),
        "Secondary #pi^{-} efficiency (MC, inelastic daughters)", "True #pi^{-} KE [MeV]", "Efficiency",
        FontStyle, TextSize, 2
    );

    ////////////////////////////////////////
    // Pion threshold: reco outcome vs KE //
    ////////////////////////////////////////

    DrawFractionCurves(
        SaveDir, "PionThresholdOutcome",
        {getMC("hDecidingPionKERecoScatter"), getMC("hDecidingPionKERecoNp"), getMC("hDecidingPionKEReco0pCandidate"), getMC("hDecidingPionKERejected")},
        getMC("hDecidingPionKEAll"),
        {"Reco scatter", "Reco abs Np", "Reco abs 0p (before cluster cut)", "Rejected"},
        {colBlue, colOrange, colGreen, colGray},
        1000 * PION_SCATTERING_ENERGY_THRESHOLD, Form("Signal threshold (%.0f MeV)", 1000 * PION_SCATTERING_ENERGY_THRESHOLD),
        "True scatter candidates by reco outcome (MC)", "True outgoing #pi^{-} KE [MeV]", "Fraction",
        FontStyle, TextSize, 2
    );

    // Print the outcome fractions as a table
    {
        TH1D* hAll = getMC("hDecidingPionKEAll");
        std::cout << std::endl << "=== Reco outcome vs true outgoing pion KE ===" << std::endl;
        std::cout << Form("%-12s %10s %8s %8s %8s %8s", "KE [MeV]", "N", "Scat", "Np", "0p", "Rej") << std::endl;
        for (int b = 1; b <= hAll->GetNbinsX(); ++b) {
            double n = hAll->GetBinContent(b);
            if (n <= 0) continue;
            std::cout << Form("[%3.0f, %3.0f]  %10.1f %8.3f %8.3f %8.3f %8.3f",
                hAll->GetXaxis()->GetBinLowEdge(b), hAll->GetXaxis()->GetBinUpEdge(b), n,
                getMC("hDecidingPionKERecoScatter")->GetBinContent(b) / n,
                getMC("hDecidingPionKERecoNp")->GetBinContent(b) / n,
                getMC("hDecidingPionKEReco0pCandidate")->GetBinContent(b) / n,
                getMC("hDecidingPionKERejected")->GetBinContent(b) / n) << std::endl;
        }
    }

    ////////////////////////////////////////////////////
    // Data vs MC for secondaries near the vertex     //
    ////////////////////////////////////////////////////

    // Normalize MC to the number of data events reaching the secondary PID stage
    double scale = getData("hSecStudyNumEvents")->Integral() / getMC("hSecStudyNumEvents")->Integral();
    std::cout << std::endl << "Events at secondary PID: data = " << getData("hSecStudyNumEvents")->Integral()
              << ", MC (weighted) = " << getMC("hSecStudyNumEvents")->Integral() << ", scale = " << scale << std::endl;

    // True protons / pions split by whether they are above or below the signal threshold
    std::vector<TString> trueNames  = {"ProtonAbove", "ProtonBelow", "PionAbove", "PionBelow", "Muon", "EM", "Other"};
    std::vector<TString> trueLabels = {
        Form("True p, KE #geq %.0f MeV", 1000 * PROTON_ENERGY_LOWER_BOUND), Form("True p, KE < %.0f MeV", 1000 * PROTON_ENERGY_LOWER_BOUND),
        Form("True #pi^{#pm}, KE #geq %.0f MeV", 1000 * PION_SCATTERING_ENERGY_THRESHOLD), Form("True #pi^{#pm}, KE < %.0f MeV", 1000 * PION_SCATTERING_ENERGY_THRESHOLD),
        "True #mu", "True e/#gamma", "Other / unmatched"
    };
    std::vector<int> mcColors = {
        colBlue, TColor::GetColor(140, 200, 240), colOrange, TColor::GetColor(245, 180, 130),
        colGreen, TColor::GetColor(200, 185, 0), TColor::GetColor(204, 121, 167)
    };

    // shapeOnly: normalize the MC stack to the data integral (shape comparison)
    auto getMCSplit = [&](const TString& base, bool shapeOnly = false, int rebin = 1) {
        std::vector<TH1*> hs;
        double total = 0;
        for (const TString& name : trueNames) {
            TH1* h = rebinned(getMC(base + name), rebin);
            total += h->Integral();
            hs.push_back(h);
        }
        double s = shapeOnly ? getData(base)->Integral() / total : scale;
        for (TH1* h : hs) h->Scale(s);
        return hs;
    };

    // Rebin factors: KE 10 MeV -> 40 MeV, length 1 cm -> 5 cm, chi^2 0.2 -> 1, cos 0.1 -> 0.4, mean dE/dx 0.25 -> 1 MeV/cm
    const int rebinKE = 4, rebinLength = 5, rebinChi2 = 5, rebinCos = 4, rebinMeanDEDX = 4;

    struct Comparison { TString base, name, title, xlabel; int rebin; };
    std::vector<Comparison> comparisons = {
        {"hSecRecoKE",           "RecoKEAllSecondaries",     "Secondary tracks near vertex",                       "Reco calorimetric KE [MeV]", rebinKE},
        {"hSecRecoKEProtonChi2", "RecoKEProtonTaggedChi2",   "Secondary tracks tagged proton (#chi^{2})",          "Reco calorimetric KE [MeV]", rebinKE},
        {"hSecRecoKEProtonDEDX", "RecoKEProtonTaggedDEDX",   "Secondary tracks tagged proton (mean dE/dx only)",   "Reco calorimetric KE [MeV]", rebinKE},
        {"hSecRecoKEPionChi2",   "RecoKEPionTaggedChi2",     "Secondary tracks tagged pion (#chi^{2})",            "Reco calorimetric KE [MeV]", rebinKE},
        {"hSecRecoKEPionDEDX",   "RecoKEPionTaggedDEDX",     "Secondary tracks tagged pion (mean dE/dx only)",     "Reco calorimetric KE [MeV]", rebinKE},
        {"hSecProtonChi2Value",  "ProtonChi2AllSecondaries", "Secondary tracks near vertex",                       "Proton #chi^{2}",            rebinChi2},
        {"hSecLength",           "LengthAllSecondaries",     "Secondary tracks near vertex",                       "Track length [cm]",          rebinLength},
        {"hSecLengthProtonChi2", "LengthProtonTaggedChi2",   "Secondary tracks tagged proton (#chi^{2})",          "Track length [cm]",          rebinLength},
        {"hSecLengthProtonDEDX", "LengthProtonTaggedDEDX",   "Secondary tracks tagged proton (mean dE/dx only)",   "Track length [cm]",          rebinLength},
        {"hSecLengthPionChi2",   "LengthPionTaggedChi2",     "Secondary tracks tagged pion (#chi^{2})",            "Track length [cm]",          rebinLength},
        {"hSecLengthPionDEDX",   "LengthPionTaggedDEDX",     "Secondary tracks tagged pion (mean dE/dx only)",     "Track length [cm]",          rebinLength},
        {"hSecMeanDEDX",           "MeanDEDXAllSecondaries",  "Secondary tracks near vertex",                         "Mean dE/dx [MeV/cm]", rebinMeanDEDX},
        {"hSecMeanDEDXPionChi2",   "MeanDEDXPionTaggedChi2",  "Secondary tracks tagged pion (#chi^{2})",              "Mean dE/dx [MeV/cm]", rebinMeanDEDX},
        {"hSecMeanDEDXProtonChi2", "MeanDEDXProtonTaggedChi2","Secondary tracks tagged proton (#chi^{2})",            "Mean dE/dx [MeV/cm]", rebinMeanDEDX},
        {"hSecMeanDEDXNoChi2",     "MeanDEDXNoChi2Tag",       "Secondary tracks not tagged by #chi^{2} (mean dE/dx PID)", "Mean dE/dx [MeV/cm]", rebinMeanDEDX},
    };

    // Each comparison twice: MC normalized to the number of events at secondary PID (total), and to the data integral (shape)
    std::cout << std::endl << "=== Tracks per event at secondary PID (data / MC) ===" << std::endl;
    for (const auto& comp : comparisons) {
        if (!getData(comp.base) || !getMC(comp.base + trueNames[0])) continue;

        PrintDataVsMCContribPlot(
            SaveDir, comp.name,
            rebinned(getData(comp.base), comp.rebin),
            getMCSplit(comp.base, false, comp.rebin),
            trueLabels, mcColors,
            comp.title + " (normalized to events)", comp.xlabel, "Tracks",
            FontStyle, TextSize, true, nullptr, true
        );
        PrintDataVsMCContribPlot(
            SaveDir, comp.name + "_Shape",
            rebinned(getData(comp.base), comp.rebin),
            getMCSplit(comp.base, true, comp.rebin),
            trueLabels, mcColors,
            comp.title + " (area normalized)", comp.xlabel, "Tracks",
            FontStyle, TextSize, true, nullptr, true
        );

        double nData = getData(comp.base)->Integral(), nMC = 0;
        for (TH1* h : getMCSplit(comp.base)) nMC += h->Integral();
        std::cout << Form("%-28s data %6.0f (%.3f/evt)   MC %8.1f (%.3f/evt)   data/MC %.2f",
            comp.name.Data(), nData, nData / getData("hSecStudyNumEvents")->Integral(), nMC, nMC / getData("hSecStudyNumEvents")->Integral(), nData / nMC) << std::endl;
    }

    // Angle relative to the z axis, MC normalized to the data integral (shape comparison)
    std::vector<Comparison> shapeComparisons = {
        {"hSecCosZ",             "CosZAllSecondaries",     "Secondary tracks near vertex (area normalized)",       "cos#theta_{z}",              rebinCos},
        {"hSecCosZProtonTag",    "CosZProtonTagged",       "Proton-tagged secondaries (area normalized)",          "cos#theta_{z}",              rebinCos},
        {"hSecLeadProtonCosZ",   "CosZLeadingProton",      "Leading proton-tagged secondary (area normalized)",    "cos#theta_{z}",              rebinCos},
        {"hSecLeadProtonRecoKE", "RecoKELeadingProton",    "Leading proton-tagged secondary (area normalized)",    "Reco calorimetric KE [MeV]", rebinKE},
    };
    for (const auto& comp : shapeComparisons) {
        PrintDataVsMCContribPlot(
            SaveDir, comp.name,
            rebinned(getData(comp.base), comp.rebin),
            getMCSplit(comp.base, true, comp.rebin),
            trueLabels, mcColors,
            comp.title, comp.xlabel, "Tracks",
            FontStyle, TextSize, true, nullptr, true
        );
    }

    // Angle of secondaries near the vertex, separately for events selected as scatter and as abs Np.
    // Built from the per-track summaries in the SelectedEvents trees (stitched primary segments not included).
    {
        TPRegexp cosZRegex("cosZ=([-0-9.]+).*?(?:true pdg=(-?\\d+) KE=([-0-9.]+) MeV)?$");
        TPRegexp dedxRegex("<dEdx>=([-0-9.]+).*?(?:true pdg=(-?\\d+) KE=([-0-9.]+) MeV)?$");
        auto fillFromTree = [&](TFile* f, int selCategory, bool isMC, TH1D* hData, std::vector<TH1*>& hMC, bool protonsOnly = false, TPRegexp* trackRegex = nullptr) {
            if (!trackRegex) trackRegex = &cosZRegex;
            TTree* t = (TTree*) f->Get("SelectedEvents");
            int category; double weight; std::string* tracks = nullptr;
            t->SetBranchAddress("category", &category);
            t->SetBranchAddress("weight", &weight);
            t->SetBranchAddress("tracks", &tracks);
            for (Long64_t i = 0; i < t->GetEntries(); ++i) {
                t->GetEntry(i);
                if (category != selCategory) continue;
                TObjArray* lines = TString(tracks->c_str()).Tokenize("\n");
                for (int j = 0; j < lines->GetEntries(); ++j) {
                    TString line = ((TObjString*) lines->At(j))->GetString();
                    if (!line.BeginsWith("trk ")) continue;
                    if (protonsOnly && !line.Contains("tag=p(")) continue;
                    TObjArray* m = trackRegex->MatchS(line);
                    if (m->GetEntries() >= 2) {
                        double value = ((TObjString*) m->At(1))->GetString().Atof();
                        if (!isMC) hData->Fill(value);
                        else {
                            int pdgCode = (m->GetEntries() >= 4) ? ((TObjString*) m->At(2))->GetString().Atoi() : 0;
                            double ke   = (m->GetEntries() >= 4) ? ((TObjString*) m->At(3))->GetString().Atof() / 1000. : 0;
                            int cat;
                            if (pdgCode == 2212)                              cat = (ke >= PROTON_ENERGY_LOWER_BOUND) ? 0 : 1;
                            else if (std::abs(pdgCode) == 211)                cat = (ke >= PION_SCATTERING_ENERGY_THRESHOLD) ? 2 : 3;
                            else if (std::abs(pdgCode) == 13)                 cat = 4;
                            else if (std::abs(pdgCode) == 11 || pdgCode == 22) cat = 5;
                            else                                              cat = 6;
                            hMC[cat]->Fill(value, weight);
                        }
                    }
                    delete m;
                }
                delete lines;
            }
            t->ResetBranchAddresses();
        };

        std::vector<std::pair<int, TString>> selections = {{2, "Scatter"}, {1, "AbsNp"}};
        std::vector<TString> selectionLabels = {"selected as scatter", "selected as abs Np"};
        for (size_t s = 0; s < selections.size(); ++s) {
            TString tag = selections[s].second;
            TH1D* hData = new TH1D("hCosZData" + tag, ";cos#theta_{z};", 5, -1, 1);
            std::vector<TH1*> hMC;
            for (const TString& name : trueNames) hMC.push_back(new TH1D("hCosZMC" + tag + name, ";cos#theta_{z};", 5, -1, 1));
            fillFromTree(fData.get(), selections[s].first, false, hData, hMC);
            fillFromTree(fMC.get(), selections[s].first, true, hData, hMC);

            double mcTotal = 0;
            for (TH1* h : hMC) mcTotal += h->Integral();
            for (TH1* h : hMC) h->Scale(hData->Integral() / mcTotal);
            std::cout << "CosZ secondaries, " << selectionLabels[s] << ": data tracks = " << hData->Integral() << std::endl;

            PrintDataVsMCContribPlot(
                SaveDir, "CosZAllSecondaries" + tag,
                hData, hMC, trueLabels, mcColors,
                "Secondary tracks near vertex, " + selectionLabels[s] + " (area normalized)", "cos#theta_{z}", "Tracks",
                FontStyle, TextSize, true, nullptr, true
            );
        }

        // Same, but only proton-tagged tracks (chi^2 or mean dE/dx) in the selected events
        for (size_t s = 0; s < selections.size(); ++s) {
            TString tag = selections[s].second;
            TH1D* hData = new TH1D("hCosZProtonsData" + tag, ";cos#theta_{z};", 5, -1, 1);
            std::vector<TH1*> hMC;
            for (const TString& name : trueNames) hMC.push_back(new TH1D("hCosZProtonsMC" + tag + name, ";cos#theta_{z};", 5, -1, 1));
            fillFromTree(fData.get(), selections[s].first, false, hData, hMC, true);
            fillFromTree(fMC.get(), selections[s].first, true, hData, hMC, true);

            double mcTotal = 0;
            for (TH1* h : hMC) mcTotal += h->Integral();
            for (TH1* h : hMC) h->Scale(hData->Integral() / mcTotal);
            std::cout << "CosZ proton-tagged secondaries, " << selectionLabels[s] << ": data tracks = " << hData->Integral()
                      << ", MC true p >= threshold fraction = " << hMC[0]->Integral() / hData->Integral() << std::endl;

            PrintDataVsMCContribPlot(
                SaveDir, "CosZSelectedProtons" + tag,
                hData, hMC, trueLabels, mcColors,
                "Proton-tagged tracks, " + selectionLabels[s] + " (area normalized)", "cos#theta_{z}", "Tracks",
                FontStyle, TextSize, true, nullptr, true
            );
        }

        // Mean dE/dx of secondaries near the vertex, separately for events selected as scatter and as abs Np
        for (size_t s = 0; s < selections.size(); ++s) {
            TString tag = selections[s].second;
            TH1D* hData = new TH1D("hMeanDEDXData" + tag, ";Mean dE/dx [MeV/cm];", 10, 0, 10);
            std::vector<TH1*> hMC;
            for (const TString& name : trueNames) hMC.push_back(new TH1D("hMeanDEDXMC" + tag + name, ";Mean dE/dx [MeV/cm];", 10, 0, 10));
            fillFromTree(fData.get(), selections[s].first, false, hData, hMC, false, &dedxRegex);
            fillFromTree(fMC.get(), selections[s].first, true, hData, hMC, false, &dedxRegex);

            double mcTotal = 0;
            for (TH1* h : hMC) mcTotal += h->Integral();
            for (TH1* h : hMC) h->Scale(hData->Integral() / mcTotal);

            PrintDataVsMCContribPlot(
                SaveDir, "MeanDEDXSecondaries" + tag,
                hData, hMC, trueLabels, mcColors,
                "Secondary tracks near vertex, " + selectionLabels[s] + " (area normalized)", "Mean dE/dx [MeV/cm]", "Tracks",
                FontStyle, TextSize, true, nullptr, true
            );
        }
    }

    ////////////////////////////////////////////////
    // Proton threshold: reco outcome vs true KE  //
    ////////////////////////////////////////////////

    DrawFractionCurves(
        SaveDir, "ProtonThresholdOutcome",
        {getMC("hLeadingProtonKERecoScatter"), getMC("hLeadingProtonKERecoNp"), getMC("hLeadingProtonKEReco0pCandidate"), getMC("hLeadingProtonKERejected")},
        getMC("hLeadingProtonKEAll"),
        {"Reco scatter", "Reco abs Np", "Reco abs 0p (before cluster cut)", "Rejected"},
        {colBlue, colOrange, colGreen, colGray},
        1000 * PROTON_ENERGY_LOWER_BOUND, Form("Signal threshold (%.0f MeV)", 1000 * PROTON_ENERGY_LOWER_BOUND),
        "True absorption by reco outcome (MC, [-10, 0) = no proton)", "True leading proton KE [MeV]", "Fraction",
        FontStyle, TextSize
    );

    DrawFractionCurves(
        SaveDir, "ProtonAngleOutcome",
        {getMC("hLeadingProtonCosRecoScatter"), getMC("hLeadingProtonCosRecoNp"), getMC("hLeadingProtonCosReco0pCandidate"), getMC("hLeadingProtonCosRejected")},
        getMC("hLeadingProtonCosAll"),
        {"Reco scatter", "Reco abs Np", "Reco abs 0p (before cluster cut)", "Rejected"},
        {colBlue, colOrange, colGreen, colGray},
        -2, "",
        Form("True absorption, leading p KE #geq %.0f MeV, by reco outcome (MC)", 1000 * PROTON_ENERGY_LOWER_BOUND), "True leading proton cos#theta_{z}", "Fraction",
        FontStyle, TextSize, 2
    );

    DrawFractionCurves(
        SaveDir, "ProtonEfficiencyVsAngle",
        {getMC("hEffProtonCosReco"), getMC("hEffProtonCosTagChi2"), getMC("hEffProtonCosTagAny")},
        getMC("hEffProtonCosAll"),
        {"Reconstructed near vertex", "Tagged proton (#chi^{2})", "Tagged proton (#chi^{2} or mean dE/dx)"},
        {colGray, colBlue, colOrange},
        -2, "",
        Form("Secondary proton efficiency (MC, KE #geq %.0f MeV)", 1000 * PROTON_ENERGY_LOWER_BOUND), "True cos#theta_{z}", "Efficiency",
        FontStyle, TextSize, 2
    );

    // Print the proton outcome fractions as a table
    {
        TH1D* hAll = getMC("hLeadingProtonKEAll");
        std::cout << std::endl << "=== Reco outcome vs true leading proton KE (true absorption) ===" << std::endl;
        std::cout << Form("%-12s %10s %8s %8s %8s %8s", "KE [MeV]", "N", "Scat", "Np", "0p", "Rej") << std::endl;
        for (int b = 1; b <= hAll->GetNbinsX(); ++b) {
            double n = hAll->GetBinContent(b);
            if (n <= 0) continue;
            std::cout << Form("[%3.0f, %3.0f]  %10.1f %8.3f %8.3f %8.3f %8.3f",
                hAll->GetXaxis()->GetBinLowEdge(b), hAll->GetXaxis()->GetBinUpEdge(b), n,
                getMC("hLeadingProtonKERecoScatter")->GetBinContent(b) / n,
                getMC("hLeadingProtonKERecoNp")->GetBinContent(b) / n,
                getMC("hLeadingProtonKEReco0pCandidate")->GetBinContent(b) / n,
                getMC("hLeadingProtonKERejected")->GetBinContent(b) / n) << std::endl;
        }
    }

    ///////////////////////////////////////
    // True vs reco KE for each species  //
    ///////////////////////////////////////

    // Also split by whether the far end of the track is contained in the TPC
    std::vector<TString> speciesNames  = {"Proton", "Pion", "Muon", "EM", "Other"};
    std::vector<TString> speciesLabels = {"True p", "True #pi^{#pm}", "True #mu", "True e/#gamma", "Other / unmatched"};
    std::vector<std::pair<TString, TString>> containmentSplits = {{"", ""}, {"Contained", ", contained"}, {"Exiting", ", exiting"}};
    for (const auto& split : containmentSplits)
    for (size_t i = 0; i < speciesNames.size(); ++i) {
        TH2D* h = (TH2D*) fMC->Get("hSecTrueVsRecoKE" + split.first + speciesNames[i]);
        if (!h) continue;
        h->Rebin2D(2, 2);
        TCanvas* c = new TCanvas("c2d_" + split.first + speciesNames[i], "", 1200, 1000);
        c->SetLeftMargin(0.14); c->SetRightMargin(0.15); c->SetBottomMargin(0.13); c->SetLogz();
        h->SetTitle("Secondary tracks near vertex" + split.second + ": " + speciesLabels[i] + " (MC)");
        for (TAxis* ax : {h->GetXaxis(), h->GetYaxis(), h->GetZaxis()}) {
            ax->SetTitleFont(FontStyle); ax->SetLabelFont(FontStyle);
            ax->SetTitleSize(TextSize * 0.9); ax->SetLabelSize(TextSize * 0.8);
        }
        h->GetXaxis()->CenterTitle(); h->GetYaxis()->CenterTitle();
        h->GetYaxis()->SetTitleOffset(1.2);
        h->Draw("COLZ");
        TLine* diag = new TLine(0, 0, 400, 400);
        diag->SetLineStyle(2); diag->SetLineWidth(2); diag->Draw();
        c->SaveAs(SaveDir + "TrueVsRecoKE" + split.first + speciesNames[i] + ".png");
        delete c;
    }

    //////////////////////////////////////////////////////////
    // dE/dx vs residual range of secondaries, with templates //
    //////////////////////////////////////////////////////////

    TGraph* gProton = new TGraph(); initializeProtonPoints(gProton);
    TGraph* gPion   = new TGraph(); initializePionPoints(gPion);
    gProton->SetLineColor(kRed);  gProton->SetLineWidth(3);
    gPion->SetLineColor(kBlack);  gPion->SetLineWidth(3); gPion->SetLineStyle(2);

    std::vector<std::pair<TH2D*, TString>> dedxPlots = {{(TH2D*) fData->Get("hSecDEDXvsResR"), "Data"}};
    for (size_t i = 0; i < speciesNames.size(); ++i) {
        dedxPlots.push_back({(TH2D*) fMC->Get("hSecDEDXvsResR" + speciesNames[i]), "MC " + speciesLabels[i]});
    }
    TH2D* hMCAll = nullptr;
    for (size_t i = 1; i < dedxPlots.size(); ++i) {
        if (!dedxPlots[i].first) continue;
        if (!hMCAll) hMCAll = (TH2D*) dedxPlots[i].first->Clone("hSecDEDXvsResRMCAll");
        else hMCAll->Add(dedxPlots[i].first);
    }
    dedxPlots.push_back({hMCAll, "MC all"});

    for (size_t i = 0; i < dedxPlots.size(); ++i) {
        TH2D* h = dedxPlots[i].first;
        if (!h) continue;
        TString name = dedxPlots[i].second; name.ReplaceAll(" ", ""); name.ReplaceAll("#", ""); name.ReplaceAll("{", ""); name.ReplaceAll("}", ""); name.ReplaceAll("/", ""); name.ReplaceAll("^", ""); name.ReplaceAll("pm", "");
        TCanvas* c = new TCanvas("cdedx_" + name, "", 1200, 1000);
        c->SetLeftMargin(0.14); c->SetRightMargin(0.15); c->SetBottomMargin(0.13); c->SetLogz();
        h->SetTitle("Secondary tracks near vertex: " + dedxPlots[i].second);
        for (TAxis* ax : {h->GetXaxis(), h->GetYaxis(), h->GetZaxis()}) {
            ax->SetTitleFont(FontStyle); ax->SetLabelFont(FontStyle);
            ax->SetTitleSize(TextSize * 0.9); ax->SetLabelSize(TextSize * 0.8);
        }
        h->GetXaxis()->CenterTitle(); h->GetYaxis()->CenterTitle();
        h->Draw("COLZ");
        gProton->Draw("L SAME");
        gPion->Draw("L SAME");
        TLegend* leg = new TLegend(0.55, 0.75, 0.83, 0.88);
        leg->SetTextFont(FontStyle); leg->SetTextSize(TextSize * 0.7);
        leg->AddEntry(gProton, "Proton template", "l");
        leg->AddEntry(gPion, "Pion template", "l");
        leg->Draw();
        c->SaveAs(SaveDir + "DEDXvsResR_" + name + ".png");
        delete c;
    }

    // Profile comparison: mean dE/dx vs residual range for data and MC (all secondaries)
    if (fData->Get("hSecDEDXvsResR") && hMCAll) {
        TProfile* pData = ((TH2D*) fData->Get("hSecDEDXvsResR"))->ProfileX("pDEDXData");
        TProfile* pMC   = hMCAll->ProfileX("pDEDXMC");
        TProfile* pMCp  = ((TH2D*) fMC->Get("hSecDEDXvsResRProton"))->ProfileX("pDEDXMCProton");
        TCanvas* c = new TCanvas("cdedxprof", "", 1300, 900);
        c->SetLeftMargin(0.13); c->SetBottomMargin(0.14); c->SetGridy();
        pMC->SetLineColor(colBlue); pMC->SetMarkerColor(colBlue); pMC->SetMarkerStyle(21);
        pMCp->SetLineColor(colOrange); pMCp->SetMarkerColor(colOrange); pMCp->SetMarkerStyle(22);
        pData->SetLineColor(kBlack); pData->SetMarkerColor(kBlack); pData->SetMarkerStyle(20);
        pMC->SetTitle("Mean dE/dx of secondaries near vertex;Residual range [cm];Mean dE/dx [MeV/cm]");
        pMC->GetYaxis()->SetRangeUser(0, 25);
        for (TAxis* ax : {pMC->GetXaxis(), pMC->GetYaxis()}) {
            ax->SetTitleFont(FontStyle); ax->SetLabelFont(FontStyle);
            ax->SetTitleSize(TextSize); ax->SetLabelSize(TextSize * 0.9); ax->CenterTitle();
        }
        pMC->Draw("E1P");
        pMCp->Draw("E1P SAME");
        pData->Draw("E1P SAME");
        gProton->Draw("L SAME");
        gPion->Draw("L SAME");
        TLegend* leg = new TLegend(0.50, 0.62, 0.93, 0.90);
        leg->SetTextFont(FontStyle); leg->SetTextSize(TextSize * 0.75); leg->SetBorderSize(0);
        leg->AddEntry(pData, "Data (all secondaries)", "lep");
        leg->AddEntry(pMC, "MC (all secondaries)", "lep");
        leg->AddEntry(pMCp, "MC (true protons)", "lep");
        leg->AddEntry(gProton, "Proton template", "l");
        leg->AddEntry(gPion, "Pion template", "l");
        leg->Draw();
        c->SaveAs(SaveDir + "DEDXvsResRProfile.png");
        delete c;
    }

    TH1* hMCNumProtons = (TH1*) getMC("hSecNumProtonTagged")->Clone("hMCNumProtonsScaled");
    hMCNumProtons->Scale(scale);
    PrintDataVsMCContribPlot(
        SaveDir, "NumProtonTagged",
        getData("hSecNumProtonTagged"),
        {hMCNumProtons}, {"MC"}, {colBlue},
        "Proton-tagged secondaries per event", "Number of proton-tagged tracks (4 = 4+)", "Events",
        FontStyle, TextSize, true, nullptr, true
    );

    ////////////////////////////////////////////
    // Containment of secondaries near vertex //
    ////////////////////////////////////////////

    // Far end of the track inside the TPC (CONTAINMENT_MARGIN from the x/y walls), MC area normalized to data
    std::vector<std::pair<TString, TString>> containmentTags = {
        {"", "Secondary tracks near vertex"}, {"ProtonTag", "Proton-tagged secondaries"}, {"PionTag", "Pion-tagged secondaries"}
    };
    std::cout << std::endl << "=== Contained fraction of secondaries near the vertex ===" << std::endl;
    for (const auto& tag : containmentTags) {
        TH1D* hData = (TH1D*) getData("hSecContained" + tag.first);
        if (!hData) continue;
        hData = (TH1D*) hData->Clone("hContainedData" + tag.first);
        hData->GetXaxis()->SetBinLabel(1, "Exiting");
        hData->GetXaxis()->SetBinLabel(2, "Contained");
        std::vector<TH1*> hMC = getMCSplit("hSecContained" + tag.first, true);

        PrintDataVsMCContribPlot(
            SaveDir, "Containment" + tag.first,
            hData, hMC, trueLabels, mcColors,
            tag.second + " (area normalized)", "", "Tracks",
            FontStyle, TextSize, true, nullptr, true
        );

        double nData = hData->Integral(), fData = hData->GetBinContent(2) / nData;
        double nMC = 0, nMCContained = 0;
        for (TH1* h : hMC) { nMC += h->Integral(); nMCContained += h->GetBinContent(2); }
        std::cout << Form("%-28s data %.3f +- %.3f (%4.0f tracks)   MC %.3f",
            tag.second.Data(), fData, std::sqrt(fData * (1 - fData) / nData), nData, nMCContained / nMC) << std::endl;
    }

    //////////////////////////////////////////////////////////////////
    // Length vs opening angle, before and after primary stitching  //
    //////////////////////////////////////////////////////////////////

    // Before = vertex at the end of the primary track; After = vertex at the stitching break point, with the stitched segment as a secondary
    std::vector<std::pair<TString, TString>> lenAngleSets = {
        {"Before",         "before stitching"},
        {"After",          "after stitching"},
        {"BeforeStitched", "before stitching, stitched events only"},
        {"AfterStitched",  "after stitching, stitched events only"}
    };
    std::vector<TString> lenAngleCats = trueNames; lenAngleCats.push_back("StitchedSegment");
    std::vector<TString> lenAngleLabels = trueLabels; lenAngleLabels.push_back("Stitched segment");
    std::vector<int> lenAngleColors = mcColors; lenAngleColors.push_back(kGray + 1);
    const int rebinAngle = 4; // 5 deg -> 20 deg

    for (const auto& set : lenAngleSets) {
        TH2D* hData2D = (TH2D*) fData->Get("hSecLenAngle" + set.first);
        if (!hData2D) continue;
        hData2D = (TH2D*) hData2D->Clone("hLenAngleData" + set.first);
        hData2D->Rebin2D(rebinAngle, rebinLength);

        std::vector<TH2D*> hMC2D;
        TH2D* hMCTotal2D = nullptr;
        for (const TString& cat : lenAngleCats) {
            TH2D* h = (TH2D*) ((TH2D*) fMC->Get("hSecLenAngle" + set.first + cat))->Clone("hLenAngleMC" + set.first + cat);
            h->Rebin2D(rebinAngle, rebinLength);
            hMC2D.push_back(h);
            if (!hMCTotal2D) hMCTotal2D = (TH2D*) h->Clone("hLenAngleMCTotal" + set.first);
            else hMCTotal2D->Add(h);
        }
        double nData = hData2D->Integral(), nMC = hMCTotal2D->Integral();
        std::cout << "Length vs angle " << set.second << ": data tracks = " << nData << ", MC (weighted) = " << nMC << std::endl;
        if (nData <= 0 || nMC <= 0) continue;

        // 2D: data and MC side by side, both normalized to unit area with a common z scale
        TH2D* hDataNorm = (TH2D*) hData2D->Clone("hLenAngleDataNorm" + set.first);
        TH2D* hMCNorm   = (TH2D*) hMCTotal2D->Clone("hLenAngleMCNorm" + set.first);
        hDataNorm->Scale(1. / nData);
        hMCNorm->Scale(1. / nMC);
        double zMax = std::max(hDataNorm->GetMaximum(), hMCNorm->GetMaximum());

        TCanvas* c = new TCanvas("cLenAngle" + set.first, "", 2200, 1000);
        c->Divide(2, 1);
        std::vector<std::pair<TH2D*, TString>> panels = {{hDataNorm, "Data"}, {hMCNorm, "MC"}};
        for (size_t iPanel = 0; iPanel < panels.size(); ++iPanel) {
            c->cd(iPanel + 1);
            gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.16); gPad->SetBottomMargin(0.13);
            TH2D* h = panels[iPanel].first;
            h->SetTitle(panels[iPanel].second + ": secondaries near vertex, " + set.second + " (area normalized)");
            for (TAxis* ax : {h->GetXaxis(), h->GetYaxis(), h->GetZaxis()}) {
                ax->SetTitleFont(FontStyle); ax->SetLabelFont(FontStyle);
                ax->SetTitleSize(TextSize * 0.9); ax->SetLabelSize(TextSize * 0.8);
            }
            h->GetXaxis()->CenterTitle(); h->GetYaxis()->CenterTitle();
            h->GetZaxis()->SetTitle("Fraction of tracks");
            h->SetMinimum(0); h->SetMaximum(zMax);
            h->Draw("COLZ");
        }
        c->SaveAs(SaveDir + "LengthVsAngle" + set.first + ".png");
        delete c;

        // 1D projections, MC stacked by true particle and area normalized to data
        for (int axis = 0; axis < 2; ++axis) {
            TString proj = (axis == 0) ? "Angle" : "Length";
            TH1D* hDataProj = (axis == 0) ? hData2D->ProjectionX("hLenAngleDataProj" + proj + set.first) : hData2D->ProjectionY("hLenAngleDataProj" + proj + set.first);
            std::vector<TH1*> hMCProj;
            for (size_t iCat = 0; iCat < hMC2D.size(); ++iCat) {
                TString name = "hLenAngleMCProj" + proj + set.first + lenAngleCats[iCat];
                TH1D* h = (axis == 0) ? hMC2D[iCat]->ProjectionX(name) : hMC2D[iCat]->ProjectionY(name);
                h->Scale(nData / nMC);
                hMCProj.push_back(h);
            }
            PrintDataVsMCContribPlot(
                SaveDir, proj + set.first,
                hDataProj, hMCProj, lenAngleLabels, lenAngleColors,
                "Secondaries near vertex, " + set.second + " (area normalized)",
                (axis == 0) ? "Opening angle w.r.t. incident track [deg]" : "Track length [cm]", "Tracks",
                FontStyle, TextSize, true, nullptr, true
            );
        }
    }
}
