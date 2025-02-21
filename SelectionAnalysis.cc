#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TEfficiency.h>

#include <vector>
#include <map>

struct EventInfo {
    int pionTruePDG;
    std::string pionTrueProcess;
    std::vector<int> trueDaughtersPDG;
    std::vector<std::string> trueDaughterProcess;
};

void SelectionAnalysis() {
    // Set defaults
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    // gStyle->SetOptStat(0);

    int FontStyle = 132;
    double TextSize = 0.06;
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis/figs/SelectionAnalysis/";

    // Load root files
    TString RootFilePath = "/exp/lariat/app/users/epelaez/files/PionAbsorptionSelection_histo.root"; 
    std::unique_ptr<TFile> File(TFile::Open(RootFilePath));
    TDirectory* Directory = (TDirectory*)File->Get("PionAbsorptionSelection");
    TTree* tree = (TTree*) Directory->Get<TTree>("PionAbsorptionSelectionTree");

    // Grab file with true selected events
    TString TrueSelectedPath = "/exp/lariat/app/users/epelaez/files/SignalPionAbsorptionOnSelectedEvents_histo.root";
    std::unique_ptr<TFile> TrueSelectedFile(TFile::Open(TrueSelectedPath));
    TDirectory* TrueSelectedDirectory = (TDirectory*)TrueSelectedFile->Get("SignalPionAbsorption");
    TTree* trueRecoTree = (TTree*) TrueSelectedDirectory->Get<TTree>("PionAbsorptionTree");

    // Grab file with all true events
    TString TruePath = "/exp/lariat/app/users/epelaez/files/SignalPionAbsorption_histo.root";
    std::unique_ptr<TFile> TrueFile(TFile::Open(TruePath));
    TDirectory* TrueDirectory = (TDirectory*)TrueFile->Get("SignalPionAbsorption");
    TTree* trueTree = (TTree*) TrueDirectory->Get<TTree>("PionAbsorptionTree");

    // Load reco branches
    int run, subrun, event;
    tree->SetBranchAddress("Run", &run); 
    tree->SetBranchAddress("Subrun", &subrun); 
    tree->SetBranchAddress("Event", &event);

    int pionTruthPDG;
    std::string* pionTruthProcess = new std::string();
    std::vector<int>* pionDaughtersPDG = nullptr;
    std::vector<std::string>* pionDaughtersProcess = nullptr;
    tree->SetBranchAddress("pionTruthPDG", &pionTruthPDG);
    tree->SetBranchAddress("pionTruthProcess", &pionTruthProcess);
    tree->SetBranchAddress("pionDaughtersPDG", &pionDaughtersPDG);
    tree->SetBranchAddress("pionDaughtersProcess", &pionDaughtersProcess);

    int protonCount;
    tree->SetBranchAddress("protonCount", &protonCount);

    std::vector<double>* protonLength = nullptr;
    tree->SetBranchAddress("protonLength", &protonLength);

    std::vector<std::string>* protonTrueProcess = nullptr;
    tree->SetBranchAddress("protonTrueProcess", &protonTrueProcess);

    // Load true reco branches
    int trueRecoRun, trueRecoSubrun, trueRecoEvent;
    trueRecoTree->SetBranchAddress("run", &trueRecoRun);
    trueRecoTree->SetBranchAddress("subrun", &trueRecoSubrun);
    trueRecoTree->SetBranchAddress("event", &trueRecoEvent);

    int trueRecoProtonCount;
    trueRecoTree->SetBranchAddress("iNumProtons", &trueRecoProtonCount);

    // Declare histograms
    TH1I* hRecoNumProtonTracks = new TH1I("hRecoNumProtonTracks", "NumProtonTracks;;", 5, 0, 5);
    TH1I* hTrueNumProtonTracks = new TH1I("hTrueNumProtonTracks", "NumProtonTracks;;", 5, 0, 5);

    // We will keep track of weird events
    std::vector<EventInfo> FlaggedEvents;
    std::map<int, int> pionPDGCount;

    // Loop over reco events
    Int_t NumEntries = (Int_t) tree->GetEntries();
    for (Int_t i = 0; i < NumEntries; ++i) {
        tree->GetEntry(i);

        // Fill histograms
        hRecoNumProtonTracks->Fill(protonCount);
 
        // Flag events
        bool flagEvent = false;

        if (pionTruthPDG != -211) flagEvent = true;
        if (*pionTruthProcess != "primary") flagEvent = true;

        int numPionDaughters = pionDaughtersPDG->size();
        if (!flagEvent) {
            for (int iDaughter = 0; iDaughter < numPionDaughters; ++iDaughter) {
                int daughterPDG = pionDaughtersPDG->at(iDaughter);
                std::string daughterProcess = pionDaughtersProcess->at(iDaughter);

                if ((daughterPDG == 11) && (daughterProcess == "hIoni")) continue;
                if ((daughterPDG == 111) || (daughterPDG == 211) || (daughterPDG == -211)) { flagEvent = true; break; }
                if ((daughterProcess == "Decay") || (daughterProcess == "hBertiniCaptureAtRest")) { flagEvent = true; break; }
                if ((daughterPDG == 2212) && (daughterProcess != "pi-Inelastic")) {flagEvent = true; break; }
            }
        }

        if (flagEvent) {
            EventInfo event;
            event.pionTruePDG = pionTruthPDG;
            event.pionTrueProcess = *pionTruthProcess;
            event.trueDaughtersPDG = *pionDaughtersPDG;
            event.trueDaughterProcess = *pionDaughtersProcess;
            FlaggedEvents.push_back(event);

            pionPDGCount[pionTruthPDG]++;
        }
    }

    std::ofstream outFile("FlaggedEventsOutput.txt");
    outFile << "Num of flagged events: " << FlaggedEvents.size() << std::endl;
    for (const auto &entry : pionPDGCount) {
        outFile << "   Pion truth matched PDG: " << entry.first << " Count: " << entry.second << std::endl;
    }
    outFile << std::endl;

    for (const auto &event : FlaggedEvents) {
        outFile << "Pion truth matched PDG: " << event.pionTruePDG << std::endl;
        outFile << "Pion truth matched process: " << event.pionTrueProcess << std::endl;
        outFile << "Daughters: " << std::endl;
        for (int n = 0; n < event.trueDaughterProcess.size(); ++n) {
            if (event.trueDaughtersPDG[n] == 11) continue;
            outFile << "    PDG: " << event.trueDaughtersPDG[n];
            outFile << " Process: " << event.trueDaughterProcess[n];
            outFile << std::endl;
        }
        outFile << std::endl;
    }

    // Loop over true events
    Int_t NumTrueRecoEntries = (Int_t) trueRecoTree->GetEntries();
    for (Int_t i = 0; i < NumTrueRecoEntries; ++i) {
        trueRecoTree->GetEntry(i);

        // Fill histograms
        hTrueNumProtonTracks->Fill(trueRecoProtonCount);
    }

    // Setup for drawing plots
    std::vector<int> Colors = {
        kBlack, kBlue, kRed, kGreen
    };

    std::vector<std::vector<TH1*>> PlotGroups = {
        {hRecoNumProtonTracks, hTrueNumProtonTracks}
    };

    std::vector<std::vector<TString>> PlotLabelGroups = {
        {"Reco", "True reco"}, 
    };

    std::vector<TString> PlotTitles = {
        "NumProtons",
    };

    std::vector<TString> XLabels = {
        "# of protons",
    };

    int numPlots = PlotGroups.size();
    for (int iPlot = 0; iPlot < numPlots; ++iPlot) {
        // Set up canvas for drawing
        TCanvas* PlotCanvas = new TCanvas("Canvas", "Canvas", 205, 34, 1024, 768);
        PlotCanvas->cd();
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetBottomMargin(0.15);

        TLegend* leg = new TLegend(0.65,0.65,0.85,0.75);
        leg->SetBorderSize(0);
        leg->SetTextSize(TextSize * 0.8);
        leg->SetTextFont(FontStyle);

        // Get histograms
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
            leg->AddEntry(Plots[iSubPlot], Labels[iSubPlot], "l");
            Plots[iSubPlot]->SetLineWidth(2);
            Plots[iSubPlot]->SetLineColor(Colors.at(iSubPlot));
            Plots[iSubPlot]->Draw("hist same");
        }

        leg->Draw();
        PlotCanvas->SaveAs(SaveDir + PlotTitles.at(iPlot) + ".png");
        delete PlotCanvas;
    }

    // Compute purity and efficiency
    Int_t NumTrueEntries = (Int_t) trueTree->GetEntries();
    std::cout << "Number of reco events that pass selection criteria: " << NumEntries << std::endl;
    std::cout << "Number of reco events that pass selection criteria that are true events: " << NumTrueRecoEntries << std::endl;
    std::cout << "Number of background events: " << NumEntries - NumTrueRecoEntries << std::endl;
    std::cout << "Number of true events: " << NumTrueEntries << std::endl;
    std::cout << "Efficiency: " << 100 * ((float)NumTrueRecoEntries / (float)NumTrueEntries) << "%" << std::endl; 
    std::cout << "Purity: " << 100 * ((float)NumTrueRecoEntries / (float)NumEntries) << "%" << std::endl;
}