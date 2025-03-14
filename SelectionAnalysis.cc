#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TEfficiency.h>

#include <vector>
#include <map>

struct EventInfo {
    int run; int subrun; int event;

    int                      wcMatchPDG;
    std::string              wcMatchProcess;
    std::vector<int>         wcMatchDaughtersPDG;
    std::vector<std::string> wcMatchDaughtersProcess;

    int                      truthPrimaryPDG;
    std::vector<int>         truthPrimaryDaughtersPDG;
    std::vector<std::string> truthPrimaryDaughtersProcess;
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

    // Truth information about the particle the wire chamber track matched to
    int                       wcMatchPDG;
    std::string*              wcMatchProcess = new std::string();
    std::vector<int>*         wcMatchDaughtersPDG = nullptr;
    std::vector<std::string>* wcMatchDaughtersProcess = nullptr;
    tree->SetBranchAddress("wcMatchPDG", &wcMatchPDG);
    tree->SetBranchAddress("wcMatchProcess", &wcMatchProcess);
    tree->SetBranchAddress("wcMatchDaughtersPDG", &wcMatchDaughtersPDG);
    tree->SetBranchAddress("wcMatchDaughtersProcess", &wcMatchDaughtersProcess);

    // Truth information about the truth primary particle
    int                       truthPrimaryPDG;
    std::vector<int>*         truthPrimaryDaughtersPDG = nullptr;
    std::vector<std::string>* truthPrimaryDaughtersProcess = nullptr;
    std::vector<float>*       truthPrimaryDaughtersKE = nullptr;
    bool                      isPionAbsorptionSignal;
    int                       numVisibleProtons;
    tree->SetBranchAddress("truthPrimaryPDG", &truthPrimaryPDG);
    tree->SetBranchAddress("truthPrimaryDaughtersPDG", &truthPrimaryDaughtersPDG);
    tree->SetBranchAddress("truthPrimaryDaughtersProcess", &truthPrimaryDaughtersProcess);
    tree->SetBranchAddress("truthPrimaryDaughtersKE", &truthPrimaryDaughtersKE);
    tree->SetBranchAddress("isPionAbsorptionSignal", &isPionAbsorptionSignal);
    tree->SetBranchAddress("numVisibleProtons", &numVisibleProtons);

    // Information about reco'ed protons
    int                       protonCount;
    std::vector<double>*      protonLength = nullptr;
    std::vector<std::string>* protonTrueProcess = nullptr;
    tree->SetBranchAddress("protonCount", &protonCount);
    tree->SetBranchAddress("protonLength", &protonLength);
    tree->SetBranchAddress("protonTrueProcess", &protonTrueProcess);

    // Load true reco branches
    int trueRecoRun, trueRecoSubrun, trueRecoEvent;
    trueRecoTree->SetBranchAddress("run", &trueRecoRun);
    trueRecoTree->SetBranchAddress("subrun", &trueRecoSubrun);
    trueRecoTree->SetBranchAddress("event", &trueRecoEvent);

    int trueRecoProtonCount;
    trueRecoTree->SetBranchAddress("iNumProtons", &trueRecoProtonCount);

    // Load true branches
    int trueProtonCount;
    trueTree->SetBranchAddress("iNumProtons", &trueProtonCount);

    // Declare histograms
    TH1I* hRecoNumProtonTracks = new TH1I("hRecoNumProtonTracks", "NumProtonTracks;;", 5, 0, 5);
    TH1I* hRecoTrueNumProtons = new TH1I("hRecoTrueNumProtons", "NumProtonTracks;;", 5, 0, 5);
    TH1I* hTruthRecoTrueNumProtons = new TH1I("hTruthRecoTrueNumProtons", "NumProtonTracks;;", 5, 0, 5);

    /////////////////
    // Selection tree
    /////////////////
    
    // Keep track of events in 0p and Np bins
    int numSelectedEvents0Protons = 0;
    int numSelectedEventsNProtons = 0;
    int numSelectedTrueEvents0Protons = 0;
    int numSelectedTrueEventsNProtons = 0;

    // We will keep track of weird events
    std::vector<EventInfo> FlaggedEvents;
    std::map<int, std::map<std::string, int>> pionBackground;
    std::map<int, int> wcMatchPDGCount;
    std::map<int, int> primaryPDGCount;

    // Loop over reco events
    Int_t NumEntries = (Int_t) tree->GetEntries();
    for (Int_t i = 0; i < NumEntries; ++i) {
        tree->GetEntry(i);
 
        if (protonCount >= 1) {
            numSelectedEventsNProtons++;
            if (isPionAbsorptionSignal && (numVisibleProtons >= 1)) numSelectedTrueEventsNProtons++;
        } else {
            numSelectedEvents0Protons++;
            if (isPionAbsorptionSignal && (numVisibleProtons == 0)) numSelectedTrueEvents0Protons++;
        }

        int numWCMatchProtonDaughters = 0;
        int numWCMatchDaughters = wcMatchDaughtersPDG->size();
        for (int iDaughter = 0; iDaughter < numWCMatchDaughters; ++iDaughter) {
            int daughterPDG = wcMatchDaughtersPDG->at(iDaughter);
            std::string daughterProcess = wcMatchDaughtersProcess->at(iDaughter);

            if (daughterPDG == 2212) numWCMatchProtonDaughters++;
        }
        
        // Flag events
        if (!isPionAbsorptionSignal) {
            EventInfo flaggedEvent;

            flaggedEvent.run    = run; 
            flaggedEvent.subrun = subrun;
            flaggedEvent.event  = event;

            flaggedEvent.wcMatchPDG              = wcMatchPDG;
            flaggedEvent.wcMatchProcess          = *wcMatchProcess;
            flaggedEvent.wcMatchDaughtersPDG     = *wcMatchDaughtersPDG;
            flaggedEvent.wcMatchDaughtersProcess = *wcMatchDaughtersProcess;

            flaggedEvent.truthPrimaryPDG              = truthPrimaryPDG;
            flaggedEvent.truthPrimaryDaughtersPDG     = *truthPrimaryDaughtersPDG;
            flaggedEvent.truthPrimaryDaughtersProcess = *truthPrimaryDaughtersProcess;

            FlaggedEvents.push_back(flaggedEvent);
            
            int numTruthPrimaryDaughters = truthPrimaryDaughtersPDG->size();
            if (truthPrimaryPDG == -211) {
                for (int iDaughter = 0; iDaughter < numTruthPrimaryDaughters; ++iDaughter) {
                    int daughterPDG = truthPrimaryDaughtersPDG->at(iDaughter);
                    std::string daughterProcess = truthPrimaryDaughtersProcess->at(iDaughter);
                    if (!((daughterPDG == 111) || (daughterPDG == 211) || (daughterPDG == -211))) continue;
                    pionBackground[daughterPDG][daughterProcess]++;
                }
            }

            wcMatchPDGCount[wcMatchPDG]++;
            primaryPDGCount[truthPrimaryPDG]++;
        }

        // Fill histograms
        hRecoTrueNumProtons->Fill(numWCMatchProtonDaughters);
        hRecoNumProtonTracks->Fill(protonCount);
    }

    std::ofstream outFile("FlaggedEventsOutput.txt");
    outFile << "Num of flagged events: " << FlaggedEvents.size() << std::endl;
    outFile << std::endl;

    outFile << "WC match PDG count: " << std::endl;
    for (const auto &entry : wcMatchPDGCount) {
        outFile << "    WC match PDG: " << entry.first << " Count: " << entry.second << std::endl;
    }
    outFile << std::endl;

    outFile << "Primary PDG count: " << std::endl;
    for (const auto &entry : primaryPDGCount) {
        outFile << "    True primary PDG: " << entry.first << " Count: " << entry.second << std::endl;
    }
    outFile << std::endl;

    outFile << "Detailed info about primary -211 events:" << std::endl;
    for (const auto& outerPair : pionBackground) {
        outFile << "    Daughter pion: " << outerPair.first << ":\n";
        for (const auto& innerPair : outerPair.second) {
            outFile << "        " << innerPair.first << " -> " << innerPair.second << "\n";
        }
    }
    outFile << std::endl;

    for (const auto &flaggedEvent : FlaggedEvents) {
        outFile << "Run: " << flaggedEvent.run << " subrun: " << flaggedEvent.subrun << " event: " << flaggedEvent.event << std::endl;
        outFile << "WC match PDG: " << flaggedEvent.wcMatchPDG << std::endl;
        outFile << "WC match process: " << flaggedEvent.wcMatchProcess << std::endl;
        outFile << "WC match daughters: " << std::endl;
        for (int n = 0; n < flaggedEvent.wcMatchDaughtersPDG.size(); ++n) {
            if (flaggedEvent.wcMatchDaughtersPDG[n] == 11 && flaggedEvent.wcMatchDaughtersProcess[n] == "hIoni") continue;
            outFile << "    PDG: " << flaggedEvent.wcMatchDaughtersPDG[n];
            outFile << " Process: " << flaggedEvent.wcMatchDaughtersProcess[n];
            outFile << std::endl;
        }
        outFile << "Primary PDG: " << flaggedEvent.truthPrimaryPDG << std::endl;
        outFile << "Primary daughters: " << std::endl;
        for (int n = 0; n < flaggedEvent.truthPrimaryDaughtersPDG.size(); ++n) {
            if (flaggedEvent.truthPrimaryDaughtersPDG[n] == 11 && flaggedEvent.truthPrimaryDaughtersProcess[n] == "hIoni") continue;
            outFile << "    PDG: " << flaggedEvent.truthPrimaryDaughtersPDG[n];
            outFile << " Process: " << flaggedEvent.truthPrimaryDaughtersProcess[n];
            outFile << std::endl;
        }
        outFile << "Tracks reco'ed as protons: " << protonCount << std::endl;
        outFile << std::endl;
    }

    ///////////////////
    // True reco tree
    ///////////////////

    // Loop over true reco events
    Int_t NumTrueRecoEntries = (Int_t) trueRecoTree->GetEntries();
    for (Int_t i = 0; i < NumTrueRecoEntries; ++i) {
        trueRecoTree->GetEntry(i);

        // Fill histograms
        hTruthRecoTrueNumProtons->Fill(trueRecoProtonCount);
    }

    ///////////////////
    // True signal tree
    ///////////////////

    // Keep track of 0p and Np true events
    int numSignalEvents0Protons = 0;
    int numSignalEventsNProtons = 0;

    // Loop over true events
    Int_t NumTrueEntries = (Int_t) trueTree->GetEntries();
    for (Int_t i = 0; i < NumTrueEntries; ++i) {
        trueTree->GetEntry(i);

        if (trueProtonCount >= 1) {
            numSignalEventsNProtons++;   
        } else {
            numSignalEvents0Protons++;
        }
    }

    /////////////////////
    // Drawing histograms
    /////////////////////

    // Setup for drawing plots
    std::vector<int> Colors = {
        kBlack, kBlue, kRed, kGreen
    };

    std::vector<std::vector<TH1*>> PlotGroups = {
        {hRecoTrueNumProtons, hTruthRecoTrueNumProtons},
        {hRecoNumProtonTracks, hRecoTrueNumProtons}
    };

    std::vector<std::vector<TString>> PlotLabelGroups = {
        {"Reco", "True reco"}, 
        {"Reco tracks", "True protons"}
    };

    std::vector<TString> PlotTitles = {
        "TrueNumProtons",
        "RecoProtons"
    };

    std::vector<TString> XLabels = {
        "# of true protons",
        "# of protons"
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
    std::cout << std::endl;
    std::cout << "Number of reco events that pass selection criteria: " << NumEntries << std::endl;
    std::cout << "Number of reco events that pass selection criteria that are true events: " << NumTrueRecoEntries << std::endl;
    std::cout << "Number of background events: " << NumEntries - NumTrueRecoEntries << std::endl;
    std::cout << "Number of true events: " << NumTrueEntries << std::endl;
    std::cout << "Efficiency: " << 100 * ((float)NumTrueRecoEntries / (float)NumTrueEntries) << "%" << std::endl; 
    std::cout << "Purity:     " << 100 * ((float)NumTrueRecoEntries / (float)NumEntries) << "%" << std::endl;
    std::cout << std::endl;
    std::cout << "Selected 0p events: " << numSelectedEvents0Protons << std::endl;
    std::cout << "Selected 0p events that are true 0p events: " << numSelectedTrueEvents0Protons << std::endl;
    std::cout << "Signal 0p events: " << numSignalEvents0Protons << std::endl;
    std::cout << "Efficiency 0p: " << 100 * ((float) numSelectedTrueEvents0Protons / (float) numSignalEvents0Protons) << "%" << std::endl; 
    std::cout << "Purity 0p:     " << 100 * ((float) numSelectedTrueEvents0Protons / (float) numSelectedEvents0Protons) << "%" << std::endl; 
    std::cout << std::endl;
    std::cout << "Selected Np events: " << numSelectedEventsNProtons << std::endl;
    std::cout << "Selected Np events that are true Np events: " << numSelectedTrueEventsNProtons << std::endl;
    std::cout << "Signal Np events: " << numSignalEventsNProtons << std::endl;
    std::cout << "Efficiency Np: " << 100 * ((float) numSelectedTrueEventsNProtons / (float) numSignalEventsNProtons) << "%" << std::endl; 
    std::cout << "Purity Np:     " << 100 * ((float) numSelectedTrueEventsNProtons / (float) numSelectedEventsNProtons) << "%" << std::endl; 
}