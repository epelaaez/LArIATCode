#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TEfficiency.h>

#include <vector>
#include <map>

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
    {11, "Other"}
};

struct EventInfo {
    int run; int subrun; int event;
    float vertexX; float vertexY; float vertexZ;

    int                      wcMatchPDG;
    std::string              wcMatchProcess;
    std::vector<int>         wcMatchDaughtersPDG;
    std::vector<std::string> wcMatchDaughtersProcess;

    int                      truthPrimaryPDG;
    std::vector<int>         truthPrimaryDaughtersPDG;
    std::vector<std::string> truthPrimaryDaughtersProcess;

    int                      visibleProtons;
    int                      recoProtonCount;
    int                      backgroundNum;
};

void printBackgroundInfo(TH1D* background_histo, std::ostream& os);
void printEventInfo(EventInfo event, std::ostream& os);

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

    // Histgrams for event counting
    TH1D* hTotalEvents                 = (TH1D*) Directory->Get<TH1D>("hTotalEvents");
    TH1D* hTotalEvents0pSignal         = (TH1D*) Directory->Get<TH1D>("hTotalEvents0pSignal");
    TH1D* hTotalEventsNpSignal         = (TH1D*) Directory->Get<TH1D>("hTotalEventsNpSignal");
    TH1D* hTotalBackground             = (TH1D*) Directory->Get<TH1D>("hTotalBackground");
    TH1D* hWCExists                    = (TH1D*) Directory->Get<TH1D>("hWCExists");
    TH1D* hWCExists0pSignal            = (TH1D*) Directory->Get<TH1D>("hWCExists0pSignal");
    TH1D* hWCExistsNpSignal            = (TH1D*) Directory->Get<TH1D>("hWCExistsNpSignal");
    TH1D* hWCExistsBackground          = (TH1D*) Directory->Get<TH1D>("hWCExistsBackground");
    TH1D* hPionInRedVolume             = (TH1D*) Directory->Get<TH1D>("hPionInRedVolume");
    TH1D* hPionInRedVolume0pSignal     = (TH1D*) Directory->Get<TH1D>("hPionInRedVolume0pSignal");
    TH1D* hPionInRedVolumeNpSignal     = (TH1D*) Directory->Get<TH1D>("hPionInRedVolumeNpSignal");
    TH1D* hPionInRedVolumeBackground   = (TH1D*) Directory->Get<TH1D>("hPionInRedVolumeBackground");
    TH1D* hNoOutgoingPion              = (TH1D*) Directory->Get<TH1D>("hNoOutgoingPion");
    TH1D* hNoOutgoingPion0pSignal      = (TH1D*) Directory->Get<TH1D>("hNoOutgoingPion0pSignal");
    TH1D* hNoOutgoingPionNpSignal      = (TH1D*) Directory->Get<TH1D>("hNoOutgoingPionNpSignal");
    TH1D* hNoOutgoingPion0pBackground  = (TH1D*) Directory->Get<TH1D>("hNoOutgoingPion0pBackground");
    TH1D* hNoOutgoingPionNpBackground  = (TH1D*) Directory->Get<TH1D>("hNoOutgoingPionNpBackground");
    TH1D* hSmallTracks                 = (TH1D*) Directory->Get<TH1D>("hSmallTracks");
    TH1D* hSmallTracks0pSignal         = (TH1D*) Directory->Get<TH1D>("hSmallTracks0pSignal");
    TH1D* hSmallTracksNpSignal         = (TH1D*) Directory->Get<TH1D>("hSmallTracksNpSignal");
    TH1D* hSmallTracks0pBackground     = (TH1D*) Directory->Get<TH1D>("hSmallTracks0pBackground");
    TH1D* hSmallTracksNpBackground     = (TH1D*) Directory->Get<TH1D>("hSmallTracksNpBackground");
    TH1D* hMeanCurvature               = (TH1D*) Directory->Get<TH1D>("hMeanCurvature");
    TH1D* hMeanCurvature0pSignal       = (TH1D*) Directory->Get<TH1D>("hMeanCurvature0pSignal");
    TH1D* hMeanCurvatureNpSignal       = (TH1D*) Directory->Get<TH1D>("hMeanCurvatureNpSignal");
    TH1D* hMeanCurvature0pBackground   = (TH1D*) Directory->Get<TH1D>("hMeanCurvature0pBackground");
    TH1D* hMeanCurvatureNpBackground   = (TH1D*) Directory->Get<TH1D>("hMeanCurvatureNpBackground");
    TH1D* hFinalRecoEvents             = (TH1D*) Directory->Get<TH1D>("hFinalRecoEvents");
    TH1D* hFinalRecoEvents0pSignal     = (TH1D*) Directory->Get<TH1D>("hFinalRecoEvents0pSignal");
    TH1D* hFinalRecoEventsNpSignal     = (TH1D*) Directory->Get<TH1D>("hFinalRecoEventsNpSignal");
    TH1D* hFinalRecoEvents0pBackground = (TH1D*) Directory->Get<TH1D>("hFinalRecoEvents0pBackground");
    TH1D* hFinalRecoEventsNpBackground = (TH1D*) Directory->Get<TH1D>("hFinalRecoEventsNpBackground");

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
    float                     truthPrimaryVertexX;
    float                     truthPrimaryVertexY;
    float                     truthPrimaryVertexZ;
    bool                      isPionAbsorptionSignal;
    int                       numVisibleProtons;
    int                       backgroundType;
    tree->SetBranchAddress("truthPrimaryPDG", &truthPrimaryPDG);
    tree->SetBranchAddress("truthPrimaryDaughtersPDG", &truthPrimaryDaughtersPDG);
    tree->SetBranchAddress("truthPrimaryDaughtersProcess", &truthPrimaryDaughtersProcess);
    tree->SetBranchAddress("truthPrimaryDaughtersKE", &truthPrimaryDaughtersKE);
    tree->SetBranchAddress("isPionAbsorptionSignal", &isPionAbsorptionSignal);
    tree->SetBranchAddress("numVisibleProtons", &numVisibleProtons);
    tree->SetBranchAddress("truthPrimaryVertexX", &truthPrimaryVertexX);
    tree->SetBranchAddress("truthPrimaryVertexY", &truthPrimaryVertexY);
    tree->SetBranchAddress("truthPrimaryVertexZ", &truthPrimaryVertexZ);
    tree->SetBranchAddress("backgroundType", &backgroundType);

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

    // We will keep track of weird events
    std::vector<EventInfo> FlaggedEvents;
    std::vector<EventInfo> Reco0pBackgroundNp;
    std::vector<EventInfo> RecoNpBackground0p;
    std::map<int, int> wcMatchPDGCount;
    std::map<int, int> primaryPDGCount;

    // Loop over reco events
    Int_t NumEntries = (Int_t) tree->GetEntries();
    for (Int_t i = 0; i < NumEntries; ++i) {
        tree->GetEntry(i);

        bool is0pRecoNpBackground = false;
        bool isNpReco0pBackground = false;

        int numWCMatchProtonDaughters = 0;
        int numWCMatchDaughters = wcMatchDaughtersPDG->size();
        for (int iDaughter = 0; iDaughter < numWCMatchDaughters; ++iDaughter) {
            int daughterPDG = wcMatchDaughtersPDG->at(iDaughter);
            std::string daughterProcess = wcMatchDaughtersProcess->at(iDaughter);

            if (daughterPDG == 2212) numWCMatchProtonDaughters++;
        }

        if ((isPionAbsorptionSignal) && (numVisibleProtons == 0) && (protonCount > 0)) isNpReco0pBackground = true;
        if ((isPionAbsorptionSignal) && (numVisibleProtons > 0) && (protonCount == 0)) is0pRecoNpBackground = true;
        
        // Flag events
        if ((!isPionAbsorptionSignal) || isNpReco0pBackground || is0pRecoNpBackground) {
            EventInfo flaggedEvent;

            flaggedEvent.run    = run; 
            flaggedEvent.subrun = subrun;
            flaggedEvent.event  = event;

            flaggedEvent.vertexX = truthPrimaryVertexX;
            flaggedEvent.vertexY = truthPrimaryVertexY;
            flaggedEvent.vertexZ = truthPrimaryVertexZ;

            flaggedEvent.wcMatchPDG              = wcMatchPDG;
            flaggedEvent.wcMatchProcess          = *wcMatchProcess;
            flaggedEvent.wcMatchDaughtersPDG     = *wcMatchDaughtersPDG;
            flaggedEvent.wcMatchDaughtersProcess = *wcMatchDaughtersProcess;

            flaggedEvent.truthPrimaryPDG              = truthPrimaryPDG;
            flaggedEvent.truthPrimaryDaughtersPDG     = *truthPrimaryDaughtersPDG;
            flaggedEvent.truthPrimaryDaughtersProcess = *truthPrimaryDaughtersProcess;

            flaggedEvent.visibleProtons  = numVisibleProtons;
            flaggedEvent.backgroundNum   = backgroundType;
            flaggedEvent.recoProtonCount = protonCount;

            if (!isPionAbsorptionSignal) FlaggedEvents.push_back(flaggedEvent);
            if (is0pRecoNpBackground)    Reco0pBackgroundNp.push_back(flaggedEvent);
            if (isNpReco0pBackground)    RecoNpBackground0p.push_back(flaggedEvent);

            int numTruthPrimaryDaughters = truthPrimaryDaughtersPDG->size();
            if (truthPrimaryPDG == -211) {
                for (int iDaughter = 0; iDaughter < numTruthPrimaryDaughters; ++iDaughter) {
                    int daughterPDG = truthPrimaryDaughtersPDG->at(iDaughter);
                    std::string daughterProcess = truthPrimaryDaughtersProcess->at(iDaughter);
                    if (!((daughterPDG == 111) || (daughterPDG == 211) || (daughterPDG == -211))) continue;
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

    for (const auto &flaggedEvent : FlaggedEvents) {
        printEventInfo(flaggedEvent, outFile);
    }

    std::ofstream out0pFile("0pRecoNpBackground.txt");
    out0pFile << "Events reconstructed as 0p absorption but Np absorption at truth level: " << Reco0pBackgroundNp.size() << std::endl;
    out0pFile << std::endl;
    for (const auto &flaggedEvent : Reco0pBackgroundNp) {
        printEventInfo(flaggedEvent, out0pFile);
    }

    std::ofstream outNpFile("NpReco0pBackground.txt");
    outNpFile << "Events reconstructed as Np absorption but 0p absorption at truth level: " << RecoNpBackground0p.size() << std::endl;
    outNpFile << std::endl;
    for (const auto &flaggedEvent : RecoNpBackground0p) {
        printEventInfo(flaggedEvent, outNpFile);
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
    int total0pSignalEvents = hTotalEvents0pSignal->Integral();
    int totalNpSignalEvents = hTotalEventsNpSignal->Integral();
    int totalSignalEvents   = total0pSignalEvents + totalNpSignalEvents;
    std::cout << "Total events: " << hTotalEvents->Integral() << std::endl;
    std::cout << "Total signal events: " << totalSignalEvents << std::endl;
    std::cout << "Total 0p signal events: " << total0pSignalEvents << std::endl;
    std::cout << "Total Np signal events: " << totalNpSignalEvents << std::endl;
    std::cout << std::endl;
    std::cout << "Cut statistics:" << std::endl;
    std::cout << "  WC to TPC match: " << std::endl;
    std::cout << "    Total reco events: " << hWCExists->Integral() << std::endl;
    std::cout << "    0p reco signal events: " << hWCExists0pSignal->Integral() << std::endl;
    std::cout << "    Np reco signal events: " << hWCExistsNpSignal->Integral() << std::endl;
    std::cout << "    Overall purity: " << (hWCExists0pSignal->Integral() + hWCExistsNpSignal->Integral()) / hWCExists->Integral() << " and efficiency: " << (hWCExists0pSignal->Integral() + hWCExistsNpSignal->Integral()) / totalSignalEvents << std::endl;
    std::cout << std::endl;
    std::cout << "  Pion vertex in reduced volume: " << std::endl;
    std::cout << "    Total reco events: " << hPionInRedVolume->Integral() << std::endl;
    std::cout << "    0p signal events: " << hPionInRedVolume0pSignal->Integral() << std::endl;
    std::cout << "    Np signal events: " << hPionInRedVolumeNpSignal->Integral() << std::endl;
    std::cout << "    Overall purity: " << (hPionInRedVolume0pSignal->Integral() + hPionInRedVolumeNpSignal->Integral()) / hPionInRedVolume->Integral() << " and efficiency: " << (hPionInRedVolume0pSignal->Integral() + hPionInRedVolumeNpSignal->Integral()) / totalSignalEvents << std::endl;
    std::cout << std::endl;
    std::cout << "  No outgoing pion: " << std::endl;
    std::cout << "    Total reco events: " << hNoOutgoingPion->Integral() << std::endl;
    std::cout << "    0p signal events: " << hNoOutgoingPion0pSignal->Integral() << std::endl;
    std::cout << "    Np signal events: " <<  hNoOutgoingPionNpSignal->Integral() << std::endl;
    std::cout << "    Overall purity: " << (hNoOutgoingPion0pSignal->Integral() + hNoOutgoingPionNpSignal->Integral()) / hNoOutgoingPion->Integral() << " and efficiency: " << (hNoOutgoingPion0pSignal->Integral() + hNoOutgoingPionNpSignal->Integral()) / totalSignalEvents << std::endl;
    std::cout << std::endl;
    std::cout << "    Reco 0p events: " << hNoOutgoingPion->Integral(1, 1) << std::endl;
    std::cout << "    Reco 0p true signal events: " << hNoOutgoingPion0pSignal->Integral(1, 1) << std::endl;
    std::cout << "    0p purity: " << hNoOutgoingPion0pSignal->Integral(1, 1) / hNoOutgoingPion->Integral(1,1) << " and efficiency: " << hNoOutgoingPion0pSignal->Integral(1, 1) / total0pSignalEvents << std::endl;
    std::cout << std::endl;
    std::cout << "    Reco Np events: " << hNoOutgoingPion->Integral(2, 2) << std::endl;
    std::cout << "    Reco Np true signal events: " << hNoOutgoingPionNpSignal->Integral(2, 2) << std::endl;
    std::cout << "    Np purity: " << hNoOutgoingPionNpSignal->Integral(2, 2) / hNoOutgoingPion->Integral(2,2) << " and efficiency: " << hNoOutgoingPionNpSignal->Integral(2, 2) / totalNpSignalEvents << std::endl;
    std::cout << std::endl;
    std::cout << "  No small tracks: " << std::endl;
    std::cout << "    Total reco events: " << hSmallTracks->Integral() << std::endl;
    std::cout << "    0p signal events: " << hSmallTracks0pSignal->Integral() << std::endl;
    std::cout << "    Np signal events: " <<  hSmallTracksNpSignal->Integral() << std::endl;
    std::cout << "    Overall purity: " << (hSmallTracks0pSignal->Integral() + hSmallTracksNpSignal->Integral()) / hSmallTracks->Integral() << " and efficiency: " << (hSmallTracks0pSignal->Integral() + hSmallTracksNpSignal->Integral()) / totalSignalEvents << std::endl;
    std::cout << std::endl;
    std::cout << "    Reco 0p events: " << hSmallTracks->Integral(1, 1) << std::endl;
    std::cout << "    Reco 0p true signal events: " << hSmallTracks0pSignal->Integral(1,1) << std::endl;
    std::cout << "    0p purity: " << hSmallTracks0pSignal->Integral(1, 1) / hSmallTracks->Integral(1, 1) << " and efficiency: " << hSmallTracks0pSignal->Integral(1, 1) / total0pSignalEvents << std::endl;
    std::cout << std::endl;
    std::cout << "    Reco Np events: " << hSmallTracks->Integral(2, 2) << std::endl;
    std::cout << "    Reco Np true signal events: " << hSmallTracksNpSignal->Integral(2,2) << std::endl;
    std::cout << "    Np purity: " << hSmallTracksNpSignal->Integral(2, 2) / hSmallTracks->Integral(2, 2) << " and efficiency: " << hSmallTracksNpSignal->Integral(2, 2) / totalNpSignalEvents << std::endl;
    std::cout << std::endl;
    std::cout << "  Curvature: " << std::endl;
    std::cout << "    Total reco events: " << hMeanCurvature->Integral() << std::endl;
    std::cout << "    0p signal events: " << hMeanCurvature0pSignal->Integral() << std::endl;
    std::cout << "    Np signal events: " <<  hMeanCurvatureNpSignal->Integral() << std::endl;
    std::cout << "    Overall purity: " << (hMeanCurvature0pSignal->Integral() + hMeanCurvatureNpSignal->Integral()) / hMeanCurvature->Integral() << " and efficiency: " << (hMeanCurvature0pSignal->Integral() + hMeanCurvatureNpSignal->Integral()) / totalSignalEvents << std::endl;
    std::cout << std::endl;
    std::cout << "    Reco 0p events: " << hMeanCurvature->Integral(1, 1) << std::endl;
    std::cout << "    Reco 0p true signal events: " << hMeanCurvature0pSignal->Integral(1,1) << std::endl;
    std::cout << "    0p purity: " << hMeanCurvature0pSignal->Integral(1, 1) / hMeanCurvature->Integral(1, 1) << " and efficiency: " << hMeanCurvature0pSignal->Integral(1, 1) / total0pSignalEvents << std::endl;
    std::cout << std::endl;
    std::cout << "    Reco Np events: " << hMeanCurvature->Integral(2, 2) << std::endl;
    std::cout << "    Reco Np true signal events: " << hMeanCurvatureNpSignal->Integral(2,2) << std::endl;
    std::cout << "    Np purity: " << hMeanCurvatureNpSignal->Integral(2, 2) / hMeanCurvature->Integral(2, 2) << " and efficiency: " << hMeanCurvatureNpSignal->Integral(2, 2) / totalNpSignalEvents << std::endl;
    std::cout << std::endl;

    int totalFinalSignal0p = hFinalRecoEvents0pSignal->Integral();
    int totalFinalSignalNp = hFinalRecoEventsNpSignal->Integral();
    int totalFinalSignal   = totalFinalSignal0p + totalFinalSignalNp;
    std::cout << "Final stats: " << std::endl;
    std::cout << "  Reco events: " << hFinalRecoEvents->Integral() << std::endl;
    std::cout << "  Overall purity: " << totalFinalSignal / hFinalRecoEvents->Integral() << " and efficiency: " << (double)totalFinalSignal / totalSignalEvents << std::endl;
    std::cout << "  Reco 0p events: " << hFinalRecoEvents->Integral(1,1) << std::endl;
    std::cout << "  0p channel purity: " << hFinalRecoEvents0pSignal->Integral(1,1) / hFinalRecoEvents->Integral(1,1) << " and efficiency: " << hFinalRecoEvents0pSignal->Integral(1,1) / total0pSignalEvents << std::endl;
    std::cout << "  Reco Np events: " << hFinalRecoEvents->Integral(2,2) << std::endl;
    std::cout << "  0p channel purity: " << hFinalRecoEventsNpSignal->Integral(2,2) / hFinalRecoEvents->Integral(2,2) << " and efficiency: " << hFinalRecoEventsNpSignal->Integral(2,2) / totalNpSignalEvents << std::endl;
    std::cout << std::endl;

    std::cout << "Reco 0p background composition: " << std::endl;
    printBackgroundInfo(hFinalRecoEvents0pBackground, std::cout);
    std::cout << std::endl;

    std::cout << "Reco Np background composition: " << std::endl;
    printBackgroundInfo(hFinalRecoEventsNpBackground, std::cout);
    std::cout << std::endl;

    std::ofstream outFileBackgroundBreakdown("BackgroundBreakdown.txt");
    outFileBackgroundBreakdown << "Overall background after WC existence cut:" << std::endl;
    printBackgroundInfo(hWCExistsBackground, outFileBackgroundBreakdown);
    outFileBackgroundBreakdown << std::endl;
    
    outFileBackgroundBreakdown << "Overall background after pion in reduced volume cut:" << std::endl;
    printBackgroundInfo(hPionInRedVolumeBackground, outFileBackgroundBreakdown);
    outFileBackgroundBreakdown << std::endl;

    outFileBackgroundBreakdown << "Reco 0p background after no outgoing pion cut:" << std::endl;
    printBackgroundInfo(hNoOutgoingPion0pBackground, outFileBackgroundBreakdown);
    outFileBackgroundBreakdown << std::endl;
    outFileBackgroundBreakdown << "Reco Np background after no outgoing pion cut:" << std::endl;
    printBackgroundInfo(hNoOutgoingPionNpBackground, outFileBackgroundBreakdown);
    outFileBackgroundBreakdown << std::endl;

    outFileBackgroundBreakdown << "Reco 0p background after small tracks cut:" << std::endl;
    printBackgroundInfo(hSmallTracks0pBackground, outFileBackgroundBreakdown);
    outFileBackgroundBreakdown << std::endl;
    outFileBackgroundBreakdown << "Reco Np background after small tracks cut:" << std::endl;
    printBackgroundInfo(hSmallTracksNpBackground, outFileBackgroundBreakdown);
    outFileBackgroundBreakdown << std::endl;

    outFileBackgroundBreakdown << "Reco 0p background after mean curvature cut:" << std::endl;
    printBackgroundInfo(hMeanCurvature0pBackground, outFileBackgroundBreakdown);
    outFileBackgroundBreakdown << std::endl;
    outFileBackgroundBreakdown << "Reco Np background after mean curvature cut:" << std::endl;
    printBackgroundInfo(hMeanCurvatureNpBackground, outFileBackgroundBreakdown);
    outFileBackgroundBreakdown << std::endl;
}

void printBackgroundInfo(TH1D* background_histo, std::ostream& os) {
    int pionAbs0p            = background_histo->GetBinContent(1);
    int pionAbsNp            = background_histo->GetBinContent(2);
    int primaryMuon          = background_histo->GetBinContent(3);
    int primaryElectron      = background_histo->GetBinContent(4);
    int otherPrimary         = background_histo->GetBinContent(5);
    int pionOutRedVol        = background_histo->GetBinContent(6);
    int pionInelScatter      = background_histo->GetBinContent(7);
    int chargeExchange       = background_histo->GetBinContent(8);
    int doubleChargeExchange = background_histo->GetBinContent(9);
    int captureAtRest        = background_histo->GetBinContent(10);
    int pionDecay            = background_histo->GetBinContent(11);
    int other                = background_histo->GetBinContent(12);

    os << "  Abs 0p: " << pionAbs0p << std::endl;;
    os << "  Abs Np: " << pionAbsNp << std::endl;;
    os << "  Primary muon: " << primaryMuon << std::endl;;
    os << "  Primary electron: " << primaryElectron << std::endl;;
    os << "  Other primary: " << otherPrimary << std::endl;;
    os << "  Outside reduced volume: " << pionOutRedVol << std::endl;;
    os << "  Inelastic scattering: " << pionInelScatter << std::endl;;
    os << "  Charge exchange: " << chargeExchange << std::endl;;
    os << "  Double charge exchange: " << doubleChargeExchange << std::endl;;
    os << "  Capture at rest: " << captureAtRest << std::endl;;
    os << "  Decay: " << pionDecay << std::endl;;
    os << "  Other: " << other << std::endl;
}

void printEventInfo(EventInfo event, std::ostream& os) {
    os << "Run: " << event.run << " subrun: " << event.subrun << " event: " << event.event << std::endl;
    os << "Background type: " << backgroundTypes[event.backgroundNum] << std::endl;
    os << "WC match PDG: " << event.wcMatchPDG << std::endl;
    os << "WC match process: " << event.wcMatchProcess << std::endl;
    os << "WC match daughters: " << std::endl;
    for (int n = 0; n < event.wcMatchDaughtersPDG.size(); ++n) {
        if (event.wcMatchDaughtersPDG[n] == 11 && event.wcMatchDaughtersProcess[n] == "hIoni") continue;
        os << "    PDG: " << event.wcMatchDaughtersPDG[n];
        os << " Process: " << event.wcMatchDaughtersProcess[n];
        os << std::endl;
    }
    os << "Primary PDG: " << event.truthPrimaryPDG << std::endl;
    os << "Primary vertex: " << event.vertexX << ", " << event.vertexY << ", " << event.vertexZ << std::endl; 
    os << "Primary daughters: " << std::endl;
    for (int n = 0; n < event.truthPrimaryDaughtersPDG.size(); ++n) {
        if (event.truthPrimaryDaughtersPDG[n] == 11 && event.truthPrimaryDaughtersProcess[n] == "hIoni") continue;
        os << "    PDG: " << event.truthPrimaryDaughtersPDG[n];
        os << " Process: " << event.truthPrimaryDaughtersProcess[n];
        os << std::endl;
    }
    os << "Tracks reco'ed as protons: " << event.recoProtonCount << std::endl;
    os << "True visible protons: " << event.visibleProtons << std::endl;
    os << std::endl;
}