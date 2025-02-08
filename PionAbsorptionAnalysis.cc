#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TTree.h>

#include <vector>

void PionAbsorptionAnalysis() {
    // Set defaults
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    // gStyle->SetOptStat(0);

    int FontStyle = 132;
    double TextSize = 0.06;

    // Load root file
    TString RootFilePath = "/exp/lariat/app/users/epelaez/test_dev/srcs/lariatsoft/LArIATFilterModule/SelectPionAbsorption_histo.root"; 
    std::unique_ptr<TFile> File(TFile::Open(RootFilePath));
    TDirectory* Directory = (TDirectory*)File->Get("SelectPionAbsorption");

    // Load tree and set branch addresses
    TTree* tree = (TTree*) Directory->Get<TTree>("PionAbsorptionTree");

    int run, subrun, event;
    tree->SetBranchAddress("run", &run); 
    tree->SetBranchAddress("subrun", &subrun); 
    tree->SetBranchAddress("event", &event);

    int iNumNucleons, iNumNeutrons, iNumProtons;
    tree->SetBranchAddress("iNumNucleons", &iNumNucleons); 
    tree->SetBranchAddress("iNumNeutrons", &iNumNeutrons);
    tree->SetBranchAddress("iNumProtons", &iNumProtons);

    float fPionInitialEnergy, fPionInitialKEnergy, fPionVertexEnergy, fPionVertexKEnergy;
    tree->SetBranchAddress("fPionInitialEnergy", &fPionInitialEnergy);
    tree->SetBranchAddress("fPionInitialKEnergy", &fPionInitialKEnergy);
    tree->SetBranchAddress("fPionVertexEnergy", &fPionVertexEnergy);
    tree->SetBranchAddress("fPionVertexKEnergy", &fPionVertexKEnergy);

    float fPionVertexMomentum, fPionVertexPx, fPionVertexPy, fPionVertexPz, fPionTrackLength;
    tree->SetBranchAddress("fPionVertexMomentum", &fPionVertexMomentum);
    tree->SetBranchAddress("fPionVertexPx", &fPionVertexPx);
    tree->SetBranchAddress("fPionVertexPy", &fPionVertexPy);
    tree->SetBranchAddress("fPionVertexPz", &fPionVertexPz);
    tree->SetBranchAddress("fPionTrackLength", &fPionTrackLength);

    float fTotalOutEnergy, fTotalOutKEnergy, fTotalOutMomentum, fTotalOutPx, fTotalOutPy, fTotalOutPz;
    tree->SetBranchAddress("fTotalOutEnergy", &fTotalOutEnergy);
    tree->SetBranchAddress("fTotalOutKEnergy", &fTotalOutKEnergy);
    tree->SetBranchAddress("fTotalOutMomentum", &fTotalOutMomentum);
    tree->SetBranchAddress("fTotalOutPx", &fTotalOutPx);
    tree->SetBranchAddress("fTotalOutPy", &fTotalOutPy);
    tree->SetBranchAddress("fTotalOutPz", &fTotalOutPz);

    float fLeadingOutMomentum, fDiffInOutMomentum, fDiffInOutPx, fDiffInOutPy, fDiffInOutPz, fDiffInEnergyOutKEnergy;
    tree->SetBranchAddress("fLeadingOutMomentum", &fLeadingOutMomentum);
    tree->SetBranchAddress("fDiffInOutMomentum", &fDiffInOutMomentum);
    tree->SetBranchAddress("fDiffInOutPx", &fDiffInOutPx);
    tree->SetBranchAddress("fDiffInOutPy", &fDiffInOutPy);
    tree->SetBranchAddress("fDiffInOutPz", &fDiffInOutPz);
    tree->SetBranchAddress("fDiffInEnergyOutKEnergy", &fDiffInEnergyOutKEnergy);

    float fCosAngleInOutMomentum;
    tree->SetBranchAddress("fCosAngleInOutMomentum", &fCosAngleInOutMomentum);

    // Declare histograms
    TH1I *hNumNucleons = new TH1I("hNumNucleons", "NumNucleons;# of nucleons;", 20, 0, 20);
    TH1I *hNumNeutrons = new TH1I("hNumNeutrons", "NumNeutrons;# of neutrons", 15, 0, 15);
    TH1I *hNumProtons = new TH1I("hNumProtons", "NumProtons;# of protons", 10, 0, 10);

    TH1F *hPionInitialEnergy = new TH1F("hPionInitialEnergy", "PionInitialEnergy;[GeV];", 20, -1, -1);
    TH1F *hPionInitialKEnergy = new TH1F("hPionInitialKEnergy", "PionInitialKEnergy;[GeV];", 20, -1, -1);
    TH1F *hPionVertexEnergy = new TH1F("hPionVertexEnergy", "PionVertexEnergy;[GeV];", 20, -1, -1);
    TH1F *hPionVertexKEnergy = new TH1F("hPionVertexKEnergy", "PionVertexKEnergy;[GeV];", 20, -1, -1);

    TH1F *hPionVertexMomentum = new TH1F("hPionVertexMomentum", "PionVertexMomentum;[GeV/c];", 20, -1, -1);
    TH1F *hPionVertexPx = new TH1F("hPionVertexPx", "PionVertexPx;[GeV/c];", 20, -1, -1);
    TH1F *hPionVertexPy = new TH1F("hPionVertexPy", "PionVertexPy;[GeV/c];", 20, -1, -1);
    TH1F *hPionVertexPz = new TH1F("hPionVertexPz", "PionVertexPz;[GeV/c];", 20, -1, -1);
    TH1F *hPionTrackLength = new TH1F("hPionTrackLength", "hPionTrackLength;[cm];", 20, -1, -1);

    TH1F *hTotalOutEnergy = new TH1F("hTotalOutEnergy", "TotalOutEnergy;[GeV];", 20, -1, -1);
    TH1F *hTotalOutKEnergy = new TH1F("hTotalOutKEnergy", "TotalOutKEnergy;[GeV];", 20, -1, -1);

    TH1F *hTotalOutMomentum = new TH1F("hTotalOutMomentum", "TotalOuMomentum;[GeV/c];", 20, -1, -1);
    TH1F *hTotalOutPx = new TH1F("hTotalOutPx", "TotalOutPx;[GeV/c];", 20, -1, -1);
    TH1F *hTotalOutPy = new TH1F("hTotalOutPy", "TotalOutPy;[GeV/c];", 20, -1, -1);
    TH1F *hTotalOutPz = new TH1F("hTotalOutPz", "TotalOutPz;[GeV/c];", 20, -1, -1);

    TH1F *hDiffInOutMomentum = new TH1F("hDiffInOutMomentum", "DiffInOutMomentum;[GeV/c];", 20, -1, -1);
    TH1F *hDiffInOutPx = new TH1F("hDiffInOutPx", "DiffInOutPx;[GeV/c];", 20, -1, -1);
    TH1F *hDiffInOutPy = new TH1F("hDiffInOutPy", "DiffInOutPy;[GeV/c];", 20, -1, -1);
    TH1F *hDiffInOutPz = new TH1F("hDiffInOutPz", "DiffInOutPz;[GeV/c];", 20, -1, -1);
    TH1F *hDiffInEnergyOutKEnergy = new TH1F("hDiffInEnergyOutKEnergy", "DiffInEnergyOutKEnergy;[GeV/c];", 20, -1, -1);

    TH1F *hCosAngleInOutMomentum = new TH1F("hCosAngleInOutMomentum", "CosAngleInOutMomentum;[unitless];",20,-1,-1);

    std::vector<TH1*> AllHistograms = {
        hNumNucleons,
        hNumNeutrons,
        hNumProtons,
        hPionInitialEnergy,
        hPionInitialKEnergy,
        hPionVertexEnergy,
        hPionVertexKEnergy,
        hPionVertexMomentum,
        hPionVertexPx,
        hPionVertexPy,
        hPionVertexPz,
        hTotalOutEnergy,
        hTotalOutKEnergy,
        hTotalOutMomentum,
        hTotalOutPx,
        hTotalOutPy,
        hTotalOutPz,
        hDiffInOutMomentum,
        hDiffInOutPx,
        hDiffInOutPy,
        hDiffInOutPz,
        hDiffInEnergyOutKEnergy,
        hCosAngleInOutMomentum
    };

    Int_t NumEntries = (Int_t) tree->GetEntries();
    for (Int_t i = 0; i < NumEntries; ++i) {
        tree->GetEntry(i);

        // Look at interesting events
        if (fPionVertexPz < 0) {
            std::cout << "run: " << run;
            std::cout << " subrun: " << subrun;
            std::cout << " event: " << event << std::endl;
        }

        // Fill histograms
        hNumNucleons->Fill(iNumNucleons);
        hNumNeutrons->Fill(iNumNeutrons);
        hNumProtons->Fill(iNumProtons);

        hPionInitialEnergy->Fill(fPionInitialEnergy);
        hPionInitialKEnergy->Fill(fPionInitialKEnergy);
        hPionVertexEnergy->Fill(fPionVertexEnergy);
        hPionVertexKEnergy->Fill(fPionVertexKEnergy);

        hPionVertexMomentum->Fill(fPionVertexMomentum);
        hPionVertexPx->Fill(fPionVertexPx);
        hPionVertexPy->Fill(fPionVertexPy);
        hPionVertexPz->Fill(fPionVertexPz);
        hPionTrackLength->Fill(fPionTrackLength);

        hTotalOutEnergy->Fill(fTotalOutEnergy);
        hTotalOutKEnergy->Fill(fTotalOutKEnergy);

        hTotalOutMomentum->Fill(fTotalOutMomentum);
        hTotalOutPx->Fill(fTotalOutPx);
        hTotalOutPy->Fill(fTotalOutPy);
        hTotalOutPz->Fill(fTotalOutPz);

        hDiffInOutMomentum->Fill(fDiffInOutMomentum);
        hDiffInOutPx->Fill(fDiffInOutPx);
        hDiffInOutPy->Fill(fDiffInOutPy);
        hDiffInOutPz->Fill(fDiffInOutPz);
        hDiffInEnergyOutKEnergy->Fill(fDiffInEnergyOutKEnergy);

        hCosAngleInOutMomentum->Fill(fCosAngleInOutMomentum);
    }    

    // TODO: Reweight needed? (Maybe ask Afro)

    // Set up canvas for drawing
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis/figs/";
    TCanvas* PlotCanvas = new TCanvas("Canvas","Canvas",205,34,1024,768);
    PlotCanvas->cd();
    PlotCanvas->SetLeftMargin(0.15);
    PlotCanvas->SetBottomMargin(0.15);

    for (int iHisto = 0; iHisto < AllHistograms.size(); ++iHisto) {
        PlotCanvas->cd();
        TH1* CurrentHisto = AllHistograms[iHisto];

        // std::cout << "Integral:" << CurrentHisto->Integral() << std::endl;

        CurrentHisto->SetLineWidth(2);
        CurrentHisto->SetLineColor(kBlack);

        CurrentHisto->GetXaxis()->SetTitleFont(FontStyle);
        CurrentHisto->GetXaxis()->SetLabelFont(FontStyle);
        CurrentHisto->GetXaxis()->SetNdivisions(8);
        CurrentHisto->GetXaxis()->SetLabelSize(TextSize);
        CurrentHisto->GetXaxis()->SetTitleSize(TextSize);	
        CurrentHisto->GetXaxis()->SetTitleOffset(1.1);					
        CurrentHisto->GetXaxis()->CenterTitle();						

        CurrentHisto->GetYaxis()->SetTitleFont(FontStyle);
        CurrentHisto->GetYaxis()->SetLabelFont(FontStyle);
        CurrentHisto->GetYaxis()->SetNdivisions(6);
        CurrentHisto->GetYaxis()->SetLabelSize(TextSize);
        CurrentHisto->GetYaxis()->SetTitleSize(TextSize);
        CurrentHisto->GetYaxis()->SetTitle("Event counts");
        CurrentHisto->GetYaxis()->SetTitleOffset(1.1);
        CurrentHisto->GetYaxis()->CenterTitle();	

        CurrentHisto->Draw("hist");
        PlotCanvas->SaveAs(SaveDir + CurrentHisto->GetTitle() + ".png");
    } 
}