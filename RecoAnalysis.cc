#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TEfficiency.h>

#include <vector>

void RecoAnalysis() {
    // Set defaults
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    // gStyle->SetOptStat(0);

    // Detector dimensions
    const double minX =  0.0;
    const double maxX = 47.0;
    const double minY =-20.0; 
    const double maxY = 20.0; 
    const double minZ =  3.0;
    const double maxZ = 87.0;

    int FontStyle = 132;
    double TextSize = 0.06;
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis/figs/RecoAnalysis/";

    // Load root file
    TString RootFilePath = "/exp/lariat/app/users/epelaez/files/RecoEval_histo.root"; 
    std::unique_ptr<TFile> File(TFile::Open(RootFilePath));
    TDirectory* Directory = (TDirectory*)File->Get("RecoEval");

    // Load tree and branches
    TTree* tree = (TTree*) Directory->Get<TTree>("RecoEvalTree");

    int run, subrun, event;
    tree->SetBranchAddress("run", &run); 
    tree->SetBranchAddress("subrun", &subrun); 
    tree->SetBranchAddress("event", &event);

    double truthPionInitialEnergy, truthPionInitialKEnergy, truthPionInitialMomentum;
    tree->SetBranchAddress("truthPionInitialEnergy", &truthPionInitialEnergy);
    tree->SetBranchAddress("truthPionInitialKEnergy", &truthPionInitialKEnergy);
    tree->SetBranchAddress("truthPionInitialMomentum", &truthPionInitialMomentum);

    double truthPionVertexEnergy, truthPionVertexKEnergy, truthPionVertexMomentum;
    tree->SetBranchAddress("truthPionVertexEnergy", &truthPionVertexEnergy);
    tree->SetBranchAddress("truthPionVertexKEnergy", &truthPionVertexKEnergy);
    tree->SetBranchAddress("truthPionVertexMomentum", &truthPionVertexMomentum);

    std::vector<double>* truthProtonsEnergy  = nullptr;
    std::vector<double>* truthProtonsKEnergy = nullptr;
    std::vector<double>* truthProtonsLength  = nullptr;
    tree->SetBranchAddress("truthProtonsEnergy", &truthProtonsEnergy);
    tree->SetBranchAddress("truthProtonsKEnergy", &truthProtonsKEnergy);
    tree->SetBranchAddress("truthProtonsLength", &truthProtonsLength);

    std::vector<double>* truthProtonsEndX = nullptr;
    std::vector<double>* truthProtonsEndY = nullptr;
    std::vector<double>* truthProtonsEndZ = nullptr;
    tree->SetBranchAddress("truthProtonsEndX", &truthProtonsEndX);
    tree->SetBranchAddress("truthProtonsEndY", &truthProtonsEndY);
    tree->SetBranchAddress("truthProtonsEndZ", &truthProtonsEndZ);

    // Reco data
    std::vector<int>* matchedIdentity = nullptr;
    tree->SetBranchAddress("matchedIdentity", &matchedIdentity);

    std::vector<double>* matchedKEnergy = nullptr;
    std::vector<double>* matchedLength  = nullptr;
    std::vector<double>* matchedEndX    = nullptr;
    std::vector<double>* matchedEndY    = nullptr;
    std::vector<double>* matchedEndZ    = nullptr;
    tree->SetBranchAddress("matchedKEnergy", &matchedKEnergy);
    tree->SetBranchAddress("matchedLength", &matchedLength);
    tree->SetBranchAddress("matchedEndX", &matchedEndX);
    tree->SetBranchAddress("matchedEndY", &matchedEndY);
    tree->SetBranchAddress("matchedEndZ", &matchedEndZ);

    std::vector<double>* recoMeanDEDX = nullptr;
    std::vector<double>* recoEndX     = nullptr;
    std::vector<double>* recoEndY     = nullptr;
    std::vector<double>* recoEndZ     = nullptr;
    tree->SetBranchAddress("recoMeanDEDX", &recoMeanDEDX);
    tree->SetBranchAddress("recoEndX", &recoEndX);
    tree->SetBranchAddress("recoEndY", &recoEndY);
    tree->SetBranchAddress("recoEndZ", &recoEndZ);

    std::vector<std::vector<double>>* recoDEDX = nullptr;
    std::vector<std::vector<double>>* recoResR = nullptr;
    tree->SetBranchAddress("recoDEDX", &recoDEDX);
    tree->SetBranchAddress("recoResR", &recoResR);

    // Declare variables for histogram creation
    std::vector<int> PartPDGCodes = {-211, 2212};
    std::vector<TString> PartNames = {"Pion", "Proton"};

    // Declare histograms
    TH1D *hPionMeanDEDX = new TH1D("hPionMeanDEDX", "PionMeanDEDX;;", 40, 1, 4.5);
    TH1D *hProtonMeanDEDX = new TH1D("hProtonMeanDEDX", "ProtonMeanDEDX;;", 40, 0, 10);

    TH1D *hTrueInitialKEnergyAllProtons = new TH1D("hTrueInitialKEnergyAllProtons", "hTrueInitialKEnergyAllProtons;;", 50, 0, 0.2);
    TH1D *hTrueInitialKEnergyRecoProtons = new TH1D("hTrueInitialKEnergyRecoProtons", "hTrueInitialKEnergyRecoProtons;;", 50, 0, 0.2);

    TH1D *hTrueLengthAllProtons = new TH1D("hTrueLengthAllProtons", "hTrueLengthAllProtons;;", 20, 0, 40);
    TH1D *hTrueLengthRecoProtons = new TH1D("hTrueLengthRecoProtons", "hTrueLengthRecoProtons;;", 20, 0, 40);

    TH2D *hEnergyLossAll = new TH2D(
        "hEnergyLossAll", 
        "hEnergyLossAll",
        45, 0, 90, // 45 bins starting from 0 to 90 for residual range
        45, 0, 50  // 45 bins starting from 0 to 50 for dE/dx
    );

    TH2D *hEnergyLossProtons = new TH2D(
        "hEnergyLossProtons", 
        "hEnergyLossProtons",
        45, 0, 90, // 45 bins starting from 0 to 90 for residual range
        45, 0, 50  // 45 bins starting from 0 to 50 for dE/dx
    );

    TH2D *hEnergyLossPions = new TH2D(
        "hEnergyLossPions", 
        "hEnergyLossPions",
        45, 0, 90, // 45 bins starting from 0 to 90 for residual range
        45, 0, 50  // 45 bins starting from 0 to 50 for dE/dx
    );

    int recoProtonsTrueUncontained = 0;
    int recoProtonsUncontained     = 0;
    int truthProtonsUncontained    = 0;

    // Loop over events
    Int_t NumEntries = (Int_t) tree->GetEntries();
    std::cout << "Num entries: " << NumEntries << std::endl;
    for (Int_t i = 0; i < NumEntries; ++i) {
        tree->GetEntry(i);

        // Loop over true daughter protons
        int numTruthProtons = truthProtonsKEnergy->size();
        for (int iTrueProton = 0; iTrueProton < numTruthProtons; ++iTrueProton) {
            hTrueInitialKEnergyAllProtons->Fill(truthProtonsKEnergy->at(iTrueProton));
            hTrueLengthAllProtons->Fill(truthProtonsLength->at(iTrueProton));

            if (
                (truthProtonsEndX->at(iTrueProton) < minX) || (truthProtonsEndX->at(iTrueProton) > maxX) ||
                (truthProtonsEndY->at(iTrueProton) < minY) || (truthProtonsEndY->at(iTrueProton) > maxY) ||
                (truthProtonsEndZ->at(iTrueProton) < minZ) || (truthProtonsEndZ->at(iTrueProton) > maxZ) 
            ) { truthProtonsUncontained++; }
        }

        // Loop over reco particles
        int numRecoParticles = matchedIdentity->size();
        for (int iParticle = 0; iParticle < numRecoParticles; ++iParticle) {
            // Get calo information for this particle
            std::vector<double> thisTrackDEDX     = recoDEDX->at(iParticle);
            std::vector<double> thisTrackRecoResR = recoResR->at(iParticle);

            // Pion
            if (matchedIdentity->at(iParticle) == -211) {
                hPionMeanDEDX->Fill(recoMeanDEDX->at(iParticle));

                int caloPoints = thisTrackDEDX.size();
                for (int iCalo = 0; iCalo < caloPoints; iCalo++) {
                    hEnergyLossPions->Fill(thisTrackRecoResR[iCalo], thisTrackDEDX[iCalo]);
                }
            } 
            // Proton
            else if (matchedIdentity->at(iParticle) == 2212) {
                hProtonMeanDEDX->Fill(recoMeanDEDX->at(iParticle));
                hTrueInitialKEnergyRecoProtons->Fill(matchedKEnergy->at(iParticle));
                hTrueLengthRecoProtons->Fill(matchedLength->at(iParticle));

                // std::cout << "Comparison: " << std::endl;
                // std::cout << recoEndX->at(iParticle) << "    " << matchedEndX->at(iParticle) << std::endl;
                // std::cout << recoEndY->at(iParticle) << "    " << matchedEndY->at(iParticle) << std::endl;
                // std::cout << recoEndZ->at(iParticle) << "    " << matchedEndZ->at(iParticle) << std::endl;

                if (
                    (recoEndX->at(iParticle) < minX) || (recoEndX->at(iParticle) > maxX) ||
                    (recoEndY->at(iParticle) < minY) || (recoEndY->at(iParticle) > maxY) ||
                    (recoEndZ->at(iParticle) < minZ) || (recoEndZ->at(iParticle) > maxZ) 
                ) { recoProtonsUncontained++; }

                if (
                    (matchedEndX->at(iParticle) < minX) || (matchedEndX->at(iParticle) > maxX) ||
                    (matchedEndY->at(iParticle) < minY) || (matchedEndY->at(iParticle) > maxY) ||
                    (matchedEndZ->at(iParticle) < minZ) || (matchedEndZ->at(iParticle) > maxZ) 
                ) { recoProtonsUncontained++; }

                int caloPoints = thisTrackDEDX.size();
                for (int iCalo = 0; iCalo < caloPoints; iCalo++) {
                    hEnergyLossProtons->Fill(thisTrackRecoResR[iCalo], thisTrackDEDX[iCalo]);
                }
            }

            if ((matchedIdentity->at(iParticle) == -211) || (matchedIdentity->at(iParticle) == 2212)) {
                int caloPoints = thisTrackDEDX.size();
                for (int iCalo = 0; iCalo < caloPoints; iCalo++) {
                    hEnergyLossAll->Fill(thisTrackRecoResR[iCalo], thisTrackDEDX[iCalo]);
                }
            }
        }
    }

    std::cout << "Number of uncontained truth protons: " << truthProtonsUncontained << std::endl;
    std::cout << "Number of uncontained reco protons: " << recoProtonsUncontained << std::endl;
    std::cout << "Number of truly uncontained reco protons: " << recoProtonsTrueUncontained << std::endl;

    // Setup for drawing plots
    std::vector<int> Colors = {
        kBlack, kBlue, kRed, kGreen
    };

    std::vector<std::vector<TH1*>> PlotGroups = {
        {hPionMeanDEDX, hProtonMeanDEDX},
        {hProtonMeanDEDX},
        {hTrueInitialKEnergyAllProtons, hTrueInitialKEnergyRecoProtons},
        {hTrueLengthAllProtons, hTrueLengthRecoProtons}
    };

    std::vector<std::vector<TString>> PlotLabelGroups = {
        {"Pion", "Proton"}, 
        {"Proton"},
        {"True protons", "Reco protons"},
        {"True protons", "Reco protons"}
    };

    std::vector<TString> PlotTitles = {
        "MeanDEDX",
        "ProtonMeanDEDX",
        "ProtonKE",
        "ProtonLength"
    };
    
    std::vector<TString> XLabels = {
        "dE/dx (MeV/cm)",
        "dE/dx (MeV/cm)",
        "Kinetic energy (GeV)",
        "Track length (cm)"
    };

    std::vector<TString> YLabels = {
        "Particle counts",
        "Particle counts",
        "Particle counts",
        "Particle counts"
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
        Plots[0]->GetYaxis()->SetTitle(YLabels.at(iPlot));
        Plots[0]->GetYaxis()->SetTitleOffset(1.1);
        Plots[0]->GetYaxis()->CenterTitle();	

        for (int iSubPlot = 0; iSubPlot < Plots.size(); ++iSubPlot) {
            leg->AddEntry(Plots[iSubPlot], Labels[iSubPlot], "l");
            Plots[iSubPlot]->SetLineWidth(2);
            Plots[iSubPlot]->SetLineColor(Colors.at(iSubPlot));
            Plots[iSubPlot]->Draw("hist same");
        }

        if (PlotTitles[iPlot] == "MeanDEDX") {
            TLine *v_line= new TLine(4,-10,4,180);
            v_line->SetLineColor(kRed);
            v_line->SetLineWidth(2);
            v_line->SetLineStyle(kDashed);
            v_line->Draw("same");
        }

        leg->Draw();
        PlotCanvas->SaveAs(SaveDir + PlotTitles.at(iPlot) + ".png");
        delete PlotCanvas;
    }

    // Set up efficiency plots
    std::vector<std::tuple<TH1*, TH1*>> EfficiencyHistoPairs = {
        {hTrueInitialKEnergyAllProtons, hTrueInitialKEnergyRecoProtons},
        {hTrueLengthAllProtons, hTrueLengthRecoProtons}
    };

    std::vector<TString> EfficiencyTitles = {
        "ProtonKEEfficiency",
        "ProtonLengthEfficiency"
    };

    std::vector<TString> EfficiencyYLabels = {
        "Truth-matched to true protons",
        "Truth-matched to true protons"
    };

    int numEffPlots = EfficiencyHistoPairs.size();
    for (int iPlot = 0; iPlot < numEffPlots; ++iPlot) {
        TCanvas* PlotCanvas = new TCanvas("Canvas","Canvas",205,34,1024,768);
        PlotCanvas->cd();
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetBottomMargin(0.15);

        std::tuple<TH1*, TH1*> EfficiencyPair = EfficiencyHistoPairs.at(iPlot);
        TH1* HistoPassed = std::get<1>(EfficiencyPair); 
        TH1* HistoTotal = std::get<0>(EfficiencyPair);

        // Change title of original histos so eff plot inherits it
        HistoPassed->GetYaxis()->SetTitle(EfficiencyYLabels.at(iPlot));
        HistoTotal->GetYaxis()->SetTitle(EfficiencyYLabels.at(iPlot));

        // Create efficiency object
        TEfficiency* Eff = new TEfficiency(*HistoPassed, *HistoTotal);

        Eff->SetMarkerStyle(21);
        Eff->SetMarkerColor(kBlack);
        Eff->Draw("AP");

        gPad->Update();
        // double max = Eff->GetPaintedGraph()->GetYaxis()->GetBinUpEdge(Eff->GetPaintedGraph()->GetYaxis()->GetNbins());
        // Eff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0., max*1.1);

        Eff->GetPaintedGraph()->GetXaxis()->CenterTitle();
        Eff->GetPaintedGraph()->GetXaxis()->SetTitleFont(FontStyle);
        Eff->GetPaintedGraph()->GetXaxis()->SetLabelFont(FontStyle);
        Eff->GetPaintedGraph()->GetXaxis()->SetNdivisions(8);
        Eff->GetPaintedGraph()->GetXaxis()->SetLabelSize(TextSize);
        Eff->GetPaintedGraph()->GetXaxis()->SetTitleSize(TextSize);
        Eff->GetPaintedGraph()->GetXaxis()->SetTitleOffset(1.1);

        Eff->GetPaintedGraph()->GetYaxis()->CenterTitle();
        Eff->GetPaintedGraph()->GetYaxis()->SetTitleFont(FontStyle);
        Eff->GetPaintedGraph()->GetYaxis()->SetLabelFont(FontStyle);
        Eff->GetPaintedGraph()->GetYaxis()->SetNdivisions(6);
        Eff->GetPaintedGraph()->GetYaxis()->SetLabelSize(TextSize);
        Eff->GetPaintedGraph()->GetYaxis()->SetTitleSize(TextSize);
        Eff->GetPaintedGraph()->GetYaxis()->SetTitleOffset(1.3);
        Eff->GetPaintedGraph()->GetYaxis()->SetTickSize(0);

        PlotCanvas->SaveAs(SaveDir + EfficiencyTitles.at(iPlot) + ".png");

        delete PlotCanvas;
    }

    TCanvas* c1 = new TCanvas("c1", "EnergyLossPlots", 800, 600);
    hEnergyLossAll->SetMinimum(0);
    hEnergyLossAll->SetMaximum(hEnergyLossAll->GetMaximum());
    hEnergyLossAll->Draw("COLZ");
    c1->SaveAs(SaveDir + "EnergyLossAll.png");

    hEnergyLossProtons->SetMinimum(0);
    hEnergyLossProtons->SetMaximum(hEnergyLossProtons->GetMaximum());
    hEnergyLossProtons->Draw("COLZ");
    c1->SaveAs(SaveDir + "EnergyLossProtons.png");

    hEnergyLossPions->SetMinimum(0);
    hEnergyLossPions->SetMaximum(hEnergyLossPions->GetMaximum());
    hEnergyLossPions->Draw("COLZ");
    c1->SaveAs(SaveDir + "EnergyLossPions.png");
}