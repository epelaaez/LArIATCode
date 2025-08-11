#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TEfficiency.h>

#include <vector>

#include "Helpers.cc"

void RecoAnalysis() {
    // Set defaults
    gStyle->SetOptStat(0); // get rid of stats box
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    gStyle->SetPalette(kRainBow);

    TGraph* gProton = new TGraph();
    TGraph* gPion   = new TGraph();
    TGraph* gMuonTG = new TGraph();

    // Initialize points
    initializeProtonPoints(gProton);
    initializePionPoints(gPion);
    initializeMuonNoBraggPoints(gMuonTG);

    // Vertex radius
    double fVertexRadius = 4.;

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

    double truthPionEndX, truthPionEndY, truthPionEndZ;
    tree->SetBranchAddress("truthPionEndX", &truthPionEndX);
    tree->SetBranchAddress("truthPionEndY", &truthPionEndY);
    tree->SetBranchAddress("truthPionEndZ", &truthPionEndZ);

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

    std::vector<double>* matchedKEnergy  = nullptr;
    std::vector<double>* matchedLength   = nullptr;
    std::vector<double>* matchedRealEndX = nullptr;
    std::vector<double>* matchedRealEndY = nullptr;
    std::vector<double>* matchedRealEndZ = nullptr;
    tree->SetBranchAddress("matchedKEnergy", &matchedKEnergy);
    tree->SetBranchAddress("matchedLength", &matchedLength);
    tree->SetBranchAddress("matchedRealEndX", &matchedRealEndX);
    tree->SetBranchAddress("matchedRealEndY", &matchedRealEndY);
    tree->SetBranchAddress("matchedRealEndZ", &matchedRealEndZ);

    std::vector<double>* recoMeanDEDX = nullptr;
    std::vector<double>* recoEndX     = nullptr;
    std::vector<double>* recoEndY     = nullptr;
    std::vector<double>* recoEndZ     = nullptr;
    std::vector<double>* recoBeginX   = nullptr;
    std::vector<double>* recoBeginY   = nullptr;
    std::vector<double>* recoBeginZ   = nullptr;
    std::vector<int>*    recoTrkID    = nullptr;
    tree->SetBranchAddress("recoMeanDEDX", &recoMeanDEDX);
    tree->SetBranchAddress("recoEndX", &recoEndX);
    tree->SetBranchAddress("recoEndY", &recoEndY);
    tree->SetBranchAddress("recoEndZ", &recoEndZ);
    tree->SetBranchAddress("recoBeginX", &recoBeginX);
    tree->SetBranchAddress("recoBeginY", &recoBeginY);
    tree->SetBranchAddress("recoBeginZ", &recoBeginZ);
    tree->SetBranchAddress("recoTrkID", &recoTrkID);

    int WC2TPCtrkID;
    tree->SetBranchAddress("WC2TPCtrkID", &WC2TPCtrkID);

    std::vector<std::vector<double>>* recoDEDX = nullptr;
    std::vector<std::vector<double>>* recoResR = nullptr;
    tree->SetBranchAddress("recoDEDX", &recoDEDX);
    tree->SetBranchAddress("recoResR", &recoResR);

    std::vector<bool>* isTrackInverted = nullptr;
    tree->SetBranchAddress("isTrackInverted", &isTrackInverted);

    // Declare variables for histogram creation
    std::vector<int> PartPDGCodes = {-211, 2212};
    std::vector<TString> PartNames = {"Pion", "Proton"};

    // Declare histograms
    TH1D *hPionMeanDEDXSecondary = new TH1D("hPionMeanDEDXSecondary", "PionMeanDEDX;;", 40, 0, 12);
    TH1D *hProtonMeanDEDXSecondary = new TH1D("hProtonMeanDEDXSecondary", "ProtonMeanDEDX;;", 40, 0, 12);

    TH1D *hPionMeanDEDXPrimary = new TH1D("hPionMeanDEDXPrimary", "PionMeanDEDX;;", 40, 0, 12);
    TH1D *hProtonMeanDEDXPrimary = new TH1D("hProtonMeanDEDXPrimary", "ProtonMeanDEDX;;", 40, 0, 12);

    TH1D* hUncontainedProtonMeanDEDX = new TH1D("hUncontainedProtonMeanDEDX", "hUncontainedProtonMeanDEDX;;", 40, 0, 12);
    TH1D* hContainedProtonMeanDEDX = new TH1D("hContainedProtonMeanDEDX", "hContainedProtonMeanDEDX;;", 40, 0, 12);

    TH1D *hTrueInitialKEnergyAllProtons = new TH1D("hTrueInitialKEnergyAllProtons", "hTrueInitialKEnergyAllProtons;;", 50, 0, 0.2);
    TH1D *hTrueInitialKEnergyRecoProtons = new TH1D("hTrueInitialKEnergyRecoProtons", "hTrueInitialKEnergyRecoProtons;;", 50, 0, 0.2);

    TH1D *hTrueLengthAllProtons = new TH1D("hTrueLengthAllProtons", "hTrueLengthAllProtons;;", 20, 0, 40);
    TH1D *hTrueLengthRecoProtons = new TH1D("hTrueLengthRecoProtons", "hTrueLengthRecoProtons;;", 20, 0, 40);

    TH1D *hProtonChi2All     = new TH1D("hProtonChi2All", "hProtonChi2All;;", 30, 0, 15);
    TH1D *hProtonChi2Protons = new TH1D("hProtonChi2Protons", "hProtonChi2Protons;;", 30, 0, 15);
    TH1D *hProtonChi2Pions   = new TH1D("hProtonChi2Pions", "hProtonChi2Pions;;", 30, 0, 15);
    TH1D *hProtonChi2Others  = new TH1D("hProtonChi2Others", "hProtonChi2Others;;", 30, 0, 15);

    TH1D *hPionChi2All     = new TH1D("hPionChi2All", "hPionChi2All;;", 20, 0, 10);
    TH1D *hPionChi2Protons = new TH1D("hPionChi2Protons", "hPionChi2Protons;;", 20, 0, 10);
    TH1D *hPionChi2Pions   = new TH1D("hPionChi2Pions", "hPionChi2Pions;;", 20, 0, 10);
    TH1D *hPionChi2Others  = new TH1D("hPionChi2Others", "hPionChi2Others;;", 20, 0, 10);

    int protonChiNumSteps = 20;
    double protonChiStart = 0;
    double protonChiEnd   = 10;
    double protonChiStep  = (protonChiEnd - protonChiStart) / ((double) protonChiNumSteps);

    int pionChiNumSteps = 10;
    double pionChiStart = 0;
    double pionChiEnd   = 5;
    double pionChiStep  = (pionChiEnd - pionChiStart) / ((double) pionChiNumSteps);

    TH2F* hProtonChi2TruePositives = new TH2F(
        "hProtonChi2TruePositives", 
        "hProtonChi2TruePositives;Proton chi squared;Pion chi squared", 
        protonChiNumSteps, protonChiStart, protonChiEnd, 
        pionChiNumSteps, pionChiStart, pionChiEnd
    );
    TH2F* hProtonChi2FalsePositives = new TH2F(
        "hProtonChi2FalsePositives", 
        "hProtonChi2FalsePositives;Proton chi squared;Pion chi squared", 
        protonChiNumSteps, protonChiStart, protonChiEnd, 
        pionChiNumSteps, pionChiStart, pionChiEnd
    );
    TH2F* hProtonChi2TrueNegatives = new TH2F(
        "hProtonChi2TrueNegatives", 
        "hProtonChi2TrueNegatives;Proton chi squared;Pion chi squared", 
        protonChiNumSteps, protonChiStart, protonChiEnd, 
        pionChiNumSteps, pionChiStart, pionChiEnd
    );
    TH2F* hProtonChi2FalseNegatives = new TH2F(
        "hProtonChi2FalseNegatives", 
        "hProtonChi2FalseNegatives;Proton chi squared;Pion chi squared", 
        protonChiNumSteps, protonChiStart, protonChiEnd, 
        pionChiNumSteps, pionChiStart, pionChiEnd
    );

    TH2F* hPionChi2TruePositives = new TH2F(
        "hPionChi2TruePositives", 
        "hPionChi2TruePositives;Proton chi squared;Pion chi squared", 
        protonChiNumSteps, protonChiStart, protonChiEnd, 
        pionChiNumSteps, pionChiStart, pionChiEnd
    );
    TH2F* hPionChi2FalsePositives = new TH2F(
        "hPionChi2FalsePositives", 
        "hPionChi2FalsePositives;Proton chi squared;Pion chi squared", 
        protonChiNumSteps, protonChiStart, protonChiEnd, 
        pionChiNumSteps, pionChiStart, pionChiEnd
    );
    TH2F* hPionChi2TrueNegatives = new TH2F(
        "hPionChi2TrueNegatives", 
        "hPionChi2TrueNegatives;Proton chi squared;Pion chi squared", 
        protonChiNumSteps, protonChiStart, protonChiEnd, 
        pionChiNumSteps, pionChiStart, pionChiEnd
    );
    TH2F* hPionChi2FalseNegatives = new TH2F(
        "hPionChi2FalseNegatives", 
        "hPionChi2FalseNegatives;Proton chi squared;Pion chi squared", 
        protonChiNumSteps, protonChiStart, protonChiEnd, 
        pionChiNumSteps, pionChiStart, pionChiEnd
    );

    double startingZBoundary = 82.0; // start at usual reduced volume
    double stepZBoundary     = 1;    // at each step, decrease by this
    int    numStepsZBoundary = 60;   // how many steps to take
    TH1D *hProtonsEscapingZReduction = new TH1D("hProtonsEscapingZReduction", "hProtonsEscapingZReduction", numStepsZBoundary, 0, numStepsZBoundary * stepZBoundary);

    double startingXBoundary = 42.0;
    double stepXBoundary     = 1;
    int    numStepsXBoundary = 30;
    TH1D *hProtonsEscapingXReduction = new TH1D("hProtonsEscapingXReduction", "hProtonsEscapingXReduction", numStepsXBoundary, 0, numStepsXBoundary * stepXBoundary);

    double startingYBoundary = 15.0;
    double stepYBoundary     = 1;
    int    numStepsYBoundary = 10;
    TH1D *hProtonsEscapingYReduction = new TH1D("hProtonsEscapingYReduction", "hProtonsEscapingYReduction", numStepsYBoundary, 0, numStepsYBoundary * stepYBoundary);

    double stepBoundary     = 1;
    int    numStepsBoundary = 10;
    TH1D *hProtonsEscapingAllReduction = new TH1D("hProtonsEscapingAllReduction", "hProtonsEscapingAllReduction", numStepsBoundary, 0, numStepsBoundary * stepBoundary);

    TH2D *hEnergyLossAll = new TH2D(
        "hEnergyLossAll", 
        "hEnergyLossAll;Residual range (cm);dE/dx (MeV/cm)",
        100, 0, 35, // 100 bins starting from 0 to 90 for residual range
        100, 0, 40  // 100 bins starting from 0 to 25 for dE/dx
    );

    TH2D *hEnergyLossProtons = new TH2D(
        "hEnergyLossProtons", 
        "hEnergyLossProtons;Residual range (cm);dE/dx (MeV/cm)",
        100, 0, 35, // 100 bins starting from 0 to 90 for residual range
        100, 0, 40  // 100 bins starting from 0 to 25 for dE/dx
    );

    TH2D *hEnergyLossUncontainedProtons = new TH2D(
        "hEnergyLossUncontainedProtons", 
        "hEnergyLossUncontainedProtons;Residual range (cm);dE/dx (MeV/cm)",
        100, 0, 35, // 100 bins starting from 0 to 90 for residual range
        100, 0, 40  // 100 bins starting from 0 to 25 for dE/dx
    );

    TH2D *hEnergyLossContainedProtons = new TH2D(
        "hEnergyLossContainedProtons", 
        "hEnergyLossContainedProtons;Residual range (cm);dE/dx (MeV/cm)",
        100, 0, 35, // 100 bins starting from 0 to 90 for residual range
        100, 0, 40  // 100 bins starting from 0 to 25 for dE/dx
    );

    TH2D *hEnergyLossPions = new TH2D(
        "hEnergyLossPions", 
        "hEnergyLossPions;Residual range (cm);dE/dx (MeV/cm)",
        100, 0, 35, // 100 bins starting from 0 to 90 for residual range
        100, 0, 40  // 100 bins starting from 0 to 25 for dE/dx
    );

    double startDEDXCut    = 3.0;
    double stepSizeDEDXCut = 0.125;
    int    numStepsDEDXCut = 60;
    TH1D *hMisIDsMeanDEDXCut       = new TH1D("hMisIDsMeanDEDXCut", "hMisIDsMeanDEDXCut", numStepsDEDXCut, startDEDXCut, startDEDXCut + numStepsDEDXCut * stepSizeDEDXCut);
    TH1D *hMisProtonIDsMeanDEDXCut = new TH1D("hMisProtonIDsMeanDEDXCut", "hMisProtonIDsMeanDEDXCut", numStepsDEDXCut, startDEDXCut, startDEDXCut + numStepsDEDXCut * stepSizeDEDXCut);
    TH1D *hMisPionIDsMeanDEDXCut   = new TH1D("hMisPionIDsMeanDEDXCut", "hMisPionIDsMeanDEDXCut", numStepsDEDXCut, startDEDXCut, startDEDXCut + numStepsDEDXCut * stepSizeDEDXCut);

    int totalRecoTracks = 0;

    int recoProtonsTrueUncontained = 0;
    int recoProtonsUncontained     = 0;
    int truthProtonsUncontained    = 0;
    int recoProtons                = 0;
    int truthProtons               = 0;
    int recoPions                  = 0;

    int reversedTracks                 = 0;
    int reversedTracksMatchedToProtons = 0;
    int reversedTracksMatchedToPions   = 0;

    std::ofstream outUncontainedProtonsFile("files/UncontainedProtons.txt");

    // Mean dE/dx cut performance
    int totalSecondaryProtons = 0;
    int totalSecondaryPions   = 0;
    int totalSecondaryMuons   = 0;
    int totalSecondaryElectrons = 0;
    int totalSecondaryOthers  = 0;
    std::vector<int> idAsProtons(numStepsDEDXCut);
    std::vector<int> idAsPions(numStepsDEDXCut);
    std::vector<int> correctIDProtons(numStepsDEDXCut);
    std::vector<int> correctIDPions(numStepsDEDXCut);

    // Loop over events
    Int_t NumEntries = (Int_t) tree->GetEntries();
    std::cout << "Num entries: " << NumEntries << std::endl;
    for (Int_t i = 0; i < NumEntries; ++i) {
        tree->GetEntry(i);

        // Get vertex position
        double primaryVertexX, primaryVertexY, primaryVertexZ;
        for (int iParticle = 0; iParticle < matchedIdentity->size(); iParticle++) {
            if (recoTrkID->at(iParticle) == WC2TPCtrkID) {
                primaryVertexX = recoEndX->at(iParticle);
                primaryVertexY = recoEndY->at(iParticle);
                primaryVertexZ = recoEndZ->at(iParticle);
                break;
            }
        }

        // Look at energy profile for events we are interested in
        if ((event == 200061) || (event == 200920) || (event == 62282)) {
            int numRecoParticles = matchedIdentity->size();
            for (int iParticle = 0; iParticle < numRecoParticles; ++iParticle) {
                std::vector<double> thisTrackDEDX     = recoDEDX->at(iParticle);
                std::vector<double> thisTrackRecoResR = recoResR->at(iParticle);
                int caloPoints = thisTrackDEDX.size();

                TH2D *hDEDXProfile = new TH2D(
                    "hDEDXProfile", 
                    "hDEDXProfile;Residual range (cm);dE/dx (MeV/cm)",
                    caloPoints, 0, 0, 
                    caloPoints, 0, 0  
                );

                for (int iCalo = 0; iCalo < caloPoints; iCalo++) {
                    hDEDXProfile->Fill(thisTrackRecoResR[iCalo], thisTrackDEDX[iCalo]);
                }

                TCanvas* c1 = new TCanvas("c1", "DEDXProfiles", 800, 600);
                hDEDXProfile->SetMinimum(0);
                hDEDXProfile->SetMaximum(1);
                hDEDXProfile->Draw("COLZ");
                Overlay_dEdx_RR_Reference_PP(gProton, gPion, gMuonTG);
                c1->SaveAs(SaveDir +  "dEdxProfiles/" + event + "_" + iParticle + "_" + matchedIdentity->at(iParticle) + ".png");
                
                delete c1;
                delete hDEDXProfile;
            }
        }

        // Loop over true daughter protons
        int numTruthProtons = truthProtonsKEnergy->size();
        for (int iTrueProton = 0; iTrueProton < numTruthProtons; ++iTrueProton) {
            hTrueInitialKEnergyAllProtons->Fill(truthProtonsKEnergy->at(iTrueProton));
            hTrueLengthAllProtons->Fill(truthProtonsLength->at(iTrueProton));

            double thisProtonTrueKE = truthProtonsKEnergy->at(iTrueProton);
            if ((thisProtonTrueKE >= PROTON_ENERGY_LOWER_BOUND) && (thisProtonTrueKE <= PROTON_ENERGY_UPPER_BOUND)) { truthProtons++; }
            else { continue; }

            if (
                (truthProtonsEndX->at(iTrueProton) < minX) || (truthProtonsEndX->at(iTrueProton) > maxX) ||
                (truthProtonsEndY->at(iTrueProton) < minY) || (truthProtonsEndY->at(iTrueProton) > maxY) ||
                (truthProtonsEndZ->at(iTrueProton) < minZ) || (truthProtonsEndZ->at(iTrueProton) > maxZ) 
            ) truthProtonsUncontained++;
        }

        // Loop over reco particles
        int numRecoParticles = matchedIdentity->size();
        int thisRecoProtonsTrueUncotained = 0;
        for (int iParticle = 0; iParticle < numRecoParticles; ++iParticle) {
            totalRecoTracks++;

            // Get calo information for this particle
            std::vector<double> thisTrackDEDX     = recoDEDX->at(iParticle);
            std::vector<double> thisTrackRecoResR = recoResR->at(iParticle);
            double thisTrackTrueKE                = matchedKEnergy->at(iParticle);
            double thisTrackMeanDEDX              = recoMeanDEDX->at(iParticle);

            int thisTrackID      = recoTrkID->at(iParticle);
            bool isPrimaryTrack  = (thisTrackID == WC2TPCtrkID);
            double startDistance = sqrt(
                pow(recoBeginX->at(iParticle) - primaryVertexX, 2) + pow(recoBeginY->at(iParticle) - primaryVertexY, 2) + pow(recoBeginZ->at(iParticle) - primaryVertexZ, 2)
            );
            double endDistance   = sqrt(
                pow(recoEndX->at(iParticle) - primaryVertexX, 2) + pow(recoEndY->at(iParticle) - primaryVertexY, 2) + pow(recoEndZ->at(iParticle) - primaryVertexZ, 2)
            );
            bool isSecondaryTrack = false;
            if ((thisTrackID != WC2TPCtrkID) && ((startDistance < fVertexRadius) || (endDistance < fVertexRadius))) isSecondaryTrack = true;

            if (isTrackInverted->at(iParticle)) reversedTracks++;

            // Get chi2 fit if secondary
            if (isSecondaryTrack) {
                // Get chi^2 fit to a pion and proton
                int caloPoints    = thisTrackDEDX.size();
                double protonChi2 = computeReducedChi2(gProton, thisTrackRecoResR, thisTrackDEDX, false, caloPoints);
                double pionChi2   = computeReducedChi2(gPion, thisTrackRecoResR, thisTrackDEDX, false, caloPoints);

                hProtonChi2All->Fill(protonChi2);
                hPionChi2All->Fill(pionChi2);

                if (matchedIdentity->at(iParticle) == -211) {
                    hProtonChi2Pions->Fill(protonChi2);
                    hPionChi2Pions->Fill(pionChi2);
                } else if (matchedIdentity->at(iParticle) == 2212) {
                    hProtonChi2Protons->Fill(protonChi2);
                    hPionChi2Protons->Fill(pionChi2);
                } else {
                    hProtonChi2Others->Fill(protonChi2);
                    hPionChi2Others->Fill(pionChi2);
                }

                // Fill chi2 fit plots
                for (int iPionChiStep = 0; iPionChiStep < pionChiNumSteps; iPionChiStep++) {
                    double currentPionChiValue = pionChiStart + (iPionChiStep * pionChiStep);
                    for (int iProtonChiStep = 0; iProtonChiStep < protonChiNumSteps; iProtonChiStep++) {
                        double currentProtonChiValue = protonChiStart + (iProtonChiStep * protonChiStep);
                        
                        if ((pionChi2 < currentPionChiValue) && (protonChi2 > currentProtonChiValue)) {
                            // Tagged as pion
                            if (matchedIdentity->at(iParticle) == -211) {
                                hPionChi2TruePositives->Fill(currentProtonChiValue, currentPionChiValue);
                            } else {
                                hPionChi2FalsePositives->Fill(currentProtonChiValue, currentPionChiValue);
                            }
                        } else if ((pionChi2 > currentPionChiValue) && (protonChi2 < currentProtonChiValue)) {
                            // Tagged as proton
                            if (matchedIdentity->at(iParticle) == 2212) {
                                hProtonChi2TruePositives->Fill(currentProtonChiValue, currentPionChiValue);
                            } else {
                                hProtonChi2FalsePositives->Fill(currentProtonChiValue, currentPionChiValue);
                            }
                        }

                        if (!((pionChi2 < currentPionChiValue) && (protonChi2 > currentProtonChiValue))) {
                            if (matchedIdentity->at(iParticle) == -211) {
                                hPionChi2FalseNegatives->Fill(currentProtonChiValue, currentPionChiValue);
                            } else {
                                hPionChi2TrueNegatives->Fill(currentProtonChiValue, currentPionChiValue);
                            }
                        }

                        if (!((pionChi2 > currentPionChiValue) && (protonChi2 < currentProtonChiValue))) {
                            if (matchedIdentity->at(iParticle) == 2212) {
                                hProtonChi2FalseNegatives->Fill(currentProtonChiValue, currentPionChiValue);
                            } else {
                                hProtonChi2TrueNegatives->Fill(currentProtonChiValue, currentPionChiValue);
                            }
                        }
                    }
                }
            }

            // Pion
            if (matchedIdentity->at(iParticle) == -211) {
                recoPions++;
                if (isTrackInverted->at(iParticle)) reversedTracksMatchedToPions++;
                
                if (isSecondaryTrack) {
                    hPionMeanDEDXSecondary->Fill(thisTrackMeanDEDX);
                } else if (isPrimaryTrack) {
                    hPionMeanDEDXPrimary->Fill(thisTrackMeanDEDX);
                }

                int caloPoints = thisTrackDEDX.size();
                for (int iCalo = 0; iCalo < caloPoints; iCalo++) {
                    hEnergyLossPions->Fill(thisTrackRecoResR[iCalo], thisTrackDEDX[iCalo]);
                }
            } 
            // Proton
            else if (matchedIdentity->at(iParticle) == 2212) {
                if (isTrackInverted->at(iParticle)) reversedTracksMatchedToProtons++;

                hTrueInitialKEnergyRecoProtons->Fill(matchedKEnergy->at(iParticle));
                hTrueLengthRecoProtons->Fill(matchedLength->at(iParticle));

                if ((thisTrackTrueKE >= PROTON_ENERGY_LOWER_BOUND) && (thisTrackTrueKE <= PROTON_ENERGY_UPPER_BOUND)) { recoProtons++; }
                else { continue; }

                if (isSecondaryTrack) {
                    hProtonMeanDEDXSecondary->Fill(thisTrackMeanDEDX);
                } else if (isPrimaryTrack) {
                    hProtonMeanDEDXPrimary->Fill(thisTrackMeanDEDX);
                }

                // if ((thisTrackMeanDEDX > 1.98) && (thisTrackMeanDEDX < 2.02)) {
                //     std::cout << event << std::endl;
                //     for (int iParticleA = 0; iParticleA < numRecoParticles; ++iParticleA) {
                //         std::cout << matchedIdentity->at(iParticleA) << "   " << recoMeanDEDX->at(iParticleA) << std::endl;
                //     }
                //     std::cout << std::endl;
                // }

                if (
                    (recoEndX->at(iParticle) < minX) || (recoEndX->at(iParticle) > maxX) ||
                    (recoEndY->at(iParticle) < minY) || (recoEndY->at(iParticle) > maxY) ||
                    (recoEndZ->at(iParticle) < minZ) || (recoEndZ->at(iParticle) > maxZ) 
                ) recoProtonsUncontained++; 

                if (
                    (matchedRealEndX->at(iParticle) < minX) || (matchedRealEndX->at(iParticle) > maxX) ||
                    (matchedRealEndY->at(iParticle) < minY) || (matchedRealEndY->at(iParticle) > maxY) ||
                    (matchedRealEndZ->at(iParticle) < minZ) || (matchedRealEndZ->at(iParticle) > maxZ) 
                ) {
                    recoProtonsTrueUncontained++; 
                    thisRecoProtonsTrueUncotained++;
                    hUncontainedProtonMeanDEDX->Fill(thisTrackMeanDEDX);
                    
                    int caloPoints = thisTrackDEDX.size();
                    for (int iCalo = 0; iCalo < caloPoints; iCalo++) {
                        hEnergyLossUncontainedProtons->Fill(thisTrackRecoResR[iCalo], thisTrackDEDX[iCalo]);
                    }
                } else {
                    hContainedProtonMeanDEDX->Fill(thisTrackMeanDEDX);

                    int caloPoints = thisTrackDEDX.size();
                    for (int iCalo = 0; iCalo < caloPoints; iCalo++) {
                        hEnergyLossContainedProtons->Fill(thisTrackRecoResR[iCalo], thisTrackDEDX[iCalo]);
                    }
                }

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

            // Find performace in cut
            if (isSecondaryTrack) {
                if (matchedIdentity->at(iParticle) == -211) totalSecondaryPions++;
                else if (matchedIdentity->at(iParticle) == 2212) totalSecondaryProtons++;
                else if ((matchedIdentity->at(iParticle) == 13) || (matchedIdentity->at(iParticle) == -13)) totalSecondaryMuons++;
                else if ((matchedIdentity->at(iParticle) == 11) || (matchedIdentity->at(iParticle) == -11)) totalSecondaryElectrons++;
                else totalSecondaryOthers++;
                
                for (int stepCut = 0; stepCut < numStepsDEDXCut; stepCut++) {
                    double cutValue = startDEDXCut + (double)(stepCut * stepSizeDEDXCut);

                    if (thisTrackMeanDEDX > cutValue) {
                        idAsProtons[stepCut]++;
                        if (matchedIdentity->at(iParticle) == -211) {
                            hMisIDsMeanDEDXCut->Fill(cutValue);
                            hMisPionIDsMeanDEDXCut->Fill(cutValue);
                        } else if (matchedIdentity->at(iParticle) == 2212) {
                            correctIDProtons[stepCut]++;
                        }
                    } else if (thisTrackMeanDEDX <= cutValue) {
                        idAsPions[stepCut]++;
                        if (matchedIdentity->at(iParticle) == -211) {
                            correctIDPions[stepCut]++;
                        } else if (matchedIdentity->at(iParticle) == 2212) {
                            hMisIDsMeanDEDXCut->Fill(cutValue);
                            hMisProtonIDsMeanDEDXCut->Fill(cutValue);
                        }
                    }
                }
            }
        }

        if (thisRecoProtonsTrueUncotained > 0) {
            outUncontainedProtonsFile << "Event number: " << event << ". Uncontained protons: " << thisRecoProtonsTrueUncotained << std::endl;
        }

        // Find bin for protons escaping plot
        int binZ = (startingZBoundary - truthPionEndZ) / stepZBoundary;
        if (binZ > numStepsZBoundary) binZ = numStepsZBoundary;
        if (binZ < 1)                 binZ = 0;
        for (int b = 0; b <= binZ; ++b) {
            hProtonsEscapingZReduction->Fill(b * stepZBoundary, thisRecoProtonsTrueUncotained);
        }

        int binX = (startingXBoundary - truthPionEndX) / stepXBoundary;
        if (binX > numStepsXBoundary) binX = numStepsXBoundary;
        if (binX < 1)                 binX = 0;
        for (int b = 0; b <= binX; ++b) {
            hProtonsEscapingXReduction->Fill(b * stepXBoundary, thisRecoProtonsTrueUncotained);
        }

        int binY = (startingYBoundary - truthPionEndY) / stepYBoundary;
        if (binY > numStepsYBoundary) binY = numStepsYBoundary;
        if (binY < 1)                 binY = 0;
        for (int b = 0; b <= binY; ++b) {
            hProtonsEscapingYReduction->Fill(b * stepYBoundary, thisRecoProtonsTrueUncotained);
        }

        int overallBinX = (startingXBoundary - truthPionEndX) / stepBoundary;
        int overallBinY = (startingYBoundary - truthPionEndY) / stepBoundary;
        int overallBinZ = (startingZBoundary - truthPionEndZ) / stepBoundary;
        int overallBin  = std::min({overallBinX, overallBinY, overallBinZ});
        for (int b = 0; b <= overallBin; ++b) {
            hProtonsEscapingAllReduction->Fill(b * stepBoundary, thisRecoProtonsTrueUncotained);
        }
    }

    TH1D *hEffProtonMeanDEDXCut  = new TH1D("hEffProtonMeanDEDXCut", "hEffProtonMeanDEDXCut", numStepsDEDXCut, startDEDXCut, startDEDXCut + numStepsDEDXCut * stepSizeDEDXCut);
    TH1D *hEffPionMeanDEDXCut    = new TH1D("hEffPionMeanDEDXCut", "hEffPionMeanDEDXCut", numStepsDEDXCut, startDEDXCut, startDEDXCut + numStepsDEDXCut * stepSizeDEDXCut);

    TH1D *hPurProtonMeanDEDXCut  = new TH1D("hPurProtonMeanDEDXCut", "hPurProtonMeanDEDXCut", numStepsDEDXCut, startDEDXCut, startDEDXCut + numStepsDEDXCut * stepSizeDEDXCut);
    TH1D *hPurPionMeanDEDXCut    = new TH1D("hPurPionMeanDEDXCut", "hPurPionMeanDEDXCut", numStepsDEDXCut, startDEDXCut, startDEDXCut + numStepsDEDXCut * stepSizeDEDXCut);

    TH1D *hEffPurProtonMeanDEDXCut = new TH1D("hEffPurProtonMeanDEDXCut", "hEffPurProtonMeanDEDXCut", numStepsDEDXCut, startDEDXCut, startDEDXCut + numStepsDEDXCut * stepSizeDEDXCut);
    TH1D *hEffPurPionMeanDEDXCut   = new TH1D("hEffPurPionMeanDEDXCut", "hEffPurPionMeanDEDXCut", numStepsDEDXCut, startDEDXCut, startDEDXCut + numStepsDEDXCut * stepSizeDEDXCut);

    std::cout << std::endl;
    std::cout << "Total secondary pions: " << totalSecondaryPions << std::endl;
    std::cout << "Total secondary protons: " << totalSecondaryProtons << std::endl;
    std::cout << "Total secondary muons: " << totalSecondaryMuons << std::endl;
    std::cout << "Total secondary electrons: " << totalSecondaryElectrons << std::endl;
    std::cout << "Total secondary others: " << totalSecondaryOthers << std::endl;
    std::cout << std::endl;

    for (int i = 0; i < numStepsDEDXCut; i++) {
        double cutValue = startDEDXCut + (double)(i * stepSizeDEDXCut);

        double efficiencyProton = correctIDProtons[i] / (double) totalSecondaryProtons;
        double purityProton     = correctIDProtons[i] / (double) idAsProtons[i];
        hEffProtonMeanDEDXCut->Fill(cutValue, efficiencyProton);
        hPurProtonMeanDEDXCut->Fill(cutValue, purityProton);
        hEffPurProtonMeanDEDXCut->Fill(cutValue, efficiencyProton * purityProton);

        double efficiencyPion = correctIDPions[i] / (double) totalSecondaryPions;
        double purityPion     = correctIDPions[i] / (double) idAsPions[i];
        // std::cout << cutValue << std::endl;
        // std::cout << "   pion purity: " << purityPion << std::endl;
        // std::cout << "   correct id as pions: " << correctIDPions[i] << std::endl;
        // std::cout << "   id as pions: " << idAsPions[i] << std::endl;
        // std::cout << std::endl;
        // std::cout << "   proton purity: " << purityProton << std::endl;
        // std::cout << "   correct id as protons: " << correctIDProtons[i] << std::endl;
        // std::cout << "   id as protons: " << idAsProtons[i] << std::endl;
        hEffPionMeanDEDXCut->Fill(cutValue, efficiencyPion);
        hPurPionMeanDEDXCut->Fill(cutValue, purityPion);
        hEffPurPionMeanDEDXCut->Fill(cutValue, efficiencyPion * purityPion);
    }

    hEnergyLossUncontainedProtons->Scale(1. / hEnergyLossUncontainedProtons->Integral());
    hEnergyLossContainedProtons->Scale(1. / hEnergyLossContainedProtons->Integral());
    hUncontainedProtonMeanDEDX->Scale(1. / hUncontainedProtonMeanDEDX->Integral());
    hContainedProtonMeanDEDX->Scale(1. / hContainedProtonMeanDEDX->Integral());

    int    minBin     = hMisIDsMeanDEDXCut->GetMinimumBin();
    double bestCutVal = hMisIDsMeanDEDXCut->GetBinLowEdge(minBin);
    double bestMisIDs = hMisIDsMeanDEDXCut->GetBinContent(minBin);

    std::cout << std::endl;
    std::cout << "Total reco tracks: " << totalRecoTracks << std::endl; 
    std::cout << "Minimum overall mis-IDs at " << bestCutVal << "MeV/cm, with " << 100 * (bestMisIDs / totalRecoTracks) << "%" << std::endl;

    // hPionMeanDEDXSecondary->Scale(1 / hPionMeanDEDXSecondary->Integral());
    // hProtonMeanDEDXSecondary->Scale(1 / hProtonMeanDEDXSecondary->Integral());

    std::cout << std::endl;
    std::cout << "Number of truth protons: " << truthProtons << std::endl;
    std::cout << "Number of uncontained truth protons: " << truthProtonsUncontained << std::endl;
    std::cout << std::endl;
    std::cout << "Number of reconstructed protons: " << recoProtons << std::endl;
    std::cout << "Number of uncontained reco protons: " << recoProtonsUncontained << std::endl;
    std::cout << "Number of truly uncontained reco protons: " << recoProtonsTrueUncontained << std::endl;
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout << "Number of reversed tracks: " << reversedTracks << std::endl;
    std::cout << "Number of reversed tracks matched to pions: " << reversedTracksMatchedToPions << std::endl;
    std::cout << "Number of reversed tracks matched to protons: " << reversedTracksMatchedToProtons << std::endl;
    std::cout << std::endl;

    // Setup for drawing plots
    std::vector<int> Colors = {
        kBlack, kBlue, kRed, kGreen
    };

    std::vector<std::vector<TH1*>> PlotGroups = {
        {hProtonMeanDEDXSecondary},
        {hPionMeanDEDXSecondary, hProtonMeanDEDXSecondary},
        {hPionMeanDEDXPrimary, hProtonMeanDEDXPrimary},
        {hTrueInitialKEnergyAllProtons, hTrueInitialKEnergyRecoProtons},
        {hTrueLengthAllProtons, hTrueLengthRecoProtons},
        {hProtonsEscapingZReduction},
        {hProtonsEscapingYReduction},
        {hProtonsEscapingXReduction},
        {hProtonsEscapingAllReduction},
        {hContainedProtonMeanDEDX, hUncontainedProtonMeanDEDX},
        {hMisIDsMeanDEDXCut, hMisPionIDsMeanDEDXCut, hMisProtonIDsMeanDEDXCut},
        {hEffPurPionMeanDEDXCut, hEffPurProtonMeanDEDXCut},
        {hEffPionMeanDEDXCut, hEffProtonMeanDEDXCut},
        {hPurPionMeanDEDXCut, hPurProtonMeanDEDXCut},
        {hProtonChi2Protons, hProtonChi2Pions},
        {hPionChi2Protons, hPionChi2Pions}
    };

    std::vector<std::vector<TString>> PlotLabelGroups = {
        {"Proton"},
        {"Pion", "Proton"}, 
        {"Pion", "Proton"}, 
        {"True protons", "Reco protons"},
        {"True protons", "Reco protons"},
        {"Protons escaping"},
        {"Protons escaping"},
        {"Protons escaping"},
        {"Protons escaping"},
        {"Contained protons", "Uncontained protons"},
        {"Total", "Pion", "Proton"},
        {"Pion", "Proton"},
        {"Pion", "Proton"},
        {"Pion", "Proton"},
        {"Protons", "Pions"},
        {"Protons", "Pions"}
    };

    std::vector<TString> PlotTitles = {
        "ProtonMeanDEDXSecondary",
        "MeanDEDXSecondary",
        "MeanDEDXPrimary",
        "ProtonKE",
        "ProtonLength",
        "ProtonsEscapingZReduction",
        "ProtonsEscapingYReduction",
        "ProtonsEscapingXReduction",
        "ProtonsEscapingAllReduction",
        "ProtonMeanDEDXUncontained",
        "MisIDsDEDXCut",
        "EfficiencyPurityDEDXCut",
        "EfficiencyDEDXCut",
        "PurityDEDXCut",
        "ProtonChi2Fit",
        "PionChi2Fit"
    };
    
    std::vector<TString> XLabels = {
        "Mean dE/dx (MeV/cm)",
        "Mean dE/dx (MeV/cm)",
        "Mean dE/dx (MeV/cm)",
        "Kinetic energy (GeV)",
        "Track length (cm)",
        "Reduced volume pushback (only z direction) (cm)",
        "Reduced volume pushback (only y direction) (cm)",
        "Reduced volume pushback (only x direction) (cm)",
        "Reduced volume pushback (all directions) (cm)",
        "Mean dE/dx (MeV/cm)",
        "Mean dE/dx cut value (MeV/cm)",
        "Mean dE/dx cut value (MeV/cm)",
        "Mean dE/dx cut value (MeV/cm)",
        "Mean dE/dx cut value (MeV/cm)",
        "Reduced chi squared",
        "Reduced chi squared"
    };

    std::vector<TString> YLabels = {
        "Particle counts",
        "Particle counts",
        "Particle counts",
        "Particle counts",
        "Particle counts",
        "Protons escaping",
        "Protons escaping",
        "Protons escaping",
        "Protons escaping",
        "Proton counts",
        "Misidentified particle counts",
        "Efficiency times purity",
        "Efficiency",
        "Purity",
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

            double imax = TMath::Max(Plots[iSubPlot]->GetMaximum(), Plots[0]->GetMaximum());

            double YAxisRange = 1.15*imax;
            Plots[iSubPlot]->GetYaxis()->SetRangeUser(0., YAxisRange);
            Plots[0]->GetYaxis()->SetRangeUser(0., YAxisRange);	
        }

        // if (PlotTitles[iPlot] == "MeanDEDX") {
        //     TLine *v_line= new TLine(4,-10,4,180);
        //     v_line->SetLineColor(kRed);
        //     v_line->SetLineWidth(2);
        //     v_line->SetLineStyle(kDashed);
        //     v_line->Draw("same");
        // }

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

    gStyle->SetOptStat(0);

    TCanvas* c1 = new TCanvas("c1", "EnergyLossPlots", 800, 600);
    hEnergyLossAll->SetMinimum(0);
    hEnergyLossAll->SetMaximum(hEnergyLossAll->GetMaximum());
    hEnergyLossAll->Draw("COLZ");
    Overlay_dEdx_RR_Reference_PP(gProton, gPion, gMuonTG);
    c1->SaveAs(SaveDir + "EnergyLossAll.png");

    hEnergyLossProtons->SetMinimum(0);
    hEnergyLossProtons->SetMaximum(hEnergyLossProtons->GetMaximum());
    hEnergyLossProtons->Draw("COLZ");
    Overlay_dEdx_RR_Reference_PP(gProton, gPion, gMuonTG);
    c1->SaveAs(SaveDir + "EnergyLossProtons.png");

    hEnergyLossPions->SetMinimum(0);
    hEnergyLossPions->SetMaximum(hEnergyLossPions->GetMaximum());
    hEnergyLossPions->Draw("COLZ");
    Overlay_dEdx_RR_Reference_PP(gProton, gPion, gMuonTG);
    c1->SaveAs(SaveDir + "EnergyLossPions.png");

    hEnergyLossUncontainedProtons->SetMinimum(0);
    hEnergyLossUncontainedProtons->SetMaximum(hEnergyLossUncontainedProtons->GetMaximum());
    hEnergyLossUncontainedProtons->Draw("COLZ");
    Overlay_dEdx_RR_Reference_PP(gProton, gPion, gMuonTG);
    c1->SaveAs(SaveDir + "EnergyLossUncontainedProtons.png");

    hEnergyLossContainedProtons->SetMinimum(0);
    hEnergyLossContainedProtons->SetMaximum(hEnergyLossContainedProtons->GetMaximum());
    hEnergyLossContainedProtons->Draw("COLZ");
    Overlay_dEdx_RR_Reference_PP(gProton, gPion, gMuonTG);
    c1->SaveAs(SaveDir + "EnergyLossContainedProtons.png");

    hProtonChi2TruePositives->SetMinimum(0);
    hProtonChi2TruePositives->SetMaximum(hProtonChi2TruePositives->GetMaximum());
    hProtonChi2TruePositives->Draw("COLZ");
    c1->SaveAs(SaveDir + "ProtonChi2TruePositives.png");

    hProtonChi2FalsePositives->SetMinimum(0);
    hProtonChi2FalsePositives->SetMaximum(hProtonChi2FalsePositives->GetMaximum());
    hProtonChi2FalsePositives->Draw("COLZ");
    c1->SaveAs(SaveDir + "ProtonChi2FalsePositives.png");

    hProtonChi2TrueNegatives->SetMinimum(0);
    hProtonChi2TrueNegatives->SetMaximum(hProtonChi2TrueNegatives->GetMaximum());
    hProtonChi2TrueNegatives->Draw("COLZ");
    c1->SaveAs(SaveDir + "ProtonChi2TrueNegatives.png");

    hProtonChi2FalseNegatives->SetMinimum(0);
    hProtonChi2FalseNegatives->SetMaximum(hProtonChi2FalseNegatives->GetMaximum());
    hProtonChi2FalseNegatives->Draw("COLZ");
    c1->SaveAs(SaveDir + "ProtonChi2FalseNegatives.png");

    hPionChi2TruePositives->SetMinimum(0);
    hPionChi2TruePositives->SetMaximum(hPionChi2TruePositives->GetMaximum());
    hPionChi2TruePositives->Draw("COLZ");
    c1->SaveAs(SaveDir + "PionChi2TruePositives.png");

    hPionChi2FalsePositives->SetMinimum(0);
    hPionChi2FalsePositives->SetMaximum(hPionChi2FalsePositives->GetMaximum());
    hPionChi2FalsePositives->Draw("COLZ");
    c1->SaveAs(SaveDir + "PionChi2FalsePositives.png");

    hPionChi2TrueNegatives->SetMinimum(0);
    hPionChi2TrueNegatives->SetMaximum(hPionChi2TrueNegatives->GetMaximum());
    hPionChi2TrueNegatives->Draw("COLZ");
    c1->SaveAs(SaveDir + "PionChi2TrueNegatives.png");

    hPionChi2FalseNegatives->SetMinimum(0);
    hPionChi2FalseNegatives->SetMaximum(hPionChi2FalseNegatives->GetMaximum());
    hPionChi2FalseNegatives->Draw("COLZ");
    c1->SaveAs(SaveDir + "PionChi2FalseNegatives.png");

    // TH2F* hPionChi2Trues = dynamic_cast<TH2F*>(hPionChi2TruePositives->Clone("hPionChi2Trues"));
    // hPionChi2Trues->Add(hPionChi2TrueNegatives);
    // hPionChi2Trues->SetTitle("hPionChi2Trues");
    // hPionChi2Trues->SetMinimum(0);
    // hPionChi2Trues->SetMaximum(hPionChi2Trues->GetMaximum());
    // hPionChi2Trues->Draw("COLZ");
    // c1->SaveAs(SaveDir + "PionChi2Trues.png");

    // Compute figure of merit for chi^2 cut
    TH2F* hProtonChi2FOM = new TH2F(
        "hProtonChi2FOM", 
        "hProtonChi2FOM;Proton chi squared;Pion chi squared", 
        protonChiNumSteps, protonChiStart, protonChiEnd, 
        pionChiNumSteps, pionChiStart, pionChiEnd
    );
    TH2F* hPionChi2FOM = new TH2F(
        "hPionChi2FOM", 
        "hPionChi2FOM;Proton chi squared;Pion chi squared", 
        protonChiNumSteps, protonChiStart, protonChiEnd, 
        pionChiNumSteps, pionChiStart, pionChiEnd
    );
    
    for (int iPionChiStep = 0; iPionChiStep < pionChiNumSteps; iPionChiStep++) {
        double currentPionChiValue = pionChiStart + (iPionChiStep * pionChiStep);
        for (int iProtonChiStep = 0; iProtonChiStep < protonChiNumSteps; iProtonChiStep++) {
            double currentProtonChiValue = protonChiStart + (iProtonChiStep * protonChiStep);

            int binX = hPionChi2TruePositives->GetXaxis()->FindFixBin(currentProtonChiValue);
            int binY = hPionChi2TruePositives->GetYaxis()->FindFixBin(currentPionChiValue);

            double pion_tp     = hPionChi2TruePositives->GetBinContent(binX, binY);
            double pion_tn     = hPionChi2TrueNegatives->GetBinContent(binX, binY);
            double pion_fp     = hPionChi2FalsePositives->GetBinContent(binX, binY);
            double pion_fn     = hPionChi2FalseNegatives->GetBinContent(binX, binY);
            double pion_eff    = pion_tp / (pion_tp + pion_tn);
            double pion_purity = pion_tp / (pion_tp + pion_fp);

            // double pion_fom = pion_tp / std::sqrt(pion_tp + pion_fp);
            double pion_fom = (2 * pion_tp) / (2 * pion_tp + pion_fp + pion_fn);
            // double pion_fom = ((pion_tp * pion_tn) - (pion_fp * pion_fn)) / std::sqrt((pion_tp + pion_fp) * (pion_tp + pion_fn) * (pion_tn * pion_fp) * (pion_tn * pion_fn));
            hPionChi2FOM->SetBinContent(binX, binY, pion_fom);

            double proton_tp     = hProtonChi2TruePositives->GetBinContent(binX, binY);
            double proton_tn     = hProtonChi2TrueNegatives->GetBinContent(binX, binY);
            double proton_fp     = hProtonChi2FalsePositives->GetBinContent(binX, binY);
            double proton_fn     = hProtonChi2FalseNegatives->GetBinContent(binX, binY);
            double proton_eff    = proton_tp / (proton_tp + proton_tn);
            double proton_purity = proton_tp / (proton_tp + proton_fp);

            // double proton_fom = proton_tp / std::sqrt(proton_tp + proton_fp);
            double proton_fom = (2 * proton_tp) / (2 * proton_tp + proton_fp + proton_fn);
            // double proton_fom = ((proton_tp * proton_tn) - (proton_fp * proton_fn)) / std::sqrt((proton_tp + proton_fp) * (proton_tp + proton_fn) * (proton_tn * proton_fp) * (proton_tn * proton_fn));
            hProtonChi2FOM->SetBinContent(binX, binY, proton_fom);
        }
    }

    hPionChi2FOM->SetMinimum(0);
    hPionChi2FOM->SetMaximum(hPionChi2FOM->GetMaximum());
    hPionChi2FOM->Draw("COLZ");
    c1->SaveAs(SaveDir + "PionChi2FOM.png");

    hProtonChi2FOM->SetMinimum(0);
    hProtonChi2FOM->SetMaximum(hProtonChi2FOM->GetMaximum());
    hProtonChi2FOM->Draw("COLZ");
    c1->SaveAs(SaveDir + "ProtonChi2FOM.png");

    Int_t binxPion, binyPion, binzPion;
    Int_t globalBinPion   = hPionChi2FOM->GetMaximumBin(binxPion, binyPion, binzPion);

    double maxContentPion = hPionChi2FOM->GetBinContent(globalBinPion);
    double xAtMaxPion     = hPionChi2FOM->GetXaxis()->GetBinCenter(binxPion);
    double yAtMaxPion     = hPionChi2FOM->GetYaxis()->GetBinCenter(binyPion);

    std::cout << std::endl;
    std::cout << "[PION]  Max‑FOM = " << maxContentPion
            << " at (χ²_p, χ²_π) = (" << xAtMaxPion << ", " << yAtMaxPion << ")\n";

    Int_t binxProt, binyProt, binzProt;
    Int_t globalBinProt   = hProtonChi2FOM->GetMaximumBin(binxProt, binyProt, binzProt);

    double maxContentProt = hProtonChi2FOM->GetBinContent(globalBinProt);
    double xAtMaxProt     = hProtonChi2FOM->GetXaxis()->GetBinCenter(binxProt);
    double yAtMaxProt     = hProtonChi2FOM->GetYaxis()->GetBinCenter(binyProt);

    std::cout << "[PROTON] Max‑FOM = " << maxContentProt
            << " at (χ²_p, χ²_π) = (" << xAtMaxProt << ", " << yAtMaxProt << ")\n";
    std::cout << std::endl;
}
