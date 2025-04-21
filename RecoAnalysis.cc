#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TEfficiency.h>

#include <vector>

void Overlay_dEdx_RR_Reference_PP(bool addLegend = true, TVirtualPad* pad = gPad);

void RecoAnalysis() {
    // Set defaults
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    // gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);

    // Proton energy bounds
    double PROTON_ENERGY_LOWER_BOUND = 0.075;
    double PROTON_ENERGY_UPPER_BOUND = 1.0;

    // Detector dimensions
    const double minX =  0.0;
    const double maxX = 47.0;
    const double minY =-20.0; 
    const double maxY = 20.0; 
    const double minZ =  3.0;
    const double maxZ = 87.0;

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
                Overlay_dEdx_RR_Reference_PP();
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
        {hPurPionMeanDEDXCut, hPurProtonMeanDEDXCut}
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
        {"Pion", "Proton"}
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
        "PurityDEDXCut"
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
        "Mean dE/dx cut value (MeV/cm)"
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
        "Purity"
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

    TCanvas* c1 = new TCanvas("c1", "EnergyLossPlots", 800, 600);
    hEnergyLossAll->SetMinimum(0);
    hEnergyLossAll->SetMaximum(hEnergyLossAll->GetMaximum());
    hEnergyLossAll->Draw("COLZ");
    Overlay_dEdx_RR_Reference_PP();
    c1->SaveAs(SaveDir + "EnergyLossAll.png");

    hEnergyLossProtons->SetMinimum(0);
    hEnergyLossProtons->SetMaximum(hEnergyLossProtons->GetMaximum());
    hEnergyLossProtons->Draw("COLZ");
    Overlay_dEdx_RR_Reference_PP();
    c1->SaveAs(SaveDir + "EnergyLossProtons.png");

    hEnergyLossPions->SetMinimum(0);
    hEnergyLossPions->SetMaximum(hEnergyLossPions->GetMaximum());
    hEnergyLossPions->Draw("COLZ");
    Overlay_dEdx_RR_Reference_PP();
    c1->SaveAs(SaveDir + "EnergyLossPions.png");

    hEnergyLossUncontainedProtons->SetMinimum(0);
    hEnergyLossUncontainedProtons->SetMaximum(hEnergyLossUncontainedProtons->GetMaximum());
    hEnergyLossUncontainedProtons->Draw("COLZ");
    Overlay_dEdx_RR_Reference_PP();
    c1->SaveAs(SaveDir + "EnergyLossUncontainedProtons.png");

    hEnergyLossContainedProtons->SetMinimum(0);
    hEnergyLossContainedProtons->SetMaximum(hEnergyLossContainedProtons->GetMaximum());
    hEnergyLossContainedProtons->Draw("COLZ");
    Overlay_dEdx_RR_Reference_PP();
    c1->SaveAs(SaveDir + "EnergyLossContainedProtons.png");
}

void Overlay_dEdx_RR_Reference_PP(bool addLegend = true, TVirtualPad* pad = gPad) {
    const int N = 107;

    static TGraph *gProton = nullptr;
    static TGraph *gPion   = nullptr;

    if (!gProton) {
        /* ---------- Proton ---------- */
        gProton = new TGraph(N);
        gProton->SetLineColor(kRed+1);
        gProton->SetLineWidth(3);
        gProton->SetName("Proton");

        gProton->SetPoint(0,   31.95, 4.14);
        gProton->SetPoint(1,   31.65, 4.16);
        gProton->SetPoint(2,   31.35, 4.17);
        gProton->SetPoint(3,   31.05, 4.18);
        gProton->SetPoint(4,   30.75, 4.20);
        gProton->SetPoint(5,   30.45, 4.21);
        gProton->SetPoint(6,   30.15, 4.23);
        gProton->SetPoint(7,   29.85, 4.25);
        gProton->SetPoint(8,   29.55, 4.26);
        gProton->SetPoint(9,   29.25, 4.28);
        gProton->SetPoint(10,  28.95, 4.29);
        gProton->SetPoint(11,  28.65, 4.31);
        gProton->SetPoint(12,  28.35, 4.33);
        gProton->SetPoint(13,  28.05, 4.34);
        gProton->SetPoint(14,  27.75, 4.36);
        gProton->SetPoint(15,  27.45, 4.38);
        gProton->SetPoint(16,  27.15, 4.40);
        gProton->SetPoint(17,  26.85, 4.42);
        gProton->SetPoint(18,  26.55, 4.43);
        gProton->SetPoint(19,  26.25, 4.45);
        gProton->SetPoint(20,  25.95, 4.47);
        gProton->SetPoint(21,  25.65, 4.49);
        gProton->SetPoint(22,  25.35, 4.51);
        gProton->SetPoint(23,  25.05, 4.53);
        gProton->SetPoint(24,  24.75, 4.55);
        gProton->SetPoint(25,  24.45, 4.57);
        gProton->SetPoint(26,  24.15, 4.60);
        gProton->SetPoint(27,  23.85, 4.62);
        gProton->SetPoint(28,  23.55, 4.64);
        gProton->SetPoint(29,  23.25, 4.66);
        gProton->SetPoint(30,  22.95, 4.69);
        gProton->SetPoint(31,  22.65, 4.71);
        gProton->SetPoint(32,  22.35, 4.73);
        gProton->SetPoint(33,  22.05, 4.76);
        gProton->SetPoint(34,  21.75, 4.78);
        gProton->SetPoint(35,  21.45, 4.81);
        gProton->SetPoint(36,  21.15, 4.83);
        gProton->SetPoint(37,  20.85, 4.86);
        gProton->SetPoint(38,  20.55, 4.89);
        gProton->SetPoint(39,  20.25, 4.92);
        gProton->SetPoint(40,  19.95, 4.94);
        gProton->SetPoint(41,  19.65, 4.97);
        gProton->SetPoint(42,  19.35, 5.00);
        gProton->SetPoint(43,  19.05, 5.03);
        gProton->SetPoint(44,  18.75, 5.07);
        gProton->SetPoint(45,  18.45, 5.10);
        gProton->SetPoint(46,  18.15, 5.13);
        gProton->SetPoint(47,  17.85, 5.16);
        gProton->SetPoint(48,  17.55, 5.20);
        gProton->SetPoint(49,  17.25, 5.23);
        gProton->SetPoint(50,  16.95, 5.27);
        gProton->SetPoint(51,  16.65, 5.31);
        gProton->SetPoint(52,  16.35, 5.35);
        gProton->SetPoint(53,  16.05, 5.39);
        gProton->SetPoint(54,  15.75, 5.43);
        gProton->SetPoint(55,  15.45, 5.47);
        gProton->SetPoint(56,  15.15, 5.51);
        gProton->SetPoint(57,  14.85, 5.56);
        gProton->SetPoint(58,  14.55, 5.60);
        gProton->SetPoint(59,  14.25, 5.65);
        gProton->SetPoint(60,  13.95, 5.70);
        gProton->SetPoint(61,  13.65, 5.75);
        gProton->SetPoint(62,  13.35, 5.80);
        gProton->SetPoint(63,  13.05, 5.85);
        gProton->SetPoint(64,  12.75, 5.91);
        gProton->SetPoint(65,  12.45, 5.97);
        gProton->SetPoint(66,  12.15, 6.03);
        gProton->SetPoint(67,  11.85, 6.09);
        gProton->SetPoint(68,  11.55, 6.15);
        gProton->SetPoint(69,  11.25, 6.22);
        gProton->SetPoint(70,  10.95, 6.29);
        gProton->SetPoint(71,  10.65, 6.36);
        gProton->SetPoint(72,  10.35, 6.44);
        gProton->SetPoint(73,  10.05, 6.52);
        gProton->SetPoint(74,   9.75, 6.60);
        gProton->SetPoint(75,   9.45, 6.68);
        gProton->SetPoint(76,   9.15, 6.77);
        gProton->SetPoint(77,   8.85, 6.87);
        gProton->SetPoint(78,   8.55, 6.97);
        gProton->SetPoint(79,   8.25, 7.08);
        gProton->SetPoint(80,   7.95, 7.19);
        gProton->SetPoint(81,   7.65, 7.30);
        gProton->SetPoint(82,   7.35, 7.43);
        gProton->SetPoint(83,   7.05, 7.56);
        gProton->SetPoint(84,   6.75, 7.70);
        gProton->SetPoint(85,   6.45, 7.85);
        gProton->SetPoint(86,   6.15, 8.02);
        gProton->SetPoint(87,   5.85, 8.19);
        gProton->SetPoint(88,   5.55, 8.38);
        gProton->SetPoint(89,   5.25, 8.58);
        gProton->SetPoint(90,   4.95, 8.81);
        gProton->SetPoint(91,   4.65, 9.05);
        gProton->SetPoint(92,   4.35, 9.32);
        gProton->SetPoint(93,   4.05, 9.61);
        gProton->SetPoint(94,   3.75, 9.94);
        gProton->SetPoint(95,   3.45,10.32);
        gProton->SetPoint(96,   3.15,10.74);
        gProton->SetPoint(97,   2.85,11.23);
        gProton->SetPoint(98,   2.55,11.80);
        gProton->SetPoint(99,   2.25,12.48);
        gProton->SetPoint(100,  1.95,13.31);
        gProton->SetPoint(101,  1.65,14.35);
        gProton->SetPoint(102,  1.35,15.71);
        gProton->SetPoint(103,  1.05,17.59);
        gProton->SetPoint(104,  0.75,20.44);
        gProton->SetPoint(105,  0.45,25.48);
        gProton->SetPoint(106,  0.15,38.12);

        /* ---------- Pion ---------- */
        gPion = new TGraph(N);
        gPion->SetLineColor(kRed+3);
        gPion->SetLineWidth(3);
        gPion->SetName("Pion");

        gPion->SetPoint(0,   31.95, 2.4);
        gPion->SetPoint(1,   31.65, 2.4);
        gPion->SetPoint(2,   31.35, 2.4);
        gPion->SetPoint(3,   31.05, 2.4);
        gPion->SetPoint(4,   30.75, 2.4);
        gPion->SetPoint(5,   30.45, 2.4);
        gPion->SetPoint(6,   30.15, 2.4);
        gPion->SetPoint(7,   29.85, 2.4);
        gPion->SetPoint(8,   29.55, 2.4);
        gPion->SetPoint(9,   29.25, 2.4);
        gPion->SetPoint(10,  28.95, 2.4);
        gPion->SetPoint(11,  28.65, 2.4);
        gPion->SetPoint(12,  28.35, 2.4);
        gPion->SetPoint(13,  28.05, 2.4);
        gPion->SetPoint(14,  27.75, 2.5);
        gPion->SetPoint(15,  27.45, 2.5);
        gPion->SetPoint(16,  27.15, 2.5);
        gPion->SetPoint(17,  26.85, 2.5);
        gPion->SetPoint(18,  26.55, 2.5);
        gPion->SetPoint(19,  26.25, 2.5);
        gPion->SetPoint(20,  25.95, 2.5);
        gPion->SetPoint(21,  25.65, 2.5);
        gPion->SetPoint(22,  25.35, 2.5);
        gPion->SetPoint(23,  25.05, 2.5);
        gPion->SetPoint(24,  24.75, 2.5);
        gPion->SetPoint(25,  24.45, 2.5);
        gPion->SetPoint(26,  24.15, 2.5);
        gPion->SetPoint(27,  23.85, 2.5);
        gPion->SetPoint(28,  23.55, 2.5);
        gPion->SetPoint(29,  23.25, 2.6);
        gPion->SetPoint(30,  22.95, 2.6);
        gPion->SetPoint(31,  22.65, 2.6);
        gPion->SetPoint(32,  22.35, 2.6);
        gPion->SetPoint(33,  22.05, 2.6);
        gPion->SetPoint(34,  21.75, 2.6);
        gPion->SetPoint(35,  21.45, 2.6);
        gPion->SetPoint(36,  21.15, 2.6);
        gPion->SetPoint(37,  20.85, 2.6);
        gPion->SetPoint(38,  20.55, 2.6);
        gPion->SetPoint(39,  20.25, 2.6);
        gPion->SetPoint(40,  19.95, 2.6);
        gPion->SetPoint(41,  19.65, 2.7);
        gPion->SetPoint(42,  19.35, 2.7);
        gPion->SetPoint(43,  19.05, 2.7);
        gPion->SetPoint(44,  18.75, 2.7);
        gPion->SetPoint(45,  18.45, 2.7);
        gPion->SetPoint(46,  18.15, 2.7);
        gPion->SetPoint(47,  17.85, 2.7);
        gPion->SetPoint(48,  17.55, 2.7);
        gPion->SetPoint(49,  17.25, 2.8);
        gPion->SetPoint(50,  16.95, 2.8);
        gPion->SetPoint(51,  16.65, 2.8);
        gPion->SetPoint(52,  16.35, 2.8);
        gPion->SetPoint(53,  16.05, 2.8);
        gPion->SetPoint(54,  15.75, 2.8);
        gPion->SetPoint(55,  15.45, 2.8);
        gPion->SetPoint(56,  15.15, 2.9);
        gPion->SetPoint(57,  14.85, 2.9);
        gPion->SetPoint(58,  14.55, 2.9);
        gPion->SetPoint(59,  14.25, 2.9);
        gPion->SetPoint(60,  13.95, 2.9);
        gPion->SetPoint(61,  13.65, 2.9);
        gPion->SetPoint(62,  13.35, 3.0);
        gPion->SetPoint(63,  13.05, 3.0);
        gPion->SetPoint(64,  12.75, 3.0);
        gPion->SetPoint(65,  12.45, 3.0);
        gPion->SetPoint(66,  12.15, 3.0);
        gPion->SetPoint(67,  11.85, 3.1);
        gPion->SetPoint(68,  11.55, 3.1);
        gPion->SetPoint(69,  11.25, 3.1);
        gPion->SetPoint(70,  10.95, 3.1);
        gPion->SetPoint(71,  10.65, 3.2);
        gPion->SetPoint(72,  10.35, 3.2);
        gPion->SetPoint(73,  10.05, 3.2);
        gPion->SetPoint(74,   9.75, 3.3);
        gPion->SetPoint(75,   9.45, 3.3);
        gPion->SetPoint(76,   9.15, 3.3);
        gPion->SetPoint(77,   8.85, 3.4);
        gPion->SetPoint(78,   8.55, 3.4);
        gPion->SetPoint(79,   8.25, 3.4);
        gPion->SetPoint(80,   7.95, 3.5);
        gPion->SetPoint(81,   7.65, 3.5);
        gPion->SetPoint(82,   7.35, 3.6);
        gPion->SetPoint(83,   7.05, 3.6);
        gPion->SetPoint(84,   6.75, 3.7);
        gPion->SetPoint(85,   6.45, 3.7);
        gPion->SetPoint(86,   6.15, 3.8);
        gPion->SetPoint(87,   5.85, 3.9);
        gPion->SetPoint(88,   5.55, 3.9);
        gPion->SetPoint(89,   5.25, 4.0);
        gPion->SetPoint(90,   4.95, 4.1);
        gPion->SetPoint(91,   4.65, 4.2);
        gPion->SetPoint(92,   4.35, 4.3);
        gPion->SetPoint(93,   4.05, 4.4);
        gPion->SetPoint(94,   3.75, 4.6);
        gPion->SetPoint(95,   3.45, 4.7);
        gPion->SetPoint(96,   3.15, 4.9);
        gPion->SetPoint(97,   2.85, 5.1);
        gPion->SetPoint(98,   2.55, 5.3);
        gPion->SetPoint(99,   2.25, 5.6);
        gPion->SetPoint(100,  1.95, 5.9);
        gPion->SetPoint(101,  1.65, 6.4);
        gPion->SetPoint(102,  1.35, 6.9);
        gPion->SetPoint(103,  1.05, 7.7);
        gPion->SetPoint(104,  0.75, 8.9);
        gPion->SetPoint(105,  0.45,11.0);
        gPion->SetPoint(106,  0.15,16.5);
    }

    if (pad) pad->cd();
    gProton->Draw("L same");
    gPion  ->Draw("L same");
    
    if (addLegend) {
        static TLegend *leg = nullptr;
        if (!leg) {
            leg = new TLegend(0.50,0.80,0.75,0.88);
            leg->AddEntry(gProton, "Proton", "l");
            leg->AddEntry(gPion,   "Pion",   "l");
            leg->SetFillStyle(1001);
            leg->SetBorderSize(0);
        }
        leg->Draw();
    }
}