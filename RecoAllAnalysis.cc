#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TEfficiency.h>

#include <vector>

void initializeProtonPoints(TGraph* gProton);
void initializePionPoints(TGraph* gPion);
void initializeMuonNoBraggPoints(TGraph* gMuonTG);
double computeReducedChi2(const TGraph* theory, std::vector<double> xData, std::vector<double> yData, int nPoints);

void RecoAllAnalysis() {
    // Set defaults
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
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis/figs/RecoAllAnalysis/";

    // Load root file
    TString RootFilePath = "/exp/lariat/app/users/epelaez/files/RecoAllEval_histo.root"; 
    std::unique_ptr<TFile> File(TFile::Open(RootFilePath));
    TDirectory* Directory = (TDirectory*)File->Get("RecoAllEval");

    ///////////////////
    // Load branches //
    ///////////////////

    // Load tree and branches
    TTree* tree = (TTree*) Directory->Get<TTree>("RecoEvalTree");

    int run, subrun, event;
    tree->SetBranchAddress("run", &run); 
    tree->SetBranchAddress("subrun", &subrun); 
    tree->SetBranchAddress("event", &event);

    // Signal information
    bool isPionAbsorptionSignal; int backgroundType, numVisibleProtons; 
    tree->SetBranchAddress("isPionAbsorptionSignal", &isPionAbsorptionSignal);
    tree->SetBranchAddress("backgroundType", &backgroundType);
    tree->SetBranchAddress("numVisibleProtons", &numVisibleProtons);

    // Cut information
    bool passesPionInRedVolume, passesNoOutgoingPion, passesSmallTracksCut, passesMeanCurvatureCut;
    tree->SetBranchAddress("passesPionInRedVolume", &passesPionInRedVolume);
    tree->SetBranchAddress("passesNoOutgoingPion", &passesNoOutgoingPion);
    tree->SetBranchAddress("passesSmallTracksCut", &passesSmallTracksCut);
    tree->SetBranchAddress("passesMeanCurvatureCut", &passesMeanCurvatureCut);

    // WC match information
    int WC2TPCtrkID;
    double WC2TPCPrimaryEndX, WC2TPCPrimaryEndY, WC2TPCPrimaryEndZ;
    std::vector<double>* wcMatchResR = nullptr;
    std::vector<double>* wcMatchDEDX = nullptr;
    tree->SetBranchAddress("WC2TPCtrkID", &WC2TPCtrkID);
    tree->SetBranchAddress("WC2TPCPrimaryEndX", &WC2TPCPrimaryEndX);
    tree->SetBranchAddress("WC2TPCPrimaryEndY", &WC2TPCPrimaryEndY);
    tree->SetBranchAddress("WC2TPCPrimaryEndZ", &WC2TPCPrimaryEndZ);
    tree->SetBranchAddress("wcMatchResR", &wcMatchResR);
    tree->SetBranchAddress("wcMatchDEDX", &wcMatchDEDX);

    // Reco information
    std::vector<double>* recoMeanDEDX      = nullptr;
    std::vector<double>* recoEndX          = nullptr;
    std::vector<double>* recoEndY          = nullptr;
    std::vector<double>* recoEndZ          = nullptr;
    std::vector<double>* recoBeginX        = nullptr;
    std::vector<double>* recoBeginY        = nullptr;
    std::vector<double>* recoBeginZ        = nullptr;
    std::vector<int>*    recoTrkID         = nullptr;
    std::vector<bool>*   isTrackNearVertex = nullptr;
    tree->SetBranchAddress("recoMeanDEDX", &recoMeanDEDX);
    tree->SetBranchAddress("recoEndX", &recoEndX);
    tree->SetBranchAddress("recoEndY", &recoEndY);
    tree->SetBranchAddress("recoEndZ", &recoEndZ);
    tree->SetBranchAddress("recoBeginX", &recoBeginX);
    tree->SetBranchAddress("recoBeginY", &recoBeginY);
    tree->SetBranchAddress("recoBeginZ", &recoBeginZ);
    tree->SetBranchAddress("recoTrkID", &recoTrkID);
    tree->SetBranchAddress("isTrackNearVertex", &isTrackNearVertex);

    // Calorimetry information
    std::vector<std::vector<double>>* recoResR = nullptr;
    std::vector<std::vector<double>>* recoDEDX = nullptr;
    tree->SetBranchAddress("recoResR", &recoResR);
    tree->SetBranchAddress("recoDEDX", &recoDEDX);

    // Out files
    std::ofstream outStitchFile("files/WCMatchStitching.txt");

    //////////////////////
    // Loop over events //
    //////////////////////

    Int_t NumEntries = (Int_t) tree->GetEntries();
    std::cout << "Num entries: " << NumEntries << std::endl;

    for (Int_t i = 0; i < NumEntries; ++i) {
        tree->GetEntry(i);

        // If no track matched to wire-chamber, skip
        if (WC2TPCtrkID == -99999) continue;

        // Perform chi^2 stitching for primary track
        // For this analysis, we get the following chi^2 values:
        //   - Pion (with Bragg peak) chi^2
        //     - Capture at rest
        //   - Pion (through-going) chi^2
        //     - 0p or Np
        //   - Proton chi^2
        //     - Interaction in detector front face
        //   - Scanning fit with rhs to through-going and lhs to proton
        //     - Pion + proton stitched as one track; repeat analysis with new vertex

        int totalCaloPoints     = wcMatchDEDX->size();

        double pionChi2         = computeReducedChi2(gPion, *wcMatchResR, *wcMatchDEDX, totalCaloPoints);
        double throughGoingChi2 = computeReducedChi2(gMuonTG, *wcMatchResR, *wcMatchDEDX, totalCaloPoints);
        double protonChi2       = computeReducedChi2(gProton, *wcMatchResR, *wcMatchDEDX, totalCaloPoints);

        double minStitchedChi2 = std::numeric_limits<double>::max();
        int bestBreakPoint = -1;
        for (int caloBreakPoint = 1; caloBreakPoint < totalCaloPoints - 1; ++caloBreakPoint) {
            std::vector<double> leftResR(wcMatchResR->begin(), wcMatchResR->begin() + caloBreakPoint);
            std::vector<double> leftDEDX(wcMatchDEDX->begin(), wcMatchDEDX->begin() + caloBreakPoint);

            std::vector<double> rightResR(wcMatchResR->begin() + caloBreakPoint, wcMatchResR->end());
            std::vector<double> rightDEDX(wcMatchDEDX->begin() + caloBreakPoint, wcMatchDEDX->end());

            double chi2LHS = computeReducedChi2(gProton, leftResR, leftDEDX, leftResR.size());
            double chi2RHS = computeReducedChi2(gMuonTG, rightResR, rightDEDX, rightResR.size());

            double totalChi2 = (chi2LHS * leftResR.size() + chi2RHS * rightResR.size()) / totalCaloPoints;
            
            if (totalChi2 < minStitchedChi2) {
                minStitchedChi2 = totalChi2;
                bestBreakPoint = caloBreakPoint;
            }
        }
        
        if (std::min({minStitchedChi2, pionChi2, throughGoingChi2, protonChi2}) != throughGoingChi2) {
            outStitchFile << "Event number: " << event << std::endl; 
            outStitchFile << "  Pion with Bragg peak chi^2: " << pionChi2 << std::endl; 
            outStitchFile << "  Through-going pion chi^2: " << throughGoingChi2 << std::endl; 
            outStitchFile << "  Proton with Bragg peak chi^2: " << protonChi2 << std::endl; 
            outStitchFile << "  Best break point: " << bestBreakPoint << " with stitched chi^2: " << minStitchedChi2 << std::endl;
            outStitchFile << std::endl;
        }
    }
}

double computeReducedChi2(const TGraph* theory, std::vector<double> xData, std::vector<double> yData, int nPoints) {
    double chi2 = 0.0;

    for (int i = 0; i < nPoints; ++i) {
        double theoryY = theory->Eval(xData[i]); // interpolate the theory at xData[i]
        double deltaY = yData[i] - theoryY;
        chi2 += (deltaY * deltaY) / theoryY;
    }

    // Currently, no fixed parameters, so dof = nPoints
    int dof = nPoints;
    return dof > 0 ? chi2 / dof : 0.0; 
}

void initializeProtonPoints(TGraph* gProton) {
    double protonData[107][2] = {
        {31.95, 4.14}, {31.65, 4.16}, {31.35, 4.17}, {31.05, 4.18}, {30.75, 4.20},
        {30.45, 4.21}, {30.15, 4.23}, {29.85, 4.25}, {29.55, 4.26}, {29.25, 4.28},
        {28.95, 4.29}, {28.65, 4.31}, {28.35, 4.33}, {28.05, 4.34}, {27.75, 4.36},
        {27.45, 4.38}, {27.15, 4.40}, {26.85, 4.42}, {26.55, 4.43}, {26.25, 4.45},
        {25.95, 4.47}, {25.65, 4.49}, {25.35, 4.51}, {25.05, 4.53}, {24.75, 4.55},
        {24.45, 4.57}, {24.15, 4.60}, {23.85, 4.62}, {23.55, 4.64}, {23.25, 4.66},
        {22.95, 4.69}, {22.65, 4.71}, {22.35, 4.73}, {22.05, 4.76}, {21.75, 4.78},
        {21.45, 4.81}, {21.15, 4.83}, {20.85, 4.86}, {20.55, 4.89}, {20.25, 4.92},
        {19.95, 4.94}, {19.65, 4.97}, {19.35, 5.00}, {19.05, 5.03}, {18.75, 5.07},
        {18.45, 5.10}, {18.15, 5.13}, {17.85, 5.16}, {17.55, 5.20}, {17.25, 5.23},
        {16.95, 5.27}, {16.65, 5.31}, {16.35, 5.35}, {16.05, 5.39}, {15.75, 5.43},
        {15.45, 5.47}, {15.15, 5.51}, {14.85, 5.56}, {14.55, 5.60}, {14.25, 5.65},
        {13.95, 5.70}, {13.65, 5.75}, {13.35, 5.80}, {13.05, 5.85}, {12.75, 5.91},
        {12.45, 5.97}, {12.15, 6.03}, {11.85, 6.09}, {11.55, 6.15}, {11.25, 6.22},
        {10.95, 6.29}, {10.65, 6.36}, {10.35, 6.44}, {10.05, 6.52}, {9.75, 6.60},
        {9.45, 6.68}, {9.15, 6.77}, {8.85, 6.87}, {8.55, 6.97}, {8.25, 7.08},
        {7.95, 7.19}, {7.65, 7.30}, {7.35, 7.43}, {7.05, 7.56}, {6.75, 7.70},
        {6.45, 7.85}, {6.15, 8.02}, {5.85, 8.19}, {5.55, 8.38}, {5.25, 8.58},
        {4.95, 8.81}, {4.65, 9.05}, {4.35, 9.32}, {4.05, 9.61}, {3.75, 9.94},
        {3.45, 10.32}, {3.15, 10.74}, {2.85, 11.23}, {2.55, 11.80}, {2.25, 12.48},
        {1.95, 13.31}, {1.65, 14.35}, {1.35, 15.71}, {1.05, 17.59}, {0.75, 20.44},
        {0.45, 25.48}, {0.15, 38.12}
    };

    for (int i = 0; i < 107; ++i) {
        gProton->SetPoint(i, protonData[i][0], protonData[i][1]);
    }
}

void initializePionPoints(TGraph* gPion) {
    double pionData[107][2] = {
        {31.95, 2.4}, {31.65, 2.4}, {31.35, 2.4}, {31.05, 2.4}, {30.75, 2.4},
        {30.45, 2.4}, {30.15, 2.4}, {29.85, 2.4}, {29.55, 2.4}, {29.25, 2.4},
        {28.95, 2.4}, {28.65, 2.4}, {28.35, 2.4}, {28.05, 2.4}, {27.75, 2.5},
        {27.45, 2.5}, {27.15, 2.5}, {26.85, 2.5}, {26.55, 2.5}, {26.25, 2.5},
        {25.95, 2.5}, {25.65, 2.5}, {25.35, 2.5}, {25.05, 2.5}, {24.75, 2.5},
        {24.45, 2.5}, {24.15, 2.5}, {23.85, 2.5}, {23.55, 2.5}, {23.25, 2.6},
        {22.95, 2.6}, {22.65, 2.6}, {22.35, 2.6}, {22.05, 2.6}, {21.75, 2.6},
        {21.45, 2.6}, {21.15, 2.6}, {20.85, 2.6}, {20.55, 2.6}, {20.25, 2.6},
        {19.95, 2.6}, {19.65, 2.7}, {19.35, 2.7}, {19.05, 2.7}, {18.75, 2.7},
        {18.45, 2.7}, {18.15, 2.7}, {17.85, 2.7}, {17.55, 2.7}, {17.25, 2.8},
        {16.95, 2.8}, {16.65, 2.8}, {16.35, 2.8}, {16.05, 2.8}, {15.75, 2.8},
        {15.45, 2.8}, {15.15, 2.9}, {14.85, 2.9}, {14.55, 2.9}, {14.25, 2.9},
        {13.95, 2.9}, {13.65, 2.9}, {13.35, 3.0}, {13.05, 3.0}, {12.75, 3.0},
        {12.45, 3.0}, {12.15, 3.0}, {11.85, 3.1}, {11.55, 3.1}, {11.25, 3.1},
        {10.95, 3.1}, {10.65, 3.2}, {10.35, 3.2}, {10.05, 3.2}, {9.75, 3.3},
        {9.45, 3.3}, {9.15, 3.3}, {8.85, 3.4}, {8.55, 3.4}, {8.25, 3.4},
        {7.95, 3.5}, {7.65, 3.5}, {7.35, 3.6}, {7.05, 3.6}, {6.75, 3.7},
        {6.45, 3.7}, {6.15, 3.8}, {5.85, 3.9}, {5.55, 3.9}, {5.25, 4.0},
        {4.95, 4.1}, {4.65, 4.2}, {4.35, 4.3}, {4.05, 4.4}, {3.75, 4.6},
        {3.45, 4.7}, {3.15, 4.9}, {2.85, 5.1}, {2.55, 5.3}, {2.25, 5.6},
        {1.95, 5.9}, {1.65, 6.4}, {1.35, 6.9}, {1.05, 7.7}, {0.75, 8.9},
        {0.45, 11.0}, {0.15, 16.5}
    };

    for (int i = 0; i < 107; ++i) {
        gPion->SetPoint(i, pionData[i][0], pionData[i][1]);
    }
}

void initializeMuonNoBraggPoints(TGraph* gMuonTG) {
    double muonTGData[107][2] = {
        {31.95, 2.3}, {31.65, 2.3}, {31.35, 2.3}, {31.05, 2.3}, {30.75, 2.3},
        {30.45, 2.3}, {30.15, 2.3}, {29.85, 2.3}, {29.55, 2.3}, {29.25, 2.3},
        {28.95, 2.3}, {28.65, 2.3}, {28.35, 2.3}, {28.05, 2.3}, {27.75, 2.3},
        {27.45, 2.3}, {27.15, 2.3}, {26.85, 2.3}, {26.55, 2.3}, {26.25, 2.3},
        {25.95, 2.3}, {25.65, 2.3}, {25.35, 2.3}, {25.05, 2.3}, {24.75, 2.3},
        {24.45, 2.3}, {24.15, 2.3}, {23.85, 2.3}, {23.55, 2.3}, {23.25, 2.3},
        {22.95, 2.3}, {22.65, 2.3}, {22.35, 2.3}, {22.05, 2.3}, {21.75, 2.3},
        {21.45, 2.3}, {21.15, 2.3}, {20.85, 2.3}, {20.55, 2.3}, {20.25, 2.3},
        {19.95, 2.3}, {19.65, 2.3}, {19.35, 2.3}, {19.05, 2.3}, {18.75, 2.3},
        {18.45, 2.3}, {18.15, 2.3}, {17.85, 2.3}, {17.55, 2.3}, {17.25, 2.3},
        {16.95, 2.3}, {16.65, 2.3}, {16.35, 2.3}, {16.05, 2.3}, {15.75, 2.3},
        {15.45, 2.3}, {15.15, 2.3}, {14.85, 2.3}, {14.55, 2.3}, {14.25, 2.3},
        {13.95, 2.3}, {13.65, 2.3}, {13.35, 2.3}, {13.05, 2.3}, {12.75, 2.3},
        {12.45, 2.3}, {12.15, 2.3}, {11.85, 2.3}, {11.55, 2.3}, {11.25, 2.3},
        {10.95, 2.3}, {10.65, 2.3}, {10.35, 2.3}, {10.05, 2.3}, {9.75, 2.3},
        {9.45, 2.3}, {9.15, 2.3}, {8.85, 2.3}, {8.55, 2.3}, {8.25, 2.3},
        {7.95, 2.3}, {7.65, 2.3}, {7.35, 2.3}, {7.05, 2.3}, {6.75, 2.3},
        {6.45, 2.3}, {6.15, 2.3}, {5.85, 2.3}, {5.55, 2.3}, {5.25, 2.3},
        {4.95, 2.3}, {4.65, 2.3}, {4.35, 2.3}, {4.05, 2.3}, {3.75, 2.3},
        {3.45, 2.3}, {3.15, 2.3}, {2.85, 2.3}, {2.55, 2.3}, {2.25, 2.3},
        {1.95, 2.3}, {1.65, 2.3}, {1.35, 2.3}, {1.05, 2.3}, {0.75, 2.3},
        {0.45, 2.3}, {0.15, 2.3}
    };

    for (int i = 0; i < 107; ++i) {
        gMuonTG->SetPoint(i, muonTGData[i][0], muonTGData[i][1]);
    }
}
