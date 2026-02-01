#include "Helpers.h"

#include <vector>
#include <string>
#include <algorithm>
#include "TGraph.h"

double energyLossCalculation() { return 40.; }

double energyLossCalculation(double x, double px, bool isData) {
    // x in cm and px in MeV
    double discriminant = 0.0733 * px + 1.3 * x - 31; 
    if (discriminant > 0) {
        // particles going through the halo hole
        if (isData) return 17.5;
        else return 24.5;
    } else {
        // particles going through the halo paddle
        if (isData) return 25.5;
        else return 32.5;
    }
}

bool isHitNearPrimary(
    std::vector<int>* primaryKey, 
    std::vector<float>* hitX, 
    std::vector<float>* hitW, 
    std::vector<int>* hitPlane,
    float thisHitX, 
    float thisHitW, 
    int thisHitPlane,
    float threshold,
    bool onlyVertex
) {
    int nPrimaryHits = primaryKey->size();

    if (onlyVertex) {
        int maxWireIndex = -1;
        float maxWire = -1e9;

        // Find primary hit with the highest wire on the same plane
        for (int iHit = 0; iHit < nPrimaryHits; ++iHit) {
            int key = primaryKey->at(iHit);
            if (hitPlane->at(key) != thisHitPlane) continue;

            float w = hitW->at(key);
            if (w > maxWire) {
                maxWire = w;
                maxWireIndex = key;
            }
        }

        if (maxWireIndex < 0) return false;

        float dX = std::abs(thisHitX - hitX->at(maxWireIndex));
        float dW = std::abs(thisHitW - hitW->at(maxWireIndex));
        float d  = std::sqrt(std::pow(dX, 2) + std::pow(dW, 2));

        return (d < threshold);
    }

    for (int iHit = 0; iHit < nPrimaryHits; ++iHit) {
        if (hitPlane->at(primaryKey->at(iHit)) != thisHitPlane) continue;

        float dX = std::abs(thisHitX - hitX->at(primaryKey->at(iHit)));
        float dW = std::abs(thisHitW - hitW->at(primaryKey->at(iHit)));
        float d  = std::sqrt(std::pow(dX, 2) + std::pow(dW, 2));

        if (d < threshold) return true;
    }
    return false;
}

double computeReducedChi2(const TGraph* theory, std::vector<double> xData, std::vector<double> yData, bool dataReversed, int nPoints, int nOutliersToDiscard, int nTrim) {
    std::vector<double> chi2Contributions;
    if (dataReversed) {
        std::reverse(xData.begin(), xData.end());
        std::reverse(yData.begin(), yData.end());
    }

    int startIdx = 0;
    int endIdx   = nPoints;
    if (2 * nTrim < nPoints) {
        startIdx = nTrim;
        endIdx   = nPoints - nTrim;
    }

    for (int i = startIdx; i < endIdx; ++i) {
        double theoryY = theory->Eval(xData[i]); // interpolate the theory at xData[i]
        double deltaY  = yData[i] - theoryY;
        double chi2    = (deltaY * deltaY) / std::abs(theoryY);
        chi2Contributions.push_back(chi2);
    }

    std::sort(chi2Contributions.begin(), chi2Contributions.end());
    int totalContrib = chi2Contributions.size();
    int nUsed        = std::max(0, totalContrib - nOutliersToDiscard);
    double chi2Sum   = 0.0;
    for (int i = 0; i < nUsed; ++i) {
        chi2Sum += chi2Contributions[i];
    }

    // Return reduced chi²
    return nUsed > 0 ? chi2Sum / nUsed : 0.0;
}

double meanDEDX(std::vector<double> trackDEDX, bool isTrackReversed, int pointsToUse) {
    double dEdx = 0.;
    if (isTrackReversed) std::reverse(trackDEDX.begin(), trackDEDX.end());

    unsigned int bound = pointsToUse;
    if (pointsToUse > trackDEDX.size()) bound = trackDEDX.size();
    for (unsigned int i = 0; i < bound; ++i) dEdx += trackDEDX.at(i);
    dEdx /= bound;
    return dEdx;
}

int isSecondaryInteractionAbsorption(std::vector<int> daughtersPDG, std::vector<string> daughtersProcess, std::vector<double> daughtersKE) {
    int numDaughters = daughtersPDG.size();
    if (numDaughters == 0) return 11;

    int tempNumProtons = 0;
    for (int i = 0; i < numDaughters; ++i) {
        if ((daughtersPDG[i] == 11) && (daughtersProcess[i] == "hIoni")) continue;
        if (daughtersPDG[i] == -211) return 6;
        if (daughtersPDG[i] == 111) return 7;
        if (daughtersPDG[i] == 211) return 8;
        if (daughtersProcess[i] == "Decay") return 10;
        if (daughtersProcess[i] == "hBertiniCaptureAtRest") return 9;

        if (daughtersProcess[i] == "pi-Inelastic") {
            if ((daughtersPDG[i] == 13) || (daughtersPDG[i] == -13)) return 11;
            else if ((daughtersPDG[i] == 321) || (daughtersPDG[i] == -321) || (daughtersPDG[i] == 311)) return 11;
            else if (daughtersPDG[i] == 2212) {
                if ((daughtersKE[i] >= PROTON_ENERGY_LOWER_BOUND) && (daughtersKE[i] <= PROTON_ENERGY_UPPER_BOUND)) {
                    tempNumProtons++;
                }
            }
        }
    }

    if (tempNumProtons == 0) return 0;
    return 1;
}

std::vector<double> calcLinearityProfile(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, int nb) {
    size_t N = vx.size();
    std::vector<double> linearity(N, 0.0);

    if (N < 3 || vx.size() != vy.size() || vx.size() != vz.size()) return linearity;

    for (size_t index = 0; index < N; ++index) {
        int k1 = std::max(int(0), int(index - nb));
        int k2 = std::min(int(N - 1), int(index + nb));
        int M  = k2 - k1 + 1;

        double mx = 0., my = 0., mz = 0.;
        for (int i = k1; i <= k2; ++i) { mx += vx[i]; my += vy[i]; mz += vz[i]; }
        mx /= M; my /= M; mz /= M;

        double Cxx = 0., Cxy = 0., Cxz = 0.;
        double Cyy = 0., Cyz = 0., Czz = 0.;

        for (int i = k1; i <= k2; ++i) {
            double dx = vx[i] - mx;
            double dy = vy[i] - my;
            double dz = vz[i] - mz;

            Cxx += dx * dx; Cxy += dx * dy; Cxz += dx * dz;
            Cyy += dy * dy; Cyz += dy * dz;
            Czz += dz * dz;
        }

        Cxx /= M; Cxy /= M; Cxz /= M;
        Cyy /= M; Cyz /= M;
        Czz /= M;

        TMatrixDSym cov(3);
        cov(0,0) = Cxx; cov(0,1) = Cxy; cov(0,2) = Cxz;
        cov(1,0) = Cxy; cov(1,1) = Cyy; cov(1,2) = Cyz;
        cov(2,0) = Cxz; cov(2,1) = Cyz; cov(2,2) = Czz;

        TMatrixDSymEigen eigenSolver(cov);
        TVectorD evals = eigenSolver.GetEigenValues();
        std::vector<double> eigs = { evals[0], evals[1], evals[2] };

        double lambda_sum = eigs[0] + eigs[1] + eigs[2];
        if (lambda_sum == 0.0f) {
            linearity[index] = 1.0;
        } else {
            linearity[index] = eigs[0] / lambda_sum;
        }
    }
    return linearity;
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

void Overlay_dEdx_RR_Reference_PP(
    TGraph* gProton,
    TGraph* gPion,
    TGraph* gMIP,
    bool truncate,
    double MIPStart,
    bool addLegend,
    TVirtualPad* pad
) {
    // Helper to interpolate (or extrapolate) y at a given x in a TGraph
    auto interpolateAt = [](TGraph* g, double x_target) -> double {
        int n = g->GetN();
        double* x = g->GetX();
        double* y = g->GetY();

        for (int i = 0; i < n - 1; ++i) {
            if ((x[i] <= x_target && x_target <= x[i + 1]) ||
                (x[i] >= x_target && x_target >= x[i + 1])) {
                // Linear interpolation
                double t = (x_target - x[i]) / (x[i + 1] - x[i]);
                return y[i] + t * (y[i + 1] - y[i]);
            }
        }

        // If out of bounds, return endpoint (flat extrapolation)
        return (x_target < x[0]) ? y[0] : y[n - 1];
    };

    auto truncateGraphAt = [&](TGraph* g, double x_max) -> TGraph* {
        TGraph* gTrunc = new TGraph();
        int n = g->GetN();
        double* x = g->GetX();
        double* y = g->GetY();

        for (int i = 0; i < n; ++i) {
            if (x[i] <= x_max) {
                gTrunc->SetPoint(gTrunc->GetN(), x[i], y[i]);
            } else {
                break;
            }
        }

        // Add interpolated point at x = x_max if needed
        if (x_max > x[0] && x_max < x[n - 1]) {
            double y_interp = interpolateAt(g, x_max);
            gTrunc->SetPoint(gTrunc->GetN(), x_max, y_interp);
        }

        return gTrunc;
    };

    // Clone and truncate Proton and Pion
    TGraph* gProtonGraph; TGraph* gPionGraph;
    if (truncate) {
        gProtonGraph = truncateGraphAt(gProton, MIPStart);
        gPionGraph   = truncateGraphAt(gPion,   MIPStart);
    } else {
        gProtonGraph = gProton;
        gPionGraph   = gPion;
    }

    gProtonGraph->SetLineColor(kRed + 1);
    gProtonGraph->SetLineWidth(3);
    gProtonGraph->SetName("Proton");

    gPionGraph->SetLineColor(kRed + 3);
    gPionGraph->SetLineWidth(3);
    gPionGraph->SetName("Pion");

    // Shift MIP to start at x = MIPStart
    int nMIP = gMIP->GetN();
    double* xMIP = gMIP->GetX();
    double* yMIP = gMIP->GetY();

    TGraph* gMIPShifted = new TGraph(nMIP);
    for (int i = 0; i < nMIP; ++i) {
        gMIPShifted->SetPoint(i, xMIP[i] + MIPStart, yMIP[i]);
    }

    gMIPShifted->SetLineColor(kOrange + 7);
    gMIPShifted->SetLineWidth(3);
    gMIPShifted->SetName("MIP");

    // Draw
    if (pad) pad->cd();
    gProtonGraph->Draw("L same");
    gPionGraph->Draw("L same");
    gMIPShifted->Draw("L same");

    // Legend
    if (addLegend) {
        static TLegend *leg = nullptr;
        if (!leg) {
            leg = new TLegend(0.50, 0.80, 0.75, 0.88);
            leg->AddEntry(gProtonGraph, "Proton", "l");
            leg->AddEntry(gPionGraph,   "Pion",   "l");
            leg->AddEntry(gMIPShifted,  "MIP",    "l");
            leg->SetFillStyle(1001);
            leg->SetBorderSize(0);
        }
        leg->Draw();
    }
}

double getClusterElongation(HitCluster cluster) {
    const std::vector<float>& X = cluster.hitX;
    const std::vector<float>& W = cluster.hitW;

    const std::size_t N = std::min(X.size(), W.size());
    if (N <= 2) return 0.0;

    // Means
    double meanX = 0.0, meanW = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        meanX += static_cast<double>(X[i]);
        meanW += static_cast<double>(W[i]);
    }
    meanX /= static_cast<double>(N);
    meanW /= static_cast<double>(N);

    // Covariance
    double Sxx = 0.0, Sww = 0.0, Sxw = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        const double dx = static_cast<double>(X[i]) - meanX;
        const double dw = static_cast<double>(W[i]) - meanW;
        Sxx += dx * dx;
        Sww += dw * dw;
        Sxw += dx * dw;
    }
    const double denom = static_cast<double>(N - 1);
    Sxx /= denom; Sww /= denom; Sxw /= denom;

    TMatrixDSym cov(2);
    cov(0, 0) = Sxx; cov(0, 1) = Sxw;
    cov(1, 0) = Sxw; cov(1, 1) = Sww;

    TMatrixDSymEigen eigenSolver(cov);
    TVectorD evals = eigenSolver.GetEigenValues();
    double lambda1 = evals[0]; double lambda2 = evals[1];

    // Ensure ordering: lambda1 >= lambda2
    if (lambda2 > lambda1) std::swap(lambda1, lambda2);

    return lambda2 / lambda1;
}

double getClusterWidth(HitCluster cluster) {
    const auto& X = cluster.hitX;
    const auto& W = cluster.hitW;
    const std::size_t N = (X.size() < W.size()) ? X.size() : W.size();
    if (N < 2) return 0.0;

    // Means
    double mx = 0.0, mw = 0.0;
    for (std::size_t i = 0; i < N; ++i) { mx += X[i]; mw += W[i]; }
    mx /= double(N); mw /= double(N);

    // Covariance (unbiased)
    double Sxx = 0.0, Sww = 0.0, Sxw = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        const double dx = double(X[i]) - mx;
        const double dw = double(W[i]) - mw;
        Sxx += dx*dx; Sww += dw*dw; Sxw += dx*dw;
    }
    const double d = double(N - 1);
    Sxx /= d; Sww /= d; Sxw /= d;

    // Eigenvalues of 2x2 symmetric matrix
    const double T = Sxx + Sww;
    const double D = Sxx*Sww - Sxw*Sxw;
    const double h = 0.5 * T;
    double disc2 = h*h - D;
    if (disc2 < 0.0 && disc2 > -1e-18) disc2 = 0.0; // clamp tiny negative
    const double disc = std::sqrt(disc2 > 0.0 ? disc2 : 0.0);

    const double lambda1 = h + disc;
    const double lambda2 = h - disc; // minor eigenvalue = variance ⟂ to principal axis

    // RMS width orthogonal to the principal axis
    return (lambda2 > 0.0) ? std::sqrt(lambda2) : 0.0;
}

double distance(double x1, double x2, double y1, double y2, double z1, double z2) {
    return sqrt(
        pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)
    );
}

bool isWithinReducedVolume(double x, double y, double z) {
    return (
        (x > RminX) && (x < RmaxX) && 
        (y > RminY) && (y < RmaxY) && 
        (z > RminZ) && (z < RmaxZ)
    );
}

bool isWithinActiveVolume(double x, double y, double z) {
    return (
        (x > minX) && (x < maxX) && 
        (y > minY) && (y < maxY) && 
        (z > minZ) && (z < maxZ)
    );
}

void printEfficiencyPlots(
    TString dir, int fontStyle, double textSize,
    std::vector<TEfficiency*> efficiencies,
    std::vector<TString> titles,
    std::vector<TString> xlabels
) {
    for (size_t iPlot = 0; iPlot < efficiencies.size(); ++iPlot) {
        TEfficiency* eff = efficiencies.at(iPlot);

        TCanvas* PlotCanvas = new TCanvas("Canvas","Canvas",205,34,1024,768);
        PlotCanvas->cd();
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetBottomMargin(0.15);

        eff->SetTitle(titles.at(iPlot) + ";" + xlabels.at(iPlot) + ";Efficiency");
        eff->SetLineWidth(2);
        eff->SetMarkerStyle(20);
        eff->SetMarkerSize(1.0);
        eff->Draw("AP");

        PlotCanvas->Update();
        TGraph* gr = eff->GetPaintedGraph();
        if (gr) {
            gr->GetXaxis()->SetTitleFont(fontStyle);
            gr->GetXaxis()->SetLabelFont(fontStyle);
            gr->GetXaxis()->SetNdivisions(8);
            gr->GetXaxis()->SetLabelSize(textSize);
            gr->GetXaxis()->SetTitleSize(textSize);
            gr->GetXaxis()->SetTitle(xlabels[iPlot]);
            gr->GetYaxis()->SetTitleOffset(1.1);
            gr->GetYaxis()->CenterTitle();

            gr->GetYaxis()->SetTitleFont(fontStyle);
            gr->GetYaxis()->SetLabelFont(fontStyle);
            gr->GetYaxis()->SetNdivisions(8);
            gr->GetYaxis()->SetLabelSize(textSize);
            gr->GetYaxis()->SetTitleSize(textSize);
            gr->GetYaxis()->SetTitle("Efficiency");
            gr->GetYaxis()->SetTitleOffset(1.1);
            gr->GetYaxis()->CenterTitle();
        }
        PlotCanvas->Update();

        // Save plot and delete canvas
        PlotCanvas->SaveAs(dir + titles.at(iPlot) + ".png");
        delete PlotCanvas;
    }
}

void printOneDPlots(
    const TString& dir, 
    int fontStyle, 
    double textSize,
    std::vector<std::vector<TH1*>>& groups,
    std::vector<int>& colors, 
    std::vector<std::vector<TString>>& labels, 
    std::vector<TString>& titles, 
    std::vector<TString>& xlabels, 
    std::vector<TString>& ylabels,
    std::vector<bool>& stack,
    std::vector<std::vector<bool>> asPoints,
    std::vector<std::vector<bool>> addText
) {
    int numPlots = groups.size();

    if (asPoints.empty()) {
        asPoints.resize(numPlots);
        for (size_t i = 0; i < numPlots; ++i) {
            asPoints[i].assign(groups[i].size(), false);
        }
    }

    if (addText.empty()) {
        addText.resize(numPlots);
        for (size_t i = 0; i < numPlots; ++i) {
            addText[i].assign(groups[i].size(), false);
        }
    }

    auto styleAxes = [&](TH1* h, int iPlot){
        h->SetTitle(titles.at(iPlot));
        h->GetXaxis()->SetTitleFont(fontStyle);
        h->GetXaxis()->SetLabelFont(fontStyle);
        h->GetXaxis()->SetNdivisions(8);
        h->GetXaxis()->SetLabelSize(textSize);
        h->GetXaxis()->SetTitleSize(textSize);
        h->GetXaxis()->SetTitle(xlabels.at(iPlot));
        h->GetXaxis()->SetTitleOffset(1.1);
        h->GetXaxis()->CenterTitle();

        h->GetYaxis()->SetTitleFont(fontStyle);
        h->GetYaxis()->SetLabelFont(fontStyle);
        h->GetYaxis()->SetNdivisions(6);
        h->GetYaxis()->SetLabelSize(textSize);
        h->GetYaxis()->SetTitleSize(textSize);
        h->GetYaxis()->SetTitle(ylabels.at(iPlot));
        h->GetYaxis()->SetTitleOffset(1.1);
        h->GetYaxis()->CenterTitle();
    };

    auto styleFilled = [&](TH1* h, int color){
        h->SetLineWidth(2);
        h->SetLineColor(color);
        h->SetFillColor(color);
        h->SetFillColorAlpha(color, 0.2);
        h->SetFillStyle(3001);
    };

    auto styleLineOnly = [&](TH1* h, int color){
        h->SetLineWidth(2);
        h->SetLineColor(color);
        h->SetFillStyle(0);     // no fill
        h->SetFillColor(0);
    };

    auto stylePoints = [&](TH1* h, int color){
        h->SetLineColor(color);
        h->SetMarkerColor(color);
        h->SetMarkerStyle(20);
        h->SetMarkerSize(1.0);
        h->Sumw2(kTRUE);
    };

    auto getGlobalMax = [&](const std::vector<TH1*>& vv){
        double mx = 0.0;
        for (auto* h : vv) if (h) mx = std::max(mx, h->GetMaximum());
        return mx;
    };

    auto getGlobalMin = [&](const std::vector<TH1*>& vv){
        double mi = 0.0;
        for (auto* h : vv) if (h) mi = std::min(mi, h->GetMinimum());
        return mi;
    };

    for (int iPlot = 0; iPlot < numPlots; ++iPlot) {
        TCanvas* PlotCanvas = new TCanvas(Form("Canvas_%d", iPlot), Form("Canvas_%d", iPlot), 205, 34, 1300, 768);

        TPad* mainPad = new TPad(Form("mainPad_%d", iPlot), "", 0.0, 0.0, 0.80, 1.0);
        mainPad->SetRightMargin(0.05);
        mainPad->SetLeftMargin(0.15);
        mainPad->SetBottomMargin(0.15);
        mainPad->Draw();
        mainPad->cd();

        TLegend* leg = new TLegend(0.0, 0.0, 1.0, 1.0);
        leg->SetTextSize(textSize * 2.5);
        leg->SetTextFont(fontStyle);
        leg->SetBorderSize(1);
        leg->SetLineColor(kBlack);

        TPad* legendPad = new TPad(Form("legendPad_%d", iPlot), "", 0.80, 0.4, 1.0, 0.8);
        legendPad->SetLeftMargin(0.05);
        legendPad->SetRightMargin(0.05);
        legendPad->SetBottomMargin(0.15);
        legendPad->SetFillStyle(0);
        legendPad->SetBorderMode(0);

        std::vector<TH1*> Plots = groups.at(iPlot);
        std::vector<TString> Labels = labels.at(iPlot);
        std::vector<bool> AsPts = asPoints.at(iPlot);
        std::vector<bool> AddText = addText.at(iPlot);

        if (stack.at(iPlot)) {
            // Stacked: non-points -> filled stack; points -> overlay
            THStack hs("hs", titles.at(iPlot));
            std::vector<TH1*> overlayPoints;

            for (size_t k = 0; k < Plots.size(); ++k) {
                TH1* h = Plots[k];
                if (!h) continue;

                if (AsPts[k]) {
                    stylePoints(h, colors.at(k));
                    overlayPoints.push_back(h);
                    leg->AddEntry(h, Labels[k], "lep");
                } else {
                    styleFilled(h, colors.at(k));
                    hs.Add(h, "HIST");
                    leg->AddEntry(h, Labels[k], "f");
                }
            }

            // Draw stack or create frame if only points
            TH1* frameForAxes = nullptr;
            if (hs.GetHists() && hs.GetHists()->GetSize() > 0) {
                hs.Draw("HIST");
                frameForAxes = (TH1*)hs.GetHistogram();
            } else if (!Plots.empty() && Plots[0]) {
                frameForAxes = (TH1*)Plots[0]->Clone(Form("frame_%d", iPlot));
                frameForAxes->Reset();
                frameForAxes->Draw("AXIS");
            }

            if (frameForAxes) styleAxes(frameForAxes, iPlot);

            // Y range
            double ymax = std::max(getGlobalMax(Plots), frameForAxes ? frameForAxes->GetMaximum() : 0.0);
            double ymin = std::min(getGlobalMin(Plots), frameForAxes ? frameForAxes->GetMinimum() : 0.0);

            // Expand the range with a 15% buffer
            double yr_high = 1.15 * (ymax > 0 ? ymax : 1.0);
            double yr_low  = (ymin < 0) ? 1.15 * ymin : 0.0;

            if (frameForAxes) frameForAxes->GetYaxis()->SetRangeUser(yr_low, yr_high);
            if (hs.GetHistogram()) {
                hs.SetMaximum(yr_high);
                hs.SetMinimum(yr_low);
            }

            // Re-draw stack (ROOT quirk) and overlay points
            if (hs.GetHists() && hs.GetHists()->GetSize() > 0) hs.Draw("HIST SAME");
            for (auto* h : overlayPoints) h->Draw("E1P SAME");

            gPad->Update();
            if (frameForAxes) frameForAxes->GetXaxis()->SetMaxDigits(3);
            gPad->Modified(); gPad->Update();

            PlotCanvas->cd();
            legendPad->Draw(); legendPad->cd(); leg->Draw();
            PlotCanvas->SaveAs(dir + titles.at(iPlot) + ".png");
        } else {
            // Unstacked: non-points -> line-only; points -> E1P
            // Find first non-null to set axes
            int firstIdx = -1;
            for (size_t k = 0; k < Plots.size(); ++k) if (Plots[k]) { firstIdx = (int)k; break; }
            if (firstIdx < 0) { delete PlotCanvas; continue; }

            TH1* h0 = Plots[firstIdx];
            styleAxes(h0, iPlot);

            if (AsPts[firstIdx]) {
                stylePoints(h0, colors.at(firstIdx));
                h0->Draw("E1P");
                leg->AddEntry(h0, Labels[firstIdx], "lep");
            } else {
                styleLineOnly(h0, colors.at(firstIdx));
                if (AddText[firstIdx]) h0->Draw("HIST TEXT0");
                else h0->Draw("HIST"); // line only (no fill because FillStyle=0)
                leg->AddEntry(h0, Labels[firstIdx], "l");
            }

            for (size_t k = 0; k < Plots.size(); ++k) {
                if ((int)k == firstIdx) continue;
                TH1* h = Plots[k];
                if (!h) continue;

                if (AsPts[k]) {
                    stylePoints(h, colors.at(k));
                    h->Draw("E1P SAME");
                    leg->AddEntry(h, Labels[k], "lep");
                } else {
                    styleLineOnly(h, colors.at(k));
                    if (AddText[k]) h->Draw("HIST TEXT0 SAME");
                    else h->Draw("HIST SAME"); 
                    leg->AddEntry(h, Labels[k], "l");
                }
            }

            // Unified Y range
            double ymax = getGlobalMax(Plots);
            double ymin = getGlobalMin(Plots);

            double origYmin = h0->GetYaxis()->GetXmin();
            double origYmax = h0->GetYaxis()->GetXmax();

            double yr_high = 1.15 * (ymax > 0 ? ymax : 1.0);
            double yr_low  = (ymin < 0) ? 1.15 * ymin : 0.0;

            h0->GetYaxis()->SetRangeUser(yr_low, yr_high);

            gPad->Update();
            h0->GetXaxis()->SetMaxDigits(3);
            gPad->Modified(); gPad->Update();

            PlotCanvas->cd();
            legendPad->Draw(); legendPad->cd(); leg->Draw();
            PlotCanvas->SaveAs(dir + titles.at(iPlot) + ".png");

            // Reset min and max for future plotting
            h0->SetMinimum(); h0->SetMaximum();
        }

        delete PlotCanvas; // cleans pads too
    }
}

std::pair<TH1*, TH1*> getBinEfficiencyAndPurity(TH1* hTrue, TH1* hReco, TH1* hRecoTrue) {
    TH1* hEfficiency = (TH1*) hTrue->Clone("Efficiency"); hEfficiency->Reset("ICE");
    TH1* hPurity     = (TH1*) hTrue->Clone("Purity");     hPurity->Reset("ICE");

    for (int iBin = 1; iBin <= hTrue->GetNbinsX(); ++iBin) {
        if (hTrue->GetBinContent(iBin) == 0 || hReco->GetBinContent(iBin) == 0) { hEfficiency->SetBinContent(iBin, 0); hPurity->SetBinContent(iBin, 0); continue; }
        hEfficiency->SetBinContent(iBin, hRecoTrue->GetBinContent(iBin) / hTrue->GetBinContent(iBin));
        hPurity->SetBinContent(iBin, hRecoTrue->GetBinContent(iBin) / hReco->GetBinContent(iBin));
    }
    return std::make_pair(hEfficiency, hPurity);
}

void printTwoDPlots(
    const TString& dir, 
    const std::vector<TH2*>& plots, 
    const std::vector<TString>& titles, 
    const std::vector<std::pair<double,double>>& zRanges,
    const std::vector<bool>& displayNumbers
) {
    int nPlots = plots.size();

    TCanvas* c1 = new TCanvas("c1", "EnergyLossPlots", 850, 700);
    c1->SetRightMargin(0.20);

    for (int iPlot = 0; iPlot < nPlots; ++iPlot) {
        TH1* hPlot = plots.at(iPlot);

        if (iPlot < (int) zRanges.size() && zRanges[iPlot].first < zRanges[iPlot].second) {
            hPlot->GetZaxis()->SetRangeUser(zRanges[iPlot].first, zRanges[iPlot].second);
        } else {
            double min = 0.;
            if (hPlot->GetMinimum() < min) min = hPlot->GetMinimum();
            hPlot->SetMinimum(min);
            hPlot->SetMaximum(hPlot->GetMaximum());
        }

        if (iPlot < (int)displayNumbers.size() && displayNumbers[iPlot]) {
            gStyle->SetPaintTextFormat(".2g");
            hPlot->Draw("COLZ TEXT");
        } else {
            hPlot->Draw("COLZ");
        }

        c1->SaveAs(dir + titles.at(iPlot) + ".png");
    }
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
    int elasticScattering    = background_histo->GetBinContent(13);

    // extra scattering bins
    int scattering0p         = background_histo->GetBinContent(14);
    int scatteringNp         = background_histo->GetBinContent(15);

    int total = background_histo->Integral(1, 15);

    os << "  Abs 0p: " << pionAbs0p << "    " << 100 * (pionAbs0p / (double) total) << "%" << std::endl;
    os << "  Abs Np: " << pionAbsNp << "    " << 100 * (pionAbsNp / (double) total) << "%" << std::endl;
    os << "  Primary muon: " << primaryMuon << "    " << 100 * (primaryMuon / (double) total) << "%" << std::endl;
    os << "  Primary electron: " << primaryElectron << "    " << 100 * (primaryElectron / (double) total) << "%" << std::endl;
    os << "  Other primary: " << otherPrimary << "    " << 100 * (otherPrimary / (double) total) << "%" << std::endl;
    os << "  Outside reduced volume: " << pionOutRedVol << "    " << 100 * (pionOutRedVol / (double) total) << "%" << std::endl;
    if (scattering0p + scatteringNp == 0) {
        os << "  Inelastic scattering: " << pionInelScatter << "    " << 100 * (pionInelScatter / (double) total) << "%" << std::endl;
        os << "  Elastic scattering: " << elasticScattering << "    " << 100 * (elasticScattering / (double) total) << "%" << std::endl;
    } else {
        os << "  Scattering 0p: " << scattering0p << "    " << 100 * (scattering0p / (double) total) << "%" << std::endl;
        os << "  Scattering Np: " << scatteringNp << "    " << 100 * (scatteringNp / (double) total) << "%" << std::endl;
    }
    os << "  Charge exchange: " << chargeExchange << "    " << 100 * (chargeExchange / (double) total) << "%" << std::endl;
    os << "  Double charge exchange: " << doubleChargeExchange << "    " << 100 * (doubleChargeExchange / (double) total) << "%" << std::endl;
    os << "  Capture at rest: " << captureAtRest << "    " << 100 * (captureAtRest / (double) total) << "%" << std::endl;
    os << "  Decay: " << pionDecay << "    " << 100 * (pionDecay / (double) total) << "%" << std::endl;
    os << "  Other: " << other << "    " << 100 * (other / (double) total) << "%" << std::endl;
}

void printEventInfo(EventInfo event, std::ostream& os) {
    os << "Run: " << event.run << " subrun: " << event.subrun << " event: " << event.event << std::endl;
    os << "  Background type: " << backgroundTypes[event.backgroundNum] << std::endl;
    os << "  WC match PDG: " << event.wcMatchPDG << std::endl;
    os << "  WC match process: " << event.wcMatchProcess << std::endl;
    os << "  WC match daughters: " << std::endl;
    for (int n = 0; n < event.wcMatchDaughtersPDG.size(); ++n) {
        if (event.wcMatchDaughtersPDG[n] == 11 && event.wcMatchDaughtersProcess[n] == "hIoni") continue;
        os << "      PDG: " << event.wcMatchDaughtersPDG[n];
        os << " Process: " << event.wcMatchDaughtersProcess[n];
        os << std::endl;
    }
    os << std::endl;
    os << "  Primary PDG: " << event.truthPrimaryPDG << std::endl;
    os << "  Primary vertex: " << event.vertexX << ", " << event.vertexY << ", " << event.vertexZ << std::endl; 
    os << "  Primary daughters: " << std::endl;
    for (int n = 0; n < event.truthPrimaryDaughtersPDG.size(); ++n) {
        if (event.truthPrimaryDaughtersPDG[n] == 11 && event.truthPrimaryDaughtersProcess[n] == "hIoni") continue;
        os << "      PDG: " << event.truthPrimaryDaughtersPDG[n];
        os << " Process: " << event.truthPrimaryDaughtersProcess[n];
        os << std::endl;
    }
    os << "  Secondary pion daughters: " << std::endl;
    for (int n = 0; n < event.truthSecondaryDaughtersPDG.size(); ++n) {
        if (event.truthSecondaryDaughtersPDG[n] == 11 && event.truthSecondaryDaughtersProcess[n] == "hIoni") continue;
        os << "      PDG: " << event.truthSecondaryDaughtersPDG[n];
        os << " Process: " << event.truthSecondaryDaughtersProcess[n];
        os << std::endl;
    }
    os << "  Secondary interaction tagged as: " << backgroundTypes[event.secondaryInteractionTag] << std::endl;
    os << std::endl;
    os << "  Tracks reco'ed near vertex: " << event.tracksNearVertex << std::endl;
    os << "  Tracks reco'ed as protons: " << event.recoProtonCount << std::endl;
    os << "  True visible protons: " << event.visibleProtons << std::endl;
    os << std::endl;
    os << "  Number of hit clusters found near primary track: " << event.numHitClusters << std::endl;
    os << std::endl;
    os << "  Minimum primary track linearity: " << event.minLocalLinearity << std::endl;
    os << "  Maximum primary track linearity derivative: " << event.maxLocalLinearityD << " above threshold: " << (event.maxLocalLinearityD >= LINEARITY_DERIVATIVE_THRESHOLD) << std::endl;
    os << std::endl;
}

int flattenIndex (int beta, int j, int S) { 
    return beta * S + j; 
}

std::pair<int,int> unflattenIndex(int f, int S) {
    int beta = f / S;
    int j    = f % S;
    return {beta, j};
}

void H2M(const TH2D* histo, TMatrixD& mat, bool rowcolumn) {
    // Fill 2D histogram into matrix
    // If TH2D(i, j) = Matrix(i, j), rowcolumn = kTRUE, else rowcolumn = kFALSE
    for (Int_t i=0; i<histo->GetNbinsX(); i++) {
        for (Int_t j=0; j<histo->GetNbinsY(); j++) {
            if (rowcolumn) { mat(i, j) = histo->GetBinContent(i+1, j+1); }
            else { mat(j, i) = histo->GetBinContent(i+1, j+1); }
        }
    }
}

void H2V(const TH1D* histo, TVectorD& vec) {
    if (!histo) return;
    for (Int_t i=0; i<histo->GetNbinsX(); i++) {
        vec(i) = histo->GetBinContent(i+1);
    }
}

void M2H(const TMatrixD mat, TH2D* histo) {
    for (Int_t i=0; i<mat.GetNrows(); i++) {
        for (Int_t j=0; j<mat.GetNcols(); j++) {
            histo->SetBinContent(i+1, j+1, mat(i, j));
        }
    }
}

void V2H(const TVectorD vec, TH1D* histo) {
    for (Int_t i=0; i<vec.GetNrows(); i++) {
        histo->SetBinContent(i+1, vec(i));
    }
}

int getBin(double data, const std::vector<double>& edges) {
    // edges should be sorted ascending and have size nBins+1
    if (edges.size() < 2) return -1;

    // Out of range
    if (data < edges.front() || data >= edges.back()) {
        return -1;
    }

    // Find bin index
    for (size_t i = 0; i + 1 < edges.size(); ++i) {
        if (data >= edges[i] && data < edges[i+1]) {
            return static_cast<int>(i); // 0-based index
        }
    }

    // Edge case: exactly equal to last edge
    if (data == edges.back()) {
        return static_cast<int>(edges.size() - 2);
    }

    return -1; // Shouldn’t happen
}

void reweightOneDHisto(TH1D* histo, double weight) {
    // Safety checks
    if (!histo || weight == 0) return;

    for (int i = 1; i <= histo->GetNbinsX(); ++i) {
        double content = histo->GetBinContent(i);
        double error   = histo->GetBinError(i);
        double width   = histo->GetBinWidth(i);

        // scaling factor = width / weight
        double factor = width / weight;

        if (factor != 0) {
            content /= factor;
            error   /= factor;  // scale errors consistently
        }

        histo->SetBinContent(i, content);
        histo->SetBinError(i, error);
    }
}

void removeRepeatedPoints(std::vector<double>* x, std::vector<double>* y, std::vector<double>* z) {
    int repeat_count = 0;
    int num_points = x->size();
    if (num_points > 1) {
        double last_x = x->at(num_points - 1);
        double last_y = y->at(num_points - 1);
        double last_z = z->at(num_points - 1);
        for (int j = num_points - 2; j >= 0; --j) {
            if (x->at(j) == last_x &&
                y->at(j) == last_y &&
                z->at(j) == last_z) {
                repeat_count++;
            } else {
                break;
            }
        }
    }
    
    if (repeat_count > 0) {
        x->erase(x->begin() + (num_points - repeat_count), x->end() - 1);
        y->erase(y->begin() + (num_points - repeat_count), y->end() - 1);
        z->erase(z->begin() + (num_points - repeat_count), z->end() - 1);
    }
}

void histoToText(const TH1D* hist, const std::string& filename) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: could not open file " << filename << std::endl;
        return;
    }

    int nbins = hist->GetNbinsX();
    for (int i = 1; i <= nbins; ++i) {
        double bin_low  = hist->GetBinLowEdge(i);
        double bin_high = hist->GetBinLowEdge(i + 1);
        double count    = hist->GetBinContent(i);
        outfile << bin_low << ", " << bin_high << ", " << count << "\n";
    }
    outfile.close();
}

////////////////////
// Geometry stuff //
////////////////////

static double pointToSegmentDist2(
    double px, double py, double pz, // point to check
    double ax, double ay, double az, // segment start
    double bx, double by, double bz, // segment end
    double* cx = nullptr, double* cy = nullptr, double* cz = nullptr, // closest point in AB
    double* t_out = nullptr // t in (C, or) Q = A + t AB
) {
    // Get AB = B - A
    const double abx = bx - ax;
    const double aby = by - ay;
    const double abz = bz - az;

    // Get AP = P - A
    const double apx = px - ax;
    const double apy = py - ay;
    const double apz = pz - az;

    // Get magntiude of AB
    const double ab2 = abx*abx + aby*aby + abz*abz;

    // Find t by projecting AP into AB
    double t = 0.0;
    if (ab2 > 0.0) {
        t = (apx*abx + apy*aby + apz*abz) / ab2; // projection parameter
        if (t < 0.0) t = 0.0;
        else if (t > 1.0) t = 1.0;
    } // else A==B (degenerate), keep t=0 => closest to A

    // Compute Q = A + t AB, point in AB closest to P
    const double qx = ax + t * abx;
    const double qy = ay + t * aby;
    const double qz = az + t * abz;

    // Save optional outputs
    if (cx) *cx = qx;
    if (cy) *cy = qy;
    if (cz) *cz = qz;
    if (t_out) *t_out = t;

    // Get distance between Q and P
    const double dx = px - qx;
    const double dy = py - qy;
    const double dz = pz - qz;
    return dx*dx + dy*dy + dz*dz;
}

bool IsPointInsideTrackCylinder (
    const std::vector<double>* primaryX,
    const std::vector<double>* primaryY,
    const std::vector<double>* primaryZ,
    double x, double y, double z,
    double radius,
    double* minDist2_out = nullptr,
    std::size_t* segIndex_out = nullptr
) {
    if (!primaryX || !primaryY || !primaryZ) {
        return false; // invalid inputs
    }
    const auto& X = *primaryX;
    const auto& Y = *primaryY;
    const auto& Z = *primaryZ;

    const std::size_t n = X.size();
    if (Y.size() != n || Z.size() != n || n == 0) {
        return false; // size mismatch or empty track
    }

    // If only one point, treat the "track" as a sphere around that point
    const double r2 = radius * radius;
    double best = std::numeric_limits<double>::infinity();
    std::size_t bestSeg = static_cast<std::size_t>(-1);

    if (n == 1) {
        // Just compute distance to point
        const double dx = x - X[0], dy = y - Y[0], dz = z - Z[0];
        best    = dx*dx + dy*dy + dz*dz;
        bestSeg = 0;
    } else {
        // Find minimum distance from point to segment in track
        for (std::size_t i = 0; i + 1 < n; ++i) {
            double d2 = pointToSegmentDist2(
                x, y, z,
                X[i], Y[i], Z[i],
                X[i+1], Y[i+1], Z[i+1]
            );
            if (d2 < best) {
                best = d2;
                bestSeg = i;
                if (best <= r2) {
                    // Exit if already inside
                    if (minDist2_out) *minDist2_out = best;
                    if (segIndex_out) *segIndex_out = bestSeg;
                    return true;
                }
            }
        }
    }

    if (minDist2_out) *minDist2_out = best;
    if (segIndex_out) *segIndex_out = bestSeg;
    return best <= r2;
}

bool IsPointInsideSlicedCone(
    double px, double py, double pz,
    double ax, double ay, double az,
    double vx, double vy, double vz,
    double height,
    double r0,
    double r1
) {
    // Compute direction vector norm and normalize
    double vnorm = std::sqrt(vx*vx + vy*vy + vz*vz);
    if (vnorm == 0.0 || height <= 0.0) return false;
    double dx = vx / vnorm;
    double dy = vy / vnorm;
    double dz = vz / vnorm;

    // Vector from base to point
    double apx = px - ax;
    double apy = py - ay;
    double apz = pz - az;

    // Project AP onto direction vector to get t (distance along axis)
    double t = apx*dx + apy*dy + apz*dz;

    // Check if point projects within [0, height] along axis
    if (t < 0.0 || t > height) return false;

    // Find axis point at t
    double cx = ax + t*dx;
    double cy = ay + t*dy;
    double cz = az + t*dz;

    // Compute radius at this height (linear interpolation)
    double r = r0 + (r1 - r0) * (t / height);

    // Distance squared from axis
    double dist2 = (px - cx)*(px - cx) + (py - cy)*(py - cy) + (pz - cz)*(pz - cz);

    return dist2 <= r*r;
}

std::tuple<double,double> azimuth_polar_from_points(
    double xs, double ys, double zs,
    double xe, double ye, double ze
) {
    const double dx = xe - xs;
    const double dy = ye - ys;
    const double dz = ze - zs;
    const double r  = std::sqrt(dx*dx + dy*dy + dz*dz);

    if (r == 0.0) {
        return {std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};
    }

    double phi = std::atan2(dy, dx); // (-pi, pi]

    // Clamp to avoid tiny numeric drift outside [-1,1]
    double c = dz / r;
    c = std::max(-1.0, std::min(1.0, c));
    double theta = std::acos(c); // [0, pi]

    return {phi, theta};
}

std::vector<double> projToZ(const std::vector<double>& hit0, const std::vector<double>& hit1, double zpos) {
  // This is used to find the point at z = zpos along the line created by hit0 and hit1. 
  // Uses the parameterized vector form of a line
  //
  // <x, y, z> = <sx, sy, sz> * t + <startx, starty, startz>
  //             sx, sy, sz are all slopes
  //
  // (z - startz) / sz = t
  // x = sx * t + startx
  // y = sx * t + srarty
      
  double sx = hit1[0] - hit0[0];
  double sy = hit1[1] - hit0[1];
  double sz = hit1[2] - hit0[2];

  double t = (zpos - hit0[2]) / sz;

  std::vector<double> result {sx * t + hit0[0], sy * t + hit0[1], zpos};

  return result;
}

std::vector<double> getAverageDir(const std::vector<std::vector<double>>& points) {
    std::vector<double> avgDir(3, 0.0);
    for (size_t i = 0; i < points.size() - 1; ++i) {
        avgDir[0] += (points[i + 1][0] - points[i][0]) / points.size();
        avgDir[1] += (points[i + 1][1] - points[i][1]) / points.size();
        avgDir[2] += (points[i + 1][2] - points[i][2]) / points.size();
    }
    return avgDir;
}

std::pair<size_t, size_t> find_unique_position(const std::vector<std::vector<int>>* v, int n) {
    for (size_t i = 0; i < v->size(); ++i) {
        const auto& row = (*v)[i];
        for (size_t j = 0; j < row.size(); ++j) {
            if (row[j] == n) {
                return {i, j};
            }
        }
    }

    // If you truly guarantee it exists, reaching here is a logic error:
    assert(false && "n not found in recoTrackHitIndices");
    return {0, 0}; // unreachable, but silences warnings in some builds
}

bool EqualApprox(const TMatrixD& A, const TMatrixD& B, double rtol, double atol) {
    if (A.GetNrows() != B.GetNrows() || A.GetNcols() != B.GetNcols()) return false;

    for (int i = 0; i < A.GetNrows(); ++i) {
        for (int j = 0; j < A.GetNcols(); ++j) {
            const double a = A(i,j);
            const double b = B(i,j);
            const double diff = std::fabs(a - b);
            const double scale = atol + rtol * std::max(std::fabs(a), std::fabs(b));
            if (diff > scale) return false;
        }
    }
    return true;
}

///////////////
// WC checks //
///////////////

// Magnet 1
bool CheckUpstreamMagnetAperture(const std::vector<double>& hit1, const std::vector<double>& hit2) {
   std::vector<double> USHit = projToZ(hit1, hit2, zcentMagnet1[0]);
   std::vector<double> DSHit = projToZ(hit1, hit2, zcentMagnet1[1]);
   
   if (USHit[0] < xboundMagnet1[0] || USHit[0] > xboundMagnet1[1]) return false;      // Upstream Aperture X Check
   else if (USHit[1] < yboundMagnet1[0] || USHit[1] > yboundMagnet1[1]) return false; // Upstream Aperture Y Check
   else if (DSHit[0] < xboundMagnet1[2] || DSHit[0] > xboundMagnet1[3]) return false; // Downstream Aperture X check
   else if (DSHit[1] < yboundMagnet1[2] || DSHit[1] > yboundMagnet1[3]) return false; // Downstream Aperture Y Check
   
   return true;
}

// Magnet 2
bool CheckDownstreamMagnetAperture(const std::vector<double>& hit1, const std::vector<double>& hit2) {
   std::vector<double> USHit = projToZ(hit1,hit2,zcentMagnet2[0]);
   std::vector<double> DSHit = projToZ(hit1,hit2,zcentMagnet2[1]);
   
    if (USHit[0] < xboundMagnet2[0] || USHit[0] > xboundMagnet2[1]) return false;      // Upstream Aperture X check
    else if (USHit[1] < yboundMagnet2[0] || USHit[1] > yboundMagnet2[1]) return false; // Upstream Aperture Y Check
    else if (DSHit[0] < xboundMagnet2[2] || DSHit[0] > xboundMagnet2[3]) return false; // Downstream Aperture X Check   
    else if (DSHit[1] < yboundMagnet2[2] || DSHit[1] > yboundMagnet2[3]) return false; // Downstream Aperture Y Check

   return true;
}

// Downstream collimator
bool CheckDownstreamCollimatorAperture(const std::vector<double>& hit1, const std::vector<double>& hit2) {
   std::vector<double> USHit = projToZ(hit1,hit2,zcentDSCol[0]);
   std::vector<double> DSHit = projToZ(hit1,hit2,zcentDSCol[1]);
   
   if (USHit[0] < xboundDSCol[0] || USHit[0] > xboundDSCol[1]) return false;      // Upstream Aperture X Check
   else if (USHit[1] < yboundDSCol[0] || USHit[1] > yboundDSCol[1]) return false; // Upstream Aperture Y Check
   else if (DSHit[0] < xboundDSCol[2] || DSHit[0] > xboundDSCol[3]) return false; // Downstream Aperture X Check   
   else if (DSHit[1] < yboundDSCol[2] || DSHit[1] > yboundDSCol[3]) return false; // Downstream Aperture Y Check

   return true;
}

//////////////////////////
// Wiener SVD unfolding //
// Author: Hanyu WEI    //
//////////////////////////

TMatrixD Matrix_C(Int_t n, Int_t type) {
    //Type 0: Unit matrix
    //Type I: First derivative matrix
    //Type II: Second derivative matrix
    Double_t epsilon = 1e-6; // needed for 2nd derivative matrix inversion
    Double_t epsilon2 = 1e-2; // needed for 3rd derivative matrix inversion
    TMatrixD C(n, n);
    for(Int_t i=0; i<n; i++) {
        for(Int_t j=0; j<n; j++) {
            C(i, j) = 0;
            if(type == 0) {
                if(i == j) C(i, j) = 1;
            } else if (type == 1) {
                if( j-i == 1 ) C(i, j) = 1;
                if(i==j) {
                    C(i, j) = -1;
                }
            } else if (type == 2) {
                if(TMath::Abs(i-j) == 1) C(i, j) = 1;
                if(i==j) {
                    C(i, j) = -2+epsilon;
                    if(i==0 || i==n-1) {
                        C(i, j) = -1+epsilon;
                    }
                }
            } else { // third derivative matrix 
                if(i-j == 2) { 
                    C(i, j) = -1;
                }
                if(i-j == 1) {
                    C(i, j) = 2;
                    if(i==1) {
                        C(i, j) = 0;
                    }
                }
                if(i-j == -1) {
                    C(i, j) = -2;
                    if (i==n-2) {
                        C(i, j) = 0;
                    }
                }
                if(i-j == -2) {
                    C(i, j) = 1;
                }
                if(i==j) {
                    C(i, j) = 0 + epsilon2;
                    if (i==0 || i==1) {
                        C(i, j) = 1 + epsilon2;
                    }
                    if(i==n-1 || i==n-2) {
                        C(i, j) = -1 + epsilon2;
                    }
                }
            }
        }
    }
    return C;
}

TVectorD WienerSVD(
    TMatrixD Response,     // response matrix 
    TVectorD Signal,       // true var values (all generated)
    TVectorD Measure,      // reco var values
    TMatrixD Covariance,   // cov matrix
    Int_t C_type,          // 0: unit, 1: first derivative, 2: second derivative
    Float_t Norm_type,     // norm type to exponentiate
    TMatrixD& AddSmear,    // 
    TVectorD& WF,          //
    TMatrixD& UnfoldCov,   //
    TMatrixD& CovRotation,  //
    TMatrixD& AddSmearInverse
) {
    Int_t m = Response.GetNrows(); // measure, M
    Int_t n = Response.GetNcols(); // signal, S

    // Decomposition of Covariance Matrix to get orthogonal Q to rotate the current frame, 
    // then make the uncertainty for each bin equal to 1
    TDecompSVD decV(Covariance);
    TMatrixD Q0 (TMatrixD::kTransposed, decV.GetV());
    TVectorD err0 = decV.GetSig();

    TMatrixD err(m, m);
    for(Int_t i=0; i<m; i++) {
        for(Int_t j=0; j<m; j++) {
            err(i, j) = 0;
            if (i == j) {
                if(err0(i)) err(i, j) = 1./TMath::Sqrt( err0(i) );
                else err(i, j) = 0;
            }
        }
    }

    TMatrixD Q = err*Q0;
    // Transform Measure and Response
    TVectorD M = Q*Measure;
    TMatrixD R = Q*Response;

    // For addtion of smoothness matrix, e.g. 2nd derivative C
    TMatrixD C0(n, n);
    C0 = Matrix_C(n, C_type);
    TMatrixD normsig(n, n);
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            normsig(i, j) = 0;
            if(i==j) normsig(i, j) = 1./TMath::Power(Signal(i), Norm_type);
        }
    }
    C0 = C0*normsig;
    if (C_type == 0 && Norm_type == 0) C0.UnitMatrix();

    TMatrixD C = C0;
    C0.Invert();
    TMatrixD C_inv = C0;
    Signal = C*Signal;
    R = R*C_inv;
  
    // SVD decomposition of R
    TDecompSVD udv(R);
    TMatrixD U = udv.GetU();
    TMatrixD U_t (TMatrixD::kTransposed, U);
    TMatrixD V = udv.GetV();
    TMatrixD V_t (TMatrixD::kTransposed, V);
    TVectorD D = udv.GetSig();
    // for matrix D transverse
    TMatrixD D_t(n, m);
    for(Int_t i=0; i<n; i++) {
        for(Int_t j=0; j<m; j++) {
            D_t(i,j) = 0;
            if(i==j) {
                D_t(i, j) = D(i);
            }
        }
    }

    TVectorD S = V_t*Signal;
    // Wiener Filter
    TMatrixD W(n, n);
    TMatrixD W0(n, n);
    for(Int_t i=0; i<n; i++) {
        for(Int_t j=0; j<n; j++) {
            W(i, j)  = 0;
            W0(i, j) = 0;
            if(i == j) {
                if (Norm_type == 0 && C_type == 0) {
                    W(i, j) = 1. / (D(i)*D(i)); // removed 2e-7 tolerance
                } else {
                    W(i, j) = S(i)*S(i) / ( D(i)*D(i)*S(i)*S(i)+1 );
                }
                WF(i)    = D(i)*D(i)*W(i, j);
                W0(i, j) = WF(i);
            }
        }
    }

    // Inverse of smearing matrix
    TMatrixD W0_inv = W0;
    for (int i = 0; i < W0_inv.GetNrows(); ++i) {
        W0_inv(i,i ) = 1.0 / W0(i,i);
    }
    
    TVectorD unfold(n);
    TMatrixD covRotation(n,m);

    AddSmearInverse = C_inv * V * W0_inv * V_t * C;
    AddSmear        = C_inv * V * W0 * V_t * C;
    unfold          = C_inv * V * W * D_t * U_t * M;
    covRotation     = C_inv * V * W * D_t * U_t * Q;

    CovRotation = covRotation;
    TMatrixD covRotation_t (TMatrixD::kTransposed, covRotation); 
    UnfoldCov = covRotation * Covariance * covRotation_t;  

    return unfold;
}

///////////////////////////
// Get variables for BDT //
///////////////////////////

void getBDTVariables(
    int WC2TPCtrkID,
    const std::vector<double>* wcLocationsX,
    const std::vector<double>* wcLocationsY,
    const std::vector<double>* wcLocationsZ,
    const std::vector<double>* wcX,
    const std::vector<double>* wcY,
    const std::vector<double>* wcZ,
    const std::vector<double>* recoBeginX,
    const std::vector<double>* recoBeginY,
    const std::vector<double>* recoBeginZ,
    const std::vector<double>* recoEndX,
    const std::vector<double>* recoEndY,
    const std::vector<double>* recoEndZ,
    const std::vector<std::vector<double>>* recoDEDX,
    const std::vector<int>*    recoTrkID,
    int maxRecoTrks,

    Float_t* bdt_WC2TPCTrackLength,
    Float_t* bdt_WC2TPCBeginX,
    Float_t* bdt_WC2TPCBeginY,
    Float_t* bdt_WC2TPCBeginZ,
    Float_t* bdt_WC2TPCEndX,
    Float_t* bdt_WC2TPCEndY,
    Float_t* bdt_WC2TPCEndZ,
    Float_t* bdt_WC2TPCdEdx,

    Float_t* bdt_recoTrkBeginX,
    Float_t* bdt_recoTrkBeginY,
    Float_t* bdt_recoTrkBeginZ,
    Float_t* bdt_recoTrkEndX,
    Float_t* bdt_recoTrkEndY,
    Float_t* bdt_recoTrkEndZ,
    Float_t* bdt_recoTrkLen,
    Float_t* bdt_recoTrkdEdx,
    Float_t* bdt_numRecoTrksInCylinder
) {
    if (
        wcLocationsX && 
        wcLocationsY && 
        wcLocationsZ &&
        wcLocationsX->size() > 0 &&
        wcLocationsY->size() > 0 &&
        wcLocationsZ->size() > 0
    ) {
        // Track length: sum of segment distances
        double length = 0.0;
        for (size_t i = wcLocationsX->size() - 1; i > 0; --i) {
            double dx = wcLocationsX->at(i) - wcLocationsX->at(i - 1);
            double dy = wcLocationsY->at(i) - wcLocationsY->at(i - 1);
            double dz = wcLocationsZ->at(i) - wcLocationsZ->at(i - 1);
            length += std::sqrt(dx * dx + dy * dy + dz * dz);
        }
        *bdt_WC2TPCTrackLength = length;

        *bdt_WC2TPCBeginX = wcLocationsX->at(0);
        *bdt_WC2TPCBeginY = wcLocationsY->at(0);
        *bdt_WC2TPCBeginZ = wcLocationsZ->at(0);

        *bdt_WC2TPCEndX = wcLocationsX->at(wcLocationsX->size() - 1);
        *bdt_WC2TPCEndY = wcLocationsY->at(wcLocationsY->size() - 1);
        *bdt_WC2TPCEndZ = wcLocationsZ->at(wcLocationsZ->size() - 1);
    } else {
        *bdt_WC2TPCTrackLength = -9999;
        *bdt_WC2TPCBeginX = -9999;
        *bdt_WC2TPCBeginY = -9999;
        *bdt_WC2TPCBeginZ = -9999;
        *bdt_WC2TPCEndX = -9999;
        *bdt_WC2TPCEndY = -9999;
        *bdt_WC2TPCEndZ = -9999;
        *bdt_WC2TPCdEdx = -9999;
    }

    // Clean up entries
    for (int i = 0; i < maxRecoTrks; ++i) {
        bdt_recoTrkBeginX[i] = -9999; bdt_recoTrkBeginY[i] = -9999; bdt_recoTrkBeginZ[i] = -9999;
        bdt_recoTrkEndX[i]   = -9999; bdt_recoTrkEndY[i]   = -9999; bdt_recoTrkEndZ[i]   = -9999;
        bdt_recoTrkLen[i]    = -9999; bdt_recoTrkdEdx[i]   = -9999;
    }

    int trks_in_cylinder = 0;
    std::vector<std::pair<int, double>> trk_idx_len;

    for (int trk_idx = 0; trk_idx < recoTrkID->size(); ++trk_idx) {
        std::vector<double> thisTrackDEDX;
        if (recoDEDX && recoDEDX->size() > trk_idx) {
            thisTrackDEDX = recoDEDX->at(trk_idx);
        }

        if (recoTrkID->at(trk_idx) == WC2TPCtrkID) {
            if (thisTrackDEDX.size() > 0) {
                *bdt_WC2TPCdEdx = (Float_t) std::accumulate(thisTrackDEDX.begin(), thisTrackDEDX.end(), 0.0) / thisTrackDEDX.size();
            } else {
                *bdt_WC2TPCdEdx = -9999;
            }
            continue;
        }

        double bx = recoBeginX->at(trk_idx);
        double by = recoBeginY->at(trk_idx);
        double bz = recoBeginZ->at(trk_idx);
        double ex = recoEndX->at(trk_idx);
        double ey = recoEndY->at(trk_idx);
        double ez = recoEndZ->at(trk_idx);

        bool startInCylinder = IsPointInsideTrackCylinder(
            wcX, wcY, wcZ,
            bx, by, bz,
            CYLINDER_RADIUS
        );
        bool endInCylinder = IsPointInsideTrackCylinder(
            wcX, wcY, wcZ,
            ex, ey, ez,
            CYLINDER_RADIUS
        );

        double trk_len = std::sqrt((bx - ex) * (bx - ex) + (by - ey) * (by - ey) + (bz - ez) * (bz - ez));

        if (startInCylinder && endInCylinder) {
            trks_in_cylinder += 1;
            trk_idx_len.emplace_back(trk_idx, trk_len);
        }
    }

    // Sort tracks by length, descending
    std::sort(trk_idx_len.begin(), trk_idx_len.end(),
    [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
        return a.second > b.second;
    });

    for (int i = 0; i < std::min<int>(maxRecoTrks, (int) trk_idx_len.size()); ++i) {
        int idx = trk_idx_len[i].first;

        std::vector<double> thisTrackDEDX;
        if (recoDEDX && recoDEDX->size() > idx) {
            thisTrackDEDX = recoDEDX->at(idx);
        }
        bdt_recoTrkBeginX[i] = (Float_t)(*recoBeginX)[idx];
        bdt_recoTrkBeginY[i] = (Float_t)(*recoBeginY)[idx];
        bdt_recoTrkBeginZ[i] = (Float_t)(*recoBeginZ)[idx];
        bdt_recoTrkEndX[i]   = (Float_t)(*recoEndX)[idx];
        bdt_recoTrkEndY[i]   = (Float_t)(*recoEndY)[idx];
        bdt_recoTrkEndZ[i]   = (Float_t)(*recoEndZ)[idx];
        bdt_recoTrkLen[i]    = (Float_t)trk_idx_len[i].second;
        bdt_recoTrkdEdx[i]   = (thisTrackDEDX.size() > 0) ? (Float_t)(std::accumulate(thisTrackDEDX.begin(), thisTrackDEDX.end(), 0.0) / thisTrackDEDX.size()) : -9999;
    }
    *bdt_numRecoTrksInCylinder = trks_in_cylinder;
}

/////////////////
// Systematics //
/////////////////

std::vector<TH1D*> MakeUniverseHists(
    const std::string& baseName,
    const std::string& title,
    int nBins,
    const double* binEdges,
    int nUniverses
) {
    std::vector<TH1D*> hists;
    hists.reserve(nUniverses);

    for (int i = 0; i < nUniverses; ++i) {
        auto h = std::make_unique<TH1D>(
            (baseName + "_" + std::to_string(i)).c_str(),
            title.c_str(),
            nBins,
            binEdges
        );
        hists.push_back(h.release());
    }

    return hists;
}

std::vector<std::vector<std::vector<TH1D*>>> MakeUniverseHistBlock(
    const std::string& baseName,        // e.g., "hTrueAbs0p"
    const std::string& titleFmt,        // e.g., "Abs0p (U=%d): type %d in bin %d"
    int nUniverses,
    int nOuter,                         // e.g., NUM_BINS_KE
    int nInner,                         // e.g., NUM_SIGNAL_TYPES
    int nbins,
    const double* binEdges
) {
    std::vector<std::vector<std::vector<TH1D*>>> block;
    block.resize(nUniverses);

    for (int u = 0; u < nUniverses; ++u) {
        block[u].resize(nOuter);

        for (int i = 0; i < nOuter; ++i) {
            block[u][i].reserve(nInner);

            for (int j = 0; j < nInner; ++j) {
                // Unique histogram name
                std::string name = Form(
                    "%s_U%d_%d_%d", 
                    baseName.c_str(), u, i, j
                );

                // Title with formatting
                std::string title = Form(
                    titleFmt.c_str(), u, j, i
                );

                auto h = std::make_unique<TH1D>(
                    name.c_str(),
                    title.c_str(),
                    nbins,
                    binEdges
                );

                block[u][i].push_back(h.release());
            }
        }
    }

    return block;
}

void GetResponseMatrix(
    int SizeOuter, int SizeInner,
    const std::vector<TH1D*>& TotalEventsHistos,
    const std::vector<std::vector<std::vector<TH1D*>>>& TrueRecoAsByBin,
    const TH1D* RecoIncidentFlux,
    const TH1D* TrueIncidentFlux,
    TH2D* ReponseMatrix
) {
    for (int iOuterSignalBin = 0; iOuterSignalBin < SizeOuter; ++iOuterSignalBin) {
        for (int iOuterEnergyBin = 0; iOuterEnergyBin < SizeInner; ++iOuterEnergyBin) {
            int row      = flattenIndex(iOuterSignalBin, iOuterEnergyBin, SizeInner);
            double denom = TotalEventsHistos.at(iOuterSignalBin)->GetBinContent(iOuterEnergyBin + 1);

            for (int iInnerSignalBin = 0; iInnerSignalBin < SizeOuter; ++iInnerSignalBin) {
                for (int iInnerEnergyBin = 0; iInnerEnergyBin < SizeInner; ++iInnerEnergyBin) {
                    int  column = flattenIndex(iInnerSignalBin, iInnerEnergyBin, SizeInner);
                    double prob = 0;
                    
                    // Get total number of true (i, \alpha) events that were reconstructed as (j, \beta)
                    if (denom > 0) prob = TrueRecoAsByBin.at(iOuterSignalBin).at(iOuterEnergyBin).at(iInnerSignalBin)->GetBinContent(iInnerEnergyBin + 1) / denom;

                    // Scale for different incident flux between energy bins
                    prob = prob * (TrueIncidentFlux->GetBinContent(iOuterEnergyBin + 1) / RecoIncidentFlux->GetBinContent(iInnerEnergyBin + 1));

                    ReponseMatrix->SetBinContent(column + 1, row + 1, prob);
                }
            }
        }
    }
}

void GetCovMatrix(
    TH1* RecoHisto,
    std::vector<TH1D*> UnivRecoHisto,
    TH2* CovMatrix
) {
    int NUniv = UnivRecoHisto.size();

    // Create covariance matrices
    int n = RecoHisto->GetNbinsX();
    for (int iUniv = 0; iUniv < NUniv; ++iUniv) {
        for (int x = 1; x < n + 1; x++) {
            double XEventRateCV = RecoHisto->GetBinContent(x);
            double XEventRateVar = UnivRecoHisto[iUniv]->GetBinContent(x);
            for (int y = 1; y < n + 1; y++) {
                double YEventRateCV = RecoHisto->GetBinContent(y);
                double YEventRateVar = UnivRecoHisto[iUniv]->GetBinContent(y);

                double Value = ((XEventRateVar - XEventRateCV) * (YEventRateVar - YEventRateCV)) / NUniv;
                if (TMath::Abs(Value) < 1e-14) Value = 1e-14;

                // if (x == y) {
                //     std::cout << "(" << x << ", " << y << ")" << std::endl;
                //     std::cout << XEventRateVar << "  " << XEventRateCV << "  " << XEventRateVar - XEventRateCV << std::endl;
                //     std::cout << YEventRateVar << "  " << YEventRateCV << "  " << YEventRateVar - YEventRateCV << std::endl;
                //     std::cout << Value << std::endl;
                //     std::cout << std::endl;
                // }

                CovMatrix->Fill(
                    RecoHisto->GetXaxis()->GetBinCenter(x),
                    RecoHisto->GetXaxis()->GetBinCenter(y),
                    Value
                );
            }
        }
    }
}

void GetFracCovAndCorrMatrix(
    TH1* RecoHisto,
    TH2* CovMatrix,
    TH2* FracCovMatrix,
    TH2* CorrMatrix
) {
    int n = RecoHisto->GetNbinsX();
    for (int x = 1; x < n + 1; x++) {
        double XEventRateCV = RecoHisto->GetBinContent(x);
        for (int y = 1; y < n + 1; y++) {
            double YEventRateCV = RecoHisto->GetBinContent(y);
            double CovBinValue = CovMatrix->GetBinContent(x,y);
            double XBinValue = CovMatrix->GetBinContent(x,x);
            double YBinValue = CovMatrix->GetBinContent(y,y);

            // Fill frac cov matrix
            double FracValue = (XEventRateCV == 0. || YEventRateCV == 0.) ? 1e-14 : CovBinValue / (XEventRateCV * YEventRateCV);
            if (TMath::Abs(FracValue) < 1e-14) FracValue = 1e-14;
            FracCovMatrix->SetBinContent(x, y, FracValue);

            // Fill corr matrix
            double CorrValue = (XBinValue == 0. || YBinValue == 0.) ? 1e-14 : CovBinValue / (TMath::Sqrt(XBinValue) * TMath::Sqrt(YBinValue));
            if (TMath::Abs(CorrValue) < 1e-14) CorrValue = 1e-14;
            CorrMatrix->SetBinContent(x, y, CorrValue);
        }
    }
}

double FindQuantile(double frac, std::vector<double>& xs_in) {
    // This turns out to be a much more fraught issue than you would naively
    // expect. This algorithm is equivalent to R-6 here:
    // https://en.wikipedia.org/wiki/Quantile#Estimating_quantiles_from_a_sample

    // In principle we could use std::nth_element(). Probably doesn't matter
    // much in practice since this is only for plotting.
    std::vector<double> xs = xs_in;
    std::sort(xs.begin(), xs.end());

    const int N = xs.size();
    // The index we would ideally be sampling at
    const double h = frac*(N+1);
    // The indices on either side where we have to actually evaluate
    const unsigned int h0 = std::floor(h);
    const unsigned int h1 = std::ceil(h);
    if(h0 == 0) return xs[0]; // Don't underflow indexing
    if(h1 > xs.size()) return xs.back(); // Don't overflow indexing
    // The values at those indices
    const double x0 = xs[h0-1]; // wikipedia is using 1-based indexing
    const double x1 = xs[h1-1];

    if(h0 == h1) return x0;

    // Linear interpolation
    return (h1-h)*x0 + (h-h0)*x1;
}

void DrawHistosWithErrorBands(
    TH1D* RecoHisto,
    std::vector<TH1D*> UnivRecoHisto,
    const TString& dir,
    const TString& SystName,
    const TString& PlotName,
    double textSize,
    int fontStyle,
    const TString& title,
    const TString& xlabel,
    const TString& ylabel
) {
    // Basic sanity
    if (!RecoHisto) return;
    const int NUniv = (int) UnivRecoHisto.size();
    if (NUniv == 0) return;

    auto styleAxes = [&](TH1* h){
        h->SetTitle(title);

        h->GetXaxis()->SetTitleFont(fontStyle);
        h->GetXaxis()->SetLabelFont(fontStyle);
        h->GetXaxis()->SetNdivisions(8);
        h->GetXaxis()->SetLabelSize(textSize);
        h->GetXaxis()->SetTitleSize(textSize);
        h->GetXaxis()->SetTitle(xlabel);
        h->GetXaxis()->SetTitleOffset(1.1);
        h->GetXaxis()->CenterTitle();

        h->GetYaxis()->SetTitleFont(fontStyle);
        h->GetYaxis()->SetLabelFont(fontStyle);
        h->GetYaxis()->SetNdivisions(6);
        h->GetYaxis()->SetLabelSize(textSize);
        h->GetYaxis()->SetTitleSize(textSize);
        h->GetYaxis()->SetTitle(ylabel);
        h->GetYaxis()->SetTitleOffset(1.1);
        h->GetYaxis()->CenterTitle();
    };

    auto styleLineOnly = [&](TH1* h, int color){
        h->SetLineWidth(2);
        h->SetLineColor(color);
        h->SetFillStyle(0);
        h->SetFillColor(0);
    };

    auto getGlobalMax = [&](TH1* h){
        return h ? h->GetMaximum() : 0.0;
    };

    auto getGlobalMin = [&](TH1* h){
        return h ? h->GetMinimum() : 0.0;
    };

    TCanvas* PlotCanvas = new TCanvas("Selection","Selection", 205, 34, 1300, 768);

    TPad* mainPad = new TPad("mainPad", "", 0.0, 0.0, 0.80, 1.0);
    mainPad->SetRightMargin(0.05);
    mainPad->SetLeftMargin(0.15);
    mainPad->SetBottomMargin(0.15);
    mainPad->Draw();
    mainPad->cd();

    TPad* legendPad = new TPad("legendPad", "", 0.80, 0.4, 1.0, 0.8);
    legendPad->SetLeftMargin(0.05);
    legendPad->SetRightMargin(0.05);
    legendPad->SetBottomMargin(0.15);
    legendPad->SetFillStyle(0);
    legendPad->SetBorderMode(0);

    TLegend* leg = new TLegend(0.0, 0.0, 1.0, 1.0);
    leg->SetTextSize(textSize * 2.5);
    leg->SetTextFont(fontStyle);
    leg->SetBorderSize(1);
    leg->SetLineColor(kBlack);

    styleAxes(RecoHisto);
    styleLineOnly(RecoHisto, kBlue + 2);

    const int n = RecoHisto->GetNbinsX();
    TGraphAsymmErrors* RecoErrorBand = new TGraphAsymmErrors;
    RecoErrorBand->Set(n);

    double ymaxBand = -1e300;
    for (int binIdx = 1; binIdx <= n; ++binIdx) {
        const double x = RecoHisto->GetXaxis()->GetBinCenter(binIdx);
        const double y = RecoHisto->GetBinContent(binIdx);
        RecoErrorBand->SetPoint(binIdx, x, y);

        const double dx = RecoHisto->GetXaxis()->GetBinWidth(std::min(std::max(binIdx, 1), n));

        std::vector<double> ys;
        ys.reserve(NUniv);
        for (int iUniv = 0; iUniv < NUniv; ++iUniv) {
            if (!UnivRecoHisto[iUniv]) continue;
            ys.push_back(UnivRecoHisto[iUniv]->GetBinContent(binIdx));
        }

        double y0 = y, y1 = y;
        if (!ys.empty()) {
            y0 = FindQuantile(.5 - 0.6827/2, ys);
            y1 = FindQuantile(.5 + 0.6827/2, ys);
        }

        ymaxBand = std::max(ymaxBand, std::max(y, y1));

        RecoErrorBand->SetPointError(
            binIdx,
            dx/2, dx/2,
            std::max(y - y0, 0.0),
            std::max(y1 - y, 0.0)
        );
    }

    // Legend entry
    leg->AddEntry(RecoHisto, "Reconstructed", "l");

    const double ymax = std::max(getGlobalMax(RecoHisto), ymaxBand);
    const double ymin = getGlobalMin(RecoHisto);

    const double yr_high = 1.15 * (ymax > 0 ? ymax : 1.0);
    const double yr_low  = (ymin < 0) ? 1.15 * ymin : 0.0;
    RecoHisto->GetYaxis()->SetRangeUser(yr_low, yr_high);

    // --- Draw ---
    RecoHisto->Draw("HIST");
    DrawErrorBand(RecoHisto, RecoErrorBand, -1, 1);

    gPad->Update();
    RecoHisto->GetXaxis()->SetMaxDigits(3);
    gPad->Modified(); gPad->Update();

    PlotCanvas->cd();
    legendPad->Draw();
    legendPad->cd();
    leg->Draw();

    PlotCanvas->SaveAs(dir + SystName + "/" + PlotName + ".png");

    delete PlotCanvas;
    delete leg;
    delete RecoErrorBand;
}

void DrawErrorBand(TH1* nom, TGraphAsymmErrors* band, int bandCol, double alpha) {
    if(bandCol == -1) bandCol = nom->GetLineColor()-10; // hopefully a lighter version

    // Check if this pad has already been drawn in
    const bool same = gPad && !gPad->GetListOfPrimitives()->IsEmpty();

    nom->Draw(same ? "hist same" : "hist");

    band->SetFillColorAlpha(bandCol, alpha);
    band->Draw("e2 same");

    nom->Draw("hist same");

    // If we are the first spectrum to draw, scale the y-axis appropriately to
    // fit the error band as well as the nominal
    if(!same){
    double maxY = 0;
    // Don't consider underflow or overflow bins when determining maximum
    for(int i = 1; i < band->GetN()-1; ++i){
        maxY = std::max(maxY, band->GetY()[i] + band->GetErrorYhigh(i));
    }

    // Use non-zero lower value so that SetLogy() still works
    nom->GetYaxis()->SetRangeUser(1e-10, 1.1 * maxY);
    }
    gPad->RedrawAxis();
}

void PrintTruthUnfoldedStatShapePlot(
    const TString& SaveDir,
    const TString& Name,
    TH1* hTrue,
    TH1* hUnfolded,
    TH1* hStat,
    TH1* hShape,
    const TString& PlotTitle,
    const TString& XTitle,
    const TString& YTitle,
    int FontStyle,
    double TextSize
) {
    if (!hTrue || !hUnfolded || !hStat || !hShape) {
        std::cerr << "PrintTruthUnfoldedStatShapePlot: null input histogram.\n";
        return;
    }

    const int nb = hUnfolded->GetNbinsX();
    if (hTrue->GetNbinsX() != nb || hStat->GetNbinsX() != nb || hShape->GetNbinsX() != nb) {
        std::cerr << "PrintTruthUnfoldedStatShapePlot: WARNING nbins differ. "
                  << "Using unfolded nbins = " << nb << " for errors.\n";
    }

    // Make "bands" as vertical error bars (no fill, just I-shaped bars)
    TGraphAsymmErrors* StatErrorBand  = new TGraphAsymmErrors(nb);
    TGraphAsymmErrors* ShapeErrorBand = new TGraphAsymmErrors(nb);

    // Optional: if you also want total stat⊕shape as a separate band
    TGraphAsymmErrors* TotErrorBand   = new TGraphAsymmErrors(nb);

    for (int iBin = 1; iBin <= nb; ++iBin) {
        const double xnom = hUnfolded->GetXaxis()->GetBinCenter(iBin);
        const double ynom = hUnfolded->GetBinContent(iBin);

        const double s_stat  = std::max(0.0, hStat->GetBinContent(iBin));   // sigma_stat
        const double s_shape = std::max(0.0, hShape->GetBinContent(iBin));  // sigma_shape
        const double s_tot   = std::sqrt(s_stat*s_stat + s_shape*s_shape);  // stat⊕shape

        const int ip = iBin - 1; // TGraph point index is 0-based

        StatErrorBand->SetPoint(ip, xnom, ynom);
        StatErrorBand->SetPointError(ip, 0.0, 0.0, s_stat, s_stat);

        // ShapeErrorBand shows the *outer* envelope; here I interpret "shape" as
        // the additional component on top of stat, so the outer bar is total.
        ShapeErrorBand->SetPoint(ip, xnom, ynom);
        ShapeErrorBand->SetPointError(ip, 0.0, 0.0, s_tot, s_tot);

        TotErrorBand->SetPoint(ip, xnom, ynom);
        TotErrorBand->SetPointError(ip, 0.0, 0.0, s_tot, s_tot);
    }

    // Canvas
    static int ic = 0;
    TCanvas* c = new TCanvas(Form("c_ts_%d", ic++), "", 1400, 1000);
    c->SetLeftMargin(0.14);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.14);
    c->SetTopMargin(0.10);
    c->SetTicks(1,1);

    gStyle->SetGridStyle(3); // 3 = dotted
    gStyle->SetGridWidth(1);
    c->SetGridx(1);
    c->SetGridy(1);

    // Style: truth (red step)
    hTrue->SetLineColor(kRed);
    hTrue->SetLineWidth(3);
    hTrue->SetFillStyle(0);

    // Style: unfolded points
    hUnfolded->SetMarkerStyle(20);
    hUnfolded->SetMarkerSize(1.0);
    hUnfolded->SetMarkerColor(kBlack);
    hUnfolded->SetLineColor(kBlack);

    // Style: error bars (vertical "I" bars)
    // Outer (stat⊕shape): thicker
    ShapeErrorBand->SetLineColor(kBlack);
    ShapeErrorBand->SetLineWidth(2);
    ShapeErrorBand->SetMarkerStyle(1);
    ShapeErrorBand->SetMarkerSize(0.0);

    // Inner (stat): thinner
    StatErrorBand->SetLineColor(kBlack);
    StatErrorBand->SetLineWidth(1);
    StatErrorBand->SetMarkerStyle(1);
    StatErrorBand->SetMarkerSize(0.0);

    gStyle->SetEndErrorSize(4);

    // Style title
    hTrue->SetTitle(PlotTitle);
    gStyle->SetTitleFont(FontStyle, "t");
    gStyle->SetTitleSize(TextSize * 1.1, "t");
    gStyle->SetTitleAlign(23); // centered, top
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleY(0.98);

    hTrue->GetXaxis()->SetTitle(XTitle);
    hTrue->GetXaxis()->CenterTitle(true);
    hTrue->GetXaxis()->SetTitleFont(FontStyle);
    hTrue->GetXaxis()->SetLabelFont(FontStyle);
    hTrue->GetXaxis()->SetTitleSize(TextSize);
    hTrue->GetXaxis()->SetLabelSize(TextSize * 0.85);
    hTrue->GetXaxis()->SetTitleOffset(1.15);

    hTrue->GetYaxis()->SetTitle(YTitle);
    hTrue->GetYaxis()->CenterTitle(true);
    hTrue->GetYaxis()->SetTitleFont(FontStyle);
    hTrue->GetYaxis()->SetLabelFont(FontStyle);
    hTrue->GetYaxis()->SetTitleSize(TextSize);
    hTrue->GetYaxis()->SetLabelSize(TextSize * 0.85);
    hTrue->GetYaxis()->SetTitleOffset(1.0);

    // Plot title positioning
    hTrue->SetTitleFont(FontStyle, "t");
    hTrue->SetTitleSize(TextSize * 1.05, "t");
    hTrue->SetTitleOffset(0.95, "t");
    
    // Major + minor tick density
    hTrue->GetXaxis()->SetNdivisions(7, 5, 0, kTRUE);
    hTrue->GetYaxis()->SetNdivisions(5, 5, 0, kTRUE);

    hTrue->GetXaxis()->SetTickLength(0.02);
    hTrue->GetYaxis()->SetTickLength(0.02);

    // Set Y range to include outer bars
    double ymax = 0.0;
    for (int iBin = 1; iBin <= nb; ++iBin) {
        const double y_true = hTrue->GetBinContent(iBin);
        const double y_unf  = hUnfolded->GetBinContent(iBin);
        const double s_stat  = std::max(0.0, hStat->GetBinContent(iBin));
        const double s_shape = std::max(0.0, hShape->GetBinContent(iBin));
        const double s_tot   = std::sqrt(s_stat*s_stat + s_shape*s_shape);
        ymax = std::max(ymax, std::max(y_true, y_unf + s_tot));
    }
    if (ymax <= 0) ymax = 1.0;
    hTrue->SetMaximum(1.30 * ymax);
    hTrue->SetMinimum(0.0);

    // Draw order: truth -> outer bars -> inner bars -> points
    hTrue->Draw("HIST");
    ShapeErrorBand->Draw("E1 SAME");
    StatErrorBand->Draw("E1 SAME");
    hUnfolded->Draw("P SAME");

    // Legend: use the graphs with option "e" so the legend glyph is a vertical error bar
    TLegend* leg = new TLegend(0.20, 0.75, 0.55, 0.85);
    leg->SetBorderSize(0);
    leg->SetNColumns(2);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->SetTextFont(FontStyle);
    leg->SetTextSize(TextSize * 0.75);

    leg->AddEntry(hTrue, "True", "l");
    leg->AddEntry(hUnfolded, "Unfolded", "p");
    leg->AddEntry(ShapeErrorBand, "Stat#oplusShape", "e");
    leg->Draw();

    c->Update();
    c->Print(SaveDir + Name + ".png");

    delete leg;
    delete c;
    delete StatErrorBand;
    delete ShapeErrorBand;
    delete TotErrorBand;
}


void PrintFDPlot(
    const TString& SaveDir,
    const TString& Name,
    TH1* hTrue,
    TH1* hNomTrue,
    TH1* hUnfolded,
    TH1* hSysts,
    const TString& PlotTitle,
    const TString& XTitle,
    const TString& YTitle,
    std::pair<double, double> chi,
    std::pair<int, int> ndof,
    std::pair<double, double> pval,
    std::pair<double, double> sigma,
    int FontStyle,
    double TextSize
) {
    const int nb = hUnfolded->GetNbinsX();

    // Optional: if you also want total stat⊕shape as a separate band
    TGraphAsymmErrors* TotErrorBand   = new TGraphAsymmErrors(nb);

    for (int iBin = 1; iBin <= nb; ++iBin) {
        const double xnom = hUnfolded->GetXaxis()->GetBinCenter(iBin);
        const double ynom = hUnfolded->GetBinContent(iBin);

        const double s_tot = hSysts->GetBinContent(iBin);

        const int ip = iBin - 1; // TGraph point index is 0-based

        TotErrorBand->SetPoint(ip, xnom, ynom);
        TotErrorBand->SetPointError(ip, 0.0, 0.0, s_tot, s_tot);
    }

    // Canvas
    static int ic = 0;
    TCanvas* c = new TCanvas(Form("c_ts_%d", ic++), "", 1400, 1000);
    c->SetLeftMargin(0.14);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.14);
    c->SetTopMargin(0.10);
    c->SetTicks(1,1);

    gStyle->SetGridStyle(3); // 3 = dotted
    gStyle->SetGridWidth(1);
    c->SetGridx(1);
    c->SetGridy(1);

    // Style: truth (red step)
    hTrue->SetLineColor(kRed);
    hTrue->SetLineWidth(3);
    hTrue->SetFillStyle(0);

    // Style: nom truth (blue step)
    hNomTrue->SetLineColor(kBlue);
    hNomTrue->SetLineWidth(3);
    hNomTrue->SetFillStyle(0);

    // Style: unfolded points
    hUnfolded->SetMarkerStyle(20);
    hUnfolded->SetMarkerSize(1.0);
    hUnfolded->SetMarkerColor(kBlack);
    hUnfolded->SetLineColor(kBlack);

    // Style: error bars (vertical "I" bars)
    TotErrorBand->SetLineColor(kBlack);
    TotErrorBand->SetLineWidth(2);
    TotErrorBand->SetMarkerStyle(1);
    TotErrorBand->SetMarkerSize(0.0);

    gStyle->SetEndErrorSize(4);

    // Style title
    hTrue->SetTitle(PlotTitle);
    gStyle->SetTitleFont(FontStyle, "t");
    gStyle->SetTitleSize(TextSize * 1.1, "t");
    gStyle->SetTitleAlign(23); // centered, top
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleY(0.98);

    hTrue->GetXaxis()->SetTitle(XTitle);
    hTrue->GetXaxis()->CenterTitle(true);
    hTrue->GetXaxis()->SetTitleFont(FontStyle);
    hTrue->GetXaxis()->SetLabelFont(FontStyle);
    hTrue->GetXaxis()->SetTitleSize(TextSize);
    hTrue->GetXaxis()->SetLabelSize(TextSize * 0.85);
    hTrue->GetXaxis()->SetTitleOffset(1.15);

    hTrue->GetYaxis()->SetTitle(YTitle);
    hTrue->GetYaxis()->CenterTitle(true);
    hTrue->GetYaxis()->SetTitleFont(FontStyle);
    hTrue->GetYaxis()->SetLabelFont(FontStyle);
    hTrue->GetYaxis()->SetTitleSize(TextSize);
    hTrue->GetYaxis()->SetLabelSize(TextSize * 0.85);
    hTrue->GetYaxis()->SetTitleOffset(1.0);

    // Plot title positioning
    hTrue->SetTitleFont(FontStyle, "t");
    hTrue->SetTitleSize(TextSize * 1.05, "t");
    hTrue->SetTitleOffset(0.95, "t");
    
    // Major + minor tick density
    hTrue->GetXaxis()->SetNdivisions(7, 5, 0, kTRUE);
    hTrue->GetYaxis()->SetNdivisions(5, 5, 0, kTRUE);

    hTrue->GetXaxis()->SetTickLength(0.02);
    hTrue->GetYaxis()->SetTickLength(0.02);

    // Set Y range to include outer bars
    double ymax = 0.0;
    for (int iBin = 1; iBin <= nb; ++iBin) {
        const double y_true = hTrue->GetBinContent(iBin);
        const double y_unf  = hUnfolded->GetBinContent(iBin);
        const double s_tot   = std::max(0.0, hSysts->GetBinContent(iBin));
        ymax = std::max(ymax, std::max(y_true, y_unf + s_tot));
    }
    if (ymax <= 0) ymax = 1.0;
    hTrue->SetMaximum(1.35 * ymax);
    hTrue->SetMinimum(0.0);

    // Draw order: truth -> outer bars -> inner bars -> points
    hTrue->Draw("HIST");
    hNomTrue->Draw("HIST SAME");
    TotErrorBand->Draw("E1 SAME");
    hUnfolded->Draw("PE1 SAME");

    // Legend: use the graphs with option "e" so the legend glyph is a vertical error bar
    TLegend* leg = new TLegend(0.15, 0.70, 0.55, 0.85);
    leg->SetBorderSize(0);
    // leg->SetNColumns(2);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->SetTextFont(FontStyle);
    leg->SetTextSize(TextSize * 0.75);

    std::ostringstream oss1, oss2;
    oss1 << std::fixed << std::setprecision(1);
    oss2 << std::fixed << std::setprecision(1);

    oss1 << "Fake True ("
        << chi.first << "/" << ndof.first
        << ", " << pval.first
        << ", " << (sigma.first >= 8 ? ">8.0" : Form("%.2g", sigma.first)) 
        << "#sigma)";

    oss2 << "Nominal True ("
        << chi.second << "/" << ndof.second
        << ", " << pval.second
        << ", " << (sigma.second >= 8 ? ">8.0" : Form("%.2g", sigma.second)) 
        << "#sigma)";

    std::string fakeTrue = oss1.str();
    std::string nominalTrue = oss2.str();

    leg->AddEntry(hTrue, fakeTrue.c_str(), "l");
    leg->AddEntry(hNomTrue, nominalTrue.c_str(), "l");
    leg->AddEntry(hUnfolded, "Unfolded", "pe");
    leg->Draw();

    c->Update();
    c->Print(SaveDir + Name + ".png");

    delete leg;
    delete c;
    delete TotErrorBand;
}

void CalcChiSquared(
    TH1D* h_model, // true
    TH1D* h_data,  // unfolded
    TH2D* cov,     // unfolded cov
    double &chi, 
    int &ndof, 
    double &pval, 
    double &sigma
) {
	// Clone them so we can scale them 
	TH1D* h_model_clone = (TH1D*)h_model->Clone();
	TH1D* h_data_clone  = (TH1D*)h_data->Clone();
	TH2D* h_cov_clone   = (TH2D*)cov->Clone();
	int NBins = h_cov_clone->GetNbinsX();

	// Getting covariance matrix in TMatrix form
	TMatrixD cov_m;
	cov_m.Clear();
	cov_m.ResizeTo(NBins,NBins);

	// loop over rows
	for (int i = 0; i < NBins; i++) {			
		// loop over columns
		for (int j = 0; j < NBins; j++) {
            cov_m[i][j] = h_cov_clone->GetBinContent(i+1, j+1);
		}
	}
	TMatrixD copy_cov_m = cov_m;

	// Eigen-decomposition
    TMatrixDEigen eig(cov_m);
    TMatrixD evals  = eig.GetEigenValues();
    TMatrixD evecs  = eig.GetEigenVectors();

    double maxEval = 0.0;
    for (int k = 0; k < NBins; ++k) {
        if (evals[k][k] > maxEval) maxEval = evals[k][k];
    }

    // Build diagonal matrix of 1/λ with cutoff
    double tol = 1e-4 * maxEval;
    TMatrixD Dinv(NBins, NBins);
    Dinv.Zero();

    int nModesKept = 0;
    for (int k = 0; k < NBins; ++k) {
        if (evals[k][k] > tol) {
            Dinv(k,k) = 1.0 / evals[k][k];
            ++nModesKept;
        } else {
            // Regularize: treat this mode as unconstrained -> 0 in inverse
            Dinv(k,k) = 0.0;
        }
    }

    // Pseudo-inverse: V * Dinv * V^T
    TMatrixD evecsT(TMatrixD::kTransposed, evecs);
    TMatrixD inverse_cov_m = evecs * Dinv * evecsT;

	// Calculating the chi2 = Summation_ij{ (x_i - mu_j)*E_ij^(-1)*(x_j - mu_j)  }
	// x = data, mu = model, E^(-1) = inverted covariance matrix 
	chi = 0.;
	
	for (int i = 0; i < NBins; i++) {
		//double XWidth = h_data_clone->GetBinWidth(i+1);
		for (int j = 0; j < NBins; j++) {
			//double YWidth = h_data_clone->GetBinWidth(i+1);
			double diffi = h_data_clone->GetBinContent(i+1) - h_model_clone->GetBinContent(i+1);
			double diffj = h_data_clone->GetBinContent(j+1) - h_model_clone->GetBinContent(j+1);
			double LocalChi = diffi * inverse_cov_m[i][j] * diffj;
			chi += LocalChi;
		}
	}
	// ndof  = h_data_clone->GetNbinsX();
    ndof = nModesKept;
    pval = TMath::Prob(chi, ndof);
    if (pval < 1e-16) {
        sigma = 8.0;
    } else {
        sigma = TMath::Sqrt(TMath::ChisquareQuantile(1 - pval, 1));
    }

	delete h_model_clone;
	delete h_data_clone;
	delete h_cov_clone;
}