#include "Helpers.h"

#include <vector>
#include <string>
#include <algorithm>
#include "TGraph.h"

int getCorrespondingBin(double value, int num_bins, double low, double high) {
    double bin_size = (high - low) / (double) num_bins;
    return std::ceil(value / bin_size); // 1 is first bin, 2 second, etc.
}

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

bool isHitNearPrimary(std::vector<int>* primaryKey, std::vector<float>* hitX, std::vector<float>* hitW, float thisHitX, float thisHitW, float xThreshold, float wThreshold) {
    int nPrimaryHits = primaryKey->size();
    for (int iHit = 0; iHit < nPrimaryHits; ++iHit) {
        float dX = std::abs(thisHitX - hitX->at(primaryKey->at(iHit)));
        float dW = std::abs(thisHitW - hitW->at(primaryKey->at(iHit)));
        if ((dX < xThreshold) && (dW < wThreshold)) return true;
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
    std::vector<std::vector<bool>> asPoints
) {
    int numPlots = groups.size();

    if (asPoints.empty()) {
        asPoints.resize(numPlots);
        for (size_t i = 0; i < numPlots; ++i) {
            asPoints[i].assign(groups[i].size(), false);
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
            double yr = 1.15 * (ymax > 0 ? ymax : 1.0);
            if (frameForAxes) frameForAxes->GetYaxis()->SetRangeUser(0., yr);
            if (hs.GetHistogram()) hs.SetMaximum(yr);

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
                h0->Draw("HIST");                 // line only (no fill because FillStyle=0)
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
                    h->Draw("HIST SAME");
                    leg->AddEntry(h, Labels[k], "l");
                }
            }

            // Unified Y range
            double ymax = getGlobalMax(Plots);
            double yr = 1.15 * (ymax > 0 ? ymax : 1.0);
            h0->GetYaxis()->SetRangeUser(0., yr);

            gPad->Update();
            h0->GetXaxis()->SetMaxDigits(3);
            gPad->Modified(); gPad->Update();

            PlotCanvas->cd();
            legendPad->Draw(); legendPad->cd(); leg->Draw();
            PlotCanvas->SaveAs(dir + titles.at(iPlot) + ".png");
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

    TCanvas* c1 = new TCanvas("c1", "EnergyLossPlots", 800, 600);
    for (int iPlot = 0; iPlot < nPlots; ++iPlot) {
        TH1* hPlot = plots.at(iPlot);

        if (iPlot < (int) zRanges.size() && zRanges[iPlot].first < zRanges[iPlot].second) {
            hPlot->GetZaxis()->SetRangeUser(zRanges[iPlot].first, zRanges[iPlot].second);
        } else {
            hPlot->SetMinimum(0);
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
    int scattering0p         = background_histo->GetBinContent(14);
    int scatteringNp         = background_histo->GetBinContent(15);

    os << "  Abs 0p: " << pionAbs0p << std::endl;
    os << "  Abs Np: " << pionAbsNp << std::endl;
    os << "  Primary muon: " << primaryMuon << std::endl;
    os << "  Primary electron: " << primaryElectron << std::endl;
    os << "  Other primary: " << otherPrimary << std::endl;
    os << "  Outside reduced volume: " << pionOutRedVol << std::endl;
    if (scattering0p + scatteringNp == 0) {
        os << "  Inelastic scattering: " << pionInelScatter << std::endl;
        os << "  Elastic scattering: " << elasticScattering << std::endl;
    } else {
        os << "  Scattering 0p: " << scattering0p << std::endl;
        os << "  Scattering Np: " << scatteringNp << std::endl;
    }
    os << "  Charge exchange: " << chargeExchange << std::endl;
    os << "  Double charge exchange: " << doubleChargeExchange << std::endl;
    os << "  Capture at rest: " << captureAtRest << std::endl;
    os << "  Decay: " << pionDecay << std::endl;
    os << "  Other: " << other << std::endl;
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

int getBin(double data, double lower, double upper, int nBins) {
    if (data < lower || data >= upper) {
        return -1;
    }

    double binWidth = (upper - lower) / nBins;
    int bin = static_cast<int>((data - lower) / binWidth);

    // Safety check in case of edge case at upper bound
    if (bin >= nBins) bin = nBins - 1;

    return bin; // 0-based index
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
    TMatrixD& CovRotation  //
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
    // transform Measure and Response
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
                //W(i, j) = 1./(D(i)*D(i)+2e-7); //S(i)*S(i) / ( D(i)*D(i)*S(i)*S(i)+1 );
                //WF(i) = D(i)*D(i)*W(i, j);//S(i)*S(i) / ( D(i)*D(i)*S(i)*S(i)+1 );
                W(i, j) = S(i)*S(i) / ( D(i)*D(i)*S(i)*S(i)+1 );
                WF(i) = D(i)*D(i)*W(i, j);
                W0(i, j) = WF(i); 
            }
        }
    }

    TVectorD unfold = C_inv*V*W*D_t*U_t*M;
    AddSmear = C_inv*V*W0*V_t*C;

    // covariance matrix of the unfolded spectrum
    TMatrixD covRotation = C_inv*V*W*D_t*U_t*Q;
    CovRotation = covRotation;
    TMatrixD covRotation_t (TMatrixD::kTransposed, covRotation); 
    UnfoldCov = covRotation*Covariance*covRotation_t;  

    return unfold;
}