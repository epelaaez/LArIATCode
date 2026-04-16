#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TVector3.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"
#include "TEfficiency.h"
#include "TMVA/Reader.h"

#include <vector>

#include "Helpers.cc"

void DataUnfold() {
    // Set defaults
    gStyle->SetOptStat(0); // get rid of stats box
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    TH1::AddDirectory(false);
    gStyle->SetPalette(kRainBow);
    gStyle->SetPaintTextFormat("4.2f");

    int FontStyle = 132;
    double TextSize = 0.06;
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis/figs/DataUnfold/";

    int N = NUM_SIGNAL_TYPES * NUM_BINS_KE;

    // Load data histograms
    std::unique_ptr<TFile> DataFile(TFile::Open("/exp/lariat/app/users/epelaez/histos/data/Reco.root"));
    TH1D* MeasureData    = (TH1D*) DataFile->Get("hMeasureVector");
    TH2D* StatCovariance = (TH2D*) DataFile->Get("hStatCovariance");

    // Get histograms from MC
    std::unique_ptr<TFile> NominalFile(TFile::Open("/exp/lariat/app/users/epelaez/histos/nominal/Reco.root"));
    TH1D* BackgroundNominal = (TH1D*) NominalFile->Get("hBackgroundVectorNominal");
    TH1D* TrueNominal       = (TH1D*) NominalFile->Get("hSignalVectorNominal");
    TH2D* ResponseNominal   = (TH2D*) NominalFile->Get("hResponseMatrix");

    // Get covariances from MC
    std::unique_ptr<TFile> GeneratorFile(TFile::Open("/exp/lariat/app/users/epelaez/histos/generator/Measure.root"));
    TH2D* GeneratorCovariance = (TH2D*) GeneratorFile->Get("hMeasureCovMatrix");

    std::unique_ptr<TFile> EnergyRecoFile(TFile::Open("/exp/lariat/app/users/epelaez/histos/energyreco/Measure.root"));
    TH2D* EnergyRecoCovariance = (TH2D*) EnergyRecoFile->Get("hMeasureCovMatrix");

    std::unique_ptr<TFile> BeamlineMuFile(TFile::Open("/exp/lariat/app/users/epelaez/histos/beamlinemu/Measure.root"));
    TH2D* BeamlineMuCovariance = (TH2D*) BeamlineMuFile->Get("hMeasureCovMatrix");

    std::unique_ptr<TFile> BeamlineElFile(TFile::Open("/exp/lariat/app/users/epelaez/histos/beamlineel/Measure.root"));
    TH2D* BeamlineElCovariance = (TH2D*) BeamlineElFile->Get("hMeasureCovMatrix");

    std::unique_ptr<TFile> MCStatFile(TFile::Open("/exp/lariat/app/users/epelaez/histos/mcstat/Measure.root"));
    TH2D* MCStatCovariance = (TH2D*) MCStatFile->Get("hMeasureCovMatrix");

    // Convert histograms into matrices/vectors
    TVectorD Measure(N); H2V(MeasureData, Measure);
    TVectorD Background(N); H2V(BackgroundNominal, Background);
    TVectorD TrueSignal(N); H2V(TrueNominal, TrueSignal);
    TMatrixD Response(N, N); H2M(ResponseNominal, Response, kTRUE);

    // ============================================================
    // DEBUG: Print input vectors for AbsNp bins
    // ============================================================
    std::cout << "\n========== DEBUG: INPUT VECTORS (AbsNp bins) ==========\n";
    for (int i = 0; i < N; ++i) {
        auto [signalBin, energyBin] = unflattenIndex(i, NUM_BINS_KE);
        if (signalBin == 1) { // AbsNp
            std::cout << "  Flat bin " << i << " (AbsNp, E-bin " << energyBin << "):"
                      << "  Measure=" << Measure(i)
                      << "  Background=" << Background(i)
                      << "  Meas-Bkg=" << Measure(i) - Background(i)
                      << "  TrueSignal=" << TrueSignal(i)
                      << std::endl;
        }
    }

    // DEBUG: Print response matrix rows for AbsNp
    std::cout << "\n========== DEBUG: RESPONSE MATRIX (rows = AbsNp bins) ==========\n";
    for (int i = 0; i < N; ++i) {
        auto [signalBin, energyBin] = unflattenIndex(i, NUM_BINS_KE);
        if (signalBin == 1) { // AbsNp row
            double rowSum = 0;
            std::cout << "  Row " << i << " (AbsNp, E-bin " << energyBin << "): ";
            for (int j = 0; j < N; ++j) {
                if (std::fabs(Response(i, j)) > 1e-12) {
                    auto [sj, ej] = unflattenIndex(j, NUM_BINS_KE);
                    std::cout << "[col " << j << " (s" << sj << ",e" << ej << ")=" << Response(i, j) << "] ";
                }
                rowSum += Response(i, j);
            }
            std::cout << "  rowSum=" << rowSum << std::endl;
        }
    }

    // DEBUG: Print response matrix columns for AbsNp (what maps INTO AbsNp true bins)
    std::cout << "\n========== DEBUG: RESPONSE MATRIX (cols = AbsNp bins) ==========\n";
    for (int j = 0; j < N; ++j) {
        auto [signalBin, energyBin] = unflattenIndex(j, NUM_BINS_KE);
        if (signalBin == 1) { // AbsNp column
            double colSum = 0;
            std::cout << "  Col " << j << " (AbsNp, E-bin " << energyBin << "): ";
            for (int i = 0; i < N; ++i) {
                if (std::fabs(Response(i, j)) > 1e-12) {
                    auto [si, ei] = unflattenIndex(i, NUM_BINS_KE);
                    std::cout << "[row " << i << " (s" << si << ",e" << ei << ")=" << Response(i, j) << "] ";
                }
                colSum += Response(i, j);
            }
            std::cout << "  colSum=" << colSum << std::endl;
        }
    }

    // DEBUG: Check diagonal dominance for AbsNp
    std::cout << "\n========== DEBUG: RESPONSE DIAGONAL vs OFF-DIAGONAL (AbsNp) ==========\n";
    for (int i = 0; i < N; ++i) {
        auto [signalBin, energyBin] = unflattenIndex(i, NUM_BINS_KE);
        if (signalBin == 1) {
            double diag = Response(i, i);
            double offDiagSum = 0;
            for (int j = 0; j < N; ++j) {
                if (j != i) offDiagSum += std::fabs(Response(i, j));
            }
            std::cout << "  Bin " << i << " (AbsNp, E-bin " << energyBin << "):"
                      << "  diag=" << diag
                      << "  |off-diag| sum=" << offDiagSum
                      << "  ratio diag/offdiag=" << (offDiagSum > 0 ? diag/offDiagSum : 999)
                      << std::endl;
        }
    }

    // Organize covariance matrices
    std::vector<TH2D*> CovarianceMatrices = {
        StatCovariance,
        MCStatCovariance, 
        GeneratorCovariance,
        EnergyRecoCovariance,
        BeamlineMuCovariance,
        BeamlineElCovariance
    };

    std::vector<std::string> CovarianceLabels = {
        "Stat",
        "MCStat",
        "Generator",
        "EnergyReco",
        "BeamlineMu",
        "BeamlineEl"
    };

    std::vector<TMatrixD> Covariances; Covariances.reserve(CovarianceMatrices.size());
    TMatrixD TotalCovariance(N, N); TotalCovariance.Zero();
    for (int iCov = 0; iCov < CovarianceMatrices.size(); ++iCov) {
        TMatrixD Cov(N, N);
        H2M(CovarianceMatrices[iCov], Cov, kTRUE);
        Covariances.push_back(Cov);
        TotalCovariance += Cov;
    }

    // DEBUG: Print total covariance diagonal for AbsNp bins
    std::cout << "\n========== DEBUG: TOTAL COVARIANCE DIAGONAL (AbsNp) ==========\n";
    for (int i = 0; i < N; ++i) {
        auto [signalBin, energyBin] = unflattenIndex(i, NUM_BINS_KE);
        if (signalBin == 1) {
            std::cout << "  Bin " << i << " (AbsNp, E-bin " << energyBin << "):"
                      << "  TotalCov(" << i << "," << i << ")=" << TotalCovariance(i, i)
                      << "  sqrt=" << std::sqrt(std::max(0.0, TotalCovariance(i, i)))
                      << std::endl;
        }
    }

    // Vectors to store output
    TMatrixD AddSmear(N, N);
    TMatrixD AddSmearInverse(N, N);
    TVectorD WF(N);
    TMatrixD UnfoldCov(N, N);
    TMatrixD CovRotation(N, N);

    // Subtract background from measured
    TVectorD MeasureMinusBackground = Measure - Background;

    // Unfold using the total covariance matrix
    TVectorD Unfolded = WienerSVD(
        Response,
        TrueSignal,
        MeasureMinusBackground,
        TotalCovariance,
        0,
        0,
        AddSmear,
        WF,
        UnfoldCov, 
        CovRotation,
        AddSmearInverse
    );

    // ============================================================
    // DEBUG: Print unfolded results for AbsNp
    // ============================================================
    std::cout << "\n========== DEBUG: UNFOLDED RESULT (AbsNp bins) ==========\n";
    for (int i = 0; i < N; ++i) {
        auto [signalBin, energyBin] = unflattenIndex(i, NUM_BINS_KE);
        if (signalBin == 1) {
            double sigma = std::sqrt(std::max(0.0, UnfoldCov(i, i)));
            std::cout << "  Bin " << i << " (AbsNp, E-bin " << energyBin << "):"
                      << "  Unfolded=" << Unfolded(i)
                      << "  sigma=" << sigma
                      << "  Unf/True=" << (TrueSignal(i) != 0 ? Unfolded(i)/TrueSignal(i) : 0)
                      << (Unfolded(i) < 0 ? "  *** NEGATIVE ***" : "")
                      << std::endl;
        }
    }

    // DEBUG: Print Wiener filter values for AbsNp
    std::cout << "\n========== DEBUG: WIENER FILTER (AbsNp bins) ==========\n";
    for (int i = 0; i < N; ++i) {
        auto [signalBin, energyBin] = unflattenIndex(i, NUM_BINS_KE);
        if (signalBin == 1) {
            std::cout << "  Bin " << i << " (AbsNp, E-bin " << energyBin << "):"
                      << "  WF=" << WF(i)
                      << std::endl;
        }
    }

    // DEBUG: Print CovRotation matrix row for AbsNp (how measured bins contribute to unfolded AbsNp)
    std::cout << "\n========== DEBUG: COV ROTATION ROW (AbsNp bins) ==========\n";
    std::cout << "  Shows how each measured bin contributes to unfolded AbsNp bins\n";
    for (int i = 0; i < N; ++i) {
        auto [signalBin, energyBin] = unflattenIndex(i, NUM_BINS_KE);
        if (signalBin == 1) {
            std::cout << "  Row " << i << " (AbsNp, E-bin " << energyBin << "): ";
            // Print the largest (positive and negative) contributions
            double posSum = 0, negSum = 0;
            int posCnt = 0, negCnt = 0;
            for (int j = 0; j < N; ++j) {
                double val = CovRotation(i, j);
                if (val > 0) { posSum += val; posCnt++; }
                else { negSum += val; negCnt++; }
            }
            std::cout << "posSum=" << posSum << " (" << posCnt << " entries)"
                      << "  negSum=" << negSum << " (" << negCnt << " entries)"
                      << std::endl;

            // Print top 5 positive and negative entries
            std::vector<std::pair<double, int>> entries;
            for (int j = 0; j < N; ++j) entries.push_back({CovRotation(i,j), j});
            std::sort(entries.begin(), entries.end());
            std::cout << "    Most negative: ";
            for (int k = 0; k < std::min(5, (int)entries.size()); ++k) {
                auto [sj, ej] = unflattenIndex(entries[k].second, NUM_BINS_KE);
                std::cout << "[col " << entries[k].second << " (s" << sj << ",e" << ej << ")=" << entries[k].first << "] ";
            }
            std::cout << std::endl;
            std::cout << "    Most positive: ";
            for (int k = entries.size()-1; k >= std::max(0, (int)entries.size()-5); --k) {
                auto [sj, ej] = unflattenIndex(entries[k].second, NUM_BINS_KE);
                std::cout << "[col " << entries[k].second << " (s" << sj << ",e" << ej << ")=" << entries[k].first << "] ";
            }
            std::cout << std::endl;
        }
    }

    // DEBUG: Decompose the unfolded value: CovRotation * MeasureMinusBackground
    std::cout << "\n========== DEBUG: DECOMPOSITION of Unfolded AbsNp ==========\n";
    std::cout << "  Unfolded(i) = sum_j CovRotation(i,j) * MeasMinusBkg(j)\n";
    for (int i = 0; i < N; ++i) {
        auto [signalBin, energyBin] = unflattenIndex(i, NUM_BINS_KE);
        if (signalBin == 1) {
            double posContrib = 0, negContrib = 0;
            std::cout << "  Bin " << i << " (AbsNp, E-bin " << energyBin << "):\n";
            for (int j = 0; j < N; ++j) {
                double contrib = CovRotation(i, j) * MeasureMinusBackground(j);
                if (std::fabs(contrib) > 1e-6) {
                    auto [sj, ej] = unflattenIndex(j, NUM_BINS_KE);
                    std::cout << "    j=" << j << " (s" << sj << ",e" << ej << ")"
                              << "  CovRot=" << CovRotation(i,j)
                              << "  MeasMinBkg=" << MeasureMinusBackground(j)
                              << "  contrib=" << contrib
                              << std::endl;
                }
                if (contrib > 0) posContrib += contrib;
                else negContrib += contrib;
            }
            std::cout << "    => posContrib=" << posContrib
                      << "  negContrib=" << negContrib
                      << "  total=" << posContrib + negContrib
                      << "  (Unfolded=" << Unfolded(i) << ")"
                      << std::endl;
        }
    }

    // DEBUG: Print AddSmear matrix for AbsNp (additional smearing from regularization)
    std::cout << "\n========== DEBUG: ADDITIONAL SMEARING MATRIX (AbsNp rows) ==========\n";
    for (int i = 0; i < N; ++i) {
        auto [signalBin, energyBin] = unflattenIndex(i, NUM_BINS_KE);
        if (signalBin == 1) {
            std::cout << "  Row " << i << " (AbsNp, E-bin " << energyBin << "): diag=" << AddSmear(i,i);
            double offDiag = 0;
            for (int j = 0; j < N; ++j) {
                if (j != i) offDiag += std::fabs(AddSmear(i,j));
            }
            std::cout << "  |off-diag| sum=" << offDiag << std::endl;
        }
    }

    std::cout << "\n========== END DEBUG ==========\n\n";

    TH1D* hUnfolded = new TH1D("hUnfolded", "hUnfolded;;", N, 0, N); V2H(Unfolded, hUnfolded);
    TMatrixD CovRotationT(TMatrixD::kTransposed, CovRotation);

    // Unfold each covariance matrix
    std::vector<TMatrixD> UnfCovSources;
    UnfCovSources.reserve(Covariances.size());

    TMatrixD UnfCovSum(N, N); UnfCovSum.Zero();
    for (int iCov = 0; iCov < (int) Covariances.size(); ++iCov) {
        TMatrixD Ck_unf = CovRotation * Covariances[iCov] * CovRotationT;
        UnfCovSources.push_back(Ck_unf);
        UnfCovSum += Ck_unf;
    }
    TMatrixD UnfCovFromTotalRotate = CovRotation * TotalCovariance * CovRotationT;

    // Sanity check
    if (!EqualApprox(UnfoldCov, UnfCovSum, 1e-6, 1e-6)) {
        std::cerr << "Warning: UnfoldCov (from WienerSVD) != sum_k(U Ck U^T) within tolerance.\n";
    }

    /////////////////////
    // Organize output //
    /////////////////////

    // Unfolded total uncertainty
    TH1D* hUncTotal     = new TH1D("hUncTotal", "hUncTotal;Bin (i, #alpha);Cross section [barn] per 100 MeV", N, 0, N);
    TH1D* hUncStat      = new TH1D("hUncStat", "hUncStat;Bin (i, #alpha);Cross section [barn] per 100 MeV", N, 0, N);
    TH1D* hUncShape     = new TH1D("hUncShape", "hUncShape;Bin (i, #alpha);Cross section [barn] per 100 MeV", N, 0, N);
    TH1D* hFracUncTotal = new TH1D("hFracUncTotal", "hFracUncTotal;Bin (i, #alpha);Cross section unc. [%]", N, 0, N);

    // Total uncertainty per interaction type
    std::vector<TH1D*> UncTotal; UncTotal.reserve(UnfCovSources.size());
    std::vector<TH1D*> FracUncTotal; FracUncTotal.reserve(UnfCovSources.size());

    // Unfolded vector by interaction type
    TH1D* hUnfAbs0p   = new TH1D("hUnfAbs0p", "hUnfAbs0p;Kinetic Energy [MeV];Cross section [barn] per 100 MeV", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hUnfAbsNp   = new TH1D("hUnfAbsNp", "hUnfAbsNp;Kinetic Energy [MeV];Cross section [barn] per 100 MeV", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hUnfScatter = new TH1D("hUnfScatter", "hUnfScatter;Kinetic Energy [MeV];Cross section [barn] per 100 MeV", NUM_BINS_KE, ARRAY_KE_BINS.data());
    std::vector<TH1D*> UnfPerInteraction = {hUnfAbs0p, hUnfAbsNp, hUnfScatter};

    TH1D* Abs0pUncTotal   = new TH1D("Abs0pUncTotal", "Abs0pUncTotal;Kinetic Energy [MeV];Cross section [barn] per 100 MeV", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* AbsNpUncTotal   = new TH1D("AbsNpUncTotal", "AbsNpUncTotal;Kinetic Energy [MeV];Cross section [barn] per 100 MeV", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* ScatterUncTotal = new TH1D("ScatterUncTotal", "ScatterUncTotal;Kinetic Energy [MeV];Cross section [barn] per 100 MeV", NUM_BINS_KE, ARRAY_KE_BINS.data());
    std::vector<TH1D*> TotalUncPerInteraction = {Abs0pUncTotal, AbsNpUncTotal, ScatterUncTotal};

    TH1D* Abs0pUncStat  = new TH1D("Abs0pUncStat", "Abs0pUncStat;Kinetic Energy [MeV];Cross section [barn] per 100 MeV", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* AbsNpUncStat  = new TH1D("AbsNpUncStat", "AbsNpUncStat;Kinetic Energy [MeV];Cross section [barn] per 100 MeV", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* ScatterUncStat  = new TH1D("ScatterUncStat", "ScatterUncStat;Kinetic Energy [MeV];Cross section [barn] per 100 MeV", NUM_BINS_KE, ARRAY_KE_BINS.data());
    std::vector<TH1D*> StatUncPerInteraction = {Abs0pUncStat, AbsNpUncStat, ScatterUncStat};

    TH1D* Abs0pUncShape = new TH1D("Abs0pUncShape", "Abs0pUncShape;Kinetic Energy [MeV];Cross section [barn] per 100 MeV", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* AbsNpUncShape = new TH1D("AbsNpUncShape", "AbsNpUncShape;Kinetic Energy [MeV];Cross section [barn] per 100 MeV", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* ScatterUncShape = new TH1D("ScatterUncShape", "ScatterUncShape;Kinetic Energy [MeV];Cross section [barn] per 100 MeV", NUM_BINS_KE, ARRAY_KE_BINS.data());
    std::vector<TH1D*> ShapeUncPerInteraction = {Abs0pUncShape, AbsNpUncShape, ScatterUncShape};

    TH1D* Abs0pFracUncTotal   = new TH1D("Abs0pFracUncTotal", "Abs0pFracUncTotal;Kinetic Energy [MeV];Cross section [barn] per 100 MeV", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* AbsNpFracUncTotal   = new TH1D("AbsNpFracUncTotal", "AbsNpFracUncTotal;Kinetic Energy [MeV];Cross section [barn] per 100 MeV", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* ScatterFracUncTotal = new TH1D("ScatterFracUncTotal", "ScatterFracUncTotal;Kinetic Energy [MeV];Cross section [barn] per 100 MeV", NUM_BINS_KE, ARRAY_KE_BINS.data());
    std::vector<TH1D*> TotalFracUncPerInteraction = {Abs0pFracUncTotal, AbsNpFracUncTotal, ScatterFracUncTotal};

    // Uncertainty per source and interaction type
    std::vector<std::vector<TH1D*>> SourceUncPerInteraction(NUM_SIGNAL_TYPES);
    std::vector<std::vector<TH1D*>> SourceFracUncPerInteraction(NUM_SIGNAL_TYPES);

    // Populate histograms with total uncertainty
    for (int i = 0; i < N; ++i) {
        auto [signalBin, energyBin] = unflattenIndex(i, NUM_BINS_KE);

        const double x   = Unfolded(i);
        const double var = UnfoldCov(i,i);

        const double sigma   = std::sqrt(std::max(0.0, var));
        const double fracPct = (std::fabs(x) > 0.0) ? (100.0 * sigma / std::fabs(x)) : 0.0;

        const double statvar  = UnfCovSources[0](i,i);
        double shapevar = 0.0;
        for (int j = 1; j < (int) UnfCovSources.size(); ++j) {
            shapevar += UnfCovSources[j](i,i);
        }
        const double statsigma  = std::sqrt(std::max(0.0, statvar));
        const double shapesigma = std::sqrt(std::max(0.0, shapevar));

        std::cout << "Bin " << i << ": x = " << x << ", sigma = " << sigma << ", fracPct = " << fracPct << std::endl;
        // std::cout << "   " << "stat sigma = " << statsigma << ", shape sigma = " << shapesigma << std::endl;
        // std::cout << "   stat + shape in quadrature = " << std::sqrt(statsigma*statsigma + shapesigma*shapesigma) << std::endl;
        // std::cout << "   stat + shape linearly = " << (statsigma + shapesigma) << std::endl;

        hUncTotal->SetBinContent(i + 1, sigma);
        hFracUncTotal->SetBinContent(i + 1, fracPct);

        hUncStat->SetBinContent(i + 1, statsigma);
        hUncShape->SetBinContent(i + 1, shapesigma);

        UnfPerInteraction[signalBin]->SetBinContent(energyBin + 1, x);
        UnfPerInteraction[signalBin]->SetBinError(energyBin + 1, sigma);

        TotalUncPerInteraction[signalBin]->SetBinContent(energyBin + 1, sigma);
        TotalFracUncPerInteraction[signalBin]->SetBinContent(energyBin + 1, fracPct);

        StatUncPerInteraction[signalBin]->SetBinContent(energyBin + 1, statsigma);
        ShapeUncPerInteraction[signalBin]->SetBinContent(energyBin + 1, shapesigma);
    }
    for (int i = 0; i < NUM_SIGNAL_TYPES; ++i) {
        reweightOneDHisto(UnfPerInteraction[i], 100.);
        reweightOneDHisto(TotalUncPerInteraction[i], 100.);
        reweightOneDHisto(StatUncPerInteraction[i], 100.);
        reweightOneDHisto(ShapeUncPerInteraction[i], 100.);
    }

    // Total-from-sources in quadrature (per bin) (store sigma, variance, and %)
    std::vector<double> varFromSources(N, 0.0);

    // Get uncertainty per source
    for (int k = 0; k < (int) UnfCovSources.size(); ++k) {
        TH1D* hUnc = new TH1D(
            Form("hUnc%s", CovarianceLabels[k].c_str()),
            Form("Uncertainty %s; bin;#sigma/x [%%]", CovarianceLabels[k].c_str()),
            N, 0, N
        );
        TH1D* hFracUnc = new TH1D(
            Form("hFracUnc%s", CovarianceLabels[k].c_str()),
            Form("Uncertainty %s; bin;#sigma/x [%%]", CovarianceLabels[k].c_str()),
            N, 0, N
        );

        TH1D* hUncAbs0p_k = new TH1D(
            Form("hUncAbs0p_%s", CovarianceLabels[k].c_str()),
            Form("Abs0p %s;Kinetic Energy [MeV];#sigma [barn] per 100 MeV", CovarianceLabels[k].c_str()),
            NUM_BINS_KE, ARRAY_KE_BINS.data()
        );
        TH1D* hFracAbs0p_k = new TH1D(
            Form("hFracUncAbs0p_%s", CovarianceLabels[k].c_str()),
            Form("Abs0p %s;Kinetic Energy [MeV];#sigma/x [%%]", CovarianceLabels[k].c_str()),
            NUM_BINS_KE, ARRAY_KE_BINS.data()
        );

        TH1D* hUncAbsNp_k = new TH1D(
            Form("hUncAbsNp_%s", CovarianceLabels[k].c_str()),
            Form("AbsNp %s;Kinetic Energy [MeV];#sigma [barn] per 100 MeV", CovarianceLabels[k].c_str()),
            NUM_BINS_KE, ARRAY_KE_BINS.data()
        );
        TH1D* hFracAbsNp_k = new TH1D(
            Form("hFracUncAbsNp_%s", CovarianceLabels[k].c_str()),
            Form("AbsNp %s;Kinetic Energy [MeV];#sigma/x [%%]", CovarianceLabels[k].c_str()),
            NUM_BINS_KE, ARRAY_KE_BINS.data()
        );

        TH1D* hUncScatter_k = new TH1D(
            Form("hUncScatter_%s", CovarianceLabels[k].c_str()),
            Form("Scatter %s;Kinetic Energy [MeV];#sigma [barn] per 100 MeV", CovarianceLabels[k].c_str()),
            NUM_BINS_KE, ARRAY_KE_BINS.data()
        );
        TH1D* hFracScatter_k = new TH1D(
            Form("hFracUncScatter_%s", CovarianceLabels[k].c_str()),
            Form("Scatter %s;Kinetic Energy [MeV];#sigma/x [%%]", CovarianceLabels[k].c_str()),
            NUM_BINS_KE, ARRAY_KE_BINS.data()
        );

        for (int i = 0; i < N; ++i) {
            auto [signalBin, energyBin] = unflattenIndex(i, NUM_BINS_KE);

            const double x   = Unfolded(i);
            const double var = UnfCovSources[k](i,i);

            const double sigma   = std::sqrt(std::max(0.0, var));
            const double fracPct = (std::fabs(x) > 0.0) ? (100.0 * sigma / std::fabs(x)) : 0.0;

            hUnc->SetBinContent(i+1, sigma);
            hFracUnc->SetBinContent(i+1, fracPct);
            varFromSources[i] += std::max(0.0, var);

            if (signalBin == 0) {
                hUncAbs0p_k->SetBinContent(energyBin + 1, sigma);
                hFracAbs0p_k->SetBinContent(energyBin + 1, fracPct);
            } else if (signalBin == 1) {
                hUncAbsNp_k->SetBinContent(energyBin + 1, sigma);
                hFracAbsNp_k->SetBinContent(energyBin + 1, fracPct);
            } else if (signalBin == 2) {
                hUncScatter_k->SetBinContent(energyBin + 1, sigma);
                hFracScatter_k->SetBinContent(energyBin + 1, fracPct);
            }
        }
        reweightOneDHisto(hUncAbs0p_k,   100.);
        reweightOneDHisto(hUncAbsNp_k,   100.);
        reweightOneDHisto(hUncScatter_k, 100.);

        // store totals
        UncTotal.push_back(hUnc);
        FracUncTotal.push_back(hFracUnc);

        // store per-interaction (source-by-source)
        SourceUncPerInteraction[0].push_back(hUncAbs0p_k);
        SourceFracUncPerInteraction[0].push_back(hFracAbs0p_k);

        SourceUncPerInteraction[1].push_back(hUncAbsNp_k);
        SourceFracUncPerInteraction[1].push_back(hFracAbsNp_k);

        SourceUncPerInteraction[2].push_back(hUncScatter_k);
        SourceFracUncPerInteraction[2].push_back(hFracScatter_k);

    }

    // Sanity check
    double maxAbsDiffPct = 0.0;
    int worstBin         = -1;
    for (int i = 0; i < N; ++i) {
        const double x = Unfolded(i);

        // total sigma from sources
        const double sigFromSources     = std::sqrt(std::max(0.0, varFromSources[i]));
        const double fracFromSourcesPct = (std::fabs(x) > 0.0) ? (100.0 * sigFromSources / std::fabs(x)) : 0.0;

        // total sigma from total unfolded covariance
        const double sigFromTotalCov     = std::sqrt(std::max(0.0, UnfoldCov(i,i)));
        const double fracFromTotalCovPct = (std::fabs(x) > 0.0) ? (100.0 * sigFromTotalCov / std::fabs(x)) : 0.0;

        const double absDiff = std::fabs(fracFromSourcesPct - fracFromTotalCovPct);
        if (absDiff > maxAbsDiffPct) { maxAbsDiffPct = absDiff; worstBin = i; }
    }

    //////////////////
    // Create plots //
    //////////////////
    
    std::vector<int> Colors = {
        kBlack,
        kBlue,
        kRed,
        kGreen,
        kOrange+1,
        kMagenta,
        kCyan+1,
        kViolet+1,
        kAzure+1,
        kPink+6
    };

    std::vector<std::vector<TH1*>> PlotGroups;
    std::vector<std::vector<TString>> PlotLabelGroups;
    std::vector<TString> PlotTitles;
    std::vector<TString> XLabels;
    std::vector<TString> YLabels;
    std::vector<bool> PlotStacked;
    std::vector<std::vector<bool>> PlotsAsPoints;
    std::vector<std::vector<bool>> PlotNumbers;

    // Add total histograms
    std::vector<TH1*> AllUnc;
    std::vector<TH1*> AllFracUnc;
    std::vector<TString> UncLabels;
    std::vector<bool> Points;
    std::vector<bool> Numbers;

    AllUnc.push_back(hUncTotal);
    AllFracUnc.push_back(hFracUncTotal);
    UncLabels.push_back("Total");
    Points.push_back(false);
    Numbers.push_back(true);
    for (int i = 0; i < CovarianceLabels.size(); ++i) {
        AllUnc.push_back(UncTotal[i]);
        AllFracUnc.push_back(FracUncTotal[i]);
        UncLabels.push_back(CovarianceLabels[i]);
        Points.push_back(false);
        Numbers.push_back(true);
    }

    PlotGroups.push_back(AllUnc);
    PlotLabelGroups.push_back(UncLabels);
    PlotTitles.push_back("UncertaintyBreakdown/TotalUnc");
    XLabels.push_back("Bin (i, #alpha)");
    YLabels.push_back("Cross section unc. [barn] per 100 MeV");
    PlotStacked.push_back(false);
    PlotsAsPoints.push_back(Points);
    PlotNumbers.push_back(Numbers);

    PlotGroups.push_back(AllFracUnc);
    PlotLabelGroups.push_back(UncLabels);
    PlotTitles.push_back("UncertaintyBreakdown/TotalFracUnc");
    XLabels.push_back("Bin (i, #alpha)");
    YLabels.push_back("Cross section unc. [%]");
    PlotStacked.push_back(false);
    PlotsAsPoints.push_back(Points);
    PlotNumbers.push_back(Numbers);

    // For each signal type, get uncertainty breakdown
    for (int iSignal = 0; iSignal < NUM_SIGNAL_TYPES; ++iSignal) {
        std::vector<TH1*> UncPlots;
        std::vector<TH1*> FracUncPlots;
        std::vector<TString> LabelPlots;
        std::vector<bool> PointsPlots;
        std::vector<bool> NumbersPlots;

        UncPlots.push_back(TotalUncPerInteraction[iSignal]);
        FracUncPlots.push_back(TotalFracUncPerInteraction[iSignal]);
        LabelPlots.push_back("Total");
        PointsPlots.push_back(false);
        NumbersPlots.push_back(true);
        for (int iSource = 0; iSource < CovarianceLabels.size(); ++iSource) {
            UncPlots.push_back(SourceUncPerInteraction[iSignal][iSource]);
            FracUncPlots.push_back(SourceFracUncPerInteraction[iSignal][iSource]);
            LabelPlots.push_back(CovarianceLabels[iSource]);
            PointsPlots.push_back(false);
            NumbersPlots.push_back(true);
        }

        PlotGroups.push_back(UncPlots);
        PlotLabelGroups.push_back(LabelPlots);
        if      (iSignal == 0) PlotTitles.push_back("UncertaintyBreakdown/Abs0pUnc");
        else if (iSignal == 1) PlotTitles.push_back("UncertaintyBreakdown/AbsNpUnc");
        else if (iSignal == 2) PlotTitles.push_back("UncertaintyBreakdown/ScatterUnc");
        XLabels.push_back("Kinetic Energy [MeV]");
        YLabels.push_back("Cross section unc. [barn] per 100 MeV");
        PlotStacked.push_back(false);
        PlotsAsPoints.push_back(PointsPlots);
        PlotNumbers.push_back(NumbersPlots);

        PlotGroups.push_back(FracUncPlots);
        PlotLabelGroups.push_back(LabelPlots);
        if      (iSignal == 0) PlotTitles.push_back("UncertaintyBreakdown/Abs0pFracUnc");
        else if (iSignal == 1) PlotTitles.push_back("UncertaintyBreakdown/AbsNpFracUnc");
        else if (iSignal == 2) PlotTitles.push_back("UncertaintyBreakdown/ScatterFracUnc");
        XLabels.push_back("Kinetic Energy [MeV]");
        YLabels.push_back("Cross section unc. [%]");
        PlotStacked.push_back(false);
        PlotsAsPoints.push_back(PointsPlots);
        PlotNumbers.push_back(NumbersPlots);
    }

    printOneDPlots(
        SaveDir, FontStyle, TextSize,
        PlotGroups,
        Colors,
        PlotLabelGroups,
        PlotTitles,
        XLabels,
        YLabels,
        PlotStacked,
        PlotsAsPoints,
        PlotNumbers
    );

    ///////////////////////////
    // Two-dimensional plots //
    ///////////////////////////

    std::vector<TH2*> TwoDPlots;
    std::vector<TString> TwoDTitles;
    std::vector<std::pair<double,double>> TwoDRanges;
    std::vector<bool> TwoDDisplayNumbers;

    // Add all unfolded covariances
    for (int i = 0; i < UnfCovSources.size(); ++i) {
        TH2D* unfcov = new TH2D(
            (TString) CovarianceLabels[i] + " Unfolded Covariance",
            (TString) CovarianceLabels[i] + " Unfolded Covariance;Reco (j, #beta);True (i, #alpha)",
            N, 0, N,
            N, 0, N
        ); M2H(UnfCovSources[i], unfcov);

        TH2D* unffraccov = new TH2D(
            (TString) CovarianceLabels[i] + " Unfolded Frac. Covariance",
            (TString) CovarianceLabels[i] + " Unfolded Frac. Covariance;Reco (j, #beta);True (i, #alpha)",
            N, 0, N,
            N, 0, N
        );
        TH2D* unfcorr = new TH2D(
            (TString) CovarianceLabels[i] + " Unfolded Correlation",
            (TString) CovarianceLabels[i] + " Unfolded Correlation;Reco (j, #beta);True (i, #alpha)",
            N, 0, N,
            N, 0, N
        );
        GetFracCovAndCorrMatrix(hUnfolded, unfcov, unffraccov, unfcorr);

        TwoDPlots.push_back(unfcov);
        TwoDTitles.push_back(CovarianceLabels[i] + "/UnfoldedCovariance");
        TwoDRanges.push_back({0,0});
        TwoDDisplayNumbers.push_back(false);

        TwoDPlots.push_back(unffraccov);
        TwoDTitles.push_back(CovarianceLabels[i] + "/UnfoldedFracCovariance");
        TwoDRanges.push_back({0,0});
        TwoDDisplayNumbers.push_back(false);

        TwoDPlots.push_back(unfcorr);
        TwoDTitles.push_back(CovarianceLabels[i] + "/UnfoldedCorrelation");
        TwoDRanges.push_back({-1,-1});
        TwoDDisplayNumbers.push_back(false);
    }

    // Get total covariance
    TH2D* totalunfcov = new TH2D(
        "Total Unfolded Covariance",
        "Total Unfolded Covariance;Reco (j, #beta);True (i, #alpha)",
        N, 0, N,
        N, 0, N
    ); M2H(UnfoldCov, totalunfcov);

    TH2D* totalunffraccov = new TH2D(
        "Total Unfolded Frac. Covariance",
        "Total Unfolded Frac. Covariance;Reco (j, #beta);True (i, #alpha)",
        N, 0, N,
        N, 0, N
    );
    TH2D* totalunfcorr = new TH2D(
        "Total Unfolded Correlation",
        "Total Unfolded Correlation;Reco (j, #beta);True (i, #alpha)",
        N, 0, N,
        N, 0, N
    );

    GetFracCovAndCorrMatrix(hUnfolded, totalunfcov, totalunffraccov, totalunfcorr);

    TwoDPlots.push_back(totalunfcov);
    TwoDTitles.push_back("TotalUnfoldedCovariance");
    TwoDRanges.push_back({0,0});
    TwoDDisplayNumbers.push_back(false);

    TwoDPlots.push_back(totalunffraccov);
    TwoDTitles.push_back("TotalUnfoldedFracCovariance");
    TwoDRanges.push_back({0,0});
    TwoDDisplayNumbers.push_back(false);

    TwoDPlots.push_back(totalunfcorr);
    TwoDTitles.push_back("TotalUnfoldedCorrelation");
    TwoDRanges.push_back({-1,-1});
    TwoDDisplayNumbers.push_back(false);

    // Total covariance before unfolding
    TH2D* totalcov = new TH2D(
        "Total Covariance",
        "Total Covariance;Reco (j, #beta);True (i, #alpha)",
        N, 0, N,
        N, 0, N
    ); M2H(TotalCovariance, totalcov);

    TH2D* totalfraccov = new TH2D(
        "Total Frac. Covariance",
        "Total Frac. Covariance;Reco (j, #beta);True (i, #alpha)",
        N, 0, N,
        N, 0, N
    );
    TH2D* totalcorr = new TH2D(
        "Total Correlation",
        "Total Correlation;Reco (j, #beta);True (i, #alpha)",
        N, 0, N,
        N, 0, N
    );

    TH1D* MeasureMinusBkg = new TH1D("MeasureMinusBkg", "MeasureMinusBkg;;", N, 0, N);
    V2H(MeasureMinusBackground, MeasureMinusBkg);
    GetFracCovAndCorrMatrix(MeasureMinusBkg, totalcov, totalfraccov, totalcorr);

    TwoDPlots.push_back(totalcov);
    TwoDTitles.push_back("TotalCovariance");
    TwoDRanges.push_back({0,0});
    TwoDDisplayNumbers.push_back(false);

    TwoDPlots.push_back(totalfraccov);
    TwoDTitles.push_back("TotalFracCovariance");
    TwoDRanges.push_back({0,0});
    TwoDDisplayNumbers.push_back(false);

    TwoDPlots.push_back(totalcorr);
    TwoDTitles.push_back("TotalCorrelation");
    TwoDRanges.push_back({-1,-1});
    TwoDDisplayNumbers.push_back(false);

    printTwoDPlots(SaveDir, TwoDPlots, TwoDTitles, TwoDRanges, TwoDDisplayNumbers);

    ////////////////
    // Nice plots //
    ////////////////

    PrintDataUnfoldedStatShapePlot(
        SaveDir,
        "TotalUnf",
        TrueNominal,
        hUnfolded,
        hUncStat,
        hUncShape,
        "Unfolded Total Vector",
        "Bin number",
        "Cross section [barn]"
    );

    TH1D* hTrueAbs0p   = new TH1D("hTrueAbs0p_plot",   ";Kinetic Energy [MeV];Cross section [barn]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueAbsNp   = new TH1D("hTrueAbsNp_plot",   ";Kinetic Energy [MeV];Cross section [barn]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    TH1D* hTrueScatter = new TH1D("hTrueScatter_plot", ";Kinetic Energy [MeV];Cross section [barn]", NUM_BINS_KE, ARRAY_KE_BINS.data());
    std::vector<TH1D*> TrueKE = {hTrueAbs0p, hTrueAbsNp, hTrueScatter};

    for (int idx = 0; idx < N; ++idx) {
        auto [signalBin, energyBin] = unflattenIndex(idx, NUM_BINS_KE);
        TrueKE[signalBin]->SetBinContent(energyBin + 1, TrueSignal(idx));
    }
    for (int i = 0; i < NUM_SIGNAL_TYPES; ++i) reweightOneDHisto(TrueKE[i],  100.);

    PrintDataUnfoldedStatShapePlot(
        SaveDir,
        "Abs0pUnf",
        TrueKE[0],
        UnfPerInteraction[0],
        StatUncPerInteraction[0],
        ShapeUncPerInteraction[0],
        "Unfolded Absorption 0p",
        "Kinetic Energy [MeV]",
        "Cross section [barn] per 100 MeV"
    );
    PrintDataUnfoldedStatShapePlot(
        SaveDir,
        "AbsNpUnf",
        TrueKE[1],
        UnfPerInteraction[1],
        StatUncPerInteraction[1],
        ShapeUncPerInteraction[1],
        "Unfolded Absorption Np",
        "Kinetic Energy [MeV]",
        "Cross section [barn] per 100 MeV"
    );
    PrintDataUnfoldedStatShapePlot(
        SaveDir,
        "ScatterUnf",
        TrueKE[2],
        UnfPerInteraction[2],
        StatUncPerInteraction[2],
        ShapeUncPerInteraction[2],
        "Unfolded Scatter",
        "Kinetic Energy [MeV]",
        "Cross section [barn] per 100 MeV"
    );

    ////////////////////////////////////
    // Save everything from later use //
    ////////////////////////////////////

    TString outPath = "/exp/lariat/app/users/epelaez/histos/data/Unfold1.root";
    TFile outAll(outPath, "RECREATE");
    outAll.cd();

    std::unordered_set<std::string> written;

    // Helper lambda: write once by name
    auto writeOnce = [&](TObject* obj) {
        if (!obj) return;
        const std::string name = obj->GetName();
        if (name.empty()) return;
        if (written.insert(name).second) {
            obj->Write(name.c_str(), TObject::kOverwrite);
        }
    };

    // 1D groups
    for (auto& group : PlotGroups) {
        for (auto* h : group) writeOnce(h);
    }

    // 1D data plots
    for (auto* h : UnfPerInteraction) writeOnce(h);
    for (auto* h : StatUncPerInteraction) writeOnce(h);
    for (auto* h : ShapeUncPerInteraction) writeOnce(h);
    for (auto* h : TrueKE) writeOnce(h);

    outAll.Write();
    outAll.Close();

    std::cout << "\nWrote " << written.size() << " objects to " << outPath << std::endl;
}