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

void UnfoldFD() {
    // Set defaults
    gStyle->SetOptStat(0); // get rid of stats box
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    TH1::AddDirectory(false);
    gStyle->SetPalette(kRainBow);
    gStyle->SetPaintTextFormat("4.2f");

    int FontStyle = 132;
    double TextSize = 0.06;
    TString SaveDir = "/exp/lariat/app/users/epelaez/analysis/figs/FakeData/";

    int N = NUM_SIGNAL_TYPES * NUM_BINS_KE;

    // Load nominal response matrix and true
    std::unique_ptr<TFile> NominalFile(TFile::Open("/exp/lariat/app/users/epelaez/histos/nominal/Reco.root"));
    TH1D* TrueNominal     = (TH1D*) NominalFile->Get("hSignalVectorNominal");
    TH2D* ResponseNominal = (TH2D*) NominalFile->Get("hResponseMatrix");
    TH1D* BackgroundNominal = (TH1D*) NominalFile->Get("hBackgroundVectorNominal");

    // Load fake data samples
    std::unique_ptr<TFile> FD1File(TFile::Open("/exp/lariat/app/users/epelaez/histos/fake_data/FD1.root"));
    TH1D* FD1True = (TH1D*) FD1File->Get("hSignal");
    TH1D* FD1Reco = (TH1D*) FD1File->Get("hMeasure");

    std::vector<TH1D*> FDTrue = {
        FD1True
    };

    std::vector<TH1D*> FDReco = {
        FD1Reco
    };

    std::vector<std::string> FDLabel = {
        "Sample1"
    };

    // MC stat covariance
    TH2D* StatCovariance = (TH2D*) NominalFile->Get("hStatCovariance");

    // Generator covariance
    std::unique_ptr<TFile> GeneratorFile(TFile::Open("/exp/lariat/app/users/epelaez/histos/generator/Measure.root"));
    TH2D* GeneratorCovariance = (TH2D*) GeneratorFile->Get("hMeasureCovMatrix");

    // Energy reconstruction covariance
    std::unique_ptr<TFile> EnergyRecoFile(TFile::Open("/exp/lariat/app/users/epelaez/histos/energyreco/Measure.root"));
    TH2D* EnergyRecoCovariance = (TH2D*) EnergyRecoFile->Get("hMeasureCovMatrix");

    // Beam line muon covariance
    std::unique_ptr<TFile> BeamlineMuFile(TFile::Open("/exp/lariat/app/users/epelaez/histos/beamlinemu/Measure.root"));
    TH2D* BeamlineMuCovariance = (TH2D*) BeamlineMuFile->Get("hMeasureCovMatrix");

    // Beam line electron covariance
    std::unique_ptr<TFile> BeamlineElFile(TFile::Open("/exp/lariat/app/users/epelaez/histos/beamlineel/Measure.root"));
    TH2D* BeamlineElCovariance = (TH2D*) BeamlineElFile->Get("hMeasureCovMatrix");

    // Convert histograms into matrices/vectors
    TVectorD Background(N); H2V(BackgroundNominal, Background);
    TMatrixD Response(N, N); H2M(ResponseNominal, Response, kTRUE);

    std::vector<TH2D*> CovarianceMatrices = {
        StatCovariance, 
        GeneratorCovariance,
        EnergyRecoCovariance,
        BeamlineMuCovariance,
        BeamlineElCovariance
    };

    std::vector<std::string> CovarianceLabels = {
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

    // Unfold fake data
    for (int iFDUniv = 0; iFDUniv < FDReco.size(); ++iFDUniv) {
        TVectorD Measure(N); H2V(FDReco[iFDUniv], Measure);
        TVectorD True(N); H2V(FDTrue[iFDUniv], True);

        // Subtract background
        TVectorD MeasureMinusBackground = Measure - Background;

        // Vectors to store output
        TMatrixD AddSmear(N, N); TMatrixD AddSmearInverse(N, N);
        TVectorD WF(N); TMatrixD UnfoldCov(N, N); TMatrixD CovRotation(N, N);

        // Unfold
        TVectorD Unfolded = WienerSVD(
            Response,
            True,
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
        TH1D* hUnfoldedUnc = new TH1D("hUnfoldedUnc", "hUnfoldedUnc;;", N, 0, N);
        TH1D* hUnfolded    = new TH1D("hUnfolded", "hUnfolded;;", N, 0, N); V2H(Unfolded, hUnfolded);

        TH2D* hUnfCovariance = new TH2D(
            "Unfolded Covariance", "Unfolded Covariance Matrix;(i, #alpha);(j, #beta)",
            NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE, 
            NUM_SIGNAL_TYPES * NUM_BINS_KE, 0, NUM_SIGNAL_TYPES * NUM_BINS_KE 
        ); M2H(UnfoldCov, hUnfCovariance);

        // Populate output
        for (int i = 0; i < N; ++i) {
            auto [signalBin, energyBin] = unflattenIndex(i, NUM_BINS_KE);
            const double var   = UnfoldCov(i,i);
            const double sigma = std::sqrt(std::max(0.0, var));
            hUnfoldedUnc->SetBinContent(i, sigma);
        }

        // Get chi2, ndof, pval, sigma
        double f_chi; int f_ndof; double f_pval; double f_sigma;
        CalcChiSquared(FDTrue[iFDUniv], hUnfolded, hUnfCovariance, f_chi, f_ndof, f_pval, f_sigma);

        double n_chi; int n_ndof; double n_pval; double n_sigma;
        CalcChiSquared(TrueNominal, hUnfolded, hUnfCovariance, n_chi, n_ndof, n_pval, n_sigma);

        // TODO: temporal stat temp
        PrintFDPlot(
            SaveDir + FDLabel[iFDUniv] + "/",
            "Unfolded",
            FDTrue[iFDUniv],
            TrueNominal,
            hUnfolded,
            hUnfoldedUnc,
            "Unfolded Total Vector",
            "Kinetic Energy [MeV]",
            "Cross section [barn]",
            {f_chi, n_chi},
            {f_ndof, n_ndof},
            {f_pval, n_pval},
            {f_sigma, n_sigma}
        );
    }
}
