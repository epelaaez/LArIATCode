#include "TFile.h"
#include "TH1D.h"

#include <iostream>
#include <iomanip>

double FullIntegral(TH1D* h)
{
    // Include underflow and overflow so this measures all
    // selected events, not just those inside the plotted KE range.
    return h->Integral(0, h->GetNbinsX() + 1);
}

void PrintBeamComposition()
{
    TFile* f = TFile::Open(
        "/exp/lariat/app/users/mking/abs4/analysis/LArIATCode/"
        "histos/RecoClassify3Cat_AllHists.root",
        "READ"
    );

    if (!f || f->IsZombie()) {
        std::cerr << "Could not open input ROOT file." << std::endl;
        return;
    }

    TH1D* hAll = (TH1D*)f->Get("hFrontFacePionKE");
    TH1D* hPi  = (TH1D*)f->Get("hFrontFacePionKEPion");
    TH1D* hMu  = (TH1D*)f->Get("hFrontFacePionKEMuon");
    TH1D* hEl  = (TH1D*)f->Get("hFrontFacePionKEElectron");

    if (!hAll || !hPi || !hMu || !hEl) {
        std::cerr << "Could not find one or more beam-composition histograms."
                  << std::endl;
        return;
    }

    const double nAll = FullIntegral(hAll);
    const double nPi  = FullIntegral(hPi);
    const double nMu  = FullIntegral(hMu);
    const double nEl  = FullIntegral(hEl);

    const double nPiMu = nPi + nMu;

    std::cout << std::fixed << std::setprecision(4);

    std::cout << "\n========================================\n";
    std::cout << "Selected beam composition\n";
    std::cout << "After TG + small-track cuts\n";
    std::cout << "Includes front-face KE reweighting\n";
    std::cout << "========================================\n";

    std::cout << "\nWeighted event counts:\n";
    std::cout << "  Total    = " << nAll << "\n";
    std::cout << "  Pions    = " << nPi  << "\n";
    std::cout << "  Muons    = " << nMu  << "\n";
    std::cout << "  Electrons= " << nEl  << "\n";

    std::cout << "\nFractions of full selected beam:\n";
    std::cout << "  pion     = " << 100.0 * nPi / nAll << " %\n";
    std::cout << "  muon     = " << 100.0 * nMu / nAll << " %\n";
    std::cout << "  electron = " << 100.0 * nEl / nAll << " %\n";

    std::cout << "\nComposition within pion + muon component:\n";
    std::cout << "  pion     = " << 100.0 * nPi / nPiMu << " %\n";
    std::cout << "  muon     = " << 100.0 * nMu / nPiMu << " %\n";

    std::cout << "\nClosure:\n";
    std::cout << "  Pi + Mu + El = " << nPi + nMu + nEl << "\n";
    std::cout << "  All          = " << nAll << "\n";
    std::cout << "  Difference   = "
              << nAll - (nPi + nMu + nEl) << "\n";

    f->Close();
}
