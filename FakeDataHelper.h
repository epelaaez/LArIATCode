#pragma once
namespace FakeDataFD {
    enum class Scenario {
        Nominal = 0,  // no reweighting (weight = 1)
        FDa,
        FDb,
        FDc,
        FDd
    };

    inline bool IsAbs0p(int id) {
        return (id == 0); // "Abs 0p"
    }

    inline bool IsAbsNp(int id) {
        return (id == 1); // "Abs Np"
    }

    inline bool IsScattering(int id) {
    switch (id) {
        case 6:  // Inelastic scattering
        case 12: // Elastic scattering
            return true;
        default:
            return false;
        }
    }

    ////////////////////
    // Energy regions //
    ////////////////////

    enum class ERegion {
        Low,
        Mid,
        High
    };

    inline ERegion GetEnergyRegion(double trueKE_GeV) {
        constexpr double kLowMax = 200;  // MeV
        constexpr double kMidMax = 400;  // MeV

        if (trueKE_GeV < kLowMax) {
            return ERegion::Low;
        } else if (trueKE_GeV < kMidMax) {
            return ERegion::Mid;
        } else {
            return ERegion::High;
        }
    }

    ////////////////
    // Get weight //
    ////////////////

    inline void GetRawChannelWeights(
        Scenario scenario,
        ERegion region,
        double &w0p_raw,
        double &wNp_raw,
        double &wScatt_raw
    ) {
        // Default: no change
        w0p_raw    = 1.0;
        wNp_raw    = 1.0;
        wScatt_raw = 1.0;

        switch (scenario) {
            case Scenario::Nominal:
                // all = 1
                return;

            case Scenario::FDa:
                // FD-1a:
                // Low:  0p x1.3, Np x0.7, scat x1.0
                // Mid:  unchanged
                // High: 0p x0.8, Np x1.3, scat x1.0
                switch (region) {
                    case ERegion::Low:
                        w0p_raw    = 1.3;
                        wNp_raw    = 0.7;
                        wScatt_raw = 1.0;
                        break;
                    case ERegion::Mid:
                        w0p_raw    = 1.0;
                        wNp_raw    = 1.0;
                        wScatt_raw = 1.0;
                        break;
                    case ERegion::High:
                        w0p_raw    = 0.8;
                        wNp_raw    = 1.3;
                        wScatt_raw = 1.0;
                        break;
                    }
                return;

            case Scenario::FDb:
                // FD-1b:
                // Low:  0p x1.4, Np x0.6, scat x1.0
                // Mid:  0p x0.9, Np x1.1, scat x1.0
                // High: 0p x0.7, Np x1.4, scat x1.0
                switch (region) {
                    case ERegion::Low:
                        w0p_raw    = 1.4;
                        wNp_raw    = 0.6;
                        wScatt_raw = 1.0;
                        break;
                    case ERegion::Mid:
                        w0p_raw    = 0.9;
                        wNp_raw    = 1.1;
                        wScatt_raw = 1.0;
                        break;
                    case ERegion::High:
                        w0p_raw    = 0.7;
                        wNp_raw    = 1.4;
                        wScatt_raw = 1.0;
                        break;
                    }
                return;

            case Scenario::FDc:
                // FD-1c:
                // Low:  unchanged
                // Mid:  0p x0.9, Np x0.9, scat x1.4 (Δ-like scattering bump)
                // High: 0p x1.1, Np x1.1, scat x0.8
                switch (region) {
                    case ERegion::Low:
                        w0p_raw    = 1.0;
                        wNp_raw    = 1.0;
                        wScatt_raw = 1.0;
                        break;
                    case ERegion::Mid:
                        w0p_raw    = 0.9;
                        wNp_raw    = 0.9;
                        wScatt_raw = 1.4;
                        break;
                    case ERegion::High:
                        w0p_raw    = 1.1;
                        wNp_raw    = 1.1;
                        wScatt_raw = 0.8;
                        break;
                    }
                return;

            case Scenario::FDd:
                // FD-1d:
                // Low:  0p x1.1, Np x0.6, scat x1.2
                // Mid:  0p x1.0, Np x0.7, scat x1.2
                // High: 0p x0.9, Np x0.6, scat x1.3
                switch (region) {
                    case ERegion::Low:
                        w0p_raw    = 1.1;
                        wNp_raw    = 0.6;
                        wScatt_raw = 1.2;
                        break;
                    case ERegion::Mid:
                        w0p_raw    = 1.0;
                        wNp_raw    = 0.7;
                        wScatt_raw = 1.2;
                        break;
                    case ERegion::High:
                        w0p_raw    = 0.9;
                        wNp_raw    = 0.6;
                        wScatt_raw = 1.3;
                        break;
                }
                return;
        }
        throw std::runtime_error("Unknown FakeData::Scenario");
    }



    inline double GetFDWeight(
        Scenario scenario,
        int id,
        double trueKE
    ) {
        // Non-signal categories: leave unchanged.
        const bool is0p    = IsAbs0p(id);
        const bool isNp    = IsAbsNp(id);
        const bool isScatt = IsScattering(id);

        if (!is0p && !isNp && !isScatt) {
            return 1.0;
        }

        const ERegion region = GetEnergyRegion(trueKE);

        double w0p_raw    = 1.0;
        double wNp_raw    = 1.0;
        double wScatt_raw = 1.0;

        GetRawChannelWeights(scenario, region, w0p_raw, wNp_raw, wScatt_raw);

        if (is0p)    return w0p_raw;
        if (isNp)    return wNp_raw;
        if (isScatt) return wScatt_raw;

        return 1.0;
    }

} // namespace FakeDataFD
