#pragma once
namespace FakeDataFD {
    enum class Scenario {
        Nominal = 0,  // no reweighting (weight = 1)
        FDa,
        FDb,
        FDc
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
                switch (region) {
                    case ERegion::Low:
                        w0p_raw    = 1.6;
                        wNp_raw    = 0.4;
                        wScatt_raw = 0.8;
                        break;
                    case ERegion::Mid:
                        w0p_raw    = 1.2;
                        wNp_raw    = 0.8;
                        wScatt_raw = 1.0;
                        break;
                    case ERegion::High:
                        w0p_raw    = 0.8;
                        wNp_raw    = 1.0;
                        wScatt_raw = 1.4;
                        break;
                    }
                return;

            case Scenario::FDb:
                // FD-1b:
                // Mid:  0p x0.9, Np x0.9, scat x1.4 (Δ-like scattering bump)
                switch (region) {
                    case ERegion::Low:
                        w0p_raw    = 0.8;
                        wNp_raw    = 0.8;
                        wScatt_raw = 0.8;
                        break;
                    case ERegion::Mid:
                        w0p_raw    = 0.8;
                        wNp_raw    = 0.8;
                        wScatt_raw = 1.6;
                        break;
                    case ERegion::High:
                        w0p_raw    = 1.2;
                        wNp_raw    = 1.2;
                        wScatt_raw = 0.8;
                        break;
                    }
                return;

            case Scenario::FDc:
                // Uniform 50% increase
                switch (region) {
                    case ERegion::Low:
                        w0p_raw    = 1.5;
                        wNp_raw    = 1.5;
                        wScatt_raw = 1.5;
                        break;
                    case ERegion::Mid:
                        w0p_raw    = 1.5;
                        wNp_raw    = 1.5;
                        wScatt_raw = 1.5;
                        break;
                    case ERegion::High:
                        w0p_raw    = 1.5;
                        wNp_raw    = 1.5;
                        wScatt_raw = 1.5;
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
