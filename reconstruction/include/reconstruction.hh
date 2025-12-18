#ifndef LOR_HH
#define LOR_HH

#include <iostream>
#include <cmath>
#include "TVector3.h"

// --------------------------
// Event inputs per detector
// --------------------------
struct Hits
{
    double   time;      // ns
    TVector3 pos;       // 1st interaction (scatter) position S
    double   energy;    // MeV (scatter energy used for Compton angle)

    Hits()
    : time(std::sqrt(-1.0))
    , pos(std::sqrt(-1.0), std::sqrt(-1.0), std::sqrt(-1.0))
    , energy(std::sqrt(-1.0))
    {}
};

// --------------------------
// Reconstruction class
// --------------------------
class Reconstruction
{
public:
    Reconstruction();
    ~Reconstruction();

    // Inputs
    void SetFirstHit (const Hits& h) { hit1 = h; }
    void SetSecondHit(const Hits& h) { hit2 = h; }

    // Detector geometry (inward normals must be unit vectors)
    void SetDetectorNormals(const TVector3& n_inward_1, const TVector3& n_inward_2);
    void SetAcceptanceHalfAngles(double psi1_rad, double psi2_rad); // radians

    // Optional tuning
    void SetTOFSigmaMM(double sigma_mm);     // std dev of TOF prior in mm (default ~15 mm)
    void SetSearchBracket(double s_min_mm, double s_max_mm); // default [-500, +500]
    void SetEnergyWindow(double e_min_mev, double e_max_mev); // default [0.4, 0.6]
    void SetTimeGate(double abs_dt_ns_max); // default 0.2 ns

    // Run reconstruction
    // - If normals/acceptances not set (NaN), falls back to classic TOF-LOR.
    void Run();

    // Result
    TVector3 GetResult() const { return annihilationPoint; }

private:
    // Physical constants / gates
    const double me;  // electron mass (MeV)
    double eMin, eMax;         // energy window
    double dtAbsMax;           // |dt| gate (ns)
    double sigma_s_mm;         // TOF prior sigma (mm)
    double sMin, sMax;         // search bracket along LOR (mm)

    // Geometry
    TVector3 n1, n2;           // inward normals (unit) for head 1,2
    double psi1, psi2;         // acceptance half-angles (rad)

    // Event
    Hits hit1, hit2;

    // Output
    TVector3 annihilationPoint;

    // Helpers
    static bool IsNaN(double x) { return x != x; } // NaN != NaN
    static bool IsFiniteVec3(const TVector3& v)
    {
        return !(IsNaN(v.X()) || IsNaN(v.Y()) || IsNaN(v.Z()));
    }
    static double Clamp(double x, double a, double b)
    {
        return (x < a) ? a : ((x > b) ? b : x);
    }

    // Core steps
    void RunBasicTOFLOR(); // your original method
    double CalculateComptonAngle(double energy); // radians

    // Internals for 1D minimization
    struct AcceptanceRelaxedResidual;
};

// --------------------------
// Residual functor declaration
// --------------------------
struct Reconstruction::AcceptanceRelaxedResidual {
    TVector3 S1, S2, n1, n2, mid, u;
    double theta1, theta2;   // radians
    double psi1, psi2;       // radians (acceptance half-angles)
    double dTOF;             // mm (c*dt/2)
    double sigma_s_mm;       // TOF sigma in mm

    AcceptanceRelaxedResidual();

    static double Clamp(double x, double a, double b)
    {
        return (x < a) ? a : ((x > b) ? b : x);
    }

    double operator()(double s);
};

// --------------------------
// Golden-section minimizer
// --------------------------
namespace recon_min {
    template <typename F>
    double GoldenSectionMin(F& f, double a, double b, int iters, double tol)
    {
        const double phi = (std::sqrt(5.0) - 1.0) / 2.0;
        double c = b - phi * (b - a);
        double d = a + phi * (b - a);
        double fc = f(c);
        double fd = f(d);

        int it = 0;
        while (it < iters && std::fabs(b - a) > tol) {
            if (fc < fd) {
                b = d; d = c; fd = fc;
                c = b - phi*(b - a); fc = f(c);
            } else {
                a = c; c = d; fc = fd;
                d = a + phi*(b - a); fd = f(d);
            }
            ++it;
        }
        return 0.5 * (a + b);
    }
}

#endif // LOR_HH
