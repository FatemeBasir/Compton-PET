#include "reconstruction.hh"

// --------------------------
// Constructor / config
// --------------------------
Reconstruction::Reconstruction()
: me(0.511)
, eMin(0.4), eMax(0.6)
, dtAbsMax(0.2)            // 200 ps
, sigma_s_mm(15.0)         // ~100 ps timing sigma → ~15 mm
, sMin(-500.0), sMax(500.0)
, n1(std::sqrt(-1.0), std::sqrt(-1.0), std::sqrt(-1.0))
, n2(std::sqrt(-1.0), std::sqrt(-1.0), std::sqrt(-1.0))
, psi1(std::sqrt(-1.0))
, psi2(std::sqrt(-1.0))
{
    annihilationPoint = TVector3(std::sqrt(-1.0), std::sqrt(-1.0), std::sqrt(-1.0));
}

Reconstruction::~Reconstruction() {}

void Reconstruction::SetDetectorNormals(const TVector3& n_inward_1, const TVector3& n_inward_2)
{
    n1 = n_inward_1.Unit();
    n2 = n_inward_2.Unit();
}

void Reconstruction::SetAcceptanceHalfAngles(double psi1_rad, double psi2_rad)
{
    psi1 = psi1_rad;
    psi2 = psi2_rad;
}

void Reconstruction::SetTOFSigmaMM(double sigma_mm) { sigma_s_mm = sigma_mm; }
void Reconstruction::SetSearchBracket(double s_min_mm, double s_max_mm) { sMin = s_min_mm; sMax = s_max_mm; }
void Reconstruction::SetEnergyWindow(double e_min_mev, double e_max_mev) { eMin = e_min_mev; eMax = e_max_mev; }
void Reconstruction::SetTimeGate(double abs_dt_ns_max) { dtAbsMax = abs_dt_ns_max; }

// --------------------------
// Public entry point
// --------------------------
void Reconstruction::Run()
{
    // If normals/acceptances are not set, fall back to classic TOF-LOR.
    if (!IsFiniteVec3(n1) || !IsFiniteVec3(n2) || IsNaN(psi1) || IsNaN(psi2)) {
        RunBasicTOFLOR();
        return;
    }

    // Gates
    if (hit1.energy < eMin || hit1.energy > eMax ||
        hit2.energy < eMin || hit2.energy > eMax) {
        return;
    }

    const double dt = hit1.time - hit2.time;  // ns
    if (std::fabs(dt) > dtAbsMax) return;

    // Compton angles from scatter energy (radians)
    const double theta1 = CalculateComptonAngle(hit1.energy);
    const double theta2 = CalculateComptonAngle(hit2.energy);

    // LOR geometry
    const TVector3 S1  = hit1.pos;
    const TVector3 S2  = hit2.pos;
    if (!IsFiniteVec3(S1) || !IsFiniteVec3(S2)) return;

    const TVector3 u   = (S2 - S1).Unit();
    const TVector3 mid = 0.5 * (S1 + S2);

    // TOF shift
    const double c_mm_per_ns = 299.792458;         // mm/ns
    const double dTOF = 0.5 * c_mm_per_ns * dt;    // mm

    // Residual functor
    AcceptanceRelaxedResidual R;
    R.S1 = S1; R.S2 = S2; R.n1 = n1; R.n2 = n2; R.mid = mid; R.u = u;
    R.theta1 = theta1; R.theta2 = theta2;
    R.psi1 = psi1; R.psi2 = psi2;
    R.dTOF = dTOF;
    R.sigma_s_mm = sigma_s_mm;

    // 1D minimize
    const double sOpt = recon_min::GoldenSectionMin(R, sMin, sMax, 80, 1e-4);
    annihilationPoint = mid + u * sOpt;
}

// --------------------------
// Classic TOF-LOR (fallback)
// --------------------------
void Reconstruction::RunBasicTOFLOR()
{
    if (hit1.energy < eMin || hit1.energy > eMax ||
        hit2.energy < eMin || hit2.energy > eMax) {
        return;
    }

    const double dt = hit1.time - hit2.time;  // ns
    if (std::fabs(dt) > dtAbsMax) return;

    const double dist = 0.5 * (299.792458 * dt);   // mm

    const TVector3 midPoint = 0.5 * (hit2.pos + hit1.pos);
    const TVector3 LOR      = (hit2.pos - hit1.pos).Unit();

    // Optional diagnostics:
    // std::cout << "dt(ns): " << dt << "  d(mm): " << dist << std::endl;

    annihilationPoint = midPoint + LOR * dist;
}

// --------------------------
// Compton angle from energy
// --------------------------
double Reconstruction::CalculateComptonAngle(double energy)
{
    // Your original simplified relation: cos(theta) = 2 - me/E
    double cosTheta = 2.0 - me / energy;
    cosTheta = Clamp(cosTheta, -1.0, 1.0);
    return std::acos(cosTheta);
}

// --------------------------
// Residual functor impl.
// --------------------------
Reconstruction::AcceptanceRelaxedResidual::AcceptanceRelaxedResidual()
: theta1(0.0), theta2(0.0), psi1(0.35), psi2(0.35), // defaults ~20 degrees
  dTOF(0.0), sigma_s_mm(15.0)
{}

double Reconstruction::AcceptanceRelaxedResidual::operator()(double s)
{
    const TVector3 X = mid + u * s;

    TVector3 v1 = S1 - X; if (v1.Mag2()==0.0) return 1e30;
    TVector3 v2 = S2 - X; if (v2.Mag2()==0.0) return 1e30;

    TVector3 v1u = v1.Unit();
    TVector3 v2u = v2.Unit();

    // angles between incoming direction and inward normals
    const double cA1 = Clamp(v1u.Dot(n1), -1.0, 1.0);
    const double cA2 = Clamp(v2u.Dot(n2), -1.0, 1.0);
    const double alpha1 = std::acos(cA1);
    const double alpha2 = std::acos(cA2);

    // acceptance-relaxed mismatch to Compton angles
    double gap1 = std::fabs(alpha1 - theta1) - psi1;
    double gap2 = std::fabs(alpha2 - theta2) - psi2;
    if (gap1 < 0.0) gap1 = 0.0;
    if (gap2 < 0.0) gap2 = 0.0;

    // TOF prior
    const double lambdaTOF = 1.0 / (sigma_s_mm * sigma_s_mm);
    const double rTOF = s - dTOF;

    return gap1*gap1 + gap2*gap2 + lambdaTOF * rTOF * rTOF;
}
