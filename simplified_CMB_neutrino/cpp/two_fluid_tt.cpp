#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double C_LIGHT_KM_S = 299792.458;
constexpr double FID_AS = 2.1e-9;
constexpr double FID_NS = 0.965;
constexpr double ELL_PIVOT = 750.0;
constexpr double T_CMB_MICROK = 2.7255e6;

struct Params {
  double As = FID_AS;
  double omega_cdm_h2 = 0.1201;
  double omega_b_h2 = 0.0223;
  double h = 0.67;
  double ns = FID_NS;
  double N_eff = 3.0;
  int n_source = 1001;
  int ell_min = -1;
  int ell_max = -1;
  int ell_step = 1;
  std::string ell_grid = "class";
  bool interpolate_output = true;
  std::string output = "two_fluid_tt_cpp.dat";
  std::string bessel_cache;
  std::string bessel_x_cache;
  double bessel_x_dx = 0.10;
  double bessel_x_max = 15000.0;
  bool bessel_class_memory = false;
  double hyper_sampling_flat = 8.0;
  double k_samples_per_period = 40.0;
  bool server = false;
  bool scalar_enabled = false;
  double scalar_m_over_H0 = 0.0;
  double scalar_f = 0.0;
  double scalar_theta_ini = 0.0;
  double scalar_theta_prime_ini = 0.0;
  int scalar_n = 1;
  int background_steps = 20000;
  int n_isw = 320;
  double isw_k_eta0_max = 900.0;
  std::string source_mode = "full";
  std::string recfast_output;
  int recfast_steps = 50000;
  double helium_y_p = 0.24;
};

struct Background {
  explicit Background(const Params &p)
      : h0(p.h), H0(100.0 * p.h / C_LIGHT_KM_S),
        omega_gamma(2.47e-5 / (p.h * p.h)),
        omega_cdm(p.omega_cdm_h2 / (p.h * p.h)),
        omega_b(p.omega_b_h2 / (p.h * p.h)) {
    omega_nu = 0.875 * std::pow(4.0 / 11.0, 4.0 / 3.0) * p.N_eff * omega_gamma;
    omega_rad = omega_gamma + omega_nu;
    omega_m = omega_cdm + omega_b;
    omega_lambda = 1.0 - (omega_cdm + omega_b + omega_rad);
  }

  double h0, H0;
  double omega_gamma, omega_nu, omega_rad;
  double omega_cdm, omega_b, omega_m, omega_lambda;
};

// >>>>>>>>>>>>>>>>>>>>>>>>> RECFAST RECOMBINATION <<<<<<<<<<<<<<<<<<<<<<<<<
// This compact background module follows the hydrogen part of RECFAST 1.5.2
// (Scott/Seager/Sasselov).  It computes x_e(z) relative to hydrogen nuclei.
// The full public RECFAST source also evolves helium with singlet/triplet
// corrections; those refinements are intentionally left outside this first
// pedagogical integration point.
struct RecfastHistory {
  std::vector<double> z;
  std::vector<double> xe;
  std::vector<double> tm;
};

struct RecombinationScales {
  double z_star = 0.0;
  double eta_star = 0.0;
  double visibility_sigma_eta = 0.0;
  double silk_length_sq = 0.0;
};

RecfastHistory compute_recfast_history(const Params &p, const Background &bg) {
  constexpr double C = 2.99792458e8;
  constexpr double k_B = 1.380658e-23;
  constexpr double h_P = 6.6260755e-34;
  constexpr double m_e = 9.1093897e-31;
  constexpr double m_H = 1.673575e-27;
  constexpr double sigma_T = 6.6524616e-29;
  constexpr double a_rad = 7.565914e-16;
  constexpr double G = 6.6742e-11;
  constexpr double Lambda_H = 8.2245809;
  constexpr double L_H_ion = 1.096787737e7;
  constexpr double L_H_alpha = 8.225916453e6;
  constexpr double Tnow = 2.7255;
  constexpr double Mpc = 3.0856775814913673e22;

  const double H0_si = 100.0 * p.h * 1000.0 / Mpc;
  const double mu_H = 1.0 / (1.0 - p.helium_y_p);
  const double nH0 = 3.0 * H0_si * H0_si * bg.omega_b /
                     (8.0 * PI * G * mu_H * m_H);
  const double Lalpha = 1.0 / L_H_alpha;
  const double DeltaB = h_P * C * (L_H_ion - L_H_alpha);
  const double CDB = DeltaB / k_B;
  const double CB1 = h_P * C * L_H_ion / k_B;
  const double CR = 2.0 * PI * (m_e / h_P) * (k_B / h_P);
  const double CK = Lalpha * Lalpha * Lalpha / (8.0 * PI);
  const double CL = C * h_P / (k_B * Lalpha);
  const double CT = (8.0 / 3.0) * (sigma_T / (m_e * C)) * a_rad;
  const double fu = 1.14;
  const double H_frac = 1.0e-3;
  const double fHe = p.helium_y_p / (3.9715 * (1.0 - p.helium_y_p));

  auto H_si = [&](double z) {
    const double zp1 = 1.0 + z;
    return H0_si * std::sqrt(bg.omega_rad * std::pow(zp1, 4) +
                             bg.omega_m * std::pow(zp1, 3) +
                             bg.omega_lambda);
  };
  auto dHdz_si = [&](double z) {
    const double zp1 = 1.0 + z;
    const double Hz = H_si(z);
    return (H0_si * H0_si / (2.0 * Hz)) *
           (4.0 * bg.omega_rad * std::pow(zp1, 3) +
            3.0 * bg.omega_m * zp1 * zp1);
  };
  auto saha_xe = [&](double z) {
    const double rhs =
        std::exp(1.5 * std::log(CR * Tnow / (1.0 + z)) -
                 CB1 / (Tnow * (1.0 + z))) /
        nH0;
    return 0.5 * (std::sqrt(rhs * rhs + 4.0 * rhs) - rhs);
  };

  auto deriv = [&](double z, const std::array<double, 2> &y) {
    const double xH = std::clamp(y[0], 1.0e-12, 1.0);
    const double Tmat = std::max(1.0, y[1]);
    const double Trad = Tnow * (1.0 + z);
    const double n = nH0 * std::pow(1.0 + z, 3);
    const double Hz = H_si(z);
    const double Rdown = 1.0e-19 * 4.309 * std::pow(Tmat / 1.0e4, -0.6166) /
                         (1.0 + 0.6703 * std::pow(Tmat / 1.0e4, 0.5300));
    const double Rup =
        Rdown * std::pow(CR * Tmat, 1.5) * std::exp(-CDB / Tmat);
    const double K = CK / Hz;

    std::array<double, 2> dy{};
    dy[0] = ((xH * xH * n * Rdown -
              Rup * (1.0 - xH) * std::exp(-CL / Tmat)) *
             (1.0 + K * Lambda_H * n * (1.0 - xH))) /
            (Hz * (1.0 + z) *
             (1.0 / fu + K * Lambda_H * n * (1.0 - xH) / fu +
              K * Rup * n * (1.0 - xH)));

    const double timeTh = (1.0 / (CT * std::pow(Trad, 4))) *
                          (1.0 + xH + fHe) / xH;
    const double timeH = 2.0 / (3.0 * H0_si * std::pow(1.0 + z, 1.5));
    if (timeTh < H_frac * timeH) {
      const double epsilon = Hz * (1.0 + xH + fHe) /
                             (CT * std::pow(Trad, 3) * xH);
      dy[1] = Tnow +
              epsilon * ((1.0 + fHe) / (1.0 + fHe + xH)) * (dy[0] / xH) -
              epsilon * dHdz_si(z) / Hz + 3.0 * epsilon / (1.0 + z);
    } else {
      dy[1] = CT * std::pow(Trad, 4) * xH / (1.0 + xH + fHe) *
                  (Tmat - Trad) / (Hz * (1.0 + z)) +
              2.0 * Tmat / (1.0 + z);
    }
    return dy;
  };

  const int n_steps = std::max(200, p.recfast_steps);
  const double z_start = 10000.0;
  const double z_end = 0.0;
  const double dz = (z_end - z_start) / static_cast<double>(n_steps);
  std::array<double, 2> y{1.0, Tnow * (1.0 + z_start)};
  RecfastHistory hist;
  hist.z.reserve(static_cast<size_t>(n_steps + 1));
  hist.xe.reserve(static_cast<size_t>(n_steps + 1));
  hist.tm.reserve(static_cast<size_t>(n_steps + 1));

  double z = z_start;
  bool ode_started = false;
  for (int i = 0; i <= n_steps; ++i) {
    if (!ode_started) {
      y[0] = (z > 2000.0) ? 1.0 : saha_xe(z);
      y[1] = Tnow * (1.0 + z);
      if (z <= 2000.0 && y[0] <= 0.985) ode_started = true;
    }
    hist.z.push_back(z);
    hist.xe.push_back(std::clamp(y[0], 1.0e-12, 1.0 + fHe));
    hist.tm.push_back(y[1]);
    if (i == n_steps) break;
    if (!ode_started) {
      z += dz;
      continue;
    }
    const auto k1 = deriv(z, y);
    std::array<double, 2> yt{};
    for (int j = 0; j < 2; ++j) yt[j] = y[j] + 0.5 * dz * k1[j];
    const auto k2 = deriv(z + 0.5 * dz, yt);
    for (int j = 0; j < 2; ++j) yt[j] = y[j] + 0.5 * dz * k2[j];
    const auto k3 = deriv(z + 0.5 * dz, yt);
    for (int j = 0; j < 2; ++j) yt[j] = y[j] + dz * k3[j];
    const auto k4 = deriv(z + dz, yt);
    for (int j = 0; j < 2; ++j) {
      y[j] += dz * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]) / 6.0;
    }
    y[0] = std::clamp(y[0], 1.0e-12, 1.0 + fHe);
    y[1] = std::max(1.0, y[1]);
    z += dz;
  }
  return hist;
}
// >>>>>>>>>>>>>>>>>>>>>>> END RECFAST RECOMBINATION <<<<<<<<<<<<<<<<<<<<<<<

constexpr int STATE_SIZE = 7;
using State = std::array<double, STATE_SIZE>;

[[maybe_unused]] double simpson_integral(const Background &bg, double amax,
                                         int n = 65536) {
  if (n % 2 != 0) n += 1;
  const double h = amax / static_cast<double>(n);
  auto f = [&](double a) {
    return 1.0 / (bg.H0 * std::sqrt(bg.omega_m * a + bg.omega_rad +
                                     bg.omega_lambda * a * a * a * a));
  };
  double sum = f(0.0) + f(amax);
  for (int i = 1; i < n; ++i) {
    sum += (i % 2 == 0 ? 2.0 : 4.0) * f(i * h);
  }
  return sum * h / 3.0;
}

// >>>>>>>>>>>>>>>>>>>>>>>>> SCALAR FIELD EQUATIONS <<<<<<<<<<<<<<<<<<<<<<<<<
// AxiCLASS evolves phi_scf and phi_prime_scf with the Klein-Gordon equation
// phi'' + 2 calH phi' + a^2 dV/dphi = 0, where prime is d/d conformal time.
// This teaching solver uses theta = phi/f and the axion-like potential
// V(theta) = (m/H0)^2 f^2 [1 - cos(theta)]^n.
// The same helpers are used by both the background ODE and perturbation ODE.
double scalar_u(const Params &p, double theta) {
  if (!p.scalar_enabled || p.scalar_f == 0.0 || p.scalar_m_over_H0 == 0.0) {
    return 0.0;
  }
  return std::pow(std::max(0.0, 1.0 - std::cos(theta)), p.scalar_n);
}

double scalar_du(const Params &p, double theta) {
  if (!p.scalar_enabled || p.scalar_f == 0.0 || p.scalar_m_over_H0 == 0.0) {
    return 0.0;
  }
  const double base = std::max(0.0, 1.0 - std::cos(theta));
  if (p.scalar_n == 1) return std::sin(theta);
  return static_cast<double>(p.scalar_n) * std::pow(base, p.scalar_n - 1) *
         std::sin(theta);
}

double scalar_ddu(const Params &p, double theta) {
  if (!p.scalar_enabled || p.scalar_f == 0.0 || p.scalar_m_over_H0 == 0.0) {
    return 0.0;
  }
  const double base = std::max(0.0, 1.0 - std::cos(theta));
  if (p.scalar_n == 1) return std::cos(theta);
  const double n = static_cast<double>(p.scalar_n);
  const double term1 = n * (n - 1.0) * std::pow(base, p.scalar_n - 2) *
                       std::sin(theta) * std::sin(theta);
  const double term2 = n * std::pow(base, p.scalar_n - 1) * std::cos(theta);
  return term1 + term2;
}

double scalar_potential(const Params &p, double theta) {
  return p.scalar_m_over_H0 * p.scalar_m_over_H0 * p.scalar_f * p.scalar_f *
         scalar_u(p, theta);
}

[[maybe_unused]] double scalar_dpotential_dtheta(const Params &p, double theta) {
  return p.scalar_m_over_H0 * p.scalar_m_over_H0 * p.scalar_f * p.scalar_f *
         scalar_du(p, theta);
}

[[maybe_unused]] double scalar_ddpotential_dtheta2(const Params &p, double theta) {
  return p.scalar_m_over_H0 * p.scalar_m_over_H0 * p.scalar_f * p.scalar_f *
         scalar_ddu(p, theta);
}

double scalar_rho(const Params &p, double a, double theta, double theta_prime) {
  if (!p.scalar_enabled) return 0.0;
  return 0.5 * p.scalar_f * p.scalar_f * theta_prime * theta_prime / (a * a) +
         scalar_potential(p, theta);
}

[[maybe_unused]] double scalar_pressure(const Params &p, double a, double theta,
                                        double theta_prime) {
  if (!p.scalar_enabled) return 0.0;
  return 0.5 * p.scalar_f * p.scalar_f * theta_prime * theta_prime / (a * a) -
         scalar_potential(p, theta);
}
// >>>>>>>>>>>>>>>>>>>>>>> END SCALAR FIELD EQUATIONS <<<<<<<<<<<<<<<<<<<<<<<

struct BackgroundSolution {
  std::vector<double> a_grid;
  std::vector<double> eta_grid;
  std::vector<double> theta_grid;
  std::vector<double> theta_prime_grid;

  double eta_of_a(double a) const {
    if (a <= a_grid.front()) return eta_grid.front();
    if (a >= a_grid.back()) return eta_grid.back();
    auto it = std::lower_bound(a_grid.begin(), a_grid.end(), a);
    const size_t hi = static_cast<size_t>(it - a_grid.begin());
    const size_t lo = hi - 1;
    const double t = (a - a_grid[lo]) / (a_grid[hi] - a_grid[lo]);
    return (1.0 - t) * eta_grid[lo] + t * eta_grid[hi];
  }

  double a_of_eta(double eta) const {
    if (eta <= eta_grid.front()) return a_grid.front();
    if (eta >= eta_grid.back()) return a_grid.back();
    auto it = std::lower_bound(eta_grid.begin(), eta_grid.end(), eta);
    const size_t hi = static_cast<size_t>(it - eta_grid.begin());
    const size_t lo = hi - 1;
    const double t = (eta - eta_grid[lo]) / (eta_grid[hi] - eta_grid[lo]);
    return (1.0 - t) * a_grid[lo] + t * a_grid[hi];
  }

  void scalar_at_a(double a, double &theta, double &theta_prime) const {
    if (a <= a_grid.front()) {
      theta = theta_grid.front();
      theta_prime = theta_prime_grid.front();
      return;
    }
    if (a >= a_grid.back()) {
      theta = theta_grid.back();
      theta_prime = theta_prime_grid.back();
      return;
    }
    auto it = std::lower_bound(a_grid.begin(), a_grid.end(), a);
    const size_t hi = static_cast<size_t>(it - a_grid.begin());
    const size_t lo = hi - 1;
    const double t = (a - a_grid[lo]) / (a_grid[hi] - a_grid[lo]);
    theta = (1.0 - t) * theta_grid[lo] + t * theta_grid[hi];
    theta_prime = (1.0 - t) * theta_prime_grid[lo] + t * theta_prime_grid[hi];
  }
};

BackgroundSolution integrate_background_ode(const Params &p, const Background &bg) {
  const int n = std::max(1000, p.background_steps);
  const double a_min = 1.0e-8;
  const double loga_min = std::log(a_min);
  const double h = -loga_min / static_cast<double>(n);

  BackgroundSolution sol;
  sol.a_grid.reserve(static_cast<size_t>(n + 1));
  sol.eta_grid.reserve(static_cast<size_t>(n + 1));
  sol.theta_grid.reserve(static_cast<size_t>(n + 1));
  sol.theta_prime_grid.reserve(static_cast<size_t>(n + 1));

  std::array<double, 3> y{
      a_min / (bg.H0 * std::sqrt(bg.omega_rad)),
      p.scalar_theta_ini,
      p.scalar_theta_prime_ini,
  };

  auto deriv = [&](double loga, const std::array<double, 3> &state) {
    const double a = std::exp(loga);
    const double theta = state[1];
    const double theta_prime = state[2];
    const double e2 = bg.omega_m / (a * a * a) +
                      bg.omega_rad / (a * a * a * a) + bg.omega_lambda +
                      scalar_rho(p, a, theta, theta_prime);
    if (e2 <= 0.0) throw std::runtime_error("background H(a)^2 became non-positive");
    const double e = std::sqrt(e2);
    const double calH = a * bg.H0 * e;

    std::array<double, 3> dy{};
    dy[0] = 1.0 / (a * bg.H0 * e);
    dy[1] = theta_prime / (a * bg.H0 * e);
    dy[2] = (-2.0 * calH * theta_prime -
             a * a * p.scalar_m_over_H0 * p.scalar_m_over_H0 *
                 scalar_du(p, theta)) /
            (a * bg.H0 * e);
    return dy;
  };

  double loga = loga_min;
  for (int i = 0; i <= n; ++i) {
    const double a = std::exp(loga);
    sol.a_grid.push_back(a);
    sol.eta_grid.push_back(y[0]);
    sol.theta_grid.push_back(y[1]);
    sol.theta_prime_grid.push_back(y[2]);
    if (i == n) break;

    const auto k1 = deriv(loga, y);
    std::array<double, 3> yt{};
    for (int j = 0; j < 3; ++j) yt[j] = y[j] + 0.5 * h * k1[j];
    const auto k2 = deriv(loga + 0.5 * h, yt);
    for (int j = 0; j < 3; ++j) yt[j] = y[j] + 0.5 * h * k2[j];
    const auto k3 = deriv(loga + 0.5 * h, yt);
    for (int j = 0; j < 3; ++j) yt[j] = y[j] + h * k3[j];
    const auto k4 = deriv(loga + h, yt);
    for (int j = 0; j < 3; ++j) {
      y[j] += h * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]) / 6.0;
    }
    loga += h;
  }
  return sol;
}

class CubicSpline {
 public:
  CubicSpline() = default;

  CubicSpline(const std::vector<double> &x_in, const std::vector<double> &y_in)
      : x(x_in), y(y_in), y2(x_in.size(), 0.0) {
    if (x.size() < 3 || x.size() != y.size()) {
      throw std::runtime_error("bad spline input");
    }
    std::vector<double> u(x.size() - 1, 0.0);
    y2[0] = 0.0;
    u[0] = 0.0;
    for (size_t i = 1; i + 1 < x.size(); ++i) {
      const double sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
      const double p = sig * y2[i - 1] + 2.0;
      y2[i] = (sig - 1.0) / p;
      const double dd =
          (y[i + 1] - y[i]) / (x[i + 1] - x[i]) -
          (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
      u[i] = (6.0 * dd / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
    }
    y2.back() = 0.0;
    for (int k = static_cast<int>(x.size()) - 2; k >= 0; --k) {
      y2[k] = y2[k] * y2[k + 1] + u[k];
    }
  }

  double operator()(double xx) const {
    auto it = std::lower_bound(x.begin(), x.end(), xx);
    size_t khi = 0;
    if (it == x.begin()) {
      khi = 1;
    } else if (it == x.end()) {
      khi = x.size() - 1;
    } else {
      khi = static_cast<size_t>(it - x.begin());
    }
    size_t klo = khi - 1;
    const double h = x[khi] - x[klo];
    const double a = (x[khi] - xx) / h;
    const double b = (xx - x[klo]) / h;
    return a * y[klo] + b * y[khi] +
           ((a * a * a - a) * y2[klo] + (b * b * b - b) * y2[khi]) *
               h * h / 6.0;
  }

 private:
  std::vector<double> x, y, y2;
};

class TwoFluidModel {
 public:
  explicit TwoFluidModel(const Params &params)
      : p(params), bg(params), bgsol(integrate_background_ode(params, bg)) {
    recfast_history = compute_recfast_history(p, bg);
    have_recfast_history = true;
    recombination_scales = compute_recombination_scales_from_recfast();
    arec = 1.0 / (1.0 + recombination_scales.z_star);
    aeq = bg.omega_rad / bg.omega_m;
    etatoday = bgsol.eta_of_a(1.0);
    etarec = recombination_scales.eta_star;
    etastar = etarec;
    eta_eq = bgsol.eta_of_a(aeq);
    xs = std::sqrt(std::max(0.0, recombination_scales.silk_length_sq));
    visibility_sigma_eta = recombination_scales.visibility_sigma_eta;
  }

  struct Spectrum {
    std::vector<int> ell;
    std::vector<double> cl;
  };

  void write_recfast_history(const std::string &path) const {
    if (!have_recfast_history) return;
    std::ofstream out(path);
    if (!out) throw std::runtime_error("could not open RECFAST output: " + path);
    out << "# z xe Tmat_K\n";
    for (size_t i = 0; i < recfast_history.z.size(); ++i) {
      out << std::setprecision(17) << recfast_history.z[i] << " "
          << recfast_history.xe[i] << " " << recfast_history.tm[i] << "\n";
    }
  }

  void write_scale_summary(std::ostream &out) const {
    const double H_eq = bg.H0 * std::sqrt(bg.omega_rad / std::pow(aeq, 4) +
                                          bg.omega_m / std::pow(aeq, 3) +
                                          bg.omega_lambda);
    const double k_eq = aeq * H_eq;
    const double cs_rec = sound_speed_rec();
    const double k_D =
        (recombination_scales.silk_length_sq > 0.0)
            ? 1.0 / std::sqrt(recombination_scales.silk_length_sq)
            : 0.0;
    const double k_width =
        (cs_rec * visibility_sigma_eta > 0.0)
            ? 1.0 / (cs_rec * visibility_sigma_eta)
            : 0.0;
    out << std::setprecision(8)
        << "Scales from RECFAST/background:\n"
        << "  z_star = " << recombination_scales.z_star << "\n"
        << "  N_eff = " << p.N_eff << "\n"
        << "  Omega_nu = " << bg.omega_nu << "\n"
        << "  eta_star = " << etarec << " Mpc\n"
        << "  eta_0 = " << etatoday << " Mpc\n"
        << "  chi_star = " << (etatoday - etarec) << " Mpc\n"
        << "  a_eq = " << aeq << "\n"
        << "  k_eq = " << k_eq << " Mpc^-1\n"
        << "  sigma_eta = " << visibility_sigma_eta << " Mpc\n"
        << "  k_width = " << k_width << " Mpc^-1\n"
        << "  k_D = " << k_D << " Mpc^-1\n"
        << "  sound_speed_star = " << cs_rec << "\n";
  }

  Spectrum compute_spectrum() const {
    std::vector<double> source_k, sw, dop;
    solve_sources(source_k, sw, dop);

    CubicSpline sw_spline(source_k, sw);
    CubicSpline dop_spline(source_k, dop);

    const double kmin = 1.0 / etatoday;
    const double kmax = 1.0e4 / etatoday;
    const double dk = 2.0 * PI / (p.k_samples_per_period * etatoday);
    std::vector<double> k_grid;
    for (double k = kmin; k <= kmax + 0.5 * dk; k += dk) k_grid.push_back(k);

    std::vector<int> ell = ell_list();
    std::vector<double> cl(ell.size(), 0.0);
    const double chi = etatoday - etastar;
    const bool sw_only = p.source_mode == "sw";
    const bool doppler_only = p.source_mode == "doppler";
    const bool isw_only = p.source_mode == "isw";
    const bool include_primary = !isw_only;
    const bool include_sw = p.source_mode == "full" || sw_only;
    const bool include_doppler = p.source_mode == "full" || doppler_only;
    const bool include_isw = p.source_mode == "full" || isw_only;
    std::vector<std::vector<double>> isw_source;
    if (include_isw) {
      isw_source = compute_isw_sources(ell, source_k);
    }
    // In CMB TT without explicit reionization, this amplitude should be read as
    // A_s.  With reionization included phenomenologically, high-ell TT measures
    // approximately A_s * exp(-2 tau_reio).
    const double amp = 4.0 * PI * p.As * T_CMB_MICROK * T_CMB_MICROK;
    const double cs_rec = sound_speed_rec();
    const double finite_width_damping =
        cs_rec * cs_rec * visibility_sigma_eta * visibility_sigma_eta;

    std::vector<double> xarg(k_grid.size()), damping(k_grid.size()), swsrc(k_grid.size()),
        dopv(k_grid.size());
    for (size_t i = 0; i < k_grid.size(); ++i) {
      xarg[i] = k_grid[i] * chi;
      damping[i] = std::exp(
          -k_grid[i] * k_grid[i] * (2.0 * xs * xs + finite_width_damping));
      swsrc[i] = include_sw ? sw_spline(k_grid[i]) : 0.0;
      dopv[i] = include_doppler ? dop_spline(k_grid[i]) : 0.0;
    }

    std::vector<int> needed_orders = required_bessel_orders(ell);
    std::vector<std::vector<double>> bessel;
    if (p.bessel_class_memory) {
      bessel = bessel_class_memory_table(needed_orders, xarg, p.hyper_sampling_flat);
    } else if (!p.bessel_x_cache.empty()) {
      if (!read_bessel_x_cache(p.bessel_x_cache, needed_orders, xarg, bessel)) {
        std::cerr << "Building fixed-x Bessel cache " << p.bessel_x_cache << "\n";
        write_bessel_x_cache(p.bessel_x_cache, needed_orders, p.bessel_x_dx,
                             p.bessel_x_max);
        if (!read_bessel_x_cache(p.bessel_x_cache, needed_orders, xarg, bessel)) {
          throw std::runtime_error("new fixed-x Bessel cache is not valid for this grid");
        }
      }
    } else if (!p.bessel_cache.empty()) {
      if (!read_bessel_cache(p.bessel_cache, needed_orders, xarg, bessel)) {
        std::cerr << "Building Bessel cache " << p.bessel_cache << "\n";
        bessel = spherical_bessel_table(needed_orders, xarg);
        write_bessel_cache(p.bessel_cache, needed_orders, xarg, bessel);
      } else {
        std::cerr << "Loaded Bessel cache " << p.bessel_cache << "\n";
      }
    } else {
      bessel = spherical_bessel_table(needed_orders, xarg);
    }
    auto order_index = [&](int order) {
      auto it = std::lower_bound(needed_orders.begin(), needed_orders.end(), order);
      if (it == needed_orders.end() || *it != order) {
        throw std::runtime_error("missing Bessel order");
      }
      return static_cast<size_t>(it - needed_orders.begin());
    };

    for (size_t j = 0; j < ell.size(); ++j) {
      const int l = ell[j];
      CubicSpline isw_spline;
      if (include_isw) {
        isw_spline = CubicSpline(source_k, isw_source[j]);
      }
      const size_t idx_lm = order_index(l - 1);
      const size_t idx_l = order_index(l);
      const size_t idx_lp = order_index(l + 1);
      std::vector<double> integrand(k_grid.size(), 0.0);
      for (size_t i = 0; i < k_grid.size(); ++i) {
        const double x = xarg[i];
        const double jl = bessel[idx_l][i];
        const double jlm = bessel[idx_lm][i];
        const double jlp = bessel[idx_lp][i];
        const double deriv = -jl / (2.0 * x) + 0.5 * (jlm - jlp);
        const double primary_src = swsrc[i] * jl + dopv[i] * deriv;
        const double isw = include_isw ? isw_spline(k_grid[i]) : 0.0;
        const double src =
            (include_primary ? std::sqrt(damping[i]) * primary_src : 0.0) + isw;
        integrand[i] = src * src / k_grid[i];
      }
      double val = 0.0;
      for (size_t i = 0; i + 1 < k_grid.size(); ++i) {
        val += 0.5 * (integrand[i] + integrand[i + 1]) *
               (k_grid[i + 1] - k_grid[i]);
      }
      const double tilt = std::pow(static_cast<double>(l) / ELL_PIVOT, p.ns - FID_NS);
      cl[j] = amp * tilt * val;
    }
    return Spectrum{ell, cl};
  }

 private:
  Params p;
  Background bg;
  BackgroundSolution bgsol;
  RecfastHistory recfast_history;
  bool have_recfast_history = false;
  RecombinationScales recombination_scales;
  double arec, aeq, etatoday, etarec, etastar, eta_eq, xs, visibility_sigma_eta;

  double sound_speed_rec() const {
    const double R_rec = 0.75 * (bg.omega_b / bg.omega_gamma) * arec;
    return 1.0 / std::sqrt(3.0 * (1.0 + R_rec));
  }

  RecombinationScales compute_recombination_scales_from_recfast() const {
    constexpr double Mpc = 3.0856775814913673e22;
    constexpr double G = 6.6742e-11;
    constexpr double m_H = 1.673575e-27;
    constexpr double sigma_T = 6.6524616e-29;

    const size_t n = recfast_history.z.size();
    if (n < 3) throw std::runtime_error("RECFAST history is too short");

    const double H0_si = 100.0 * p.h * 1000.0 / Mpc;
    const double mu_H = 1.0 / (1.0 - p.helium_y_p);
    const double nH0 = 3.0 * H0_si * H0_si * bg.omega_b /
                       (8.0 * PI * G * mu_H * m_H);

    std::vector<double> eta(n), kappa_dot(n), visibility(n), tau(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
      const double z = recfast_history.z[i];
      const double a = 1.0 / (1.0 + z);
      eta[i] = bgsol.eta_of_a(a);
      const double ne = recfast_history.xe[i] * nH0 * std::pow(1.0 + z, 3);
      kappa_dot[i] = a * ne * sigma_T * Mpc;
    }

    for (int i = static_cast<int>(n) - 2; i >= 0; --i) {
      const double deta = eta[static_cast<size_t>(i + 1)] - eta[static_cast<size_t>(i)];
      tau[static_cast<size_t>(i)] =
          tau[static_cast<size_t>(i + 1)] +
          0.5 * (kappa_dot[static_cast<size_t>(i)] +
                 kappa_dot[static_cast<size_t>(i + 1)]) *
              deta;
    }

    size_t peak = 0;
    double gmax = -1.0;
    for (size_t i = 0; i < n; ++i) {
      visibility[i] = kappa_dot[i] * std::exp(-tau[i]);
      const double z = recfast_history.z[i];
      if (z > 500.0 && z < 2000.0 && visibility[i] > gmax) {
        gmax = visibility[i];
        peak = i;
      }
    }
    if (gmax <= 0.0) throw std::runtime_error("could not find RECFAST visibility peak");

    double norm = 0.0;
    double mean = 0.0;
    double second = 0.0;
    for (size_t i = 0; i + 1 < n; ++i) {
      const double zmid = 0.5 * (recfast_history.z[i] + recfast_history.z[i + 1]);
      if (zmid < 500.0 || zmid > 2000.0) continue;
      const double deta = eta[i + 1] - eta[i];
      const double g0 = visibility[i];
      const double g1 = visibility[i + 1];
      const double e0 = eta[i];
      const double e1 = eta[i + 1];
      norm += 0.5 * (g0 + g1) * deta;
      mean += 0.5 * (g0 * e0 + g1 * e1) * deta;
      second += 0.5 * (g0 * e0 * e0 + g1 * e1 * e1) * deta;
    }
    if (norm <= 0.0) throw std::runtime_error("RECFAST visibility normalization failed");
    mean /= norm;
    second /= norm;
    const double sigma_eta =
        std::sqrt(std::max(0.0, second - mean * mean));

    const double eta_star = eta[peak];
    double silk_length_sq = 0.0;
    for (size_t i = 0; i + 1 < peak; ++i) {
      const double a0 = 1.0 / (1.0 + recfast_history.z[i]);
      const double a1 = 1.0 / (1.0 + recfast_history.z[i + 1]);
      const double R0 = 0.75 * (bg.omega_b / bg.omega_gamma) * a0;
      const double R1 = 0.75 * (bg.omega_b / bg.omega_gamma) * a1;
      const double fac0 =
          (R0 * R0 + 16.0 * (1.0 + R0) / 15.0) /
          (6.0 * kappa_dot[i] * (1.0 + R0) * (1.0 + R0));
      const double fac1 =
          (R1 * R1 + 16.0 * (1.0 + R1) / 15.0) /
          (6.0 * kappa_dot[i + 1] * (1.0 + R1) * (1.0 + R1));
      silk_length_sq += 0.5 * (fac0 + fac1) * (eta[i + 1] - eta[i]);
    }

    return RecombinationScales{
        recfast_history.z[peak],
        eta_star,
        sigma_eta,
        silk_length_sq,
    };
  }

  State rhs(double eta, const State &u, double k) const {
    const double delta_gamma = u[0];
    const double delta_nu = u[1];
    const double delta_c = u[2];
    const double theta_gammab = u[3];
    const double theta_nu = u[4];
    const double theta_c = u[5];
    const double phi = u[6];
    // This is approximately true during tight coupling.
    const double delta_b = 0.75 * delta_gamma;
    const double a = bgsol.a_of_eta(eta);
    const double rho_gamma = bg.omega_gamma / (a * a * a * a);
    const double rho_nu = bg.omega_nu / (a * a * a * a);
    const double rho_b = bg.omega_b / (a * a * a);
    const double rho_c = bg.omega_cdm / (a * a * a);
    const double rho_lambda = bg.omega_lambda;
    const double rho_tot = rho_gamma + rho_nu + rho_b + rho_c + rho_lambda;
    const double H_of_a = bg.H0 * std::sqrt(rho_tot);
    const double Hconf = a * H_of_a;
    const double R = (4.0 / 3.0) * rho_gamma / rho_b;
    const double delta_tot =
        (rho_gamma * delta_gamma + rho_nu * delta_nu + rho_b * delta_b +
         rho_c * delta_c) /
        rho_tot;
    const double phi_p =
        (-k * k * phi - 1.5 * Hconf * Hconf * delta_tot) /
            (3.0 * Hconf) -
        Hconf * phi;
    State du{};
    du[0] = -(4.0 / 3.0) * theta_gammab + 4.0 * phi_p;
    du[1] = -(4.0 / 3.0) * theta_nu + 4.0 * phi_p;
    du[2] = -theta_c + 3.0 * phi_p;
    du[3] = -Hconf * theta_gammab / (1.0 + R) +
            k * k * phi + 0.25 * R * k * k * delta_gamma / (1.0 + R);
    du[4] = k * k * phi + 0.25 * k * k * delta_nu;
    du[5] = -Hconf * theta_c + k * k * phi;
    du[6] = phi_p;
    return du;
  }

  std::pair<double, State> initial_conditions(double k) const {
    const double eta_i = std::min(1.0e-3 / k, 1.0e-3 * eta_eq);
    const double phi_i = -2.0 / 3.0;
    const double dg_i = -2.0 * phi_i;
    const double dnu_i = dg_i;
    const double dc_i = 0.75 * dg_i;
    const double theta_i = 0.5 * eta_i * k * k * phi_i;
    return {eta_i, State{dg_i, dnu_i, dc_i, theta_i, theta_i, theta_i, phi_i}};
  }

  State rk45_step(double eta, const State &u, double h, double k,
                  State &err) const {
    static constexpr double c2 = 1.0 / 5.0, c3 = 3.0 / 10.0, c4 = 4.0 / 5.0,
                            c5 = 8.0 / 9.0;
    static constexpr double a21 = 1.0 / 5.0;
    static constexpr double a31 = 3.0 / 40.0, a32 = 9.0 / 40.0;
    static constexpr double a41 = 44.0 / 45.0, a42 = -56.0 / 15.0,
                            a43 = 32.0 / 9.0;
    static constexpr double a51 = 19372.0 / 6561.0,
                            a52 = -25360.0 / 2187.0,
                            a53 = 64448.0 / 6561.0,
                            a54 = -212.0 / 729.0;
    static constexpr double a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0,
                            a63 = 46732.0 / 5247.0,
                            a64 = 49.0 / 176.0,
                            a65 = -5103.0 / 18656.0;
    static constexpr double a71 = 35.0 / 384.0, a73 = 500.0 / 1113.0,
                            a74 = 125.0 / 192.0,
                            a75 = -2187.0 / 6784.0,
                            a76 = 11.0 / 84.0;
    static constexpr double b1 = 35.0 / 384.0, b3 = 500.0 / 1113.0,
                            b4 = 125.0 / 192.0,
                            b5 = -2187.0 / 6784.0, b6 = 11.0 / 84.0;
    static constexpr double bs1 = 5179.0 / 57600.0,
                            bs3 = 7571.0 / 16695.0,
                            bs4 = 393.0 / 640.0,
                            bs5 = -92097.0 / 339200.0,
                            bs6 = 187.0 / 2100.0, bs7 = 1.0 / 40.0;

    auto add = [](const State &u0, double hstep,
                  const std::vector<std::pair<double, State>> &terms) {
      State out = u0;
      for (int i = 0; i < STATE_SIZE; ++i) {
        for (const auto &term : terms) out[i] += hstep * term.first * term.second[i];
      }
      return out;
    };

    const State k1 = rhs(eta, u, k);
    const State k2 = rhs(eta + c2 * h, add(u, h, {{a21, k1}}), k);
    const State k3 = rhs(eta + c3 * h, add(u, h, {{a31, k1}, {a32, k2}}), k);
    const State k4 =
        rhs(eta + c4 * h, add(u, h, {{a41, k1}, {a42, k2}, {a43, k3}}), k);
    const State k5 = rhs(eta + c5 * h,
                         add(u, h, {{a51, k1}, {a52, k2}, {a53, k3}, {a54, k4}}),
                         k);
    const State k6 = rhs(eta + h,
                         add(u, h, {{a61, k1}, {a62, k2}, {a63, k3}, {a64, k4},
                                    {a65, k5}}),
                         k);
    const State k7 = rhs(eta + h,
                         add(u, h, {{a71, k1}, {a73, k3}, {a74, k4}, {a75, k5},
                                    {a76, k6}}),
                         k);

    State y5{}, y4{};
    for (int i = 0; i < STATE_SIZE; ++i) {
      y5[i] = u[i] + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] +
                          b5 * k5[i] + b6 * k6[i]);
      y4[i] = u[i] + h * (bs1 * k1[i] + bs3 * k3[i] + bs4 * k4[i] +
                          bs5 * k5[i] + bs6 * k6[i] + bs7 * k7[i]);
      err[i] = y5[i] - y4[i];
    }
    return y5;
  }

  State integrate_source(double k) const {
    auto init = initial_conditions(k);
    double eta = init.first;
    State u = init.second;
    double h = std::min(1.0, (etarec - eta) / 20.0);
    const double rtol = 2.0e-8;
    const double atol = 2.0e-10;
    int steps = 0;
    while (eta < etarec && steps < 200000) {
      if (eta + h > etarec) h = etarec - eta;
      State err{};
      State trial = rk45_step(eta, u, h, k, err);
      double err_norm = 0.0;
      for (int i = 0; i < STATE_SIZE; ++i) {
        const double scale = atol + rtol * std::max(std::abs(u[i]), std::abs(trial[i]));
        err_norm = std::max(err_norm, std::abs(err[i]) / scale);
      }
      if (err_norm <= 1.0) {
        eta += h;
        u = trial;
      }
      const double factor =
          (err_norm == 0.0) ? 5.0 : std::clamp(0.9 * std::pow(err_norm, -0.2), 0.2, 5.0);
      h *= factor;
      ++steps;
    }
    if (steps >= 200000) throw std::runtime_error("ODE integration did not converge");
    return u;
  }

  std::vector<double> phi_prime_history(double k,
                                        const std::vector<double> &eta_nodes) const {
    auto init = initial_conditions(k);
    double eta = init.first;
    State u = init.second;
    double h = std::max(1.0e-6, std::min(2.0, (eta_nodes.front() - eta) / 20.0));
    const double rtol = 1.0e-6;
    const double atol = 1.0e-9;
    int steps = 0;
    const int max_steps = 3000000;
    std::vector<double> history(eta_nodes.size(), 0.0);

    for (size_t target_index = 0; target_index < eta_nodes.size(); ++target_index) {
      const double target = eta_nodes[target_index];
      while (eta < target && steps < max_steps) {
        if (eta + h > target) h = target - eta;
        State err{};
        State trial = rk45_step(eta, u, h, k, err);
        double err_norm = 0.0;
        for (int i = 0; i < STATE_SIZE; ++i) {
          const double scale = atol + rtol * std::max(std::abs(u[i]), std::abs(trial[i]));
          err_norm = std::max(err_norm, std::abs(err[i]) / scale);
        }
        if (err_norm <= 1.0) {
          eta += h;
          u = trial;
        }
        const double factor =
            (err_norm == 0.0) ? 5.0
                              : std::clamp(0.9 * std::pow(err_norm, -0.2), 0.2, 5.0);
        h *= factor;
        ++steps;
      }
      if (steps >= max_steps) {
        throw std::runtime_error("ISW ODE integration did not converge");
      }
      history[target_index] = rhs(target, u, k)[6];
      h = std::min(h, std::max(1.0e-6, (target_index + 1 < eta_nodes.size()
                                            ? eta_nodes[target_index + 1] - target
                                            : h)));
    }
    return history;
  }

  std::vector<std::vector<double>> compute_isw_sources(
      const std::vector<int> &ell,
      const std::vector<double> &source_k) const {
    std::vector<double> eta_nodes(static_cast<size_t>(std::max(2, p.n_isw)), 0.0);
    for (size_t i = 0; i < eta_nodes.size(); ++i) {
      const double t = static_cast<double>(i) / static_cast<double>(eta_nodes.size() - 1);
      eta_nodes[i] = etarec + t * (etatoday - etarec);
    }

    std::vector<std::vector<double>> phi_p(
        source_k.size(), std::vector<double>(eta_nodes.size(), 0.0));
    for (size_t ik = 0; ik < source_k.size(); ++ik) {
      if (source_k[ik] * etatoday <= p.isw_k_eta0_max) {
        phi_p[ik] = phi_prime_history(source_k[ik], eta_nodes);
      }
    }

    std::vector<double> xarg;
    xarg.reserve(source_k.size() * eta_nodes.size());
    for (double k : source_k) {
      for (double eta_node : eta_nodes) {
        xarg.push_back(k * (etatoday - eta_node));
      }
    }

    std::vector<int> orders(ell.begin(), ell.end());
    std::vector<std::vector<double>> bessel = spherical_bessel_table(orders, xarg);
    std::vector<std::vector<double>> isw(
        ell.size(), std::vector<double>(source_k.size(), 0.0));

    for (size_t il = 0; il < ell.size(); ++il) {
      for (size_t ik = 0; ik < source_k.size(); ++ik) {
        double val = 0.0;
        const size_t offset = ik * eta_nodes.size();
        for (size_t ix = 0; ix + 1 < eta_nodes.size(); ++ix) {
          const double f0 = 2.0 * phi_p[ik][ix] * bessel[il][offset + ix];
          const double f1 = 2.0 * phi_p[ik][ix + 1] * bessel[il][offset + ix + 1];
          val += 0.5 * (f0 + f1) * (eta_nodes[ix + 1] - eta_nodes[ix]);
        }
        isw[il][ik] = val;
      }
    }
    return isw;
  }

  void solve_sources(std::vector<double> &k_values, std::vector<double> &sw,
                     std::vector<double> &dop) const {
    k_values.resize(p.n_source);
    sw.resize(p.n_source);
    dop.resize(p.n_source);

    for (int i = 0; i < p.n_source; ++i) {
      const double grid = 4.0 * static_cast<double>(i) / static_cast<double>(p.n_source - 1);
      k_values[i] = std::pow(10.0, grid) / etatoday;
      State u = integrate_source(k_values[i]);
      sw[i] = u[6] + u[0] / 4.0;
      dop[i] = -u[3] / k_values[i];
    }
  }

  std::vector<int> ell_list() const {
    if (p.ell_grid == "class") {
      return class_ell_list((p.ell_max > 0) ? p.ell_max : 2500);
    }

    if (p.ell_grid != "sparse") {
      throw std::runtime_error("ell-grid must be 'class' or 'sparse'");
    }

    if (p.ell_min > 0 || p.ell_max > 0) {
      const int ell_min = (p.ell_min > 0) ? p.ell_min : 2;
      const int ell_max = (p.ell_max > 0) ? p.ell_max : 2000;
      const int ell_step = std::max(1, p.ell_step);
      if (ell_max < ell_min) throw std::runtime_error("ell-max must be >= ell-min");
      std::vector<int> ell;
      for (int l = ell_min; l <= ell_max; l += ell_step) ell.push_back(l);
      if (ell.back() != ell_max) ell.push_back(ell_max);
      return ell;
    }

    std::vector<int> ell{2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20};
    for (int l = 30; l <= 100; l += 10) ell.push_back(l);
    for (int l = 120; l <= 500; l += 20) ell.push_back(l);
    for (int l = 550; l <= 2000; l += 50) ell.push_back(l);
    std::sort(ell.begin(), ell.end());
    ell.erase(std::unique(ell.begin(), ell.end()), ell.end());
    return ell;
  }

  static std::vector<int> class_ell_list(int ell_max) {
    constexpr double l_logstep = 1.12;
    constexpr int l_linstep = 40;
    constexpr double angular_rescaling = 1.0;

    if (ell_max < 2) throw std::runtime_error("CLASS ell grid needs ell-max >= 2");

    std::vector<int> ell{2};
    int current_l = 2;
    int increment = std::max(
        static_cast<int>(current_l *
                         (std::pow(l_logstep, angular_rescaling) - 1.0)),
        1);

    while (((current_l + increment) < ell_max) &&
           (increment < static_cast<int>(l_linstep * angular_rescaling))) {
      current_l += increment;
      ell.push_back(current_l);
      increment = std::max(
          static_cast<int>(current_l *
                           (std::pow(l_logstep, angular_rescaling) - 1.0)),
          1);
    }

    increment = static_cast<int>(l_linstep * angular_rescaling);
    while ((current_l + increment) <= ell_max) {
      current_l += increment;
      ell.push_back(current_l);
    }

    if (current_l != ell_max) ell.push_back(ell_max);
    return ell;
  }

  static std::vector<int> required_bessel_orders(const std::vector<int> &ell) {
    std::set<int> orders;
    for (int l : ell) {
      orders.insert(l - 1);
      orders.insert(l);
      orders.insert(l + 1);
    }
    return std::vector<int>(orders.begin(), orders.end());
  }

  static std::vector<std::vector<double>> spherical_bessel_table(
      const std::vector<int> &orders, const std::vector<double> &xarg) {
    const int max_order = orders.back();
    std::vector<std::vector<double>> table(orders.size(),
                                           std::vector<double>(xarg.size(), 0.0));

    for (size_t ix = 0; ix < xarg.size(); ++ix) {
      const double x = xarg[ix];
      if (x == 0.0) {
        for (size_t io = 0; io < orders.size(); ++io) {
          table[io][ix] = (orders[io] == 0) ? 1.0 : 0.0;
        }
        continue;
      }

      const int L = std::max(max_order + 80, static_cast<int>(std::ceil(x)) + 80);
      std::vector<double> j(static_cast<size_t>(L + 2), 0.0);
      j[L + 1] = 0.0;
      j[L] = 1.0;
      for (int n = L; n >= 1; --n) {
        j[n - 1] = (2.0 * n + 1.0) * j[n] / x - j[n + 1];
        const double mag = std::max(std::abs(j[n - 1]), std::abs(j[n]));
        if (mag > 1.0e100) {
          for (int m = n - 1; m <= L + 1; ++m) j[m] *= 1.0e-100;
        }
      }

      const double true_j0 = std::sin(x) / x;
      const double true_j1 = std::sin(x) / (x * x) - std::cos(x) / x;
      double norm = 0.0;
      if (std::abs(j[0]) > std::abs(j[1]) && std::abs(j[0]) > 0.0) {
        norm = true_j0 / j[0];
      } else if (std::abs(j[1]) > 0.0) {
        norm = true_j1 / j[1];
      } else {
        norm = true_j0 / j[0];
      }

      for (size_t io = 0; io < orders.size(); ++io) {
        table[io][ix] = norm * j[orders[io]];
      }
    }
    return table;
  }

  struct BesselHermiteGrid {
    std::vector<std::vector<double>> value;
    std::vector<std::vector<double>> slope;
  };

  static std::vector<std::vector<double>> interpolate_uniform_bessel_hermite(
      const BesselHermiteGrid &grid,
      const std::vector<double> &xarg,
      double dx) {
    const size_t n_x = grid.value.front().size();
    std::vector<std::vector<double>> table(
        grid.value.size(), std::vector<double>(xarg.size(), 0.0));

    for (size_t io = 0; io < grid.value.size(); ++io) {
      const std::vector<double> &row = grid.value[io];
      const std::vector<double> &slope = grid.slope[io];
      for (size_t ix = 0; ix < xarg.size(); ++ix) {
        const double u = xarg[ix] / dx;
        size_t i = static_cast<size_t>(std::floor(u));
        if (i + 1 >= n_x) i = n_x - 2;
        const double t = u - static_cast<double>(i);
        const double t2 = t * t;
        const double t3 = t2 * t;

        const double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
        const double h10 = t3 - 2.0 * t2 + t;
        const double h01 = -2.0 * t3 + 3.0 * t2;
        const double h11 = t3 - t2;
        table[io][ix] = h00 * row[i] + h10 * dx * slope[i] +
                        h01 * row[i + 1] + h11 * dx * slope[i + 1];
      }
    }
    return table;
  }

  static BesselHermiteGrid spherical_bessel_hermite_grid(
      const std::vector<int> &orders, const std::vector<double> &xarg) {
    const int max_order = orders.back();
    BesselHermiteGrid grid;
    grid.value.assign(orders.size(), std::vector<double>(xarg.size(), 0.0));
    grid.slope.assign(orders.size(), std::vector<double>(xarg.size(), 0.0));

    for (size_t ix = 0; ix < xarg.size(); ++ix) {
      const double x = xarg[ix];
      if (x == 0.0) {
        for (size_t io = 0; io < orders.size(); ++io) {
          const int order = orders[io];
          grid.value[io][ix] = (order == 0) ? 1.0 : 0.0;
          grid.slope[io][ix] = (order == 1) ? 1.0 / 3.0 : 0.0;
        }
        continue;
      }

      const int L = std::max(max_order + 81, static_cast<int>(std::ceil(x)) + 81);
      std::vector<double> j(static_cast<size_t>(L + 2), 0.0);
      j[L + 1] = 0.0;
      j[L] = 1.0;
      for (int n = L; n >= 1; --n) {
        j[n - 1] = (2.0 * n + 1.0) * j[n] / x - j[n + 1];
        const double mag = std::max(std::abs(j[n - 1]), std::abs(j[n]));
        if (mag > 1.0e100) {
          for (int m = n - 1; m <= L + 1; ++m) j[m] *= 1.0e-100;
        }
      }

      const double true_j0 = std::sin(x) / x;
      const double true_j1 = std::sin(x) / (x * x) - std::cos(x) / x;
      double norm = 0.0;
      if (std::abs(j[0]) > std::abs(j[1]) && std::abs(j[0]) > 0.0) {
        norm = true_j0 / j[0];
      } else if (std::abs(j[1]) > 0.0) {
        norm = true_j1 / j[1];
      } else {
        norm = true_j0 / j[0];
      }

      for (size_t io = 0; io < orders.size(); ++io) {
        const int order = orders[io];
        const double jl = norm * j[order];
        const double jlp = norm * j[order + 1];
        grid.value[io][ix] = jl;
        grid.slope[io][ix] = static_cast<double>(order) * jl / x - jlp;
      }
    }
    return grid;
  }

  static std::vector<std::vector<double>> interpolate_uniform_bessel_table(
      const std::vector<std::vector<double>> &grid_table,
      const std::vector<double> &xarg,
      double dx) {
    const size_t n_x = grid_table.front().size();
    std::vector<std::vector<double>> table(
        grid_table.size(), std::vector<double>(xarg.size(), 0.0));

    for (size_t io = 0; io < grid_table.size(); ++io) {
      const std::vector<double> &row = grid_table[io];
      for (size_t ix = 0; ix < xarg.size(); ++ix) {
        const double u = xarg[ix] / dx;
        size_t i = static_cast<size_t>(std::floor(u));
        if (i + 1 >= n_x) i = n_x - 2;
        const double t = u - static_cast<double>(i);

        if (i >= 1 && i + 2 < n_x) {
          const double y0 = row[i - 1];
          const double y1 = row[i];
          const double y2 = row[i + 1];
          const double y3 = row[i + 2];
          const double t2 = t * t;
          const double t3 = t2 * t;
          table[io][ix] =
              0.5 * ((2.0 * y1) + (-y0 + y2) * t +
                     (2.0 * y0 - 5.0 * y1 + 4.0 * y2 - y3) * t2 +
                     (-y0 + 3.0 * y1 - 3.0 * y2 + y3) * t3);
        } else {
          table[io][ix] = (1.0 - t) * row[i] + t * row[i + 1];
        }
      }
    }
    return table;
  }

  static std::vector<std::vector<double>> bessel_class_memory_table(
      const std::vector<int> &orders,
      const std::vector<double> &xarg,
      double hyper_sampling_flat) {
    if (hyper_sampling_flat <= 0.0) {
      throw std::runtime_error("--hyper-sampling-flat must be positive");
    }
    const double dx = 2.0 * PI / hyper_sampling_flat;
    const double x_max = *std::max_element(xarg.begin(), xarg.end()) + 2.0 * dx;
    const size_t n_x = static_cast<size_t>(std::floor(x_max / dx)) + 1;
    std::vector<double> xgrid(n_x, 0.0);
    for (size_t i = 0; i < n_x; ++i) {
      xgrid[i] = dx * static_cast<double>(i);
    }
    BesselHermiteGrid grid = spherical_bessel_hermite_grid(orders, xgrid);
    return interpolate_uniform_bessel_hermite(grid, xarg, dx);
  }

  static bool read_bessel_cache(const std::string &path,
                                const std::vector<int> &expected_orders,
                                const std::vector<double> &expected_xarg,
                                std::vector<std::vector<double>> &table) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;

    char magic[16] = {};
    in.read(magic, sizeof(magic));
    if (!in || std::strncmp(magic, "TFBESSEL1", 9) != 0) return false;

    std::uint64_t n_orders = 0;
    std::uint64_t n_x = 0;
    in.read(reinterpret_cast<char *>(&n_orders), sizeof(n_orders));
    in.read(reinterpret_cast<char *>(&n_x), sizeof(n_x));
    if (!in || n_orders != expected_orders.size() || n_x != expected_xarg.size()) {
      return false;
    }

    std::vector<std::int32_t> orders(n_orders);
    in.read(reinterpret_cast<char *>(orders.data()),
            static_cast<std::streamsize>(orders.size() * sizeof(std::int32_t)));
    if (!in) return false;
    for (size_t i = 0; i < orders.size(); ++i) {
      if (orders[i] != expected_orders[i]) return false;
    }

    std::vector<double> xarg(n_x);
    in.read(reinterpret_cast<char *>(xarg.data()),
            static_cast<std::streamsize>(xarg.size() * sizeof(double)));
    if (!in) return false;
    double max_rel = 0.0;
    double max_abs = 0.0;
    for (size_t i = 0; i < xarg.size(); ++i) {
      const double abs_diff = std::abs(xarg[i] - expected_xarg[i]);
      const double scale = std::max(1.0, std::abs(xarg[i]));
      max_abs = std::max(max_abs, abs_diff);
      max_rel = std::max(max_rel, abs_diff / scale);
      if (abs_diff > 1.0e-10 * scale) {
        std::cerr << "Bessel cache grid mismatch: max_abs_so_far=" << max_abs
                  << " max_scaled_so_far=" << max_rel << "\n";
        return false;
      }
    }

    table.assign(n_orders, std::vector<double>(n_x, 0.0));
    for (size_t i = 0; i < n_orders; ++i) {
      in.read(reinterpret_cast<char *>(table[i].data()),
              static_cast<std::streamsize>(n_x * sizeof(double)));
      if (!in) return false;
    }
    return true;
  }

  static void write_bessel_cache(const std::string &path,
                                 const std::vector<int> &orders,
                                 const std::vector<double> &xarg,
                                 const std::vector<std::vector<double>> &table) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("could not open Bessel cache for writing: " + path);

    char magic[16] = {};
    std::memcpy(magic, "TFBESSEL1", 9);
    out.write(magic, sizeof(magic));

    const std::uint64_t n_orders = static_cast<std::uint64_t>(orders.size());
    const std::uint64_t n_x = static_cast<std::uint64_t>(xarg.size());
    out.write(reinterpret_cast<const char *>(&n_orders), sizeof(n_orders));
    out.write(reinterpret_cast<const char *>(&n_x), sizeof(n_x));

    for (int order : orders) {
      const std::int32_t value = static_cast<std::int32_t>(order);
      out.write(reinterpret_cast<const char *>(&value), sizeof(value));
    }
    out.write(reinterpret_cast<const char *>(xarg.data()),
              static_cast<std::streamsize>(xarg.size() * sizeof(double)));
    for (const auto &row : table) {
      if (row.size() != xarg.size()) throw std::runtime_error("bad Bessel table shape");
      out.write(reinterpret_cast<const char *>(row.data()),
                static_cast<std::streamsize>(row.size() * sizeof(double)));
    }
    if (!out) throw std::runtime_error("failed while writing Bessel cache: " + path);
    std::cerr << "Wrote Bessel cache " << path << "\n";
  }

  static bool read_bessel_x_cache(const std::string &path,
                                  const std::vector<int> &expected_orders,
                                  const std::vector<double> &xarg,
                                  std::vector<std::vector<double>> &table) {
    static std::string cached_path;
    static std::vector<std::int32_t> cached_orders;
    static std::vector<std::vector<double>> cached_rows;
    static std::uint64_t cached_n_x = 0;
    static double cached_dx = 0.0;
    static double cached_x_max = 0.0;

    if (cached_path != path || cached_rows.empty()) {
      std::ifstream in(path, std::ios::binary);
      if (!in) return false;

      char magic[16] = {};
      in.read(magic, sizeof(magic));
      if (!in || std::strncmp(magic, "TFBESSELX1", 10) != 0) return false;

      std::uint64_t n_orders = 0;
      std::uint64_t n_x = 0;
      double dx = 0.0;
      double x_max = 0.0;
      in.read(reinterpret_cast<char *>(&n_orders), sizeof(n_orders));
      in.read(reinterpret_cast<char *>(&n_x), sizeof(n_x));
      in.read(reinterpret_cast<char *>(&dx), sizeof(dx));
      in.read(reinterpret_cast<char *>(&x_max), sizeof(x_max));
      if (!in || n_orders == 0 || n_x < 2 || dx <= 0.0) {
        return false;
      }

      std::vector<std::int32_t> orders(n_orders);
      in.read(reinterpret_cast<char *>(orders.data()),
              static_cast<std::streamsize>(orders.size() * sizeof(std::int32_t)));
      if (!in) return false;

      std::vector<std::vector<double>> rows(
          n_orders, std::vector<double>(static_cast<size_t>(n_x), 0.0));
      for (size_t io = 0; io < n_orders; ++io) {
        in.read(reinterpret_cast<char *>(rows[io].data()),
                static_cast<std::streamsize>(rows[io].size() * sizeof(double)));
        if (!in) return false;
      }

      cached_path = path;
      cached_orders = std::move(orders);
      cached_rows = std::move(rows);
      cached_n_x = n_x;
      cached_dx = dx;
      cached_x_max = x_max;
      std::cerr << "Loaded fixed-x Bessel cache into memory " << path << "\n";
    }

    if (cached_orders.size() != expected_orders.size()) return false;
    for (size_t i = 0; i < cached_orders.size(); ++i) {
      if (cached_orders[i] != expected_orders[i]) return false;
    }
    if (*std::max_element(xarg.begin(), xarg.end()) > cached_x_max) return false;

    table.assign(cached_orders.size(), std::vector<double>(xarg.size(), 0.0));
    for (size_t io = 0; io < cached_orders.size(); ++io) {
      const std::vector<double> &row = cached_rows[io];
      for (size_t ix = 0; ix < xarg.size(); ++ix) {
        const double u = xarg[ix] / cached_dx;
        size_t i0 = static_cast<size_t>(std::floor(u));
        if (i0 + 1 >= cached_n_x) i0 = static_cast<size_t>(cached_n_x - 2);
        const double frac = u - static_cast<double>(i0);
        table[io][ix] = (1.0 - frac) * row[i0] + frac * row[i0 + 1];
      }
    }
    return true;
  }

  static void write_bessel_x_cache(const std::string &path,
                                   const std::vector<int> &orders,
                                   double dx,
                                   double x_max) {
    if (dx <= 0.0 || x_max <= 0.0) {
      throw std::runtime_error("bad fixed-x Bessel cache grid");
    }
    const std::uint64_t n_x =
        static_cast<std::uint64_t>(std::floor(x_max / dx)) + 1;
    std::vector<double> xgrid(static_cast<size_t>(n_x), 0.0);
    for (size_t i = 0; i < xgrid.size(); ++i) {
      xgrid[i] = dx * static_cast<double>(i);
    }
    const double stored_x_max = xgrid.back();
    std::vector<std::vector<double>> table = spherical_bessel_table(orders, xgrid);

    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("could not open fixed-x Bessel cache: " + path);

    char magic[16] = {};
    std::memcpy(magic, "TFBESSELX1", 10);
    out.write(magic, sizeof(magic));

    const std::uint64_t n_orders = static_cast<std::uint64_t>(orders.size());
    out.write(reinterpret_cast<const char *>(&n_orders), sizeof(n_orders));
    out.write(reinterpret_cast<const char *>(&n_x), sizeof(n_x));
    out.write(reinterpret_cast<const char *>(&dx), sizeof(dx));
    out.write(reinterpret_cast<const char *>(&stored_x_max), sizeof(stored_x_max));

    for (int order : orders) {
      const std::int32_t value = static_cast<std::int32_t>(order);
      out.write(reinterpret_cast<const char *>(&value), sizeof(value));
    }
    for (const auto &row : table) {
      out.write(reinterpret_cast<const char *>(row.data()),
                static_cast<std::streamsize>(row.size() * sizeof(double)));
    }
    if (!out) throw std::runtime_error("failed while writing fixed-x Bessel cache: " + path);
    const double mib =
        static_cast<double>(n_orders * n_x * sizeof(double)) / (1024.0 * 1024.0);
    std::cerr << "Wrote fixed-x Bessel cache " << path << " ("
              << n_orders << " orders, " << n_x << " x-points, "
              << mib << " MiB values)\n";
  }
};

Params parse_args(int argc, char **argv) {
  Params p;
  for (int i = 1; i < argc; ++i) {
    std::string key = argv[i];
    auto need_value = [&](const std::string &name) {
      if (i + 1 >= argc) throw std::runtime_error("missing value for " + name);
      return std::string(argv[++i]);
    };
    if (key == "--As") p.As = std::stod(need_value(key));
    else if (key == "--omega-cdm") p.omega_cdm_h2 = std::stod(need_value(key));
    else if (key == "--omega-b") p.omega_b_h2 = std::stod(need_value(key));
    else if (key == "--h") p.h = std::stod(need_value(key));
    else if (key == "--ns") p.ns = std::stod(need_value(key));
    else if (key == "--N-eff" || key == "--n-eff") p.N_eff = std::stod(need_value(key));
    else if (key == "--n-source") p.n_source = std::stoi(need_value(key));
    else if (key == "--ell-grid") p.ell_grid = need_value(key);
    else if (key == "--ell-min") p.ell_min = std::stoi(need_value(key));
    else if (key == "--ell-max") p.ell_max = std::stoi(need_value(key));
    else if (key == "--ell-step") p.ell_step = std::stoi(need_value(key));
    else if (key == "--sampled-output") p.interpolate_output = false;
    else if (key == "--interpolated-output") p.interpolate_output = true;
    else if (key == "--output") p.output = need_value(key);
    else if (key == "--bessel-cache") p.bessel_cache = need_value(key);
    else if (key == "--bessel-x-cache") p.bessel_x_cache = need_value(key);
    else if (key == "--bessel-x-dx") p.bessel_x_dx = std::stod(need_value(key));
    else if (key == "--bessel-x-max") p.bessel_x_max = std::stod(need_value(key));
    else if (key == "--bessel-class-memory") p.bessel_class_memory = true;
    else if (key == "--hyper-sampling-flat") p.hyper_sampling_flat = std::stod(need_value(key));
    else if (key == "--k-samples-per-period") p.k_samples_per_period = std::stod(need_value(key));
    else if (key == "--server") p.server = true;
    else if (key == "--scalar") p.scalar_enabled = true;
    else if (key == "--scalar-m-over-H0") p.scalar_m_over_H0 = std::stod(need_value(key));
    else if (key == "--scalar-f") p.scalar_f = std::stod(need_value(key));
    else if (key == "--scalar-theta-ini") p.scalar_theta_ini = std::stod(need_value(key));
    else if (key == "--scalar-theta-prime-ini") p.scalar_theta_prime_ini = std::stod(need_value(key));
    else if (key == "--scalar-n") p.scalar_n = std::stoi(need_value(key));
    else if (key == "--background-steps") p.background_steps = std::stoi(need_value(key));
    else if (key == "--n-isw") p.n_isw = std::stoi(need_value(key));
    else if (key == "--isw-k-eta0-max") p.isw_k_eta0_max = std::stod(need_value(key));
    else if (key == "--source-mode") p.source_mode = need_value(key);
    else if (key == "--recfast-output") p.recfast_output = need_value(key);
    else if (key == "--recfast-steps") p.recfast_steps = std::stoi(need_value(key));
    else if (key == "--helium-y-p") p.helium_y_p = std::stod(need_value(key));
    else if (key == "--help") {
      std::cout << "Usage: two_fluid_tt_neutrino [--As val] [--omega-cdm val] [--omega-b val]\n"
                   "                    [--h val] [--ns val] [--N-eff val] [--n-source N]\n"
                   "                    [--ell-grid sparse|class]\n"
                   "                    [--ell-min L] [--ell-max L] [--ell-step N]\n"
                   "                    [--interpolated-output|--sampled-output]\n"
                   "                    [--output path]\n"
                   "                    [--bessel-cache path]\n"
                   "                    [--bessel-x-cache path] [--bessel-x-dx dx]\n"
                   "                    [--bessel-x-max xmax]\n"
                   "                    [--bessel-class-memory]\n"
                   "                    [--hyper-sampling-flat samples]\n"
                   "                    [--k-samples-per-period samples]\n"
                   "                    [--scalar]\n"
                   "                    [--scalar-m-over-H0 val] [--scalar-f val]\n"
                   "                    [--scalar-theta-ini val]\n"
                   "                    [--scalar-theta-prime-ini val]\n"
                   "                    [--scalar-n N] [--background-steps N]\n"
                   "                    [--n-isw N] [--isw-k-eta0-max K]\n"
                   "                    [--source-mode full|sw|doppler|isw]\n"
                   "                    [--recfast-output path]\n"
                   "                    [--recfast-steps N] [--helium-y-p Yp]\n"
                   "                    [--server]\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown argument: " + key);
    }
  }
  if (p.n_source < 11) throw std::runtime_error("--n-source must be at least 11");
  if (p.N_eff < 0.0) throw std::runtime_error("--N-eff must be non-negative");
  if (p.scalar_n < 1) throw std::runtime_error("--scalar-n must be at least 1");
  if (p.background_steps < 1000) throw std::runtime_error("--background-steps must be at least 1000");
  if (p.k_samples_per_period <= 0.0) {
    throw std::runtime_error("--k-samples-per-period must be positive");
  }
  if (p.n_isw < 2) throw std::runtime_error("--n-isw must be at least 2");
  if (p.isw_k_eta0_max <= 0.0) throw std::runtime_error("--isw-k-eta0-max must be positive");
  if (p.source_mode != "full" && p.source_mode != "sw" &&
      p.source_mode != "doppler" && p.source_mode != "isw") {
    throw std::runtime_error("--source-mode must be full, sw, doppler, or isw");
  }
  if (p.recfast_steps < 200) throw std::runtime_error("--recfast-steps must be at least 200");
  if (p.helium_y_p < 0.0 || p.helium_y_p > 0.5) {
    throw std::runtime_error("--helium-y-p must be between 0 and 0.5");
  }
  return p;
}

std::vector<double> dense_dell_over_2pi(const TwoFluidModel::Spectrum &spectrum) {
  std::vector<double> ell_x(spectrum.ell.size(), 0.0);
  for (size_t i = 0; i < spectrum.ell.size(); ++i) {
    ell_x[i] = static_cast<double>(spectrum.ell[i]);
  }
  CubicSpline cl_spline(ell_x, spectrum.cl);
  const int ell_min = spectrum.ell.front();
  const int ell_max = spectrum.ell.back();
  std::vector<double> values;
  values.reserve(static_cast<size_t>(ell_max - ell_min + 1));
  for (int ell = ell_min; ell <= ell_max; ++ell) {
    const double l = static_cast<double>(ell);
    const double cl = cl_spline(l);
    values.push_back(l * (l + 1.0) * cl / (2.0 * PI));
  }
  return values;
}

void write_spectrum_file(const Params &p, const TwoFluidModel::Spectrum &spectrum) {
  std::ofstream out(p.output);
  out << "# ell Cl ell_ellplus1_Cl_over_2pi quadrupole_normalized\n";
  const double c2 = spectrum.cl.front();

  if (p.interpolate_output) {
    std::vector<double> ell_x(spectrum.ell.size(), 0.0);
    for (size_t i = 0; i < spectrum.ell.size(); ++i) {
      ell_x[i] = static_cast<double>(spectrum.ell[i]);
    }
    CubicSpline cl_spline(ell_x, spectrum.cl);
    const int ell_min = spectrum.ell.front();
    const int ell_max = spectrum.ell.back();
    for (int ell = ell_min; ell <= ell_max; ++ell) {
      const double l = static_cast<double>(ell);
      const double cl = cl_spline(l);
      const double dell_over_2pi = l * (l + 1.0) * cl / (2.0 * PI);
      const double qnorm = l * (l + 1.0) * cl / (6.0 * c2);
      out << ell << " " << std::setprecision(17) << cl << " "
          << dell_over_2pi << " " << qnorm << "\n";
    }
  } else {
    for (size_t i = 0; i < spectrum.ell.size(); ++i) {
      const double l = static_cast<double>(spectrum.ell[i]);
      const double dell_over_2pi = l * (l + 1.0) * spectrum.cl[i] / (2.0 * PI);
      const double qnorm = l * (l + 1.0) * spectrum.cl[i] / (6.0 * c2);
      out << spectrum.ell[i] << " " << std::setprecision(17) << spectrum.cl[i] << " "
          << dell_over_2pi << " " << qnorm << "\n";
    }
  }
}

void run_server(const Params &base) {
  std::cerr << "SERVER_READY\n";
  std::string line;
  while (std::getline(std::cin, line)) {
    if (line.empty()) continue;
    if (line == "quit" || line == "QUIT") break;

    std::istringstream in(line);
    Params p = base;
    in >> p.As >> p.omega_cdm_h2 >> p.omega_b_h2 >> p.h >> p.ns;
    if (!in) {
      std::cout << "ERR bad_parameter_line\n" << std::flush;
      continue;
    }

    try {
      TwoFluidModel model(p);
      TwoFluidModel::Spectrum spectrum = model.compute_spectrum();
      std::vector<double> values = dense_dell_over_2pi(spectrum);

      std::cout << "OK " << values.size();
      for (double value : values) {
        std::cout << " " << std::setprecision(17) << value;
      }
      std::cout << "\n" << std::flush;
    } catch (const std::exception &e) {
      std::cout << "ERR " << e.what() << "\n" << std::flush;
    }
  }
}

}  // namespace

int main(int argc, char **argv) {
  try {
    Params p = parse_args(argc, argv);
    if (p.server) {
      run_server(p);
      return 0;
    }

    auto t0 = std::chrono::steady_clock::now();
    TwoFluidModel model(p);
    TwoFluidModel::Spectrum spectrum = model.compute_spectrum();
    auto t1 = std::chrono::steady_clock::now();
    const double seconds = std::chrono::duration<double>(t1 - t0).count();
    model.write_scale_summary(std::cerr);
    write_spectrum_file(p, spectrum);
    if (!p.recfast_output.empty()) model.write_recfast_history(p.recfast_output);
    std::cerr << "Computed " << spectrum.ell.size() << " multipoles in " << seconds
              << " s\n";
    std::cerr << "Wrote " << p.output << "\n";
    if (!p.recfast_output.empty()) std::cerr << "Wrote " << p.recfast_output << "\n";
  } catch (const std::exception &e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
