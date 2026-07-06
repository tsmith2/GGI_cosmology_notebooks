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
constexpr double FID_AS = 2.1e-9;
constexpr double FID_NS = 0.965;
constexpr double ELL_PIVOT = 750.0;
constexpr double T_CMB_MICROK = 2.7255e6;

struct Params {
  double As = FID_AS;
  double omega_cdm_h2 = 0.1201;
  double omega_b_h2 = 0.0223;
  // H0 is in km/s/Mpc.  The dimensionless h is computed internally as H0/100.
  double H0 = 67.0;
  double ns = FID_NS;
  int n_source = 301;
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
  bool server = false;
  int n_isw = 96;
  std::string isw_mode = "early-late";
  double late_isw_chi_eff = 3.1;
  int late_isw_ell_max = 120;
  // Gaussian visibility width from the recombination history, sigma_eta / eta_rec.
  double visibility_width_ratio = 0.05386821332646303;
};

struct Background {
  explicit Background(const Params &p)
      : h0(p.H0 / 100.0), H0(p.H0),
        omega_gamma(2.47e-5 / (h0 * h0)),
        omega_cdm(p.omega_cdm_h2 / (h0 * h0)),
        omega_b(p.omega_b_h2 / (h0 * h0)) {
    omega_rad = omega_gamma;
    omega_m = omega_cdm + omega_b;
    omega_lambda = 1.0 - (omega_cdm + omega_b + omega_rad);
  }

  double h0, H0;
  double omega_gamma, omega_rad;
  double omega_cdm, omega_b, omega_m, omega_lambda;
};

using State = std::array<double, 5>;

struct LateISWApprox {
  double growth_decay = 0.0;
  double chi_eff = 0.0;
};

double simpson_integral(const Background &bg, double amax, int n = 65536) {
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
  explicit TwoFluidModel(const Params &params) : p(params), bg(params) {
    arec = 1.0 / (1.0 + zrec(bg.omega_b));
    aeq = bg.omega_gamma / bg.omega_m;
    etatoday = simpson_integral(bg, 1.0);
    etarec = simpson_integral(bg, arec);
    etastar = etarec;
    eta_eq = simpson_integral(bg, aeq);
    tau_r = 1.0 / (std::sqrt(bg.omega_m / arec) * bg.H0 / 2.0);
    alpha = std::sqrt(arec / aeq);
    alpha2 = alpha * alpha;
    yb_prefactor = bg.omega_b / bg.omega_m;
    yc_prefactor = bg.omega_cdm / bg.omega_m;
    xrec = (std::sqrt(alpha2 + 1.0) - 1.0) / alpha;
    xeq = eta_eq / tau_r;
    xs = 0.6 * std::pow(bg.omega_m, 0.25) * std::pow(bg.omega_b, -0.5) *
         std::pow(arec, 0.75) * std::pow(bg.h0, -0.5);
  }

  struct Spectrum {
    std::vector<int> ell;
    std::vector<double> cl;
  };

  Spectrum compute_spectrum() const {
    std::vector<double> source_kappa, sw, dop, phi_rec, radiation_fraction;
    solve_sources(source_kappa, sw, dop, phi_rec, radiation_fraction);

    CubicSpline sw_spline(source_kappa, sw);
    CubicSpline dop_spline(source_kappa, dop);
    CubicSpline phi_spline(source_kappa, phi_rec);
    CubicSpline fr_spline(source_kappa, radiation_fraction);

    const double kmin = tau_r / etatoday;
    const double kmax = tau_r * 1.0e4 / etatoday;
    const double dk = tau_r * 2.0 * PI / (10.0 * etatoday);
    std::vector<double> kappa;
    for (double k = kmin; k <= kmax + 0.5 * dk; k += dk) kappa.push_back(k);

    std::vector<int> ell = ell_list();
    std::vector<double> cl(ell.size(), 0.0);
    const double chi = (etatoday - etastar) / tau_r;
    const bool use_los_isw = p.isw_mode == "los";
    const bool use_seljak_isw = p.isw_mode == "seljak";
    const bool use_early_isw = p.isw_mode == "early" || p.isw_mode == "early-late";
    const bool use_late_isw = p.isw_mode == "early-late";
    const double seljak_isw = use_seljak_isw ? 2.0 * delta_phi() : 0.0;
    const LateISWApprox late_isw = use_late_isw ? late_isw_approximation() : LateISWApprox{};
    std::vector<std::vector<double>> isw_source;
    if (use_los_isw) {
      isw_source = compute_isw_sources(ell, source_kappa);
    }
    // In CMB TT without explicit reionization, this amplitude should be read as
    // A_s.  With reionization included phenomenologically, high-ell TT measures
    // approximately A_s * exp(-2 tau_reio).
    const double amp = 4.0 * PI * p.As * T_CMB_MICROK * T_CMB_MICROK;
    const double cs_rec = sound_speed_rec();
    const double sigma_x = p.visibility_width_ratio * xrec;
    const double finite_width_damping = cs_rec * cs_rec * sigma_x * sigma_x;

    std::vector<double> xarg(kappa.size()), damping(kappa.size()), swisw(kappa.size()),
        dopv(kappa.size());
    for (size_t i = 0; i < kappa.size(); ++i) {
      xarg[i] = kappa[i] * chi;
      damping[i] = std::exp(
          -kappa[i] * kappa[i] * (2.0 * xs * xs + finite_width_damping));
      swisw[i] = sw_spline(kappa[i]) + seljak_isw;
      dopv[i] = dop_spline(kappa[i]);
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

    std::vector<std::vector<double>> late_bessel;
    std::vector<int> late_bessel_row(ell.size(), -1);
    if (use_late_isw) {
      std::vector<double> late_xarg(kappa.size(), 0.0);
      for (size_t i = 0; i < kappa.size(); ++i) {
        late_xarg[i] = kappa[i] * late_isw.chi_eff;
      }
      std::vector<int> late_orders;
      for (size_t j = 0; j < ell.size(); ++j) {
        if (p.late_isw_ell_max < 0 || ell[j] <= p.late_isw_ell_max) {
          late_bessel_row[j] = static_cast<int>(late_orders.size());
          late_orders.push_back(ell[j]);
        }
      }
      if (!late_orders.empty()) {
        if (!p.bessel_x_cache.empty()) {
          if (!read_bessel_x_cache(p.bessel_x_cache, late_orders, late_xarg,
                                   late_bessel)) {
            late_bessel = spherical_bessel_table(late_orders, late_xarg);
          }
        } else {
          late_bessel = spherical_bessel_table(late_orders, late_xarg);
        }
      }
    }

    for (size_t j = 0; j < ell.size(); ++j) {
      const int l = ell[j];
      CubicSpline isw_spline;
      if (use_los_isw) {
        isw_spline = CubicSpline(source_kappa, isw_source[j]);
      }
      const size_t idx_lm = order_index(l - 1);
      const size_t idx_l = order_index(l);
      const size_t idx_lp = order_index(l + 1);
      std::vector<double> integrand(kappa.size(), 0.0);
      for (size_t i = 0; i < kappa.size(); ++i) {
        const double x = xarg[i];
        const double jl = bessel[idx_l][i];
        const double jlm = bessel[idx_lm][i];
        const double jlp = bessel[idx_lp][i];
        const double deriv = -jl / (2.0 * x) + 0.5 * (jlm - jlp);
        const double primary_src = swisw[i] * jl + dopv[i] * deriv;
        double src = primary_src;
        if (use_los_isw || use_early_isw) {
          double isw = 0.0;
          if (use_los_isw) {
            isw = isw_spline(kappa[i]);
          } else {
            const double phi_star = phi_spline(kappa[i]);
            const double fr = fr_spline(kappa[i]);
            isw = -2.0 * phi_star * fr *
                  (jl + (kappa[i] * xrec / 3.0) * deriv);
            if (use_late_isw && late_bessel_row[j] >= 0) {
              const double fm = 1.0 - fr;
              const size_t late_row = static_cast<size_t>(late_bessel_row[j]);
              isw += 2.0 * phi_star * fm * late_isw.growth_decay *
                     late_bessel[late_row][i];
            }
          }
          src = std::sqrt(damping[i]) * primary_src + isw;
          integrand[i] = src * src / kappa[i];
        } else {
          integrand[i] = damping[i] * src * src / kappa[i];
        }
      }
      double val = 0.0;
      for (size_t i = 0; i + 1 < kappa.size(); ++i) {
        val += 0.5 * (integrand[i] + integrand[i + 1]) *
               (kappa[i + 1] - kappa[i]);
      }
      const double tilt = std::pow(static_cast<double>(l) / ELL_PIVOT, p.ns - FID_NS);
      cl[j] = amp * tilt * val;
    }
    return Spectrum{ell, cl};
  }

 private:
  Params p;
  Background bg;
  double arec, aeq, etatoday, etarec, etastar, eta_eq, tau_r;
  double alpha, alpha2, yb_prefactor, yc_prefactor, xrec, xeq, xs;

  static double zrec(double omega_b) {
    return 1000.0 * std::pow(omega_b, -0.027 / (1.0 + 0.11 * std::log(omega_b)));
  }

  double y(double x) const { return alpha2 * x * x + 2.0 * alpha * x; }

  double delta_phi() const {
    const double yrec = y(xrec);
    return (2.0 - 8.0 / yrec + 16.0 * xrec / (yrec * yrec * yrec)) /
           (10.0 * yrec);
  }

  double sound_speed_rec() const {
    const double yb_rec = yb_prefactor * y(xrec);
    return 1.0 / std::sqrt(3.0 * (1.0 + 0.75 * yb_rec));
  }

  double e_of_a(double a) const {
    return std::sqrt(bg.omega_rad / (a * a * a * a) +
                     bg.omega_m / (a * a * a) + bg.omega_lambda);
  }

  double omega_m_of_a(double a) const {
    const double e2 = bg.omega_rad / (a * a * a * a) +
                      bg.omega_m / (a * a * a) + bg.omega_lambda;
    return (bg.omega_m / (a * a * a)) / e2;
  }

  LateISWApprox late_isw_approximation() const {
    const int n = 4096;
    std::vector<double> a(n + 1, 0.0), eta(n + 1, 0.0), ln_d(n + 1, 0.0);
    const double log_amin = std::log(arec);
    for (int i = 0; i <= n; ++i) {
      const double t = static_cast<double>(i) / static_cast<double>(n);
      a[i] = std::exp(log_amin * (1.0 - t));
    }

    eta[0] = etarec;
    for (int i = 0; i < n; ++i) {
      const double amid = std::sqrt(a[i] * a[i + 1]);
      const double deta_da = 1.0 / (bg.H0 * amid * amid * e_of_a(amid));
      eta[i + 1] = eta[i] + deta_da * (a[i + 1] - a[i]);
    }

    // D(1)=1 and d ln D / d ln a = Omega_m(a)^0.55.
    ln_d[n] = 0.0;
    for (int i = n - 1; i >= 0; --i) {
      const double amid = std::sqrt(a[i] * a[i + 1]);
      const double f = std::pow(omega_m_of_a(amid), 0.55);
      ln_d[i] = ln_d[i + 1] - f * (std::log(a[i + 1]) - std::log(a[i]));
    }

    const double g_star = std::exp(ln_d[0]) / a[0];
    std::vector<double> q(n + 1, 0.0);
    for (int i = 0; i <= n; ++i) {
      q[i] = (std::exp(ln_d[i]) / a[i]) / g_star;
    }

    const double dq_total = q[n] - q[0];
    double eta_weighted = 0.0;
    for (int i = 0; i < n; ++i) {
      const double dq = q[i + 1] - q[i];
      eta_weighted += 0.5 * (eta[i] + eta[i + 1]) * dq;
    }
    const double eta_eff = std::abs(dq_total) > 1.0e-14
                               ? eta_weighted / dq_total
                               : 0.5 * (etarec + etatoday);
    double chi_eff = (etatoday - eta_eff) / tau_r;
    if (p.late_isw_chi_eff >= 0.0) chi_eff = p.late_isw_chi_eff;
    return LateISWApprox{dq_total, chi_eff};
  }

  State rhs(double x, const State &u, double kappa) const {
    const double delta_gamma = u[0];
    const double v_gamma = u[2];
    const double v_c = u[3];
    const double phi = u[4];
    const double yy = y(x);
    const double eta = 2.0 * alpha * (alpha * x + 1.0) / yy;
    const double yb = yb_prefactor * yy;
    const double yc = yc_prefactor * yy;
    const double phi_p =
        -eta * phi +
        3.0 * eta * eta *
            (v_gamma * (4.0 / 3.0 + yy - yc) + yc * v_c) /
            (2.0 * (1.0 + yy) * kappa);
    State du{};
    du[0] = -(4.0 / 3.0) * kappa * v_gamma + 4.0 * phi_p;
    du[1] = -kappa * v_c + 3.0 * phi_p;
    du[2] = (-eta * yb * v_gamma + kappa * delta_gamma / 3.0) /
                (4.0 / 3.0 + yb) +
            kappa * phi;
    du[3] = -eta * v_c + kappa * phi;
    du[4] = phi_p;
    return du;
  }

  std::pair<double, State> initial_conditions(double kappa) const {
    const double xi = std::min(1.0e-4 / kappa, xeq / 10000.0);
    // Ma & Bertschinger Eq. (98), conformal Newtonian growing mode.
    // CLASS uses the primordial-curvature convention corresponding to C=1/2.
    const double c_mb = 0.5;
    const double phi_i = 20.0 * c_mb / 15.0;
    const double yi = y(xi);
    const double eta_i = 2.0 * alpha * (alpha * xi + 1.0) / yi;
    const double dg_i = -2.0 * phi_i * (1.0 + 3.0 * yi / 16.0);
    const double dc_i = 0.75 * dg_i;
    const double vg_i =
        (-kappa / eta_i) *
        (dg_i / 4.0 +
         2.0 * kappa * kappa * (1.0 + yi) * phi_i /
             (9.0 * eta_i * eta_i * (4.0 / 3.0 + yi)));
    return {xi, State{dg_i, dc_i, vg_i, vg_i, phi_i}};
  }

  State rk45_step(double x, const State &u, double h, double kappa,
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
      for (size_t i = 0; i < out.size(); ++i) {
        for (const auto &term : terms) out[i] += hstep * term.first * term.second[i];
      }
      return out;
    };

    const State k1 = rhs(x, u, kappa);
    const State k2 = rhs(x + c2 * h, add(u, h, {{a21, k1}}), kappa);
    const State k3 = rhs(x + c3 * h, add(u, h, {{a31, k1}, {a32, k2}}), kappa);
    const State k4 =
        rhs(x + c4 * h, add(u, h, {{a41, k1}, {a42, k2}, {a43, k3}}), kappa);
    const State k5 = rhs(x + c5 * h,
                         add(u, h, {{a51, k1}, {a52, k2}, {a53, k3}, {a54, k4}}),
                         kappa);
    const State k6 = rhs(x + h,
                         add(u, h, {{a61, k1}, {a62, k2}, {a63, k3}, {a64, k4},
                                    {a65, k5}}),
                         kappa);
    const State k7 = rhs(x + h,
                         add(u, h, {{a71, k1}, {a73, k3}, {a74, k4}, {a75, k5},
                                    {a76, k6}}),
                         kappa);

    State y5{}, y4{};
    for (size_t i = 0; i < y5.size(); ++i) {
      y5[i] = u[i] + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] +
                          b5 * k5[i] + b6 * k6[i]);
      y4[i] = u[i] + h * (bs1 * k1[i] + bs3 * k3[i] + bs4 * k4[i] +
                          bs5 * k5[i] + bs6 * k6[i] + bs7 * k7[i]);
      err[i] = y5[i] - y4[i];
    }
    return y5;
  }

  State integrate_source(double kappa) const {
    auto init = initial_conditions(kappa);
    double x = init.first;
    State u = init.second;
    double h = std::min(1.0e-3, (xrec - x) / 20.0);
    const double rtol = 2.0e-8;
    const double atol = 2.0e-10;
    int steps = 0;
    while (x < xrec && steps < 200000) {
      if (x + h > xrec) h = xrec - x;
      State err{};
      State trial = rk45_step(x, u, h, kappa, err);
      double err_norm = 0.0;
      for (size_t i = 0; i < u.size(); ++i) {
        const double scale = atol + rtol * std::max(std::abs(u[i]), std::abs(trial[i]));
        err_norm = std::max(err_norm, std::abs(err[i]) / scale);
      }
      if (err_norm <= 1.0) {
        x += h;
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

  std::vector<double> phi_prime_history(double kappa,
                                        const std::vector<double> &x_nodes) const {
    auto init = initial_conditions(kappa);
    double x = init.first;
    State u = init.second;
    double h = std::min(1.0e-3, (x_nodes.front() - x) / 20.0);
    const double rtol = 2.0e-8;
    const double atol = 2.0e-10;
    int steps = 0;
    std::vector<double> history(x_nodes.size(), 0.0);

    for (size_t target_index = 0; target_index < x_nodes.size(); ++target_index) {
      const double target = x_nodes[target_index];
      while (x < target && steps < 400000) {
        if (x + h > target) h = target - x;
        State err{};
        State trial = rk45_step(x, u, h, kappa, err);
        double err_norm = 0.0;
        for (size_t i = 0; i < u.size(); ++i) {
          const double scale = atol + rtol * std::max(std::abs(u[i]), std::abs(trial[i]));
          err_norm = std::max(err_norm, std::abs(err[i]) / scale);
        }
        if (err_norm <= 1.0) {
          x += h;
          u = trial;
        }
        const double factor =
            (err_norm == 0.0) ? 5.0
                              : std::clamp(0.9 * std::pow(err_norm, -0.2), 0.2, 5.0);
        h *= factor;
        ++steps;
      }
      if (steps >= 400000) throw std::runtime_error("ISW ODE integration did not converge");
      history[target_index] = rhs(target, u, kappa)[4];
      h = std::min(h, std::max(1.0e-6, (target_index + 1 < x_nodes.size()
                                            ? x_nodes[target_index + 1] - target
                                            : h)));
    }
    return history;
  }

  std::vector<std::vector<double>> compute_isw_sources(
      const std::vector<int> &ell,
      const std::vector<double> &source_kappa) const {
    const double x_today = etatoday / tau_r;
    std::vector<double> x_nodes(static_cast<size_t>(std::max(2, p.n_isw)), 0.0);
    for (size_t i = 0; i < x_nodes.size(); ++i) {
      const double t = static_cast<double>(i) / static_cast<double>(x_nodes.size() - 1);
      x_nodes[i] = xrec + t * (x_today - xrec);
    }

    std::vector<std::vector<double>> phi_p(
        source_kappa.size(), std::vector<double>(x_nodes.size(), 0.0));
    for (size_t ik = 0; ik < source_kappa.size(); ++ik) {
      phi_p[ik] = phi_prime_history(source_kappa[ik], x_nodes);
    }

    std::vector<double> xarg;
    xarg.reserve(source_kappa.size() * x_nodes.size());
    for (double kappa : source_kappa) {
      for (double xnode : x_nodes) {
        xarg.push_back(kappa * (x_today - xnode));
      }
    }

    std::vector<int> orders(ell.begin(), ell.end());
    std::vector<std::vector<double>> bessel = spherical_bessel_table(orders, xarg);
    std::vector<std::vector<double>> isw(
        ell.size(), std::vector<double>(source_kappa.size(), 0.0));

    for (size_t il = 0; il < ell.size(); ++il) {
      for (size_t ik = 0; ik < source_kappa.size(); ++ik) {
        double val = 0.0;
        const size_t offset = ik * x_nodes.size();
        for (size_t ix = 0; ix + 1 < x_nodes.size(); ++ix) {
          const double f0 = 2.0 * phi_p[ik][ix] * bessel[il][offset + ix];
          const double f1 = 2.0 * phi_p[ik][ix + 1] * bessel[il][offset + ix + 1];
          val += 0.5 * (f0 + f1) * (x_nodes[ix + 1] - x_nodes[ix]);
        }
        isw[il][ik] = val;
      }
    }
    return isw;
  }

  void solve_sources(std::vector<double> &kappas, std::vector<double> &sw,
                     std::vector<double> &dop, std::vector<double> &phi_rec,
                     std::vector<double> &radiation_fraction) const {
    kappas.resize(p.n_source);
    sw.resize(p.n_source);
    dop.resize(p.n_source);
    phi_rec.resize(p.n_source);
    radiation_fraction.resize(p.n_source);

    for (int i = 0; i < p.n_source; ++i) {
      const double grid = 4.0 * static_cast<double>(i) / static_cast<double>(p.n_source - 1);
      kappas[i] = tau_r * std::pow(10.0, grid) / etatoday;
      State u = integrate_source(kappas[i]);
      sw[i] = u[4] + u[0] / 4.0;
      dop[i] = u[2];
      phi_rec[i] = u[4];

      const double yy = y(xrec);
      const double radiation_source = u[0];
      const double baryon_delta = 0.75 * u[0];
      const double matter_source =
          yb_prefactor * yy * baryon_delta + yc_prefactor * yy * u[1];
      const double total_source = radiation_source + matter_source;
      radiation_fraction[i] =
          std::abs(total_source) > 1.0e-30 ? radiation_source / total_source : 0.0;
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
    else if (key == "--H0" || key == "--H_0") p.H0 = std::stod(need_value(key));
    else if (key == "--h") p.H0 = 100.0 * std::stod(need_value(key));
    else if (key == "--ns") p.ns = std::stod(need_value(key));
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
    else if (key == "--n-isw") p.n_isw = std::stoi(need_value(key));
    else if (key == "--isw-mode") p.isw_mode = need_value(key);
    else if (key == "--late-isw-chi-eff") p.late_isw_chi_eff = std::stod(need_value(key));
    else if (key == "--late-isw-ell-max") p.late_isw_ell_max = std::stoi(need_value(key));
    else if (key == "--visibility-width-ratio") p.visibility_width_ratio = std::stod(need_value(key));
    else if (key == "--server") p.server = true;
    else if (key == "--help") {
      std::cout << "Usage: two_fluid_tt [--As val] [--omega-cdm val] [--omega-b val]\n"
                   "                    [--H0 val] [--ns val] [--n-source N]\n"
                   "                    [--ell-grid sparse|class]\n"
                   "                    [--ell-min L] [--ell-max L] [--ell-step N]\n"
                   "                    [--interpolated-output|--sampled-output]\n"
                   "                    [--output path]\n"
                   "                    [--bessel-cache path]\n"
                   "                    [--bessel-x-cache path] [--bessel-x-dx dx]\n"
                   "                    [--bessel-x-max xmax]\n"
                   "                    [--bessel-class-memory]\n"
                   "                    [--hyper-sampling-flat samples]\n"
                   "                    [--n-isw N]\n"
                   "                    [--isw-mode none|seljak|early|early-late|los]\n"
                   "                    [--late-isw-chi-eff chi]\n"
                   "                    [--late-isw-ell-max L]\n"
                   "                    [--visibility-width-ratio sigma_eta_over_eta_rec]\n"
                   "                    [--server]\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown argument: " + key);
    }
  }
  if (p.n_source < 11) throw std::runtime_error("--n-source must be at least 11");
  if (p.n_isw < 2) throw std::runtime_error("--n-isw must be at least 2");
  if (p.isw_mode != "none" && p.isw_mode != "seljak" &&
      p.isw_mode != "early" && p.isw_mode != "early-late" &&
      p.isw_mode != "los") {
    throw std::runtime_error("--isw-mode must be none, seljak, early, early-late, or los");
  }
  if (p.visibility_width_ratio < 0.0) {
    throw std::runtime_error("--visibility-width-ratio must be non-negative");
  }
  if (p.late_isw_chi_eff < -1.0) {
    throw std::runtime_error("--late-isw-chi-eff must be non-negative, or -1 for automatic");
  }
  if (p.late_isw_ell_max != -1 && p.late_isw_ell_max < 2) {
    throw std::runtime_error("--late-isw-ell-max must be >= 2, or -1 for no cutoff");
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
    in >> p.As >> p.omega_cdm_h2 >> p.omega_b_h2 >> p.H0 >> p.ns;
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
    write_spectrum_file(p, spectrum);
    std::cerr << "Computed " << spectrum.ell.size() << " multipoles in " << seconds
              << " s\n";
    std::cerr << "Wrote " << p.output << "\n";
  } catch (const std::exception &e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
