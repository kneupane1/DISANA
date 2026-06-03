#ifndef DISANAMATH_H
#define DISANAMATH_H

// ROOT + std
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TVector2.h>
#include <TVector3.h>

#include <ROOT/RDataFrame.hxx>
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <cmath>
#include <limits>

double NormalizeHCorrByMiddle(TH1* hCorr)
{
    if (!hCorr) return 1.0;

    const int nb = hCorr->GetNbinsX();
    if (nb <= 0) return 1.0;

    double mid = 1.0;

    if (nb % 2 == 1) {
        const int bmid = (nb + 1) / 2;     // e.g. nb=5 -> 3
        mid = hCorr->GetBinContent(bmid);
    } else {
        const int b1 = nb / 2;             // e.g. nb=6 -> 3
        const int b2 = b1 + 1;             // -> 4
        const double v1 = hCorr->GetBinContent(b1);
        const double v2 = hCorr->GetBinContent(b2);
        mid = 0.5 * (v1 + v2);
    }

    if (!std::isfinite(mid) || mid == 0.0) {
        // std::cerr << "Warning: mid is invalid (" << mid << "), skip scaling.\n";
        return mid;
    }

    hCorr->Scale(1.0 / mid);
    return mid;
}


// -----------------------------------------------------------------------------
// Constants
// -----------------------------------------------------------------------------
constexpr double pi = 3.14159265358979323846;
const double m_e = 0.000511;  // GeV
// const double m_e = 0.0;       // GeV
const double m_p = 0.938272;  // GeV
// const double m_p = 0.938;       // GeV
const double m_kMinus = 0.493677;  // GeV
const double m_kPlus = 0.493677;   // GeV
const double m_phi = 1.019461;     // GeV

struct EqualStatBinningResult {
  std::vector<double> q2Edges;      // x-axis (Q^2)
  std::vector<double> tprimeEdges;  // y-axis (t')
};

struct DVCSWeightInput {
  int pho_det_region = -1;
  int pro_det_region = -1;
  double recel_p = 0.0;
  double recel_theta = 0.0;
  double recel_phi = 0.0;
  double recpho_p = 0.0;
  double recpho_theta = 0.0;
  double recpho_phi = 0.0;
  double recpro_p = 0.0;
  double recpro_theta = 0.0;
  double recpro_phi = 0.0;
};

using DVCSWeightFunction = std::function<double(const DVCSWeightInput&)>;

// -----------------------------------------------------------------------------
// (p,theta,phi) helpers
// -----------------------------------------------------------------------------
namespace {
TVector3 SphericalToCartesian(double p, double theta, double phi) {
  const double px = p * std::sin(theta) * std::cos(phi);
  const double py = p * std::sin(theta) * std::sin(phi);
  const double pz = p * std::cos(theta);
  return TVector3(px, py, pz);
}
TLorentzVector Build4Vector(double p, double theta, double phi, double mass) {
  TVector3 v = SphericalToCartesian(p, theta, phi);
  const double E = std::sqrt(v.Mag2() + mass * mass);
  return TLorentzVector(v, E);
}
}  // namespace

// -----------------------------------------------------------------------------
// Bin manager
// -----------------------------------------------------------------------------
class BinManager {
 public:
  BinManager() {
    q2_bins_ = {1.0, 2.0, 4.0, 6.0};
    t_bins_ = {0.1, 0.3, 0.6, 1.0};
    tprime_bins_ = {0.0, 10.0};
    xb_bins_ = {0.1, 0.2, 0.4, 0.6};
    W_bins_ = {0.1, 10.0};
    cos_thetaKK_bins_ = {-1.0, 1.0};
    trento_phi_bins_ = {-180.0, 180.0};
    z_phi_bins_      = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};  // example, change as you like
  }

  const std::vector<double> &GetQ2Bins() const { return q2_bins_; }
  const std::vector<double> &GetTBins() const { return t_bins_; }
  const std::vector<double> &GetXBBins() const { return xb_bins_; }
  const std::vector<double> &GetWBins() const { return W_bins_; }
  const std::vector<double> &GetTprimeBins() const { return tprime_bins_; }
  const std::vector<double>& GetCosThetaKKBins() const { return cos_thetaKK_bins_; }
  const std::vector<double>& GetTrentoPhiBins() const { return trento_phi_bins_; }
  const std::vector<double>& GetZPhiBins() const { return z_phi_bins_; }
  void SetQ2Bins(const std::vector<double> &v) { q2_bins_ = v; }
  void SetTBins(const std::vector<double> &v) { t_bins_ = v; }
  void SetXBBins(const std::vector<double> &v) { xb_bins_ = v; }
  void SetWBins(const std::vector<double> &v) { W_bins_ = v; }
  void SetTprimeBins(const std::vector<double> &v) { tprime_bins_ = v; }
  void SetCosThetaKKBins(const std::vector<double>& v) { cos_thetaKK_bins_ = v; }
  void SetTrentoPhiBins(const std::vector<double>& v) { trento_phi_bins_ = v; }
  void SetZPhiBins(const std::vector<double>& v) { z_phi_bins_ = v; }

  std::vector<double> ComputeAxisEdges2D(const TH2 *h, bool alongX, int nDesiredBins) {
    std::vector<double> edges;

    if (!h) {
      std::cerr << "[ComputeAxisEdges2D] ERROR: null histogram pointer\n";
      return edges;
    }
    if (nDesiredBins <= 0) {
      std::cerr << "[ComputeAxisEdges2D] ERROR: nDesiredBins <= 0 (" << nDesiredBins << ")\n";
      return edges;
    }

    const int nx = h->GetNbinsX();
    const int ny = h->GetNbinsY();

    const int nAxisBins = alongX ? nx : ny;
    std::vector<double> proj(nAxisBins, 0.0);

    // Build 1D projection along chosen axis
    if (alongX) {
      for (int ix = 1; ix <= nx; ++ix) {
        double sum = 0.0;
        for (int iy = 1; iy <= ny; ++iy) {
          sum += h->GetBinContent(ix, iy);
        }
        proj[ix - 1] = sum;
      }
    } else {
      for (int iy = 1; iy <= ny; ++iy) {
        double sum = 0.0;
        for (int ix = 1; ix <= nx; ++ix) {
          sum += h->GetBinContent(ix, iy);
        }
        proj[iy - 1] = sum;
      }
    }

    double total = 0.0;
    for (double v : proj) total += v;

    std::cout << "[ComputeAxisEdges2D] alongX=" << alongX << " nDesiredBins=" << nDesiredBins << " nAxisBins=" << nAxisBins << " total=" << total << std::endl;

    const TAxis *ax = alongX ? h->GetXaxis() : h->GetYaxis();

    if (total <= 0.0) {
      std::cerr << "[ComputeAxisEdges2D] WARNING: histogram appears empty; "
                << "returning original axis edges." << std::endl;

      edges.reserve(ax->GetNbins() + 1);
      edges.push_back(ax->GetBinLowEdge(1));
      for (int i = 1; i <= ax->GetNbins(); ++i) {
        edges.push_back(ax->GetBinUpEdge(i));
      }
      return edges;
    }

    edges.reserve(nDesiredBins + 1);
    edges.push_back(ax->GetBinLowEdge(1));  // first edge

    double cumulative = 0.0;
    double target = total / nDesiredBins;  // target entries per bin

    for (int i = 0; i < nAxisBins; ++i) {
      cumulative += proj[i];

      // Whenever we pass a multiple of target, define a new edge
      while (cumulative >= target && static_cast<int>(edges.size()) < nDesiredBins) {
        double upper = ax->GetBinUpEdge(i + 1);
        edges.push_back(upper);

        // Next target is (number_of_edges_so_far) * total / nDesiredBins
        target = total * static_cast<double>(edges.size()) / nDesiredBins;
      }
    }

    // Ensure last edge at axis maximum
    double lastEdge = ax->GetBinUpEdge(ax->GetNbins());
    if (edges.empty() || std::fabs(edges.back() - lastEdge) > 1e-12) {
      edges.push_back(lastEdge);
    }

    // Remove duplicates (can happen if some regions are empty)
    edges.erase(std::unique(edges.begin(), edges.end(), [](double a, double b) { return std::fabs(a - b) < 1e-10; }), edges.end());

    return edges;
  }
  // Public interface: build both Q^2 and t' edges
  EqualStatBinningResult MakeEqualStatBinning(const TH2 *hQ2t, int nQ2Bins, int nTprimeBins) {
    EqualStatBinningResult res;
    res.tprimeEdges = ComputeAxisEdges2D(hQ2t, /*alongX=*/true, nTprimeBins);  // x = t'
    res.q2Edges = ComputeAxisEdges2D(hQ2t, /*alongX=*/false, nQ2Bins);         // y = Q^2
    return res;
  }

  void DrawQ2TprimeWithGrid(TH2 *h, const std::vector<double> &q2Edges, const std::vector<double> &tprimeEdges, const std::string &canvasName, const std::string &outFileName) {
    if (!h) return;

    TCanvas *c = new TCanvas(canvasName.c_str(), canvasName.c_str(), 900, 700);
    c->SetRightMargin(0.15);
    c->SetBottomMargin(0.12);
    c->SetLeftMargin(0.12);

    gPad->SetRightMargin(0.16);
    h->SetStats(0);
    h->SetTitle("");
    h->GetYaxis()->SetNoExponent(true);
    h->GetYaxis()->SetLabelFont(42);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetYaxis()->SetTitleOffset(1.0);
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetNdivisions(410);
    h->GetXaxis()->SetTitleSize(0.065);
    h->GetXaxis()->SetLabelFont(42);
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitleOffset(0.9);
    h->GetXaxis()->SetNdivisions(205);
    h->GetZaxis()->SetNdivisions(410);
    h->GetZaxis()->SetLabelSize(0.06);
    h->GetZaxis()->SetTitleOffset(1.5);
    h->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h->DrawCopy("COLZ");
    gPad->SetLogz();  // optional

    double ymin = h->GetYaxis()->GetBinLowEdge(1);
    double ymax = h->GetYaxis()->GetBinUpEdge(h->GetYaxis()->GetNbins());
    double xmin = h->GetXaxis()->GetBinLowEdge(1);
    double xmax = h->GetXaxis()->GetBinUpEdge(h->GetXaxis()->GetNbins());

    // Vertical lines at Q^2 edges
    for (double q : tprimeEdges) {
      TLine *lv = new TLine(q, ymin, q, ymax);
      lv->SetLineColor(kRed + 1);
      lv->SetLineWidth(2);
      lv->SetLineStyle(2);
      lv->Draw("same");
    }

    // Horizontal lines at t' edges
    for (double t : q2Edges) {
      TLine *lh = new TLine(xmin, t, xmax, t);
      lh->SetLineColor(kRed + 1);
      lh->SetLineWidth(2);
      lh->SetLineStyle(2);
      lh->Draw("same");
    }

    c->Update();
    c->SaveAs(outFileName.c_str());
  };

 private:
  std::vector<double> q2_bins_, t_bins_, xb_bins_, W_bins_, tprime_bins_, cos_thetaKK_bins_,trento_phi_bins_,  z_phi_bins_;  
};

// -----------------------------------------------------------------------------
// DISANAMath
// -----------------------------------------------------------------------------
struct Pi0Tag {};
class DISANAMath {
 private:
  // Kinematics
  double Q2_{}, xB_{}, t_{}, phi_deg_{}, W_{}, nu_{}, y_{}, cosTheta_KK_{}, cosPhi_KK_{}, z_phi_{};
  double xi_{};
  double costheta_SDHEP_{};
  double phi_SDHEP_deg_{};

  // Exclusivity
  double mx2_ep_{};
  double emiss_{};
  double ptmiss_{};
  double mx2_epg_{};
  double mx2_epKpKm_{};
  double delta_phi_{};
  double theta_gg_{};
  double theta_gphi_{};
  double mx2_egamma_{};
  double mx2_eKp_{};
  double mx2_eKpKm_{};
  double mx2_epKp_{-1.};
  double mx2_epKm_{};
  double Theta_e_gamma_{};
  double Theta_e_phimeson_{};
  double Cone_Kp_{};
  double Cone_Km_{};
  double Cone_p_{};
  double tmin_{-999.0};
  double coplanarity_had_normals_deg_{};

  double Mass_pi0_{};
  double Mx2_eppi0_{};
  double emiss_pi0_{};
  double Mx2_ep_pi0_{};
  double Mx2_epi0_{};
  double ptmiss_pi0_{};
  double theta_pi0pi0_{};
  double delta_phi_pi0_{};
  double theta_epho1_{};
  double theta_epho2_{};
  double theta_pho1pho2_{};

  double DeltaE_{};

  double EpsContract(const TLorentzVector& a,
                   const TLorentzVector& b,
                   const TLorentzVector& c,
                   const TLorentzVector& d)
  {
    TMatrixD M(4,4);
    M(0,0)=a.E();  M(0,1)=a.Px();  M(0,2)=a.Py();  M(0,3)=a.Pz();
    M(1,0)=b.E();  M(1,1)=b.Px();  M(1,2)=b.Py();  M(1,3)=b.Pz();
    M(2,0)=c.E();  M(2,1)=c.Px();  M(2,2)=c.Py();  M(2,3)=c.Pz();
    M(3,0)=d.E();  M(3,1)=d.Px();  M(3,2)=d.Py();  M(3,3)=d.Pz();

    return M.Determinant();
  }
  // Inelasticity / radiative-photon quantities
  double vm_cut_{-999.0};      // maximal inelasticity  [GeV^2]
  double Egamma_star_{-999.0}; // radiative-photon energy in CM(X) [GeV]

 public:
  DISANAMath() = default;

  // DVCS-style (photon provided as p,theta,phi)
  DISANAMath(double e_in_E, double e_out_p, double e_out_theta, double e_out_phi, double p_out_p, double p_out_theta, double p_out_phi, double g_p, double g_theta, double g_phi) {
    TLorentzVector e_in(0, 0, e_in_E, e_in_E);
    TLorentzVector e_out = Build4Vector(e_out_p, e_out_theta, e_out_phi, m_e);
    TLorentzVector p_in(0, 0, 0, m_p);
    TLorentzVector p_out = Build4Vector(p_out_p, p_out_theta, p_out_phi, m_p);
    TLorentzVector photon = Build4Vector(g_p, g_theta, g_phi, 0.0);
    ComputeKinematics(e_in, e_out, p_in, p_out, photon);
  }

  // DVPi0-style (photon provided as p,theta,phi)
  DISANAMath(Pi0Tag, double e_in_E, double e_out_p, double e_out_theta, double e_out_phi, double p_out_p, double p_out_theta, double p_out_phi, double g_p, double g_theta,
             double g_phi, double g2_p, double g2_theta, double g2_phi) {
    TLorentzVector e_in(0, 0, e_in_E, e_in_E);
    TLorentzVector e_out = Build4Vector(e_out_p, e_out_theta, e_out_phi, m_e);
    TLorentzVector p_in(0, 0, 0, m_p);
    TLorentzVector p_out = Build4Vector(p_out_p, p_out_theta, p_out_phi, m_p);
    TLorentzVector photon = Build4Vector(g_p, g_theta, g_phi, 0.0);
    TLorentzVector photon2 = Build4Vector(g2_p, g2_theta, g2_phi, 0.0);
    ComputePi0Kinematics(e_in, e_out, p_in, p_out, photon, photon2);
  }

  // φ→K+K− analysis (both kaons given)
  DISANAMath(double e_in_E, double e_out_p, double e_out_theta, double e_out_phi, double p_out_p, double p_out_theta, double p_out_phi, double kMinus_p, double kMinus_theta,
             double kMinus_phi, double kPlus_p, double kPlus_theta, double kPlus_phi) {
    TLorentzVector e_in(0, 0, e_in_E, e_in_E);
    TLorentzVector e_out = Build4Vector(e_out_p, e_out_theta, e_out_phi, m_e);
    TLorentzVector p_in(0, 0, 0, m_p);
    TLorentzVector p_out = Build4Vector(p_out_p, p_out_theta, p_out_phi, m_p);
    TLorentzVector kMinus = Build4Vector(kMinus_p, kMinus_theta, kMinus_phi, m_kMinus);
    TLorentzVector kPlus = Build4Vector(kPlus_p, kPlus_theta, kPlus_phi, m_kPlus);
    ComputeKinematics(e_in, e_out, p_in, p_out, kPlus, kMinus);
  }

  // Mixed: one kaon measured, the other missing (set IsMissing=true if k is the missing one)
  DISANAMath(double e_in_E, double e_out_p, double e_out_theta, double e_out_phi, double p_out_p, double p_out_theta, double p_out_phi, double k_p, double k_theta, double k_phi,
             bool IsMissing) {
    TLorentzVector e_in(0, 0, e_in_E, e_in_E);
    TLorentzVector e_out = Build4Vector(e_out_p, e_out_theta, e_out_phi, m_e);
    TLorentzVector p_in(0, 0, 0, m_p);
    TLorentzVector p_out = Build4Vector(p_out_p, p_out_theta, p_out_phi, m_p);
    TLorentzVector k_rec = Build4Vector(k_p, k_theta, k_phi, m_kMinus);  // treat as kaon mass
    TLorentzVector k_miss = e_in + p_in - e_out - p_out - k_rec;
    TLorentzVector kMinus, kPlus;
    if (IsMissing) {
      kMinus = k_miss;
      kPlus = k_rec;
    } else {
      kMinus = k_rec;
      kPlus = k_miss;
    }
    ComputeKinematics(e_in, e_out, p_in, p_out, kPlus, kMinus);
  }

  // Accessors
  double GetQ2() const { return Q2_; }
  double GetxB() const { return xB_; }
  double GetT() const { return t_; }
  double GetTmin() const { return tmin_; }
  double GetPhi() const { return phi_deg_; }
  double GetW() const { return W_; }
  double GetNu() const { return nu_; }
  double Gety() const { return y_; }
  double GetXi() const { return xi_; }
  double GetCostheta_SDHEP() const { return costheta_SDHEP_; }
  double GetPhi_SDHEP() const { return phi_SDHEP_deg_; }
  double GetCosTheta_KK() const { return cosTheta_KK_; }
  double GetCosPhi_KK() const { return cosPhi_KK_; }
  double GetZ_phi() const { return z_phi_; }

  double GetMx2_ep() const { return mx2_ep_; }
  double GetEmiss() const { return emiss_; }
  double GetPTmiss() const { return ptmiss_; }
  double GetMx2_epg() const { return mx2_epg_; }
  double GetMx2_epKpKm() const { return mx2_epKpKm_; }
  double GetDeltaPhi() const { return delta_phi_; }
  double GetTheta_gamma_gamma() const { return theta_gg_; }
  double GetMx2_egamma() const { return mx2_egamma_; }
  double GetMx2_eKpKm() const { return mx2_eKpKm_; }
  double GetMx2_eKp() const { return mx2_eKp_; }
  double GetMx_eKp() const { return std::sqrt(mx2_eKp_); }
  double GetMx2_epKp() const { return mx2_epKp_; }
  double GetMx2_epKm() const { return mx2_epKm_; }
  double GetTheta_e_gamma() const { return Theta_e_gamma_; }
  double GetTheta_g_phimeson() const { return theta_gphi_; }
  double GetTheta_e_phimeson() const { return Theta_e_phimeson_; }
  double GetDeltaE() const { return DeltaE_; }
  double GetMx_epKp() const { return (mx2_epKp_ > 0) ? std::sqrt(mx2_epKp_) : -999.0; }
  double GetCone_Kp() const { return Cone_Kp_; }
  double GetCone_Km() const { return Cone_Km_; }
  double GetCone_p() const { return Cone_p_; }

  double GetMass_pi0() const { return Mass_pi0_; }
  double GetMx2_eppi0() const { return Mx2_eppi0_; }
  double GetEmiss_pi0() const { return emiss_pi0_; }
  double GetMx2_ep_pi0() const { return Mx2_ep_pi0_; }
  double GetMx2_epi0() const { return Mx2_epi0_; }
  double GetPTmiss_pi0() const { return ptmiss_pi0_; }
  double GetTheta_pi0pi0() const { return theta_pi0pi0_; }
  double GetTheta_epho1() const { return theta_epho1_; }
  double GetTheta_epho2() const { return theta_epho2_; }
  double GetTheta_pho1pho2() const { return theta_pho1pho2_; }
  double GetDeltaPhi_pi0() const { return delta_phi_pi0_; }

  double GetCoplanarity_had_normals_deg() const { return coplanarity_had_normals_deg_; }

  double GetVmCut()       const { return vm_cut_; }
  double GetEgammaStar()  const { return Egamma_star_; }

  // -------------------------------------------------------------------------
  // ComputeVmApprox  —  approximate maximal inelasticity for diffractive
  // vector-meson electroproduction  (see Vcut.pdf notes, Section 6)
  //
  //   vm(t) = (m_v^2 * W^2) / (M_p * Q^2 * Sx)
  //            * |t| * ( (W^2 - m_v^2)/2 - |t| )
  //
  // with  S  = 2 M_p E_beam
  //       y  = Q^2 / (xB * S)
  //       Sx = y * S
  //       W^2 = M_p^2 + Sx - Q^2
  //
  // Returns max(vm, 0).  Returns -999 for unphysical kinematics.
  // -------------------------------------------------------------------------
  static double ComputeVmApprox(double Q2, double xB, double t,
                                 double mv, double beam_energy)
  {
    if (Q2 <= 0.0 || xB <= 0.0 || beam_energy <= 0.0) return -999.0;

    const double S   = 2.0 * m_p * beam_energy;
    const double y   = Q2 / (xB * S);
    if (y <= 0.0 || y >= 1.0) return 0.0;

    const double Sx  = y * S;
    const double W2  = m_p * m_p + Sx - Q2;
    if (W2 <= 0.0) return 0.0;

    const double mv2     = mv * mv;
    const double t_abs   = std::abs(t);
    const double bracket = (W2 - mv2) / 2.0 - t_abs;
    if (bracket <= 0.0) return 0.0;          // past the kinematic endpoint

    const double pref = (mv2 * W2) / (m_p * Q2 * Sx);
    const double vm   = pref * t_abs * bracket;
    return std::max(vm, 0.0);
  }

  // Helpers
  // got from the PAC39 need for tmin, Källén function λ(x,y,z) = x^2 + y^2 + z^2 − 2xy − 2xz − 2yz
  double kallen(double x, double y, double z) { return std::max(0.0, x * x + y * y + z * z - 2.0 * (x * y + x * z + y * z)); }
  // t_min for γ* p → φ p given Q2 and xB  (proton at rest DIS vars)
  // Returns t_min (note: usually negative); also useful: -t_min.
  double tmin_phi_from_Q2_xB(double Q2, double xB) {
    // physical masses (GeV)
    const double m_p2 = m_p * m_p;
    const double m_phi2 = m_phi * m_phi;

    // W^2 = M^2 + Q^2 * (1/xB - 1)
    const double s = m_p2 + Q2 * (1.0 / xB - 1.0);  // γ* p c.m. energy squared. :contentReference[oaicite:1]{index=1}

    // masses-squared of a+b→c+d: a=γ* (m^2=-Q^2), b=p, c=φ, d=p
    const double m1sq = -Q2;     // γ*
    const double m2sq = m_p2;    // p
    const double m3sq = m_phi2;  // φ
    const double m4sq = m_p2;    // p

    // handy pieces for PDG 2→2 t-limits at θ=0 (forward)
    const double A = (s + m1sq - m2sq) * (s + m3sq - m4sq);
    const double B = std::sqrt(kallen(s, m1sq, m2sq)) * std::sqrt(kallen(s, m3sq, m4sq));

    // forward-angle limit (a.k.a. "t0" in PDG) — this is the experimental t_min (smallest |t|)
    const double t_min = (m1sq + m3sq) - (A - B) / (2.0 * s);
    return t_min;
  }

  double ComputePhiH(const TVector3 &q1_, const TVector3 &k1_, const TVector3 &q2_) const {
    const double t1 = ((q1_.Cross(k1_)).Dot(q2_)) / std::abs((q1_.Cross(k1_)).Dot(q2_));
    const TVector3 t2 = q1_.Cross(k1_);
    const TVector3 t3 = q1_.Cross(q2_);
    const double n2 = t2.Mag(), n3 = t3.Mag();
    const double t4 = (t2.Dot(t3)) / (n2 * n3);
    return t1 * std::acos(t4) * 180. / pi + 180.;
  }

  void ComputePi0Kinematics(const TLorentzVector &electron_in, const TLorentzVector &electron_out, const TLorentzVector &proton_in, const TLorentzVector &proton_out,
                            const TLorentzVector &photon1, const TLorentzVector &photon2) {
    TLorentzVector q = electron_in - electron_out;
    TLorentzVector pi0 = photon1 + photon2;
    Mass_pi0_ = pi0.Mag();
    Mx2_eppi0_ = (electron_in + proton_in - electron_out - proton_out - pi0).Mag2();
    emiss_pi0_ = (electron_in + proton_in - electron_out - proton_out - pi0).E();
    Mx2_ep_pi0_ = (electron_in + proton_in - electron_out - proton_out).Mag2();
    Mx2_epi0_ = (electron_in + proton_in - electron_out - pi0).Mag2();
    ptmiss_pi0_ = (electron_in + proton_in - electron_out - proton_out - pi0).Vect().Perp();
    theta_pi0pi0_ = pi0.Angle((electron_in + proton_in - (electron_out + proton_out)).Vect()) * 180. / pi;
    delta_phi_pi0_ = std::abs(ComputePhiH(q.Vect(), electron_in.Vect(), pi0.Vect()) - ComputePhiH(q.Vect(), electron_in.Vect(), -proton_out.Vect()));
    theta_epho1_ = electron_out.Angle(photon1.Vect()) * 180. / pi;
    theta_epho2_ = electron_out.Angle(photon2.Vect()) * 180. / pi;
    theta_pho1pho2_ = photon1.Angle(photon2.Vect()) * 180. / pi;
  }

  // Core (DVCS-like)
  void ComputeKinematics(const TLorentzVector &electron_in, const TLorentzVector &electron_out, const TLorentzVector &proton_in, const TLorentzVector &proton_out,
                         const TLorentzVector &photon) {
    TLorentzVector q = electron_in - electron_out;

    Q2_ = -q.Mag2();
    nu_ = q.E();
    y_ = nu_ / electron_in.E();
    W_ = (proton_in + q).Mag();
    xB_ = Q2_ / (2.0 * proton_in.Dot(q));
    t_ = std::abs((proton_in - proton_out).Mag2());

    TVector3 nL = electron_in.Vect().Cross(electron_out.Vect()).Unit();
    TVector3 nH = photon.Vect().Cross(q.Vect()).Unit();
    const double cos_phi = nL.Dot(nH);
    const double sin_phi = (nL.Cross(nH)).Dot(q.Vect().Unit());
    phi_deg_ = (std::atan2(sin_phi, cos_phi) + pi) * 180. / pi;

    const TLorentzVector totI = electron_in + proton_in;
    const TLorentzVector totF = electron_out + proton_out + photon;
    const TLorentzVector miss = totI - totF;

    mx2_ep_ = (totI - electron_out - proton_out).Mag2();
    emiss_ = miss.E();
    ptmiss_ = miss.Vect().Perp();
    mx2_epg_ = miss.Mag2();

    // Coplanarity (hadron vs q)
    TVector3 qv = q.Vect();
    TVector3 ev = electron_in.Vect();
    TVector3 gv = photon.Vect();
    TVector3 pv = proton_out.Vect();
    delta_phi_ = std::abs(ComputePhiH(qv, ev, gv) - ComputePhiH(qv, ev, -pv));
    theta_gg_ = photon.Angle((totI - (electron_out + proton_out)).Vect()) * 180. / pi;

    mx2_egamma_ = (electron_in + proton_in - electron_out - photon).Mag2();
    Theta_e_gamma_ = electron_out.Angle(photon.Vect()) * 180. / pi;
    DeltaE_ = (electron_in.E() + proton_in.E()) - (electron_out.E() + proton_out.E() + photon.E());


    TLorentzVector Delta =  proton_in - proton_out;
    double xi = (proton_in-proton_out).Dot(electron_in)/((proton_in + proton_out).Dot(electron_in));
    double t0 = -4.0 * m_p * m_p * xi * xi / (1.0 - xi * xi);
    double s_hat = (electron_in + Delta).Mag2();
    double tt = Delta.Mag2();
    double s = (electron_in + proton_in).Mag2();

    double deltalprime = Delta.Dot(electron_out);
    double llprime = electron_in.Dot(electron_out);
    double twoPlprime = (proton_in + proton_out).Dot(electron_out);
    double twoepsilonM =  EpsContract(electron_out,electron_in,proton_in+proton_out,Delta);
    double nom = (1-xi*xi)*s_hat*(t0-tt)*(deltalprime-tt/2);

    double costheta = -2*((s_hat-tt)*deltalprime-(s_hat+tt)*llprime)/(s_hat*(s_hat-tt));
    double cosphi = xi*(s-m_p*m_p)*(twoPlprime*xi-deltalprime)+(1+xi)*llprime*tt;
    cosphi = cosphi/((1+xi)*sqrt(nom));
    double sinphi = -xi*twoepsilonM/sqrt(nom);
    double phi_SDHEP = std::atan2(sinphi, cosphi);
    if (phi_SDHEP<0) phi_SDHEP = phi_SDHEP +2*pi;
    double phi_SDHEP_deg = phi_SDHEP * 180/pi;

    double qT = 0.5*sqrt(s_hat*(1-costheta*costheta));

    phi_SDHEP_deg_ = phi_SDHEP_deg;
    xi_ = xi;
    costheta_SDHEP_ = costheta;
    //s_hat_ = s_hat;
    //Q2_ = qT;
  }

  // Core (phi)
  void ComputeKinematics(const TLorentzVector &electron_in, const TLorentzVector &electron_out, const TLorentzVector &proton_in, const TLorentzVector &proton_out,
                         const TLorentzVector &kPlus, const TLorentzVector &kMinus) {
    TLorentzVector q = electron_in - electron_out;
    TLorentzVector phi = kPlus + kMinus;

    Q2_ = -q.Mag2();
    nu_ = q.E();
    y_ = nu_ / electron_in.E();
    W_ = (proton_in + q).Mag();
    xB_ = Q2_ / (2.0 * proton_in.Dot(q));
    t_ = std::abs((proton_in - proton_out).Mag2());
    tmin_ = tmin_phi_from_Q2_xB(Q2_, xB_);
    z_phi_=1-t_*xB_/Q2_;

    TVector3 nL = electron_in.Vect().Cross(electron_out.Vect()).Unit();
    TVector3 nH = q.Vect().Cross(proton_out.Vect()).Unit();
    const double cos_phi = nL.Dot(nH);
    const double sin_phi = (nL.Cross(nH)).Dot(q.Vect().Unit());
    phi_deg_ = (std::atan2(sin_phi, cos_phi) + pi) * 180. / pi;

    const TLorentzVector totI = electron_in + proton_in;
    const TLorentzVector totF = electron_out + proton_out + phi;
    const TLorentzVector miss = totI - totF;

    mx2_ep_ = (totI - electron_out - proton_out).Mag2();
    emiss_ = miss.E();
    ptmiss_ = miss.Vect().Perp();
    mx2_epKpKm_ = miss.Mag2();

    TVector3 qv = q.Vect();
    TVector3 ev = electron_in.Vect();
    TVector3 phiv = phi.Vect();
    TVector3 pv = proton_out.Vect();
    TVector3 kpv = kPlus.Vect();
    TVector3 kmv = kMinus.Vect();
    TVector3 n_qphi = qv.Cross(phiv);
    TVector3 n_pphi = pv.Cross(phiv);

    delta_phi_ = std::abs(ComputePhiH(qv, ev, phiv) - ComputePhiH(qv, ev, -pv));
    theta_gphi_ = phiv.Angle((totI - (electron_out + proton_out)).Vect()) * 180. / pi;
    Cone_p_ = (electron_in + proton_in - electron_out - phi).Angle(pv) * 180. / pi;                    // p vs K+
    Cone_Km_ = (electron_in + proton_in - electron_out - kPlus - proton_out).Angle(kmv) * 180. / pi;   // p vs K−
    Cone_Kp_ = (electron_in + proton_in - electron_out - kMinus - proton_out).Angle(kpv) * 180. / pi;  // p vs K−

    mx2_eKpKm_ = (electron_in + proton_in - electron_out - phi).Mag2();
    mx2_epKp_ = (electron_in + proton_in - electron_out - kPlus - proton_out).Mag2();
    mx2_epKm_ = (electron_in + proton_in - electron_out - kMinus - proton_out).Mag2();
    mx2_eKp_ = (electron_in + proton_in - electron_out - kPlus).Mag2();

    Theta_e_phimeson_ = electron_out.Angle(phi.Vect()) * 180. / pi;
    DeltaE_ = (electron_in.E() + proton_in.E()) - (electron_out.E() + proton_out.E() + phi.E());

    // ------------------------------------------------------------------
    // Radiative-photon energy in the CM frame of the missing system X
    //   E_gamma* = (M_X^2 - M_p^2) / (2 * M_X)
    // where M_X^2 = mx2_eKpKm_ = (p + q - phi)^2
    // vm_cut is computed per-event as a column in DISANAMMUtils (needs
    // beam_energy), but E_gamma* can be computed here.
    // ------------------------------------------------------------------
    {
      const double Mx2 = mx2_eKpKm_;
      if (Mx2 > 0.0) {
        const double Mx      = std::sqrt(Mx2);
        Egamma_star_ = (Mx2 - m_p * m_p) / (2.0 * Mx);
      } else {
        Egamma_star_ = -999.0;
      }
    }
    double coplanarity_had_normals_deg_ = std::numeric_limits<double>::quiet_NaN();
    if (n_qphi.Mag() > 0 && n_pphi.Mag() > 0) {
      n_qphi = n_qphi.Unit();
      n_pphi = n_pphi.Unit();

      double c = n_qphi.Dot(n_pphi);
      // numerical safety
      c = std::max(-1.0, std::min(1.0, c));
      double ang = std::acos(c) * 180. / pi;  // in [0, 180]

      // treat ±normals as equivalent (optional but common)
      coplanarity_had_normals_deg_ = std::min(ang, 180.0 - ang);
    }

    // --------------------------------------------------------------------
    // NEW: cos(theta) and cos(phi) in the K+K- COM (phi rest frame)
    //      Following the rho->pi pi construction:
    //      lab --> (gamma* p) c.m. --> phi (K+K-) c.m.
    // --------------------------------------------------------------------
    cosTheta_KK_ = std::numeric_limits<double>::quiet_NaN();
    cosPhi_KK_ = std::numeric_limits<double>::quiet_NaN();

    // --------------------- 1) gamma* p c.m. frame -----------------------
    TLorentzVector gp_sys = q + proton_in;  // total gamma* + p_in
    TVector3 beta_gp = gp_sys.BoostVector();

    TLorentzVector phi_gp = phi;
    TLorentzVector p_out_gp = proton_out;
    TLorentzVector kPlus_gp = kPlus;
    TLorentzVector kMinus_gp = kMinus;
    TLorentzVector q_gp = q;

    phi_gp.Boost(-beta_gp);
    p_out_gp.Boost(-beta_gp);
    kPlus_gp.Boost(-beta_gp);
    kMinus_gp.Boost(-beta_gp);
    q_gp.Boost(-beta_gp);

    // helicity axis in gamma* p c.m.: along phi momentum (same as -p_out_gp)
    TVector3 zprime(0., 0., 0.);
    if (phi_gp.Vect().Mag2() > 0.0) {
      zprime = phi_gp.Vect().Unit();
    } else if (p_out_gp.Vect().Mag2() > 0.0) {
      zprime = (-p_out_gp.Vect()).Unit();
    }

    // --------------------- 2) phi (K+K-) rest frame ---------------------
    TVector3 beta_phi_gp = phi_gp.BoostVector();

    TLorentzVector kPlus_phi = kPlus_gp;
    TLorentzVector p_out_phi = p_out_gp; // Recoil proton
    TLorentzVector q_phi     = q_gp;     // Virtual photon

    kPlus_phi.Boost(-beta_phi_gp);
    p_out_phi.Boost(-beta_phi_gp);
    q_phi.Boost(-beta_phi_gp);

    // --------------------------------------------------------------------
    // Construct the Helicity Coordinate System (Right-Handed)
    // Reference: Diehl, arXiv:0704.1565, Page 5
    // --------------------------------------------------------------------

    // 1. Z-axis: Opposite to the recoil proton in the meson rest frame.
    //    (Diehl p.5: "p' points in the negative z direction")
    TVector3 z_axis = -1.0 * p_out_phi.Vect().Unit();

    // 2. Y-axis: Normal to the production plane.
    //    Defined by q x p' (virtual photon cross recoil proton).
    //    Note: This is the same in Lab, CM, or Rest frame if q and p' are collinear.
    //    Usually calculated using the vectors available in the rest frame.
    TVector3 y_axis = q_phi.Vect().Cross(p_out_phi.Vect()).Unit();

    // 3. X-axis: Completes the right-handed triad (y cross z).
    TVector3 x_axis = y_axis.Cross(z_axis).Unit();

    // ------------------------ Calculate Angles --------------------------
    // Get the K+ momentum vector in the rest frame
    TVector3 k_vec = kPlus_phi.Vect();

    // Polar Angle (Theta): Angle between K+ and Z-axis
    // Range: [0, pi], cosTheta in [-1, 1]
    cosTheta_KK_ = std::numeric_limits<double>::quiet_NaN();
    if(k_vec.Mag() > 0){
        cosTheta_KK_ = z_axis.Dot(k_vec.Unit());
    }

    // Azimuthal Angle (Phi): Angle of K+ in the XY plane
    // Range: [-pi, pi] usually mapped to degrees [-180, 180] or [0, 360]
    // Calculated using projections onto X and Y axes.
    cosPhi_KK_ = std::numeric_limits<double>::quiet_NaN();
    // Use atan2(y, x) for the full angle, or just calculate cosPhi if that's all you need.
    // However, for consistency with interference terms (sin(phi)), it is better to compute phi first.
    
    double phi_rad = std::atan2(k_vec.Dot(y_axis), k_vec.Dot(x_axis)); 
    cosPhi_KK_ = std::cos(phi_rad);
  }

  // ---------------------------------------------------------------------------
  // dσ/dφ (DVCS) and helpers (unchanged except includes)
  // ---------------------------------------------------------------------------
  std::vector<std::vector<std::vector<TH1D *>>> ComputeDVCS_CrossSection(ROOT::RDF::RNode df, const BinManager &bins, double luminosity) {
    TStopwatch timer;
    timer.Start();
    constexpr double phi_min = 0.0, phi_max = 360.0;
    constexpr int n_phi_bins = 18;

    const auto &q2_bins = bins.GetQ2Bins();
    const auto &t_bins = bins.GetTBins();
    const auto &xb_bins = bins.GetXBBins();

    const size_t n_q2 = q2_bins.size() - 1;
    const size_t n_t = t_bins.size() - 1;
    const size_t n_xb = xb_bins.size() - 1;

    std::vector<std::vector<std::vector<TH1D *>>> hist(n_xb, std::vector<std::vector<TH1D *>>(n_q2, std::vector<TH1D *>(n_t, nullptr)));
    std::vector<std::vector<std::vector<double>>> q2xBtbins(n_xb, std::vector<std::vector<double>>(n_q2, std::vector<double>(n_t, 0.0)));

    for (size_t ix = 0; ix < n_xb; ++ix)
      for (size_t iq = 0; iq < n_q2; ++iq)
        for (size_t it = 0; it < n_t; ++it) {
          const double qmin = q2_bins[iq], qmax = q2_bins[iq + 1];
          const double tmin = t_bins[it], tmax = t_bins[it + 1];
          const double xbmin = xb_bins[ix], xbmax = xb_bins[ix + 1];
          auto name = Form("hphi_q%.1f_t%.1f_xb%.2f", qmin, tmin, xbmin);
          auto title = Form("d#sigma/d#phi (Q^{2}=[%.1f,%.1f], t=[%.1f,%.1f], x_{B}=[%.2f,%.2f])", qmin, qmax, tmin, tmax, xbmin, xbmax);
          hist[ix][iq][it] = new TH1D(name, title, n_phi_bins, phi_min, phi_max);
          hist[ix][iq][it]->SetDirectory(nullptr);
          q2xBtbins[ix][iq][it] = (qmax - qmin) * (tmax - tmin) * (xbmax - xbmin);
        }

    auto findBin = [](double v, const std::vector<double> &e) -> int {
      auto it = std::upper_bound(e.begin(), e.end(), v);
      if (it == e.begin() || it == e.end()) return -1;
      return int(it - e.begin()) - 1;
    };

    auto fill = [&](double Q2, double t, double xB, double phi) {
      int iq = findBin(Q2, q2_bins), it = findBin(t, t_bins), ix = findBin(xB, xb_bins);
      if (iq >= 0 && it >= 0 && ix >= 0) hist[ix][iq][it]->Fill(phi);
    };

    df.Foreach(fill, {"Q2", "t", "xB", "phi"});

    const double bin_width = (pi / 180.) * (phi_max - phi_min) / n_phi_bins;
    for (size_t ix = 0; ix < n_xb; ++ix)
      for (size_t iq = 0; iq < n_q2; ++iq)
        for (size_t it = 0; it < n_t; ++it) {
          TH1D *h = hist[ix][iq][it];
          if (!h) continue;
          const double vol = q2xBtbins[ix][iq][it];
          for (int b = 1; b <= h->GetNbinsX(); ++b) {
            const double N = h->GetBinContent(b);
            const double ds = N / (luminosity * bin_width * vol);
            const double es = std::sqrt(std::max(0.0, N)) / (luminosity * bin_width * vol);
            h->SetBinContent(b, ds);
            h->SetBinError(b, es);
          }
        }

    std::cout << "DVCS cross-sections computed in a single pass.\n";
    timer.Stop();
    std::cout << "Time elapsed: " << timer.RealTime() << " s (real), " << timer.CpuTime() << " s (CPU)\n";
    return hist;
  }

    std::vector<std::vector<std::vector<TH1D *>>> ComputeDVCS_CrossSection_Weighted(ROOT::RDF::RNode df, const BinManager &bins, double luminosity,
                                                                                   DVCSWeightFunction weightFunc = nullptr)
  {
    TStopwatch timer;
    timer.Start();
    constexpr double phi_min = 0.0, phi_max = 360.0;
    constexpr int n_phi_bins = 18;

    const auto &q2_bins = bins.GetQ2Bins();
    const auto &t_bins  = bins.GetTBins();
    const auto &xb_bins = bins.GetXBBins();

    const size_t n_q2 = q2_bins.size() - 1;
    const size_t n_t  = t_bins.size()  - 1;
    const size_t n_xb = xb_bins.size() - 1;

    // -----------------------
    // helper: find bin index
    // -----------------------
    auto findBin = [](double v, const std::vector<double> &e) -> int {
      auto it = std::upper_bound(e.begin(), e.end(), v);
      if (it == e.begin() || it == e.end()) return -1;
      return int(it - e.begin()) - 1;
    };

    // -----------------------
    // allocate final hist + volume
    // -----------------------
    std::vector<std::vector<std::vector<TH1D *>>> hist(
        n_xb, std::vector<std::vector<TH1D *>>(n_q2, std::vector<TH1D *>(n_t, nullptr)));

    std::vector<std::vector<std::vector<double>>> q2xBtbins(
        n_xb, std::vector<std::vector<double>>(n_q2, std::vector<double>(n_t, 0.0)));

    for (size_t ix = 0; ix < n_xb; ++ix)
      for (size_t iq = 0; iq < n_q2; ++iq)
        for (size_t it = 0; it < n_t; ++it) {
          const double qmin = q2_bins[iq], qmax = q2_bins[iq + 1];
          const double tmin = t_bins[it],  tmax = t_bins[it + 1];
          const double xbmin = xb_bins[ix], xbmax = xb_bins[ix + 1];

          auto name  = Form("hphi_q%.2f_t%.2f_xb%.3f", qmin, tmin, xbmin);
          auto title = Form("d#sigma/d#phi (Q^{2}=[%.2f,%.2f], t=[%.2f,%.2f], x_{B}=[%.3f,%.3f])",
                            qmin, qmax, tmin, tmax, xbmin, xbmax);

          hist[ix][iq][it] = new TH1D(name, title, n_phi_bins, phi_min, phi_max);
          hist[ix][iq][it]->SetDirectory(nullptr);
          hist[ix][iq][it]->Sumw2(); // IMPORTANT: weighted errors

          q2xBtbins[ix][iq][it] = (qmax - qmin) * (tmax - tmin) * (xbmax - xbmin);
        }

    // ------------------------------------------------------------
    // thread-safe filling: one histogram set per slot, then merge
    // ------------------------------------------------------------
    const unsigned int nSlots = df.GetNSlots();

    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> hslot(
        nSlots,
        std::vector<std::vector<std::vector<TH1D*>>>(
            n_xb, std::vector<std::vector<TH1D*>>(n_q2, std::vector<TH1D*>(n_t, nullptr))));

    for (unsigned int s = 0; s < nSlots; ++s) {
      for (size_t ix = 0; ix < n_xb; ++ix)
        for (size_t iq = 0; iq < n_q2; ++iq)
          for (size_t it = 0; it < n_t; ++it) {
            // clone "shape" (bins), unique name per slot
            hslot[s][ix][iq][it] = (TH1D*)hist[ix][iq][it]->Clone(
                Form("%s_slot%u", hist[ix][iq][it]->GetName(), s));
            hslot[s][ix][iq][it]->Reset();
            hslot[s][ix][iq][it]->SetDirectory(nullptr);
            //hslot[s][ix][iq][it]->Sumw2();
          }
    }

    auto fillSlot = [&](unsigned int slot, double Q2, double t, double xB, double phi, int pho_det_region, int pro_det_region,
                        double recel_p, double recel_theta, double recel_phi, double recpho_p, double recpho_theta, double recpho_phi,
                        double recpro_p, double recpro_theta, double recpro_phi) {
      const int iq = findBin(Q2, q2_bins);
      const int it = findBin(t,  t_bins);
      const int ix = findBin(xB, xb_bins);
      if (iq < 0 || it < 0 || ix < 0) return;

      const DVCSWeightInput weightInput{pho_det_region, pro_det_region, recel_p, recel_theta, recel_phi,
                                        recpho_p, recpho_theta, recpho_phi, recpro_p, recpro_theta, recpro_phi};
      const double w = weightFunc ? weightFunc(weightInput) : 1.0;
      hslot[slot][ix][iq][it]->Fill(phi, w);
    };

    df.ForeachSlot(fillSlot, {"Q2", "t", "xB", "phi", "pho_det_region", "pro_det_region",
                              "recel_p", "recel_theta", "recel_phi", "recpho_p", "recpho_theta", "recpho_phi",
                              "recpro_p", "recpro_theta", "recpro_phi"});

    // merge slot hists into final
    for (unsigned int s = 0; s < nSlots; ++s) {
      for (size_t ix = 0; ix < n_xb; ++ix)
        for (size_t iq = 0; iq < n_q2; ++iq)
          for (size_t it = 0; it < n_t; ++it) {
            hist[ix][iq][it]->Add(hslot[s][ix][iq][it]);
            delete hslot[s][ix][iq][it];
            hslot[s][ix][iq][it] = nullptr;
          }
    }

    // ------------------------------------------------------------
    // convert yields -> cross section
    // Weighted hist: bin error is sqrt(sum w^2) already.
    // Use Scale so errors scale correctly.
    // ------------------------------------------------------------
    const double bin_width = (pi / 180.) * (phi_max - phi_min) / n_phi_bins;

    for (size_t ix = 0; ix < n_xb; ++ix)
      for (size_t iq = 0; iq < n_q2; ++iq)
        for (size_t it = 0; it < n_t; ++it) {
          TH1D *h = hist[ix][iq][it];
          if (!h) continue;
          const double vol = q2xBtbins[ix][iq][it];
          const double scale = 1.0 / (luminosity * bin_width * vol);
          h->Scale(scale); // scales both content and errors
        }

    std::cout << "DVCS cross-sections computed " << (weightFunc ? "(with custom event weights)" : "(unweighted)") << " in a single pass.\n";
    timer.Stop();
    std::cout << "Time elapsed: " << timer.RealTime() << " s (real), " << timer.CpuTime() << " s (CPU)\n";
    return hist;
  }


  // φ dσ/dt (unchanged)
  std::vector<TH1D *> ComputePhi_CrossSection(ROOT::RDF::RNode df, const BinManager &bins, double luminosity) {
    TStopwatch timer;
    timer.Start();
    const auto &q2_bins = bins.GetQ2Bins();
    const auto &t_bins = bins.GetTBins();
    const size_t n_q2 = q2_bins.size() - 1, n_t = t_bins.size() - 1;
    std::vector<TH1D *> hist(n_q2, nullptr);

    for (size_t iq = 0; iq < n_q2; ++iq) {
      const double qmin = q2_bins[iq], qmax = q2_bins[iq + 1];
      auto hname = Form("phi_cs_q2bin%zu", iq);
      auto htitle = Form("d#sigma/dt for Q^{2}=[%.2f, %.2f]", qmin, qmax);
      hist[iq] = new TH1D(hname, htitle, n_t, &t_bins[0]);
      hist[iq]->SetDirectory(nullptr);
      hist[iq]->GetXaxis()->SetTitle("-t [GeV^{2}]");
      hist[iq]->GetYaxis()->SetTitle("d#sigma/dt [nb/GeV^{2}]");
    }

    auto findBin = [](double v, const std::vector<double> &e) -> int {
      auto it = std::upper_bound(e.begin(), e.end(), v);
      if (it == e.begin() || it == e.end()) return -1;
      return int(it - e.begin()) - 1;
    };
    auto fill = [&](double Q2, double t) {
      int iq = findBin(Q2, q2_bins), it = findBin(t, t_bins);
      if (iq >= 0 && it >= 0) hist[iq]->Fill(t);
    };
    df.Foreach(fill, {"Q2", "t"});

    for (auto *h : hist) {
      if (!h) continue;
      for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double w = h->GetBinWidth(b), N = h->GetBinContent(b);
        const double ds = N / (luminosity * w), es = std::sqrt(std::max(0.0, N)) / (luminosity * w);
        h->SetBinContent(b, ds);
        h->SetBinError(b, es);
      }
    }
    timer.Stop();
    std::cout << "[ComputePhi_CrossSection] Time elapsed: " << timer.RealTime() << " s (real), " << timer.CpuTime() << " s (CPU)\n";
    return hist;
  }
  // === NEW: 3D Beam-Spin Asymmetry ============================================Add commentMore actions
  std::vector<std::vector<std::vector<TH1D *>>> ComputeBeamSpinAsymmetry(const std::vector<std::vector<std::vector<TH1D *>>> &sigma_pos,
                                                                         const std::vector<std::vector<std::vector<TH1D *>>> &sigma_neg, double pol = 1.0) {
    if (sigma_pos.size() != sigma_neg.size()) {
      std::cerr << "ERROR: xB dim mismatch!\n";
      return {};
    }
    const size_t n_xb = sigma_pos.size();
    std::vector<std::vector<std::vector<TH1D *>>> asym(n_xb);
    for (size_t ix = 0; ix < n_xb; ++ix) {
      if (sigma_pos[ix].size() != sigma_neg[ix].size()) {
        std::cerr << "ERROR: Q² dim mismatch at xB " << ix << "!\n";
        return {};
      }

      const size_t n_q2 = sigma_pos[ix].size();
      asym[ix].resize(n_q2);
      for (size_t iq = 0; iq < n_q2; ++iq) {
        if (sigma_pos[ix][iq].size() != sigma_neg[ix][iq].size()) {
          std::cerr << "ERROR: t dim mismatch at (xB,q2)=(" << ix << ',' << iq << ")!\n";
          return {};
        }

        const size_t n_t = sigma_pos[ix][iq].size();
        asym[ix][iq].resize(n_t, nullptr);
        for (size_t it = 0; it < n_t; ++it) {
          TH1D *hp = sigma_pos[ix][iq][it];
          TH1D *hm = sigma_neg[ix][iq][it];
          if (!hp || !hm) {
            std::cerr << "WARNING: null hist @(" << ix << ',' << iq << ',' << it << ")\n";
            continue;
          }

          std::string nm = std::string(hp->GetName()) + "_BSA";
          TH1D *ha = dynamic_cast<TH1D *>(hp->Clone(nm.c_str()));
          ha->Reset();
          ha->SetTitle(("Beam Spin Asymmetry of " + std::string(hp->GetTitle())).c_str());
          for (int b = 1; b <= hp->GetNbinsX(); ++b) {
            const double Np = hp->GetBinContent(b), Nm = hm->GetBinContent(b);
            const double Ep = hp->GetBinError(b), Em = hm->GetBinError(b);
            const double den = Np + Nm, num = Np - Nm;
            double A = (den != 0) ? num / den : 0.0;
            double E = (den != 0) ? 2.0 / (den * den) * std::sqrt(std::pow(Nm * Ep, 2) + std::pow(Np * Em, 2)) : 0.0;
            ha->SetBinContent(b, A / pol);
            ha->SetBinError(b, E / pol);
          }
          asym[ix][iq][it] = ha;
        }
      }
    }
    std::cout << "Beam-Spin Asymmetries (3D) computed.\n";
    return asym;
  }

  std::vector<std::vector<std::vector<TH1D *>>> CalcPi0Corr(ROOT::RDF::RNode df_dvcs_pi0mc, ROOT::RDF::RNode df_pi0_pi0mc, ROOT::RDF::RNode df_dvcs_data,
                                                            ROOT::RDF::RNode df_pi0_data, const BinManager &xBins) {
    const size_t n_t = xBins.GetTBins().size() - 1;
    const size_t n_q2 = xBins.GetQ2Bins().size() - 1;
    const size_t n_xb = xBins.GetXBBins().size() - 1;

    DISANAMath pi0Corr;
    auto df_dvcs_pi0mc_CrossSection = pi0Corr.ComputeDVCS_CrossSection(df_dvcs_pi0mc, xBins, 1);
    auto df_pi0_pi0mc_CrossSection = pi0Corr.ComputeDVCS_CrossSection(df_pi0_pi0mc, xBins, 1);
    auto df_dvcs_data_CrossSection = pi0Corr.ComputeDVCS_CrossSection(df_dvcs_data, xBins, 1);
    auto df_pi0_data_CrossSection = pi0Corr.ComputeDVCS_CrossSection(df_pi0_data, xBins, 1);

    std::vector<std::vector<std::vector<TH1D *>>> hCorr(n_xb, std::vector<std::vector<TH1D *>>(n_q2, std::vector<TH1D *>(n_t, nullptr)));

    for (size_t t_bin = 0; t_bin < n_t; ++t_bin) {
      for (size_t q2_bin = 0; q2_bin < n_q2; ++q2_bin) {
        for (size_t xb_bin = 0; xb_bin < n_xb; ++xb_bin) {
          TH1D *h_dvcs_pi0mc = df_dvcs_pi0mc_CrossSection[xb_bin][q2_bin][t_bin];
          TH1D *h_pi0_pi0mc = df_pi0_pi0mc_CrossSection[xb_bin][q2_bin][t_bin];
          TH1D *h_dvcs_data = df_dvcs_data_CrossSection[xb_bin][q2_bin][t_bin];
          TH1D *h_pi0_data = df_pi0_data_CrossSection[xb_bin][q2_bin][t_bin];

          if (!h_dvcs_pi0mc || !h_pi0_pi0mc || !h_dvcs_data || !h_pi0_data) {
            std::cerr << "Pi0 corr: Missing histogram for Q² bin " << q2_bin << ", xB bin " << xb_bin << ", t bin " << t_bin << "\n";
            continue;
          }
          TH1D *hRatio = static_cast<TH1D *>(h_dvcs_pi0mc->Clone(Form("hPi0Corr_xb%zu_q2%zu_t%zu", xb_bin, q2_bin, t_bin)));
          hRatio->Reset();
          hRatio->Divide(h_dvcs_pi0mc, h_pi0_pi0mc);
          hRatio->Multiply(hRatio, h_pi0_data);
          hRatio->Divide(hRatio, h_dvcs_data);
          hCorr[xb_bin][q2_bin][t_bin] = hRatio;
        }
      }
    }
    return hCorr;
  }

  std::vector<std::vector<std::vector<TH1D *>>> CalcAcceptanceCorr(ROOT::RDF::RNode df_gen_dvcsmc, ROOT::RDF::RNode df_accept_dvcsmc, const BinManager &xBins) {
    const size_t n_t = xBins.GetTBins().size() - 1;
    const size_t n_q2 = xBins.GetQ2Bins().size() - 1;
    const size_t n_xb = xBins.GetXBBins().size() - 1;

    DISANAMath AccCorr;
    auto df_gen_dvcsmc_CrossSection = AccCorr.ComputeDVCS_CrossSection(df_gen_dvcsmc, xBins, 1);
    auto df_accept_dvcsmc_CrossSection = AccCorr.ComputeDVCS_CrossSection(df_accept_dvcsmc, xBins, 1);

    std::vector<std::vector<std::vector<TH1D *>>> hCorr(n_xb, std::vector<std::vector<TH1D *>>(n_q2, std::vector<TH1D *>(n_t, nullptr)));

    for (size_t t_bin = 0; t_bin < n_t; ++t_bin) {
      for (size_t q2_bin = 0; q2_bin < n_q2; ++q2_bin) {
        for (size_t xb_bin = 0; xb_bin < n_xb; ++xb_bin) {
          TH1D *h_gen_dvcsmc = df_gen_dvcsmc_CrossSection[xb_bin][q2_bin][t_bin];
          TH1D *h_accept_dvcsmc = df_accept_dvcsmc_CrossSection[xb_bin][q2_bin][t_bin];

          if (!h_gen_dvcsmc || !h_accept_dvcsmc) {
            std::cerr << "Acceptance Corr: Missing histogram for Q² bin " << q2_bin << ", xB bin " << xb_bin << ", t bin " << t_bin << "\n";
            continue;
          }
          TH1D *hRatio = static_cast<TH1D *>(h_gen_dvcsmc->Clone(Form("hAccCorr_xb%zu_q2%zu_t%zu", xb_bin, q2_bin, t_bin)));
          hRatio->Reset();
          hRatio->Divide(h_accept_dvcsmc, h_gen_dvcsmc);
          hCorr[xb_bin][q2_bin][t_bin] = hRatio;
        }
      }
    }
    return hCorr;
  }

  std::vector<std::vector<std::vector<TH1D *>>> CalcEfficiencyCorr_old(ROOT::RDF::RNode df_dvcsmc_bkg, ROOT::RDF::RNode df_dvcsmc_nobkg, const BinManager &xBins) {
    const size_t n_t = xBins.GetTBins().size() - 1;
    const size_t n_q2 = xBins.GetQ2Bins().size() - 1;
    const size_t n_xb = xBins.GetXBBins().size() - 1;

    DISANAMath EffCorr;
    auto df_dvcsmc_bkg_CrossSection = EffCorr.ComputeDVCS_CrossSection(df_dvcsmc_bkg, xBins, 1);
    auto df_dvcsmc_nobkg_CrossSection = EffCorr.ComputeDVCS_CrossSection(df_dvcsmc_nobkg, xBins, 1);
    std::vector<std::vector<std::vector<TH1D *>>> hCorr(n_xb, std::vector<std::vector<TH1D *>>(n_q2, std::vector<TH1D *>(n_t, nullptr)));
    for (size_t t_bin = 0; t_bin < n_t; ++t_bin) {
      for (size_t q2_bin = 0; q2_bin < n_q2; ++q2_bin) {
        for (size_t xb_bin = 0; xb_bin < n_xb; ++xb_bin) {
          TH1D *h_dvcsmc_bkg = df_dvcsmc_bkg_CrossSection[xb_bin][q2_bin][t_bin];
          TH1D *h_dvcsmc_nobkg = df_dvcsmc_nobkg_CrossSection[xb_bin][q2_bin][t_bin];

          if (!h_dvcsmc_bkg || !h_dvcsmc_nobkg) {
            std::cerr << "Missing histogram for Q² bin " << q2_bin << ", xB bin " << xb_bin << ", t bin " << t_bin << "\n";
            continue;
          }
          TH1D *hRatio = static_cast<TH1D *>(h_dvcsmc_bkg->Clone(Form("hEffCorr_xb%zu_q2%zu_t%zu", xb_bin, q2_bin, t_bin)));
          hRatio->Reset();
          hRatio->Divide(h_dvcsmc_bkg, h_dvcsmc_nobkg);
          hCorr[xb_bin][q2_bin][t_bin] = hRatio;
        }
      }
    }
    return hCorr;
  }

  std::vector<std::vector<std::vector<TH1D *>>> CalcEfficiencyCorr(
    ROOT::RDF::RNode df_dvcsmc_bkg,
    ROOT::RDF::RNode df_dvcsmc_nobkg,
    const BinManager &xBins,
    double I_avg = 60.0,
    double I_mc = 60.0,
    double eff_corr = 1.0) 
  {
    const size_t n_t  = xBins.GetTBins().size()  - 1;
    const size_t n_q2 = xBins.GetQ2Bins().size() - 1;
    const size_t n_xb = xBins.GetXBBins().size() - 1;

    DISANAMath EffCorr;
    auto df_dvcsmc_bkg_CrossSection   = EffCorr.ComputeDVCS_CrossSection(df_dvcsmc_bkg,   xBins, 1);
    auto df_dvcsmc_nobkg_CrossSection = EffCorr.ComputeDVCS_CrossSection(df_dvcsmc_nobkg, xBins, 1);

    std::vector<std::vector<std::vector<TH1D *>>> hCorr(
        n_xb,
        std::vector<std::vector<TH1D *>>(
            n_q2, std::vector<TH1D *>(n_t, nullptr)));

    if (I_mc == 0.0) {
      std::cerr << "[CalcEfficiencyCorr] Error: I_mc = 0, cannot compute I_avg/I_mc.\n";
      return hCorr;
    }

    const double r = I_avg / I_mc;
    const double c = eff_corr;

    std::cout<<"[CalcEfficiencyCorr] I_avg = " << I_avg << ", I_mc = " << I_mc << ", r = " << r << ", eff_corr = " << c << "\n";

    for (size_t t_bin = 0; t_bin < n_t; ++t_bin) {
      for (size_t q2_bin = 0; q2_bin < n_q2; ++q2_bin) {
        for (size_t xb_bin = 0; xb_bin < n_xb; ++xb_bin) {

          TH1D *h_dvcsmc_bkg   = df_dvcsmc_bkg_CrossSection[xb_bin][q2_bin][t_bin];
          TH1D *h_dvcsmc_nobkg = df_dvcsmc_nobkg_CrossSection[xb_bin][q2_bin][t_bin];

          if (!h_dvcsmc_bkg || !h_dvcsmc_nobkg) {
            std::cerr << "[CalcEfficiencyCorr] Missing histogram for "
                      << "xB bin " << xb_bin
                      << ", Q2 bin " << q2_bin
                      << ", t bin "  << t_bin << "\n";
            continue;
          }

          TH1D *hRatio = static_cast<TH1D *>(
              h_dvcsmc_bkg->Clone(Form("hEffCorr_xb%zu_q2%zu_t%zu", xb_bin, q2_bin, t_bin)));
          hRatio->Reset();
          //hRatio->Sumw2();

          const int nBinsX = h_dvcsmc_bkg->GetNbinsX();

          for (int ibin = 1; ibin <= nBinsX; ++ibin) {
            const double N_bkg   = h_dvcsmc_bkg->GetBinContent(ibin);
            const double N_nobkg = h_dvcsmc_nobkg->GetBinContent(ibin);

            const double err_bkg   = h_dvcsmc_bkg->GetBinError(ibin);
            const double err_nobkg = h_dvcsmc_nobkg->GetBinError(ibin);

            double corr_val = 0.0;
            double corr_err = 0.0;

            if (N_nobkg > 0.0) {
              const double eff = N_bkg / N_nobkg;

              // sigma_eff^2 = (sigma_bkg / N_nobkg)^2
              //             + (N_bkg * sigma_nobkg / N_nobkg^2)^2
              const double eff_err2 =
                  std::pow(err_bkg / N_nobkg, 2) +
                  std::pow(N_bkg * err_nobkg / (N_nobkg * N_nobkg), 2);

              const double eff_err = (eff_err2 > 0.0) ? std::sqrt(eff_err2) : 0.0;

              corr_val = (1.0 - (1.0 - eff) * r) * c;
              // corr = (1 - r + r*eff) * c
              // d(corr)/d(eff) = r*c
              corr_err = std::abs(r * c) * eff_err;

              if (corr_val < 0.0) corr_val = 0.0;
              if (corr_val > 1.0) corr_val = 1.0;
            } 
            else {
              corr_val = 0.0;
              corr_err = 0.0;
            }
            hRatio->SetBinContent(ibin, corr_val);
            hRatio->SetBinError(ibin, corr_err);
          }

          hCorr[xb_bin][q2_bin][t_bin] = hRatio;
        }
      }
    }

    return hCorr;
  }


  std::vector<std::vector<std::vector<TH1D *>>> CalcRadiativeCorr(ROOT::RDF::RNode df_dvcs_rad, ROOT::RDF::RNode df_dvcs_norad, const BinManager &xBins) {
    const size_t n_t = xBins.GetTBins().size() - 1;
    const size_t n_q2 = xBins.GetQ2Bins().size() - 1;
    const size_t n_xb = xBins.GetXBBins().size() - 1;

    DISANAMath RadCorr;
    auto df_dvcs_rad_CrossSection = RadCorr.ComputeDVCS_CrossSection(df_dvcs_rad, xBins, 1);
    auto df_dvcs_norad_CrossSection = RadCorr.ComputeDVCS_CrossSection(df_dvcs_norad, xBins, 1);
    std::vector<std::vector<std::vector<TH1D *>>> hCorr(n_xb, std::vector<std::vector<TH1D *>>(n_q2, std::vector<TH1D *>(n_t, nullptr)));
    for (size_t t_bin = 0; t_bin < n_t; ++t_bin) {
      for (size_t q2_bin = 0; q2_bin < n_q2; ++q2_bin) {
        for (size_t xb_bin = 0; xb_bin < n_xb; ++xb_bin) {
          TH1D *h_dvcs_rad = df_dvcs_rad_CrossSection[xb_bin][q2_bin][t_bin];
          TH1D *h_dvcs_norad = df_dvcs_norad_CrossSection[xb_bin][q2_bin][t_bin];
          if (!h_dvcs_rad || !h_dvcs_norad) {
            std::cerr << " Rad Corr: Missing histogram for Q² bin " << q2_bin << ", xB bin " << xb_bin << ", t bin " << t_bin << "\n";
            continue;
          }
          TH1D *hRatio = static_cast<TH1D *>(h_dvcs_rad->Clone(Form("hRadCorr_xb%zu_q2%zu_t%zu", xb_bin, q2_bin, t_bin)));
          hRatio->SetDirectory(nullptr);  // avoid ROOT ownership surprises
          hRatio->Divide(h_dvcs_norad);   // uses default error propagation
          hCorr[xb_bin][q2_bin][t_bin] = hRatio;
        }
      }
    }
    return hCorr;
  }

  std::vector<std::vector<std::vector<TH1D *>>> CalcP1Cut(ROOT::RDF::RNode df_dvcsmc_p1cut, ROOT::RDF::RNode df_dvcsmc_norad, const BinManager &xBins) {
    const size_t n_t = xBins.GetTBins().size() - 1;
    const size_t n_q2 = xBins.GetQ2Bins().size() - 1;
    const size_t n_xb = xBins.GetXBBins().size() - 1;

    DISANAMath P1Corr;
    auto df_dvcsmc_p1cut_CrossSection = P1Corr.ComputeDVCS_CrossSection(df_dvcsmc_p1cut, xBins, 1);
    auto df_dvcsmc_norad_CrossSection = P1Corr.ComputeDVCS_CrossSection(df_dvcsmc_norad, xBins, 1);
    std::vector<std::vector<std::vector<TH1D *>>> hCorr(n_xb, std::vector<std::vector<TH1D *>>(n_q2, std::vector<TH1D *>(n_t, nullptr)));
    for (size_t t_bin = 0; t_bin < n_t; ++t_bin) {
      for (size_t q2_bin = 0; q2_bin < n_q2; ++q2_bin) {
        for (size_t xb_bin = 0; xb_bin < n_xb; ++xb_bin) {
          TH1D *h_dvcsmc_p1cut = df_dvcsmc_p1cut_CrossSection[xb_bin][q2_bin][t_bin];
          TH1D *h_dvcsmc_norad = df_dvcsmc_norad_CrossSection[xb_bin][q2_bin][t_bin];
          if (!h_dvcsmc_p1cut || !h_dvcsmc_norad) {
            std::cerr << " P1 Cut Corr: Missing histogram for Q² bin " << q2_bin << ", xB bin " << xb_bin << ", t bin " << t_bin << "\n";
            continue;
          }
          TH1D *hRatio = static_cast<TH1D *>(h_dvcsmc_p1cut->Clone(Form("hP1CutCorr_xb%zu_q2%zu_t%zu", xb_bin, q2_bin, t_bin)));
          hRatio->Reset();
          hRatio->Divide(h_dvcsmc_norad, h_dvcsmc_p1cut);
          TH1D* hRatioN = hRatio ? (TH1D*)hRatio->Clone("hCorr_norm") : nullptr;
          if (hRatioN) hRatioN->SetDirectory(nullptr);

          NormalizeHCorrByMiddle(hRatioN);

          hCorr[xb_bin][q2_bin][t_bin] = hRatioN;
        }
      }
    }
    return hCorr;
  }
};

#endif  // DISANAMATH_H
