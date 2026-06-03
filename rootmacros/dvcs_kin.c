// dvcs_kin.C
// Run: root -l -q dvcs_kin.C
// Or:  root -l
//      .L dvcs_kin.C
//      dvcs_kin();

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"

#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TString.h"
#include "TBox.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <cstring>

static const double MP = 0.9382720813;

static void PrintP(const char* name, const TLorentzVector& p4){
  TVector3 v = p4.Vect();
  double p  = v.Mag();
  double th = v.Theta() * TMath::RadToDeg();
  double ph = v.Phi()   * TMath::RadToDeg(); // [-180,180]
  std::cout << name
            << "  p=" << p
            << "  th=" << th << " deg"
            << "  ph=" << ph << " deg"
            << "  (px,py,pz)=(" << v.X() << "," << v.Y() << "," << v.Z() << ")"
            << "  E=" << p4.E()
            << std::endl;
}

// Källén lambda
static double Kallen(double a, double b, double c){
  return a*a + b*b + c*c - 2*a*b - 2*a*c - 2*b*c;
}

// Return 0 if ok; negative if impossible
int DVCS_Kin(double xB, double Q2, double t, double phi_deg, double EB,
             TLorentzVector &eprime_lab,
             TLorentzVector &pprime_lab,
             TLorentzVector &gamma_lab)
{
  // (keep your current convention; if you want +180 shift, do it here)
  // phi_deg = fmod(phi_deg + 180.0, 360.0);

  if (xB<=0 || Q2<=0 || EB<=0) return -1;
  if (t>=0) return -2; // expect negative t

  // ---- 1) electron kinematics in lab ----
  double nu = Q2/(2.0*MP*xB);
  double Eprime = EB - nu;
  if (nu<=0 || Eprime<=0) return -3;

  // Q2 = 4 E E' sin^2(theta/2)
  double s2 = Q2/(4.0*EB*Eprime);
  if (s2<0 || s2>1) return -4;
  double theta_e = 2.0*asin(sqrt(s2)); // rad

  TLorentzVector k(0,0,EB,EB); // incoming e along +z, massless
  TLorentzVector kp(Eprime*sin(theta_e), 0, Eprime*cos(theta_e), Eprime);
  eprime_lab = kp;

  TLorentzVector p_target(0,0,0,MP);
  TLorentzVector q = k - kp;
  TVector3 q3 = q.Vect();
  double qmag = q3.Mag();
  if (qmag<=0) return -5;

  // ---- Build Trento-like basis in lab: zhat || q, yhat normal to lepton plane, xhat = y x z
  TVector3 yhat = k.Vect().Cross(kp.Vect());
  if (yhat.Mag()<=0) return -6;
  yhat = yhat.Unit();

  TVector3 zhat = q3.Unit();
  TVector3 xhat = yhat.Cross(zhat);
  if (xhat.Mag()<=0) return -7;
  xhat = xhat.Unit();
  // re-orthogonalize yhat
  yhat = zhat.Cross(xhat).Unit();

  // ---- 2) invariants for gamma* p system ----
  double W2 = MP*MP + 2.0*MP*nu - Q2;
  if (W2<=0) return -8;
  double W = sqrt(W2);

  double lam_in  = Kallen(W2, -Q2, MP*MP);
  double lam_out = Kallen(W2, 0.0,  MP*MP);
  if (lam_in<0 || lam_out<0) return -9;

  double p_in  = sqrt(lam_in)  / (2.0*W);
  double p_out = sqrt(lam_out) / (2.0*W);

  double E1_cm = (W2 + (-Q2) - MP*MP) / (2.0*W); // gamma*
  double E3_cm = (W2 + 0.0   - MP*MP) / (2.0*W); // real gamma (=p_out)

  // ---- 3) solve theta_cm from t ----
  double denom = 2.0*p_in*p_out;
  if (fabs(denom)<1e-15) return -10;

  double cos_th = (t + Q2 + 2.0*E1_cm*E3_cm)/denom;
  if (cos_th < -1.0-1e-9 || cos_th > 1.0+1e-9) return -11;
  if (cos_th<-1) cos_th=-1;
  if (cos_th> 1) cos_th= 1;

  double th_cm  = acos(cos_th);
  double sin_th = sin(th_cm);
  double phi    = phi_deg*TMath::DegToRad();

  // ---- 4) outgoing gamma and proton in CM (axes: z along q) ----
  TVector3 g_cm_3( p_out*sin_th*cos(phi),
                   p_out*sin_th*sin(phi),
                   p_out*cos_th );

  TLorentzVector g_cm(g_cm_3.X(), g_cm_3.Y(), g_cm_3.Z(), p_out); // massless
  TLorentzVector p_cm(-g_cm_3.X(), -g_cm_3.Y(), -g_cm_3.Z(), sqrt(MP*MP + p_out*p_out));

  // ---- 5) boost CM -> lab along +q direction ----
  double beta = qmag/(MP + nu);
  if (beta>=1.0) return -12;

  TVector3 boost_vec(0,0,beta);
  TLorentzVector g_lab_basis = g_cm; g_lab_basis.Boost(boost_vec);
  TLorentzVector p_lab_basis = p_cm; p_lab_basis.Boost(boost_vec);

  // rotate basis -> lab xyz using (xhat,yhat,zhat)
  TVector3 g3b = g_lab_basis.Vect();
  TVector3 p3b = p_lab_basis.Vect();

  TVector3 g3_lab = g3b.X()*xhat + g3b.Y()*yhat + g3b.Z()*zhat;
  TVector3 p3_lab = p3b.X()*xhat + p3b.Y()*yhat + p3b.Z()*zhat;

  gamma_lab.SetPxPyPzE(g3_lab.X(), g3_lab.Y(), g3_lab.Z(), g_lab_basis.E());
  pprime_lab.SetPxPyPzE(p3_lab.X(), p3_lab.Y(), p3_lab.Z(), p_lab_basis.E());

  return 0;
}

// Convenience wrapper that prints
void DVCS_Kin(double xB, double Q2, double t, double phi_deg, double EB){
  TLorentzVector epr, ppr, gam;
  int rc = DVCS_Kin(xB,Q2,t,phi_deg,EB, epr,ppr,gam);
  if (rc!=0){
    std::cout << "DVCS_Kin: kinematics not feasible, rc=" << rc << std::endl;
    return;
  }
  PrintP("e'  :", epr);
  PrintP("p'  :", ppr);
  PrintP("gam :", gam);
}

// -------- helper: Q2(xB) curve for fixed theta_e ----------
static double Q2_of_xB_for_thetae(double xB, double EB, double thetae_deg){
  double s = pow(sin(0.5 * thetae_deg * TMath::DegToRad()), 2);
  double num = 4.0 * EB * EB * s;
  double den = 1.0 + (2.0 * EB * s) / (MP * xB);
  return num / den;
}

// -------- helper: Q2(xB) curve for fixed W ----------
static double Q2_of_xB_for_W(double xB, double W){
  if (xB >= 1.0) return -1.0;
  double W2 = W*W;
  double num = (W2 - MP*MP) * xB;
  double den = (1.0 - xB);
  return num / den;
}

// -------- helper: Q2(xB) curve for fixed y ----------
static double Q2_of_xB_for_y(double xB, double EB, double y0){
  return 2.0 * MP * EB * y0 * xB;
}

// -------- helper: Q2(xB) curve for fixed p_e' ----------
static double Q2_of_xB_for_pe(double xB, double EB, double pe0){
  return 2.0 * MP * xB * (EB - pe0);
}

// Draw user-defined grid lines (your specified numbers; NOT hist bins)
void DrawCustomGrid(const std::vector<double>& xB_lines,
                    const std::vector<double>& Q2_lines,
                    double xBmin, double xBmax,
                    double Q2min, double Q2max,
                    int lineStyle=3, int lineWidth=2, int lineColor=kBlack)
{
  for (double x : xB_lines){
    if (x < xBmin || x > xBmax) continue;
    TLine* lv = new TLine(x, Q2min, x, Q2max);
    lv->SetLineStyle(lineStyle);
    lv->SetLineWidth(lineWidth);
    lv->SetLineColor(lineColor);
    lv->Draw("SAME");
  }
  for (double q2 : Q2_lines){
    if (q2 < Q2min || q2 > Q2max) continue;
    TLine* lh = new TLine(xBmin, q2, xBmax, q2);
    lh->SetLineStyle(lineStyle);
    lh->SetLineWidth(lineWidth);
    lh->SetLineColor(lineColor);
    lh->Draw("SAME");
  }
}

// ===================== MAIN PLOT FUNCTION (with theta_gamma band + theta_p band) =====================
// Band rule you requested:
//   - gamma band drawn ONLY if (draw_gam_band && gam_lo != gam_hi)
//   - proton band drawn ONLY if (draw_p_band   && p_lo   != p_hi)
// You can set each band's color/alpha separately.
void PlotThetaMap(double EB=7.546, double t=-0.55, double phi_deg=0.0,
                  const char* which="gamma",
                  double xBmin=0.125, double xBmax=0.30, int nXB=120,
                  double Q2min=1.0,   double Q2max=2.0,  int nQ2=140,
                  double thetae0_deg=10.0,
                  double W0=2.0,
                  double y0=0.0,     // y=const ; y0>0 draws
                  double pe0=0.0,    // p_e'=const ; pe0>0 draws
                  bool   draw_custom_grid=true,
                  double zmin=0.0, double zmax=40.0,
                  // ---- gamma band ----
                  bool   draw_gam_band=true,
                  double gam_lo_deg=2.5,
                  double gam_hi_deg=6.0,
                  int    gam_band_color=kRed,
                  double gam_band_alpha=0.35,
                  // ---- proton band ----
                  bool   draw_p_band=true,
                  double p_lo_deg=20.0,
                  double p_hi_deg=30.0,
                  int    p_band_color=kAzure+7,
                  double p_band_alpha=0.25,
                  const char* outname="theta_map.png")
{
  gStyle->SetOptStat(0);

  // 2D hist (value = theta of chosen particle)
  TString hname  = TString::Format("hTheta_%s", which);
  TString htitle = TString::Format("#theta_{%s} in lab; x_{B}; Q^{2} [GeV^{2}]", which);
  TH2D* h = new TH2D(hname, htitle, nXB, xBmin, xBmax, nQ2, Q2min, Q2max);

  // flags for band drawing
  bool doGamBand = draw_gam_band && (fabs(gam_hi_deg - gam_lo_deg) > 1e-12);
  bool doPBand   = draw_p_band   && (fabs(p_hi_deg   - p_lo_deg)   > 1e-12);

  // store whether each bin is inside each band
  std::vector<char> inGamBand((nXB+1)*(nQ2+1), 0);
  std::vector<char> inPBand  ((nXB+1)*(nQ2+1), 0);
  auto idx = [&](int ix, int iy){ return ix*(nQ2+1) + iy; };

  for (int ix=1; ix<=nXB; ix++){
    double xB = h->GetXaxis()->GetBinCenter(ix);
    for (int iy=1; iy<=nQ2; iy++){
      double Q2 = h->GetYaxis()->GetBinCenter(iy);

      TLorentzVector epr, ppr, gam;
      int rc = DVCS_Kin(xB, Q2, t, phi_deg, EB, epr, ppr, gam);
      if (rc != 0) continue;

      // value to plot
      double theta_deg = 0.0;
      if (strcmp(which,"gamma")==0){
        theta_deg = gam.Theta()*TMath::RadToDeg();
      } else if (strcmp(which,"e")==0 || strcmp(which,"eprime")==0){
        theta_deg = epr.Theta()*TMath::RadToDeg();
      } else if (strcmp(which,"p")==0 || strcmp(which,"proton")==0){
        theta_deg = ppr.Theta()*TMath::RadToDeg();
      } else {
        std::cout << "[ERROR] which must be gamma/e/p, got: " << which << std::endl;
        return;
      }
      h->SetBinContent(ix, iy, theta_deg);

      // band decisions use theta_gamma and theta_p in lab
      if (doGamBand){
        double thg = gam.Theta()*TMath::RadToDeg();
        if (thg >= gam_lo_deg && thg <= gam_hi_deg) inGamBand[idx(ix,iy)] = 1;
      }
      if (doPBand){
        double thp = ppr.Theta()*TMath::RadToDeg();
        if (thp >= p_lo_deg && thp <= p_hi_deg) inPBand[idx(ix,iy)] = 1;
      }
    }
  }

  // Curves: theta_e=thetae0, W=W0, y=y0, pe'=pe0
  const int N = 800;
  TGraph* gThetaE = new TGraph(); gThetaE->SetName("gThetaE");
  TGraph* gW      = new TGraph(); gW->SetName("gW");
  TGraph* gY      = new TGraph(); gY->SetName("gY");
  TGraph* gPe     = new TGraph(); gPe->SetName("gPe");

  int p1=0, p2=0, p3=0, p4=0;
  for (int i=0; i<N; i++){
    double xB = xBmin + (xBmax-xBmin) * (i/(double)(N-1));

    double Q2t = Q2_of_xB_for_thetae(xB, EB, thetae0_deg);
    if (Q2t>=Q2min && Q2t<=Q2max && std::isfinite(Q2t)) gThetaE->SetPoint(p1++, xB, Q2t);

    double Q2w = Q2_of_xB_for_W(xB, W0);
    if (Q2w>=Q2min && Q2w<=Q2max && std::isfinite(Q2w)) gW->SetPoint(p2++, xB, Q2w);

    if (y0 > 0.0){
      double Q2y = Q2_of_xB_for_y(xB, EB, y0);
      if (Q2y>=Q2min && Q2y<=Q2max && std::isfinite(Q2y)) gY->SetPoint(p3++, xB, Q2y);
    }

    if (pe0 > 0.0){
      double Q2p = Q2_of_xB_for_pe(xB, EB, pe0);
      if (Q2p>=Q2min && Q2p<=Q2max && std::isfinite(Q2p)) gPe->SetPoint(p4++, xB, Q2p);
    }
  }

  // ============ draw ============
  TCanvas* c = new TCanvas("cThetaMap","cThetaMap",1100,850);

  h->GetZaxis()->SetRangeUser(zmin, zmax);
  h->Draw("COLZ");

  // --- gamma band overlay (semi-transparent boxes per bin) ---
  if (doGamBand){
    for (int ix=1; ix<=nXB; ix++){
      double x1 = h->GetXaxis()->GetBinLowEdge(ix);
      double x2 = h->GetXaxis()->GetBinUpEdge(ix);
      for (int iy=1; iy<=nQ2; iy++){
        if (!inGamBand[idx(ix,iy)]) continue;
        double y1 = h->GetYaxis()->GetBinLowEdge(iy);
        double y2 = h->GetYaxis()->GetBinUpEdge(iy);

        TBox* b = new TBox(x1, y1, x2, y2);
        b->SetFillColorAlpha(gam_band_color, gam_band_alpha);
        b->SetLineColorAlpha(gam_band_color, 0.0);
        b->Draw("SAME");
      }
    }
  }

  // --- proton band overlay ---
  if (doPBand){
    for (int ix=1; ix<=nXB; ix++){
      double x1 = h->GetXaxis()->GetBinLowEdge(ix);
      double x2 = h->GetXaxis()->GetBinUpEdge(ix);
      for (int iy=1; iy<=nQ2; iy++){
        if (!inPBand[idx(ix,iy)]) continue;
        double y1 = h->GetYaxis()->GetBinLowEdge(iy);
        double y2 = h->GetYaxis()->GetBinUpEdge(iy);

        TBox* b = new TBox(x1, y1, x2, y2);
        b->SetFillColorAlpha(p_band_color, p_band_alpha);
        b->SetLineColorAlpha(p_band_color, 0.0);
        b->Draw("SAME");
      }
    }
  }

  // --- your specified "grid" lines ---
  if (draw_custom_grid){
    std::vector<double> xB_lines = {0.125, 0.150, 0.180, 0.210, 0.240, 0.285, 0.350, 0.430};
    std::vector<double> Q2_lines = {1.00, 1.25, 1.50, 1.75, 2.00, 2.40, 2.90, 3.50};
    DrawCustomGrid(xB_lines, Q2_lines, xBmin, xBmax, Q2min, Q2max, 3, 2, kBlack);
  }

  // --- style curves ---
  gThetaE->SetLineColor(kRed);
  gThetaE->SetLineStyle(1);
  gThetaE->SetLineWidth(3);

  gW->SetLineColor(kBlue);
  gW->SetLineStyle(2);
  gW->SetLineWidth(3);

  gY->SetLineColor(kGreen+2);
  gY->SetLineStyle(3);
  gY->SetLineWidth(3);

  gPe->SetLineColor(kMagenta);
  gPe->SetLineStyle(7);
  gPe->SetLineWidth(3);

  gThetaE->Draw("L SAME");
  gW->Draw("L SAME");
  if (y0  > 0.0) gY->Draw("L SAME");
  if (pe0 > 0.0) gPe->Draw("L SAME");

  // --- legend ---
  TLegend* leg = new TLegend(0.14, 0.70, 0.78, 0.92);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  if (strcmp(which,"e")!=0){
    leg->AddEntry((TObject*)0,
                  TString::Format("E_{b}=%.3f GeV, t=%.3f GeV^{2}, #phi=%.1f^{#circ}", EB, t, phi_deg),
                  "");
  } else {
    leg->AddEntry((TObject*)0,
                  TString::Format("E_{b}=%.3f GeV", EB),
                  "");
  }

  // band legend entries (use dummy boxes)
  if (doGamBand){
    TBox* bleg = new TBox(0,0,1,1);
    bleg->SetFillColorAlpha(gam_band_color, gam_band_alpha);
    bleg->SetLineColorAlpha(gam_band_color, 0.0);
    leg->AddEntry(bleg, TString::Format("#theta_{#gamma} in [%.1f, %.1f] deg band", gam_lo_deg, gam_hi_deg), "f");
  }
  if (doPBand){
    TBox* blegp = new TBox(0,0,1,1);
    blegp->SetFillColorAlpha(p_band_color, p_band_alpha);
    blegp->SetLineColorAlpha(p_band_color, 0.0);
    leg->AddEntry(blegp, TString::Format("#theta_{p} in [%.1f, %.1f] deg band", p_lo_deg, p_hi_deg), "f");
  }

  leg->AddEntry(gThetaE, TString::Format("#theta_{e}=%.1f^{#circ}", thetae0_deg), "l");
  leg->AddEntry(gW,      TString::Format("W=%.2f GeV", W0), "l");
  if (y0  > 0.0) leg->AddEntry(gY,  TString::Format("y=%.2f", y0), "l");
  if (pe0 > 0.0) leg->AddEntry(gPe, TString::Format("p_{e'}=%.2f GeV", pe0), "l");
  leg->Draw();

  c->SaveAs(outname);
  std::cout << "[Saved] " << outname << std::endl;
}

// Auto-run entry point
void dvcs_kin(){
  double t   = -0.150;
  double phi = 0.0;
  double EB  = 7.546;

  double xB  = 0.23;
  double Q2  = 1.3;

  DVCS_Kin(xB, Q2, t, phi, EB);
  TLorentzVector epr, ppr, gam;
  int rc = DVCS_Kin(xB, Q2, t, phi, EB, epr, ppr, gam);
  std::cout << "proton theta = " << ppr.Theta()*TMath::RadToDeg() << " deg" << std::endl;

  // Example 1: gamma theta map, with BOTH bands
  PlotThetaMap(EB, t, phi, "gamma",
               0.1, 0.6, 120,
               1.0,   5.0,  140,
               9.3,   // theta_e0
               2.0,    // W0
               0.0,    // y0
               2.0,    // pe0
               true,   // custom grid
               2.5, 40.0,   // z range
               true,  5, 35.0,  kRed, 0.50,     // gamma band (lo,hi,color,alpha)
               true,  0.0,  0.0,  kAzure+7, 0.50, // proton band (lo,hi,color,alpha)
               "theta_gamma_map.pdf");

  // Example 2: e theta map, band_lo==band_hi => band NOT drawn
  PlotThetaMap(EB, t, phi, "e",
               0.1, 0.60, 120,
               1.0,   5.0,  140,
               9.3,   // theta_e0
               2.0,    // W0
               0.0,    // y0
               2.0,    // pe0
               true,   // custom grid
               2.5, 40.0,   // z range
               true,  2.5, 2.5, kRed, 0.50,         // gamma band OFF (lo==hi)
               true,  25.0, 25.0, kAzure+7, 0.50,   // proton band OFF (lo==hi)
               "theta_e_map.pdf");
  PlotThetaMap(EB, t, phi, "p",
               0.1, 0.43, 120,
               1.0,   3.5,  140,
               9.3,   // theta_e0
               2.0,    // W0
               0.0,    // y0
               2.0,    // pe0
               true,   // custom grid
               5, 80.0,   // z range
               true,  25.0, 25.0,  kRed, 0.50,     // gamma band (lo,hi,color,alpha)
               true,  57,  62.0,  kBlack, 0.50, // proton band (lo,hi,color,alpha)
               "theta_p_map.pdf");
}
