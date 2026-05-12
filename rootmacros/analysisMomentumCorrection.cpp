#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <map>
#include <cmath>
#include <TF1.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TSpline.h>
#include <memory>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <array>
#include <limits>


using namespace ROOT::VecOps;



//================ Utility =================
// ===== Helpers =====
static inline std::string DetNorm(const std::string& det) {
    if (det == "DC") return "CD";   // allow "DC" alias
    return det;
}
static inline std::string ParticleNameSimple(int pid) {
    if (pid == 11) return "electron";
    if (pid == 22) return "photon";
    if (pid == 2212) return "proton";
    return "pid" + std::to_string(pid);
}

static inline bool IsProtonFD(int pid, const std::string& det) {
    return (pid == 2212 && DetNorm(det) == "FD");
}
static inline double EvalA(double t, const std::array<double,3>& p){ return p[0] + p[1]*t + p[2]*t*t; }
static inline double EvalB(double t, const std::array<double,3>& p){ return p[0] + p[1]*t + p[2]*t*t; }
static inline double EvalC(double t, const std::array<double,3>& p){ return p[0] + p[1]*t + p[2]*t*t; }

// Reads coefficients (p0,p1,p2) from ParamFits text like:
//
// # Fit expression: [0] + [1]*x+ [2]*x*x
// p0 = <val> ± <err>
// p1 = <val> ± <err>
// p2 = <val> ± <err>
//
static bool ReadQuadParams(const std::string& txtPath, std::array<double,3>& pars) {
    std::ifstream fin(txtPath);
    if (!fin) return false;
    std::string line;
    bool found[3] = {false, false, false};
    while (std::getline(fin, line)) {
        // tolerate both "±" and "+/-"
        auto eat = [&](const std::string& key)->bool {
            size_t k = line.find(key);
            if (k == std::string::npos) return false;
            size_t eq = line.find('=', k);
            if (eq == std::string::npos) return false;
            std::string rhs = line.substr(eq+1);
            // strip after ± or +/-
            size_t pm = rhs.find('\xC2'); // start of UTF-8 ±
            size_t slash = rhs.find("+/-");
            if (pm != std::string::npos) rhs = rhs.substr(0, pm);
            else if (slash != std::string::npos) rhs = rhs.substr(0, slash);
            try {
                double val = std::stod(rhs);
                int idx = (key=="p0")?0:(key=="p1")?1:2;
                pars[idx] = val;
                found[idx] = true;
                return true;
            } catch (...) { return false; }
        };
        eat("p0");
        eat("p1");
        eat("p2");
        if (found[0] && found[1] && found[2]) break;
    }
    return (found[0] && found[1] && found[2]);
}

struct XYErrPoints {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> err;
};

static bool ReadPointFile3Cols(const std::string& txtPath, XYErrPoints& pts) {
    std::ifstream fin(txtPath);
    if (!fin) return false;

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        double x = 0.0, y = 0.0, ey = 0.0;
        if (!(iss >> x >> y >> ey)) continue;

        pts.x.push_back(x);
        pts.y.push_back(y);
        pts.err.push_back(ey);
    }

    return !pts.x.empty();
}

static std::string ThetaDirToLabel(const std::string& thetaDirName) {
    std::size_t us = thetaDirName.find('_');
    if (us == std::string::npos || us + 1 >= thetaDirName.size()) return thetaDirName;

    std::string body = thetaDirName.substr(us + 1);
    if (body.rfind("lt", 0) == 0) {
        return "#theta < " + body.substr(2) + "#circ";
    }

    std::size_t toPos = body.find("to");
    if (toPos != std::string::npos) {
        return body.substr(0, toPos) + "#circ < #theta < " + body.substr(toPos + 2) + "#circ";
    }

    if (body.rfind("ge", 0) == 0) {
        return "#theta #geq " + body.substr(2) + "#circ";
    }

    return thetaDirName;
}

static bool ThetaDirToBoundsDeg(const std::string& thetaDirName, double& thetaMinDeg, double& thetaMaxDeg) {
    std::size_t us = thetaDirName.find('_');
    if (us == std::string::npos || us + 1 >= thetaDirName.size()) return false;

    std::string body = thetaDirName.substr(us + 1);
    if (body.rfind("lt", 0) == 0) {
        thetaMinDeg = 0.0;
        thetaMaxDeg = std::stod(body.substr(2));
        return true;
    }

    std::size_t toPos = body.find("to");
    if (toPos != std::string::npos) {
        thetaMinDeg = std::stod(body.substr(0, toPos));
        thetaMaxDeg = std::stod(body.substr(toPos + 2));
        return true;
    }

    if (body.rfind("ge", 0) == 0) {
        thetaMinDeg = std::stod(body.substr(2));
        thetaMaxDeg = 180.0;
        return true;
    }

    return false;
}

static std::vector<std::pair<double,double>> MidpointsToBinEdges(const std::vector<double>& mids) {
    std::vector<std::pair<double,double>> edges;
    if (mids.empty()) return edges;

    edges.resize(mids.size());
    for (std::size_t i = 0; i < mids.size(); ++i) {
        double lo = 0.0;
        double hi = 0.0;

        if (mids.size() == 1) {
            lo = mids[i] - 0.5;
            hi = mids[i] + 0.5;
        } else if (i == 0) {
            const double halfStep = 0.5 * (mids[1] - mids[0]);
            lo = mids[0] - halfStep;
            hi = 0.5 * (mids[0] + mids[1]);
        } else if (i + 1 == mids.size()) {
            const double halfStep = 0.5 * (mids[i] - mids[i - 1]);
            lo = 0.5 * (mids[i - 1] + mids[i]);
            hi = mids[i] + halfStep;
        } else {
            lo = 0.5 * (mids[i - 1] + mids[i]);
            hi = 0.5 * (mids[i] + mids[i + 1]);
        }

        edges[i] = {lo, hi};
    }
    return edges;
}

// Core evaluator for Δp given p, θ and quadratic-in-θ parameter sets
static inline double DeltaP_of_p_theta(int pid,
                                       const std::string& det,
                                       double p,
                                       double thetaDeg,
                                       const std::array<double,3>& Acoef,
                                       const std::array<double,3>& Bcoef,
                                       const std::array<double,3>& Ccoef)
{
    double A = EvalA(thetaDeg, Acoef);
    double B = EvalB(thetaDeg, Bcoef);
    double C = EvalC(thetaDeg, Ccoef);

    if (IsProtonFD(pid, det)) {
        if (p <= 0) return std::numeric_limits<double>::quiet_NaN();
        return A + B/p + C/(p*p);
    } else {
        return A + B*p + C*p*p;
    }
}

// ====== OPTION 1: from ParamFits text files ======
void PlotMomentumCorrectionVsTheta_FromParamFits(const int selectedPid,
                                                 const std::string& selecteddetector, // "FD","CD"/"DC","ALL"
                                                 const std::vector<double>& pValues,  // e.g. {0.3,0.5,1.0,1.5,2.0}
                                                 double thetaMinDeg,
                                                 double thetaMaxDeg,
                                                 int nThetaPoints = 200,
                                                 const std::string& baseDir = "ParticleDeltaPPlots/ParamFits", bool isOutBend = true)
{
    // ---- inputs & files -----------------------------------------------------
    std::string det = DetNorm(selecteddetector);
    std::string prefix = ParticleNameSimple(selectedPid) + "_" + det;

    std::string fA = baseDir + "/" + prefix + "_A_p_vs_theta.txt";
    std::string fB = baseDir + "/" + prefix + "_B_p_vs_theta.txt";
    std::string fC = baseDir + "/" + prefix + "_C_p_vs_theta.txt";
    std::array<double,3> Acoef{}, Bcoef{}, Ccoef{};
    bool okA = ReadQuadParams(fA, Acoef);
    bool okB = ReadQuadParams(fB, Bcoef);
    bool okC = ReadQuadParams(fC, Ccoef);

    if (!(okA && okB && okC)) {
        std::cerr << "[PlotMomentumCorrectionVsTheta_FromParamFits] Could not read parameter files for "
                  << prefix << ". Expected:\n  " << fA << "\n  " << fB << "\n  " << fC << "\n";
        return;
    }

    if (thetaMaxDeg <= thetaMinDeg) std::swap(thetaMinDeg, thetaMaxDeg);
    nThetaPoints = std::max(nThetaPoints, 10);
    if (pValues.empty()) {
        std::cerr << "[PlotMomentumCorrectionVsTheta_FromParamFits] pValues is empty.\n";
        return;
    }

    // ---- canvas & style (local tweaks, no global gStyle upheaval) -----------
    auto *c = new TCanvas(Form("c_deltaP_vs_theta_%s_%s",
                               ParticleNameSimple(selectedPid).c_str(), det.c_str()),
                          "", 1400, 950);
    c->SetMargin(0.12, 0.04, 0.12, 0.07);  // L, R, B, T
    c->SetTicks(1,1);
    c->SetGridx();
    c->SetGridy();

    auto *mg = new TMultiGraph();
    mg->SetTitle(Form("#delta p vs #theta  (%s, %s);#theta [deg];#delta p [GeV]",
                      ParticleNameSimple(selectedPid).c_str(), det.c_str()));

    // A pleasant, high-contrast color cycle
    const std::vector<int> colors = {
        kAzure+2, kOrange+7, kTeal+3, kMagenta+2, kGreen+2,
        kRed+1, kViolet+5, kBlue+1, kPink+7, kCyan+2
    };
    const int baseMarker = 20;

    // Legend: adapt columns for many lines
    /*int ncol = (int)pValues.size() > 8 ? 3 : ((int)pValues.size() > 4 ? 2 : 1);
    auto *leg = new TLegend(0.62, 0.62, 0.90, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.028);
    leg->SetNColumns(ncol);
    leg->SetHeader("p (GeV)", "C");
    */
   TLegend *leg = new TLegend(0.58, 0.68, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  
  leg->Draw();

    // ---- build curves; track global y-range --------------------------------
    double yMin = std::numeric_limits<double>::infinity();
    double yMax = -std::numeric_limits<double>::infinity();

    std::vector<TGraph*> keepAlive; keepAlive.reserve(pValues.size());

    for (size_t ip = 0; ip < pValues.size(); ++ip) {
        double p = pValues[ip];
        std::vector<double> x(nThetaPoints), y(nThetaPoints);
        for (int i = 0; i < nThetaPoints; ++i) {
            double t = thetaMinDeg + (thetaMaxDeg - thetaMinDeg) * (double(i) / (nThetaPoints - 1));
            x[i] = t;
            y[i] = DeltaP_of_p_theta(selectedPid, det, p, t, Acoef, Bcoef, Ccoef);
            yMin = std::min(yMin, y[i]);
            yMax = std::max(yMax, y[i]);
        }

        auto *gr = new TGraph(nThetaPoints, x.data(), y.data());
        gr->SetLineWidth(3);
        gr->SetLineColor(colors[ip % colors.size()]);
        gr->SetMarkerStyle(baseMarker + int(ip) % 5);
        gr->SetMarkerSize(0.8);
        gr->SetMarkerColor(gr->GetLineColor());

        // draw points subtly on top of lines for readability
        mg->Add(gr, "L");
        keepAlive.push_back(gr);

        //leg->AddEntry(gr, Form("%.3g", p), "l");
        leg->AddEntry(gr, Form("p = %.2f GeV", pValues[ip]), "l");
    }

    // ---- draw, format axes, set ranges -------------------------------------
    mg->Draw("A");

    // Axis cosmetics
    auto *xax = mg->GetXaxis();
    auto *yax = mg->GetYaxis();
    xax->SetTitleOffset(1.2);
    yax->SetTitleOffset(1.25);
    xax->SetTitleSize(0.045);
    yax->SetTitleSize(0.045);
    xax->SetLabelSize(0.038);
    yax->SetLabelSize(0.038);
    xax->SetMoreLogLabels(false);
    xax->SetNdivisions(510, /*optimize=*/true);
    yax->SetNdivisions(506, /*optimize=*/true);

    // Sensible default ranges (your original special cases) with fallback to autoscale
    bool appliedSpecial = false;
    if (det == "FD" && selectedPid == 2212) {
        yax->SetRangeUser(0.0, 0.051);
        appliedSpecial = true;
    } else if ((det == "CD" || det == "DC") && selectedPid == 2212) {
        yax->SetRangeUser(-0.03, 0.03);
        appliedSpecial = true;
    } else if (det == "ALL") {
        yax->SetRangeUser(-0.15, 0.15);
        appliedSpecial = true;
    }
    if (!appliedSpecial) {
        // Auto-pad the y-range by 15%
        if (std::isfinite(yMin) && std::isfinite(yMax) && yMax > yMin) {
            double pad = 0.15 * std::max(std::abs(yMax), std::abs(yMin));
            yax->SetRangeUser(yMin - pad, yMax + pad);
        }
    }

    // Zero line for reference
    TLine *z = new TLine(thetaMinDeg, 0., thetaMaxDeg, 0.);
    z->SetLineStyle(2);
    z->SetLineColor(kGray+2);
    z->Draw("same");

    // Legend & an unobtrusive subtitle
    leg->Draw();

    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextSize(0.032);
    lat.SetTextColor(kGray+2);
    //lat.DrawLatex(0.13, 0.94,
      //            Form("#font[42]{%s, %s   (%g#circ to %g#circ)}",
       ///                ParticleNameSimple(selectedPid).c_str(), det.c_str(),
        //               thetaMinDeg, thetaMaxDeg));

    gPad->Modified();
    gPad->Update();

    // ---- save outputs -------------------------------------------------------
    gSystem->Exec(("mkdir -p " + baseDir).c_str());
    std::string outBase = baseDir + "/" + prefix + "_deltaP_vs_theta_multiP";
    std::string outPng  = outBase + ".png";
    std::string outPdf  = outBase + ".pdf";

    c->SaveAs(outPng.c_str());
    c->SaveAs(outPdf.c_str());
    std::cout << "Saved: " << outPng << "\nSaved: " << outPdf << std::endl;
}

void PlotMomentumCorrection_AllDetectors_FromFiles(const int selectedPid,
                                                   const std::vector<double>& pValues,
                                                   // theta ranges to sample for each detector panel:
                                                   double thetaMinFD,  double thetaMaxFD,
                                                   double thetaMinCD,  double thetaMaxCD,
                                                   double thetaMinALL, double thetaMaxALL,
                                                   int nThetaPoints = 200,
                                                   const std::string& baseDir = "ParticleDeltaPPlots/ParamFits")
{
    PlotMomentumCorrectionVsTheta_FromParamFits(selectedPid, "FD",  pValues, thetaMinFD,  thetaMaxFD,  nThetaPoints, baseDir);
    PlotMomentumCorrectionVsTheta_FromParamFits(selectedPid, "CD",  pValues, thetaMinCD,  thetaMaxCD,  nThetaPoints, baseDir);
    PlotMomentumCorrectionVsTheta_FromParamFits(selectedPid, "ALL", pValues, thetaMinALL, thetaMaxALL, nThetaPoints, baseDir);
}

void PlotMomentumCorrectionVsPhi_MultiP_FromParamVsPhiPoints(
    const int selectedPid,
    const std::string& selecteddetector,
    const std::vector<double>& pValues,
    const std::string& baseDir = "ProtonMomCorrPhi",
    const std::string& outTag = "")
{
    std::string det = DetNorm(selecteddetector);
    std::string prefix = ParticleNameSimple(selectedPid) + "_" + det;
    std::string paramVsPhiDir = baseDir + "/ParamVsPhi";

    void* dirp = gSystem->OpenDirectory(paramVsPhiDir.c_str());
    if (!dirp) {
        std::cerr << "[PlotMomentumCorrectionVsPhi_MultiP_FromParamVsPhiPoints] Cannot open directory: "
                  << paramVsPhiDir << "\n";
        return;
    }

    std::vector<std::string> thetaDirs;
    while (const char* entry = gSystem->GetDirEntry(dirp)) {
        std::string name(entry);
        if (name == "." || name == "..") continue;
        if (name.rfind("T", 0) != 0) continue;
        thetaDirs.push_back(name);
    }
    gSystem->FreeDirectory(dirp);
    std::sort(thetaDirs.begin(), thetaDirs.end());

    if (thetaDirs.empty()) {
        std::cerr << "[PlotMomentumCorrectionVsPhi_MultiP_FromParamVsPhiPoints] No theta directories found under "
                  << paramVsPhiDir << "\n";
        return;
    }

    struct CurveInfo {
        std::string thetaTag;
        std::string thetaLabel;
        std::vector<double> phi;
        std::vector<double> deltaP;
        std::vector<double> deltaPErr;
    };

    std::vector<std::vector<CurveInfo>> curvesByTheta(thetaDirs.size(),
                                                      std::vector<CurveInfo>(pValues.size()));
    double phiMin =  std::numeric_limits<double>::infinity();
    double phiMax = -std::numeric_limits<double>::infinity();

    for (std::size_t itheta = 0; itheta < thetaDirs.size(); ++itheta) {
        const std::string& thetaTag = thetaDirs[itheta];
        const std::string thetaDir = paramVsPhiDir + "/" + thetaTag;
        const std::string fileBase = thetaDir + "/" + prefix + "_" + thetaTag;

        XYErrPoints aPts, bPts, cPts;
        bool okA = ReadPointFile3Cols(fileBase + "_A_p_vs_phi_points.txt", aPts);
        bool okB = ReadPointFile3Cols(fileBase + "_B_p_vs_phi_points.txt", bPts);
        bool okC = ReadPointFile3Cols(fileBase + "_C_p_vs_phi_points.txt", cPts);
        if (!(okA && okB && okC)) continue;

        const std::size_t nPts = std::min({aPts.x.size(), bPts.x.size(), cPts.x.size(),
                                           aPts.y.size(), bPts.y.size(), cPts.y.size(),
                                           aPts.err.size(), bPts.err.size(), cPts.err.size()});
        if (nPts < 2) continue;

        for (std::size_t ip = 0; ip < pValues.size(); ++ip) {
            const double p = pValues[ip];
            CurveInfo& curve = curvesByTheta[itheta][ip];
            curve.thetaTag = thetaTag;
            curve.thetaLabel = ThetaDirToLabel(thetaTag);

            for (std::size_t i = 0; i < nPts; ++i) {
                const double phi = aPts.x[i];
                const double A = aPts.y[i];
                const double B = bPts.y[i];
                const double C = cPts.y[i];
                const double AErr = aPts.err[i];
                const double BErr = bPts.err[i];
                const double CErr = cPts.err[i];

                double deltaP = 0.0;
                double deltaPErr = 0.0;
                if (IsProtonFD(selectedPid, det)) {
                    if (p <= 0.0) continue;
                    deltaP = A + B / p + C / (p * p);
                    deltaPErr = std::sqrt(AErr * AErr +
                                          (BErr / p) * (BErr / p) +
                                          (CErr / (p * p)) * (CErr / (p * p)));
                } else {
                    deltaP = A + B * p + C * p * p;
                    deltaPErr = std::sqrt(AErr * AErr +
                                          (BErr * p) * (BErr * p) +
                                          (CErr * p * p) * (CErr * p * p));
                }

                curve.phi.push_back(phi);
                curve.deltaP.push_back(deltaP);
                curve.deltaPErr.push_back(deltaPErr);

                phiMin = std::min(phiMin, phi);
                phiMax = std::max(phiMax, phi);
            }
        }
    }

    if (!std::isfinite(phiMin) || !std::isfinite(phiMax)) {
        std::cerr << "[PlotMomentumCorrectionVsPhi_MultiP_FromParamVsPhiPoints] No valid phi points found.\n";
        return;
    }

    gSystem->Exec(("mkdir -p " + baseDir).c_str());
    std::string safeTag = outTag.empty() ? "" : ("_" + outTag);

    for (std::size_t itheta = 0; itheta < thetaDirs.size(); ++itheta) {
        bool hasAnyCurve = false;
        for (std::size_t ip = 0; ip < pValues.size(); ++ip) {
            if (!curvesByTheta[itheta][ip].phi.empty()) {
                hasAnyCurve = true;
                break;
            }
        }
        if (!hasAnyCurve) continue;

        auto *c = new TCanvas(
            Form("c_deltaP_vs_phi_multiP_%s_%s_%s",
                 ParticleNameSimple(selectedPid).c_str(),
                 det.c_str(),
                 thetaDirs[itheta].c_str()),
            "",
            2200, 1600
        );
        c->Divide(2, 2, 0.002, 0.002);

        std::vector<TGraphErrors*> keepAliveGraphs;
        std::vector<TH1F*> frames;
        std::vector<TLine*> zeroLines;

        for (std::size_t ip = 0; ip < pValues.size(); ++ip) {
            c->cd(ip + 1);
            gPad->SetMargin(0.10, 0.04, 0.11, 0.08);
            gPad->SetTicks(1, 1);
            gPad->SetGridx();
            gPad->SetGridy();

            auto *frame = new TH1F(
                Form("frame_deltaP_phi_%s_%zu", thetaDirs[itheta].c_str(), ip),
                "",
                100,
                phiMin - 5.0,
                phiMax + 5.0
            );
            frames.push_back(frame);
            frame->SetMinimum(-0.05);
            frame->SetMaximum(0.05);
            frame->SetTitle(Form("%s, p = %.2f GeV;#phi [deg];#delta p [GeV]",
                                 ThetaDirToLabel(thetaDirs[itheta]).c_str(),
                                 pValues[ip]));
            frame->GetXaxis()->SetTitleSize(0.050);
            frame->GetYaxis()->SetTitleSize(0.050);
            frame->GetXaxis()->SetLabelSize(0.040);
            frame->GetYaxis()->SetLabelSize(0.040);
            frame->GetYaxis()->SetTitleOffset(1.0);
            frame->Draw();

            TLine *z = new TLine(phiMin - 5.0, 0.0, phiMax + 5.0, 0.0);
            zeroLines.push_back(z);
            z->SetLineStyle(2);
            z->SetLineColor(kGray + 2);
            z->Draw("SAME");

            const auto& curve = curvesByTheta[itheta][ip];
            if (!curve.phi.empty()) {
                auto *g = new TGraphErrors(curve.phi.size());
                keepAliveGraphs.push_back(g);

                for (std::size_t i = 0; i < curve.phi.size(); ++i) {
                    g->SetPoint(i, curve.phi[i], curve.deltaP[i]);
                    g->SetPointError(i, 0.0, curve.deltaPErr[i]);
                }

                g->SetLineColor(kBlue + 2);
                g->SetMarkerColor(kBlue + 2);
                g->SetLineWidth(2);
                g->SetMarkerStyle(20);
                g->SetMarkerSize(0.9);
                g->Draw("PL SAME");
            }

            TLatex lat;
            lat.SetNDC(true);
            lat.SetTextSize(0.035);
            lat.SetTextColor(kGray + 2);
            lat.DrawLatex(0.12, 0.92, Form("%s, %s",
                                           ParticleNameSimple(selectedPid).c_str(),
                                           det.c_str()));
        }

        c->cd();
        c->Modified();
        c->Update();

        std::string outBase = baseDir + "/" + prefix + "_" + thetaDirs[itheta] +
                              "_deltaP_vs_phi_multiP_4panel" + safeTag;
        c->SaveAs((outBase + ".png").c_str());
        c->SaveAs((outBase + ".pdf").c_str());
        std::cout << "Saved: " << outBase << ".png\nSaved: " << outBase << ".pdf" << std::endl;

        for (auto *g : keepAliveGraphs) delete g;
        for (auto *f : frames) delete f;
        for (auto *z : zeroLines) delete z;
        delete c;
    }
}

void GeneratePiecewiseProtonCorrectionFunctionFromParamVsPhi(
    const std::string& baseDir = "ProtonMomCorrPhi36",
    const std::string& outputFile = "GeneratedPiecewiseProtonCorrection_RunDVCSAnalysis.inc",
    const std::string& functionName = "AddProtonFDPiecewiseCorrection_FromProtonMomCorrPhi36",
    const std::string& dataconfigGuard = "rgkfa18_7546")
{
    const int selectedPid = 2212;
    const std::string det = "FD";
    const std::string prefix = ParticleNameSimple(selectedPid) + "_" + det;
    const std::string paramVsPhiDir = baseDir + "/ParamVsPhi";

    void* dirp = gSystem->OpenDirectory(paramVsPhiDir.c_str());
    if (!dirp) {
        std::cerr << "[GeneratePiecewiseProtonCorrectionFunctionFromParamVsPhi] Cannot open directory: "
                  << paramVsPhiDir << "\n";
        return;
    }

    std::vector<std::string> thetaDirs;
    while (const char* entry = gSystem->GetDirEntry(dirp)) {
        std::string name(entry);
        if (name == "." || name == "..") continue;
        if (name.rfind("T", 0) != 0) continue;
        thetaDirs.push_back(name);
    }
    gSystem->FreeDirectory(dirp);
    std::sort(thetaDirs.begin(), thetaDirs.end());

    if (thetaDirs.empty()) {
        std::cerr << "[GeneratePiecewiseProtonCorrectionFunctionFromParamVsPhi] No theta directories found under "
                  << paramVsPhiDir << "\n";
        return;
    }

    std::ofstream fout(outputFile);
    if (!fout) {
        std::cerr << "[GeneratePiecewiseProtonCorrectionFunctionFromParamVsPhi] Cannot write file: "
                  << outputFile << "\n";
        return;
    }

    fout << "#include <iostream>\n";
    fout << "#include <memory>\n";
    fout << "#include <cmath>\n\n";
    fout << "// Auto-generated from " << baseDir << "/ParamVsPhi\n";
    fout << "// One AddPiecewiseCorrection per (phi, theta) bin; no theta/phi fit is used.\n\n";
    fout << "void " << functionName << "(const std::shared_ptr<MomentumCorrection>& corr) {\n";
    if (!dataconfigGuard.empty()) {
        fout << "  std::cout << \"Applying piecewise proton momentum correction from " << baseDir << "\" << std::endl;\n";
    }

    fout << std::fixed << std::setprecision(6);

    int nRegions = 0;
    for (const std::string& thetaTag : thetaDirs) {
        const std::string thetaDir = paramVsPhiDir + "/" + thetaTag;
        const std::string fileBase = thetaDir + "/" + prefix + "_" + thetaTag;

        XYErrPoints aPts, bPts, cPts;
        bool okA = ReadPointFile3Cols(fileBase + "_A_p_vs_phi_points.txt", aPts);
        bool okB = ReadPointFile3Cols(fileBase + "_B_p_vs_phi_points.txt", bPts);
        bool okC = ReadPointFile3Cols(fileBase + "_C_p_vs_phi_points.txt", cPts);
        if (!(okA && okB && okC)) continue;

        const std::size_t nPts = std::min({aPts.x.size(), bPts.x.size(), cPts.x.size(),
                                           aPts.y.size(), bPts.y.size(), cPts.y.size()});
        if (nPts == 0) continue;

        double thetaMinDeg = 0.0;
        double thetaMaxDeg = 180.0;
        if (!ThetaDirToBoundsDeg(thetaTag, thetaMinDeg, thetaMaxDeg)) continue;

        std::vector<double> phiMids(aPts.x.begin(), aPts.x.begin() + nPts);
        const auto phiEdges = MidpointsToBinEdges(phiMids);

        for (std::size_t i = 0; i < nPts; ++i) {
            const double phiMinDeg = phiEdges[i].first;
            const double phiMaxDeg = phiEdges[i].second;
            const double A = aPts.y[i];
            const double B = bPts.y[i];
            const double C = cPts.y[i];

            fout << "  corr->AddPiecewiseCorrection(  // " << thetaTag
                 << ", phi [" << phiMinDeg << ", " << phiMaxDeg << ") deg\n";
            fout << "      2212, MomentumCorrection::RegionWithDetector{0.0, 10.0, "
                 << thetaMinDeg << " * M_PI / 180, " << thetaMaxDeg << " * M_PI / 180, "
                 << phiMinDeg << " * M_PI / 180, " << phiMaxDeg << " * M_PI / 180, MomentumCorrection::FD}, "
                 << "[](double p, double theta, double phi) {\n";
            fout << "        const double A_p = " << A << ";\n";
            fout << "        const double B_p = " << B << ";\n";
            fout << "        const double C_p = " << C << ";\n";
            fout << "        return p + (A_p + B_p / p + C_p / (p * p));\n";
            fout << "      });\n";
            ++nRegions;
        }
        fout << "\n";
    }

    fout << "}\n";
    fout.close();

    std::cout << "[GeneratePiecewiseProtonCorrectionFunctionFromParamVsPhi] Wrote " << nRegions
              << " piecewise regions to " << outputFile << std::endl;
}

int GetThetaRegionIndex(float thetaDeg, const std::vector<float> &thetaCuts) {
    for (size_t i = 0; i < thetaCuts.size(); ++i)
        if (thetaDeg <= thetaCuts[i]) return i;
    return thetaCuts.size();
}

std::string GetThetaBinLabel(int ti, const std::vector<float>& thetaCuts) {
    char buf[64];
    if (ti == 0)
        sprintf(buf, "#theta #leq %.1f#circ", thetaCuts[0]);
    else if (ti == (int)thetaCuts.size())
        sprintf(buf, "#theta > %.1f#circ", thetaCuts.back());
    else
        sprintf(buf, "%.1f#circ < #theta #leq %.1f#circ", thetaCuts[ti-1], thetaCuts[ti]);
    return std::string(buf);
}

std::string GetDetectorPart(short status) {
    short absStatus = std::abs(status);
    if (absStatus >= 1000 && absStatus < 2000) return "FT";
    else if (absStatus >= 2000 && absStatus < 3000)  return "FD";
    else if (absStatus >= 4000 && absStatus < 5000) return "CD";
    else return "Unknown";
}

TGraph* MakePeakGraph(TH2D* hist, const std::string& outTxtPath)
{
    const int minEntries = 100;
    const double maxChi2Ndf = 100.0;
    const double maxerrPeak = 0.2;

    std::vector<double> xVec, yVec, yerrVec, sigmaVec, chi2ndfVec;

    int nxbins = hist->GetNbinsX();
    for (int ix = 1; ix <= nxbins; ++ix) {
        std::unique_ptr<TH1D> proj(hist->ProjectionY(Form("proj_y_%d", ix), ix, ix, "e"));
        if (!proj || proj->GetEntries() < minEntries) continue;

        int maxBin = proj->GetMaximumBin();
        double mu0 = proj->GetBinCenter(maxBin);
        double sigma0 = proj->GetRMS();
        if (sigma0 == 0) continue;

        TF1 gaus("gaus", "gaus", mu0 - 2*sigma0, mu0 + 2*sigma0);
        gaus.SetParameters(proj->GetMaximum(), mu0, sigma0);
        int fitStatus = proj->Fit(&gaus, "QNR");

        bool fitOK = (fitStatus == 0) &&
                     (gaus.GetParameter(2) > 0) &&
                     (std::abs(gaus.GetParError(1)) < maxerrPeak) &&
                     (gaus.GetChisquare() / std::max(1.0, static_cast<double>(gaus.GetNDF())) < maxChi2Ndf);

        if (!fitOK) continue;

        double xC = hist->GetXaxis()->GetBinCenter(ix);
        double yPeak = gaus.GetParameter(1);
        double peakError = gaus.GetParError(1);
        double sigma = gaus.GetParameter(2);
        double chi2ndf = gaus.GetChisquare() / std::max(1.0, static_cast<double>(gaus.GetNDF()));

        xVec.push_back(xC);
        yVec.push_back(yPeak);
        yerrVec.push_back(peakError);
        sigmaVec.push_back(sigma);
        chi2ndfVec.push_back(chi2ndf);
    }

    // 删除边界点
    if (xVec.size() > 3) {
        xVec.erase(xVec.begin());        xVec.erase(xVec.begin());
        yVec.erase(yVec.begin());        yVec.erase(yVec.begin());
        yerrVec.erase(yerrVec.begin()); yerrVec.erase(yerrVec.begin());
        sigmaVec.erase(sigmaVec.begin()); sigmaVec.erase(sigmaVec.begin());
        chi2ndfVec.erase(chi2ndfVec.begin()); chi2ndfVec.erase(chi2ndfVec.begin());

        xVec.pop_back();    yVec.pop_back();
        yerrVec.pop_back();
        sigmaVec.pop_back(); chi2ndfVec.pop_back();
    }


    // 写入 TXT 文件
    std::ofstream fout(outTxtPath);
    fout << "# xCenter  yPeak yerr sigma  chi2/ndf\n";
    for (size_t i = 0; i < xVec.size(); ++i)
        fout << xVec[i] << "\t" << yVec[i] << "\t" << yerrVec[i]<<"\t" << sigmaVec[i] << "\t" << chi2ndfVec[i] << "\n";
    fout.close();
/*
    // 生成 TGraph
    TGraph* g = new TGraph(xVec.size(), xVec.data(), yVec.data());
    g->SetLineColor(kRed);
    g->SetLineWidth(3);
    return g;
*/
    auto *gErr = new TGraphErrors(xVec.size());
    for (size_t i = 0; i < xVec.size(); ++i) {
        double ex = 0.5 * hist->GetXaxis()->GetBinWidth(               // Bin half-width
                         hist->GetXaxis()->FindBin(xVec[i]));
        double ey = yerrVec[i];                                       // σ from fit
        gErr->SetPoint(i, xVec[i], yVec[i]);           // (xi , yi)
        gErr->SetPointError(i, ex, ey);                // (Δx , Δy)
    }

    // 样式：点 + 误差条
    gErr->SetMarkerStyle(20);
    gErr->SetMarkerSize(0.9);
    gErr->SetLineWidth(1);
    gErr->SetLineColor(kRed);
    gErr->SetMarkerColor(kRed);

    return gErr;
}


//================ Core =================
void DrawParticleKinematicsByThetaBins(
    const int &selectedPid,
    const std::vector<float> &thetaCuts,
    const std::string selecteddetector,                                        // in degree
    const std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> &plotVars,     // {name, nbins, xmin, xmax}
    const std::string &filename,
    const std::string &treename)
    {
    TStopwatch timer;
    timer.Start();

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df(treename, filename);
    //gStyle->SetOptStat(1110);

    struct VarInfo {
        std::string saveName;
        std::string title;  // for future use
        std::string name;
        int nbins;
        double xmin, xmax;
        TH1D* overall;
        std::vector<TH1D*> binHists;   // size = thetaCuts.size()+1
    };

    //==== create histograms ====
    std::vector<VarInfo> vars;
    for (auto &cfg : plotVars) {
        VarInfo v;
        v.saveName = std::get<0>(cfg);
        v.title    = std::get<1>(cfg);
        v.name  = std::get<2>(cfg);
        v.nbins = std::get<3>(cfg);
        v.xmin  = std::get<4>(cfg);
        v.xmax  = std::get<5>(cfg);
        v.overall = new TH1D((v.name+"_overall").c_str(), "", v.nbins, v.xmin, v.xmax);
        v.overall->SetDirectory(nullptr);
        v.binHists.resize(thetaCuts.size()+1);
        for (size_t ti=0; ti<=thetaCuts.size(); ++ti) {
            std::string hname = v.name + Form("_T%zu", ti);
            v.binHists[ti] = new TH1D(hname.c_str(), "", v.nbins, v.xmin, v.xmax);
            v.binHists[ti]->SetDirectory(nullptr);
        }
        vars.push_back(v);
    }

    //==== fill histograms ====
    std::vector<std::string> colNames = {"REC_Particle_pid","REC_Particle_theta","REC_Particle_phi","REC_Particle_p","REC_Particle_status","REC_Particle_pass","REC_Particle_vz","REC_Particle_beta"};

    df.Foreach([&](const RVec<int> &pid,
                   const RVec<float> &theta,
                   const RVec<float> &phi,
                   const RVec<float> &p,
                   const RVec<short> &status,
                   const RVec<int> &passPar,
                   const RVec<float> &vz,
                   const RVec<float> &beta
                   )  // passFid is not used but included for consistency
    {
        for (size_t i=0;i<pid.size();++i) {
            if (pid[i]!=selectedPid) continue;
            if (!passPar[i]) continue;  // only fill if the particle passes the fiducial cuts
            if (GetDetectorPart(status[i])!=selecteddetector && selecteddetector!="ALL") continue;  // check detector part
            float thetaDeg = theta[i]*180.0/M_PI;
            int ti = GetThetaRegionIndex(thetaDeg, thetaCuts);
            for (auto &v: vars) {
                double value = NAN;
                if (v.name=="theta") value = theta[i]*180.0/M_PI;
                else if (v.name=="phi"){value = phi[i]*180.0/M_PI; if (value < 0) value += 360.0;} // ensure phi is in [0, 360)
                else if (v.name=="p") value = p[i];
                else if (v.name=="vz") value = vz[i];  // assuming vz is the same as p for this example
                else if (v.name=="beta") value = beta[i];
                if (std::isnan(value)) continue;
                v.binHists[ti]->Fill(value);
                v.overall->Fill(value);
            }
        }
    }, colNames);

    //==== output directory ====
    std::string outDir = "ParticleKinematicPlots";
    gSystem->Exec(("mkdir -p "+outDir).c_str());

    //==== draw and save ====
    for (auto &v: vars) {
        // overall
        {
            TCanvas *c = new TCanvas(("c_"+v.name+"_overall").c_str(),"",1600,1200);
            v.overall->SetTitle((selecteddetector + " " +v.title + " " + v.name).c_str());
            if (v.name == "phi" || v.name == "theta") v.overall->GetXaxis()->SetTitle((v.name+" [deg]").c_str());
            else v.overall->GetXaxis()->SetTitle((v.name+" GeV").c_str());
            v.overall->GetYaxis()->SetTitle("Counts");
            v.overall->Draw();
            std::string out = outDir + "/" + v.saveName + "_" +selecteddetector+ "_" + v.name + "_overall.png";
            c->SaveAs(out.c_str());
            delete c;
            std::cout << "Saved: " << out << std::endl;
        }
        // theta bins
        for (size_t ti=0; ti<=thetaCuts.size(); ++ti) {
            TCanvas *c = new TCanvas(("c_"+v.name+Form("_T%zu",ti)).c_str(),"",1600,1200);
            if (v.name == "phi" || v.name == "theta") v.binHists[ti]->GetXaxis()->SetTitle((v.name+ " [deg]").c_str());
            else v.binHists[ti]->GetXaxis()->SetTitle((v.name+" [GeV]").c_str()); 
            v.binHists[ti]->SetTitle((selecteddetector + " " + v.title+" "+v.name + " in " + GetThetaBinLabel(ti, thetaCuts)).c_str());
            v.binHists[ti]->GetYaxis()->SetTitle("Counts");
            v.binHists[ti]->Draw();
            std::string out = outDir + "/" + v.saveName + "_" + selecteddetector + "_" + v.name + Form("_T%zu.png",ti);
            c->SaveAs(out.c_str());
            delete c;
            std::cout << "Saved: " << out << std::endl;
        }
    }

    timer.Stop();
    std::cout << "Time for DrawParticleKinematicsByThetaBins: ";
    timer.Print();
}

void Draw2DParticleKinematicsByThetaBins(
    const int &selectedPid,
    const std::vector<float> &thetaCuts,
    const std::string selecteddetector,
    const std::vector<std::tuple<std::string,std::string,std::string,int,double,double,int,double,double>> &plotVars,
    const std::string &filename,
    const std::string &treename) {
    TStopwatch timer;
    timer.Start();

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df(treename, filename);

    struct Var2DInfo {
        std::string saveName;
        std::string title;
        std::string name; // format: x:y (e.g., "p:theta")
        int nxbins;
        double xmin, xmax;
        int nybins;
        double ymin, ymax;
        TH2D* overall;
        std::vector<TH2D*> binHists;
    };

    std::vector<Var2DInfo> vars;
    for (auto &cfg : plotVars) {
        Var2DInfo v;
        v.saveName = std::get<0>(cfg);
        v.title    = std::get<1>(cfg);
        v.name     = std::get<2>(cfg);
        v.nxbins   = std::get<3>(cfg);
        v.xmin     = std::get<4>(cfg);
        v.xmax     = std::get<5>(cfg);
        v.nybins   = std::get<6>(cfg);
        v.ymin     = std::get<7>(cfg);
        v.ymax     = std::get<8>(cfg);

        v.overall = new TH2D((v.saveName+"_overall").c_str(), "", v.nxbins, v.xmin, v.xmax, v.nybins, v.ymin, v.ymax);
        v.overall->SetDirectory(nullptr);

        v.binHists.resize(thetaCuts.size() + 1);
        for (size_t ti = 0; ti <= thetaCuts.size(); ++ti) {
            std::string hname = v.saveName + Form("_T%zu", ti);
            v.binHists[ti] = new TH2D(hname.c_str(), "", v.nxbins, v.xmin, v.xmax, v.nybins, v.ymin, v.ymax);
            v.binHists[ti]->SetDirectory(nullptr);
        }
        vars.push_back(v);
    }

    std::vector<std::string> colNames = {"REC_Particle_pid", "REC_Particle_theta", "REC_Particle_phi", "REC_Particle_p", "REC_Particle_status", "REC_Particle_pass"};

    df.Foreach([&](const RVec<int> &pid,
                   const RVec<float> &theta,
                   const RVec<float> &phi,
                   const RVec<float> &p,
                   const RVec<short> &status,
                   const RVec<int> &passPar)
    {
        for (size_t i = 0; i < pid.size(); ++i) {
            if (pid[i] != selectedPid) continue;
            if (!passPar[i]) continue;
            if (GetDetectorPart(status[i]) != selecteddetector && selecteddetector != "ALL") continue;

            float thetaDeg = theta[i] * 180.0 / M_PI;
            float phiDeg = phi[i] * 180.0 / M_PI;
            if (phiDeg < 0) phiDeg += 360.0;
            int ti = GetThetaRegionIndex(thetaDeg, thetaCuts);

            for (auto &v : vars) {
                double xval = NAN, yval = NAN;
                if (v.name == "p:theta") { xval = p[i]; yval = thetaDeg; }
                else if (v.name == "theta:p") { xval = thetaDeg; yval = p[i]; }
                else if (v.name == "phi:theta") { xval = phiDeg; yval = thetaDeg; }
                else if (v.name == "theta:phi") { xval = thetaDeg; yval = phiDeg; }
                else continue;

                if (std::isnan(xval) || std::isnan(yval)) continue;
                v.binHists[ti]->Fill(xval, yval);
                v.overall->Fill(xval, yval);
            }
        }
    }, colNames);

    std::string outDir = "Particle2DKinematicPlots";
    gSystem->Exec(("mkdir -p " + outDir).c_str());

    for (auto &v : vars) {
        {
            TCanvas *c = new TCanvas(("c_" + v.saveName + "_overall").c_str(), "", 1600, 1200);
            v.overall->SetTitle((selecteddetector + " " + v.title).c_str());
            v.overall->Draw("COLZ");
            if (v.name == "p:theta") {
                v.overall->GetXaxis()->SetTitle("p [GeV]");
                v.overall->GetYaxis()->SetTitle("#theta [deg]");
            } else if (v.name == "theta:p") {
                v.overall->GetXaxis()->SetTitle("#theta [deg]");
                v.overall->GetYaxis()->SetTitle("p [GeV]");
            } else if (v.name == "phi:theta") {
                v.overall->GetXaxis()->SetTitle("#phi [deg]");
                v.overall->GetYaxis()->SetTitle("#theta [deg]");
            } else if (v.name == "theta:phi") {
                v.overall->GetXaxis()->SetTitle("#theta [deg]");
                v.overall->GetYaxis()->SetTitle("#phi [deg]");
            }
            std::string out = outDir + "/" + v.saveName + "_" + selecteddetector + "_overall.png";
            c->SaveAs(out.c_str());
            delete c;
            std::cout << "Saved: " << out << std::endl;
        }

        for (size_t ti = 0; ti <= thetaCuts.size(); ++ti) {
            TCanvas *c = new TCanvas(("c_" + v.saveName + Form("_T%zu", ti)).c_str(), "", 1600, 1200);
            v.binHists[ti]->SetTitle((selecteddetector + " " + v.title + " in " + GetThetaBinLabel(ti, thetaCuts)).c_str());
            v.binHists[ti]->Draw("COLZ");
            if (v.name == "p:theta"){
                v.binHists[ti]->GetXaxis()->SetTitle("p [GeV]");
                v.binHists[ti]->GetYaxis()->SetTitle("#theta [deg]");
            } else if (v.name == "theta:p") {
                v.binHists[ti]->GetXaxis()->SetTitle("#theta [deg]");
                v.binHists[ti]->GetYaxis()->SetTitle("p [GeV]");
            } else if (v.name == "phi:theta") {
                v.binHists[ti]->GetXaxis()->SetTitle("#phi [deg]");
                v.binHists[ti]->GetYaxis()->SetTitle("#theta [deg]");
            } else if (v.name == "theta:phi") {
                v.binHists[ti]->GetXaxis()->SetTitle("#theta [deg]");
                v.binHists[ti]->GetYaxis()->SetTitle("#phi [deg]");
            }
            std::string out = outDir + "/" + v.saveName + "_" + selecteddetector + Form("_T%zu.png", ti);
            c->SaveAs(out.c_str());
            delete c;
            std::cout << "Saved: " << out << std::endl;
        }
    }

    timer.Stop();
    std::cout << "Time for Draw2DParticleKinematicsByThetaBins: ";
    timer.Print();
}

void DrawMCParticleKinematicsByThetaBins(
    const int &selectedPid,
    const std::vector<float> &thetaCuts,                                     // in degree
    const std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> &plotVars,     // {name, nbins, xmin, xmax}
    const std::string &filename,
    const std::string &treename)
    {
    TStopwatch timer;
    timer.Start();

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df(treename, filename);
    //gStyle->SetOptStat(1110);

    struct VarInfo {
        std::string saveName;
        std::string title;  // for future use
        std::string name;
        int nbins;
        double xmin, xmax;
        TH1D* overall;
        std::vector<TH1D*> binHists;   // size = thetaCuts.size()+1
    };

    //==== create histograms ====
    std::vector<VarInfo> vars;
    for (auto &cfg : plotVars) {
        VarInfo v;
        v.saveName = std::get<0>(cfg);
        v.title    = std::get<1>(cfg);
        v.name  = std::get<2>(cfg);
        v.nbins = std::get<3>(cfg);
        v.xmin  = std::get<4>(cfg);
        v.xmax  = std::get<5>(cfg);
        v.overall = new TH1D((v.name+"_overall").c_str(), "", v.nbins, v.xmin, v.xmax);
        v.overall->SetDirectory(nullptr);
        v.binHists.resize(thetaCuts.size()+1);
        for (size_t ti=0; ti<=thetaCuts.size(); ++ti) {
            std::string hname = v.name + Form("_T%zu", ti);
            v.binHists[ti] = new TH1D(hname.c_str(), "", v.nbins, v.xmin, v.xmax);
            v.binHists[ti]->SetDirectory(nullptr);
        }
        vars.push_back(v);
    }

    //==== fill histograms ====
    std::vector<std::string> colNames = {"MC_Lund_pid","MC_Lund_px","MC_Lund_py","MC_Lund_pz"};

    df.Foreach([&](const RVec<int> &pid,
                   const RVec<float> &px,
                   const RVec<float> &py,
                   const RVec<float> &pz)  // passFid is not used but included for consistency
    {
        for (size_t i=0;i<pid.size();++i) {
            if (pid[i]!=selectedPid) continue;
            float p = std::sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
            float phi = std::atan2(py[i], px[i]);
            if (phi < 0) phi += 2 * M_PI;
            float theta = std::atan2(std::sqrt(px[i]*px[i] + py[i]*py[i]), pz[i]);
            float thetaDeg = theta*180.0/M_PI;
            int ti = GetThetaRegionIndex(thetaDeg, thetaCuts);
            for (auto &v: vars) {
                double value = NAN;
                if (v.name=="theta") value = theta*180.0/M_PI;
                else if (v.name=="phi") value = phi*180.0/M_PI; // ensure phi is in [0, 360)
                else if (v.name=="p") value = p;
                if (std::isnan(value)) continue;
                v.binHists[ti]->Fill(value);
                v.overall->Fill(value);
            }
        }
    }, colNames);

    //==== output directory ====
    std::string outDir = "MCParticleKinematicPlots";
    gSystem->Exec(("mkdir -p "+outDir).c_str());

    //==== draw and save ====
    for (auto &v: vars) {
        // overall
        {
            TCanvas *c = new TCanvas(("c_"+v.name+"_overall").c_str(),"",1600,1200);
            v.overall->SetTitle((v.title + " " + v.name).c_str());
            if (v.name == "phi" || v.name == "theta") v.overall->GetXaxis()->SetTitle((v.name+" [deg]").c_str());
            else v.overall->GetXaxis()->SetTitle((v.name+" GeV").c_str());
            v.overall->GetYaxis()->SetTitle("Counts");
            v.overall->Draw();
            std::string out = outDir + "/" + v.saveName + "_" + v.name + "_overall.png";
            c->SaveAs(out.c_str());
            delete c;
            std::cout << "Saved: " << out << std::endl;
        }
        // theta bins
        for (size_t ti=0; ti<=thetaCuts.size(); ++ti) {
            TCanvas *c = new TCanvas(("c_"+v.name+Form("_T%zu",ti)).c_str(),"",1600,1200);
            if (v.name == "phi" || v.name == "theta") v.binHists[ti]->GetXaxis()->SetTitle((v.name+ " [deg]").c_str());
            else v.binHists[ti]->GetXaxis()->SetTitle((v.name+" [GeV]").c_str()); 
            v.binHists[ti]->SetTitle((v.title+" "+v.name + " in " + GetThetaBinLabel(ti, thetaCuts)).c_str());
            v.binHists[ti]->GetYaxis()->SetTitle("Counts");
            v.binHists[ti]->Draw();
            std::string out = outDir + "/" + v.saveName + "_" + v.name + Form("_T%zu.png",ti);
            c->SaveAs(out.c_str());
            delete c;
            std::cout << "Saved: " << out << std::endl;
        }
    }

    timer.Stop();
    std::cout << "Time for DrawParticleKinematicsByThetaBins: ";
    timer.Print();
}

void DrawDeltaPByThetaBins(
    const int &selectedPid,
    const std::vector<float> &thetaCuts,
    const std::string selecteddetector,
    const std::vector<std::tuple<std::string,std::string,std::string,int,double,double,int,double,double>> &plotVars,
    const std::string &filename,
    const std::string &treename, const std::string& outDir ="ProtonMomCorr", bool isOutBend = false) {
    TStopwatch timer;
    timer.Start();

    ROOT::EnableImplicitMT();

    std::vector<double> thetaMidVec, aVec, bVec, cVec;
    std::vector<double> aErrVec, bErrVec, cErrVec;

    auto GetParticleName = [](int pid) -> std::string {
        if (pid == 11) return "electron";
        if (pid == 22) return "photon";
        if (pid == 2212) return "proton";
        return "pid" + std::to_string(pid);
    };
    std::string prefix = GetParticleName(selectedPid) + "_" + selecteddetector;


    ROOT::RDataFrame df(treename, filename);
    gStyle->SetOptStat(0);

    struct Var2DInfo {
        std::string saveName;
        std::string title;
        std::string name; // format: x:y (e.g., "deltaP:p")
        int nxbins;
        double xmin, xmax;
        int nybins;
        double ymin, ymax;
        TH2D* overall;
        std::vector<TH2D*> binHists;
    };

    std::vector<Var2DInfo> vars;
    for (auto &cfg : plotVars) {
        Var2DInfo v;
        v.saveName = std::get<0>(cfg);
        v.title    = std::get<1>(cfg);
        v.name     = std::get<2>(cfg);
        v.nxbins   = std::get<3>(cfg);
        v.xmin     = std::get<4>(cfg);
        v.xmax     = std::get<5>(cfg);
        v.nybins   = std::get<6>(cfg);
        v.ymin     = std::get<7>(cfg);
        v.ymax     = std::get<8>(cfg);

        v.overall = new TH2D((v.saveName+"_overall").c_str(), "", v.nxbins, v.xmin, v.xmax, v.nybins, v.ymin, v.ymax);
        v.overall->SetDirectory(nullptr);

        v.binHists.resize(thetaCuts.size() + 1);
        for (size_t ti = 0; ti <= thetaCuts.size(); ++ti) {
            std::string hname = v.saveName + Form("_T%zu", ti);
            v.binHists[ti] = new TH2D(hname.c_str(), "", v.nxbins, v.xmin, v.xmax, v.nybins, v.ymin, v.ymax);
            v.binHists[ti]->SetDirectory(nullptr);
        }
        vars.push_back(v);
    }

    std::vector<std::string> colNames = {"REC_Particle_pid", "REC_Particle_theta", "REC_Particle_phi", "REC_Particle_p", "REC_Particle_status", "REC_Particle_pass",
                                         "MC_Particle_pid", "MC_Particle_px", "MC_Particle_py", "MC_Particle_pz"};

    df.Foreach([&](const RVec<int> &pid,
                   const RVec<float> &theta,
                   const RVec<float> &phi,
                   const RVec<float> &p,
                   const RVec<short> &status,
                   const RVec<int> &passPar,
                   const RVec<int> &mcpid,
                   const RVec<float> &mcpx,
                   const RVec<float> &mcpy,
                   const RVec<float> &mcpz)
    {
        std::vector<double> thetaCenters;
        std::vector<double> fitParamA, fitParamB, fitParamC;

        for (size_t i = 0; i < pid.size(); ++i) {
            if (pid[i] != selectedPid) continue;
            if (!passPar[i]) continue;
            if (GetDetectorPart(status[i]) != selecteddetector && selecteddetector != "ALL") continue;

            float thetaDeg = theta[i] * 180.0 / M_PI;
            int ti = GetThetaRegionIndex(thetaDeg, thetaCuts);

            float p_rec = p[i];
            float p_mc = NAN;
            for (size_t j = 0; j < mcpid.size(); ++j) {
                if (mcpid[j] == selectedPid) {
                    p_mc = std::sqrt(mcpx[j]*mcpx[j] + mcpy[j]*mcpy[j] + mcpz[j]*mcpz[j]);
                    break;
                }
            }
            if (std::isnan(p_mc)) continue;
            float deltaP = p_mc - p_rec;

            for (auto &v : vars) {
                double xval = NAN, yval = NAN;
                if (v.name == "deltaP:p") { xval = p_rec; yval = deltaP; }
                else if (v.name == "p:deltaP") { xval = deltaP; yval = p_rec; }
                else continue;

                if (std::isnan(xval) || std::isnan(yval)) continue;
                v.binHists[ti]->Fill(xval, yval);
                v.overall->Fill(xval, yval);
            }
        }
    }, colNames);

    //std::string outDir = "ProtonMomCorr_Fall2018_inb_DeltaPPlots";
    gSystem->Exec(("mkdir -p " + outDir).c_str());

    for (auto &v : vars) {
        {
            TCanvas *c = new TCanvas(("c_" + v.saveName + "_overall").c_str(), "", 3000, 1500);
            v.overall->SetTitle((selecteddetector + " " + v.title).c_str());
            gPad->SetLogz();
            v.overall->Draw("COLZ");
            //TGraph* gPeak = MakePeakGraph(v.overall,outDir + "/" + v.saveName + "_" + selecteddetector + "_overall.txt");
            //gPeak->Draw("PEZ SAME");
            //drawHist(v.overall, selecteddetector + " " + v.title, "p [GeV]", "#delta p [GeV]");
            if (v.name == "deltaP:p") {
                v.overall->GetXaxis()->SetTitle("p [GeV]");
                v.overall->GetYaxis()->SetTitle("#delta p [GeV]");
            } else if (v.name == "p:deltaP") {
                v.overall->GetXaxis()->SetTitle("#delta p [GeV]");
                v.overall->GetYaxis()->SetTitle("p [GeV]");
            }
            std::string out = outDir + "/" + v.saveName + "_" + selecteddetector + "_overall.png";
            c->SaveAs(out.c_str());
            delete c;
            std::cout << "Saved: " << out << std::endl;
        }

        for (size_t ti = 0; ti <= thetaCuts.size(); ++ti) {
            TCanvas *c = new TCanvas(("c_" + v.saveName + Form("_T%zu", ti)).c_str(), "", 3000, 1500);
            v.binHists[ti]->SetTitle((selecteddetector + " " + v.title + " in " + GetThetaBinLabel(ti, thetaCuts)).c_str());
            gPad->SetLogz();
            v.binHists[ti]->Draw("COLZ");
            TGraph* gPeak2 = MakePeakGraph(v.binHists[ti],outDir + "/" + v.saveName + "_" + selecteddetector + Form("_T%zu.txt", ti));
            gPeak2->Draw("PEZ SAME");
  
            TF1* fitFunc = nullptr;
            std::string fitExpr;

            if (selectedPid == 2212) { // Proton
                if (selecteddetector == "FD") {
                    fitExpr = "[0] + [1]/x+[2]/(x*x)";
                } else if (selecteddetector == "CD") {
                    fitExpr = "[0] + [1]*x + [2]*x*x";  // a + bx + cx^2
                } else {
                    std::cerr << "Unknown detector for proton: " << selecteddetector << "\n";
                    continue;
                }
            } else if (selectedPid == 11 || selectedPid == 22) { // Electron or Photon
                fitExpr = "[0] + [1]*x + [2]*x*x"; 
            } else {
                std::cerr << "Unknown pid for fitting: " << selectedPid << "\n";
                fitExpr = "[0] + [1]*x + [2]*x*x"; // fallback
            }

            fitFunc = new TF1("fitFunc", fitExpr.c_str(), v.binHists[ti]->GetXaxis()->GetXmin(), v.binHists[ti]->GetXaxis()->GetXmax());
             if (selecteddetector == "FD"&& isOutBend) {
                    //fitFunc->FixParameter(2, 0.0); // For outbending FD protons, fix the C constant term to 0
            }
            //fitFunc->SetParameters(0.01, -0.01, 0.01); // Initial guesses
            //fitFunc->SetParameters(0.01, -0.01, 0.01); // Initial guesses
            //fitFunc->SetParameters(0.01, -0.01, 0.01); // Initial guesses


            if (gPeak2->GetN() > 2) {
                gPeak2->Fit(fitFunc, "Q");  // Quiet mode

                std::ofstream foutFit(outDir + "/" + v.saveName + "_" + selecteddetector + Form("_T%zu_fit.txt", ti));
                foutFit << "# Fit expression: " << fitExpr << "\n";
                foutFit << "# Fit parameters:\n";
                for (int i = 0; i < fitFunc->GetNpar(); ++i)
                    foutFit << "p" << i << " = " << fitFunc->GetParameter(i)
                            << " ± " << fitFunc->GetParError(i) << "\n";
                foutFit.close();

                fitFunc->SetLineColor(kBlue + 2);
                fitFunc->SetLineWidth(4);
                fitFunc->Draw("SAME");

                // === 记录每个 bin 的 theta center 和 拟合参数 ===
                double thetaMid;
                if (ti == 0)
                    thetaMid = thetaCuts[0];
                else if (ti == thetaCuts.size())
                    thetaMid = thetaCuts.back();
                else
                    thetaMid = 0.5 * (thetaCuts[ti] + thetaCuts[ti - 1]);

                if (fitFunc->GetNpar() >= 3 && ti < thetaCuts.size()-1) {
                    thetaMidVec.push_back(thetaMid);
                    aVec.push_back(fitFunc->GetParameter(0));
                    bVec.push_back(fitFunc->GetParameter(1));
                    cVec.push_back(fitFunc->GetParameter(2));
                    aErrVec.push_back(fitFunc->GetParError(0));
                    bErrVec.push_back(fitFunc->GetParError(1));
                    cErrVec.push_back(fitFunc->GetParError(2));

                }
            }

            //drawHist(v.binHists[ti],selecteddetector + " " + v.title + " in " + GetThetaBinLabel(ti, thetaCuts),"p [GeV]", "#delta p [GeV]");

            if (v.name == "deltaP:p") {
                v.binHists[ti]->GetXaxis()->SetTitle("p [GeV]");
                v.binHists[ti]->GetYaxis()->SetTitle("#delta p [GeV]");
            } else if (v.name == "p:deltaP") {
                v.binHists[ti]->GetXaxis()->SetTitle("#delta p [GeV]");
                v.binHists[ti]->GetYaxis()->SetTitle("p [GeV]");
            }
            std::string out = outDir + "/" + v.saveName + "_" + selecteddetector + Form("_T%zu.png", ti);
            c->SaveAs(out.c_str());
            delete c;
            std::cout << "Saved: " << out << std::endl;
        }
    }
    std::string fitParamOutDir = outDir + "/ParamFits";
    gSystem->Exec(("mkdir -p " + fitParamOutDir).c_str());

    auto PlotParamVsTheta = [&](const std::vector<double>& theta,
                             const std::vector<double>& param,
                             const std::vector<double>& errors,
                             const std::string& pname,
                             const std::string& fitExpr) {
        if (theta.size() < 3) return;

        TCanvas *c = new TCanvas(("c_" + prefix + "_" + pname + "_vs_theta").c_str(), "", 1800, 1200);
        TGraphErrors* g = new TGraphErrors(theta.size());
        for (size_t i = 0; i < theta.size(); ++i) {
            g->SetPoint(i, theta[i], param[i]);
            g->SetPointError(i, 0.0, errors[i]); // X误差设为0，Y误差来自拟合
        }

        TF1* f = new TF1(("fit_" + pname).c_str(), fitExpr.c_str(),
                     *std::min_element(theta.begin(), theta.end()),
                     *std::max_element(theta.begin(), theta.end()));
        g->Fit(f, "Q");

        g->SetTitle((prefix + ": " + pname + " vs #theta").c_str());
        g->GetXaxis()->SetTitle("#theta [deg]");
        g->GetYaxis()->SetTitle((pname + " value").c_str());
         g->GetYaxis()->SetTitleOffset(1.1);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.2);
        g->Draw("AP");

        f->SetLineColor(kRed + 1);
        f->SetLineWidth(2);
        f->Draw("SAME");
        g->GetYaxis()->SetRangeUser(-0.05, 0.05);
        if(prefix =="proton_FD")g->GetYaxis()->SetRangeUser(-0.05, 0.05);
        if(prefix =="proton_CD")g->GetYaxis()->SetRangeUser(-0.05, 0.05);


        std::string imgPath = fitParamOutDir + "/" + prefix + "_" + pname + "_vs_theta.png";
        c->SaveAs(imgPath.c_str());
        delete c;
        std::cout << "Saved: " << imgPath << std::endl;

        std::ofstream fout(fitParamOutDir + "/" + prefix + "_" + pname + "_vs_theta.txt");
        fout << "# Fit expression: " << fitExpr << "\n";
        for (int i = 0; i < f->GetNpar(); ++i)
        fout << "p" << i << " = " << f->GetParameter(i)
             << " ± " << f->GetParError(i) << "\n";
        fout.close();
    };

    if (selecteddetector == "FD") {
         if (isOutBend) {
            PlotParamVsTheta(thetaMidVec, aVec, aErrVec, "A_p", "[0] + [1]*x");
            PlotParamVsTheta(thetaMidVec, bVec, bErrVec, "B_p", "[0] + [1]*x");
            PlotParamVsTheta(thetaMidVec, cVec, cErrVec, "C_p", "[0] + [1]*x");
        }
        else{
            PlotParamVsTheta(thetaMidVec, aVec, aErrVec, "A_p", "[0] + [1]*x");
            PlotParamVsTheta(thetaMidVec, bVec, bErrVec, "B_p", "[0] + [1]*x");
            PlotParamVsTheta(thetaMidVec, cVec, cErrVec, "C_p", "[0] + [1]*x");
        }
   
    } else if (selecteddetector == "CD") {
        PlotParamVsTheta(thetaMidVec, aVec, aErrVec, "A_p", "[0] + [1]*x");
        PlotParamVsTheta(thetaMidVec, bVec, bErrVec, "B_p", "[0] + [1]*x");
        PlotParamVsTheta(thetaMidVec, cVec, cErrVec, "C_p", "[0] + [1]*x");
    }



    timer.Stop();
    std::cout << "Time for DrawDeltaPByThetaBins: ";
    timer.Print();
}

void DrawDeltaPByThetaPhiBinsold(
    const int &selectedPid,
    const std::vector<float> &thetaCuts,
    const std::string selecteddetector,
    const int nPhiBins,
    double phiStartDeg,
    const std::vector<std::tuple<std::string,std::string,std::string,int,double,double,int,double,double>> &plotVars,
    const std::string &filename,
    const std::string &treename,
    const std::string& outDir = "ProtonMomCorrPhi",
    bool isOutBend = false)
{
    TStopwatch timer;
    timer.Start();

    ROOT::EnableImplicitMT();

    auto GetParticleName = [](int pid) -> std::string {
        if (pid == 11)   return "electron";
        if (pid == 22)   return "photon";
        if (pid == 2212) return "proton";
        return "pid" + std::to_string(pid);
    };

    auto WrapToRange = [](double x, double start) -> double {
        double y = std::fmod(x - start, 360.0);
        if (y < 0) y += 360.0;
        return y + start;
    };

    auto GetPhiBinIndex = [&](double phiDeg) -> int {
        double phiWrapped = WrapToRange(phiDeg, phiStartDeg);
        double dphi = 360.0 / nPhiBins;
        int idx = static_cast<int>((phiWrapped - phiStartDeg) / dphi);

        if (idx < 0) idx = 0;
        if (idx >= nPhiBins) idx = nPhiBins - 1;
        return idx;
    };

    auto GetPhiBinLabel = [&](int phiBinIdx) -> std::string {
        double dphi = 360.0 / nPhiBins;
        double phiMin = phiStartDeg + phiBinIdx * dphi;
        double phiMax = phiStartDeg + (phiBinIdx + 1) * dphi;
        char buf[128];
        sprintf(buf, "%.1f#circ #leq #phi < %.1f#circ", phiMin, phiMax);
        return std::string(buf);
    };

    auto GetPhiBinTag = [&](int phiBinIdx) -> std::string {
        double dphi = 360.0 / nPhiBins;
        double phiMin = phiStartDeg + phiBinIdx * dphi;
        double phiMax = phiStartDeg + (phiBinIdx + 1) * dphi;
        char buf[128];
        sprintf(buf, "P%02d_%.0fto%.0f", phiBinIdx, phiMin, phiMax);
        return std::string(buf);
    };

    auto GetParamYRange = [&](const std::string& prefix, const std::string& paramName) -> std::pair<double,double> {
        if (prefix == "proton_FD") {
            if (paramName == "A_p") return {-0.1, 0.1};
            if (paramName == "B_p") return {-0.1, 0.1};
            if (paramName == "C_p") return {-0.1, 0.1};
            return {-0.1, 0.1};
        }
        if (prefix == "proton_CD") {
            if (paramName == "A_p") return {-0.1, 0.1};
            if (paramName == "B_p") return {-0.1, 0.1};
            if (paramName == "C_p") return {-0.1, 0.1};
            return {-0.1, 0.1};
        }
        return {-0.1, 0.1};
    };

    std::string particleDetectorPrefix = GetParticleName(selectedPid) + "_" + selecteddetector;

    struct ParamThetaPlotInfo {
        std::string phiTag;
        std::string phiLabel;
        std::vector<double> theta;
        std::vector<double> param;
        std::vector<double> errors;
        std::string fitFormula;
        bool hasData = false;
    };

    // 改成 map<paramName, vector fixed-size = nPhiBins>
    std::map<std::string, std::vector<ParamThetaPlotInfo>> summaryPadPlotsByParam;

    ROOT::RDataFrame df(treename, filename);
    gStyle->SetOptStat(0);

    struct Var3DInfo {
        std::string saveName;
        std::string title;
        std::string name; // "deltaP:p" or "p:deltaP"
        int nxbins;
        double xmin, xmax;
        int nybins;
        double ymin, ymax;

        TH2D* overall;
        std::vector<TH2D*> phiOverall;
        std::vector<std::vector<TH2D*>> binHists;
    };

    std::vector<Var3DInfo> vars;
    for (auto &cfg : plotVars) {
        Var3DInfo v;
        v.saveName = std::get<0>(cfg);
        v.title    = std::get<1>(cfg);
        v.name     = std::get<2>(cfg);
        v.nxbins   = std::get<3>(cfg);
        v.xmin     = std::get<4>(cfg);
        v.xmax     = std::get<5>(cfg);
        v.nybins   = std::get<6>(cfg);
        v.ymin     = std::get<7>(cfg);
        v.ymax     = std::get<8>(cfg);

        v.overall = new TH2D((v.saveName + "_overall").c_str(), "",
                             v.nxbins, v.xmin, v.xmax, v.nybins, v.ymin, v.ymax);
        v.overall->SetDirectory(nullptr);

        v.phiOverall.resize(nPhiBins, nullptr);
        for (int phiBinIdx = 0; phiBinIdx < nPhiBins; ++phiBinIdx) {
            std::string histName = v.saveName + Form("_P%02d_overall", phiBinIdx);
            v.phiOverall[phiBinIdx] = new TH2D(histName.c_str(), "",
                                               v.nxbins, v.xmin, v.xmax, v.nybins, v.ymin, v.ymax);
            v.phiOverall[phiBinIdx]->SetDirectory(nullptr);
        }

        v.binHists.resize(nPhiBins);
        for (int phiBinIdx = 0; phiBinIdx < nPhiBins; ++phiBinIdx) {
            v.binHists[phiBinIdx].resize(thetaCuts.size() + 1, nullptr);
            for (size_t thetaBinIdx = 0; thetaBinIdx <= thetaCuts.size(); ++thetaBinIdx) {
                std::string histName = v.saveName + Form("_P%02d_T%zu", phiBinIdx, thetaBinIdx);
                v.binHists[phiBinIdx][thetaBinIdx] = new TH2D(histName.c_str(), "",
                                                              v.nxbins, v.xmin, v.xmax, v.nybins, v.ymin, v.ymax);
                v.binHists[phiBinIdx][thetaBinIdx]->SetDirectory(nullptr);
            }
        }

        vars.push_back(v);
    }

    std::vector<std::string> colNames = {
        "REC_Particle_pid", "REC_Particle_theta", "REC_Particle_phi", "REC_Particle_p",
        "REC_Particle_status", "REC_Particle_pass",
        "MC_Particle_pid", "MC_Particle_px", "MC_Particle_py", "MC_Particle_pz"
    };

    df.Foreach([&](const RVec<int> &pid,
                   const RVec<float> &theta,
                   const RVec<float> &phi,
                   const RVec<float> &p,
                   const RVec<short> &status,
                   const RVec<int> &passPar,
                   const RVec<int> &mcpid,
                   const RVec<float> &mcpx,
                   const RVec<float> &mcpy,
                   const RVec<float> &mcpz)
    {
        for (size_t i = 0; i < pid.size(); ++i) {
            if (pid[i] != selectedPid) continue;
            if (!passPar[i]) continue;
            if (GetDetectorPart(status[i]) != selecteddetector && selecteddetector != "ALL") continue;

            double thetaDeg = theta[i] * 180.0 / M_PI;
            double phiDeg = phi[i] * 180.0 / M_PI;
            double phiWrapped = WrapToRange(phiDeg, phiStartDeg);

            int thetaBinIdx = GetThetaRegionIndex(thetaDeg, thetaCuts);
            int phiBinIdx = GetPhiBinIndex(phiDeg);

            double pRec = p[i];
            double pMC  = NAN;

            for (size_t j = 0; j < mcpid.size(); ++j) {
                if (mcpid[j] == selectedPid) {
                    pMC = std::sqrt(mcpx[j]*mcpx[j] + mcpy[j]*mcpy[j] + mcpz[j]*mcpz[j]);
                    break;
                }
            }
            if (std::isnan(pMC)) continue;

            double deltaP = pMC - pRec;

            for (auto &v : vars) {
                double xval = NAN, yval = NAN;

                if (v.name == "deltaP:p") {
                    xval = pRec;
                    yval = deltaP;
                } else if (v.name == "p:deltaP") {
                    xval = deltaP;
                    yval = pRec;
                } else {
                    continue;
                }

                if (std::isnan(xval) || std::isnan(yval)) continue;

                v.overall->Fill(xval, yval);
                v.phiOverall[phiBinIdx]->Fill(xval, yval);
                v.binHists[phiBinIdx][thetaBinIdx]->Fill(xval, yval);
            }
        }
    }, colNames);

    gSystem->Exec(("mkdir -p " + outDir).c_str());

    for (auto &v : vars) {
        {
            TCanvas *c = new TCanvas(("c_" + v.saveName + "_overall").c_str(), "", 3000, 1500);
            v.overall->SetTitle((selecteddetector + " " + v.title + " (all #phi)").c_str());
            gPad->SetLogz();
            v.overall->Draw("COLZ");

            if (v.name == "deltaP:p") {
                v.overall->GetXaxis()->SetTitle("p [GeV]");
                v.overall->GetYaxis()->SetTitle("#delta p [GeV]");
            } else if (v.name == "p:deltaP") {
                v.overall->GetXaxis()->SetTitle("#delta p [GeV]");
                v.overall->GetYaxis()->SetTitle("p [GeV]");
            }

            std::string out = outDir + "/" + v.saveName + "_" + selecteddetector + "_overall.png";
            c->SaveAs(out.c_str());
            delete c;
            std::cout << "Saved: " << out << std::endl;
        }

        for (int phiBinIdx = 0; phiBinIdx < nPhiBins; ++phiBinIdx) {
            std::string phiTag   = GetPhiBinTag(phiBinIdx);
            std::string phiLabel = GetPhiBinLabel(phiBinIdx);

            std::string phiDir = outDir + "/" + phiTag;
            std::string fitParamOutDir = phiDir + "/ParamFits";
            gSystem->Exec(("mkdir -p " + phiDir).c_str());
            gSystem->Exec(("mkdir -p " + fitParamOutDir).c_str());

            std::vector<double> thetaMidVec, aVec, bVec, cVec;
            std::vector<double> aErrVec, bErrVec, cErrVec;

            {
                TCanvas *c = new TCanvas(("c_" + v.saveName + "_" + phiTag + "_overall").c_str(), "", 3000, 1500);
                v.phiOverall[phiBinIdx]->SetTitle((selecteddetector + " " + v.title + " in " + phiLabel).c_str());
                gPad->SetLogz();
                v.phiOverall[phiBinIdx]->Draw("COLZ");

                if (v.name == "deltaP:p") {
                    v.phiOverall[phiBinIdx]->GetXaxis()->SetTitle("p [GeV]");
                    v.phiOverall[phiBinIdx]->GetYaxis()->SetTitle("#delta p [GeV]");
                } else if (v.name == "p:deltaP") {
                    v.phiOverall[phiBinIdx]->GetXaxis()->SetTitle("#delta p [GeV]");
                    v.phiOverall[phiBinIdx]->GetYaxis()->SetTitle("p [GeV]");
                }

                std::string out = phiDir + "/" + v.saveName + "_" + selecteddetector + "_" + phiTag + "_overall.png";
                c->SaveAs(out.c_str());
                delete c;
                std::cout << "Saved: " << out << std::endl;
            }

            for (size_t thetaBinIdx = 0; thetaBinIdx <= thetaCuts.size(); ++thetaBinIdx) {
                TCanvas *c = new TCanvas(
                    ("c_" + v.saveName + "_" + phiTag + Form("_T%zu", thetaBinIdx)).c_str(),
                    "", 3000, 1500
                );

                TH2D* h2 = v.binHists[phiBinIdx][thetaBinIdx];
                h2->SetTitle((selecteddetector + " " + v.title + " in " + phiLabel + ", " +
                              GetThetaBinLabel(thetaBinIdx, thetaCuts)).c_str());
                gPad->SetLogz();
                h2->Draw("COLZ");

                TGraph* gPeak2 = MakePeakGraph(
                    h2,
                    phiDir + "/" + v.saveName + "_" + selecteddetector + "_" + phiTag + Form("_T%zu.txt", thetaBinIdx)
                );
                if (gPeak2 && gPeak2->GetN() > 0) gPeak2->Draw("PEZ SAME");

                TF1* deltaPFitModel = nullptr;
                std::string deltaPFitFormula;

                if (selectedPid == 2212) {
                    if (selecteddetector == "FD") {
                        deltaPFitFormula = "[0] + [1]/x + [2]/(x*x)";
                    } else if (selecteddetector == "CD" || selecteddetector == "ALL") {
                        deltaPFitFormula = "[0] + [1]*x + [2]*x*x";
                    } else {
                        std::cerr << "Unknown detector for proton: " << selecteddetector << "\n";
                        if (gPeak2) delete gPeak2;
                        delete c;
                        continue;
                    }
                } else if (selectedPid == 11 || selectedPid == 22) {
                    deltaPFitFormula = "[0] + [1]*x + [2]*x*x";
                } else {
                    std::cerr << "Unknown pid for fitting: " << selectedPid << "\n";
                    deltaPFitFormula = "[0] + [1]*x + [2]*x*x";
                }

                std::string fitName = "deltaPFitModel_" + phiTag + Form("_T%zu", thetaBinIdx);
                deltaPFitModel = new TF1(
                    fitName.c_str(),
                    deltaPFitFormula.c_str(),
                    h2->GetXaxis()->GetXmin(),
                    h2->GetXaxis()->GetXmax()
                );

                //deltaPFitModel->SetParameters(0.01, -0.01, 0.01);

                if (gPeak2 && gPeak2->GetN() > 2) {
                    gPeak2->Fit(deltaPFitModel, "Q0N");

                    std::ofstream foutFit(
                        phiDir + "/" + v.saveName + "_" + selecteddetector + "_" + phiTag +
                        Form("_T%zu_fit.txt", thetaBinIdx)
                    );
                    foutFit << "# Fit expression: " << deltaPFitFormula << "\n";
                    foutFit << "# Fit parameters:\n";
                    for (int ip = 0; ip < deltaPFitModel->GetNpar(); ++ip) {
                        foutFit << "p" << ip << " = " << deltaPFitModel->GetParameter(ip)
                                << " ± " << deltaPFitModel->GetParError(ip) << "\n";
                    }
                    foutFit.close();

                    deltaPFitModel->SetLineColor(kBlue + 2);
                    deltaPFitModel->SetLineWidth(4);
                    deltaPFitModel->Draw("SAME");

                    double thetaMid = 0.0;
                    if (thetaBinIdx == 0) {
                        thetaMid = thetaCuts[0];
                    } else if (thetaBinIdx == thetaCuts.size()) {
                        thetaMid = thetaCuts.back();
                    } else {
                        thetaMid = 0.5 * (thetaCuts[thetaBinIdx] + thetaCuts[thetaBinIdx - 1]);
                    }

                    if (deltaPFitModel->GetNpar() >= 3 && thetaBinIdx < thetaCuts.size()) {
                        thetaMidVec.push_back(thetaMid);
                        aVec.push_back(deltaPFitModel->GetParameter(0));
                        bVec.push_back(deltaPFitModel->GetParameter(1));
                        cVec.push_back(deltaPFitModel->GetParameter(2));
                        aErrVec.push_back(deltaPFitModel->GetParError(0));
                        bErrVec.push_back(deltaPFitModel->GetParError(1));
                        cErrVec.push_back(deltaPFitModel->GetParError(2));
                    }
                }

                if (v.name == "deltaP:p") {
                    h2->GetXaxis()->SetTitle("p [GeV]");
                    h2->GetYaxis()->SetTitle("#delta p [GeV]");
                } else if (v.name == "p:deltaP") {
                    h2->GetXaxis()->SetTitle("#delta p [GeV]");
                    h2->GetYaxis()->SetTitle("p [GeV]");
                }

                std::string out = phiDir + "/" + v.saveName + "_" + selecteddetector + "_" + phiTag +
                                  Form("_T%zu.png", thetaBinIdx);
                c->SaveAs(out.c_str());

                delete deltaPFitModel;
                if (gPeak2) delete gPeak2;
                delete c;

                std::cout << "Saved: " << out << std::endl;
            }

            auto PlotParamVsTheta = [&](const std::vector<double>& theta,
                                        const std::vector<double>& param,
                                        const std::vector<double>& errors,
                                        const std::string& paramName,
                                        const std::string& paramVsThetaFormula)
            {
                std::string phiPrefix = particleDetectorPrefix + "_" + phiTag;

                // 先确保 summary 结构按 phi 固定占位
                if (summaryPadPlotsByParam.find(paramName) == summaryPadPlotsByParam.end()) {
                    summaryPadPlotsByParam[paramName].resize(nPhiBins);
                    for (int iphi = 0; iphi < nPhiBins; ++iphi) {
                        summaryPadPlotsByParam[paramName][iphi].phiTag   = GetPhiBinTag(iphi);
                        summaryPadPlotsByParam[paramName][iphi].phiLabel = GetPhiBinLabel(iphi);
                        summaryPadPlotsByParam[paramName][iphi].fitFormula = paramVsThetaFormula;
                        summaryPadPlotsByParam[paramName][iphi].hasData = false;
                    }
                }

                // 无点也保留空位，不跳过
                if (theta.empty()) {
                    summaryPadPlotsByParam[paramName][phiBinIdx].fitFormula = paramVsThetaFormula;
                    summaryPadPlotsByParam[paramName][phiBinIdx].hasData = false;
                    return;
                }

                TCanvas *c = new TCanvas(
                    ("c_" + phiPrefix + "_" + paramName + "_vs_theta").c_str(),
                    "", 1800, 1200
                );

                TGraphErrors* g = new TGraphErrors(theta.size());
                for (size_t i = 0; i < theta.size(); ++i) {
                    g->SetPoint(i, theta[i], param[i]);
                    g->SetPointError(i, 0.0, errors[i]);
                }

                g->SetTitle((phiPrefix + ": " + paramName + " vs #theta").c_str());
                g->GetXaxis()->SetTitle("#theta [deg]");
                g->GetYaxis()->SetTitle((paramName + " value").c_str());
                g->GetYaxis()->SetTitleOffset(1.1);
                g->SetMarkerStyle(20);
                g->SetMarkerSize(1.2);
                g->SetMarkerColor(kBlack);
                g->SetLineColor(kBlack);
                g->Draw("AP");

                TF1* paramVsThetaFitModel = nullptr;
                if (theta.size() >= 3) {
                    paramVsThetaFitModel = new TF1(
                        ("paramVsThetaFit_" + phiPrefix + "_" + paramName).c_str(),
                        paramVsThetaFormula.c_str(),
                        *std::min_element(theta.begin(), theta.end()),
                        *std::max_element(theta.begin(), theta.end())
                    );

                    g->Fit(paramVsThetaFitModel, "Q");
                    paramVsThetaFitModel->SetLineColor(kRed + 1);
                    paramVsThetaFitModel->SetLineWidth(2);
                    paramVsThetaFitModel->Draw("SAME");

                    std::ofstream fout(fitParamOutDir + "/" + phiPrefix + "_" + paramName + "_vs_theta.txt");
                    fout << "# Fit expression: " << paramVsThetaFormula << "\n";
                    for (int ip = 0; ip < paramVsThetaFitModel->GetNpar(); ++ip) {
                        fout << "p" << ip << " = " << paramVsThetaFitModel->GetParameter(ip)
                             << " ± " << paramVsThetaFitModel->GetParError(ip) << "\n";
                    }
                    fout.close();
                }

                auto yRange = GetParamYRange(particleDetectorPrefix, paramName);
                g->GetYaxis()->SetRangeUser(yRange.first, yRange.second);

                std::string imgPath = fitParamOutDir + "/" + phiPrefix + "_" + paramName + "_vs_theta.png";
                c->SaveAs(imgPath.c_str());
                std::cout << "Saved: " << imgPath << std::endl;

                ParamThetaPlotInfo info;
                info.phiTag = phiTag;
                info.phiLabel = phiLabel;
                info.theta = theta;
                info.param = param;
                info.errors = errors;
                info.fitFormula = paramVsThetaFormula;
                info.hasData = true;
                summaryPadPlotsByParam[paramName][phiBinIdx] = info;

                if (paramVsThetaFitModel) delete paramVsThetaFitModel;
                delete g;
                delete c;
            };

            if (selecteddetector == "FD") {
                if (isOutBend) {
                    PlotParamVsTheta(thetaMidVec, aVec, aErrVec, "A_p", "[0] + [1]*x + [2]*x*x");
                    PlotParamVsTheta(thetaMidVec, bVec, bErrVec, "B_p", "[0] + [1]*x + [2]*x*x");
                    PlotParamVsTheta(thetaMidVec, cVec, cErrVec, "C_p", "[0] + [1]*x + [2]*x*x");
                } else {
                    PlotParamVsTheta(thetaMidVec, aVec, aErrVec, "A_p", "[0] + [1]*x + [2]*x*x");
                    PlotParamVsTheta(thetaMidVec, bVec, bErrVec, "B_p", "[0] + [1]*x + [2]*x*x");
                    PlotParamVsTheta(thetaMidVec, cVec, cErrVec, "C_p", "[0] + [1]*x + [2]*x*x");
                }
            } else if (selecteddetector == "CD" || selecteddetector == "ALL") {
                PlotParamVsTheta(thetaMidVec, aVec, aErrVec, "A_p", "[0]");
                PlotParamVsTheta(thetaMidVec, bVec, bErrVec, "B_p", "[0]");
                PlotParamVsTheta(thetaMidVec, cVec, cErrVec, "C_p", "[0]");
            }
        }
    }

    {
        std::string summaryDir = outDir + "/SummaryCanvas";
        gSystem->Exec(("mkdir -p " + summaryDir).c_str());

        for (const auto &kv : summaryPadPlotsByParam) {
            const std::string &paramName = kv.first;
            const auto &plotInfos = kv.second;

            if (plotInfos.empty()) continue;

            int nPads = nPhiBins;
            int nCols = std::ceil(std::sqrt((double)nPads));
            int nRows = std::ceil((double)nPads / nCols);

            TCanvas *cSummary = new TCanvas(
                ("c_summary_grid_" + particleDetectorPrefix + "_" + paramName).c_str(),
                "",
                700 * nCols,
                500 * nRows
            );
            cSummary->Divide(nCols, nRows, 0.001, 0.001);

            // 关键：保存到最后再 delete，避免 subplot 变空
            std::vector<TGraphErrors*> summaryGraphs;
            std::vector<TF1*> summaryFits;
            std::vector<TH1F*> summaryFrames;

            for (int iPad = 0; iPad < nPads; ++iPad) {
                cSummary->cd(iPad + 1);
                gPad->SetMargin(0.14, 0.05, 0.13, 0.10);
                gPad->SetTicks(1, 1);
                gPad->SetGridx();
                gPad->SetGridy();

                const auto &info = plotInfos[iPad];
                auto yRange = GetParamYRange(particleDetectorPrefix, paramName);

                double xMin = 0.0;
                double xMax = 0.0;
                if (!thetaCuts.empty()) {
                    xMin = thetaCuts.front() - 2.0;
                    xMax = thetaCuts.back() + 2.0;
                } else {
                    xMin = 0.0;
                    xMax = 50.0;
                }

                TH1F *frame = new TH1F(
                    Form("frame_%s_%s_%d", particleDetectorPrefix.c_str(), paramName.c_str(), iPad),
                    "",
                    100, xMin, xMax
                );
                summaryFrames.push_back(frame);

                frame->SetTitle((info.phiTag + ": " + paramName).c_str());
                frame->GetXaxis()->SetTitle("#theta [deg]");
                frame->GetYaxis()->SetTitle((paramName + " value").c_str());
                frame->GetXaxis()->SetTitleSize(0.055);
                frame->GetYaxis()->SetTitleSize(0.055);
                frame->GetXaxis()->SetLabelSize(0.045);
                frame->GetYaxis()->SetLabelSize(0.045);
                frame->GetYaxis()->SetTitleOffset(1.15);
                frame->SetMinimum(yRange.first);
                frame->SetMaximum(yRange.second);
                frame->Draw();

                TLatex lat;
                lat.SetNDC(true);
                lat.SetTextSize(0.05);
                lat.SetTextColor(kBlue + 2);
                lat.DrawLatex(0.17, 0.88, info.phiLabel.c_str());

                if (info.hasData && !info.theta.empty()) {
                    TGraphErrors *g = new TGraphErrors(info.theta.size());
                    summaryGraphs.push_back(g);

                    for (size_t i = 0; i < info.theta.size(); ++i) {
                        g->SetPoint(i, info.theta[i], info.param[i]);
                        g->SetPointError(i, 0.0, info.errors[i]);
                    }

                    g->SetMarkerStyle(20);
                    g->SetMarkerSize(0.9);
                    g->SetMarkerColor(kBlack);
                    g->SetLineColor(kBlack);
                    g->SetLineWidth(1);
                    g->Draw("P SAME");

                    if (info.theta.size() >= 3) {
                        TF1 *paramVsThetaFitModel = new TF1(
                            Form("fit_summary_%s_%s_%d", particleDetectorPrefix.c_str(), paramName.c_str(), iPad),
                            info.fitFormula.c_str(),
                            *std::min_element(info.theta.begin(), info.theta.end()),
                            *std::max_element(info.theta.begin(), info.theta.end())
                        );
                        summaryFits.push_back(paramVsThetaFitModel);

                        g->Fit(paramVsThetaFitModel, "Q");
                        paramVsThetaFitModel->SetLineColor(kRed + 1);
                        paramVsThetaFitModel->SetLineWidth(2);
                        paramVsThetaFitModel->Draw("SAME");
                    }
                } else {
                    TLatex lat2;
                    lat2.SetNDC(true);
                    lat2.SetTextSize(0.055);
                    lat2.SetTextColor(kGray + 2);
                    lat2.DrawLatex(0.30, 0.50, "No valid points");
                }
            }

            cSummary->Modified();
            cSummary->Update();

            std::string outPng = summaryDir + "/" + particleDetectorPrefix + "_" + paramName + "_summary_grid.png";
            std::string outPdf = summaryDir + "/" + particleDetectorPrefix + "_" + paramName + "_summary_grid.pdf";
            cSummary->SaveAs(outPng.c_str());
            cSummary->SaveAs(outPdf.c_str());

            std::cout << "Saved: " << outPng << std::endl;
            std::cout << "Saved: " << outPdf << std::endl;

            for (auto *f : summaryFits)   delete f;
            for (auto *g : summaryGraphs) delete g;
            for (auto *fr : summaryFrames) delete fr;
            delete cSummary;
        }
    }

    timer.Stop();
    std::cout << "Time for DrawDeltaPByThetaPhiBins: ";
    timer.Print();
}


void DrawDeltaPByThetaPhiBins(
    const int &selectedPid,
    const std::vector<float> &thetaCuts,
    const std::string selecteddetector,
    const int nPhiBins,
    double phiStartDeg,
    const std::vector<std::tuple<std::string,std::string,std::string,int,double,double,int,double,double>> &plotVars,
    const std::string &filename,
    const std::string &treename,
    const std::string& outDir = "ProtonMomCorrPhi",
    bool isOutBend = false,
    bool drawDeltaPFitCurve = true,
    const std::map<std::string, std::string>& paramVsThetaFormulaMap = {},
    const std::map<std::string, std::string>& paramVsPhiFormulaMap   = {})
{
    TStopwatch timer;
    timer.Start();

    ROOT::EnableImplicitMT();

    auto GetParticleName = [](int pid) -> std::string {
        if (pid == 11)   return "electron";
        if (pid == 22)   return "photon";
        if (pid == 2212) return "proton";
        return "pid" + std::to_string(pid);
    };

    auto WrapToRange = [](double x, double start) -> double {
        double y = std::fmod(x - start, 360.0);
        if (y < 0) y += 360.0;
        return y + start;
    };

    auto GetPhiBinIndex = [&](double phiDeg) -> int {
        double phiWrapped = WrapToRange(phiDeg, phiStartDeg);
        double dphi = 360.0 / nPhiBins;
        int idx = static_cast<int>((phiWrapped - phiStartDeg) / dphi);

        if (idx < 0) idx = 0;
        if (idx >= nPhiBins) idx = nPhiBins - 1;
        return idx;
    };

    auto GetPhiBinLabel = [&](int phiBinIdx) -> std::string {
        double dphi = 360.0 / nPhiBins;
        double phiMin = phiStartDeg + phiBinIdx * dphi;
        double phiMax = phiStartDeg + (phiBinIdx + 1) * dphi;
        char buf[128];
        sprintf(buf, "%.1f#circ #leq #phi < %.1f#circ", phiMin, phiMax);
        return std::string(buf);
    };

    auto GetPhiBinTag = [&](int phiBinIdx) -> std::string {
        double dphi = 360.0 / nPhiBins;
        double phiMin = phiStartDeg + phiBinIdx * dphi;
        double phiMax = phiStartDeg + (phiBinIdx + 1) * dphi;
        char buf[128];
        sprintf(buf, "P%02d_%.0fto%.0f", phiBinIdx, phiMin, phiMax);
        return std::string(buf);
    };

    auto GetPhiMid = [&](int phiBinIdx) -> double {
        double dphi = 360.0 / nPhiBins;
        double phiMin = phiStartDeg + phiBinIdx * dphi;
        double phiMax = phiStartDeg + (phiBinIdx + 1) * dphi;
        return 0.5 * (phiMin + phiMax);
    };

    auto GetThetaBinTag = [&](size_t thetaBinIdx) -> std::string {
        char buf[128];
        if (thetaCuts.empty()) {
            sprintf(buf, "T%02zu_all", thetaBinIdx);
            return std::string(buf);
        }

        if (thetaBinIdx == 0) {
            sprintf(buf, "T%02zu_lt%.1f", thetaBinIdx, thetaCuts[0]);
        } else if (thetaBinIdx == thetaCuts.size()) {
            sprintf(buf, "T%02zu_ge%.1f", thetaBinIdx, thetaCuts.back());
        } else {
            sprintf(buf, "T%02zu_%.1fto%.1f", thetaBinIdx, thetaCuts[thetaBinIdx - 1], thetaCuts[thetaBinIdx]);
        }
        return std::string(buf);
    };

    auto GetThetaMid = [&](size_t thetaBinIdx) -> double {
        if (thetaCuts.empty()) return 0.0;
        if (thetaBinIdx == 0) return thetaCuts[0];
        if (thetaBinIdx == thetaCuts.size()) return thetaCuts.back();
        return 0.5 * (thetaCuts[thetaBinIdx] + thetaCuts[thetaBinIdx - 1]);
    };

    auto GetParamYRange = [&](const std::string& prefix, const std::string& paramName) -> std::pair<double,double> {
        if (prefix == "proton_FD") {
            if (paramName == "A_p") return {-0.1, 0.1};
            if (paramName == "B_p") return {-0.1, 0.1};
            if (paramName == "C_p") return {-0.1, 0.1};
            return {-0.1, 0.1};
        }
        if (prefix == "proton_CD") {
            if (paramName == "A_p") return {-0.1, 0.1};
            if (paramName == "B_p") return {-0.1, 0.1};
            if (paramName == "C_p") return {-0.1, 0.1};
            return {-0.1, 0.1};
        }
        return {-0.1, 0.1};
    };

    auto GetDefaultParamVsThetaFormula = [&](const std::string& paramName) -> std::string {
        auto it = paramVsThetaFormulaMap.find(paramName);
        if (it != paramVsThetaFormulaMap.end()) return it->second;

        if (selecteddetector == "FD") {
            return "[0] + [1]*x + [2]*x*x";
        } else if (selecteddetector == "CD" || selecteddetector == "ALL") {
            return "[0]";
        }
        return "[0]";
    };

    auto GetDefaultParamVsPhiFormula = [&](const std::string& paramName) -> std::string {
        auto it = paramVsPhiFormulaMap.find(paramName);
        if (it != paramVsPhiFormulaMap.end()) return it->second;

        // 默认给一个可拟合的二次函数；你也可以改成周期函数，比如 [0]+[1]*sin((x-[2])*TMath::DegToRad())
        return "[0]";
    };

    auto HasPhysicalLowPMomentumTrend = [](TF1* fitModel, TGraph* peakGraph) -> bool {
        if (!fitModel || !peakGraph) return true;
        if (peakGraph->GetN() < 3 || fitModel->GetNpar() < 3) return true;

        double x0 = 0.0, y0 = 0.0;
        double x1 = 0.0, y1 = 0.0;
        peakGraph->GetPoint(0, x0, y0);
        peakGraph->GetPoint(1, x1, y1);
        if (x1 <= x0) return true;

        const double lowPSlope = fitModel->Derivative(x0);
        const double lowPEndDrop = fitModel->Eval(x0) - fitModel->Eval(x1);

        // FD proton low-p region should bend upward as p decreases:
        // when moving to larger p from the low edge, the curve should fall or stay flat.
        if (lowPSlope > 0.0 || lowPEndDrop < 0.0) return false;

        if (peakGraph->GetN() >= 3) {
            double x2 = 0.0, y2 = 0.0;
            peakGraph->GetPoint(2, x2, y2);
            if (x2 > x1 && fitModel->Eval(x1) < fitModel->Eval(x2)) return false;
        }

        return true;
    };

    auto WriteGraphPoints = [](const std::string& outFile,
                               const std::vector<double>& x,
                               const std::vector<double>& y,
                               const std::vector<double>& yerr,
                               const std::string& xTitle,
                               const std::string& yTitle)
    {
        std::ofstream fout(outFile);
        fout << "# " << xTitle << " " << yTitle << " " << yTitle << "_err\n";
        for (size_t i = 0; i < x.size(); ++i) {
            double ey = (i < yerr.size() ? yerr[i] : 0.0);
            fout << x[i] << " " << y[i] << " " << ey << "\n";
        }
        fout.close();
    };

    auto WriteFitCurve = [](const std::string& outFile,
                            TF1* f,
                            double xMin,
                            double xMax,
                            int nPoints,
                            const std::string& xTitle,
                            const std::string& yTitle)
    {
        if (!f) return;
        std::ofstream fout(outFile);
        fout << "# " << xTitle << " " << yTitle << "_fit\n";
        for (int i = 0; i < nPoints; ++i) {
            double x = xMin + (xMax - xMin) * i / std::max(1, nPoints - 1);
            fout << x << " " << f->Eval(x) << "\n";
        }
        fout.close();
    };

    auto WriteFitParameters = [](const std::string& outFile,
                                 const std::string& formula,
                                 TF1* f)
    {
        std::ofstream fout(outFile);
        fout << "# Fit expression: " << formula << "\n";
        if (f) {
            fout << "# Chi2 = " << f->GetChisquare() << "\n";
            fout << "# NDF  = " << f->GetNDF() << "\n";
            for (int ip = 0; ip < f->GetNpar(); ++ip) {
                fout << "p" << ip << " = " << f->GetParameter(ip)
                     << " +- " << f->GetParError(ip) << "\n";
            }
        } else {
            fout << "# No fit performed.\n";
        }
        fout.close();
    };

    std::string particleDetectorPrefix = GetParticleName(selectedPid) + "_" + selecteddetector;
    const size_t nThetaFitBins = thetaCuts.size(); // 只保存真正用于拟合的 theta bin，不含最后 overflow

    struct ParamThetaPlotInfo {
        std::string phiTag;
        std::string phiLabel;
        std::vector<double> theta;
        std::vector<double> param;
        std::vector<double> errors;
        std::string fitFormula;
        bool hasData = false;
    };

    struct ParamPhiPlotInfo {
        std::string thetaTag;
        std::string thetaLabel;
        std::vector<double> phi;
        std::vector<double> param;
        std::vector<double> errors;
        std::string fitFormula;
        bool hasData = false;
    };

    std::map<std::string, std::vector<ParamThetaPlotInfo>> summaryPadPlotsByParamVsTheta;
    std::map<std::string, std::vector<ParamPhiPlotInfo>>   summaryPadPlotsByParamVsPhi;

    // paramName -> thetaBinIdx -> info
    std::map<std::string, std::vector<ParamPhiPlotInfo>> paramVsPhiStorage;

    ROOT::RDataFrame df(treename, filename);
    gStyle->SetOptStat(0);

    struct Var3DInfo {
        std::string saveName;
        std::string title;
        std::string name;
        int nxbins;
        double xmin, xmax;
        int nybins;
        double ymin, ymax;

        TH2D* overall;
        std::vector<TH2D*> phiOverall;
        std::vector<std::vector<TH2D*>> binHists;
    };

    std::vector<Var3DInfo> vars;
    for (auto &cfg : plotVars) {
        Var3DInfo v;
        v.saveName = std::get<0>(cfg);
        v.title    = std::get<1>(cfg);
        v.name     = std::get<2>(cfg);
        v.nxbins   = std::get<3>(cfg);
        v.xmin     = std::get<4>(cfg);
        v.xmax     = std::get<5>(cfg);
        v.nybins   = std::get<6>(cfg);
        v.ymin     = std::get<7>(cfg);
        v.ymax     = std::get<8>(cfg);

        v.overall = new TH2D((v.saveName + "_overall").c_str(), "",
                             v.nxbins, v.xmin, v.xmax, v.nybins, v.ymin, v.ymax);
        v.overall->SetDirectory(nullptr);

        v.phiOverall.resize(nPhiBins, nullptr);
        for (int phiBinIdx = 0; phiBinIdx < nPhiBins; ++phiBinIdx) {
            std::string histName = v.saveName + Form("_P%02d_overall", phiBinIdx);
            v.phiOverall[phiBinIdx] = new TH2D(histName.c_str(), "",
                                               v.nxbins, v.xmin, v.xmax, v.nybins, v.ymin, v.ymax);
            v.phiOverall[phiBinIdx]->SetDirectory(nullptr);
        }

        v.binHists.resize(nPhiBins);
        for (int phiBinIdx = 0; phiBinIdx < nPhiBins; ++phiBinIdx) {
            v.binHists[phiBinIdx].resize(thetaCuts.size() + 1, nullptr);
            for (size_t thetaBinIdx = 0; thetaBinIdx <= thetaCuts.size(); ++thetaBinIdx) {
                std::string histName = v.saveName + Form("_P%02d_T%zu", phiBinIdx, thetaBinIdx);
                v.binHists[phiBinIdx][thetaBinIdx] = new TH2D(histName.c_str(), "",
                                                              v.nxbins, v.xmin, v.xmax, v.nybins, v.ymin, v.ymax);
                v.binHists[phiBinIdx][thetaBinIdx]->SetDirectory(nullptr);
            }
        }

        vars.push_back(v);
    }

    for (const std::string& paramName : {"A_p", "B_p", "C_p"}) {
        summaryPadPlotsByParamVsTheta[paramName].resize(nPhiBins);
        for (int iphi = 0; iphi < nPhiBins; ++iphi) {
            summaryPadPlotsByParamVsTheta[paramName][iphi].phiTag      = GetPhiBinTag(iphi);
            summaryPadPlotsByParamVsTheta[paramName][iphi].phiLabel    = GetPhiBinLabel(iphi);
            summaryPadPlotsByParamVsTheta[paramName][iphi].fitFormula  = GetDefaultParamVsThetaFormula(paramName);
            summaryPadPlotsByParamVsTheta[paramName][iphi].hasData     = false;
        }

        paramVsPhiStorage[paramName].resize(nThetaFitBins);
        summaryPadPlotsByParamVsPhi[paramName].resize(nThetaFitBins);
        for (size_t itheta = 0; itheta < nThetaFitBins; ++itheta) {
            paramVsPhiStorage[paramName][itheta].thetaTag     = GetThetaBinTag(itheta);
            paramVsPhiStorage[paramName][itheta].thetaLabel   = GetThetaBinLabel(itheta, thetaCuts);
            paramVsPhiStorage[paramName][itheta].fitFormula   = GetDefaultParamVsPhiFormula(paramName);
            paramVsPhiStorage[paramName][itheta].hasData      = false;

            summaryPadPlotsByParamVsPhi[paramName][itheta]    = paramVsPhiStorage[paramName][itheta];
        }
    }

    std::vector<std::string> colNames = {
        "REC_Particle_pid", "REC_Particle_theta", "REC_Particle_phi", "REC_Particle_p",
        "REC_Particle_status", "REC_Particle_pass",
        "MC_Particle_pid", "MC_Particle_px", "MC_Particle_py", "MC_Particle_pz"
    };

    df.Foreach([&](const ROOT::VecOps::RVec<int> &pid,
                   const ROOT::VecOps::RVec<float> &theta,
                   const ROOT::VecOps::RVec<float> &phi,
                   const ROOT::VecOps::RVec<float> &p,
                   const ROOT::VecOps::RVec<short> &status,
                   const ROOT::VecOps::RVec<int> &passPar,
                   const ROOT::VecOps::RVec<int> &mcpid,
                   const ROOT::VecOps::RVec<float> &mcpx,
                   const ROOT::VecOps::RVec<float> &mcpy,
                   const ROOT::VecOps::RVec<float> &mcpz)
    {
        for (size_t i = 0; i < pid.size(); ++i) {
            if (pid[i] != selectedPid) continue;
            if (!passPar[i]) continue;
            if (GetDetectorPart(status[i]) != selecteddetector && selecteddetector != "ALL") continue;

            double thetaDeg = theta[i] * 180.0 / M_PI;
            double phiDeg = phi[i] * 180.0 / M_PI;

            int thetaBinIdx = GetThetaRegionIndex(thetaDeg, thetaCuts);
            int phiBinIdx = GetPhiBinIndex(phiDeg);

            double pRec = p[i];
            double pMC  = NAN;

            for (size_t j = 0; j < mcpid.size(); ++j) {
                if (mcpid[j] == selectedPid) {
                    pMC = std::sqrt(mcpx[j]*mcpx[j] + mcpy[j]*mcpy[j] + mcpz[j]*mcpz[j]);
                    break;
                }
            }
            if (std::isnan(pMC)) continue;

            double deltaP = pMC - pRec;

            for (auto &v : vars) {
                double xval = NAN, yval = NAN;

                if (v.name == "deltaP:p") {
                    xval = pRec;
                    yval = deltaP;
                } else if (v.name == "p:deltaP") {
                    xval = deltaP;
                    yval = pRec;
                } else {
                    continue;
                }

                if (std::isnan(xval) || std::isnan(yval)) continue;

                v.overall->Fill(xval, yval);
                v.phiOverall[phiBinIdx]->Fill(xval, yval);
                v.binHists[phiBinIdx][thetaBinIdx]->Fill(xval, yval);
            }
        }
    }, colNames);

    gSystem->Exec(("mkdir -p " + outDir).c_str());

    for (auto &v : vars) {
        {
            TCanvas *c = new TCanvas(("c_" + v.saveName + "_overall").c_str(), "", 3000, 1500);
            v.overall->SetTitle((selecteddetector + " " + v.title + " (all #phi)").c_str());
            gPad->SetLogz();
            v.overall->Draw("COLZ");

            if (v.name == "deltaP:p") {
                v.overall->GetXaxis()->SetTitle("p [GeV]");
                v.overall->GetYaxis()->SetTitle("#delta p [GeV]");
            } else if (v.name == "p:deltaP") {
                v.overall->GetXaxis()->SetTitle("#delta p [GeV]");
                v.overall->GetYaxis()->SetTitle("p [GeV]");
            }

            std::string out = outDir + "/" + v.saveName + "_" + selecteddetector + "_overall.png";
            c->SaveAs(out.c_str());
            delete c;
            std::cout << "Saved: " << out << std::endl;
        }

        for (int phiBinIdx = 0; phiBinIdx < nPhiBins; ++phiBinIdx) {
            std::string phiTag   = GetPhiBinTag(phiBinIdx);
            std::string phiLabel = GetPhiBinLabel(phiBinIdx);

            std::string phiDir = outDir + "/" + phiTag;
            std::string fitParamOutDir = phiDir + "/ParamFits";
            gSystem->Exec(("mkdir -p " + phiDir).c_str());
            gSystem->Exec(("mkdir -p " + fitParamOutDir).c_str());

            std::vector<double> thetaMidVec;
            std::map<std::string, std::vector<double>> paramThetaValues;
            std::map<std::string, std::vector<double>> paramThetaErrors;

            {
                TCanvas *c = new TCanvas(("c_" + v.saveName + "_" + phiTag + "_overall").c_str(), "", 3000, 1500);
                v.phiOverall[phiBinIdx]->SetTitle((selecteddetector + " " + v.title + " in " + phiLabel).c_str());
                gPad->SetLogz();
                v.phiOverall[phiBinIdx]->Draw("COLZ");

                if (v.name == "deltaP:p") {
                    v.phiOverall[phiBinIdx]->GetXaxis()->SetTitle("p [GeV]");
                    v.phiOverall[phiBinIdx]->GetYaxis()->SetTitle("#delta p [GeV]");
                } else if (v.name == "p:deltaP") {
                    v.phiOverall[phiBinIdx]->GetXaxis()->SetTitle("#delta p [GeV]");
                    v.phiOverall[phiBinIdx]->GetYaxis()->SetTitle("p [GeV]");
                }

                std::string out = phiDir + "/" + v.saveName + "_" + selecteddetector + "_" + phiTag + "_overall.png";
                c->SaveAs(out.c_str());
                delete c;
                std::cout << "Saved: " << out << std::endl;
            }

            for (size_t thetaBinIdx = 0; thetaBinIdx <= thetaCuts.size(); ++thetaBinIdx) {
                TCanvas *c = new TCanvas(
                    ("c_" + v.saveName + "_" + phiTag + Form("_T%zu", thetaBinIdx)).c_str(),
                    "", 3000, 1500
                );

                TH2D* h2 = v.binHists[phiBinIdx][thetaBinIdx];
                h2->SetTitle((selecteddetector + " " + v.title + " in " + phiLabel + ", " +
                              GetThetaBinLabel(thetaBinIdx, thetaCuts)).c_str());
                gPad->SetLogz();
                h2->Draw("COLZ");

                TGraph* gPeak2 = MakePeakGraph(
                    h2,
                    phiDir + "/" + v.saveName + "_" + selecteddetector + "_" + phiTag + Form("_T%zu.txt", thetaBinIdx)
                );
                if (gPeak2 && gPeak2->GetN() > 0) gPeak2->Draw("PEZ SAME");

                TF1* deltaPFitModel = nullptr;
                std::string deltaPFitFormula;

                if (selectedPid == 2212) {
                    if (selecteddetector == "FD") {
                        deltaPFitFormula = "[0] + [1]/x + [2]/(x*x)";
                    } else if (selecteddetector == "CD" || selecteddetector == "ALL") {
                        deltaPFitFormula = "[0] + [1]*x + [2]*x*x";
                    } else {
                        std::cerr << "Unknown detector for proton: " << selecteddetector << "\n";
                        if (gPeak2) delete gPeak2;
                        delete c;
                        continue;
                    }
                } else if (selectedPid == 11 || selectedPid == 22) {
                    deltaPFitFormula = "[0] + [1]*x + [2]*x*x";
                } else {
                    std::cerr << "Unknown pid for fitting: " << selectedPid << "\n";
                    deltaPFitFormula = "[0] + [1]*x + [2]*x*x";
                }

                std::string fitName = "deltaPFitModel_" + phiTag + Form("_T%zu", thetaBinIdx);
                deltaPFitModel = new TF1(
                    fitName.c_str(),
                    deltaPFitFormula.c_str(),
                    h2->GetXaxis()->GetXmin(),
                    h2->GetXaxis()->GetXmax()
                );

                if (gPeak2 && gPeak2->GetN() > 2) {
                    gPeak2->Fit(deltaPFitModel, "Q0N");

                    if (selectedPid == 2212 && selecteddetector == "FD" &&
                        !HasPhysicalLowPMomentumTrend(deltaPFitModel, gPeak2)) {
                        TF1* constrainedFitModel = new TF1(
                            (fitName + "_constrained").c_str(),
                            deltaPFitFormula.c_str(),
                            h2->GetXaxis()->GetXmin(),
                            h2->GetXaxis()->GetXmax()
                        );

                        const double seedA = deltaPFitModel->GetParameter(0);
                        const double seedB = std::min(deltaPFitModel->GetParameter(1), -1e-6);
                        const double seedC = std::max(deltaPFitModel->GetParameter(2), 1e-6);
                        constrainedFitModel->SetParameters(seedA, seedB, seedC);

                        // Enforce the low-p upturn shape for FD protons.
                        constrainedFitModel->SetParLimits(1, -10.0, 0.0);
                        constrainedFitModel->SetParLimits(2, 0.0, 10.0);

                        const int constrainedStatus = gPeak2->Fit(constrainedFitModel, "Q0N");
                        if (constrainedStatus == 0) {
                            delete deltaPFitModel;
                            deltaPFitModel = constrainedFitModel;
                        } else {
                            delete constrainedFitModel;
                        }
                    }

                    std::string fitParFile = phiDir + "/" + v.saveName + "_" + selecteddetector + "_" + phiTag +
                                             Form("_T%zu_fit.txt", thetaBinIdx);
                    WriteFitParameters(fitParFile, deltaPFitFormula, deltaPFitModel);

                    if (drawDeltaPFitCurve) {
                        deltaPFitModel->SetLineColor(kBlue + 2);
                        deltaPFitModel->SetLineWidth(4);
                        deltaPFitModel->Draw("SAME");
                    }

                    double thetaMid = GetThetaMid(thetaBinIdx);

                    if (deltaPFitModel->GetNpar() >= 3 && thetaBinIdx < nThetaFitBins) {
                        double A = deltaPFitModel->GetParameter(0);
                        double B = deltaPFitModel->GetParameter(1);
                        double C = deltaPFitModel->GetParameter(2);

                        double AErr = deltaPFitModel->GetParError(0);
                        double BErr = deltaPFitModel->GetParError(1);
                        double CErr = deltaPFitModel->GetParError(2);

                        thetaMidVec.push_back(thetaMid);
                        paramThetaValues["A_p"].push_back(A);
                        paramThetaValues["B_p"].push_back(B);
                        paramThetaValues["C_p"].push_back(C);
                        paramThetaErrors["A_p"].push_back(AErr);
                        paramThetaErrors["B_p"].push_back(BErr);
                        paramThetaErrors["C_p"].push_back(CErr);

                        double phiMid = GetPhiMid(phiBinIdx);

                        paramVsPhiStorage["A_p"][thetaBinIdx].phi.push_back(phiMid);
                        paramVsPhiStorage["A_p"][thetaBinIdx].param.push_back(A);
                        paramVsPhiStorage["A_p"][thetaBinIdx].errors.push_back(AErr);
                        paramVsPhiStorage["A_p"][thetaBinIdx].hasData = true;

                        paramVsPhiStorage["B_p"][thetaBinIdx].phi.push_back(phiMid);
                        paramVsPhiStorage["B_p"][thetaBinIdx].param.push_back(B);
                        paramVsPhiStorage["B_p"][thetaBinIdx].errors.push_back(BErr);
                        paramVsPhiStorage["B_p"][thetaBinIdx].hasData = true;

                        paramVsPhiStorage["C_p"][thetaBinIdx].phi.push_back(phiMid);
                        paramVsPhiStorage["C_p"][thetaBinIdx].param.push_back(C);
                        paramVsPhiStorage["C_p"][thetaBinIdx].errors.push_back(CErr);
                        paramVsPhiStorage["C_p"][thetaBinIdx].hasData = true;

                    }
                }

                if (v.name == "deltaP:p") {
                    h2->GetXaxis()->SetTitle("p [GeV]");
                    h2->GetYaxis()->SetTitle("#delta p [GeV]");
                } else if (v.name == "p:deltaP") {
                    h2->GetXaxis()->SetTitle("#delta p [GeV]");
                    h2->GetYaxis()->SetTitle("p [GeV]");
                }

                std::string out = phiDir + "/" + v.saveName + "_" + selecteddetector + "_" + phiTag +
                                  Form("_T%zu.png", thetaBinIdx);
                c->SaveAs(out.c_str());

                delete deltaPFitModel;
                if (gPeak2) delete gPeak2;
                delete c;

                std::cout << "Saved: " << out << std::endl;
            }

            auto PlotParamVsTheta = [&](const std::vector<double>& theta,
                                        const std::vector<double>& param,
                                        const std::vector<double>& errors,
                                        const std::string& paramName,
                                        const std::string& paramVsThetaFormula)
            {
                std::string phiPrefix = particleDetectorPrefix + "_" + phiTag;

                if (theta.empty()) {
                    summaryPadPlotsByParamVsTheta[paramName][phiBinIdx].fitFormula = paramVsThetaFormula;
                    summaryPadPlotsByParamVsTheta[paramName][phiBinIdx].hasData = false;
                    return;
                }

                TCanvas *c = new TCanvas(
                    ("c_" + phiPrefix + "_" + paramName + "_vs_theta").c_str(),
                    "", 1800, 1200
                );

                TGraphErrors* g = new TGraphErrors(theta.size());
                for (size_t i = 0; i < theta.size(); ++i) {
                    g->SetPoint(i, theta[i], param[i]);
                    g->SetPointError(i, 0.0, errors[i]);
                }

                g->SetTitle((phiPrefix + ": " + paramName + " vs #theta").c_str());
                g->GetXaxis()->SetTitle("#theta [deg]");
                g->GetYaxis()->SetTitle((paramName + " value").c_str());
                g->GetYaxis()->SetTitleOffset(1.1);
                g->SetMarkerStyle(20);
                g->SetMarkerSize(1.2);
                g->SetMarkerColor(kBlack);
                g->SetLineColor(kBlack);
                g->Draw("AP");

                TF1* paramVsThetaFitModel = nullptr;
                if (theta.size() >= 3) {
                    paramVsThetaFitModel = new TF1(
                        ("paramVsThetaFit_" + phiPrefix + "_" + paramName).c_str(),
                        paramVsThetaFormula.c_str(),
                        *std::min_element(theta.begin(), theta.end()),
                        *std::max_element(theta.begin(), theta.end())
                    );

                    g->Fit(paramVsThetaFitModel, "Q");
                    paramVsThetaFitModel->SetLineColor(kRed + 1);
                    paramVsThetaFitModel->SetLineWidth(2);
                    paramVsThetaFitModel->Draw("SAME");
                }

                auto yRange = GetParamYRange(particleDetectorPrefix, paramName);
                g->GetYaxis()->SetRangeUser(yRange.first, yRange.second);

                std::string base = fitParamOutDir + "/" + phiPrefix + "_" + paramName + "_vs_theta";
                c->SaveAs((base + ".png").c_str());
                c->SaveAs((base + ".pdf").c_str());

                WriteFitParameters(base + "_fit.txt", paramVsThetaFormula, paramVsThetaFitModel);
                WriteGraphPoints(base + "_points.txt", theta, param, errors, "theta_deg", paramName);
                if (paramVsThetaFitModel) {
                    WriteFitCurve(base + "_curve.txt",
                                  paramVsThetaFitModel,
                                  *std::min_element(theta.begin(), theta.end()),
                                  *std::max_element(theta.begin(), theta.end()),
                                  300,
                                  "theta_deg",
                                  paramName);
                }

                std::cout << "Saved: " << base << ".png" << std::endl;

                ParamThetaPlotInfo info;
                info.phiTag = phiTag;
                info.phiLabel = phiLabel;
                info.theta = theta;
                info.param = param;
                info.errors = errors;
                info.fitFormula = paramVsThetaFormula;
                info.hasData = true;
                summaryPadPlotsByParamVsTheta[paramName][phiBinIdx] = info;

                if (paramVsThetaFitModel) delete paramVsThetaFitModel;
                delete g;
                delete c;
            };

            for (const std::string& paramName : {"A_p", "B_p", "C_p"}) {
                PlotParamVsTheta(thetaMidVec,
                                 paramThetaValues[paramName],
                                 paramThetaErrors[paramName],
                                 paramName,
                                 GetDefaultParamVsThetaFormula(paramName));
            }
        }
    }

    // =========================
    // Plot param vs phi
    // =========================
    {
        std::string phiSummaryBaseDir = outDir + "/ParamVsPhi";
        gSystem->Exec(("mkdir -p " + phiSummaryBaseDir).c_str());

        auto PlotParamVsPhi = [&](const std::vector<double>& phi,
                                  const std::vector<double>& param,
                                  const std::vector<double>& errors,
                                  const std::string& paramName,
                                  size_t thetaBinIdx,
                                  const std::string& paramVsPhiFormula)
        {
            std::string thetaTag   = GetThetaBinTag(thetaBinIdx);
            std::string thetaLabel = GetThetaBinLabel(thetaBinIdx, thetaCuts);
            std::string thetaDir   = phiSummaryBaseDir + "/" + thetaTag;
            gSystem->Exec(("mkdir -p " + thetaDir).c_str());

            if (phi.empty()) {
                summaryPadPlotsByParamVsPhi[paramName][thetaBinIdx].thetaTag    = thetaTag;
                summaryPadPlotsByParamVsPhi[paramName][thetaBinIdx].thetaLabel  = thetaLabel;
                summaryPadPlotsByParamVsPhi[paramName][thetaBinIdx].fitFormula  = paramVsPhiFormula;
                summaryPadPlotsByParamVsPhi[paramName][thetaBinIdx].hasData     = false;
                return;
            }

            TCanvas *c = new TCanvas(
                ("c_" + particleDetectorPrefix + "_" + thetaTag + "_" + paramName + "_vs_phi").c_str(),
                "",
                1800, 1200
            );

            TGraphErrors* g = new TGraphErrors(phi.size());
            for (size_t i = 0; i < phi.size(); ++i) {
                g->SetPoint(i, phi[i], param[i]);
                g->SetPointError(i, 0.0, errors[i]);
            }

            g->SetTitle((particleDetectorPrefix + " " + thetaTag + ": " + paramName + " vs #phi").c_str());
            g->GetXaxis()->SetTitle("#phi [deg]");
            g->GetYaxis()->SetTitle((paramName + " value").c_str());
            g->GetYaxis()->SetTitleOffset(1.1);
            g->SetMarkerStyle(20);
            g->SetMarkerSize(1.2);
            g->SetMarkerColor(kBlack);
            g->SetLineColor(kBlack);
            g->Draw("AP");

            TF1* fitModel = nullptr;
            if (phi.size() >= 3) {
                fitModel = new TF1(
                    ("fit_" + particleDetectorPrefix + "_" + thetaTag + "_" + paramName + "_vs_phi").c_str(),
                    paramVsPhiFormula.c_str(),
                    *std::min_element(phi.begin(), phi.end()),
                    *std::max_element(phi.begin(), phi.end())
                );
                g->Fit(fitModel, "Q");
                fitModel->SetLineColor(kRed + 1);
                fitModel->SetLineWidth(2);
                fitModel->Draw("SAME");
            }

            auto yRange = GetParamYRange(particleDetectorPrefix, paramName);
            g->GetYaxis()->SetRangeUser(yRange.first, yRange.second);

            TLatex lat;
            lat.SetNDC(true);
            lat.SetTextSize(0.04);
            lat.DrawLatex(0.15, 0.92, thetaLabel.c_str());

            std::string base = thetaDir + "/" + particleDetectorPrefix + "_" + thetaTag + "_" + paramName + "_vs_phi";
            c->SaveAs((base + ".png").c_str());
            c->SaveAs((base + ".pdf").c_str());

            WriteFitParameters(base + "_fit.txt", paramVsPhiFormula, fitModel);
            WriteGraphPoints(base + "_points.txt", phi, param, errors, "phi_deg", paramName);
            if (fitModel) {
                WriteFitCurve(base + "_curve.txt",
                              fitModel,
                              *std::min_element(phi.begin(), phi.end()),
                              *std::max_element(phi.begin(), phi.end()),
                              400,
                              "phi_deg",
                              paramName);
            }

            ParamPhiPlotInfo info;
            info.thetaTag = thetaTag;
            info.thetaLabel = thetaLabel;
            info.phi = phi;
            info.param = param;
            info.errors = errors;
            info.fitFormula = paramVsPhiFormula;
            info.hasData = true;
            summaryPadPlotsByParamVsPhi[paramName][thetaBinIdx] = info;

            if (fitModel) delete fitModel;
            delete g;
            delete c;
        };

        for (const std::string& paramName : {"A_p", "B_p", "C_p"}) {
            for (size_t thetaBinIdx = 0; thetaBinIdx < nThetaFitBins; ++thetaBinIdx) {
                const auto& info = paramVsPhiStorage[paramName][thetaBinIdx];
                PlotParamVsPhi(info.phi,
                               info.param,
                               info.errors,
                               paramName,
                               thetaBinIdx,
                               GetDefaultParamVsPhiFormula(paramName));
            }
        }
    }

    // =========================
    // summary canvas: param vs theta
    // =========================
    {
        std::string summaryDir = outDir + "/SummaryCanvas";
        gSystem->Exec(("mkdir -p " + summaryDir).c_str());

        for (const auto &kv : summaryPadPlotsByParamVsTheta) {
            const std::string &paramName = kv.first;
            const auto &plotInfos = kv.second;

            if (plotInfos.empty()) continue;

            int nPads = nPhiBins;
            int nCols = 6;
            int nRows = std::ceil((double)nPads / nCols);

            TCanvas *cSummary = new TCanvas(
                ("c_summary_grid_theta_" + particleDetectorPrefix + "_" + paramName).c_str(),
                "",
                700 * nCols,
                500 * nRows
            );
            cSummary->Divide(nCols, nRows, 0.001, 0.001);

            std::vector<TGraphErrors*> summaryGraphs;
            std::vector<TF1*> summaryFits;
            std::vector<TH1F*> summaryFrames;

            for (int iPad = 0; iPad < nPads; ++iPad) {
                cSummary->cd(iPad + 1);
                gPad->SetMargin(0.14, 0.05, 0.13, 0.10);
                gPad->SetTicks(1, 1);
                gPad->SetGridx();
                gPad->SetGridy();

                const auto &info = plotInfos[iPad];
                auto yRange = GetParamYRange(particleDetectorPrefix, paramName);

                double xMin = 0.0;
                double xMax = 50.0;
                if (!thetaCuts.empty()) {
                    xMin = thetaCuts.front() - 2.0;
                    xMax = thetaCuts.back() + 2.0;
                }

                TH1F *frame = new TH1F(
                    Form("frame_theta_%s_%s_%d", particleDetectorPrefix.c_str(), paramName.c_str(), iPad),
                    "",
                    100, xMin, xMax
                );
                summaryFrames.push_back(frame);

                frame->SetTitle((info.phiTag + ": " + paramName + " vs #theta").c_str());
                frame->GetXaxis()->SetTitle("#theta [deg]");
                frame->GetYaxis()->SetTitle((paramName + " value").c_str());
                frame->GetXaxis()->SetTitleSize(0.055);
                frame->GetYaxis()->SetTitleSize(0.055);
                frame->GetXaxis()->SetLabelSize(0.045);
                frame->GetYaxis()->SetLabelSize(0.045);
                frame->GetYaxis()->SetTitleOffset(1.15);
                frame->SetMinimum(yRange.first);
                frame->SetMaximum(yRange.second);
                frame->Draw();

                TLatex lat;
                lat.SetNDC(true);
                lat.SetTextSize(0.05);
                lat.SetTextColor(kBlue + 2);
                lat.DrawLatex(0.17, 0.88, info.phiLabel.c_str());

                if (info.hasData && !info.theta.empty()) {
                    TGraphErrors *g = new TGraphErrors(info.theta.size());
                    summaryGraphs.push_back(g);

                    for (size_t i = 0; i < info.theta.size(); ++i) {
                        g->SetPoint(i, info.theta[i], info.param[i]);
                        g->SetPointError(i, 0.0, info.errors[i]);
                    }

                    g->SetMarkerStyle(20);
                    g->SetMarkerSize(0.9);
                    g->SetMarkerColor(kBlack);
                    g->SetLineColor(kBlack);
                    g->SetLineWidth(1);
                    g->Draw("P SAME");

                    if (info.theta.size() >= 3) {
                        TF1 *fit = new TF1(
                            Form("fit_summary_theta_%s_%s_%d", particleDetectorPrefix.c_str(), paramName.c_str(), iPad),
                            info.fitFormula.c_str(),
                            *std::min_element(info.theta.begin(), info.theta.end()),
                            *std::max_element(info.theta.begin(), info.theta.end())
                        );
                        summaryFits.push_back(fit);

                        g->Fit(fit, "Q");
                        fit->SetLineColor(kRed + 1);
                        fit->SetLineWidth(2);
                        fit->Draw("SAME");
                    }
                } else {
                    TLatex lat2;
                    lat2.SetNDC(true);
                    lat2.SetTextSize(0.055);
                    lat2.SetTextColor(kGray + 2);
                    lat2.DrawLatex(0.30, 0.50, "No valid points");
                }
            }

            cSummary->Modified();
            cSummary->Update();

            std::string outPng = summaryDir + "/" + particleDetectorPrefix + "_" + paramName + "_vs_theta_summary_grid.png";
            std::string outPdf = summaryDir + "/" + particleDetectorPrefix + "_" + paramName + "_vs_theta_summary_grid.pdf";
            cSummary->SaveAs(outPng.c_str());
            cSummary->SaveAs(outPdf.c_str());

            std::cout << "Saved: " << outPng << std::endl;
            std::cout << "Saved: " << outPdf << std::endl;

            for (auto *f : summaryFits) delete f;
            for (auto *g : summaryGraphs) delete g;
            for (auto *fr : summaryFrames) delete fr;
            delete cSummary;
        }
    }

    // =========================
    // summary canvas: param vs phi
    // =========================
    {
        std::string summaryDir = outDir + "/SummaryCanvas";
        gSystem->Exec(("mkdir -p " + summaryDir).c_str());

        for (const auto &kv : summaryPadPlotsByParamVsPhi) {
            const std::string &paramName = kv.first;
            const auto &plotInfos = kv.second;

            if (plotInfos.empty()) continue;

            int nPads = (int)nThetaFitBins;
            if (nPads <= 0) continue;

            int nCols = std::ceil(std::sqrt((double)nPads));
            int nRows = std::ceil((double)nPads / nCols);

            TCanvas *cSummary = new TCanvas(
                ("c_summary_grid_phi_" + particleDetectorPrefix + "_" + paramName).c_str(),
                "",
                700 * nCols,
                500 * nRows
            );
            cSummary->Divide(nCols, nRows, 0.001, 0.001);

            std::vector<TGraphErrors*> summaryGraphs;
            std::vector<TF1*> summaryFits;
            std::vector<TH1F*> summaryFrames;

            double phiMinFrame = phiStartDeg;
            double phiMaxFrame = phiStartDeg + 360.0;

            for (int iPad = 0; iPad < nPads; ++iPad) {
                cSummary->cd(iPad + 1);
                gPad->SetMargin(0.14, 0.05, 0.13, 0.10);
                gPad->SetTicks(1, 1);
                gPad->SetGridx();
                gPad->SetGridy();

                const auto &info = plotInfos[iPad];
                auto yRange = GetParamYRange(particleDetectorPrefix, paramName);

                TH1F *frame = new TH1F(
                    Form("frame_phi_%s_%s_%d", particleDetectorPrefix.c_str(), paramName.c_str(), iPad),
                    "",
                    100, phiMinFrame, phiMaxFrame
                );
                summaryFrames.push_back(frame);

                frame->SetTitle((info.thetaTag + ": " + paramName + " vs #phi").c_str());
                frame->GetXaxis()->SetTitle("#phi [deg]");
                frame->GetYaxis()->SetTitle((paramName + " value").c_str());
                frame->GetXaxis()->SetTitleSize(0.055);
                frame->GetYaxis()->SetTitleSize(0.055);
                frame->GetXaxis()->SetLabelSize(0.045);
                frame->GetYaxis()->SetLabelSize(0.045);
                frame->GetYaxis()->SetTitleOffset(1.15);
                frame->SetMinimum(yRange.first);
                frame->SetMaximum(yRange.second);
                frame->Draw();

                TLatex lat;
                lat.SetNDC(true);
                lat.SetTextSize(0.05);
                lat.SetTextColor(kBlue + 2);
                lat.DrawLatex(0.17, 0.88, info.thetaLabel.c_str());

                if (info.hasData && !info.phi.empty()) {
                    TGraphErrors *g = new TGraphErrors(info.phi.size());
                    summaryGraphs.push_back(g);

                    for (size_t i = 0; i < info.phi.size(); ++i) {
                        g->SetPoint(i, info.phi[i], info.param[i]);
                        g->SetPointError(i, 0.0, info.errors[i]);
                    }

                    g->SetMarkerStyle(20);
                    g->SetMarkerSize(0.9);
                    g->SetMarkerColor(kBlack);
                    g->SetLineColor(kBlack);
                    g->SetLineWidth(1);
                    g->Draw("P SAME");

                    if (info.phi.size() >= 3) {
                        TF1 *fit = new TF1(
                            Form("fit_summary_phi_%s_%s_%d", particleDetectorPrefix.c_str(), paramName.c_str(), iPad),
                            info.fitFormula.c_str(),
                            *std::min_element(info.phi.begin(), info.phi.end()),
                            *std::max_element(info.phi.begin(), info.phi.end())
                        );
                        summaryFits.push_back(fit);

                        g->Fit(fit, "Q");
                        fit->SetLineColor(kRed + 1);
                        fit->SetLineWidth(2);
                        fit->Draw("SAME");
                    }
                } else {
                    TLatex lat2;
                    lat2.SetNDC(true);
                    lat2.SetTextSize(0.055);
                    lat2.SetTextColor(kGray + 2);
                    lat2.DrawLatex(0.30, 0.50, "No valid points");
                }
            }

            cSummary->Modified();
            cSummary->Update();

            std::string outPng = summaryDir + "/" + particleDetectorPrefix + "_" + paramName + "_vs_phi_summary_grid.png";
            std::string outPdf = summaryDir + "/" + particleDetectorPrefix + "_" + paramName + "_vs_phi_summary_grid.pdf";
            cSummary->SaveAs(outPng.c_str());
            cSummary->SaveAs(outPdf.c_str());

            std::cout << "Saved: " << outPng << std::endl;
            std::cout << "Saved: " << outPdf << std::endl;

            for (auto *f : summaryFits) delete f;
            for (auto *g : summaryGraphs) delete g;
            for (auto *fr : summaryFrames) delete fr;
            delete cSummary;
        }
    }

    timer.Stop();
    std::cout << "Time for DrawDeltaPByThetaPhiBins: ";
    timer.Print();
}

//================ example driver =================
void analysisMomentumCorrection() {
     ROOT::EnableImplicitMT(6); 
    //std::string path = "/work/clas12/yijie/clas12ana/analysis801/DISANA/build/nobkg/";
    std::string path = "../build/File7546/data/";
    //std::string path = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/sims/clasdis/outb/";
    //std::string path = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/sims/clasdis/spring2018/outb/";
    //std::string path = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/sims/clasdis/spring2018/inb/";
    //std::string path = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/fall2018/sims/DVCS/inb/class_dis/";
    std::string filename = path + "dfSelected_afterFid.root";
    std::string filenameCorrected = path + "dfSelected_afterFid_afterCorr.root";
    std::string treename = "dfSelected_afterFid";
    std::string treenameCorrected = "dfSelected_afterFid_afterCorr";
    bool isOutBend = true; // set to true for outbending data, false for inbending data
    const std::string& outDir ="ProtonMomCorr/";

    std::vector<float> thetaCutsFDelectron = {10,15,20,25};
    std::vector<float> thetaCutsFDproton = {10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,
                                            26,27,28,29,30,31,32,33,34,35};
    //std::vector<float> thetaCutsFDproton = {10,12,14,16,18,20,22,24,26,28,30,32,34};
    std::vector<float> thetaCutsFDphoton = {10,15,20,25};

    std::vector<float> thetaCutsCDproton = {33.0,36.1,39.2,42.3,45.3,48.4,51.5,54.6,57.7,60.7,63.8,66.9,70.0};

    std::vector<float> thetaCutsFTelectron = {3,3.5,4};
    std::vector<float> thetaCutsFTphoton = {3,3.5,4};

    std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> plotVarsFDelectron = {
        {"electron","electron","theta", 500, 0, 60},
        {"electron","electron","phi",   500, 0, 360},
        {"electron","electron","p",     500, 0, 8.0}
    };
    std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> plotVarsFDproton = {
        {"proton","proton","theta", 500, 0, 60},
        {"proton","proton","phi",   500, 0, 360},
        {"proton","proton","p",     500, 0, 3.0}
    };
    std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> plotVarsFDphoton = {
        {"photon","photon","theta", 500, 0, 60},
        {"photon","photon","phi",   500, 0, 360},
        {"photon","photon","p",     500, 0, 8.0}
    };
    std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> plotVarsCDproton = {
        {"proton","proton","theta", 500, 20, 150},
        {"proton","proton","phi",   500, 0, 360},
        {"proton","proton","p",     500, 0, 2.0}
    };
    std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> plotVarsFTelectron = {
        {"electron","electron","theta", 500, 0, 10},
        {"electron","electron","phi",   500, 0, 360},
        {"electron","electron","p",     500, 0, 8.0}
    };
    std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> plotVarsFTphoton = {
        {"photon","photon","theta", 500, 0, 10},
        {"photon","photon","phi",   500, 0, 360},
        {"photon","photon","p",     500, 0, 8.0}
    };
    std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> plotVarsAllelectron = {
        {"electron","electron","theta", 500, 0, 60},
        {"electron","electron","phi",   500, 0, 360},
        {"electron","electron","p",     500, 0, 8.0},
        {"electron","electron","beta",  500, 0, 2.0}
    };
    std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> plotVarsAllphoton = {
        {"photon","photon","theta", 500, 0, 60},
        {"photon","photon","phi",   500, 0, 360},
        {"photon","photon","p",     500, 0, 8.0},
        {"photon","photon","beta",  500, 0, 2.0}
    };
    std::vector<std::tuple<std::string,std::string,std::string,int,double,double>> plotVarsAllproton = {
        {"proton","proton","theta", 500, 0, 150},
        {"proton","proton","phi",   500, 0, 360},
        {"proton","proton","p",     500, 0, 3.0},
        {"proton","proton","beta",  500, 0, 2.0}
    };
/*
    DrawParticleKinematicsByThetaBins(22, {0}, "ALL", plotVarsAllphoton, filename, treename); // photon
    DrawParticleKinematicsByThetaBins(11, {0}, "ALL", plotVarsAllelectron, filenameCorrected, treenameCorrected); // electron
    DrawParticleKinematicsByThetaBins(2212, {0}, "ALL", plotVarsAllproton, filenameCorrected, treenameCorrected); // proton
*/
/*

    DrawParticleKinematicsByThetaBins(11, thetaCutsFDelectron, "FD", plotVarsFDelectron, filename, treename); // electron
    DrawParticleKinematicsByThetaBins(2212, thetaCutsFDproton, "FD", plotVarsFDproton, filename, treename); // proton
    DrawParticleKinematicsByThetaBins(22, thetaCutsFDphoton, "FD", plotVarsFDphoton, filename, treename); // photon

    DrawParticleKinematicsByThetaBins(2212, thetaCutsCDproton, "CD", plotVarsCDproton, filename, treename); // proton

    DrawParticleKinematicsByThetaBins(11, thetaCutsFTelectron, "FT", plotVarsFTelectron, filename, treename); // electron
    DrawParticleKinematicsByThetaBins(22, thetaCutsFTphoton, "FT", plotVarsFTphoton, filename, treename); // photon

    DrawParticleKinematicsByThetaBins(11, {0}, "ALL", plotVarsAllelectron, filename, treename); // electron
    DrawParticleKinematicsByThetaBins(22, {0}, "ALL", plotVarsAllphoton, filename, treename); // photon
    DrawParticleKinematicsByThetaBins(2212, {0}, "ALL", plotVarsAllproton, filename, treename); // proton

    DrawMCParticleKinematicsByThetaBins(11, {0}, plotVarsAllelectron, filename, treename); // electron
    DrawMCParticleKinematicsByThetaBins(2212, {0}, plotVarsAllproton, filename, treename); // proton
    DrawMCParticleKinematicsByThetaBins(22, {0}, plotVarsAllphoton, filename, treename); // photon
*/
 /* Draw2DParticleKinematicsByThetaBins(11,{0},"FD",{{"electron_p_vs_theta","electron p vs theta","p:theta",500,0,10,500,0,40}},filenameCorrected,treenameCorrected);
    Draw2DParticleKinematicsByThetaBins(11,{0},"FD",{{"electron_phi_vs_theta","electron phi vs theta","phi:theta",500,0,360,500,0,40}},filenameCorrected,treenameCorrected);
    Draw2DParticleKinematicsByThetaBins(11,{0},"FT",{{"electron_p_vs_theta","electron p vs theta","p:theta",500,0,10,500,0,10}},filenameCorrected,treenameCorrected);
    Draw2DParticleKinematicsByThetaBins(11,{0},"FT",{{"electron_phi_vs_theta","electron phi vs theta","phi:theta",500,0,360,500,0,10}},filenameCorrected,treenameCorrected);
    Draw2DParticleKinematicsByThetaBins(11,{0},"ALL",{{"electron_p_vs_theta","electron p vs theta","p:theta",500,0,10,500,0,40}},filenameCorrected,treenameCorrected);


    Draw2DParticleKinematicsByThetaBins(2212,{0},"FD",{{"proton_p_vs_theta","proton p vs theta","p:theta",500,0,3,500,0,60}},filenameCorrected,treenameCorrected);
    Draw2DParticleKinematicsByThetaBins(2212,{0},"FD",{{"proton_phi_vs_theta","proton phi vs theta","phi:theta",500,0,360,500,0,60}},filenameCorrected,treenameCorrected);
    Draw2DParticleKinematicsByThetaBins(2212,{0},"CD",{{"proton_p_vs_theta","proton p vs theta","p:theta",500,0,2,500,20,150}},filenameCorrected,treenameCorrected);
    Draw2DParticleKinematicsByThetaBins(2212,{0},"CD",{{"proton_phi_vs_theta","proton phi vs theta","phi:theta",500,0,360,500,20,150}},filenameCorrected,treenameCorrected);
    Draw2DParticleKinematicsByThetaBins(2212,{0},"ALL",{{"proton_p_vs_theta","proton p vs theta","p:theta",500,0,3,500,0,150}},filenameCorrected,treenameCorrected);
  

    Draw2DParticleKinematicsByThetaBins(22,{0},"FD",{{"photon_p_vs_theta","photon p vs theta","p:theta",500,0,8,500,0,40}},filenameCorrected,treenameCorrected);
    Draw2DParticleKinematicsByThetaBins(22,{0},"FD",{{"photon_phi_vs_theta","photon phi vs theta","phi:theta",500,0,360,500,0,40}},filenameCorrected,treenameCorrected);
    Draw2DParticleKinematicsByThetaBins(22,{0},"FT",{{"photon_p_vs_theta","photon p vs theta","p:theta",500,0,8,500,0,10}},filenameCorrected,treenameCorrected);
    Draw2DParticleKinematicsByThetaBins(22,{0},"FT",{{"photon_phi_vs_theta","photon phi vs theta","phi:theta",500,0,360,500,0,10}},filenameCorrected,treenameCorrected);
    Draw2DParticleKinematicsByThetaBins(22,{0},"ALL",{{"photon_p_vs_theta","photon p vs theta","p:theta",500,0,8,500,0,40}},filenameCorrected,treenameCorrected);
    Draw2DParticleKinematicsByThetaBins(22,{0},"ALL",{{"photon_phi_vs_theta","photon phi vs theta","phi:theta",500,0,360,500,2,8}},filenameCorrected,treenameCorrected);
*/
/*
    DrawDeltaPByThetaBins(11,thetaCutsFDelectron,"FD",{{"electron_deltaP_vs_p","electron #delta p vs p","deltaP:p",500,0,8,500,-0.03,0.03}},filename,treename);
    DrawDeltaPByThetaBins(11,thetaCutsFTelectron,"FT",{{"electron_deltaP_vs_p","electron #delta p vs p","deltaP:p",500,0,8,500,-1.0,1.0}},filename,treename);
    DrawDeltaPByThetaBins(11,{0},"ALL",{{"electron_deltaP_vs_p","electron #delta p vs p","deltaP:p",500,0,8,500,-0.1,0.1}},filename,treename);
  */

    //Draw2DParticleKinematicsByThetaBins(2212,{0},"FD",{{"proton_p_vs_theta","proton p vs theta","p:theta",500,0,5,500,0,60}},filename,treename);
    Draw2DParticleKinematicsByThetaBins(2212,{0},"FD",{{"proton_phi_vs_theta","proton phi vs theta","phi:theta",500,0,360,500,0,60}},filenameCorrected,treenameCorrected);
    //Draw2DParticleKinematicsByThetaBins(2212,{0},"CD",{{"proton_p_vs_theta","proton p vs theta","p:theta",500,0,5,500,20,150}},filename,treename);
    Draw2DParticleKinematicsByThetaBins(2212,{0},"CD",{{"proton_phi_vs_theta","proton phi vs theta","phi:theta",500,0,360,500,20,150}},filenameCorrected,treenameCorrected);
    //Draw2DParticleKinematicsByThetaBins(2212,{0},"FD",{{"proton_p_vs_theta","proton p vs theta","p:theta",200,0.4,2,200,10,50}},filenameCorrected,treenameCorrected);
 
   
    //DrawDeltaPByThetaBins(2212,thetaCutsFDproton,"FD",{{"proton_deltaP_vs_p","proton #delta p vs p","deltaP:p",100,0.2,3,100,-0.05,0.05}},filename,treename,outDir, isOutBend);
    //DrawDeltaPByThetaBins(2212,thetaCutsCDproton,"CD",{{"proton_deltaP_vs_p","proton #delta p vs p","deltaP:p",100,0.01,2.5,100,-0.2,0.2}},filename,treename,outDir, isOutBend);
    //DrawDeltaPByThetaBins(2212,{0},"ALL",{{"proton_deltaP_vs_p","proton #delta p vs p","deltaP:p",100,0,8,100,-0.1,0.1}},filename,treename,outDir, isOutBend);
    //DrawDeltaPByThetaBins(2212,thetaCutsFDproton,"FD",{{"proton_deltaP_vs_pcorr","corrected proton #delta p vs p","deltaP:p",100,0,2.5,100,-0.1,0.1}},filenameCorrected,treenameCorrected);
    //DrawDeltaPByThetaBins(2212,thetaCutsCDproton,"CD",{{"proton_deltaP_vs_pcorr","corrected proton #delta p vs p","deltaP:p",100,0,2.5,100,-0.2,0.2}},filenameCorrected,treenameCorrected);
/*
    DrawDeltaPByThetaPhiBins(2212,thetaCutsFDproton,"FD",1,-37,{{"proton_deltaP_vs_p","proton #delta p vs p","deltaP:p",100,0.2,5,100,-0.1,0.1}},filename,treename,"ProtonMomCorr",isOutBend,true);
    DrawDeltaPByThetaPhiBins(2212,thetaCutsFDproton,"FD",1,-37,{{"proton_deltaP_vs_p","proton #delta p vs p","deltaP:p",100,0.2,5,100,-0.1,0.1}},filenameCorrected,treenameCorrected,"ProtonMomCorrPhiCorrected",isOutBend,false);
    DrawDeltaPByThetaPhiBins(2212,thetaCutsFDproton,"FD",36,-37,{{"proton_deltaP_vs_p","proton #delta p vs p","deltaP:p",100,0.2,5,100,-0.1,0.1}},filename,treename,"ProtonMomCorrPhi36",isOutBend,true);
    DrawDeltaPByThetaPhiBins(2212,thetaCutsFDproton,"FD",36,-37,{{"proton_deltaP_vs_p","proton #delta p vs p","deltaP:p",100,0.2,5,100,-0.1,0.1}},filenameCorrected,treenameCorrected,"ProtonMomCorrPhiCorrected36",isOutBend,false);
    PlotMomentumCorrectionVsPhi_MultiP_FromParamVsPhiPoints(2212,"FD",{0.8, 1.0, 1.5, 2.0},"ProtonMomCorrPhi36","from_txt");
    PlotMomentumCorrectionVsPhi_MultiP_FromParamVsPhiPoints(2212,"FD",{0.8, 1.0, 1.5, 2.0},"ProtonMomCorrPhiCorrected36","from_txt");
*/
    //GeneratePiecewiseProtonCorrectionFunctionFromParamVsPhi(
    //    "ProtonMomCorrPhi36",
    //    "ProtonMomCorrPhi36/GeneratedPiecewiseProtonCorrection_RunDVCSAnalysis.inc",
    //    "AddProtonFDPiecewiseCorrection_FromProtonMomCorrPhi36"
    //);


    //DrawDeltaPByThetaPhiBins(2212,thetaCutsCDproton,"CD",6,-37,{{"proton_deltaP_vs_p","proton #delta p vs p","deltaP:p",100,0.01,2.5,100,-0.2,0.2}},filename,treename,"ProtonMomCorrPhi",isOutBend);
    //DrawDeltaPByThetaPhiBins(2212,thetaCutsFDproton,"FD",18,{{"proton_deltaP_vs_p","proton #delta p vs p","deltaP:p",100,0.2,3,100,-0.05,0.05}},filenameCorrected,treenameCorrected,"ProtonMomCorrPhiCorr",isOutBend);
    //DrawDeltaPByThetaPhiBins(2212,thetaCutsCDproton,"CD",6,{{"proton_deltaP_vs_p","proton #delta p vs p","deltaP:p",100,0.01,2.5,100,-0.2,0.2}},filenameCorrected,treenameCorrected,"ProtonMomCorrPhiCorr",isOutBend);
    
    //DrawDeltaPByThetaBins(22,thetaCutsFDphoton,"FD",{{"photon_deltaP_vs_p","photon #delta p vs p","deltaP:p",500,0,8,500,-0.5,0.5}},filename,treename);
    //DrawDeltaPByThetaBins(22,thetaCutsFTphoton,"FT",{{"photon_deltaP_vs_p","photon #delta p vs p","deltaP:p",500,0,8,500,-0.2,0.2}},filename,treename);
    //DrawDeltaPByThetaBins(22,{0},"ALL",{{"photon_deltaP_vs_p","photon #delta p vs p","deltaP:p",500,0,8,500,-0.5,0.5}},filename,treename);

   

//    std::vector<double> pgridP = {0.4, 0.75, 1.10, 1.75, 2.75};
//PlotMomentumCorrection_AllDetectors_FromFiles(
  //  /*selectedPid=*/2212,
  //  /*pValues=*/pgridP,
  //  /*thetaMinFD, thetaMaxFD=*/ 5.0, 40.0,
  //  /*thetaMinCD, thetaMaxCD=*/ 20.0, 70.0,
  //  /*thetaMinALL,thetaMaxALL=*/5.0, 70.0,
  //  /*nThetaPoints=*/300,
  //  /*baseDir=*/outDir+"ParamFits"
//);

    gApplication->Terminate(0);
}
