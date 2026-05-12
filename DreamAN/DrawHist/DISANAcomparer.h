#ifndef DISANA_COMPARER_H
#define DISANA_COMPARER_H

// ROOT headers
#include <TCanvas.h>
#include <TLegend.h>
#include <sys/stat.h>
// STL headers
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>

// Project-specific headers
#include <chrono>

#include "DISANAplotter.h"
#include "DrawStyle.h"

namespace fs = std::filesystem;
// color palettes for different models
std::vector<std::tuple<double, double, double>> modelShades = {
    {0.20, 0.30, 0.85},  // Blue
    {0.90, 0.45, 0.10},  // Orange
    {0.00, 0.60, 0.60},  // Teal green
    {0.00, 0.70, 0.00},  // Green
    {0.60, 0.30, 0.80},  // Purple
    {0.85, 0.10, 0.25},  // Red
    {0.40, 0.40, 0.40}   // Gray (fallback)
};

class DISANAcomparer {
 public:
  // Set the bin ranges used for cross-section calculations and plotting
  void SetXBinsRanges(BinManager bins) { fXbins = bins; }

  void NormalizeHistogram(TH1* hist) {
    if (!hist) return;
    double integral = hist->Integral();
    if (integral > 0) hist->Scale(1.0 / integral);
  }
  // Add a new model with its DataFrame, label, and beam energy
  void AddModelwithPi0Corr(ROOT::RDF::RNode df_dvcs_data, ROOT::RDF::RNode df_pi0_data, ROOT::RDF::RNode df_dvcs_pi0mc, ROOT::RDF::RNode df_pi0_pi0mc,
                           ROOT::RDF::RNode df_gen_dvcsmc, ROOT::RDF::RNode df_accept_dvcsmc, ROOT::RDF::RNode df_dvcsmc_bkg, ROOT::RDF::RNode df_dvcsmc_nobkg,
                           ROOT::RDF::RNode df_dvcsmc_rad, ROOT::RDF::RNode df_dvcsmc_norad, ROOT::RDF::RNode df_dvcsmc_p1cut, const std::string& label, double beamEnergy,
                           bool fPi0Correction = false, bool fAcceptanceCorrection = false, bool fEfficiencyCorrection = false, bool fRadiativeCorrection = false,
                           bool fP1cut = false, double luminosity = 1.0, double I_avg = 60.0, double I_mc = 60.0, double eff_corr = 1.0) {
    auto plotter = std::make_unique<DISANAplotter>(DVCSModeTag{}, df_dvcs_data, beamEnergy, luminosity, I_avg, I_mc, eff_corr, df_pi0_data, df_dvcs_pi0mc, df_pi0_pi0mc, df_gen_dvcsmc, df_accept_dvcsmc,
                                                   df_dvcsmc_bkg, df_dvcsmc_nobkg, df_dvcsmc_rad, df_dvcsmc_norad, df_dvcsmc_p1cut);
    std::cout << "Adding model: " << label << " with beam energy: " << beamEnergy << " GeV with Pi0 Correction: " << fPi0Correction
              << ", Acceptance Correction: " << fAcceptanceCorrection << ", Background Merging efficiency: " << fEfficiencyCorrection
              << ", Radiative Correction: " << fRadiativeCorrection << ", P1 cut: " << fP1cut << std::endl;
    plotter->SetPlotApplyCorrection(fPi0Correction);
    plotter->SetPlotApplyAcceptanceCorrection(fAcceptanceCorrection);
    plotter->SetPlotApplyEfficiencyCorrection(fEfficiencyCorrection);
    plotter->SetPlotApplyRadiativeCorrection(fRadiativeCorrection);
    plotter->SetPlotApplyP1Cut(fP1cut);
    plotter->GenerateKinematicHistos("el");
    plotter->GenerateKinematicHistos("pro");
    plotter->GenerateKinematicHistos("pho");
    plotter->GeneratePi0KinematicHistos("el");
    plotter->GeneratePi0KinematicHistos("pro");
    plotter->GeneratePi0KinematicHistos("pho");
    plotter->GeneratePi0KinematicHistos("pho2");
    labels.push_back(label);
    plotters.push_back(std::move(plotter));
  }

  void AddModel(ROOT::RDF::RNode df, const std::string& label, double beamEnergy, double luminosity = 1.0) {
    auto plotter = std::make_unique<DISANAplotter>(DVCSModeTag{}, df, beamEnergy, luminosity);
    std::cout << "Adding model: " << label << " with beam energy: " << beamEnergy << " GeV without Pi0 Correction" << std::endl;
    plotter->GenerateKinematicHistos("el");
    plotter->GenerateKinematicHistos("pro");
    plotter->GenerateKinematicHistos("pho");
    labels.push_back(label);
    plotters.push_back(std::move(plotter));
  }

  void AddModelPhi(ROOT::RDF::RNode df_data, const std::string& label, double beamEnergy, double luminosity) {
    auto plotter = std::make_unique<DISANAplotter>(PhiModeTag{}, df_data, beamEnergy, luminosity);
    labels.push_back(label);
    plotter->GeneratePhiKinematicHistos("el");
    plotter->GeneratePhiKinematicHistos("pro");
    plotter->GeneratePhiKinematicHistos("kMinus");
    plotter->GeneratePhiKinematicHistos("kPlus");
    plotters.push_back(std::move(plotter));
  }

  void AddModelPhi(ROOT::RDF::RNode df_data, const std::string& label, double beamEnergy, double luminosity, ROOT::RDF::RNode df_gen, ROOT::RDF::RNode df_rec,
                   ROOT::RDF::RNode df_radRatio, bool fAcc = false, bool fEff = false, bool fRad = false) {
    auto plotter = std::make_unique<DISANAplotter>(PhiModeTag{}, df_data, beamEnergy, luminosity, df_gen, df_rec, std::nullopt, std::nullopt, df_radRatio);
    plotter->GeneratePhiKinematicHistos("el");
    plotter->GeneratePhiKinematicHistos("pro");
    plotter->GeneratePhiKinematicHistos("kMinus");
    plotter->GeneratePhiKinematicHistos("kPlus");
    plotter->SetPlotApplyAcceptanceCorrection(fAcc);
    plotter->SetPlotApplyEfficiencyCorrection(fEff);
    plotter->SetPlotApplyRadiativeCorrection(fRad);

    labels.push_back(label);
    plotters.push_back(std::move(plotter));
  }

  void AddModelPhi(ROOT::RDF::RNode df_data, const std::string& label, double beamEnergy, double luminosity, ROOT::RDF::RNode df_gen, ROOT::RDF::RNode df_rec,
                   ROOT::RDF::RNode df_bkg, ROOT::RDF::RNode df_nobkg, ROOT::RDF::RNode df_radRatio, bool fAcc = false, bool fEff = false, bool fRad = false) {
    auto plotter = std::make_unique<DISANAplotter>(PhiModeTag{}, df_data, beamEnergy, luminosity, df_gen, df_rec, df_bkg, df_nobkg, df_radRatio);
    plotter->GeneratePhiKinematicHistos("el");
    plotter->GeneratePhiKinematicHistos("pro");
    plotter->GeneratePhiKinematicHistos("kMinus");
    plotter->GeneratePhiKinematicHistos("kPlus");
    plotter->SetPlotApplyAcceptanceCorrection(fAcc);
    plotter->SetPlotApplyEfficiencyCorrection(fEff);
    plotter->SetPlotApplyRadiativeCorrection(fRad);

    labels.push_back(label);
    plotters.push_back(std::move(plotter));
  }

  // Set the output directory for saving plots
  void SetOutputDir(const std::string& outdir) {
    outputDir = outdir;
    if (!fs::exists(outputDir)) {
      fs::create_directories(outputDir);
    }
  }

  // Enable or disable individual variable plotting
  void PlotIndividual(bool plotInd) { plotIndividual = plotInd; }

  // Set plot styles for various plot types
  void SetKinStyle(const DrawStyle& style) { styleKin_ = style; }
  void SetDVCSStyle(const DrawStyle& style) { styleDVCS_ = style; }
  void SetCrossSectionStyle(const DrawStyle& style) { styleCrossSection_ = style; }
  void SetBSAStyle(const DrawStyle& style) { styleBSA_ = style; }
  void UseFittedPhiYields(bool on = true) { useFittedYields_ = on; }

  // Enable or disable correctio
  void SetApplyCorrection(bool apply) { applyCorrection = apply; }

  // Load correction histogram from ROOT file
  void LoadCorrectionHistogram(const std::string& filename, const std::string& histoname = "h_correction") {
    correctionHist = nullptr;  // Reset after applying to avoid reusing the same histogram
    TFile* f = TFile::Open(filename.c_str(), "READ");

    if (!f || f->IsZombie()) {
      std::cerr << "Error: Cannot open correction file: " << filename << "\n";
      return;
    }

    correctionHist = dynamic_cast<THnSparseD*>(f->Get(histoname.c_str()));
    if (!correctionHist) {
      std::cerr << "Error: Correction histogram '" << histoname << "' not found in file: " << filename << "\n";
      return;
    }

    // correctionHist->SetDirectory(0);  // Detach from file
    f->Close();
    delete f;
    std::cout << "✅ Correction histogram loaded: " << histoname << "\n";
  }
  /// get mean values of Q^2 and x_B
  std::vector<std::vector<std::vector<std::tuple<double, double, double>>>> getMeanQ2xBt_old(const BinManager& bins, std::unique_ptr<DISANAplotter>& plotter) {
    const auto& xb_bins = bins.GetXBBins();
    const auto& q2_bins = bins.GetQ2Bins();
    const auto& t_bins = bins.GetTBins();

    size_t n_xb = xb_bins.size() - 1;
    size_t n_q2 = q2_bins.size() - 1;
    size_t n_t = t_bins.size() - 1;

    auto rdf = plotter->GetRDF();

    std::vector<std::vector<std::vector<std::tuple<double, double, double>>>> result(
        n_xb, std::vector<std::vector<std::tuple<double, double, double>>>(n_q2, std::vector<std::tuple<double, double, double>>(n_t)));

    for (size_t ix = 0; ix < n_xb; ++ix) {
      for (size_t iq = 0; iq < n_q2; ++iq) {
        for (size_t it = 0; it < n_t; ++it) {
          double xb_lo = xb_bins[ix], xb_hi = xb_bins[ix + 1];
          double q2_lo = q2_bins[iq], q2_hi = q2_bins[iq + 1];
          double t_lo = t_bins[it], t_hi = t_bins[it + 1];

          // Apply filter
          auto rdf_cut = rdf.Filter(Form("xB >= %f && xB < %f", xb_lo, xb_hi)).Filter(Form("Q2 >= %f && Q2 < %f", q2_lo, q2_hi)).Filter(Form("t >= %f && t < %f", t_lo, t_hi));

          // Compute means
          double mean_xB = rdf_cut.Mean("xB").GetValue();
          double mean_Q2 = rdf_cut.Mean("Q2").GetValue();
          double mean_t = rdf_cut.Mean("t").GetValue();

          result[ix][iq][it] = std::make_tuple(mean_xB, mean_Q2, mean_t);
        }
      }
    }

    return result;
  }

  std::vector<std::vector<std::vector<std::tuple<double, double, double>>>>
  getMeanQ2xBt(const BinManager& bins, std::unique_ptr<DISANAplotter>& plotter) {
    const auto& xb_bins = bins.GetXBBins();
    const auto& q2_bins = bins.GetQ2Bins();
    const auto& t_bins  = bins.GetTBins();
  
    const size_t n_xb = xb_bins.size() - 1;
    const size_t n_q2 = q2_bins.size() - 1;
    const size_t n_t  = t_bins.size() - 1;
  
    auto rdf = plotter->GetRDF();
  
    // 保持原来的返回结构
    std::vector<std::vector<std::vector<std::tuple<double, double, double>>>> result(
        n_xb,
        std::vector<std::vector<std::tuple<double, double, double>>>(
            n_q2,
            std::vector<std::tuple<double, double, double>>(n_t, std::make_tuple(0.0, 0.0, 0.0))
        )
    );
  
    struct Accumulator {
      double sum_xB = 0.0;
      double sum_Q2 = 0.0;
      double sum_t  = 0.0;
      unsigned long long count = 0;
    };
  
    auto findBin = [](double x, const std::vector<double>& edges) -> int {
      // 与你原来的条件一致：lo <= x < hi
      auto it = std::upper_bound(edges.begin(), edges.end(), x);
      if (it == edges.begin() || it == edges.end()) return -1;
      return static_cast<int>(std::distance(edges.begin(), it) - 1);
    };
  
    const size_t nCells = n_xb * n_q2 * n_t;
    const unsigned int nSlots = rdf.GetNSlots();
  
    // 每个 slot 一份局部缓存，避免多线程竞争
    std::vector<std::vector<Accumulator>> slotAccs(
        nSlots, std::vector<Accumulator>(nCells)
    );
  
    auto rdf_binned = rdf
      .Define("ix_bin", [xb_bins, findBin](double xB) {
        return findBin(xB, xb_bins);
      }, {"xB"})
      .Define("iq_bin", [q2_bins, findBin](double Q2) {
        return findBin(Q2, q2_bins);
      }, {"Q2"})
      .Define("it_bin", [t_bins, findBin](double t) {
        return findBin(t, t_bins);
      }, {"t"})
      .Filter([](int ix, int iq, int it) {
        return ix >= 0 && iq >= 0 && it >= 0;
      }, {"ix_bin", "iq_bin", "it_bin"});
  
    rdf_binned.ForeachSlot(
        [&slotAccs, n_q2, n_t](unsigned int slot, double xB, double Q2, double t,
                               int ix, int iq, int it) {
          const size_t idx = static_cast<size_t>(ix) * n_q2 * n_t
                           + static_cast<size_t>(iq) * n_t
                           + static_cast<size_t>(it);
  
          auto& acc = slotAccs[slot][idx];
          acc.sum_xB += xB;
          acc.sum_Q2 += Q2;
          acc.sum_t  += t;
          ++acc.count;
        },
        {"xB", "Q2", "t", "ix_bin", "iq_bin", "it_bin"}
    );
  
    // 合并各个 slot 的结果
    std::vector<Accumulator> totalAccs(nCells);
    for (unsigned int slot = 0; slot < nSlots; ++slot) {
      for (size_t idx = 0; idx < nCells; ++idx) {
        totalAccs[idx].sum_xB += slotAccs[slot][idx].sum_xB;
        totalAccs[idx].sum_Q2 += slotAccs[slot][idx].sum_Q2;
        totalAccs[idx].sum_t  += slotAccs[slot][idx].sum_t;
        totalAccs[idx].count  += slotAccs[slot][idx].count;
      }
    }
  
    // 写回原来的 3D tuple 结构
    for (size_t ix = 0; ix < n_xb; ++ix) {
      for (size_t iq = 0; iq < n_q2; ++iq) {
        for (size_t it = 0; it < n_t; ++it) {
          const size_t idx = ix * n_q2 * n_t + iq * n_t + it;
          const auto& acc = totalAccs[idx];
  
          if (acc.count > 0) {
            result[ix][iq][it] = std::make_tuple(
                acc.sum_xB / acc.count,
                acc.sum_Q2 / acc.count,
                acc.sum_t  / acc.count
            );
          } else {
            // 空 bin 时保持默认值
            result[ix][iq][it] = std::make_tuple(0.0, 0.0, 0.0);
            // 如果你更想看 NaN，可以改成：
            // const double nan = std::numeric_limits<double>::quiet_NaN();
            // result[ix][iq][it] = std::make_tuple(nan, nan, nan);
          }
        }
      }
    }
  
    return result;
  }

  // Plot all basic kinematic distributions (p, theta, phi) for all particle types
  void PlotKinematicComparison_phiAna() {
    TCanvas* canvas = new TCanvas("KinematicComparison", "Kinematic Comparison", 1800, 1200);
    canvas->Divide(4, 4);

    std::vector<std::string> types = {"el", "pro", "kPlus", "kMinus"};
    std::vector<std::string> vars = {"p", "theta", "phi", "vz"};

    int pad = 1;
    for (const auto& type : types) {
      for (const auto& var : vars) {
        PlotVariableComparison(type, var, pad++, canvas);
      }
    }

    canvas->Update();
    canvas->SaveAs((outputDir + "KinematicComparison_phiAna.pdf").c_str());

    // Optionally save individual plots
    if (plotIndividual) {
      for (const auto& type : types) {
        for (const auto& var : vars) {
          PlotSingleVariableComparison(type, var);
        }
      }
    }

    std::cout << "Saved kinematic comparison plots to: " << outputDir + "/KinematicComparison_phiAna.pdf" << std::endl;
    delete canvas;
  }

  // Plot all basic kinematic distributions (p, theta, phi) for all particle types
  void PlotKinematicComparison() {
    TCanvas* canvas = new TCanvas("KinematicComparison", "Kinematic Comparison", 1800, 1200);
    canvas->Divide(3, 3);

    std::vector<std::string> types = {"el", "pro", "pho"};
    std::vector<std::string> vars = {"p", "theta", "phi"};

    int pad = 1;
    for (const auto& type : types) {
      for (const auto& var : vars) {
        PlotVariableComparison(type, var, pad++, canvas);
      }
    }

    canvas->Update();
    canvas->SaveAs((outputDir + "KinematicComparison.pdf").c_str());

    // Optionally save individual plots
    if (plotIndividual) {
      for (const auto& type : types) {
        for (const auto& var : vars) {
          PlotSingleVariableComparison(type, var);
        }
      }
    }

    std::cout << "Saved kinematic comparison plots to: " << outputDir + "/KinematicComparison.pdf" << std::endl;
    delete canvas;
  }

  // Plot one specific variable (e.g., p) for a given particle type (e.g., "el")
  void PlotVariableComparison(const std::string& type, const std::string& var, int padIndex, TCanvas* canvas) {
    canvas->cd(padIndex);
    std::string hname_target = "rec" + type + "_" + var;

    TLegend* legend = new TLegend(0.6, 0.7, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    bool first = true;
    styleKin_.StylePad((TPad*)gPad);

    for (size_t i = 0; i < plotters.size(); ++i) {
      const auto& histograms = plotters[i]->GetAllHistograms();
      TH1* target = nullptr;

      for (TH1* h : histograms) {
        if (std::string(h->GetName()) == hname_target) {
          target = h;
          break;
        }
      }

      if (!target) {
        std::cerr << "[PlotVariableComparison]: Histogram " << hname_target << " not found for model [" << labels[i] << "]\n";
        continue;
      }
      NormalizeHistogram(target);
      styleKin_.StyleTH1(target);
      // styleKin_.StyleTH1(target, typeToParticle[type].c_str());
      auto [cr, cg, cb] = modelShades[i % modelShades.size()];
      const int colorIdx = 4000 + int(i) * 20;
      if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
      target->SetMarkerColor(colorIdx);
      target->SetLineColorAlpha(colorIdx, 0.8);
      target->SetLineWidth(1);
      // target->SetTitleSize(0.09);
      target->SetTitle(typeToParticle[type].c_str());
      target->GetXaxis()->SetTitle(VarName[var].c_str());
      target->GetYaxis()->SetTitle("Count");
      gStyle->SetTitleFillColor(0);
      gStyle->SetTitleBorderSize(0);
      gStyle->SetOptTitle(1);
      // After the first draw, remove the title box
      gStyle->SetTitleX(0.5);
      if (first) {
        target->Draw("HIST");
        first = false;
      } else {
        target->Draw("HIST SAME");
      }

      legend->AddEntry(target, labels[i].c_str(), "l");
    }

    legend->Draw("SAME");
  }

  void PlotPi0KinematicComparison() {
    TCanvas* canvas = new TCanvas("KinematicComparison", "Kinematic Comparison", 1800, 1200);
    canvas->Divide(3, 4);

    std::vector<std::string> types = {"el", "pro", "pho", "pho2"};
    std::vector<std::string> vars = {"p", "theta", "phi"};

    int pad = 1;
    for (const auto& type : types) {
      for (const auto& var : vars) {
        PlotPi0VariableComparison(type, var, pad++, canvas);
      }
    }

    canvas->Update();
    canvas->SaveAs((outputDir + "Pi0KinematicComparison.pdf").c_str());

    std::cout << "Saved pi0 kinematic comparison plots to: " << outputDir + "/Pi0KinematicComparison.pdf" << std::endl;
    delete canvas;
  }

  // Plot one specific variable (e.g., p) for a given particle type (e.g., "el")
  void PlotPi0VariableComparison(const std::string& type, const std::string& var, int padIndex, TCanvas* canvas) {
    canvas->cd(padIndex);
    std::string hname_target = "pi0_rec" + type + "_" + var;

    TLegend* legend = new TLegend(0.6, 0.7, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    bool first = true;
    styleKin_.StylePad((TPad*)gPad);

    for (size_t i = 0; i < plotters.size(); ++i) {
      const auto& histograms = plotters[i]->GetAllHistograms();
      TH1* target = nullptr;

      for (TH1* h : histograms) {
        if (std::string(h->GetName()) == hname_target) {
          target = h;
          break;
        }
      }

      if (!target) {
        std::cerr << "[PlotVariableComparison]: Histogram " << hname_target << " not found for model [" << labels[i] << "]\n";
        continue;
      }
      NormalizeHistogram(target);
      styleKin_.StyleTH1(target);
      auto [cr, cg, cb] = modelShades[i % modelShades.size()];
      const int colorIdx = 4000 + int(i) * 20;
      if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
      target->SetMarkerColor(colorIdx);
      target->SetLineColorAlpha(colorIdx, 0.8);
      target->SetLineWidth(1);
      target->SetTitle(Form("%s;%s;Count", typeToParticle[type].c_str(), VarName[var].c_str()));

      if (first) {
        target->Draw("HIST");
        first = false;
      } else {
        target->Draw("HIST SAME");
      }

      legend->AddEntry(target, labels[i].c_str(), "l");
    }

    legend->Draw();
  }

  // Save an individual variable comparison plot as PNG
  void PlotSingleVariableComparison(const std::string& type, const std::string& var) {
    TCanvas* canvas = new TCanvas(("c_" + type + "_" + var).c_str(), ("Comparison " + type + " " + var).c_str(), 800, 600);
    gPad->SetGrid();

    std::string hname_target = "rec" + type + "_" + var;

    TLegend* legend = new TLegend(0.6, 0.7, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    bool first = true;

    for (size_t i = 0; i < plotters.size(); ++i) {
      const auto& histograms = plotters[i]->GetAllHistograms();
      TH1* target = nullptr;

      for (TH1* h : histograms) {
        if (std::string(h->GetName()) == hname_target) {
          target = h;
          break;
        }
      }

      if (!target) {
        std::cerr << "[PlotSingleVariableComparison]: Histogram " << hname_target << " not found for model [" << labels[i] << "]\n";
        continue;
      }

      target->SetLineColor(i + 1);
      // target->SetTitle(Form("%s;%s;Count", typeToParticle[type].c_str()), VarName[var].c_str());

      if (first) {
        target->Draw("HIST");
        first = false;
      } else {
        target->Draw("HIST SAME");
      }

      legend->AddEntry(target, labels[i].c_str(), "l");
    }

    legend->Draw();
    canvas->Update();
    canvas->SaveAs((outputDir + "/compare_" + type + "_" + var + ".pdf").c_str());
    delete canvas;
  }

  void PlotPhiDVEPKinematicsPlots(bool plotIndividual = false) {
    // Store current global TGaxis state
    int oldMaxDigits = TGaxis::GetMaxDigits();
    // variables to draw (9 one-dimensional plots)
    std::vector<std::string> variables = {"Q2", "xB", "t", "W", "phi", "mtprime", "tmin", "cos_thetaKK", "cos_phiKK"};

    std::map<std::string, std::string> titles = {{"Q2", "Q^{2} [GeV^{2}]"},
                                                 {"xB", "x_{B}"},
                                                 {"t", "-t [GeV^{2}]"},
                                                 {"mtprime", "-t' \\equiv (t_{min}-t) [GeV^{2}]"},
                                                 {"tmin", "t_{min} [GeV^{2}]"},
                                                 {"W", "W [GeV]"},
                                                 {"phi", "#phi [deg]"},
                                                 {"cos_thetaKK", "cos(#theta_{K^{+}K^{-}})"},
                                                 {"cos_phiKK", "cos(#phi_{K^{+}K^{-}})"}};

    // First pass: compute global [min,max] per variable across all plotters for consistent binning
    struct Range {
      double lo{+1e300}, hi{-1e300};
    };
    std::map<std::string, Range> globalRange;

    for (const auto& var : variables) {
      Range R;
      for (size_t i = 0; i < plotters.size(); ++i) {
        auto r = plotters[i]->GetRDF();
        const std::string col = var;
        if (!r.HasColumn(col)) continue;

        double lo = *(r.Min(col));
        double hi = *(r.Max(col));
        if (std::isfinite(lo) && std::isfinite(hi)) {
          R.lo = std::min(R.lo, lo);
          R.hi = std::max(R.hi, hi);
        }
      }
      if (R.lo > R.hi) {  // fallback if empty
        R.lo = 0.0;
        R.hi = 1.0;
      }
      // pad a margin
      double margin = std::max(1e-3, 0.05 * (R.hi - R.lo));
      globalRange[var].lo = R.lo - margin;
      globalRange[var].hi = R.hi + margin;
    }

    // Canvas: 3x4 = 12 pads → 9×1D + 3×2D
    TCanvas* canvas = new TCanvas("PhiEPVars", "Phi electroproduction — kinematic comparison", 1800, 1400);
    canvas->Divide(3, 4);

    int pad = 1;

    // ---------------------------------------------------------------------------
    //  9 × 1D plots
    // ---------------------------------------------------------------------------
    for (const auto& var : variables) {
      canvas->cd(pad++);
      styleDVCS_.StylePad((TPad*)gPad);

      TLegend* legend = new TLegend(0.60, 0.55, 0.88, 0.88);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      legend->SetTextSize(0.04);

      std::vector<TH1D*> histos_to_draw;

      for (size_t i = 0; i < plotters.size(); ++i) {
        auto r = plotters[i]->GetRDF();
        if (!r.HasColumn(var)) {
          std::cerr << "[WARN] Column " << var << " not present after derivation for model " << labels[i] << "\n";
          continue;
        }

        const auto& R = globalRange[var];
        const int nbins = 100;
        const std::string hname = Form("h_%s_%zu_%u", var.c_str(), i, gRandom->Integer(1u << 30));

        auto htmp = r.Histo1D({hname.c_str(), titles[var].c_str(), nbins, R.lo, R.hi}, var);
        auto* h = (TH1D*)htmp->Clone((hname + "_clone").c_str());
        if (!h) continue;

        h->SetDirectory(0);
        NormalizeHistogram(h);
        styleDVCS_.StyleTH1(h);

        auto [cr, cg, cb] = modelShades[i % modelShades.size()];
        const int colorIdx = 4000 + int(i) * 20;
        if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
        h->SetMarkerColor(colorIdx);
        h->SetLineColorAlpha(colorIdx, 0.95);
        h->SetLineWidth(1);
        h->GetXaxis()->SetTitle(titles[var].c_str());
        h->GetYaxis()->SetTitle("Normalized counts");

        histos_to_draw.push_back(h);
        legend->AddEntry(h, labels[i].c_str(), "l");
      }

      for (size_t j = 0; j < histos_to_draw.size(); ++j) {
        histos_to_draw[j]->Draw(j == 0 ? "HIST" : "HIST SAME");
      }
      if (!histos_to_draw.empty()) legend->Draw();

      if (plotIndividual && (var == "xB" || var == "Q2" || var == "t" || var == "W" || var == "phi" || var == "mtprime" || var == "tmin")) {  // <- fixed name
        PlotSingleVariableComparison("el", var);
      }
    }

    // For 2D plots we just use the first plotter (change if you need overlays)
    auto r2d = plotters.front()->GetRDF();

    // Common style lambda
    auto style_h2d = [&]() {
      styleDVCS_.StylePad((TPad*)gPad);
      gPad->SetRightMargin(0.16);
      TH2* h = (TH2*)gPad->GetListOfPrimitives()->At(0);
      if (!h) return;
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
    };

    // ---------------------------------------------------------------------------
    //  3 × 2D plots (pads 10, 11, 12)
    // ---------------------------------------------------------------------------

    // 1) Q2 vs -t' (mtprime)
    if (pad <= 12) {
      canvas->cd(pad++);
      auto h2d = r2d.Histo2D({"h_Q2_vs_tprime", "Q^{2} vs -t';-t' [GeV^{2}];Q^{2} [GeV^{2}]", 60, std::max(0.0, globalRange["mtprime"].lo), globalRange["mtprime"].hi, 60,
                              globalRange["Q2"].lo, globalRange["Q2"].hi},
                             "mtprime", "Q2");
      h2d->DrawCopy("COLZ");
      style_h2d();
    }

    // 2) Q2 vs W
    if (pad <= 12) {
      canvas->cd(pad++);
      auto h2d2 =
          r2d.Histo2D({"h_Q2_vs_W", "Q^{2} vs W;W [GeV];Q^{2} [GeV^{2}]", 60, globalRange["W"].lo, globalRange["W"].hi, 60, globalRange["Q2"].lo, globalRange["Q2"].hi}, "W", "Q2");
      h2d2->DrawCopy("COLZ");
      style_h2d();
    }

    // 3) Q2 vs tmin  (requested last plot)
    if (pad <= 12) {
      canvas->cd(pad++);
      auto h2d3 = r2d.Histo2D({"h_Q2_vs_tmin", "Q^{2} vs t_{min};t_{min} [GeV^{2}];Q^{2} [GeV^{2}]", 60, globalRange["tmin"].lo, globalRange["tmin"].hi, 60, globalRange["Q2"].lo,
                               globalRange["Q2"].hi},
                              "tmin", "Q2");
      h2d3->DrawCopy("COLZ");
      style_h2d();
    }

    // Final save and cleanup
    canvas->SaveAs((outputDir + "/PhiAna_Kinematics_Comparison.pdf").c_str());
    std::cout << "Saved phi electroproduction kinematics comparison to: " << (outputDir + "/PhiAna_Kinematics_Comparison.pdf") << std::endl;

    delete canvas;
    TGaxis::SetMaxDigits(oldMaxDigits);
  }

  void DrawCustomGrid(const std::vector<double>& xB_lines,
                    const std::vector<double>& Q2_lines,
                    double xBmin, double xBmax,
                    double Q2min, double Q2max,
                    int lineStyle=2, int lineWidth=1, int lineColor=kRed)
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


  void PlotDVCSKinematicsComparison(bool plotIndividual = false) {
    // Store current global TGaxis state
    int oldMaxDigits = TGaxis::GetMaxDigits();

    std::vector<std::string> variables = {"Q2", "xB", "t", "W", "phi"};
    std::map<std::string, std::string> titles = {{"Q2", "Q^{2} [GeV^{2}]"}, {"xB", "x_{B}"}, {"t", "-t [GeV^{2}]"}, {"W", "W [GeV]"}, {"phi", "#phi [deg]"}};

    TCanvas* canvas = new TCanvas("DVCSVars", "DVCS Kinematic Comparison", 1800, 1400);
    canvas->Divide(3, 2);

    int pad = 1;
    for (const auto& var : variables) {
      canvas->cd(pad++);
      styleDVCS_.StylePad((TPad*)gPad);

      TLegend* legend = new TLegend(0.6, 0.7, 0.88, 0.88);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);

      bool first = true;
      std::vector<TH1D*> histos_to_draw;

      for (size_t i = 0; i < plotters.size(); ++i) {
        auto rdf = plotters[i]->GetRDF();
        if (!rdf.HasColumn(var)) {
          std::cerr << "[ERROR] Column " << var << " not found in RDF for model " << labels[i] << "\n";
          continue;
        }

        double min = *(rdf.Min(var));
        double max = *(rdf.Max(var));
        if (min == max) {
          min -= 0.1;
          max += 0.1;
        }
        double margin = std::max(1e-3, 0.05 * (max - min));

        // Get histogram (RResultPtr) and clone it
        auto htmp = rdf.Histo1D({Form("h_%s_%zu", var.c_str(), i), titles[var].c_str(), 100, min - margin, max + margin}, var);
        auto h = (TH1D*)htmp->Clone(Form("h_%s_%zu_clone", var.c_str(), i));

        if (!h) continue;  // guard against failed clone

        h->SetDirectory(0);  // prevent ROOT from managing ownership
        NormalizeHistogram(h);
        styleDVCS_.StyleTH1(h);
        h->SetLineColor(i + 2);
        h->SetLineWidth(1);
        h->GetXaxis()->SetTitle(titles[var].c_str());
        h->GetYaxis()->SetTitle("Counts");

        histos_to_draw.push_back(h);
        legend->AddEntry(h, labels[i].c_str(), "l");
      }

      for (size_t j = 0; j < histos_to_draw.size(); ++j) {
        histos_to_draw[j]->Draw(j == 0 ? "HIST" : "HIST SAME");
      }

      if (!histos_to_draw.empty()) {
        legend->Draw();
      }

      if (plotIndividual && (var == "xB" || var == "Q2" || var == "t" || var == "W" || var == "phi")) {
        PlotSingleVariableComparison("el", var);
      }
    }
    canvas->cd(pad);
    auto rdf = plotters.front()->GetRDF();
    auto h2d = rdf.Histo2D({"h_Q2_vs_xB", "Q^{2} vs x_{B};x_{B};Q^{2} [GeV^{2}]", 60, 0, 1.0, 60, 0, 10.0}, "xB", "Q2");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d->GetYaxis()->SetNoExponent(true);
    h2d->SetStats(0);
    h2d->SetTitle("");
    h2d->GetYaxis()->SetLabelFont(42);
    h2d->GetYaxis()->SetLabelSize(0.06);
    h2d->GetYaxis()->SetTitleOffset(1.0);
    h2d->GetYaxis()->SetTitleSize(0.06);
    h2d->GetYaxis()->SetNdivisions(410);

    h2d->GetXaxis()->SetTitleSize(0.065);
    h2d->GetXaxis()->SetLabelFont(42);
    h2d->GetXaxis()->SetLabelSize(0.06);
    h2d->GetXaxis()->SetTitleOffset(0.9);
    h2d->GetXaxis()->SetNdivisions(205);

    h2d->GetZaxis()->SetNdivisions(410);
    h2d->GetZaxis()->SetLabelSize(0.06);
    h2d->GetZaxis()->SetTitleOffset(1.5);
    h2d->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d->DrawCopy("COLZ");
    // Final save and cleanup
    canvas->SaveAs((outputDir + "/DVCS_Kinematics_Comparison.pdf").c_str());
    std::cout << "Saved DVCS kinematics comparison to: " << outputDir + "/DVCS_Kinematics_Comparison.pdf" << std::endl;
    delete canvas;
    TGaxis::SetMaxDigits(oldMaxDigits);
  }

  void PlotDVPi0KinematicsComparison(bool plotIndividual = false) {
    // Store current global TGaxis state
    int oldMaxDigits = TGaxis::GetMaxDigits();

    std::vector<std::string> variables = {"Q2", "xB", "t", "W", "phi"};
    std::map<std::string, std::string> titles = {{"Q2", "Q^{2} [GeV^{2}]"}, {"xB", "x_{B}"}, {"t", "-t [GeV^{2}]"}, {"W", "W [GeV]"}, {"phi", "#phi [deg]"}};

    TCanvas* canvas = new TCanvas("DVCSVars", "DVCS Kinematic Comparison", 1800, 1400);
    canvas->Divide(3, 2);

    int pad = 1;
    for (const auto& var : variables) {
      canvas->cd(pad++);
      styleDVCS_.StylePad((TPad*)gPad);

      TLegend* legend = new TLegend(0.6, 0.7, 0.88, 0.88);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);

      bool first = true;
      std::vector<TH1D*> histos_to_draw;

      for (size_t i = 0; i < plotters.size(); ++i) {
        auto rdf = plotters[i]->GetRDF_Pi0Data();
        if (!rdf.HasColumn(var)) {
          std::cerr << "[ERROR] Column " << var << " not found in RDF for model " << labels[i] << "\n";
          continue;
        }

        double min = *(rdf.Min(var));
        double max = *(rdf.Max(var));
        if (min == max) {
          min -= 0.1;
          max += 0.1;
        }
        double margin = std::max(1e-3, 0.05 * (max - min));

        // Get histogram (RResultPtr) and clone it
        auto htmp = rdf.Histo1D({Form("h_%s_%zu", var.c_str(), i), titles[var].c_str(), 100, min - margin, max + margin}, var);
        auto h = (TH1D*)htmp->Clone(Form("h_%s_%zu_clone", var.c_str(), i));

        if (!h) continue;  // guard against failed clone

        h->SetDirectory(0);  // prevent ROOT from managing ownership
        NormalizeHistogram(h);
        styleDVCS_.StyleTH1(h);
        h->SetLineColor(i + 2);
        h->SetLineWidth(1);
        h->GetXaxis()->SetTitle(titles[var].c_str());
        h->GetYaxis()->SetTitle("Counts");

        histos_to_draw.push_back(h);
        legend->AddEntry(h, labels[i].c_str(), "l");
      }

      for (size_t j = 0; j < histos_to_draw.size(); ++j) {
        histos_to_draw[j]->Draw(j == 0 ? "HIST" : "HIST SAME");
      }

      if (!histos_to_draw.empty()) {
        legend->Draw();
      }

      if (plotIndividual && (var == "xB" || var == "Q2" || var == "t" || var == "W" || var == "phi")) {
        PlotSingleVariableComparison("el", var);
      }
    }
    canvas->cd(pad);
    auto rdf = plotters.front()->GetRDF_Pi0Data();
    auto h2d = rdf.Histo2D({"h_Q2_vs_xB", "Q^{2} vs x_{B};x_{B};Q^{2} [GeV^{2}]", 60, 0, 1.0, 60, 0, 10.0}, "xB", "Q2");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d->GetYaxis()->SetNoExponent(true);
    h2d->SetStats(0);
    h2d->SetTitle("");
    h2d->GetYaxis()->SetLabelFont(42);
    h2d->GetYaxis()->SetLabelSize(0.06);
    h2d->GetYaxis()->SetTitleOffset(1.0);
    h2d->GetYaxis()->SetTitleSize(0.06);
    h2d->GetYaxis()->SetNdivisions(410);

    h2d->GetXaxis()->SetTitleSize(0.065);
    h2d->GetXaxis()->SetLabelFont(42);
    h2d->GetXaxis()->SetLabelSize(0.06);
    h2d->GetXaxis()->SetTitleOffset(0.9);
    h2d->GetXaxis()->SetNdivisions(205);

    h2d->GetZaxis()->SetNdivisions(410);
    h2d->GetZaxis()->SetLabelSize(0.06);
    h2d->GetZaxis()->SetTitleOffset(1.5);
    h2d->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d->DrawCopy("COLZ");
    // Final save and cleanup
    canvas->SaveAs((outputDir + "/DVPi0_Kinematics_Comparison.pdf").c_str());
    std::cout << "Saved DVPi0 kinematics comparison to: " << outputDir + "/DVPi0_Kinematics_Comparison.pdf" << std::endl;
    delete canvas;
    TGaxis::SetMaxDigits(oldMaxDigits);
  }

  void PlotxBQ2tBin(bool plotIndividual = false) {
    // Store current global TGaxis state
    int oldMaxDigits = TGaxis::GetMaxDigits();

    std::vector<std::string> variables = {"Q2", "xB", "t", "W", "phi"};
    std::map<std::string, std::string> titles = {{"Q2", "Q^{2} [GeV^{2}]"}, {"xB", "x_{B}"}, {"t", "-t [GeV^{2}]"}, {"W", "W [GeV]"}, {"phi", "#phi [deg]"}};

    TCanvas* canvas = new TCanvas("xBQ2tBin", "xB-Q2-t-Bin Set", 5400, 1800);
    canvas->Divide(3, 1);

    canvas->cd(1);
    auto rdf = plotters.front()->GetRDF();
    double xBmin = 0.1, xBmax = 0.6;
    double Q2min = 1.0, Q2max = 5.0;
    double tmin = 0.0, tmax = 1.0;
    auto h2d = rdf.Histo2D({"h_Q2_vs_xB", "Q^{2} vs x_{B};x_{B};Q^{2} [GeV^{2}]", 500, xBmin, xBmax, 500, Q2min, Q2max}, "xB", "Q2");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d->GetYaxis()->SetNoExponent(true);
    h2d->SetStats(0);
    h2d->SetTitle("");
    h2d->GetYaxis()->SetLabelFont(42);
    h2d->GetYaxis()->SetLabelSize(0.06);
    h2d->GetYaxis()->SetTitleOffset(1.0);
    h2d->GetYaxis()->SetTitleSize(0.06);
    h2d->GetYaxis()->SetNdivisions(410);

    h2d->GetXaxis()->SetTitleSize(0.065);
    h2d->GetXaxis()->SetLabelFont(42);
    h2d->GetXaxis()->SetLabelSize(0.06);
    h2d->GetXaxis()->SetTitleOffset(0.9);
    h2d->GetXaxis()->SetNdivisions(205);

    h2d->GetZaxis()->SetNdivisions(410);
    h2d->GetZaxis()->SetLabelSize(0.06);
    h2d->GetZaxis()->SetTitleOffset(1.5);
    h2d->GetZaxis()->SetTitleSize(0.06);
    //std::vector<double> xB_lines = {0.150, 0.180, 0.210, 0.240, 0.285, 0.350, 0.430};
    //std::vector<double> Q2_lines = {1.00, 1.25, 1.50, 1.75, 2.00, 2.40, 2.90};
    std::vector<double> xB_lines = {0.125, 0.150, 0.180, 0.210, 0.240, 0.285, 0.350, 0.43};
    std::vector<double> Q2_lines = {1.00, 1.25, 1.50, 1.75, 2.00, 2.40, 2.90, 3.50};
    std::vector<double> t_lines = {0.13, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0};
    TGaxis::SetMaxDigits(3);
    h2d->DrawCopy("COLZ");
    DrawCustomGrid(xB_lines, Q2_lines, xBmin, xBmax, Q2min, Q2max, 1, 1, kRed);
    gPad->RedrawAxis();

    canvas->cd(2);
    auto h2d2 = rdf.Histo2D({"h_Q2_vs_t", "Q^{2} vs -t;-t[GeV^{2}];Q^{2} [GeV^{2}]", 500, tmin, tmax, 500, Q2min, Q2max}, "t", "Q2");
    
    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d2->GetYaxis()->SetNoExponent(true);
    h2d2->SetStats(0);
    h2d2->SetTitle("");
    h2d2->GetYaxis()->SetLabelFont(42);
    h2d2->GetYaxis()->SetLabelSize(0.06);
    h2d2->GetYaxis()->SetTitleOffset(1.0);
    h2d2->GetYaxis()->SetTitleSize(0.06);
    h2d2->GetYaxis()->SetNdivisions(410);
    h2d2->GetXaxis()->SetTitleSize(0.065);
    h2d2->GetXaxis()->SetLabelFont(42);
    h2d2->GetXaxis()->SetLabelSize(0.06);
    h2d2->GetXaxis()->SetTitleOffset(0.9);
    h2d2->GetXaxis()->SetNdivisions(205);
    h2d2->GetZaxis()->SetNdivisions(410);
    h2d2->GetZaxis()->SetLabelSize(0.06);
    h2d2->GetZaxis()->SetTitleOffset(1.5);
    h2d2->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d2->DrawCopy("COLZ");
    DrawCustomGrid(t_lines, Q2_lines, tmin, tmax, Q2min, Q2max, 1, 1, kRed);
    gPad->RedrawAxis();

    canvas->cd(3);
    auto h2d3 = rdf.Histo2D({"h_xB_vs_t", "x_{B} vs -t;-t[GeV^{2}];x_{B}", 500, tmin, tmax, 500, xBmin, xBmax}, "t", "xB");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d3->GetYaxis()->SetNoExponent(true);
    h2d3->SetStats(0);
    h2d3->SetTitle("");
    h2d3->GetYaxis()->SetLabelFont(42);
    h2d3->GetYaxis()->SetLabelSize(0.06);
    h2d3->GetYaxis()->SetTitleOffset(1.0);
    h2d3->GetYaxis()->SetTitleSize(0.06);
    h2d3->GetYaxis()->SetNdivisions(410);
    h2d3->GetXaxis()->SetTitleSize(0.065);
    h2d3->GetXaxis()->SetLabelFont(42);
    h2d3->GetXaxis()->SetLabelSize(0.06);
    h2d3->GetXaxis()->SetTitleOffset(0.9);
    h2d3->GetXaxis()->SetNdivisions(205);
    h2d3->GetZaxis()->SetNdivisions(410);
    h2d3->GetZaxis()->SetLabelSize(0.06);
    h2d3->GetZaxis()->SetTitleOffset(1.5);
    h2d3->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d3->DrawCopy("COLZ");
    DrawCustomGrid(t_lines, xB_lines, tmin, tmax, xBmin, xBmax, 1, 1, kRed);
    gPad->RedrawAxis();
    // Final save and cleanup
    canvas->SaveAs((outputDir + "/xBQ2tBin.pdf").c_str());
    std::cout << "Saved xBQ2tBin kinematics to: " << outputDir + "/xBQ2tBin.pdf" << std::endl;
    delete canvas;
    TGaxis::SetMaxDigits(oldMaxDigits);
  }

  void PlotxBQ2tBinMC(bool plotIndividual = false) {
    // Store current global TGaxis state
    int oldMaxDigits = TGaxis::GetMaxDigits();

    std::vector<std::string> variables = {"Q2", "xB", "t", "W", "phi"};
    std::map<std::string, std::string> titles = {{"Q2", "Q^{2} [GeV^{2}]"}, {"xB", "x_{B}"}, {"t", "-t [GeV^{2}]"}, {"W", "W [GeV]"}, {"phi", "#phi [deg]"}};

    TCanvas* canvas = new TCanvas("xBQ2tBin", "xB-Q2-t-Bin Set", 5400, 1800);
    canvas->Divide(3, 1);

    canvas->cd(1);
    auto rdf = plotters.front()->GetRDF_DVCSMC();
    double xBmin = 0.1, xBmax = 0.6;
    double Q2min = 1.0, Q2max = 5.0;
    double tmin = 0.0, tmax = 1.0;
    auto h2d = rdf.Histo2D({"h_Q2_vs_xB", "Q^{2} vs x_{B};x_{B};Q^{2} [GeV^{2}]", 500, xBmin, xBmax, 500, Q2min, Q2max}, "xB", "Q2");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d->GetYaxis()->SetNoExponent(true);
    h2d->SetStats(0);
    h2d->SetTitle("");
    h2d->GetYaxis()->SetLabelFont(42);
    h2d->GetYaxis()->SetLabelSize(0.06);
    h2d->GetYaxis()->SetTitleOffset(1.0);
    h2d->GetYaxis()->SetTitleSize(0.06);
    h2d->GetYaxis()->SetNdivisions(410);

    h2d->GetXaxis()->SetTitleSize(0.065);
    h2d->GetXaxis()->SetLabelFont(42);
    h2d->GetXaxis()->SetLabelSize(0.06);
    h2d->GetXaxis()->SetTitleOffset(0.9);
    h2d->GetXaxis()->SetNdivisions(205);

    h2d->GetZaxis()->SetNdivisions(410);
    h2d->GetZaxis()->SetLabelSize(0.06);
    h2d->GetZaxis()->SetTitleOffset(1.5);
    h2d->GetZaxis()->SetTitleSize(0.06);
    //std::vector<double> xB_lines = {0.150, 0.180, 0.210, 0.240, 0.285, 0.350, 0.430};
    //std::vector<double> Q2_lines = {1.00, 1.25, 1.50, 1.75, 2.00, 2.40, 2.90};
    std::vector<double> xB_lines = {0.125, 0.150, 0.180, 0.210, 0.240, 0.285, 0.350, 0.43};
    std::vector<double> Q2_lines = {1.00, 1.25, 1.50, 1.75, 2.00, 2.40, 2.90, 3.50};
    std::vector<double> t_lines = {0.13, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0};
    TGaxis::SetMaxDigits(3);
    h2d->DrawCopy("COLZ");
    DrawCustomGrid(xB_lines, Q2_lines, xBmin, xBmax, Q2min, Q2max, 1, 1, kRed);
    gPad->RedrawAxis();

    canvas->cd(2);
    auto h2d2 = rdf.Histo2D({"h_Q2_vs_t", "Q^{2} vs -t;-t[GeV^{2}];Q^{2} [GeV^{2}]", 500, tmin, tmax, 500, Q2min, Q2max}, "t", "Q2");
    
    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d2->GetYaxis()->SetNoExponent(true);
    h2d2->SetStats(0);
    h2d2->SetTitle("");
    h2d2->GetYaxis()->SetLabelFont(42);
    h2d2->GetYaxis()->SetLabelSize(0.06);
    h2d2->GetYaxis()->SetTitleOffset(1.0);
    h2d2->GetYaxis()->SetTitleSize(0.06);
    h2d2->GetYaxis()->SetNdivisions(410);
    h2d2->GetXaxis()->SetTitleSize(0.065);
    h2d2->GetXaxis()->SetLabelFont(42);
    h2d2->GetXaxis()->SetLabelSize(0.06);
    h2d2->GetXaxis()->SetTitleOffset(0.9);
    h2d2->GetXaxis()->SetNdivisions(205);
    h2d2->GetZaxis()->SetNdivisions(410);
    h2d2->GetZaxis()->SetLabelSize(0.06);
    h2d2->GetZaxis()->SetTitleOffset(1.5);
    h2d2->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d2->DrawCopy("COLZ");
    DrawCustomGrid(t_lines, Q2_lines, tmin, tmax, Q2min, Q2max, 1, 1, kRed);
    gPad->RedrawAxis();

    canvas->cd(3);
    auto h2d3 = rdf.Histo2D({"h_xB_vs_t", "x_{B} vs -t;-t[GeV^{2}];x_{B}", 500, tmin, tmax, 500, xBmin, xBmax}, "t", "xB");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d3->GetYaxis()->SetNoExponent(true);
    h2d3->SetStats(0);
    h2d3->SetTitle("");
    h2d3->GetYaxis()->SetLabelFont(42);
    h2d3->GetYaxis()->SetLabelSize(0.06);
    h2d3->GetYaxis()->SetTitleOffset(1.0);
    h2d3->GetYaxis()->SetTitleSize(0.06);
    h2d3->GetYaxis()->SetNdivisions(410);
    h2d3->GetXaxis()->SetTitleSize(0.065);
    h2d3->GetXaxis()->SetLabelFont(42);
    h2d3->GetXaxis()->SetLabelSize(0.06);
    h2d3->GetXaxis()->SetTitleOffset(0.9);
    h2d3->GetXaxis()->SetNdivisions(205);
    h2d3->GetZaxis()->SetNdivisions(410);
    h2d3->GetZaxis()->SetLabelSize(0.06);
    h2d3->GetZaxis()->SetTitleOffset(1.5);
    h2d3->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d3->DrawCopy("COLZ");
    DrawCustomGrid(t_lines, xB_lines, tmin, tmax, xBmin, xBmax, 1, 1, kRed);
    gPad->RedrawAxis();
    // Final save and cleanup
    canvas->SaveAs((outputDir + "/xBQ2tBinMC.pdf").c_str());
    std::cout << "Saved xBQ2tBinMC kinematics to: " << outputDir + "/xBQ2tBinMC.pdf" << std::endl;
    delete canvas;
    TGaxis::SetMaxDigits(oldMaxDigits);
  }

  void PlotxBQ2tBinPi0(bool plotIndividual = false) {
    // Store current global TGaxis state
    int oldMaxDigits = TGaxis::GetMaxDigits();

    std::vector<std::string> variables = {"Q2", "xB", "t", "W", "phi"};
    std::map<std::string, std::string> titles = {{"Q2", "Q^{2} [GeV^{2}]"}, {"xB", "x_{B}"}, {"t", "-t [GeV^{2}]"}, {"W", "W [GeV]"}, {"phi", "#phi [deg]"}};

    TCanvas* canvas = new TCanvas("xBQ2tBin", "xB-Q2-t-Bin Set", 5400, 1800);
    canvas->Divide(3, 1);

    canvas->cd(1);
    auto rdf = plotters.front()->GetRDF_Pi0Data();
    double xBmin = 0.1, xBmax = 0.6;
    double Q2min = 1.0, Q2max = 5.0;
    double tmin = 0.0, tmax = 1.0;
    auto h2d = rdf.Histo2D({"h_Q2_vs_xB", "Q^{2} vs x_{B};x_{B};Q^{2} [GeV^{2}]", 500, xBmin, xBmax, 500, Q2min, Q2max}, "xB", "Q2");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d->GetYaxis()->SetNoExponent(true);
    h2d->SetStats(0);
    h2d->SetTitle("");
    h2d->GetYaxis()->SetLabelFont(42);
    h2d->GetYaxis()->SetLabelSize(0.06);
    h2d->GetYaxis()->SetTitleOffset(1.0);
    h2d->GetYaxis()->SetTitleSize(0.06);
    h2d->GetYaxis()->SetNdivisions(410);

    h2d->GetXaxis()->SetTitleSize(0.065);
    h2d->GetXaxis()->SetLabelFont(42);
    h2d->GetXaxis()->SetLabelSize(0.06);
    h2d->GetXaxis()->SetTitleOffset(0.9);
    h2d->GetXaxis()->SetNdivisions(205);

    h2d->GetZaxis()->SetNdivisions(410);
    h2d->GetZaxis()->SetLabelSize(0.06);
    h2d->GetZaxis()->SetTitleOffset(1.5);
    h2d->GetZaxis()->SetTitleSize(0.06);
    //std::vector<double> xB_lines = {0.150, 0.180, 0.210, 0.240, 0.285, 0.350, 0.430};
    //std::vector<double> Q2_lines = {1.00, 1.25, 1.50, 1.75, 2.00, 2.40, 2.90};
    std::vector<double> xB_lines = {0.125, 0.150, 0.180, 0.210, 0.240, 0.285, 0.350, 0.43};
    std::vector<double> Q2_lines = {1.00, 1.25, 1.50, 1.75, 2.00, 2.40, 2.90, 3.50};
    std::vector<double> t_lines = {0.13, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0};
    TGaxis::SetMaxDigits(3);
    h2d->DrawCopy("COLZ");
    DrawCustomGrid(xB_lines, Q2_lines, xBmin, xBmax, Q2min, Q2max, 1, 1, kRed);
    gPad->RedrawAxis();

    canvas->cd(2);
    auto h2d2 = rdf.Histo2D({"h_Q2_vs_t", "Q^{2} vs -t;-t[GeV^{2}];Q^{2} [GeV^{2}]", 500, tmin, tmax, 500, Q2min, Q2max}, "t", "Q2");
    
    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d2->GetYaxis()->SetNoExponent(true);
    h2d2->SetStats(0);
    h2d2->SetTitle("");
    h2d2->GetYaxis()->SetLabelFont(42);
    h2d2->GetYaxis()->SetLabelSize(0.06);
    h2d2->GetYaxis()->SetTitleOffset(1.0);
    h2d2->GetYaxis()->SetTitleSize(0.06);
    h2d2->GetYaxis()->SetNdivisions(410);
    h2d2->GetXaxis()->SetTitleSize(0.065);
    h2d2->GetXaxis()->SetLabelFont(42);
    h2d2->GetXaxis()->SetLabelSize(0.06);
    h2d2->GetXaxis()->SetTitleOffset(0.9);
    h2d2->GetXaxis()->SetNdivisions(205);
    h2d2->GetZaxis()->SetNdivisions(410);
    h2d2->GetZaxis()->SetLabelSize(0.06);
    h2d2->GetZaxis()->SetTitleOffset(1.5);
    h2d2->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d2->DrawCopy("COLZ");
    DrawCustomGrid(t_lines, Q2_lines, tmin, tmax, Q2min, Q2max, 1, 1, kRed);
    gPad->RedrawAxis();

    canvas->cd(3);
    auto h2d3 = rdf.Histo2D({"h_xB_vs_t", "x_{B} vs -t;-t[GeV^{2}];x_{B}", 500, tmin, tmax, 500, xBmin, xBmax}, "t", "xB");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d3->GetYaxis()->SetNoExponent(true);
    h2d3->SetStats(0);
    h2d3->SetTitle("");
    h2d3->GetYaxis()->SetLabelFont(42);
    h2d3->GetYaxis()->SetLabelSize(0.06);
    h2d3->GetYaxis()->SetTitleOffset(1.0);
    h2d3->GetYaxis()->SetTitleSize(0.06);
    h2d3->GetYaxis()->SetNdivisions(410);
    h2d3->GetXaxis()->SetTitleSize(0.065);
    h2d3->GetXaxis()->SetLabelFont(42);
    h2d3->GetXaxis()->SetLabelSize(0.06);
    h2d3->GetXaxis()->SetTitleOffset(0.9);
    h2d3->GetXaxis()->SetNdivisions(205);
    h2d3->GetZaxis()->SetNdivisions(410);
    h2d3->GetZaxis()->SetLabelSize(0.06);
    h2d3->GetZaxis()->SetTitleOffset(1.5);
    h2d3->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d3->DrawCopy("COLZ");
    DrawCustomGrid(t_lines, xB_lines, tmin, tmax, xBmin, xBmax, 1, 1, kRed);
    gPad->RedrawAxis();
    // Final save and cleanup
    canvas->SaveAs((outputDir + "/xBQ2tBinPi0.pdf").c_str());
    std::cout << "Saved xBQ2tBinPi0 kinematics to: " << outputDir + "/xBQ2tBinPi0.pdf" << std::endl;
    delete canvas;
    TGaxis::SetMaxDigits(oldMaxDigits);
  }


  void PlotxBQ2tBinPi0MC(bool plotIndividual = false) {
    // Store current global TGaxis state
    int oldMaxDigits = TGaxis::GetMaxDigits();

    std::vector<std::string> variables = {"Q2", "xB", "t", "W", "phi"};
    std::map<std::string, std::string> titles = {{"Q2", "Q^{2} [GeV^{2}]"}, {"xB", "x_{B}"}, {"t", "-t [GeV^{2}]"}, {"W", "W [GeV]"}, {"phi", "#phi [deg]"}};

    TCanvas* canvas = new TCanvas("xBQ2tBin", "xB-Q2-t-Bin Set", 5400, 1800);
    canvas->Divide(3, 1);

    canvas->cd(1);
    auto rdf = plotters.front()->GetRDF_Pi0MC();
    double xBmin = 0.1, xBmax = 0.6;
    double Q2min = 1.0, Q2max = 5.0;
    double tmin = 0.0, tmax = 1.0;
    auto h2d = rdf.Histo2D({"h_Q2_vs_xB", "Q^{2} vs x_{B};x_{B};Q^{2} [GeV^{2}]", 500, xBmin, xBmax, 500, Q2min, Q2max}, "xB", "Q2");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d->GetYaxis()->SetNoExponent(true);
    h2d->SetStats(0);
    h2d->SetTitle("");
    h2d->GetYaxis()->SetLabelFont(42);
    h2d->GetYaxis()->SetLabelSize(0.06);
    h2d->GetYaxis()->SetTitleOffset(1.0);
    h2d->GetYaxis()->SetTitleSize(0.06);
    h2d->GetYaxis()->SetNdivisions(410);

    h2d->GetXaxis()->SetTitleSize(0.065);
    h2d->GetXaxis()->SetLabelFont(42);
    h2d->GetXaxis()->SetLabelSize(0.06);
    h2d->GetXaxis()->SetTitleOffset(0.9);
    h2d->GetXaxis()->SetNdivisions(205);

    h2d->GetZaxis()->SetNdivisions(410);
    h2d->GetZaxis()->SetLabelSize(0.06);
    h2d->GetZaxis()->SetTitleOffset(1.5);
    h2d->GetZaxis()->SetTitleSize(0.06);
    //std::vector<double> xB_lines = {0.150, 0.180, 0.210, 0.240, 0.285, 0.350, 0.430};
    //std::vector<double> Q2_lines = {1.00, 1.25, 1.50, 1.75, 2.00, 2.40, 2.90};
    std::vector<double> xB_lines = {0.125, 0.150, 0.180, 0.210, 0.240, 0.285, 0.350, 0.43};
    std::vector<double> Q2_lines = {1.00, 1.25, 1.50, 1.75, 2.00, 2.40, 2.90, 3.50};
    std::vector<double> t_lines = {0.13, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0};
    TGaxis::SetMaxDigits(3);
    h2d->DrawCopy("COLZ");
    DrawCustomGrid(xB_lines, Q2_lines, xBmin, xBmax, Q2min, Q2max, 1, 1, kRed);
    gPad->RedrawAxis();

    canvas->cd(2);
    auto h2d2 = rdf.Histo2D({"h_Q2_vs_t", "Q^{2} vs -t;-t[GeV^{2}];Q^{2} [GeV^{2}]", 500, tmin, tmax, 500, Q2min, Q2max}, "t", "Q2");
    
    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d2->GetYaxis()->SetNoExponent(true);
    h2d2->SetStats(0);
    h2d2->SetTitle("");
    h2d2->GetYaxis()->SetLabelFont(42);
    h2d2->GetYaxis()->SetLabelSize(0.06);
    h2d2->GetYaxis()->SetTitleOffset(1.0);
    h2d2->GetYaxis()->SetTitleSize(0.06);
    h2d2->GetYaxis()->SetNdivisions(410);
    h2d2->GetXaxis()->SetTitleSize(0.065);
    h2d2->GetXaxis()->SetLabelFont(42);
    h2d2->GetXaxis()->SetLabelSize(0.06);
    h2d2->GetXaxis()->SetTitleOffset(0.9);
    h2d2->GetXaxis()->SetNdivisions(205);
    h2d2->GetZaxis()->SetNdivisions(410);
    h2d2->GetZaxis()->SetLabelSize(0.06);
    h2d2->GetZaxis()->SetTitleOffset(1.5);
    h2d2->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d2->DrawCopy("COLZ");
    DrawCustomGrid(t_lines, Q2_lines, tmin, tmax, Q2min, Q2max, 1, 1, kRed);
    gPad->RedrawAxis();

    canvas->cd(3);
    auto h2d3 = rdf.Histo2D({"h_xB_vs_t", "x_{B} vs -t;-t[GeV^{2}];x_{B}", 500, tmin, tmax, 500, xBmin, xBmax}, "t", "xB");

    styleDVCS_.StylePad((TPad*)gPad);
    gPad->SetRightMargin(0.16);
    h2d3->GetYaxis()->SetNoExponent(true);
    h2d3->SetStats(0);
    h2d3->SetTitle("");
    h2d3->GetYaxis()->SetLabelFont(42);
    h2d3->GetYaxis()->SetLabelSize(0.06);
    h2d3->GetYaxis()->SetTitleOffset(1.0);
    h2d3->GetYaxis()->SetTitleSize(0.06);
    h2d3->GetYaxis()->SetNdivisions(410);
    h2d3->GetXaxis()->SetTitleSize(0.065);
    h2d3->GetXaxis()->SetLabelFont(42);
    h2d3->GetXaxis()->SetLabelSize(0.06);
    h2d3->GetXaxis()->SetTitleOffset(0.9);
    h2d3->GetXaxis()->SetNdivisions(205);
    h2d3->GetZaxis()->SetNdivisions(410);
    h2d3->GetZaxis()->SetLabelSize(0.06);
    h2d3->GetZaxis()->SetTitleOffset(1.5);
    h2d3->GetZaxis()->SetTitleSize(0.06);
    TGaxis::SetMaxDigits(3);
    h2d3->DrawCopy("COLZ");
    DrawCustomGrid(t_lines, xB_lines, tmin, tmax, xBmin, xBmax, 1, 1, kRed);
    gPad->RedrawAxis();
    // Final save and cleanup
    canvas->SaveAs((outputDir + "/xBQ2tBinPi0MC.pdf").c_str());
    std::cout << "Saved xBQ2tBinPi0MC kinematics to: " << outputDir + "/xBQ2tBinPi0MC.pdf" << std::endl;
    delete canvas;
    TGaxis::SetMaxDigits(oldMaxDigits);
  }
  /// For exclusivity cuts, you can use the following function to select one triplet
  void PlotExclusivityComparisonByDetectorCases(const std::vector<std::pair<std::string, std::string>>& detectorCuts) {
    std::vector<std::tuple<std::string, std::string, std::string, double, double>> vars = {
        {"Mx2_ep", "Missing Mass Squared (ep)", "MM^{2}(ep) [GeV^{2}]", -0.6, 0.6},
        {"Emiss", "Missing Energy", "E_{miss} [GeV]", -1.0, 2.0},
        {"PTmiss", "Transverse Missing Momentum", "P_{T}^{miss} [GeV/c]", -0.1, 0.4},
        {"Theta_gamma_gamma", "#theta(#gamma, #vec{q})", "#theta_{#gamma#gamma'} [deg]", -2.0, 4.0},
        {"DeltaPhi", "Coplanarity Angle", "#Delta#phi [deg]", 0, 20},
        {"Mx2_epg", "Missing Mass Squared (ep#gamma)", "MM^{2}(ep#gamma) [GeV^{2}]", -0.05, 0.05},
        {"Mx2_eg", "Invariant Mass (e#gamma)", "M^{2}(e#gamma) [GeV^{2}]", -0.5, 3},
        {"Theta_e_gamma", "Angle: e-#gamma", "#theta(e, #gamma) [deg]", 0.0, 60.0},
        {"DeltaE", "Energy Balance", "#DeltaE [GeV]", -1.0, 2.0}};

    for (const auto& [cutExpr, cutLabel] : detectorCuts) {
      std::string cleanName = cutLabel;
      std::replace(cleanName.begin(), cleanName.end(), ' ', '_');
      std::replace(cleanName.begin(), cleanName.end(), ',', '_');

      TCanvas* canvas = new TCanvas(("c_" + cleanName).c_str(), cutLabel.c_str(), 1800, 1200);
      int cols = 3;
      int rows = (vars.size() + cols - 1) / cols;
      canvas->Divide(cols, rows);

      for (size_t i = 0; i < vars.size(); ++i) {
        canvas->cd(i + 1);
        const auto& [var, title, xlabel, xmin, xmax] = vars[i];
        gPad->SetTicks();
        styleKin_.StylePad((TPad*)gPad);

        TLegend* legend = new TLegend(0.6, 0.7, 0.88, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.04);

        bool first = true;

        for (size_t m = 0; m < plotters.size(); ++m) {
          auto rdf_cut = plotters[m]->GetRDF().Filter(cutExpr, cutLabel);
          if (!rdf_cut.HasColumn(var)) continue;

          auto h = rdf_cut.Histo1D({Form("h_%s_%s_%zu", var.c_str(), cleanName.c_str(), m), (title + ";" + xlabel + ";Counts").c_str(), 100, xmin, xmax}, var);
          h.GetValue();

          TH1D* h_clone = (TH1D*)h.GetPtr()->Clone();
          h_clone->SetDirectory(0);
          NormalizeHistogram(h_clone);

          styleKin_.StyleTH1(h_clone);
          h_clone->SetLineColor(m + 2);
          h_clone->SetLineWidth(2);

          double mean = h_clone->GetMean();
          double sigma = h_clone->GetStdDev();
          double x1 = mean - 3 * sigma;
          double x2 = mean + 3 * sigma;

          TLine* line1 = new TLine(x1, 0, x1, h_clone->GetMaximum() * 0.5);
          TLine* line2 = new TLine(x2, 0, x2, h_clone->GetMaximum() * 0.5);
          line1->SetLineColor(m + 2);
          line2->SetLineColor(m + 2);
          line1->SetLineStyle(2);  // Dashed
          line2->SetLineStyle(2);

          if (first) {
            h_clone->Draw("HIST");
            first = false;
          } else {
            h_clone->Draw("HIST SAME");
          }

          legend->AddEntry(h_clone, labels[m].c_str(), "l");
          std::ostringstream stats;
          stats << "#mu = " << std::fixed << std::setprecision(2) << mean << ", #sigma = " << std::fixed << std::setprecision(2) << sigma;
          legend->AddEntry((TObject*)0, stats.str().c_str(), "");
          line1->Draw("SAME");
          line2->Draw("SAME");
        }

        legend->Draw();
      }

      std::string outpath = outputDir + "/Exclusivity_" + cleanName + ".pdf";
      canvas->SaveAs(outpath.c_str());
      std::cout << "Saved detector-specific comparison to: " << outpath << "\n";
      delete canvas;
    }
  };

  void PlotPi0ExclusivityComparisonByDetectorCases(const std::vector<std::pair<std::string, std::string>>& detectorCuts) {
    std::vector<std::tuple<std::string, std::string, std::string, double, double>> vars = {
        {"Mass_pi0", "Pi0 Mass", "M_{#pi0} [GeV]", 0.07, 0.2},
        {"Emiss_pi0", "Missing Energy", "E_{miss} [GeV]", -1.0, 1.0},
        {"PTmiss_pi0", "Transverse Missing Momentum", "P_{T}^{miss} [GeV/c]", 0.0, 0.3},
        {"Theta_pi0pi0", "#theta(#pi0, #pi')", "#theta_{#pi#pi'} [deg]", 0.0, 3.0},
        {"DeltaPhi_pi0", "Coplanarity Angle", "#Delta#phi [deg]", 0.0, 15.0},
        {"Mx2_eppi0", "Missing Mass Squared (ep#pi)", "MM^{2}(ep#pi) [GeV^{2}]", -0.03, 0.03},
        {"Mx2_ep_pi0", "Missing Mass Squared (ep#pi)", "MM^{2}(ep) [GeV^{2}]", -0.4, 0.6},
        {"Mx2_epi0", "Missing Mass Squared (e#pi)", "MM^{2}(e#pi) [GeV^{2}]", -0.2, 2.0},
        {"Theta_epho1", "Angle: e-#gamma_{1}", "#theta(e, #gamma_{1}) [deg]", 0.0, 60.0},
        {"Theta_epho2", "Angle: e-#gamma_{2}", "#theta(e, #gamma_{2}) [deg]", 0.0, 60.0},
        {"Theta_pho1pho2", "Angle: #gamma_{1}-#gamma_{2}", "#theta(#gamma_{1}, #gamma_{2}) [deg]", 0.0, 15.0}};

    for (const auto& [cutExpr, cutLabel] : detectorCuts) {
      std::string cleanName = cutLabel;
      std::replace(cleanName.begin(), cleanName.end(), ' ', '_');
      std::replace(cleanName.begin(), cleanName.end(), ',', '_');

      TCanvas* canvas = new TCanvas(("c_" + cleanName).c_str(), cutLabel.c_str(), 4800, 3200);
      int cols = 4;
      int rows = (vars.size() + cols - 1) / cols;
      canvas->Divide(cols, rows);

      for (size_t i = 0; i < vars.size(); ++i) {
        canvas->cd(i + 1);
        const auto& [var, title, xlabel, xmin, xmax] = vars[i];

        gPad->SetTicks();
        styleKin_.StylePad((TPad*)gPad);

        TLegend* legend = new TLegend(0.50, 0.70, 0.88, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.040);

        bool first = true;

        for (size_t m = 0; m < plotters.size(); ++m) {
          auto rdf_data = plotters[m]->GetRDF_Pi0Data().Filter(cutExpr, cutLabel);
          auto rdf_mc = plotters[m]->GetRDF_DVCSPi0MC().Filter(cutExpr, cutLabel);

          if (!rdf_data.HasColumn(var)) continue;

          double xlo = xmin, xhi = xmax;

          auto h_data = rdf_data.Histo1D({Form("h_data_%s_%s_%zu", var.c_str(), cleanName.c_str(), m), (title + ";" + xlabel + ";Counts").c_str(), 100, xlo, xhi}, var);
          h_data.GetValue();  // materialize
          TH1D* hD = (TH1D*)h_data.GetPtr()->Clone();
          hD->SetDirectory(0);
          NormalizeHistogram(hD);
          styleKin_.StyleTH1(hD);
          hD->SetLineColor(m + 2);
          hD->SetLineWidth(2);
          hD->SetLineStyle(1);

          TH1D* hM = nullptr;
          if (rdf_mc.HasColumn(var)) {
            auto h_mc = rdf_mc.Histo1D({Form("h_mc_%s_%s_%zu", var.c_str(), cleanName.c_str(), m), (title + ";" + xlabel + ";Counts").c_str(), 100, xmin, xmax}, var);
            h_mc.GetValue();
            hM = (TH1D*)h_mc.GetPtr()->Clone();
            hM->SetDirectory(0);
            NormalizeHistogram(hM);
            styleKin_.StyleTH1(hM);
            hM->SetLineColor(m + 2);
            hM->SetLineWidth(2);
            hM->SetLineStyle(2);
          }

          if (first) {
            hD->Draw("HIST");
            if (hM) hM->Draw("HIST SAME");
            first = false;
          } else {
            hD->Draw("HIST SAME");
            if (hM) hM->Draw("HIST SAME");
          }
          double maxY = hD->GetMaximum();
          if (hM) maxY = std::max(maxY, hM->GetMaximum());
          hD->SetMaximum(maxY * 1.2);
          gPad->Modified();
          gPad->Update();

          legend->AddEntry(hD, (labels[m] + " Data").c_str(), "l");
          if (hM) legend->AddEntry(hM, (labels[m] + " MC").c_str(), "l");

          double mean = hD->GetMean();
          double sigma = hD->GetStdDev();
          double x1 = mean - 3.0 * sigma;
          double x2 = mean + 3.0 * sigma;

          double meanM = hM->GetMean();
          double sigmaM = hM->GetStdDev();
          double x1M = meanM - 3.0 * sigmaM;
          double x2M = meanM + 3.0 * sigmaM;

          double ymax = std::max(hD->GetMaximum(), hM ? hM->GetMaximum() : 0.0);
          TLine* line1 = new TLine(x1, 0.0, x1, ymax * 0.5);
          TLine* line2 = new TLine(x2, 0.0, x2, ymax * 0.5);
          line1->SetLineColor(m + 2);
          line2->SetLineColor(m + 2);
          line1->SetLineStyle(3);
          line2->SetLineStyle(3);
          line1->Draw("SAME");
          line2->Draw("SAME");

          TLine* line1M = new TLine(x1M, 0.0, x1M, ymax * 0.5);
          TLine* line2M = new TLine(x2M, 0.0, x2M, ymax * 0.5);
          line1M->SetLineColor(m + 2);
          line2M->SetLineColor(m + 2);
          line1M->SetLineStyle(3);
          line2M->SetLineStyle(3);
          line1M->Draw("SAME");
          line2M->Draw("SAME");

          std::ostringstream stats;
          stats << "Data #mu = " << std::fixed << std::setprecision(3) << mean << ", #sigma = " << std::fixed << std::setprecision(3) << sigma;
          legend->AddEntry((TObject*)nullptr, stats.str().c_str(), "");

          std::ostringstream statsM;
          statsM << "MC #mu = " << std::fixed << std::setprecision(3) << meanM << ", #sigma = " << std::fixed << std::setprecision(3) << sigmaM;
          legend->AddEntry((TObject*)nullptr, statsM.str().c_str(), "");
        }

        legend->Draw();
      }

      std::string outpath = outputDir + "/Pi0Exclusivity_" + cleanName + ".pdf";
      canvas->SaveAs(outpath.c_str());
      std::cout << "Saved detector-specific comparison to: " << outpath << "\n";
      delete canvas;
    }
  }

  void PlotPi0ExclusivityComparisonByDetectorCases_old(const std::vector<std::pair<std::string, std::string>>& detectorCuts) {
    std::vector<std::tuple<std::string, std::string, std::string, double, double>> vars = {{"Mass_pi0", "Pi0 Mass", "M_{#pi0} [GeV]", -0.6, 0.6},
                                                                                           {"Emiss_pi0", "Missing Energy", "E_{miss} [GeV]", -1.0, 2.0},
                                                                                           {"PTmiss_pi0", "Transverse Missing Momentum", "P_{T}^{miss} [GeV/c]", -0.1, 0.4},
                                                                                           {"Theta_pi0pi0", "#theta(#pi0, #pi')", "#theta_{#pi#pi'} [deg]", -2.0, 4.0},
                                                                                           {"DeltaPhi_pi0", "Coplanarity Angle", "#Delta#phi [deg]", 0, 20},
                                                                                           {"Mx2_eppi0", "Missing Mass Squared (ep#pi)", "MM^{2}(ep#pi) [GeV^{2}]", -0.05, 0.05},
                                                                                           {"Mx2_ep_pi0", "Missing Mass Squared (ep#pi)", "MM^{2}(ep) [GeV^{2}]", -0.5, 1},
                                                                                           {"Mx2_epi0", "Missing Mass Squared (e#pi)", "MM^{2}(e#pi) [GeV^{2}]", 0.0, 2.0},
                                                                                           {"Theta_epho1", "Angle: e-#gamma_{1}", "#theta(e, #gamma_{1}) [deg]", 0.0, 60.0},
                                                                                           {"Theta_epho2", "Angle: e-#gamma_{2}", "#theta(e, #gamma_{2}) [deg]", 0.0, 60.0}};

    for (const auto& [cutExpr, cutLabel] : detectorCuts) {
      std::string cleanName = cutLabel;
      std::replace(cleanName.begin(), cleanName.end(), ' ', '_');
      std::replace(cleanName.begin(), cleanName.end(), ',', '_');

      TCanvas* canvas = new TCanvas(("c_" + cleanName).c_str(), cutLabel.c_str(), 1800, 1200);
      int cols = 4;
      int rows = (vars.size() + cols - 1) / cols;
      canvas->Divide(cols, rows);

      for (size_t i = 0; i < vars.size(); ++i) {
        canvas->cd(i + 1);
        const auto& [var, title, xlabel, xmin, xmax] = vars[i];
        gPad->SetTicks();
        styleKin_.StylePad((TPad*)gPad);

        TLegend* legend = new TLegend(0.6, 0.7, 0.88, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.04);

        bool first = true;

        for (size_t m = 0; m < plotters.size(); ++m) {
          auto rdf_cutdata = plotters[m]->GetRDF_Pi0Data().Filter(cutExpr, cutLabel);
          auto rdf_cut = plotters[m]->GetRDF_DVCSPi0MC().Filter(cutExpr, cutLabel);
          if (!rdf_cut.HasColumn(var)) continue;

          auto h = rdf_cut.Histo1D({Form("h_%s_%s_%zu", var.c_str(), cleanName.c_str(), m), (title + ";" + xlabel + ";Counts").c_str(), 100, xmin, xmax}, var);
          h.GetValue();

          TH1D* h_clone = (TH1D*)h.GetPtr()->Clone();
          h_clone->SetDirectory(0);
          NormalizeHistogram(h_clone);

          styleKin_.StyleTH1(h_clone);
          h_clone->SetLineColor(m + 2);
          h_clone->SetLineWidth(2);

          double mean = h_clone->GetMean();
          double sigma = h_clone->GetStdDev();
          double x1 = mean - 3 * sigma;
          double x2 = mean + 3 * sigma;

          TLine* line1 = new TLine(x1, 0, x1, h_clone->GetMaximum() * 0.5);
          TLine* line2 = new TLine(x2, 0, x2, h_clone->GetMaximum() * 0.5);
          line1->SetLineColor(m + 2);
          line2->SetLineColor(m + 2);
          line1->SetLineStyle(2);  // Dashed
          line2->SetLineStyle(2);

          if (first) {
            h_clone->Draw("HIST");
            first = false;
          } else {
            h_clone->Draw("HIST SAME");
          }

          legend->AddEntry(h_clone, labels[m].c_str(), "l");
          std::ostringstream stats;
          stats << "#mu = " << std::fixed << std::setprecision(2) << mean << ", #sigma = " << std::fixed << std::setprecision(2) << sigma;
          legend->AddEntry((TObject*)0, stats.str().c_str(), "");
          line1->Draw("SAME");
          line2->Draw("SAME");
        }

        legend->Draw();
      }

      std::string outpath = outputDir + "/Pi0Exclusivity_" + cleanName + ".pdf";
      canvas->SaveAs(outpath.c_str());
      std::cout << "Saved detector-specific comparison to: " << outpath << "\n";
      delete canvas;
    }
  };

  // phi analysis
  /// For exclusivity cuts, you can use the following function to select one triplet
  /*void PlotPhiAnaExclusivityComparisonByDetectorCases(const std::vector<std::pair<std::string, std::string>>& detectorCuts) {
    std::vector<std::tuple<std::string, std::string, std::string, double, double>> vars = {
        {"Mx2_ep", "Missing Mass Squared (ep)", "MM^{2}(ep) [GeV^{2}]", 0.8, 1.3},
        {"Mx2_epKpKm", "Missing Mass Squared (epK^{+}K^{-})", "MM^{2}(epK^{+}K^{-}) [GeV^{2}]", -0.07, 0.07},
        {"Mx2_eKpKm", "Invariant Mass (eK^{+}K^{-})", "M^{2}(eK^{+}K^{-}) [GeV^{2}]", -0.5, 3},
        {"Mx2_epKp", "Missing Mass Squared (epK^{+})", "MM^{2}(epK^{+}) [GeV^{2}]", -0.5, 1.5},
        {"Mx2_epKm", "Missing Mass Squared (epK^{-})", "MM^{2}(epK^{-}) [GeV^{2}]", -0.5, 1.5},
        {"Emiss", "Missing Energy", "E_{miss} [GeV]", -1.0, 2.0},
        {"PTmiss", "Transverse Missing Momentum", "P_{T}^{miss} [GeV/c]", -0.1, 0.5},
        {"DeltaPhi", "Coplanarity Angle", "#Delta#phi [deg]", 0, 20},
        {"Theta_e_phimeson", "Angle: e-#phi", "#theta(e, #phi) [deg]", 0.0, 60.0}};

    for (const auto& [cutExpr, cutLabel] : detectorCuts) {
      std::string cleanName = cutLabel;
      std::replace(cleanName.begin(), cleanName.end(), ' ', '_');
      std::replace(cleanName.begin(), cleanName.end(), ',', '_');

      TCanvas* canvas = new TCanvas(("c_" + cleanName).c_str(), cutLabel.c_str(), 1800, 1200);
      int cols = 3;
      int rows = (vars.size() + cols - 1) / cols;
      canvas->Divide(cols, rows);

      for (size_t i = 0; i < vars.size(); ++i) {
        canvas->cd(i + 1);
        const auto& [var, title, xlabel, xmin, xmax] = vars[i];
        gPad->SetTicks();
        styleKin_.StylePad((TPad*)gPad);

        TLegend* legend = new TLegend(0.6, 0.55, 0.88, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.04);

        bool first = true;

        for (size_t m = 0; m < plotters.size(); ++m) {
          auto rdf_cut = plotters[m]->GetRDF().Filter(cutExpr, cutLabel);
          if (!rdf_cut.HasColumn(var)) continue;

          auto h = rdf_cut.Histo1D({Form("h_%s_%s_%zu", var.c_str(), cleanName.c_str(), m), (title + ";" + xlabel + ";Counts").c_str(), 100, xmin, xmax}, var);
          h.GetValue();

          TH1D* h_clone = (TH1D*)h.GetPtr()->Clone();
          h_clone->SetDirectory(0);
          NormalizeHistogram(h_clone);

          styleKin_.StyleTH1(h_clone);
          h_clone->SetLineColorAlpha(m + 4,0.8);
          auto [cr, cg, cb] = modelShades[m % modelShades.size()];
          const int colorIdx = 4000 + int(m) * 20;
          if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
          h_clone->SetMarkerColor(colorIdx);
          h_clone->SetLineColorAlpha(colorIdx, 0.8);
          h_clone->SetLineWidth(1);

          double mean = h_clone->GetMean();
          double sigma = h_clone->GetStdDev();
          double x1 = mean - 3 * sigma;
          double x2 = mean + 3 * sigma;

          TLine* line1 = new TLine(x1, 0, x1, h_clone->GetMaximum() * 0.5);
          TLine* line2 = new TLine(x2, 0, x2, h_clone->GetMaximum() * 0.5);
          line1->SetLineColorAlpha(colorIdx, 0.8);
          line2->SetLineColorAlpha(colorIdx, 0.8);
          line1->SetLineStyle(2);  // Dashed
          line2->SetLineStyle(2);

          if (first) {
            h_clone->Draw("HIST");
            first = false;
          } else {
            h_clone->Draw("HIST SAME");
          }

          legend->AddEntry(h_clone, labels[m].c_str(), "l");
          std::ostringstream stats;
          stats << "#mu = " << std::fixed << std::setprecision(2) << mean << ", #sigma = " << std::fixed << std::setprecision(2) << sigma;
          legend->AddEntry((TObject*)0, stats.str().c_str(), "");
          line1->Draw("SAME");
          line2->Draw("SAME");
        }

        legend->Draw();
      }

      std::string outpath = outputDir + "/Exclusivity_Phi_Ana" + cleanName + ".pdf";
      canvas->SaveAs(outpath.c_str());
      std::cout << "Saved detector-specific comparison to: " << outpath << "\n";
      delete canvas;
    }
  };*/

  void PlotPhiAnaExclusivityComparisonByDetectorCases(const std::vector<std::pair<std::string, std::string>>& detectorCuts) {
    std::vector<std::tuple<std::string, std::string, std::string, double, double>> vars = {
        {"Mx2_ep", "Missing Mass Squared (ep)", "MM^{2}(ep) [GeV^{2}]", 0.8, 1.3},
        {"Emiss", "Missing Energy", "E_{miss} [GeV]", -1.0, 2.0},
        {"PTmiss", "Transverse Missing Momentum", "P_{T}^{miss} [GeV/c]", -0.1, 0.5},
        {"Mx2_epKpKm", "Missing Mass Squared (epK^{+}K^{-})", "MM^{2}(epK^{+}K^{-}) [GeV^{2}]", -0.07, 0.07},
        {"Mx2_eKpKm", "Invariant Mass Squared (eK^{+}K^{-})", "M^{2}(eK^{+}K^{-}) [GeV^{2}]", -0.5, 3.0},
        {"Mx2_epKm", "Missing Mass Squared (epK^{-})", "MM^{2}(epK^{-}) [GeV^{2}]", -0.5, 1.5},
        {"Mx2_epKp", "Missing Mass Squared (epK^{+})", "MM^{2}(epK^{+}) [GeV^{2}]", -0.5, 1.5},
        {"DeltaPhi", "Coplanarity Angle", "#Delta#phi [deg]", 0.0, 20.0},
        {"Theta_g_phimeson", "Angle: #gamma - #phi", "#theta(#gamma, #phi) [deg]", 0.0, 10.0},
        {"Theta_e_phimeson", "Angle: e - #phi", "#theta(e, #phi) [deg]", 0.0, 60.0},
        {"DeltaE", "Energy Difference", "#DeltaE [GeV]", -1.0, 1.0},
        {"Cone_p", "Cone Angle (p)", "Cone(p) [deg]", 0.0, 20.0},
        {"Cone_Kp", "Cone Angle (K^{+})", "Cone(K^{+}) [deg]", 0.0, 20.0},
        {"Cone_Km", "Cone Angle (K^{-})", "Cone(K^{-}) [deg]", 0.0, 20.0},
        {"Coplanarity_had_normals_deg", "Coplanarity of Hadronic Normals", "Coplanarity_{had} [deg]", 0.0, 20.0},

        // New: Missing mass of eK+
        {"Mx_eKp", "Missing Mass (eK^{+})", "MM(eK^{+}) [GeV]", 0.0, 4.0}};

    for (const auto& [cutExpr, cutLabel] : detectorCuts) {
      std::string cleanName = cutLabel;
      std::replace(cleanName.begin(), cleanName.end(), ' ', '_');
      std::replace(cleanName.begin(), cleanName.end(), ',', '_');

      // 4 columns now
      int cols = 4;
      int rows = (vars.size() + cols - 1) / cols;
      TCanvas* canvas = new TCanvas(("c_" + cleanName).c_str(), cutLabel.c_str(), 2000, 1200);
      canvas->Divide(cols, rows);

      for (size_t i = 0; i < vars.size(); ++i) {
        canvas->cd(i + 1);
        const auto& [var, title, xlabel, xmin, xmax] = vars[i];
        gPad->SetTicks();
        styleKin_.StylePad((TPad*)gPad);

        TLegend* legend = new TLegend(0.6, 0.55, 0.88, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.04);

        bool first = true;

        for (size_t m = 0; m < plotters.size(); ++m) {
          auto rdf_cut = plotters[m]->GetRDF().Filter(cutExpr, cutLabel);
          if (!rdf_cut.HasColumn(var)) continue;

          auto h = rdf_cut.Histo1D({Form("h_%s_%s_%zu", var.c_str(), cleanName.c_str(), m), (title + ";" + xlabel + ";Counts").c_str(), 100, xmin, xmax}, var);
          h.GetValue();

          TH1D* h_clone = (TH1D*)h.GetPtr()->Clone();
          h_clone->SetDirectory(0);
          NormalizeHistogram(h_clone);

          styleKin_.StyleTH1(h_clone);
          h_clone->SetLineColorAlpha(m + 4, 0.8);
          auto [cr, cg, cb] = modelShades[m % modelShades.size()];
          const int colorIdx = 4000 + int(m) * 20;
          if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
          h_clone->SetMarkerColor(colorIdx);
          h_clone->SetLineColorAlpha(colorIdx, 0.8);
          h_clone->SetLineWidth(1);

          double mean = h_clone->GetMean();
          double sigma = h_clone->GetStdDev();
          double x1 = mean - 3 * sigma;
          double x2 = mean + 3 * sigma;

          TLine* line1 = new TLine(x1, 0, x1, h_clone->GetMaximum() * 0.5);
          TLine* line2 = new TLine(x2, 0, x2, h_clone->GetMaximum() * 0.5);
          line1->SetLineColorAlpha(colorIdx, 0.8);
          line2->SetLineColorAlpha(colorIdx, 0.8);
          line1->SetLineStyle(2);  // Dashed
          line2->SetLineStyle(2);

          if (first) {
            h_clone->Draw("HIST");
            first = false;
          } else {
            h_clone->Draw("HIST SAME");
          }

          legend->AddEntry(h_clone, labels[m].c_str(), "l");
          std::ostringstream stats;
          stats << "#mu = " << std::fixed << std::setprecision(2) << mean << ", #sigma = " << std::fixed << std::setprecision(2) << sigma;
          legend->AddEntry((TObject*)0, stats.str().c_str(), "");
          line1->Draw("SAME");
          line2->Draw("SAME");
        }

        legend->Draw();
      }

      std::string outpath = outputDir + "/Exclusivity_Phi_Ana" + cleanName + ".pdf";
      canvas->SaveAs(outpath.c_str());
      std::cout << "Saved detector-specific comparison to: " << outpath << "\n";
      delete canvas;
    }
  }

  void PlotDIS_BSA_Cross_Section_AndCorr_Comparison(double pol = 1.0, bool plotBSA = true, bool plotDVCSCross = false, bool plotPi0Corr = false, bool plotAccCorr = false,
                                                    bool plotEffCorr = false, bool plotRadCorr = false, bool plotP1Cut = false, bool meanKinVar = false) {
    if (plotters.empty()) {
      std::cerr << "No models loaded to compare.\n";
      return;
    }
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allBSA;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allDVCSCross;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allDVCSPolCross_postive;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allDVCSPolCross_negative;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allPi0Corr;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allPi0DVCSdiffmc;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allPi0DVCSdiffexp;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allAccCorr;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allEffCorr;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allRadCorr;
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> allP1Cut;
    // job for chatgpt
    std::vector<std::vector<std::vector<std::vector<std::tuple<double, double, double>>>>> allBSAmeans;

    for (auto& p : plotters) {
      if (plotBSA) {
        auto h = p->ComputeBSA(fXbins, pol);
        allBSA.push_back(std::move(h));
      }
      if (plotDVCSCross) {
        auto hists = p->ComputeDVCS_CrossSection(fXbins);
        allDVCSCross.push_back(std::move(hists));
        auto hists_pol_postive = p->ComputePolDVCS_CrossSection(1,fXbins);
        auto hists_pol_negative = p->ComputePolDVCS_CrossSection(-1,fXbins);
        allDVCSPolCross_postive.push_back(std::move(hists_pol_postive));
        allDVCSPolCross_negative.push_back(std::move(hists_pol_negative));
      }
      if (plotPi0Corr) {
        auto hcorr = p->ComputePi0Corr(fXbins);
        auto hpi0dvcsdiffmc = p->ComputePi0DVCSdiffmc(fXbins);
        auto hpi0dvcsdiffexp = p->ComputePi0DVCSdiffexp(fXbins);
        allPi0Corr.push_back(std::move(hcorr));
        allPi0DVCSdiffmc.push_back(std::move(hpi0dvcsdiffmc));
        allPi0DVCSdiffexp.push_back(std::move(hpi0dvcsdiffexp));
      }
      if (plotAccCorr) {
        auto hacc = p->ComputeAccCorr(fXbins);
        allAccCorr.push_back(std::move(hacc));
      }
      if (plotEffCorr) {
        auto heff = p->ComputeEffCorr(fXbins);
        allEffCorr.push_back(std::move(heff));
      }
      if (plotRadCorr) {
        auto hrad = p->ComputeRadCorr(fXbins);
        allRadCorr.push_back(std::move(hrad));
      }
      if (plotP1Cut) {
        auto hp1 = p->ComputeP1CutEffect(fXbins);
        allP1Cut.push_back(std::move(hp1));
      }
      if (meanKinVar) {
        allBSAmeans.push_back(getMeanQ2xBt(fXbins, p));
      }
    }

    if (plotBSA) MakeTiledGridComparison("DIS_BSA", "A_{LU}", allBSA, &allBSAmeans, -0.65, 0.65, "pdf", true, true, false, false, meanKinVar);
    if (plotDVCSCross){
        MakeTiledGridComparison("DIS_Cross_Section", "d#sigma/d#phi [nb/GeV^4]", allDVCSCross, &allBSAmeans, 0.0001, 1, "pdf", false, false, true, true, meanKinVar);
        MakeTiledGridComparison("DIS_PolCross_Section_Positive", "d#sigma^{+}/d#phi [nb/GeV^4]", allDVCSPolCross_postive, &allBSAmeans, 0.0001, 1, "pdf", false, false, true, true, meanKinVar);
        MakeTiledGridComparison("DIS_PolCross_Section_Negative", "d#sigma^{-}/d#phi [nb/GeV^4]", allDVCSPolCross_negative, &allBSAmeans, 0.0001, 1, "pdf", false, false, true, true, meanKinVar);
      }
    if (plotPi0Corr) MakeTiledGridComparison("DIS_pi0Corr", "#eta^{#pi^{0}}", allPi0Corr, &allBSAmeans, 0.0, 1, "pdf", false, true, false, false, meanKinVar);
    if (plotPi0Corr) MakeTiledGridComparison("DIS_pi0DVCSdiffmc", "d_{mc}", allPi0DVCSdiffmc, &allBSAmeans, 0.0, 2, "pdf", false, true, false, false, meanKinVar);
    if (plotPi0Corr) MakeTiledGridComparison("DIS_pi0DVCSdiffexp", "d_{exp}", allPi0DVCSdiffexp, &allBSAmeans, 0.0, 2, "pdf", false, true, false, false, meanKinVar);
    if (plotAccCorr) MakeTiledGridComparison("DIS_accCorr", "A_{acc}", allAccCorr, &allBSAmeans, 0.01, 1.0, "pdf", false, true, true, false, meanKinVar);
    if (plotEffCorr) MakeTiledGridComparison("DIS_effCorr", "A_{eff}", allEffCorr, &allBSAmeans, 0.1, 1.1, "pdf", false, true, false, false, meanKinVar);
    if (plotRadCorr) MakeTiledGridComparison("DIS_radCorr", "C_{rad}", allRadCorr, &allBSAmeans, 0.5, 1.5, "pdf", false, true, false, false, meanKinVar);
    if (plotP1Cut) MakeTiledGridComparison("DIS_P1Cut", "C_{P1}", allP1Cut, &allBSAmeans, 0.0, 1.2, "pdf", false, true, false, false, meanKinVar);
  }

  bool file_exists(const char* name) {
    struct stat buffer;
    return (stat(name, &buffer) == 0);
  }

  void dumpHistogram(TH1D* h,
                     double xB, double Q2, double t,
                     double xBmin, double xBmax,
                     double Q2min, double Q2max,
                     double tmin, double tmax,
                     const std::string& label,
                     const char* filename = "h_data.txt") {

    bool exists = file_exists(filename);

    std::ofstream fout(filename, std::ios::out | std::ios::app);
    if (!fout.is_open()) {
      std::cerr << "cannot open " << filename << " to write!\n";
      return;
    }

    if (!exists) {
      fout << "# xB\tQ2\t-t\tphi\tvalue\terror\t"
           << "xBmin\txBmax\tQ2min\tQ2max\ttmin\ttmax\tlabel\n";
    }

    for (int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
      double phi   = h->GetBinCenter(ibin);
      double value = h->GetBinContent(ibin);
      double err   = h->GetBinError(ibin);

      fout << xB << "\t"
           << Q2 << "\t"
           << t << "\t"
           << phi << "\t"
           << value << "\t"
           << err << "\t"
           << xBmin << "\t"
           << xBmax << "\t"
           << Q2min << "\t"
           << Q2max << "\t"
           << tmin << "\t"
           << tmax << "\t"
           << label << "\n";
    }

    fout.close();
    std::cout << "Data " << h->GetName()
              << " written into " << filename
              << (exists ? " (appended)" : "") << "\n";
  }

  void MakeTiledGridComparison(const std::string& observableName, const std::string& yAxisTitle, const std::vector<std::vector<std::vector<std::vector<TH1D*>>>>& histograms,
                               const std::vector<std::vector<std::vector<std::vector<std::tuple<double, double, double>>>>>* meanValues, double yMin, double yMax,
                               const std::string& suffix = "png", bool fitSinusoid = false, bool setManualYrange = false, bool setLogY = false, bool fixlineYrang = true,
                               bool showMeanKin = false) {
    if (histograms.empty() || histograms[0].empty() || histograms[0][0].empty() || histograms[0][0][0].empty()) {
      std::cerr << "No histograms to compare.\n";
      return;
    }
    const auto& q2_edges = fXbins.GetQ2Bins();
    const auto& t_edges = fXbins.GetTBins();
    const auto& xb_edges = fXbins.GetXBBins();

    const size_t n_q2 = q2_edges.size() - 1;
    const size_t n_t = t_edges.size() - 1;
    const size_t n_xb = xb_edges.size() - 1;

    const int rows = n_q2;
    const int cols = n_xb;

    bool Doplot = true;
    int first_perbin_xb = 0;
    bool first_perbin_q2 = true;
    bool first_first_perbin_q2 = true;
    double fixlineYMin = 0.00001;
    double fixlineYMax = 1.0;

    for (size_t t_bin = 0; t_bin < n_t; ++t_bin) {
      TString cname = Form("DIS_BSA_t[%zu]", t_bin);
      TCanvas* c = new TCanvas(cname, cname, 2200, 1600);
      first_perbin_xb = 0;
      first_first_perbin_q2 = true;
      double canvasBorderX = 0.06;
      double canvasBorderY = 0.08;
      double gpad_margin_ratio = 0.2;

      double cellW = (1 - 2 * canvasBorderX) / cols, cellH = (1 - 2 * canvasBorderY) / rows;

      for (size_t q2_bin = 0; q2_bin < n_q2; ++q2_bin) {
        first_perbin_q2 = true;

        /// xbin loop
        for (size_t xb_bin = 0; xb_bin < n_xb; ++xb_bin) {
          int visualRow = rows - 1 - q2_bin;
          int pad = visualRow * cols + xb_bin + 1;
          c->cd();

          bool first = true;
          gStyle->SetCanvasPreferGL(true);

          TLegend* leg = new TLegend(0.35, 0.85, 0.85, 0.95);
          leg->SetBorderSize(0);
          leg->SetFillStyle(0);
          leg->SetTextSize(0.08);

          TLegend* legParams = new TLegend(0.35, 0.16, 0.85, 0.32);  // Bottom legend for a₁
          legParams->SetBorderSize(0);
          legParams->SetFillStyle(0);
          legParams->SetTextSize(0.08);

          TPad* thisPad = new TPad(Form("%zu_%zu", xb_bin, q2_bin), Form("%zu_%zu", xb_bin, q2_bin), cellW * xb_bin + canvasBorderX, cellH * (q2_bin) + canvasBorderY,
                                   cellW * (xb_bin + 1) + canvasBorderX, cellH * (q2_bin + 1) + canvasBorderY);
          double l = 0.00, r = 0.00, b = 0.00, t = 0.00;
          Doplot = false;

          for (size_t m = 0; m < histograms.size(); ++m) {
            // if(q2_bin == 2 && xb_bin == 2) continue; // save first per bin xb
            //  Pad margins
            /*
            double l = (first_perbin_q2) ? 0.2 : 0.00;
            double r = (!first_perbin_q2) ? 0.00 : 0.00;
            double b = (xb_bin == first_perbin_xb) ? 0.16 : 0.00;
            double t = (visualRow == 0) ? 0.000 : 0.00;
            */

            // const int idx = q2_bin * (n_t * n_xb) + t_bin * n_xb + xb_bin;

            TH1D* h = histograms[m][xb_bin][q2_bin][t_bin];

            styleBSA_.StyleTH1(h);

            auto [r, g, b] = modelShades[m % modelShades.size()];

            int colorIdx = 3000 + m * 20;  // Avoid low TColor indices

            if (!gROOT->GetColor(colorIdx)) {
              new TColor(colorIdx, r, g, b);
            }

            h->SetLineColor(colorIdx);
            h->SetMarkerColor(colorIdx);
            h->SetFillColorAlpha(colorIdx, 1.0);

            h->SetLineColor(colorIdx);
            h->SetMarkerColor(colorIdx);
            h->SetFillColorAlpha(colorIdx, 1.0);
            h->SetLineWidth(1);
            h->SetMarkerStyle(20);
            h->SetMarkerSize(1.0);
            h->SetStats(0);

            if (first) {
              l = (first_perbin_q2) ? (gpad_margin_ratio) / (1 + gpad_margin_ratio) : 0.00;
              r = (xb_bin == first_perbin_xb) ? (gpad_margin_ratio) / (1 + gpad_margin_ratio) : 0.00;
              b = (xb_bin == first_perbin_xb) ? (gpad_margin_ratio) / (1 + gpad_margin_ratio) : 0.00;
              t = (visualRow == 0) ? 0.000 : 0.00;

              l = (first_perbin_q2 && xb_bin == first_perbin_xb) ? (gpad_margin_ratio) / (1 + 2 * gpad_margin_ratio) : l;
              r = (xb_bin == first_perbin_xb && first_perbin_q2) ? (gpad_margin_ratio) / (1 + 2 * gpad_margin_ratio) : r;

              styleBSA_.StylePad(thisPad, l, r, b, t);
              thisPad->SetTicks(1, 0);
              thisPad->SetFillStyle(4000);
              // std::cout << "xb_bin: " << xb_bin << ", q2_bin: " << q2_bin << ", first_perbin_xb: " << first_perbin_xb << ", first_perbin_q2: " << first_perbin_q2 << ",
              // first_first_perbin_q2: " <<first_first_perbin_q2<< std::endl;

              h->GetXaxis()->SetTitle((xb_bin == first_perbin_xb) ? "#phi [deg]" : "");
              h->GetYaxis()->SetTitle((first_perbin_q2) ? yAxisTitle.c_str() : "");
              h->GetXaxis()->SetLabelSize((xb_bin == first_perbin_xb) ? 0.085 : 0.0);
              h->GetXaxis()->SetTitleSize((xb_bin == first_perbin_xb) ? 0.095 : 0.0);
              h->GetYaxis()->SetLabelSize((first_perbin_q2) ? 0.085 : 0.0);
              h->GetYaxis()->SetTitleSize((first_perbin_q2) ? 0.1 : 0.0);
              if (xb_bin == first_perbin_xb && first_perbin_q2) {
                h->GetYaxis()->SetLabelSize(0.085 * (1 + gpad_margin_ratio) / (1 + 2 * gpad_margin_ratio));
                h->GetYaxis()->SetTitleSize(0.1 * (1 + gpad_margin_ratio) / (1 + 2 * gpad_margin_ratio));
              }
            }

            h->GetXaxis()->SetTitleOffset((xb_bin == first_perbin_xb) ? 0.82 : 0.0);
            h->GetYaxis()->SetTitleOffset((first_perbin_q2) ? 0.82 : 0.0);

            h->GetXaxis()->SetNdivisions(4, false);
            h->GetYaxis()->SetNdivisions(6, true);
            if (setManualYrange) h->GetYaxis()->SetRangeUser(yMin, yMax);
            if (fixlineYrang && first_perbin_q2 && h->GetMinimum() > 0) {
              fixlineYMin = h->GetMinimum() * 0.3;
              fixlineYMax = h->GetMaximum() * 3;
              h->GetYaxis()->SetRangeUser(fixlineYMin, fixlineYMax);
            }
            if (fixlineYrang && !first_perbin_q2) {
              h->GetYaxis()->SetRangeUser(fixlineYMin, fixlineYMax);
            }
            h->GetXaxis()->SetRangeUser(0, 360);

            h->GetXaxis()->CenterTitle(true);
            h->GetYaxis()->CenterTitle(true);
            Doplot = !(!h || h->GetBinContent(5) == 0) || Doplot;
            if (!h || h->GetBinContent(5) == 0) {
              continue;
            }

            if (first) {
              if (xb_bin == first_perbin_xb && first_perbin_q2) {
                thisPad->SetPad(cellW * (xb_bin - gpad_margin_ratio) + canvasBorderX, cellH * (q2_bin - gpad_margin_ratio) + canvasBorderY,
                                cellW * (xb_bin + 1 + gpad_margin_ratio) + canvasBorderX, cellH * (q2_bin + 1) + canvasBorderY);
              } else if (xb_bin == first_perbin_xb && !first_perbin_q2) {
                h->GetXaxis()->ChangeLabel(1, -1, 0, -1, -1, -1, "");  // blank it out
                thisPad->SetPad(cellW * (xb_bin) + canvasBorderX, cellH * (q2_bin - gpad_margin_ratio) + canvasBorderY, cellW * (xb_bin + 1 + gpad_margin_ratio) + canvasBorderX,
                                cellH * (q2_bin + 1) + canvasBorderY);
              } else if (first_perbin_q2 && xb_bin != first_perbin_xb) {
                thisPad->SetPad(cellW * (xb_bin - gpad_margin_ratio) + canvasBorderX, cellH * (q2_bin) + canvasBorderY, cellW * (xb_bin + 1) + canvasBorderX,
                                cellH * (q2_bin + 1) + canvasBorderY);
              } else if (!first_perbin_q2 && xb_bin != first_perbin_xb) {
                thisPad->SetPad(cellW * (xb_bin) + canvasBorderX, cellH * (q2_bin) + canvasBorderY, cellW * (xb_bin + 1) + canvasBorderX, cellH * (q2_bin + 1) + canvasBorderY);
              }
            }

            if (first) thisPad->Draw();
            thisPad->cd();
            if (setLogY) thisPad->SetLogy();
            h->Draw(first ? "E1X0" : "E1X0 SAME");
            first = false;
            first_perbin_q2 = false;

            // Fit function and extract a₁
            if (fitSinusoid) {
              TF1* fitFunc = new TF1(Form("fit_%zu_%zu_%zu_%zu", m, t_bin, q2_bin, xb_bin), "[0] + ([1]*sin(x*TMath::DegToRad())) / (1 + [2]*cos(x*TMath::DegToRad()))", 0, 360);
              fitFunc->SetParameters(0.0, 0.2, 0.1);
              fitFunc->SetFillColorAlpha(colorIdx, 0.5);
              fitFunc->SetLineColorAlpha(colorIdx, 0.5);
              fitFunc->SetLineStyle(2);
              fitFunc->SetLineWidth(1);
              h->Fit(fitFunc, "Q0");
              fitFunc->Draw("SAME");

              double a1 = fitFunc->GetParameter(1);
              double a1e = fitFunc->GetParError(1);
              TString a1label = Form("a_{1} = %.2f #pm %.2f", a1, a1e);
              legParams->AddEntry(fitFunc, a1label, "l");
            }
            leg->AddEntry(h, labels[m].c_str(), "p");
            // auto [mean_xB, mean_Q2, mean_t] = meanValues[m][xb_bin][q2_bin][t_bin];
            if (showMeanKin) {
              auto [mean_xB, mean_Q2, mean_t] = (*meanValues)[m][xb_bin][q2_bin][t_bin];
              TString meanText = Form("<x_{B}> = %.2f, <Q^{2}> = %.2f, <|t|> = %.2f", mean_xB, mean_Q2, mean_t);
              TLatex* meanLatex = new TLatex(0.25, 0.78 - m * 0.10, meanText.Data());
              meanLatex->SetTextSize(0.05);
              meanLatex->SetNDC();
              meanLatex->SetTextFont(42);
              meanLatex->Draw();
              dumpHistogram(h, mean_xB, mean_Q2, mean_t, xb_edges[xb_bin], xb_edges[xb_bin + 1], q2_edges[q2_bin], q2_edges[q2_bin + 1], t_edges[t_bin], t_edges[t_bin + 1], labels[m],
                            Form("datapoint_%s.txt", observableName.c_str()));
            }
          }
          if (!Doplot) {
            if (first_first_perbin_q2) first_perbin_xb++;
            std::cout << "No data for this bin combination, skipping...\n";
            continue;
          }
          /*
                    // Annotate bin ranges
                    double xB_low = xb_edges[xb_bin], xB_high = xb_edges[xb_bin + 1];
                    double Q2_low = q2_edges[q2_bin], Q2_high = q2_edges[q2_bin + 1];
                    TString labelText = Form("x_{B} #in [%.2f, %.2f], Q^{2} #in [%.1f, %.1f]", xB_low, xB_high, Q2_low, Q2_high);
                    TLatex* latex = new TLatex(0.25, 0.82, labelText.Data());
                    latex->SetTextSize(0.055);
                    latex->SetNDC();
                    latex->SetTextFont(42);
                    latex->Draw();
          */
          leg->Draw();
          if (fitSinusoid) legParams->Draw();

          thisPad->Modified();
          thisPad->Update();
          c->Modified();
          c->Update();

          if (first_first_perbin_q2 && xb_bin == first_perbin_xb) {
            first_perbin_xb++;
            first_first_perbin_q2 = false;
            // std::cout<< "1"<<std::endl;
          } else if (!first_first_perbin_q2 && xb_bin == first_perbin_xb) {
            first_perbin_xb++;
            // std::cout<< "2"<<std::endl;
          } else if (first_first_perbin_q2 && xb_bin != first_perbin_xb) {
            first_perbin_xb++;
            // std::cout<< "3"<<std::endl;
          } else if (!first_first_perbin_q2 && xb_bin != first_perbin_xb) {
            // std::cout<< "4"<<std::endl;
          }
        }
      }
      TString outfile = Form("%s/%s_t_%.2f-%.2f.%s", outputDir.c_str(), observableName.c_str(), t_edges[t_bin], t_edges[t_bin + 1], suffix.c_str());
      c->SaveAs(outfile);

      // std::cout << "Saved: " << outfile << '\n';

      delete c;
    }
  }
  // === Add to DISANAcomparer (public): =========================
  void PlotPhiInvMassPerBin_AllModels(const std::string& baseOutDir = "PhiInvMassFits", int nBins = 120, double mMin = 0.98, double mMax = 1.08, bool constrainSigma = true,
                                      double sigmaRef = 0.004, double sigmaFrac = 0.25, double branching = 1.0, bool doAcceptanceCorr = false, bool doEfficiencyCorr = false,
                                      bool doRadCorr = false) {
    if (plotters.empty()) {
      std::cerr << "[PlotPhiInvMassPerBin_AllModels] no models.\n";
      return;
    }

    gSystem->Exec(Form("mkdir -p %s", baseOutDir.c_str()));
    for (size_t i = 0; i < plotters.size(); ++i) {
      DISANAplotter* p = plotters[i].get();
      if (!p) continue;

      const std::string subdir = baseOutDir + "/" + labels[i];
      gSystem->Exec(Form("mkdir -p %s", subdir.c_str()));

      std::cout << "→ Fitting/drawing per-bin K^{+}K^{-} mass for model: " << labels[i] << " → " << subdir << "  (AccCorr=" << (doAcceptanceCorr ? "ON" : "OFF")
                << ", RadCorr=" << (doRadCorr ? "ON" : "OFF") << ")\n";

      if (doAcceptanceCorr) {
        p->MakePhiAcceptanceCorrection3D(fXbins, subdir + "/acc/", "mtprime", "mtprime", /*minGenEntries=*/10, constrainSigma, sigmaRef, sigmaFrac);
      }

      if (doRadCorr) {
        p->MakePhiRadiativeCorrection3D_FromRatio(fXbins, subdir + "/rad/", "rad_corr", "mtprime");
      }
      if (doEfficiencyCorr) {
        p->MakePhiEfficiencyCorrection3D(fXbins, subdir + "/eff/");
      }

      p->SetPlotApplyAcceptanceCorrection(doAcceptanceCorr);
      p->SetPlotApplyRadiativeCorrection(doRadCorr);

      (void)p->MakePhiMassFitCanvases3D(fXbins, subdir, nBins, mMin, mMax, constrainSigma, sigmaRef, sigmaFrac, branching);
    }
  }

  // =========================================================================
  // PlotPhiInvMassPerBin_XBBins_AllModels
  // =========================================================================
  // Drives MakePhiMassFitCanvases_XBBins() for every loaded model.
  // Bins dσ/dt' in (x_B, t') — the correct binning for gluon radius
  // extraction because at fixed x_B the value of W (and hence |t_min|)
  // is kinematically determined event-by-event instead of smeared across
  // a wide W range.
  // =========================================================================
  void PlotPhiInvMassPerBin_XBBins_AllModels(const std::string& baseOutDir = "PhiInvMassFits_xBBins", int nBins = 120, double mMin = 0.98, double mMax = 1.08,
                                             bool constrainSigma = true, double sigmaRef = 0.004, double sigmaFrac = 0.25, double branching = 1.0, bool doAcceptanceCorr = false,
                                             bool doEfficiencyCorr = false, bool doRadCorr = false) {
    if (plotters.empty()) {
      std::cerr << "[PlotPhiInvMassPerBin_XBBins_AllModels] no models.\n";
      return;
    }

    // Require x_B binning to be set
    const auto& xbEdges = fXbins.GetXBBins();
    if (xbEdges.size() < 2) {
      std::cerr << "[PlotPhiInvMassPerBin_XBBins_AllModels] "
                   "No x_B bins found in BinManager. "
                   "Call SetXBinsRanges() with a BinManager that has GetXBBins() set.\n";
      return;
    }

    gSystem->Exec(Form("mkdir -p %s", baseOutDir.c_str()));

    for (size_t i = 0; i < plotters.size(); ++i) {
      DISANAplotter* p = plotters[i].get();
      if (!p) continue;

      const std::string subdir = baseOutDir + "/" + labels[i];
      gSystem->Exec(Form("mkdir -p %s", subdir.c_str()));

      std::cout << "→ [xB-binned] Fitting K^{+}K^{-} mass for model: " << labels[i] << " → " << subdir << "  (AccCorr=" << (doAcceptanceCorr ? "ON" : "OFF")
                << ", RadCorr=" << (doRadCorr ? "ON" : "OFF") << ")\n";

      p->SetPlotApplyAcceptanceCorrection(doAcceptanceCorr);
      p->SetPlotApplyRadiativeCorrection(doRadCorr);

      p->MakePhiMassFitCanvases_XBBins(fXbins, subdir, nBins, mMin, mMax, constrainSigma, sigmaRef, sigmaFrac, branching);
    }
  }

  // =========================================================================
  // PlotPhiDSigmaDt_XBBins_FromCache
  // =========================================================================
  // Plots and exports dσ/dt' vs t' in each x_B bin for all models.
  // The CSV written per (model, x_B bin) has columns designed so that the
  // Python script can compute a physically correct |t_min| from the
  // data-weighted mean W (stored as W_mean) rather than from bin-boundary
  // arithmetic. The key columns written are:
  //
  //   tprime_center, tprime_lo, tprime_hi
  //   xB_lo, xB_hi
  //   W_mean        ← data-weighted mean W  (use this for tmin)
  //   Q2_mean       ← data-weighted mean Q²
  //   xB_mean       ← data-weighted mean x_B
  //   tmin_mean     ← mean |t_min|(Q²_i, W_i) averaged event-by-event
  //   GammaV_mean
  //   CrossSection, CrossSection_Err
  //   RawCounts, RawCounts_Err
  // =========================================================================
  void PlotPhiDSigmaDt_XBBins_FromCache(bool writeCSV = true, bool plotData = true) {
    if (plotters.empty()) return;

    const auto& xbEdges = fXbins.GetXBBins();
    const auto& tpEdges = fXbins.GetTprimeBins();
    if (xbEdges.size() < 2 || tpEdges.size() < 2) {
      std::cerr << "[PlotPhiDSigmaDt_XBBins_FromCache] "
                   "Missing x_B or t' bin edges in BinManager.\n";
      return;
    }

    const size_t nXB = xbEdges.size() - 1;
    const size_t nT = tpEdges.size() - 1;

    // ---- CSV writer --------------------------------------------------------
    auto writeCSVForXBBin = [&](size_t im, size_t ixB) {
      TH1D* hXS = nullptr;
      TH1D* hAcc = nullptr;
      TH1D* hEff = nullptr;
      TH1D* hRad = nullptr;
      TH1D* hNsig = nullptr;

      {
        auto& H = plotters[im]->GetPhiDSigmaDt_XBT();
        if (ixB < H.size()) hXS = H[ixB];
      }
      {
        auto& A = plotters[im]->GetPhiAcceptance_XBT();
        if (ixB < A.size()) hAcc = A[ixB];
      }
      {
        auto& E = plotters[im]->GetPhiEfficiency_XBT();
        if (ixB < E.size()) hEff = E[ixB];
      }
      {
        auto& R = plotters[im]->GetPhiRadCorr_XBT();
        if (ixB < R.size()) hRad = R[ixB];
      }
      {
        auto& N = plotters[im]->GetPhiRawCounts_XBT();
        if (ixB < N.size()) hNsig = N[ixB];
      }
      if (!hXS) return;

      // Sanitise label for directory name
      std::string safeLabel = labels[im];
      for (char& ch : safeLabel)
        if (ch == ' ' || ch == '/' || ch == '\\') ch = '_';

      const std::string csvDir = outputDir + "/CSVs_xBBins/" + safeLabel;
      gSystem->Exec(Form("mkdir -p \"%s\"", csvDir.c_str()));

      const std::string fname = Form("%s/dsdt_xB%zu.csv", csvDir.c_str(), ixB);

      std::ofstream csv(fname);
      if (!csv.is_open()) {
        std::cerr << "[writeCSV_xB] Cannot open " << fname << "\n";
        return;
      }

      // Header — all quantities needed for tmin and xB extraction downstream
      csv << "tprime_center,tprime_lo,tprime_hi"
          << ",xB_lo,xB_hi"
          << ",xB_mean"    // data-weighted mean x_B in this (xB, t') cell
          << ",W_mean"     // data-weighted mean W  — USE THIS for tmin/xB
          << ",Q2_mean"    // data-weighted mean Q²
          << ",tmin_mean"  // mean |t_min|(Q²_i,W_i) averaged event-by-event
          << ",GammaV_mean"
          << ",CrossSection,CrossSection_Err"
          << ",Acceptance,Efficiency,RadCorr"
          << ",RawCounts,RawCounts_Err\n";
      csv << std::setprecision(6) << std::scientific;

      const double xbLo = xbEdges[ixB];
      const double xbHi = xbEdges[ixB + 1];

      for (int ib = 1; ib <= hXS->GetNbinsX(); ++ib) {
        const size_t it = static_cast<size_t>(ib - 1);
        const double tp_cen = hXS->GetBinCenter(ib);
        const double tp_lo = hXS->GetBinLowEdge(ib);
        const double tp_hi = hXS->GetBinLowEdge(ib + 1);

        const double xs = hXS->GetBinContent(ib);
        const double xs_err = hXS->GetBinError(ib);
        const double acc = hAcc ? hAcc->GetBinContent(ib) : std::numeric_limits<double>::quiet_NaN();
        const double eff = hEff ? hEff->GetBinContent(ib) : std::numeric_limits<double>::quiet_NaN();
        const double rad = hRad ? hRad->GetBinContent(ib) : std::numeric_limits<double>::quiet_NaN();
        const double nsig = hNsig ? hNsig->GetBinContent(ib) : std::numeric_limits<double>::quiet_NaN();
        const double nsig_err = hNsig ? hNsig->GetBinError(ib) : std::numeric_limits<double>::quiet_NaN();

        // Per-(xB,t') data-weighted means from new getters
        const double xb_mean = plotters[im]->GetPhiMeanXB_XBT(ixB, it);
        const double w_mean = plotters[im]->GetPhiMeanW_XBT(ixB, it);
        const double q2_mean = plotters[im]->GetPhiMeanQ2_XBT(ixB, it);
        const double tmin_mean = plotters[im]->GetPhiMeanTmin_XBT(ixB, it);
        const double gv_mean = plotters[im]->GetPhiMeanGammaV_XBT(ixB, it);

        csv << tp_cen << "," << tp_lo << "," << tp_hi << "," << xbLo << "," << xbHi << "," << xb_mean << "," << w_mean << "," << q2_mean << "," << tmin_mean << "," << gv_mean
            << "," << xs << "," << xs_err << "," << acc << "," << eff << "," << rad << "," << nsig << "," << nsig_err << "\n";
      }
      csv.close();
      std::cout << "[CSV-xB] Written → " << fname << "\n";
    };

    // ---- Plotting ----------------------------------------------------------
    auto doOnePlot = [&](size_t ixB) {
      const double xbLo = xbEdges[ixB];
      const double xbHi = xbEdges[ixB + 1];
      const TString head = Form("x_{B} #in [%.3f, %.3f]", xbLo, xbHi);

      TCanvas* c = new TCanvas(Form("c_phi_dsdt_xB%zu", ixB), "", 1200, 900);
      styleCrossSection_.StylePad((TPad*)gPad);
      gPad->SetTicks(1, 1);
      gPad->SetLogy();

      TLegend* leg = new TLegend(0.60, 0.72, 0.92, 0.90);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.035);

      bool first = true;
      for (size_t im = 0; im < plotters.size(); ++im) {
        auto& H = plotters[im]->GetPhiDSigmaDt_XBT();
        if (ixB >= H.size() || !H[ixB]) continue;
        TH1D* h = H[ixB];

        styleCrossSection_.StyleTH1(h);
        auto [cr, cg, cb] = modelShades[im % modelShades.size()];
        const int colorIdx = 5500 + (int)im * 20;
        if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
        h->SetLineColor(colorIdx);
        h->SetMarkerColor(colorIdx);
        h->SetMarkerStyle(20);
        h->SetMarkerSize(1.0);
        h->SetTitle("");
        h->GetXaxis()->SetTitle("-t' [GeV^{2}]");
        h->GetYaxis()->SetTitle("d#sigma/dt' [nb/GeV^{2}]");

        if (first) {
          h->Draw("E1X0");
          TLatex latex;
          latex.SetNDC();
          latex.SetTextFont(42);
          latex.SetTextSize(0.040);
          latex.DrawLatex(0.14, 0.93, head);

          // Annotate with mean kinematics from the lowest-t' bin (bin 0)
          const double wmean = plotters[im]->GetPhiMeanW_XBT(ixB, 0);
          const double q2mean = plotters[im]->GetPhiMeanQ2_XBT(ixB, 0);
          if (std::isfinite(wmean) && std::isfinite(q2mean)) latex.DrawLatex(0.14, 0.88, Form("<W> = %.2f GeV,  <Q^{2}> = %.2f GeV^{2}", wmean, q2mean));
          first = false;
        } else {
          h->Draw("E1X0 SAME");
        }
        leg->AddEntry(h, labels[im].c_str(), "lep");
      }
      leg->Draw();

      const TString out = Form("%s/phi_dsdt_xBBin_%zu.pdf", outputDir.c_str(), ixB);
      c->SaveAs(out);
      delete leg;
      delete c;
    };

    // ---- Main loop ---------------------------------------------------------
    for (size_t ixB = 0; ixB < nXB; ++ixB) {
      if (plotData) doOnePlot(ixB);
      if (writeCSV) {
        for (size_t im = 0; im < plotters.size(); ++im) writeCSVForXBBin(im, ixB);
      }
    }
  }

  void PlotPhiDSigmaDt_FromCache(bool plotData = true, bool plotAcc = false, bool plotRadCorr = false, bool plotEff = false, bool writeCSV = true) {
    if (plotters.empty()) return;

    const auto& q2 = fXbins.GetQ2Bins();
    const auto& w = fXbins.GetWBins();
    const auto& t_edges = fXbins.GetTBins();  // |t| bin edges (may be empty)
    const bool hasW = !w.empty();
    const bool hasT = (t_edges.size() > 1);

    const size_t nQ = (q2.size() > 1) ? (q2.size() - 1) : 0;
    const size_t nW = hasW ? (w.size() - 1) : 1;

    // -----------------------------------------------------------------------
    //  CSV helper – writes one CSV per (model, Q2-bin, W-bin)
    //  Columns: tprime_center, tprime_lo, tprime_hi,
    //           Q2_center, Q2_lo, Q2_hi,
    //           t_center, t_lo, t_hi,       (only when t-bins exist)
    //           W_center, W_lo, W_hi,        (only when W-bins exist)
    //           CrossSection, Mean, Err,
    //           Acceptance, RadCorr
    // -----------------------------------------------------------------------
    auto writeCSVForBin = [&](size_t im, size_t iq, size_t iw) {
      // --- gather the three histograms (xs mandatory, acc/rad optional) ---
      TH1D* hXS = nullptr;
      TH1D* hAcc = nullptr;
      TH1D* hEff = nullptr;
      TH1D* hRad = nullptr;
      TH1D* hNsig = nullptr;

      {
        auto& H = plotters[im]->GetPhiDSigmaDt3D();
        if (iq < H.size() && iw < H[iq].size()) hXS = H[iq][iw];
      }
      {
        auto& A = plotters[im]->GetPhiAcceptance3D();
        if (iq < A.size() && iw < A[iq].size()) hAcc = A[iq][iw];
      }
      {
        auto& E = plotters[im]->GetPhiEfficiency3D();
        if (iq < E.size() && iw < E[iq].size()) hEff = E[iq][iw];
      }
      {
        auto& R = plotters[im]->GetPhiRadCorr3D();
        if (iq < R.size() && iw < R[iq].size()) hRad = R[iq][iw];
      }
      {
        auto& N = plotters[im]->GetPhiRawCounts3D();
        if (iq < N.size() && iw < N[iq].size()) hNsig = N[iq][iw];
      }

      if (!hXS) return;  // nothing to write

      // --- build output path: CSVs/<ModelName>/ ---
      // Sanitise label for use as directory name (replace spaces / slashes)
      std::string safeLabel = labels[im];
      for (char& ch : safeLabel)
        if (ch == ' ' || ch == '/' || ch == '\\') ch = '_';

      const std::string csvDir = outputDir + "/CSVs/" + safeLabel;
      gSystem->Exec(Form("mkdir -p \"%s\"", csvDir.c_str()));

      // One CSV per (Q2-bin, W-bin)
      std::string fname;
      if (hasW)
        fname = Form("%s/dsdt_Q%zu_W%zu.csv", csvDir.c_str(), iq, iw);
      else
        fname = Form("%s/dsdt_Q%zu.csv", csvDir.c_str(), iq);

      std::ofstream csv(fname);
      if (!csv.is_open()) {
        std::cerr << "[writeCSV] Cannot open " << fname << "\n";
        return;
      }

      // --- convenience bin-centre/edge values for outer axes ---
      const double q2_lo = q2[iq];
      const double q2_hi = q2[iq + 1];
      const double q2_cen = 0.5 * (q2_lo + q2_hi);

      const double w_lo = hasW ? w[iw] : std::numeric_limits<double>::quiet_NaN();
      const double w_hi = hasW ? w[iw + 1] : std::numeric_limits<double>::quiet_NaN();
      // NOTE: w_cen is the arithmetic midpoint of the W bin boundaries (e.g. 5.9 GeV
      // for a 1.8–10 GeV bin).  It is only used for W2 columns; the per-row W_mean
      // below comes from the actual data-weighted mean and is the physically meaningful
      // quantity for tmin and xB computation.
      const double w_cen = hasW ? 0.5 * (w_lo + w_hi) : std::numeric_limits<double>::quiet_NaN();

      // --- write header ---
      csv << "tprime_center,tprime_lo,tprime_hi"
          << ",t_center,t_lo,t_hi"
          << ",Q2_center,Q2_lo,Q2_hi";
      if (hasW)
        csv << ",W2_center,W2_lo,W2_hi"
            << ",W_lo,W_hi"  // bin boundaries kept for reference
            << ",W_center";  // data-weighted mean W (correct for tmin/xB)

      // NEW mean-kin columns (per (Q2,W,t') bin)
      csv << ",xB_mean,W_mean,GammaV_mean";

      csv << ",CrossSection,CrossSection_Err,tprime_mean"
          << ",Acceptance,Efficiency,RadCorr"
          << ",RawCounts,RawCounts_Err\n";
      csv << std::setprecision(6) << std::scientific;

      const int nBins = hXS->GetNbinsX();
      for (int ib = 1; ib <= nBins; ++ib) {
        const double tp_cen = hXS->GetBinCenter(ib);
        const double tp_lo = hXS->GetBinLowEdge(ib);
        const double tp_hi = hXS->GetBinLowEdge(ib + 1);

        // |t| from t-edges array if available, else fall back to tprime bin
        double t_cen, t_lo_val, t_hi_val;
        if (hasT && (ib - 1) < (int)(t_edges.size() - 1)) {
          t_lo_val = t_edges[ib - 1];
          t_hi_val = t_edges[ib];
          t_cen = 0.5 * (t_lo_val + t_hi_val);
        } else {
          t_cen = tp_cen;
          t_lo_val = tp_lo;
          t_hi_val = tp_hi;
        }

        const double xs = hXS->GetBinContent(ib);
        const double err = hXS->GetBinError(ib);
        const double mean_val = tp_cen;

        const double acc = hAcc ? hAcc->GetBinContent(ib) : std::numeric_limits<double>::quiet_NaN();
        const double eff = hEff ? hEff->GetBinContent(ib) : std::numeric_limits<double>::quiet_NaN();
        const double rad = hRad ? hRad->GetBinContent(ib) : std::numeric_limits<double>::quiet_NaN();
        const double nsig = hNsig ? hNsig->GetBinContent(ib) : std::numeric_limits<double>::quiet_NaN();
        const double nsig_err = hNsig ? hNsig->GetBinError(ib) : std::numeric_limits<double>::quiet_NaN();
        const int ib0 = ib - 1;

        const double xBmean = plotters[im]->GetPhiMeanXB(iq, iw, ib0);
        const double Wmean = plotters[im]->GetPhiMeanW(iq, iw, ib0);
        const double Gvmean = plotters[im]->GetPhiMeanGammaV(iq, iw, ib0);

        // W_center: use data-weighted mean W when available (finite and physical),
        // fall back to bin-boundary midpoint only as last resort.
        const double W_center_out = (std::isfinite(Wmean) && Wmean > 0.5) ? Wmean : w_cen;

        csv << tp_cen << "," << tp_lo << "," << tp_hi << "," << t_cen << "," << t_lo_val << "," << t_hi_val << "," << q2_cen << "," << q2_lo << "," << q2_hi;

        if (hasW) csv << "," << w_cen * w_cen << "," << w_lo * w_lo << "," << w_hi * w_hi << "," << w_lo << "," << w_hi << "," << W_center_out;

        // NEW means
        csv << "," << xBmean << "," << Wmean << "," << Gvmean;

        csv << "," << xs << "," << err << "," << mean_val << "," << acc << "," << eff << "," << rad << "," << nsig << "," << nsig_err << "\n";
      }

      csv.close();
      std::cout << "[CSV] Written → " << fname << "\n";
    };

    // -----------------------------------------------------------------------
    //  Main plotting lambda (unchanged from original)
    // -----------------------------------------------------------------------
    auto doOne = [&](const TString& tag, const TString& yTitle, double yMin, double yMax, std::function<TH1D*(DISANAplotter*, size_t, size_t)> getHist) {
      for (size_t iq = 0; iq < nQ; ++iq)
        for (size_t iw = 0; iw < nW; ++iw) {
          TCanvas* c = new TCanvas(Form("c_%s_Q%zu_W%zu", tag.Data(), iq, iw), "", 1200, 900);
          styleCrossSection_.StylePad((TPad*)gPad);
          gPad->SetTicks(1, 1);

          TLegend* leg = new TLegend(0.60, 0.72, 0.92, 0.90);
          leg->SetBorderSize(0);
          leg->SetFillStyle(0);
          leg->SetTextSize(0.035);

          TString head = hasW ? Form("Q^{2}[%.2f, %.2f]   W[%.1f, %.1f]", q2[iq], q2[iq + 1], w[iw], w[iw + 1]) : Form("Q^{2}[%.2f, %.2f]", q2[iq], q2[iq + 1]);

          bool first = true;

          for (size_t im = 0; im < plotters.size(); ++im) {
            TH1D* h = getHist(plotters[im].get(), iq, iw);
            if (!h) continue;

            styleCrossSection_.StyleTH1(h);

            auto [cr, cg, cb] = modelShades[im % modelShades.size()];
            const int colorIdx = 5000 + int(im) * 20;
            if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);

            h->SetLineColor(colorIdx);
            h->SetMarkerColor(colorIdx);
            h->SetMarkerStyle(20);
            h->SetMarkerSize(1.0);

            h->SetTitle("");
            h->GetXaxis()->SetTitle("-t' [GeV^{2}]");
            h->GetYaxis()->SetTitle(yTitle);
            if (yTitle == "d#sigma/dt' [nb/GeV^{2}]") {
              gPad->SetLogy();
            }

            if (first) {
              h->Draw("E1X0");
              TLatex latex;
              latex.SetNDC();
              latex.SetTextFont(42);
              latex.SetTextSize(0.040);
              latex.DrawLatex(0.14, 0.93, head);
              first = false;
            } else {
              h->Draw("E1X0 SAME");
            }

            leg->AddEntry(h, labels[im].c_str(), "lep");
          }

          leg->Draw();

          TString out = hasW ? Form("%s/%s_Q%zu_W%zu.pdf", outputDir.c_str(), tag.Data(), iq, iw) : Form("%s/%s_Q%zu.pdf", outputDir.c_str(), tag.Data(), iq);

          c->SaveAs(out);

          delete leg;
          delete c;
        }
    };

    // ---- Run enabled plot types ----
    if (plotData) {
      doOne("phi_dsdt", "d#sigma/dt' [nb/GeV^{2}]", 0.0, 10.0, [](DISANAplotter* P, size_t iq, size_t iw) -> TH1D* {
        auto& H = P->GetPhiDSigmaDt3D();
        if (iq >= H.size() || iw >= H[iq].size()) return nullptr;
        return H[iq][iw];
      });
    }

    if (plotAcc) {
      doOne("phi_accept", "Acceptance", 0.0, .05, [](DISANAplotter* P, size_t iq, size_t iw) -> TH1D* {
        auto& A = P->GetPhiAcceptance3D();
        if (iq >= A.size() || iw >= A[iq].size()) return nullptr;
        return A[iq][iw];
      });
    }
    if (plotEff) {
      doOne("phi_efficiency", "Efficiency", 0.0, .05, [](DISANAplotter* P, size_t iq, size_t iw) -> TH1D* {
        auto& E = P->GetPhiEfficiency3D();
        if (iq >= E.size() || iw >= E[iq].size()) return nullptr;
        return E[iq][iw];
      });
    }

    if (plotRadCorr) {
      doOne("phi_radCorr", "C_{rad}", 0.5, 1.5, [](DISANAplotter* P, size_t iq, size_t iw) -> TH1D* {
        auto& R = P->GetPhiRadCorr3D();
        if (iq >= R.size() || iw >= R[iq].size()) return nullptr;
        return R[iq][iw];
      });
    }

    // ---- CSV output: one CSV per (model × Q2-bin × W-bin) ----
    // Requires #include <fstream> and <iomanip> (add at top if not already present)
    if (writeCSV) {
      for (size_t im = 0; im < plotters.size(); ++im)
        for (size_t iq = 0; iq < nQ; ++iq)
          for (size_t iw = 0; iw < nW; ++iw) writeCSVForBin(im, iq, iw);
    }
  }

  // Generic per-(Q2,W) grid plotter of TH1 vs t'
  void PlotPhiPerBin_FromCache(const std::string& tag,                                       // used in output filename + canvas name
                               const std::string& yTitle,                                    // axis title
                               double yMin, double yMax,                                     // y range
                               std::function<TH1D*(DISANAplotter*, size_t, size_t)> getHist  // how to fetch TH1
  ) {
    if (plotters.empty()) return;

    const auto& q2 = fXbins.GetQ2Bins();
    const auto& w = fXbins.GetWBins();
    const bool hasW = !w.empty();

    const size_t nQ = (q2.size() > 1) ? (q2.size() - 1) : 0;
    const size_t nW = hasW ? (w.size() - 1) : 1;

    for (size_t iq = 0; iq < nQ; ++iq) {
      for (size_t iw = 0; iw < nW; ++iw) {
        auto c = new TCanvas(Form("c_%s_Q%zu_W%zu", tag.c_str(), iq, iw), "", 1200, 900);
        styleCrossSection_.StylePad((TPad*)gPad);
        gPad->SetTicks(1, 1);

        TLegend* leg = new TLegend(0.60, 0.72, 0.92, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);

        TString head = hasW ? Form("Q^{2}[%.2f, %.2f]   W[%.1f, %.1f]", q2[iq], q2[iq + 1], w[iw], w[iw + 1]) : Form("Q^{2}[%.2f, %.2f]", q2[iq], q2[iq + 1]);

        bool first = true;

        for (size_t im = 0; im < plotters.size(); ++im) {
          DISANAplotter* P = plotters[im].get();
          TH1D* h = getHist(P, iq, iw);
          if (!h) continue;

          // Style (same logic as your cross-section plot)
          styleCrossSection_.StyleTH1(h);

          auto [cr, cg, cb] = modelShades[im % modelShades.size()];
          const int colorIdx = 5000 + int(im) * 20;
          if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
          h->SetLineColor(colorIdx);
          h->SetMarkerColor(colorIdx);
          h->SetMarkerStyle(20);
          h->SetMarkerSize(1.0);

          h->SetTitle("");
          h->GetXaxis()->SetTitle("-t' [GeV^{2}]");
          h->GetYaxis()->SetTitle(yTitle.c_str());
          h->GetYaxis()->SetRangeUser(yMin, yMax);

          if (first) {
            h->Draw("E1X0");
            TLatex latex;
            latex.SetNDC();
            latex.SetTextFont(42);
            latex.SetTextSize(0.040);
            latex.DrawLatex(0.14, 0.93, head);
            first = false;
          } else {
            h->Draw("E1X0 SAME");
          }

          leg->AddEntry(h, labels[im].c_str(), "lep");
        }

        leg->Draw();

        TString out = hasW ? Form("%s/%s_Q%zu_W%zu.pdf", outputDir.c_str(), tag.c_str(), iq, iw) : Form("%s/%s_Q%zu.pdf", outputDir.c_str(), tag.c_str(), iq);

        c->SaveAs(out);

        delete leg;
        delete c;
      }
    }
  }

  // === NEW: A_LU(cos(theta_KK)) workflow (mirrors inv-mass fits + cached dσ/dt pattern) ===
  void PlotPhiALUCosThetaPerBin_AllModels(const std::string& baseOutDir = "PhiALUCosThetaFits", int nMassBins = 200, double mMin = 0.9874, double mMax = 1.120,
                                          bool constrainSigma = true, double sigmaRef = 0.004, double sigmaFrac = 0.30, double beamPol = 1.0) {
    if (plotters.empty()) {
      std::cerr << "[PlotPhiALUCosThetaPerBin] no models.\n";
      return;
    }
    gSystem->Exec(Form("mkdir -p %s", baseOutDir.c_str()));
    for (size_t i = 0; i < plotters.size(); ++i) {
      const std::string subdir = baseOutDir + "/" + labels[i];
      std::cout << "→ Fitting/drawing helicity-separated K^{+}K^{-} mass per cos(theta_KK) bin for model: " << labels[i] << " → " << subdir << std::endl;
      plotters[i]->MakePhiBSAMassFitCanvases3D(fXbins, subdir, nMassBins, mMin, mMax, constrainSigma, sigmaRef, sigmaFrac, beamPol);
    }
  }
  void PlotPhiALUCosTheta_FromCache(const std::string& outDir = "PhiALUCosTheta") {
    if (plotters.empty()) return;

    const std::string outBase = outputDir + "/" + outDir;
    gSystem->Exec(Form("mkdir -p \"%s\"", outBase.c_str()));

    const auto& q2 = fXbins.GetQ2Bins();
    const auto& w = fXbins.GetWBins();
    const bool hasW = !w.empty();

    const size_t nQ = q2.size() ? q2.size() - 1 : 0;
    const size_t nW = hasW ? (w.size() - 1) : 1;

    for (size_t iq = 0; iq < nQ; ++iq) {
      for (size_t iw = 0; iw < nW; ++iw) {
        auto c = new TCanvas(Form("c_phi_alu_costh_Q%zu_W%zu", iq, iw), "", 1200, 900);
        styleCrossSection_.StylePad((TPad*)gPad);
        gPad->SetFillStyle(4000);
        gPad->SetTicks(1, 1);

        TLegend* leg = new TLegend(0.60, 0.72, 0.92, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);

        double yAbsMax = 0.0;
        for (size_t im = 0; im < plotters.size(); ++im) {
          const auto& bsa3D = plotters[im]->GetPhiALUCosTheta3D();
          if (iq >= bsa3D.size() || iw >= bsa3D[iq].size()) continue;
          TH1D* h = bsa3D[iq][iw];
          if (!h) continue;
          yAbsMax = std::max(yAbsMax, std::fabs(h->GetMaximum()));
          yAbsMax = std::max(yAbsMax, std::fabs(h->GetMinimum()));
        }
        yAbsMax = std::max(0.20, 1.2 * yAbsMax);

        const TString head = hasW ? Form("Q^{2}[%.2f, %.2f]   W[%.1f, %.1f]", q2[iq], q2[iq + 1], w[iw], w[iw + 1]) : Form("Q^{2}[%.2f, %.2f]", q2[iq], q2[iq + 1]);

        bool first = true;
        for (size_t im = 0; im < plotters.size(); ++im) {
          const auto& bsa3D = plotters[im]->GetPhiALUCosTheta3D();
          if (iq >= bsa3D.size() || iw >= bsa3D[iq].size()) continue;
          TH1D* h = bsa3D[iq][iw];
          if (!h) continue;

          styleCrossSection_.StyleTH1(h);
          auto [cr, cg, cb] = modelShades[im % modelShades.size()];
          const int colorIdx = 4000 + int(im) * 20;
          if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
          h->SetLineColor(colorIdx);
          h->SetMarkerColor(colorIdx);
          h->SetMarkerStyle(20);
          h->SetMarkerSize(1.0);
          h->SetLineWidth(1);
          h->SetTitle("");
          h->GetXaxis()->SetTitle("cos#theta_{K^{+}K^{-}}");
          h->GetYaxis()->SetTitle("A_{LU}");
          h->GetXaxis()->CenterTitle(true);
          h->GetYaxis()->CenterTitle(true);
          h->GetXaxis()->SetNdivisions(505);
          h->GetYaxis()->SetNdivisions(510);
          h->GetYaxis()->SetRangeUser(-yAbsMax * 2.0, yAbsMax * 2.0);

          if (first) {
            h->Draw("E1X0");
            TLine z(h->GetXaxis()->GetXmin(), 0.0, h->GetXaxis()->GetXmax(), 0.0);
            z.SetLineStyle(2);
            z.SetLineWidth(2);
            z.Draw("SAME");

            TLatex latex;
            latex.SetNDC();
            latex.SetTextFont(42);
            latex.SetTextSize(0.040);
            latex.DrawLatex(0.14, 0.93, head);
          } else {
            h->Draw("E1X0 SAME");
          }

          leg->AddEntry(h, labels[im].c_str(), "lep");
          first = false;
        }

        leg->Draw();
        c->Update();

        TString out = hasW ? Form("%s/BSA_vs_CosKK_Q%zu_W%zu.pdf", outBase.c_str(), iq, iw) : Form("%s/BSA_vs_CosKK_Q%zu.pdf", outBase.c_str(), iq);
        c->SaveAs(out);

        delete leg;
        delete c;
      }
    }
  }

  void PlotPhiBSATrentoPhiPerBin_AllModels(const std::string& baseOutDir = "PhiBSATrentoPhiFits", int nMassBins = 200, double mMin = 0.9874, double mMax = 1.120,
                                           bool constrainSigma = true, double sigmaRef = 0.004, double sigmaFrac = 0.30, double beamPol = 1.0) {
    if (plotters.empty()) {
      std::cerr << "[PlotPhiBSATrentoPhiPerBin] no models.\n";
      return;
    }

    gSystem->Exec(Form("mkdir -p %s", baseOutDir.c_str()));
    for (size_t i = 0; i < plotters.size(); ++i) {
      const std::string subdir = baseOutDir + "/" + labels[i];
      std::cout << "→ Fitting/drawing helicity-separated K^{+}K^{-} mass per Trento-phi bin for model: " << labels[i] << " → " << subdir << std::endl;

      plotters[i]->MakePhiBSATrentoPhiMassFitCanvases3D(fXbins, subdir, nMassBins, mMin, mMax, constrainSigma, sigmaRef, sigmaFrac, beamPol);
    }
  }

  void PlotPhiBSATrentoPhi_FromCache(const std::string& outDir = "PhiBSATrentoPhi") {
    if (plotters.empty()) return;

    const std::string outBase = outputDir + "/" + outDir;
    gSystem->Exec(Form("mkdir -p \"%s\"", outBase.c_str()));

    const auto& q2 = fXbins.GetQ2Bins();
    const auto& w = fXbins.GetWBins();
    const bool hasW = !w.empty();

    const size_t nQ = q2.size() ? q2.size() - 1 : 0;
    const size_t nW = hasW ? (w.size() - 1) : 1;

    for (size_t iq = 0; iq < nQ; ++iq) {
      for (size_t iw = 0; iw < nW; ++iw) {
        auto c = new TCanvas(Form("c_phi_bsa_trentophi_Q%zu_W%zu", iq, iw), "", 1200, 900);
        styleCrossSection_.StylePad((TPad*)gPad);
        gPad->SetFillStyle(4000);
        gPad->SetTicks(1, 1);

        TLegend* leg = new TLegend(0.60, 0.72, 0.92, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);

        double yAbsMax = 0.0;
        for (size_t im = 0; im < plotters.size(); ++im) {
          const auto& bsa3D = plotters[im]->GetPhiBSATrentoPhi3D();
          if (iq >= bsa3D.size() || iw >= bsa3D[iq].size()) continue;
          TH1D* h = bsa3D[iq][iw];
          if (!h) continue;
          yAbsMax = std::max(yAbsMax, std::fabs(h->GetMaximum()));
          yAbsMax = std::max(yAbsMax, std::fabs(h->GetMinimum()));
        }
        yAbsMax = std::max(0.20, 1.2 * yAbsMax);

        const TString head = hasW ? Form("Q^{2}[%.2f, %.2f]   W[%.1f, %.1f]", q2[iq], q2[iq + 1], w[iw], w[iw + 1]) : Form("Q^{2}[%.2f, %.2f]", q2[iq], q2[iq + 1]);

        bool first = true;
        for (size_t im = 0; im < plotters.size(); ++im) {
          const auto& bsa3D = plotters[im]->GetPhiBSATrentoPhi3D();
          if (iq >= bsa3D.size() || iw >= bsa3D[iq].size()) continue;
          TH1D* h = bsa3D[iq][iw];
          if (!h) continue;

          styleCrossSection_.StyleTH1(h);
          auto [cr, cg, cb] = modelShades[im % modelShades.size()];
          const int colorIdx = 4000 + int(im) * 20;
          if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
          h->SetLineColor(colorIdx);
          h->SetMarkerColor(colorIdx);
          h->SetMarkerStyle(20);
          h->SetMarkerSize(1.0);
          h->SetLineWidth(1);

          h->SetTitle("");
          h->GetXaxis()->SetTitle("#phi_{Trento} [deg]");
          h->GetYaxis()->SetTitle("A_{LU}");
          h->GetXaxis()->CenterTitle(true);
          h->GetYaxis()->CenterTitle(true);
          h->GetXaxis()->SetNdivisions(505);
          h->GetYaxis()->SetNdivisions(510);
          h->GetYaxis()->SetRangeUser(-0.750, 0.750);

          if (first) {
            h->Draw("E1X0");
            TLine z(h->GetXaxis()->GetXmin(), 0.0, h->GetXaxis()->GetXmax(), 0.0);
            z.SetLineStyle(2);
            z.SetLineWidth(2);
            z.Draw("SAME");

            TLatex latex;
            latex.SetNDC();
            latex.SetTextFont(42);
            latex.SetTextSize(0.040);
            latex.DrawLatex(0.14, 0.93, head);
          } else {
            h->Draw("E1X0 SAME");
          }

          leg->AddEntry(h, labels[im].c_str(), "lep");
          first = false;
        }

        leg->Draw();
        c->Update();

        TString out = hasW ? Form("%s/BSA_vs_trentoPhi_Q%zu_W%zu.pdf", outBase.c_str(), iq, iw) : Form("%s/BSA_vs_trentoPhi_Q%zu.pdf", outBase.c_str(), iq);
        c->SaveAs(out);

        delete leg;
        delete c;
      }
    }
  }

  // === NEW: A_LU(z_phi) workflow, analogous to cos(theta_KK) ===
  void PlotPhiALUZPhiPerBin_AllModels(const std::string& baseOutDir = "PhiALUZPhiFits", int nMassBins = 200, double mMin = 0.9874, double mMax = 1.120, bool constrainSigma = true,
                                      double sigmaRef = 0.004, double sigmaFrac = 0.30, double beamPol = 1.0) {
    if (plotters.empty()) {
      std::cerr << "[PlotPhiALUZPhiPerBin_AllModels] no models.\n";
      return;
    }

    gSystem->Exec(Form("mkdir -p %s", baseOutDir.c_str()));

    for (size_t i = 0; i < plotters.size(); ++i) {
      const std::string subdir = baseOutDir + "/" + labels[i];
      std::cout << "→ Fitting/drawing helicity-separated K^{+}K^{-} mass per z_phi bin for model: " << labels[i] << " → " << subdir << std::endl;

      plotters[i]->MakePhiALUZPhiMassFitCanvases3D(fXbins, subdir, nMassBins, mMin, mMax, constrainSigma, sigmaRef, sigmaFrac, beamPol);
    }
  }

  void PlotPhiALUZPhi_FromCache(const std::string& outDir = "PhiALUZPhi") {
    if (plotters.empty()) return;

    const std::string outBase = outputDir + "/" + outDir;
    gSystem->Exec(Form("mkdir -p \"%s\"", outBase.c_str()));

    const auto& q2 = fXbins.GetQ2Bins();
    const auto& w = fXbins.GetWBins();
    const bool hasW = !w.empty();

    const size_t nQ = q2.size() ? q2.size() - 1 : 0;
    const size_t nW = hasW ? (w.size() - 1) : 1;

    for (size_t iq = 0; iq < nQ; ++iq) {
      for (size_t iw = 0; iw < nW; ++iw) {
        auto c = new TCanvas(Form("c_phi_alu_zphi_Q%zu_W%zu", iq, iw), "", 1200, 900);
        styleCrossSection_.StylePad((TPad*)gPad);
        gPad->SetFillStyle(4000);
        gPad->SetTicks(1, 1);

        TLegend* leg = new TLegend(0.60, 0.72, 0.92, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);

        double yAbsMax = 0.0;
        for (size_t im = 0; im < plotters.size(); ++im) {
          const auto& alu3D = plotters[im]->GetPhiALUZPhi3D();
          if (iq >= alu3D.size() || iw >= alu3D[iq].size()) continue;
          TH1D* h = alu3D[iq][iw];
          if (!h) continue;
          yAbsMax = std::max(yAbsMax, std::fabs(h->GetMaximum()));
          yAbsMax = std::max(yAbsMax, std::fabs(h->GetMinimum()));
        }
        yAbsMax = std::max(0.20, 1.2 * yAbsMax);

        const TString head = hasW ? Form("Q^{2}[%.2f, %.2f]   W[%.1f, %.1f]", q2[iq], q2[iq + 1], w[iw], w[iw + 1]) : Form("Q^{2}[%.2f, %.2f]", q2[iq], q2[iq + 1]);

        bool first = true;

        for (size_t im = 0; im < plotters.size(); ++im) {
          const auto& alu3D = plotters[im]->GetPhiALUZPhi3D();
          if (iq >= alu3D.size() || iw >= alu3D[iq].size()) continue;
          TH1D* h = alu3D[iq][iw];
          if (!h) continue;

          styleCrossSection_.StyleTH1(h);
          auto [cr, cg, cb] = modelShades[im % modelShades.size()];
          const int colorIdx = 4000 + int(im) * 20;
          if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
          h->SetLineColor(colorIdx);
          h->SetMarkerColor(colorIdx);
          h->SetMarkerStyle(20);
          h->SetMarkerSize(1.0);
          h->SetLineWidth(1);

          h->SetTitle("");
          h->GetXaxis()->SetTitle("z_{#phi}");
          h->GetYaxis()->SetTitle("A_{LU}");
          h->GetXaxis()->CenterTitle(true);
          h->GetYaxis()->CenterTitle(true);
          h->GetXaxis()->SetNdivisions(505);
          h->GetYaxis()->SetNdivisions(510);
          h->GetYaxis()->SetRangeUser(-yAbsMax, yAbsMax);

          if (first) {
            h->Draw("E1X0");
            TLine z(h->GetXaxis()->GetXmin(), 0.0, h->GetXaxis()->GetXmax(), 0.0);
            z.SetLineStyle(2);
            z.SetLineWidth(2);
            z.Draw("SAME");

            TLatex latex;
            latex.SetNDC();
            latex.SetTextFont(42);
            latex.SetTextSize(0.040);
            latex.DrawLatex(0.14, 0.93, head);
          } else {
            h->Draw("E1X0 SAME");
          }

          leg->AddEntry(h, labels[im].c_str(), "lep");
          first = false;
        }

        leg->Draw();
        c->Update();

        TString out = hasW ? Form("%s/ALU_vs_zphi_Q%zu_W%zu.pdf", outBase.c_str(), iq, iw) : Form("%s/ALU_vs_zphi_Q%zu.pdf", outBase.c_str(), iq);
        c->SaveAs(out);

        delete leg;
        delete c;
      }
    }
  }

  // ================================================================
  //  R = sigma_L / sigma_T  extraction from cos(theta_H) angular fit
  // ================================================================

  /// Step 1 — run mass fits per (Q2, W, t', cos θ) bin and fill the
  ///          r04_00 / R cache in each plotter.
  void PlotPhiRLTPerBin_AllModels(const std::string& baseOutDir = "PhiRLTFits", int nMassBins = 40, double mMin = 0.988, double mMax = 1.15, bool constrainSigma = true,
                                  double sigmaRef = 0.004, double sigmaFrac = 0.25, double beamEnergyGeV = 10.6) {
    if (plotters.empty()) {
      std::cerr << "[PlotPhiRLTPerBin_AllModels] no models.\n";
      return;
    }
    gSystem->Exec(Form("mkdir -p %s", baseOutDir.c_str()));
    for (size_t i = 0; i < plotters.size(); ++i) {
      const std::string subdir = baseOutDir + "/" + labels[i];
      std::cout << "→ Extracting R=sigma_L/sigma_T from cos(theta_KK) angular fit for model: " << labels[i] << " → " << subdir << std::endl;
      plotters[i]->MakePhiRLTFromCosTheta3D(fXbins, subdir, nMassBins, mMin, mMax, constrainSigma, sigmaRef, sigmaFrac, beamEnergyGeV);
    }
  }

  /// Step 2 — plot R and r04_00 vs t' (one canvas per Q2 bin) and
  ///          write a summary CSV.
  void PlotPhiRLT_FromCache(const std::string& outDir = "PhiRLT") {
    if (plotters.empty()) return;

    const std::string outBase = outputDir + "/" + outDir;
    gSystem->Exec(Form("mkdir -p \"%s\"", outBase.c_str()));

    const auto& q2 = fXbins.GetQ2Bins();
    const auto& w = fXbins.GetWBins();
    const bool hasW = !w.empty();
    const size_t nQ = q2.size() > 1 ? q2.size() - 1 : 0;
    const size_t nW = hasW ? (w.size() - 1) : 1;

    // ---- Helper lambda: draw one observable (R or r04_00) ----
    auto drawObs = [&](const std::string& tag, const std::string& yTitle, double yLo, double yHi, std::function<TH1D*(DISANAplotter*, size_t, size_t)> getter, bool drawZeroLine) {
      for (size_t iq = 0; iq < nQ; ++iq) {
        for (size_t iw = 0; iw < nW; ++iw) {
          TCanvas* c = new TCanvas(Form("c_%s_Q%zu_W%zu", tag.c_str(), iq, iw), "", 1200, 900);
          styleCrossSection_.StylePad((TPad*)gPad);
          gPad->SetTicks(1, 1);

          TLegend* leg = new TLegend(0.60, 0.72, 0.92, 0.90);
          leg->SetBorderSize(0);
          leg->SetFillStyle(0);
          leg->SetTextSize(0.035);

          const TString head = hasW ? Form("Q^{2}[%.2f, %.2f]   W[%.1f, %.1f]", q2[iq], q2[iq + 1], w[iw], w[iw + 1]) : Form("Q^{2}[%.2f, %.2f]", q2[iq], q2[iq + 1]);

          bool first = true;
          for (size_t im = 0; im < plotters.size(); ++im) {
            TH1D* h = getter(plotters[im].get(), iq, iw);
            if (!h) continue;
            styleCrossSection_.StyleTH1(h);
            auto [cr, cg, cb] = modelShades[im % modelShades.size()];
            const int colorIdx = 6000 + (int)im * 20;
            if (!gROOT->GetColor(colorIdx)) new TColor(colorIdx, cr, cg, cb);
            h->SetLineColor(colorIdx);
            h->SetMarkerColor(colorIdx);
            h->SetMarkerStyle(20);
            h->SetMarkerSize(1.0);
            h->SetTitle("");
            h->GetXaxis()->SetTitle("-t' [GeV^{2}]");
            h->GetYaxis()->SetTitle(yTitle.c_str());
            h->GetXaxis()->CenterTitle(true);
            h->GetYaxis()->CenterTitle(true);
            h->GetYaxis()->SetRangeUser(yLo, yHi);
            if (first) {
              h->Draw("E1X0");
              if (drawZeroLine) {
                TLine zl(h->GetXaxis()->GetXmin(), 0.0, h->GetXaxis()->GetXmax(), 0.0);
                zl.SetLineStyle(2);
                zl.SetLineWidth(2);
                zl.Draw("SAME");
              }
              TLatex latex;
              latex.SetNDC();
              latex.SetTextFont(42);
              latex.SetTextSize(0.040);
              latex.DrawLatex(0.14, 0.93, head);
              first = false;
            } else {
              h->Draw("E1X0 SAME");
            }
            leg->AddEntry(h, labels[im].c_str(), "lep");
          }
          leg->Draw();
          TString out = hasW ? Form("%s/%s_Q%zu_W%zu.pdf", outBase.c_str(), tag.c_str(), iq, iw) : Form("%s/%s_Q%zu.pdf", outBase.c_str(), tag.c_str(), iq);
          c->SaveAs(out);
          delete leg;
          delete c;
        }
      }
    };

    // ---- Draw R vs t' ----
    drawObs(
        "phi_RLT", "R = #sigma_{L}/#sigma_{T}", 0.0, 5.0,
        [](DISANAplotter* P, size_t iq, size_t iw) -> TH1D* {
          auto& H = P->GetPhiRLT3D();
          if (iq >= H.size() || iw >= H[iq].size()) return nullptr;
          return H[iq][iw];
        },
        false);

    // ---- Draw r04_00 vs t' ----
    drawObs(
        "phi_r04_00", "r^{04}_{00}", 0.0, 1.0,
        [](DISANAplotter* P, size_t iq, size_t iw) -> TH1D* {
          auto& H = P->GetPhiR04_003D();
          if (iq >= H.size() || iw >= H[iq].size()) return nullptr;
          return H[iq][iw];
        },
        false);

    // ---- Write summary CSV ----
    const std::string csvPath = outBase + "/RLT_summary.csv";
    std::ofstream csv(csvPath);
    if (csv.is_open()) {
      csv << "model,Q2_lo,Q2_hi,W_lo,W_hi,tprime_center,r04_00,r04_00_err,R,R_err\n";
      csv << std::setprecision(6) << std::scientific;
      for (size_t im = 0; im < plotters.size(); ++im) {
        auto& Hrlt = plotters[im]->GetPhiRLT3D();
        auto& Hr04 = plotters[im]->GetPhiR04_003D();
        for (size_t iq = 0; iq < nQ; ++iq) {
          for (size_t iw = 0; iw < nW; ++iw) {
            TH1D* hR = (iq < Hrlt.size() && iw < Hrlt[iq].size()) ? Hrlt[iq][iw] : nullptr;
            TH1D* hr04 = (iq < Hr04.size() && iw < Hr04[iq].size()) ? Hr04[iq][iw] : nullptr;
            if (!hR || !hr04) continue;
            for (int b = 1; b <= hR->GetNbinsX(); ++b) {
              csv << labels[im] << "," << q2[iq] << "," << q2[iq + 1] << "," << (hasW ? w[iw] : 0.0) << "," << (hasW ? w[iw + 1] : 0.0) << "," << hR->GetBinCenter(b) << ","
                  << hr04->GetBinContent(b) << "," << hr04->GetBinError(b) << "," << hR->GetBinContent(b) << "," << hR->GetBinError(b) << "\n";
            }
          }
        }
      }
      csv.close();
      std::cout << "[CSV] R=sigma_L/sigma_T summary → " << csvPath << "\n";
    }
  }

  // === A_LU^{sin(phi_Trento)}(cos(theta_KK)) using the moment method ===
  // Uses the invariant-mass window [mMin, mMax] instead of mass fits.
  void PlotPhiALUCosThetaPerBin_AllModels_SinPhiMoment(const std::string& baseOutDir = "PhiALUCosTheta_SinPhiMoment", double mMin = 0.9874, double mMax = 1.120,
                                                       double beamPol = 1.0) {
    if (plotters.empty()) {
      std::cerr << "[PlotPhiALUCosThetaPerBin_AllModels_SinPhiMoment] no models.\n";
      return;
    }

    gSystem->Exec(Form("mkdir -p %s", baseOutDir.c_str()));

    for (size_t i = 0; i < plotters.size(); ++i) {
      const std::string subdir = baseOutDir + "/" + labels[i];

      std::cout << "→ Computing A_{LU}^{sin#phi}(cos#theta_{KK}) via sin(phi_Trento) moment "
                << "in m(K^{+}K^{-}) ∈ [" << mMin << ", " << mMax << "]"
                << " for model: " << labels[i] << " → " << subdir << std::endl;

      plotters[i]->MakePhiALUCosThetaSinPhiMoment3D(fXbins, subdir, mMin, mMax, beamPol);
    }
  }

  // A_LU(cos(theta_KK)) using a fit A*sin(phi)/(1 + b cos(phi))
  // in each cos(theta_KK) bin. Uses invMass_KpKm window [mMin, mMax].
  void PlotPhiALUCosThetaPerBin_AllModels_SinOver1PlusbCosFit(const std::string& baseOutDir = "PhiALUCosTheta_SinOver1PlusbCosFit", double mMin = 1.01, double mMax = 1.03,
                                                              double beamPol = 1.0) {
    if (plotters.empty()) {
      std::cerr << "[PlotPhiALUCosThetaPerBin_AllModels_SinOver1PlusbCosFit] no models.\n";
      return;
    }

    gSystem->Exec(Form("mkdir -p %s", baseOutDir.c_str()));

    for (size_t i = 0; i < plotters.size(); ++i) {
      const std::string subdir = baseOutDir + "/" + labels[i];

      std::cout << "→ Computing A_{LU}(cos#theta_{KK}) via fit "
                << "A*sin(#phi)/(1 + b cos(#phi)) "
                << "in m(K^{+}K^{-}) ∈ [" << mMin << ", " << mMax << "]"
                << " for model: " << labels[i] << " → " << subdir << std::endl;

      plotters[i]->MakePhiALUCosTheta_SinOver1PlusbCosFit3D(fXbins, subdir, mMin, mMax, beamPol);
    }
  }

 private:
  BinManager fXbins;
  bool plotIndividual = false;
  bool useFittedYields_ = true;
  bool applyCorrection = false;

  DrawStyle style_;              // Default style
  DrawStyle styleKin_;           // Kin plot style
  DrawStyle styleDVCS_;          // DVCS plot style
  DrawStyle styleCrossSection_;  // Cross-section plot style
  DrawStyle styleBSA_;           // BSA plot style

  THnSparseD* correctionHist = nullptr;

  std::unique_ptr<ROOT::RDF::RNode> rdf;
  std::string outputDir = ".";

  std::vector<std::unique_ptr<DISANAplotter>> plotters;
  std::vector<std::string> labels;

  std::vector<std::string> particleName = {"e", "p", "#gamma"};
  std::map<std::string, std::string> typeToParticle = {{"el", "electron"},     {"pro", "proton"},   {"pho", "#gamma_{1}"},
                                                       {"pho2", "#gamma_{2}"}, {"kMinus", "K^{-}"}, {"kPlus", "K^{+}"}};
  std::map<std::string, std::string> VarName = {{"p", "p (GeV/#it{c})"}, {"theta", "#theta (rad)"}, {"phi", "#phi(rad)"}, {"vz", "v_{z}(cm)"}};
};
#endif  // DISANA_COMPARER_H