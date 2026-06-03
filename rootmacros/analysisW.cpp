#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVirtualFitter.h>

#include <algorithm>
#include <atomic>
#include <array>
#include <chrono>
#include <cmath>
#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

constexpr double kElectronMassGeV = 0.00051099895;
constexpr double kProtonMassGeV = 0.9382720813;

struct ElectronWRow {
  double W = 0.0;
  double p = 0.0;
  double phi = 0.0;    // degree, expected [0, 360)
  double theta = 0.0;  // degree
  int sector = -1;
};

struct HistPack {
  std::string yVar;
  std::string xVar;
  std::string label;
  int sector = 0;
  int phiThird = -2;
  int thetaBin = -1;
  double thetaMin = 0.0;
  double thetaMax = 0.0;
  double yTarget = 0.0;
  std::string yTitle;
  std::string peakLabel;
  TH2D* displayHist = nullptr;
  TH2D* fitHist = nullptr;
};

std::string Trim(const std::string& s) {
  const auto first = s.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return "";
  const auto last = s.find_last_not_of(" \t\r\n");
  return s.substr(first, last - first + 1);
}

bool ParseElectronWRow(const std::string& line, ElectronWRow& row) {
  std::stringstream ss(line);
  std::string field;
  std::vector<std::string> fields;
  while (std::getline(ss, field, ',')) fields.push_back(Trim(field));
  if (fields.size() < 5) return false;

  try {
    row.W = std::stod(fields[0]);
    row.p = std::stod(fields[1]);
    row.phi = std::stod(fields[2]);
    row.theta = std::stod(fields[3]);
    row.sector = std::stoi(fields[4]);
  } catch (...) {
    return false;
  }

  while (row.phi < 0.0) row.phi += 360.0;
  while (row.phi >= 360.0) row.phi -= 360.0;

  return std::isfinite(row.W) && std::isfinite(row.p) && std::isfinite(row.phi) &&
         std::isfinite(row.theta) && row.W > 0.0 && row.sector >= 1 && row.sector <= 6;
}

std::vector<ElectronWRow> ReadElectronWCSV(const std::string& csvPath) {
  const std::string cleanPath = Trim(csvPath);
  std::ifstream in(cleanPath);
  if (!in.is_open()) {
    throw std::runtime_error("analysisW: cannot open CSV file: " + cleanPath);
  }

  std::vector<ElectronWRow> rows;
  std::string line;
  bool firstLine = true;
  size_t skipped = 0;
  while (std::getline(in, line)) {
    line = Trim(line);
    if (line.empty()) continue;
    if (firstLine) {
      firstLine = false;
      if (line.find("W") != std::string::npos && line.find("e_p") != std::string::npos) continue;
    }

    ElectronWRow row;
    if (ParseElectronWRow(line, row)) {
      rows.push_back(row);
    } else {
      ++skipped;
    }
  }

  std::cout << "[analysisW] Loaded " << rows.size() << " rows from " << cleanPath;
  if (skipped > 0) std::cout << " (skipped " << skipped << " malformed rows)";
  std::cout << std::endl;
  return rows;
}

double Wrap360(double phiDeg) {
  while (phiDeg < 0.0) phiDeg += 360.0;
  while (phiDeg >= 360.0) phiDeg -= 360.0;
  return phiDeg;
}

struct SectorPhiRange {
  bool valid = false;
  double start = 0.0;
  double width = 60.0;
};

double PhiOffsetFromRangeStart(double phiDeg, const SectorPhiRange& range) {
  double offset = Wrap360(phiDeg) - range.start;
  if (offset < 0.0) offset += 360.0;
  return offset;
}

double PhiForSectorAxis(double phiDeg, const SectorPhiRange& range) {
  const double plotStart = range.start > 180.0 ? range.start - 360.0 : range.start;
  return plotStart + PhiOffsetFromRangeStart(phiDeg, range);
}

double SectorPhiPlotStart(const SectorPhiRange& range) {
  return range.start > 180.0 ? range.start - 360.0 : range.start;
}

double ElasticElectronPTrue(double beamEnergyGeV, double thetaDeg) {
  if (beamEnergyGeV <= 0.0) return std::numeric_limits<double>::quiet_NaN();
  const double thetaRad = thetaDeg * M_PI / 180.0;
  const double sinHalfTheta = std::sin(0.5 * thetaRad);
  const double scatteredE = kProtonMassGeV * beamEnergyGeV /
                            (kProtonMassGeV + 2.0 * beamEnergyGeV * sinHalfTheta * sinHalfTheta);
  const double p2 = scatteredE * scatteredE - kElectronMassGeV * kElectronMassGeV;
  return p2 > 0.0 ? std::sqrt(p2) : std::numeric_limits<double>::quiet_NaN();
}

double DeltaPTrueMinusExp(const ElectronWRow& row, double beamEnergyGeV) {
  const double pTrue = ElasticElectronPTrue(beamEnergyGeV, row.theta);
  return std::isfinite(pTrue) ? pTrue - row.p : std::numeric_limits<double>::quiet_NaN();
}

std::streamoff CSVFileSize(const std::string& csvPath) {
  std::ifstream in(csvPath, std::ios::binary | std::ios::ate);
  if (!in.is_open()) {
    throw std::runtime_error("analysisW: cannot open CSV file: " + csvPath);
  }
  return in.tellg();
}

int ResolveThreadCount(int requestedThreads, std::streamoff fileSize) {
  int nThreads = requestedThreads;
  if (nThreads <= 0) nThreads = static_cast<int>(std::thread::hardware_concurrency());
  if (nThreads <= 0) nThreads = 1;
  if (fileSize > 0) {
    const int maxUsefulThreads = std::max<int>(1, static_cast<int>(fileSize / (1024 * 1024)));
    nThreads = std::min(nThreads, maxUsefulThreads);
  }
  return std::max(1, nThreads);
}

std::vector<std::pair<std::streamoff, std::streamoff>> MakeCSVRanges(std::streamoff fileSize, int nThreads) {
  std::vector<std::pair<std::streamoff, std::streamoff>> ranges;
  ranges.reserve(nThreads);
  for (int i = 0; i < nThreads; ++i) {
    const std::streamoff begin = fileSize * i / nThreads;
    const std::streamoff end = fileSize * (i + 1) / nThreads;
    if (begin < end) ranges.push_back({begin, end});
  }
  return ranges;
}

void PrintProgressBar(const std::string& label, long long done, long long total) {
  const int barWidth = 40;
  const double frac = total > 0 ? std::min(1.0, std::max(0.0, static_cast<double>(done) / static_cast<double>(total))) : 1.0;
  const int filled = static_cast<int>(frac * barWidth);

  std::cout << "\r[analysisW] " << label << " [";
  for (int i = 0; i < barWidth; ++i) std::cout << (i < filled ? "=" : " ");
  std::cout << "] " << std::fixed << std::setprecision(1) << (100.0 * frac) << "%" << std::flush;
}

void ProgressMonitor(const std::string& label,
                     std::atomic<long long>& doneBytes,
                     long long totalBytes,
                     std::atomic<bool>& finished) {
  while (!finished.load()) {
    PrintProgressBar(label, doneBytes.load(), totalBytes);
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
  }
  PrintProgressBar(label, totalBytes, totalBytes);
  std::cout << std::endl;
}

template <typename Handler>
size_t ForEachCSVRowInRange(const std::string& csvPath,
                            std::streamoff begin,
                            std::streamoff end,
                            Handler&& handler,
                            std::atomic<long long>* progressBytes = nullptr) {
  std::ifstream in(csvPath);
  if (!in.is_open()) {
    throw std::runtime_error("analysisW: cannot open CSV file: " + csvPath);
  }

  std::streamoff lastProgressPos = begin;
  if (begin > 0) {
    in.seekg(begin);
    std::string partialLine;
    std::getline(in, partialLine);
    const std::streampos pos = in.tellg();
    if (progressBytes && pos != std::streampos(-1)) {
      const std::streamoff current = std::min(end, static_cast<std::streamoff>(pos));
      if (current > lastProgressPos) {
        progressBytes->fetch_add(static_cast<long long>(current - lastProgressPos));
        lastProgressPos = current;
      }
    }
  } else {
    in.seekg(0);
  }

  size_t accepted = 0;
  std::string line;
  while (true) {
    const std::streampos pos = in.tellg();
    if (pos == std::streampos(-1)) break;
    const std::streamoff lineStart = static_cast<std::streamoff>(pos);
    if (lineStart >= end) break;
    if (!std::getline(in, line)) break;

    if (progressBytes) {
      const std::streampos afterPos = in.tellg();
      const std::streamoff current = afterPos == std::streampos(-1) ? end : std::min(end, static_cast<std::streamoff>(afterPos));
      if (current > lastProgressPos) {
        progressBytes->fetch_add(static_cast<long long>(current - lastProgressPos));
        lastProgressPos = current;
      }
    }

    ElectronWRow row;
    if (!ParseElectronWRow(Trim(line), row)) continue;
    handler(row);
    ++accepted;
  }
  if (progressBytes && end > lastProgressPos) {
    progressBytes->fetch_add(static_cast<long long>(end - lastProgressPos));
  }
  return accepted;
}

std::array<SectorPhiRange, 7> MakeFixedSectorPhiRanges(double sector1PhiStartDeg) {
  std::array<SectorPhiRange, 7> ranges;
  for (int sector = 1; sector <= 6; ++sector) {
    ranges[sector].valid = true;
    ranges[sector].start = Wrap360(sector1PhiStartDeg + 60.0 * (sector - 1));
    ranges[sector].width = 60.0;

    const double plotStart = SectorPhiPlotStart(ranges[sector]);
    std::cout << "[analysisW] sector " << sector
              << " fixed phi range: start=" << ranges[sector].start
              << " width=60 deg"
              << " plot=" << plotStart << " to " << plotStart + 60.0 << " deg" << std::endl;
  }
  return ranges;
}

int PhiThirdInSector(double phiDeg, const SectorPhiRange& range) {
  if (!range.valid || range.width <= 0.0) return -1;
  const double offset = PhiOffsetFromRangeStart(phiDeg, range);
  if (offset > range.width) return -1;
  const double thirdWidth = range.width / 3.0;
  if (offset < thirdWidth) return 0;
  if (offset < 2.0 * thirdWidth) return 1;
  return 2;
}

std::string PhiThirdTag(int third) {
  if (third == 0) return "phiLow";
  if (third == 1) return "phiMid";
  return "phiHigh";
}

std::string PhiThirdLabel(int third, const SectorPhiRange& range) {
  if (!range.valid) return Form("#phi third %d", third + 1);
  const double thirdWidth = range.width / 3.0;
  const double lo = Wrap360(range.start + thirdWidth * third);
  const double hi = Wrap360(range.start + thirdWidth * (third + 1));
  return Form("#phi %.1f-%.1f#circ", lo, hi);
}

void PrepareHistPacks(std::vector<HistPack>& hists) {
  for (auto& item : hists) {
    item.displayHist->SetDirectory(nullptr);
    item.fitHist->SetDirectory(nullptr);
  }
}

std::vector<HistPack> CloneHistPacks(const std::vector<HistPack>& hists, const std::string& suffix) {
  std::vector<HistPack> clones;
  clones.reserve(hists.size());
  for (const auto& item : hists) {
    auto* display = dynamic_cast<TH2D*>(item.displayHist->Clone((std::string(item.displayHist->GetName()) + suffix).c_str()));
    auto* fit = dynamic_cast<TH2D*>(item.fitHist->Clone((std::string(item.fitHist->GetName()) + suffix).c_str()));
    display->Reset("ICES");
    fit->Reset("ICES");
    display->SetDirectory(nullptr);
    fit->SetDirectory(nullptr);
    clones.push_back({item.yVar,
                      item.xVar,
                      item.label,
                      item.sector,
                      item.phiThird,
                      item.thetaBin,
                      item.thetaMin,
                      item.thetaMax,
                      item.yTarget,
                      item.yTitle,
                      item.peakLabel,
                      display,
                      fit});
  }
  return clones;
}

void DeleteHistPacks(std::vector<HistPack>& hists) {
  for (auto& item : hists) {
    delete item.displayHist;
    delete item.fitHist;
    item.displayHist = nullptr;
    item.fitHist = nullptr;
  }
}

bool FillHistPacks(const ElectronWRow& row,
                   const std::array<SectorPhiRange, 7>& sectorPhiRanges,
                   std::vector<HistPack>& hists,
                   double beamEnergyGeV) {
  const auto& phiRange = sectorPhiRanges[row.sector];
  if (!phiRange.valid) return false;
  const double phiAxis = PhiForSectorAxis(row.phi, phiRange);
  const double phiPlotStart = SectorPhiPlotStart(phiRange);
  if (phiAxis < phiPlotStart || phiAxis >= phiPlotStart + phiRange.width) return false;

  for (auto& item : hists) {
    if (item.xVar == "phi") {
      double fillPhi = row.phi;
      if (item.sector > 0) {
        if (item.sector != row.sector) continue;
        fillPhi = phiAxis;
      }
      if (item.thetaBin >= 0 && (row.theta < item.thetaMin || row.theta >= item.thetaMax)) continue;
      const double yValue = item.yVar == "deltaP" ? DeltaPTrueMinusExp(row, beamEnergyGeV) : row.W;
      if (!std::isfinite(yValue)) continue;
      item.displayHist->Fill(fillPhi, yValue);
      item.fitHist->Fill(fillPhi, yValue);
    }
  }
  return true;
}

size_t FillHistPacksFromCSVParallel(const std::string& csvPath,
                                    const std::vector<std::pair<std::streamoff, std::streamoff>>& ranges,
                                    const std::array<SectorPhiRange, 7>& sectorPhiRanges,
                                    std::vector<HistPack>& hists,
                                    double beamEnergyGeV) {
  std::vector<std::vector<HistPack>> localHists;
  localHists.reserve(ranges.size());
  for (size_t i = 0; i < ranges.size(); ++i) {
    localHists.push_back(CloneHistPacks(hists, Form("_thread%zu", i)));
  }

  std::vector<size_t> accepted(ranges.size(), 0);
  std::vector<std::thread> workers;
  workers.reserve(ranges.size());

  long long totalBytes = 0;
  for (const auto& range : ranges) totalBytes += static_cast<long long>(range.second - range.first);
  std::atomic<long long> doneBytes{0};
  std::atomic<bool> finished{false};
  std::thread monitor(ProgressMonitor, "fill histograms", std::ref(doneBytes), totalBytes, std::ref(finished));

  for (size_t i = 0; i < ranges.size(); ++i) {
    workers.emplace_back([&, i]() {
      accepted[i] = ForEachCSVRowInRange(csvPath, ranges[i].first, ranges[i].second, [&](const ElectronWRow& row) {
        FillHistPacks(row, sectorPhiRanges, localHists[i], beamEnergyGeV);
      }, &doneBytes);
    });
  }
  for (auto& worker : workers) worker.join();
  finished.store(true);
  monitor.join();

  for (auto& item : hists) {
    item.displayHist->Reset("ICES");
    item.fitHist->Reset("ICES");
  }

  for (auto& oneThreadHists : localHists) {
    for (size_t ih = 0; ih < hists.size(); ++ih) {
      hists[ih].displayHist->Add(oneThreadHists[ih].displayHist);
      hists[ih].fitHist->Add(oneThreadHists[ih].fitHist);
    }
    DeleteHistPacks(oneThreadHists);
  }

  size_t totalAccepted = 0;
  for (const auto count : accepted) totalAccepted += count;
  return totalAccepted;
}

std::string AxisTitle(const std::string& xVar) {
  if (xVar == "theta") return "electron #theta [deg]";
  if (xVar == "phi") return "electron #phi [deg]";
  return xVar;
}

std::string FormatTargetValue(double value) {
  if (std::abs(value) < 5e-7) return "0";
  return Form("%.6f", value);
}

std::string ShiftAxisTitle(const std::string& peakLabel, double yTarget) {
  if (peakLabel.find("#Delta p") != std::string::npos && std::abs(yTarget) < 5e-7) return "#Delta p [GeV]";
  if (peakLabel.find("W") != std::string::npos && std::abs(yTarget - kProtonMassGeV) < 1e-4) return "W - m_{p} [GeV]";
  return Form("%s - %s [GeV]", peakLabel.c_str(), FormatTargetValue(yTarget).c_str());
}

void Style2D(TH2D* h, const std::string& xVar, const std::string& yTitle) {
  h->SetStats(0);
  h->SetTitle("");
  h->GetXaxis()->SetTitle(AxisTitle(xVar).c_str());
  h->GetYaxis()->SetTitle(yTitle.c_str());
  h->GetXaxis()->SetTitleSize(0.052);
  h->GetYaxis()->SetTitleSize(0.052);
  h->GetXaxis()->SetLabelSize(0.044);
  h->GetYaxis()->SetLabelSize(0.044);
  h->GetXaxis()->SetTitleOffset(1.05);
  h->GetYaxis()->SetTitleOffset(1.20);
}

void StyleGraph(TGraphErrors* g, const std::string& xVar, double yTarget, const std::string& peakLabel) {
  g->SetTitle("");
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.85);
  g->SetMarkerColor(kBlack);
  g->SetLineColor(kBlack);
  g->SetLineWidth(2);
  g->GetXaxis()->SetTitle(AxisTitle(xVar).c_str());
  g->GetYaxis()->SetTitle(ShiftAxisTitle(peakLabel, yTarget).c_str());
  g->GetXaxis()->SetTitleSize(0.052);
  g->GetYaxis()->SetTitleSize(0.052);
  g->GetXaxis()->SetLabelSize(0.044);
  g->GetYaxis()->SetLabelSize(0.044);
  g->GetXaxis()->SetTitleOffset(1.05);
  g->GetYaxis()->SetTitleOffset(1.22);
}

struct PeakConfig {
  int minEntries = 80;
  double fitWindowSigma = 1.6;
  double maxPeakError = 0.08;
};

TF1 MakeAsymGaussianFit(const std::string& name, double fitLo, double fitHi) {
  TF1 fit(name.c_str(),
          [](double* x, double* p) {
            const double xx = x[0];
            const double amp = p[0];
            const double mu = p[1];
            const double sigmaL = std::max(std::abs(p[2]), 1e-6);
            const double sigmaR = std::max(std::abs(p[3]), 1e-6);
            const double sigma = (xx < mu) ? sigmaL : sigmaR;
            const double signal = amp * std::exp(-0.5 * std::pow((xx - mu) / sigma, 2));
            const double bkg = p[4] + p[5] * xx;
            return signal + bkg;
          },
          fitLo, fitHi, 6);
  fit.SetParNames("A", "mu", "sigmaL", "sigmaR", "c0", "c1");
  return fit;
}

std::pair<double, double> EstimatePeakFitRange(TH1D* proj, int maxBin, double fallbackSigma, double windowSigma) {
  const double xmin = proj->GetXaxis()->GetXmin();
  const double xmax = proj->GetXaxis()->GetXmax();
  const double peakX = proj->GetBinCenter(maxBin);
  const double halfMax = 0.5 * proj->GetBinContent(maxBin);

  int loBin = maxBin;
  while (loBin > 1 && proj->GetBinContent(loBin) > halfMax) --loBin;
  int hiBin = maxBin;
  while (hiBin < proj->GetNbinsX() && proj->GetBinContent(hiBin) > halfMax) ++hiBin;

  double lo = proj->GetBinCenter(loBin);
  double hi = proj->GetBinCenter(hiBin);

  const double minWidth = std::max(0.012, 1.2 * proj->GetXaxis()->GetBinWidth(1));
  if ((hi - lo) < minWidth) {
    const double sigma = fallbackSigma > 0.0 ? fallbackSigma : 0.025;
    lo = peakX - windowSigma * sigma;
    hi = peakX + windowSigma * sigma;
  } else {
    const double pad = 0.75 * (hi - lo);
    lo -= pad;
    hi += pad;
  }

  return {std::max(xmin, lo), std::min(xmax, hi)};
}

TGraphErrors* MakeWPeakGraph(TH2D* hist,
                             const std::string& outTxtPath,
                             const std::string& fitCheckDir,
                             const std::string& outputName,
                             double yTarget,
                             const std::string& yTitle,
                             const std::string& peakLabel,
                             const PeakConfig& cfg,
                             bool saveDiagnostics) {
  std::unique_ptr<std::ofstream> fout;
  if (saveDiagnostics) {
    fout = std::make_unique<std::ofstream>(outTxtPath);
    *fout << "# x_center peak delta_peak peak_error sigmaL sigmaR entries chi2_ndf\n";
    gSystem->Exec(("mkdir -p " + fitCheckDir).c_str());
  }

  auto* g = new TGraphErrors();
  g->SetName((outputName + "_peakGraph").c_str());

  const bool skipPhiEdgeBins = outputName.find("_vs_phi") != std::string::npos;
  int ip = 0;
  for (int ix = 1; ix <= hist->GetNbinsX(); ++ix) {
    if (skipPhiEdgeBins && (ix == 1 || ix == hist->GetNbinsX())) continue;
    std::unique_ptr<TH1D> proj(hist->ProjectionY(Form("%s_py_%d", outputName.c_str(), ix), ix, ix, "e"));
    if (!proj || proj->GetEntries() < cfg.minEntries) continue;

    const int maxBin = proj->GetMaximumBin();
    const double mu0 = proj->GetBinCenter(maxBin);
    const double rms = proj->GetRMS();
    if (rms <= 0.0) continue;

    const auto [fitLo, fitHi] = EstimatePeakFitRange(proj.get(), maxBin, rms, cfg.fitWindowSigma);
    if (fitHi <= fitLo) continue;

    TF1 fit = MakeAsymGaussianFit(Form("%s_asymGaus_%d", outputName.c_str(), ix), fitLo, fitHi);
    const double sigmaSeed = std::max(0.01, std::min(rms, 0.08));
    fit.SetParameters(proj->GetMaximum(), mu0, sigmaSeed, sigmaSeed, 0.0, 0.0);
    fit.SetParLimits(0, 0.0, 10.0 * std::max(1.0, proj->GetMaximum()));
    fit.SetParLimits(1, fitLo, fitHi);
    fit.SetParLimits(2, 0.003, 0.20);
    fit.SetParLimits(3, 0.003, 0.20);
    const int fitStatus = proj->Fit(&fit, "QNR");
    if (fitStatus != 0) continue;

    const double peak = fit.GetParameter(1);
    const double peakErr = fit.GetParError(1);
    const double sigmaL = std::abs(fit.GetParameter(2));
    const double sigmaR = std::abs(fit.GetParameter(3));
    const double chi2ndf = fit.GetNDF() > 0 ? fit.GetChisquare() / fit.GetNDF() : 0.0;
    if (!std::isfinite(peak) || !std::isfinite(peakErr) || peakErr <= 0.0 || peakErr > cfg.maxPeakError) continue;

    const double x = hist->GetXaxis()->GetBinCenter(ix);
    const double ex = 0.5 * hist->GetXaxis()->GetBinWidth(ix);
    const double deltaPeak = peak - yTarget;

    if (saveDiagnostics) {
      auto* cfit = new TCanvas(Form("cfit_%s_%d", outputName.c_str(), ix), "", 900, 700);
      cfit->SetMargin(0.14, 0.05, 0.13, 0.06);
      cfit->SetTicks(1, 1);
      proj->SetStats(0);
      proj->SetTitle("");
      proj->SetLineColor(kAzure + 2);
      proj->SetLineWidth(3);
      proj->GetXaxis()->SetTitle(yTitle.c_str());
      proj->GetYaxis()->SetTitle("Counts");
      proj->GetXaxis()->SetRangeUser(hist->GetYaxis()->GetXmin(), hist->GetYaxis()->GetXmax());
      proj->SetMaximum(proj->GetMaximum() > 0.0 ? proj->GetMaximum() * 1.30 : 1.0);
      proj->Draw("HIST");

      fit.SetLineColor(kRed + 1);
      fit.SetLineWidth(3);
      fit.DrawCopy("SAME");

      TLine target(yTarget, 0.0, yTarget, proj->GetMaximum() * 0.85);
      target.SetLineColor(kGray + 2);
      target.SetLineStyle(2);
      target.SetLineWidth(2);
      target.Draw("SAME");

      TLatex fitLatex;
      fitLatex.SetNDC();
      fitLatex.SetTextFont(42);
      fitLatex.SetTextSize(0.038);
      fitLatex.DrawLatex(0.18, 0.86, Form("%s", outputName.c_str()));
      fitLatex.DrawLatex(0.18, 0.80, Form("%s bin %.3f #pm %.3f", hist->GetXaxis()->GetTitle(), x, ex));
      fitLatex.DrawLatex(0.18, 0.74, Form("%s=%.5f #pm %.5f GeV", peakLabel.c_str(), peak, peakErr));
      fitLatex.DrawLatex(0.18, 0.68, Form("#sigma_{L}=%.5f, #sigma_{R}=%.5f", sigmaL, sigmaR));
      fitLatex.DrawLatex(0.18, 0.62, Form("N=%.0f, #chi^{2}/ndf=%.2f", proj->GetEntries(), chi2ndf));

      cfit->SaveAs(Form("%s/%s_bin%03d.png", fitCheckDir.c_str(), outputName.c_str(), ix));
      delete cfit;
    }

    g->SetPoint(ip, x, deltaPeak);
    g->SetPointError(ip, ex, peakErr);

    if (fout) {
      *fout << std::fixed << std::setprecision(8)
            << x << " " << peak << " " << deltaPeak << " " << peakErr << " "
            << sigmaL << " " << sigmaR << " " << proj->GetEntries() << " " << chi2ndf << "\n";
    }
    ++ip;
  }

  return g;
}

void WriteFitParameters(const std::string& outPath, const std::string& formula, TF1* fit) {
  std::ofstream fout(outPath);
  fout << "# Fit expression: " << formula << "\n";
  if (!fit) {
    fout << "# Fit was not performed.\n";
    return;
  }
  fout << std::setprecision(12);
  for (int i = 0; i < fit->GetNpar(); ++i) {
    fout << "p" << i << " = " << fit->GetParameter(i)
         << " +/- " << fit->GetParError(i) << "\n";
  }
}

TGraphErrors* BuildFitGraphWithErrorFloor(TGraphErrors* g, double yErrorFloor) {
  if (!g) return nullptr;
  auto* out = new TGraphErrors();
  int ip = 0;
  for (int i = 0; i < g->GetN(); ++i) {
    double x = 0.0, y = 0.0;
    g->GetPoint(i, x, y);
    const double ey = std::max(g->GetErrorY(i), yErrorFloor);
    out->SetPoint(ip, x, y);
    out->SetPointError(ip, 0.0, ey);
    ++ip;
  }
  return out;
}

TGraphErrors* BuildRobustQuadraticFitGraph(TGraphErrors* g,
                                           const std::string& name,
                                           const std::string& fitFormula,
                                           double xFitMin,
                                           double xFitMax,
                                           double yErrorFloor) {
  std::unique_ptr<TGraphErrors> floored(BuildFitGraphWithErrorFloor(g, yErrorFloor));
  if (!floored || floored->GetN() < 2) return floored.release();

  TF1 firstFit((name + "_prefit").c_str(), fitFormula.c_str(), xFitMin, xFitMax);
  floored->Fit(&firstFit, "QEX0");

  std::vector<double> absResiduals;
  absResiduals.reserve(floored->GetN());
  for (int i = 0; i < floored->GetN(); ++i) {
    double x = 0.0, y = 0.0;
    floored->GetPoint(i, x, y);
    absResiduals.push_back(std::abs(y - firstFit.Eval(x)));
  }
  if (absResiduals.empty()) return floored.release();

  std::nth_element(absResiduals.begin(), absResiduals.begin() + absResiduals.size() / 2, absResiduals.end());
  const double medianAbsResidual = absResiduals[absResiduals.size() / 2];
  const double minRobustSigma = yErrorFloor > 0.0 ? 2.0 * yErrorFloor : 0.0;
  const double robustSigma = std::max(1.4826 * medianAbsResidual, minRobustSigma);
  const double residualCut = 3.5 * robustSigma;

  auto* robust = new TGraphErrors();
  int ip = 0;
  for (int i = 0; i < floored->GetN(); ++i) {
    double x = 0.0, y = 0.0;
    floored->GetPoint(i, x, y);
    if (std::abs(y - firstFit.Eval(x)) > residualCut) continue;
    robust->SetPoint(ip, x, y);
    robust->SetPointError(ip, 0.0, floored->GetErrorY(i));
    ++ip;
  }

  if (robust->GetN() < 2) {
    delete robust;
    return floored.release();
  }
  return robust;
}

TF1* DrawQuadraticFitWithBand(TGraphErrors* g,
                              const std::string& name,
                              const std::string& fitFormula,
                              double xFitMin,
                              double xFitMax,
                              double yErrorFloor = 0.003) {
  if (!g || g->GetN() < 2) return nullptr;

  std::unique_ptr<TGraphErrors> fitGraph(BuildRobustQuadraticFitGraph(g, name, fitFormula, xFitMin, xFitMax, yErrorFloor));
  if (!fitGraph || fitGraph->GetN() < 2) return nullptr;

  auto* fit = new TF1(name.c_str(), fitFormula.c_str(), xFitMin, xFitMax);
  fitGraph->Fit(fit, "QEX0");

  auto* fitBand = new TGraphErrors(200);
  fitBand->SetName((name + "_band").c_str());
  for (int i = 0; i < 200; ++i) {
    const double x = xFitMin + (xFitMax - xFitMin) * static_cast<double>(i) / 199.0;
    fitBand->SetPoint(i, x, fit->Eval(x));
    fitBand->SetPointError(i, 0.0, 0.0);
  }

  TVirtualFitter* fitter = TVirtualFitter::GetFitter();
  if (fitter) {
      fitter->GetConfidenceIntervals(fitBand, 0.68);
      fitBand->SetFillColorAlpha(kRed, 0.22);
      fitBand->SetLineColor(kRed);
      fitBand->SetBit(TObject::kCanDelete);
      fitBand->Draw("E3 SAME");
      g->Draw("PE SAME");
  }

  fit->SetLineColor(kRed + 1);
  fit->SetLineWidth(3);
  fit->Draw("SAME");
  return fit;
}

void DrawWMapOnly(TH2D* displayHist,
                  const std::string& outBase,
                  const std::string& xVar,
                  const std::string& label,
                  const std::string& yTitle) {
  gSystem->Exec(("mkdir -p " + outBase).c_str());

  const bool isPhi = (xVar == "phi");
  auto* c2 = new TCanvas((std::string("c2_") + displayHist->GetName()).c_str(), "", isPhi ? 1800 : 1100, isPhi ? 650 : 850);
  c2->SetMargin(0.13, 0.14, 0.13, 0.06);
  c2->SetTicks(1, 1);
  Style2D(displayHist, xVar, yTitle);
  displayHist->Draw("COLZ");
  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(42);
  latex.SetTextSize(0.04);
  latex.DrawLatex(0.16, 0.88, Form("%s, %s vs %s", label.c_str(), yTitle.c_str(), xVar.c_str()));
  c2->SaveAs((outBase + "/" + std::string(displayHist->GetName()) + ".png").c_str());
  delete c2;
}

TGraphErrors* DrawOneWMapAndPeaks(TH2D* displayHist,
                                  TH2D* fitHist,
                                  const std::string& outBase,
                                  const std::string& xVar,
                                  const std::string& label,
                                  double yTarget,
                                  const std::string& yTitle,
                                  const std::string& peakLabel,
                                  const PeakConfig& peakCfg,
                                  bool savePlots) {
  gSystem->Exec(("mkdir -p " + outBase).c_str());

  const bool isPhi = (xVar == "phi");
  if (savePlots) {
    auto* c2 = new TCanvas((std::string("c2_") + displayHist->GetName()).c_str(), "", isPhi ? 1800 : 1100, isPhi ? 650 : 850);
    c2->SetMargin(0.13, 0.14, 0.13, 0.06);
    c2->SetTicks(1, 1);
    Style2D(displayHist, xVar, yTitle);
    displayHist->Draw("COLZ");
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.16, 0.88, Form("%s, %s vs %s", label.c_str(), yTitle.c_str(), xVar.c_str()));
    c2->SaveAs((outBase + "/" + std::string(displayHist->GetName()) + ".png").c_str());
    delete c2;
  }

  const std::string pointFile = outBase + "/" + std::string(displayHist->GetName()) + "_points.txt";
  const std::string fitCheckDir = outBase + "/FitChecks_" + std::string(displayHist->GetName());
  Style2D(fitHist, xVar, yTitle);
  std::unique_ptr<TGraphErrors> g(MakeWPeakGraph(fitHist,
                                                 pointFile,
                                                 fitCheckDir,
                                                 displayHist->GetName(),
                                                 yTarget,
                                                 yTitle,
                                                 peakLabel,
                                                 peakCfg,
                                                 savePlots));
  StyleGraph(g.get(), xVar, yTarget, peakLabel);

  if (!savePlots) {
    return dynamic_cast<TGraphErrors*>(g->Clone((std::string(displayHist->GetName()) + "_peakGraphClone").c_str()));
  }

  auto* c1 = new TCanvas((std::string("c1_") + displayHist->GetName()).c_str(), "", isPhi ? 1800 : 1100, isPhi ? 650 : 850);
  c1->SetMargin(0.14, 0.05, 0.13, 0.06);
  c1->SetTicks(1, 1);
  g->Draw("APE");

  const std::string fitFormula = "[0] + [1]*x";
  TF1* fit = DrawQuadraticFitWithBand(g.get(),
                                      std::string("fit_") + displayHist->GetName(),
                                      fitFormula,
                                      fitHist->GetXaxis()->GetXmin(),
                                      fitHist->GetXaxis()->GetXmax());

  TLine zero(g->GetXaxis()->GetXmin(), 0.0, g->GetXaxis()->GetXmax(), 0.0);
  zero.SetLineColor(kGray + 2);
  zero.SetLineStyle(2);
  zero.SetLineWidth(2);
  zero.Draw("SAME");

  TLatex peakLatex;
  peakLatex.SetNDC();
  peakLatex.SetTextFont(42);
  peakLatex.SetTextSize(0.04);
  peakLatex.DrawLatex(0.17, 0.88, Form("%s, %s shift", label.c_str(), peakLabel.c_str()));

  c1->SaveAs((outBase + "/" + std::string(displayHist->GetName()) + "_peakShift.png").c_str());

  WriteFitParameters(outBase + "/" + std::string(displayHist->GetName()) + "_fit.txt", fitFormula, fit);

  auto* graphOut = dynamic_cast<TGraphErrors*>(g->Clone((std::string(displayHist->GetName()) + "_peakGraphClone").c_str()));
  delete fit;
  delete c1;
  return graphOut;
}

void DrawSectorPhiThirdComparison(const std::array<TGraphErrors*, 3>& graphs,
                                  int sector,
                                  const SectorPhiRange& phiRange,
                                  const std::string& outBase,
                                  double targetW) {
  TGraphErrors* first = nullptr;
  for (auto* g : graphs) {
    if (g && g->GetN() > 0) {
      first = g;
      break;
    }
  }
  if (!first) return;

  auto* c = new TCanvas(Form("c_csv_S%d_phiThird_theta_compare", sector), "", 1100, 850);
  c->SetMargin(0.14, 0.05, 0.13, 0.06);
  c->SetTicks(1, 1);

  TH1F* frame = c->DrawFrame(0.0, -0.08, 40.0, 0.08);
  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("electron #theta [deg]");
  frame->GetYaxis()->SetTitle(Form("W_{peak} - %.6f [GeV]", targetW));
  frame->GetXaxis()->SetTitleSize(0.052);
  frame->GetYaxis()->SetTitleSize(0.052);
  frame->GetXaxis()->SetLabelSize(0.044);
  frame->GetYaxis()->SetLabelSize(0.044);
  frame->GetXaxis()->SetTitleOffset(1.05);
  frame->GetYaxis()->SetTitleOffset(1.22);

  const int colors[3] = {kAzure + 2, kOrange + 7, kGreen + 2};
  const int markers[3] = {20, 21, 22};

  TLegend* leg = new TLegend(0.62, 0.68, 0.90, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.038);

  for (int i = 0; i < 3; ++i) {
    if (!graphs[i] || graphs[i]->GetN() == 0) continue;
    graphs[i]->SetMarkerColor(colors[i]);
    graphs[i]->SetLineColor(colors[i]);
    graphs[i]->SetMarkerStyle(markers[i]);
    graphs[i]->SetMarkerSize(0.95);
    graphs[i]->SetLineWidth(2);
    graphs[i]->Draw("PE SAME");
    leg->AddEntry(graphs[i], PhiThirdLabel(i, phiRange).c_str(), "p");
  }

  TLine zero(0.0, 0.0, 40.0, 0.0);
  zero.SetLineColor(kGray + 2);
  zero.SetLineStyle(2);
  zero.SetLineWidth(2);
  zero.Draw("SAME");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(42);
  latex.SetTextSize(0.04);
  latex.DrawLatex(0.17, 0.88, Form("S%d, #phi ranges", sector));
  leg->Draw();

  const std::string out = outBase + Form("/electron_csv_S%d_phiThird_WPeak_vs_theta", sector);
  c->SaveAs((out + ".png").c_str());
  delete c;
}

void DrawThetaBinSectorSummary(const std::array<TGraphErrors*, 7>& graphs,
                               const std::array<SectorPhiRange, 7>& sectorPhiRanges,
                               int thetaBin,
                               double thetaMin,
                               double thetaMax,
                               const std::string& outBase,
                               double yTarget,
                               const std::string& peakLabel,
                               const std::string& quantityTag,
                               std::array<TF1*, 7>* fitClones = nullptr) {
  bool hasAny = false;
  for (int sector = 1; sector <= 6; ++sector) {
    if (graphs[sector] && graphs[sector]->GetN() > 0) {
      hasAny = true;
      break;
    }
  }
  if (!hasAny) return;

  auto* c = new TCanvas(Form("c_%s_thetaBin%02d_allSectors_phi", quantityTag.c_str(), thetaBin), "", 1800, 1100);
  c->Divide(3, 2);
  std::vector<TF1*> fitsToDelete;

  for (int sector = 1; sector <= 6; ++sector) {
    c->cd(sector);
    gPad->SetMargin(0.16, 0.06, 0.15, 0.10);
    gPad->SetTicks(1, 1);

    const double phiStart = SectorPhiPlotStart(sectorPhiRanges[sector]);
    const double phiEnd = phiStart + sectorPhiRanges[sector].width;
    TH1F* frame = gPad->DrawFrame(phiStart, -0.08, phiEnd, 0.08);
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("electron #phi [deg]");
    frame->GetYaxis()->SetTitle(ShiftAxisTitle(peakLabel, yTarget).c_str());
    frame->GetXaxis()->SetTitleSize(0.055);
    frame->GetYaxis()->SetTitleSize(0.055);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetXaxis()->SetTitleOffset(1.05);
    frame->GetYaxis()->SetTitleOffset(1.25);

    TLine zero(phiStart, 0.0, phiEnd, 0.0);
    zero.SetLineColor(kGray + 2);
    zero.SetLineStyle(2);
    zero.SetLineWidth(2);
    zero.Draw("SAME");

    if (graphs[sector] && graphs[sector]->GetN() > 0) {
      graphs[sector]->SetMarkerStyle(20);
      graphs[sector]->SetMarkerSize(0.9);
      graphs[sector]->SetMarkerColor(kAzure + 2);
      graphs[sector]->SetLineColor(kAzure + 2);
      graphs[sector]->Draw("PE SAME");
      fitsToDelete.push_back(DrawQuadraticFitWithBand(graphs[sector],
                                                      Form("fit_summary_thetaBin%02d_S%d", thetaBin, sector),
                                                      "[0] + [1]*x",
                                                      phiStart,
                                                      phiEnd));
      if (fitClones && fitsToDelete.back()) {
        (*fitClones)[sector] = dynamic_cast<TF1*>(fitsToDelete.back()->Clone(Form("fitClone_%s_thetaBin%02d_S%d", quantityTag.c_str(), thetaBin, sector)));
      }
    }

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.055);
    latex.DrawLatex(0.20, 0.80, Form("S%d", sector));
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.20, 0.73, Form("#theta %.1f-%.1f#circ", thetaMin, thetaMax));
  }

  const std::string out = outBase + Form("/thetaBin%02d_allSectors_%sPeak_vs_phi", thetaBin, quantityTag.c_str());
  c->SaveAs((out + ".png").c_str());
  for (auto* fit : fitsToDelete) delete fit;
  delete c;
}

void DrawSectorThetaStrip(const std::vector<std::array<TGraphErrors*, 7>>& summaryGraphs,
                          const std::array<SectorPhiRange, 7>& sectorPhiRanges,
                          const std::vector<std::pair<double, double>>& thetaRanges,
                          int sector,
                          const std::string& outBase,
                          double targetW) {
  bool hasAny = false;
  for (const auto& graphs : summaryGraphs) {
    if (graphs[sector] && graphs[sector]->GetN() > 0) {
      hasAny = true;
      break;
    }
  }
  if (!hasAny) return;

  const int nThetaBins = static_cast<int>(summaryGraphs.size());
  const int canvasWidth = std::max(1200, 320 * nThetaBins);
  auto* c = new TCanvas(Form("c_S%d_allThetaBins_phi_strip", sector), "", canvasWidth, 360);
  c->Divide(nThetaBins, 1, 0.001, 0.001);

  const double phiStart = SectorPhiPlotStart(sectorPhiRanges[sector]);
  const double phiEnd = phiStart + sectorPhiRanges[sector].width;

  for (int thetaBin = 0; thetaBin < nThetaBins; ++thetaBin) {
    c->cd(thetaBin + 1);
    gPad->SetMargin(0.18, 0.04, 0.20, 0.12);
    gPad->SetTicks(1, 1);

    TH1F* frame = gPad->DrawFrame(phiStart, -0.08, phiEnd, 0.08);
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("#phi [deg]");
    frame->GetYaxis()->SetTitle(Form("W_{peak} - %.6f", targetW));
    frame->GetXaxis()->SetTitleSize(0.070);
    frame->GetYaxis()->SetTitleSize(0.070);
    frame->GetXaxis()->SetLabelSize(0.060);
    frame->GetYaxis()->SetLabelSize(0.060);
    frame->GetXaxis()->SetTitleOffset(1.05);
    frame->GetYaxis()->SetTitleOffset(1.25);

    TLine zero(phiStart, 0.0, phiEnd, 0.0);
    zero.SetLineColor(kGray + 2);
    zero.SetLineStyle(2);
    zero.SetLineWidth(2);
    zero.Draw("SAME");

    auto* graph = summaryGraphs[thetaBin][sector];
    if (graph && graph->GetN() > 0) {
      graph->SetMarkerStyle(20);
      graph->SetMarkerSize(0.85);
      graph->SetMarkerColor(kAzure + 2);
      graph->SetLineColor(kAzure + 2);
      graph->Draw("PE SAME");
    }

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.070);
    latex.DrawLatex(0.22, 0.86, Form("#theta %.1f-%.1f#circ", thetaRanges[thetaBin].first, thetaRanges[thetaBin].second));
  }

  const std::string out = outBase + Form("/Sector%d_allThetaBins_WPeak_vs_phi_strip", sector);
  c->SaveAs((out + ".png").c_str());
  delete c;
}

void DrawFitParameterThetaSummary(const std::array<TGraphErrors*, 7>& graphs,
                                  const std::string& paramName,
                                  const std::string& yTitle,
                                  const std::string& outBase,
                                  int polyOrder) {
  bool hasAny = false;
  double xMin = 1e9;
  double xMax = -1e9;
  double yMin = 1e9;
  double yMax = -1e9;
  for (int sector = 1; sector <= 6; ++sector) {
    if (graphs[sector] && graphs[sector]->GetN() > 0) {
      hasAny = true;
      for (int i = 0; i < graphs[sector]->GetN(); ++i) {
        double x = 0.0;
        double y = 0.0;
        graphs[sector]->GetPoint(i, x, y);
        const double ex = graphs[sector]->GetErrorX(i);
        const double ey = graphs[sector]->GetErrorY(i);
        xMin = std::min(xMin, x - ex);
        xMax = std::max(xMax, x + ex);
        yMin = std::min(yMin, y - ey);
        yMax = std::max(yMax, y + ey);
      }
    }
  }
  if (!hasAny) return;

  if (!(xMax > xMin)) {
    xMin = 5.5;
    xMax = 25.5;
  } else {
    const double xPad = 0.04 * (xMax - xMin);
    xMin -= xPad;
    xMax += xPad;
  }
  if (!(yMax > yMin)) {
    yMin -= 0.01;
    yMax += 0.01;
  } else {
    const double yPad = 0.18 * (yMax - yMin);
    yMin -= yPad;
    yMax += yPad;
  }
  if (yMin > 0.0) yMin = 0.0;
  if (yMax < 0.0) yMax = 0.0;

  auto* c = new TCanvas(Form("c_deltaP_%s_vs_theta_allSectors", paramName.c_str()), "", 1800, 1100);
  c->Divide(3, 2);
  std::vector<TF1*> fitsToDelete;
  const int requestedOrder = std::max(0, polyOrder);

  for (int sector = 1; sector <= 6; ++sector) {
    c->cd(sector);
    gPad->SetMargin(0.16, 0.06, 0.15, 0.10);
    gPad->SetTicks(1, 1);

    TH1F* frame = gPad->DrawFrame(xMin, yMin, xMax, yMax);
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("electron #theta [deg]");
    frame->GetYaxis()->SetTitle(yTitle.c_str());
    frame->GetXaxis()->SetTitleSize(0.055);
    frame->GetYaxis()->SetTitleSize(0.055);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetXaxis()->SetTitleOffset(1.05);
    frame->GetYaxis()->SetTitleOffset(paramName == "a_1" ? 1.55 : 1.25);

    TLine zero(xMin, 0.0, xMax, 0.0);
    zero.SetLineColor(kGray + 2);
    zero.SetLineStyle(2);
    zero.SetLineWidth(2);
    zero.Draw("SAME");

    if (graphs[sector] && graphs[sector]->GetN() > 0) {
      graphs[sector]->SetMarkerStyle(20);
      graphs[sector]->SetMarkerSize(0.9);
      graphs[sector]->SetMarkerColor(kAzure + 2);
      graphs[sector]->SetLineColor(kAzure + 2);
      graphs[sector]->Draw("PE SAME");
      if (graphs[sector]->GetN() >= 2) {
        const int fitOrder = std::min(requestedOrder, graphs[sector]->GetN() - 1);
        const std::string fitFormula = Form("pol%d", fitOrder);
        TF1* fit = DrawQuadraticFitWithBand(graphs[sector],
                                            Form("fit_deltaP_%s_vs_theta_S%d_pol%d", paramName.c_str(), sector, fitOrder),
                                            fitFormula,
                                            xMin,
                                            xMax,
                                            0.0);
        fitsToDelete.push_back(fit);
        WriteFitParameters(outBase + Form("/deltaP_%s_vs_theta_S%d_pol%d_fit.txt", paramName.c_str(), sector, fitOrder),
                           fitFormula,
                           fit);
      }
    }

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.055);
    latex.DrawLatex(0.20, 0.80, Form("S%d", sector));
  }

  const std::string out = outBase + Form("/deltaP_%s_vs_theta_allSectors", paramName.c_str());
  c->SaveAs((out + ".png").c_str());
  for (auto* fit : fitsToDelete) delete fit;
  delete c;
}

void DrawGlobalPhiThetaStrip(const std::vector<TGraphErrors*>& graphs,
                             const std::vector<std::pair<double, double>>& thetaRanges,
                             const std::string& outBase,
                             double targetW) {
  bool hasAny = false;
  for (auto* graph : graphs) {
    if (graph && graph->GetN() > 0) {
      hasAny = true;
      break;
    }
  }
  if (!hasAny) return;

  const int nThetaBins = static_cast<int>(graphs.size());
  const int canvasWidth = std::max(1400, 360 * nThetaBins);
  auto* c = new TCanvas("c_allThetaBins_fullPhi_strip", "", canvasWidth, 380);
  c->Divide(nThetaBins, 1, 0.001, 0.001);

  for (int thetaBin = 0; thetaBin < nThetaBins; ++thetaBin) {
    c->cd(thetaBin + 1);
    gPad->SetMargin(0.16, 0.04, 0.20, 0.12);
    gPad->SetTicks(1, 1);

    TH1F* frame = gPad->DrawFrame(0.0, -0.08, 360.0, 0.08);
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("#phi [deg]");
    frame->GetYaxis()->SetTitle(Form("W_{peak} - %.6f", targetW));
    frame->GetXaxis()->SetTitleSize(0.070);
    frame->GetYaxis()->SetTitleSize(0.070);
    frame->GetXaxis()->SetLabelSize(0.060);
    frame->GetYaxis()->SetLabelSize(0.060);
    frame->GetXaxis()->SetTitleOffset(1.05);
    frame->GetYaxis()->SetTitleOffset(1.20);

    TLine zero(0.0, 0.0, 360.0, 0.0);
    zero.SetLineColor(kGray + 2);
    zero.SetLineStyle(2);
    zero.SetLineWidth(2);
    zero.Draw("SAME");

    auto* graph = graphs[thetaBin];
    if (graph && graph->GetN() > 0) {
      graph->SetMarkerStyle(20);
      graph->SetMarkerSize(0.80);
      graph->SetMarkerColor(kAzure + 2);
      graph->SetLineColor(kAzure + 2);
      graph->Draw("PE SAME");
    }

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.070);
    latex.DrawLatex(0.20, 0.86, Form("#theta %.1f-%.1f#circ", thetaRanges[thetaBin].first, thetaRanges[thetaBin].second));
  }

  const std::string out = outBase + "/allThetaBins_fullPhi_WPeak_vs_phi_strip";
  c->SaveAs((out + ".png").c_str());
  delete c;
}

}  // namespace

void analysisW(const std::string& csvPath = "../build/Elastic7546Corr/electron_w_afterCorr_merged_raw.csv",
               const std::string& outDir = "ElectronWCorrectionCSV",
               double targetW = kProtonMassGeV,
               int minEntries = 80,
               int nThreads = 400,
               double sector1PhiStartDeg = 335.0,
               int nPhiBinsPerSector = 6,
               const std::vector<double>& thetaBinEdges = {6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 13.0, 15.0, 25.0},
               double beamEnergyGeV = 7.546,
               int fitParameterPolyOrder = 2) {
  TStopwatch timer;
  timer.Start();

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gSystem->Exec(("mkdir -p " + outDir).c_str());

  const std::string cleanCsvPath = Trim(csvPath);
  const std::streamoff fileSize = CSVFileSize(cleanCsvPath);
  const int actualThreads = ResolveThreadCount(nThreads, fileSize);
  const auto csvRanges = MakeCSVRanges(fileSize, actualThreads);
  std::cout << "[analysisW] CSV: " << cleanCsvPath << std::endl;
  std::cout << "[analysisW] Using " << csvRanges.size() << " reader/fill thread(s)" << std::endl;
  std::cout << "[analysisW] Beam energy for p_true = " << beamEnergyGeV << " GeV" << std::endl;
  std::cout << "[analysisW] Sector 1 phi start = " << Wrap360(sector1PhiStartDeg) << " deg" << std::endl;
  if (thetaBinEdges.size() < 2) {
    throw std::runtime_error("analysisW: thetaBinEdges must contain at least two values.");
  }
  for (size_t i = 1; i < thetaBinEdges.size(); ++i) {
    if (thetaBinEdges[i] <= thetaBinEdges[i - 1]) {
      throw std::runtime_error("analysisW: thetaBinEdges must be strictly increasing.");
    }
  }
  const int nThetaBins = static_cast<int>(thetaBinEdges.size()) - 1;
  std::cout << "[analysisW] nPhiBinsPerSector = " << nPhiBinsPerSector
            << ", nThetaBins = " << nThetaBins << std::endl;

  const auto sectorPhiRanges = MakeFixedSectorPhiRanges(sector1PhiStartDeg);

  std::vector<HistPack> hists;
  const int displayPhiBins = std::max(2, 6 * nPhiBinsPerSector);
  const int fitPhiBins = std::max(1, nPhiBinsPerSector);
  std::vector<std::pair<double, double>> thetaRanges(nThetaBins);
  auto addHistPack = [&](const std::string& yVar,
                         const std::string& yTag,
                         const std::string& yTitle,
                         const std::string& peakLabel,
                         double yTarget,
                         double yMin,
                         double yMax,
                         int sector,
                         int thetaBin,
                         double thetaLo,
                         double thetaHi,
                         int nDisplayPhiBins,
                         int nFitPhiBins,
                         double phiMin,
                         double phiMax) {
    const std::string sectorLabel = sector > 0 ? Form("S%d", sector) : "all sectors";
    const std::string objectSectorTag = sector > 0 ? Form("S%d", sector) : "allSectors";
    const std::string displayName = Form("electron_csv_%s_thetaBin%02d_%s_vs_phi", objectSectorTag.c_str(), thetaBin, yTag.c_str());
    const std::string fitName = Form("electron_csv_%s_thetaBin%02d_%s_vs_phi_fitBins", objectSectorTag.c_str(), thetaBin, yTag.c_str());
    hists.push_back({yVar,
                     "phi",
                     Form("%s #theta %.1f-%.1f#circ", sectorLabel.c_str(), thetaLo, thetaHi),
                     sector,
                     -1,
                     thetaBin,
                     thetaLo,
                     thetaHi,
                     yTarget,
                     yTitle,
                     peakLabel,
                     new TH2D(displayName.c_str(),
                              "",
                              nDisplayPhiBins,
                              phiMin,
                              phiMax,
                              180,
                              yMin,
                              yMax),
                     new TH2D(fitName.c_str(),
                              "",
                              nFitPhiBins,
                              phiMin,
                              phiMax,
                              180,
                              yMin,
                              yMax)});
  };
  for (int thetaBin = 0; thetaBin < nThetaBins; ++thetaBin) {
    const double thetaLo = thetaBinEdges[thetaBin];
    const double thetaHi = thetaBinEdges[thetaBin + 1];
    thetaRanges[thetaBin] = {thetaLo, thetaHi};
    addHistPack("W", "W", "W [GeV]", "W_{peak}", targetW, 0.80, 1.10, 0, thetaBin, thetaLo, thetaHi, std::max(2, 6 * 6 * nPhiBinsPerSector), std::max(1, 6 * nPhiBinsPerSector), 0.0, 360.0);
    addHistPack("deltaP", "deltaP", "#Delta p = p_{true} - p_{exp} [GeV]", "#Delta p_{peak}", 0.0, -0.10, 0.10, 0, thetaBin, thetaLo, thetaHi, std::max(2, 6 * 6 * nPhiBinsPerSector), std::max(1, 6 * nPhiBinsPerSector), 0.0, 360.0);
  }
  for (int sector = 1; sector <= 6; ++sector) {
    const double phiStart = SectorPhiPlotStart(sectorPhiRanges[sector]);
    const double phiEnd = phiStart + sectorPhiRanges[sector].width;
    for (int thetaBin = 0; thetaBin < nThetaBins; ++thetaBin) {
      const double thetaLo = thetaBinEdges[thetaBin];
      const double thetaHi = thetaBinEdges[thetaBin + 1];
      addHistPack("W", "W", "W [GeV]", "W_{peak}", targetW, 0.80, 1.10, sector, thetaBin, thetaLo, thetaHi, displayPhiBins, fitPhiBins, phiStart, phiEnd);
      addHistPack("deltaP", "deltaP", "#Delta p = p_{true} - p_{exp} [GeV]", "#Delta p_{peak}", 0.0, -0.10, 0.10, sector, thetaBin, thetaLo, thetaHi, displayPhiBins, fitPhiBins, phiStart, phiEnd);
    }
  }

  PrepareHistPacks(hists);

  const size_t filledRows = FillHistPacksFromCSVParallel(cleanCsvPath, csvRanges, sectorPhiRanges, hists, beamEnergyGeV);
  std::cout << "[analysisW] Filled histograms from " << filledRows << " valid rows" << std::endl;
  if (filledRows == 0) {
    std::cerr << "[analysisW] No valid rows. Stop." << std::endl;
    DeleteHistPacks(hists);
    return;
  }

  const PeakConfig peakCfg{minEntries, 1.6, 0.08};
  std::vector<std::array<TGraphErrors*, 7>> summaryGraphsW(nThetaBins);
  std::vector<std::array<TGraphErrors*, 7>> summaryGraphsDeltaP(nThetaBins);
  std::array<TGraphErrors*, 7> deltaPA0Graphs = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  std::array<TGraphErrors*, 7> deltaPA1Graphs = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  for (auto& arr : summaryGraphsW) arr = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  for (auto& arr : summaryGraphsDeltaP) arr = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  for (int sector = 1; sector <= 6; ++sector) {
    deltaPA0Graphs[sector] = new TGraphErrors();
    deltaPA0Graphs[sector]->SetName(Form("deltaP_a_0_vs_theta_S%d", sector));
    deltaPA1Graphs[sector] = new TGraphErrors();
    deltaPA1Graphs[sector]->SetName(Form("deltaP_a_1_vs_theta_S%d", sector));
  }
  auto addGraphPoint = [](TGraphErrors* graph, double x, double y, double ex, double ey) {
    if (!graph) return;
    const int idx = graph->GetN();
    graph->SetPoint(idx, x, y);
    graph->SetPointError(idx, ex, ey);
  };

  for (auto& item : hists) {
    if (item.sector == 0) {
      DrawWMapOnly(item.displayHist, outDir, item.xVar, item.label, item.yTitle);
      delete item.displayHist;
      delete item.fitHist;
      continue;
    }

    const bool savePlots = true;
    const std::string sectorOutDir = outDir + Form("/Sector%d", item.sector);
    TGraphErrors* graph = DrawOneWMapAndPeaks(item.displayHist, item.fitHist, sectorOutDir, item.xVar, item.label, item.yTarget, item.yTitle, item.peakLabel, peakCfg, savePlots);
    if (item.thetaBin >= 0 && item.thetaBin < static_cast<int>(summaryGraphsW.size()) && item.sector >= 1 && item.sector <= 6) {
      if (item.yVar == "deltaP") {
        summaryGraphsDeltaP[item.thetaBin][item.sector] = graph;
      } else {
        summaryGraphsW[item.thetaBin][item.sector] = graph;
      }
    } else {
      delete graph;
    }
    delete item.displayHist;
    delete item.fitHist;
  }

  for (int thetaBin = 0; thetaBin < nThetaBins; ++thetaBin) {
    DrawThetaBinSectorSummary(summaryGraphsW[thetaBin],
                              sectorPhiRanges,
                              thetaBin,
                              thetaRanges[thetaBin].first,
                              thetaRanges[thetaBin].second,
                              outDir,
                              targetW,
                              "W_{peak}",
                              "W");
    std::array<TF1*, 7> deltaPFits = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    DrawThetaBinSectorSummary(summaryGraphsDeltaP[thetaBin],
                              sectorPhiRanges,
                              thetaBin,
                              thetaRanges[thetaBin].first,
                              thetaRanges[thetaBin].second,
                              outDir,
                              0.0,
                              "#Delta p_{peak}",
                              "deltaP",
                              &deltaPFits);
    const double thetaCenter = 0.5 * (thetaRanges[thetaBin].first + thetaRanges[thetaBin].second);
    const double thetaError = 0.5 * (thetaRanges[thetaBin].second - thetaRanges[thetaBin].first);
    for (int sector = 1; sector <= 6; ++sector) {
      if (!deltaPFits[sector]) continue;
      addGraphPoint(deltaPA0Graphs[sector],
                    thetaCenter,
                    deltaPFits[sector]->GetParameter(0),
                    thetaError,
                    deltaPFits[sector]->GetParError(0));
      addGraphPoint(deltaPA1Graphs[sector],
                    thetaCenter,
                    deltaPFits[sector]->GetParameter(1),
                    thetaError,
                    deltaPFits[sector]->GetParError(1));
      delete deltaPFits[sector];
    }
  }

  DrawFitParameterThetaSummary(deltaPA0Graphs, "a_0", "a_{0} [GeV]", outDir, fitParameterPolyOrder);
  DrawFitParameterThetaSummary(deltaPA1Graphs, "a_1", "a_{1} [GeV/deg]", outDir, fitParameterPolyOrder);

  for (int thetaBin = 0; thetaBin < nThetaBins; ++thetaBin) {
    for (auto* graph : summaryGraphsW[thetaBin]) delete graph;
    for (auto* graph : summaryGraphsDeltaP[thetaBin]) delete graph;
  }
  for (int sector = 1; sector <= 6; ++sector) {
    delete deltaPA0Graphs[sector];
    delete deltaPA1Graphs[sector];
  }

  timer.Stop();
  std::cout << "[analysisW] finished in " << timer.RealTime() << " s" << std::endl;
}
