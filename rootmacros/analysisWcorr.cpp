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
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <atomic>
#include <array>
#include <chrono>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
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
  double phi = 0.0;    // degree, wrapped to [0, 360)
  double theta = 0.0;  // degree
  int sector = -1;
};

struct SectorPhiRange {
  bool valid = false;
  double start = 0.0;
  double width = 60.0;
};

struct HistPack {
  std::string yVar;
  std::string yTag;
  std::string yTitle;
  std::string peakLabel;
  double yTarget = 0.0;
  int dataset = -1;
  int sector = 0;
  int thetaBin = -1;
  double thetaMin = 0.0;
  double thetaMax = 0.0;
  TH2D* fitHist = nullptr;
};

struct DatasetGraphs {
  std::string legend;
  int color = kBlack;
  int marker = 20;
  std::vector<std::array<TGraphErrors*, 7>> graphsW;
  std::vector<std::array<TGraphErrors*, 7>> graphsDeltaP;
};

struct DatasetConfig {
  std::string csvPath;
  std::string legend;
  int color = kBlack;
  int marker = 20;
};

struct PeakConfig {
  int minEntries = 80;
  double fitWindowSigma = 1.6;
  double maxPeakError = 0.08;
};

std::string Trim(const std::string& s) {
  const auto first = s.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return "";
  const auto last = s.find_last_not_of(" \t\r\n");
  return s.substr(first, last - first + 1);
}

std::string CleanObjectTag(std::string s) {
  for (auto& ch : s) {
    if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '_' && ch != '-') ch = '_';
  }
  return s.empty() ? "dataset" : s;
}

double Wrap360(double phiDeg) {
  while (phiDeg < 0.0) phiDeg += 360.0;
  while (phiDeg >= 360.0) phiDeg -= 360.0;
  return phiDeg;
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
    row.phi = Wrap360(std::stod(fields[2]));
    row.theta = std::stod(fields[3]);
    row.sector = std::stoi(fields[4]);
  } catch (...) {
    return false;
  }

  return std::isfinite(row.W) && std::isfinite(row.p) && std::isfinite(row.phi) &&
         std::isfinite(row.theta) && row.W > 0.0 && row.sector >= 1 && row.sector <= 6;
}

std::streamoff CSVFileSize(const std::string& csvPath) {
  std::ifstream in(csvPath, std::ios::binary | std::ios::ate);
  if (!in.is_open()) {
    throw std::runtime_error("analysisWcorr: cannot open CSV file: " + csvPath);
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

  std::cout << "\r[analysisWcorr] " << label << " [";
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
    throw std::runtime_error("analysisWcorr: cannot open CSV file: " + csvPath);
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

double SectorPhiPlotStart(const SectorPhiRange& range) {
  return range.start > 180.0 ? range.start - 360.0 : range.start;
}

std::array<SectorPhiRange, 7> MakeFixedSectorPhiRanges(double sector1PhiStartDeg) {
  std::array<SectorPhiRange, 7> ranges;
  for (int sector = 1; sector <= 6; ++sector) {
    ranges[sector].valid = true;
    ranges[sector].start = Wrap360(sector1PhiStartDeg + 60.0 * (sector - 1));
    ranges[sector].width = 60.0;

    const double plotStart = SectorPhiPlotStart(ranges[sector]);
    std::cout << "[analysisWcorr] sector " << sector
              << " fixed phi range: start=" << ranges[sector].start
              << " width=60 deg"
              << " plot=" << plotStart << " to " << plotStart + 60.0 << " deg" << std::endl;
  }
  return ranges;
}

double PhiOffsetFromRangeStart(double phiDeg, const SectorPhiRange& range) {
  double offset = Wrap360(phiDeg) - range.start;
  if (offset < 0.0) offset += 360.0;
  return offset;
}

double PhiForSectorAxis(double phiDeg, const SectorPhiRange& range) {
  const double plotStart = SectorPhiPlotStart(range);
  return plotStart + PhiOffsetFromRangeStart(phiDeg, range);
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

std::string ShiftAxisTitle(const std::string& peakLabel, double yTarget) {
  if (peakLabel.find("#Delta p") != std::string::npos && std::abs(yTarget) < 5e-7) return "#Delta p [GeV]";
  if (peakLabel.find("W") != std::string::npos && std::abs(yTarget - kProtonMassGeV) < 1e-4) return "W - m_{p} [GeV]";
  return Form("%s - %.6f [GeV]", peakLabel.c_str(), yTarget);
}

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

TGraphErrors* MakePeakGraph(TH2D* hist,
                            double yTarget,
                            const std::string& peakLabel,
                            const PeakConfig& cfg) {
  auto* g = new TGraphErrors();
  g->SetName((std::string(hist->GetName()) + "_peakGraph").c_str());

  int ip = 0;
  for (int ix = 1; ix <= hist->GetNbinsX(); ++ix) {
    if (ix == 1 || ix == hist->GetNbinsX()) continue;
    std::unique_ptr<TH1D> proj(hist->ProjectionY(Form("%s_py_%d", hist->GetName(), ix), ix, ix, "e"));
    if (!proj || proj->GetEntries() < cfg.minEntries) continue;

    const int maxBin = proj->GetMaximumBin();
    const double mu0 = proj->GetBinCenter(maxBin);
    const double rms = proj->GetRMS();
    if (rms <= 0.0) continue;

    const auto [fitLo, fitHi] = EstimatePeakFitRange(proj.get(), maxBin, rms, cfg.fitWindowSigma);
    if (fitHi <= fitLo) continue;

    TF1 fit = MakeAsymGaussianFit(Form("%s_asymGaus_%d", hist->GetName(), ix), fitLo, fitHi);
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
    if (!std::isfinite(peak) || !std::isfinite(peakErr) || peakErr <= 0.0 || peakErr > cfg.maxPeakError) continue;

    const double x = hist->GetXaxis()->GetBinCenter(ix);
    const double ex = 0.5 * hist->GetXaxis()->GetBinWidth(ix);
    const double deltaPeak = peak - yTarget;
    g->SetPoint(ip, x, deltaPeak);
    g->SetPointError(ip, ex, peakErr);
    ++ip;
  }

  g->SetTitle(Form(";%s;%s", "electron #phi [deg]", ShiftAxisTitle(peakLabel, yTarget).c_str()));
  return g;
}

void StyleGraph(TGraphErrors* g, int color, int marker) {
  if (!g) return;
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(0.85);
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->SetLineWidth(2);
}

void PrepareHistPacks(std::vector<HistPack>& hists) {
  for (auto& item : hists) {
    item.fitHist->SetDirectory(nullptr);
  }
}

std::vector<HistPack> CloneHistPacks(const std::vector<HistPack>& hists, const std::string& suffix) {
  std::vector<HistPack> clones;
  clones.reserve(hists.size());
  for (const auto& item : hists) {
    auto* fit = dynamic_cast<TH2D*>(item.fitHist->Clone((std::string(item.fitHist->GetName()) + suffix).c_str()));
    fit->Reset("ICES");
    fit->SetDirectory(nullptr);
    clones.push_back({item.yVar,
                      item.yTag,
                      item.yTitle,
                      item.peakLabel,
                      item.yTarget,
                      item.dataset,
                      item.sector,
                      item.thetaBin,
                      item.thetaMin,
                      item.thetaMax,
                      fit});
  }
  return clones;
}

void DeleteHistPacks(std::vector<HistPack>& hists) {
  for (auto& item : hists) {
    delete item.fitHist;
    item.fitHist = nullptr;
  }
}

bool FillHistPacks(const ElectronWRow& row,
                   int dataset,
                   const std::array<SectorPhiRange, 7>& sectorPhiRanges,
                   std::vector<HistPack>& hists,
                   double beamEnergyGeV) {
  const auto& phiRange = sectorPhiRanges[row.sector];
  if (!phiRange.valid) return false;
  const double phiAxis = PhiForSectorAxis(row.phi, phiRange);
  const double phiPlotStart = SectorPhiPlotStart(phiRange);
  if (phiAxis < phiPlotStart || phiAxis >= phiPlotStart + phiRange.width) return false;

  for (auto& item : hists) {
    if (item.dataset != dataset) continue;
    if (item.sector != row.sector) continue;
    if (item.thetaBin >= 0 && (row.theta < item.thetaMin || row.theta >= item.thetaMax)) continue;
    const double yValue = item.yVar == "deltaP" ? DeltaPTrueMinusExp(row, beamEnergyGeV) : row.W;
    if (!std::isfinite(yValue)) continue;
    item.fitHist->Fill(phiAxis, yValue);
  }
  return true;
}

size_t FillHistPacksFromCSVParallel(const std::string& csvPath,
                                    int dataset,
                                    const std::vector<std::pair<std::streamoff, std::streamoff>>& ranges,
                                    const std::array<SectorPhiRange, 7>& sectorPhiRanges,
                                    std::vector<HistPack>& hists,
                                    double beamEnergyGeV,
                                    const std::string& progressLabel) {
  std::vector<std::vector<HistPack>> localHists;
  localHists.reserve(ranges.size());
  for (size_t i = 0; i < ranges.size(); ++i) {
    localHists.push_back(CloneHistPacks(hists, Form("_d%d_thread%zu", dataset, i)));
  }

  std::vector<size_t> accepted(ranges.size(), 0);
  std::vector<std::thread> workers;
  workers.reserve(ranges.size());

  long long totalBytes = 0;
  for (const auto& range : ranges) totalBytes += static_cast<long long>(range.second - range.first);
  std::atomic<long long> doneBytes{0};
  std::atomic<bool> finished{false};
  std::thread monitor(ProgressMonitor, progressLabel, std::ref(doneBytes), totalBytes, std::ref(finished));

  for (size_t i = 0; i < ranges.size(); ++i) {
    workers.emplace_back([&, i]() {
      accepted[i] = ForEachCSVRowInRange(csvPath, ranges[i].first, ranges[i].second, [&](const ElectronWRow& row) {
        FillHistPacks(row, dataset, sectorPhiRanges, localHists[i], beamEnergyGeV);
      }, &doneBytes);
    });
  }
  for (auto& worker : workers) worker.join();
  finished.store(true);
  monitor.join();

  for (auto& item : hists) {
    if (item.dataset == dataset) item.fitHist->Reset("ICES");
  }

  for (auto& oneThreadHists : localHists) {
    for (size_t ih = 0; ih < hists.size(); ++ih) {
      if (hists[ih].dataset == dataset) hists[ih].fitHist->Add(oneThreadHists[ih].fitHist);
    }
    DeleteHistPacks(oneThreadHists);
  }

  size_t totalAccepted = 0;
  for (const auto count : accepted) totalAccepted += count;
  return totalAccepted;
}

void DrawThetaBinSectorSummary(const std::vector<DatasetGraphs>& datasets,
                               const std::array<SectorPhiRange, 7>& sectorPhiRanges,
                               int thetaBin,
                               double thetaMin,
                               double thetaMax,
                               const std::string& outBase,
                               double yTarget,
                               const std::string& peakLabel,
                               const std::string& quantityTag) {
  bool hasAny = false;
  for (const auto& dataset : datasets) {
    const auto& graphs = quantityTag == "W" ? dataset.graphsW : dataset.graphsDeltaP;
    for (int sector = 1; sector <= 6; ++sector) {
      if (graphs[thetaBin][sector] && graphs[thetaBin][sector]->GetN() > 0) {
        hasAny = true;
        break;
      }
    }
  }
  if (!hasAny) return;

  auto* c = new TCanvas(Form("c_%s_thetaBin%02d_allSectors_phi", quantityTag.c_str(), thetaBin), "", 1800, 1100);
  c->Divide(3, 2);

  for (int sector = 1; sector <= 6; ++sector) {
    c->cd(sector);
    gPad->SetMargin(0.16, 0.06, 0.15, 0.10);
    gPad->SetTicks(1, 1);

    const double phiStart = SectorPhiPlotStart(sectorPhiRanges[sector]);
    const double phiEnd = phiStart + sectorPhiRanges[sector].width;
    const double yMin = quantityTag == "W" ? -0.08 : -0.08;
    const double yMax = quantityTag == "W" ? 0.08 : 0.08;
    TH1F* frame = gPad->DrawFrame(phiStart, yMin, phiEnd, yMax);
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("electron #phi [deg]");
    frame->GetYaxis()->SetTitle(ShiftAxisTitle(peakLabel, yTarget).c_str());
    frame->GetXaxis()->SetTitleSize(0.055);
    frame->GetYaxis()->SetTitleSize(0.055);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetXaxis()->SetTitleOffset(1.05);
    frame->GetYaxis()->SetTitleOffset(1.25);

    auto* legend = new TLegend(0.20, 0.70, 0.62, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.040);

    for (const auto& dataset : datasets) {
      const auto& graphs = quantityTag == "W" ? dataset.graphsW : dataset.graphsDeltaP;
      TGraphErrors* graph = graphs[thetaBin][sector];
      if (!graph || graph->GetN() <= 0) continue;
      graph->Draw("PE SAME");
      legend->AddEntry(graph, dataset.legend.c_str(), "p");
    }
    auto* zero = new TLine(phiStart, 0.0, phiEnd, 0.0);
    zero->SetLineColor(kRed + 1);
    zero->SetLineStyle(2);
    zero->SetLineWidth(2);
    zero->SetBit(TObject::kCanDelete);
    zero->Draw("SAME");
    legend->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextAlign(31);
    latex.SetTextSize(0.055);
    latex.DrawLatex(0.92, 0.86, Form("S%d", sector));
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.92, 0.79, Form("#theta %.1f-%.1f#circ", thetaMin, thetaMax));
  }

  const std::string out = outBase + Form("/thetaBin%02d_allSectors_%sPeak_vs_phi.png", thetaBin, quantityTag.c_str());
  c->SaveAs(out.c_str());
  delete c;
}

void DeleteDatasetGraphs(DatasetGraphs& dataset) {
  for (auto& thetaArr : dataset.graphsW) {
    for (int sector = 1; sector <= 6; ++sector) {
      delete thetaArr[sector];
      thetaArr[sector] = nullptr;
    }
  }
  for (auto& thetaArr : dataset.graphsDeltaP) {
    for (int sector = 1; sector <= 6; ++sector) {
      delete thetaArr[sector];
      thetaArr[sector] = nullptr;
    }
  }
}

}  // namespace

void analysisWcorr(const std::string& csvPath1 = "../build/Elastic6535/electron_w_afterCorr_merged_raw.csv",
                   const std::string& csvPath2 = "../build/Elastic6535Corr/electron_w_afterCorr_merged_raw.csv",
                   const std::string& outDir = "ElectronWCorrectionCompare",
                   const std::string& legend1 = "before",
                   const std::string& legend2 = "after",
                   double beamEnergyGeV = 6.535,
                   int minEntries = 80,
                   int nThreads = 800,
                   double sector1PhiStartDeg = 335.0,
                   int nPhiBinsPerSector = 6,
                   const std::vector<double>& thetaBinEdges = {6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 13.0, 15.0, 25.0}) {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gSystem->Exec(("mkdir -p " + outDir).c_str());

  if (thetaBinEdges.size() < 2) {
    throw std::runtime_error("analysisWcorr: thetaBinEdges must contain at least two values.");
  }
  for (size_t i = 1; i < thetaBinEdges.size(); ++i) {
    if (thetaBinEdges[i] <= thetaBinEdges[i - 1]) {
      throw std::runtime_error("analysisWcorr: thetaBinEdges must be strictly increasing.");
    }
  }

  const std::vector<DatasetConfig> configs = {
      {Trim(csvPath1), legend1, kBlack, 20},
      {Trim(csvPath2), legend2, kRed + 1, 24},
  };
  const int nThetaBins = static_cast<int>(thetaBinEdges.size()) - 1;
  const auto sectorPhiRanges = MakeFixedSectorPhiRanges(sector1PhiStartDeg);
  const PeakConfig peakCfg{minEntries, 1.6, 0.08};

  std::cout << "[analysisWcorr] Beam energy for p_true = " << beamEnergyGeV << " GeV" << std::endl;
  std::cout << "[analysisWcorr] nPhiBinsPerSector = " << nPhiBinsPerSector
            << ", nThetaBins = " << nThetaBins
            << ", minEntries = " << minEntries << std::endl;

  std::vector<HistPack> hists;
  hists.reserve(configs.size() * nThetaBins * 6 * 2);
  std::vector<DatasetGraphs> datasets(configs.size());
  for (size_t id = 0; id < configs.size(); ++id) {
    datasets[id].legend = configs[id].legend;
    datasets[id].color = configs[id].color;
    datasets[id].marker = configs[id].marker;
    datasets[id].graphsW.resize(nThetaBins);
    datasets[id].graphsDeltaP.resize(nThetaBins);
    for (auto& arr : datasets[id].graphsW) arr = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    for (auto& arr : datasets[id].graphsDeltaP) arr = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};

    const std::string tag = CleanObjectTag(configs[id].legend);
    for (int sector = 1; sector <= 6; ++sector) {
      const double phiStart = SectorPhiPlotStart(sectorPhiRanges[sector]);
      const double phiEnd = phiStart + sectorPhiRanges[sector].width;
      for (int thetaBin = 0; thetaBin < nThetaBins; ++thetaBin) {
        const double thetaLo = thetaBinEdges[thetaBin];
        const double thetaHi = thetaBinEdges[thetaBin + 1];
        hists.push_back({"W",
                         "W",
                         "W [GeV]",
                         "W_{peak}",
                         kProtonMassGeV,
                         static_cast<int>(id),
                         sector,
                         thetaBin,
                         thetaLo,
                         thetaHi,
                         new TH2D(Form("h_%s_S%d_thetaBin%02d_W_vs_phi_fitBins", tag.c_str(), sector, thetaBin),
                                  "",
                                  std::max(1, nPhiBinsPerSector),
                                  phiStart,
                                  phiEnd,
                                  180,
                                  0.80,
                                  1.10)});
        hists.push_back({"deltaP",
                         "deltaP",
                         "#Delta p = p_{true} - p_{exp} [GeV]",
                         "#Delta p_{peak}",
                         0.0,
                         static_cast<int>(id),
                         sector,
                         thetaBin,
                         thetaLo,
                         thetaHi,
                         new TH2D(Form("h_%s_S%d_thetaBin%02d_deltaP_vs_phi_fitBins", tag.c_str(), sector, thetaBin),
                                  "",
                                  std::max(1, nPhiBinsPerSector),
                                  phiStart,
                                  phiEnd,
                                  180,
                                  -0.10,
                                  0.10)});
      }
    }
  }
  PrepareHistPacks(hists);

  for (size_t id = 0; id < configs.size(); ++id) {
    const std::streamoff fileSize = CSVFileSize(configs[id].csvPath);
    const int actualThreads = ResolveThreadCount(nThreads, fileSize);
    const auto ranges = MakeCSVRanges(fileSize, actualThreads);
    std::cout << "[analysisWcorr] CSV " << id + 1 << ": " << configs[id].csvPath << std::endl;
    std::cout << "[analysisWcorr] Using " << ranges.size() << " reader/fill thread(s) for " << configs[id].legend << std::endl;
    const size_t filledRows = FillHistPacksFromCSVParallel(configs[id].csvPath,
                                                           static_cast<int>(id),
                                                           ranges,
                                                           sectorPhiRanges,
                                                           hists,
                                                           beamEnergyGeV,
                                                           "fill " + configs[id].legend);
    std::cout << "[analysisWcorr] Filled histograms from " << filledRows
              << " valid rows for " << configs[id].legend << std::endl;
  }

  for (auto& item : hists) {
    TGraphErrors* graph = MakePeakGraph(item.fitHist, item.yTarget, item.peakLabel, peakCfg);
    StyleGraph(graph, datasets[item.dataset].color, datasets[item.dataset].marker);
    if (item.yVar == "W") {
      datasets[item.dataset].graphsW[item.thetaBin][item.sector] = graph;
    } else {
      datasets[item.dataset].graphsDeltaP[item.thetaBin][item.sector] = graph;
    }
  }

  for (int thetaBin = 0; thetaBin < nThetaBins; ++thetaBin) {
    DrawThetaBinSectorSummary(datasets,
                              sectorPhiRanges,
                              thetaBin,
                              thetaBinEdges[thetaBin],
                              thetaBinEdges[thetaBin + 1],
                              outDir,
                              kProtonMassGeV,
                              "W_{peak}",
                              "W");
    DrawThetaBinSectorSummary(datasets,
                              sectorPhiRanges,
                              thetaBin,
                              thetaBinEdges[thetaBin],
                              thetaBinEdges[thetaBin + 1],
                              outDir,
                              0.0,
                              "#Delta p_{peak}",
                              "deltaP");
  }

  DeleteHistPacks(hists);
  for (auto& dataset : datasets) DeleteDatasetGraphs(dataset);
  std::cout << "[analysisWcorr] Saved theta-bin comparison plots to " << outDir << std::endl;
}
