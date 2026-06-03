#include "DVCSAnalysis.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <stdexcept>

#include "AnalysisTaskManager.h"
#include "PerRunCounter.h"
#include "ROOT/RVec.hxx"

namespace {
constexpr double kElectronMassGeV = 0.00051099895;
constexpr double kProtonMassGeV = 0.9382720813;

static inline int SelectedElectronIndex(const std::vector<int>& pid, const std::vector<int>& pass) {
  const auto n = std::min(pid.size(), pass.size());
  for (size_t i = 0; i < n; ++i) {
    if (pid[i] == 11 && pass[i]) return static_cast<int>(i);
  }
  return -1;
}

static inline float FloatAtIndex(const std::vector<float>& values, int idx) {
  return (idx >= 0 && static_cast<size_t>(idx) < values.size()) ? values[idx] : -999.0f;
}

static inline int TrackSectorForParticle(int particleIndex,
                                         const std::vector<int16_t>& trackPindex,
                                         const std::vector<int16_t>& trackSector) {
  if (particleIndex < 0) return -1;
  for (size_t i = 0; i < trackPindex.size(); ++i) {
    if (i >= trackSector.size()) continue;
    if (trackPindex[i] != particleIndex) continue;
    const int sector = static_cast<int>(trackSector[i]);
    if (sector >= 1 && sector <= 6) return sector;
  }
  return -1;
}

static inline float InclusiveWFromElectron(float beamEnergyGeV, float electronPGeV, float electronThetaRad) {
  if (beamEnergyGeV <= 0.0f || electronPGeV <= 0.0f) return -999.0f;
  const double scatteredE = std::sqrt(electronPGeV * electronPGeV + kElectronMassGeV * kElectronMassGeV);
  const double sinHalfTheta = std::sin(0.5 * electronThetaRad);
  const double q2 = 4.0 * beamEnergyGeV * scatteredE * sinHalfTheta * sinHalfTheta;
  const double nu = beamEnergyGeV - scatteredE;
  const double w2 = kProtonMassGeV * kProtonMassGeV + 2.0 * kProtonMassGeV * nu - q2;
  return w2 > 0.0 ? static_cast<float>(std::sqrt(w2)) : -999.0f;
}

static inline float PhiRadToWrappedDeg(float phiRad) {
  float phiDeg = phiRad * 180.0f / static_cast<float>(M_PI);
  while (phiDeg < 0.0f) phiDeg += 360.0f;
  while (phiDeg >= 360.0f) phiDeg -= 360.0f;
  return phiDeg;
}
}  // namespace

static inline std::pair<std::string, std::string> PickRunEventCols(ROOT::RDF::RNode df) {
  auto cols = df.GetColumnNames();
  auto has = [&](const std::string& n) { return std::find(cols.begin(), cols.end(), n) != cols.end(); };

  if (has("RUN_config_run") && has("RUN_config_event")) return {"RUN_config_run", "RUN_config_event"};

  if (has("RUN::config.run") && has("RUN::config.event")) return {"RUN::config.run", "RUN::config.event"};

  throw std::runtime_error("QADB: cannot find run/event columns");
}
#include "ROOT/RVec.hxx"

// helper: detect substring
static inline bool has_substr(const std::string& s, const std::string& sub) {
  return s.find(sub) != std::string::npos;
}

static inline ROOT::RDF::RNode
DefineRunEventScalars(ROOT::RDF::RNode df, const std::string& runCol, const std::string& evCol) {
  const auto runType = df.GetColumnType(runCol);
  const auto evType  = df.GetColumnType(evCol);

  // ----- scalar cases -----
  auto is_scalar_int = [&](const std::string& t){
    return t=="int" || t=="Int_t" || t=="unsigned int" || t=="UInt_t";
  };
  auto is_scalar_long = [&](const std::string& t){
    return t=="Long64_t" || t=="long" || t=="long long" || t=="ULong64_t" || t=="unsigned long long";
  };
  auto is_scalar_short = [&](const std::string& t){
    return t=="short" || t=="Short_t" || t=="unsigned short" || t=="UShort_t";
  };

  if (is_scalar_int(runType) && is_scalar_int(evType)) {
    return df.Define("RUN_run",  [](int r){ return r; }, {runCol})
             .Define("RUN_event",[](int e){ return e; }, {evCol});
  }
  if (is_scalar_long(runType) && is_scalar_long(evType)) {
    return df.Define("RUN_run",  [](Long64_t r){ return (int)r; }, {runCol})
             .Define("RUN_event",[](Long64_t e){ return (int)e; }, {evCol});
  }
  if (is_scalar_short(runType) && is_scalar_short(evType)) {
    return df.Define("RUN_run",  [](Short_t r){ return (int)r; }, {runCol})
             .Define("RUN_event",[](Short_t e){ return (int)e; }, {evCol});
  }

  // ----- RVec cases -----
  auto is_rvec = [&](const std::string& t){ return has_substr(t, "ROOT::VecOps::RVec"); };

  // RVec<int> / RVec<Int_t> / etc (match by substring)
  if (is_rvec(runType) && is_rvec(evType) && (has_substr(runType,"<int") || has_substr(runType,"<Int_t"))) {
    return df.Define("RUN_run",  [](const ROOT::VecOps::RVec<int>& v){ return v.empty() ? -1 : (int)v[0]; }, {runCol})
             .Define("RUN_event",[](const ROOT::VecOps::RVec<int>& v){ return v.empty() ? -1 : (int)v[0]; }, {evCol});
  }
  if (is_rvec(runType) && is_rvec(evType) && (has_substr(runType,"<Long64_t") || has_substr(runType,"<long long") || has_substr(runType,"<long"))) {
    return df.Define("RUN_run",  [](const ROOT::VecOps::RVec<Long64_t>& v){ return v.empty() ? -1 : (int)v[0]; }, {runCol})
             .Define("RUN_event",[](const ROOT::VecOps::RVec<Long64_t>& v){ return v.empty() ? -1 : (int)v[0]; }, {evCol});
  }
  if (is_rvec(runType) && is_rvec(evType) && (has_substr(runType,"<Short_t") || has_substr(runType,"<short"))) {
    return df.Define("RUN_run",  [](const ROOT::VecOps::RVec<Short_t>& v){ return v.empty() ? -1 : (int)v[0]; }, {runCol})
             .Define("RUN_event",[](const ROOT::VecOps::RVec<Short_t>& v){ return v.empty() ? -1 : (int)v[0]; }, {evCol});
  }

  // ----- std::vector cases -----
  auto is_stdvec = [&](const std::string& t){ return has_substr(t, "std::vector") || has_substr(t, "vector<"); };

  if (is_stdvec(runType) && is_stdvec(evType) && (has_substr(runType,"<int") || has_substr(runType,"<Int_t"))) {
    return df.Define("RUN_run",  [](const std::vector<int>& v){ return v.empty() ? -1 : (int)v[0]; }, {runCol})
             .Define("RUN_event",[](const std::vector<int>& v){ return v.empty() ? -1 : (int)v[0]; }, {evCol});
  }
  if (is_stdvec(runType) && is_stdvec(evType) && (has_substr(runType,"<Long64_t") || has_substr(runType,"<long long") || has_substr(runType,"<long"))) {
    return df.Define("RUN_run",  [](const std::vector<Long64_t>& v){ return v.empty() ? -1 : (int)v[0]; }, {runCol})
             .Define("RUN_event",[](const std::vector<Long64_t>& v){ return v.empty() ? -1 : (int)v[0]; }, {evCol});
  }
  if (is_stdvec(runType) && is_stdvec(evType) && (has_substr(runType,"<Short_t") || has_substr(runType,"<short"))) {
    return df.Define("RUN_run",  [](const std::vector<Short_t>& v){ return v.empty() ? -1 : (int)v[0]; }, {runCol})
             .Define("RUN_event",[](const std::vector<Short_t>& v){ return v.empty() ? -1 : (int)v[0]; }, {evCol});
  }

  throw std::runtime_error("QADB: unsupported run/event types: run=" + runType + " ev=" + evType);
}

DVCSAnalysis::DVCSAnalysis(bool IsMC, bool IsReproc, bool IsMinBook) : IsMC(IsMC), IsReproc(IsReproc), IsMinBooking(IsMinBook), fHistPhotonP(nullptr) {}
DVCSAnalysis::~DVCSAnalysis() {}

void DVCSAnalysis::UserCreateOutputObjects() {}

void DVCSAnalysis::UserExec(ROOT::RDF::RNode& df) {
  using namespace std;

  if (fMaxEvents > 0) {
    df = df.Range(0, fMaxEvents);  // only process the first fMaxEvents
  }
  if (!fTrackCuts || !fEventCuts) throw std::runtime_error("DVCSAnalysis: One or more cut not set.");

  fTrackCutsNoFid = std::make_shared<TrackCut>(*fTrackCuts);
  fTrackCutsWithFid = std::make_shared<TrackCut>(*fTrackCuts);
  fTrackCutsWithFid->SetDoFiducialCut(true);
  fTrackCutsWithFid->SetFiducialCutOptions(true, true);  // apply both DC and ECAL cuts

  // Cache column names
  // auto colnames = df.GetColumnNames();
  auto dfDefs = df;
  dfDefs = DefineOrRedefine(dfDefs, "REC_Particle_num", [](const std::vector<int>& pid) { return static_cast<int>(pid.size()); }, {"REC_Particle_pid"});
  dfDefs = DefineOrRedefine(dfDefs, "REC_Particle_theta", RECParticletheta(), RECParticle::All());
  dfDefs = DefineOrRedefine(dfDefs, "REC_Particle_phi", RECParticlephi(), RECParticle::All());
  dfDefs = DefineOrRedefine(dfDefs, "REC_Particle_p", RECParticleP(), RECParticle::All());
  dforginal = dfDefs;
  // Fiducial cuts
  auto dfDefsWithTraj = dfDefs;
  auto [runCol, evCol] = PickRunEventCols(dfDefsWithTraj);
  dfDefsWithTraj = DefineRunEventScalars(dfDefsWithTraj, runCol, evCol);
  // QADB cuts should be place in the first to reduce the computation load
  if (fIsQADBCut && fQADBCuts) {
    std::cout << "Applying QADB cut..." << std::endl;
    auto qadb = *fQADBCuts;
    dfDefsWithTraj =
        dfDefsWithTraj.Define("REC_QADB_pass", [qadb](int run, int ev) mutable { return qadb(run, ev); }, {"RUN_run", "RUN_event"}).Filter("REC_QADB_pass", "QADB pass");
  }

  auto trajCols = CombineColumns(RECTraj::ForFiducialCut(), std::vector<std::string>{"REC_Particle_pid"}, std::vector<std::string>{"REC_Particle_num"});
  auto caloCols =
      CombineColumns(RECCalorimeter::ForFiducialCut(), std::vector<std::string>{"REC_Particle_pid"}, std::vector<std::string>{"REC_Particle_p"}, std::vector<std::string>{"REC_Particle_num"});
  auto fwdtagCols = CombineColumns(RECForwardTagger::ForFiducialCut(), std::vector<std::string>{"REC_Particle_pid"}, std::vector<std::string>{"REC_Particle_num"});

  dfDefsWithTraj = DefineOrRedefine(dfDefsWithTraj, "REC_Track_pass_nofid", fTrackCutsNoFid->RECTrajPass(), trajCols);
  dfDefsWithTraj = DefineOrRedefine(dfDefsWithTraj, "REC_Traj_pass_fid", fTrackCutsWithFid->RECTrajPass(), trajCols);
  dfDefsWithTraj = DefineOrRedefine(dfDefsWithTraj, "REC_Calorimeter_pass_fid", fTrackCutsWithFid->RECCalorimeterPass(), caloCols);

  if (fFTonConfig) {
    dfDefsWithTraj = DefineOrRedefine(dfDefsWithTraj, "REC_ForwardTagger_pass_fid", fTrackCutsWithFid->RECForwardTaggerPass(), fwdtagCols);
  }
  dfDefsWithTraj = DefineOrRedefine(dfDefsWithTraj, "REC_Track_pass_fid", Columns::LogicalAND2(),
                                    CombineColumns(std::vector<std::string>{"REC_Traj_pass_fid"}, std::vector<std::string>{"REC_Calorimeter_pass_fid"}));
  if (fFTonConfig) {
    dfDefsWithTraj = DefineOrRedefine(dfDefsWithTraj, "REC_Track_pass_fid", Columns::LogicalAND2(),
                                      CombineColumns(std::vector<std::string>{"REC_Track_pass_fid"}, std::vector<std::string>{"REC_ForwardTagger_pass_fid"}));
  }
  auto AllCols = CombineColumns(trajCols, caloCols);

  auto cols_track_fid = CombineColumns(RECParticle::All(), std::vector<std::string>{"REC_Track_pass_fid"});
  auto cols_track_nofid = CombineColumns(RECParticle::All(), std::vector<std::string>{"REC_Track_pass_nofid"});

  if (fAcceptAll) {
    fEventCuts->AcceptEverything(true);
  }

  dfSelected = dfDefsWithTraj;
  dfSelected = DefineOrRedefine(*dfSelected, "EventCutResult", *fEventCuts, cols_track_nofid);
  dfSelected = DefineOrRedefine(*dfSelected, "REC_Event_pass", [](const EventCutResult& result) { return result.eventPass; }, {"EventCutResult"});
  dfSelected = DefineOrRedefine(*dfSelected, "REC_Particle_pass", [](const EventCutResult& result) { return result.particlePass; }, {"EventCutResult"});
  dfSelected = DefineOrRedefine(*dfSelected, "REC_Photon_MaxE", [](const EventCutResult& result) { return result.MaxPhotonEnergyPass; }, {"EventCutResult"});

  if (fDoInvMassCut) {
    fEventCuts->SetDoCutMotherInvMass(true);
    dfSelected = DefineOrRedefine(*dfSelected, "REC_DaughterParticle_pass", [](const EventCutResult& result) { return result.particleDaughterPass; }, {"EventCutResult"});
    dfSelected = DefineOrRedefine(*dfSelected, "REC_MotherMass", [](const EventCutResult& result) { return result.MotherMass; }, {"EventCutResult"});
  }
  dfSelected = dfSelected->Filter("REC_Event_pass");

  // After fiducial cut
  if (fFiducialCut) {
    dfSelected_afterFid = dfDefsWithTraj;
    dfSelected_afterFid = DefineOrRedefine(*dfSelected_afterFid, "EventCutResult", *fEventCuts, cols_track_fid);
    dfSelected_afterFid = DefineOrRedefine(*dfSelected_afterFid, "REC_Event_pass", [](const EventCutResult& result) { return result.eventPass; }, {"EventCutResult"});
    dfSelected_afterFid = DefineOrRedefine(*dfSelected_afterFid, "REC_Particle_pass", [](const EventCutResult& result) { return result.particlePass; }, {"EventCutResult"});
    dfSelected_afterFid = DefineOrRedefine(*dfSelected_afterFid, "REC_Photon_MaxE", [](const EventCutResult& result) { return result.MaxPhotonEnergyPass; }, {"EventCutResult"});
    if (fDoInvMassCut) {
      fEventCuts->SetDoCutMotherInvMass(true);
      dfSelected_afterFid =
          DefineOrRedefine(*dfSelected_afterFid, "REC_DaughterParticle_pass", [](const EventCutResult& result) { return result.particleDaughterPass; }, {"EventCutResult"});
      dfSelected_afterFid = DefineOrRedefine(*dfSelected_afterFid, "REC_MotherMass", [](const EventCutResult& result) { return result.MotherMass; }, {"EventCutResult"});
    }
    dfSelected_afterFid = dfSelected_afterFid->Filter("REC_Event_pass");
  }

  dfSelected_afterFid_afterCorr = dfSelected_afterFid;

  if (fMomCorr && fDoMomentumCorrection) {
    std::cout << "Applying momentum correction..." << std::endl;
    dfSelected_afterFid_afterCorr = DefineOrRedefine(*dfSelected_afterFid_afterCorr, "REC_Particle_px", fMomCorr->RECParticlePxCorrected(), RECParticle::Extend());
    dfSelected_afterFid_afterCorr = DefineOrRedefine(*dfSelected_afterFid_afterCorr, "REC_Particle_py", fMomCorr->RECParticlePyCorrected(), RECParticle::Extend());
    dfSelected_afterFid_afterCorr = DefineOrRedefine(*dfSelected_afterFid_afterCorr, "REC_Particle_pz", fMomCorr->RECParticlePzCorrected(), RECParticle::Extend());
    dfSelected_afterFid_afterCorr = DefineOrRedefine(*dfSelected_afterFid_afterCorr, "REC_Particle_theta", RECParticletheta(), RECParticle::All());
    dfSelected_afterFid_afterCorr = DefineOrRedefine(*dfSelected_afterFid_afterCorr, "REC_Particle_phi", RECParticlephi(), RECParticle::All());
    dfSelected_afterFid_afterCorr = DefineOrRedefine(*dfSelected_afterFid_afterCorr, "REC_Particle_p", RECParticleP(), RECParticle::All());
  }
}

std::vector<std::string> DVCSAnalysis::MinimalColumns() const {
  using V = std::vector<std::string>;

  auto cols = CombineColumns(RECParticle::All());
  cols = CombineColumns(cols, RECTraj::ForFiducialCut());
  cols = CombineColumns(cols, RECCalorimeter::ForFiducialCut());
  cols = CombineColumns(cols, RECTrack::ForFiducialCut());
  if (fFTonConfig)
    cols = CombineColumns(cols, RECForwardTagger::ForFiducialCut());
    cols.push_back("REC_Event_helicity");

  for (const auto& c : V{"RUN_config_run", "RUN_config_event",
                            "RUN::config.run", "RUN::config.event",
                            "RUN_run", "RUN_event",
                            "REC_Photon_MaxE"})
    cols.push_back(c);

  // Pre-computed kinematics and analysis decision columns.
  for (const auto& c : V{
      "num_events",
      "REC_Particle_num",
      "REC_Particle_theta",
      "REC_Particle_phi",
      "REC_Particle_p",
      "REC_Particle_pass",
      "REC_Event_pass",
      "REC_MotherMass",
      "REC_DaughterParticle_pass"
  })
    cols.push_back(c);
  
  if (IsMC){
    for (const auto& c : V{"MC_Particle_pid", "MC_Particle_px", "MC_Particle_py", "MC_Particle_pz", "MC_Particle_vx", "MC_Particle_vy", "MC_Particle_vz", "MC_Particle_vt", "MC_Event_weight",
         "MC_Event_pbeam",  // include if this exists
         "MC_Event_ptarget", "MC_Event_ebeam"})
    cols.push_back(c);
  }
  return cols;
}

void DVCSAnalysis::WriteFinalElectronWCSV(ROOT::RDF::RNode df, const std::string& csvPath) {
  auto out = df;
  out = DefineOrRedefine(out, "WCSV_e_idx", SelectedElectronIndex, {"REC_Particle_pid", "REC_Particle_pass"});
  out = DefineOrRedefine(out, "WCSV_e_p", FloatAtIndex, {"REC_Particle_p", "WCSV_e_idx"});
  out = DefineOrRedefine(out, "WCSV_e_theta_rad", FloatAtIndex, {"REC_Particle_theta", "WCSV_e_idx"});
  out = DefineOrRedefine(out, "WCSV_e_phi_rad", FloatAtIndex, {"REC_Particle_phi", "WCSV_e_idx"});
  out = DefineOrRedefine(out, "WCSV_e_sector", TrackSectorForParticle, {"WCSV_e_idx", "REC_Track_pindex", "REC_Track_sector"});
  out = DefineOrRedefine(out, "WCSV_W", [beam = fbeam_energy](float p, float theta) { return InclusiveWFromElectron(beam, p, theta); }, {"WCSV_e_p", "WCSV_e_theta_rad"});
  out = DefineOrRedefine(out, "WCSV_e_phi", PhiRadToWrappedDeg, {"WCSV_e_phi_rad"});
  out = DefineOrRedefine(out, "WCSV_e_theta", [](float theta) { return theta * 180.0f / static_cast<float>(M_PI); }, {"WCSV_e_theta_rad"});
  out = out.Filter("WCSV_e_idx >= 0 && WCSV_W >= 0.8 && WCSV_W <= 1.1 && WCSV_e_sector >= 1 && WCSV_e_sector <= 6",
                   "final electron W CSV rows");

  std::ofstream csv(csvPath);
  if (!csv.is_open()) {
    throw std::runtime_error("DVCSAnalysis: cannot open W CSV output: " + csvPath);
  }
  csv << "W,e_p,e_phi,e_theta,e_sector\n";
  csv << std::fixed << std::setprecision(8);

  std::mutex csvMutex;
  size_t rows = 0;
  out.Foreach([&](float w, float p, float phi, float theta, int sector) {
    std::lock_guard<std::mutex> lock(csvMutex);
    csv << w << "," << p << "," << phi << "," << theta << "," << sector << "\n";
    ++rows;
  }, {"WCSV_W", "WCSV_e_p", "WCSV_e_phi", "WCSV_e_theta", "WCSV_e_sector"});

  std::cout << "[SaveOutput] Wrote final corrected electron W CSV: " << csvPath
            << " (rows = " << rows << ")\n";
}

void DVCSAnalysis::SaveOutput() {
  if (fOutputWCSV) {
    if (!dfSelected.has_value()) {
      std::cerr << "DVCSAnalysis::SaveOutput: dfSelected not set!" << std::endl;
      return;
    }
    ROOT::RDF::RNode finalDf = dfSelected_afterFid_afterCorr.has_value()
                                   ? *dfSelected_afterFid_afterCorr
                                   : (dfSelected_afterFid.has_value() ? *dfSelected_afterFid : *dfSelected);
    const std::string csvPath = fOutputDir + "/" + fOutputWCSVName;
    std::cout << "[SaveOutput] Output-W mode ON — writing CSV only, no ROOT snapshots.\n";
    WriteFinalElectronWCSV(finalDf, csvPath);
    if (fIsQADBCut) {
      std::cout << "\n[QADB] total accumulated charge analyzed: " << fQADBCuts->GetAccumulatedCharge() / 1e6 << " mC (Do NOT use this number if you enable MT)\n";
    }
    return;
  }

  if (IsMC) {
    // snapshot of the MC bank for efficiency and other studies
    dforginal->Snapshot(
        "dfSelectedMC", Form("%s/%s", fOutputDir.c_str(), "dfSelectedMC.root"),
        {"MC_Particle_pid", "MC_Particle_px", "MC_Particle_py", "MC_Particle_pz", "MC_Particle_vx", "MC_Particle_vy", "MC_Particle_vz", "MC_Particle_vt", "MC_Event_weight",
         "MC_Event_pbeam",  // include if this exists
         "MC_Event_ptarget", "MC_Event_ebeam"});
  }

  if (!dfSelected.has_value()) {
    std::cerr << "DVCSAnalysis::SaveOutput: dfSelected not set!" << std::endl;
    return;
  }

  if (fOptimizeColumns) {
    std::cout << "[SaveOutput] Column optimisation ON — writing only analysis-used columns.\n";
  } else {
    std::cout << "[SaveOutput] Column optimisation OFF — writing all columns.\n";
  }

  auto resolveColumns = [this](ROOT::RDF::RNode& node) -> std::vector<std::string> {
    if (fOptimizeColumns) {
      return ResolveSnapshotColumns(node, MinimalColumns());
    } else {
      return SafeSnapshotColumns(node, {"EventCutResult"});
    }
  };

  if (!IsReproc) {
    auto cols = resolveColumns(*dfSelected);
    auto cnt = dfSelected->Count();
    dfSelected->Snapshot("dfSelected",
                    Form("%s/%s", fOutputDir.c_str(), "dfSelected.root"), cols);  // triggers loop
    std::cout << "Events selected: " << *cnt << std::endl;
  }
  if (IsReproc) SafeSnapshot(*dfSelected, "dfSelected_reproc", Form("%s/%s", fOutputDir.c_str(), "dfSelected_reproc.root"));
  if (fFiducialCut && dfSelected_afterFid.has_value()) {
    std::cout << "output directory is : " << fOutputDir.c_str() << std::endl;

    if (IsReproc) {SafeSnapshot(*dfSelected_afterFid,"dfSelected_afterFid_reprocessed",
                                Form("%s/%s", fOutputDir.c_str(),"dfSelected_afterFid_reprocessed.root"));
    } else {
      if (!IsMinBooking) {
        const std::string root_afterFid = Form("%s/%s", fOutputDir.c_str(), "dfSelected_afterFid.root");
        auto cnt_afterFid = dfSelected_afterFid->Count();
        dfSelected_afterFid->Snapshot("dfSelected_afterFid", root_afterFid, resolveColumns(*dfSelected_afterFid));
        std::cout << "Events after fiducial selected: " << *cnt_afterFid << std::endl;
      }
    }
  }
  if (fDoMomentumCorrection && dfSelected_afterFid_afterCorr.has_value()) {
    auto cnt_afterFid_afterCorr = dfSelected_afterFid_afterCorr->Count();
    dfSelected_afterFid_afterCorr->Snapshot("dfSelected_afterFid_afterCorr", Form("%s/%s", fOutputDir.c_str(), "dfSelected_afterFid_afterCorr.root"), resolveColumns(*dfSelected_afterFid_afterCorr));
    std::cout << "Events after fiducial and momentum correction selected: " << *cnt_afterFid_afterCorr << std::endl;
    const std::string root_afterFid_afterCorr = Form("%s/%s", fOutputDir.c_str(), "dfSelected_afterFid_afterCorr.root");
    const std::string csvpath = fOutputDir + "/events_per_run_afterFid.csv";
    try {
      ROOT::RDataFrame rdf_afterFid_afterCorr("dfSelected_afterFid_afterCorr", root_afterFid_afterCorr);
      auto items = CountPerRunAndWriteCSV<int>(rdf_afterFid_afterCorr, "RUN_run", csvpath);

      std::cout << "[INFO] Wrote per-run counts to " << csvpath << " (unique runs = " << items.size() << ")\n";
    } catch (const std::exception& e) {
      std::cerr << "[WARN] Failed to write per-run CSV after fiducial snapshot: " << e.what() << std::endl;
    }
  }
  if (fIsQADBCut) {
    std::cout << "\n[QADB] total accumulated charge analyzed: " << fQADBCuts->GetAccumulatedCharge() / 1e6 << " mC (Do NOT use this number if you enable MT)\n";
  }
}

void DVCSAnalysis::SetOutputDir(const std::string& dir) { fOutputDir = dir; }
