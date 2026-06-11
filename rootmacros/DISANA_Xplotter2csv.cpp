#include <THnSparse.h>

#include "../DreamAN/DrawHist/DISANAcomparer.h"
#include "../DreamAN/DrawHist/DrawStyle.h"
#include "../DreamAN/DrawHist/DISANAMath.h"

#include <ROOT/RDataFrame.hxx>

#include <array>
#include <cctype>
#include <cmath>
#include <iostream>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>

bool BranchIsBool(const std::string& file,
                  const std::string& tree,
                  const std::string& branch) {
  TFile f(file.c_str(), "READ");
  auto *t = dynamic_cast<TTree*>(f.Get(tree.c_str()));
  if (!t) throw std::runtime_error("Tree not found: " + tree);
  auto *leaf = t->GetLeaf(branch.c_str());
  if (!leaf) throw std::runtime_error("Leaf not found: " + branch);
  std::string ty = leaf->GetTypeName();
  std::cout << "[BranchIsBool] " << file << "  " << branch
            << " type = " << ty << std::endl;
  return ty == "Bool_t" || ty == "bool" || ty == "vector<bool>";
}

ROOT::RDF::RNode NormalizePassColumns(const std::string& file,
                                      const std::string& tree) {
  ROOT::RDataFrame rdf(tree, file);

  const bool pass_is_bool = BranchIsBool(file, tree, "REC_Particle_pass");
  const bool maxe_is_bool = BranchIsBool(file, tree, "REC_Photon_MaxE");

  ROOT::RDF::RNode df = rdf;

  if (pass_is_bool) {
    df = df.Define("REC_Particle_pass_std",
                   [](const ROOT::VecOps::RVec<bool>& v) {
                     return v; 
                   },
                   {"REC_Particle_pass"});
  } else {
    std::cout << "coverting int pass to bool pass"<<std::endl;
    df = df.Define("REC_Particle_pass_std",
                   [](const ROOT::VecOps::RVec<int>& v) {
                     ROOT::VecOps::RVec<bool> out(v.size());
                     for (size_t i = 0; i < v.size(); ++i)
                       out[i] = (v[i] != 0);
                     return out;
                   },
                   {"REC_Particle_pass"});
  }

  if (maxe_is_bool) {
    df = df.Define("REC_Photon_MaxE_std",
                   [](const ROOT::VecOps::RVec<bool>& v) {
                     return v;
                   },
                   {"REC_Photon_MaxE"});
  } else {
    std::cout << "coverting int maxE to bool pass"<<std::endl;
    df = df.Define("REC_Photon_MaxE_std",
                   [](const ROOT::VecOps::RVec<int>& v) {
                     ROOT::VecOps::RVec<bool> out(v.size());
                     for (size_t i = 0; i < v.size(); ++i)
                       out[i] = (v[i] != 0);
                     return out;
                   },
                   {"REC_Photon_MaxE"});
  }

  return df;
}


ROOT::RDF::RNode RejectPi0TwoPhoton(ROOT::RDF::RNode df_, float beam_energy);
ROOT::RDF::RNode SelectPi0Event(ROOT::RDF::RNode df);

ROOT::RDF::RNode ApplyFinalDVCSSelections(
    ROOT::RDF::RNode df,
    const std::string& rec_csv = "dvcs_cuts_rec_corr.csv",
    double windowScaleFactor = 1.0);
ROOT::RDF::RNode ApplyFinalDVCSRadSelections(ROOT::RDF::RNode df);
ROOT::RDF::RNode ApplyFinalGenDVCSSelections(
    ROOT::RDF::RNode df,
    const std::string& rec_csv = "dvcs_cuts_gen_corr.csv",
    double windowScaleFactor = 1.0);

ROOT::RDF::RNode DefineDVPi0Pass(ROOT::RDF::RNode df);
ROOT::RDF::RNode DefineGenDVPi0Pass(ROOT::RDF::RNode df);
ROOT::RDF::RNode ApplyFinalDVPi0Selections(ROOT::RDF::RNode df);
ROOT::RDF::RNode ApplyFinalGenDVPi0Selections(ROOT::RDF::RNode df);

ROOT::RDF::RNode InitKinematics(const std::string& filename_ = "", const std::string& treename_ = "", float beam_energy = 0);
ROOT::RDF::RNode momentumcorr(ROOT::RDF::RNode df_, bool applyCorrection = true,
                              double correctionScaleFactor = 1.0);
ROOT::RDF::RNode Init2PhotonKinematics(ROOT::RDF::RNode df_, float beam_energy = 0);
ROOT::RDF::RNode InitGenKinematics(const std::string& filename_ = "", const std::string& treename_ = "", float beam_energy = 0);

ROOT::RDF::RNode GetSlim_exclusive(ROOT::RDF::RNode df_, const std::string& filename_slim, const std::string& treename_slim, bool isGen = false);
ROOT::RDF::RNode WriteSlimAndReload_exclusive(ROOT::RDF::RNode df_, const std::string& filename_slim, const std::string& treename_slim, bool isGen = false);

void PlotAllRecoDistributions(ROOT::RDF::RNode df, const std::string& out = "reco_kinematics_grid.png",
                              int bins_p = 120, int bins_theta = 120, int bins_phi = 120);
void WriteRunEventCounts(ROOT::RDF::RNode df,
                         const std::string& out_csv,
                         const std::string& label = "");

static double MomentumFunc(float px, float py, float pz) { return std::sqrt(px * px + py * py + pz * pz); }
static double ThetaFunc(float px, float py, float pz) { return std::acos(pz / std::sqrt(px * px + py * py + pz * pz)); }
static double PhiFunc(float px, float py) {
  double phi = std::atan2(py, px);
  return phi < 0 ? phi + 2 * M_PI : phi;
}

static double gDVPi0WindowScaleFactor = 1.0;

static inline void SetDVPi0WindowScaleFactor(double factor) {
  if (!std::isfinite(factor) || factor < 0.0) {
    throw std::invalid_argument(
        "DVpi0 window scale factor must be finite and non-negative");
  }
  gDVPi0WindowScaleFactor = factor;
}

bool Inrange(double var, double min, double max) {
  return (var >= min && var < max);
}

bool InrangeVaried(double var, double min, double max) {
  const double center = 0.5 * (min + max);
  const double halfWidth =
      0.5 * gDVPi0WindowScaleFactor * (max - min);
  return (var >= center - halfWidth && var < center + halfWidth);
}

static double MomentumCorrectionScale7p5(bool isPhoton, int det, double p, double theta) {
  if (!std::isfinite(p) || !std::isfinite(theta) || p <= 0.0) return 1.0;

  double a0_p0 = 1.0;
  double a0_p1 = 0.0;
  double a1_p0 = 0.0;
  double a1_p1 = 0.0;
  bool hasCorrection = false;

  if (isPhoton && det == 0) {  // FT photon energy
    a0_p0 = 1.01384962371;
    a0_p1 = -0.00345633005311;
    a1_p0 = 0.00392055912406;
    a1_p1 = 0.000101895000551;
    hasCorrection = true;
  } else if (isPhoton && det == 1) {  // FD photon energy
    a0_p0 = 1.07754923622;
    a0_p1 = -0.00196102164162;
    a1_p0 = 0.00448549463246;
    a1_p1 = -0.00179558909233;
    hasCorrection = true;
  } else if (!isPhoton && det == 1) {  // FD proton momentum
    a0_p0 = 1.03324429434;
    a0_p1 = -0.0627598761543;
    a1_p0 = 0.0019642331647;
    a1_p1 = -0.000215812664027;
    hasCorrection = true;
  } else if (!isPhoton && det == 2) {  // CD proton momentum
    a0_p0 = 0.966806323053;
    a0_p1 = 0.0711872437488;
    a1_p0 = 0.000536019455192;
    a1_p1 = -0.00140324153035;
    hasCorrection = true;
  }

  if (!hasCorrection) return 1.0;

  const double thetaDeg = theta * 180.0 / M_PI;
  const double a0 = a0_p0 + a0_p1 * p;
  const double a1 = a1_p0 + a1_p1 * p;
  const double scale = a0 + a1 * thetaDeg;
  return (std::isfinite(scale) && scale > 0.0) ? scale : 1.0;
}

template <typename Method>
ROOT::RDF::RNode redefine_DISCAT(ROOT::RDF::RNode node, const std::string& name, const Method method, float beam_energy) {
  return node.Redefine(name,
                       [method, beam_energy](double recel_p, double recel_theta, double recel_phi, double recpro_p, double recpro_theta, double recpro_phi,
                                             double recpho_p, double recpho_theta, double recpho_phi) {
                         return (DISANAMath(beam_energy, recel_p, recel_theta, recel_phi, recpro_p, recpro_theta, recpro_phi, recpho_p, recpho_theta,
                                            recpho_phi).*method)();
                       },
                       {"recel_p", "recel_theta", "recel_phi", "recpro_p", "recpro_theta", "recpro_phi", "recpho_p", "recpho_theta", "recpho_phi"});
}

ROOT::RDF::RNode momentumcorr(ROOT::RDF::RNode df_, bool applyCorrection,
                              double correctionScaleFactor) {
  if (!applyCorrection) return df_;
  if (!std::isfinite(correctionScaleFactor) || correctionScaleFactor < 0.0) {
    throw std::invalid_argument(
        "momentumcorr correctionScaleFactor must be finite and non-negative");
  }

  // Scale the correction relative to no correction:
  // factor=0 -> no correction, factor=1 -> nominal correction.
  const auto varyCorrection = [correctionScaleFactor](double nominalScale) {
    const double variedScale =
        1.0 + correctionScaleFactor * (nominalScale - 1.0);
    return (std::isfinite(variedScale) && variedScale > 0.0)
               ? variedScale
               : 1.0;
  };

  df_ = df_
      .Define("momcorr_pho_scale",
              [varyCorrection](double p, double theta, int det) {
                return varyCorrection(
                    MomentumCorrectionScale7p5(true, det, p, theta));
              },
              {"recpho_p", "recpho_theta", "pho_det_region"})
      .Define("momcorr_pro_scale",
              [varyCorrection](double p, double theta, int det) {
                return varyCorrection(
                    MomentumCorrectionScale7p5(false, det, p, theta));
              },
              {"recpro_p", "recpro_theta", "pro_det_region"})
      .Redefine("pho_px",
                [](float px, double scale) -> float { return static_cast<float>(px * scale); },
                {"pho_px", "momcorr_pho_scale"})
      .Redefine("pho_py",
                [](float py, double scale) -> float { return static_cast<float>(py * scale); },
                {"pho_py", "momcorr_pho_scale"})
      .Redefine("pho_pz",
                [](float pz, double scale) -> float { return static_cast<float>(pz * scale); },
                {"pho_pz", "momcorr_pho_scale"})
      .Redefine("recpho_p",
                [](double p, double scale) { return p * scale; },
                {"recpho_p", "momcorr_pho_scale"})
      .Redefine("pro_px",
                [](float px, double scale) -> float { return static_cast<float>(px * scale); },
                {"pro_px", "momcorr_pro_scale"})
      .Redefine("pro_py",
                [](float py, double scale) -> float { return static_cast<float>(py * scale); },
                {"pro_py", "momcorr_pro_scale"})
      .Redefine("pro_pz",
                [](float pz, double scale) -> float { return static_cast<float>(pz * scale); },
                {"pro_pz", "momcorr_pro_scale"})
      .Redefine("recpro_p",
                [](double p, double scale) { return p * scale; },
                {"recpro_p", "momcorr_pro_scale"});

  const float beam_energy = 7.546;
  df_ = redefine_DISCAT(df_, "Q2", &DISANAMath::GetQ2, beam_energy);
  df_ = redefine_DISCAT(df_, "xB", &DISANAMath::GetxB, beam_energy);
  df_ = redefine_DISCAT(df_, "t", &DISANAMath::GetT, beam_energy);
  df_ = redefine_DISCAT(df_, "phi", &DISANAMath::GetPhi, beam_energy);
  df_ = redefine_DISCAT(df_, "W", &DISANAMath::GetW, beam_energy);
  df_ = redefine_DISCAT(df_, "nu", &DISANAMath::GetNu, beam_energy);
  df_ = redefine_DISCAT(df_, "y", &DISANAMath::Gety, beam_energy);
  df_ = redefine_DISCAT(df_, "Mx2_ep", &DISANAMath::GetMx2_ep, beam_energy);
  df_ = redefine_DISCAT(df_, "Emiss", &DISANAMath::GetEmiss, beam_energy);
  df_ = redefine_DISCAT(df_, "PTmiss", &DISANAMath::GetPTmiss, beam_energy);
  df_ = redefine_DISCAT(df_, "Mx2_epg", &DISANAMath::GetMx2_epg, beam_energy);
  df_ = redefine_DISCAT(df_, "Mx2_eg", &DISANAMath::GetMx2_egamma, beam_energy);
  df_ = redefine_DISCAT(df_, "Theta_e_gamma", &DISANAMath::GetTheta_e_gamma, beam_energy);
  df_ = redefine_DISCAT(df_, "DeltaE", &DISANAMath::GetDeltaE, beam_energy);
  df_ = redefine_DISCAT(df_, "DeltaPhi", &DISANAMath::GetDeltaPhi, beam_energy);
  df_ = redefine_DISCAT(df_, "Theta_gamma_gamma", &DISANAMath::GetTheta_gamma_gamma, beam_energy);

  return df_;
}

/// styling plots
// double double titleSize = 0.05, double labelSize = 0.04,double xTitleOffset = 1.1, double yTitleOffset = 1.6, int font = 42, int maxDigits = 5, int nDivisions = 510, double
// leftMargin = 0.16, double rightMargin = 0.07, double bottomMargin = 0.13, double topMargin = 0.06
DrawStyle KinStyle(0.07, 0.06, 0.9, 1.2);                                       // For Kin plots
DrawStyle dvcsStyle(0.06, 0.06, 1.2, 1.4, 42, 5, 510, 0.14, 0.07, 0.13, 0.06);  // For DVCS plots
DrawStyle csStyle(0.05, 0.04, 1.0, 1.3);                                        // For Cross-Sections
DrawStyle bsaStyle(0.06, 0.045, 1.0, 1.2);                                      // For BSA

// for exclusivity plots
std::vector<std::pair<std::string, std::string>> detCuts = {
    {"pho_det_region == 0 && pro_det_region == 2", "FT-CD"},
    {"pho_det_region == 1 && pro_det_region == 2", "FD-CD"},
    {"pho_det_region == 1 && pro_det_region == 1", "FD-FD"},
};

std::vector<std::pair<std::string, std::string>> detCutsPi0 = {
    {"pho_det_region == 0 && pho2_det_region == 0 && pro_det_region == 2", "FT-FT-CD"},
    {"pho_det_region == 1 && pho2_det_region == 1 && pro_det_region == 2", "FD-FD-CD"},
    {"pho_det_region == 1 && pho2_det_region == 1 && pro_det_region == 1", "FD-FD-FD"},
};

template <typename Method>
ROOT::RDF::RNode define_DISCAT(ROOT::RDF::RNode node, const std::string& name, const Method method, float beam_energy) {
  return node.Define(name,
                     [method, beam_energy](double recel_p, double recel_theta, double recel_phi, double recpro_p, double recpro_theta, double recpro_phi, double recpho_p,
                                           double recpho_theta, double recpho_phi) {
                       return (DISANAMath(beam_energy, recel_p, recel_theta, recel_phi, recpro_p, recpro_theta, recpro_phi, recpho_p, recpho_theta, recpho_phi).*method)();
                     },
                     {"recel_p", "recel_theta", "recel_phi", "recpro_p", "recpro_theta", "recpro_phi", "recpho_p", "recpho_theta", "recpho_phi"});
}

template <typename Method>
ROOT::RDF::RNode define_DISCAT_pi0(ROOT::RDF::RNode node, const std::string& name, const Method method, float beam_energy) {
  return node.Define(name,
                     [method, beam_energy](double recel_p, double recel_theta, double recel_phi, double recpro_p, double recpro_theta, double recpro_phi,
                                           double recpho_p, double recpho_theta, double recpho_phi, double recpho2_p, double recpho2_theta,
                                           double recpho2_phi) {
                       return (DISANAMath(Pi0Tag{}, beam_energy, recel_p, recel_theta, recel_phi, recpro_p, recpro_theta, recpro_phi, recpho_p, recpho_theta,
                                          recpho_phi, recpho2_p, recpho2_theta, recpho2_phi).*method)();
                     },
                     {"recel_p", "recel_theta", "recel_phi", "recpro_p", "recpro_theta", "recpro_phi", "recpho_p", "recpho_theta", "recpho_phi",
                      "recpho2_p", "recpho2_theta", "recpho2_phi"});
}

void WriteRunEventCounts(ROOT::RDF::RNode df,
                         const std::string& out_csv,
                         const std::string& label) {
  auto runs = df.Take<int>("RUN_config_run");

  std::map<int, unsigned long long> counts_by_run;
  unsigned long long total = 0;
  for (const auto run : *runs) {
    ++counts_by_run[run];
    ++total;
  }

  std::ofstream fout(out_csv);
  if (!fout) {
    std::cerr << "[WriteRunEventCounts] ERROR: cannot open output file: "
              << out_csv << std::endl;
    return;
  }

  fout << "RUN_config_run,event_count\n";
  for (const auto& item : counts_by_run) {
    fout << item.first << "," << item.second << "\n";
  }

  std::cout << "[WriteRunEventCounts] ";
  if (!label.empty()) std::cout << label << ": ";
  std::cout << "saved " << counts_by_run.size()
            << " runs, total events = " << total
            << " to " << out_csv << std::endl;
}

void DISANA_Xplotter2csv() {
  bool ComputeBgk_core = false;  // Set to true if you want to compute background
  bool DoBkgCorr = true;       // Set to true if you want to apply background correction

  ROOT::EnableImplicitMT(40);

  std::string input_path_from_analysisRun_7546_data = "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File7546BSA/";

  std::string input_path_from_analysisRun_7546_pi0MC = "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File7546BSA/Pi0MC/";

  std::string input_path_from_analysisRun_7546_dvcsmc_gen = "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File7546/nobkgall";
  std::string input_path_from_analysisRun_7546_dvcsmc_rec = "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File7546/nobkg"; //(no bkg merged)

  std::string input_path_from_analysisRun_7546_dvcsmc_bkg = "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File7546/bkg";
  std::string input_path_from_analysisRun_7546_dvcsmc_nobkg = "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File7546/nobkg";

  std::string filename_afterFid_7546_data = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_7546_data.c_str());

  std::string filename_afterFid_7546_pi0MC = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_7546_pi0MC.c_str());

  std::string filename_afterFid_7546_dvcsmc_gen = Form("%s/dfSelected.root", input_path_from_analysisRun_7546_dvcsmc_gen.c_str());
  std::string filename_afterFid_7546_dvcsmc_rec = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_7546_dvcsmc_rec.c_str());

  std::string filename_afterFid_7546_dvcsmc_bkg = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_7546_dvcsmc_bkg.c_str());
  std::string filename_afterFid_7546_dvcsmc_nobkg = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_7546_dvcsmc_nobkg.c_str());

  std::string filename_afterFid_7546_dvcsmc_rad = "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File7546/radandP1/raddelta0p1v0p6Max.root";
  std::string filename_afterFid_7546_dvcsmc_norad = "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File7546/radandP1/nor.root";
  std::string filename_afterFid_7546_dvcsmc_p1cut = "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File7546/radandP1/norP1_2.root";

  std::string filename_afterFid_7546_dvcsmc_km15_rad = "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File7546/radandP1/radkm15_7546.root";
  std::string filename_afterFid_7546_dvcsmc_km15_norad = "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File7546/radandP1/noradkm15_7546.root";
  

  float beam_energy = 7.546;
  bool momentumcorr_sys = false;  // Set to true to apply momentum correction systematic variation
  float momentumcorr_scale_factor = 0.9;  // Set to 1.0 for nominal, >1.0 for up variation, <1.0 for down variation
  bool dvcs_selection_sys = false;  // Vary DVCS CSV exclusivity windows
  float dvcs_selection_scale_factor = 0.9;
  bool dvpi0_selection_sys = false;  // Vary only DVpi0 InrangeVaried windows
  float dvpi0_selection_scale_factor = 0.9;

  const double momentumcorr_factor =
      momentumcorr_sys ? momentumcorr_scale_factor : 1.0;
  std::cout << "[momentumcorr] mode="
            << (momentumcorr_sys ? "systematic" : "nominal")
            << ", correction scale factor=" << momentumcorr_factor
            << std::endl;
  const double dvcs_selection_factor =
      dvcs_selection_sys ? dvcs_selection_scale_factor : 1.0;
  const double dvpi0_selection_factor =
      dvpi0_selection_sys ? dvpi0_selection_scale_factor : 1.0;
  std::cout << "[DVCS selection] mode="
            << (dvcs_selection_sys ? "systematic" : "nominal")
            << ", window half-width scale factor=" << dvcs_selection_factor
            << std::endl;
  std::cout << "[DVpi0 selection] mode="
            << (dvpi0_selection_sys ? "systematic" : "nominal")
            << ", window half-width scale factor=" << dvpi0_selection_factor
            << std::endl;
  SetDVPi0WindowScaleFactor(dvpi0_selection_factor);

  ROOT::RDF::RNode df_afterFid_7546_data_init = momentumcorr(InitKinematics(filename_afterFid_7546_data,"dfSelected_afterFid_afterCorr", beam_energy),true, momentumcorr_factor);

  ROOT::RDF::RNode df_afterFid_7546_pi0MC_init = InitKinematics(filename_afterFid_7546_pi0MC, "dfSelected_afterFid_afterCorr", beam_energy);

  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_gen_init = InitGenKinematics(filename_afterFid_7546_dvcsmc_gen, "dfSelected", beam_energy);
  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_rec_init = InitKinematics(filename_afterFid_7546_dvcsmc_rec, "dfSelected_afterFid_afterCorr", beam_energy);

  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_bkg_init = InitKinematics(filename_afterFid_7546_dvcsmc_bkg, "dfSelected_afterFid_afterCorr", beam_energy);
  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_nobkg_init = InitKinematics(filename_afterFid_7546_dvcsmc_nobkg, "dfSelected_afterFid_afterCorr", beam_energy);

  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_rad_init = InitGenKinematics(filename_afterFid_7546_dvcsmc_rad, "MC", beam_energy);
  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_norad_init = InitGenKinematics(filename_afterFid_7546_dvcsmc_norad, "MC", beam_energy);
  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_p1cut_init = InitGenKinematics(filename_afterFid_7546_dvcsmc_p1cut, "MC", beam_energy);

  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_km15_rad_init = InitGenKinematics(filename_afterFid_7546_dvcsmc_km15_rad, "MC", beam_energy);
  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_km15_norad_init = InitGenKinematics(filename_afterFid_7546_dvcsmc_km15_norad, "MC", beam_energy);

  ROOT::RDF::RNode df_afterFid_7546_data = GetSlim_exclusive(df_afterFid_7546_data_init, "dfSlim_7546_data.root", "dfSlim_7546_data", false);
  ROOT::RDF::RNode df_afterFid_7546_pi0MC = GetSlim_exclusive(df_afterFid_7546_pi0MC_init, "dfSlim_7546_pi0MC.root", "dfSlim_7546_pi0MC", false);
  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_gen = GetSlim_exclusive(df_afterFid_7546_dvcsmc_gen_init, "dfSlim_7546_dvcsmc_gen.root", "dfSlim_7546_dvcsmc_gen", true);
  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_rec = GetSlim_exclusive(df_afterFid_7546_dvcsmc_rec_init, "dfSlim_7546_dvcsmc_rec.root", "dfSlim_7546_dvcsmc_rec", false);
  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_bkg = GetSlim_exclusive(df_afterFid_7546_dvcsmc_bkg_init, "dfSlim_7546_dvcsmc_bkg.root", "dfSlim_7546_dvcsmc_bkg", false);
  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_nobkg = GetSlim_exclusive(df_afterFid_7546_dvcsmc_nobkg_init, "dfSlim_7546_dvcsmc_nobkg.root", "dfSlim_7546_dvcsmc_nobkg", false);
  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_rad_temp = GetSlim_exclusive(df_afterFid_7546_dvcsmc_rad_init, "dfSlim_7546_dvcsmc_rad.root", "dfSlim_7546_dvcsmc_rad", true);
  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_norad = GetSlim_exclusive(df_afterFid_7546_dvcsmc_norad_init, "dfSlim_7546_dvcsmc_norad.root", "dfSlim_7546_dvcsmc_norad", true);
  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_p1cut = GetSlim_exclusive(df_afterFid_7546_dvcsmc_p1cut_init, "dfSlim_7546_dvcsmc_p1cut.root", "dfSlim_7546_dvcsmc_p1cut", true);

  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_km15_rad = GetSlim_exclusive(df_afterFid_7546_dvcsmc_km15_rad_init, "dfSlim_7546_dvcsmc_km15_rad.root", "dfSlim_7546_dvcsmc_km15_rad", true);
  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_km15_norad = GetSlim_exclusive(df_afterFid_7546_dvcsmc_km15_norad_init, "dfSlim_7546_dvcsmc_km15_norad.root", "dfSlim_7546_dvcsmc_km15_norad", true);

  DrawStyle fitStyle(0.06, 0.05, 1.0, 1.3);  // You can tweak this

  // Apply final DVCS cuts
  auto df_final_dvcs_7546_data = ApplyFinalDVCSSelections(df_afterFid_7546_data, "dvcs_cuts_rec_corr.csv",dvcs_selection_factor);
  auto df_final_dvcsPi_rejected_7546_data = RejectPi0TwoPhoton(df_final_dvcs_7546_data, beam_energy);

  auto df_final_dvcs_7546_pi0MC = ApplyFinalGenDVCSSelections(df_afterFid_7546_pi0MC, "dvcs_cuts_gen_corr.csv",dvcs_selection_factor);
  auto df_final_dvcsPi_rejected_7546_pi0MC = RejectPi0TwoPhoton(df_final_dvcs_7546_pi0MC, beam_energy);

  auto df_final_OnlPi0_7546_data = ApplyFinalDVPi0Selections(Init2PhotonKinematics(SelectPi0Event(df_afterFid_7546_data), beam_energy));
  auto df_final_OnlPi0_7546_pi0MC = ApplyFinalGenDVPi0Selections(Init2PhotonKinematics(SelectPi0Event(df_afterFid_7546_pi0MC), beam_energy));

  auto df_final_dvcs_7546_dvcsmc_rec = ApplyFinalGenDVCSSelections(df_afterFid_7546_dvcsmc_rec, "dvcs_cuts_gen_corr.csv",dvcs_selection_factor);
  auto df_final_dvcsPi_rejected_7546_dvcsmc_rec = RejectPi0TwoPhoton(df_final_dvcs_7546_dvcsmc_rec, beam_energy);

  auto df_final_dvcs_7546_dvcsmc_bkg = ApplyFinalGenDVCSSelections(df_afterFid_7546_dvcsmc_bkg, "dvcs_cuts_gen_corr.csv",dvcs_selection_factor);
  auto df_final_dvcsPi_rejected_7546_dvcsmc_bkg = RejectPi0TwoPhoton(df_final_dvcs_7546_dvcsmc_bkg, beam_energy);

  auto df_final_dvcs_7546_dvcsmc_nobkg = ApplyFinalGenDVCSSelections(df_afterFid_7546_dvcsmc_nobkg, "dvcs_cuts_gen_corr.csv",dvcs_selection_factor);
  auto df_final_dvcsPi_rejected_7546_dvcsmc_nobkg = RejectPi0TwoPhoton(df_final_dvcs_7546_dvcsmc_nobkg, beam_energy);

  auto df_afterFid_7546_dvcsmc_rad = df_afterFid_7546_dvcsmc_rad_temp;

  DISANAcomparer comparer;
  comparer.SetOutputDir("./");
  comparer.SetKinStyle(KinStyle);
  comparer.SetDVCSStyle(dvcsStyle);
  comparer.SetCrossSectionStyle(csStyle);
  comparer.SetBSAStyle(bsaStyle);

  comparer.SetDVCSWeightFunction([](const DVCSWeightInput& v) -> double {
    // Configure DVCS event weights here. Available variables:
    // v.pho_det_region, v.pro_det_region,
    // v.recel_p, v.recel_theta, v.recel_phi,
    // v.recpho_p, v.recpho_theta, v.recpho_phi,
    // v.recpro_p, v.recpro_theta, v.recpro_phi.
    //if (v.pho_det_region == 1) return 1.0; // FD, e.g. 0.7828
    //if (v.pho_det_region == 0) return 1.0; // FT, e.g. 1.2079
    if (v.pro_det_region == 2 ){
      double pro_theta = v.recpro_theta*180.0/M_PI;
      double weight2  = 174.874782895;
      weight2 = weight2 - 13.4093933552 * pro_theta;
      weight2 = weight2 + 0.386737126334 * pro_theta * pro_theta;
      weight2 = weight2 - 0.00494497943272 * pro_theta * pro_theta * pro_theta;
      weight2 = weight2 + 0.0000236221306856 * pro_theta * pro_theta * pro_theta * pro_theta;
      //return 1/(weight2);
      return 1;
    }
    if (v.pro_det_region == 1 ){
      double pro_theta = v.recpro_theta*180.0/M_PI;
      double weight1  =   7.62332107254;
      weight1 = weight1 - 0.804406387543 * pro_theta;
      weight1 = weight1 + 0.0327055624737 * pro_theta * pro_theta;
      weight1 = weight1 - 0.000521842763762 * pro_theta * pro_theta * pro_theta;
      weight1 = weight1 + 0.00000241692488631 * pro_theta * pro_theta * pro_theta * pro_theta;
      //return 1/(weight1);
      return 1;
    }
    return 1.0;
  });

  comparer.PlotIndividual(false);
  /// bins for cross-section plots
  BinManager xBins;

  //xBins.SetQ2Bins({1.0, 1.5, 2.0});
  //xBins.SetTBins({0.1, 0.2, 0.6});
  //xBins.SetXBBins({0.15, 0.2, 0.3});

  //xBins.SetQ2Bins({1.0, 1.25, 1.5, 1.7, 2.0});
  //xBins.SetQ2Bins({1.0, 1.25});
  //xBins.SetTBins({0.11, 0.15, 0.25, 0.40, 0.60, 0.80, 1.00});
  //xBins.SetTBins({0.40, 0.60});
  //xBins.SetXBBins({0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3});
  //xBins.SetXBBins({0.15, 0.175});

  //*******xBins.SetQ2Bins({1.00, 1.25, 1.50, 1.75, 2.00, 2.40, 2.90, 3.50});
  //*******xBins.SetTBins({0.13, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0});
  //*******xBins.SetXBBins({0.125, 0.150, 0.180, 0.210, 0.240, 0.285, 0.350, 0.430});

  //xBins.SetQ2Bins({1.0, 1.25, 1.5, 2.0});
  //xBins.SetTBins({0.2, 0.4, 0.6, 1.0});
  //xBins.SetXBBins({0.125, 0.15, 0.175, 0.2, 0.25, 0.3});

  xBins.SetQ2Bins({1.00, 1.20, 1.50, 2.00, 2.80, 6.00});
  xBins.SetTBins({0.13, 0.23, 0.43, 1.0, 2.0});
  // One xB edge vector for each Q2 interval above. For example, with
  // Q2 edges {1.0, 1.5, 2.0}, pass {{...for 1.0-1.5...}, {...for 1.5-2.0...}}.
  xBins.SetXBBinsByQ2({
      {0.090, 0.160, 0.285}, // for 1.00-1.25
      {0.120, 0.200, 0.330}, // for 1.25-1.50
      {0.140, 0.250, 0.400}, // for 1.50-2.00
      {0.190, 0.480}, // for 2.00-2.80
      {0.270, 0.630} // for 2.80-6.00
  });

  comparer.SetXBinsRanges(xBins);

  /*df_final_dvcsPi_rejected_7546_data.Count();
  df_final_OnlPi0_7546_data.Count();
  df_final_dvcsPi_rejected_7546_pi0MC.Count();
  df_final_OnlPi0_7546_pi0MC.Count();
  df_afterFid_7546_dvcsmc_gen.Count();
  df_final_dvcsPi_rejected_7546_dvcsmc_rec.Count();
  df_final_dvcsPi_rejected_7546_dvcsmc_bkg.Count();
  df_final_dvcsPi_rejected_7546_dvcsmc_nobkg.Count();*/

  double charge=8.2054; // mC (5681-5757, 5757-5870(trigger prescale 2), 5694, 5699, 5702, 5707, 5725, 5733, 5758, 5759, 5767, 5850, 5855 removed)
  std::cout<<"Total effective charge (mC): "<<charge<<std::endl;
  double luminosity = (charge)*1.33*pow(10,6);  // Set your desired luminosity here nb^-1
  std::cout<<"Luminosity (nb^-1): "<<luminosity<<std::endl;
  double polarisation = 0.8617;  // Set your desired polarisation here

  //df_final_dvcsPi_rejected_7546_dvcsmc_rec = df_final_dvcsPi_rejected_7546_dvcsmc_rec.Filter("t < 0.3 && t >0.2", "Cut: t > 0.5 GeV^2");

  /*comparer.AddModelwithPi0Corr(df_final_dvcsPi_rejected_7546_dvcsmc_rec,
                              //df_afterFid_7546_dvcsmc_gen,
                              df_final_OnlPi0_7546_pi0MC,
                              df_final_dvcsPi_rejected_7546_pi0MC,
                              df_final_OnlPi0_7546_pi0MC,
                              df_afterFid_7546_dvcsmc_gen,
                              df_final_dvcsPi_rejected_7546_dvcsmc_rec,
                              df_final_dvcsPi_rejected_7546_dvcsmc_bkg,
                              df_final_dvcsPi_rejected_7546_dvcsmc_nobkg,
                              df_afterFid_7546_dvcsmc_rad,
                              df_afterFid_7546_dvcsmc_norad,
                              df_afterFid_7546_dvcsmc_p1cut,
                              "RGK 7.5GeV mc", beam_energy, true, true, true, true, true);*/

  //df_final_dvcsPi_rejected_7546_data = df_final_dvcsPi_rejected_7546_data.Filter("t < 1.0 && t > 0.8", "Cut: t > 0.5 GeV^2");
  //df_final_OnlPi0_7546_data = df_final_OnlPi0_7546_data.Filter("t < 1.0 && t > 0.8", "Cut: t > 0.5 GeV^2");
  df_final_dvcsPi_rejected_7546_data = df_final_dvcsPi_rejected_7546_data.Filter(
                                        "RUN_config_run!=5694 && "
                                        "RUN_config_run!=5699 && "
                                        "RUN_config_run!=5702 && "
                                        "RUN_config_run!=5707 && "
                                        "RUN_config_run!=5725 && "
                                        "RUN_config_run!=5733 && "
                                        "RUN_config_run!=5758 && "
                                        "RUN_config_run!=5759 && "
                                        "RUN_config_run!=5767 && "
                                        "RUN_config_run!=5850 && "
                                        "RUN_config_run!=5855",
                                        "Cut: runs with bad beam conditions");
  df_final_OnlPi0_7546_data = df_final_OnlPi0_7546_data.Filter(
                                        "RUN_config_run!=5694 && "
                                        "RUN_config_run!=5699 && "
                                        "RUN_config_run!=5702 && "
                                        "RUN_config_run!=5707 && "
                                        "RUN_config_run!=5725 && "
                                        "RUN_config_run!=5733 && "
                                        "RUN_config_run!=5758 && "
                                        "RUN_config_run!=5759 && "
                                        "RUN_config_run!=5767 && "
                                        "RUN_config_run!=5850 && "
                                        "RUN_config_run!=5855",
                                        "Cut: runs with bad beam conditions");         

  // Uncomment or copy one of these lines to count events by RUN_config_run
  // for any specific dataframe you want to inspect.
  
  WriteRunEventCounts(df_final_dvcsPi_rejected_7546_data,
                       "run7p5_counts_dvcs_data.csv",
                       "DVCS data after bad-run removal");
  WriteRunEventCounts(df_final_dvcsPi_rejected_7546_dvcsmc_bkg,
                       "run7p5_counts_dvcs_bkg.csv",
                        "DVCS background MC");
  WriteRunEventCounts(df_final_dvcsPi_rejected_7546_dvcsmc_nobkg,
                       "run7p5_counts_dvcs_nobkg.csv",
                        "DVCS no-background MC");
  
  /*comparer.AddModelwithPi0Corr(//df_final_dvcsPi_rejected_7546_data,
                              df_final_dvcsPi_rejected_7546_dvcsmc_rec,
                              df_final_OnlPi0_7546_data,
                              df_final_dvcsPi_rejected_7546_pi0MC,
                              df_final_OnlPi0_7546_pi0MC,
                              df_afterFid_7546_dvcsmc_gen,
                              df_final_dvcsPi_rejected_7546_dvcsmc_rec,
                              df_final_dvcsPi_rejected_7546_dvcsmc_bkg,
                              df_final_dvcsPi_rejected_7546_dvcsmc_nobkg,
                              df_afterFid_7546_dvcsmc_km15_rad,
                              df_afterFid_7546_dvcsmc_km15_norad,
                              df_afterFid_7546_dvcsmc_km15_norad,
                              "RGK 7.5 GeV MC", beam_energy, true, true, true, true, true, luminosity, 39.32, 45, 0.9837);*/

  comparer.AddModelwithPi0Corr(df_final_dvcsPi_rejected_7546_data,
                              //df_afterFid_7546_dvcsmc_gen,
                              df_final_OnlPi0_7546_data,
                              df_final_dvcsPi_rejected_7546_pi0MC,
                              df_final_OnlPi0_7546_pi0MC,
                              df_afterFid_7546_dvcsmc_gen,
                              df_final_dvcsPi_rejected_7546_dvcsmc_rec,
                              df_final_dvcsPi_rejected_7546_dvcsmc_bkg,
                              df_final_dvcsPi_rejected_7546_dvcsmc_nobkg,
                              df_afterFid_7546_dvcsmc_rad,
                              df_afterFid_7546_dvcsmc_norad,
                              df_afterFid_7546_dvcsmc_p1cut,
                              "RGK 7.5GeV", beam_energy, true, true, true, true, true, luminosity, 39.32, 45, 0.9837/*0.9069/0.9647*/);

  //comparer.PlotKinematicComparison();
  //comparer.PlotPi0KinematicComparison();
  //comparer.PlotxBQ2tBin();
  //comparer.PlotxBQ2tBinPi0();
  //comparer.PlotxBQ2tBinMC();
  //comparer.PlotxBQ2tBinPi0MC();
  //comparer.PlotDVCSKinematicsComparison();
  //comparer.PlotDVPi0KinematicsComparison();
  //comparer.PlotDIS_BSA_Cross_Section_AndCorr_Comparison(polarisation, true, true, true, true, true, true, true, true);   
  //comparer.PlotDISCrossSectionComparison(luminosity);  // argument is Luminosity, polarisation
  //comparer.PlotDIS_BSA_Comparison(luminosity, polarisation);         // argument is Luminosity
  //comparer.PlotDIS_Pi0CorrComparison();
  //comparer.PlotMomentumCorrection();
  //comparer.PlotExclusivityComparisonByDetectorCases(detCuts);
  comparer.PlotExclusivityComparisonByDetectorCaseswithPi0(detCuts);
  bool draw_dvpi0_mc = true;
  bool output_wide_mpi0 = false;
  //comparer.PlotPi0ExclusivityComparisonByDetectorCases(detCutsPi0, draw_dvpi0_mc, output_wide_mpi0);
  gApplication->Terminate(0);
}

ROOT::RDF::RNode InitKinematics(const std::string& filename_, const std::string& treename_, float beam_energy) {
  //ROOT::RDataFrame rdf(treename_, filename_);
  //auto df_ = std::make_unique<ROOT::RDF::RNode>(rdf);
  auto df_ = NormalizePassColumns(filename_, treename_);
  df_ = df_.Define("ele_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass_std"})
             .Define("ele_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass_std"})
             .Define("ele_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass_std"})
             .Define("recel_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& vz, const ROOT::VecOps::RVec<bool>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && pass[i]) return vz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass_std"})
             .Define("pho_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& trackpass, const ROOT::VecOps::RVec<bool>& maxEpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && trackpass[i] && maxEpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass_std", "REC_Photon_MaxE_std"})
             .Define("pho_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<bool>& trackpass, const ROOT::VecOps::RVec<bool>& maxEpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && trackpass[i] && maxEpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass_std", "REC_Photon_MaxE_std"})
             .Define("pho_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<bool>& trackpass, const ROOT::VecOps::RVec<bool>& maxEpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && trackpass[i] && maxEpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass_std", "REC_Photon_MaxE_std"})
             .Define("recpho_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& vz, const ROOT::VecOps::RVec<bool>& pass, const ROOT::VecOps::RVec<bool>& maxEpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && pass[i] && maxEpass[i]) return vz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass_std", "REC_Photon_MaxE_std"})
             .Define("recpho_beta",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& beta, const ROOT::VecOps::RVec<bool>& trackpass, const ROOT::VecOps::RVec<bool>& maxEpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && trackpass[i] && maxEpass[i]) return beta[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_beta", "REC_Particle_pass_std", "REC_Photon_MaxE_std"})
             .Define("pro_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass_std"})
             .Define("pro_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass_std"})
             .Define("pro_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass_std"})
             .Define("recpro_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& vz, const ROOT::VecOps::RVec<bool>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && pass[i]) return vz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass_std"})
             .Filter([](float ex, float gx, float px) { return ex != -999 && gx != -999 && px != -999; }, {"ele_px", "pho_px", "pro_px"})
             .Define("recel_p", MomentumFunc, {"ele_px", "ele_py", "ele_pz"})
             .Define("recel_theta", ThetaFunc, {"ele_px", "ele_py", "ele_pz"})
             .Define("recel_phi", PhiFunc, {"ele_px", "ele_py"})
             .Define("recpho_p", MomentumFunc, {"pho_px", "pho_py", "pho_pz"})
             .Define("recpho_theta", ThetaFunc, {"pho_px", "pho_py", "pho_pz"})
             .Define("recpho_phi", PhiFunc, {"pho_px", "pho_py"})
             .Define("recpro_p", MomentumFunc, {"pro_px", "pro_py", "pro_pz"})
             .Define("recpro_theta", ThetaFunc, {"pro_px", "pro_py", "pro_pz"})
             .Define("recpro_phi", PhiFunc, {"pro_px", "pro_py"})
             .Define("pho_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<bool>& pass, const ROOT::VecOps::RVec<bool>& maxEpass) {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 22 && pass[i] && maxEpass[i]) {
                           int abs_status = std::abs(status[i]);
                           if (abs_status >= 1000 && abs_status < 2000)
                             return 0;  // FT
                           else if (abs_status >= 2000 && abs_status < 3000)
                             return 1;  // FD
                           else if (abs_status >= 4000 && abs_status < 5000)
                             return 2;  // CD
                           else
                             return -1;  // Unknown/Other
                         }
                       }
                       return -1;
                     },
                     {"REC_Particle_pid", "REC_Particle_status", "REC_Particle_pass_std", "REC_Photon_MaxE_std"})

             .Define("pro_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<bool>& pass) {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 2212 && pass[i]) {
                           int abs_status = std::abs(status[i]);
                           if (abs_status >= 1000 && abs_status < 2000)
                             return 0;  // FT (probably rare for protons)
                           else if (abs_status >= 2000 && abs_status < 3000)
                             return 1;  // FD
                           else if (abs_status >= 4000 && abs_status < 5000)
                             return 2;  // CD
                           else
                             return -1;
                         }
                       }
                       return -1;
                     },
                     {"REC_Particle_pid", "REC_Particle_status", "REC_Particle_pass_std"})
             .Define("ele_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<bool>& pass) {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 11 && pass[i]) {
                           int abs_status = std::abs(status[i]);
                           if (abs_status >= 1000 && abs_status < 2000)
                             return 0;  // FT (probably rare for protons)
                           else if (abs_status >= 2000 && abs_status < 3000)
                             return 1;  // FD
                           else if (abs_status >= 4000 && abs_status < 5000)
                             return 2;  // CD
                           else
                             return -1;
                         }
                       }
                       return -1;
                     },
                     {"REC_Particle_pid", "REC_Particle_status", "REC_Particle_pass_std"});

  df_ = define_DISCAT(df_, "Q2", &DISANAMath::GetQ2, beam_energy);
  df_ = define_DISCAT(df_, "xB", &DISANAMath::GetxB, beam_energy);
  df_ = define_DISCAT(df_, "t", &DISANAMath::GetT, beam_energy);
  df_ = define_DISCAT(df_, "phi", &DISANAMath::GetPhi, beam_energy);
  df_ = define_DISCAT(df_, "W", &DISANAMath::GetW, beam_energy);
  df_ = define_DISCAT(df_, "nu", &DISANAMath::GetNu, beam_energy);
  df_ = define_DISCAT(df_, "y", &DISANAMath::Gety, beam_energy);
  df_ = define_DISCAT(df_, "Mx2_ep", &DISANAMath::GetMx2_ep, beam_energy);
  df_ = define_DISCAT(df_, "Emiss", &DISANAMath::GetEmiss, beam_energy);
  df_ = define_DISCAT(df_, "PTmiss", &DISANAMath::GetPTmiss, beam_energy);
  df_ = define_DISCAT(df_, "Mx2_epg", &DISANAMath::GetMx2_epg, beam_energy);
  df_ = define_DISCAT(df_, "Mx2_eg", &DISANAMath::GetMx2_egamma, beam_energy);
  df_ = define_DISCAT(df_, "Theta_e_gamma", &DISANAMath::GetTheta_e_gamma, beam_energy);
  df_ = define_DISCAT(df_, "DeltaE", &DISANAMath::GetDeltaE, beam_energy);
  df_ = define_DISCAT(df_, "DeltaPhi", &DISANAMath::GetDeltaPhi, beam_energy);
  df_ = define_DISCAT(df_, "Theta_gamma_gamma", &DISANAMath::GetTheta_gamma_gamma, beam_energy);

  return df_;
}

ROOT::RDF::RNode Init2PhotonKinematics(ROOT::RDF::RNode df_, float beam_energy) {
  df_ = df_.Define("pho2_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& trackpass, const ROOT::VecOps::RVec<bool>& maxEpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && trackpass[i] && !maxEpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass_std", "REC_Photon_MaxE_std"})
             .Define("pho2_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<bool>& trackpass, const ROOT::VecOps::RVec<bool>& maxEpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && trackpass[i] && !maxEpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass_std", "REC_Photon_MaxE_std"})
             .Define("pho2_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<bool>& trackpass, const ROOT::VecOps::RVec<bool>& maxEpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && trackpass[i] && !maxEpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass_std", "REC_Photon_MaxE_std"})
             .Define("recpho2_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& vz, const ROOT::VecOps::RVec<bool>& pass, const ROOT::VecOps::RVec<bool>& maxEpass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && pass[i] && !maxEpass[i]) return vz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass_std", "REC_Photon_MaxE_std"})
             .Define("recpho2_beta",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& beta, const ROOT::VecOps::RVec<bool>& trackpass, const ROOT::VecOps::RVec<bool>& maxEpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && trackpass[i] && !maxEpass[i]) return beta[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_beta", "REC_Particle_pass_std", "REC_Photon_MaxE_std"})
             .Define("hasrecpho2", [](float px) { return (px != -999.0f); }, {"pho2_px"})
             .Define("recpho2_p", MomentumFunc, {"pho2_px", "pho2_py", "pho2_pz"})
             .Define("recpho2_theta", ThetaFunc, {"pho2_px", "pho2_py", "pho2_pz"})
             .Define("recpho2_phi", PhiFunc, {"pho2_px", "pho2_py"})
             .Define("pho2_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<bool>& pass, const ROOT::VecOps::RVec<bool>& maxEpass) {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 22 && pass[i] && !maxEpass[i]) {
                           int abs_status = std::abs(status[i]);
                           if (abs_status >= 1000 && abs_status < 2000)
                             return 0;  // FT
                           else if (abs_status >= 2000 && abs_status < 3000)
                             return 1;  // FD
                           else if (abs_status >= 4000 && abs_status < 5000)
                             return 2;  // CD
                           else
                             return -1;  // Unknown/Other
                         }
                       }
                       return -1;
                     },
                     {"REC_Particle_pid", "REC_Particle_status", "REC_Particle_pass_std", "REC_Photon_MaxE_std"});
  df_ = define_DISCAT_pi0(df_, "Mass_pi0", &DISANAMath::GetMass_pi0, beam_energy);
  df_ = define_DISCAT_pi0(df_, "Mx2_eppi0", &DISANAMath::GetMx2_eppi0, beam_energy);
  df_ = define_DISCAT_pi0(df_, "Emiss_pi0", &DISANAMath::GetEmiss_pi0, beam_energy);
  df_ = define_DISCAT_pi0(df_, "Mx2_ep_pi0", &DISANAMath::GetMx2_ep_pi0, beam_energy);
  df_ = define_DISCAT_pi0(df_, "Mx2_epi0", &DISANAMath::GetMx2_epi0, beam_energy);
  df_ = define_DISCAT_pi0(df_, "PTmiss_pi0", &DISANAMath::GetPTmiss_pi0, beam_energy);
  df_ = define_DISCAT_pi0(df_, "Theta_pi0pi0", &DISANAMath::GetTheta_pi0pi0, beam_energy);
  df_ = define_DISCAT_pi0(df_, "DeltaPhi_pi0", &DISANAMath::GetDeltaPhi_pi0, beam_energy);
  df_ = define_DISCAT_pi0(df_, "Theta_epho1", &DISANAMath::GetTheta_epho1, beam_energy);
  df_ = define_DISCAT_pi0(df_, "Theta_epho2", &DISANAMath::GetTheta_epho2, beam_energy);
  df_ = define_DISCAT_pi0(df_, "Theta_pho1pho2", &DISANAMath::GetTheta_pho1pho2, beam_energy);
  return df_;
}

//
ROOT::RDF::RNode RejectPi0TwoPhoton(ROOT::RDF::RNode df_, float beam_energy) {
  //df_ = Init2PhotonKinematics(df_, beam_energy);
  //df_ = DefineDVPi0Pass(df_);
  //df_ = df_.Filter("!DVPi0_pass", "Cut: reject pi0");
  return df_.Filter(
      [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<bool>& pass) {
        int e = 0, g = 0, p = 0;
        for (size_t i = 0; i < pid.size(); ++i) {
          if (!pass[i]) continue;
          if (pid[i] == 11)
            e++;
          else if (pid[i] == 22)
            g++;  // photon must NOT be from pi0
          else if (pid[i] == 2212)
            p++;
        }
        return (e == 1 && g == 1 && p == 1);
      },
      {"REC_Particle_pid", "REC_Particle_pass_std"}, "Cut: one good e, γ , p");
}
// pi-0 event selection cuts for single photon contaminations
ROOT::RDF::RNode SelectPi0Event(ROOT::RDF::RNode df_) {
  return df_.Filter(
      [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<bool>& pass) {
        int e = 0, g = 0, p = 0;
        bool result = false;
        for (size_t i = 0; i < pid.size(); ++i) {
          //if (!pass[i]) continue;
          if (pid[i] == 11 && pass[i]){
            e++;
          }
          else if (pid[i] == 22 && pass[i]){
            g++;
          }
          else if (pid[i] == 2212 && pass[i]){
            p++;
          }
        }
        //if(g==2) std::cout << "e = " << e << " g = " << g << " p = " << p << std::endl;
        result = (e == 1 && g >= 2 && p == 1 );  // at least one photon
        return result;  // at least one photon, 
      },
      {"REC_Particle_pid", "REC_Particle_pass_std"}, "Cut: one good e, γ (not π⁰-like), p");
}
// exclusivity cuts

// ------------------------
// small utils
// ------------------------
static inline std::string Trim(std::string s) {
  auto is_space = [](unsigned char c){ return std::isspace(c); };
  while (!s.empty() && is_space((unsigned char)s.front())) s.erase(s.begin());
  while (!s.empty() && is_space((unsigned char)s.back()))  s.pop_back();
  return s;
}

static inline std::vector<std::string> SplitCSVLine(const std::string& line) {
  // simple csv split (no quoted commas)
  std::vector<std::string> out;
  std::stringstream ss(line);
  std::string item;
  while (std::getline(ss, item, ',')) out.push_back(Trim(item));
  return out;
}

struct Win { double lo=0.0, hi=0.0; };
static inline bool InWin(double x, const Win& w) {
  return (x > w.lo && x < w.hi);  // keep strict like your original
}

// key = (var,tbin,cfg)
struct Key3 {
  std::string var;
  int tbin;
  int cfg;
  bool operator==(const Key3& o) const {
    return var==o.var && tbin==o.tbin && cfg==o.cfg;
  }
};
struct Key3Hash {
  std::size_t operator()(const Key3& k) const noexcept {
    std::size_t h1 = std::hash<std::string>{}(k.var);
    std::size_t h2 = std::hash<int>{}(k.tbin);
    std::size_t h3 = std::hash<int>{}(k.cfg);
    // combine
    std::size_t h = h1;
    h ^= (h2 + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2));
    h ^= (h3 + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2));
    return h;
  }
};

struct CutTableTDep {
  std::vector<double> tEdges;  // size = ntbin+1
  std::unordered_map<Key3, Win, Key3Hash> winmap;

  int NTbin() const { return (int)tEdges.size() - 1; }
};

static inline int FindBin(double x, const std::vector<double>& edges) {
  const int nb = (int)edges.size() - 1;
  if (nb <= 0) return -1;
  if (x < edges.front() || x >= edges.back()) return -1;
  for (int i = 0; i < nb; i++) {
    if (x >= edges[i] && x < edges[i+1]) return i;
  }
  return -1;
}

// your 3 configs
static inline int ConfigId(int pho_det_region, int pro_det_region) {
  if (pho_det_region==0 && pro_det_region==2) return 0;
  if (pho_det_region==1 && pro_det_region==1) return 1;
  if (pho_det_region==1 && pro_det_region==2) return 2;
  return -1;
}

static CutTableTDep LoadCutsTDepCSV(const std::string& path) {
  std::ifstream fin(path);
  if (!fin) throw std::runtime_error("LoadCutsTDepCSV: cannot open: " + path);

  CutTableTDep tab;
  std::string line;
  int lineno = 0;

  while (std::getline(fin, line)) {
    lineno++;
    line = Trim(line);
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    auto cols = SplitCSVLine(line);
    if (cols.empty()) continue;

    if (cols[0] == "tEdges") {
      if (cols.size() < 3) {
        throw std::runtime_error("CSV parse error: tEdges needs >=2 numbers at line " + std::to_string(lineno));
      }
      tab.tEdges.clear();
      for (size_t i = 1; i < cols.size(); i++) tab.tEdges.push_back(std::stod(cols[i]));
      for (size_t i = 1; i < tab.tEdges.size(); i++) {
        if (!(tab.tEdges[i] > tab.tEdges[i-1])) {
          throw std::runtime_error("CSV parse error: tEdges must be strictly increasing");
        }
      }
      continue;
    }

    // data row: var,tbin,cfg,lo,hi
    if (cols.size() != 5) {
      throw std::runtime_error("CSV parse error at line " + std::to_string(lineno) +
                               ": expect 5 columns var,tbin,cfg,lo,hi");
    }
    if (tab.tEdges.empty()) {
      throw std::runtime_error("CSV parse error: tEdges must appear before data rows");
    }

    Key3 k{cols[0], std::stoi(cols[1]), std::stoi(cols[2])};
    Win  w{std::stod(cols[3]), std::stod(cols[4])};

    if (k.cfg < 0 || k.cfg >= 3) {
      throw std::runtime_error("CSV parse error: cfg must be 0/1/2");
    }
    if (k.tbin < 0 || k.tbin >= tab.NTbin()) {
      throw std::runtime_error("CSV parse error: tbin out of range (0.." + std::to_string(tab.NTbin()-1) + ")");
    }
    if (!(w.hi > w.lo)) {
      throw std::runtime_error("CSV parse error: invalid window hi<=lo for var=" + k.var);
    }

    tab.winmap[k] = w;
  }

  if (tab.tEdges.empty()) throw std::runtime_error("CSV missing tEdges line");

  return tab;
}

static inline void ScaleCutWindows(CutTableTDep& tab, double factor) {
  if (!std::isfinite(factor) || factor < 0.0) {
    throw std::invalid_argument(
        "event-selection window scale factor must be finite and non-negative");
  }
  for (auto& entry : tab.winmap) {
    Win& window = entry.second;
    const double center = 0.5 * (window.lo + window.hi);
    const double halfWidth = 0.5 * factor * (window.hi - window.lo);
    window.lo = center - halfWidth;
    window.hi = center + halfWidth;
  }
}

static inline bool PassVarTDep(const CutTableTDep& tab,
                               const std::string& var,
                               double x, int tbin, int cfg) {
  auto it = tab.winmap.find(Key3{var,tbin,cfg});
  if (it == tab.winmap.end()) return false;
  return InWin(x, it->second);
}

// ------------------------
// core selection: shared by REC/GEN
// ------------------------
static inline ROOT::RDF::RNode ApplyFinalDVCSSelections_TDep(ROOT::RDF::RNode df,
                                                            std::shared_ptr<CutTableTDep> cuts) {
  auto d0 = df
    // base kinematics
    //.Filter([](double p_out_theta, int pro_det) { return (pro_det == 1 && p_out_theta < 40.00 * M_PI/180) || (p_out_theta > 46.50 * M_PI/180 && p_out_theta < 64.23 * M_PI/180); }, {"recpro_theta", "pro_det_region"})
    .Filter("Q2 > 1.0",  "Cut: Q2 > 1 GeV^2")
    //.Filter("t < 1.0",   "Cut: t < 1 GeV^2")     // 你也可以删掉它，完全依赖 tEdges
    .Filter("W > 2.0",   "Cut: W > 2.0 GeV")
    .Filter("recpho_p > 2.0", "Cut: recpho_p > 2.0 GeV")
    .Filter("ele_det_region == 1", "Cut: electron in FD")

    // loose pre-cuts (same as your original)
    .Filter("Mx2_ep > -1.5 && Mx2_ep < 1.5", "Cut: MM^2(ep) in 3sigma (loose)")
    .Filter("Emiss < 1.0",                  "Cut: Missing energy (loose)")
    .Filter("PTmiss < 0.25",                "Cut: Transverse missing momentum (loose)")
    .Filter("Theta_e_gamma > 5",            "Cut: Theta_e_gamma (loose)")
    .Filter("Theta_gamma_gamma < 3.0",      "Cut: photon-missing angle (loose)")

    // labels
    .Define("cfg",  [](int pho, int pro){ return ConfigId(pho, pro); },
            {"pho_det_region","pro_det_region"})
    .Define("tbin", [cuts](double tt){ return FindBin(tt, cuts->tEdges); }, {"t"})
    .Filter("cfg >= 0",  "Cut: three config")
    .Filter("tbin >= 0", "Cut: t in tEdges");

  // cfg+tbin dependent windows
  auto d1 = d0
    .Define("pass_Mx2_ep", [cuts](double v, int tb, int cg){
        return PassVarTDep(*cuts, "Mx2_ep", v, tb, cg);
      }, {"Mx2_ep","tbin","cfg"})
    .Filter("pass_Mx2_ep", "Cut: Mx2_ep in 3sigma (tbin+cfg)")

    .Define("pass_Emiss", [cuts](double v, int tb, int cg){
        return PassVarTDep(*cuts, "Emiss", v, tb, cg);
      }, {"Emiss","tbin","cfg"})
    .Filter("pass_Emiss", "Cut: Emiss in 3sigma (tbin+cfg)")

    .Define("pass_PTmiss", [cuts](double v, int tb, int cg){
        return PassVarTDep(*cuts, "PTmiss", v, tb, cg);
      }, {"PTmiss","tbin","cfg"})
    .Filter("pass_PTmiss", "Cut: PTmiss in 3sigma (tbin+cfg)")

    .Define("pass_Theta_gamma_gamma", [cuts](double v, int tb, int cg){
        return PassVarTDep(*cuts, "Theta_gamma_gamma", v, tb, cg);
      }, {"Theta_gamma_gamma","tbin","cfg"})
    .Filter("pass_Theta_gamma_gamma", "Cut: Theta_gamma_gamma in 3sigma (tbin+cfg)")

    .Define("pass_DeltaPhi", [cuts](double v, int tb, int cg){
        return PassVarTDep(*cuts, "DeltaPhi", v, tb, cg);
      }, {"DeltaPhi","tbin","cfg"})
    .Filter("pass_DeltaPhi", "Cut: DeltaPhi in 3sigma (tbin+cfg)")

    .Define("pass_Mx2_epg", [cuts](double v, int tb, int cg){
        return PassVarTDep(*cuts, "Mx2_epg", v, tb, cg);
      }, {"Mx2_epg","tbin","cfg"})
    .Filter("pass_Mx2_epg", "Cut: Mx2_epg in 3sigma (tbin+cfg)")

    .Define("pass_Mx2_eg", [cuts](double v, int tb, int cg){
        return PassVarTDep(*cuts, "Mx2_eg", v, tb, cg);
      }, {"Mx2_eg","tbin","cfg"})
    .Filter("pass_Mx2_eg", "Cut: Mx2_eg in 3sigma (tbin+cfg)");

  return d1;
}

// ------------------------
// public APIs: REC/GEN
// ------------------------
ROOT::RDF::RNode ApplyFinalDVCSSelections(ROOT::RDF::RNode df,
                                         const std::string& rec_csv,
                                         double windowScaleFactor) {
  auto cuts = std::make_shared<CutTableTDep>(LoadCutsTDepCSV(rec_csv));
  ScaleCutWindows(*cuts, windowScaleFactor);
  return ApplyFinalDVCSSelections_TDep(df, cuts);
}

ROOT::RDF::RNode ApplyFinalGenDVCSSelections(ROOT::RDF::RNode df,
                                            const std::string& gen_csv,
                                            double windowScaleFactor) {
  auto cuts = std::make_shared<CutTableTDep>(LoadCutsTDepCSV(gen_csv));
  ScaleCutWindows(*cuts, windowScaleFactor);
  return ApplyFinalDVCSSelections_TDep(df, cuts);
}


ROOT::RDF::RNode ApplyFinalDVCSRadSelections(ROOT::RDF::RNode df) {
  return df;
  //return df
  //  .Filter("Emiss < 1.0", "Cut: Missing energy")
  //  .Filter("Mx2_ep > -0.20 && Mx2_ep < 0.20", "Cut: MM^2(ep)")
  //  .Filter("Theta_e_gamma > 5 ", "Cut: Theta_e_gamma");
}

ROOT::RDF::RNode DefineDVPi0Pass(ROOT::RDF::RNode df){
  return df.Define("DVPi0_pass",
      [](bool& haspho2, double& mass_pi0, double& mx2_eppi0, double& emiss_pi0, double& mx2_ep_pi0, double& mx2_epi0, double& ptmiss_pi0, double& theta_pi0pi0, double& deltaphi_pi0,
              int& pho_det_region, int& pho2_det_region, double& recpho_p, double& recpho2_p, int& pro_det_region,
              double& Q2, double& t, double& W,
              double& Theta_epho1, double& Theta_epho2, double& Theta_pho1pho2) {
        bool pass = false;
        if (haspho2 && recpho_p > 2.0 && recpho2_p > 0.4 && Q2 > 1.0 && t < 2.0 && W > 2.0) {
          if (pho_det_region == 0 && pho2_det_region == 0 && pro_det_region ==2) {
            pass = Inrange(emiss_pi0, -0.4, 0.4);
            pass = pass && Inrange(mx2_epi0, 0.5, 2.0);
            pass = pass && Inrange(ptmiss_pi0, 0.0, 0.2);
            //uncorr ver.
            //pass = pass && Inrange(mx2_eppi0, -0.014, 0.010);
            //pass = pass && Inrange(mx2_ep_pi0, -0.185, 0.299);
            //pass = pass && Inrange(mx2_epi0, 0.522, 1.434);
            //pass = pass && Inrange(emiss_pi0, -0.219, 0.321);
            //pass = pass && Inrange(deltaphi_pi0, 0.0, 3.969);
            //pass = pass && Inrange(ptmiss_pi0, 0.0, 0.091);
            //pass = pass && Inrange(theta_pi0pi0, 0.0, 1.472);
            //pass = pass && Inrange(mass_pi0, 0.09, 0.178);
            //corr ver.
            pass = pass && InrangeVaried(mx2_eppi0, -0.014, 0.010);
            pass = pass && InrangeVaried(mx2_ep_pi0, -0.177, 0.323);
            pass = pass && InrangeVaried(mx2_epi0, 0.459, 1.367);
            pass = pass && InrangeVaried(emiss_pi0, -0.251, 0.285);
            pass = pass && InrangeVaried(deltaphi_pi0, 0.0, 4.094);
            pass = pass && InrangeVaried(ptmiss_pi0, 0.0, 0.088);
            pass = pass && InrangeVaried(theta_pi0pi0, 0.0, 1.288);
            pass = pass && InrangeVaried(mass_pi0, 0.091, 0.179);
          } else if (pho_det_region == 1 && pho2_det_region == 1 && pro_det_region ==1) {
            pass = Inrange(Theta_epho1, 20.0, 999.0);
            pass = pass && Inrange(Theta_epho2, 10.0, 999.0);
            pass = pass && Inrange(Theta_pho1pho2, 2.0, 999.0);
            pass = pass && Inrange(ptmiss_pi0, 0.0, 0.2);
            //uncorr ver.
            //pass = pass && Inrange(mx2_eppi0, -0.019, 0.017);
            //pass = pass && Inrange(mx2_ep_pi0, -0.203, 0.401);
            //pass = pass && Inrange(mx2_epi0, 0.458, 1.598);
            //pass = pass && Inrange(emiss_pi0, -0.3, 0.544);
            //pass = pass && Inrange(deltaphi_pi0, 0.0, 8.796);
            //pass = pass && Inrange(ptmiss_pi0, 0.0, 0.149);
            //pass = pass && Inrange(theta_pi0pi0, 0.0, 2.107);
            //pass = pass && Inrange(mass_pi0, 0.113, 0.153);
            //corr ver.
            pass = pass && InrangeVaried(mx2_eppi0, -0.021, 0.015);
            pass = pass && InrangeVaried(mx2_ep_pi0, -0.227, 0.449);
            pass = pass && InrangeVaried(mx2_epi0, 0.181, 1.385);
            pass = pass && InrangeVaried(emiss_pi0, -0.516, 0.392);
            pass = pass && InrangeVaried(deltaphi_pi0, 0.0, 9.051);
            pass = pass && InrangeVaried(ptmiss_pi0, 0.0, 0.165);
            pass = pass && InrangeVaried(theta_pi0pi0, 0.0, 2.132);
            pass = pass && InrangeVaried(mass_pi0, 0.116, 0.160);
          } else if (pho_det_region == 1 && pho2_det_region == 1 && pro_det_region ==2) {
            pass = Inrange(Theta_epho1, 10.0, 999.0);
            pass = pass && Inrange(Theta_epho2, 10.0, 999.0);
            pass = pass && Inrange(Theta_pho1pho2, 3.0, 999.0);
            pass = pass && Inrange(ptmiss_pi0, 0.0, 0.2);
            //uncorr ver.
            //pass = pass && Inrange(mx2_eppi0, -0.017, 0.015);
            //pass = pass && Inrange(mx2_ep_pi0, -0.165, 0.351);
            //pass = pass && Inrange(mx2_epi0, 0.324, 1.696);
            //pass = pass && Inrange(emiss_pi0, -0.383, 0.565);
            //pass = pass && Inrange(deltaphi_pi0, 0.0, 8.065);
            //pass = pass && Inrange(ptmiss_pi0, 0.0, 0.120);
            //pass = pass && Inrange(theta_pi0pi0, 0.0, 2.051);
            //pass = pass && Inrange(mass_pi0, 0.114, 0.154);
            //corr ver.
            pass = pass && InrangeVaried(mx2_eppi0, -0.019, 0.013);
            pass = pass && InrangeVaried(mx2_ep_pi0, -0.187, 0.393);
            pass = pass && InrangeVaried(mx2_epi0, -0.012, 1.464);
            pass = pass && InrangeVaried(emiss_pi0, -0.598, 0.406);
            pass = pass && InrangeVaried(deltaphi_pi0, 0.0, 7.925);
            pass = pass && InrangeVaried(ptmiss_pi0, 0.0, 0.142);
            pass = pass && InrangeVaried(theta_pi0pi0, 0.0, 1.745);
            pass = pass && InrangeVaried(mass_pi0, 0.117, 0.161);
          }
        }
        return pass;
      },
      {"hasrecpho2", "Mass_pi0", "Mx2_eppi0", "Emiss_pi0", "Mx2_ep_pi0", "Mx2_epi0", "PTmiss_pi0", "Theta_pi0pi0", "DeltaPhi_pi0","pho_det_region","pho2_det_region", "recpho_p", "recpho2_p", "pro_det_region","Q2","t","W", "Theta_epho1", "Theta_epho2", "Theta_pho1pho2"});
}

ROOT::RDF::RNode DefineGenDVPi0Pass(ROOT::RDF::RNode df){
  return df.Define("DVPi0_pass",
      [](bool& haspho2, double& mass_pi0, double& mx2_eppi0, double& emiss_pi0, double& mx2_ep_pi0, double& mx2_epi0, double& ptmiss_pi0, double& theta_pi0pi0, double& deltaphi_pi0,
              int& pho_det_region, int& pho2_det_region, double& recpho_p, double& recpho2_p, int& pro_det_region,
              double& Q2, double& t, double& W,
              double& Theta_epho1, double& Theta_epho2, double& Theta_pho1pho2) {
        bool pass = false;
        if (haspho2 && recpho_p > 2.0 && recpho2_p > 0.4 && Q2 > 1.0 && t < 2.0 && W > 2.0) {
          if (pho_det_region == 0 && pho2_det_region == 0 && pro_det_region ==2) {
            pass = Inrange(emiss_pi0, -0.4, 0.4);
            pass = pass && Inrange(mx2_epi0, 0.5, 2.0);
            pass = pass && Inrange(ptmiss_pi0, 0.0, 0.2);
            //uncorr ver.
            //pass = pass && Inrange(mx2_eppi0, -0.007, 0.005);
            //pass = pass && Inrange(mx2_ep_pi0, -0.126, 0.174);
            //pass = pass && Inrange(mx2_epi0, 0.578, 1.214);
            //pass = pass && Inrange(emiss_pi0, -0.180, 0.192);
            //pass = pass && Inrange(deltaphi_pi0, 0.0, 1.590);
            //pass = pass && Inrange(ptmiss_pi0, 0.0, 0.044);
            //pass = pass && Inrange(theta_pi0pi0, 0.0, 0.689);
            //pass = pass && Inrange(mass_pi0, 0.130, 0.142);
            //corr ver.
            pass = pass && InrangeVaried(mx2_eppi0, -0.007, 0.005);
            pass = pass && InrangeVaried(mx2_ep_pi0, -0.134, 0.194);
            pass = pass && InrangeVaried(mx2_epi0, 0.585, 1.205);
            pass = pass && InrangeVaried(emiss_pi0, -0.172, 0.188);
            pass = pass && InrangeVaried(deltaphi_pi0, 0.0, 2.137);
            pass = pass && InrangeVaried(ptmiss_pi0, 0.0, 0.047);
            pass = pass && InrangeVaried(theta_pi0pi0, 0.0, 0.63);
            pass = pass && InrangeVaried(mass_pi0, 0.128, 0.144);
          } else if (pho_det_region == 1 && pho2_det_region == 1 && pro_det_region ==1) {
            pass = Inrange(Theta_epho1, 20.0, 999.0);
            pass = pass && Inrange(Theta_epho2, 10.0, 999.0);
            pass = pass && Inrange(Theta_pho1pho2, 2.0, 999.0);
            pass = pass && Inrange(ptmiss_pi0, 0.0, 0.2);
            //uncorr ver.
            //pass = pass && Inrange(mx2_eppi0, -0.016, 0.012);
            //pass = pass && Inrange(mx2_ep_pi0, -0.190, 0.306);
            //pass = pass && Inrange(mx2_epi0, 0.369, 1.417);
            //pass = pass && Inrange(emiss_pi0, -0.365, 0.399);
            //pass = pass && Inrange(deltaphi_pi0, 0.0, 6.577);
            //pass = pass && Inrange(ptmiss_pi0, 0.0, 0.129);
            //pass = pass && Inrange(theta_pi0pi0, 0.0, 1.538);
            //pass = pass && Inrange(mass_pi0, 0.117, 0.153);
            //corr ver.
            pass = pass && InrangeVaried(mx2_eppi0, -0.018, 0.014);
            pass = pass && InrangeVaried(mx2_ep_pi0, -0.229, 0.347);
            pass = pass && InrangeVaried(mx2_epi0, 0.340, 1.372);
            pass = pass && InrangeVaried(emiss_pi0, -0.412, 0.392);
            pass = pass && InrangeVaried(deltaphi_pi0, 0.0, 7.743);
            pass = pass && InrangeVaried(ptmiss_pi0, 0.0, 0.130);
            pass = pass && InrangeVaried(theta_pi0pi0, 0.0, 1.438);
            pass = pass && InrangeVaried(mass_pi0, 0.116, 0.156);
          } else if (pho_det_region == 1 && pho2_det_region == 1 && pro_det_region ==2) {
            pass = Inrange(Theta_epho1, 10.0, 999.0);
            pass = pass && Inrange(Theta_epho2, 10.0, 999.0);
            pass = pass && Inrange(Theta_pho1pho2, 3.0, 999.0);
            pass = pass && Inrange(ptmiss_pi0, 0.0, 0.2);
            //uncorr ver.
            //pass = pass && Inrange(mx2_eppi0, -0.011, 0.009);
            //pass = pass && Inrange(mx2_ep_pi0, -0.131, 0.209);
            //pass = pass && Inrange(mx2_epi0, 0.315, 1.515);
            //pass = pass && Inrange(emiss_pi0, -0.374, 0.422);
            //pass = pass && Inrange(deltaphi_pi0, 0.0, 5.012);
            //pass = pass && Inrange(ptmiss_pi0, 0.0, 0.093);
            //pass = pass && Inrange(theta_pi0pi0, 0.0, 1.241);
            //pass = pass && Inrange(mass_pi0, 0.117, 0.153);
            //corr ver.
            pass = pass && InrangeVaried(mx2_eppi0, -0.011, 0.009);
            pass = pass && InrangeVaried(mx2_ep_pi0, -0.149, 0.231);
            pass = pass && InrangeVaried(mx2_epi0, 0.228, 1.532);
            pass = pass && InrangeVaried(emiss_pi0, -0.427, 0.433);
            pass = pass && InrangeVaried(deltaphi_pi0, 0.0, 5.471);
            pass = pass && InrangeVaried(ptmiss_pi0, 0.0, 0.097);
            pass = pass && InrangeVaried(theta_pi0pi0, 0.0, 1.128);
            pass = pass && InrangeVaried(mass_pi0, 0.116, 0.156);
          }
        }
        return pass;
      },
      {"hasrecpho2", "Mass_pi0", "Mx2_eppi0", "Emiss_pi0", "Mx2_ep_pi0", "Mx2_epi0", "PTmiss_pi0", "Theta_pi0pi0", "DeltaPhi_pi0","pho_det_region","pho2_det_region", "recpho_p", "recpho2_p", "pro_det_region","Q2","t","W", "Theta_epho1", "Theta_epho2", "Theta_pho1pho2"});
}

ROOT::RDF::RNode ApplyFinalDVPi0Selections(ROOT::RDF::RNode df) {
  df = df.Filter("Q2 > 1.0", "Cut: Q2 > 1 GeV^2")
         //.Filter("t < 1.0", "Cut: t < 1 GeV^2")
      //.Filter("recel_p > 6.0", "Cut: recel_p > 0.6")

      // 5. W > 2
         .Filter("W > 2.0", "Cut: W > 1.8 GeV");
      //.Filter("phi > 100.0 && phi < 300 ", "Cut: phi")
  df = DefineDVPi0Pass(df);
  return df//.Filter([](double p_out_theta, int pro_det) { return (pro_det == 1 && p_out_theta < 40.00 * M_PI/180) || (p_out_theta > 46.50 * M_PI/180 && p_out_theta < 64.23 * M_PI/180); }, {"recpro_theta", "pro_det_region"})
           .Filter("DVPi0_pass", "Cut: DVPi0 event selection");
}

ROOT::RDF::RNode ApplyFinalGenDVPi0Selections(ROOT::RDF::RNode df) {
  df = df.Filter("Q2 > 1.0", "Cut: Q2 > 1 GeV^2")
         //.Filter("t < 1.0", "Cut: t < 1 GeV^2")
      //.Filter("recel_p > 6.0", "Cut: recel_p > 0.6")

      // 5. W > 2
         .Filter("W > 2.0", "Cut: W > 1.8 GeV");
      //.Filter("phi > 100.0 && phi < 300 ", "Cut: phi")
  df = DefineGenDVPi0Pass(df);
  return df//.Filter([](double p_out_theta, int pro_det) { return (pro_det == 1 && p_out_theta < 40.00 * M_PI/180) || (p_out_theta > 46.50 * M_PI/180 && p_out_theta < 64.23 * M_PI/180); }, {"recpro_theta", "pro_det_region"})
           .Filter("DVPi0_pass", "Cut: DVPi0 event selection");
}


ROOT::RDF::RNode InitGenKinematics(const std::string& filename_, const std::string& treename_, float beam_energy) {
  ROOT::RDataFrame rdf(treename_, filename_);
  auto df_ = std::make_unique<ROOT::RDF::RNode>(rdf);
  *df_ = df_->Define("ele_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11) return px[i];
                       return -999.0f;
                     },
                     {"MC_Particle_pid", "MC_Particle_px"})
             .Define("ele_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11) return py[i];
                       return -999.0f;
                     },
                     {"MC_Particle_pid", "MC_Particle_py"})
             .Define("ele_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11) return pz[i];
                       return -999.0f;
                     },
                     {"MC_Particle_pid", "MC_Particle_pz"})
             .Define("pho_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<float>& pz) {
                      for (size_t i = 0; i < pid.size(); ++i){
                         if (pid[i]==22) return px[i];
                      }
                         //if (pid[i] == 22 /*&& maxEpass[i]*/) return px[i];
                       return -999.0f;
                     },
                     {"MC_Particle_pid", "MC_Particle_px","MC_Particle_py","MC_Particle_pz"})
             .Define("pho_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<float>& pz) {
                      for (size_t i = 0; i < pid.size(); ++i){
                         if (pid[i]==22) return py[i];
                       }
                         //if (pid[i] == 22 /*&& maxEpass[i]*/) return px[i];
                       return -999.0f;
                     },
                     {"MC_Particle_pid", "MC_Particle_px","MC_Particle_py","MC_Particle_pz"})
             .Define("pho_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<float>& pz) {
                      float temp_p = 0.0; 
                      float result = -999.0f;
                      for (size_t i = 0; i < pid.size(); ++i){
                         if (pid[i]==22) return pz[i];
                       }
                         //if (pid[i] == 22 /*&& maxEpass[i]*/) return px[i];
                       return -999.0f;
                     },
                     {"MC_Particle_pid", "MC_Particle_px","MC_Particle_py","MC_Particle_pz"})
             .Define("pro_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 ) return px[i];
                       return -999.0f;
                     },
                     {"MC_Particle_pid", "MC_Particle_px"})
             .Define("pro_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 ) return py[i];
                       return -999.0f;
                     },
                     {"MC_Particle_pid", "MC_Particle_py"})
             .Define("pro_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 ) return pz[i];
                       return -999.0f;
                     },
                     {"MC_Particle_pid", "MC_Particle_pz"})
             .Filter([](float ex, float gx, float px) { return ex != -999 && gx != -999 && px != -999; }, {"ele_px", "pho_px", "pro_px"})
             .Define("recel_p", MomentumFunc, {"ele_px", "ele_py", "ele_pz"})
             .Define("recel_theta", ThetaFunc, {"ele_px", "ele_py", "ele_pz"})
             .Define("recel_phi", PhiFunc, {"ele_px", "ele_py"})
             .Define("recpho_p", MomentumFunc, {"pho_px", "pho_py", "pho_pz"})
             .Define("recpho_theta", ThetaFunc, {"pho_px", "pho_py", "pho_pz"})
             .Define("recpho_phi", PhiFunc, {"pho_px", "pho_py"})
             .Define("recpro_p", MomentumFunc, {"pro_px", "pro_py", "pro_pz"})
             .Define("recpro_theta", ThetaFunc, {"pro_px", "pro_py", "pro_pz"})
             .Define("recpro_phi", PhiFunc, {"pro_px", "pro_py"})
             .Define("recel_vz", [](){ return 0.0; })
             .Define("recpho_vz", [](){ return 0.0; })
             .Define("recpro_vz", [](){ return 0.0; })
             .Define("pho_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid) {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 22  /*&& maxEpass[i]*/) {
                           return 1;
                         }
                       }
                       return -1;
                     },
                     {"MC_Particle_pid"})

             .Define("pro_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid) {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 2212 ) {
                           return 1;
                         }
                       }
                       return -1;
                     },
                     {"MC_Particle_pid"})
             .Define("ele_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid) {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 11 ) {
                           return 1;
                         }
                       }
                       return -1;
                     },
                     {"MC_Particle_pid"});

  *df_ = define_DISCAT(*df_, "Q2", &DISANAMath::GetQ2, beam_energy);
  *df_ = define_DISCAT(*df_, "xB", &DISANAMath::GetxB, beam_energy);
  *df_ = define_DISCAT(*df_, "t", &DISANAMath::GetT, beam_energy);
  *df_ = define_DISCAT(*df_, "phi", &DISANAMath::GetPhi, beam_energy);
  *df_ = define_DISCAT(*df_, "W", &DISANAMath::GetW, beam_energy);
  *df_ = define_DISCAT(*df_, "nu", &DISANAMath::GetNu, beam_energy);
  *df_ = define_DISCAT(*df_, "y", &DISANAMath::Gety, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_ep", &DISANAMath::GetMx2_ep, beam_energy);
  *df_ = define_DISCAT(*df_, "Emiss", &DISANAMath::GetEmiss, beam_energy);
  *df_ = define_DISCAT(*df_, "PTmiss", &DISANAMath::GetPTmiss, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epg", &DISANAMath::GetMx2_epg, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_eg", &DISANAMath::GetMx2_egamma, beam_energy);
  *df_ = define_DISCAT(*df_, "Theta_e_gamma", &DISANAMath::GetTheta_e_gamma, beam_energy);
  *df_ = define_DISCAT(*df_, "DeltaE", &DISANAMath::GetDeltaE, beam_energy);
  *df_ = define_DISCAT(*df_, "DeltaPhi", &DISANAMath::GetDeltaPhi, beam_energy);
  *df_ = define_DISCAT(*df_, "Theta_gamma_gamma", &DISANAMath::GetTheta_gamma_gamma, beam_energy);

  return *df_;
}

ROOT::RDF::RNode GetSlim_exclusive(ROOT::RDF::RNode src, const std::string& f, const std::string& t, bool InitGen) {
  const bool fileExists = !gSystem->AccessPathName(f.c_str());  // note the '!' (exists == true)
  if (fileExists) {
    std::cout << "Slim file " << f << " exists, loading it." << std::endl;
    return ROOT::RDataFrame(t, f);
  } else {
    std::cout << "Trimming the file " << std::endl;
    return WriteSlimAndReload_exclusive(src, f, t, InitGen);
  }
}


ROOT::RDF::RNode WriteSlimAndReload_exclusive(ROOT::RDF::RNode df, const std::string& outFile, const std::string& outTree, bool InitGen) {
  // Keep EXACTLY these columns (update this list if you add/remove defs)
  const std::vector<std::string> keep_InitKinematics = {
    // ===== RAW REC BANKS =====
    "REC_Particle_pid",
    "REC_Particle_px",
    "REC_Particle_py",
    "REC_Particle_pz",
    "REC_Particle_vz",
    "REC_Particle_beta",
    "REC_Particle_status",
    "REC_Particle_pass_std",
    "REC_Photon_MaxE_std",
    "REC_Event_helicity",
    "RUN_config_run",
    //"REC_Track_pindex",
    //"REC_Track_sector",

    // ===== Picked particles =====
    "ele_px","ele_py","ele_pz",
    "pho_px","pho_py","pho_pz",
    "pro_px","pro_py","pro_pz",

    "recel_vz","recpho_vz","recpro_vz",
    "recpho_beta",

    // ===== Derived kinematics =====
    "recel_p","recel_theta","recel_phi",
    "recpho_p","recpho_theta","recpho_phi",
    "recpro_p","recpro_theta","recpro_phi",

    "ele_det_region","pho_det_region","pro_det_region",

    // ===== DISANAMath DVCS =====
    "Q2","xB","t","phi","W","nu","y",
    "Mx2_ep","Emiss","PTmiss",
    "Mx2_epg","Mx2_eg",
    "Theta_e_gamma","DeltaE","DeltaPhi","Theta_gamma_gamma"
  };

  const std::vector<std::string> keep_InitGenKinematics = {

    // ===== RAW GEN BANKS =====
    "MC_Particle_pid",
    "MC_Particle_px",
    "MC_Particle_py",
    "MC_Particle_pz",

    // ===== picked particles =====
    "ele_px","ele_py","ele_pz",
    "pho_px","pho_py","pho_pz",
    "pro_px","pro_py","pro_pz",

    "recel_vz","recpho_vz","recpro_vz",

    // ===== derived kinematics =====
    "recel_p","recel_theta","recel_phi",
    "recpho_p","recpho_theta","recpho_phi",
    "recpro_p","recpro_theta","recpro_phi",

    "ele_det_region","pho_det_region","pro_det_region",

    // ===== DISANAMath DVCS =====
    "Q2","xB","t","phi","W","nu","y",
    "Mx2_ep","Emiss","PTmiss",
    "Mx2_epg","Mx2_eg",
    "Theta_e_gamma","DeltaE","DeltaPhi","Theta_gamma_gamma"
  };

  if (InitGen){
    df.Snapshot(outTree, outFile, keep_InitGenKinematics);
  }
  else{
    df.Snapshot(outTree, outFile, keep_InitKinematics);
  }

  // Reload a much lighter dataframe
  ROOT::RDataFrame slim(outTree, outFile);
  return slim;  // implicitly converts to RNode
}
