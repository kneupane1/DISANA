#include <THnSparse.h>

#include "../DreamAN/DrawHist/DISANAcomparer.h"
#include "../DreamAN/DrawHist/DrawStyle.h"
#include "../DreamAN/DrawHist/DISANAMath.h"

#include <ROOT/RDataFrame.hxx>

#include <array>
#include <cctype>
#include <fstream>
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

ROOT::RDF::RNode ApplyFinalDVCSSelections(ROOT::RDF::RNode df, const std::string& rec_csv = "dvcs_cuts_rec_6p5.csv");
ROOT::RDF::RNode ApplyFinalDVCSRadSelections(ROOT::RDF::RNode df);
ROOT::RDF::RNode ApplyFinalGenDVCSSelections(ROOT::RDF::RNode df, const std::string& rec_csv = "dvcs_cuts_gen_6p5.csv");

ROOT::RDF::RNode DefineDVPi0Pass(ROOT::RDF::RNode df);
ROOT::RDF::RNode DefineGenDVPi0Pass(ROOT::RDF::RNode df);
ROOT::RDF::RNode ApplyFinalDVPi0Selections(ROOT::RDF::RNode df);
ROOT::RDF::RNode ApplyFinalGenDVPi0Selections(ROOT::RDF::RNode df);

ROOT::RDF::RNode InitKinematics(const std::string& filename_ = "", const std::string& treename_ = "", float beam_energy = 0);
ROOT::RDF::RNode Init2PhotonKinematics(ROOT::RDF::RNode df_, float beam_energy = 0);
ROOT::RDF::RNode InitGenKinematics(const std::string& filename_ = "", const std::string& treename_ = "", float beam_energy = 0);

ROOT::RDF::RNode GetSlim_exclusive(ROOT::RDF::RNode df_, const std::string& filename_slim, const std::string& treename_slim, bool isGen = false);
ROOT::RDF::RNode WriteSlimAndReload_exclusive(ROOT::RDF::RNode df_, const std::string& filename_slim, const std::string& treename_slim, bool isGen = false);

void PlotAllRecoDistributions(ROOT::RDF::RNode df, const std::string& out = "reco_kinematics_grid.png",
                              int bins_p = 120, int bins_theta = 120, int bins_phi = 120);

static double MomentumFunc(float px, float py, float pz) { return std::sqrt(px * px + py * py + pz * pz); }
static double ThetaFunc(float px, float py, float pz) { return std::acos(pz / std::sqrt(px * px + py * py + pz * pz)); }
static double PhiFunc(float px, float py) {
  double phi = std::atan2(py, px);
  return phi < 0 ? phi + 2 * M_PI : phi;
}

bool Inrange(double var, double min, double max) { return (var >= min && var < max); }

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

void DISANA_Xplotter2csv6p5() {
  bool ComputeBgk_core = false;  // Set to true if you want to compute background
  bool DoBkgCorr = true;       // Set to true if you want to apply background correction

  ROOT::EnableImplicitMT(40);
 
  std::string input_path_from_analysisRun_6535_data = "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File6535/data";
  std::string input_path_from_analysisRun_6535_pi0MC = "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File6535/pi0/";
  std::string input_path_from_analysisRun_6535_dvcsmc_rec = "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File6535/nobkg/";
  std::string input_path_from_analysisRun_6535_dvcsmc_gen = "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File6535/nobkgall/";
  std::string input_path_from_analysisRun_6535_bkg = "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File6535/bkg/";
  std::string input_path_from_analysisRun_6535_nobkg = "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File6535/nobkg/";

  std::string filename_afterFid_6535_data = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_6535_data.c_str());
  std::string filename_afterFid_6535_pi0MC = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_6535_pi0MC.c_str());
  std::string filename_afterFid_6535_dvcsmc_gen = Form("%s/dfSelected.root", input_path_from_analysisRun_6535_dvcsmc_gen.c_str());
  std::string filename_afterFid_6535_dvcsmc_rec = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_6535_dvcsmc_rec.c_str());
  std::string filename_afterFid_6535_bkg = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_6535_bkg.c_str());
  std::string filename_afterFid_6535_nobkg = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_6535_nobkg.c_str());
  std::string filename_afterFid_6535_dvcsmc_rad = "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File6535/radandP1gen/rad6535v0p4.root";
  std::string filename_afterFid_6535_dvcsmc_norad = "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File6535/radandP1gen/nor6535.root";
  std::string filename_afterFid_6535_dvcsmc_p1cut = "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File6535/radandP1gen/nor6535P1.root";


  float beam_energy = 6.535;

  ROOT::RDF::RNode df_afterFid_6535_data_init = InitKinematics(filename_afterFid_6535_data, "dfSelected_afterFid_afterCorr", beam_energy);
  ROOT::RDF::RNode df_afterFid_6535_pi0MC_init = InitKinematics(filename_afterFid_6535_pi0MC, "dfSelected_afterFid_afterCorr", beam_energy);
  ROOT::RDF::RNode df_afterFid_6535_dvcsmc_gen_init = InitGenKinematics(filename_afterFid_6535_dvcsmc_gen, "dfSelected", beam_energy);
  ROOT::RDF::RNode df_afterFid_6535_dvcsmc_rec_init = InitKinematics(filename_afterFid_6535_dvcsmc_rec, "dfSelected_afterFid_afterCorr", beam_energy);
  ROOT::RDF::RNode df_afterFid_6535_bkg_init = InitKinematics(filename_afterFid_6535_bkg, "dfSelected_afterFid_afterCorr", beam_energy);
  ROOT::RDF::RNode df_afterFid_6535_nobkg_init = InitKinematics(filename_afterFid_6535_nobkg, "dfSelected_afterFid_afterCorr", beam_energy);
  ROOT::RDF::RNode df_afterFid_6535_dvcsmc_rad_init = InitGenKinematics(filename_afterFid_6535_dvcsmc_rad, "MC", beam_energy);
  ROOT::RDF::RNode df_afterFid_6535_dvcsmc_norad_init = InitGenKinematics(filename_afterFid_6535_dvcsmc_norad, "MC", beam_energy);
  ROOT::RDF::RNode df_afterFid_6535_dvcsmc_p1cut_init = InitGenKinematics(filename_afterFid_6535_dvcsmc_p1cut, "MC", beam_energy);

  ROOT::RDF::RNode df_afterFid_6535_data = GetSlim_exclusive(df_afterFid_6535_data_init, "dfSlim_6535_data.root", "dfSlim_6535_data", false);
  ROOT::RDF::RNode df_afterFid_6535_pi0MC = GetSlim_exclusive(df_afterFid_6535_pi0MC_init, "dfSlim_6535_pi0MC.root", "dfSlim_6535_pi0MC", false);
  ROOT::RDF::RNode df_afterFid_6535_dvcsmc_gen = GetSlim_exclusive(df_afterFid_6535_dvcsmc_gen_init, "dfSlim_6535_dvcsmc_gen.root", "dfSlim_6535_dvcsmc_gen", true);
  ROOT::RDF::RNode df_afterFid_6535_dvcsmc_rec = GetSlim_exclusive(df_afterFid_6535_dvcsmc_rec_init, "dfSlim_6535_dvcsmc_rec.root", "dfSlim_6535_dvcsmc_rec", false);
  ROOT::RDF::RNode df_afterFid_6535_bkg = GetSlim_exclusive(df_afterFid_6535_bkg_init, "dfSlim_6535_bkg.root", "dfSlim_6535_bkg", false);
  ROOT::RDF::RNode df_afterFid_6535_nobkg = GetSlim_exclusive(df_afterFid_6535_nobkg_init, "dfSlim_6535_nobkg.root", "dfSlim_6535_nobkg", false);
  ROOT::RDF::RNode df_afterFid_6535_dvcsmc_rad = GetSlim_exclusive(df_afterFid_6535_dvcsmc_rad_init, "dfSlim_6535_dvcsmc_rad.root", "dfSlim_6535_dvcsmc_rad", true);
  ROOT::RDF::RNode df_afterFid_6535_dvcsmc_norad = GetSlim_exclusive(df_afterFid_6535_dvcsmc_norad_init, "dfSlim_6535_dvcsmc_norad.root", "dfSlim_6535_dvcsmc_norad", true);
  ROOT::RDF::RNode df_afterFid_6535_dvcsmc_p1cut = GetSlim_exclusive(df_afterFid_6535_dvcsmc_p1cut_init, "dfSlim_6535_dvcsmc_p1cut.root", "dfSlim_6535_dvcsmc_p1cut", true);

  DrawStyle fitStyle(0.06, 0.05, 1.0, 1.3);  // You can tweak this


  // Apply final DVCS cuts
  auto df_final_dvcs_6535_data = ApplyFinalDVCSSelections(df_afterFid_6535_data);
  auto df_final_dvcsPi_rejected_6535_data = RejectPi0TwoPhoton(df_final_dvcs_6535_data, beam_energy);

  auto df_final_dvcs_6535_pi0MC = ApplyFinalGenDVCSSelections(df_afterFid_6535_pi0MC);
  auto df_final_dvcsPi_rejected_6535_pi0MC = RejectPi0TwoPhoton(df_final_dvcs_6535_pi0MC, beam_energy);

  auto df_final_dvcs_6535_dvcsmc_rec = ApplyFinalGenDVCSSelections(df_afterFid_6535_dvcsmc_rec);
  auto df_final_dvcsPi_rejected_6535_dvcsmc_rec = RejectPi0TwoPhoton(df_final_dvcs_6535_dvcsmc_rec, beam_energy);

  auto df_final_dvcs_6535_dvcsmc_bkg = ApplyFinalGenDVCSSelections(df_afterFid_6535_bkg);
  auto df_final_dvcsPi_rejected_6535_dvcsmc_bkg = RejectPi0TwoPhoton(df_final_dvcs_6535_dvcsmc_bkg, beam_energy);

  auto df_final_dvcs_6535_dvcsmc_nobkg = ApplyFinalGenDVCSSelections(df_afterFid_6535_nobkg);
  auto df_final_dvcsPi_rejected_6535_dvcsmc_nobkg = RejectPi0TwoPhoton(df_final_dvcs_6535_dvcsmc_nobkg, beam_energy);

  auto df_final_OnlPi0_6535_data = ApplyFinalDVPi0Selections(Init2PhotonKinematics(SelectPi0Event(df_afterFid_6535_data), beam_energy));
  auto df_final_OnlPi0_6535_pi0MC = ApplyFinalGenDVPi0Selections(Init2PhotonKinematics(SelectPi0Event(df_afterFid_6535_pi0MC), beam_energy));

  DISANAcomparer comparer;
  comparer.SetOutputDir("./");
  comparer.SetKinStyle(KinStyle);
  comparer.SetDVCSStyle(dvcsStyle);
  comparer.SetCrossSectionStyle(csStyle);
  comparer.SetBSAStyle(bsaStyle);

  comparer.PlotIndividual(false);
  /// bins for cross-section plots
  BinManager xBins;

  xBins.SetQ2Bins({1.00, 1.25, 1.50, 1.75, 2.00, 2.40, 2.90, 3.50});
  xBins.SetTBins({0.13, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0});
  xBins.SetXBBins({0.150, 0.180, 0.210, 0.240, 0.285, 0.350, 0.430});

  comparer.SetXBinsRanges(xBins);

  //double charge=4.815525219658029+(8.88177914805192-0.2128897513862203)*0.5; // mC (5681-5757, 5757-5870(trigger prescale 2), 5758removed) 7.5GeV
  double charge=17.644; // mC 6.5GeV 5919 removed
  std::cout<<"Total effective charge (mC): "<<charge<<std::endl;
  double luminosity = (charge)*1.33*pow(10,6);  // Set your desired luminosity here nb^-1
  std::cout<<"Luminosity (nb^-1): "<<luminosity<<std::endl;
  double polarisation = 0.85;  // Set your desired polarisation here
  std::cout<<"Polarisation: "<<polarisation<<std::endl;

  /*comparer.AddModelwithPi0Corr(//df_final_dvcsPi_rejected_6535_data,
                              df_final_dvcsPi_rejected_6535_dvcsmc_rec,
                              df_final_OnlPi0_6535_pi0MC,
                              df_final_dvcsPi_rejected_6535_pi0MC,
                              df_final_OnlPi0_6535_pi0MC,
                              df_afterFid_6535_dvcsmc_gen,
                              df_final_dvcsPi_rejected_6535_dvcsmc_rec,
                              df_final_dvcsPi_rejected_6535_dvcsmc_bkg,
                              df_final_dvcsPi_rejected_6535_dvcsmc_nobkg,
                              df_afterFid_6535_dvcsmc_rad,
                              df_afterFid_6535_dvcsmc_norad,
                              df_afterFid_6535_dvcsmc_p1cut,
                              "RGK 6.5GeV GEMC", beam_energy, true, true, true, true, false, luminosity);*/
  comparer.AddModelwithPi0Corr(df_final_dvcsPi_rejected_6535_data,
                              //df_final_dvcsPi_rejected_6535_dvcsmc_rec,
                              df_final_OnlPi0_6535_data,
                              df_final_dvcsPi_rejected_6535_pi0MC,
                              df_final_OnlPi0_6535_pi0MC,
                              df_afterFid_6535_dvcsmc_gen,
                              df_final_dvcsPi_rejected_6535_dvcsmc_rec,
                              df_final_dvcsPi_rejected_6535_dvcsmc_bkg,
                              df_final_dvcsPi_rejected_6535_dvcsmc_nobkg,
                              df_afterFid_6535_dvcsmc_rad,
                              df_afterFid_6535_dvcsmc_norad,
                              df_afterFid_6535_dvcsmc_p1cut,
                              "RGK 6.5GeV", beam_energy, true, true, true, true, true, luminosity, 55.09, 60, 0.8280/0.9495);


  //comparer.PlotKinematicComparison();
  //comparer.PlotPi0KinematicComparison();
  //comparer.PlotxBQ2tBin();
  //comparer.PlotxBQ2tBinPi0();
  //comparer.PlotxBQ2tBinMC();
  //comparer.PlotxBQ2tBinPi0MC();
  //comparer.PlotDVCSKinematicsComparison();
  //comparer.PlotDVPi0KinematicsComparison();
  comparer.PlotDIS_BSA_Cross_Section_AndCorr_Comparison(polarisation, true, true, true, true, true, true, true, true);   
  //comparer.PlotDISCrossSectionComparison(luminosity);  // argument is Luminosity, polarisation
  //comparer.PlotDIS_BSA_Comparison(luminosity, polarisation);         // argument is Luminosity
  //comparer.PlotDIS_Pi0CorrComparison();
  //comparer.PlotPi0ExclusivityComparisonByDetectorCases(detCutsPi0);
  //comparer.PlotExclusivityComparisonByDetectorCases(detCuts);

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

static inline bool PassVarTDep(const CutTableTDep& tab,
                               const std::string& var,
                               double x, int tbin, int cfg) {
  auto it = tab.winmap.find(Key3{var,tbin,cfg});
  if (it == tab.winmap.end()) return false;
  return InWin(x, it->second);
}

static inline bool PassVarTDep_Debug(const CutTableTDep& tab,
                                     const std::string& var,
                                     double v, int tbin, int cfg) 
{
  auto it = tab.winmap.find(Key3{var, tbin, cfg});

  if (it == tab.winmap.end()) {
    std::cout << "[CutDebug] " << var
              << " cfg=" << cfg << " tbin=" << tbin
              << "  --> NO WINDOW FOUND (FAIL)" << std::endl;
    return false;
  }

  double lo = it->second.lo;
  double hi = it->second.hi;

  bool pass = (v >= lo && v <= hi);

  std::cout << "[ApplyCut] " << var
            << " (cfg=" << cfg << ", tbin=" << tbin << ") : "
            << lo << " < " << var << " < " << hi
            << "  | value=" << v
            << " --> " << (pass ? "PASS" : "FAIL")
            << std::endl;

  return pass;
}
// ------------------------
// core selection: shared by REC/GEN
// ------------------------
static inline ROOT::RDF::RNode ApplyFinalDVCSSelections_TDep(ROOT::RDF::RNode df,
                                                            std::shared_ptr<CutTableTDep> cuts) {
  auto d0 = df
    // base kinematics
    .Filter("Q2 > 1.0",  "Cut: Q2 > 1 GeV^2")
    .Filter("t < 1.0",   "Cut: t < 1 GeV^2")     // 你也可以删掉它，完全依赖 tEdges
    .Filter("W > 2.0",   "Cut: W > 2.0 GeV")
    .Filter("recpho_p > 2.0", "Cut: recpho_p > 2.0 GeV")
    .Filter("ele_det_region == 1", "Cut: electron in FD")
    .Filter("pho_det_region != 0", "Cut: photon not in FT")

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
                                         const std::string& rec_csv = "dvcs_cuts_rec_6p5.csv") {
  auto cuts = std::make_shared<CutTableTDep>(LoadCutsTDepCSV(rec_csv));
  return ApplyFinalDVCSSelections_TDep(df, cuts);
}

ROOT::RDF::RNode ApplyFinalGenDVCSSelections(ROOT::RDF::RNode df,
                                            const std::string& gen_csv = "dvcs_cuts_gen_6p5.csv") {
  auto cuts = std::make_shared<CutTableTDep>(LoadCutsTDepCSV(gen_csv));
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
        if (haspho2 && recpho_p > 2.0 && recpho2_p > 0.4 && Q2 > 1.0 && t < 1.0 && W > 2.0) {
          if (pho_det_region == 1 && pho2_det_region == 1 && pro_det_region ==1) {
            pass = Inrange(Theta_epho1, 20.0, 999.0);
            pass = pass && Inrange(Theta_epho2, 10.0, 999.0);
            pass = pass && Inrange(Theta_pho1pho2, 2.0, 999.0);
            pass = pass && Inrange(ptmiss_pi0, 0.0, 0.2);
            pass = pass && Inrange(mx2_eppi0, -0.019,0.017);
            pass = pass && Inrange(mx2_ep_pi0, -0.237,0.339);
            pass = pass && Inrange(mx2_epi0, -0.468,1.596);
            pass = pass && Inrange(emiss_pi0, -0.278,0.542);
            pass = pass && Inrange(deltaphi_pi0, 0,8.120);
            pass = pass && Inrange(ptmiss_pi0, 0,0.149);
            pass = pass && Inrange(theta_pi0pi0, 0,2.283);
            pass = pass && Inrange(mass_pi0, 0.11,154);
            //pass = Inrange(mass_pi0, 0.1, 0.16);
            //pass = pass && Inrange(emiss_pi0, -0.4, 0.4);
            //pass = pass && Inrange(ptmiss_pi0, 0.0, 0.2);
            //pass = pass && Inrange(theta_pi0pi0, 0.0, 1.5);
            //pass = pass && Inrange(deltaphi_pi0, 0.0, 8.0);
            //pass = pass && Inrange(mx2_eppi0, -0.02, 0.02);
            //pass = pass && Inrange(mx2_ep_pi0, -0.2, 0.2);
            //pass = pass && Inrange(mx2_epi0, 0.5, 1.5);
          } else if (pho_det_region == 1 && pho2_det_region == 1 && pro_det_region ==2) {
            pass = Inrange(Theta_epho1, 10.0, 999.0);
            pass = pass && Inrange(Theta_epho2, 10.0, 999.0);
            pass = pass && Inrange(Theta_pho1pho2, 3.0, 999.0);
            pass = pass && Inrange(ptmiss_pi0, 0.0, 0.2);
            pass = pass && Inrange(mx2_eppi0, -0.017,0.015);
            pass = pass && Inrange(mx2_ep_pi0, -0.155,0.329);
            pass = pass && Inrange(mx2_epi0, -0.360,1.684);
            pass = pass && Inrange(emiss_pi0, -0.357,0.559);
            pass = pass && Inrange(deltaphi_pi0, 0,8.543);
            pass = pass && Inrange(ptmiss_pi0, 0,0.12);
            pass = pass && Inrange(theta_pi0pi0, 0,2.12);
            pass = pass && Inrange(mass_pi0, 0.112,0.156);
            //pass = Inrange(mass_pi0, 0.11, 0.15);
            //pass = pass && Inrange(emiss_pi0, -0.4, 0.4);
            //pass = pass && Inrange(ptmiss_pi0, 0.0, 0.2);
            //pass = pass && Inrange(theta_pi0pi0, 0.0, 1.5);
            //pass = pass && Inrange(deltaphi_pi0, 0.0, 8.0);
            //pass = pass && Inrange(mx2_eppi0, -0.02, 0.02);
            //pass = pass && Inrange(mx2_ep_pi0, -0.2, 0.2);
            //pass = pass && Inrange(mx2_epi0, 0.5, 1.5);
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
        if (haspho2 && recpho_p > 2.0 && recpho2_p > 0.4 && Q2 > 1.0 && t < 1.0 && W > 2.0) {
          if (pho_det_region == 1 && pho2_det_region == 1 && pro_det_region ==1) {
            pass = Inrange(Theta_epho1, 20.0, 999.0);
            pass = pass && Inrange(Theta_epho2, 10.0, 999.0);
            pass = pass && Inrange(Theta_pho1pho2, 2.0, 999.0);
            pass = pass && Inrange(ptmiss_pi0, 0.0, 0.2);
            pass = pass && Inrange(mx2_eppi0, -0.016,0.012);
            pass = pass && Inrange(mx2_ep_pi0, -0.209,0.311);
            pass = pass && Inrange(mx2_epi0, -0.396,1.396);
            pass = pass && Inrange(emiss_pi0, -0.333,0.387);
            pass = pass && Inrange(deltaphi_pi0, 0,6.872);
            pass = pass && Inrange(ptmiss_pi0, 0,0.128);
            pass = pass && Inrange(theta_pi0pi0, 0,1.522);
            pass = pass && Inrange(mass_pi0, 0.115,155);
            //pass = Inrange(mass_pi0, 0.1, 0.16);
            //pass = pass && Inrange(emiss_pi0, -0.4, 0.4);
            //pass = pass && Inrange(ptmiss_pi0, 0.0, 0.2);
            //pass = pass && Inrange(theta_pi0pi0, 0.0, 1.5);
            //pass = pass && Inrange(deltaphi_pi0, 0.0, 8.0);
            //pass = pass && Inrange(mx2_eppi0, -0.02, 0.02);
            //pass = pass && Inrange(mx2_ep_pi0, -0.2, 0.2);
            //pass = pass && Inrange(mx2_epi0, 0.5, 1.5);
          } else if (pho_det_region == 1 && pho2_det_region == 1 && pro_det_region ==2) {
            pass = Inrange(Theta_epho1, 10.0, 999.0);
            pass = pass && Inrange(Theta_epho2, 10.0, 999.0);
            pass = pass && Inrange(Theta_pho1pho2, 3.0, 999.0);
            pass = pass && Inrange(ptmiss_pi0, 0.0, 0.2);
            pass = pass && Inrange(mx2_eppi0, -0.011,0.009);
            pass = pass && Inrange(mx2_ep_pi0, -0.135,0.213);
            pass = pass && Inrange(mx2_epi0, -0.339,1.499);
            pass = pass && Inrange(emiss_pi0, -0.358,0.414);
            pass = pass && Inrange(deltaphi_pi0, 0,5.65);
            pass = pass && Inrange(ptmiss_pi0, 0,0.093);
            pass = pass && Inrange(theta_pi0pi0, 0,1.328);
            pass = pass && Inrange(mass_pi0, 0.115,0.155);
            //pass = Inrange(mass_pi0, 0.11, 0.15);
            //pass = pass && Inrange(emiss_pi0, -0.4, 0.4);
            //pass = pass && Inrange(ptmiss_pi0, 0.0, 0.2);
            //pass = pass && Inrange(theta_pi0pi0, 0.0, 1.5);
            //pass = pass && Inrange(deltaphi_pi0, 0.0, 8.0);
            //pass = pass && Inrange(mx2_eppi0, -0.02, 0.02);
            //pass = pass && Inrange(mx2_ep_pi0, -0.2, 0.2);
            //pass = pass && Inrange(mx2_epi0, 0.5, 1.5);
          }
        }
        return pass;
      },
      {"hasrecpho2", "Mass_pi0", "Mx2_eppi0", "Emiss_pi0", "Mx2_ep_pi0", "Mx2_epi0", "PTmiss_pi0", "Theta_pi0pi0", "DeltaPhi_pi0","pho_det_region","pho2_det_region", "recpho_p", "recpho2_p", "pro_det_region","Q2","t","W", "Theta_epho1", "Theta_epho2", "Theta_pho1pho2"});
}

ROOT::RDF::RNode ApplyFinalDVPi0Selections(ROOT::RDF::RNode df) {
  //df = df.Filter("Q2 > 1.0", "Cut: Q2 > 1 GeV^2")
  //    .Filter("t < 1.0", "Cut: t < 1 GeV^2")
      //.Filter("recel_p > 6.0", "Cut: recel_p > 0.6")

      // 5. W > 2
  //    .Filter("W > 2.0", "Cut: W > 1.8 GeV");
      //.Filter("phi > 100.0 && phi < 300 ", "Cut: phi")
  df = DefineDVPi0Pass(df);
  return df.Filter("DVPi0_pass", "Cut: DVPi0 event selection");
}

ROOT::RDF::RNode ApplyFinalGenDVPi0Selections(ROOT::RDF::RNode df) {
  //df = df.Filter("Q2 > 1.0", "Cut: Q2 > 1 GeV^2")
  //    .Filter("t < 1.0", "Cut: t < 1 GeV^2")
      //.Filter("recel_p > 6.0", "Cut: recel_p > 0.6")

      // 5. W > 2
  //    .Filter("W > 2.0", "Cut: W > 1.8 GeV");
      //.Filter("phi > 100.0 && phi < 300 ", "Cut: phi")
  df = DefineGenDVPi0Pass(df);
  return df.Filter("DVPi0_pass", "Cut: DVPi0 event selection");
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
  //return src;
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
    "REC_Track_pindex",
    "REC_Track_sector",

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
