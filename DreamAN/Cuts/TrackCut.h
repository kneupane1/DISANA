#ifndef TRACKCUT_H_
#define TRACKCUT_H_

#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <TMath.h>

struct FiducialAxisCut {
  std::vector<std::pair<float, float>> excludedRanges;  // e.g., {{100, 120}, {240, 260}}
  std::set<float> excludedStrips;                       // e.g., {128.0, 256.0}
};

struct FiducialCut3D {
  FiducialAxisCut luCut;
  FiducialAxisCut lvCut;
  FiducialAxisCut lwCut;
};

struct SFCutABC {
  float A0;
  float Bm1;
  float Cm2;
};

class TrackCut {
 public:
  TrackCut();
  virtual ~TrackCut();
  // it is requited when copying the object itself
  TrackCut(const TrackCut& other);

  void SetSectorCut(int SSector, int selectpid, int selectdetector, bool selectSector);
  void SetPositionCut(float minX, float maxX, float minY, float maxY, float minZ, float maxZ);
  void SetDirectionCut(float minCX, float maxCX, float minCY, float maxCY, float minCZ, float maxCZ);
  void SetPathLengthCut(float minPath, float maxPath);
  void SetDCEdgeCut(float minEdge, float maxEdge);
  void SetECALEdgeCut(float minEdge, float maxEdge);
  void SetSectorCut_Bhawani(const std::vector<int>& sectors, int selectpid, int selectdetector, bool selectSector);
  void SetDoFiducialCut(bool isFiducial) { fDoFiducialCut = isFiducial; }
  void SetSFCut(bool doSFCut,int SFpid, float SFmin, float SFminP) { fdoSFCut = doSFCut; fSFpid = SFpid; fSFmin = SFmin; fSFminP = SFminP; }
  void SetFiducialCutOptions(bool doDC, bool doECAL) {
    fDoDCFiducial = doDC;
    fDoECALFiducial = doECAL;
  }
  void SetDCEdgeCut(float minTheta, float maxTheta, int dcRegion, float edgeCut);
  void SetThetaBins(const std::vector<std::pair<float, float>>& thetaBins);
  void SetDCEdgeCuts(int pid, const std::vector<float>& edgeCutsPerRegion);
  float GetEdgeCut(int pid, int region) const;
  void SetCVTEdgeCuts(int pid, const std::vector<float>& edgeCutsPerLayer);
  float GetCVTEdgeCut(int pid, int layer) const;
  const std::map<int, std::vector<float>>& GetEdgeCuts() const;
  const std::map<int, std::vector<float>>& GetCVTEdgeCuts() const;


  bool operator()(const std::vector<int16_t>& pindex, const std::vector<int16_t>& index, const std::vector<int>& detector, const std::vector<int>& layer,
                  const std::vector<float>& x, const std::vector<float>& y, const std::vector<float>& z, const std::vector<float>& cx, const std::vector<float>& cy,
                  const std::vector<float>& cz, const std::vector<float>& path, const std::vector<float>& edge) const;

  // DC / CVT filter — ForFiducialCut() columns only
  std::function<std::vector<int>(const std::vector<int16_t>& pindex,
                                 const std::vector<int16_t>& detector,
                                 const std::vector<int16_t>& layer,
                                 const std::vector<float>& x,
                                 const std::vector<float>& y,
                                 const std::vector<float>& z,
                                 const std::vector<float>& edge,
                                 const std::vector<int>& pid,
                                 const int& REC_Particle_num)>
  RECTrajPass() const;

  // Calorimeter filter — ForFiducialCut() columns only
  std::function<std::vector<int>(const std::vector<int16_t>&,  // pindex
                                 const std::vector<int16_t>&,  // detector
                                 const std::vector<int16_t>&,  // sector
                                 const std::vector<int16_t>&,  // layer
                                 const std::vector<float>&,    // energy
                                 const std::vector<float>&,    // lu
                                 const std::vector<float>&,    // lv
                                 const std::vector<float>&,    // lw
                                 const std::vector<int>&,      // pid
                                 const std::vector<float>&,    // p
                                 const int& REC_Particle_num)>
  RECCalorimeterPass() const;

  // Forward Tagger filter — ForFiducialCut() columns only
  std::function<std::vector<int>(const std::vector<short>&,    // pindex
                                 const std::vector<int16_t>&,  // detector
                                 const std::vector<int16_t>&,  // layer
                                 const std::vector<float>&,    // x
                                 const std::vector<float>&,    // y
                                 const std::vector<int>&,      // pid
                                 const int& REC_Particle_num)>
  RECForwardTaggerPass() const;


// Fiducial filter function combined DC and Calorimeter
  // This function will apply the fiducial cuts for both DC and Calorimeter
  // It will return a vector of integers indicating whether each track passes the fiducial cuts
  // The vector will have the same size as the number of particles in the event
  // A value of 1 indicates that the track passes the fiducial cuts, while a value of 0 indicates that it does not pass
std::function<std::vector<int>(
    // RECTraj
    const std::vector<int16_t>& traj_pindex,
    const std::vector<int16_t>& traj_index,
    const std::vector<int16_t>& traj_detector,
    const std::vector<int16_t>& traj_layer,
    const std::vector<float>& x, 
    const std::vector<float>& y, 
    const std::vector<float>& z, 
    const std::vector<float>& cx, 
    const std::vector<float>& cy,
    const std::vector<float>& cz, 
    const std::vector<float>& path,
    const std::vector<float>& traj_edge,
    // RECCalorimeter
    const std::vector<int16_t>& calo_pindex,
    const std::vector<int16_t>& calo_index,
    const std::vector<int16_t>& calo_detector,
    const std::vector<int16_t>& calo_sector,
    const std::vector<int16_t>& calo_layer,
    const std::vector<float>& calo_energy,
    const std::vector<float>& calo_time,
    const std::vector<float>& calo_path,
    const std::vector<float>& calo_chi2,
    const std::vector<float>& calo_x,
    const std::vector<float>& calo_y,
    const std::vector<float>& calo_z,
    const std::vector<float>& calo_hx,
    const std::vector<float>& calo_hy,
    const std::vector<float>& calo_hz,
    const std::vector<float>& calo_lu,
    const std::vector<float>& calo_lv,
    const std::vector<float>& calo_lw,
    const std::vector<float>& calo_du,
    const std::vector<float>& calo_dv,
    const std::vector<float>& calo_dw,
    const std::vector<float>& calo_m2u,
    const std::vector<float>& calo_m2v,
    const std::vector<float>& calo_m2w,
    const std::vector<float>& calo_m3u,
    const std::vector<float>& calo_m3v,
    const std::vector<float>& calo_m3w,
    const std::vector<short>& calo_status,
    const std::vector<int>& pid,
    const int& REC_Particle_num)>
    RECFiducialPass() const;
  // ===== New sector-specific fiducial setters =====
void AddCVTFiducialRange(int pid, int layer, const std::string& axis, float min, float max);
void AddFTCalFiducialRange(int pid, int layer, float x, float y, float rmin, float rmax);
void AddPCalFiducialRange(int pid, int sector, const std::string& axis, float min, float max);
void AddECinFiducialRange(int pid, int sector, const std::string& axis, float min, float max);
void AddECoutFiducialRange(int pid, int sector, const std::string& axis, float min, float max);
void SetMinECALEnergyCut(int pid, int layer, float minEnergy);

void AddSamplingFractionMinCut(int pid, int sector, float A0, float Bm1, float Cm2);
void AddSamplingFractionMaxCut(int pid, int sector, float A0, float Bm1, float Cm2);
void SetSamplingFractionSystematics(bool enabled, float threeSigmaFactor) {
  fPIDSystematics = enabled;
  fSFthreeSigmaFactor = threeSigmaFactor;
}



 private:
  bool fselectSector = false;
  bool fDoFiducialCut = false;
  bool fDoDCFiducial = false;
  bool fDoECALFiducial = false;
  bool fdoSFCut = false; 

  int fSFpid = 11;  // Default PID for SF cut
  float fSFmin = 0.0;  // Default minimum SF value
  float fSFminP = 100.0;  // Default minimum SF value for positive particles
  bool fPIDSystematics = false;
  float fSFthreeSigmaFactor = 1.0f;

  struct FiducialAxisCut {
    std::vector<std::pair<float, float>> excludedRanges;
  };

  struct FiducialRingCut {
    std::vector<std::tuple<float, float, float, float>> excludedRanges;
  };

  struct FiducialCut3D {
    FiducialAxisCut luCut;
    FiducialAxisCut lvCut;
    FiducialAxisCut lwCut;
  };

  struct FiducialCut2D_CVT {
    FiducialAxisCut thetaCut;
    FiducialAxisCut phiCut;
  };

  struct FiducialCutRing_FTCal {
    FiducialRingCut ringCut;
  };

  float fSector = -1;
  int fselectPID = -1;
  int fselectdetector = 1000;
  std::set<int> fSectors;  // allowed sectors for DC
  float fMinX = -999999, fMaxX = 999999;
  float fMinY = -999999, fMaxY = 999999;
  float fMinZ = -999999, fMaxZ = 999999;
  float fMinCX = -999999, fMaxCX = 999999;
  float fMinCY = -999999, fMaxCY = 999999;
  float fMinCZ = -999999, fMaxCZ = 999999;
  float fMinPath = -999999, fMaxPath = 999999;
  float fDCMinEdge = -999999, fDCMaxEdge = 999999;
  float fECALMinEdge = -999999, fECALMaxEdge = 999999;
  std::vector<std::pair<float, float>> fThetaBins;                   // Still used for reference or binning
  std::map<int, std::vector<float>> fDCEdgeCutsPerPID;
  std::map<int, std::vector<float>> fCVTEdgeCutsPerPID;

  /// ECin ECout and PCal Fiducial cuts
  std::map<int, std::map<int, FiducialCut2D_CVT>> fFiducialCutsCVT;
  std::map<int, std::map<int, FiducialCut2D_CVT>> fFiducialCutsCVT_Bhawani;
  std::map<int, std::map<int, FiducialCutRing_FTCal>> fFiducialCutsFTCal;
  std::map<int, std::map<int, FiducialCut3D>> fFiducialCutsPCal;
  std::map<int, std::map<int, FiducialCut3D>> fFiducialCutsECin;
  std::map<int, std::map<int, FiducialCut3D>> fFiducialCutsECout;

  std::map<int, std::map<int, SFCutABC>> fSFCutsMinCut;  // Sampling Fraction cuts per PID and sector
  std::map<int, std::map<int, SFCutABC>> fSFCutsMaxCut;  // Sampling Fraction cuts per PID and sector

  ///ECAL min energy cuts
  std::map<int, std::map<int, float>> fMinECALEnergyCutPerPIDLayer;

  template <typename T>
  bool IsInRange(T value, T min, T max) const {
    return value >= min && value <= max;
  }
};

#endif  // TRACKCUT_H_
