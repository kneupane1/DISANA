#include "TrackCut.h"

#include <cmath>
#include <iostream>
#include <map>

TrackCut::TrackCut() = default;
TrackCut::~TrackCut() = default;
TrackCut::TrackCut(const TrackCut& other) {
  this->fDCEdgeCutsPerPID = other.fDCEdgeCutsPerPID;
  this->fCVTEdgeCutsPerPID = other.fCVTEdgeCutsPerPID;
  this->fThetaBins = other.fThetaBins;
  this->fselectPID = other.fselectPID;
  this->fselectdetector = other.fselectdetector;
  this->fselectSector = other.fselectSector;
  this->fSector = other.fSector;
  this->fSectors = other.fSectors;
  this->fDoFiducialCut = other.fDoFiducialCut;
  this->fDoDCFiducial = other.fDoDCFiducial;
  this->fDoECALFiducial = other.fDoECALFiducial;

  this->fdoSFCut = other.fdoSFCut;
  this->fSFpid = other.fSFpid;
  this->fSFmin = other.fSFmin;
  this->fSFminP = other.fSFminP;
  this->fPIDSystematics = other.fPIDSystematics;
  this->fSFthreeSigmaFactor = other.fSFthreeSigmaFactor;

  this->fSFCutsMinCut = other.fSFCutsMinCut;
  this->fSFCutsMaxCut = other.fSFCutsMaxCut;

  // Also copy other necessary cuts if needed
  this->fMinX = other.fMinX;
  this->fMaxX = other.fMaxX;
  this->fMinY = other.fMinY;
  this->fMaxY = other.fMaxY;
  this->fMinZ = other.fMinZ;
  this->fMaxZ = other.fMaxZ;

  this->fMinCX = other.fMinCX;
  this->fMaxCX = other.fMaxCX;
  this->fMinCY = other.fMinCY;
  this->fMaxCY = other.fMaxCY;
  this->fMinCZ = other.fMinCZ;
  this->fMaxCZ = other.fMaxCZ;

  this->fMinPath = other.fMinPath;
  this->fMaxPath = other.fMaxPath;

  this->fDCMinEdge = other.fDCMinEdge;
  this->fDCMaxEdge = other.fDCMaxEdge;
  this->fECALMinEdge = other.fECALMinEdge;
  this->fECALMaxEdge = other.fECALMaxEdge;

  this->fFiducialCutsPCal = other.fFiducialCutsPCal;
  this->fFiducialCutsECin = other.fFiducialCutsECin;
  this->fFiducialCutsECout = other.fFiducialCutsECout;

  this->fFiducialCutsCVT = other.fFiducialCutsCVT;
  this->fFiducialCutsFTCal = other.fFiducialCutsFTCal;
}

void TrackCut::SetSectorCut(int SSector, int selectpid, int selectdetector, bool selectSector) {
  fSector = SSector;
  fselectSector = selectSector;
  fselectPID = selectpid;
  fselectdetector = selectdetector;
}

void TrackCut::SetSectorCut_Bhawani(const std::vector<int>& sectors, int selectpid, int selectdetector, bool selectSector) {
  fSectors = std::set<int>(sectors.begin(), sectors.end());
  fselectSector = selectSector;
  fselectPID = selectpid;
  fselectdetector = selectdetector;
}

void TrackCut::SetPositionCut(float minX, float maxX, float minY, float maxY, float minZ, float maxZ) {
  fMinX = minX;
  fMaxX = maxX;
  fMinY = minY;
  fMaxY = maxY;
  fMinZ = minZ;
  fMaxZ = maxZ;
}

void TrackCut::SetDirectionCut(float minCX, float maxCX, float minCY, float maxCY, float minCZ, float maxCZ) {
  fMinCX = minCX;
  fMaxCX = maxCX;
  fMinCY = minCY;
  fMaxCY = maxCY;
  fMinCZ = minCZ;
  fMaxCZ = maxCZ;
}

void TrackCut::SetPathLengthCut(float minPath, float maxPath) {
  fMinPath = minPath;
  fMaxPath = maxPath;
}

void TrackCut::SetDCEdgeCut(float minEdge, float maxEdge) {
  fDCMinEdge = minEdge;
  fDCMaxEdge = maxEdge;
}

void TrackCut::SetECALEdgeCut(float minEdge, float maxEdge) {
  fECALMinEdge = minEdge;
  fECALMaxEdge = maxEdge;
}

void TrackCut::SetThetaBins(const std::vector<std::pair<float, float>>& thetaBins) { fThetaBins = thetaBins; }

void TrackCut::SetDCEdgeCuts(int pid, const std::vector<float>& edgeCutsPerRegion) {
  if (edgeCutsPerRegion.size() != 3) {
    throw std::runtime_error("DC Edge cuts must have 3 values (for regions 1, 2, 3)");
  }
  fDCEdgeCutsPerPID[pid] = edgeCutsPerRegion;
  std::cout << "[Info] DC edge cuts for PID " << pid << ": ";
  for (auto e : edgeCutsPerRegion) std::cout << e << " ";
  std::cout << std::endl;
}

float TrackCut::GetEdgeCut(int pid, int region) const {
  auto it = fDCEdgeCutsPerPID.find(pid);
  if (it == fDCEdgeCutsPerPID.end()) {
    throw std::runtime_error("DC Edge cuts not defined for PID: " + std::to_string(pid));
  }
  return it->second[region - 1];
}

void TrackCut::SetCVTEdgeCuts(int pid, const std::vector<float>& edgeCutsPerLayer) {
  if (edgeCutsPerLayer.size() != 5) {
    throw std::runtime_error("CVT Edge cuts must have 5 values (for layers 1, 3, 5, 7, 12)");
  }
  fCVTEdgeCutsPerPID[pid] = edgeCutsPerLayer;
  std::cout << "[Info] CVT edge cuts for PID " << pid << ": ";
  for (auto e : edgeCutsPerLayer) std::cout << e << " ";
  std::cout << std::endl;
}

float TrackCut::GetCVTEdgeCut(int pid, int layer) const {
  auto it = fCVTEdgeCutsPerPID.find(pid);
  if (it == fCVTEdgeCutsPerPID.end()) {
    throw std::runtime_error("CVT Edge cuts not defined for PID: " + std::to_string(pid));
  }
  return it->second[layer - 1];
}

void TrackCut::AddCVTFiducialRange(int pid, int layer, const std::string& axis, float min, float max) {
  if (axis == "theta")
    fFiducialCutsCVT[pid][layer].thetaCut.excludedRanges.emplace_back(min, max);
  else if (axis == "phi")
    fFiducialCutsCVT[pid][layer].phiCut.excludedRanges.emplace_back(min, max);
}

void TrackCut::AddFTCalFiducialRange(int pid, int layer, float x, float y, float rmin, float rmax) {
  fFiducialCutsFTCal[pid][layer].ringCut.excludedRanges.emplace_back(x, y, rmin, rmax);
}

void TrackCut::AddPCalFiducialRange(int pid, int sector, const std::string& axis, float min, float max) {
  if (axis == "lu")
    fFiducialCutsPCal[pid][sector].luCut.excludedRanges.emplace_back(min, max);
  else if (axis == "lv")
    fFiducialCutsPCal[pid][sector].lvCut.excludedRanges.emplace_back(min, max);
  else if (axis == "lw")
    fFiducialCutsPCal[pid][sector].lwCut.excludedRanges.emplace_back(min, max);
}

void TrackCut::AddECinFiducialRange(int pid, int sector, const std::string& axis, float min, float max) {
  if (axis == "lu")
    fFiducialCutsECin[pid][sector].luCut.excludedRanges.emplace_back(min, max);
  else if (axis == "lv")
    fFiducialCutsECin[pid][sector].lvCut.excludedRanges.emplace_back(min, max);
  else if (axis == "lw")
    fFiducialCutsECin[pid][sector].lwCut.excludedRanges.emplace_back(min, max);
}

void TrackCut::AddECoutFiducialRange(int pid, int sector, const std::string& axis, float min, float max) {
  if (axis == "lu")
    fFiducialCutsECout[pid][sector].luCut.excludedRanges.emplace_back(min, max);
  else if (axis == "lv")
    fFiducialCutsECout[pid][sector].lvCut.excludedRanges.emplace_back(min, max);
  else if (axis == "lw")
    fFiducialCutsECout[pid][sector].lwCut.excludedRanges.emplace_back(min, max);
}

void TrackCut::SetMinECALEnergyCut(int pid, int layer, float minEnergy) {
  fMinECALEnergyCutPerPIDLayer[pid][layer] = minEnergy;
}

void TrackCut::AddSamplingFractionMinCut(int pid, int sector, float A0, float Bm1, float Cm2) {
  fSFCutsMinCut[pid][sector] = {A0, Bm1, Cm2};
}
void TrackCut::AddSamplingFractionMaxCut(int pid, int sector, float A0, float Bm1, float Cm2) {
  fSFCutsMaxCut[pid][sector] = {A0, Bm1, Cm2};
}

const std::map<int, std::vector<float>>& TrackCut::GetEdgeCuts() const { return fDCEdgeCutsPerPID; }
const std::map<int, std::vector<float>>& TrackCut::GetCVTEdgeCuts() const { return fCVTEdgeCutsPerPID; }

bool TrackCut::operator()(const std::vector<int16_t>& pindex, const std::vector<int16_t>& index, const std::vector<int>& detector, const std::vector<int>& layer,
                          const std::vector<float>& x, const std::vector<float>& y, const std::vector<float>& z, const std::vector<float>& cx, const std::vector<float>& cy,
                          const std::vector<float>& cz, const std::vector<float>& path, const std::vector<float>& edge) const {
  return true;  //
}

/// copy the fiducical cuts here and should be used it from the EventFilte
std::function<std::vector<int>(const std::vector<int16_t>& pindex,
                               const std::vector<int16_t>& detector,
                               const std::vector<int16_t>& layer,
                               const std::vector<float>& x,
                               const std::vector<float>& y,
                               const std::vector<float>& z,
                               const std::vector<float>& edge,
                               const std::vector<int>& pid,
                               const int& REC_Particle_num)>
TrackCut::RECTrajPass() const {
  return [this](const std::vector<int16_t>& pindex, const std::vector<int16_t>& detector, const std::vector<int16_t>& layer,
                const std::vector<float>& x, const std::vector<float>& y, const std::vector<float>& z,
                const std::vector<float>& edge, const std::vector<int>& pid,
                const int& REC_Particle_num) -> std::vector<int> {
    std::vector<int> pass_values(REC_Particle_num, 1);
    auto isExcluded = [](float value, const FiducialAxisCut& cut) -> bool {
      for (const auto& range : cut.excludedRanges) {
        if (value >= range.first && value <= range.second) return true;
      }
      return false;
    };
    for (size_t i = 0; i < pindex.size(); ++i) {
      if (detector[i] == 6) {
          if (fDoFiducialCut) {
              int region = 0;
              int absLayer = std::abs(layer[i]);
              if (absLayer == 6)       region = 1;
              else if (absLayer == 18) region = 2;
              else if (absLayer == 36) region = 3;
              if (region == 0) continue;
              if (pindex[i] < 0 || pindex[i] >= static_cast<int>(pid.size())) continue;
              int cur_pid = pid[pindex[i]];
              auto pidCuts = fDCEdgeCutsPerPID.find(cur_pid);
              if (pidCuts == fDCEdgeCutsPerPID.end()) continue;

              float edgeCut = pidCuts->second[region - 1];
              if (edge[i] <= edgeCut) {
                  if (pindex[i] < static_cast<int>(pass_values.size()))
                      pass_values[pindex[i]] = 0;
                  continue;
              }
          }
      }


      if (detector[i] == 5) {  // CVT
        if (fDoFiducialCut) {
          int region = 0;
          int absLayer = std::abs(layer[i]);
          // AFTER (safe):
          if (absLayer == 1)  region = 1;
          else if (absLayer == 3)  region = 2;
          else if (absLayer == 5)  region = 3;
          else if (absLayer == 7)  region = 4;
          else if (absLayer == 12) region = 5;
          if (region == 0) continue;
          if (pindex[i] < 0 || pindex[i] >= static_cast<int>(pid.size())) continue;

          int cur_pid = pid[pindex[i]];
          auto pidCuts = fCVTEdgeCutsPerPID.find(cur_pid);
          if (pidCuts == fCVTEdgeCutsPerPID.end()) continue;
          float edgeCut = pidCuts->second[region - 1];  // now safe
          if (edge[i] <= edgeCut) {
              if (pindex[i] < static_cast<int>(pass_values.size()))  // ADD THIS
                  pass_values[pindex[i]] = 0;
              continue;
          }

          const std::map<int, std::map<int, FiducialCut2D_CVT>>* cutMap = nullptr;
          cutMap = &fFiducialCutsCVT;

          if (cutMap) {
            int cur_pid = pid[pindex[i]];
            auto pidMapIt = cutMap->find(cur_pid);
            if (pidMapIt != cutMap->end()) {
              const auto& layerMap = pidMapIt->second;
              auto it = layerMap.find(layer[i]);
              if (it != layerMap.end()) {
                const FiducialCut2D_CVT& cut = it->second;
                float CVTtheta = 180.0 / TMath::Pi() * TMath::ACos(z[i] / sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]));
                float CVTphi = 180.0 / TMath::Pi() * TMath::ATan2(y[i], x[i]);
                if (isExcluded(CVTtheta, cut.thetaCut) || isExcluded(CVTphi, cut.phiCut)) {
                  if (pindex[i] < static_cast<int>(pass_values.size()))  // ADD THIS
                      pass_values[pindex[i]] = 0;
                  continue;
              }
              }
            }
          }
        }
      }
    }

    return pass_values;
  };
}

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
TrackCut::RECCalorimeterPass() const {
  return [this](const std::vector<int16_t>& pindex, const std::vector<int16_t>& detector, const std::vector<int16_t>& sector,
                const std::vector<int16_t>& layer, const std::vector<float>& energy,
                const std::vector<float>& lu, const std::vector<float>& lv, const std::vector<float>& lw,
                const std::vector<int>& pid, const std::vector<float>& p,
                const int& REC_Particle_num) -> std::vector<int> {
    // Initialize return_values with size REC_Particle_num and default value 9999.0
    std::vector<int> return_values(REC_Particle_num, 1);
    std::vector<float> SF14(REC_Particle_num, 0.0);
    std::vector<float> SF(REC_Particle_num, 0.0);
    std::vector<int> REC_Particle_Sector(REC_Particle_num, -1);// Later save it in the RDF
    auto isExcluded = [](float value, const FiducialAxisCut& cut) -> bool {
      for (const auto& range : cut.excludedRanges) {
        if (value >= range.first && value <= range.second) return true;
      }
      return false;
    };

    for (size_t i = 0; i < pindex.size(); ++i) {
      if (detector[i] == 7) {
        if (fDoFiducialCut) {
          // Guard: REC::Calorimeter pindex is -1 (or otherwise invalid) for
          // calorimeter hits that have no matched REC::Particle row.  Without
          // this check the SF14/SF/REC_Particle_Sector writes below would
          // index before the start of the heap buffer, corrupting memory and
          // eventually triggering std::length_error in a later vector resize.
          if (pindex[i] < 0 || pindex[i] >= REC_Particle_num) continue;

          const std::map<int, std::map<int, FiducialCut3D>>* cutMap = nullptr;
          if (layer[i] == 1)
            cutMap = &fFiducialCutsPCal;
          else if (layer[i] == 4)
            cutMap = &fFiducialCutsECin;
          else if (layer[i] == 7)
            cutMap = &fFiducialCutsECout;
          
          if (layer[i] == 1 || layer[i] == 4) SF14[pindex[i]] = SF14[pindex[i]] + energy[i];
          SF[pindex[i]] = SF[pindex[i]] + energy[i];
          REC_Particle_Sector[pindex[i]] = sector[i];

          if (cutMap) {
            int cur_pid = pid[pindex[i]];
            auto pidMapIt = cutMap->find(cur_pid);
            if (pidMapIt != cutMap->end()) {
              const auto& sectorMap = pidMapIt->second;
              auto it = sectorMap.find(sector[i]);
              if (it != sectorMap.end()) {
                const FiducialCut3D& cut = it->second;
                if (isExcluded(lu[i], cut.luCut) || isExcluded(lv[i], cut.lvCut) || isExcluded(lw[i], cut.lwCut)) {
                  return_values[pindex[i]] = 0;
                  continue;
                }
              }
            }
          }
        }
      }
    }
    for (size_t i = 0; i < SF.size(); ++i) {
      if (!fdoSFCut) continue;
      
      int cur_pid = pid[i];
      if (p[i] <= 0) continue;
      SF[i] = SF[i] / p[i];
      SF14[i] = SF14[i] / p[i];

      // Sector-dependent minimum cut
      auto minIt = fSFCutsMinCut.find(cur_pid);
      if (minIt != fSFCutsMinCut.end()) {
        auto& sectorMap = minIt->second;
        if (sectorMap.count(REC_Particle_Sector[i])) {
          const auto& abc = sectorMap.at(REC_Particle_Sector[i]);
          float minCut = abc.A0 + abc.Bm1 * p[i] + abc.Cm2 * (p[i] * p[i]);
          if (fPIDSystematics) {
            auto maxPidIt = fSFCutsMaxCut.find(cur_pid);
            if (maxPidIt != fSFCutsMaxCut.end()) {
              auto maxSectorIt = maxPidIt->second.find(REC_Particle_Sector[i]);
              if (maxSectorIt != maxPidIt->second.end()) {
                const auto& maxABC = maxSectorIt->second;
                const float maxCut =
                    maxABC.A0 + maxABC.Bm1 * p[i] + maxABC.Cm2 * (p[i] * p[i]);
                const float center = 0.5f * (minCut + maxCut);
                const float halfWidth = 0.5f * (maxCut - minCut);
                minCut = center - halfWidth * fSFthreeSigmaFactor;
              }
            }
          }
          if (SF[i] < minCut) {
            return_values[i] = 0;
          }
        }
      }

      // Sector-dependent maximum cut
      auto maxIt = fSFCutsMaxCut.find(cur_pid);
      if (maxIt != fSFCutsMaxCut.end()) {
        auto& sectorMap = maxIt->second;
        if (sectorMap.count(REC_Particle_Sector[i])) {
          const auto& abc = sectorMap.at(REC_Particle_Sector[i]);
          float maxCut = abc.A0 + abc.Bm1 * p[i] + abc.Cm2 * (p[i] * p[i]);
          if (fPIDSystematics) {
            auto minPidIt = fSFCutsMinCut.find(cur_pid);
            if (minPidIt != fSFCutsMinCut.end()) {
              auto minSectorIt = minPidIt->second.find(REC_Particle_Sector[i]);
              if (minSectorIt != minPidIt->second.end()) {
                const auto& minABC = minSectorIt->second;
                const float minCut =
                    minABC.A0 + minABC.Bm1 * p[i] + minABC.Cm2 * (p[i] * p[i]);
                const float center = 0.5f * (minCut + maxCut);
                const float halfWidth = 0.5f * (maxCut - minCut);
                maxCut = center + halfWidth * fSFthreeSigmaFactor;
              }
            }
          }
          if (SF[i] > maxCut) {
            return_values[i] = 0;
          }
        }
      }

      if (cur_pid == fSFpid && SF14[i] < fSFmin && p[i] > fSFminP) {
        return_values[i] = 0;
      }
    }
    return return_values;
  };
}

std::function<std::vector<int>(const std::vector<short>&,    // pindex
                               const std::vector<int16_t>&,  // detector
                               const std::vector<int16_t>&,  // layer
                               const std::vector<float>&,    // x
                               const std::vector<float>&,    // y
                               const std::vector<int>&,      // pid
                               const int& REC_Particle_num)>
TrackCut::RECForwardTaggerPass() const {
  return [this](const std::vector<short>& pindex, const std::vector<int16_t>& detector,
                const std::vector<int16_t>& layer,
                const std::vector<float>& x, const std::vector<float>& y,
                const std::vector<int>& pid,
                const int& REC_Particle_num) -> std::vector<int> {
    // Initialize return_values with size REC_Particle_num and default value 9999.0

    auto isExcluded = [](float x, float y, const FiducialRingCut& cut) -> bool {
      for (const auto& range : cut.excludedRanges) {
        float x0 = std::get<0>(range);
        float y0 = std::get<1>(range);
        float rmin = std::get<2>(range);
        float rmax = std::get<3>(range);
        if ((x0-x)*(x0-x)+(y0-y)*(y0-y) >= rmin*rmin && (x0-x)*(x0-x)+(y0-y)*(y0-y) <= rmax*rmax) return true;
      }
      return false;
    };

    std::vector<int> return_values(REC_Particle_num, 1);

    for (size_t i = 0; i < pindex.size(); ++i) {
      if (detector[i] == 10) {
        if (fDoFiducialCut) {
          // Guard: REC::ForwardTagger pindex can be -1 for unmatched hits.
          if (pindex[i] < 0 || pindex[i] >= REC_Particle_num) continue;

          const std::map<int, std::map<int, FiducialCutRing_FTCal>>* cutMap = nullptr;
          cutMap = &fFiducialCutsFTCal;

          if (cutMap) {
            int cur_pid = pid[pindex[i]];
            auto pidMapIt = cutMap->find(cur_pid);
            if (pidMapIt != cutMap->end()) {
              const auto& layerMap = pidMapIt->second;
              auto it = layerMap.find(layer[i]);
              if (it != layerMap.end()) {
                const FiducialCutRing_FTCal& cut = it->second;
                if (isExcluded(x[i], y[i], cut.ringCut) ) {
                  return_values[pindex[i]] = 0;
                  continue;
                }
              }
            }
          }
        }
        
      }
    }
    return return_values;
  };
}

std::function<std::vector<int>(
    const std::vector<int16_t>& traj_pindex,
    const std::vector<int16_t>& traj_index,
    const std::vector<int16_t>& traj_detector,
    const std::vector<int16_t>& traj_layer,
    const std::vector<float>& x, const std::vector<float>& y,
    const std::vector<float>& z, const std::vector<float>& cx,
    const std::vector<float>& cy, const std::vector<float>& cz,
    const std::vector<float>& path, const std::vector<float>& traj_edge,
    const std::vector<int16_t>& calo_pindex, const std::vector<int16_t>& calo_index,
    const std::vector<int16_t>& calo_detector, const std::vector<int16_t>& calo_sector,
    const std::vector<int16_t>& calo_layer, const std::vector<float>& calo_energy,
    const std::vector<float>& calo_time, const std::vector<float>& calo_path,
    const std::vector<float>& calo_chi2, const std::vector<float>& calo_x,
    const std::vector<float>& calo_y, const std::vector<float>& calo_z,
    const std::vector<float>& calo_hx, const std::vector<float>& calo_hy,
    const std::vector<float>& calo_hz, const std::vector<float>& calo_lu,
    const std::vector<float>& calo_lv, const std::vector<float>& calo_lw,
    const std::vector<float>& calo_du, const std::vector<float>& calo_dv,
    const std::vector<float>& calo_dw, const std::vector<float>& calo_m2u,
    const std::vector<float>& calo_m2v, const std::vector<float>& calo_m2w,
    const std::vector<float>& calo_m3u, const std::vector<float>& calo_m3v,
    const std::vector<float>& calo_m3w, const std::vector<short>& calo_status,
    const std::vector<int>& pid, const int& REC_Particle_num)>
TrackCut::RECFiducialPass() const {
  return [this](
      const std::vector<int16_t>& traj_pindex,
      const std::vector<int16_t>& traj_index,
      const std::vector<int16_t>& traj_detector,
      const std::vector<int16_t>& traj_layer,
      const std::vector<float>& x, const std::vector<float>& y,
      const std::vector<float>& z, const std::vector<float>& cx,
      const std::vector<float>& cy, const std::vector<float>& cz,
      const std::vector<float>& path, const std::vector<float>& traj_edge,
      const std::vector<int16_t>& calo_pindex, const std::vector<int16_t>& calo_index,
      const std::vector<int16_t>& calo_detector, const std::vector<int16_t>& calo_sector,
      const std::vector<int16_t>& calo_layer, const std::vector<float>& calo_energy,
      const std::vector<float>& calo_time, const std::vector<float>& calo_path,
      const std::vector<float>& calo_chi2, const std::vector<float>& calo_x,
      const std::vector<float>& calo_y, const std::vector<float>& calo_z,
      const std::vector<float>& calo_hx, const std::vector<float>& calo_hy,
      const std::vector<float>& calo_hz, const std::vector<float>& calo_lu,
      const std::vector<float>& calo_lv, const std::vector<float>& calo_lw,
      const std::vector<float>& calo_du, const std::vector<float>& calo_dv,
      const std::vector<float>& calo_dw, const std::vector<float>& calo_m2u,
      const std::vector<float>& calo_m2v, const std::vector<float>& calo_m2w,
      const std::vector<float>& calo_m3u, const std::vector<float>& calo_m3v,
      const std::vector<float>& calo_m3w, const std::vector<short>& calo_status,
      const std::vector<int>& pid, const int& REC_Particle_num) -> std::vector<int> {

    std::vector<int> result(REC_Particle_num, 1);
    //std::cout << "traj_pindex.size() " << result.size() << "calo_pindex.size() " << calo_pindex.size() << std::endl;

    auto isExcluded = [](float val, const FiducialAxisCut& cut) -> bool {
      for (const auto& range : cut.excludedRanges) {
        if (val >= range.first && val <= range.second) return true;
      }
      return false;
    };

    // === DC (TRAJ) cuts ===
    for (size_t i = 0; i < traj_pindex.size(); ++i) {
      if (!fDoDCFiducial || traj_detector[i] != 6) continue;
      if (traj_pindex[i] < 0 || traj_pindex[i] >= static_cast<int>(pid.size())) continue;

      int region = 0;
      int absLayer = std::abs(traj_layer[i]);
      if (absLayer == 6) region = 1;
      else if (absLayer == 18) region = 2;
      else if (absLayer == 36) region = 3;
      else continue;

      int pid_i = pid[traj_pindex[i]];
      auto it = fDCEdgeCutsPerPID.find(pid_i);
      if (it != fDCEdgeCutsPerPID.end()) {
        float edgeCut = it->second[region - 1];
        if (traj_edge[i] <= edgeCut) {
          result[traj_pindex[i]] = 0;
        }
      }
    }
 
    // === ECAL (CALO) cuts ===
    for (size_t i = 0; i < calo_pindex.size(); ++i) {
      if (!fDoECALFiducial || calo_detector[i] != 7) continue;
      if (calo_pindex[i] < 0 || calo_pindex[i] >= static_cast<int>(pid.size())) continue;

      const std::map<int, std::map<int, FiducialCut3D>>* cutMap = nullptr;
      if (calo_layer[i] == 1)
        cutMap = &fFiducialCutsPCal;
      else if (calo_layer[i] == 4)
        cutMap = &fFiducialCutsECin;
      else if (calo_layer[i] == 7)
        cutMap = &fFiducialCutsECout;
      else
        continue;

      int pid_i = pid[calo_pindex[i]];
      auto pidMapIt = cutMap->find(pid_i);
      if (pidMapIt != cutMap->end()) {
        auto sectorIt = pidMapIt->second.find(calo_sector[i]);
        if (sectorIt != pidMapIt->second.end()) {
          const auto& cut = sectorIt->second;
          if (isExcluded(calo_lu[i], cut.luCut) ||
              isExcluded(calo_lv[i], cut.lvCut) ||
              isExcluded(calo_lw[i], cut.lwCut)) {
            result[calo_pindex[i]] = 0;
          }
        }
      }
    }

    return result;
  };
}
