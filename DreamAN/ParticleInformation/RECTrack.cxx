#include "RECTrack.h"
#include <iostream>
#include <cmath>
#include <functional>

const std::vector<std::string>& RECTrack::All() {
    static const std::vector<std::string> names = {
        "REC_Track_pindex",
        "REC_Track_detector",
        "REC_Track_sector",
        "REC_Track_chi2",
        "REC_Track_NDF",
    };
    return names;
}

const std::vector<std::string>& RECTrack::ForFiducialCut() {
    static const std::vector<std::string> names = {
        "REC_Track_pindex",
        "REC_Track_detector",
        "REC_Track_sector",
    };
    return names;
}

std::function<std::vector<float>(const std::vector<int16_t>& pindex,
                                 const std::vector<int16_t>& detector,
                                 const std::vector<int16_t>& sector,
                                 const std::vector<float>& chi2,
                                 const std::vector<int16_t>& NDF,
                                 const int& REC_Particle_num)> RECTrackchi2perndf(int target_detector) {
    return [target_detector](
                  const std::vector<int16_t>& pindex,
                  const std::vector<int16_t>& detector,
                  const std::vector<int16_t>& sector,
                  const std::vector<float>& chi2,
                  const std::vector<int16_t>& NDF,
                  const int& REC_Particle_num) -> std::vector<float> {
        // 初始化 pass_values，大小为 REC_Particle_num，默认所有粒子通过
        std::vector<float> return_values(REC_Particle_num, 9999);
        for (size_t i = 0; i < pindex.size(); ++i) {
            if (detector[i] == target_detector){
                return_values[pindex[i]] = chi2[i]/NDF[i];
            }
        }

        return return_values;
    };
}