#ifndef RECTRACK_H
#define RECTRACK_H

#include <vector>
#include <string>
#include <functional>


struct RECTrack {
    static const std::vector<std::string>& All();
    static const std::vector<std::string>& Extend();
    static const std::vector<std::string>& ForFiducialCut();

    using AllTypes = std::tuple<
        const std::vector<int>&,      // pindex
        const std::vector<int>&,      // detector
        const std::vector<int>&,      // sector
        const std::vector<float>&,    // chi2
        const std::vector<int>&    // NDF
    >;
};

std::function<std::vector<float>(const std::vector<int16_t>& pindex,
                                 const std::vector<int16_t>& detector,
                                 const std::vector<int16_t>& sector,
                                 const std::vector<float>& chi2,
                                 const std::vector<int16_t>& NDF,
                                 const int& REC_Particle_num)> RECTrackchi2perndf(int target_detector);

#endif // RECTRACK_H