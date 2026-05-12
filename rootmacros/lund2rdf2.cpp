// Run with: root -l 'lund2rdf.cpp("dvcsgen*.dat","out.root")'
#include <TFile.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <glob.h>   // POSIX glob for pattern matching
#include <cmath>

// Split a line by whitespace
static std::vector<std::string> split_ws(const std::string &line) {
    std::istringstream iss(line);
    std::vector<std::string> toks;
    std::string t;
    while (iss >> t) toks.push_back(t);
    return toks;
}

// Trim leading/trailing spaces
static std::string trim(const std::string &s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

// Expand a glob pattern into a list of file names
static std::vector<std::string> glob_files(const char* pattern) {
    glob_t results;
    std::vector<std::string> files;
    if (glob(pattern, 0, nullptr, &results) == 0) {
        for (size_t i = 0; i < results.gl_pathc; ++i) {
            files.emplace_back(results.gl_pathv[i]);
        }
    }
    globfree(&results);
    return files;
}

// Main function: read all files matching pattern, write a single ROOT file
void lund2rdf2(const char* inPattern = "/work/clas12/yijie/Simulation/DVCS/km15gen/km15gen/output_7p5_1/*.txt",
               const char* outPath   = "noradkm15_7546.root")
{
    // Prepare output ROOT file
    TFile fout(outPath, "RECREATE");
    if (fout.IsZombie()) {
        std::cerr << "ERROR: cannot create output file: " << outPath << std::endl;
        return;
    }

    // Branch containers
    std::vector<int>   v_pid;
    std::vector<float> v_px, v_py, v_pz;

    TTree tree("MC", "MC particles from LUND (dvcsgen)");
    tree.Branch("MC_Particle_pid", &v_pid);
    tree.Branch("MC_Particle_px",  &v_px);
    tree.Branch("MC_Particle_py",  &v_py);
    tree.Branch("MC_Particle_pz",  &v_pz);

    long long evCount = 0;

    // Expand glob pattern
    auto files = glob_files(inPattern);
    if (files.empty()) {
        std::cerr << "No files match pattern: " << inPattern << std::endl;
        return;
    }

    // Loop over files
    for (auto& fname : files) {
        std::ifstream fin(fname);
        if (!fin.is_open()) {
            std::cerr << "WARNING: cannot open file: " << fname << std::endl;
            continue;
        }
        std::cout << "Processing " << fname << " ..." << std::endl;

        std::string line;
        while (std::getline(fin, line)) {
            line = trim(line);
            if (line.empty()) continue;

            auto toks = split_ws(line);
            if (toks.empty()) continue;

            int nPart = -1;
            try {
                nPart = std::stoi(toks[0]); // Npart
            } catch (...) {
                continue; // not an event header
            }
            if (nPart <= 0) continue;

            v_pid.clear(); v_px.clear(); v_py.clear(); v_pz.clear();
            v_pid.reserve(nPart);
            v_px .reserve(nPart);
            v_py .reserve(nPart);
            v_pz .reserve(nPart);

            // === 新增：记录能量最大的 photon ===
            int   bestPhotonIdx = -1;
            float bestPhotonE   = -1.0f;

            for (int i = 0; i < nPart; ++i) {
                std::string pline;
                if (!std::getline(fin, pline)) break;
                pline = trim(pline);
                if (pline.empty()) { --i; continue; }

                auto ptoks = split_ws(pline);
                if (ptoks.size() < 9) continue;

                try {
                    int   pid = std::stoi(ptoks[3]);
                    float px  = std::stof(ptoks[6]);
                    float py  = std::stof(ptoks[7]);
                    float pz  = std::stof(ptoks[8]);

                    // 先把粒子统一 push 进容器
                    int idx = (int)v_pid.size();
                    v_pid.push_back(pid);
                    v_px .push_back(px);
                    v_py .push_back(py);
                    v_pz .push_back(pz);

                    // 如果是 photon，更新“最大能量 photon”的 index
                    if (pid == 22) {
                        float E = std::sqrt(px*px + py*py + pz*pz); // 对 photon: E ≈ |p|
                        if (E > bestPhotonE) {
                            bestPhotonE   = E;
                            bestPhotonIdx = idx;
                        }
                    }
                } catch (...) {
                    // skip bad line
                }
            }

            // === 第二步：只保留最大能量的 photon，其余 photon 删掉 ===
            if (bestPhotonIdx >= 0) {
                std::vector<int>   new_pid;
                std::vector<float> new_px, new_py, new_pz;

                new_pid.reserve(v_pid.size());
                new_px .reserve(v_px.size());
                new_py .reserve(v_py.size());
                new_pz .reserve(v_pz.size());

                for (size_t i = 0; i < v_pid.size(); ++i) {
                    int pid = v_pid[i];
                    // 非 photon 全留下；photon 只留下 bestPhotonIdx 那一个
                    if (pid != 22 || (int)i == bestPhotonIdx) {
                        new_pid.push_back(pid);
                        new_px .push_back(v_px[i]);
                        new_py .push_back(v_py[i]);
                        new_pz .push_back(v_pz[i]);
                    }
                }

                v_pid.swap(new_pid);
                v_px .swap(new_px);
                v_py .swap(new_py);
                v_pz .swap(new_pz);
            }

            tree.Fill();
            ++evCount;
        }
    }

    fout.Write();
    fout.Close();

    std::cout << "Done. Wrote " << evCount << " events to " << outPath << std::endl;
}
