#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TStopwatch.h>
#include <iostream>
#include <map>
#include <vector>
#include <tuple>

using namespace ROOT::VecOps;

int thetaRegionIndex(float thetaRad, const std::vector<float> &thetaCuts) {
    float deg = thetaRad * 180.0 / TMath::Pi();
    for (size_t i = 0; i < thetaCuts.size(); ++i) {
        if (deg <= thetaCuts[i]) return i;
    }
    return thetaCuts.size();
}

void DrawCVTChi2ndf_Optimized(const int &selectedPid, const int &selecteddetector,
                               const double ymin, const double ymax, const int nbinsx,
                               const std::vector<int> &layers,
                               const std::vector<float> &xmins,
                               const std::vector<float> &xmaxs,
                               const std::vector<float> &thetaCuts,
                               const std::string &filename, const std::string &treename, const bool doFid) {
    TStopwatch timer;
    timer.Start();
    ROOT::EnableImplicitMT();  // Enable multi-threading for RDataFrame
    ROOT::RDataFrame df(treename, filename);
    gStyle->SetOptStat(0);

    std::map<std::string, TH2F*> histos;

    for (size_t li = 0; li < layers.size(); ++li) {
        for (size_t ti = 0; ti <= thetaCuts.size(); ++ti) {
            std::string key = "L" + std::to_string(layers[li]) + "_T" + std::to_string(ti);
            histos[key] = new TH2F(key.c_str(), ";edge [cm];<chi2/ndf>", nbinsx, xmins[li], xmaxs[li], 500, ymin, ymax);
            histos[key]->SetDirectory(nullptr);
        }
    }

    df.Foreach([&](const RVec<short> &det, const RVec<short> &layer,
                const RVec<float> &edge, const RVec<float> &theta,
                const RVec<float> &chi2, const RVec<int16_t> &ndf,
                const RVec<int16_t> &trackpindex,
                const RVec<int16_t> &pindex, const RVec<int> &pid,
                const RVec<int> &passFid) {

    for (size_t i = 0; i < layer.size(); ++i) {

        if (i >= det.size() || i >= edge.size() ||
            i >= theta.size() || i >= pindex.size()) {
            continue;
        }

        int idx = pindex[i];

        if (idx < 0 || idx >= (int)pid.size()) {
            continue;
        }

        if (doFid && idx >= (int)passFid.size()) {
            continue;
        }

        int fidpass = doFid ? passFid[idx] : 1;

        if (det[i] != selecteddetector || pid[idx] != selectedPid || !fidpass) {
            continue;
        }

        bool foundTrack = false;
        float trackchi2 = 0.0;
        float trackndf  = 0.0;

        for (size_t j = 0; j < trackpindex.size(); ++j) {
            if (j >= chi2.size() || j >= ndf.size()) {
                continue;
            }

            if (trackpindex[j] == idx) {
                trackchi2 = chi2[j];
                trackndf  = ndf[j];
                foundTrack = true;
                break;
            }
        }

        if (!foundTrack) {
            continue;
        }

        if (trackndf <= 0) {
            continue;
        }

        float chi2ndf = trackchi2 / trackndf;

        if (!std::isfinite(chi2ndf)) {
            continue;
        }

        if (chi2ndf >= ymax / 2.0) {
            continue;
        }

        int lay = layer[i];

        auto it = std::find(layers.begin(), layers.end(), lay);
        if (it == layers.end()) {
            continue;
        }

        int ti = thetaRegionIndex(theta[i], thetaCuts);

        std::string key = "L" + std::to_string(lay)
                        + "_T" + std::to_string(ti);

        auto hIt = histos.find(key);
        if (hIt == histos.end() || hIt->second == nullptr) {
            continue;
        }

        hIt->second->Fill(edge[i], chi2ndf);
    }

}, {"REC_Traj_detector", "REC_Traj_layer", "REC_Traj_edge",
    "REC_Particle_theta", "REC_Track_chi2", "REC_Track_NDF",
    "REC_Track_pindex", "REC_Traj_pindex", "REC_Particle_pid",
    "REC_Track_pass_fid"});

    for (size_t li = 0; li < layers.size(); ++li) {
        TCanvas *c = new TCanvas(("c_layer" + std::to_string(layers[li])).c_str(), "", 2000, 1500);
        TLegend *leg = new TLegend(0.65, 0.7, 0.88, 0.88);
        leg->SetTextFont(42);
        leg->SetTextSize(0.03);
        int color = 1;

        for (int ti = 0; ti <= (int)thetaCuts.size(); ++ti) {
            std::string key = "L" + std::to_string(layers[li]) + "_T" + std::to_string(ti);
            if (!histos.count(key)) continue;

            TProfile *pfx = histos[key]->ProfileX(("prof_" + key).c_str(), 1, -1, "s");
            pfx->SetBins(nbinsx, xmins[li], xmaxs[li]);
            pfx->GetYaxis()->SetRangeUser(ymin, ymax/2);
            pfx->SetLineColor(color);
            pfx->SetMarkerColor(color);
            pfx->SetMarkerStyle(20);
            pfx->SetMarkerSize(1.0);
            pfx->SetLineWidth(2);
            if (selectedPid == 2212){
                if (doFid){
                    pfx->SetTitle(Form("proton CVT layer %d;edge [cm];#LT#chi^{2}/ndf#GT", layers[li]));
                }
                else {
                    pfx->SetTitle(Form("proton CVT layer %d (no fiducial cuts);edge [cm];#LT#chi^{2}/ndf#GT", layers[li]));
                }
            }
            else {
                if (doFid)
                    pfx->SetTitle(Form("pid%d CVT layer %d;edge [cm];#LT#chi^{2}/ndf#GT", selectedPid, layers[li]));
                else
                    pfx->SetTitle(Form("pid%d CVT layer %d (no fiducial cuts);edge [cm];#LT#chi^{2}/ndf#GT", selectedPid, layers[li]));
            }

            //pfx->SetTitle(Form("proton CVT layer %d", layers[li]));
            pfx->Draw(ti == 0 ? "PE" : "PE SAME");
            color++;

            std::string label;
            if (ti == 0)
                label = Form("#scale[0.9]{#theta/deg #leq %.0f}", thetaCuts[0]);
            else if (ti == (int)thetaCuts.size())
                label = Form("#scale[0.9]{#theta/deg > %.0f}", thetaCuts.back());
            else
                label = Form("#scale[0.9]{%.0f < #theta/deg #leq %.0f}", thetaCuts[ti - 1], thetaCuts[ti]);

            leg->AddEntry(pfx, label.c_str(), "lp");
        }
        leg->Draw();

        std::string outdir = "CVTChi2ndfPlots";
        gSystem->Exec(("mkdir -p " + outdir).c_str());

        std::string outname;
        if (selectedPid == 2212) {
            if (doFid)
                outname = outdir + "/proton_CVT_layer" + std::to_string(layers[li]) + ".png";
            else
                outname = outdir + "/proton_CVT_layer" + std::to_string(layers[li]) + "_noFid.png";
        } else {
            if (doFid)
                outname = outdir + "/pid" + std::to_string(selectedPid) + "_CVT_layer" + std::to_string(layers[li]) + ".png";
            else
                outname = outdir + "/pid" + std::to_string(selectedPid) + "_CVT_layer" + std::to_string(layers[li]) + "_noFid.png";
        }

        c->SaveAs(outname.c_str());
        std::cout << "Saved: " << outname << std::endl;
        delete c;

    }
    timer.Stop();
    std::cout << "Time for DrawCVTChi2ndf_Optimized: ";
    timer.Print();
}

void DrawCVTHitResponse(const int &selectedPid, const int &selecteddetector,
                        const std::vector<int> &layers,
                        const std::string &filename, const std::string &treename,
                        const bool doFid) {
    TStopwatch timer;
    timer.Start();
    ROOT::EnableImplicitMT();  // Enable multi-threading for RDataFrame
    ROOT::RDataFrame df(treename, filename);
    gStyle->SetOptStat(0);

    auto dfWithAngles = df.Define("theta_deg", [](const RVec<float> &x, const RVec<float> &y, const RVec<float> &z) {
                RVec<float> theta;
                for (size_t i = 0; i < x.size(); ++i)
                    theta.push_back(180.0 / TMath::Pi() * TMath::ACos(z[i] / sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i])));
                return theta;
            }, {"REC_Traj_x", "REC_Traj_y", "REC_Traj_z"})
            .Define("phi_deg", [](const RVec<float> &x, const RVec<float> &y) {
                RVec<float> phi;
                for (size_t i = 0; i < x.size(); ++i)
                    phi.push_back(180.0 / TMath::Pi() * TMath::ATan2(y[i], x[i]));
                return phi;
            }, {"REC_Traj_x", "REC_Traj_y"});

    // 初始化所有层的TH2F
    std::map<int, TH2F*> histos;
    for (int layer : layers) {
        std::string name;
        std::string title;
        if (selectedPid == 2212) {
            if (doFid) {
                name = Form("proton_theta_vs_phi_L%d", layer);
                title = Form("proton CVT layer %d;#phi [deg];#theta [deg]", layer);
            } else {
                name = Form("proton_theta_vs_phi_L%d_noFid", layer);
                title = Form("proton CVT layer %d (no fiducial cuts);#phi [deg];#theta [deg]", layer);
            }
        } else {
            if (doFid) {
                name = Form("pid%d_theta_vs_phi_L%d", selectedPid, layer);
                title = Form("pid%d CVT layer %d;#phi [deg];#theta [deg]", selectedPid, layer);
            } else {
                name = Form("pid%d_theta_vs_phi_L%d_noFid", selectedPid, layer);
                title = Form("pid%d CVT layer %d (no fiducial cuts);#phi [deg];#theta [deg]", selectedPid, layer);
            }
        }
        histos[layer] = new TH2F(name.c_str(), title.c_str(), 360, -180, 180, 360, 20, 140);
        histos[layer]->SetDirectory(nullptr);
    }

    dfWithAngles.Foreach([&](const RVec<short> &det, const RVec<short> &layer,
                             const RVec<float> &theta_deg, const RVec<float> &phi_deg,
                             const RVec<int16_t> &pindex, const RVec<int> &pid,
                             const RVec<int> &passFid) {
        for (size_t i = 0; i < det.size(); ++i) {
            if (det[i] != selecteddetector) continue;
            int lay = layer[i];
            int idx = pindex[i];
            int fidpass = doFid ? passFid[idx] : 1;  // 如果不需要fiducial cuts，则始终通过
            if (pid[idx] != selectedPid || !fidpass) continue;
            if (histos.count(lay))
                histos[lay]->Fill(phi_deg[i], theta_deg[i]);
        }
    }, {"REC_Traj_detector", "REC_Traj_layer", "theta_deg", "phi_deg",  "REC_Traj_pindex", "REC_Particle_pid", "REC_Track_pass_fid"});

    std::string outdir = "CVTHitPlots";
    gSystem->Exec(("mkdir -p " + outdir).c_str());

    for (int layer : layers) {
        TCanvas *c = new TCanvas(Form("c_hit_layer_%d", layer), "", 2400, 2000);
        histos[layer]->Draw("COLZ");

        std::string outname;
        if (selectedPid == 2212) {
            outname = outdir + "/" + Form(doFid ? "proton_CVT_hit_layer%d.png" : "proton_CVT_hit_layer%d_noFid.png", layer);
        } else {
            outname = outdir + "/" + Form(doFid ? "pid%d_CVT_hit_layer%d.png" : "pid%d_CVT_hit_layer%d_noFid.png", selectedPid, layer);
        }

        c->SaveAs(outname.c_str());
        std::cout << "Saved: " << outname << std::endl;
        delete c;
    }


    timer.Stop();
    std::cout << "Time for DrawCVTHitResponse: ";
    timer.Print();
}


void analysisCVTFid() {
    ROOT::EnableImplicitMT();

    std::string path = "../build/";

    std::vector<int> layers = {1, 3, 5, 7, 12};
    std::vector<float> xmins = {-0.5, -0.5, -0.5, -4.0, -4.0};
    std::vector<float> xmaxs = {2.5, 2.5, 2.5, 20.0, 20.0};
    std::vector<float> thetaCuts = {55, 85, 115};
    //DrawCVTChi2ndf_Optimized(2212, 5, 0, 60, 50, layers, xmins, xmaxs, thetaCuts, path + "dfSelected_afterFid_afterCorr.root", "dfSelected_afterFid_afterCorr", true);
    DrawCVTHitResponse(2212, 5, layers, path + "dfSelected_afterFid_afterCorr.root", "dfSelected_afterFid_afterCorr", true);
    //DrawCVTChi2ndf_Optimized(2212, 5, 0, 60, 50, layers, xmins, xmaxs, thetaCuts, path + "dfSelected.root", "dfSelected", false);
    DrawCVTHitResponse(2212, 5, layers, path + "dfSelected.root", "dfSelected", false);
    gApplication->Terminate(0);
}
