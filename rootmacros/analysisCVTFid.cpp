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
                               const std::string &filename, const std::string &treename) {
    TStopwatch timer;
    timer.Start();
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
                    const RVec<int16_t> &pindex, const RVec<int> &pid, const RVec<int> &passFid) {
        for (size_t i = 0; i < layer.size(); ++i) {
            int idx = pindex[i];
            if (det[i] != selecteddetector || pid[idx] != selectedPid || !passFid[idx]) continue;

            float chi2ndf = ndf[idx] > 0 ? chi2[idx] / ndf[idx] : 0;
            int lay = layer[i];
            float edg = edge[i];
            float th = theta[i];

            auto it = std::find(layers.begin(), layers.end(), lay);
            if (it == layers.end()) continue;
            size_t li = std::distance(layers.begin(), it);
            int ti = thetaRegionIndex(th, thetaCuts);

            std::string key = "L" + std::to_string(lay) + "_T" + std::to_string(ti);
            histos[key]->Fill(edg, chi2ndf);
        }
    }, {"REC_Traj_detector", "REC_Traj_layer", "REC_Traj_edge", "REC_Particle_theta", "REC_Track_chi2",
        "REC_Track_NDF", "REC_Traj_pindex", "REC_Particle_pid", "REC_Track_pass_fid"});

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
            pfx->SetTitle(Form("proton CVT layer %d;edge [cm];#LT#chi^{2}/ndf#GT", layers[li]));
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

        std::string outname = "proton_CVT_layer" + std::to_string(layers[li]) + ".png";
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
                        const std::string &filename, const std::string &treename) {
    TStopwatch timer;
    timer.Start();

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
        std::string name = Form("theta_vs_phi_L%d", layer);
        std::string title = Form("proton CVT layer %d;#phi [deg];#theta [deg]", layer);
        histos[layer] = new TH2F(name.c_str(), title.c_str(), 360, -180, 180, 360, 0, 180);
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
            if (pid[idx] != selectedPid || !passFid[idx]) continue;
            if (histos.count(lay))
                histos[lay]->Fill(phi_deg[i], theta_deg[i]);
        }
    }, {"REC_Traj_detector", "REC_Traj_layer", "theta_deg", "phi_deg",  "REC_Traj_pindex", "REC_Particle_pid", "REC_Track_pass_fid"});

    for (int layer : layers) {
        TCanvas *c = new TCanvas(Form("c_hit_layer_%d", layer), "", 2400, 2000);
        histos[layer]->Draw("COLZ");
        //gPad->SetLogz();

        c->SaveAs(Form("proton_CVT_hit_layer%d.png", layer));
        delete c;
    }

    timer.Stop();
    std::cout << "Time for DrawCVTHitResponse: ";
    timer.Print();
}


void analysisCVTFid() {
    //std::string path = "/work/clas12/yijie/clas12ana/analysis203/DISANA/build/bbbs/";
    std::string path = "./../build/";
    std::vector<int> layers = {1, 3, 5, 7, 12};
    std::vector<float> xmins = {-0.5, -0.5, -0.5, -4.0, -5.0};
    std::vector<float> xmaxs = {2.5, 2.5, 2.5, 20.0, 25.0};
    std::vector<float> thetaCuts = {55, 85, 115};
    DrawCVTChi2ndf_Optimized(2212, 5, 0, 400, 50, layers, xmins, xmaxs, thetaCuts, path + "dfSelected_before_fiducialCuts.root", "dfSelected_before");
    DrawCVTHitResponse(2212, 5, layers, path + "dfSelected_before_fiducialCuts.root", "dfSelected_before");
    gApplication->Terminate(0);
}
