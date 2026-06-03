#include "DISANA_Xplotter2csv.cpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <utility>
#include <vector>
#include <TF1.h>
#include <TGraphErrors.h>

// ================================================================
// Define reconstructed proton sector
// ================================================================
ROOT::RDF::RNode DefineRecoProtonSector(ROOT::RDF::RNode df) {
    return df.Define("recpro_sector",
        [](const ROOT::VecOps::RVec<int>& pid,
           const ROOT::VecOps::RVec<bool>& pass,
           const ROOT::VecOps::RVec<int16_t>& track_pindex,
           const ROOT::VecOps::RVec<short>& track_sector) {

            int proton_idx = -1;

            for (size_t i = 0; i < pid.size(); ++i) {
                if (pid[i] == 2212 && pass[i]) {
                    proton_idx = static_cast<int>(i);
                    break;
                }
            }

            if (proton_idx < 0) return -1;

            for (size_t j = 0; j < track_pindex.size(); ++j) {
                if (j >= track_sector.size()) continue;

                if (track_pindex[j] == proton_idx) {
                    int sec = track_sector[j];
                    if (sec >= 1 && sec <= 6) return sec;
                }
            }

            return -1;
        },
        {"REC_Particle_pid",
         "REC_Particle_pass_std",
         "REC_Track_pindex",
         "REC_Track_sector"});
}

double ProtonEfficiencyDataWeight(int pro_det_region,
                                  double recpro_theta,
                                  double recpro_p,
                                  int recpro_sector) {
    double pro_theta = recpro_theta * 180.0 / std::acos(-1.0);

    auto EvalPol4 = [](double x,
                       double p0,
                       double p1,
                       double p2,
                       double p3,
                       double p4) {
        return p0 + x * (p1 + x * (p2 + x * (p3 + x * p4)));
    };

    if (pro_det_region == 2) {
        // CD proton weights can also be split by recpro_p here.
        // Example:
        // if (recpro_p < 0.6) { ... }
        double weight2  = 194.113650031;
        weight2 = weight2 - 14.9235275433 * pro_theta;
        weight2 = weight2 + 0.431130425403 * pro_theta * pro_theta;
        weight2 = weight2 - 0.00551949772084 * pro_theta * pro_theta * pro_theta;
        weight2 = weight2 + 0.0000263906217339 * pro_theta * pro_theta * pro_theta * pro_theta;
        return 1.0/weight2;
    }

    if (pro_det_region == 1) {
        if (recpro_p >= 0.0 && recpro_p < 1.2) {
            double weight1 = 1.0;

            if (recpro_sector == 1) {
                weight1 = EvalPol4(pro_theta,
                                    -19.7716898552,
                                    3.01402168055,
                                    -0.163630873818,
                                    0.0038696341387,
                                    -3.35266065754e-05);
            } else if (recpro_sector == 2) {
                weight1 = EvalPol4(pro_theta,
                                    18.7984939041,
                                    -2.35340457263,
                                    0.110872026708,
                                    -0.00222906760896,
                                    1.60904132746e-05);
            } else if (recpro_sector == 3) {
                weight1 = EvalPol4(pro_theta,
                                    5.63566735157,
                                    -0.522093175844,
                                    0.0173336373725,
                                    -0.000143110588775,
                                    -1.10839774731e-06);
            } else if (recpro_sector == 4) {
                weight1 = EvalPol4(pro_theta,
                                    -30.1303875783,
                                    4.2743643373,
                                    -0.219263641405,
                                    0.00495202478552,
                                    -4.16003363644e-05);
            } else if (recpro_sector == 5) {
                weight1 = EvalPol4(pro_theta,
                                    24.284454148,
                                    -3.04287347804,
                                    0.14139967477,
                                    -0.00278355773718,
                                    1.93900421411e-05);
            } else if (recpro_sector == 6) {
                weight1 = EvalPol4(pro_theta,
                                    -7.37655860381,
                                    1.11750148989,
                                    -0.058354913403,
                                    0.0013760829902,
                                    -1.23555511421e-05);
            }

            return 1.0 / weight1;
        }

        return 1.0;
    }

    return 1.0;
}


// ================================================================
// Plot proton theta, p, phi comparison + Data/MC ratio
// theta, p: 10 bins
// phi: 30 bins
// ================================================================
void PlotProtonThetaPComparison(ROOT::RDF::RNode df1,
                                ROOT::RDF::RNode df2,
                                const std::string& label1 = "Pi0 MC",
                                const std::string& label2 = "Pi0 Data",
                                const std::string& outname = "ProtonThetaPComparison.png",
                                double theta_min = 0.0,
                                double theta_max = 80.0,
                                double p_min = 0.0,
                                double p_max = 4.0,
                                double luminosity1 = 1.0,
                                double luminosity2 = 1.0,
                                int theta_ratio_fit_order = 0) {

    if (luminosity1 <= 0 || luminosity2 <= 0) {
        std::cerr << "ERROR: luminosities must be > 0" << std::endl;
        return;
    }

    std::string tag = outname;
    for (auto& ch : tag) {
        if (!std::isalnum(static_cast<unsigned char>(ch))) ch = '_';
    }

    auto df1_deg = df1
        .Define("recpro_theta_deg", "recpro_theta * 180.0 / TMath::Pi()")
        .Define("recpro_phi_deg",   "recpro_phi   * 180.0 / TMath::Pi()");

    auto df2_weighted = df2
        .Define("data_pro_eff_weight",
                "ProtonEfficiencyDataWeight(pro_det_region, recpro_theta, recpro_p, recpro_sector)");

    auto df2_deg = df2_weighted
        .Define("recpro_theta_deg", "recpro_theta * 180.0 / TMath::Pi()")
        .Define("recpro_phi_deg",   "recpro_phi   * 180.0 / TMath::Pi()");

    auto h_theta_mc = df1_deg.Histo1D(
        {Form("h_theta_mc_%s", tag.c_str()),
         ";Proton #theta [deg];Counts / fb^{-1}",
         30, theta_min, theta_max},
        "recpro_theta_deg");

    auto h_theta_data = df2_deg.Histo1D(
        {Form("h_theta_data_%s", tag.c_str()),
         ";Proton #theta [deg];Counts / fb^{-1}",
         30, theta_min, theta_max},
        "recpro_theta_deg",
        "data_pro_eff_weight");

    auto h_p_mc = df1.Histo1D(
        {Form("h_p_mc_%s", tag.c_str()),
         ";Proton Momentum p [GeV/c];Counts / fb^{-1}",
         30, p_min, p_max},
        "recpro_p");

    auto h_p_data = df2_weighted.Histo1D(
        {Form("h_p_data_%s", tag.c_str()),
         ";Proton Momentum p [GeV/c];Counts / fb^{-1}",
         30, p_min, p_max},
        "recpro_p",
        "data_pro_eff_weight");

    auto h_phi_mc = df1_deg.Histo1D(
        {Form("h_phi_mc_%s", tag.c_str()),
         ";Proton #phi [deg];Counts / fb^{-1}",
         30, 0.0, 360.0},
        "recpro_phi_deg");

    auto h_phi_data = df2_deg.Histo1D(
        {Form("h_phi_data_%s", tag.c_str()),
         ";Proton #phi [deg];Counts / fb^{-1}",
         30, 0.0, 360.0},
        "recpro_phi_deg",
        "data_pro_eff_weight");

    h_theta_mc->Scale(1.0 / luminosity1);
    h_p_mc->Scale(1.0 / luminosity1);
    h_phi_mc->Scale(1.0 / luminosity1);

    h_theta_data->Scale(1.0 / luminosity2);
    h_p_data->Scale(1.0 / luminosity2);
    h_phi_data->Scale(1.0 / luminosity2);

    auto StyleMCData = [](TH1D* h_mc, TH1D* h_data) {
        h_mc->SetLineColor(kRed + 1);
        h_mc->SetLineWidth(3);

        h_data->SetLineColor(kBlue + 1);
        h_data->SetLineWidth(3);
    };

    StyleMCData(h_theta_mc.GetPtr(), h_theta_data.GetPtr());
    StyleMCData(h_p_mc.GetPtr(),     h_p_data.GetPtr());
    StyleMCData(h_phi_mc.GetPtr(),   h_phi_data.GetPtr());

    TCanvas* c = new TCanvas(
        Form("c_proton_compare_%s", tag.c_str()),
        "Proton Kinematics Comparison",
        1800, 600);

    c->Divide(3, 1);

    auto DrawCompare = [&](int pad, TH1D* h_mc, TH1D* h_data) {
        c->cd(pad);
        gPad->SetGrid();

        double ymax = std::max(h_mc->GetMaximum(), h_data->GetMaximum());
        h_mc->SetMaximum(ymax * 1.25);

        h_mc->Draw("hist");
        h_data->Draw("hist same");

        TLegend* leg = new TLegend(0.58, 0.74, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(h_mc, label1.c_str(), "l");
        leg->AddEntry(h_data, label2.c_str(), "l");
        leg->Draw();
    };

    DrawCompare(1, h_theta_mc.GetPtr(), h_theta_data.GetPtr());
    DrawCompare(2, h_p_mc.GetPtr(),     h_p_data.GetPtr());
    DrawCompare(3, h_phi_mc.GetPtr(),   h_phi_data.GetPtr());

    c->SaveAs(outname.c_str());

    TH1D* r_theta = (TH1D*)h_theta_data->Clone(Form("r_theta_%s", tag.c_str()));
    TH1D* r_p     = (TH1D*)h_p_data->Clone(Form("r_p_%s", tag.c_str()));
    TH1D* r_phi   = (TH1D*)h_phi_data->Clone(Form("r_phi_%s", tag.c_str()));

    r_theta->Divide(h_theta_mc.GetPtr());
    r_p->Divide(h_p_mc.GetPtr());
    r_phi->Divide(h_phi_mc.GetPtr());

    r_theta->SetTitle(";Proton #theta [deg];Data / MC");
    r_p->SetTitle(";Proton Momentum p [GeV/c];Data / MC");
    r_phi->SetTitle(";Proton #phi [deg];Data / MC");

    for (auto* r : {r_theta, r_p, r_phi}) {
        r->SetLineColor(kBlack);
        r->SetMarkerStyle(20);
        r->SetMarkerSize(0.8);
        r->SetMinimum(0.0);
        r->SetMaximum(2.0);
    }

    if (theta_ratio_fit_order < 0) theta_ratio_fit_order = 0;
    if (theta_ratio_fit_order > 9) theta_ratio_fit_order = 9;

    TF1* f_theta_ratio = new TF1(
        Form("f_theta_ratio_%s", tag.c_str()),
        Form("pol%d", theta_ratio_fit_order),
        theta_min,
        theta_max);
    f_theta_ratio->SetLineColor(kGreen + 2);
    f_theta_ratio->SetLineWidth(3);
    f_theta_ratio->SetLineStyle(1);
    auto fit_result = r_theta->Fit(f_theta_ratio, "Q0RS");

    TGraphErrors* g_theta_ratio_uncertainty = new TGraphErrors(200);
    auto fit_covariance = fit_result->GetCovarianceMatrix();
    for (int i = 0; i < g_theta_ratio_uncertainty->GetN(); ++i) {
        double x = theta_min;
        if (g_theta_ratio_uncertainty->GetN() > 1) {
            x += (theta_max - theta_min) * i / (g_theta_ratio_uncertainty->GetN() - 1);
        }

        double variance = 0.0;
        for (int ipar = 0; ipar <= theta_ratio_fit_order; ++ipar) {
            const double xi = std::pow(x, ipar);
            for (int jpar = 0; jpar <= theta_ratio_fit_order; ++jpar) {
                const double xj = std::pow(x, jpar);
                variance += xi * fit_covariance(ipar, jpar) * xj;
            }
        }

        g_theta_ratio_uncertainty->SetPoint(i, x, f_theta_ratio->Eval(x));
        g_theta_ratio_uncertainty->SetPointError(
            i,
            0.0,
            std::sqrt(std::max(0.0, variance)));
    }
    g_theta_ratio_uncertainty->SetFillColorAlpha(kGreen + 2, 0.25);
    g_theta_ratio_uncertainty->SetFillStyle(1001);
    g_theta_ratio_uncertainty->SetLineColor(kGreen + 2);
    g_theta_ratio_uncertainty->SetLineWidth(0);

    std::string ratio_outname = outname;
    size_t pos = ratio_outname.rfind(".png");
    if (pos != std::string::npos) {
        ratio_outname.replace(pos, 4, "_ratio.png");
    } else {
        ratio_outname += "_ratio.png";
    }

    std::string fit_txt_outname = outname;
    pos = fit_txt_outname.rfind(".png");
    if (pos != std::string::npos) {
        fit_txt_outname.replace(
            pos,
            4,
            Form("_theta_ratio_pol%d_fit.txt", theta_ratio_fit_order));
    } else {
        fit_txt_outname += Form("_theta_ratio_pol%d_fit.txt", theta_ratio_fit_order);
    }

    std::ofstream fit_txt(fit_txt_outname);
    fit_txt << std::setprecision(12);
    fit_txt << "# theta ratio Data/MC fit\n";
    fit_txt << "# function: pol" << theta_ratio_fit_order << "\n";
    fit_txt << "# range: [" << theta_min << ", " << theta_max << "] deg\n";
    fit_txt << "# status: " << static_cast<int>(fit_result) << "\n";
    fit_txt << "# chi2: " << f_theta_ratio->GetChisquare() << "\n";
    fit_txt << "# ndf: " << f_theta_ratio->GetNDF() << "\n";
    fit_txt << "# index parameter error\n";
    for (int ipar = 0; ipar <= theta_ratio_fit_order; ++ipar) {
        fit_txt << ipar << " "
                << f_theta_ratio->GetParameter(ipar) << " "
                << f_theta_ratio->GetParError(ipar) << "\n";
    }

    TCanvas* c_ratio = new TCanvas(
        Form("c_proton_ratio_%s", tag.c_str()),
        "Data / MC Proton Kinematics Ratio",
        1800, 600);

    c_ratio->Divide(3, 1);

    auto DrawRatio = [&](int pad,
                         TH1D* r,
                         TF1* fit = nullptr,
                         TGraphErrors* uncertainty = nullptr) {
        c_ratio->cd(pad);
        gPad->SetGrid();

        r->Draw("E1");

        double xmin = r->GetXaxis()->GetXmin();
        double xmax = r->GetXaxis()->GetXmax();

        TLine* line = new TLine(xmin, 1.0, xmax, 1.0);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->SetLineStyle(2);
        line->Draw("same");

        if (uncertainty) {
            uncertainty->Draw("3 same");
        }

        if (fit) {
            fit->Draw("same");
            r->Draw("E1 same");

            TLegend* leg = new TLegend(0.52, 0.76, 0.88, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(r, "Data / MC", "lep");
            leg->AddEntry(fit, Form("pol%d fit", theta_ratio_fit_order), "l");
            leg->Draw();
        }
    };

    DrawRatio(1, r_theta, f_theta_ratio, g_theta_ratio_uncertainty);
    DrawRatio(2, r_p);
    DrawRatio(3, r_phi);

    c_ratio->SaveAs(ratio_outname.c_str());

    std::cout << "Saved comparison plot to: " << outname << std::endl;
    std::cout << "Saved ratio plot to: " << ratio_outname << std::endl;
    std::cout << "Saved theta ratio fit parameters to: " << fit_txt_outname << std::endl;
}

std::string MomentumRangeTag(double p_min, double p_max) {
    auto format = [](double v) {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(2) << v;
        std::string s = ss.str();
        std::replace(s.begin(), s.end(), '.', 'p');
        return s;
    };
    return "p" + format(p_min) + "_" + format(p_max);
}

void PlotProtonThetaPComparisonByPRange(
        ROOT::RDF::RNode df1,
        ROOT::RDF::RNode df2,
        const std::vector<std::pair<double, double>>& p_ranges,
        const std::string& label1,
        const std::string& label2,
        const std::string& outname_base,
        double theta_min,
        double theta_max,
        double luminosity1,
        double luminosity2,
        int theta_ratio_fit_order = 0) {

    for (const auto& p_range : p_ranges) {
        const double p_min = p_range.first;
        const double p_max = p_range.second;

        if (p_max <= p_min) {
            std::cerr << "Skip invalid proton p range: ["
                      << p_min << ", " << p_max << "]" << std::endl;
            continue;
        }

        const std::string pcut = Form("recpro_p >= %.12g && recpro_p < %.12g", p_min, p_max);
        const std::string ptag = MomentumRangeTag(p_min, p_max);
        const std::string mc_filter_label = "MC proton " + ptag;
        const std::string data_filter_label = "Data proton " + ptag;

        std::string outname = outname_base;
        const size_t pos = outname.rfind(".png");
        if (pos != std::string::npos) {
            outname.replace(pos, 4, "_" + ptag + ".png");
        } else {
            outname += "_" + ptag + ".png";
        }

        PlotProtonThetaPComparison(
            df1.Filter(pcut, mc_filter_label),
            df2.Filter(pcut, data_filter_label),
            label1 + " " + ptag,
            label2 + " " + ptag,
            outname,
            theta_min,
            theta_max,
            p_min,
            p_max,
            luminosity1,
            luminosity2,
            theta_ratio_fit_order);
    }
}


// ================================================================
// Main
// ================================================================
void analysisProEff() {
    ROOT::EnableImplicitMT(40);
    float beam_energy = 7.546;

    std::string input_path_from_analysisRun_7546_pi0MC =
        "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File7546Pi0/mc/";

    std::string filename_afterFid_7546_pi0MC =
        Form("%s/dfSelected_afterFid_afterCorr.root",
             input_path_from_analysisRun_7546_pi0MC.c_str());

    ROOT::RDF::RNode df_afterFid_7546_pi0MC_init =
        InitKinematics(filename_afterFid_7546_pi0MC,
                       "dfSelected_afterFid_afterCorr",
                       beam_energy);

    ROOT::RDF::RNode df_afterFid_7546_pi0MC =
        GetSlim_exclusive(df_afterFid_7546_pi0MC_init,
                          "ProtonEff/dfSlim_7546_pi0MC.root",
                          "dfSlim_7546_pi0MC",
                          false);

    auto df_final_OnlPi0_7546_pi0MC =
        ApplyFinalGenDVPi0Selections(
            Init2PhotonKinematics(
                SelectPi0Event(df_afterFid_7546_pi0MC),
                beam_energy));


    std::string input_path_from_analysisRun_7546_pi0data =
        "/work/clas12/yijie/clas12ana/analysis1301/DISANA/build/File7546Pi0/data/";

    std::string filename_afterFid_7546_pi0data =
        Form("%s/dfSelected_afterFid_afterCorr.root",
             input_path_from_analysisRun_7546_pi0data.c_str());

    ROOT::RDF::RNode df_afterFid_7546_pi0data_init =
        InitKinematics(filename_afterFid_7546_pi0data,
                       "dfSelected_afterFid_afterCorr",
                       beam_energy);

    ROOT::RDF::RNode df_afterFid_7546_pi0data =
        GetSlim_exclusive(df_afterFid_7546_pi0data_init,
                          "ProtonEff/dfSlim_7546_pi0data.root",
                          "dfSlim_7546_pi0data",
                          false);

    auto df_final_OnlPi0_7546_pi0data =
        ApplyFinalDVPi0Selections(
            Init2PhotonKinematics(
                SelectPi0Event(df_afterFid_7546_pi0data),
                beam_energy));

    double luminosity_pi0MC   = 16.4 * (4.3301-0.1083) * 0.5;
    double luminosity_pi0data = 10.913182; // fb^-1

    df_final_OnlPi0_7546_pi0data =
        df_final_OnlPi0_7546_pi0data.Filter(
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

    df_final_OnlPi0_7546_pi0MC =
        DefineRecoProtonSector(df_final_OnlPi0_7546_pi0MC);

    df_final_OnlPi0_7546_pi0data =
        DefineRecoProtonSector(df_final_OnlPi0_7546_pi0data);

    // Configure proton momentum slices here.
    // Each range is [p_min, p_max) in GeV/c.
    std::vector<std::pair<double, double>> fd_proton_p_ranges = {
        //{0.367, 0.458},
        //{0.458, 0.567},
        //{0.567, 0.672},
        //{0.672, 0.775},
        //{0.775, 0.876},
        //{0.876, 1.074},
        //{1.074, 1.270},
        {0.0, 1.2},
    };

    std::vector<std::pair<double, double>> cd_proton_p_ranges = {
        //{0.2, 0.3},
        //{0.367, 0.458},
        //{0.458, 0.567},
        //{0.567, 0.672},
        //{0.672, 0.775},
        //{0.775, 0.876},
        //{0.876, 1.074},
        //{1.074, 1.270},
        {0.0, 1.2},
    };

    // ============================================================
    // FD proton: overall
    // ============================================================
    auto df_pi0MC_FD =
        df_final_OnlPi0_7546_pi0MC.Filter(
            "pro_det_region == 1",
            "Proton in FD");

    auto df_pi0Data_FD =
        df_final_OnlPi0_7546_pi0data.Filter(
            "pro_det_region == 1",
            "Proton in FD");

    PlotProtonThetaPComparisonByPRange(
        df_pi0MC_FD,
        df_pi0Data_FD,
        fd_proton_p_ranges,
        "Pi0 MC FD",
        "Pi0 Data FD",
        "ProtonThetaPPhiComparison_FD.png",
        10.0, 50.0,
        luminosity_pi0MC,
        luminosity_pi0data,
        4);

    // ============================================================
    // FD proton: sector-by-sector
    // ============================================================
    
    for (int sec = 1; sec <= 6; ++sec) {
        auto df_pi0MC_FD_sec =
            df_final_OnlPi0_7546_pi0MC.Filter(
                Form("pro_det_region == 1 && recpro_sector == %d", sec),
                Form("FD proton sector %d", sec));

        auto df_pi0Data_FD_sec =
            df_final_OnlPi0_7546_pi0data.Filter(
                Form("pro_det_region == 1 && recpro_sector == %d", sec),
                Form("FD proton sector %d", sec));

        PlotProtonThetaPComparisonByPRange(
            df_pi0MC_FD_sec,
            df_pi0Data_FD_sec,
            fd_proton_p_ranges,
            Form("Pi0 MC FD S%d", sec),
            Form("Pi0 Data FD S%d", sec),
            Form("ProtonThetaPPhiComparison_FD_S%d.png", sec),
            10.0, 50.0,
            luminosity_pi0MC,
            luminosity_pi0data,
            4);
    }
    
    // ============================================================
    // CD proton: overall
    // ============================================================
    auto df_pi0MC_CD =
        df_final_OnlPi0_7546_pi0MC.Filter(
            "pro_det_region == 2",
            "Proton in CD");

    auto df_pi0Data_CD =
        df_final_OnlPi0_7546_pi0data.Filter(
            "pro_det_region == 2",
            "Proton in CD");

    PlotProtonThetaPComparisonByPRange(
        df_pi0MC_CD,
        df_pi0Data_CD,
        cd_proton_p_ranges,
        "Pi0 MC (CD Proton)",
        "Pi0 Data (CD Proton)",
        "ProtonThetaPPhiComparison_CD.png",
        35.0, 75.0,
        luminosity_pi0MC,
        luminosity_pi0data,
        4);
}
