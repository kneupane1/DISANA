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




void DrawECALSF(const int &selectedPid, const int &selecteddetector,
                                             //const std::vector<int> &layers,
                                             const std::vector<int> &sectors,
                                             const std::string &filename, const std::string &treename,
                                             const bool doFid, double nSigma, double fitpMin, double fitpMax) {
    TStopwatch timer;
    timer.Start();
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df(treename, filename);
    gStyle->SetOptStat(0);

    std::map<int, TH2F*> hist_SF, hist_Triangle;

    for (int sector : sectors) {
        hist_SF[sector] = new TH2F(Form("lw_S%d", sector), ";PCAL lw [cm];E/p", 200, 0, 10, 200, 0, 0.4);
        hist_Triangle[sector] = new TH2F(Form("lv_S%d", sector), ";PCAL lv [cm];E/p", 200, 0, 0.4, 200, 0, 0.4);
        hist_SF[sector]->SetDirectory(nullptr);
        hist_Triangle[sector]->SetDirectory(nullptr);
    }

    df.Foreach([&](const RVec<short> &det, const RVec<short> &layer,
                   const RVec<short> &sector, const RVec<int16_t> &pindex,
                   const RVec<int> &pid, const RVec<int> &passFid,
                   const RVec<float> &lv, const RVec<float> &lw,
                   const RVec<float> &energy,
                   const RVec<float> &px, const RVec<float> &py, const RVec<float> &pz, const RVec<int> &parPass) {
        struct HitInfo {
            float totalE = 0;
            float PCAL_E = 0;
            float ECIN_E = 0;
            float ECOUT_E = 0;
            float parE = 0 ;
            bool hasL1 = false;
            bool hasL4 = false;
            bool hasL7 = false;
        };

        std::map<std::pair<int16_t, int>, HitInfo> hitMap;  // key: (pindex, sector)

        for (size_t i = 0; i < det.size(); ++i) {
            if (det[i] != selecteddetector) continue;
            int lay = layer[i];
            int sec = sector[i];
            int16_t idx = pindex[i];
            if (pid[idx] != selectedPid || (doFid && !passFid[idx])) continue;
            //if (std::find(layers.begin(), layers.end(), lay) == layers.end()) continue;

            auto key = std::make_pair(idx, sec);
            hitMap[key].totalE += energy[i];
            hitMap[key].parE = std::sqrt(px[idx]*px[idx] + py[idx]*py[idx] + pz[idx]*pz[idx]);
            if (lay == 1) {
                hitMap[key].PCAL_E = energy[i];
                hitMap[key].hasL1 = true;
            }
            if (lay == 4) {
                hitMap[key].ECIN_E = energy[i];
                hitMap[key].hasL4 = true;
            }
            if (lay == 7) {
                hitMap[key].ECOUT_E = energy[i];
                hitMap[key].hasL7 = true;
            }
        }

        for (const auto &[key, info] : hitMap) {
            int16_t idx = key.first;
            int sec = key.second;

            float p = std::sqrt(px[idx]*px[idx] + py[idx]*py[idx] + pz[idx]*pz[idx]);
            if (p < 1e-3) continue;

            float eOverP = info.totalE / p;
            float PCALeOverP = info.PCAL_E / p;
            float ECINeOverP = info.ECIN_E / p;
            if (parPass[idx] && eOverP!=0) hist_SF[sec]->Fill(p, eOverP);
            if (parPass[idx] && info.hasL1 && info.hasL4) hist_Triangle[sec]->Fill(ECINeOverP, PCALeOverP);
        }
    }, {"REC_Calorimeter_detector", "REC_Calorimeter_layer", "REC_Calorimeter_sector",
        "REC_Calorimeter_pindex", "REC_Particle_pid", "REC_Track_pass_fid",
        "REC_Calorimeter_lv", "REC_Calorimeter_lw", "REC_Calorimeter_energy",
        "REC_Particle_px", "REC_Particle_py", "REC_Particle_pz","REC_Particle_pass"});

    std::string outdir = "ECALSFPlots";
    gSystem->Exec(("mkdir -p " + outdir).c_str());

    for (int sector : sectors) {
        if (!hist_SF.count(sector)) continue;

        auto grUpper = new TGraph();
        auto grLower = new TGraph();
        grUpper->SetLineColor(kRed);   grUpper->SetLineStyle(2);    grUpper->SetLineWidth(5);
        grLower->SetLineColor(kRed);   grLower->SetLineStyle(2);    grLower->SetLineWidth(5);
        int ip = 0;

        for (int ix = 1; ix <= hist_SF[sector]->GetNbinsX(); ++ix) {

            std::unique_ptr<TH1D> proj(hist_SF[sector]->ProjectionY(
                       Form("py_%d_%d",sector,ix), ix, ix, "e")); 

            if (proj->GetEntries() < 30) continue;        
            TF1 fgaus("fgaus", "gaus", 0.05, 0.35);       
            if (proj->Fit(&fgaus, "QNR") != 0) continue;  

            double mu    = fgaus.GetParameter(1);
            double sigma = fgaus.GetParameter(2);
            double pCen  = hist_SF[sector]->GetXaxis()->GetBinCenter(ix);

            grUpper->SetPoint(ip, pCen, mu + nSigma*sigma);
            grLower->SetPoint(ip, pCen, mu - nSigma*sigma);
            ++ip;
        }

        double pMin = fitpMin;
        double pMax = fitpMax;
  
        TF1 fUpFit(Form("fUpFit_S%d",sector),"[0]+[1]*x+[2]*(x*x)",pMin,pMax);
        TF1 fLoFit(Form("fLoFit_S%d",sector),"[0]+[1]*x+[2]*(x*x)",pMin,pMax);
        std::cout << "Sector " << sector << " fit range: " << pMin << " to " << pMax << std::endl;
        fUpFit.SetLineColor(kBlack);   fUpFit.SetLineWidth(5); fUpFit.SetLineStyle(7);
        fLoFit.SetLineColor(kBlack);   fLoFit.SetLineWidth(5); fLoFit.SetLineStyle(7);

        grUpper->Fit(&fUpFit,"QR");
        grLower->Fit(&fLoFit,"QR");

        std::ofstream fout(Form("ECALSFPlots/fitParams_sector%d.txt",sector));
        auto dump = [&](const char* tag, TF1& f){
            std::cout << "p[0]+p[1]*x+p[2]*x^2 "<<std::endl;
            std::cout << Form("Sector %d  %s  a=%.4g±%.3g  b=%.4g±%.3g  c=%.4g±%.3g\n",
                    sector,tag,
                    f.GetParameter(0),f.GetParError(0),
                    f.GetParameter(1),f.GetParError(1),
                    f.GetParameter(2),f.GetParError(2));
            fout << tag << "p[0]+p[1]*x+p[2]*x^2 \n "
                 << f.GetParameter(0) << ", "
                 << f.GetParameter(1) << ", "
                 << f.GetParameter(2) << " "<<"\n";
        };
        dump("Upper",fUpFit);
        dump("Lower",fLoFit);
        fout.close();


        TCanvas *c = new TCanvas(Form("c_ecal_EoverP_S%d", sector), "", 2000, 1500);
        if (selectedPid == 11) {
            hist_SF[sector]->SetTitle(Form("Electron E_{dep}/p in Sector %d; p [GeV];E_{dep}/p", sector));
        }
        else if (selectedPid == 22) {
            hist_SF[sector]->SetTitle(Form("Photon E_{dep}/p in Sector %d; p [GeV];E_{dep}/p", sector));
        } else {
            hist_SF[sector]->SetTitle(Form("PID %d E_{dep}/p in Sector %d; p [GeV];E_{dep}/p", selectedPid, sector));
        }
        hist_SF[sector]->GetYaxis()->SetRangeUser(0.08, 0.352);
        hist_SF[sector]->GetXaxis()->SetRangeUser(fitpMin, fitpMax);
        gPad->SetLogz();
        hist_SF[sector]->Draw("COLZ");
        grUpper->SetLineWidth(10);
        grLower->SetLineWidth(10);
        fUpFit.SetLineWidth(5);
        fLoFit.SetLineWidth(5);
        grUpper->Draw("L SAME");
        grLower->Draw("L SAME");
        fUpFit.Draw("L SAME");      // --- NEW ---
        fLoFit.Draw("L SAME");

        std::string outname;
        if (selectedPid == 11) {
            outname = Form("%s/electron_EoverP_P_sector%d%s.png", outdir.c_str(), sector, doFid ? "" : "_noFid");
        } else if (selectedPid == 22) {
            outname = Form("%s/photon_EoverP_P_sector%d%s.png", outdir.c_str(), sector, doFid ? "" : "_noFid");
        } else {
            outname = Form("%s/pid%d_EoverP_P_sector%d%s.png", outdir.c_str(), selectedPid, sector, doFid ? "" : "_noFid");
        }

        c->SaveAs(outname.c_str());
        std::cout << "Saved: " << outname << std::endl;
        delete c;

        if (!hist_Triangle.count(sector)) continue;

        TCanvas *c2 = new TCanvas(Form("c_ecal_triangle_S%d", sector), "", 2000, 1500);
        if (selectedPid == 11) {
            hist_Triangle[sector]->SetTitle(Form("Electron E_{dep} in Sector %d; ECIN E_{dep}/p; PCAL E_{dep}/p", sector));
        }
        else if (selectedPid == 22) {
            hist_Triangle[sector]->SetTitle(Form("Photon E_{dep} in Sector %d; ECIN E_{dep}/p; PCAL E_{dep}/p", sector));
        } else {
            hist_Triangle[sector]->SetTitle(Form("PID %d E_{dep} in Sector %d; ECIN E_{dep}/p; PCAL E_{dep}/p", selectedPid, sector));
        }
        hist_Triangle[sector]->GetYaxis()->SetRangeUser(0.00, 0.3);
        hist_Triangle[sector]->GetXaxis()->SetRangeUser(0.00, 0.3);
        gPad->SetLogz();
        hist_Triangle[sector]->Draw("COLZ");

        std::string outname2;
        if (selectedPid == 11) {
            outname2 = Form("%s/electron_Edep_triangle_sector%d%s.png", outdir.c_str(), sector, doFid ? "" : "_noFid");
        } else if (selectedPid == 22) {
            outname2 = Form("%s/photon_Edep_triangle_sector%d%s.png", outdir.c_str(), sector, doFid ? "" : "_noFid");
        } else {
            outname2 = Form("%s/pid%d_Edep_triangle_sector%d%s.png", outdir.c_str(), selectedPid, sector, doFid ? "" : "_noFid");
        }
        c2->SaveAs(outname2.c_str());
        std::cout << "Saved: " << outname2 << std::endl;
        delete c2;
    }

    timer.Stop();
    std::cout << "Time for DrawECALSF: ";
    timer.Print();
}




void analysisECALSF() {
    ROOT::EnableImplicitMT();

    std::string path = "../build/File8477/"; 
    //std::string path = "../../../analysis1101/DISANA/build/data/";
    
    //==================// DVCS analysis paths==========================
    std::vector<int> layers = {1, 4, 7};
    std::vector<int> sectors = {1, 2, 3, 4, 5, 6};
    DrawECALSF(11, 7, sectors, path + "dfSelected_afterFid_afterCorr.root", "dfSelected_afterFid_afterCorr", true, 3, 2, 7.8);
    //DrawECALSF(11, 7, sectors, path + "dfSelected.root", "dfSelected", false, 3, 2, 7.8);
    gApplication->Terminate(0);
}
