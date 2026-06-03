#ifndef DVCSANALYSIS_H
#define DVCSANALYSIS_H

#include <TH1F.h>

#include <ROOT/RDF/RInterface.hxx>
#include <optional>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include "../Cuts/EventCut.h"
#include "../Cuts/TrackCut.h"
#include "../Cuts/QADBCuts.h"
#include "../Correction/MomentumCorrection.h"
#include "../Math/ParticleMassTable.h"
#include "../Math/RECParticleKinematic.h"
#include "../ParticleInformation/RECCalorimeter.h"
#include "../ParticleInformation/RECParticle.h"
#include "../ParticleInformation/RECTrack.h"
#include "../ParticleInformation/RECTraj.h"
#include "../ParticleInformation/RECForwardTagger.h"
#include "../core/Columns.h"
#include "AnalysisTask.h"

class DVCSAnalysis : public AnalysisTask {
 public:
  DVCSAnalysis(bool IsMC = false, bool IsReproc = false, bool IsMinBook = false);
  virtual ~DVCSAnalysis();

  void UserCreateOutputObjects() override;
  void UserExec(ROOT::RDF::RNode &df) override;
  void SaveOutput() override;
  void SetMaxEvents(size_t n) { fMaxEvents = n; }
  // Bad for raw pointer setup
  void SetTrackCuts(std::shared_ptr<TrackCut> cuts) { fTrackCuts = std::move(cuts); };

  void SetEventCuts(EventCut *evtCuts) { fEventCuts = evtCuts; };
  void SetDoInvMassCut(bool cut) { fDoInvMassCut = cut; };
  void SetAcceptEverything(bool accept) { fAcceptAll = accept; };
 
  void SetDoFiducialCut(bool cut) { fFiducialCut = cut; };

  void SetBeamEnergy(float beam_energy) { fbeam_energy = beam_energy; };
  void SetOutputDir(const std::string &dir) override;

  void SetFTonConfig(bool config) { fFTonConfig = config; }

  void SetDoMomentumCorrection(bool do_correction) { fDoMomentumCorrection = do_correction; }
  void SetMomentumCorrection(std::shared_ptr<MomentumCorrection> corr) { fMomCorr = std::move(corr); }

  void SetQADBCuts(std::shared_ptr<QADBCuts> qadbcuts){ fQADBCuts = std::move(qadbcuts); };
  void SetDoQADBCuts(bool charge_output) { fIsQADBCut = charge_output; }

  void SetOptimizeColumns(bool optimize) { fOptimizeColumns = optimize; }
  void SetOutputWCSV(bool output) { fOutputWCSV = output; }
  void SetOutputWCSVName(const std::string& name) { fOutputWCSVName = name; }

 private:
  std::vector<std::string> MinimalColumns() const;
  void WriteFinalElectronWCSV(ROOT::RDF::RNode df, const std::string& csvPath);
  bool IsMC = false;
  bool fDoInvMassCut = false;  // Flag to indicate if invMass cut is applied
  bool fAcceptAll = false;  // Flag to indicate if all events are accepted without cuts
  bool IsMinBooking = false;  // reduces the output to minimum only after fiducial 
  bool IsReproc = false;  // Flag to indicate if fiducial cut is applied
  bool fFiducialCut = false;  // Flag to indicate if fiducial cut is applied
  bool fFTonConfig = true;
  bool fDoMomentumCorrection = false;  // Flag to indicate if momentum correction is applied
  bool fOptimizeColumns = false;  // Flag to indicate if column optimization is enabled
  bool fOutputWCSV = false;  // Output final corrected electron W CSV and skip ROOT snapshots
  std::string fOutputWCSVName = "electron_w_afterCorr.csv";

  bool fIsQADBCut = false;  // Flag to indicate if QADB cut is applied
  bool fChargeOutput = false; // Flag to indicate if output the accumulated charge from QADB

  size_t fMaxEvents{0}; // Maximum number of events to process, 0 means no limit

  std::optional<ROOT::RDF::RNode> dforginal;

  std::optional<ROOT::RDF::RNode> dfSelected;
  std::optional<ROOT::RDF::RNode> dfSelected_afterFid;  // DataFrame after fiducial cuts
  std::optional<ROOT::RDF::RNode> dfSelected_afterFid_afterCorr;  // DataFrame after fiducial cuts and momentum correction
  std::string fOutputDir;
  
  float fbeam_energy = 10.6;
  
  TH1F *fHistPhotonP = nullptr;

  std::shared_ptr<TrackCut> fTrackCuts;
  std::shared_ptr<QADBCuts> fQADBCuts;
  EventCut *fEventCuts = nullptr;
  std::shared_ptr<TrackCut> fTrackCutsNoFid;
  std::shared_ptr<TrackCut> fTrackCutsWithFid;

  std::shared_ptr<MomentumCorrection> fMomCorr = nullptr;  // Pointer to momentum correction object
  

  TFile *fOutFile = nullptr;  // Output file pointer set by manager
};

#endif
