// system include files
#include <iostream>
#include <memory>
#include "TMath.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"  
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include<algorithm>
#define Pi 3.141593
#include "Math/VectorUtil.h"
#include "TMath.h"
#include <TFormula.h>

struct sortPt
{
   bool operator()(TLorentzVector* s1, TLorentzVector* s2) const
   {
      return s1->Pt() >= s2->Pt();
   }
} mysortPt;
//
// class declaration
//

class PKUTreeMaker : public edm::EDAnalyzer {
public:
  explicit PKUTreeMaker(const edm::ParameterSet&);
  ~PKUTreeMaker();
  //static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
    
    enum PhotonMatchType {UNMATCHED = 0,
        MATCHED_FROM_GUDSCB,
        MATCHED_FROM_PI0,
        MATCHED_FROM_OTHER_SOURCES};
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  virtual void beginRun(const edm::Run&, const edm::EventSetup&) override;
  virtual void endRun(const edm::Run&, const edm::EventSetup&) override;
  virtual void addTypeICorr( edm::Event const & event );
  virtual double getJEC( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ );
  virtual double getJECOffset( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ );
  math::XYZTLorentzVector getNeutrinoP4(double& MetPt, double& MetPhi, TLorentzVector& lep, int lepType);
  bool hasMatchedPromptElectron(const reco::SuperClusterRef &sc, const edm::Handle<edm::View<pat::Electron> > &eleCol,
                                  const edm::Handle<reco::ConversionCollection> &convCol, const math::XYZPoint &beamspot,
                                  float lxyMin=2.0, float probMin=1e-6, unsigned int nHitsBeforeVtxMax=0);
  int matchToTruth(const reco::Photon &pho,
                     const edm::Handle<edm::View<reco::GenParticle>>  &genParticles, bool &ISRPho, double &dR);
    
  void findFirstNonPhotonMother(const reco::Candidate *particle,
                                  int &ancestorPID, int &ancestorStatus);
  float EAch(float x); 
  float EAnh(float x);
  float EApho(float x);
  std::vector<std::string> offsetCorrLabel_;
  FactorizedJetCorrector* jecOffset_;
  std::vector<std::string> jetCorrLabel_;
  edm::Handle< double >  rho_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<pat::METCollection>  metInputToken_;
  std::vector<edm::EDGetTokenT<pat::METCollection>> mettokens;
  edm::EDGetTokenT<pat::METCollection> metToken_;
  edm::EDGetTokenT<edm::View<pat::Electron> > electronToken_ ;
  edm::EDGetTokenT<edm::View<pat::Photon> > photonToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<std::vector<reco::Conversion> > conversionsToken_;
  edm::EDGetTokenT<edm::View<pat::Electron> > looseelectronToken_ ; 
  edm::EDGetTokenT<edm::View<pat::Muon> > loosemuonToken_; 

// Filter
  edm::EDGetTokenT<edm::TriggerResults> 		     noiseFilterToken_;
  edm::Handle< edm::TriggerResults> 			     noiseFilterBits_;
  std::string HBHENoiseFilter_Selector_;
  edm::EDGetTokenT<bool> EarlyRunsHBHENoiseFilter_Selector_;
  edm::EDGetTokenT<bool> HBHENoiseFilterResultToken_;
  std::string CSCHaloNoiseFilter_Selector_;
  std::string HCALlaserNoiseFilter_Selector_;
  std::string ECALDeadCellNoiseFilter_Selector_;
  std::string GoodVtxNoiseFilter_Selector_;
  std::string TrkFailureNoiseFilter_Selector_;
  std::string EEBadScNoiseFilter_Selector_;
  std::string ECALlaserNoiseFilter_Selector_;
  std::string TrkPOGNoiseFilter_Selector_;
  std::string TrkPOG_manystrip_NoiseFilter_Selector_;
  std::string TrkPOG_toomanystrip_NoiseFilter_Selector_;
  std::string TrkPOG_logError_NoiseFilter_Selector_;
  std::string METFilters_Selector_;


  // ----------member data ---------------------------
  TTree* outTree_;

  double MW_; //Jing
  int nevent, run, ls;
  int nVtx;
  double triggerWeight, lumiWeight, pileupWeight;
  double theWeight;
  double  nump=0.;
  double  numm=0.;

  double ptVlep, yVlep, phiVlep, massVlep, mtVlep;
  double ptVlepJEC, yVlepJEC, phiVlepJEC, massVlepJEC, mtVlepJEC;
  double ptlep1, etalep1, philep1;
  int  lep, nlooseeles,nloosemus;
  double met, metPhi, j1metPhi, j2metPhi;
  //Met JEC
  double METraw_et, METraw_phi, METraw_sumEt;
  double MET_et, MET_phi, MET_sumEt, MET_corrPx, MET_corrPy;
  double useless;
  // AK4 Jets
  double ak4jet_pt[6],ak4jet_eta[6],ak4jet_phi[6],ak4jet_e[6];
  double ak4jet_pt_jer[6];
  double ak4jet_csv[6],ak4jet_icsv[6];
  double drjetlep[6], drjetphoton[6];
  //Photon
  double photon_pt[6],photon_eta[6],photon_phi[6],photon_e[6];
  double drphotonlep[6];
  double photonet, photoneta, photonphi, photone;
  double photonsieie, photonphoiso, photonchiso, photonnhiso;
  int iphoton;
  double drla;
  bool passEleVeto;
  //Photon gen match
  std::vector<Int_t>   isTrue_;
  bool ISRPho;
  double dR_;
  //Jets
  double jet1pt, jet1eta, jet1phi, jet1e, jet1csv, jet1icsv;
  double jet2pt, jet2eta, jet2phi, jet2e, jet2csv, jet2icsv;
  double drj1a, drj2a, drj1l, drj2l;
  double Mjj, deltaeta, zepp;
   
  edm::InputTag electronIdTag_;
  void setDummyValues();
    
  /// Parameters to steer the treeDumper
  int originalNEvents_;
  double crossSectionPb_;
  double targetLumiInvPb_;
  std::string PKUChannel_;
  bool isGen_ , RunOnMC_;
  std::string leptonicVSrc_;
  std::string ak4jetsSrc_;
  std::vector<std::string> jecAK4Labels_;
  //correction jet
  FactorizedJetCorrector* jecAK4_;
  std::string photonSrc_;
  std::string gravitonSrc_, metSrc_;
  std::map<std::string,double>  TypeICorrMap_;
  edm::InputTag mets_;

  //High Level Trigger
  HLTConfigProvider hltConfig;
  edm::EDGetTokenT<edm::TriggerResults> hltToken_;
  std::vector<std::string> elPaths_;
  std::vector<std::string> muPaths_;
  std::vector<std::string> elPaths;
  std::vector<std::string> muPaths;
  int  HLT_Ele;
  int  HLT_Mu;

// filter
  bool passFilter_HBHE_                   ;
  bool passFilter_CSCHalo_                ;
  bool passFilter_HCALlaser_              ;
  bool passFilter_ECALDeadCell_           ;
  bool passFilter_GoodVtx_                ;
  bool passFilter_TrkFailure_             ;
  bool passFilter_EEBadSc_                ;
  bool passFilter_ECALlaser_              ;
  bool passFilter_TrkPOG_                 ;
  bool passFilter_TrkPOG_manystrip_       ;
  bool passFilter_TrkPOG_toomanystrip_    ;
  bool passFilter_TrkPOG_logError_        ;
  bool passFilter_METFilters_             ;

};

float PKUTreeMaker::EAch( float x){
 float EA = 0.0158;
 if(x>1.0)   EA = 0.0143;
 if(x>1.479) EA = 0.0115;
 if(x>2.0)   EA = 0.0094;
 if(x>2.2)   EA = 0.0095;
 if(x>2.3)   EA = 0.0068;
 if(x>2.4)   EA = 0.0053;
 return EA;
}

float PKUTreeMaker::EAnh( float x){
 float EA = 0.0143;
 if(x>1.0)   EA = 0.0210;
 if(x>1.479) EA = 0.0147;
 if(x>2.0)   EA = 0.0082;
 if(x>2.2)   EA = 0.0124;
 if(x>2.3)   EA = 0.0186;
 if(x>2.4)   EA = 0.0320;
 return EA;
}

float PKUTreeMaker::EApho( float x){
 float EA = 0.0725;
 if(x>1.0)   EA = 0.0604;
 if(x>1.479) EA = 0.0320;
 if(x>2.0)   EA = 0.0512;
 if(x>2.2)   EA = 0.0766;
 if(x>2.3)   EA = 0.0949;
 if(x>2.4)   EA = 0.1160;
 return EA;
}


//
// constructors and destructor
//
PKUTreeMaker::PKUTreeMaker(const edm::ParameterSet& iConfig):
  hltToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("hltToken"))),
  elPaths_(iConfig.getParameter<std::vector<std::string>>("elPaths")),
  muPaths_(iConfig.getParameter<std::vector<std::string>>("muPaths"))
{
  originalNEvents_ = iConfig.getParameter<int>("originalNEvents");
  crossSectionPb_  = iConfig.getParameter<double>("crossSectionPb");
  targetLumiInvPb_ = iConfig.getParameter<double>("targetLumiInvPb");
  PKUChannel_     = iConfig.getParameter<std::string>("PKUChannel");
  isGen_           = iConfig.getParameter<bool>("isGen");
  RunOnMC_           = iConfig.getParameter<bool>("RunOnMC");
  leptonicVSrc_ = iConfig.getParameter<std::string>("leptonicVSrc");
  rhoToken_  = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));
  ak4jetsSrc_      = iConfig.getParameter<std::string>("ak4jetsSrc");
  jecAK4Labels_   =  iConfig.getParameter<std::vector<std::string>>("jecAK4chsPayloadNames");
  metSrc_          = iConfig.getParameter<std::string>("metSrc");
  metToken_ = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metSrc"));
  mettokens.push_back( metToken_ );
  metInputToken_ = mettokens[0];
  electronIdTag_   = iConfig.getParameter<edm::InputTag>("electronIDs");
  electronToken_    = (consumes<edm::View<pat::Electron> > (iConfig.getParameter<edm::InputTag>("electrons")))            ;
  looseelectronToken_    = (consumes<edm::View<pat::Electron> > (iConfig.getParameter<edm::InputTag>("looseelectronSrc"))) ;
  loosemuonToken_    = (consumes<edm::View<pat::Muon> > (iConfig.getParameter<edm::InputTag>("loosemuonSrc")))              ;
  photonSrc_      = iConfig.getParameter<std::string>("photonSrc");
  beamSpotToken_    = (consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))) ;
  conversionsToken_ = (consumes<std::vector<reco::Conversion> >(iConfig.getParameter<edm::InputTag>("conversions"))) ;
    
   
   jetCorrLabel_ = jecAK4Labels_;
   offsetCorrLabel_.push_back(jetCorrLabel_[0]);

// filter
   noiseFilterToken_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("noiseFilter"));
   HBHENoiseFilterResultToken_ = consumes<bool>(iConfig.getParameter<edm::InputTag>("noiseFilterSelection_EarlyRunsHBHENoiseFilter"));
   HBHENoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_HBHENoiseFilter");
   CSCHaloNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_CSCTightHaloFilter");
   HCALlaserNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_hcalLaserEventFilter");
   ECALDeadCellNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter");
   GoodVtxNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_goodVertices");
   TrkFailureNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_trackingFailureFilter");
   EEBadScNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_eeBadScFilter");
   ECALlaserNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_ecalLaserCorrFilter");
   TrkPOGNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_trkPOGFilters");
   TrkPOG_manystrip_NoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_trkPOG_manystripclus53X");
   TrkPOG_toomanystrip_NoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_trkPOG_toomanystripclus53X");
   TrkPOG_logError_NoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_trkPOG_logErrorTooManyClusters");
   METFilters_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_metFilters");


   MW_=80.385; //Jing 
  //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  outTree_ = fs->make<TTree>("PKUCandidates","PKU Candidates");

  /// Basic event quantities
  outTree_->Branch("event"           ,&nevent         ,"event/I"          );
  outTree_->Branch("nVtx"            ,&nVtx           ,"nVtx/I"           );
  outTree_->Branch("theWeight"           ,&theWeight         ,"theWeight/D"          );
  outTree_->Branch("nump"           ,&nump         ,"nump/D"          );
  outTree_->Branch("numm"           ,&numm         ,"numm/D"          );
  outTree_->Branch("lep"             ,&lep            ,"lep/I"            );
  outTree_->Branch("ptVlep"          ,&ptVlep         ,"ptVlep/D"         );
  outTree_->Branch("yVlep"           ,&yVlep          ,"yVlep/D"          );
  outTree_->Branch("phiVlep"         ,&phiVlep        ,"phiVlep/D"        );
  outTree_->Branch("massVlep"        ,&massVlep       ,"massVlep/D"       );
  outTree_->Branch("mtVlep"          ,&mtVlep         ,"mtVlep/D"         );
  outTree_->Branch("ptVlepJEC"          ,&ptVlepJEC         ,"ptVlepJEC/D"         );
  outTree_->Branch("yVlepJEC"           ,&yVlepJEC          ,"yVlepJEC/D"          );
  outTree_->Branch("phiVlepJEC"         ,&phiVlepJEC        ,"phiVlepJEC/D"        );
  outTree_->Branch("massVlepJEC"        ,&massVlepJEC       ,"massVlepJEC/D"       );
  outTree_->Branch("mtVlepJEC"          ,&mtVlepJEC         ,"mtVlepJEC/D"         );
  outTree_->Branch("nlooseeles"          ,&nlooseeles         ,"nlooseeles/I"         );
  outTree_->Branch("nloosemus"          ,&nloosemus         ,"nloosemus/I"         );

  /// Photon
  outTree_->Branch("photon_pt"        , photon_pt       ,"photon_pt[6]/D"       );
  outTree_->Branch("photon_eta"        , photon_eta       ,"photon_eta[6]/D"       );
  outTree_->Branch("photon_phi"        , photon_phi       ,"photon_phi[6]/D"       );
  outTree_->Branch("photon_e"        , photon_e       ,"photon_e[6]/D"       );
  outTree_->Branch("passEleVeto"        , &passEleVeto       ,"passEleVeto/O"       );
  outTree_->Branch("photonet"          ,&photonet         ,"photonet/D"         );
  outTree_->Branch("photoneta"          ,&photoneta         ,"photoneta/D"         );
  outTree_->Branch("photonphi"          ,&photonphi         ,"photonphi/D"         );
  outTree_->Branch("photone"          ,&photone         ,"photone/D"         );
  outTree_->Branch("photonsieie"          ,&photonsieie         ,"photonsieie/D"         );
  outTree_->Branch("photonphoiso"          ,&photonphoiso         ,"photonphoiso/D"         );
  outTree_->Branch("photonchiso"          ,&photonchiso         ,"photonchiso/D"         );
  outTree_->Branch("photonnhiso"          ,&photonnhiso         ,"photonnhiso/D"         );
  outTree_->Branch("iphoton"             ,&iphoton            ,"iphoton/I"            );
  outTree_->Branch("drla"          ,&drla         ,"drla/D"         );
    //photon gen match
    outTree_->Branch("dR"    , &dR_, "dR/D");
    outTree_->Branch("ISRPho"        , &ISRPho       ,"ISRPho/O"       );
    outTree_->Branch("isTrue", &isTrue_);
//jets
  outTree_->Branch("jet1pt"          ,&jet1pt         ,"jet1pt/D"         );
  outTree_->Branch("jet1eta"          ,&jet1eta         ,"jet1eta/D"         );
  outTree_->Branch("jet1phi"          ,&jet1phi         ,"jet1phi/D"         );
  outTree_->Branch("jet1e"          ,&jet1e         ,"jet1e/D"         );
  outTree_->Branch("jet1csv"          ,&jet1csv         ,"jet1csv/D"         );
  outTree_->Branch("jet1icsv"          ,&jet1icsv         ,"jet1icsv/D"         );
  outTree_->Branch("jet2pt"          ,&jet2pt         ,"jet2pt/D"         );
  outTree_->Branch("jet2eta"          ,&jet2eta         ,"jet2eta/D"         );
  outTree_->Branch("jet2phi"          ,&jet2phi         ,"jet2phi/D"         );
  outTree_->Branch("jet2e"          ,&jet2e         ,"jet2e/D"         );
  outTree_->Branch("jet2csv"          ,&jet2csv         ,"jet2csv/D"         );
  outTree_->Branch("jet2icsv"          ,&jet2icsv         ,"jet2icsv/D"         );
  outTree_->Branch("drj1a"          ,&drj1a         ,"drj1a/D"         );
  outTree_->Branch("drj2a"          ,&drj2a         ,"drj2a/D"         );
  outTree_->Branch("drj1l"          ,&drj1l         ,"drj1l/D"         );
  outTree_->Branch("drj2l"          ,&drj2l         ,"drj2l/D"         );
  outTree_->Branch("Mjj"          ,&Mjj         ,"Mjj/D"         );
  outTree_->Branch("deltaeta"          ,&deltaeta         ,"deltaeta/D"         );
  outTree_->Branch("zepp"          ,&zepp         ,"zepp/D"         );
  /// Generic kinematic quantities
  outTree_->Branch("ptlep1"          ,&ptlep1         ,"ptlep1/D"         );
  outTree_->Branch("etalep1"         ,&etalep1        ,"etalep1/D"        );
  outTree_->Branch("philep1"         ,&philep1        ,"philep1/D"        );
  outTree_->Branch("met"             ,&met            ,"met/D"            );
  outTree_->Branch("metPhi"          ,&metPhi         ,"metPhi/D"         );
  outTree_->Branch("j1metPhi"          ,&j1metPhi         ,"j1metPhi/D"         );
  outTree_->Branch("j2metPhi"          ,&j2metPhi         ,"j2metPhi/D"         );
   
  outTree_->Branch("METraw_et",&METraw_et,"METraw_et/D");
  outTree_->Branch("METraw_phi",&METraw_phi,"METraw_phi/D");
  outTree_->Branch("METraw_sumEt",&METraw_sumEt,"METraw_sumEt/D");
  outTree_->Branch("MET_et",&MET_et,"MET_et/D");
  outTree_->Branch("MET_phi",&MET_phi,"MET_phi/D");
  outTree_->Branch("MET_sumEt",&MET_sumEt,"MET_sumEt/D");
  outTree_->Branch("MET_corrPx",&MET_corrPx,"MET_corrPx/D");
  outTree_->Branch("MET_corrPy",&MET_corrPy,"MET_corrPy/D");

  //HLT bits
  outTree_->Branch("HLT_Ele"  ,&HLT_Ele ,"HLT_Ele/I" );
  outTree_->Branch("HLT_Mu"   ,&HLT_Mu  ,"HLT_Mu/I"  );
  
  // filter
  outTree_->Branch("passFilter_HBHE"                 ,&passFilter_HBHE_                ,"passFilter_HBHE_/O");
  outTree_->Branch("passFilter_CSCHalo"              ,&passFilter_CSCHalo_             ,"passFilter_CSCHalo_/O");
  outTree_->Branch("passFilter_HCALlaser"            ,&passFilter_HCALlaser_           ,"passFilter_HCALlaser_/O");
  outTree_->Branch("passFilter_ECALDeadCell"         ,&passFilter_ECALDeadCell_        ,"passFilter_ECALDeadCell_/O");
  outTree_->Branch("passFilter_GoodVtx"              ,&passFilter_GoodVtx_             ,"passFilter_GoodVtx_/O");
  outTree_->Branch("passFilter_TrkFailure"           ,&passFilter_TrkFailure_          ,"passFilter_TrkFailure_/O");
  outTree_->Branch("passFilter_EEBadSc"              ,&passFilter_EEBadSc_             ,"passFilter_EEBadSc_/O");
  outTree_->Branch("passFilter_ECALlaser"            ,&passFilter_ECALlaser_           ,"passFilter_ECALlaser_/O");
  outTree_->Branch("passFilter_TrkPOG"               ,&passFilter_TrkPOG_              ,"passFilter_TrkPOG_/O");
  outTree_->Branch("passFilter_TrkPOG_manystrip"     ,&passFilter_TrkPOG_manystrip_    ,"passFilter_TrkPOG_manystrip_/O");
  outTree_->Branch("passFilter_TrkPOG_toomanystrip"  ,&passFilter_TrkPOG_toomanystrip_ ,"passFilter_TrkPOG_toomanystrip_/O");
  outTree_->Branch("passFilter_TrkPOG_logError"      ,&passFilter_TrkPOG_logError_     ,"passFilter_TrkPOG_logError_/O");
  outTree_->Branch("passFilter_METFilters"           ,&passFilter_METFilters_          ,"passFilter_METFilters_/O");

  /// Other quantities
  outTree_->Branch("triggerWeight"   ,&triggerWeight  ,"triggerWeight/D"  );
  outTree_->Branch("lumiWeight"      ,&lumiWeight     ,"lumiWeight/D"     );
  outTree_->Branch("pileupWeight"    ,&pileupWeight   ,"pileupWeight/D"   );
}


double PKUTreeMaker::getJEC( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ ){
    std::vector<JetCorrectorParameters> vPar;
    //         vPar.clear();
    for ( std::vector<std::string>::const_iterator payloadBegin = jecAK4Labels_.begin(), payloadEnd = jecAK4Labels_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
        JetCorrectorParameters pars(*ipayload);
        vPar.push_back(pars);
    }
    jecAK4_ = new FactorizedJetCorrector(vPar);
    double jetCorrFactor = 1.;
    if ( fabs(rawJetP4.eta()) < jetCorrEtaMax ){
        jecAK4_->setJetEta( rawJetP4.eta() );
        jecAK4_->setJetPt ( rawJetP4.pt() );
        jecAK4_->setJetE  ( rawJetP4.energy() );
        jecAK4_->setJetPhi( rawJetP4.phi()    );
        jecAK4_->setJetA  ( jet.jetArea() );
        jecAK4_->setRho   ( *(rho_.product()) );
        jecAK4_->setNPV   ( nVtx );
        jetCorrFactor = jecAK4_->getCorrection();
    }
    reco::Candidate::LorentzVector corrJetP4 = rawJetP4;
    corrJetP4 *= jetCorrFactor;
    return jetCorrFactor;
}

double PKUTreeMaker::getJECOffset( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ ){
    std::vector<JetCorrectorParameters> vPar;
    //         vPar.clear();
    for ( std::vector<std::string>::const_iterator payloadBegin = offsetCorrLabel_.begin(), payloadEnd = offsetCorrLabel_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
        JetCorrectorParameters pars(*ipayload);
        vPar.push_back(pars);
    }
    jecOffset_ = new FactorizedJetCorrector(vPar);
    double jetCorrFactor = 1.;
    if ( fabs(rawJetP4.eta()) < jetCorrEtaMax ){
        jecOffset_->setJetEta( rawJetP4.eta()     );
        jecOffset_->setJetPt ( rawJetP4.pt()      );
        jecOffset_->setJetE  ( rawJetP4.energy()  );
        jecOffset_->setJetPhi( rawJetP4.phi()     );
        jecOffset_->setJetA  ( jet.jetArea()      );
        jecOffset_->setRho   ( *(rho_.product())  );
        jecOffset_->setNPV   ( nVtx  );
        jetCorrFactor = jecOffset_->getCorrection();
    }
    reco::Candidate::LorentzVector corrJetP4 = rawJetP4;
    corrJetP4 *= jetCorrFactor;
    return jetCorrFactor;
}


void PKUTreeMaker::addTypeICorr( edm::Event const & event ){
    TypeICorrMap_.clear();
    edm::Handle<pat::JetCollection> jets_;
    event.getByLabel("slimmedJets", jets_);
    event.getByToken(rhoToken_      , rho_     );
    edm::Handle<reco::VertexCollection> vertices_;
    event.getByLabel("offlineSlimmedPrimaryVertices", vertices_);
    edm::Handle<edm::View<pat::Muon>> muons_;
    event.getByLabel("slimmedMuons",muons_);

    bool skipEM_                    = true;
    double skipEMfractionThreshold_ = 0.9;
    bool skipMuons_                 = true;
    double jetCorrEtaMax_           = 9.9;
    double type1JetPtThreshold_     = 10.0;
    double corrEx    = 0;
    double corrEy    = 0;
    double corrSumEt = 0;
    
    for (const pat::Jet &jet : *jets_) {
        double emEnergyFraction = jet.chargedEmEnergyFraction() + jet.neutralEmEnergyFraction();
        if ( skipEM_ && emEnergyFraction > skipEMfractionThreshold_ ) continue;
        
        reco::Candidate::LorentzVector rawJetP4 = jet.correctedP4(0);
        double corr = getJEC(rawJetP4, jet, jetCorrEtaMax_, jetCorrLabel_);
        if ( skipMuons_ && jet.muonMultiplicity() != 0 ) {
            for (const pat::Muon &muon : *muons_) {
                if( !muon.isGlobalMuon() && !muon.isStandAloneMuon() ) continue;
                TLorentzVector muonV; muonV.SetPtEtaPhiE(muon.p4().pt(),muon.p4().eta(),muon.p4().phi(),muon.p4().e());
                TLorentzVector jetV; jetV.SetPtEtaPhiE(jet.p4().pt(),jet.p4().eta(),jet.p4().phi(),jet.p4().e());
                if( muonV.DeltaR(jetV) < 0.5 ){
                    reco::Candidate::LorentzVector muonP4 = muon.p4();
                    rawJetP4 -= muonP4;
                }
            }
        }
        reco::Candidate::LorentzVector corrJetP4 = corr*rawJetP4;
        if ( corrJetP4.pt() > type1JetPtThreshold_ ) {
            reco::Candidate::LorentzVector tmpP4 = jet.correctedP4(0);
            corr = getJECOffset(tmpP4, jet, jetCorrEtaMax_, offsetCorrLabel_);
            reco::Candidate::LorentzVector rawJetP4offsetCorr = corr*rawJetP4;
            corrEx    -= (corrJetP4.px() - rawJetP4offsetCorr.px());
            corrEy    -= (corrJetP4.py() - rawJetP4offsetCorr.py());
            corrSumEt += (corrJetP4.Et() - rawJetP4offsetCorr.Et());
        }
    }
    TypeICorrMap_["corrEx"]    = corrEx;
    TypeICorrMap_["corrEy"]    = corrEy;
    TypeICorrMap_["corrSumEt"] = corrSumEt;
}

//-------------------------------------------------------------------------------------------------------------------------------------//
math::XYZTLorentzVector
PKUTreeMaker::getNeutrinoP4(double& MetPt, double& MetPhi, TLorentzVector& lep, int lepType){ // Jing
    double leppt = lep.Pt();
    double lepphi = lep.Phi();
    double lepeta = lep.Eta();
    double lepenergy = lep.Energy();
    
    double metpt = MetPt;
    double metphi = MetPhi;
    
    double  px = metpt*cos(metphi);
    double  py = metpt*sin(metphi);
    double  pz = 0;
    double  pxl= leppt*cos(lepphi);
    double  pyl= leppt*sin(lepphi);
    double  pzl= leppt*sinh(lepeta);
    double  El = lepenergy;
    double  a = pow(MW_,2) + pow(px+pxl,2) + pow(py+pyl,2) - px*px - py*py - El*El + pzl*pzl;
    double  b = 2.*pzl;
    double  A = b*b -4.*El*El;
    double  B = 2.*a*b;
    double  C = a*a-4.*(px*px+py*py)*El*El;
    
    ///////////////////////////pz for fnal
    double M_mu =  0;
    
    //if(lepType==1)M_mu=0.105658367;//mu
    //if(lepType==0)M_mu=0.00051099891;//electron
    
    int type=2; // use the small abs real root
    
    a = MW_*MW_ - M_mu*M_mu + 2.0*pxl*px + 2.0*pyl*py;
    A = 4.0*(El*El - pzl*pzl);
    B = -4.0*a*pzl;
    C = 4.0*El*El*(px*px + py*py) - a*a;
    
    
    double tmproot = B*B - 4.0*A*C;
    
    if (tmproot<0) {
        //std::cout << "Complex root detected, taking real part..." << std::endl;
        pz = - B/(2*A); // take real part of complex roots
    }
    else {
        
        double tmpsol1 = (-B + sqrt(tmproot))/(2.0*A);
        double tmpsol2 = (-B - sqrt(tmproot))/(2.0*A);
        
        //std::cout << " Neutrino Solutions: " << tmpsol1 << ", " << tmpsol2 << std::endl;
        
        if (type == 0 ) {
            // two real roots, pick the one closest to pz of muon
            if (TMath::Abs(tmpsol2-pzl) < TMath::Abs(tmpsol1-pzl)) { pz = tmpsol2; }
            else { pz = tmpsol1; }
            // if pz is > 300 pick the most central root
            if ( abs(pz) > 300. ) {
                if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) { pz = tmpsol1; }
                else { pz = tmpsol2; }
            }
        }
        if (type == 1 ) {
            // two real roots, pick the one closest to pz of muon
            if (TMath::Abs(tmpsol2-pzl) < TMath::Abs(tmpsol1-pzl)) { pz = tmpsol2; }
            else {pz = tmpsol1; }
        }
        if (type == 2 ) {
            // pick the most central root.
            if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) { pz = tmpsol1; }
            else { pz = tmpsol2; }
        }
        /*if (type == 3 ) {
         // pick the largest value of the cosine
         TVector3 p3w, p3mu;
         p3w.SetXYZ(pxl+px, pyl+py, pzl+ tmpsol1);
         p3mu.SetXYZ(pxl, pyl, pzl );
         
         double sinthcm1 = 2.*(p3mu.Perp(p3w))/MW_;
         p3w.SetXYZ(pxl+px, pyl+py, pzl+ tmpsol2);
         double sinthcm2 = 2.*(p3mu.Perp(p3w))/MW_;
         
         double costhcm1 = sqrt(1. - sinthcm1*sinthcm1);
         double costhcm2 = sqrt(1. - sinthcm2*sinthcm2);
         
         if ( costhcm1 > costhcm2 ) { pz = tmpsol1; otherSol_ = tmpsol2; }
         else { pz = tmpsol2;otherSol_ = tmpsol1; }
         
         }*///end of type3
        
    }//endl of if real root
    
    //dont correct pt neutrino
    math::XYZTLorentzVector outP4(px,py,pz,sqrt(px*px+py*py+pz*pz));
    return outP4;
    
}//end neutrinoP4
//---------------------------------------------------


bool PKUTreeMaker::hasMatchedPromptElectron(const reco::SuperClusterRef &sc, const edm::Handle<edm::View<pat::Electron> > &eleCol, const edm::Handle<reco::ConversionCollection> &convCol, const math::XYZPoint &beamspot,  float lxyMin, float probMin, unsigned int nHitsBeforeVtxMax) {
    //check if a given SuperCluster matches to at least one GsfElectron having zero expected inner hits
    //and not matching any conversion in the collection passing the quality cuts
    if (sc.isNull()) return false;
    for (edm::View<pat::Electron>::const_iterator it = eleCol->begin(); it!=eleCol->end(); ++it) {
        //match electron to supercluster
        if (it->superCluster()!=sc) continue;
        //check expected inner hits
        if (it->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0) continue;
        //check if electron is matching to a conversion
        if (ConversionTools::hasMatchedConversion(*it,convCol,beamspot)) continue;
        return true;
    }
    return false;
}

int PKUTreeMaker::matchToTruth(const reco::Photon &pho,
                                      const edm::Handle<edm::View<reco::GenParticle>>
                                      &genParticles, bool &ISRPho, double &dR)
{
    //
    // Explicit loop and geometric matching method
    //
    
    // Find the closest status 1 gen photon to the reco photon
  //  double dR = 999;
    const reco::Candidate *closestPhoton = 0;
  //std::cout<<"genParticles->size() = "<<genParticles->size()<<std::endl;
    for(size_t i=0; i<genParticles->size();i++){
        const reco::Candidate *particle = &(*genParticles)[i];
        // Drop everything that is not photon or not status 1
        if( abs(particle->pdgId()) != 22 || particle->status() != 1 )
        continue;
//    std::cout<<"particle->pdgId()) = "<<particle->pdgId()<<std::endl;    
//        double dRtmp = ROOT::Math::VectorUtil::DeltaR( pho.p4(), particle->p4() );
//std::cout<<"pho->superCluster()->eta(),pho->superCluster()->phi(),particle->eta(),particle->phi()"<<pho.eta()<<" "<<pho.phi()<<" "<<particle->eta()<<" "<<particle->phi()<<std::endl;  
    double dRtmp = deltaR(pho.eta(),pho.phi(),particle->eta(),particle->phi());  
      if( dRtmp < dR ){
            dR = dRtmp;
            closestPhoton = particle;
        }
//std::cout<<"closestPhoton ="<<closestPhoton->pdgId()<<std::endl;    
  }
    // See if the closest photon (if it exists) is close enough.
    // If not, no match found.
    if( !(closestPhoton != 0 && dR < 0.3) ) {
        return UNMATCHED;
       // ISRPho = false;
    }
    // Find ID of the parent of the found generator level photon match
    int ancestorPID = -999;
    int ancestorStatus = -999;
    findFirstNonPhotonMother(closestPhoton, ancestorPID, ancestorStatus);
    // Allowed parens: quarks pdgId 1-5, or a gluon 21
    std::vector<int> allowedParents { -1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -21, 21,-11,11,-13,13,-23,23,-24,24 };
    if( !(std::find(allowedParents.begin(),
                    allowedParents.end(), ancestorPID)
          != allowedParents.end()) ){
        // So it is not from g, u, d, s, c, b. Check if it is from pi0 or not.
        if( abs(ancestorPID) == 111 )
        return MATCHED_FROM_PI0;
       // ISRPho =true;
        else
      std::cout<<"Mother = "<<abs(ancestorPID)<<std::endl;
        return MATCHED_FROM_OTHER_SOURCES;
      //  ISRPho =true;
    }
    return MATCHED_FROM_GUDSCB;
     //   ISRPho =true;
}



void PKUTreeMaker::findFirstNonPhotonMother(const reco::Candidate *particle,
                                                   int &ancestorPID, int &ancestorStatus){
    if( particle == 0 ){
        printf("SimplePhotonNtupler: ERROR! null candidate pointer, this should never happen\n");
        return;
    }
    // Is this the first non-photon parent? If yes, return, otherwise
    // go deeper into recursion
    if( abs(particle->pdgId()) == 22 ){
        findFirstNonPhotonMother(particle->mother(0), ancestorPID, ancestorStatus);
    }else{
        ancestorPID = particle->pdgId();
        ancestorStatus = particle->status();
    }
    return;
}

//----------------------------------------------------------------

PKUTreeMaker::~PKUTreeMaker()
{
    // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//-------------------------------------------------------------------------------------------------------------------------------------//
void
PKUTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   isTrue_.clear();
   setDummyValues(); //Initalize variables with dummy values
   nevent = iEvent.eventAuxiliary().event();
   run    = iEvent.eventAuxiliary().run();
   ls     = iEvent.eventAuxiliary().luminosityBlock();
//events weight
   if (RunOnMC_){
        edm::Handle<GenEventInfoProduct> genEvtInfo;
        iEvent.getByLabel( "generator", genEvtInfo );
        theWeight = genEvtInfo->weight();
        if(theWeight>0) nump = nump+1;
        if(theWeight<0) numm = numm+1;
   }


   Handle<TriggerResults> trigRes;
   iEvent.getByToken(hltToken_, trigRes);
   HLT_Ele = (int)trigRes->accept(hltConfig.triggerIndex(elPaths[0]));
   HLT_Mu  = (int)trigRes->accept(hltConfig.triggerIndex(muPaths[0]));

   edm::Handle<edm::View<reco::Candidate> > leptonicVs;
   iEvent.getByLabel(leptonicVSrc_.c_str(), leptonicVs);

   if (leptonicVs->empty()) {  outTree_->Fill(); return;  }


 
   iEvent.getByToken(rhoToken_      , rho_     );
   double fastJetRho = *(rho_.product());
   useless = fastJetRho;
  
   edm::Handle<edm::View<pat::Jet> > ak4jets;
   iEvent.getByLabel(ak4jetsSrc_.c_str(), ak4jets);

   edm::Handle<edm::View<pat::Photon> > photons;
   iEvent.getByLabel(photonSrc_.c_str(), photons);

   edm::Handle<edm::View<reco::GenParticle> > genParticles;//define genParticle
   iEvent.getByLabel(InputTag("prunedGenParticles"), genParticles);

   edm::Handle<edm::View<pat::Muon>> loosemus;
   iEvent.getByToken(loosemuonToken_,loosemus); 

   edm::Handle<edm::View<pat::Electron>> looseeles;
   iEvent.getByToken(looseelectronToken_,looseeles); 

   edm::Handle<edm::View<reco::Candidate> > metHandle; 
   iEvent.getByLabel(metSrc_.c_str(), metHandle); 



//filter
   iEvent.getByToken(noiseFilterToken_, noiseFilterBits_);
   const edm::TriggerNames &names = iEvent.triggerNames(*noiseFilterBits_);
   bool HcalNoiseFilter = false;
   for (unsigned int i = 0, n = noiseFilterBits_->size(); i < n; ++i) {
   if (names.triggerName(i) == HBHENoiseFilter_Selector_)
            HcalNoiseFilter = noiseFilterBits_->accept(i); // TO BE USED
   if (names.triggerName(i) == CSCHaloNoiseFilter_Selector_)
            passFilter_CSCHalo_ = noiseFilterBits_->accept(i); // TO BE USED
   if (names.triggerName(i) == HCALlaserNoiseFilter_Selector_)
            passFilter_HCALlaser_ = noiseFilterBits_->accept(i); // DEPRECATED
   if (names.triggerName(i) == ECALDeadCellNoiseFilter_Selector_)
            passFilter_ECALDeadCell_ = noiseFilterBits_->accept(i); // under scrutiny
   if (names.triggerName(i) == GoodVtxNoiseFilter_Selector_)
            passFilter_GoodVtx_ = noiseFilterBits_->accept(i); // TO BE USED
   if (names.triggerName(i) == TrkFailureNoiseFilter_Selector_)
            passFilter_TrkFailure_ = noiseFilterBits_->accept(i); // DEPRECATED
   if (names.triggerName(i) == EEBadScNoiseFilter_Selector_)
            passFilter_EEBadSc_ = noiseFilterBits_->accept(i); // under scrutiny
   if (names.triggerName(i) == ECALlaserNoiseFilter_Selector_)
            passFilter_ECALlaser_ = noiseFilterBits_->accept(i); // DEPRECATED
   if (names.triggerName(i) == TrkPOGNoiseFilter_Selector_)
            passFilter_TrkPOG_ = noiseFilterBits_->accept(i); // DEPRECATED
   if (names.triggerName(i) == TrkPOG_manystrip_NoiseFilter_Selector_)
            passFilter_TrkPOG_manystrip_ = noiseFilterBits_->accept(i); // DEPRECATED
   if (names.triggerName(i) == TrkPOG_toomanystrip_NoiseFilter_Selector_)
            passFilter_TrkPOG_toomanystrip_ = noiseFilterBits_->accept(i); // DEPRECATED
   if (names.triggerName(i) == TrkPOG_logError_NoiseFilter_Selector_)
            passFilter_TrkPOG_logError_ = noiseFilterBits_->accept(i); // DEPRECATED
   if (names.triggerName(i) == METFilters_Selector_)
            passFilter_METFilters_ = noiseFilterBits_->accept(i); // DEPRECATED
   }
   passFilter_HBHE_ = HcalNoiseFilter;

  
   const reco::Candidate& leptonicV = leptonicVs->at(0);
   const reco::Candidate& metCand = metHandle->at(0);
   const reco::Candidate& lepton = (*leptonicV.daughter(0));
       
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByLabel("offlineSlimmedPrimaryVertices", vertices);
   if (vertices->empty()) { outTree_->Fill(); return;} // skip the event if no PV found
   nVtx = vertices->size();
   reco::VertexCollection::const_iterator firstGoodVertex = vertices->end();
   for (reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx) {
     // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
    // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
       if (  /*!vtx->isFake() &&*/ 
           !(vtx->chi2()==0 && vtx->ndof()==0) 
           &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
           && fabs(vtx->position().Z())<=24.0) {
           firstGoodVertex = vtx;
           break;
          }           
      }
   if ( firstGoodVertex==vertices->end() ) {outTree_->Fill();  return;} // skip event if there are no good PVs


    //************************* MET **********************//
    //*****************************************************************//
    edm::Handle<pat::METCollection>  METs_;
    bool defaultMET = iEvent.getByToken(metInputToken_ , METs_ );
    if(defaultMET){
        addTypeICorr(iEvent);
        for (const pat::MET &met : *METs_) {
         const float  rawPt    = met.shiftedPt(pat::MET::METUncertainty::NoShift, pat::MET::METUncertaintyLevel::Raw);
         const float  rawPhi   = met.shiftedPhi(pat::MET::METUncertainty::NoShift, pat::MET::METUncertaintyLevel::Raw);
         const float  rawSumEt = met.shiftedSumEt(pat::MET::METUncertainty::NoShift, pat::MET::METUncertaintyLevel::Raw);
//              const float rawPt = met.uncorPt();
//              const float rawPhi = met.uncorPhi();
//              const float rawSumEt = met.uncorSumEt();


            TVector2 rawMET_;
            rawMET_.SetMagPhi (rawPt, rawPhi );
            Double_t rawPx = rawMET_.Px();
            Double_t rawPy = rawMET_.Py();
            Double_t rawEt = std::hypot(rawPx,rawPy);
            METraw_et = rawEt;
            METraw_phi = rawPhi;
            METraw_sumEt = rawSumEt;
            double pxcorr = rawPx+TypeICorrMap_["corrEx"];
            double pycorr = rawPy+TypeICorrMap_["corrEy"];
            double et     = std::hypot(pxcorr,pycorr);
            double sumEtcorr = rawSumEt+TypeICorrMap_["corrSumEt"];
            TLorentzVector corrmet; corrmet.SetPxPyPzE(pxcorr,pycorr,0.,et);
            useless = sumEtcorr;
            useless = rawEt;
            MET_et = et;
            MET_phi = corrmet.Phi();
            MET_sumEt = sumEtcorr;
            MET_corrPx = TypeICorrMap_["corrEx"];
            MET_corrPy = TypeICorrMap_["corrEy"];
        }
    }
//---------------------------
       /// For the time being, set these to 1
       triggerWeight=1.0;
       pileupWeight=1.0;
       double targetEvents = targetLumiInvPb_*crossSectionPb_;
       lumiWeight = targetEvents/originalNEvents_;

       lep          = std::max(abs(leptonicV.daughter(0)->pdgId()), abs(leptonicV.daughter(1)->pdgId()));
       ptVlep       = leptonicV.pt();
       yVlep        = leptonicV.eta();
       phiVlep      = leptonicV.phi();
       massVlep     = leptonicV.mass();
       mtVlep       = leptonicV.mt();
       ptlep1       = leptonicV.daughter(1)->pt();
       etalep1      = leptonicV.daughter(1)->eta();
       philep1      = leptonicV.daughter(1)->phi();
       double energylep1     = leptonicV.daughter(1)->energy();
       if(leptonicV.daughter(0)->isElectron()||leptonicV.daughter(0)->isMuon() ) {
       ptlep1       = leptonicV.daughter(0)->pt();
       etalep1      = leptonicV.daughter(0)->eta();
       philep1      = leptonicV.daughter(0)->phi(); 
       energylep1     = leptonicV.daughter(0)->energy(); }


       met          = metCand.pt();
       metPhi       = metCand.phi();

       nlooseeles = looseeles->size(); 
       nloosemus = loosemus->size(); 

       TLorentzVector  glepton;
       glepton.SetPtEtaPhiE(ptlep1, etalep1, philep1, energylep1);
       math::XYZTLorentzVector neutrinoP4 = getNeutrinoP4(MET_et, MET_phi, glepton, 1);
       reco::CandidateBaseRef METBaseRef = metHandle->refAt(0);  //?????
       reco::ShallowCloneCandidate neutrino(METBaseRef, 0 , neutrinoP4);
       reco::CompositeCandidate WLeptonic;
       WLeptonic.addDaughter(lepton);
       WLeptonic.addDaughter(neutrino); 
       AddFourMomenta addP4;
       addP4.set(WLeptonic);
       ptVlepJEC       = WLeptonic.pt();
       yVlepJEC        = WLeptonic.eta();
       phiVlepJEC      = WLeptonic.phi();
       massVlepJEC     = WLeptonic.mass();
       mtVlepJEC       = WLeptonic.mt();

     //******************************************************************//
     //************************* Photon Jets Information******************//
     //*******************************************************************//
         double rhoVal_;
         rhoVal_=-99.;
         rhoVal_ = *rho_;

         photonet=-100.;  iphoton=-1;
            for (size_t ip=0; ip<photons->size();ip++)
         {
            if(ip<6)  {
                photon_pt[ip] = (*photons)[ip].pt();
                photon_eta[ip] = (*photons)[ip].eta();
                photon_phi[ip] = (*photons)[ip].phi();
                photon_e[ip] = (*photons)[ip].energy();
            }

            int istightphoton=0;

             edm::Handle<edm::View<pat::Electron> > electrons;
             iEvent.getByToken(electronToken_, electrons);
             edm::Handle<reco::BeamSpot> beamSpot;
             iEvent.getByToken(beamSpotToken_,beamSpot);
             edm::Handle<std::vector<reco::Conversion> > conversions;
             iEvent.getByToken(conversionsToken_,conversions);
             
             passEleVeto = (!hasMatchedPromptElectron((*photons)[ip].superCluster(),electrons, conversions, beamSpot->position() ) );

            double phoiso=std::max((*photons)[ip].photonIso()-rhoVal_*EApho(fabs((*photons)[ip].eta())),0.0);
            double chiso=std::max((*photons)[ip].chargedHadronIso()-rhoVal_*EAch(fabs((*photons)[ip].eta())),0.0);
            double nhiso=std::max((*photons)[ip].neutralHadronIso()-rhoVal_*EAnh(fabs((*photons)[ip].eta())),0.0);

            if(passEleVeto && (*photons)[ip].isEB() && (*photons)[ip].hadTowOverEm()<0.050 && (*photons)[ip].sigmaIetaIeta()<0.01 && chiso<0.91 && nhiso<(0.33 + exp(0.0044*(*photons)[ip].pt()+0.5809)) && phoiso<(0.61+0.0043*(*photons)[ip].pt())) {istightphoton=1;}
            if(passEleVeto && (*photons)[ip].isEE() && (*photons)[ip].hadTowOverEm()<0.050 && (*photons)[ip].sigmaIetaIeta()<0.0267 && chiso<0.65 && nhiso<(0.93 + exp(0.0040*(*photons)[ip].pt()+0.9402)) && phoiso<(0.54+0.0041*(*photons)[ip].pt())) {istightphoton=1;}

             if(istightphoton==1 && deltaR(photon_eta[ip],photon_phi[ip],etalep1,philep1) > 0.5) {
                 if(photon_pt[ip]>photonet) 
                         {
                      iphoton=ip;
        
          }
         }
         }
             //Gen photon matching
    if(RunOnMC_ && iphoton>-1){
             const auto pho = photons->ptrAt(iphoton);
             isTrue_.push_back( matchToTruth(*pho, genParticles, ISRPho, dR_));
      //  std::cout << "dR " << dR_ << std::endl;
    }

         if(iphoton>-1) {
               photonet=(*photons)[iphoton].pt();
               photoneta=(*photons)[iphoton].eta();
               photonphi=(*photons)[iphoton].phi();
               photone=(*photons)[iphoton].energy();
               photonsieie=(*photons)[iphoton].sigmaIetaIeta();
               photonphoiso=std::max((*photons)[iphoton].photonIso()-rhoVal_*EApho(fabs((*photons)[iphoton].eta())),0.0);
               photonchiso=std::max((*photons)[iphoton].chargedHadronIso()-rhoVal_*EAch(fabs((*photons)[iphoton].eta())),0.0);
               photonnhiso=std::max((*photons)[iphoton].neutralHadronIso()-rhoVal_*EAnh(fabs((*photons)[iphoton].eta())),0.0);
               drla=deltaR(photon_eta[iphoton],photon_phi[iphoton],etalep1,philep1);
         }  
//std::cout<<iphoton<<" "<<photonet<<std::endl;
    
             //******************************************************************//
            //************************* AK4 Jets Information******************//
             //******************************************************************//
    Int_t jetindexphoton12[2] = {-1,-1}; 

    std::vector<JetCorrectorParameters> vPar;
    for ( std::vector<std::string>::const_iterator payloadBegin = jecAK4Labels_.begin(), payloadEnd = jecAK4Labels_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
         JetCorrectorParameters pars(*ipayload);
         vPar.push_back(pars);    }
    jecAK4_ = new FactorizedJetCorrector(vPar);
    vPar.clear();

    int nujets=0 ;
    double tmpjetptcut=20.0;
    std::vector<TLorentzVector*> jets;
    std::vector<TLorentzVector*> ak4jet_p4_jer;
   
//################Jet Correction##########################
        for (size_t ik=0; ik<ak4jets->size();ik++)
         {
            reco::Candidate::LorentzVector uncorrJet = (*ak4jets)[ik].correctedP4(0);
            jecAK4_->setJetEta( uncorrJet.eta() );
            jecAK4_->setJetPt ( uncorrJet.pt() );
            jecAK4_->setJetE ( uncorrJet.energy() );
            jecAK4_->setRho ( rhoVal_ );
            jecAK4_->setNPV ( vertices->size() );
            jecAK4_->setJetA ( (*ak4jets)[ik].jetArea() );
            double corr = jecAK4_->getCorrection();

            if(corr*uncorrJet.pt()>tmpjetptcut) {
            TLorentzVector *dummy = new TLorentzVector(0,0,0,0);    
            dummy->SetPtEtaPhiE(corr*uncorrJet.pt(), uncorrJet.eta(), uncorrJet.phi(), corr*uncorrJet.energy());
            jets.push_back(dummy);
            ++nujets;
            }   

            if(ik<6)  {   
                ak4jet_pt[ik] =  corr*uncorrJet.pt();
                ak4jet_eta[ik] = (*ak4jets)[ik].eta();
                ak4jet_phi[ik] = (*ak4jets)[ik].phi();
                ak4jet_e[ik] =   corr*uncorrJet.energy();
                ak4jet_csv[ik] = (*ak4jets)[ik].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
                ak4jet_icsv[ik] = (*ak4jets)[ik].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");   }
          }
    
    sort (jets.begin (), jets.end (), mysortPt);
           for (size_t i=0;i<jets.size();i++) {
             if(iphoton>0) {
               double drtmp1=deltaR(jets.at(i)->Eta(), jets.at(i)->Phi(), photoneta,photonphi);
              if(drtmp1>0.5 && jetindexphoton12[0]==-1&&jetindexphoton12[1]==-1) {
                     jetindexphoton12[0] = i;
                     continue;
              }
              if(drtmp1>0.5 && jetindexphoton12[0]!=-1&&jetindexphoton12[1]==-1) {
                     jetindexphoton12[1] = i;
                     continue;
              }
            }
         }


         if(jetindexphoton12[0]>-1 && jetindexphoton12[1]>-1) {
            jet1pt=jets[jetindexphoton12[0]]->Pt();
            jet1eta=jets[jetindexphoton12[0]]->Eta();
            jet1phi=jets[jetindexphoton12[0]]->Phi();
            jet1e=jets[jetindexphoton12[0]]->E();
            jet2pt=jets[jetindexphoton12[1]]->Pt();
            jet2eta=jets[jetindexphoton12[1]]->Eta();
            jet2phi=jets[jetindexphoton12[1]]->Phi();
            jet2e=jets[jetindexphoton12[1]]->E();

            drj1a=deltaR(jet1eta,jet1phi,photoneta,photonphi);
            drj2a=deltaR(jet2eta,jet2phi,photoneta,photonphi);
            drj1l=deltaR(jet1eta,jet1phi,etalep1,philep1);
            drj2l=deltaR(jet2eta,jet2phi,etalep1,philep1);
            
            TLorentzVector j1p4;
            j1p4.SetPtEtaPhiE(jet1pt, jet1eta, jet1phi, jet1e);
            TLorentzVector j2p4;
            j2p4.SetPtEtaPhiE(jet2pt, jet2eta, jet2phi, jet2e);
            TLorentzVector photonp4;
            photonp4.SetPtEtaPhiE(photonet, photoneta, photonphi, photone);
            TLorentzVector vp4;
            vp4.SetPtEtaPhiE(leptonicV.pt(), leptonicV.eta(), leptonicV.phi(), leptonicV.energy());

            j1metPhi=fabs(jet1phi-metPhi);
            if(j1metPhi>Pi) {j1metPhi=2.0*Pi-j1metPhi;}
            

            j2metPhi=fabs(jet2phi-metPhi);
            if(j2metPhi>Pi) {j2metPhi=2.0*Pi-j2metPhi;}

 
            Mjj=(j1p4 + j2p4).M();
            deltaeta = fabs(jet1eta - jet2eta);
            zepp = fabs((vp4+photonp4).Rapidity() - (j1p4.Rapidity() + j2p4.Rapidity())/ 2.0); 
         }
   
       outTree_->Fill();
   }
   

//-------------------------------------------------------------------------------------------------------------------------------------//


void PKUTreeMaker::setDummyValues() {
     nVtx           = -1e1;
     triggerWeight  = -1e1;
     pileupWeight   = -1e1;
     lumiWeight     = -1e1;
     theWeight = -99;
     lep            = -1e1;
     nlooseeles=-1e1;
     nloosemus=-1e1;
     ptVlep         = -1e1;
     yVlep          = -1e1;
     phiVlep        = -1e1;
     massVlep       = -1e1;
     mtVlep         = -1e1;
     ptVlepJEC         = -1e1;
     yVlepJEC          = -1e1;
     phiVlepJEC        = -1e1;
     massVlepJEC       = -1e1;
     mtVlepJEC         = -1e1;
     ptlep1         = -1e1;
     etalep1        = -1e1;
     philep1        = -1e1;
     met            = -1e1;
     metPhi         = -1e1;
     j1metPhi         = -1e1;
     j2metPhi         = -1e1;
     METraw_et = -99;
     METraw_phi = -99;
     METraw_sumEt = -99;
     MET_et = -99;
     MET_phi = -99;
     MET_sumEt = -99;
     MET_corrPx = -99;
     MET_corrPy = -99;
     photon_pt[0] = -1e1;
     photon_pt[1] = -1e1;
     photon_pt[2] = -1e1;
     photon_pt[3] = -1e1;
     photon_pt[4] = -1e1;
     photon_pt[5] = -1e1;
     photon_eta[0] = -1e1;
     photon_eta[1] = -1e1;
     photon_eta[2] = -1e1;
     photon_eta[3] = -1e1;
     photon_eta[4] = -1e1;
     photon_eta[5] = -1e1;
     photon_phi[0] = -1e1;
     photon_phi[1] = -1e1;
     photon_phi[2] = -1e1;
     photon_phi[3] = -1e1;
     photon_phi[4] = -1e1;
     photon_phi[5] = -1e1;
     photon_e[0] = -1e1;
     photon_e[1] = -1e1;
     photon_e[2] = -1e1;
     photon_e[3] = -1e1;
     photon_e[4] = -1e1;
     photon_e[5] = -1e1;

     photonet=-1e1;
     photoneta=-1e1;
     photonphi=-1e1;
     photone=-1e1;
     photonsieie=-1e1;
     photonphoiso=-1e1;
     photonchiso=-1e1;
     photonnhiso=-1e1;
     iphoton=-1;
     drla=1e1;
     passEleVeto=false;
    
     ISRPho = false;
     dR_ = 999;
 
     jet1pt=-1e1;
     jet1eta=-1e1;
     jet1phi=-1e1;
     jet1e=-1e1;
     jet1csv=-1e1;
     jet1icsv=-1e1;
     jet2pt=-1e1;
     jet2eta=-1e1;
     jet2phi=-1e1;
     jet2e=-1e1;
     jet2csv=-1e1;
     jet2icsv=-1e1;
     drj1a=1e1;
     drj2a=1e1;
     drj1l=1e1;
     drj2l=1e1;
     Mjj=-1e1;
     deltaeta=-1e1;
     zepp=-1e1;
  
     HLT_Ele=-99;
     HLT_Mu=-99;
 
     passFilter_HBHE_                  = false;
     passFilter_CSCHalo_               = false;
     passFilter_HCALlaser_             = false;
     passFilter_ECALDeadCell_          = false;
     passFilter_GoodVtx_               = false;
     passFilter_TrkFailure_            = false;
     passFilter_EEBadSc_               = false;
     passFilter_ECALlaser_             = false;
     passFilter_TrkPOG_                = false;
     passFilter_TrkPOG_manystrip_      = false;
     passFilter_TrkPOG_toomanystrip_   = false;
     passFilter_TrkPOG_logError_       = false;
     passFilter_METFilters_            = false; 
}

// ------------ method called once each job just before starting event loop  ------------
void 
PKUTreeMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void PKUTreeMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
 {
   bool changed;
   if ( !hltConfig.init(iRun, iSetup, "HLT", changed) ) {
        edm::LogError("HltAnalysis") << "Initialization of HLTConfigProvider failed!!";
       return;
      }

   for (size_t i = 0; i < elPaths_.size(); i++) {
         std::vector<std::string> foundPaths = hltConfig.matched( hltConfig.triggerNames(), elPaths_[i] );
         while ( !foundPaths.empty() ){
               elPaths.push_back( foundPaths.back() );
               foundPaths.pop_back();
                                      }
                                                }
   for (size_t i = 0; i < muPaths_.size(); i++) {
         std::vector<std::string> foundPaths = hltConfig.matched( hltConfig.triggerNames(), muPaths_[i] );
         while ( !foundPaths.empty() ){
               muPaths.push_back( foundPaths.back() );
               foundPaths.pop_back();
                                      }
                                                }
   std::cout<<"\n************** HLT Information **************\n";
   for (size_t i=0; i < elPaths.size(); i++) std::cout << "\n Electron paths: " << elPaths[i].c_str() <<"\t"<< std::endl;
   for (size_t i=0; i < muPaths.size(); i++) std::cout << "\n Muon paths    : " << muPaths[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

}

void PKUTreeMaker::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
PKUTreeMaker::endJob() {
  std::cout << "PKUTreeMaker endJob()..." << std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(PKUTreeMaker);
