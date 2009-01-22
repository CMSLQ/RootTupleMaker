// -*- C++ -*-
//
// Package:    RootTupleMaker
// Class:      RootTupleMaker
// 
/**\class RootTupleMaker RootTupleMaker.cc RootTuple/RootTupleMaker/src/RootTupleMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Ellie Lockner
//         Created:  Tue Oct 21 13:56:04 CEST 2008
// $Id: RootTupleMaker.cc,v 1.11 2008/12/01 17:59:28 lockner Exp $
//
//

#include "Leptoquarks/RootTupleMaker/interface/RootTupleMaker.h"

//namesapecs
using namespace std;
using namespace edm;
using namespace reco;

typedef std::pair<const PixelMatchGsfElectron*, int> my_pair;

bool ComparePtPair(const my_pair& left , const my_pair& right)
{
  return left.first->pt() > right.first->pt();
}

bool sortByET(const reco::Muon &x, const reco::Muon &y)
    {
      return x.et()>y.et();
    }


//
// class decleration
//

class RootTupleMaker : public edm::EDAnalyzer {
   public:
      explicit RootTupleMaker(const edm::ParameterSet&);
      ~RootTupleMaker();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  void CreateParticleTree(const edm::Handle<reco::GenParticleCollection> collection);

  void SetTriggers(const edm::Event& iEvent);

  bool Skim1st2ndGenLQ(const edm::View<reco::Candidate> *emObjectHandle, const reco::MuonCollection muonCollection);

  int singleEleRelHLTCounter;
  int muonHLTCounter;

      // ----------member data ---------------------------
 // read from cfg file
  std::string          rootfile_;
  int                  maxgenparticles_;
  int                  maxgenjets_;
  int                  maxelectrons_;
  int                  maxcalojets_;
  int                  maxmuons_;
  bool                 aodsim_;
  bool                 fastSim_;
  bool                 PAT_;
  bool                 debug_;
  double               luminosity_;
  int                  numEvents_;              
  bool                 saveTrigger_;
  int                  prescaleSingleEleRel_;
  int                  prescaleMuon_;
  bool                 useSkim1st2ndGenLQ_;
  double               skim1st2ndGenLQpTEle_;
  double               skim1st2ndGenLQpTMu_; 


 //Output RootNtuple
  TTree *              m_tree;
  TFile *              m_file;

  TTree *              m_tree2;

  //Global info
  int                  Nstart;

  // Event info
  int                  event;
  int                  runnum;

  // Gen Event Quantities  
  float                m_cross_section;
  float                m_auto_cross_section;
  int                  m_processID;
  float                m_filter_eff;
  float                m_pthat;
  float                m_weight;              

  // GenParticles
  Int_t                m_GenParticleCount;
  Float_t              m_GenParticleVX[MAXGENPARTICLES];              
  Float_t              m_GenParticleVY[MAXGENPARTICLES];              
  Float_t              m_GenParticleVZ[MAXGENPARTICLES];              
  Float_t              m_GenParticleP[MAXGENPARTICLES];              
  Float_t              m_GenParticlePt[MAXGENPARTICLES];              
  Float_t              m_GenParticlePx[MAXGENPARTICLES];              
  Float_t              m_GenParticlePy[MAXGENPARTICLES];              
  Float_t              m_GenParticlePz[MAXGENPARTICLES];              
  Float_t              m_GenParticleE[MAXGENPARTICLES];              
  Float_t              m_GenParticleEta[MAXGENPARTICLES];              
  Float_t              m_GenParticlePhi[MAXGENPARTICLES];              
  Int_t                m_GenParticlePdgId[MAXGENPARTICLES];              
  Int_t                m_GenParticleMotherIndex[MAXGENPARTICLES];
  Int_t                m_GenParticleNumDaught[MAXGENPARTICLES];  

  // Trigger
  TString              aNames[MAXHLTBITS];
  char                 aHLTNames[6000];
  Int_t                hltNamesLen;
  Int_t                hltCount;
  bool                 aHLTResults[MAXHLTBITS];

  // Electrons
  Int_t                eleCount;
  Float_t              eleEta[MAXELECTRONS];
  Float_t              elePhi[MAXELECTRONS];
  Float_t              elePt[MAXELECTRONS];
  Float_t              eleEnergy[MAXELECTRONS];
  Float_t              eleCaloEnergy[MAXELECTRONS];

  Float_t              eleHoE[MAXELECTRONS];
  Float_t              eleSigmaEE[MAXELECTRONS];
  Float_t              eleDeltaPhiTrkSC[MAXELECTRONS];
  Float_t              eleDeltaEtaTrkSC[MAXELECTRONS];

  Float_t              eleTrkIso[MAXELECTRONS];
  Float_t              eleNumTrkIso[MAXELECTRONS];
  Float_t              eleReducedEcalIso[MAXELECTRONS];
  Float_t              eleHcalRecHitIso[MAXELECTRONS];
  Float_t              eleEcalRecHitIso[MAXELECTRONS];
  Float_t              eleHcalTowerIso[MAXELECTRONS];
  Int_t                eleClassif[MAXELECTRONS];

  // GenJets
  Int_t                genJetCount;
  Float_t              genJetEta[MAXGENJETS];
  Float_t              genJetPhi[MAXGENJETS];
  Float_t              genJetPt[MAXGENJETS];
  Float_t              genJetEnergy[MAXGENJETS];
  Float_t              genJetEMF[MAXGENJETS];
  Float_t              genJetHADF[MAXGENJETS];

  // CaloJets
  Int_t                caloJetIC5Count;
  Float_t              caloJetIC5Eta[MAXCALOJETS];
  Float_t              caloJetIC5Phi[MAXCALOJETS];
  Float_t              caloJetIC5EMF[MAXCALOJETS];
  Float_t              caloJetIC5HADF[MAXCALOJETS];
  Float_t              caloJetIC5Pt_raw[MAXCALOJETS];
  Float_t              caloJetIC5Energy_raw[MAXCALOJETS];
  Float_t              caloJetIC5Pt_L23[MAXCALOJETS];
  Float_t              caloJetIC5Energy_L23[MAXCALOJETS];
  Float_t              caloJetIC5Pt[MAXCALOJETS];
  Float_t              caloJetIC5Energy[MAXCALOJETS];


  // Muons
  Int_t                muonCount;
  Int_t                muonGlobalCount;
  Int_t                muonStandAloneCount;
  Int_t                muonTrackerCount;
  Float_t              muonEta[MAXMUONS];
  Float_t              muonPhi[MAXMUONS];
  Float_t              muonPt[MAXMUONS];
  Float_t              muonEnergy[MAXMUONS];
  Int_t                muonCharge[MAXMUONS];
  Float_t              muonEt[MAXMUONS];
  Float_t              muonTrkHits[MAXMUONS];
  Float_t              muonTrkD0[MAXMUONS];
  Float_t              muonTrkDz[MAXMUONS];
  Float_t              muonEcalIso[MAXMUONS];
  Float_t              muonTrkIso[MAXMUONS];
  Float_t              muonHcalIso[MAXMUONS];
  Float_t              muonHOIso[MAXMUONS];
  Float_t              muonGlobalChi2[MAXMUONS];
  Float_t              muonStandAloneChi2[MAXMUONS];
  Float_t              muonTrackerChi2[MAXMUONS];

  // MET 
  Float_t              genMET;
  Float_t              MET;

  // Parameters from cfg
  std::string mUDSCorrectorName;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
RootTupleMaker::RootTupleMaker(const edm::ParameterSet& iConfig)

{
  //get parameters from cfg file
  rootfile_          = iConfig.getUntrackedParameter<std::string>("rootfile","LQRootNtuple.root");
  maxgenparticles_   = iConfig.getUntrackedParameter<int>("maxgenparticles",100); 
  maxgenjets_        = iConfig.getUntrackedParameter<int>("maxgenjets",10); 
  maxelectrons_      = iConfig.getUntrackedParameter<int>("maxelectrons",5); 
  maxcalojets_       = iConfig.getUntrackedParameter<int>("maxcalojets",10); 
  maxmuons_          = iConfig.getUntrackedParameter<int>("maxmuons",5); 

  aodsim_            = iConfig.getUntrackedParameter<bool>("aodsim",0);
  fastSim_           = iConfig.getUntrackedParameter<bool>("fastSim",0);
  PAT_                 = iConfig.getUntrackedParameter<bool>("PAT");
  debug_             = iConfig.getUntrackedParameter<bool>("debug",0);
  luminosity_        = iConfig.getUntrackedParameter<double>("luminosity",100); // pb -1
  numEvents_         = iConfig.getUntrackedParameter<int>("numEvents",100); 

  saveTrigger_           = iConfig.getUntrackedParameter<bool>("saveTrigger",1);
  prescaleSingleEleRel_  = iConfig.getUntrackedParameter<int>("prescaleSingleEleRel",30); 
  prescaleMuon_          = iConfig.getUntrackedParameter<int>("prescaleMuon",30);

  useSkim1st2ndGenLQ_       = iConfig.getUntrackedParameter<bool>("useSkim1st2ndGenLQ",0);
  skim1st2ndGenLQpTEle_  = iConfig.getUntrackedParameter<double>("skim1st2ndGenLQpTEle",30);
  skim1st2ndGenLQpTMu_  = iConfig.getUntrackedParameter<double>("skim1st2ndGenLQpTMu",30);


  //Initialize some variables
  singleEleRelHLTCounter=0;
  muonHLTCounter=0;

  event=-999;
  runnum=-999;
  Nstart=0;

  m_cross_section=-999.;
  m_auto_cross_section=-999.;
  m_processID=-999;
  m_filter_eff=-999.;
  m_pthat=-999.;
  m_weight=-999.;              

  genMET=-999.;
  MET=-999.;

}


RootTupleMaker::~RootTupleMaker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
RootTupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  if(debug_==true)
    cout << "Analyze " << endl;
  
  // Count number of events analyzed
  Nstart++;

  // Fill Event info

  event = iEvent.id().event();
  runnum = iEvent.id().run();

  // Fill gen info

  //using HepMCProduct

  // work in CSA08 RECO //////////////////////////////////
    Handle<HepMCProduct> mcHandle;
    iEvent.getByLabel("source", mcHandle );
    const HepMCProduct* mcCollection = mcHandle.failedToGet () ? 0 : &*mcHandle;
    
    int processID = -99;
    double pthat = -99;
    if (mcCollection) {
      const HepMC::GenEvent *genEvt = mcCollection->GetEvent();
      processID = genEvt->signal_process_id();
//       pthat = genEvt->event_scale(); 
    }

    
  edm::Handle<reco::PdfInfo> pdfHandle;
  iEvent.getByLabel("genEventPdfInfo",pdfHandle);
  const PdfInfo* pdfCollection = pdfHandle.failedToGet () ? 0 : &*pdfHandle;
  if (pdfCollection){
    pthat = pdfCollection->scalePDF;
  }

  m_processID = processID;
  m_pthat = pthat;

//     cout << processID << endl;
//     cout << pthat << endl;

  //using genEvent products

  // ## don't work in CSA08 RECO
  //   Handle<int> genProcessID; 
  //   iEvent.getByLabel( "genEventProcID", genProcessID ); 
  //   double processID = *genProcessID;
  //   cout << processID << endl; 
  // ##

  // ## work in CSA08 RECO
  //   Handle<double> genEventScale; 
  //   iEvent.getByLabel( "genEventScale", genEventScale );
  //   double pthat = *genEventScale; 
  //   cout << pthat << endl;
  // ##

  // ## don't work in CSA08 RECO
  // double filter_eff = -99.; 
  //   Handle<double> genFilterEff; 
  //   iEvent.getByLabel( "genEventRunInfo", "FilterEfficiency", genFilterEff); 
  //   filter_eff = *genFilterEff;
  //   cout << filter_eff << endl; 
  // ## 

  // ## don't work in CSA08 RECO
  //   double cross_section = -99.; 
  //   Handle<double> genCrossSect; 
  //   iEvent.getByLabel( "genEventRunInfo", "PreCalculatedCrossSection", genCrossSect); 
  //   cross_section = *genCrossSect;
  //   cout << cross_section << endl;
  // ## 

  if(debug_==true)
    cout << "gen event info filled" << endl;

  
  /////////// Electrons
  ///////////////////////////////////////////////////////////////////////////////////////
  // Get the objects that were fed into the isolation producer (not necessary for method 2)
  edm::Handle< edm::View<reco::Candidate> > emObjectHandle_;
  iEvent.getByLabel("pixelMatchGsfElectrons",emObjectHandle_);
  const edm::View<reco::Candidate> *emObjectHandle = emObjectHandle_.product();

  //now get the CaloTopology and rec hits for ID
  //note in practice you wouldnt hardcode the hit InputTags
  edm::ESHandle<CaloTopology> caloTopologyHandle;
  iSetup.get<CaloTopologyRecord>().get(caloTopologyHandle);
  edm::Handle<EcalRecHitCollection> ebReducedRecHitsHandle;
  iEvent.getByLabel("reducedEcalRecHitsEB",ebReducedRecHitsHandle);
  edm::Handle<EcalRecHitCollection> eeReducedRecHitsHandle;
  iEvent.getByLabel("reducedEcalRecHitsEE",eeReducedRecHitsHandle);
 
  const CaloTopology* caloTopology = caloTopologyHandle.product();
  const EcalRecHitCollection* ebRecHits = ebReducedRecHitsHandle.product();
  const EcalRecHitCollection* eeRecHits = eeReducedRecHitsHandle.product();

  //Isolation collections
  //trkiso ( EgammaElectronTkIsolationProducer )
  edm::Handle< reco::CandViewDoubleAssociations > trkIsolationHandle;
  iEvent.getByLabel("egammaElectronTkIsolation",trkIsolationHandle);
  const CandViewDoubleAssociations* trkIsolation = trkIsolationHandle.failedToGet () ? 0 : &*trkIsolationHandle;

  //numtrksio ( EgammaElectronTkNumIsolationProducer )
  edm::Handle< reco::CandViewDoubleAssociations > trkNumIsolationHandle;
  iEvent.getByLabel("egammaElectronTkNumIsolation",trkNumIsolationHandle);
  const CandViewDoubleAssociations* trkNumIsolation = trkNumIsolationHandle.failedToGet () ? 0 : &*trkNumIsolationHandle;

  //ecaliso ( EgammaEcalRecHitIsolationProducer )
  edm::Handle< reco::CandViewDoubleAssociations > ecalRecHitIsolationHandle;
  iEvent.getByLabel("egammaEcalRecHitIsolation",ecalRecHitIsolationHandle);
  const CandViewDoubleAssociations* ecalRecHitIsolation = ecalRecHitIsolationHandle.failedToGet () ? 0 : &*ecalRecHitIsolationHandle;

  //hcaliso ( EgammaHcalIsolationProducer )
  edm::Handle< reco::CandViewDoubleAssociations > hcalIsolationHandle;
  iEvent.getByLabel("egammaHcalIsolation",hcalIsolationHandle);
  const CandViewDoubleAssociations* hcalIsolation = hcalIsolationHandle.failedToGet () ? 0 : &*hcalIsolationHandle;

  //ecaliso ( EgammareducedEcalRecHitIsolationProducer ) for AOD
  edm::Handle< reco::CandViewDoubleAssociations > reducedEcalIsolationHandle;
  iEvent.getByLabel("reducedEcalRecHitIsolation",reducedEcalIsolationHandle);
  const CandViewDoubleAssociations* reducedEcalIsolation = reducedEcalIsolationHandle.failedToGet () ? 0 : &*reducedEcalIsolationHandle;

  //hcaliso ( EgammaTowerIsolationProducer ) for AOD
  edm::Handle< reco::CandViewDoubleAssociations > towerIsolationHandle;
  iEvent.getByLabel("egammaTowerIsolation",towerIsolationHandle);
  const CandViewDoubleAssociations* towerIsolation = towerIsolationHandle.failedToGet () ? 0 : &*towerIsolationHandle;


  /// sort electrons
  std::list<my_pair> electronRefListPair;
  for(int elecand_idx = 0; elecand_idx < (int)emObjectHandle->size(); elecand_idx++) 
    {
      const PixelMatchGsfElectronRef electron = emObjectHandle->refAt(elecand_idx).castTo<PixelMatchGsfElectronRef>();
      electronRefListPair.push_front( my_pair(&*electron,electron.key()) );

      //cout << "pT: " << electron->pt() << " " << "key: " << electron.key() << endl;

    }
  electronRefListPair.sort(ComparePtPair);  

  //Loop over electrons
  eleCount = 0;
  int eleidx=-1;
  std::list<my_pair>::const_iterator electron;
  for(electron=electronRefListPair.begin(); electron!=electronRefListPair.end(); electron++) 
    {
      //cout << (*electron).first->pt() << endl;

      eleidx++;

      const reco::SuperClusterRef& SCref = (*electron).first->superCluster();  

      //## Remove electrons associated to the same SC ##
      bool IsCopy=false;      

      int eleidx_1=-1;
      std::list<my_pair>::const_iterator electron_1;
      for(electron_1=electronRefListPair.begin(); electron_1!=electronRefListPair.end(); electron_1++) 
	{
	  eleidx_1++;	  
	  if(eleidx_1<=eleidx)
	    continue;
	  
	  const reco::SuperClusterRef& SCref_1 = (*electron_1).first->superCluster();  	  
	  
  	  if(SCref==SCref_1)
  	    {
  	      IsCopy=true;
  	      break;
  	    }
  	}
      //skip this electron 
      if(IsCopy==true)
	continue;
      //## end remove electrons associated to the same SC

      //## Start counting electrons from here ##
      
      //exit from loop when you reach the required number of electrons
      if(eleCount > maxelectrons_)
	break;

      //now we need to get the basic cluster and decide if it is barrel or endcap
      const reco::BasicCluster& seedClus = *(SCref->seed());
      const DetId firstDetId = seedClus.getHitsByDetId()[0];// this is NOT the seed but all hits will be either endcap or barrel

      float sigmaee = -999;

      if(firstDetId.subdetId()==EcalBarrel){
	std::vector<float> localCov = EcalClusterTools::localCovariances(seedClus,ebRecHits,caloTopology);
	sigmaee =  sqrt(localCov[0]);
      }else if(firstDetId.subdetId()==EcalEndcap){
	//identical to sigmaEtaEta which would be EcalClusterTools::covariances(seedClus,ebRecHits,caloGeometry,caloTopology)
	std::vector<float> localCov = EcalClusterTools::localCovariances(seedClus,eeRecHits,caloTopology);
	sigmaee =  sqrt(localCov[0]);
      }

      //////////ID variables
      float hOverE = (*electron).first->hadronicOverEm();
      float deltaPhiIn = (*electron).first->deltaPhiSuperClusterTrackAtVtx();
      float deltaEtaIn = (*electron).first->deltaEtaSuperClusterTrackAtVtx();
      int classif =  (*electron).first->classification();
 
      //////////Iso variables

      double trkIso = -99;
      double trkNumIso = -99;
      double ecalRecHitIso = -99;
      double hcalIso = -99;
      double reducedEcalRecHitIso = -99;
      double towerIso = -99;

      // this retrieves the index in the original collection associated to the reference to electron
      int index = (*electron).second;

      if (trkIsolation) trkIso =(*trkIsolation)[index].second;
      if (trkNumIsolation) trkNumIso = (*trkNumIsolation)[index].second; 
      if (ecalRecHitIsolation) ecalRecHitIso =(*ecalRecHitIsolation)[index].second;
      if (hcalIsolation) hcalIso =(*hcalIsolation)[index].second;
      if (reducedEcalIsolation) reducedEcalRecHitIso =(*reducedEcalIsolation)[index].second;
      if (towerIsolation) towerIso =(*towerIsolation)[index].second;

     // Set variables in RootNtuple
      eleEta[eleCount]=(*electron).first->eta();
      elePhi[eleCount]=(*electron).first->phi();
      elePt[eleCount]=(*electron).first->pt();
      eleEnergy[eleCount]=(*electron).first->energy();
      eleCaloEnergy[eleCount]=(*electron).first->caloEnergy();

      eleHoE[eleCount]=hOverE;
      eleSigmaEE[eleCount]=sigmaee;
      eleDeltaPhiTrkSC[eleCount]=deltaPhiIn;
      eleDeltaEtaTrkSC[eleCount]=deltaEtaIn;
      eleClassif[eleCount]=classif;

      eleTrkIso[eleCount]=trkIso;
      eleNumTrkIso[eleCount]=trkNumIso;
      eleEcalRecHitIso[eleCount]=ecalRecHitIso;
      eleHcalRecHitIso[eleCount]=hcalIso;
      eleReducedEcalIso[eleCount]=reducedEcalRecHitIso;
      eleHcalTowerIso[eleCount]=towerIso;

      //go to next electron
      eleCount++;
 
    }


  //////////// Gen Particles
  //////////////////////////////////////////////////////////////////////////////////////////////
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel ("genParticles", genParticles);
  CreateParticleTree( genParticles );     
  if(debug_==true)
    cout << "gen particles filled" << endl;

//    reco::GenParticleCollection::const_iterator genParts_it;
//   for (genParts_it = genParticles->begin();
//   	genParts_it != genParticles->end(); ++genParts_it)
//        {
// 	 if (abs(genParts_it->pdgId()) == 13) cout << "Found muon" << endl;
//        }

  //////////// GenJets
  ////////////////////////////////////////////////////////////////////////////////////////
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByLabel ("iterativeCone5GenJets", genJets);
  
  genJetCount=0;
  for( GenJetCollection::const_iterator genjet = genJets->begin(); genjet != genJets->end(); genjet++ ) 
    {

      //exit from loop when you reach the required number of electrons
      if(genJetCount > maxgenjets_)
	break;

      float EMF = genjet->emEnergy() / genjet->energy();
      float HADF = genjet->hadEnergy() / genjet->energy();

      genJetPt[genJetCount]=genjet->pt();
      genJetPhi[genJetCount]=genjet->phi();
      genJetEta[genJetCount]=genjet->eta();
      genJetEnergy[genJetCount]=genjet->energy();
      genJetEMF[genJetCount]=EMF;
      genJetHADF[genJetCount]=HADF;
      
      genJetCount++;
    }

  if(debug_==true)
    cout << "genJets filled" << endl;

  ////////////// CaloJets
  ////////////////////////////////////////////////////////////////////////////////////////////

  if (aodsim_ || fastSim_){
    //Jets need L2, L3, L5 corrections

  //raw
  edm::Handle<reco::CaloJetCollection> caloJetsIC5_raw;
  iEvent.getByLabel ("iterativeCone5CaloJets", caloJetsIC5_raw); 

  //L2+L3 correction
  edm::Handle<reco::CaloJetCollection> caloJetsIC5_L23;
  iEvent.getByLabel ("L2L3CorJet", caloJetsIC5_L23); 

  //L2+L3+L5 correction
  edm::Handle<reco::CaloJetCollection> caloJetsIC5;
  iEvent.getByLabel ("L2L3L5CorJet", caloJetsIC5); 

//   // get Jet Flavor corrector (L5)
//   const JetCorrector* udsJetCorrector = JetCorrector::getJetCorrector (mUDSCorrectorName, iSetup);

  caloJetIC5Count=0;
  for( CaloJetCollection::const_iterator calojet = caloJetsIC5->begin(); calojet != caloJetsIC5->end(); calojet++ ) 
    {
      //exit from loop when you reach the required number of electrons
      if(caloJetIC5Count > maxcalojets_)
	break;
     float EMF = calojet->emEnergyFraction();
      float HADF = calojet->energyFractionHadronic();

      caloJetIC5Pt[caloJetIC5Count]=calojet->pt();
      caloJetIC5Energy[caloJetIC5Count]=calojet->energy();
      caloJetIC5Pt[caloJetIC5Count]=calojet->pt();
      caloJetIC5Energy[caloJetIC5Count]=calojet->energy();
      caloJetIC5Pt_raw[caloJetIC5Count]=(*caloJetsIC5_raw)[caloJetIC5Count].pt();
      caloJetIC5Energy_raw[caloJetIC5Count]=(*caloJetsIC5_raw)[caloJetIC5Count].energy();
      caloJetIC5Pt_L23[caloJetIC5Count]=(*caloJetsIC5_L23)[caloJetIC5Count].pt();
      caloJetIC5Energy_L23[caloJetIC5Count]=(*caloJetsIC5_L23)[caloJetIC5Count].energy();
      caloJetIC5Phi[caloJetIC5Count]=calojet->phi();
      caloJetIC5Eta[caloJetIC5Count]=calojet->eta();
      caloJetIC5EMF[caloJetIC5Count]=EMF;
      caloJetIC5HADF[caloJetIC5Count]=HADF;
      
      caloJetIC5Count++;
    }

  } // end of if aod or fastSim

  if (PAT_){
    //probably have already L2,3 corrections, so need to grab correction factor
}

  if(debug_==true)
    cout << "CaloJets filled" << endl;

  ////////// MET and GenMET
  //////////////////////////////////////////////////////////////////////////////////////////

  // GENMET
  Handle<GenMETCollection> genMETColl;
  iEvent.getByLabel( "genMet" , genMETColl );

  const GenMET & genmet = (*genMETColl)[0];
  genMET = genmet.et();

  // MET
  Handle<CaloMETCollection> recoMETColl;
  iEvent.getByLabel ( "met" , recoMETColl);

  const reco::CaloMET & recomet = (*recoMETColl)[0];
  MET = recomet.et();

  if(debug_==true)
    cout << "MET filled" << endl;



  /////////////////// HLT info
  ///////////////////////////////////////////////////////////////////////////////////////////////

  if(saveTrigger_==true)  
    {
      SetTriggers(iEvent);
      if(debug_==true)
	cout << "HLT bits filled " << endl;
    }
  else 
    {
      // reset variables in ntuple
      for(unsigned int iHLT=0; iHLT<MAXHLTBITS; ++iHLT) {
	aHLTResults[iHLT] = false;
      }

      hltCount=-999;
      hltNamesLen=-999;
      strcpy(aHLTNames,"");

      if(debug_==true)
	cout << "HLT bits not filled" << endl;
    }

  ///////////////// Muons
  /////////////////////////////////////////////////////////////////////////////////////////////
  edm::Handle<reco::MuonCollection> MuonObjectHandle_;
  std::string m_muonLabel  = "muons";
  iEvent.getByLabel(m_muonLabel,MuonObjectHandle_);

  const reco::MuonCollection* MuonObjectHandle_tmp;
  MuonObjectHandle_tmp = MuonObjectHandle_.product();

  /// sort muons in decending ET
  reco::MuonCollection muonColl = *MuonObjectHandle_tmp;
  std::sort(muonColl.begin(), muonColl.end(), sortByET);

  muonCount = 0;  
  for(reco::MuonCollection::const_iterator muon = muonColl.begin(); muon != muonColl.end(); muon++)
    {
      //exit from loop when you reach the required number of muons
      if(muonCount > maxmuons_)
	break;

      unsigned int trkhits = 0;
      float trkd0 = 0.;
      float trkdz = 0.;
      float muonGlobalChi2=0.;
      float muonStandAloneChi2 =0.;
      float muonTrackerChi2 =0.;

      /*
      if(fastSim_==1)
	{
	  //FastSim
// 	  trkhits  = muon->track()->numberOfValidHits();
// 	  trkd0    = muon->track()->d0();
// 	  trkdz    = muon->track()->dz();
	}
      */

      if(muon->isGlobalMuon()){
	
	trkhits  = muon->globalTrack()->numberOfValidHits();
	trkd0    = muon->globalTrack()->d0();
	trkdz    = muon->globalTrack()->dz();
	muonGlobalChi2 = muon->globalTrack()->normalizedChi2();
	muonGlobalCount++;
      }
      
      if(muon->isStandAloneMuon()){
	trkhits  = muon->outerTrack()->numberOfValidHits();
	trkd0    = muon->outerTrack()->d0();
	trkdz    = muon->outerTrack()->dz();
	muonStandAloneChi2 = muon->outerTrack()->normalizedChi2();
	muonStandAloneCount++;
      }

      if (muon->isTrackerMuon()) {
	trkhits  = muon->innerTrack()->numberOfValidHits();
	trkd0    = muon->innerTrack()->d0();
	trkdz    = muon->innerTrack()->dz();
	muonTrackerChi2 = muon->innerTrack()->normalizedChi2();
	muonTrackerCount++;
      }

      float ptInCone       = 0.;
      float coneEMenergy   = 0.;
      float coneHADenergy  = 0.;
      float coneHOenergy   = 0.;

      bool m_useTrackConeSize03 = true;
      bool m_useTrackConeSize05 = false;

      if(muon->isGlobalMuon()){
	if (m_useTrackConeSize03)
	  {
	    ptInCone      = muon->isolationR03().sumPt;
	    coneEMenergy  = muon->isolationR03().emEt;
	    coneHADenergy = muon->isolationR03().hadEt;
	    coneHOenergy  = muon->isolationR03().hoEt;
	  }
	else if (m_useTrackConeSize05)
	  {
	    ptInCone      = muon->isolationR05().sumPt;
	    coneEMenergy  = muon->isolationR03().emEt;
	    coneHADenergy = muon->isolationR03().hadEt;
	    coneHOenergy  = muon->isolationR03().hoEt;
	  }
	
	muonEta[muonCount] = muon->eta();
	muonPhi[muonCount] = muon->phi();
	muonPt[muonCount]= muon->pt();
	muonEnergy[muonCount] = muon->energy();
	muonCharge[muonCount] = muon->charge();
	muonEt[muonCount] = muon->et();
	muonTrkHits[muonCount] = trkhits;
	muonTrkD0[muonCount] = trkd0;
	muonTrkDz[muonCount] = trkdz;
	muonEcalIso[muonCount] = coneEMenergy;
	muonTrkIso[muonCount] = ptInCone;
	muonHcalIso[muonCount] = coneHADenergy;
	muonHOIso[muonCount] = coneHOenergy;
	
	muonCount++;
      }
    } 

  //////////////////////////////////////////////////////////////////////////////////////////////

  //  Fill tree for each event
  //  ********************************************************
  //  ********************************************************
  if(debug_==true)
    cout << "About to fill tree" << endl;
  
  if( useSkim1st2ndGenLQ_ )
    {//fill only events passing skim
      if( Skim1st2ndGenLQ(emObjectHandle,muonColl) )
	m_tree->Fill();
    }
  else if( !useSkim1st2ndGenLQ_)
    {//fill all
      m_tree->Fill();
    }
  
  //  ********************************************************
  //  ********************************************************
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
RootTupleMaker::beginJob(const edm::EventSetup&)
{
  m_file = new TFile(rootfile_.c_str(),"RECREATE");
  m_tree = NULL;
  m_tree = new TTree ("RootTupleMaker","RootTupleMaker") ;

  m_tree2 = NULL;
  m_tree2 = new TTree ("RootTupleMaker_globalInfo","RootTupleMaker_globalInfo") ;

  m_tree2->Branch("Nstart",&Nstart,"Nstart/I");

  m_tree->Branch("event",&event,"event/I");
  m_tree->Branch("run",&runnum,"runnum/I");

   m_tree->Branch("processID",&m_processID,"processID/I");
   m_tree->Branch("pthat",&m_pthat,"pthat/F");
//   m_tree->Branch("cross_section",&m_cross_section,"cross_section/F");
//   m_tree->Branch("auto_cross_section",&m_auto_cross_section,"auto_cross_section/F");
//   m_tree->Branch("filter_eff",&m_filter_eff,"filter_eff/F"); 
//   m_tree->Branch("weight",&m_weight,"weight/F");
  
  m_tree->Branch("GenParticleCount",&m_GenParticleCount,"GenParticleCount/I");
  m_tree->Branch("GenParticleE",&m_GenParticleE,"GenParticleE[GenParticleCount]/F");
  m_tree->Branch("GenParticleP",&m_GenParticleP,"GenParticleP[GenParticleCount]/F");
  m_tree->Branch("GenParticlePt",&m_GenParticlePt,"GenParticlePt[GenParticleCount]/F");
  m_tree->Branch("GenParticlePx",&m_GenParticlePx,"GenParticlePz[GenParticleCount]/F");
  m_tree->Branch("GenParticlePy",&m_GenParticlePy,"GenParticlePy[GenParticleCount]/F");
  m_tree->Branch("GenParticlePz",&m_GenParticlePz,"GenParticlePz[GenParticleCount]/F");
  m_tree->Branch("GenParticlePdgId",&m_GenParticlePdgId,"GenParticlePdgId[GenParticleCount]/I");
  m_tree->Branch("GenParticleEta",&m_GenParticleEta,"GenParticleEta[GenParticleCount]/F");
  m_tree->Branch("GenParticlePhi",&m_GenParticlePhi,"GenParticlePhi[GenParticleCount]/F");
  m_tree->Branch("GenParticleVX",&m_GenParticleVX,"GenParticleVX[GenParticleCount]/F");
  m_tree->Branch("GenParticleVY",&m_GenParticleVY,"GenParticleVY[GenParticleCount]/F");
  m_tree->Branch("GenParticleVZ",&m_GenParticleVZ,"GenParticleVZ[GenParticleCount]/F");
  m_tree->Branch("GenParticleMotherIndex",&m_GenParticleMotherIndex,"GenParticleMotherIndex[GenParticleCount]/I");
  m_tree->Branch("GenParticleNumDaught",&m_GenParticleNumDaught,"GenParticleNumDaught[GenParticleCount]/I");

  m_tree->Branch("hltCount",&hltCount,"hltCount/I");
  m_tree->Branch("hltNamesLen",&hltNamesLen,"hltNamesLen/I");
  m_tree->Branch("HLTNames",&aHLTNames,"HLTNames[hltNamesLen]/C",6000);
  m_tree->Branch("HLTResults",&aHLTResults,"HLTResults[hltCount]/O");

  m_tree->Branch("eleCount",&eleCount,"eleCount/I");
  m_tree->Branch("eleEta",&eleEta,"eleEta[eleCount]/F");
  m_tree->Branch("elePhi",&elePhi,"elePhi[eleCount]/F");
  m_tree->Branch("elePt",&elePt,"elePt[eleCount]/F");
  m_tree->Branch("eleEnergy",&eleEnergy,"eleEnergy[eleCount]/F");
  m_tree->Branch("eleCaloEnergy",&eleCaloEnergy,"eleCaloEnergy[eleCount]/F");

  m_tree->Branch("eleHoE",&eleHoE,"eleHoE[eleCount]/F");
  m_tree->Branch("eleSigmaEE",&eleSigmaEE,"eleSigmaEE[eleCount]/F");
  m_tree->Branch("eleDeltaPhiTrkSC",&eleDeltaPhiTrkSC,"eleDeltaPhiTrkSC[eleCount]/F");
  m_tree->Branch("eleDeltaEtaTrkSC",&eleDeltaEtaTrkSC,"eleDeltaEtaTrkSC[eleCount]/F");

  m_tree->Branch("eleTrkIso",&eleTrkIso,"eleTrkIso[eleCount]/F");
  m_tree->Branch("eleNumTrkIso",&eleNumTrkIso,"eleNumTrkIso[eleCount]/F");
  m_tree->Branch("eleReducedEcalIso",&eleReducedEcalIso,"eleReducedEcalIso[eleCount]/F");
  m_tree->Branch("eleHcalTowerIso",&eleHcalTowerIso,"eleHcalTowerIso[eleCount]/F");
  m_tree->Branch("eleEcalRecHitIso",&eleEcalRecHitIso,"eleEcalRecHitIso[eleCount]/F");
  m_tree->Branch("eleHcalRecHitIso",&eleHcalRecHitIso,"eleHcalRecHitIso[eleCount]/F");
  m_tree->Branch("eleClassif",&eleClassif,"eleClassif[eleCount]/I");

  m_tree->Branch("genJetCount",&genJetCount,"genJetCount/I");
  m_tree->Branch("genJetEta",&genJetEta,"genJetEta[genJetCount]/F");
  m_tree->Branch("genJetPhi",&genJetPhi,"genJetPhi[genJetCount]/F");
  m_tree->Branch("genJetPt",&genJetPt,"genJetPt[genJetCount]/F");
  m_tree->Branch("genJetEnergy",&genJetEnergy,"genJetEnergy[genJetCount]/F");
  m_tree->Branch("genJetEMF",&genJetEMF,"genJetEMF[genJetCount]/F");
  m_tree->Branch("genJetHADF",&genJetHADF,"genJetHADF[genJetCount]/F");

  m_tree->Branch("caloJetIC5Count",&caloJetIC5Count,"caloJetIC5Count/I");
  m_tree->Branch("caloJetIC5Eta",&caloJetIC5Eta,"caloJetIC5Eta[caloJetIC5Count]/F");
  m_tree->Branch("caloJetIC5Phi",&caloJetIC5Phi,"caloJetIC5Phi[caloJetIC5Count]/F");
  m_tree->Branch("caloJetIC5Pt",&caloJetIC5Pt,"caloJetIC5Pt[caloJetIC5Count]/F");
  m_tree->Branch("caloJetIC5Energy",&caloJetIC5Energy,"caloJetIC5Energy[caloJetIC5Count]/F");
  m_tree->Branch("caloJetIC5Pt_raw",&caloJetIC5Pt_raw,"caloJetIC5Pt_raw[caloJetIC5Count]/F");
  m_tree->Branch("caloJetIC5Energy_raw",&caloJetIC5Energy_raw,"caloJetIC5Energy_raw[caloJetIC5Count]/F");
  m_tree->Branch("caloJetIC5Pt_L23",&caloJetIC5Pt_L23,"caloJetIC5Pt_L23[caloJetIC5Count]/F");
  m_tree->Branch("caloJetIC5Energy_L23",&caloJetIC5Energy_L23,"caloJetIC5Energy_L23[caloJetIC5Count]/F");
  m_tree->Branch("caloJetIC5EMF",&caloJetIC5EMF,"caloJetIC5EMF[caloJetIC5Count]/F");
  m_tree->Branch("caloJetIC5HADF",&caloJetIC5HADF,"caloJetIC5HADF[caloJetIC5Count]/F");

  m_tree->Branch("muonCount",&muonCount,"muonCount/I");
  m_tree->Branch("muonStandAloneCount",&muonStandAloneCount,"muonStandAloneCount/I");
  m_tree->Branch("muonGlobalCount",&muonGlobalCount,"muonGlobalCount/I");
  m_tree->Branch("muonTrackerCount",&muonTrackerCount,"muonTrackerCount/I");
  m_tree->Branch("muonEta",&muonEta,"muonEta[muonCount]/F");
  m_tree->Branch("muonPhi",&muonPhi,"muonPhi[muonCount]/F");
  m_tree->Branch("muonPt",&muonPt,"muonPt[muonCount]/F");
  m_tree->Branch("muonEnergy",&muonEnergy,"muonEnergy[muonCount]/F");
  m_tree->Branch("muonCharge",&muonCharge,"muonCharge[muonCount]/I");
  m_tree->Branch("muonEt",&muonEt,"muonEt[muonCount]/F");
  m_tree->Branch("muonTrkHits",&muonTrkHits,"muonTrkHits[muonCount]/F");
  m_tree->Branch("muonTrkD0",&muonTrkD0,"muonTrkD0[muonCount]/F");
  m_tree->Branch("muonTrkDz",&muonTrkDz,"muonTrkDz[muonCount]/F");
  m_tree->Branch("muonEcalIso",&muonEcalIso,"muonEcalIso[muonCount]/F");
  m_tree->Branch("muonTrkIso",&muonTrkIso,"muonTrkIso[muonCount]/F");
  m_tree->Branch("muonHcalIso",&muonHcalIso,"muonHcalIso[muonCount]/F");
  m_tree->Branch("muonHOIso",&muonHOIso,"muonHOIso[muonCount]/F");
  m_tree->Branch("muonGlobalChi2",&muonGlobalChi2,"muonGlobalChi2[muonCount]/F");
  m_tree->Branch("muonStandAloneChi2",&muonStandAloneChi2,"muonStandAloneChi2[muonCount]/F");
  m_tree->Branch("muonTrackerChi2",&muonTrackerChi2,"muonTrackerChi2[muonCount]/F");

  m_tree->Branch("genMET",&genMET,"genMET/F");
  m_tree->Branch("MET",&MET,"MET/F");

}

////====================================================================
void RootTupleMaker::CreateParticleTree(edm::Handle<reco::GenParticleCollection> collection)
{
////====================================================================
  m_GenParticleCount = 0;
  reco::GenParticleCollection::const_iterator cand;
  for(cand=collection->begin(); cand!=collection->end(); cand++) 
  {

    if(m_GenParticleCount  > maxgenparticles_)
        break;

    m_GenParticleP[m_GenParticleCount] = cand->p();
    m_GenParticlePx[m_GenParticleCount] = cand->px();
    m_GenParticlePy[m_GenParticleCount] = cand->py();
    m_GenParticlePz[m_GenParticleCount] = cand->pz();
    m_GenParticlePt[m_GenParticleCount] = cand->pt();
    m_GenParticleEta[m_GenParticleCount] = cand->eta();
    m_GenParticlePhi[m_GenParticleCount] = cand->phi();
    m_GenParticleE[m_GenParticleCount] = cand->energy();
    m_GenParticlePdgId[m_GenParticleCount] = cand->pdgId();
    m_GenParticleNumDaught[m_GenParticleCount] = cand->numberOfDaughters();
   
    m_GenParticleVX[m_GenParticleCount] = cand->vx();
    m_GenParticleVY[m_GenParticleCount] = cand->vy();
    m_GenParticleVZ[m_GenParticleCount] = cand->vz();
 
    m_GenParticleMotherIndex[m_GenParticleCount] = -1;
    int idx=0;
    reco::GenParticleCollection::const_iterator candIter;
    for(candIter = collection->begin(); candIter!=collection->end(); candIter++) {
      if(&(*candIter)==cand->mother()) {
       m_GenParticleMotherIndex[m_GenParticleCount] = idx;
       break;
      }
      idx++;
    }
  m_GenParticleCount++;
  }    
}


////====================================================================
void RootTupleMaker::SetTriggers(const edm::Event& iEvent)
{
  
  ////====================================================================
  // reset variables in ntuple
  for(unsigned int iHLT=0; iHLT<MAXHLTBITS; ++iHLT) {
    aHLTResults[iHLT] = false;
  }

  strcpy(aHLTNames,"");
  hltNamesLen = 0;

  edm::InputTag tag("TriggerResults::HLT");
  edm::Handle<edm::TriggerResults> hltTriggerResults;
  iEvent.getByLabel(tag,hltTriggerResults);


  if(!hltTriggerResults.isValid()) {
    std::cout << "invalid handle for HLT TriggerResults" << std::endl;
  } else {
    
    edm::TriggerNames HLTNames;
    HLTNames.init(*hltTriggerResults);
    
    std::string tempnames;
    
    hltCount = hltTriggerResults->size();

    for(int i = 0 ; i < hltCount ; i++) {
      //cout << "HLTTrigger: " << HLTNames.triggerName(i).c_str() << "  " << hltTriggerResults->accept(i) << endl;
      tempnames += HLTNames.triggerName(i) + ":";
      aHLTResults[i] = hltTriggerResults->accept(i);
    }
    
    hltNamesLen = tempnames.length();
    strcpy(aHLTNames,tempnames.c_str());
  
  }


}

////====================================================================
bool RootTupleMaker::Skim1st2ndGenLQ(const edm::View<reco::Candidate> *eleCollection, const reco::MuonCollection muonCollection)
{
////====================================================================

//   cout << "skim1st2ndGenLQpTMu_ : " << skim1st2ndGenLQpTMu_  << endl;
//   cout << "skim1st2ndGenLQpTEle_ : " << skim1st2ndGenLQpTEle_  << endl;

  int ne = 0;
  int nmu = 0;
  int nl = 0;

  reco::MuonCollection::const_iterator muon;
  for(muon=muonCollection.begin(); muon!=muonCollection.end(); muon++) 
    {
      //skip if muon is not global
      if(!muon->isGlobalMuon())
	continue;

      //      cout << "muon->pt(): " << muon->pt() << endl;

      if( muon->pt() > skim1st2ndGenLQpTMu_ 
	  || ( muon->pt() > skim1st2ndGenLQpTMu_ && nmu>=1 ) )
	nmu++;
    }
  
  for(int elecand_idx = 0; elecand_idx < (int)eleCollection->size(); elecand_idx++) 
    {
      const PixelMatchGsfElectronRef electron = eleCollection->refAt(elecand_idx).castTo<PixelMatchGsfElectronRef>();
      //      cout << "electron->pt(): " << electron->pt() << endl;
      const reco::SuperClusterRef& SCref = electron->superCluster();  	        

      //## Remove electrons associated to the same SC ##
      bool IsCopy=false;      
      
      for(int elecand_idx1 = 0; elecand_idx1 < (int)eleCollection->size(); elecand_idx1++) 
	{
	  const PixelMatchGsfElectronRef electron1 = eleCollection->refAt(elecand_idx1).castTo<PixelMatchGsfElectronRef>();

	  if(elecand_idx1<=elecand_idx)
	    continue;
	  
	  const reco::SuperClusterRef& SCref1 = electron1->superCluster();  	  
	  
  	  if( SCref == SCref1 )
  	    {
  	      IsCopy=true;
  	      break;
  	    }
  	}

      //skip this electron 
      if(IsCopy==true)
	continue;

      //## end remove electrons associated to the same SC
      
      if( electron->pt() > skim1st2ndGenLQpTEle_ 
	  || ( electron->pt() > skim1st2ndGenLQpTEle_ && ne>=1 ) )
	ne++;
    }

  nl = ne + nmu;
  
  if(nl>=2)
    return true;
  else
    return false;

}



// ------------ method called once each job just after ending the event loop  ------------
void 
RootTupleMaker::endJob() {
  if(debug_==true)
    cout << "Before write tree" << endl;
  
  m_file = m_tree->GetCurrentFile();
  
  if(debug_==true)
    cout << "Before fill tree with global info" << endl;

  //  cout << "Nstart: " << Nstart << endl;

  m_tree2->Fill();

  if(debug_==true)
    cout << "Tree with global info filled" << endl;

  Int_t ret = m_file->Write();
  if(debug_==true)
    cout << "Both trees saved." << ret <<  endl;

  m_file->Close();

}


//define this as a plug-in
DEFINE_FWK_MODULE(RootTupleMaker);
