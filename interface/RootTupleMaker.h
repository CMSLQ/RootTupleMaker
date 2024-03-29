/* RooTupleMaker.h
  by: E. Twedt
  Nov. 7th, 2008
*/

// system include files
#include <memory>
#include <typeinfo>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "SimDataFormats/HepMCProduct/interface/GenInfoProduct.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/CachedProducts.h"

// Use for stringstream
#include <iostream>
#include <iomanip>

// root include statements
#include "TROOT.h"
#include "TVector3.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

//Event info
#include "DataFormats/HepMCCandidate/interface/PdfInfo.h"
#include <SimDataFormats/HepMCProduct/interface/HepMCProduct.h>
#include "DataFormats/HepMCCandidate/interface/PdfInfo.h"
//LHAPDF stuff

extern "C" {
  void initpdfset_(char*, int len);
  void initpdfsetm_(int &, char*);
  void initpdf_(int &);
  void evolvepdf_(double &, double &, double *);
  void numberpdf_(int &);
}


// GenParticle stuff
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// MET stuff
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/MET.h"

#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/METCollection.h"

// Jet stuff
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h" 

// Muon stuff
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
//#include "RecoMuon/MuonIdentification/interface/IdGlobalFunctions.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"

//Electron stuff
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

//Electron ID
#include "RecoEgamma/ElectronIdentification/interface/CutBasedElectronID.h"
#include "RecoEgamma/ElectronIdentification/interface/ElectronIDAlgo.h"
#include "AnalysisDataFormats/Egamma/interface/ElectronID.h"
#include "AnalysisDataFormats/Egamma/interface/ElectronIDAssociation.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"


//Electron Isolation
#include "DataFormats/Candidate/interface/CandAssociation.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronIsoCollection.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronIsoNumCollection.h"


// PAT stuff
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

// HLT
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include <FWCore/Framework/interface/TriggerNames.h> 
#include <DataFormats/Common/interface/TriggerResults.h> 

// NEW, possibly unnecessary
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

// About tree
#define MAXGENPARTICLES  2000
#define MAXHLTBITS    200
#define MAXELECTRONS  100
#define MAXGENJETS    100
#define MAXCALOJETS   100
#define MAXMUONS      100
