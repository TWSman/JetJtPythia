#include "TH1D.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "TMath.h"
#include "TFile.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Pythia8/Pythia.h"
#include <TStopwatch.h>

#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
/*
#include "fastjet/config.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include <fastjet/SISConePlugin.hh>
#include <fastjet/CDFMidPointPlugin.hh>
#ifdef FASTJET_VERSION
#include <fastjet/Selector.hh>
#include <fastjet/FunctionOfPseudoJet.hh>
#include <fastjet/tools/BackgroundEstimatorBase.hh>
#include <fastjet/tools/Subtractor.hh>
*/
#include <iostream>

#include <TClonesArray.h>
#include "src/AliJCDijetHistos.h"
#include "src/AliJBaseTrack.h"
//#include "src/JHistos.h"
#include "src/AliJCard.h"
#include "src/AliJJet.h"
#include "src/AliJHistManager.h"

#define DEBUG 0

using namespace fastjet;
using namespace std;
using namespace Pythia8; 
Double_t getDiffR(double phi1, double phi2, double eta1, double eta2);
int GetBin(TVector *array, double val){
      int iBin=-1;
      for(int i=1; i< array->GetNoElements(); i++){
        if((*array)[i] <= val && val<(*array)[i+1]){
          iBin=i-1;
          break;
        }
      }
                  
      return iBin;
    }
void CalculateJetsJt(TClonesArray *inList,
    vector<fastjet::PseudoJet> chparticles,
    int    lDebug,
    int    lCBin,
    double lParticleEtaCut,
    double lParticlePtCut,
    double lJetCone,
    double lConstituentCut,
    double lJetEtaCut,
    int iContainer,
    int doLeadingRef);

class AliJCDijetHistos;

/*
   class MyUserInfo : public PseudoJet::UserInfoBase{
   public:
// default ctor
MyUserInfo(const int & pdg_id_in,const vector<bool> &pType) :
_pdg_id(pdg_id_in){ ispid = pType;}

/// access to the PDG id
int pdg_id() const { return _pdg_id;}
void PrintIsPid() const { 
for(unsigned int i=0;i<ispid.size();i++) {
cout << ispid.at(i)<<":";
}
cout << endl;
}
bool IsType(int i) const { return ispid[i];}

protected:
int _pdg_id;         // the associated pdg id
vector<bool> ispid;
};
*/
AliJCDijetHistos *fhistos;


int main(int argc, char **argv) {

  if(argc<4){
    cout<<"usage: " << argv[0] << " pythia.config cardname <output.root> [random_seed]"<<endl;exit(1);
  }
  TStopwatch timer;
  timer.Start();

  char* pythiaconfig  = argv[1];
  char *cardName = argv[2];
  TString outputs = argv[3];
  Int_t random_seed = argc>4 ? atoi(argv[4]) : 0;//placing the inputs into variables
  cout<<"card_name input"<<endl;
  AliJCard *fCard = new AliJCard(cardName);
  fCard->PrintOut();
  fCard->ReCompile();
  fCard->PrintOut();
  TVector *fJetTriggPtBorders = fCard->GetVector("JetTriggPtBorders");
  Double_t fEff = fCard->Get("fEff");


  TFile *fout = new TFile(outputs.Data(),"RECREATE");
  fout->cd();//opening of the output file
  TDirectoryFile *fdir = new TDirectoryFile( "JCDijetBaseTask","JCDijetBaseTask" );
  fdir->cd();

  cout << pythiaconfig << endl;
  //---------------------
  //Pythia initialization 
  //---------------------
  Pythia pythia;   // Generator.
  //Event& event      = pythia.event;
  ParticleData& pdt = pythia.particleData;

  cout << pythiaconfig << endl;
  // Read in commands from external file.
  pythia.readFile(pythiaconfig);

  // Extract settings to be used in the main program.
  int    nEvent  = pythia.mode("Main:numberOfEvents");
  bool   showCS  = pythia.flag("Init:showChangedSettings");
  bool   showCPD = pythia.flag("Init:showChangedParticleData");
  double energy  = pythia.parm("Beams:eCM");

  //pythia.readString(Form("PhaseSpace:pTHatMin ==%f",3));
  //pythia.readString(Form("PhaseSpace:pTHatMax ==%f",-1));
  cout<<"Events="<<nEvent <<" RNDM seed "<< random_seed << endl;

  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed=%02d",random_seed));

  // Initialize. Beam parameters set in .cmnd file.
  pythia.init();

  // List changed data. 
  //if (showCS)  pythia.settings.listChanged();
  //if (showCPD) pdt.listChanged();

  //-------------------------------------------------------
  // Histograms and tools
  //-------------------------------------------------------
  fhistos = new AliJCDijetHistos(fCard);
  vector<double> centbins = {0.0, 100.0};
  fhistos->SetCentralityBinsHistos(centbins);
  fhistos->CreateEventTrackHistos();

  fhistos->fHMG->Print();

  TH1D *hCrossSectionInfo = new TH1D("hCrossSection","CrossSectionInfo",8,0,8);

  //------------------------------------------------------------------
  // Define jet reconstruction
  //------------------------------------------------------------------
  TClonesArray *inputList = new TClonesArray("AliJBaseTrack",1500);

  double partMinPtCut         = 0.15;// atlas 0.5 cms/alice 0.15
  double partMinEtaCut        = 1.2;
  double coneR                = 0.4; // atlas 0.6, cms 0.7 alice 0.4
  double ktconeR              = 0.4;
  double fusePionMassInktjets = false;
  double fuseDeltaPhiBGSubtr  = false;
  double jetConstituentCut    = 5.0;
  double dijetSubleadingPt    = 20.0;
  double dijetDeltaPhiCut     = 2.0; // Cut is pi/dijetDeltaPhiCut
  double jetEtaCut            = 0.25;
  int fktScheme               = 1;
  int centBin=0;
  int doLeadingRef=0;

  TString sktScheme;
  switch (fktScheme) {
    case 0:  sktScheme = "E_scheme";
             break;
    case 1:  sktScheme = "pt_scheme";
             break;
    case 2:  sktScheme = "pt2_scheme";
             break;
    case 3:  sktScheme = "Et_scheme";
             break;
    case 4:  sktScheme = "Et2_scheme";
             break;
    case 5:  sktScheme = "BIpt_scheme";
             break;
    case 6:  sktScheme = "BIpt2_scheme";
             break;
    default: sktScheme = "Unknown, check macro arguments!";
             break;
  }

  cout << endl;
  cout << "============= Settings =============" << endl;
  cout << "cent bin:                   " << centBin << endl;
  cout << "particle eta cut:           " << partMinEtaCut << endl;
  cout << "particle pt cut:            " << Form("%.2f",partMinPtCut) << endl;
  cout << "jet cone size:              " << coneR << endl;
  cout << "kt-jet cone size:           " << ktconeR << endl;
  cout << "Using pion mass in kt-jets: " << fusePionMassInktjets << endl;
  cout << "Using DeltaPhi in BG subtr: " << fuseDeltaPhiBGSubtr << endl;
  cout << "Using kt-jet scheme:        " << sktScheme.Data() << endl;
  cout << "jet costituent cut:         " << jetConstituentCut << endl;
  cout << "dijet subleading pt cut:    " << dijetSubleadingPt << endl;
  cout << "dijet DeltaPhi cut:         pi/" << dijetDeltaPhiCut << endl;
  cout << endl;

  if(fusePionMassInktjets && fktScheme!=0) {
    cout << "Warning: Using pion mass for kt-jets but not using E_scheme!" << endl;
    cout << endl;
  }

  // Save information about the settings used.
  fhistos->fh_info->Fill("Count", 1.0);
  fhistos->fh_info->Fill("MC", 1.0);
  fhistos->fh_info->Fill("Cent bin border 00", 0.0);
  fhistos->fh_info->Fill("Cent bin border 01", 100.0);
  fhistos->fh_info->Fill("Jet cone", coneR);
  fhistos->fh_info->Fill("kt-jet cone", ktconeR);
  fhistos->fh_info->Fill("Scheme", fktScheme);
  fhistos->fh_info->Fill("Use pion mass", fusePionMassInktjets);
  fhistos->fh_info->Fill("Particle eta cut", partMinEtaCut);
  fhistos->fh_info->Fill("Particle pt cut", partMinPtCut);
  fhistos->fh_info->Fill("Subleading jet cut", dijetSubleadingPt);
  fhistos->fh_info->Fill("Const. cut", jetConstituentCut);
  fhistos->fh_info->Fill("Delta phi cut pi/",dijetDeltaPhiCut);

  // Initialize fh_events so that the bin order is correct
  fhistos->fh_events[0]->Fill("events",nEvent);
  fhistos->fh_events[0]->Fill("particles",0.0);
  fhistos->fh_events[0]->Fill("acc. particles",0.0);
  fhistos->fh_events[0]->Fill("no rho calc. events",0.0);
  fhistos->fh_events[0]->Fill("rho calc. events",0.0);
  fhistos->fh_events[0]->Fill("jets",0.0);
  fhistos->fh_events[0]->Fill("acc. jets",0.0);
  fhistos->fh_events[0]->Fill("const. cut jets",0.0);
  fhistos->fh_events[0]->Fill("bg. subtr. jets",0.0);
  fhistos->fh_events[0]->Fill("bg. subtr. const. cut jets",0.0);
  fhistos->fh_events[0]->Fill("kt-jets",0.0);
  fhistos->fh_events[0]->Fill("acc. kt-jets",0.0);
  fhistos->fh_events[0]->Fill("leading jet drop",0.0);
  fhistos->fh_events[0]->Fill("subleading jet drop",0.0);
  fhistos->fh_events[0]->Fill("raw dijets",0.0);
  fhistos->fh_events[0]->Fill("raw dijets leading cut",0.0);
  fhistos->fh_events[0]->Fill("raw acc. dijets",0.0);
  fhistos->fh_events[0]->Fill("raw deltaphi cut dijets",0.0);
  fhistos->fh_events[0]->Fill("bg. subtr. dijets",0.0);
  fhistos->fh_events[0]->Fill("bg. subtr. dijets leading cut",0.0);
  fhistos->fh_events[0]->Fill("bg. subtr. acc. dijets",0.0);
  fhistos->fh_events[0]->Fill("bg. subtr. deltaphi cut dijets",0.0);
  fhistos->fh_events[0]->Fill("bg. subtr. const. cut dijets",0.0);
  fhistos->fh_events[0]->Fill("bg. subtr. const. cut dijets leading cut",0.0);
  fhistos->fh_events[0]->Fill("bg. subtr. const. cut acc. dijets",0.0);
  fhistos->fh_events[0]->Fill("bg. subtr. const. cut deltaphi cut dijets",0.0);
  fhistos->fh_events[0]->Fill("const. cut dijets",0.0);
  fhistos->fh_events[0]->Fill("const. cut dijets leading cut",0.0);
  fhistos->fh_events[0]->Fill("const. cut acc. dijets",0.0);
  fhistos->fh_events[0]->Fill("const. cut deltaphi cut dijets",0.0);
  fhistos->fh_events[0]->Fill("kt dijets",0.0);
  fhistos->fh_events[0]->Fill("kt dijets leading cut",0.0);
  fhistos->fh_events[0]->Fill("kt acc. dijets",0.0);
  fhistos->fh_events[0]->Fill("kt deltaphi cut dijets",0.0);


  //--------------------------------------------------------
  //         B e g i n    e v e n t    l o o p.
  //--------------------------------------------------------
  cout<<"Let's start" <<endl; 
  int ieout = nEvent/20;
  if (ieout<1) ieout=1;
  int EventCounter = 0;
  Int_t nTried = 0; 
  Int_t prev_nTried = 0;
  Int_t nTrial = 0;
  Int_t nAccepted = 0;
  Float_t sigmaGen = 0.0;
  Float_t ebeweight = 1.0;

  TRandom3 *rand = new TRandom3();

  TVector fR = *fCard->GetVector("fR");
  int nR = fR.GetNoElements();
  vector<fastjet::PseudoJet> chparticles;

  for(int iEvent = 0; iEvent < nEvent; ++iEvent) {//begin event loop

    if (!pythia.next()) continue;
    inputList->Clear("C");
    nTried = pythia.info.nTried();
    nTrial = nTried - prev_nTried;
    prev_nTried = nTried;
    sigmaGen = pythia.info.sigmaGen();
    ebeweight = 1.0; //no event-by-event weight at all. //sigmaGen/nTrial;
    hCrossSectionInfo->Fill(7.5,ebeweight);
    if(iEvent % ieout == 0) cout << iEvent << "\t" << int(float(iEvent)/nEvent*100) << "%, nTried:" << nTried << ", nTrial:" << nTrial << ", sigma:" << sigmaGen << endl;

    for (int i = 0; i < pythia.event.size(); ++i) {//loop over all the particles in the event
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged() && pythia.event[i].isHadron() ) { // Only check if it is final, charged and hadron since the acceptance is checked in the CalculateJetsDijets
        if(rand->Uniform() > fEff) continue;
        TLorentzVector lParticle(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
        AliJBaseTrack track( lParticle );
        track.SetID(pythia.event[i].id());
        track.SetTrackEff(1.);
        new ((*inputList)[inputList->GetEntriesFast()]) AliJBaseTrack(track);
      }
    } // end of finalparticles
    int noTracks = inputList->GetEntries();
    double pt,eta,phi;

    chparticles.clear();
    fastjet::PseudoJet fTrack;
    for (int itrack = 0; itrack < noTracks; ++itrack) {//loop over all the particles in the event
      // Building input particle list for the jet reconstruction
      AliJBaseTrack *trk = (AliJBaseTrack*)inputList->At(itrack);
      pt = trk->Pt();
      eta = trk->Eta();
      fhistos->fh_events[centBin]->Fill("particles",1.0);
      if (pt>partMinPtCut && TMath::Abs(eta) < partMinEtaCut){
        fhistos->fh_events[centBin]->Fill("acc. particles",1.0);
        phi = trk->Phi();
        fhistos->fh_eta[centBin]->Fill(eta);
        fhistos->fh_phi[centBin]->Fill(phi);
        fhistos->fh_etaPhi[centBin]->Fill(eta,phi);
        fhistos->fh_pt[centBin]->Fill(pt);
        fTrack = fastjet::PseudoJet(trk->Px(),trk->Py(),trk->Pz(),trk->E());
        fTrack.set_user_index(itrack);
        chparticles.push_back(fTrack);
      }
    }
    if(chparticles.size()==0) continue; // We are not intereted in empty events.

    // Here I call my function
    for(int iR = 1 ; iR < nR+1 ; iR++){
      CalculateJetsJt(inputList,
          chparticles,
          5, // Debug
          centBin, // Cent bin
          partMinEtaCut, // Particle eta cut
          partMinPtCut, // Particle pt cut
          fR[iR], // Jet cone size
          jetConstituentCut, // Jet constituent cut
          jetEtaCut,     // Jet eta cut
          iR - 1, //iContainer
          doLeadingRef);
    }

    EventCounter++;
    if(iEvent == nEvent-1) cout << nEvent << "\t" << "100%, nTried:" << pythia.info.nTried() << ", sigma:" << pythia.info.sigmaGen() << endl ;
  }//event loop

  nTried = pythia.info.nTried();
  nAccepted = pythia.info.nAccepted();
  sigmaGen = pythia.info.sigmaGen();
  //double sigmaErr = pythia.info.sigmaErr();
  hCrossSectionInfo->Fill(0.5,nTried);
  cout << "nTried after loop:" << nTried << endl;// print also inside event loop and in the macro.
  hCrossSectionInfo->Fill(1.5,nAccepted);
  //cout << "nAccepted after loop:" << nAccepted << endl;
  hCrossSectionInfo->Fill(2.5,sigmaGen);
  cout << "sigma after loop:" << sigmaGen << endl;
  hCrossSectionInfo->Fill(3.5,EventCounter);
  hCrossSectionInfo->Fill(4.5,energy);
  hCrossSectionInfo->Fill(5.5,1); // for counting # of merged
  hCrossSectionInfo->Fill(6.5,pythia.info.weightSum()); // for counting # of merged

  fout->Write();
  fout->Close();
  cout << EventCounter << " events are analyzed successfully."<< endl;
  timer.Print(); 
  return 0;
}

//______________________________________________________________________________
void CalculateJetsJt(TClonesArray *inList,
    vector<fastjet::PseudoJet> chparticles,
    int    lDebug,
    int    lCBin,
    double lParticleEtaCut,
    double lParticlePtCut,
    double lJetCone,
    double lConstituentCut,
    double lJetEtaCut,
    int iContainer,
    int doLeadingRef) {

  double const etaMaxCutForJet = lParticleEtaCut-lJetCone;
  double const MinJetPt = 5.0; // Min Jet Pt cut to disregard low pt jets
  double const ghost_maxrap = lParticleEtaCut;
  unsigned int const repeat = 1; // default
  double const ghost_area   = 0.005; // ALICE=0.005 // default=0.01
  //double const pionmass = 0.139570; // From PDG


  double phi, eta, pt, rho, rhom, area;
  bool leadingTrackOverThreshold = false;
  vector<fastjet::PseudoJet> ktchparticles;
  vector<fastjet::PseudoJet> jets;
  vector<fastjet::PseudoJet> rhoEstJets;
  vector<fastjet::PseudoJet> constituents;
  fastjet::PseudoJet jetAreaVector;
  fastjet::PseudoJet dijet;
  fastjet::PseudoJet jet1;


  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Run the clustering, Reconstruct jets
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  fastjet::GhostedAreaSpec const area_spec(ghost_maxrap, repeat, ghost_area);
  fastjet::AreaDefinition const area_def(fastjet::active_area, area_spec);
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm,lJetCone,fastjet::pt_scheme);
  fastjet::ClusterSequenceArea cs(chparticles,jet_def,area_def);
  jets    = fastjet::sorted_by_pt(cs.inclusive_jets(MinJetPt)); // APPLY Min pt cut for jet
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Loop over jets and fill various histos 
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  fhistos->fh_rho[lCBin]->Fill(rho);
  fhistos->fh_rhom[lCBin]->Fill(rhom);
  double z; double jt; double jtleading; double zleading;
  TLorentzVector  vOrtho;
  TLorentzVector summedJet;
  int moveJet = 1;
  int iBin2=0;
  int iBin = 0;
  int leadingTrackIndex;
  int noTracks = inList->GetEntries();

  // anti-kt jets:
  for (unsigned ijet = 0; ijet < jets.size(); ijet++) {
    eta = jets[ijet].eta();
    if(TMath::Abs(eta) > lJetEtaCut) continue;
    AliJJet *jet =  new AliJJet(jets[ijet].px(),jets[ijet].py(), jets[ijet].pz(), jets[ijet].E(), ijet,0,0);
    jet->SetArea( jet->Area() );

    //== TRACK or Particle
    int nTrack = jets[ijet].constituents().size();
    for(unsigned it=0;it<nTrack; it++){
      int itrack = jets[ijet].constituents()[it].user_index();
      AliJBaseTrack *trk = (AliJBaseTrack*)inList->At(itrack);
      jet->AddConstituent(trk);
    }

    double conPtMax =0;
    fhistos->fh_events[lCBin]->Fill("jets",1.0);
    // anti-kt-jet eta cut
    int doBkg = 1;
    int counter = 0;
    vOrtho.SetVect(jet->Vect());
    vOrtho.SetE(jet->E());
    vOrtho.SetPhi(jet->Phi()+TMath::Pi()/2);
    if(TMath::Abs(eta) < etaMaxCutForJet) {
      fhistos->fh_events[lCBin]->Fill("acc. jets",1.0);
      pt = jets[ijet].pt();
      iBin = GetBin(fhistos->fJetTriggPtBorders,pt); // fill jetPt histos
      phi = jets[ijet].phi();
      area = jets[ijet].area();
      for (int icon = 0; icon<jet->GetConstituents()->GetEntries(); icon++){
        AliJBaseTrack *con = jet->GetConstituent(icon);
        if (con->Pt()>conPtMax) conPtMax = con->Pt();
        leadingTrackIndex = icon;
      }
      AliJBaseTrack *leadingTrack = jet->GetConstituent(leadingTrackIndex);
      double leadingTrackPt = jet->LeadingParticlePt(); //FIXME? For MC tracks this is possibly a track with no charge
      int jetMult = jet->GetConstituents()->GetEntries();
      iBin2 = GetBin(fhistos->fJetLeadPtBorders,leadingTrackPt);
      int iBin3 = GetBin(fhistos->fJetMultBorders,jetMult);
      jetAreaVector = jets[ijet].area_4vector();
      fhistos->fh_jetEta[lCBin][iContainer]->Fill(eta);
      fhistos->fh_jetPhi[lCBin][iContainer]->Fill(phi - TMath::Pi()); //Pseudojet.phi range 0-2pi
      fhistos->fh_jetEtaPhi[lCBin][iContainer]->Fill(eta,phi - TMath::Pi());
      fhistos->fh_jetPt[lCBin][iContainer]->Fill(pt);
      fhistos->fhJetPtBin[iContainer][iBin]->Fill(pt);
      fhistos->fh_jetArea[lCBin][iContainer]->Fill(area);
      fhistos->fh_jetAreaRho[lCBin][iContainer]->Fill(area*rho);
      leadingTrackOverThreshold=true;
      if(lDebug > 9) cout << "Jet i=" << ijet << ", jet pt=" << pt << endl;
      for(unsigned iconst=0;iconst<jets[ijet].constituents().size(); iconst++) {
        if(lDebug > 9) cout << "Constituent i=" << iconst << ", constituent pt=" << jets[ijet].constituents()[iconst].pt() << endl;
        if(jets[ijet].constituents()[iconst].pt() > lConstituentCut) { // Jet leading constituent cut.
          leadingTrackOverThreshold=true;
          break;
        }
      }


      double maxconpt = 0;

      for (int itrack = 0; itrack < noTracks; ++itrack) {//loop over all the particles in the event
        //cout << "itrack: " << itrack << endl;
        AliJBaseTrack *track = (AliJBaseTrack*)inList->At(itrack);
        double pta = track->Pt();
        eta = track->Eta();
        if (!track){
          cout << "No track " << endl;
          continue;
        }
        /*if(track->GetCharge() == 0){
          cout << "No Charge" << endl;
          continue;
          }*/
        if (pt<lParticlePtCut || TMath::Abs(eta) > lParticleEtaCut) continue;
        phi = track->Phi();
        if (pta > maxconpt) maxconpt = pta;
        int iptaBin = GetBin(fhistos->fJetAssocPtBorders, pta);
        if( iptaBin < 0 ) continue;
        z = (track->Vect()*jet->Vect().Unit())/jet->P();
        jt = (track->Vect()-z*jet->Vect()).Mag();
        //jT for all tracks in the event
        double deltaR   = getDiffR(jet->Phi(),track->Phi(),jet->Eta(),track->Eta());
        //cout << "jet Phi: " << jet->Phi() << "Eta: " << jet->Eta() << " Track Phi: " << track->Phi() << " Eta: " << track->Eta() << endl;
        //cout << "deltaR: " << deltaR << endl;
        if(deltaR < TMath::Pi()/2){
          fhistos->fhEventJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]->Fill(jt, 1.0/jt);
          fhistos->fhEventJtWeightBin[iContainer][iBin]->Fill(jt, 1.0/jt);
          fhistos->fhEventJtBin[iContainer][iBin]->Fill(jt,1);
        }
        //Jet Cone Jt here
        if ( deltaR < lJetCone){
          fhistos->fhJetConeTrkPt[iContainer]->Fill(pta,1);
          fhistos->fhJetConeTrkPtBin[iContainer][iBin]->Fill(pta,1);
          fhistos->fhJetConeTrkPtWeightBin[iContainer][iBin]->Fill(pta,1/pta);
          fhistos->fhJetConeZ[iContainer]->Fill( z , 1.0);
          fhistos->fhJetConeZBin[iContainer][iBin]->Fill( z , 1);
          fhistos->fhJetConeJt[iContainer]->Fill( jt , 1);
          fhistos->fhJetConeJtBin[iContainer][iBin]->Fill( jt , 1);
          fhistos->fhJetConeJtWeightBin[iContainer][iBin]->Fill( jt, 1.0/jt * 1);
          fhistos->fhJetConeJtWeight2D[iContainer]->Fill(jt,jet->Pt(),1.0/jt * 1);
          if(iBin2 > -1){
            fhistos->fhJetConeJtWeightWithTrackCutBinBin[iContainer][iBin][iBin2]->Fill( jt, 1.0/jt * 1);
          }
          if(iBin3 > -1){
            fhistos->fhJetConeJtWeightWithMultiplicityCutBinBin[iContainer][iBin][iBin3]->Fill( jt, 1.0/jt * 1);
          }
          if(pta < 0.99*leadingTrackPt && doLeadingRef){
            int xlongBin = GetBin(fhistos->fXlongBorders, pta/leadingTrackPt);
            if( xlongBin < 0 ) {
              continue;
            }
            zleading = (track->Vect()*leadingTrack->Vect().Unit())/leadingTrack->P();
            jtleading =  (track->Vect()-zleading*leadingTrack->Vect()).Mag();
            fhistos->fhJetConeJtLeadingRefBin[iContainer][iBin][xlongBin]->Fill(jtleading,1);
            fhistos->fhJetConeJtWeightLeadingRefBin[iContainer][iBin][xlongBin]->Fill(jtleading,1.0/jtleading);
            if(iBin2 > -1){
              fhistos->fhJetConeJtWeightLeadingRefWithTrackCutBinBin[iContainer][iBin][iBin2][xlongBin]->Fill(jtleading,1.0/jtleading);
            }
          }

          if (iptaBin < 0) continue;
          fhistos->fhJetConeJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
            ->Fill( jt, 1.0/jt);
        }

        if ( doBkg ){
          //Background jt
          deltaR   = getDiffR(vOrtho.Phi(),track->Phi(),vOrtho.Eta(),track->Eta());
          if ( deltaR < lJetCone){
            counter++;
            fhistos->fhBgTrkPt[iContainer]->Fill(pta,1);
            fhistos->fhBgTrkPtBin[iContainer][iBin]->Fill(pta,1);
            fhistos->fhBgTrkPtWeightBin[iContainer][iBin]->Fill(pta,1.0/pta);
            if(moveJet){
              summedJet = track->GetLorentzVector() + vOrtho;
              fhistos->fhBgRBin[iContainer][iBin]->Fill(getDiffR(summedJet.Phi(), track->Phi(), summedJet.Eta(), track->Eta()),1);
              z = (track->Vect()*summedJet.Vect().Unit())/summedJet.P();
              jt = (track->Vect()-z*summedJet.Vect()).Mag();
            }else{
              fhistos->fhBgRBin[iContainer][iBin]->Fill(deltaR,1);
              z = (track->Vect()*vOrtho.Vect().Unit())/vOrtho.P();
              jt = (track->Vect()-z*vOrtho.Vect()).Mag();
            }
            fhistos->fhBgZ[iContainer]->Fill( z , 1);
            fhistos->fhBgZBin[iContainer][iBin]->Fill( z , 1);
            fhistos->fhBgJt[iContainer]->Fill( jt , 1);
            fhistos->fhBgJtBin[iContainer][iBin]->Fill( jt , 1);
            fhistos->fhBgJtWeightBin[iContainer][iBin]->Fill( jt, 1.0/jt);

            if (iptaBin < 0) continue;
            fhistos->fhBgJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]
              ->Fill( jt, 1.0/jt);
          } //End of If
        } //End of background
      } //End of track loop
      if(doBkg){
        fhistos->fhBgTrkNumber[iContainer]->Fill(counter);
        fhistos->fhBgTrkNumberBin[iContainer][iBin]->Fill(counter);
      }
    }
  }//end of the anti-kt-jet loop
}

Double_t getDiffR(double phi1, double phi2, double eta1, double eta2){
  Double_t diffPhi = TMath::Abs(phi1-phi2);
  if(diffPhi > TMath::Pi()){
    diffPhi = 2*TMath::Pi() - diffPhi;
  }
  return TMath::Sqrt(TMath::Power(diffPhi,2)+TMath::Power(eta1-eta2,2));
}
