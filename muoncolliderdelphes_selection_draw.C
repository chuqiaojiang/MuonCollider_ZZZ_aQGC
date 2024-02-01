#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TObjArray.h"

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/types.h>
#include <errno.h>
#include <assert.h>

//#ifdef __CLING__
R__ADD_INCLUDE_PATH(/home/pku/jiangcq/trial_condor_madgraph/MG5_aMC_v3_5_0/Delphes);
R__ADD_INCLUDE_PATH(/home/pku/jiangcq/trial_condor_madgraph/MG5_aMC_v3_5_0/Delphes/external/);
R__ADD_LIBRARY_PATH(/home/pku/jiangcq/trial_condor_madgraph/MG5_aMC_v3_5_0/Delphes/);
#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
R__LOAD_LIBRARY(Delphes);
//#endif

double convert_weight_index_to_paravalue(int ind, double RWstart, double RWstep, int RWnstep)
{
	std::map<int,int> map_true_ind;
	map_true_ind[0] =  0;
	map_true_ind[1] =  1;
	map_true_ind[2] =  10;
	map_true_ind[3] =  11;
	map_true_ind[4] =  12;
	map_true_ind[5] =  13;
	map_true_ind[6] =  14;
	map_true_ind[7] =  15;
	map_true_ind[8] =  16;
	map_true_ind[9] =  17;
	map_true_ind[10] = 18;
	map_true_ind[11] = 19;
	map_true_ind[12] = 2;
	map_true_ind[13] = 20;
	map_true_ind[14] = 3;
	map_true_ind[15] = 4;
	map_true_ind[16] = 5;
	map_true_ind[17] = 6;
	map_true_ind[18] = 7;
	map_true_ind[19] = 8;
	map_true_ind[20] = 9;
	double paravalue = 0;
	paravalue = (RWstart+RWstep*(map_true_ind[ind]));//*(1e-12);
	return paravalue;
}




bool AnalyseEvents(ExRootTreeReader *treeReader, TString filename, int RWnstep)
{

	if( !(access(filename,0)) ) { 
		TFile file(filename, "read");
		TTree *tree = (TTree*)(file.Get("tree"));
		bool event_exist = ( (tree!=0 && (tree->GetEntries() > 0)) ? true : false );
		if(event_exist) cout << filename << " contains " << tree->GetEntries() << " entries " << endl;
		delete tree;
		file.Close();
		return event_exist;
	}


	Long64_t numberOfEntries = treeReader->GetEntries();

	// create a root file to store the variables
	TFile file(filename, "recreate");
	TTree *tree = new TTree("tree", Form("tree_from_%s", filename.Data()));

	double sum_weights = 0;
	double CX = 0;
	double SF = 1.0;

	int inputNum=0;
	int event;
	string lep_fst="";
	double M4l; 
	double M4l2j; 
	double deltaM; 
	double deltaM1234; 
	double deltaM1423; 
	double Ptll1; 
	double Ptll2; 
	double Pll1; 
	double Pll2; 
	double Etall1; 
	double Etall2; 
	double Mll1; 
	double Mll2; 
	double Pt4l;
	double P4l;
	double Eta4l;
	double E4l;
	double Pt4l2j;
	double P4l2j;
	double Eta4l2j;
	double E4l2j;
	double met; 
	double met_byhand; 
	double meta; 
	double Mjj; 
	double Ptjj; 
	double deltaRjj; 
	double costheta; 
	double Mrecoil; 
	double ptMu_ptlead; 
	double ptMu_pttail; 
	double etaMu_ptlead; 
	double etaMu_pttail; 
	double phiMu_ptlead; 
	double phiMu_pttail;
	double ptEl_ptlead; 
	double ptEl_pttail; 
	double etaEl_ptlead; 
	double etaEl_pttail; 
	double phiEl_ptlead; 
	double phiEl_pttail;
	double ptLep_ptlead; 
	double ptLep_pttail; 
	double etaLep_ptlead; 
	double etaLep_pttail; 
	double phiLep_ptlead; 
	double phiLep_pttail;
	double ptJ_ptlead; 
	double ptJ_pttail; 
	double etaJ_ptlead; 
	double etaJ_pttail; 
	double phiJ_ptlead; 
	double phiJ_pttail;
	double MassJ_ptlead;
	double MassJ_pttail;
	double deltaR12; 
	double deltaR34;
	double deltaR13;
	double deltaR24; 
	double deltaR_ll1;
	double deltaR_ll2; 
	int numMu; 
	int numLep; 
	int numGen; 
	int numEl;
	int numJ;
	std::vector<double>* weights = 0;
	double event_weight = 0;

	//std::vector<double> ptMu;
	//std::vector<double> etaMu;
	//std::vector<double> phiMu;
	//std::vector<double> chargeMu;
	//std::vector<double> px;
	//std::vector<double> py;
	//std::vector<double> pz;
	//std::vector<double> E;
	//std::vector<double> ptEl;
	//std::vector<double> etaEl;
	//std::vector<double> phiEl;
	//std::vector<double> chargeEl;

	tree->Branch("event", &event);
	tree->Branch("event_weight", &event_weight);
	tree->Branch("weights", &weights);
	tree->Branch("lep_fst", &lep_fst);
	tree->Branch("SF", &SF);
	tree->Branch("numMu", &numMu);
	tree->Branch("numEl", &numEl);
	tree->Branch("numJ", &numJ);
	tree->Branch("numLep", &numLep);
	tree->Branch("met",&met);
	tree->Branch("met_byhand",&met_byhand);
	tree->Branch("meta",&meta);
	tree->Branch("Mjj", &Mjj);
	tree->Branch("Ptjj",&Ptjj);
	tree->Branch("deltaRjj",&deltaRjj);
	tree->Branch("M4l", &M4l);
	tree->Branch("Pt4l",&Pt4l);
	tree->Branch("P4l",&P4l);
	tree->Branch("Eta4l",&Eta4l);
	tree->Branch("E4l",&E4l);
	tree->Branch("M4l2j", &M4l2j);
	tree->Branch("Pt4l2j",&Pt4l2j);
	tree->Branch("P4l2j",&P4l2j);
	tree->Branch("Eta4l2j",&Eta4l2j);
	tree->Branch("E4l2j",&E4l2j);
	tree->Branch("deltaM", &deltaM);
	tree->Branch("deltaM1234", &deltaM1234);
	tree->Branch("deltaM1423", &deltaM1423);
	tree->Branch("Mll1",&Mll1);
	tree->Branch("Ptll1",&Ptll1);
	tree->Branch("Pll1",&Pll1);
	tree->Branch("Etall1",&Etall1);
	tree->Branch("Mll2",&Mll2);
	tree->Branch("Ptll2",&Ptll2);
	tree->Branch("Pll2",&Pll2);
	tree->Branch("Etall2",&Etall2);
	tree->Branch("Mrecoil",&Mrecoil);
	tree->Branch("ptMu_ptlead",&ptMu_ptlead); 
	tree->Branch("ptMu_pttail",&ptMu_pttail); 
	tree->Branch("etaMu_ptlead",&etaMu_ptlead); 
	tree->Branch("etaMu_pttail",&etaMu_pttail); 
	tree->Branch("phiMu_ptlead",&phiMu_ptlead); 
	tree->Branch("phiMu_pttail",&phiMu_pttail);
	tree->Branch("ptEl_ptlead",&ptEl_ptlead); 
	tree->Branch("ptEl_pttail",&ptEl_pttail); 
	tree->Branch("etaEl_ptlead",&etaEl_ptlead); 
	tree->Branch("etaEl_pttail",&etaEl_pttail); 
	tree->Branch("phiEl_ptlead",&phiEl_ptlead); 
	tree->Branch("phiEl_pttail",&phiEl_pttail);
	tree->Branch("ptLep_ptlead",&ptLep_ptlead); 
	tree->Branch("ptLep_pttail",&ptLep_pttail); 
	tree->Branch("etaLep_ptlead",&etaLep_ptlead); 
	tree->Branch("etaLep_pttail",&etaLep_pttail); 
	tree->Branch("phiLep_ptlead",&phiLep_ptlead); 
	tree->Branch("phiLep_pttail",&phiLep_pttail);
	tree->Branch("ptJ_ptlead",&ptJ_ptlead); 
	tree->Branch("ptJ_pttail",&ptJ_pttail); 
	tree->Branch("etaJ_ptlead",&etaJ_ptlead); 
	tree->Branch("etaJ_pttail",&etaJ_pttail); 
	tree->Branch("phiJ_ptlead",&phiJ_ptlead); 
	tree->Branch("phiJ_pttail",&phiJ_pttail);
	tree->Branch("MassJ_ptlead",&MassJ_ptlead);
	tree->Branch("MassJ_pttail",&MassJ_pttail);
	tree->Branch("deltaR12",&deltaR12); 
	tree->Branch("deltaR34",&deltaR34);
	tree->Branch("deltaR13",&deltaR13);
	tree->Branch("deltaR24",&deltaR24); 
	tree->Branch("deltaR_ll1",&deltaR_ll1);
	tree->Branch("deltaR_ll2",&deltaR_ll2); 

	SF = 1.0/((double)numberOfEntries);

	event_weight = 0;






	TClonesArray *branchEvent = treeReader->UseBranch("Event");
	TClonesArray *branchWeight = treeReader->UseBranch("Weight");
	TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	TClonesArray *branchMuon = treeReader->UseBranch("Muon");
	TClonesArray *branchMET = treeReader->UseBranch("MissingET");
	//TClonesArray *branchJet = treeReader->UseBranch("Jet");
	TClonesArray *branchJet = treeReader->UseBranch("KTjet");

	inputNum = numberOfEntries;
	cout << filename << "'s input file contains " << numberOfEntries << " events" << endl;
	TLorentzVector l1;
	TLorentzVector l2;
	TLorentzVector l3;
	TLorentzVector l4;
	TLorentzVector j1;
	TLorentzVector j2;
	TLorentzVector lep_lv;
	TLorentzVector jet_lv;
	TLorentzVector met_lv_byhand;

	double Mz = 91.0 ;
	double MuMass=1.056600e-01; 
	double ElMass=5.110000e-04;
	int mup_ind=0;
	int mum_ind=0;
	int elp_ind=0;
	int elm_ind=0;

	for(int count=0; count < numberOfEntries; count++) {
		//cout<<"Processing event #"<<count<<endl;
		event=-1;
		event_weight=1;
		if(weights && !(weights->empty())) weights->clear();
		lep_fst="";
		met=-1000.0;
		meta=-1000.0;
		Mjj=-1000.0;
		Ptjj=-1000.0;
		deltaRjj=-1000.0;
		numMu=0;
		numEl=0;
		numJ=0;
		costheta=-1000.0;

		ptMu_ptlead=-1000.0; 
		ptMu_pttail=-1000.0; 
		etaMu_ptlead=-1000.0; 
		etaMu_pttail=-1000.0; 
		phiMu_ptlead=-1000.0; 
		phiMu_pttail=-1000.0;
		ptEl_ptlead=-1000.0; 
		ptEl_pttail=-1000.0; 
		etaEl_ptlead=-1000.0; 
		etaEl_pttail=-1000.0; 
		phiEl_ptlead=-1000.0; 
		phiEl_pttail=-1000.0;
		ptJ_ptlead=-1000.0; 
		ptJ_pttail=-1000.0; 
		etaJ_ptlead=-1000.0; 
		etaJ_pttail=-1000.0; 
		phiJ_ptlead=-1000.0; 
		phiJ_pttail=-1000.0;
		MassJ_ptlead=-1000.0; 
		MassJ_pttail=-1000.0; 
		deltaR12=-1000.0; 
		deltaR34=-1000.0;
		deltaR13=-1000.0;
		deltaR24=-1000.0; 
		deltaR_ll1=-1000.0;
		deltaR_ll2=-1000.0; 

		numGen=0;
		numLep=0;
		M4l=-1000.0;
		Pt4l=-1000.0;
		P4l=-1000.0;
		Eta4l=-1000.0;
		M4l2j=-1000.0;
		Pt4l2j=-1000.0;
		P4l2j=-1000.0;
		Eta4l2j=-1000.0;
		deltaM=-1000.0;
		deltaM1234=-1000.0;
		deltaM1423=-1000.0;
		Ptll1=-1000.0;
		Ptll2=-1000.0;
		Pll1=-1000.0;
		Pll2=-1000.0;
		Etall1=-1000.0;
		Etall2=-1000.0;
		Mll1=-1000.0;
		Mll2=-1000.0;
		Mrecoil=-1000.0;
		E4l=-1000.0; 
		E4l2j=-1000.0; 


		//if(!ptMu.empty())ptMu.clear(); 
		//if(!etaMu.empty())etaMu.clear(); 
		//if(!phiMu.empty())phiMu.clear(); 
		//if(!chargeMu.empty())chargeMu.clear(); 
		//if(!ptEl.empty())ptEl.clear(); 
		//if(!etaEl.empty())etaEl.clear(); 
		//if(!phiEl.empty())phiEl.clear(); 
		//if(!chargeEl.empty())chargeEl.clear(); 
		//if(!px.empty())px.clear(); 
		//if(!py.empty())py.clear(); 
		//if(!pz.empty())pz.clear(); 
		//if(!E.empty())E.clear();           

		//********************************************************************
		treeReader->ReadEntry(count) ;

		Event* ev = (Event*)branchEvent->At(0);
		LHEFEvent* lhefev = (LHEFEvent*)branchEvent->At(0);
		event=ev->Number;
		event_weight=lhefev->Weight;
		Weight* ew = 0;

		// !!!!!!!!!!!! CAUTION !!!!!!!!!!!!
		// The last 3 weights in the weights vector in branchWeight of a reweighted sample are merged weights?
		// The last third one is the merged weight of weights series named by paramerters setting
		// The last second is the merged weight of weights series named by rwgt name
		// The last first is the merged weight of original weight from default initial paramerters setting
		// !!!!!!!!!!!! CAUTION !!!!!!!!!!!!
		if( ((filename).Contains("reweightscan")) ) {
			//for(int k=(RWnstep+1); k<=((2*RWnstep)+1); k++) 
			for(int k=(RWnstep); k<=((2*RWnstep)); k++) 
			{
				ew = (Weight*)branchWeight->At(k);
				weights->push_back(ew->Weight);
			}
			ew = (Weight*)branchWeight->At(42);
			event_weight=ew->Weight;
		}

		if( !((filename).Contains("reweightscan")) ) {
			ew = (Weight*)branchWeight->At(0);
			weights->push_back(ew->Weight);
			ew = (Weight*)branchWeight->At(0);
			event_weight=ew->Weight;
		}


		MissingET* Met = (MissingET *)branchMET->At(0);
		met = Met->MET;
		meta = Met->Eta;
		numMu = branchMuon->GetEntries();
		numEl = branchElectron->GetEntries();
		numGen = branchGenParticle->GetEntries();
		numLep = numMu + numEl;
		numJ = branchJet->GetEntries();



		if( !( ((numMu==2)&&(numEl==2)) || ((numMu==4)&&(numEl==0)) || ((numMu==0)&&(numEl==4)) ) ) { 
			continue; 
		}



		Jet *jet;
		Muon* Mu;
		Electron* El;


		met_lv_byhand.SetPtEtaPhiM(0,0,0,0);
		j1.SetPtEtaPhiM(0,0,0,0);
		j2.SetPtEtaPhiM(0,0,0,0);
		int j_count = 0;
		bool jet_is_lepton = false;
		for(int j=0; j<numJ; j++){
			jet = (Jet*)branchJet->At(j);
			jet_lv.SetPtEtaPhiM(jet->PT,jet->Eta,jet->Phi,jet->Mass);
			met_lv_byhand = met_lv_byhand + jet_lv;
			jet_is_lepton = false;
			for(int i=0; i<numEl; i++){
				El = (Electron*)branchElectron->At(i);
				lep_lv.SetPtEtaPhiM(El->PT,El->Eta,El->Phi,ElMass);
				if( ( (jet_lv.DeltaR(lep_lv)) < 0.1 ) && ( fabs(jet_lv.Pt() - lep_lv.Pt()) < 1 ) ) {
					jet_is_lepton = true;
					break;
				}
			}
			for(int i=0; i<numMu; i++){
				Mu = (Muon*)branchMuon->At(i);
				lep_lv.SetPtEtaPhiM(Mu->PT,Mu->Eta,Mu->Phi,MuMass);
				if( ( (jet_lv.DeltaR(lep_lv)) < 0.1 ) && ( fabs(jet_lv.Pt() - lep_lv.Pt()) < 1 ) ) {
					jet_is_lepton = true;
					break;
				}
			}
			if(jet_is_lepton) continue;
			if(j_count==0) {
				j2.SetPtEtaPhiM(jet->PT,jet->Eta,jet->Phi,jet->Mass);
				j1.SetPtEtaPhiM(jet->PT,jet->Eta,jet->Phi,jet->Mass);
				j_count++;
				continue;
			}
			if(j_count==1) {
				if( (jet->PT >= j1.Pt()) ){
					j1.SetPtEtaPhiM(jet->PT,jet->Eta,jet->Phi,jet->Mass);
				}
				if( (jet->PT < j1.Pt()) ){
					j2.SetPtEtaPhiM(jet->PT,jet->Eta,jet->Phi,jet->Mass);
				}
				j_count++;
				continue;
			}
			if( j_count>=2 && (jet->PT > j2.Pt()) ) { 
				j2.SetPtEtaPhiM(jet->PT,jet->Eta,jet->Phi,jet->Mass); 
				if( jet->PT > j1.Pt() ) { 
					j2.SetPtEtaPhiM(j1.Pt(),j1.Eta(),j1.Phi(),j1.M());
					j1.SetPtEtaPhiM(jet->PT,jet->Eta,jet->Phi,jet->Mass); 
				}
				j_count++;
				continue;
			}
		}


		met_byhand = met_lv_byhand.P();


		ptJ_ptlead = j1.Pt();
		etaJ_ptlead = fabs(j1.Eta());
		phiJ_ptlead = j1.Phi();
		MassJ_ptlead = j1.M();
		ptJ_pttail = j2.Pt();
		etaJ_pttail = fabs(j2.Eta());
		phiJ_pttail = j2.Phi();
		MassJ_pttail = j2.M();
		Mjj = (j1+j2).M();
		Ptjj = (j1+j2).Pt();
		deltaRjj = j1.DeltaR(j2);




		lep_fst = string( Form("%dMu_%dEl", numMu, numEl) );
		mup_ind=0;
		mum_ind=0;
		elp_ind=0;
		elm_ind=0;
		if((numMu==2)&&(numEl==2)){
			//cout<<"lep_fst="<<lep_fst<<endl;
			for(int j=0; j<2; j++){
				Mu = (Muon*)branchMuon->At(j);
				if(mup_ind==0 && Mu->Charge == 1) { l1.SetPtEtaPhiM(Mu->PT,Mu->Eta,Mu->Phi,MuMass); mup_ind+=1; continue; }
				if(mum_ind==0 && Mu->Charge ==-1) { l2.SetPtEtaPhiM(Mu->PT,Mu->Eta,Mu->Phi,MuMass); mum_ind+=1; continue; }

			}
			for(int j=0; j<2; j++){
				El = (Electron*)branchElectron->At(j);
				if(elp_ind==0 && El->Charge == 1) { l3.SetPtEtaPhiM(El->PT,El->Eta,El->Phi,ElMass); elp_ind+=1; continue; }
				if(elm_ind==0 && El->Charge ==-1) { l4.SetPtEtaPhiM(El->PT,El->Eta,El->Phi,ElMass); elm_ind+=1; continue; }

			}
			if( (elp_ind-elm_ind)!=0 || (mup_ind-mum_ind)!=0 ) 
			{ 
				//cout<<"Charges do not conserve!"<<endl;
				//cout<<"count = "<<count<<endl;
				continue;
			}
			deltaM = fabs((l1+l2).M()-Mz) + fabs((l3+l4).M()-Mz);
			deltaM1234=deltaM;
			deltaM1423=deltaM;
			Mll1=(l1+l2).M();
			Ptll1=(l1+l2).Pt();
			Pll1=(l1+l2).P();
			Etall1=(l1+l2).Eta();
			deltaR_ll1=l1.DeltaR(l2);
			Mll2=(l3+l4).M();
			Ptll2=(l3+l4).Pt();
			Pll2=(l3+l4).P();
			Etall2=(l3+l4).Eta();
			deltaR_ll2=l3.DeltaR(l4);
			ptMu_ptlead = l1.Pt();
			etaMu_ptlead = fabs(l1.Eta());
			phiMu_ptlead = l1.Phi();
			ptMu_pttail = l3.Pt();
			etaMu_pttail = fabs(l3.Eta());
			phiMu_pttail = l3.Phi();
			ptEl_ptlead = l2.Pt();
			etaEl_ptlead = fabs(l2.Eta());
			phiEl_ptlead = l2.Phi();
			ptEl_pttail = l4.Pt();
			etaEl_pttail = fabs(l4.Eta());
			phiEl_pttail = l4.Phi();
			if(l1.Pt() < l3.Pt()) {
				ptMu_ptlead = l3.Pt();
				etaMu_ptlead = fabs(l3.Eta());
				phiMu_ptlead = l3.Phi();
				ptMu_pttail = l1.Pt();
				etaMu_pttail = fabs(l1.Eta());
				phiMu_pttail = l1.Phi();
			}
			if(l2.Pt() < l4.Pt()) {
				ptEl_ptlead = l4.Pt();
				etaEl_ptlead = fabs(l4.Eta());
				phiEl_ptlead = l4.Phi();
				ptEl_pttail = l2.Pt();
				etaEl_pttail = fabs(l2.Eta());
				phiEl_pttail = l2.Phi();
			}

			if(ptEl_ptlead > ptMu_ptlead) {
				ptLep_ptlead  = ptEl_ptlead ; 
				etaLep_ptlead = etaEl_ptlead;
				phiLep_ptlead = phiEl_ptlead;
			}
			else {
				ptLep_ptlead  = ptMu_ptlead ; 
				etaLep_ptlead = etaMu_ptlead;
				phiLep_ptlead = phiMu_ptlead;
			}

			if(ptEl_pttail < ptMu_pttail) {
				ptLep_pttail  = ptEl_pttail ;
				etaLep_pttail = etaEl_pttail;
				phiLep_pttail = phiEl_pttail;
			}
			else {
				ptLep_ptlead  = ptMu_pttail ; 
				etaLep_ptlead = etaMu_pttail;
				phiLep_ptlead = phiMu_pttail;
			}

		}
		if((numMu==4)&&(numEl==0)){
			for(int j=0; j<4; j++){
				Mu = (Muon*)branchMuon->At(j);
				if(mup_ind==0 && Mu->Charge == 1) { l1.SetPtEtaPhiM(Mu->PT,Mu->Eta,Mu->Phi,MuMass); mup_ind+=1; continue; }
				if(mum_ind==0 && Mu->Charge ==-1) { l2.SetPtEtaPhiM(Mu->PT,Mu->Eta,Mu->Phi,MuMass); mum_ind+=1; continue; }
				if(mup_ind==1 && Mu->Charge == 1) { l3.SetPtEtaPhiM(Mu->PT,Mu->Eta,Mu->Phi,MuMass); mup_ind+=1; continue; }
				if(mum_ind==1 && Mu->Charge ==-1) { l4.SetPtEtaPhiM(Mu->PT,Mu->Eta,Mu->Phi,MuMass); mum_ind+=1; continue; }
			}
			if( (mup_ind-mum_ind)!=0 ) 
			{ 
				continue;
			}
			deltaM1234 = fabs((l1+l2).M()-Mz) + fabs((l3+l4).M()-Mz);
			deltaM1423 = fabs((l1+l4).M()-Mz) + fabs((l2+l3).M()-Mz);
			if(deltaM1234>=deltaM1423){
				deltaM=deltaM1423;
				Mll1=(l1+l4).M();
				Ptll1=(l1+l4).Pt();
				Pll1=(l1+l4).P();
				Etall1=(l1+l4).Eta();
				deltaR_ll1=l1.DeltaR(l4);
				Mll2=(l2+l3).M();
				Ptll2=(l2+l3).Pt();
				Pll2=(l2+l3).P();
				Etall2=(l2+l3).Eta();
				deltaR_ll2=l2.DeltaR(l3);
			}
			if(deltaM1234<deltaM1423){
				deltaM=deltaM1234;
				Mll1=(l1+l2).M();
				Ptll1=(l1+l2).Pt();
				Pll1=(l1+l2).P();
				Etall1=(l1+l2).Eta();
				deltaR_ll1=l1.DeltaR(l2);
				Mll2=(l3+l4).M();
				Ptll2=(l3+l4).Pt();
				Pll2=(l3+l4).P();
				Etall2=(l3+l4).Eta();
				deltaR_ll2=l3.DeltaR(l4);
			}
			ptMu_ptlead = l1.Pt();
			etaMu_ptlead = fabs(l1.Eta());
			phiMu_ptlead = l1.Phi();
			ptMu_pttail = l2.Pt();
			etaMu_pttail = fabs(l2.Eta());
			phiMu_pttail = l2.Phi();
			if(l2.Pt() > l1.Pt()) {
				ptMu_ptlead = l2.Pt();
				etaMu_ptlead = fabs(l2.Eta());
				phiMu_ptlead = l2.Phi();
				ptMu_pttail = l1.Pt();
				etaMu_pttail = fabs(l1.Eta());
				phiMu_pttail = l1.Phi();
			}
			if(l3.Pt() > ptMu_ptlead) {
				ptMu_ptlead = l3.Pt();
				etaMu_ptlead = fabs(l3.Eta());
				phiMu_ptlead = l3.Phi();
			}
			if(l3.Pt() < ptMu_pttail) {
				ptMu_pttail = l3.Pt();
				etaMu_pttail = fabs(l3.Eta());
				phiMu_pttail = l3.Phi();
			}
			if(l4.Pt() > ptMu_ptlead) {
				ptMu_ptlead = l4.Pt();
				etaMu_ptlead = fabs(l4.Eta());
				phiMu_ptlead = l4.Phi();
			}
			if(l4.Pt() < ptMu_pttail) {
				ptMu_pttail = l4.Pt();
				etaMu_pttail = fabs(l4.Eta());
				phiMu_pttail = l4.Phi();
			}

			{
				ptLep_ptlead  = ptMu_ptlead ; 
				etaLep_ptlead = etaMu_ptlead;
				phiLep_ptlead = phiMu_ptlead;
				ptLep_ptlead  = ptMu_pttail ; 
				etaLep_ptlead = etaMu_pttail;
				phiLep_ptlead = phiMu_pttail;
			}
		}

		if((numMu==0)&&(numEl==4)){
			//cout<<"lep_fst="<<lep_fst<<endl;
			for(int j=0; j<4; j++){
				El = (Electron*)branchElectron->At(j);
				if(elp_ind==0 && El->Charge == 1) { l1.SetPtEtaPhiM(El->PT,El->Eta,El->Phi,ElMass); elp_ind+=1; continue; }
				if(elm_ind==0 && El->Charge ==-1) { l2.SetPtEtaPhiM(El->PT,El->Eta,El->Phi,ElMass); elm_ind+=1; continue; }
				if(elp_ind==1 && El->Charge == 1) { l3.SetPtEtaPhiM(El->PT,El->Eta,El->Phi,ElMass); elp_ind+=1; continue; }
				if(elm_ind==1 && El->Charge ==-1) { l4.SetPtEtaPhiM(El->PT,El->Eta,El->Phi,ElMass); elm_ind+=1; continue; }
			}
			if( (elp_ind-elm_ind)!=0 )  
			{
				continue;
			}
			deltaM1234 = fabs((l1+l2).M()-Mz) + fabs((l3+l4).M()-Mz);
			deltaM1423 = fabs((l1+l4).M()-Mz) + fabs((l2+l3).M()-Mz);
			if(deltaM1234>=deltaM1423){
				deltaM=deltaM1423;
				Mll1=(l1+l4).M();
				Ptll1=(l1+l4).Pt();
				Pll1=(l1+l4).P();
				Etall1=(l1+l4).Eta();
				deltaR_ll1=l1.DeltaR(l4);
				Mll2=(l2+l3).M();
				Ptll2=(l2+l3).Pt();
				Pll2=(l2+l3).P();
				Etall2=(l2+l3).Eta();
				deltaR_ll2=l2.DeltaR(l3);
			}
			if(deltaM1234<deltaM1423){
				deltaM=deltaM1234;
				Mll1=(l1+l2).M();
				Ptll1=(l1+l2).Pt();
				Pll1=(l1+l2).P();
				Etall1=(l1+l2).Eta();
				deltaR_ll1=l1.DeltaR(l2);
				Mll2=(l3+l4).M();
				Ptll2=(l3+l4).Pt();
				Pll2=(l3+l4).P();
				Etall2=(l3+l4).Eta();
				deltaR_ll2=l3.DeltaR(l4);
			}
			ptEl_ptlead = l1.Pt();
			etaEl_ptlead = fabs(l1.Eta());
			phiEl_ptlead = l1.Phi();
			ptEl_pttail = l2.Pt();
			etaEl_pttail = fabs(l2.Eta());
			phiEl_pttail = l2.Phi();
			if(l2.Pt() > l1.Pt()) {
				ptEl_ptlead = l2.Pt();
				etaEl_ptlead = fabs(l2.Eta());
				phiEl_ptlead = l2.Phi();
				ptEl_pttail = l1.Pt();
				etaEl_pttail = fabs(l1.Eta());
				phiEl_pttail = l1.Phi();
			}
			if(l3.Pt() > ptEl_ptlead) {
				ptEl_ptlead = l3.Pt();
				etaEl_ptlead = fabs(l3.Eta());
				phiEl_ptlead = l3.Phi();
			}
			if(l3.Pt() < ptEl_pttail) {
				ptEl_pttail = l3.Pt();
				etaEl_pttail = fabs(l3.Eta());
				phiEl_pttail = l3.Phi();
			}
			if(l4.Pt() > ptEl_ptlead) {
				ptEl_ptlead = l4.Pt();
				etaEl_ptlead = fabs(l4.Eta());
				phiEl_ptlead = l4.Phi();
			}
			if(l4.Pt() < ptEl_pttail) {
				ptEl_pttail = l4.Pt();
				etaEl_pttail = fabs(l4.Eta());
				phiEl_pttail = l4.Phi();
			}

			{
				ptLep_ptlead  = ptEl_ptlead ; 
				etaLep_ptlead = etaEl_ptlead;
				phiLep_ptlead = phiEl_ptlead;
				ptLep_ptlead  = ptEl_pttail ; 
				etaLep_ptlead = etaEl_pttail;
				phiLep_ptlead = phiEl_pttail;
			}
		}


		if(Mll1 < Mll2) {
			double temp_Mll=Mll1;
			double temp_Ptll=Ptll1;
			double temp_Pll=Pll1;
			double temp_Etall=Etall1;
			double temp_deltaR_ll=deltaR_ll1;
			Mll1=Mll2;
			Ptll1=Ptll2;
			Pll1=Pll2;
			Etall1=Etall2;
			deltaR_ll1=deltaR_ll2;
			Mll2=temp_Mll;
			Ptll2=temp_Ptll;
			Pll2=temp_Pll;
			Etall2=temp_Etall;
			deltaR_ll2=temp_deltaR_ll;
		}

		Pt4l=(l1+l2+l3+l4).Pt();
		P4l=(l1+l2+l3+l4).P();
		Eta4l=(l1+l2+l3+l4).Eta();
		M4l=(l1+l2+l3+l4).M();
		E4l=(l1+l2+l3+l4).E();
		Mrecoil=sqrt(pow(1000.0-E4l,2)-pow(P4l,2));
		deltaR12=l1.DeltaR(l2);
		deltaR34=l3.DeltaR(l4);
		deltaR13=l1.DeltaR(l3);
		deltaR24=l2.DeltaR(l4);
		Pt4l2j=(l1+l2+l3+l4+j1+j2).Pt();
		P4l2j=(l1+l2+l3+l4+j1+j2).P();
		Eta4l2j=(l1+l2+l3+l4+j1+j2).Eta();
		M4l2j=(l1+l2+l3+l4+j1+j2).M();
		E4l2j=(l1+l2+l3+l4+j1+j2).E();



		//////////////////////////////////////////////////////////////////////////////
		tree->Fill();
	}
	bool event_exist = ( (tree!=0 && (tree->GetEntries() > 0)) ? true : false );
	if(event_exist) cout << filename << " contains " << tree->GetEntries() << " entries " << endl;
	file.Write();
	file.Close();

	return event_exist;
}



std::map<string, std::pair<std::pair<std::pair<double, double>, std::pair<double, double>>, double>> OptimizeCuts(std::vector<TString> filesin, std::vector<int> sample_types, std::vector<string> variables_for_cut, string& optcut, double ndivisions, string para_name, double RWstart, double RWstep, int RWnstep) {

	gStyle->SetOptStat(0000);

	int nfiles = filesin.size();
	int signfiles=0;
	for(int fi=0; fi<nfiles; fi++){
		if(sample_types.at(fi)==0 || sample_types.at(fi)==2) signfiles++;
	}
	int bkgnfiles=0;
	for(int fi=0; fi<nfiles; fi++){
		if(sample_types.at(fi)==1) bkgnfiles++;
	}

	TFile fin1(filesin.at(0), "READ");
	TTree* tin1 = (TTree*)fin1.Get("tree");
	TObjArray* objarr = tin1->GetListOfBranches();
	TObjArray* objarr1 = (TObjArray*)objarr->Clone();
	int objnum = objarr1->GetEntries();
	delete tin1;
	fin1.Close();

	std::map<string, std::pair<std::pair<std::pair<double, double>, std::pair<double, double>>, double>> output;
	std::map<string, std::pair<std::pair<std::pair<double, double>, std::pair<double, double>>, double>> output_sorted;
	std::pair<double, double> cut_limits;
	std::pair<double, double> n_sig_and_n_bkg;
	double significance=0;
	double n_sig_for_significance=0;
	double n_bkg_for_significance=0;
	string var_name="";
	double left_cut=0;
	double right_cut=0;
	double xmin=0;
	double xmax=100;
	double xrange=0;
	double n_sig=0;
	double n_bkg=0;
	double integral_left=0;
	double integral_right=100;
	double width_step=0;
	for(int k=0; k<objnum; k++){
		bool is_variable_for_cut=false;
		for(int icv=0; icv<variables_for_cut.size(); icv++) {
			if(string(variables_for_cut.at(icv)) == string("*")) {
				is_variable_for_cut = true;
				break;
			}
			if(string((objarr1->At(k))->GetName()) == string(variables_for_cut.at(icv))) {
				is_variable_for_cut = true;
			}
		}
		if(!is_variable_for_cut) continue;
		std::vector<TH1D*> opt_temps;
		xmin = 0;
		xmax = 100;
		xrange = xmax - xmin;
		width_step = xrange/ndivisions;
		for(int fid=0; fid<nfiles; fid++){ 
			if( !((TString(filesin.at(fid))).Contains("mumuTozzzTo")) ) continue;
			if( sample_types.at(fid) != 2 ) continue;
			TChain chin("tree");
			chin.Add(filesin.at(fid));
			//chin.SetBranchStatus("*", 0);
			//chin.SetBranchStatus(((objarr1->At(k)))->GetName(), 1);
			//chin.SetBranchStatus("SF", 1);
			//chin.SetBranchStatus("event_weight", 1);
			chin.SetBranchStatus("*", 1);
			if( (string(((objarr1->At(k)))->GetName()) != string("lep_fst")) ) {
				chin.Draw(Form("%s>>opt_temp%d", (objarr1->At(k))->GetName(), k), Form("10000000*SF*(%s > -99)*(weights[20])*(%s)",(objarr1->At(k))->GetName(),optcut.c_str()), "goff");
				//double xmean = ((TH1F*)(gDirectory->Get(Form("th_temp%d", k))))->GetMean();
				//double xstaddev = ((TH1F*)(gDirectory->Get(Form("th_temp%d", k))))->GetStdDev();
				//xmin = chin.GetMinimum(((objarr1->At(k)))->GetName());
				//xmax = chin.GetMaximum(((objarr1->At(k)))->GetName());
				if( ((TH1F*)(gDirectory->Get(Form("opt_temp%d", k)))) != 0 ) {
					xmin = ((TH1F*)(gDirectory->Get(Form("opt_temp%d", k))))->GetBinLowEdge(1);
					xmax = ((TH1F*)(gDirectory->Get(Form("opt_temp%d", k))))->GetBinLowEdge( ((TH1F*)(gDirectory->Get(Form("opt_temp%d", k))))->GetNcells() );
					xrange = xmax - xmin;
					width_step = xrange/ndivisions;
					xmin -= 2.*width_step;
					xmax += 2.*width_step;
				}
				else{
					xmin = 0;
					xmax = 100;
					xrange = xmax - xmin;
					width_step = xrange/ndivisions;
					xmin -= 2.0*width_step;
					xmax += 2.0*width_step;
				}
			}
			else if( (string(((objarr1->At(k)))->GetName()) == string("lep_fst")) ) {
				xmin = 0;
				xmax = 100;
				xrange = xmax - xmin;
				width_step = xrange/ndivisions;
				xmin -= 2.0*width_step;
				xmax += 2.0*width_step;
			}
		}
		for(int fid=0; fid<nfiles; fid++){ 
			if( sample_types.at(fid) < 0 || sample_types.at(fid) > 2 ) continue;
			TChain chin("tree");
			chin.Add(filesin.at(fid));
			//chin.SetBranchStatus("*", 0);
			//chin.SetBranchStatus(((objarr1->At(k)))->GetName(), 1);
			//chin.SetBranchStatus("SF", 1);
			//chin.SetBranchStatus("event_weight", 1);
			chin.SetBranchStatus("*", 1);
			opt_temps.push_back(new TH1D(Form("opt_file%d_%s",fid,((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
			if( sample_types.at(fid) != 2 ) {
				if(string(((objarr1->At(k)))->GetName()) != string("lep_fst")) {
					if( ((TH1F*)(gDirectory->Get(Form("opt_temp%d", k)))) != 0 ) {
						chin.Draw(Form("%s>>%s", ((objarr1->At(k)))->GetName(), Form("opt_file%d_%s",fid,((objarr1->At(k)))->GetName())), Form("10000000*SF*(%s > -99)*(event_weight)*(%s)",((objarr1->At(k)))->GetName(),optcut.c_str()), "goff");
					}
					else {
						((TH1F*)(gDirectory->Get(Form("opt_file%d_%s",fid,((objarr1->At(k)))->GetName()))))->Fill(1E7);
					}
				}
				else { 
					chin.Draw(Form("%s>>%s", ((objarr1->At(k)))->GetName(), Form("opt_file%d_%s",fid,((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",optcut.c_str()), "goff");
				}
			}
			if( sample_types.at(fid) == 2 ) {
				if(string(((objarr1->At(k)))->GetName()) != string("lep_fst")) {
					if( ((TH1F*)(gDirectory->Get(Form("opt_temp%d", k)))) != 0 ) {
						chin.Draw(Form("%s>>%s", ((objarr1->At(k)))->GetName(), Form("opt_file%d_%s",fid,((objarr1->At(k)))->GetName())), Form("10000000*SF*(%s > -99)*(weights[20])*(%s)",((objarr1->At(k)))->GetName(),optcut.c_str()), "goff");
					}
					else {
						((TH1F*)(gDirectory->Get(Form("opt_file%d_%s",fid,((objarr1->At(k)))->GetName()))))->Fill(1E7);
					}
				}
				else { 
					chin.Draw(Form("%s>>%s", ((objarr1->At(k)))->GetName(), Form("opt_file%d_%s",fid,((objarr1->At(k)))->GetName())), Form("10000000*SF*(weights[20])*(%s)",optcut.c_str()), "goff");
				}
			}

		}
		significance = 0;
		n_sig_for_significance = 0;
		n_bkg_for_significance = 0;
		var_name = string( ((objarr1->At(k)))->GetName() );
		//for(int nwidth=1; nwidth<=ndivisions; nwidth++) 
		//for(int nstart=1; (nstart+nwidth-1)<=(ndivisions+4); nstart++) 

		for(int nstart=3; nstart<=(ndivisions+2); nstart++) {
			for(int nwidth=1; (nstart+nwidth-1)<=(ndivisions+4); nwidth++) {
				if( ( ((TString(var_name)).Contains("M4l")) || ((TString(var_name)).Contains("Pt4l")) || ((TString(var_name)).Contains("P4l")) || ((TString(var_name)).Contains("deltaR")) ) ) {
					nwidth=(ndivisions+4+1-nstart);
				}
				n_sig = 0;
				n_bkg = 0;
				for(int fid=0; fid<nfiles; fid++) { 
					if( sample_types.at(fid) == 2 ) { 
						n_sig += ( (opt_temps.at(fid))->Integral(nstart,(nstart+nwidth-1)) );
					}
					if( sample_types.at(fid) == 0 || sample_types.at(fid) == 1 ) { 
						n_bkg += ( (opt_temps.at(fid))->Integral(nstart,(nstart+nwidth-1)) );
					}
				}
				if( (n_sig>0 && n_bkg>0) && ((n_sig/sqrt(n_bkg)) > significance) ) { 
					significance = (n_sig/sqrt(n_bkg));
					n_sig_for_significance = n_sig;
					n_bkg_for_significance = n_bkg;
					integral_left = (xmin+((nstart-1)*width_step));
					integral_right = (xmin+((nstart-1)*width_step)+(nwidth*width_step));
				}
				if( ( ((TString(var_name)).Contains("M4l")) || ((TString(var_name)).Contains("Pt4l")) || ((TString(var_name)).Contains("P4l")) || ((TString(var_name)).Contains("deltaR")) ) ) {
					break;
				}
			}
		}
		cut_limits.first  = integral_left;
		cut_limits.second = integral_right;
		n_sig_and_n_bkg.first  = n_sig_for_significance;
		n_sig_and_n_bkg.second = n_bkg_for_significance;
		(output[var_name]).first.first  = cut_limits;
		(output[var_name]).first.second = n_sig_and_n_bkg;
		(output[var_name]).second = significance;
	}

	//=======================================================================
	string cut_series = "( ";
	//-----------------------------------------------------------------------
	unsigned int n_var = output.size();
	std::vector<int> index(n_var);
	std::vector<double> significances;
	std::vector<string> var_names;
	std::vector<string> opt_cut_var_names;
	for (const auto& [key, value] : output) {
		significances.push_back(double(value.second));
		var_names.push_back(string(key));
	}
	TMath::SortItr(significances.begin(), significances.end(), index.begin(), true );

	cout<<"****************************************************"<<endl;
	for(int i=0; i<n_var; i++) {
		unsigned int j = index[i];
		(output_sorted[string(var_names.at(j))]).first.first.first    = (output[var_names.at(j)]).first.first.first  ;
		(output_sorted[string(var_names.at(j))]).first.first.second   = (output[var_names.at(j)]).first.first.second ;
		(output_sorted[string(var_names.at(j))]).first.second.first   = (output[var_names.at(j)]).first.second.first  ;
		(output_sorted[string(var_names.at(j))]).first.second.second  = (output[var_names.at(j)]).first.second.second ;
		(output_sorted[string(var_names.at(j))]).second = (output[var_names.at(j)]).second;
		bool is_variable_for_cut=false;
		for(int icv=0; icv<variables_for_cut.size(); icv++) {
			if(string(variables_for_cut.at(icv)) == string("*")) {
				is_variable_for_cut = true;
				break;
			}
			if(string(var_names.at(j)) == string(variables_for_cut.at(icv))) {
				is_variable_for_cut = true;
			}
		}
		if(!is_variable_for_cut) continue;
		if( ( ((TString(var_name)).Contains("M4l")) || ((TString(var_name)).Contains("Pt4l")) || ((TString(var_name)).Contains("P4l")) || ((TString(var_name)).Contains("deltaR")) ) ) {
			if( (output_sorted[string(var_names.at(j))]).first.first.first   < xmin+2*width_step ) (output_sorted[string(var_names.at(j))]).first.first.first  = -1E7 ;
			if( (output_sorted[string(var_names.at(j))]).first.first.second  > xmax-2*width_step ) (output_sorted[string(var_names.at(j))]).first.first.second = 1E7 ;
		}
		cut_series += string(Form("(%s >= %f && %s <= %f)", var_names.at(j).c_str(), (output[var_names.at(j)]).first.first.first, var_names.at(j).c_str(), (output[var_names.at(j)]).first.first.second));
		if( !(cut_series.empty()) ) cut_series += " && ";
		std::cout << '{' << var_names.at(j) << "}: cuts = [" << (output[var_names.at(j)]).first.first.first << ", " << (output[var_names.at(j)]).first.first.second << "] ; significance = " << (output[var_names.at(j)]).second << endl;
		opt_cut_var_names.push_back(var_names.at(j).c_str());
	}
	//-----------------------------------------------------------------------
	//cut_series += " (1==1) ";
	//for (const auto& [key, value] : output) {
	//	std::cout << '{' << key << "}: cuts = [" << value.first.first.first << ", " << value.first.first.second << "] ; significance = " << value.second << endl;
	//	if( ((TString(key)).Contains("deltaM1")) ) continue;
	//	//if( (value.second > 0.01) && (!((TString(key)).Contains("num")) && !((TString(key)).Contains("SF")) && !((TString(key)).Contains("event")) && !((TString(key)).Contains("lepton_fst")) && !((TString(key)).Contains("deltaM1"))) ) 
	//	if( ((TString(key)).Contains("M4l")) || ((TString(key)).Contains("Pt4l")) || ((TString(key)).Contains("deltaM")) || ((TString(key)).Contains("met")) || ((TString(key)).Contains("Mrecoil")) || ((TString(key)).Contains("deltaR_ll1")) || ((TString(key)).Contains("deltaR_ll2")) || ((TString(key)).Contains("Mll1")) || ((TString(key)).Contains("Mll2")) || ((TString(key)).Contains("Ptll1")) || ((TString(key)).Contains("Ptll2")) ) 
	//	{
	//		cut_series += string(Form(" && (%s >= %f && %s <= %f)", key.c_str(), value.first.first.first, key.c_str(), value.first.first.second));
	//	}
	//}
	//-----------------------------------------------------------------------
	cut_series += " ( " + optcut + " ) ";
	cut_series += " )";
	//cout<<"****************************************************"<<endl;
	//cout<<"cut_series is :"<<endl<<endl;
	//cout<<cut_series<<endl<<endl;
	//cout<<"****************************************************"<<endl;

	//=======================================================================

	std::map<string, std::pair<std::pair<std::pair<double, double>, std::pair<double, double>>, double>> series_output;
	std::pair<double, double> series_cut_limits;
	std::pair<double, double> series_n_sig_and_n_bkg;
	//string opt_cut_series = " (1==1) ";
	string opt_cut_series = " ( " + optcut + " ) ";
	for(int k=0; k<opt_cut_var_names.size(); k++){
		std::vector<TH1D*> hist_for_cuts;
		xmin = 0;
		xmax = 100;
		xrange = xmax - xmin;
		width_step = xrange/ndivisions;
		for(int fid=0; fid<nfiles; fid++){ 
			if( sample_types.at(fid) != 2 ) continue;
			TChain chin("tree");
			chin.Add(filesin.at(fid));
			chin.SetBranchStatus("*", 1);
			//chin.SetBranchStatus("SF", 1);
			//chin.SetBranchStatus("event_weight", 1);
			chin.Draw(Form("%s>>hist_for_cut%d", (opt_cut_var_names.at(k)).c_str(), k), Form("10000000*SF*(%s > -99)*(weights[20])", (opt_cut_var_names.at(k)).c_str()), "goff");
			//double xmean = ((TH1F*)(gDirectory->Get(Form("th_temp%d", k))))->GetMean();
			//double xstaddev = ((TH1F*)(gDirectory->Get(Form("th_temp%d", k))))->GetStdDev();
			//xmin = chin.GetMinimum((opt_cut_var_names.at(k)).c_str());
			//xmax = chin.GetMaximum((opt_cut_var_names.at(k)).c_str());
			if( ((TH1F*)(gDirectory->Get(Form("hist_for_cut%d", k)))) != 0 ) {
				xmin = ((TH1F*)(gDirectory->Get(Form("hist_for_cut%d", k))))->GetBinLowEdge(1);
				xmax = ((TH1F*)(gDirectory->Get(Form("hist_for_cut%d", k))))->GetBinLowEdge( ((TH1F*)(gDirectory->Get(Form("hist_for_cut%d", k))))->GetNcells() );
				xrange = xmax - xmin;
				width_step = xrange/ndivisions;
				xmin -= 2.*width_step;
				xmax += 2.*width_step;
			}
			else{
				xmin = 0;
				xmax = 100;
				xrange = xmax - xmin;
				width_step = xrange/ndivisions;
				xmin -= 2.0*width_step;
				xmax += 2.0*width_step;
			}
		}
		for(int fid=0; fid<nfiles; fid++){ 
			if( sample_types.at(fid) < 0 || sample_types.at(fid) > 2 ) continue;
			TChain chin("tree");
			chin.Add(filesin.at(fid));
			chin.SetBranchStatus("*", 1);
			//chin.SetBranchStatus("SF", 1);
			//chin.SetBranchStatus("event_weight", 1);
			hist_for_cuts.push_back(new TH1D(Form("opt_hist_for_cuts%d_%s",fid,(opt_cut_var_names.at(k)).c_str()), Form("%s",(opt_cut_var_names.at(k)).c_str()), ndivisions+4, xmin, xmax));
			if( sample_types.at(fid) != 2 ) {
				if( ((TH1F*)(gDirectory->Get(Form("hist_for_cut%d", k)))) != 0 ) {
					chin.Draw(Form("%s>>%s", (opt_cut_var_names.at(k)).c_str(), Form("opt_hist_for_cuts%d_%s",fid,(opt_cut_var_names.at(k)).c_str())), Form("10000000*SF*(%s > -99)*(event_weight)*(%s)",(opt_cut_var_names.at(k)).c_str(),opt_cut_series.c_str()), "goff");
				}
				else {
					((TH1F*)(gDirectory->Get(Form("opt_hist_for_cuts%d_%s",fid,(opt_cut_var_names.at(k)).c_str()))))->Fill(1E7);
				}
			}
			if( sample_types.at(fid) == 2 ) {
				if( ((TH1F*)(gDirectory->Get(Form("hist_for_cut%d", k)))) != 0 ) {
					chin.Draw(Form("%s>>%s", (opt_cut_var_names.at(k)).c_str(), Form("opt_hist_for_cuts%d_%s",fid,(opt_cut_var_names.at(k)).c_str())), Form("10000000*SF*(%s > -99)*(weights[20])*(%s)",(opt_cut_var_names.at(k)).c_str(),opt_cut_series.c_str()), "goff");
				}
				else {
					((TH1F*)(gDirectory->Get(Form("opt_hist_for_cuts%d_%s",fid,(opt_cut_var_names.at(k)).c_str()))))->Fill(1E7);
				}
			}
		}
		significance = 0;
		n_sig_for_significance = 0;
		n_bkg_for_significance = 0;
		var_name = string( opt_cut_var_names.at(k) );
		for(int nstart=3; nstart<=(ndivisions+2); nstart++) {
			for(int nwidth=1; (nstart+nwidth-1)<=(ndivisions+4); nwidth++) {
				if( ( ((TString(var_name)).Contains("M4l")) || ((TString(var_name)).Contains("Pt4l")) || ((TString(var_name)).Contains("met")) || ((TString(var_name)).Contains("deltaR_ll1")) || ((TString(var_name)).Contains("deltaR_ll2")) ) && !( ((TString(var_name)).Contains("met_byhand")) || ((TString(var_name)).Contains("meta")) ) ) {
					nwidth=(ndivisions+4+1-nstart);
				}
				n_sig = 0;
				n_bkg = 0;
				for(int fid=0; fid<nfiles; fid++) { 
					if( sample_types.at(fid) == 2 ) { 
						n_sig += ( (hist_for_cuts.at(fid))->Integral(nstart,(nstart+nwidth-1)) );
					}
					if( sample_types.at(fid) == 0 || sample_types.at(fid) == 1 ) { 
						n_bkg += ( (hist_for_cuts.at(fid))->Integral(nstart,(nstart+nwidth-1)) );
					}
				}
				if( (n_sig>0 && n_bkg>0) && ((n_sig/sqrt(n_bkg)) > significance) ) { 
					significance = (n_sig/sqrt(n_bkg));
					n_sig_for_significance = n_sig;
					n_bkg_for_significance = n_bkg;
					integral_left = (xmin+((nstart-1)*width_step));
					integral_right = (xmin+((nstart-1)*width_step)+(nwidth*width_step));
				}
				if( ( ((TString(var_name)).Contains("M4l")) || ((TString(var_name)).Contains("Pt4l")) || ((TString(var_name)).Contains("met")) || ((TString(var_name)).Contains("deltaR_ll1")) || ((TString(var_name)).Contains("deltaR_ll2")) ) && !( ((TString(var_name)).Contains("met_byhand")) || ((TString(var_name)).Contains("meta")) ) ) {
					break;
				}
			}
		}
		if( integral_left < xmin+2*width_step )  integral_left = -1E7 ;
		if( integral_right > xmax-2*width_step ) integral_right = 1E7 ;
		series_cut_limits.first  = integral_left;
		series_cut_limits.second = integral_right;
		series_n_sig_and_n_bkg.first  = n_sig_for_significance;
		series_n_sig_and_n_bkg.second = n_bkg_for_significance;
		(series_output[var_name]).first.first  = cut_limits;
		(series_output[var_name]).first.second = n_sig_and_n_bkg;
		(series_output[var_name]).second = significance;
		opt_cut_series += string(Form(" && (%s >= %f && %s <= %f)", var_name.c_str(), integral_left, var_name.c_str(), integral_right));
	}




	cout<<"****************************************************"<<endl;
	cout<<"optimized cut_series is :"<<endl<<endl;
	cout<<opt_cut_series<<endl<<endl;
	cout<<"****************************************************"<<endl;

	//optcut = cut_series;
	optcut = opt_cut_series;
	return output;
}












double ComputeSignificance(std::vector<TString> filesin, std::vector<int> sample_types, std::vector<int> sig_types, std::vector<int> bkg_types, string cut="") {

	gStyle->SetOptStat(0000);

	int nfiles = filesin.size();
	int signfiles=0;
	for(int fi=0; fi<nfiles; fi++){
		if(sample_types.at(fi)==2) signfiles++;
	}
	int bkgnfiles=0;
	for(int fi=0; fi<nfiles; fi++){
		if(sample_types.at(fi)==0 || sample_types.at(fi)==1) bkgnfiles++;
	}

	bool is_sig=false;
	bool is_bkg=false;
	double significance=0;
	double n_sig=0;
	double n_bkg=0;
	double total_n_bkg=0;
	std::map<string, double> samples_yield_after_cut;
	std::map<string, double> sm_background_samples_yield_after_cut;
	sm_background_samples_yield_after_cut["SMsignal"] = 0;
	sm_background_samples_yield_after_cut["mumuTomumuX"] = 0;
	sm_background_samples_yield_after_cut["mumuTovmvmX"] = 0;
	sm_background_samples_yield_after_cut["mumuTottX"] = 0;
	sm_background_samples_yield_after_cut["mumuTomultibosons"] = 0;
	cout<<cut.c_str()<<endl;
	for(int fid=0; fid<nfiles; fid++) { 
		n_sig=0;
		n_bkg=0;
		is_sig=false;
		is_bkg=false;
		for(int k=0; k<sig_types.size(); k++) {
			if(sample_types.at(fid)==sig_types.at(k)) {
				is_sig=true;
			}
		}
		for(int k=0; k<bkg_types.size(); k++) {
			if(sample_types.at(fid)==bkg_types.at(k)) {
				is_bkg=true;
			}
		}
		if( !is_sig && !is_bkg ) continue;
		TChain chin("tree");
		chin.Add(filesin.at(fid));
		chin.SetBranchStatus("*", 1);
		cout<<filesin.at(fid)<<" (after cut) = "<<chin.GetEntries(cut.c_str())<<endl;
		string process((filesin.at(fid)).Data());
		//process.erase(process.find("_delphes_preselected.root",0), string("_delphes_preselected.root").size());
		process.erase(process.find(".root",0), string(".root").size());
		while(process.find("/")>0 && process.find("/")<(process.size())) {
			 process.erase(0, (process.find("/") + 1));
		}
		//if(sample_types.at(fi)==0)  
		//if((TString(process)).Contains("mumuTozzzTo")) 
		if(is_sig)
		{
			chin.Draw(Form("%s>>event_th_file%d", "event", fid), Form("10000000*SF*(event_weight)*(%s)", cut.c_str()), "goff");
			if( ((TH1F*)(gDirectory->Get(Form("event_th_file%d", fid)))) != NULL ) {
				n_sig += ((TH1F*)(gDirectory->Get(Form("event_th_file%d", fid))))->Integral();
			}
			cout<<"n_sig = "<<n_sig<<endl;
			samples_yield_after_cut[(process+"{weights[20]@(Fx=-10E-12)}")] = n_sig;
		}
		//if(sample_types.at(fi)==1)  
		//if( !((TString(process)).Contains("mumuTozzzTo")) ) 
		if(is_bkg)
		{
			chin.Draw(Form("%s>>event_th_file%d", "event", fid), Form("10000000*SF*(event_weight)*(%s)", cut.c_str()), "goff");
			if( ((TH1F*)(gDirectory->Get(Form("event_th_file%d", fid)))) != NULL ) {
				n_bkg = ((TH1F*)(gDirectory->Get(Form("event_th_file%d", fid))))->Integral();
				samples_yield_after_cut[process] = n_bkg;
				total_n_bkg += n_bkg;
				cout<<"total_n_bkg = "<<total_n_bkg<<endl;
				{
					if(sample_types.at(fid)==0) {
						sm_background_samples_yield_after_cut["SMsignal"] += n_bkg;
					}
					else if((TString(process)).Contains("mumuTomumu")) {
						sm_background_samples_yield_after_cut["mumuTomumuX"] += n_bkg;
					}
					else if((TString(process)).Contains("mumuTovmvm")) {
						sm_background_samples_yield_after_cut["mumuTovmvmX"] += n_bkg;
					}
					else if((TString(process)).Contains("tt")) {
						sm_background_samples_yield_after_cut["mumuTottX"] += n_bkg;
					}
					else if((TString(process)).Contains("mumuTow") || (TString(process)).Contains("mumuToz") || (TString(process)).Contains("mumuToh") ) {
						sm_background_samples_yield_after_cut["mumuTomultibosons"] += n_bkg;
					}
				}
			}
		}
	}
	cout<<"Summary of yields: ==============================================="<<endl;
	int n_sample = samples_yield_after_cut.size();
	std::vector<int> eveindex(n_sample);
	std::vector<double> events_yield;
	std::vector<string> sample_name;
	for (const auto& [key, value] : samples_yield_after_cut) {
		events_yield.push_back(double(value));
		sample_name.push_back(string(key));
	}
	TMath::SortItr(events_yield.begin(), events_yield.end(), eveindex.begin(), true );

	//----------------------- Print yields summary with sorted order
	for(int k=0; k<events_yield.size(); k++){
		if((TString(sample_name.at(eveindex.at(k)))).Contains("Fx=")) {
			cout<<sample_name.at(eveindex.at(k))<<"+++ : signal yields after cut = "<<events_yield.at(eveindex.at(k))<<endl;
		}
		if( !((TString(sample_name.at(eveindex.at(k)))).Contains("Fx=")) ) {
			cout<<sample_name.at(eveindex.at(k))<<"--- : background yields after cut = "<<events_yield.at(eveindex.at(k))<<endl;
		}
	}
	//cout<<"=================================================================="<<endl;
	//cout<<"Summary of background yields: ===================================="<<endl;
	//for (const auto& [key, value] : sm_background_samples_yield_after_cut) {
	//	std::cout << '{' << key << "}: yields = [" << value << "] ;" << endl;
	//}
	cout<<"=================================================================="<<endl;
	if( (n_sig>0 && total_n_bkg>0) ) { 
		//significance = (n_sig/sqrt(total_n_bkg));
		significance = (sqrt(2.0*((n_sig+n_bkg)*(TMath::Log(1+(n_sig/n_bkg)))-n_sig)));
	}
	cout<<"significance = "<<significance<<endl;
	cout<<"=================================================================="<<endl;

	return significance;
}





using namespace RooStats;

class tools_for_aqgc_constraints {

	public:
		void clear_objects(){
			assert((paranames.size == 1) || (paranames.size == 2));
			if(!aqgc_para_names.empty())aqgc_para_names.clear();
			else aqgc_para_names = std::vector<string>();
			if(!fvectors_histsamples.empty())fvectors_histsamples.clear();
			else fvectors_histsamples = std::vector<TH1D*>();
			if(!fvectors_histsamplestypes.empty())fvectors_histsamplestypes.clear();
			else fvectors_histsamplestypes = std::vector<int>();
			if(!fvectors_histsamples2.empty())fvectors_histsamples2.clear();
			else fvectors_histsamples2 = std::vector<TH1D*>();
			if(!fvectors_histsamplestypes2.empty())fvectors_histsamplestypes2.clear();
			else fvectors_histsamplestypes2 = std::vector<int>();
			sm_signal_yields = 0;
			func_RatioToSM_VS_aqgc_1D_minX = 0;
			sm_background_yields = 0;
			func_RatioToSM_VS_aqgc_1D = 0;
			formula_func_RatioToSM_VS_aqgc_1D = "";
			func_RatioToSM_VS_aqgc_2D = 0;
			formula_func_RatioToSM_VS_aqgc_2D = "";
			graph_RatioToSM_VS_aqgc_1D = 0;
			graph_RatioToSM_VS_aqgc_2D = 0;
			RatioToSM_VS_aqgc_fitresult = 0;
		}

		void set_confidence_level(double CLv=0.95){
			confidence_level = CLv;
		}

		void set_Models_ntoys(int alt_ntoys=1000, int null_ntoys=500){
			altModel_ntoys = alt_ntoys;
			nullModel_ntoys = null_ntoys;
		}

		void set_ts_type(int tstype=3){
			ts_type = tstype;
		}

		void set_fit_for_ts(bool doFit=false){
			fit_for_ts = doFit;
		}

		void set_scan_range(double left=-50, double right=50){
			scan_range.first = left;
			scan_range.second = right;
		}

		tools_for_aqgc_constraints(std::vector<string> paranames){
			assert(paranames.size == 1);
			clear_objects();
			for(int i=0; i<paranames.size(); i++){
				aqgc_para_names.push_back(paranames[i]);
			}
			func_is_interpolated = false;
			if(func_RatioToSM_VS_aqgc_2D!=0 && func_RatioToSM_VS_aqgc_2D->IsValid()) func_is_interpolated = true;
			fvectors_histintegrals = std::vector<double>();
			fvectors_aqgcparavalues_x = std::vector<double>();
			fvectors_aqgcparavalues_y = std::vector<double>();
			wspace = new RooWorkspace("wspace");
			aqgc_modelconfig = 0; //new ModelConfig("aqgc_modelconfig", wspace);
			sm_modelconfig = 0;
			set_confidence_level();
			set_Models_ntoys();
			set_ts_type();
			set_fit_for_ts();
			set_scan_range();
			cout<<" Tools for aQGC  constraints are initiated !"<<endl;
		}

		tools_for_aqgc_constraints(std::vector<string> paranames, std::vector<TH1D*> histsamples,  std::vector<int>histsamplestypes){
			assert(paranames.size == 1);
			clear_objects();
			for(int i=0; i<paranames.size(); i++){
				aqgc_para_names.push_back(paranames[i]);
			}
			aqgc_para_cons = std::vector<std::pair<double,double>>();
			for(int i=0; i<histsamples.size(); i++){
				fvectors_histsamples.push_back(histsamples[i]);
			}
			for(int i=0; i<histsamplestypes.size(); i++){
				fvectors_histsamplestypes.push_back(histsamplestypes[i]);
			}
			func_is_interpolated = false;
			if(func_RatioToSM_VS_aqgc_2D!=0 && func_RatioToSM_VS_aqgc_2D->IsValid()) func_is_interpolated = true;
			fvectors_histintegrals = std::vector<double>();
			fvectors_aqgcparavalues_x = std::vector<double>();
			fvectors_aqgcparavalues_y = std::vector<double>();
			wspace = new RooWorkspace("wspace");
			aqgc_modelconfig = 0; //new ModelConfig("aqgc_modelconfig", wspace);
			sm_modelconfig = 0;
			set_confidence_level();
			set_Models_ntoys();
			set_ts_type();
			set_fit_for_ts();
			set_scan_range();
			cout<<" Tools for aQGC  constraints are initiated !"<<endl;
		}

		tools_for_aqgc_constraints(std::vector<string> paranames, std::vector<TH1D*> histsamples,  std::vector<int>histsamplestypes, std::vector<TH1D*> histsamples2,  std::vector<int>histsamplestypes2){
			assert(paranames.size == 2);
			clear_objects();
			for(int i=0; i<paranames.size(); i++){
				aqgc_para_names.push_back(paranames[i]);
			}
			aqgc_para_cons = std::vector<std::pair<double,double>>();
			for(int i=0; i<histsamples.size(); i++){
				fvectors_histsamples.push_back(histsamples[i]);
			}
			for(int i=0; i<histsamplestypes.size(); i++){
				fvectors_histsamplestypes.push_back(histsamplestypes[i]);
			}
			for(int i=0; i<histsamples2.size(); i++){
				fvectors_histsamples2.push_back(histsamples2[i]);
			}
			for(int i=0; i<histsamplestypes2.size(); i++){
				fvectors_histsamplestypes2.push_back(histsamplestypes2[i]);
			}
			func_is_interpolated = false;
			if(func_RatioToSM_VS_aqgc_2D!=0 && func_RatioToSM_VS_aqgc_2D->IsValid()) func_is_interpolated = true;
			fvectors_histintegrals = std::vector<double>();
			fvectors_aqgcparavalues_x = std::vector<double>();
			fvectors_aqgcparavalues_y = std::vector<double>();
			wspace = new RooWorkspace("wspace");
			aqgc_modelconfig = 0; //new ModelConfig("aqgc_modelconfig", wspace);
			sm_modelconfig = 0;
			set_confidence_level();
			set_Models_ntoys();
			set_ts_type();
			set_fit_for_ts();
			set_scan_range();
			cout<<" Tools for aQGC  constraints are initiated !"<<endl;
		}

		~tools_for_aqgc_constraints(){};

		void add_samples(std::vector<TH1D*> histsamples, std::vector<int>histsamplestypes){
			for(int i=0; i<histsamples.size(); i++){
				fvectors_histsamples.push_back(histsamples[i]);
			}
			for(int i=0; i<histsamplestypes.size(); i++){
				fvectors_histsamplestypes.push_back(histsamplestypes[i]);
			}
		}

		void calculate_aqgcparavalues_vs_histintegrals(){
			assert((aqgc_para_names.size() == 1) || (aqgc_para_names.size() == 2));
			if(aqgc_para_names.size() == 2) { cout<<"!!! Should construct 2D function, code is not prepared yet !!!"<<endl; return; }
			string value;
			double dvalue;
			std::vector<double> temp_fvectors_histintegrals;
			std::vector<double> temp_fvectors_aqgcparavalues;
			for (int i=0; i<fvectors_histsamples.size(); i++) {
				if((TString(fvectors_histsamples[i]->GetName())).Contains("th_aqgcRW") && (fvectors_histsamplestypes[i]==2)) {
					temp_fvectors_histintegrals.push_back((fvectors_histsamples[i]->Integral(1,fvectors_histsamples[i]->GetNcells()-2)));
					value = fvectors_histsamples[i]->GetName();
					value.erase(0, value.find(Form("="),0)+1);
					//value.erase(value.find(Form("e-12"),0), value.size()-1-value.find(Form("e-12")));
					value.erase(value.find(Form("e-12"),0));
					dvalue = atof(value.c_str());
					temp_fvectors_aqgcparavalues.push_back(dvalue);
				}

			}
			std::vector<int> index(temp_fvectors_aqgcparavalues.size());
			TMath::SortItr(temp_fvectors_aqgcparavalues.begin(), temp_fvectors_aqgcparavalues.end(), index.begin(), false );
			for(int i=0; i<index.size(); i++){
				fvectors_histintegrals.push_back(temp_fvectors_histintegrals[index[i]]);
				fvectors_aqgcparavalues_x.push_back(temp_fvectors_aqgcparavalues[index[i]]);
			}
			cout<<" vectors of aqgcparavalues_vs_histintegrals is prepared !"<<endl;
		}


		void construct_func_RatioToSM_VS_aqgc_1D(){
			assert(fvectors_aqgcparavalues_x.size() > 2);
			graph_RatioToSM_VS_aqgc_1D = new TGraph(fvectors_aqgcparavalues_x.size());
			//graph_RatioToSM_VS_aqgc_1D->SetTitle("graph_RatioToSM_VS_aqgc_1D");
			for(int i=0; i<fvectors_aqgcparavalues_x.size(); i++){
				graph_RatioToSM_VS_aqgc_1D->SetPoint(i, fvectors_aqgcparavalues_x[i], fvectors_histintegrals[i]);
			}
			if(!func_is_interpolated){
				func_RatioToSM_VS_aqgc_1D = new TF1("func_RatioToSM_VS_aqgc_1D","[0]+[1]*x+[2]*x*x",-150,150);
				RatioToSM_VS_aqgc_fitresult = graph_RatioToSM_VS_aqgc_1D->Fit(func_RatioToSM_VS_aqgc_1D, "SRE");
				assert(RatioToSM_VS_aqgc_fitresult->IsValid());
			}
			if(func_is_interpolated){
				cout<<"No available method founded!"<<endl;
				//func_RatioToSM_VS_aqgc_1D = new TF1("func_RatioToSM_VS_aqgc_1D",[&](double*x){ return graph_RatioToSM_VS_aqgc_1D->Eval(x[0]); }, graph_RatioToSM_VS_aqgc_1D->GetXaxis()->GetXmin(), graph_RatioToSM_VS_aqgc_1D->GetXaxis()->GetXmax(), 0);
				assert(func_RatioToSM_VS_aqgc_1D);
				assert(RatioToSM_VS_aqgc_fitresult->IsValid());
			}
			//sm_signal_yields = func_RatioToSM_VS_aqgc_1D->Eval(0); 
			//sm_signal_yields = func_RatioToSM_VS_aqgc_1D->GetMinimum(graph_RatioToSM_VS_aqgc_1D->GetXaxis()->GetXmin(), graph_RatioToSM_VS_aqgc_1D->GetXaxis()->GetXmax());

			func_RatioToSM_VS_aqgc_1D_minX = (-0.5*(func_RatioToSM_VS_aqgc_1D->GetParameter(1))/(func_RatioToSM_VS_aqgc_1D->GetParameter(2))); sm_signal_yields = func_RatioToSM_VS_aqgc_1D->Eval(func_RatioToSM_VS_aqgc_1D_minX);
			//for(int is=0; is<fvectors_histsamplestypes.size(); is++){
			//	if(fvectors_histsamplestypes.at(is)==0) sm_signal_yields += (fvectors_histsamples[is])->Integral(1, (fvectors_histsamples[is])->GetNcells()-2);
			//}
			for(int is=0; is<fvectors_histsamplestypes.size(); is++){
				if(fvectors_histsamplestypes.at(is)==1) sm_background_yields += (fvectors_histsamples[is])->Integral(1, (fvectors_histsamples[is])->GetNcells()-2);
			}
			formula_func_RatioToSM_VS_aqgc_1D = Form("( %e + %e * poi1 + %e * poi1 * poi1 )", func_RatioToSM_VS_aqgc_1D->GetParameter(0), func_RatioToSM_VS_aqgc_1D->GetParameter(1), func_RatioToSM_VS_aqgc_1D->GetParameter(2)); // formula_func_RatioToSM_VS_aqgc_1D = func_RatioToSM_VS_aqgc_1D->GetFormula();
			cout<<" func_RatioToSM_VS_aqgc_1D is builded !"<<endl;
		}


		void construct_func_RatioToSM_VS_aqgc_2D(){
			assert(fvectors_aqgcparavalues_x.size() > 2);
			// TO BE ACCOMPLISHED ...
		}


		void draw_RatioToSM_VS_aqgc_1D_plot(){
			TCanvas tempc("temp_draw_RatioToSM_VS_aqgc_1D_plot","");
			tempc.SetRightMargin(0.25);
			tempc.cd();
			tempc.SetLeftMargin(0.16);
			tempc.SetBottomMargin(0.13);
			graph_RatioToSM_VS_aqgc_1D->SetMarkerStyle(20);
			graph_RatioToSM_VS_aqgc_1D->SetMarkerSize(1.5);
			graph_RatioToSM_VS_aqgc_1D->SetMarkerColor(1);
			graph_RatioToSM_VS_aqgc_1D->SetLineColor(1);
			//graph_RatioToSM_VS_aqgc_1D->SetTitle(Form("; %s; #frac{#sigma^{aQGC}_{sig}}{#sigma^{SM}_{sig}}",aqgc_para_names[0].c_str()));
			graph_RatioToSM_VS_aqgc_1D->SetTitle(Form("; %s; N^{aQGC}_{sig}",aqgc_para_names[0].c_str()));
			graph_RatioToSM_VS_aqgc_1D->GetXaxis()->SetTitleSize(0.05);
			graph_RatioToSM_VS_aqgc_1D->GetXaxis()->SetLabelSize(0.05);
			graph_RatioToSM_VS_aqgc_1D->GetYaxis()->SetTitleSize(0.05);
			graph_RatioToSM_VS_aqgc_1D->GetYaxis()->SetLabelSize(0.05);
			graph_RatioToSM_VS_aqgc_1D->Draw("AP");
			func_RatioToSM_VS_aqgc_1D->SetLineColor(4);
			func_RatioToSM_VS_aqgc_1D->Draw("SAME");
			TLegend lg(0.75,0.80,0.99,0.99);
			lg.AddEntry(graph_RatioToSM_VS_aqgc_1D, "Data points",   "p");
			lg.AddEntry(func_RatioToSM_VS_aqgc_1D,  "Quadratic fit", "l");
			lg.AddEntry(func_RatioToSM_VS_aqgc_1D, Form("C0 = %.3e",func_RatioToSM_VS_aqgc_1D->GetParameter(0)), "p");
			lg.AddEntry(func_RatioToSM_VS_aqgc_1D, Form("C1 = %.3e",func_RatioToSM_VS_aqgc_1D->GetParameter(1)), "p");
			lg.AddEntry(func_RatioToSM_VS_aqgc_1D, Form("C2 = %.3e",func_RatioToSM_VS_aqgc_1D->GetParameter(2)), "p");
			lg.Draw();
			tempc.SaveAs(Form("RatioToSM_VS_aqgc_1D_plot_%s.png",aqgc_para_names[0].c_str()));
			cout<<" RatioToSM_VS_aqgc_1D plot is drawn !"<<endl;
		}

		double eval_sig_expect_1D(double x){
			double sig_expect = 0;
			if(!func_is_interpolated){
				assert(RatioToSM_VS_aqgc_fitresult->IsValid());
				sig_expect = func_RatioToSM_VS_aqgc_1D->Eval(x);
			}
			else{
				sig_expect = graph_RatioToSM_VS_aqgc_1D->Eval(x);
			}
			return sig_expect;
		}

		double obtain_sm_signal_yields(){
			return sm_signal_yields;
		}

		double obtain_func_RatioToSM_VS_aqgc_1D_minX(){
			return func_RatioToSM_VS_aqgc_1D_minX;
		}

		double obtain_sm_background_yields(){
			return sm_background_yields;
		}

		TF1* obtain_func_1D() {
			assert(func_RatioToSM_VS_aqgc_1D);
			return func_RatioToSM_VS_aqgc_1D;
		}

		TF1* obtain_func_2D() {
			assert(func_RatioToSM_VS_aqgc_2D);
			return func_RatioToSM_VS_aqgc_2D;
		}

		double obtain_significance(double par_value) {
			double aqgc_sig = (func_RatioToSM_VS_aqgc_1D->Eval(par_value));
			double sm_bkg = sm_background_yields;
			double significance = sqrt(2.0*((aqgc_sig+sm_bkg)*(TMath::Log(1+(aqgc_sig/sm_bkg)))-aqgc_sig));
			return significance;
		}

		void print_par_value_of_specific_significance(double target_significance) {
			ofstream output(Form("./HypoTestInvLog_%s.log", aqgc_para_names[0].c_str()), ios::app);
			cout << " +++ --------- Print AQGC "<< aqgc_para_names[0] << " values for specific significance at [significance = " << target_significance<<"] --------- +++ " << endl;
			output << " +++ --------- Print AQGC "<< aqgc_para_names[0] << " values for specific significance at [significance = " << target_significance<<"] --------- +++ " << endl;
			double temp_significance = 0;
			double nsteps = 100;
			double step = (20 - (-20))/nsteps;
			for(int istep=0; istep<nsteps; istep++) {
				if( ((obtain_significance(-20+(istep*step)) - target_significance) * (obtain_significance(-20+((istep+1)*step)) - target_significance)) <= 0 ) {
					cout << " --- AQGC "<< aqgc_para_names[0] << " value is " << (-20+((istep+1)*step)) << " --- " << endl;
					output << " --- AQGC "<< aqgc_para_names[0] << " value is " << (-20+((istep+1)*step)) << " --- " << endl;
				}
			}
			cout << " --- --------- Print ended --------- --- " <<endl;
			output << " --- --------- Print ended --------- --- " <<endl;
			output.close();
		}
				


		//------------------------
		// Manners used for hypothesis test
		//------------------------
		void setup_aqgc_modelconfigs(double exp_s=5.0, double exp_b=15.0, RooDataSet* datasample=0) {
			wspace->factory( Form("obs[10.0,0.0,300.0]") ); // observable: counting number
			//wspace->factory( Form("poi1[0.0,%e,%e]",-100.0,100.0) ); // poi: the first POI
			wspace->factory( Form("poi1[%e,%e,%e]",((scan_range.first+scan_range.second)/2.0),(scan_range.first),(scan_range.second)) ); // poi: the first POI
			//wspace->factory( Form("poi1[%e,%e,%e]",((scan_range.first+scan_range.second)/2.0),(scan_range.first-0.1*(scan_range.second-scan_range.first)),(scan_range.second+0.1*(scan_range.second-scan_range.first))) ); // poi: the first POI
			// --- terms for background uncertainty
			//wspace->factory( Form("expectB[30.0]") ); // nuisance: expected background
			//wspace->factory( Form("ratioBkg[1.0,0.0,3.0]") ); // nuisance: ratio on expected background to manifest its uncertainty
			//wspace->factory( Form("globalBkgObsNum[10.0,0.0,300.0]") ); // global observable for nuisisance: ratioBkg
			//wspace->factory( Form("prod::globalBkg(globalBkgExpNum[10.0],ratioBkg)") ); // expected global observable for nuisisance: ratioBkg
			//wspace->factory( Form("Gaussian::bkgConstraint(globalBkgObsNum,globalBkg,10.0)") );  // nuisance gaussian constraint: background --- the sigma 10.0 is set somehow arbitially
			//wspace->factory( Form("prod::Bkg(ratioBkg,expectB") ); 

			// --- 
			wspace->factory( Form("Bkg[%e]",exp_b) );
			wspace->factory( Form("a0[%e]",func_RatioToSM_VS_aqgc_1D->GetParameter(0)) );
			wspace->factory( Form("a1[%e]",func_RatioToSM_VS_aqgc_1D->GetParameter(1)) );
			wspace->factory( Form("a2[%e]",func_RatioToSM_VS_aqgc_1D->GetParameter(2)) );
			//wspace->factory( Form("expectS[%e]",exp_s) );
			//wspace->factory( Form("RooPolyVar::ratioSig(poi1,{a0,a1,a2})") );
			//wspace->factory( Form("prod::Sig(ratioSig,expectS]") );
			wspace->factory( Form("RooPolyVar::Sig(poi1,{a0,a1,a2})") ); // signal counting expected in poisson distribution
			wspace->factory( Form("RooPoisson::countingModel(obs,sum(Sig,Bkg))") ); // counting model
			wspace->Print();

			RooAbsPdf  *countingModel = (RooAbsPdf *) wspace->pdf("countingModel"); 
			RooRealVar *obs           = (RooRealVar*) wspace->var("obs"); 
			RooRealVar *poi1          = (RooRealVar*) wspace->var("poi1"); 
			//RooRealVar *ratioBkg      = (RooRealVar*) wspace->var("ratioBkg");
			RooPolyVar *s             = (RooPolyVar*) wspace->function("Sig");
			RooRealVar *b             = (RooRealVar*) wspace->var("Bkg");
			b->setConstant();

			RooArgSet paramOfInterest(*poi1);
			/// --- modelconfig for aQGC (null hypothesis)
			aqgc_modelconfig = new ModelConfig(wspace);
			aqgc_modelconfig->SetName("aqgc_modelconfig");
			aqgc_modelconfig->SetPdf(*countingModel);
			aqgc_modelconfig->SetParametersOfInterest(paramOfInterest);
			aqgc_modelconfig->SetObservables(*obs);
			aqgc_modelconfig->SetSnapshot(*(aqgc_modelconfig->GetParametersOfInterest()));
			//RooArgSet constrainedParams;
			//constrainedParams.add(*ratioBkg);
			//aqgc_modelconfig->SetNuisanceParameters(constrainedParams);
			//aqgc_modelconfig->SetGlobalObservables(RooArgSet(*globalBkgObsNum));
			/// ---
			/// --- modelconfig for sm (alternative hypothesis)
			sm_modelconfig = (ModelConfig*)(aqgc_modelconfig->Clone());
			sm_modelconfig->SetName("sm_modelconfig");      
			RooRealVar* testpoi = (RooRealVar*)(aqgc_modelconfig->GetParametersOfInterest()->first()); // for case of 1 POI
			double oldval = testpoi->getVal();
			testpoi->setVal(0); // Not sure if this should be force to 0
			//testpoi->setVal(func_RatioToSM_VS_aqgc_1D_minX);
			sm_modelconfig->SetSnapshot(RooArgSet(*testpoi));
			testpoi->setVal(oldval);
			/// ---
			wspace->import(*aqgc_modelconfig);
			wspace->import(*sm_modelconfig);
			if(datasample==0 || datasample->numEntries()<1) {
				TRandom r(0);
				//obs->setVal(r.Poisson(sm_signal_yields+sm_background_yields));
				obs->setVal((sm_signal_yields+sm_background_yields));
				datasample = new RooDataSet("exampleData", "exampleData", *obs);
				datasample->add(*obs);
				datasample->Print("V");
			}
			wspace->import(*datasample);
			//wspace->SetName("wspace");
			//wspace->writeToFile(Form("wspace_%s.root",aqgc_para_names[0].c_str()));
			cout<<" Model config is setup !"<<endl;
		}


		void buildup_hypotest(RooDataSet* data=0) {
			ofstream outputlog(Form("./HypoTestInvLog_%s.log", aqgc_para_names[0].c_str()), ios::ate|ios::out);
			// A check on data input, pseudo data is generated for zero input data, the data should not affect the expected limited.
			RooRealVar *obs = (RooRealVar*) wspace->var("obs"); 
			RooRealVar *poi1 = (RooRealVar*) wspace->var("poi1"); 
			if(data==0 || data->numEntries()<1) {
				TRandom r(0);
				//obs->setVal(r.Poisson(sm_signal_yields+sm_background_yields));
				obs->setVal((sm_signal_yields+sm_background_yields));
				data = new RooDataSet("exampleData", "exampleData", *obs);
				data->add(*obs);
			}

			// First, set up the SM prediction model (Null hypothesis): clone from the generic aQGC model (alternative hypothesis) with aQGC parameter set to be 0
			assert(aqgc_modelconfig);
			assert(sm_modelconfig);
			// Second, set up the hypothsis test calculator: here we use the frequentist one. And will pass it to hypothesis inversion calculator for upper limit calculation
			FrequentistCalculator fc(*data, *sm_modelconfig, *aqgc_modelconfig);
			fc.SetToys(altModel_ntoys,nullModel_ntoys); 
			fc.UseSameAltToys();
			fc.StoreFitInfo(true);
			// Third, choose the test statistic: here we use the ProfileLikelihoodTestStat that calculates the profile likelihood ratio
			/// --- Test statistic 1: SimpleLikelihoodRatioTestStat
			SimpleLikelihoodRatioTestStat teststat1(*(aqgc_modelconfig->GetPdf()), *(sm_modelconfig->GetPdf()));
			////// ------ null parameters must includes snapshot of poi plus the nuisance values
			if(fit_for_ts && ts_type==1) {
				RooArgSet constrainParams;
				if(aqgc_modelconfig->GetNuisanceParameters()) constrainParams.add(*aqgc_modelconfig->GetNuisanceParameters());
				RooStats::RemoveConstantParameters(&constrainParams);
				RooFitResult* fitres = 0;
				int itenum = 0;
				while( (fitres==0 || fitres->status()!=0) && (itenum<=5) ) {
					itenum+=1;
					fitres = aqgc_modelconfig->GetPdf()->fitTo(*data, RooFit::InitialHesse(false), RooFit::Hesse(false), RooFit::Minimizer("Minuit2", "Migrad"), RooFit::Strategy(1), RooFit::PrintLevel(-1), RooFit::Constrain(constrainParams), RooFit::Save(true), RooFit::Offset(true));
				}
				aqgc_modelconfig->SetSnapshot(*aqgc_modelconfig->GetParametersOfInterest());
				cout<<"The snapshot of the model for aQGC (null hypothesis) is set to the fit result of data!"<<endl;
				outputlog<<"The snapshot of the model for aQGC (null hypothesis) is set to the fit result of data!"<<endl;
			}
			RooArgSet nullParams(*aqgc_modelconfig->GetSnapshot());
			if(aqgc_modelconfig->GetNuisanceParameters()) nullParams.add(*aqgc_modelconfig->GetNuisanceParameters());
			if(aqgc_modelconfig->GetSnapshot()) teststat1.SetNullParameters(nullParams);
			RooArgSet altParams(*sm_modelconfig->GetSnapshot());
			if(sm_modelconfig->GetNuisanceParameters()) altParams.add(*sm_modelconfig->GetNuisanceParameters());
			if(sm_modelconfig->GetSnapshot()) teststat1.SetAltParameters(altParams);
			/// --- Test statistic 2: RatioOfProfiledLikelihoodsTestStat
			RatioOfProfiledLikelihoodsTestStat teststat2(*aqgc_modelconfig->GetPdf(), *sm_modelconfig->GetPdf(), sm_modelconfig->GetSnapshot());
			teststat2.SetSubtractMLE(false);
			/// --- Test statistic 3: ProfileLikelihoodTestStat
			ProfileLikelihoodTestStat teststat3(*(aqgc_modelconfig->GetPdf()));
			teststat3.SetOneSided(false);
			teststat3.SetSigned(false);
			// Fourth, configure the toyMC sampler
			TestStatistic* teststat = 0;
			if(ts_type == 1) teststat = &teststat1;
			if(ts_type == 2) teststat = &teststat2;
			if(ts_type == 3) teststat = &teststat3;
			ToyMCSampler* toymcs = (ToyMCSampler*)fc.GetTestStatSampler();
			toymcs->SetNEventsPerToy(1);
			toymcs->SetTestStatistic(teststat);
			int npoints = 50;
			double poimin = poi1->getMin();
			double poimax = poi1->getMax();
			cout << "Doing a fixed scan in interval: [" << poimin << " , " << poimax << "]" << endl;
			outputlog << "Doing a fixed scan in interval: [" << poimin << " , " << poimax << "]" << endl;
			// Fifth, draw plots and print results
			TStopwatch tw;
			/// For the lower ( < 0 ) aQGC constraint
			double temp_func_RatioToSM_VS_aqgc_1D_minX = func_RatioToSM_VS_aqgc_1D_minX;
			func_RatioToSM_VS_aqgc_1D_minX = 0; // Not sure if this should be force to 0
			cout <<"=== Estimate the lower ( < "<<func_RatioToSM_VS_aqgc_1D_minX<<" ) aQGC constraint ==="<<endl;
			outputlog <<"=== Estimate the lower ( < "<<func_RatioToSM_VS_aqgc_1D_minX<<" ) aQGC constraint ==="<<endl;
			tw.Start();
			HypoTestInverter calc(fc);
			calc.SetConfidenceLevel(confidence_level);
			//calc.UseCLs(false);
			calc.UseCLs(true);
			calc.SetVerbose(true);
			//poi1->setMin(poimin);
			poi1->setMin(scan_range.first);
			poi1->setMax(func_RatioToSM_VS_aqgc_1D_minX);
			//calc.SetFixedScan(npoints,poimin,0.0);
			calc.SetFixedScan(npoints,poi1->getMin(),poi1->getMax());
			HypoTestInverterResult* htresult_lc = calc.GetInterval();
			assert(htresult_lc!=0);
			cout << "Hypothesis test finished!"<<endl;
			cout << "Time consuming: " << endl;
			outputlog << "Hypothesis test finished!"<<endl;
			outputlog << "Time consuming: " << endl;
			tw.Stop();
			tw.Print();
			double upperLimit_lc = htresult_lc->UpperLimit();
			double ulError_lc = htresult_lc->UpperLimitEstimatedError();
			std::cout << "The computed lower constraint is: " << upperLimit_lc << " +/- " << ulError_lc << std::endl;
			std::cout << " Expected lower constraints, using the alternate model : " << std::endl;
			std::cout << " expected lower constraint (median) " << htresult_lc->GetExpectedUpperLimit(0) << std::endl;
			std::cout << " expected lower constraint (-1 sig) " << htresult_lc->GetExpectedUpperLimit(-1) << std::endl;
			std::cout << " expected lower constraint (+1 sig) " << htresult_lc->GetExpectedUpperLimit(1) << std::endl;
			std::cout << " expected lower constraint (-2 sig) " << htresult_lc->GetExpectedUpperLimit(-2) << std::endl;
			std::cout << " expected lower constraint (+2 sig) " << htresult_lc->GetExpectedUpperLimit(2) << std::endl;
			outputlog << "The computed lower constraint is: " << upperLimit_lc << " +/- " << ulError_lc << std::endl;
			outputlog << " Expected lower constraints, using the alternate model : " << std::endl;
			outputlog << " expected lower constraint (median) " << htresult_lc->GetExpectedUpperLimit(0) << std::endl;
			outputlog << " expected lower constraint (-1 sig) " << htresult_lc->GetExpectedUpperLimit(-1) << std::endl;
			outputlog << " expected lower constraint (+1 sig) " << htresult_lc->GetExpectedUpperLimit(1) << std::endl;
			outputlog << " expected lower constraint (-2 sig) " << htresult_lc->GetExpectedUpperLimit(-2) << std::endl;
			outputlog << " expected lower constraint (+2 sig) " << htresult_lc->GetExpectedUpperLimit(2) << std::endl;
			HypoTestInverterPlot *plot_lc = new HypoTestInverterPlot("HTI_Result_Plot_LowerConstraint","",htresult_lc);
			poi1->setMin(poimin);
			poi1->setMax(poimax);
			/// For the upper ( > func_RatioToSM_VS_aqgc_1D_minX ) aQGC constraint
			cout<<"=== Estimate the upper ( > "<<func_RatioToSM_VS_aqgc_1D_minX<<" ) aQGC constraint ==="<<endl;
			outputlog<<"=== Estimate the upper ( > "<<func_RatioToSM_VS_aqgc_1D_minX<<" ) aQGC constraint ==="<<endl;
			tw.Start();
			calc.Clear(); // Necessary. otherwise the HypoTestInverterResult will always use the old HypoTestInverterResult
			poi1->setMin(func_RatioToSM_VS_aqgc_1D_minX);
			//poi1->setMax(poimax);
			poi1->setMax(scan_range.second);
			calc.SetFixedScan(npoints,poi1->getMin(),poi1->getMax());
			HypoTestInverterResult* htresult_uc = calc.GetInterval();
			assert(htresult_uc!=0);
			cout << "Hypothesis test finished!"<<endl;
			cout << "Time consuming: " << endl;
			tw.Stop();
			tw.Print();
			outputlog << "Hypothesis test finished!"<<endl;
			outputlog << "Time consuming: " << tw.RealTime() << endl;
			double upperLimit_uc = htresult_uc->UpperLimit();
			double ulError_uc = htresult_uc->UpperLimitEstimatedError();
			std::cout << "The computed upper constraint is: " << upperLimit_uc << " +/- " << ulError_uc << std::endl;
			std::cout << " Expected upper constraints, using the alternate model : " << std::endl;
			std::cout << " expected upper constraint (median) " << htresult_uc->GetExpectedUpperLimit(0) << std::endl;
			std::cout << " expected upper constraint (-1 sig) " << htresult_uc->GetExpectedUpperLimit(-1) << std::endl;
			std::cout << " expected upper constraint (+1 sig) " << htresult_uc->GetExpectedUpperLimit(1) << std::endl;
			std::cout << " expected upper constraint (-2 sig) " << htresult_uc->GetExpectedUpperLimit(-2) << std::endl;
			std::cout << " expected upper constraint (+2 sig) " << htresult_uc->GetExpectedUpperLimit(2) << std::endl;
			outputlog << "The computed upper constraint is: " << upperLimit_uc << " +/- " << ulError_uc << std::endl;
			outputlog << " Expected upper constraints, using the alternate model : " << std::endl;
			outputlog << " expected upper constraint (median) " << htresult_uc->GetExpectedUpperLimit(0) << std::endl;
			outputlog << " expected upper constraint (-1 sig) " << htresult_uc->GetExpectedUpperLimit(-1) << std::endl;
			outputlog << " expected upper constraint (+1 sig) " << htresult_uc->GetExpectedUpperLimit(1) << std::endl;
			outputlog << " expected upper constraint (-2 sig) " << htresult_uc->GetExpectedUpperLimit(-2) << std::endl;
			outputlog << " expected upper constraint (+2 sig) " << htresult_uc->GetExpectedUpperLimit(2) << std::endl;
			HypoTestInverterPlot *plot_uc = new HypoTestInverterPlot("HTI_Result_Plot_LowerConstraint","",htresult_uc);
			poi1->setMin(poimin);
			poi1->setMax(poimax);
			func_RatioToSM_VS_aqgc_1D_minX = temp_func_RatioToSM_VS_aqgc_1D_minX;

			// Try to draw the two scan results in one plot
			//TCanvas c_test(Form("Hypothesis tests of lower and upprt constraints: %s for %s",htresult_lc->GetName(),aqgc_para_names[0].c_str()), "", 1600, 600);
			//c_test.SetLogy(false);
			//c_test.Divide(2,1);
			//c_test.cd(1);
			//gPad->SetRightMargin(0);
			//plot_lc->Draw("CLb 2CL");
			//c_test.cd(2);
			//gPad->SetLeftMargin(0);
			//plot_uc->Draw("CLb 2CL");
			//c_test.SaveAs(Form("HypothesisTest_SeparatelyDraw_%s_TS%i.png",aqgc_para_names[0].c_str(),ts_type));

			TCanvas c_all(Form("Hypothesis test: %s for %s",htresult_lc->GetName(),aqgc_para_names[0].c_str()));
			c_all.cd();
			c_all.SetLeftMargin(0.13);
			c_all.SetBottomMargin(0.13);
			TGraph c_all_frame(2);
			c_all_frame.SetTitle(Form(" ; %s [TeV^{-4}]; P value",aqgc_para_names[0].c_str()));
			c_all_frame.GetXaxis()->SetTitleSize(0.05);
			c_all_frame.GetXaxis()->SetLabelSize(0.05);
			c_all_frame.GetYaxis()->SetTitleSize(0.05);
			c_all_frame.GetYaxis()->SetLabelSize(0.05);
			c_all_frame.SetPoint(0,poimin,0);
			c_all_frame.SetPoint(1,poimax,1);
			c_all_frame.GetXaxis()->SetLimits( (scan_range.first-0.3*(scan_range.second-scan_range.first)), (scan_range.second+0.3*(scan_range.second-scan_range.first)) );
			c_all_frame.GetXaxis()->SetRangeUser( (scan_range.first-0.3*(scan_range.second-scan_range.first)), (scan_range.second+0.3*(scan_range.second-scan_range.first)) );
			c_all_frame.GetYaxis()->SetLimits(0,1.9);
			c_all_frame.GetYaxis()->SetRangeUser(0,1.9);
			c_all_frame.SetLineColor(0);
			c_all_frame.SetLineWidth(0);
			c_all_frame.SetMarkerColor(0);
			c_all_frame.SetMarkerSize(0);
			c_all_frame.Draw("AP");
			plot_lc->Draw("EXP SAME");
			plot_uc->Draw("EXP SAME");
			TLine pc(c_all_frame.GetXaxis()->GetXmin(),(1-confidence_level),c_all_frame.GetXaxis()->GetXmax(),(1-confidence_level));
			pc.SetLineColor(kRed);
			pc.SetLineStyle(2);
			pc.SetLineWidth(3);
			pc.Draw();
			c_all.SaveAs(Form("HypothesisTest_%s_TS%i.png",aqgc_para_names[0].c_str(),ts_type));
			outputlog.close();

		}


		//---------------------------------------------
	private:

		std::vector<string> aqgc_para_names; // size < 2
		std::vector<std::pair<double,double>> aqgc_para_cons;

		std::vector<TH1D*>  fvectors_histsamples;
		std::vector<int>    fvectors_histsamplestypes; // 0 for SM signal, 1 for SM background, 2 for AQGC signal
		std::vector<TH1D*>  fvectors_histsamples2;
		std::vector<int>    fvectors_histsamplestypes2; // 0 for SM signal, 1 for SM background, 2 for AQGC signal

		std::vector<double> fvectors_histintegrals;
		std::vector<double> fvectors_aqgcparavalues_x;
		std::vector<double> fvectors_aqgcparavalues_y;

		bool func_is_interpolated;
		double sm_signal_yields;
		double func_RatioToSM_VS_aqgc_1D_minX;
		double sm_background_yields;
		double confidence_level;
		std::pair<double,double> scan_range;
		int altModel_ntoys;
		int nullModel_ntoys;
		int ts_type;
		bool fit_for_ts;
		TF1* func_RatioToSM_VS_aqgc_1D;
		string formula_func_RatioToSM_VS_aqgc_1D;
		TF2* func_RatioToSM_VS_aqgc_2D;
		string formula_func_RatioToSM_VS_aqgc_2D;
		TGraph*  graph_RatioToSM_VS_aqgc_1D;
		TGraph2D* graph_RatioToSM_VS_aqgc_2D;
		TFitResultPtr RatioToSM_VS_aqgc_fitresult;
		RooWorkspace* wspace;
		ModelConfig* aqgc_modelconfig;
		ModelConfig* sm_modelconfig;

};




void calculate_aQGC_constraints(std::vector<string> paranames, std::vector<TH1D*> histsamples, std::vector<int> histsamplestypes, RooDataSet* datasample=0, double CLv=0.95, int TStype=1, int ntoy_alt=1000, int ntoy_null=1000) {
	tools_for_aqgc_constraints hypo_test_model_tool(paranames, histsamples, histsamplestypes);

	hypo_test_model_tool.set_confidence_level(CLv);
	hypo_test_model_tool.set_ts_type(TStype);
	hypo_test_model_tool.set_Models_ntoys(ntoy_alt, ntoy_null);
	// Necessary to set proper range for scanning of different aQGC coefficients, otherwise the program may stack (the underlying reason is unknown)
	if(TString(paranames[0]).Contains("S")) {
		hypo_test_model_tool.set_scan_range(-500, 500);
	}
	else if(TString(paranames[0]).Contains("M")) {
		hypo_test_model_tool.set_scan_range(-60, 60);
	}
	else if(TString(paranames[0]).Contains("M7")) {
		hypo_test_model_tool.set_scan_range(-130, 130);
	}
	else if(TString(paranames[0]).Contains("T")) {
		hypo_test_model_tool.set_scan_range(-25, 25);
	}

	hypo_test_model_tool.calculate_aqgcparavalues_vs_histintegrals();

	hypo_test_model_tool.construct_func_RatioToSM_VS_aqgc_1D();

	hypo_test_model_tool.draw_RatioToSM_VS_aqgc_1D_plot();

	double exp_S_sm = hypo_test_model_tool.obtain_sm_signal_yields();
	cout<<"SM expected signal: "<<exp_S_sm<<endl;
	double exp_B_sm = hypo_test_model_tool.obtain_sm_background_yields();
	cout<<"SM expected background: "<<exp_B_sm<<endl;

	hypo_test_model_tool.setup_aqgc_modelconfigs(exp_S_sm, exp_B_sm, datasample);

	hypo_test_model_tool.buildup_hypotest(datasample);
}










//////////////////////////////////////////////////////////////////////////////////
//==============================================================================//
//              Mutiple Channels Combine to Get Constraonts                     //





class tools_for_aqgc_constraints_subchannel {

	public:
		void clear_objects(){
			assert((paranames.size == 1) || (paranames.size == 2));
			if(!aqgc_para_names.empty())aqgc_para_names.clear();
			else aqgc_para_names = std::vector<string>();
			if(!fvectors_histsamples.empty())fvectors_histsamples.clear();
			else fvectors_histsamples = std::vector<TH1D*>();
			if(!fvectors_histsamplestypes.empty())fvectors_histsamplestypes.clear();
			else fvectors_histsamplestypes = std::vector<int>();
			if(!fvectors_histsamples2.empty())fvectors_histsamples2.clear();
			else fvectors_histsamples2 = std::vector<TH1D*>();
			if(!fvectors_histsamplestypes2.empty())fvectors_histsamplestypes2.clear();
			else fvectors_histsamplestypes2 = std::vector<int>();
			set_channel_label();
			set_poi1();
			set_poi2();
			set_obs();
			sm_signal_yields = 0;
			func_RatioToSM_VS_aqgc_1D_minX = 0;
			sm_background_yields = 0;
			func_RatioToSM_VS_aqgc_1D = 0;
			formula_func_RatioToSM_VS_aqgc_1D = "";
			func_RatioToSM_VS_aqgc_2D = 0;
			formula_func_RatioToSM_VS_aqgc_2D = "";
			graph_RatioToSM_VS_aqgc_1D = 0;
			graph_RatioToSM_VS_aqgc_2D = 0;
			RatioToSM_VS_aqgc_fitresult = 0;
		}

		void set_channel_label(string channellabel=""){
			channel_label = channellabel;
		}

		void set_poi1(RooRealVar* common_p1=0){
			poi1 = common_p1;
		}

		void set_poi2(RooRealVar* common_p2=0){
			poi2 = common_p2;
		}

		void set_obs(RooRealVar* common_obs=0){
			obs = common_obs;
		}


		double obtain_significance(double par_value) {
			double aqgc_sig = (func_RatioToSM_VS_aqgc_1D->Eval(par_value));
			double sm_bkg = sm_background_yields;
			double significance = sqrt(2.0*((aqgc_sig+sm_bkg)*(TMath::Log(1+(aqgc_sig/sm_bkg)))-aqgc_sig));
			return significance;
		}

		void print_par_value_of_specific_significance(double target_significance) {
			ofstream output(Form("./CombinedHypoTestInvLog_%s.log", aqgc_para_names[0].c_str()), ios::app);
			cout << " +++ --------- Print AQGC "<< aqgc_para_names[0] << " values for specific significance at [significance = " << target_significance<<"] --------- +++ " << endl;
			output << " +++ --------- Print AQGC "<< aqgc_para_names[0] << " values for specific significance at [significance = " << target_significance<<"] --------- +++ " << endl;
			double temp_significance = 0;
			double nsteps = 100;
			double step = (20 - (-20))/nsteps;
			for(int istep=0; istep<nsteps; istep++) {
				if( ((obtain_significance(-20+(istep*step)) - target_significance) * (obtain_significance(-20+((istep+1)*step)) - target_significance)) <= 0 ) {
					cout << " --- AQGC "<< aqgc_para_names[0] << " value is " << (-20+((istep+1)*step)) << " --- " << endl;
					output << " --- AQGC "<< aqgc_para_names[0] << " value is " << (-20+((istep+1)*step)) << " --- " << endl;
				}
			}
			cout << " --- --------- Print ended --------- --- " <<endl;
			output << " --- --------- Print ended --------- --- " <<endl;
			output.close();
		}
				

		tools_for_aqgc_constraints_subchannel(std::vector<string> paranames, string channellabel, RooRealVar* common_p1, RooRealVar* common_obs){
			assert(paranames.size == 1);
			clear_objects();
			for(int i=0; i<paranames.size(); i++){
				aqgc_para_names.push_back(paranames[i]);
			}
			set_channel_label(channellabel);
			set_poi1(common_p1);
			set_obs(common_obs);
			aqgc_para_cons = std::vector<std::pair<double,double>>();
			func_is_interpolated = false;
			if(func_RatioToSM_VS_aqgc_2D!=0 && func_RatioToSM_VS_aqgc_2D->IsValid()) func_is_interpolated = true;
			fvectors_histintegrals = std::vector<double>();
			fvectors_aqgcparavalues_x = std::vector<double>();
			fvectors_aqgcparavalues_y = std::vector<double>();
			wspace = new RooWorkspace(Form("wspace_%s",channel_label.c_str()));
			cout<<" Tools for aQGC of"<<channel_label<<" are initiated !"<<endl;
		}

		tools_for_aqgc_constraints_subchannel(std::vector<string> paranames, string channellabel, RooRealVar* common_p1, RooRealVar* common_obs, std::vector<TH1D*> histsamples,  std::vector<int>histsamplestypes){
			assert(paranames.size == 1);
			clear_objects();
			for(int i=0; i<paranames.size(); i++){
				aqgc_para_names.push_back(paranames[i]);
			}
			set_channel_label(channellabel);
			set_poi1(common_p1);
			set_obs(common_obs);
			aqgc_para_cons = std::vector<std::pair<double,double>>();
			for(int i=0; i<histsamples.size(); i++){
				fvectors_histsamples.push_back(histsamples[i]);
			}
			for(int i=0; i<histsamplestypes.size(); i++){
				fvectors_histsamplestypes.push_back(histsamplestypes[i]);
			}
			func_is_interpolated = false;
			if(func_RatioToSM_VS_aqgc_2D!=0 && func_RatioToSM_VS_aqgc_2D->IsValid()) func_is_interpolated = true;
			fvectors_histintegrals = std::vector<double>();
			fvectors_aqgcparavalues_x = std::vector<double>();
			fvectors_aqgcparavalues_y = std::vector<double>();
			wspace = new RooWorkspace(Form("wspace_%s",channel_label.c_str()));
			cout<<" Tools for aQGC of"<<channel_label<<" are initiated !"<<endl;
		}

		tools_for_aqgc_constraints_subchannel(std::vector<string> paranames, string channellabel, RooRealVar* common_p1, RooRealVar* common_p2, RooRealVar* common_obs, std::vector<TH1D*> histsamples,  std::vector<int>histsamplestypes, std::vector<TH1D*> histsamples2,  std::vector<int>histsamplestypes2){
			assert(paranames.size == 2);
			clear_objects();
			for(int i=0; i<paranames.size(); i++){
				aqgc_para_names.push_back(paranames[i]);
			}
			set_channel_label(channellabel);
			set_poi1(common_p1);
			set_poi2(common_p2);
			set_obs(common_obs);
			aqgc_para_cons = std::vector<std::pair<double,double>>();
			for(int i=0; i<histsamples.size(); i++){
				fvectors_histsamples.push_back(histsamples[i]);
			}
			for(int i=0; i<histsamplestypes.size(); i++){
				fvectors_histsamplestypes.push_back(histsamplestypes[i]);
			}
			for(int i=0; i<histsamples2.size(); i++){
				fvectors_histsamples2.push_back(histsamples2[i]);
			}
			for(int i=0; i<histsamplestypes2.size(); i++){
				fvectors_histsamplestypes2.push_back(histsamplestypes2[i]);
			}
			func_is_interpolated = false;
			if(func_RatioToSM_VS_aqgc_2D!=0 && func_RatioToSM_VS_aqgc_2D->IsValid()) func_is_interpolated = true;
			fvectors_histintegrals = std::vector<double>();
			fvectors_aqgcparavalues_x = std::vector<double>();
			fvectors_aqgcparavalues_y = std::vector<double>();
			wspace = new RooWorkspace(Form("wspace_%s",channel_label.c_str()));
			cout<<" Tools for aQGC of"<<channel_label<<" are initiated !"<<endl;
		}

		~tools_for_aqgc_constraints_subchannel(){};

		void add_samples(std::vector<TH1D*> histsamples, std::vector<int>histsamplestypes){
			for(int i=0; i<histsamples.size(); i++){
				fvectors_histsamples.push_back(histsamples[i]);
			}
			for(int i=0; i<histsamplestypes.size(); i++){
				fvectors_histsamplestypes.push_back(histsamplestypes[i]);
			}
		}

		void calculate_aqgcparavalues_vs_histintegrals(){
			assert((aqgc_para_names.size() == 1) || (aqgc_para_names.size() == 2));
			if(aqgc_para_names.size() == 2) { cout<<"!!! Should construct 2D function, code is not prepared yet !!!"<<endl; return; }
			string value;
			double dvalue;
			std::vector<double> temp_fvectors_histintegrals;
			std::vector<double> temp_fvectors_aqgcparavalues;
			for (int i=0; i<fvectors_histsamples.size(); i++) {
				if((TString(fvectors_histsamples[i]->GetName())).Contains("th_aqgcRW") && (fvectors_histsamplestypes[i]==2)) {
					temp_fvectors_histintegrals.push_back((fvectors_histsamples[i]->Integral(1,fvectors_histsamples[i]->GetNcells()-2)));
					value = fvectors_histsamples[i]->GetName();
					value.erase(0, value.find(Form("="),0)+1);
					//value.erase(value.find(Form("e-12"),0), value.size()-1-value.find(Form("e-12")));
					value.erase(value.find(Form("e-12"),0));
					dvalue = atof(value.c_str());
					temp_fvectors_aqgcparavalues.push_back(dvalue);
				}

			}
			std::vector<int> index(temp_fvectors_aqgcparavalues.size());
			TMath::SortItr(temp_fvectors_aqgcparavalues.begin(), temp_fvectors_aqgcparavalues.end(), index.begin(), false );
			for(int i=0; i<index.size(); i++){
				fvectors_histintegrals.push_back(temp_fvectors_histintegrals[index[i]]);
				fvectors_aqgcparavalues_x.push_back(temp_fvectors_aqgcparavalues[index[i]]);
			}
			cout<<" vectors of aqgcparavalues_vs_histintegrals is prepared !"<<endl;
		}


		void construct_func_RatioToSM_VS_aqgc_1D(){
			assert(fvectors_aqgcparavalues_x.size() > 2);
			graph_RatioToSM_VS_aqgc_1D = new TGraph(fvectors_aqgcparavalues_x.size());
			//graph_RatioToSM_VS_aqgc_1D->SetTitle("graph_RatioToSM_VS_aqgc_1D");
			for(int i=0; i<fvectors_aqgcparavalues_x.size(); i++){
				graph_RatioToSM_VS_aqgc_1D->SetPoint(i, fvectors_aqgcparavalues_x[i], fvectors_histintegrals[i]);
			}
			if(!func_is_interpolated){
				func_RatioToSM_VS_aqgc_1D = new TF1("func_RatioToSM_VS_aqgc_1D","[0]+[1]*x+[2]*x*x",-150,150);
				RatioToSM_VS_aqgc_fitresult = graph_RatioToSM_VS_aqgc_1D->Fit(func_RatioToSM_VS_aqgc_1D, "SRE");
				assert(RatioToSM_VS_aqgc_fitresult->IsValid());
			}
			if(func_is_interpolated){
				cout<<"No available method founded!"<<endl;
				//func_RatioToSM_VS_aqgc_1D = new TF1("func_RatioToSM_VS_aqgc_1D",[&](double*x){ return graph_RatioToSM_VS_aqgc_1D->Eval(x[0]); }, graph_RatioToSM_VS_aqgc_1D->GetXaxis()->GetXmin(), graph_RatioToSM_VS_aqgc_1D->GetXaxis()->GetXmax(), 0);
				assert(func_RatioToSM_VS_aqgc_1D);
				assert(RatioToSM_VS_aqgc_fitresult->IsValid());
			}
			//sm_signal_yields = func_RatioToSM_VS_aqgc_1D->Eval(0); 
			//sm_signal_yields = func_RatioToSM_VS_aqgc_1D->GetMinimum(graph_RatioToSM_VS_aqgc_1D->GetXaxis()->GetXmin(), graph_RatioToSM_VS_aqgc_1D->GetXaxis()->GetXmax());

			func_RatioToSM_VS_aqgc_1D_minX = (-0.5*(func_RatioToSM_VS_aqgc_1D->GetParameter(1))/(func_RatioToSM_VS_aqgc_1D->GetParameter(2))); sm_signal_yields = func_RatioToSM_VS_aqgc_1D->Eval(func_RatioToSM_VS_aqgc_1D_minX);
			//for(int is=0; is<fvectors_histsamplestypes.size(); is++){
			//	if(fvectors_histsamplestypes.at(is)==0) sm_signal_yields += (fvectors_histsamples[is])->Integral(1, (fvectors_histsamples[is])->GetNcells()-2);
			//}
			for(int is=0; is<fvectors_histsamplestypes.size(); is++){
				if(fvectors_histsamplestypes.at(is)==1) sm_background_yields += (fvectors_histsamples[is])->Integral(1, (fvectors_histsamples[is])->GetNcells()-2);
			}
			formula_func_RatioToSM_VS_aqgc_1D = Form("( %e + %e * poi1 + %e * poi1 * poi1 )", func_RatioToSM_VS_aqgc_1D->GetParameter(0), func_RatioToSM_VS_aqgc_1D->GetParameter(1), func_RatioToSM_VS_aqgc_1D->GetParameter(2)); // formula_func_RatioToSM_VS_aqgc_1D = func_RatioToSM_VS_aqgc_1D->GetFormula();
			cout<<" func_RatioToSM_VS_aqgc_1D is builded !"<<endl;
		}


		void construct_func_RatioToSM_VS_aqgc_2D(){
			assert(fvectors_aqgcparavalues_x.size() > 2);
			// TO BE ACCOMPLISHED ...
		}


		void draw_RatioToSM_VS_aqgc_1D_plot(){
			TCanvas tempc("temp_draw_RatioToSM_VS_aqgc_1D_plot","");
			tempc.SetRightMargin(0.25);
			tempc.cd();
			tempc.SetLeftMargin(0.16);
			tempc.SetBottomMargin(0.13);
			graph_RatioToSM_VS_aqgc_1D->SetMarkerStyle(20);
			graph_RatioToSM_VS_aqgc_1D->SetMarkerSize(1.5);
			graph_RatioToSM_VS_aqgc_1D->SetMarkerColor(1);
			graph_RatioToSM_VS_aqgc_1D->SetLineColor(1);
			//graph_RatioToSM_VS_aqgc_1D->SetTitle(Form("; %s; #frac{#sigma^{aQGC}_{sig}}{#sigma^{SM}_{sig}}",aqgc_para_names[0].c_str()));
			graph_RatioToSM_VS_aqgc_1D->SetTitle(Form("; %s; N^{aQGC}_{sig}",aqgc_para_names[0].c_str()));
			graph_RatioToSM_VS_aqgc_1D->GetXaxis()->SetTitleSize(0.05);
			graph_RatioToSM_VS_aqgc_1D->GetXaxis()->SetLabelSize(0.05);
			graph_RatioToSM_VS_aqgc_1D->GetYaxis()->SetTitleSize(0.05);
			graph_RatioToSM_VS_aqgc_1D->GetYaxis()->SetLabelSize(0.05);
			graph_RatioToSM_VS_aqgc_1D->Draw("AP");
			func_RatioToSM_VS_aqgc_1D->SetLineColor(4);
			func_RatioToSM_VS_aqgc_1D->Draw("SAME");
			TLegend lg(0.75,0.80,0.99,0.99);
			lg.AddEntry(graph_RatioToSM_VS_aqgc_1D, "Data points",   "p");
			lg.AddEntry(func_RatioToSM_VS_aqgc_1D,  "Quadratic fit", "l");
			lg.AddEntry(func_RatioToSM_VS_aqgc_1D, Form("C0 = %.3e",func_RatioToSM_VS_aqgc_1D->GetParameter(0)), "p");
			lg.AddEntry(func_RatioToSM_VS_aqgc_1D, Form("C1 = %.3e",func_RatioToSM_VS_aqgc_1D->GetParameter(1)), "p");
			lg.AddEntry(func_RatioToSM_VS_aqgc_1D, Form("C2 = %.3e",func_RatioToSM_VS_aqgc_1D->GetParameter(2)), "p");
			lg.Draw();
			tempc.SaveAs(Form("RatioToSM_VS_aqgc_1D_plot_%s_%s.png",aqgc_para_names[0].c_str(),channel_label.c_str()));
			cout<<" RatioToSM_VS_aqgc_1D plot is drawn !"<<endl;
		}

		double eval_sig_expect_1D(double x){
			double sig_expect = 0;
			if(!func_is_interpolated){
				assert(RatioToSM_VS_aqgc_fitresult->IsValid());
				sig_expect = func_RatioToSM_VS_aqgc_1D->Eval(x);
			}
			else{
				sig_expect = graph_RatioToSM_VS_aqgc_1D->Eval(x);
			}
			return sig_expect;
		}

		double obtain_sm_signal_yields(){
			return sm_signal_yields;
		}

		double obtain_func_RatioToSM_VS_aqgc_1D_minX(){
			return func_RatioToSM_VS_aqgc_1D_minX;
		}

		double obtain_sm_background_yields(){
			return sm_background_yields;
		}

		TF1* obtain_func_1D() {
			assert(func_RatioToSM_VS_aqgc_1D);
			return func_RatioToSM_VS_aqgc_1D;
		}

		TF1* obtain_func_2D() {
			assert(func_RatioToSM_VS_aqgc_2D);
			return func_RatioToSM_VS_aqgc_2D;
		}

		RooWorkspace* obtain_wspace(){
			return wspace;
		}

		double obtain_func_RatioToSM_VS_aqgc_1D_a0(){
			return func_RatioToSM_VS_aqgc_1D->GetParameter(0);
		}

		double obtain_func_RatioToSM_VS_aqgc_1D_a1(){
			return func_RatioToSM_VS_aqgc_1D->GetParameter(1);
		}

		double obtain_func_RatioToSM_VS_aqgc_1D_a2(){
			return func_RatioToSM_VS_aqgc_1D->GetParameter(2);
		}

		//------------------------
		// Manners used for hypothesis test
		//------------------------
		void setup_aqgc_worksapce(double exp_s=5.0, double exp_b=15.0) {
			// ---------------------------------------------------------------------------------------------------------------------------- 
			wspace->import(*poi1);
			wspace->import(*obs);
			wspace->factory( Form("Bkg_%s[%e]", channel_label.c_str(), exp_b) );
			RooRealVar *b = (RooRealVar*) wspace->var(Form("Bkg_%s", channel_label.c_str()));
			b->setConstant();
			wspace->factory( Form("a0_%s[%e]" , channel_label.c_str(), func_RatioToSM_VS_aqgc_1D->GetParameter(0)) );
			wspace->factory( Form("a1_%s[%e]" , channel_label.c_str(), func_RatioToSM_VS_aqgc_1D->GetParameter(1)) );
			wspace->factory( Form("a2_%s[%e]" , channel_label.c_str(), func_RatioToSM_VS_aqgc_1D->GetParameter(2)) );
			wspace->factory( Form("RooPolyVar::Sig_%s(%s,{a0_%s,a1_%s,a2_%s})", channel_label.c_str(), Form("common_poi_%s",aqgc_para_names[0].c_str()), channel_label.c_str(), channel_label.c_str(), channel_label.c_str()) ); // signal counting expected in poisson distribution
			wspace->factory( Form("sum::Total_%s(Sig_%s,Bkg_%s)", channel_label.c_str(), channel_label.c_str(), channel_label.c_str()) );
			wspace->factory( Form("RooPoisson::countingModel_%s(%s,Total_%s)", channel_label.c_str(), obs->GetName(), channel_label.c_str()) ); // counting model
			//wspace->factory( Form("sum::Total(1.0,2.0)");
			//wspace->factory( Form("RooPoisson::countingModel(%s,Total)", obs->GetName()) ); // counting model
			// --------------------------- terms for background uncertainty ---------------------------
			//wspace->factory( Form("expectB[30.0]") ); // nuisance: expected background
			//wspace->factory( Form("ratioBkg[1.0,0.0,3.0]") ); // nuisance: ratio on expected background to manifest its uncertainty
			//wspace->factory( Form("globalBkgObsNum[10.0,0.0,300.0]") ); // global obs for nuisisance: ratioBkg
			//wspace->factory( Form("prod::globalBkg(globalBkgExpNum[10.0],ratioBkg)") ); // expected global obs for nuisisance: ratioBkg
			//wspace->factory( Form("Gaussian::bkgConstraint(globalBkgObsNum,globalBkg,10.0)") );  // nuisance gaussian constraint: background --- the sigma 10.0 is set somehow arbitially
			//wspace->factory( Form("prod::Bkg(ratioBkg,expectB") ); 
			//wspace->factory( Form("expectS[%e]",exp_s) );
			//wspace->factory( Form("RooPolyVar::ratioSig(%s,{a0,a1,a2})",Form("common_poi_%s",aqgc_para_names[0].c_str())) );
			//wspace->factory( Form("prod::Sig(ratioSig,expectS]") );
			/// =========================================================================================================================== 
			///			wspace->factory( Form("obs[10.0,0.0,300.0]") ); // obs: counting number
			///			wspace->factory( Form("poi1[%e,%e,%e]",((scan_range.first+scan_range.second)/2.0),(scan_range.first),(scan_range.second)) ); // poi: the first POI
			///			wspace->import(*obs);
			///			wspace->import(*poi1);
			///			if(poi2!=0 && poi2->isValid()) wspace->import(*poi2);
			///			RooRealVar* Bkg = new RooRealVar(Form("Bkg_%s",channel_label.c_str()), exp_b);
			///			RooRealVar* a0  = new RooRealVar(Form("a0_%s" ,channel_label.c_str()), func_RatioToSM_VS_aqgc_1D->GetParameter(0));
			///			RooRealVar* a1  = new RooRealVar(Form("a1_%s" ,channel_label.c_str()), func_RatioToSM_VS_aqgc_1D->GetParameter(1));
			///			RooRealVar* a2  = new RooRealVar(Form("a2_%s" ,channel_label.c_str()), func_RatioToSM_VS_aqgc_1D->GetParameter(2));
			///			wspace->import(*Bkg);
			///			wspace->import(*a0 );
			///			wspace->import(*a1 );
			///			wspace->import(*a2 );
			///			RooPolyVar* Sig = new RooPolyVar(Form("Sig_%s",channel_label.c_str()), Form("Sig_%s",channel_label.c_str()), *poi1, RooArgList(*a0,*a1,*a2));
			///			RooAddition* sum = new RooAddition(Form("sum_%s",channel_label.c_str()), Form("sum_%s",channel_label.c_str()), RooArgList(*Sig,*Bkg));
			///			RooPoisson* countingModel = new RooPoisson(Form("countingModel_%s",channel_label.c_str()), Form("countingModel_%s",channel_label.c_str()), *obs, *sum);
			///			wspace->import(*Sig);
			///			wspace->import(*sum);
			///			wspace->import(*countingModel);


			//wspace->SetName(Form("wspace_%s",channel_label.c_str()));
			//wspace->writeToFile(Form("wspace_%s.root",aqgc_para_names[0].c_str()));
			wspace->Print();
			cout<<" Model config is setup !"<<endl;
		}



		//---------------------------------------------
	private:

		std::vector<string> aqgc_para_names; // size < 2
		string channel_label;
		RooRealVar* poi1;
		RooRealVar* poi2;
		RooRealVar* obs;
		std::vector<std::pair<double,double>> aqgc_para_cons;

		std::vector<TH1D*>  fvectors_histsamples;
		std::vector<int>    fvectors_histsamplestypes; // 0 for SM signal, 1 for SM background, 2 for AQGC signal
		std::vector<TH1D*>  fvectors_histsamples2;
		std::vector<int>    fvectors_histsamplestypes2; // 0 for SM signal, 1 for SM background, 2 for AQGC signal

		std::vector<double> fvectors_histintegrals;
		std::vector<double> fvectors_aqgcparavalues_x;
		std::vector<double> fvectors_aqgcparavalues_y;

		bool func_is_interpolated;
		double sm_signal_yields;
		double func_RatioToSM_VS_aqgc_1D_minX;
		double sm_background_yields;
		TF1* func_RatioToSM_VS_aqgc_1D;
		string formula_func_RatioToSM_VS_aqgc_1D;
		TF2* func_RatioToSM_VS_aqgc_2D;
		string formula_func_RatioToSM_VS_aqgc_2D;
		TGraph*  graph_RatioToSM_VS_aqgc_1D;
		TGraph2D* graph_RatioToSM_VS_aqgc_2D;
		TFitResultPtr RatioToSM_VS_aqgc_fitresult;
		RooWorkspace* wspace;

};



void CombinedHypoTest(std::vector<string> aqgc_para_names, RooArgSet *obs, RooRealVar *poi, std::pair<double,double> scan_range, ModelConfig* aqgc_modelconfig, ModelConfig* sm_modelconfig, RooDataSet* data, double confidence_level=0.95, int ts_type=1, bool fit_for_ts=false, int altModel_ntoys=1000, int nullModel_ntoys=1000, double func_RatioToSM_VS_aqgc_1D_minX=0) {

	ofstream outputlog(Form("./CombinedHypoTestInvLog_%s.log", aqgc_para_names[0].c_str()), ios::ate|ios::out);
	// First, set up the SM prediction model (Null hypothesis): clone from the generic aQGC model (alternative hypothesis) with aQGC parameter set to be 0
	assert(aqgc_modelconfig);
		
		
	assert(sm_modelconfig);
	// Second, set up the hypothsis test calculator: here we use the frequentist one. And will pass it to hypothesis inversion calculator for upper limit calculation
	FrequentistCalculator fc(*data, *sm_modelconfig, *aqgc_modelconfig);
	fc.SetToys(altModel_ntoys,nullModel_ntoys); 
	fc.UseSameAltToys();
	fc.StoreFitInfo(true);
	// Third, choose the test statistic: here we use the ProfileLikelihoodTestStat that calculates the profile likelihood ratio
	/// --- Test statistic 1: SimpleLikelihoodRatioTestStat
	SimpleLikelihoodRatioTestStat teststat1(*(aqgc_modelconfig->GetPdf()), *(sm_modelconfig->GetPdf()));
	////// ------ null parameters must includes snapshot of poi plus the nuisance values
	if(fit_for_ts && ts_type==1) {
		RooArgSet constrainParams;
		if(aqgc_modelconfig->GetNuisanceParameters()) constrainParams.add(*aqgc_modelconfig->GetNuisanceParameters());
		RooStats::RemoveConstantParameters(&constrainParams);
		RooFitResult* fitres = 0;
		int itenum = 0;
		while( (fitres==0 || fitres->status()!=0) && (itenum<=5) ) {
			itenum+=1;
			fitres = aqgc_modelconfig->GetPdf()->fitTo(*data, RooFit::InitialHesse(false), RooFit::Hesse(false), RooFit::Minimizer("Minuit2", "Migrad"), RooFit::Strategy(1), RooFit::PrintLevel(-1), RooFit::Constrain(constrainParams), RooFit::Save(true), RooFit::Offset(true));
		}
		aqgc_modelconfig->SetSnapshot(*aqgc_modelconfig->GetParametersOfInterest());
		cout<<"The snapshot of the model for aQGC (null hypothesis) is set to the fit result of data!"<<endl;
		outputlog<<"The snapshot of the model for aQGC (null hypothesis) is set to the fit result of data!"<<endl;
	}
	RooArgSet nullParams(*aqgc_modelconfig->GetSnapshot());
	if(aqgc_modelconfig->GetNuisanceParameters()) nullParams.add(*aqgc_modelconfig->GetNuisanceParameters());
	if(aqgc_modelconfig->GetSnapshot()) teststat1.SetNullParameters(nullParams);
	RooArgSet altParams(*sm_modelconfig->GetSnapshot());
	if(sm_modelconfig->GetNuisanceParameters()) altParams.add(*sm_modelconfig->GetNuisanceParameters());
	if(sm_modelconfig->GetSnapshot()) teststat1.SetAltParameters(altParams);
	/// --- Test statistic 2: RatioOfProfiledLikelihoodsTestStat
	RatioOfProfiledLikelihoodsTestStat teststat2(*aqgc_modelconfig->GetPdf(), *sm_modelconfig->GetPdf(), sm_modelconfig->GetSnapshot());
	teststat2.SetSubtractMLE(false);
	/// --- Test statistic 3: ProfileLikelihoodTestStat
	ProfileLikelihoodTestStat teststat3(*(aqgc_modelconfig->GetPdf()));
	teststat3.SetOneSided(false);
	teststat3.SetSigned(false);
	// Fourth, configure the toyMC sampler
	TestStatistic* teststat = 0;
	if(ts_type == 1) teststat = &teststat1;
	if(ts_type == 2) teststat = &teststat2;
	if(ts_type == 3) teststat = &teststat3;
	ToyMCSampler* toymcs = (ToyMCSampler*)fc.GetTestStatSampler();
	toymcs->SetNEventsPerToy(1);
	toymcs->SetTestStatistic(teststat);
	int npoints = 50;
	double poimin = poi->getMin();
	double poimax = poi->getMax();
	cout << "Doing a fixed scan in interval: [" << poimin << " , " << poimax << "]" << endl;
	outputlog << "Doing a fixed scan in interval: [" << poimin << " , " << poimax << "]" << endl;
	// Fifth, draw plots and print results
	TStopwatch tw;
	/// For the lower ( < 0 ) aQGC constraint
	double temp_func_RatioToSM_VS_aqgc_1D_minX = func_RatioToSM_VS_aqgc_1D_minX;
	func_RatioToSM_VS_aqgc_1D_minX = 0; // Not sure if this should be force to 0
	cout <<"=== Estimate the lower ( < "<<func_RatioToSM_VS_aqgc_1D_minX<<" ) aQGC constraint ==="<<endl;
	outputlog <<"=== Estimate the lower ( < "<<func_RatioToSM_VS_aqgc_1D_minX<<" ) aQGC constraint ==="<<endl;
	tw.Start();
	HypoTestInverter calc(fc);
	calc.SetConfidenceLevel(confidence_level);
	//calc.UseCLs(false);
	calc.UseCLs(true);
	calc.SetVerbose(true);
	//poi->setMin(poimin);
	poi->setMin(scan_range.first);
	poi->setMax(func_RatioToSM_VS_aqgc_1D_minX);
	//calc.SetFixedScan(npoints,poimin,0.0);
	calc.SetFixedScan(npoints,poi->getMin(),poi->getMax());
	HypoTestInverterResult* htresult_lc = calc.GetInterval();
	assert(htresult_lc!=0);
	cout << "CombinedHypothesis test finished!"<<endl;
	cout << "Time consuming: " << endl;
	outputlog << "CombinedHypothesis test finished!"<<endl;
	outputlog << "Time consuming: " << endl;
	tw.Stop();
	tw.Print();
	double upperLimit_lc = htresult_lc->UpperLimit();
	double ulError_lc = htresult_lc->UpperLimitEstimatedError();
	std::cout << "The computed lower constraint is: " << upperLimit_lc << " +/- " << ulError_lc << std::endl;
	std::cout << " Expected lower constraints, using the alternate model : " << std::endl;
	std::cout << " expected lower constraint (median) " << htresult_lc->GetExpectedUpperLimit(0) << std::endl;
	std::cout << " expected lower constraint (-1 sig) " << htresult_lc->GetExpectedUpperLimit(-1) << std::endl;
	std::cout << " expected lower constraint (+1 sig) " << htresult_lc->GetExpectedUpperLimit(1) << std::endl;
	std::cout << " expected lower constraint (-2 sig) " << htresult_lc->GetExpectedUpperLimit(-2) << std::endl;
	std::cout << " expected lower constraint (+2 sig) " << htresult_lc->GetExpectedUpperLimit(2) << std::endl;
	outputlog << "The computed lower constraint is: " << upperLimit_lc << " +/- " << ulError_lc << std::endl;
	outputlog << " Expected lower constraints, using the alternate model : " << std::endl;
	outputlog << " expected lower constraint (median) " << htresult_lc->GetExpectedUpperLimit(0) << std::endl;
	outputlog << " expected lower constraint (-1 sig) " << htresult_lc->GetExpectedUpperLimit(-1) << std::endl;
	outputlog << " expected lower constraint (+1 sig) " << htresult_lc->GetExpectedUpperLimit(1) << std::endl;
	outputlog << " expected lower constraint (-2 sig) " << htresult_lc->GetExpectedUpperLimit(-2) << std::endl;
	outputlog << " expected lower constraint (+2 sig) " << htresult_lc->GetExpectedUpperLimit(2) << std::endl;
	HypoTestInverterPlot *plot_lc = new HypoTestInverterPlot("HTI_Result_Plot_LowerConstraint","",htresult_lc);
	poi->setMin(poimin);
	poi->setMax(poimax);
	/// For the upper ( > func_RatioToSM_VS_aqgc_1D_minX ) aQGC constraint
	cout<<"=== Estimate the upper ( > "<<func_RatioToSM_VS_aqgc_1D_minX<<" ) aQGC constraint ==="<<endl;
	outputlog<<"=== Estimate the upper ( > "<<func_RatioToSM_VS_aqgc_1D_minX<<" ) aQGC constraint ==="<<endl;
	tw.Start();
	calc.Clear(); // Necessary. otherwise the HypoTestInverterResult will always use the old HypoTestInverterResult
	poi->setMin(func_RatioToSM_VS_aqgc_1D_minX);
	//poi->setMax(poimax);
	poi->setMax(scan_range.second);
	calc.SetFixedScan(npoints,poi->getMin(),poi->getMax());
	HypoTestInverterResult* htresult_uc = calc.GetInterval();
	assert(htresult_uc!=0);
	cout << "CombinedHypothesis test finished!"<<endl;
	cout << "Time consuming: " << endl;
	tw.Stop();
	tw.Print();
	outputlog << "CombinedHypothesis test finished!"<<endl;
	outputlog << "Time consuming: " << tw.RealTime() << endl;
	double upperLimit_uc = htresult_uc->UpperLimit();
	double ulError_uc = htresult_uc->UpperLimitEstimatedError();
	std::cout << "The computed upper constraint is: " << upperLimit_uc << " +/- " << ulError_uc << std::endl;
	std::cout << " Expected upper constraints, using the alternate model : " << std::endl;
	std::cout << " expected upper constraint (median) " << htresult_uc->GetExpectedUpperLimit(0) << std::endl;
	std::cout << " expected upper constraint (-1 sig) " << htresult_uc->GetExpectedUpperLimit(-1) << std::endl;
	std::cout << " expected upper constraint (+1 sig) " << htresult_uc->GetExpectedUpperLimit(1) << std::endl;
	std::cout << " expected upper constraint (-2 sig) " << htresult_uc->GetExpectedUpperLimit(-2) << std::endl;
	std::cout << " expected upper constraint (+2 sig) " << htresult_uc->GetExpectedUpperLimit(2) << std::endl;
	outputlog << "The computed upper constraint is: " << upperLimit_uc << " +/- " << ulError_uc << std::endl;
	outputlog << " Expected upper constraints, using the alternate model : " << std::endl;
	outputlog << " expected upper constraint (median) " << htresult_uc->GetExpectedUpperLimit(0) << std::endl;
	outputlog << " expected upper constraint (-1 sig) " << htresult_uc->GetExpectedUpperLimit(-1) << std::endl;
	outputlog << " expected upper constraint (+1 sig) " << htresult_uc->GetExpectedUpperLimit(1) << std::endl;
	outputlog << " expected upper constraint (-2 sig) " << htresult_uc->GetExpectedUpperLimit(-2) << std::endl;
	outputlog << " expected upper constraint (+2 sig) " << htresult_uc->GetExpectedUpperLimit(2) << std::endl;
	HypoTestInverterPlot *plot_uc = new HypoTestInverterPlot("HTI_Result_Plot_LowerConstraint","",htresult_uc);
	poi->setMin(poimin);
	poi->setMax(poimax);
	func_RatioToSM_VS_aqgc_1D_minX = temp_func_RatioToSM_VS_aqgc_1D_minX;

	// Try to draw the two scan results in one plot
	//TCanvas c_test(Form("CombinedHypothesis tests of lower and upprt constraints: %s for %s",htresult_lc->GetName(),aqgc_para_names[0].c_str()), "", 1600, 600);
	//c_test.SetLogy(false);
	//c_test.Divide(2,1);
	//c_test.cd(1);
	//gPad->SetRightMargin(0);
	//plot_lc->Draw("CLb 2CL");
	//c_test.cd(2);
	//gPad->SetLeftMargin(0);
	//plot_uc->Draw("CLb 2CL");
	//c_test.SaveAs(Form("CombinedHypothesisTest_SeparatelyDraw_%s_TS%i.png",aqgc_para_names[0].c_str(),ts_type));

	TCanvas c_all(Form("CombinedHypothesis test: %s for %s",htresult_lc->GetName(),aqgc_para_names[0].c_str()));
	c_all.cd();
	c_all.SetLeftMargin(0.13);
	c_all.SetBottomMargin(0.13);
	TGraph c_all_frame(2);
	c_all_frame.SetTitle(Form(" ; %s [TeV^{-4}]; P value",aqgc_para_names[0].c_str()));
	c_all_frame.GetXaxis()->SetTitleSize(0.05);
	c_all_frame.GetXaxis()->SetLabelSize(0.05);
	c_all_frame.GetYaxis()->SetTitleSize(0.05);
	c_all_frame.GetYaxis()->SetLabelSize(0.05);
	c_all_frame.SetPoint(0,poimin,0);
	c_all_frame.SetPoint(1,poimax,1);
	c_all_frame.GetXaxis()->SetLimits( (scan_range.first-0.3*(scan_range.second-scan_range.first)), (scan_range.second+0.3*(scan_range.second-scan_range.first)) );
	c_all_frame.GetXaxis()->SetRangeUser( (scan_range.first-0.3*(scan_range.second-scan_range.first)), (scan_range.second+0.3*(scan_range.second-scan_range.first)) );
	c_all_frame.GetYaxis()->SetLimits(0,1.9);
	c_all_frame.GetYaxis()->SetRangeUser(0,1.9);
	c_all_frame.SetLineColor(0);
	c_all_frame.SetLineWidth(0);
	c_all_frame.SetMarkerColor(0);
	c_all_frame.SetMarkerSize(0);
	c_all_frame.Draw("AP");
	plot_lc->Draw("EXP SAME");
	plot_uc->Draw("EXP SAME");
	TLine pc(c_all_frame.GetXaxis()->GetXmin(),(1-confidence_level),c_all_frame.GetXaxis()->GetXmax(),(1-confidence_level));
	pc.SetLineColor(kRed);
	pc.SetLineStyle(2);
	pc.SetLineWidth(3);
	pc.Draw();
	c_all.SaveAs(Form("CombinedHypothesisTest_%s_TS%i.png",aqgc_para_names[0].c_str(),ts_type));
	outputlog.close();

}









void calculate_aQGC_constraints_combined(std::vector<string> paranames, std::vector<string> channellabels, std::vector<std::vector<TH1D*>> histsamples_channels, std::vector<std::vector<int>> histsamplestypes_channels, std::vector<RooDataSet*> datasamples_channels,  double confidence_level=0.95, int ts_type=1, bool fit_for_ts=false, int altModel_ntoys=1000, int nullModel_ntoys=1000) {


	std::vector<double> exp_S_sm;
	std::vector<double> exp_B_sm;
	std::vector<double> exp_Total;
	std::vector<double> a0;
	std::vector<double> a1;
	std::vector<double> a2;
	std::vector<RooRealVar*> vec_common_observables;
	RooArgSet* common_observables = new RooArgSet();
	std::vector<tools_for_aqgc_constraints_subchannel*> subchannel_tools;

	////// Codes to make up poi for all channels
	RooWorkspace* combine_wspace = new RooWorkspace("combine_wspace");
	combine_wspace->factory( Form("common_poi_%s[0.0, -100.0, 100.0]",paranames[0].c_str()) );
	RooRealVar* common_poi = (RooRealVar*)(combine_wspace->var(Form("common_poi_%s",paranames[0].c_str())));
	if(TString(paranames[0]).Contains("S")) {
		common_poi->setRange(-500, 500);
	}
	else if(TString(paranames[0]).Contains("M")) {
		common_poi->setRange(-60, 60);
	}
	else if(TString(paranames[0]).Contains("M7")) {
		common_poi->setRange(-130, 130);
	}
	else if(TString(paranames[0]).Contains("T")) {
		common_poi->setRange(-25, 25);
	}
	////// Codes to make up observables for each channel
	for(int ichannel=0; ichannel<channellabels.size(); ichannel++) {
		combine_wspace->factory( Form("common_obs_%s_%s[5.0, 0.0, 1E5]",paranames[0].c_str(), channellabels[ichannel].c_str()));
	}
	for(int ichannel=0; ichannel<channellabels.size(); ichannel++) {
		( (RooRealVar*)(combine_wspace->var( Form("common_obs_%s_%s",paranames[0].c_str(), channellabels[ichannel].c_str()) )) )->Print("V");
		vec_common_observables.push_back( (RooRealVar*)(combine_wspace->var( Form("common_obs_%s_%s",paranames[0].c_str(), channellabels[ichannel].c_str()) )) );
		common_observables->add(*(vec_common_observables[ichannel]));
	}
	std::vector<std::vector<double>> temp_obs_values(channellabels.size());
	///// Codes to setup tools and PDFs for each channel
	for(int ichannel=0; ichannel<channellabels.size(); ichannel++) {
		subchannel_tools.push_back( new tools_for_aqgc_constraints_subchannel(paranames, channellabels[ichannel], common_poi, vec_common_observables[ichannel], histsamples_channels[ichannel], histsamplestypes_channels[ichannel]) );
		(subchannel_tools.back())->calculate_aqgcparavalues_vs_histintegrals();
		(subchannel_tools.back())->construct_func_RatioToSM_VS_aqgc_1D();
		(subchannel_tools.back())->draw_RatioToSM_VS_aqgc_1D_plot();
		exp_S_sm.push_back((subchannel_tools.back())->obtain_sm_signal_yields());
		exp_B_sm.push_back((subchannel_tools.back())->obtain_sm_background_yields());
		exp_Total.push_back(exp_S_sm.back()+exp_B_sm.back());
		cout<<"Subchannel [ "<<(channellabels[ichannel])<<" ] === "<<endl;
		cout<<" --- SM expected signal: "<<exp_S_sm.back()<<endl;
		cout<<" --- SM expected background: "<<exp_B_sm.back()<<endl;
		(subchannel_tools.back())->print_par_value_of_specific_significance(2.0);
		(subchannel_tools.back())->print_par_value_of_specific_significance(3.0);
		(subchannel_tools.back())->print_par_value_of_specific_significance(5.0);
		a0.push_back((subchannel_tools.back())->obtain_func_RatioToSM_VS_aqgc_1D_a0());
		a1.push_back((subchannel_tools.back())->obtain_func_RatioToSM_VS_aqgc_1D_a1());
		a2.push_back((subchannel_tools.back())->obtain_func_RatioToSM_VS_aqgc_1D_a2());
		combine_wspace->factory( Form("Bkg_%s[%e]", channellabels[ichannel].c_str(), exp_S_sm.back()) );
		RooRealVar *b = (RooRealVar*) combine_wspace->var(Form("Bkg_%s", channellabels[ichannel].c_str()));
		b->setConstant();
		combine_wspace->factory( Form("a0_%s[%e]" , channellabels[ichannel].c_str(), a0.back()) );
		combine_wspace->factory( Form("a1_%s[%e]" , channellabels[ichannel].c_str(), a1.back()) );
		combine_wspace->factory( Form("a2_%s[%e]" , channellabels[ichannel].c_str(), a2.back()) );
		combine_wspace->factory( Form("RooPolyVar::Sig_%s(%s,{a0_%s,a1_%s,a2_%s})", channellabels[ichannel].c_str(), Form("common_poi_%s", paranames[0].c_str()), channellabels[ichannel].c_str(), channellabels[ichannel].c_str(), channellabels[ichannel].c_str()) ); // signal counting expected in poisson distribution
		combine_wspace->factory( Form("sum::Total_%s(Sig_%s,Bkg_%s)", channellabels[ichannel].c_str(), channellabels[ichannel].c_str(), channellabels[ichannel].c_str()) );
		combine_wspace->factory( Form("RooPoisson::countingModel_%s(common_obs_%s_%s,Total_%s)", channellabels[ichannel].c_str(), paranames[0].c_str(), channellabels[ichannel].c_str(), channellabels[ichannel].c_str()) ); // counting model

		///// Prepare datasets (step 1)
		if(datasamples_channels[ichannel]==0 || (datasamples_channels[ichannel])->numEntries()<1) {
			//TRandom r(0);
			//(vec_common_observables[ichannel])->setVal(r.Poisson(exp_Total));
			//(vec_common_observables[ichannel])->setVal(exp_Total);
			temp_obs_values[ichannel].push_back(exp_Total.back());
		}
		else {
			for(int ie=0; ie<(datasamples_channels[ichannel])->numEntries(); ie++){
				if(((datasamples_channels[ichannel])->get(ie))!=0) {
					TIter itobs = ((datasamples_channels[ichannel])->get(ie))->createIterator();
					RooRealVar *riter;
					while((riter = dynamic_cast<RooRealVar*>(itobs.Next()))) {
						temp_obs_values[ichannel].push_back(riter->getVal());
					}
				}
			}
		}
	}
	combine_wspace->Print("V");

	///// Prepare datasets (step 2)
	RooDataSet *combine_datasample = new RooDataSet("combine_datasample", "combine_datasample", RooArgSet(*common_observables));
	cout<<"================== Combined Dataset"<<endl;
	for(int de=0; de<temp_obs_values[0].size(); de++) {
		for(int ichannel=0; ichannel<channellabels.size(); ichannel++) {
			((RooRealVar*)(common_observables->find(Form("common_obs_%s_%s",paranames[0].c_str(),channellabels[ichannel].c_str()))))->setVal(temp_obs_values[ichannel].at(de));
		}
		cout<<" --- entry: "<<de<<endl;
		common_observables->Print("V");
		combine_datasample->add(*common_observables);
	}
	combine_datasample->Print("V");
		

	///// Prepare PDFs (step 2)
	cout<<"================== Combined Workspace"<<endl;
	string combine_PDF_string = "PROD::countingModel(";
	for(int ichannel=0; ichannel<channellabels.size(); ichannel++) {
		if(ichannel!=(channellabels.size()-1)) combine_PDF_string += Form("countingModel_%s , ",channellabels[ichannel].c_str());
		else combine_PDF_string += Form("countingModel_%s",channellabels[ichannel].c_str());
	}
	combine_PDF_string += ")";
	//combine_PDF_string = Form("PROD::countingModel( RooPoisson(%s,Total_%s), RooPoisson(%s,Total_%s) )", vec_common_observables[0]->GetName(), "channel1", vec_common_observables[1]->GetName(), "channel2");
	combine_wspace->factory(Form("%s", combine_PDF_string.c_str()));
	combine_wspace->Print("V");


	cout<<"================== Data: "<<endl;
	combine_datasample->Print("V");
	cout<<"================== PDF:"<<endl;
	RooAbsPdf* countingModel = (RooAbsPdf*) combine_wspace->pdf("countingModel"); 
	countingModel->Print("V");
	cout<<"=================="<<endl;


	ModelConfig* aqgc_modelconfig;
	ModelConfig* sm_modelconfig;

	/// --- modelconfig for aQGC (null hypothesis)
	RooRealVar* POI = ((RooRealVar*)(combine_wspace->var(common_poi->GetName())));
	RooArgSet paramOfInterest(*POI);
	RooArgSet Observables;
	for(int ichannel=0; ichannel<channellabels.size(); ichannel++) {
		Observables.add(*(vec_common_observables[ichannel]));
	}
	aqgc_modelconfig = new ModelConfig(combine_wspace);
	aqgc_modelconfig->SetName(Form("aqgc_modelconfig"));
	aqgc_modelconfig->SetPdf(*countingModel);
	aqgc_modelconfig->SetParametersOfInterest(paramOfInterest);
	aqgc_modelconfig->SetObservables(Observables);
	aqgc_modelconfig->SetSnapshot(*(aqgc_modelconfig->GetParametersOfInterest()));
	//RooArgSet constrainedParams;
	//constrainedParams.add(*ratioBkg);
	//aqgc_modelconfig->SetNuisanceParameters(constrainedParams);
	//aqgc_modelconfig->SetGlobalObservables(RooArgSet(*globalBkgObsNum));
	/// ---
	/// --- modelconfig for sm (alternative hypothesis)
	sm_modelconfig = (ModelConfig*)(aqgc_modelconfig->Clone());
	sm_modelconfig->SetName(Form("sm_modelconfig"));      
	RooRealVar* testpoi = (RooRealVar*)(aqgc_modelconfig->GetParametersOfInterest()->first()); // for case of 1 POI
	double oldval = testpoi->getVal();
	testpoi->setVal(0); // Not sure if this should be force to 0
	sm_modelconfig->SetSnapshot(RooArgSet(*testpoi));
	testpoi->setVal(oldval);
	/// ---



	// Necessary to set proper range for scanning of different aQGC coefficients, otherwise the program may stack (the underlying reason is unknown)
	std::pair<double,double> scan_range;
	if(TString(paranames[0]).Contains("S")) {
		scan_range.first = -500;
		scan_range.second = 500;
	}
	else if(TString(paranames[0]).Contains("M")) {
		scan_range.first = -60;
		scan_range.second = 60;
	}
	else if(TString(paranames[0]).Contains("M7")) {
		scan_range.first = -130;
		scan_range.second = 130;
	}
	else if(TString(paranames[0]).Contains("T")) {
		scan_range.first = -25;
		scan_range.second = 25;
	}

	CombinedHypoTest(paranames, common_observables, common_poi, scan_range, aqgc_modelconfig, sm_modelconfig, combine_datasample, confidence_level, ts_type, fit_for_ts, altModel_ntoys, nullModel_ntoys, 0);

}



//==============================================================================//
//////////////////////////////////////////////////////////////////////////////////









void PlotEvents(std::vector<TString> filesin, std::vector<bool> if_event_exist, std::vector<int> sample_types, string plot_prefix, string cut, int ndivisions, string para_name, double RWstart, double RWstep, int RWnstep) {

	gStyle->SetOptStat(0000);
	for(int ifi=0; ifi<if_event_exist.size(); ifi++) {
		if(filesin.empty()) return;
		if( !(if_event_exist.at(ifi)) ) {
			filesin.erase(filesin.begin()+ifi);
			sample_types.erase(sample_types.begin()+ifi);
		}
	}
	if(filesin.empty()) return;
	int nfiles = filesin.size();
	TFile fin1(filesin.at(0), "READ");
	TTree* tin1 = (TTree*)fin1.Get("tree");
	TObjArray* objarr = tin1->GetListOfBranches();
	TObjArray* objarr1 = (TObjArray*)objarr->Clone();
	int objnum = objarr1->GetEntries();
	delete tin1;
	fin1.Close();

	TCanvas* c_hstacks = new TCanvas("c_hstacks","c_hstacks",2500,1900);
	double xmin=0;
	double xmax=100;
	double xrange=0;
	double width_step=0;
	double signalhistmaxy = 0;
	for(int k=0; k<objnum; k++){
		if(k>0) {
			cout<<"Clear Pad----------------------------------"<<endl;
			c_hstacks->Clear();
		}
		c_hstacks->cd();
		c_hstacks->SetRightMargin(0.25);
		TLegend lg(0.75,0.70,0.99,0.99);
		//lg.SetNColumns(2);
		cout<<"Draw "<<k<<" Hists--------------------------"<<endl;
		cout<<"Draw Hists---------------------------------- Branch Name: "<< ((objarr1->At(k)))->GetName() <<endl;

		std::vector<TH1D*> ewh_temps;
		std::vector<TH1D*> rwh_temps;
		std::vector<TH1D*> bmh_temps;
		std::vector<TH1D*> smh_temps;
		std::vector<TH1D*> rw_sub_sm_h_temps;
		std::map<int, string> benchmarks_map;
		std::map<int, double> benchmarks_values_map;
		xmin = 0;
		xmax = 100;
		xrange = xmax - xmin;
		width_step = xrange/ndivisions;
		signalhistmaxy = 0;
		for(int fid=0; fid<nfiles; fid++){ 
			if(sample_types.at(fid) == 1)continue;
			TChain chin("tree");
			chin.Add(filesin.at(fid));
			chin.SetBranchStatus("*", 1);
			chin.SetBranchStatus("SF", 1);
			if( (string(((objarr1->At(k)))->GetName()) != string("weights")) && (string(((objarr1->At(k)))->GetName()) != string("lep_fst")) ) {
				chin.Draw(Form("%s>>opt_temp%d", (objarr1->At(k))->GetName(), k), Form("10000000*SF*(event_weight)*((%s > -99) && %s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
				xmin = ((TH1F*)(gDirectory->Get(Form("opt_temp%d", k))))->GetBinLowEdge(1);
				xmax = ((TH1F*)(gDirectory->Get(Form("opt_temp%d", k))))->GetBinLowEdge( ((TH1F*)(gDirectory->Get(Form("opt_temp%d", k))))->GetNcells() );
				signalhistmaxy += ((TH1F*)(gDirectory->Get(Form("opt_temp%d", k))))->GetMaximum();
				xrange = xmax - xmin;
				width_step = xrange/ndivisions;
				xmin -= 2.0*width_step;
				xmax += 2.0*width_step;
			}
			else if( (string(((objarr1->At(k)))->GetName()) == string("weights")) && (string(((objarr1->At(k)))->GetName()) != string("lep_fst")) ) {
				chin.Draw(Form("%s>>opt_temp%d", (objarr1->At(k))->GetName(), k), Form("10000000*SF*(event_weight)*((%s > -99) && %s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
				xmin = ((TH1F*)(gDirectory->Get(Form("opt_temp%d", k))))->GetBinLowEdge(1);
				xmax = ((TH1F*)(gDirectory->Get(Form("opt_temp%d", k))))->GetBinLowEdge( ((TH1F*)(gDirectory->Get(Form("opt_temp%d", k))))->GetNcells() );
				signalhistmaxy += (chin.GetEntries());
				xrange = xmax - xmin;
				width_step = xrange/ndivisions;
				xmin -= 2.0*width_step;
				xmax += 2.0*width_step;
			}
			else if( (string(((objarr1->At(k)))->GetName()) == string("lep_fst")) ) {
				xmin = 0;
				xmax = 100;
				signalhistmaxy += (chin.GetEntries()/2);
				xrange = xmax - xmin;
				width_step = xrange/ndivisions;
				xmin -= 2.0*width_step;
				xmax += 2.0*width_step;
			}
		}
		cout<<"======================================================================== Decide histos ranges"<<endl;

		//-------------------------------------- histograms for aQGC sample with original event weights
		ewh_temps.push_back(new TH1D(Form("th_aqgcEW_%s_%s_%s","signal", para_name.c_str(), ((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
		//-------------------------------------- histograms for aQGC sample with reweighting weights
		for(int nw=0; nw<=RWnstep; nw++) {
			rwh_temps.push_back(new TH1D(Form("th_aqgcRW_%s_%s=%.3fe-12_%s","signal", para_name.c_str(), convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep), ((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
		}
		//-------------------------------------- histograms for aQGC benchmarks samples
		for(int fid=0; fid<nfiles; fid++){ 
			if( !((filesin.at(fid)).Contains("benchmark")) )continue;
			string benchmark((filesin.at(fid)).Data());
			benchmark.erase(0, (benchmark.find("_benchmark_",0) + 11));
			benchmark.erase(benchmark.find(Form("00.root"),0), string(Form("00.root")).size());
			benchmarks_map[fid] = benchmark;
			string tempbm = benchmark;
			tempbm.erase(0,3);
			benchmarks_values_map[fid] = atof(tempbm.c_str());
			bmh_temps.push_back(new TH1D(Form("th_benchmark%d_%s_%se-12_%s", fid, "signal", benchmarks_map[fid].c_str(), ((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
		}
		//-------------------------------------- histograms for SM samples
		{
			smh_temps.push_back(new TH1D(Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
			smh_temps.push_back(new TH1D(Form("th_SM_%s_%s","mumuTomumuX",((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
			smh_temps.push_back(new TH1D(Form("th_SM_%s_%s","mumuTovmvmX",((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
			smh_temps.push_back(new TH1D(Form("th_SM_%s_%s","mumuTottX",((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
			smh_temps.push_back(new TH1D(Form("th_SM_%s_%s","mumuTomultibosons",((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
			cout<<"histo vectors prepared"<<endl;
		}


		//============================================================================== Fill histograms
		for(int fid=0; fid<nfiles; fid++){ 
			TChain Chain("tree");
			Chain.Add(filesin.at(fid));
			Chain.SetBranchStatus("*", 1);
			Chain.SetBranchStatus("SF", 1);

			//------------------------------------------------- aQGC signal with reweighting
			if((filesin.at(fid)).Contains("reweightscan")) {
				if(string(((objarr1->At(k)))->GetName()) == string("weights")) { 
					Chain.Draw(Form("event_weight>>+%s", Form("th_aqgcEW_%s_%s_%s","signal", para_name.c_str(), ((objarr1->At(k)))->GetName())), Form("10000000*SF*(%s)",cut.c_str()), "goff");
					for(int nw=0; nw<=RWnstep; nw++) {
						Chain.Draw(Form("weights[%d]>>+%s", nw, Form("th_aqgcRW_%s_%s=%.3fe-12_%s","signal", para_name.c_str(), convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep), ((objarr1->At(k)))->GetName())), Form("10000000*SF*(%s)",cut.c_str()), "goff");
					}
				}
				if(string(((objarr1->At(k)))->GetName()) == string("lep_fst")) { 
					Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_aqgcEW_%s_%s_%s","signal", para_name.c_str(), ((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
					for(int nw=0; nw<=RWnstep; nw++) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_aqgcRW_%s_%s=%.3fe-12_%s","signal", para_name.c_str(), convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep), ((objarr1->At(k)))->GetName())), Form("10000000*SF*(weights[%d])*(%s)",nw,cut.c_str()), "goff");
					}
				}
				if( (string(((objarr1->At(k)))->GetName()) != string("lep_fst")) && (string(((objarr1->At(k)))->GetName()) != string("weights")) ) {
					Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_aqgcEW_%s_%s_%s","signal", para_name.c_str(), ((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
					for(int nw=0; nw<=RWnstep; nw++) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_aqgcRW_%s_%s=%.3fe-12_%s","signal", para_name.c_str(), convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep), ((objarr1->At(k)))->GetName())), Form("10000000*SF*(weights[%d])*(%s > -99)*(%s)",nw,((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
					}
				}
			}


			//------------------------------------------------- aQGC signal bench mark
			if( !((filesin.at(fid)).Contains("reweightscan")) && ((filesin.at(fid)).Contains("benchmark")) ) {
				if(string(((objarr1->At(k)))->GetName()) == string("weights")) { 
					Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_benchmark%d_%s_%se-12_%s", fid, "signal", benchmarks_map[fid].c_str(), ((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)", cut.c_str()), "goff");
				}
				if(string(((objarr1->At(k)))->GetName()) == string("lep_fst")) { 
					Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_benchmark%d_%s_%se-12_%s", fid, "signal", benchmarks_map[fid].c_str(), ((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)", cut.c_str()), "goff");
				}
				if( (string(((objarr1->At(k)))->GetName()) != string("lep_fst")) && (string(((objarr1->At(k)))->GetName()) != string("weights")) ) {
					Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_benchmark%d_%s_%se-12_%s", fid, "signal", benchmarks_map[fid].c_str(), ((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)", ((objarr1->At(k)))->GetName(), cut.c_str()), "goff");
				}
			}


			//------------------------------------------------- SM signal and backgrounds
			if( !((filesin.at(fid)).Contains("reweightscan")) && !((filesin.at(fid)).Contains("benchmark")) ) {
				if(string(((objarr1->At(k)))->GetName()) == string("lep_fst")) {
					if((filesin.at(fid)).Contains("mumuTozzzTo")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("mumuTomumu")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTomumuX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("mumuTovmvm")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTovmvmX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("tt")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTottX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("mumuTow") || (filesin.at(fid)).Contains("mumuToz") || (filesin.at(fid)).Contains("mumuToh") ) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTomultibosons",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
					}
				}
				if(string(((objarr1->At(k)))->GetName()) != string("lep_fst")) {
					if((filesin.at(fid)).Contains("mumuTozzzTo")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("mumuTomumu")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTomumuX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("mumuTovmvm")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTovmvmX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("tt")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTottX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("mumuTow") || (filesin.at(fid)).Contains("mumuToz") || (filesin.at(fid)).Contains("mumuToh") ) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTomultibosons",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
					}
				}
			}
		}

		//-------------------------------------- Set histograms for aQGC sample with reweighting weights subtracting the SM distribution
		for(int nw=0; nw<=RWnstep; nw++) {
			rw_sub_sm_h_temps.push_back( (TH1D*)( (rwh_temps.at(nw))->Clone(Form("th_aqgcRW-SM_%s_%s=%.3fe-12_%s","signal", para_name.c_str(), convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep), ((objarr1->At(k)))->GetName())) ) );
			(rw_sub_sm_h_temps.back())->Add(smh_temps.at(0), -1);
			for(int bini=0; bini<((rwh_temps.at(nw))->GetNcells()); bini++) {
				(rw_sub_sm_h_temps.back())->SetBinError(bini, ((rwh_temps.at(nw))->GetBinError(bini)));
			}
		}
		cout<<"======================================================================== Preapare hist plots"<<endl;


		c_hstacks->SetTopMargin(3.0);
		c_hstacks->cd();
		for(int hindex=0; hindex<int(ewh_temps.size()); hindex++) {
			if( !(string((ewh_temps.at(hindex))->GetName())==string(Form("th_aqgcEW_%s_%s_%s","signal", para_name.c_str(), ((objarr1->At(k)))->GetName()))) ) continue;
			cout<<"-------------------------------------------- Draw original AQGC signal original weighted plots"<<endl;
			cout<<"Integral = "<<(ewh_temps.at(hindex))->Integral()<<endl;
			cout<<"Mean = "<<(ewh_temps.at(hindex))->GetMean()<<endl;
			cout<<"RMS = "<<(ewh_temps.at(hindex))->GetRMS()<<endl;
			(ewh_temps.at(hindex))->SetTitle(Form(" ; %s; Events",((objarr1->At(k)))->GetName()));
			(ewh_temps.at(hindex))->SetTitleOffset(1.2);
			(ewh_temps.at(hindex))->SetLineWidth(2*3);
			(ewh_temps.at(hindex))->SetLineStyle(1);
			(ewh_temps.at(hindex))->SetMarkerStyle(20);
			(ewh_temps.at(hindex))->SetMarkerSize(0.9*3);
			(ewh_temps.at(hindex))->SetLineColor(2);
			(ewh_temps.at(hindex))->SetMarkerColor(2);
			//(ewh_temps.at(hindex))->SetFillColorAlpha(kBlue, 30);
			//(ewh_temps.at(hindex))->SetFillStyle(3001);
			(ewh_temps.at(hindex))->Draw("PL");
			(ewh_temps.at(hindex))->SetMaximum( signalhistmaxy*300.0 );
			string process((ewh_temps.at(hindex))->GetName());
			process.erase(process.find("th_",0), string("th_").size());
			process.erase(process.find(Form("_%s",((objarr1->At(k)))->GetName()),0), string(Form("_%s",((objarr1->At(k)))->GetName())).size());
			lg.AddEntry(ewh_temps.at(hindex), process.c_str(), "LP");
		}

		for(int hindex=0; hindex<int(bmh_temps.size()); hindex++) {
			cout<<"-------------------------------------------- Draw benchmark AQGC benchmark plots"<<endl;
			cout<<"Integral = "<<(bmh_temps.at(hindex))->Integral()<<endl;
			cout<<"Mean = "<<(bmh_temps.at(hindex))->GetMean()<<endl;
			cout<<"RMS = "<<(bmh_temps.at(hindex))->GetRMS()<<endl;
			(bmh_temps.at(hindex))->SetTitle(Form(" ; %s; Events",((objarr1->At(k)))->GetName()));
			(bmh_temps.at(hindex))->SetTitleOffset(1.2);
			(bmh_temps.at(hindex))->SetLineWidth(2*3);
			(bmh_temps.at(hindex))->SetLineStyle(2);
			(bmh_temps.at(hindex))->SetMarkerStyle(20+hindex);
			(bmh_temps.at(hindex))->SetMarkerSize(0.9*3);
			(bmh_temps.at(hindex))->SetLineColor(9+2*hindex);
			(bmh_temps.at(hindex))->SetMarkerColor(9+2*hindex);
			if( (9+2*hindex)==5 ){
				(bmh_temps.at(hindex))->SetLineColor(41+(hindex%3));
				(bmh_temps.at(hindex))->SetMarkerColor(41+(hindex%3));
			}
			(bmh_temps.at(hindex))->Draw("PL SAME");
			(bmh_temps.at(hindex))->SetMaximum( signalhistmaxy*300.0 );
			string process((bmh_temps.at(hindex))->GetName());
			process.erase(process.find("th_",0), string("th_").size());
			process.erase(process.find(Form("_%s",((objarr1->At(k)))->GetName()),0), string(Form("_%s",((objarr1->At(k)))->GetName())).size());
			lg.AddEntry(bmh_temps.at(hindex), process.c_str(), "LP");
		}

		for(int hindex=0; hindex<int(rwh_temps.size()); hindex++) {
			bool skip = true;
			for(int fid=0; fid<nfiles; fid++){
				if(fabs(convert_weight_index_to_paravalue(hindex, RWstart, RWstep, RWnstep)) == fabs(benchmarks_values_map[fid])) {
					skip = false;
				}
				//if((abs(int(/*1e12**/(convert_weight_index_to_paravalue(hindex, RWstart, RWstep, RWnstep)))-int(benchmarks_values_map[fid])) % abs(4*int(RWstep)) == 0)) {
				//	skip = false;
				//}
				//if(!skip){
				//	cout<<"(abs(int(/*1e12**/(convert_weight_index_to_paravalue(hindex, RWstart, RWstep, RWnstep)))-int(benchmarks_values_map[fid])) // abs(4*int(RWstep)))"<<(abs(int((convert_weight_index_to_paravalue(hindex, RWstart, RWstep, RWnstep)))-int(benchmarks_values_map[fid])) % abs(4*int(RWstep)))<<endl;
				//	cout<<"int(convert_weight_index_to_paravalue(hindex, RWstart, RWstep, RWnstep) = "<<int(convert_weight_index_to_paravalue(hindex, RWstart, RWstep, RWnstep))<<endl;
				//	cout<<"int(benchmarks_values_map[fid])) = "<<int(benchmarks_values_map[fid])<<endl;
				//	cout<<"abs(4*int(RWstep) = "<<abs(4*int(RWstep))<<endl;
				//}
			}
			if(fabs(convert_weight_index_to_paravalue(hindex, RWstart, RWstep, RWnstep)) == fabs(0) || fabs(convert_weight_index_to_paravalue(hindex, RWstart, RWstep, RWnstep)) == fabs(50)) {
				skip = false;
			}
			if(skip) continue;
			cout<<"-------------------------------------------- Draw reweighted AQGC signal reweighted plots"<<endl;
			cout<<"REWEIGHT convert_weight_index_to_paravalue = "<<fabs(convert_weight_index_to_paravalue(hindex, RWstart, RWstep, RWnstep))<<endl;
			cout<<"Integral = "<<(rwh_temps.at(hindex))->Integral()<<endl;
			cout<<"Mean = "<<(rwh_temps.at(hindex))->GetMean()<<endl;
			cout<<"RMS = "<<(rwh_temps.at(hindex))->GetRMS()<<endl;
			if(string((rwh_temps.at(hindex))->GetName())!=string(Form("th_aqgcEW_%s_%s_%s","signal", para_name.c_str(), ((objarr1->At(k)))->GetName()))) {
				(rwh_temps.at(hindex))->SetTitle(Form(" ; %s; Events",((objarr1->At(k)))->GetName()));
				(rwh_temps.at(hindex))->SetTitleOffset(1.2);
				(rwh_temps.at(hindex))->SetLineWidth(2*3);
				(rwh_temps.at(hindex))->SetLineStyle(7);
				(rwh_temps.at(hindex))->SetMarkerStyle(21+(hindex%11));//(68+hindex);
				(rwh_temps.at(hindex))->SetMarkerSize(0.9*3);
				(rwh_temps.at(hindex))->SetLineColor(3+(hindex%7));
				(rwh_temps.at(hindex))->SetMarkerColor(3+(hindex%7));
				if( (3+(hindex%7))==5 ){
					(rwh_temps.at(hindex))->SetLineColor(41+(hindex%3));
					(rwh_temps.at(hindex))->SetMarkerColor(41+(hindex%3));
				}
				(rwh_temps.at(hindex))->Draw("PL SAME");
				(rwh_temps.at(hindex))->SetMaximum( signalhistmaxy*300.0 );
				string process((rwh_temps.at(hindex))->GetName());
				process.erase(process.find("th_",0), string("th_").size());
				process.erase(process.find(Form("_%s",((objarr1->At(k)))->GetName()),0), string(Form("_%s",((objarr1->At(k)))->GetName())).size());
				lg.AddEntry(rwh_temps.at(hindex), process.c_str(), "LP");
			}
		}

		for(int hindex=0; hindex<int(smh_temps.size()); hindex++) {
			if(string((smh_temps.at(hindex))->GetName())==string(Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName()))) {
				string process((smh_temps.at(hindex))->GetName());
				process.erase(process.find("th_",0), string("th_").size());
				process.erase(process.find(Form("_%s",((objarr1->At(k)))->GetName()),0), string(Form("_%s",((objarr1->At(k)))->GetName())).size());
				cout<<"-------------------------------------------- Draw SM signal["<<process<<"] plots"<<endl;
				cout<<"Integral = "<<(smh_temps.at(hindex))->Integral()<<endl;
				cout<<"Mean = "<<(smh_temps.at(hindex))->GetMean()<<endl;
				cout<<"RMS = "<<(smh_temps.at(hindex))->GetRMS()<<endl;
				(smh_temps.at(hindex))->SetTitle(Form(" ; %s; Events",((objarr1->At(k)))->GetName()));
				(smh_temps.at(hindex))->SetTitleOffset(1.2);
				(smh_temps.at(hindex))->SetLineWidth(2*3);
				(smh_temps.at(hindex))->SetLineStyle(1);
				(smh_temps.at(hindex))->SetMarkerStyle(20);
				(smh_temps.at(hindex))->SetMarkerSize(0.9*3);
				(smh_temps.at(hindex))->SetLineColor(kBlack);
				(smh_temps.at(hindex))->SetMarkerColor(kBlack);
				//(smh_temps.at(hindex))->SetFillColorAlpha(kBlack, 10);
				//(smh_temps.at(hindex))->SetFillStyle(3001);
				(smh_temps.at(hindex))->Draw("HIST SAME");
				(smh_temps.at(hindex))->SetMaximum( signalhistmaxy*10.0 );
				lg.AddEntry(smh_temps.at(hindex), process.c_str(), "L");
			}
		}
		for(int hindex=0; hindex<int(smh_temps.size()); hindex++) {
			if(string((smh_temps.at(hindex))->GetName())!=string(Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName()))) {
				string process((smh_temps.at(hindex))->GetName());
				process.erase(process.find("th_",0), string("th_").size());
				process.erase(process.find(Form("_%s",((objarr1->At(k)))->GetName()),0), string(Form("_%s",((objarr1->At(k)))->GetName())).size());
				cout<<"-------------------------------------------- Draw SM background["<<process<<"] plots"<<endl;
				cout<<"Integral = "<<(smh_temps.at(hindex))->Integral()<<endl;
				cout<<"Mean = "<<(smh_temps.at(hindex))->GetMean()<<endl;
				cout<<"RMS = "<<(smh_temps.at(hindex))->GetRMS()<<endl;
				(smh_temps.at(hindex))->SetTitle(Form(" ; %s; Events",((objarr1->At(k)))->GetName()));
				(smh_temps.at(hindex))->SetTitleOffset(1.2);
				(smh_temps.at(hindex))->SetLineWidth(2*3);
				(smh_temps.at(hindex))->SetLineStyle(hindex);
				(smh_temps.at(hindex))->SetMarkerStyle(20+hindex);
				(smh_temps.at(hindex))->SetMarkerSize(0.9*3);
				(smh_temps.at(hindex))->SetLineColor(1+hindex);
				(smh_temps.at(hindex))->SetMarkerColor(1+hindex);
				if( (1+hindex)==5 ){
					(smh_temps.at(hindex))->SetLineColor(41+(hindex%3));
					(smh_temps.at(hindex))->SetMarkerColor(41+(hindex%3));
				}
				//(smh_temps.at(hindex))->SetFillColorAlpha(1+hindex, 10);
				//(smh_temps.at(hindex))->SetFillStyle(3001);
				(smh_temps.at(hindex))->Draw("HIST SAME");
				(smh_temps.at(hindex))->SetMaximum( signalhistmaxy*10.0 );
				lg.AddEntry(smh_temps.at(hindex), process.c_str(), "LF");
			}
		}
		lg.Draw("SAME");

		for(int hindex=0; hindex<int(rw_sub_sm_h_temps.size()); hindex++) {
			if(string((rw_sub_sm_h_temps.at(hindex))->GetName())==string(Form("th_aqgcRW-SM_%s_%s=%.3fe-12_%s","signal", para_name.c_str(), convert_weight_index_to_paravalue(hindex, RWstart, RWstep, RWnstep), ((objarr1->At(k)))->GetName()))) {
				string process((rw_sub_sm_h_temps.at(hindex))->GetName());
				process.erase(process.find("th_",0), string("th_").size());
				process.erase(process.find(Form("_%s",((objarr1->At(k)))->GetName()),0), string(Form("_%s",((objarr1->At(k)))->GetName())).size());
				cout<<"-------------------------------------------- Draw AQGC-SM["<<process<<"] plots"<<endl;
				cout<<"Integral = "<<(rw_sub_sm_h_temps.at(hindex))->Integral()<<endl;
				cout<<"Mean = "<<(rw_sub_sm_h_temps.at(hindex))->GetMean()<<endl;
				cout<<"RMS = "<<(rw_sub_sm_h_temps.at(hindex))->GetRMS()<<endl;
			}
		}

		cout<<"======================================================================== Draw histos"<<endl;
		c_hstacks->SetLogy();
		c_hstacks->Update();
		cout<<"Save drawing----------------------------------"<<endl;
		if( access(Form("./plots/"),0) ) mkdir(Form("./plots/"),0755);
		if( access(Form("./plots/%s/", para_name.c_str()),0) ) mkdir(Form("./plots/%s/", para_name.c_str()),0755);
		c_hstacks->SaveAs(Form("./plots/%s/%s_%s_comparison_%s.png", para_name.c_str(), para_name.c_str(), plot_prefix.c_str(), ((objarr1->At(k)))->GetName()));
		c_hstacks->SaveAs(Form("./plots/%s/%s_%s_comparison_%s.pdf", para_name.c_str(), para_name.c_str(), plot_prefix.c_str(), ((objarr1->At(k)))->GetName()));


		cout<<"======================================================================== Prepare histos stacks"<<endl;
		cout<<"Clear pad for hstack----------------------------------"<<endl;
		c_hstacks->Clear();
		c_hstacks->cd();
		//c_hstacks->SetFrameFillColor(19);
		c_hstacks->SetTopMargin(3.0);
		c_hstacks->SetRightMargin(0.25);
		THStack *hs_sm = new THStack("hs_sm",Form(" ; %s; ",((objarr1->At(k)))->GetName())); 
		THStack *hs_aqgc = new THStack("hs_aqgc",Form(" ; %s; ",((objarr1->At(k)))->GetName()));
		std::vector<TH1D*> rw_add_sm_h_temps;
		for(int nw=0; nw<=RWnstep; nw++) {
			rw_add_sm_h_temps.push_back( (TH1D*)( (rwh_temps.at(nw))->Clone(Form("th_aqgcRW+SM_%s_%s=%.3fe-12_%s","background", para_name.c_str(), convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep), ((objarr1->At(k)))->GetName())) ) );
			for(int hindex=0; hindex<int(smh_temps.size()); hindex++) {
				if(string((smh_temps.at(hindex))->GetName())!=string(Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName()))) {
					(rw_add_sm_h_temps.back())->Add((smh_temps.at(hindex)), 1.0);
				}
			}
			for(int bini=0; bini<((rwh_temps.at(nw))->GetNcells()); bini++) {
				(rw_add_sm_h_temps.back())->SetBinError(bini, ((rwh_temps.at(nw))->GetBinError(bini)));
			}
		}
		TLegend lg_aqgc(0.75,0.60,0.99,0.80,"AQGC");
		for(int nw=0; nw<=RWnstep; nw++) {
			if(string((rw_add_sm_h_temps.at(nw))->GetName())==string(Form("th_aqgcRW+SM_%s_%s=%.3fe-12_%s","background", para_name.c_str(), convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep), ((objarr1->At(k)))->GetName()))) {
				bool skip = true;
				//if(fabs(convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep)) == fabs(30) || fabs(convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep)) == fabs(70) || fabs(convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep)) == fabs(100)) 
				if(fabs(convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep)) == fabs(0) || fabs(convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep)) == fabs(10) || fabs(convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep)) == fabs(30)) {
					skip = false;
				}
				if(skip) continue;
				if(string((rw_add_sm_h_temps.at(nw))->GetName())==string(Form("th_aqgcRW+SM_%s_%s=%.3fe-12_%s","background", para_name.c_str(), convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep), ((objarr1->At(k)))->GetName()))) {
					(rw_add_sm_h_temps.at(nw))->SetTitle(Form(" ; %s; Events",((objarr1->At(k)))->GetName()));
					(rw_add_sm_h_temps.at(nw))->SetTitleOffset(1.2);
					(rw_add_sm_h_temps.at(nw))->SetLineWidth(2*3);
					(rw_add_sm_h_temps.at(nw))->SetLineStyle(7);
					(rw_add_sm_h_temps.at(nw))->SetMarkerStyle(20+(nw));
					(rw_add_sm_h_temps.at(nw))->SetMarkerSize(0.9*3);
					(rw_add_sm_h_temps.at(nw))->SetLineColor(5+(nw%5));
					(rw_add_sm_h_temps.at(nw))->SetMarkerColor(5+(nw%5));
					if( (5+(nw%5))==5 ){
						(rw_add_sm_h_temps.at(nw))->SetLineColor(41+(nw%3));
						(rw_add_sm_h_temps.at(nw))->SetMarkerColor(41+(nw%3));
					}
					(rw_add_sm_h_temps.at(nw))->Draw("PL SAME");
					(rw_add_sm_h_temps.at(nw))->SetMaximum( signalhistmaxy*300.0 );
					string process((rw_add_sm_h_temps.at(nw))->GetName());
					process.erase(process.find("th_aqgcRW+SM_",0), string("th_aqgcRW+SM_").size());
					process.erase(process.find(Form("_%s",((objarr1->At(k)))->GetName()),0), string(Form("_%s",((objarr1->At(k)))->GetName())).size());
					lg_aqgc.AddEntry(rw_add_sm_h_temps.at(nw), process.c_str(), "LP");
				}
			}
		}
		TLegend lg_sm(0.75,0.80,0.99,0.99,"SM");
		// Setup SM signal legend
		for(int hindex=0; hindex<int(smh_temps.size()); hindex++) {
			if(string((smh_temps.at(hindex))->GetName())==string(Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName()))) {
				(smh_temps.at(hindex))->SetLineWidth(2*3);
				(smh_temps.at(hindex))->SetLineStyle(1);
				(smh_temps.at(hindex))->SetMarkerStyle(0);
				(smh_temps.at(hindex))->SetMarkerSize(0);
				(smh_temps.at(hindex))->SetLineColor(kBlack);
				(smh_temps.at(hindex))->SetMarkerColor(0);
				//(smh_temps.at(hindex))->SetFillColorAlpha(0, 10);
				//(smh_temps.at(hindex))->SetFillStyle(3001);
				////// Only for legen, not for stack				hs_sm->Add((smh_temps.at(hindex)));
				string process((smh_temps.at(hindex))->GetName());
				process.erase(process.find("th_",0), string("th_").size());
				process.erase(process.find(Form("_%s",((objarr1->At(k)))->GetName()),0), string(Form("_%s",((objarr1->At(k)))->GetName())).size());
				lg_sm.AddEntry(smh_temps.at(hindex), process.c_str(), "L");
			}
		}
		// Setup SM background legend
		for(int hindex=0; hindex<int(smh_temps.size()); hindex++) {
			if(string((smh_temps.at(hindex))->GetName())!=string(Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName()))) {
				(smh_temps.at(hindex))->SetLineWidth(1*3);
				(smh_temps.at(hindex))->SetLineStyle(1);
				(smh_temps.at(hindex))->SetMarkerStyle(0);
				(smh_temps.at(hindex))->SetMarkerSize(0);
				(smh_temps.at(hindex))->SetLineColor(1+hindex);
				(smh_temps.at(hindex))->SetMarkerColor(0);
				(smh_temps.at(hindex))->SetFillColorAlpha(1+hindex, 50);
				if( (1+hindex)==5 ){
					(smh_temps.at(hindex))->SetLineColor(41+(hindex%3));
					(smh_temps.at(hindex))->SetFillColorAlpha(41+(hindex%3), 50);
				}
				(smh_temps.at(hindex))->SetFillStyle(3001);
				////// Only for legen, not for stack				hs_sm->Add((smh_temps.at(hindex)));
				string process((smh_temps.at(hindex))->GetName());
				process.erase(process.find("th_",0), string("th_").size());
				process.erase(process.find(Form("_%s",((objarr1->At(k)))->GetName()),0), string(Form("_%s",((objarr1->At(k)))->GetName())).size());
				lg_sm.AddEntry(smh_temps.at(hindex), process.c_str(), "LF");
			}
		}
		// Add SM background into stack
		for(int hindex=0; hindex<int(smh_temps.size()); hindex++) {
			if(string((smh_temps.at(hindex))->GetName())!=string(Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName()))) {
				hs_sm->Add((smh_temps.at(hindex)));
			}
		}
		// Add SM signal into stack
		for(int hindex=0; hindex<int(smh_temps.size()); hindex++) {
			if(string((smh_temps.at(hindex))->GetName())==string(Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName()))) {
				hs_sm->Add((smh_temps.at(hindex)));
			}
		}
		// Add AQGC+SM_background into stack
		for(int nw=0; nw<=RWnstep; nw++) {
			bool skip = true;
			if(fabs(convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep)) == fabs(0) || fabs(convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep)) == fabs(10) || fabs(convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep)) == fabs(30)) {
				skip = false;
			}
			if(skip) continue;
			if(string((rw_add_sm_h_temps.at(nw))->GetName())==string(Form("th_aqgcRW+SM_%s_%s=%.3fe-12_%s","background", para_name.c_str(), convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep), ((objarr1->At(k)))->GetName()))) {
				hs_aqgc->Add((rw_add_sm_h_temps.at(nw)));
			}
		}
		cout<<"======================================================================== Draw histos stacks"<<endl;
		hs_aqgc->Draw("ELP NOSTACK");
		hs_sm->Draw("HIST SAME");
		lg_sm.Draw("SAME");
		lg_aqgc.Draw("SAME");
		c_hstacks->SetLogy();
		c_hstacks->Update();
		cout<<"Save drawing----------------------------------"<<endl;
		if( access(Form("./plots/"),0) ) mkdir(Form("./plots/"),0755);
		if( access(Form("./plots/%s/", para_name.c_str()),0) ) mkdir(Form("./plots/%s/", para_name.c_str()),0755);
		c_hstacks->SaveAs(Form("./plots/%s/%s_%s_stacks_%s.png", para_name.c_str(), para_name.c_str(), plot_prefix.c_str(), ((objarr1->At(k)))->GetName()));
		c_hstacks->SaveAs(Form("./plots/%s/%s_%s_stacks_%s.pdf", para_name.c_str(), para_name.c_str(), plot_prefix.c_str(), ((objarr1->At(k)))->GetName()));



	}
}








std::pair<std::vector<TH1D*>, std::vector<int>> obtain_histograms(std::vector<TString> filesin, std::vector<bool> if_event_exist, std::vector<int> sample_types, string plot_prefix, string cut, int ndivisions, string para_name, double RWstart, double RWstep, int RWnstep, string Xvar="M4l") {

	std::vector<TH1D*> output_histograms;
	std::vector<int> output_histograms_types;
	std::pair<std::vector<TH1D*>, std::vector<int>> output;


	assert(!(filesin.empty()));
	gStyle->SetOptStat(0000);
	for(int ifi=0; ifi<if_event_exist.size(); ifi++) {
		if( !(if_event_exist.at(ifi)) ) {
			filesin.erase(filesin.begin()+ifi);
			sample_types.erase(sample_types.begin()+ifi);
		}
	}
	int nfiles = filesin.size();
	TFile fin1(filesin.at(0), "READ");
	TTree* tin1 = (TTree*)fin1.Get("tree");
	TObjArray* objarr = tin1->GetListOfBranches();
	TObjArray* objarr1 = (TObjArray*)objarr->Clone();
	int objnum = objarr1->GetEntries();
	delete tin1;
	fin1.Close();

	TCanvas* c_hstacks = new TCanvas("c_hstacks","c_hstacks",2500,1900);
	double xmin=0;
	double xmax=100;
	double xrange=0;
	double width_step=0;
	double signalhistmaxy = 0;
	for(int k=0; k<objnum; k++){
		if(string(((objarr1->At(k)))->GetName()) != string(Xvar)) continue;
		if(k>0) {
			cout<<"Clear Pad----------------------------------"<<endl;
			c_hstacks->Clear();
		}
		c_hstacks->cd();
		c_hstacks->SetRightMargin(0.25);
		TLegend lg(0.75,0.70,0.99,0.99);
		//lg.SetNColumns(2);
		cout<<"Draw "<<k<<" Hists--------------------------"<<endl;
		cout<<"Draw Hists---------------------------------- Branch Name: "<< ((objarr1->At(k)))->GetName() <<endl;

		std::vector<TH1D*> ewh_temps;
		std::vector<TH1D*> rwh_temps;
		std::vector<TH1D*> bmh_temps;
		std::vector<TH1D*> smh_temps;
		std::vector<TH1D*> rw_sub_sm_h_temps;
		std::map<int, string> benchmarks_map;
		std::map<int, double> benchmarks_values_map;
		xmin = 0;
		xmax = 100;
		xrange = xmax - xmin;
		width_step = xrange/ndivisions;
		signalhistmaxy = 0;
		for(int fid=0; fid<nfiles; fid++){ 
			if(sample_types.at(fid) == 1)continue;
			TChain chin("tree");
			chin.Add(filesin.at(fid));
			chin.SetBranchStatus("*", 1);
			chin.SetBranchStatus("SF", 1);
			if( (string(((objarr1->At(k)))->GetName()) != string("weights")) && (string(((objarr1->At(k)))->GetName()) != string("lep_fst")) ) {
				chin.Draw(Form("%s>>opt_temp%d", (objarr1->At(k))->GetName(), k), Form("10000000*SF*(event_weight)*((%s > -99) && %s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
				xmin = ((TH1F*)(gDirectory->Get(Form("opt_temp%d", k))))->GetBinLowEdge(1);
				xmax = ((TH1F*)(gDirectory->Get(Form("opt_temp%d", k))))->GetBinLowEdge( ((TH1F*)(gDirectory->Get(Form("opt_temp%d", k))))->GetNcells() );
				signalhistmaxy += ((TH1F*)(gDirectory->Get(Form("opt_temp%d", k))))->GetMaximum();
				xrange = xmax - xmin;
				width_step = xrange/ndivisions;
				xmin -= 2.0*width_step;
				xmax += 2.0*width_step;
			}
			else if( (string(((objarr1->At(k)))->GetName()) == string("weights")) && (string(((objarr1->At(k)))->GetName()) != string("lep_fst")) ) {
				chin.Draw(Form("%s>>opt_temp%d", (objarr1->At(k))->GetName(), k), Form("10000000*SF*(event_weight)*((%s > -99) && %s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
				xmin = ((TH1F*)(gDirectory->Get(Form("opt_temp%d", k))))->GetBinLowEdge(1);
				xmax = ((TH1F*)(gDirectory->Get(Form("opt_temp%d", k))))->GetBinLowEdge( ((TH1F*)(gDirectory->Get(Form("opt_temp%d", k))))->GetNcells() );
				signalhistmaxy += (chin.GetEntries());
				xrange = xmax - xmin;
				width_step = xrange/ndivisions;
				xmin -= 2.0*width_step;
				xmax += 2.0*width_step;
			}
			else if( (string(((objarr1->At(k)))->GetName()) == string("lep_fst")) ) {
				xmin = 0;
				xmax = 100;
				signalhistmaxy += (chin.GetEntries()/2);
				xrange = xmax - xmin;
				width_step = xrange/ndivisions;
				xmin -= 2.0*width_step;
				xmax += 2.0*width_step;
			}
		}
		cout<<"======================================================================== Decide histos ranges"<<endl;

		//-------------------------------------- histograms for aQGC sample with original event weights
		ewh_temps.push_back(new TH1D(Form("th_aqgcEW_%s_%s_%s","signal", para_name.c_str(), ((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
		//-------------------------------------- histograms for aQGC sample with reweighting weights
		for(int nw=0; nw<=RWnstep; nw++) {
			rwh_temps.push_back(new TH1D(Form("th_aqgcRW_%s_%s=%.3fe-12_%s","signal", para_name.c_str(), convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep), ((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
		}
		//-------------------------------------- histograms for aQGC benchmarks samples
		for(int fid=0; fid<nfiles; fid++){ 
			if( !((filesin.at(fid)).Contains("benchmark")) )continue;
			string benchmark((filesin.at(fid)).Data());
			benchmark.erase(0, (benchmark.find("_benchmark_",0) + 11));
			benchmark.erase(benchmark.find(Form("00.root"),0), string(Form("00.root")).size());
			benchmarks_map[fid] = benchmark;
			string tempbm = benchmark;
			tempbm.erase(0,3);
			benchmarks_values_map[fid] = atof(tempbm.c_str());
			bmh_temps.push_back(new TH1D(Form("th_benchmark%d_%s_%se-12_%s", fid, "signal", benchmarks_map[fid].c_str(), ((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
		}
		//-------------------------------------- histograms for SM samples
		{
			smh_temps.push_back(new TH1D(Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
			smh_temps.push_back(new TH1D(Form("th_SM_%s_%s","mumuTomumuX",((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
			smh_temps.push_back(new TH1D(Form("th_SM_%s_%s","mumuTovmvmX",((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
			smh_temps.push_back(new TH1D(Form("th_SM_%s_%s","mumuTottX",((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
			smh_temps.push_back(new TH1D(Form("th_SM_%s_%s","mumuTomultibosons",((objarr1->At(k)))->GetName()), Form("%s",((objarr1->At(k)))->GetName()), ndivisions+4, xmin, xmax));
			cout<<"histo vectors prepared"<<endl;
		}


		//============================================================================== Fill histograms
		for(int fid=0; fid<nfiles; fid++){ 
			TChain Chain("tree");
			Chain.Add(filesin.at(fid));
			Chain.SetBranchStatus("*", 1);
			Chain.SetBranchStatus("SF", 1);

			//------------------------------------------------- aQGC signal with reweighting
			if((filesin.at(fid)).Contains("reweightscan")) {
				if(string(((objarr1->At(k)))->GetName()) == string("weights")) { 
					Chain.Draw(Form("event_weight>>+%s", Form("th_aqgcEW_%s_%s_%s","signal", para_name.c_str(), ((objarr1->At(k)))->GetName())), Form("10000000*SF*(%s)",cut.c_str()), "goff");
					for(int nw=0; nw<=RWnstep; nw++) {
						Chain.Draw(Form("weights[%d]>>+%s", nw, Form("th_aqgcRW_%s_%s=%.3fe-12_%s","signal", para_name.c_str(), convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep), ((objarr1->At(k)))->GetName())), Form("10000000*SF*(%s)",cut.c_str()), "goff");
					}
				}
				if(string(((objarr1->At(k)))->GetName()) == string("lep_fst")) { 
					Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_aqgcEW_%s_%s_%s","signal", para_name.c_str(), ((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
					for(int nw=0; nw<=RWnstep; nw++) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_aqgcRW_%s_%s=%.3fe-12_%s","signal", para_name.c_str(), convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep), ((objarr1->At(k)))->GetName())), Form("10000000*SF*(weights[%d])*(%s)",nw,cut.c_str()), "goff");
					}
				}
				if( (string(((objarr1->At(k)))->GetName()) != string("lep_fst")) && (string(((objarr1->At(k)))->GetName()) != string("weights")) ) {
					Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_aqgcEW_%s_%s_%s","signal", para_name.c_str(), ((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
					for(int nw=0; nw<=RWnstep; nw++) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_aqgcRW_%s_%s=%.3fe-12_%s","signal", para_name.c_str(), convert_weight_index_to_paravalue(nw, RWstart, RWstep, RWnstep), ((objarr1->At(k)))->GetName())), Form("10000000*SF*(weights[%d])*(%s > -99)*(%s)",nw,((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
					}
				}
			}


			//------------------------------------------------- aQGC signal bench mark
			if( !((filesin.at(fid)).Contains("reweightscan")) && ((filesin.at(fid)).Contains("benchmark")) ) {
				if(string(((objarr1->At(k)))->GetName()) == string("weights")) { 
					Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_benchmark%d_%s_%se-12_%s", fid, "signal", benchmarks_map[fid].c_str(), ((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)", cut.c_str()), "goff");
				}
				if(string(((objarr1->At(k)))->GetName()) == string("lep_fst")) { 
					Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_benchmark%d_%s_%se-12_%s", fid, "signal", benchmarks_map[fid].c_str(), ((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)", cut.c_str()), "goff");
				}
				if( (string(((objarr1->At(k)))->GetName()) != string("lep_fst")) && (string(((objarr1->At(k)))->GetName()) != string("weights")) ) {
					Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_benchmark%d_%s_%se-12_%s", fid, "signal", benchmarks_map[fid].c_str(), ((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)", ((objarr1->At(k)))->GetName(), cut.c_str()), "goff");
				}
			}


			//------------------------------------------------- SM signal and backgrounds
			if( !((filesin.at(fid)).Contains("reweightscan")) && !((filesin.at(fid)).Contains("benchmark")) ) {
				if(string(((objarr1->At(k)))->GetName()) == string("lep_fst")) {
					if((filesin.at(fid)).Contains("mumuTozzzTo")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("mumuTomumu")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTomumuX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("mumuTovmvm")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTovmvmX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("tt")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTottX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("mumuTow") || (filesin.at(fid)).Contains("mumuToz") || (filesin.at(fid)).Contains("mumuToh") ) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTomultibosons",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s)",cut.c_str()), "goff");
					}
				}
				if(string(((objarr1->At(k)))->GetName()) != string("lep_fst")) {
					if((filesin.at(fid)).Contains("mumuTozzzTo")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","signal",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("mumuTomumu")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTomumuX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("mumuTovmvm")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTovmvmX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("tt")) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTottX",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
					}
					else if((filesin.at(fid)).Contains("mumuTow") || (filesin.at(fid)).Contains("mumuToz") || (filesin.at(fid)).Contains("mumuToh") ) {
						Chain.Draw(Form("%s>>+%s", ((objarr1->At(k)))->GetName(), Form("th_SM_%s_%s","mumuTomultibosons",((objarr1->At(k)))->GetName())), Form("10000000*SF*(event_weight)*(%s > -99)*(%s)",((objarr1->At(k)))->GetName(),cut.c_str()), "goff");
					}
				}
			}
		}

		for(int hi=0; hi<rwh_temps.size(); hi++){
			output_histograms.push_back(rwh_temps[hi]);
			output_histograms_types.push_back(2);
		}

		for(int hi=0; hi<smh_temps.size(); hi++){
			output_histograms.push_back(smh_temps[hi]);
			if((TString((output_histograms.back())->GetName())).Contains("signal")) output_histograms_types.push_back(0);
			else output_histograms_types.push_back(1);
		}

		output.first = output_histograms;
		output.second = output_histograms_types;

	}

	//cout<<"===== Summary of output histograms obtained ====="<<endl;
	//for(int hi=0; hi<(output.first).size(); hi++){
	//	cout<<"Histogram name: "<<(output.first)[hi]->GetName()<<endl;
	//	cout<<"Histogram integral: "<<(output.first)[hi]->Integral(1, (output.first)[hi]->GetNcells()-2)<<endl;
	//}
	//cout<<"===== Summary end ==============================="<<endl;

	return output;
}







void muoncolliderdelphes_selection_draw(string para_name, double RWstart, double RWstep, int RWnstep, string hardcut="", double CLv=0.95, int TStype=3){

	gROOT->SetBatch();

	cout<<"================================= Run Start ================================="<<endl;
	cout<<"AQGC dim-8 term: "<<para_name<<endl;
	cout<<"Reweighting start from "<<RWstart<<"E-12"<<endl;
	cout<<"Reweighting step length is "<<RWstep<<"E-12"<<endl;
	cout<<"Reweighting iterate "<<RWnstep<<" times"<<endl;
	cout<<"================================== Running =================================="<<endl;

	std::vector<string> para_names;
	para_names.push_back(para_name);	




	//--------------------------------------------------------------------------
	//------------------------ Analysis for Channel 1st ------------------------
	//--------------------------------------------------------------------------

	std::vector<TString> Infiles_names_channel1;
	std::vector<TString> Outfiles_names_preselected_channel1;
	std::vector<bool>    If_Event_exist_channel1;
	std::vector<int>     Sample_Types_channel1; // 0 for SM signal, 1 for SM background, 2 for AQGC signal

	cout<<"--------------------------------------------------------------------------"<<endl;
	cout<<"------------------------ Analysis for Channel 1st ------------------------"<<endl;
	cout<<"--------------------------------------------------------------------------"<<endl;

	string hardcut_channel1 = "( (Mrecoil>0 && Mrecoil<200) && (Mll1>80 && Mll1<100) && (Mll2>80 && Mll2<100) && (deltaM<20) )";

	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTohtt_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTomumutt_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTomumuwwTomumulvllvl_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTomumuwwz_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTomumuzTomumull_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTomumuzhTomumullh_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTomumuzhTomumuvlvlh_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTomumuzh_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTomumuzzTomumullvlvl_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTovmvmtt_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTovmvmwwh_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTovmvmwwzTovmvmlvllvll_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTovmvmzhh_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTovmvmzzh_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTovmvmzzz_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTovmvmzzTovmvmllll_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTowwhTolvllvlh_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTowwttTolvllvltt_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTowwzTolvllvlll_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTozhhTollhh_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTozhhTovlvlX_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTozttTolltt_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTozzhTollllh_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTozzhTovlvlllh_1TeV_run_01.root"))); Sample_Types_channel1.push_back(1);

	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTozzzTollllvlvl_1TeV_run_01.root"))); Sample_Types_channel1.push_back(0);

	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTozzzTollllvlvl_1_QCKM_5_reweightscan_%s_run_*.root", para_name.c_str()))); Sample_Types_channel1.push_back(2);
	Infiles_names_channel1.push_back(TString(Form("samples_Channel1/mumuTozzzTollllvlvl_1_QCKM_5_benchmark_%s_30.00000_run_*.root", para_name.c_str()))); Sample_Types_channel1.push_back(2);

	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTohtt_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTomumutt_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTomumuwwTomumulvllvl_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTomumuwwz_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTomumuzTomumull_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTomumuzhTomumullh_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTomumuzhTomumuvlvlh_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTomumuzh_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTomumuzzTomumullvlvl_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTovmvmtt_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTovmvmwwh_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTovmvmwwzTovmvmlvllvll_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTovmvmzhh_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTovmvmzzh_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTovmvmzzz_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTovmvmzzTovmvmllll_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTowwhTolvllvlh_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTowwttTolvllvltt_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTowwzTolvllvlll_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTozhhTollhh_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTozhhTovlvlX_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTozttTolltt_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTozzhTollllh_1TeV.root")));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTozzhTovlvlllh_1TeV.root")));

	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTozzzTollllvlvl_1TeV.root")));

	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTozzzTollllvlvl_1_QCKM_5_reweightscan_%s.root", para_name.c_str())));
	Outfiles_names_preselected_channel1.push_back(TString(Form("samples_Channel1/output_mumuTozzzTollllvlvl_1_QCKM_5_benchmark_%s_30.00000.root", para_name.c_str())));

	if(Infiles_names_channel1.size() != Outfiles_names_preselected_channel1.size()) {
		cout<<"Files number does not match! Exit!"<<endl;
		return;
	}

	for(size_t ifile=0; ifile<Infiles_names_channel1.size(); ifile++) {
		TChain* chain_channel1 = new TChain("Delphes");
		chain_channel1->Add(Infiles_names_channel1.at(ifile));
		ExRootTreeReader *treeReader_channel1 = new ExRootTreeReader(chain_channel1);
		Long64_t numberOfEntries_channel1 = treeReader_channel1->GetEntries();
		If_Event_exist_channel1.push_back(AnalyseEvents(treeReader_channel1, Outfiles_names_preselected_channel1.at(ifile), RWnstep));

	}

	std::vector<string> variables_for_cut_channel1;
	variables_for_cut_channel1.push_back("met");
	variables_for_cut_channel1.push_back("Pt4l");
	variables_for_cut_channel1.push_back("Ptll1");
	variables_for_cut_channel1.push_back("Ptll2");
	variables_for_cut_channel1.push_back("deltaR_ll1");
	variables_for_cut_channel1.push_back("deltaR_ll2");
	variables_for_cut_channel1.push_back("M4l");
	string cut_channel1 = hardcut_channel1;
	std::map<string, std::pair<std::pair<std::pair<double, double>, std::pair<double, double>>, double>> cut_information_channel1;
	cut_information_channel1 = OptimizeCuts(Outfiles_names_preselected_channel1, Sample_Types_channel1, variables_for_cut_channel1, cut_channel1, 20, para_name, RWstart, RWstep, RWnstep);

	cout<<"Computing significance: ----------------------------"<<endl;
	std::vector<int> Sig_Types_channel1;
	std::vector<int> Bkg_Types_channel1;
	double significance_aftersel_channel1 = 0;
	cout<<"++++++++++++++++++++++++++++++++++++++++"<<endl;
	Sig_Types_channel1.push_back(2);
	Bkg_Types_channel1.push_back(0);
	Bkg_Types_channel1.push_back(1);
	//Sig_Types_channel1.push_back(0);
	//Bkg_Types_channel1.push_back(1);
	cout<<"Signal types: ";
	for(int k=0; k<Sig_Types_channel1.size(); k++) {
		cout<<Sig_Types_channel1.at(k)<<"; ";
	}
	cout<<endl;
	cout<<"Bkg types: ";
	for(int k=0; k<Bkg_Types_channel1.size(); k++) {
		cout<<Bkg_Types_channel1.at(k)<<"; ";
	}
	cout<<endl;
	//significance_aftersel_channel1 = ComputeSignificance(Outfiles_names_preselected_channel1, Sample_Types_channel1, Sig_Types_channel1, Bkg_Types_channel1, "(1==1)");
	significance_aftersel_channel1 = ComputeSignificance(Outfiles_names_preselected_channel1, Sample_Types_channel1, Sig_Types_channel1, Bkg_Types_channel1, cut_channel1);

	cout<<endl;
	cout<<"=================================================================="<<endl;
	cout<<"Optimized cut-based selection: "<<cut_channel1<<endl;
	cout<<"----------------------------------------"<<endl;
	cout<<"Signal significance S/sqrt(B) after selection is "<<significance_aftersel_channel1<<endl;
	cout<<"=================================================================="<<endl;

	PlotEvents(Outfiles_names_preselected_channel1, If_Event_exist_channel1, Sample_Types_channel1, "bfcut_channel1", "(1==1)", 20, para_name, RWstart, RWstep, RWnstep);
	PlotEvents(Outfiles_names_preselected_channel1, If_Event_exist_channel1, Sample_Types_channel1, "afcut_channel1", cut_channel1, 20, para_name, RWstart, RWstep, RWnstep);

	std::pair<std::vector<TH1D*>, std::vector<int>> histos_and_types_channel1;
	histos_and_types_channel1 = obtain_histograms(Outfiles_names_preselected_channel1, If_Event_exist_channel1, Sample_Types_channel1, "afcut_channel1", cut_channel1, 20, para_name, RWstart, RWstep, RWnstep, "M4l");

	cout<<"--------------------------------------------------------------------------"<<endl;
	cout<<"-------------------------- Channel 1st Finished --------------------------"<<endl;
	cout<<"--------------------------------------------------------------------------"<<endl;







	//--------------------------------------------------------------------------
	//------------------------ Analysis for Channel 2nd ------------------------
	//--------------------------------------------------------------------------

	std::vector<TString> Infiles_names_channel2;
	std::vector<TString> Outfiles_names_preselected_channel2;
	std::vector<bool>    If_Event_exist_channel2;
	std::vector<int>     Sample_Types_channel2; // 0 for SM signal, 1 for SM background, 2 for AQGC signal

	cout<<"--------------------------------------------------------------------------"<<endl;
	cout<<"------------------------ Analysis for Channel 2nd ------------------------"<<endl;
	cout<<"--------------------------------------------------------------------------"<<endl;

	string hardcut_channel2 = "( (Mjj>70 && Mjj<110) && (Mll1>80 && Mll1<100) && (Mll2>80 && Mll2<100) && (deltaM<40) )";

	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTohtt_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTomumutt_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTomumuwwTomumujjjj_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTomumuwwTomumujjlvl_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTomumuwwTomumulvljj_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTomumuwwTomumulvllvl_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTomumuwwz_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTomumuzTomumujj_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTomumuzTomumull_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTomumuzhTomumujjh_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTomumuzhTomumullh_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTomumuzzTomumujjjj_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTomumuzzTomumulljj_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTomumuzzTomumullll_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTowwhTojjjjh_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTowwhTojjlvlh_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTowwhTolvljjh_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTowwhTolvllvlh_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTowwttTolvllvltt_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTowwzTojjlvlll_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTowwzTolvljjll_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTowwzTolvllvljj_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTowwzTolvllvlll_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTozh_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTozhhTojjhh_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTozhhTollhh_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTozttTolltt_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTozzhTojjllh_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTozzhTollllh_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(1);

	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTozzzTolllljj_1_jjcutTeV_run_01.root"))); Sample_Types_channel2.push_back(0);

	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTozzzTolllljj_1_QCKM_5_reweightscan_%s_run_*.root", para_name.c_str()))); Sample_Types_channel2.push_back(2);
	Infiles_names_channel2.push_back(TString(Form("samples_Channel2/mumuTozzzTolllljj_1_QCKM_5_benchmark_%s_30.00000_run_*.root", para_name.c_str()))); Sample_Types_channel2.push_back(2);

	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTohtt_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTomumutt_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTomumuwwTomumujjjj_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTomumuwwTomumujjlvl_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTomumuwwTomumulvljj_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTomumuwwTomumulvllvl_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTomumuwwz_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTomumuzTomumujj_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTomumuzTomumull_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTomumuzhTomumujjh_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTomumuzhTomumullh_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTomumuzzTomumujjjj_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTomumuzzTomumulljj_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTomumuzzTomumullll_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTowwhTojjjjh_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTowwhTojjlvlh_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTowwhTolvljjh_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTowwhTolvllvlh_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTowwttTolvllvltt_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTowwzTojjlvlll_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTowwzTolvljjll_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTowwzTolvllvljj_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTowwzTolvllvlll_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTozh_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTozhhTojjhh_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTozhhTollhh_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTozttTolltt_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTozzhTojjllh_1_jjcutTeV_run_01.root")));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTozzhTollllh_1_jjcutTeV_run_01.root")));

	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTozzzTolllljj_1_jjcutTeV_run_01.root")));

	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTozzzTolllljj_1_QCKM_5_reweightscan_%s.root", para_name.c_str())));
	Outfiles_names_preselected_channel2.push_back(TString(Form("samples_Channel2/output_mumuTozzzTolllljj_1_QCKM_5_benchmark_%s_30.00000.root", para_name.c_str())));

	if(Infiles_names_channel2.size() != Outfiles_names_preselected_channel2.size()) {
		cout<<"Files number does not match! Exit!"<<endl;
		return;
	}

	for(size_t ifile=0; ifile<Infiles_names_channel2.size(); ifile++) {
		TChain* chain_channel2 = new TChain("Delphes");
		chain_channel2->Add(Infiles_names_channel2.at(ifile));
		ExRootTreeReader *treeReader_channel2 = new ExRootTreeReader(chain_channel2);
		Long64_t numberOfEntries_channel2 = treeReader_channel2->GetEntries();
		If_Event_exist_channel2.push_back(AnalyseEvents(treeReader_channel2, Outfiles_names_preselected_channel2.at(ifile), RWnstep));

	}

	std::vector<string> variables_for_cut_channel2;
	variables_for_cut_channel2.push_back("met");
	variables_for_cut_channel2.push_back("Pt4l2j");
	variables_for_cut_channel2.push_back("deltaR_ll1");
	variables_for_cut_channel2.push_back("deltaR_ll2");
	variables_for_cut_channel2.push_back("deltaR_jj");
	variables_for_cut_channel2.push_back("Ptll1");
	variables_for_cut_channel2.push_back("Ptll2");
	variables_for_cut_channel2.push_back("Ptjj");
	variables_for_cut_channel2.push_back("M4l2j");
	string cut_channel2 = hardcut_channel2;
	std::map<string, std::pair<std::pair<std::pair<double, double>, std::pair<double, double>>, double>> cut_information_channel2;
	cut_information_channel2 = OptimizeCuts(Outfiles_names_preselected_channel2, Sample_Types_channel2, variables_for_cut_channel2, cut_channel2, 20, para_name, RWstart, RWstep, RWnstep);

	cout<<"Computing significance: ----------------------------"<<endl;
	std::vector<int> Sig_Types_channel2;
	std::vector<int> Bkg_Types_channel2;
	double significance_aftersel_channel2 = 0;
	cout<<"++++++++++++++++++++++++++++++++++++++++"<<endl;
	Sig_Types_channel2.push_back(2);
	Bkg_Types_channel2.push_back(0);
	Bkg_Types_channel2.push_back(1);
	//Sig_Types_channel2.push_back(0);
	//Bkg_Types_channel2.push_back(1);
	cout<<"Signal types: ";
	for(int k=0; k<Sig_Types_channel2.size(); k++) {
		cout<<Sig_Types_channel2.at(k)<<"; ";
	}
	cout<<endl;
	cout<<"Bkg types: ";
	for(int k=0; k<Bkg_Types_channel2.size(); k++) {
		cout<<Bkg_Types_channel2.at(k)<<"; ";
	}
	cout<<endl;
	//significance_aftersel_channel2 = ComputeSignificance(Outfiles_names_preselected_channel2, Sample_Types_channel2, Sig_Types_channel2, Bkg_Types_channel2, "(1==1)");
	significance_aftersel_channel2 = ComputeSignificance(Outfiles_names_preselected_channel2, Sample_Types_channel2, Sig_Types_channel2, Bkg_Types_channel2, cut_channel2);

	cout<<endl;
	cout<<"=================================================================="<<endl;
	cout<<"Optimized cut-based selection: "<<cut_channel2<<endl;
	cout<<"----------------------------------------"<<endl;
	cout<<"Signal significance S/sqrt(B) after selection is "<<significance_aftersel_channel2<<endl;
	cout<<"=================================================================="<<endl;

	PlotEvents(Outfiles_names_preselected_channel2, If_Event_exist_channel2, Sample_Types_channel2, "bfcut_channel2", "(1==1)", 20, para_name, RWstart, RWstep, RWnstep);
	PlotEvents(Outfiles_names_preselected_channel2, If_Event_exist_channel2, Sample_Types_channel2, "afcut_channel2", cut_channel2, 20, para_name, RWstart, RWstep, RWnstep);

	std::pair<std::vector<TH1D*>, std::vector<int>> histos_and_types_channel2;
	histos_and_types_channel2 = obtain_histograms(Outfiles_names_preselected_channel2, If_Event_exist_channel2, Sample_Types_channel2, "afcut_channel2", cut_channel2, 20, para_name, RWstart, RWstep, RWnstep, "M4l");

	cout<<"--------------------------------------------------------------------------"<<endl;
	cout<<"-------------------------- Channel 1st Finished --------------------------"<<endl;
	cout<<"--------------------------------------------------------------------------"<<endl;





	std::vector<string> channels;
	channels.push_back("channel1");
	channels.push_back("channel2");
	std::vector<std::vector<TH1D*>> histsamples_channels;
	histsamples_channels.push_back(histos_and_types_channel1.first);
	histsamples_channels.push_back(histos_and_types_channel2.first);
	std::vector<std::vector<int>> histsamplestypes_channels;
	histsamplestypes_channels.push_back(histos_and_types_channel1.second);
	histsamplestypes_channels.push_back(histos_and_types_channel2.second);
	std::vector<RooDataSet*> datasamples_channels;
	datasamples_channels.push_back(0);
	datasamples_channels.push_back(0);

	calculate_aQGC_constraints_combined(para_names, channels, histsamples_channels, histsamplestypes_channels, datasamples_channels, 0.95, 1, 1000, 1000);


}
