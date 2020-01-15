// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/AnalysisLoader.hh"

#include <numeric>
#include <functional>

namespace Rivet {


  /// @brief Add a short analysis description here
  class CMS_2019_TTH_TTZBCKG : public Analysis {

  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2019_TTH_TTZBCKG);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      FinalState lepfs;

      // Define lepton cuts
      Cut lep_cuts = (Cuts::abseta < 2.5) && (Cuts::pT > 20*GeV);

      // Initialise and register projections
      //declare(FinalState(lep_cuts, "FS"));

      // Get photons to dress leptons
      IdentifiedFinalState all_photons(lepfs);
      all_photons.acceptIdPair(PID::PHOTON);
      IdentifiedFinalState ph_id(lepfs);
      ph_id.acceptIdPair(PID::PHOTON);
      PromptFinalState photons(ph_id);
      photons.acceptTauDecays(false);
      declareProjection(photons, "photons");

      //Projection to find prompt electrons
      IdentifiedFinalState el_id(lepfs);
      el_id.acceptIdPair(PID::ELECTRON);
      PromptFinalState electrons(el_id);
      electrons.acceptTauDecays(true);
      declareProjection(electrons,"electrons");

      // Switching to new definition with the last boolean considered deprecated
      DressedLeptons dressedelectrons(photons, electrons, 0.1, lep_cuts, true);
      declareProjection(dressedelectrons, "dressedelectrons");

      //Projection to find prompt muons
      IdentifiedFinalState mu_id(lepfs);
      mu_id.acceptIdPair(PID::MUON);
      PromptFinalState muons(mu_id);
      muons.acceptTauDecays(true);
      declareProjection(muons,"muons");

      DressedLeptons dressedmuons(photons, muons, 0.1, lep_cuts, true);
      declareProjection(dressedmuons, "dressedmuons");

      // Find hadronic taus
      //TauFinder tauhadronic(TauFinder::HADRONIC);
      TauFinder tauhadronic(TauFinder::DecayMode::HADRONIC);
      declareProjection(tauhadronic,"TauHadronic");

      // Projection to find neutrinos and produce MET
      FinalState fs_neutrino(Cuts::abseta < 2.5);
      IdentifiedFinalState nu_id(fs_neutrino);
      nu_id.acceptNeutrinos();
      declareProjection(nu_id, "neutrinos");
      PromptFinalState neutrinos(nu_id);
      neutrinos.acceptTauDecays(true);

      // Jet clustering.
      VetoedFinalState vfs(lepfs);
      vfs.addVetoOnThisFinalState(dressedelectrons);
      vfs.addVetoOnThisFinalState(dressedmuons);
      vfs.addVetoOnThisFinalState(neutrinos);
      FastJets jets(vfs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles();
      declareProjection(jets, "Jets");



      // Projection to find B-jets
      declareProjection(HeavyHadrons(Cuts::abseta < 5 && Cuts::pT > 5*GeV), "BCHadrons");
      // Projection to find MET
      declareProjection(MissingMomentum(FinalState(Cuts::abseta < 5)),"MissingET");

      ht_bins ={0,120,180,240,300,360,440,540,680,900,1500};
      ht_j_bins={0,90,140,180,240,300,380,460,540,650,850,1500};
      ht_l_bins={0,20,50,80,110,150,200,300,400,550,800};
      met_bins={0,20,50,80,120,180,300,500,1200};
      lep_bins={0,20,25,33,45,60,80,110,160,500};
      jet_bins={0,20,25,33,45,60,80,110,200};
      bjet_bins={0,20,25,33,45,60,80,110,150,200,300,500};

      //presel, two light-leptons , same sign  ,lepton0 pT, geq1b,geq3jet,regions
      vector<string> s_cutDescs =  {  "Preselections","Nleps","TotCharge","lepPt0>25","mll","1b","2j",
      "0t 1b 2j", "0t 2b 2j", "1t 1b 2j"};

      int Ncuts = s_cutDescs.size();
      h_cutflow_3l[0] = book(h_cutflow_3l[0],(const std::string&)"cf3l", (size_t)Ncuts, 0., (double)Ncuts);
      h_cutflow_3l[1] = book(h_cutflow_3l[1],(const std::string&)"cf3l_raw", (size_t)Ncuts, 0., (double)Ncuts);

      vector<int> lep_bins_v={0,20,25,33,45,60,80,110,160,500};
      vector<string> region_names={"0t 1b 2j", "0t 2b 2j", "1t 1b 2j"};

      for(int i=0; i<(int)region_names.size();i++){
        _h_hist_nJets[i] = book(_h_hist_nJets[i],("nJets_"+to_string(i)), (size_t)7, 2.5, 9.5);
        _h_hist_lep_Pt_0[i] = book(_h_hist_lep_Pt_0[i],("lep_Pt_0_"+to_string(i)),{0,20,25,33,45,60,80,110,160,500});
        _h_hist_lep_Pt_1[i] = book(_h_hist_lep_Pt_1[i],("lep_Pt_1_"+to_string(i)), lep_bins);
        _h_hist_lep_Pt_2[i] = book(_h_hist_lep_Pt_2[i],("lep_Pt_2_"+to_string(i)), lep_bins);
        //jets:
        _h_hist_jet_Pt_1[i] = book(_h_hist_jet_Pt_1[i],("jet_Pt_1_"+to_string(i)), jet_bins);
        _h_hist_jet_Pt_2[i] = book(_h_hist_jet_Pt_2[i],("jet_Pt_2_"+to_string(i)), jet_bins);
        _h_hist_jet_Pt_3[i] = book(_h_hist_jet_Pt_3[i],("jet_Pt_3_"+to_string(i)), jet_bins);
        _h_hist_jet_Pt_4[i] = book(_h_hist_jet_Pt_4[i],("jet_Pt_4_"+to_string(i)), jet_bins);
        //bjets
        _h_hist_Bjet_Pt_0[i] = book(_h_hist_Bjet_Pt_0[i],("Bjet_Pt_0_"+to_string(i)),bjet_bins);
        _h_hist_Bjet_Pt_1[i] = book(_h_hist_Bjet_Pt_1[i],("Bjet_Pt_1_"+to_string(i)), jet_bins);

        _h_hist_min_DRl0j[i] = book(_h_hist_min_DRl0j[i],("min_DRl0j_"+to_string(i)), dr_bins, 0., dr_max);
        _h_hist_min_DRl1j[i] = book(_h_hist_min_DRl1j[i],("min_DRl1j_"+to_string(i)), dr_bins, 0., dr_max);
        _h_hist_min_DRl2j[i] = book(_h_hist_min_DRl2j[i],("min_DRl2j_"+to_string(i)), dr_bins, 0., dr_max);
        _h_hist_maxEta_lll[i] = book(_h_hist_maxEta_lll[i],("maxEta_lll_"+to_string(i)), (size_t)13, 0, 2.6); // maxEta = max( fabs( lep_Eta_0 ), fabs( lep_Eta_1 ) );
        _h_hist_HT_jets[i] = book(_h_hist_HT_jets[i],("HT_jets_"+to_string(i)), ht_j_bins);
        _h_hist_HT_leps[i] = book(_h_hist_HT_leps[i],("HT_leps_"+to_string(i)), ht_l_bins);
        _h_hist_HT[i] = book(_h_hist_HT[i],("HT_"+to_string(i)),ht_bins);// 100, 0., 1000.
        _h_hist_nBtagJets[i] = book(_h_hist_nBtagJets[i],("nBtagJets_"+to_string(i)), (size_t)3, 0.5, 3.5);
        _h_hist_MET[i] = book(_h_hist_MET[i],("MET_"+to_string(i)),met_bins);//100, 0., 1000.
        //
        _h_hist_lep_Eta_0[i] = book(_h_hist_lep_Eta_0[i],("lep_Eta_0_"+to_string(i)), (size_t)13, -2.6, 2.6);
        _h_hist_lep_Eta_1[i] = book(_h_hist_lep_Eta_1[i],("lep_Eta_1_"+to_string(i)),(size_t)13, -2.6, 2.6);
        _h_hist_lep_Eta_2[i] = book(_h_hist_lep_Eta_2[i],("lep_Eta_1_"+to_string(i)),(size_t)13, -2.6, 2.6);
        _h_hist_lep_Phi_0[i] = book(_h_hist_lep_Phi_0[i],("lep_Phi_0_"+to_string(i)),(size_t)16, -3.2, 3.2);
        _h_hist_lep_Phi_1[i] = book(_h_hist_lep_Phi_1[i],("lep_Phi_1_"+to_string(i)),(size_t)16, -3.2, 3.2);
        _h_hist_lep_Phi_2[i] = book(_h_hist_lep_Phi_2[i],("lep_Phi_1_"+to_string(i)),(size_t)16, -3.2, 3.2);

        _h_hist_lep_dPhi[i] = book(_h_hist_lep_dPhi[i],("lep_dPhi_"+to_string(i)),  (size_t)16, 0, 6.4);

      }


    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      /// @todo Do the event by event analysis here
      //const double weight = event.weight();
      const MissingMomentum& met = applyProjection<MissingMomentum>(event, "MissingET");
      const double event_met	 = met.vectorEt().mod();

      //Jets alljets;
      Jets alljets = applyProjection<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);
      alljets    =   sortByPt(alljets);

      Particles tauVec,lepVec;
      bool orl_lf=false;
      const Particles& electrons = applyProjection<DressedLeptons>(event,"dressedelectrons").particlesByPt(Cuts::pT > 20*GeV && Cuts::abseta < 2.5);
      for(const Particle & el : electrons ){
        orl_lf=false;
        for(const Jet& jet : alljets){
          if( fabs(deltaR(jet,el))<0.4 ) orl_lf=true;
        }
        if(!orl_lf) lepVec.push_back(el);
      }
      const Particles& muons = applyProjection<DressedLeptons>(event,"dressedmuons").particlesByPt(Cuts::pT > 20*GeV && Cuts::abseta < 2.5);
      for(const Particle &mu : muons ){
        orl_lf=false;
        for(const Jet& jet : alljets){
          if( fabs(deltaR(jet,mu))<0.4 ) orl_lf=true;
        }

        if(!orl_lf) lepVec.push_back(mu);

      }

      const TauFinder &tauhad = applyProjection<TauFinder>(event,"TauHadronic");
      for(const Particle &tau : tauhad.taus()){
        if(tau.pT()>25*GeV && fabs(tau.eta()) < 2.5 ){
          int nProng = 0;
          for(Particle p : tau.children()){
            if (p.charge3()!=0) nProng++;
          }
          if(nProng ==1 || nProng ==3){
            tauVec.push_back(tau);
          }
        }
      }

      tauVec= sortByPt(tauVec);
      lepVec= sortByPt(lepVec);

      int nLep = lepVec.size();

      double ht_jet = 0.0;
      double ht_lep = 0.0;
      for(const Jet& j : alljets) { ht_jet += j.pT(); }
      for(const Particle & part: lepVec) { ht_lep += part.pT(); }

      double ht = ht_jet+ht_lep;

      /* Identify b-jets */
      const Particles bhadrons = sortByPt(applyProjection<HeavyHadrons>(event, "BCHadrons").bHadrons());

      Jets bjets, ljets;
      for(const Jet& jet : alljets) {
        if(jet.bTagged())   bjets.push_back(jet);
        else ljets.push_back(jet);

      }

      bjets      =   sortByPt(bjets);
      ljets      =   sortByPt(ljets);

      int Njets=0;       int Nbjets=0;      int Nhtaus=0;
      Njets = alljets.size();
      Nbjets= bjets.size();
      Nhtaus= tauVec.size();

      float max_eta=999;
      if(nLep>2){
        if ( max( fabs( lepVec.at(0).eta() ) , fabs( lepVec.at(1).eta() ) ) <  fabs( lepVec.at(2).eta() )){
          max_eta = fabs( lepVec.at(2).eta() );
        }
        else{
          max_eta = max( fabs( lepVec.at(0).eta() ) , fabs( lepVec.at(1).eta() ) );
        }
      }

      float min_lj_deltaR=100; float min_l0j_deltaR=100;   float min_l1j_deltaR=100; float min_l2j_deltaR=100;
      for(const Jet& jet : alljets){
        if(nLep>0)
        if(min_l0j_deltaR > fabs(deltaR(jet,lepVec.at(0)) ) ) {min_l0j_deltaR = fabs(deltaR(jet,lepVec.at(0))); }

        if(nLep>1)
        if(min_l1j_deltaR > fabs(deltaR(jet,lepVec.at(1)))) {min_l1j_deltaR = fabs(deltaR(jet,lepVec.at(1))); }

        if(nLep>2)
        if(min_l2j_deltaR > fabs(deltaR(jet,lepVec.at(2)))) {min_l2j_deltaR = fabs(deltaR(jet,lepVec.at(2))); }

        for(const Particle & part : lepVec){
          if(min_lj_deltaR > fabs(deltaR(jet,part))) {min_lj_deltaR = fabs(deltaR(jet,part)); }
        }
      }


      bool sel_array[10];
      vector<string> region_names={"0t 1b 2j", "0t 2b 2j","1t 1b 2j"};

      int cf_counter =0;

      //presel
      h_cutflow_3l[0]->fill(cf_counter);  h_cutflow_3l[1]->fill(cf_counter,1);
      cf_counter++;

      //3 light-leptons
      if(nLep!=3) vetoEvent;
      h_cutflow_3l[0]->fill(cf_counter);  h_cutflow_3l[1]->fill(cf_counter,1);
      cf_counter++;

      //total charge
      if(fabs(lepVec.at(0).charge()+lepVec.at(1).charge()+lepVec.at(2).charge() ) != 1) vetoEvent;
      h_cutflow_3l[0]->fill(cf_counter);  h_cutflow_3l[1]->fill(cf_counter,1);
      cf_counter++;

      //leading lepton pT
      if( lepVec.at(0).pT()/GeV <25) vetoEvent;
      h_cutflow_3l[0]->fill(cf_counter);  h_cutflow_3l[1]->fill(cf_counter,1);
      cf_counter++;

      //no same flavour pair with dilepton inv. mass below 12 GeV or in Z mass window
      float l0l1_mll = (lepVec.at(0).mom() + lepVec.at(1).mom()).mass();
      float l1l2_mll = (lepVec.at(1).mom() + lepVec.at(2).mom()).mass();
      float l0l2_mll = (lepVec.at(0).mom() + lepVec.at(2).mom()).mass();
      if( lepVec.at(0).abspid() == lepVec.at(1).abspid() ){ if( fabs(l0l1_mll-91) < 10 || fabs(l0l1_mll) < 12){ vetoEvent; } }
      if( lepVec.at(1).abspid() == lepVec.at(2).abspid() ){ if( fabs(l1l2_mll-91) < 10 || fabs(l1l2_mll) < 12){ vetoEvent; } }
      if( lepVec.at(0).abspid() == lepVec.at(2).abspid() ){ if( fabs(l0l2_mll-91) < 10 || fabs(l0l2_mll) < 12){ vetoEvent; } }
      h_cutflow_3l[0]->fill(cf_counter);  h_cutflow_3l[1]->fill(cf_counter,1);
      cf_counter++;

      //geq1b
      if( Nbjets < 1 ) vetoEvent;
      h_cutflow_3l[0]->fill(cf_counter);  h_cutflow_3l[1]->fill(cf_counter,1);
      cf_counter++;

      //geq2jet
      if( Njets < 2) vetoEvent;
      h_cutflow_3l[0]->fill(cf_counter);  h_cutflow_3l[1]->fill(cf_counter,1);
      cf_counter++;

      sel_array[0]=(Nhtaus == 0 && Nbjets == 1 && Njets >= 2 );  // Region 1
      sel_array[1]=(Nhtaus == 0 && Nbjets >= 2 && Njets >= 2 );  // Region 2
      sel_array[2]=(Nhtaus == 1 && Nbjets == 1 && Njets >= 2 );  // Region 3


      for(int i=0; i<(int)region_names.size();i++){
        if(sel_array[i]){

          h_cutflow_3l[0]->fill(cf_counter);  h_cutflow_3l[1]->fill(cf_counter,1);
          cf_counter++;

          _h_hist_nJets[i]->fill(Njets);
          _h_hist_nBtagJets[i]->fill(Nbjets);
          _h_hist_lep_Pt_0[i]->fill(lepVec.at(0).pT()/GeV);
          _h_hist_lep_Pt_1[i]->fill(lepVec.at(1).pT()/GeV);
          _h_hist_lep_Pt_2[i]->fill(lepVec.at(2).pT()/GeV);
          _h_hist_HT[i]->fill(ht);
          _h_hist_HT_jets[i]->fill(ht_jet);
          _h_hist_HT_leps[i]->fill(ht_lep);
          _h_hist_jet_Pt_1[i]->fill(alljets.at(0).pT()/GeV);
          _h_hist_jet_Pt_2[i]->fill(alljets.at(1).pT()/GeV);
          if( Njets>2) _h_hist_jet_Pt_3[i]->fill(alljets.at(2).pT()/GeV);
          if( Njets>3) _h_hist_jet_Pt_4[i]->fill(alljets.at(3).pT()/GeV);
          _h_hist_Bjet_Pt_0[i]->fill(bjets.at(0).pT()/GeV);
          if( Nbjets>1)_h_hist_Bjet_Pt_1[i]->fill(bjets.at(1).pT()/GeV);
          _h_hist_min_DRl0j[i]->fill(min_l0j_deltaR);
          _h_hist_min_DRl1j[i]->fill(min_l1j_deltaR);
          _h_hist_min_DRl2j[i]->fill(min_l2j_deltaR);
          _h_hist_MET[i]->fill(event_met/GeV);
          _h_hist_maxEta_lll[i]->fill(max_eta);
          //https://rivet.hepforge.org/code/2.4.0/a00431.html#a20b7db1dfdcaeedded0f6a81c24ba6ab
          // phi: [0,2pi]
          // in atlas phi: [-pi,pi]
          //apply 1pi shift
          _h_hist_lep_dPhi[i]->fill( lepVec.at(0).phi()-lepVec.at(1).phi());
          _h_hist_lep_Phi_0[i]->fill(lepVec.at(0).phi()-pi);
          _h_hist_lep_Phi_1[i]->fill(lepVec.at(1).phi()-pi);
          _h_hist_lep_Phi_2[i]->fill(lepVec.at(2).phi()-pi);
          _h_hist_lep_Eta_0[i]->fill(lepVec.at(0).eta());
          _h_hist_lep_Eta_1[i]->fill(lepVec.at(1).eta());
          _h_hist_lep_Eta_2[i]->fill(lepVec.at(2).eta());

        }//region selection
      } //region loop
    }//analyze loop



    /// Normalise histograms etc., after the run
    void finalize() {

      //normalize(_h_YYYY); // normalize to unity
      //scale(_h_ZZZZ, crossSection()/picobarn/sumOfWeights()); // norm to cross section

    }


    private:

      // @name Histogram data members
      //@{
      //
      Histo1DPtr h_cutflow_3l[2];
      //Histo1DPtr _h_hist_DRll01[10];
      Histo1DPtr _h_hist_jet_Pt_1[10];
      Histo1DPtr _h_hist_jet_Pt_2[10];
      Histo1DPtr _h_hist_jet_Pt_3[10];
      Histo1DPtr _h_hist_jet_Pt_4[10];
      Histo1DPtr _h_hist_Bjet_Pt_0[10];
      Histo1DPtr _h_hist_Bjet_Pt_1[10];
      Histo1DPtr _h_hist_lep_Pt_0[10];
      Histo1DPtr _h_hist_lep_Pt_1[10];
      Histo1DPtr _h_hist_lep_Pt_2[10];
      Histo1DPtr _h_hist_min_DRl0j[10];
      Histo1DPtr _h_hist_min_DRl1j[10];
      Histo1DPtr _h_hist_min_DRl2j[10];
      Histo1DPtr _h_hist_maxEta_lll[10];
      Histo1DPtr _h_hist_HT_jets[10];
      Histo1DPtr _h_hist_HT_leps[10];
      Histo1DPtr _h_hist_HT[10];
      Histo1DPtr _h_hist_nJets[10];
      Histo1DPtr _h_hist_nBtagJets[10];
      Histo1DPtr _h_hist_MET[10];
      Histo1DPtr _h_hist_lep_Eta_0[10];
      Histo1DPtr _h_hist_lep_Eta_1[10];
      Histo1DPtr _h_hist_lep_Eta_2[10];
      Histo1DPtr _h_hist_lep_Phi_0[10];
      Histo1DPtr _h_hist_lep_Phi_1[10];
      Histo1DPtr _h_hist_lep_Phi_2[10];
      Histo1DPtr _h_hist_lep_dPhi[10];

      vector<double> lep_bins;
      vector<double> ht_bins;
      vector<double> ht_j_bins;
      vector<double> ht_l_bins;
      vector<double> met_bins;
      vector<double> bjet_bins;
      vector<double> jet_bins;
      float dr_max=4.8; int dr_bins=12;

};


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2019_TTH_TTZBCKG);
}
