////////////////////////////////////////////////////////////////////////
// Class:       WireDump
// Module Type: analyzer
// File:        WireDump_module.cc
//
// Generated at Thu Nov  5 08:25:45 2015 by Jonathan Asaadi using artmod
// from cetpkgsupport v1_08_06.
////////////////////////////////////////////////////////////////////////
//
// ##########################
// ### Framework includes ###
// ##########################
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
//#include "cetlib/maybe_ref.h"

// ########################
// ### LArSoft includes ###
// ########################
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/RecoBaseArt/TrackUtils.h" // lar::util::TrackPitchInView()
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"

//#include "RawData/ExternalTrigger.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larevt/Filters/ChannelFilter.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "LArIATDataProducts/WCTrack.h"
#include "LArIATDataProducts/TOF.h"
#include "LArIATDataProducts/AGCounter.h"
#include "RawDataUtilities/TriggerDigitUtility.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCStep.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"


// #####################
// ### ROOT includes ###
// #####################
#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"
#include "TTimeStamp.h"

const int kMaxTicks      = 4500;  //maximum number of tracks
const int kMaxWires	 = 500;
const int kMaxPrimaries  = 20000;  //maximum number of primary particles
const int kMaxPrimaryPart = 10;	    //maximum number of true primary particles
const int kMaxTruePrimaryPts = 5000; //maximum number of points in the true primary trajectory

   // === Storing Geant4 MC Truth Information ===
   int no_primaries;                             //<---Number of primary Geant4 particles in the event
   int geant_list_size;                          //<---Number of Geant4 particles tracked
   int pdg[kMaxPrimaries];                       //<---PDG Code number of this particle
   Float_t Eng[kMaxPrimaries];                   //<---Energy of the particle
   Float_t Px[kMaxPrimaries];                    //<---Px momentum of the particle
   Float_t Py[kMaxPrimaries];                    //<---Py momentum of the particle
   Float_t Pz[kMaxPrimaries];                    //<---Pz momentum of the particle
   Float_t StartPointx[kMaxPrimaries];           //<---X position that this Geant4 particle started at
   Float_t StartPointy[kMaxPrimaries];           //<---Y position that this Geant4 particle started at
   Float_t StartPointz[kMaxPrimaries];           //<---Z position that this Geant4 particle started at
   Float_t EndPointx[kMaxPrimaries];             //<---X position that this Geant4 particle ended at
   Float_t EndPointy[kMaxPrimaries];             //<---Y position that this Geant4 particle ended at
   Float_t EndPointz[kMaxPrimaries];             //<---Z position that this Geant4 particle ended at
   Float_t Startdcosx[kMaxPrimaries];            //<---X direction cosine that Geant4 particle started at
   Float_t Startdcosy[kMaxPrimaries];            //<---Y direction cosine that Geant4 particle started at
   Float_t Startdcosz[kMaxPrimaries];            //<---Z direction cosine that Geant4 particle started at
   int NumberDaughters[kMaxPrimaries];           //<---Number of Daughters this particle has
   int TrackId[kMaxPrimaries];                   //<---Geant4 TrackID number
   int Mother[kMaxPrimaries];                    //<---TrackID of the mother of this particle
   int process_primary[kMaxPrimaries];           //<---Is this particle primary (primary = 1, non-primary = 1)
   //std::vector<std::string> G4Process;         //<---The process which created this particle
   //std::vector<std::string> G4FinalProcess;    //<---The last process which this particle went under
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //+++++ adding more information for geant4 particle Track_length ++++++++++
   Float_t xi[kMaxPrimaries];                    //+++++<--- StartPointx - EndPointx
   Float_t yi[kMaxPrimaries];                    //+++++<--- StartPointy - EndPointy
   Float_t zi[kMaxPrimaries];                    //+++++<--- StartPointz - EndPointz
   Float_t x_i[kMaxPrimaries];                    //+++++<--- using to comput track_length "if the track hit the boundry"
   Float_t y_i[kMaxPrimaries];                    //+++++<--- using to comput track_length "if the track hit the boundry"
   Float_t z_i[kMaxPrimaries];                    //+++++<--- using to comput track_length "if the track hit the boundry"
   Float_t Track_length[kMaxPrimaries];          //+++++<---Track_length of the Geant4 particle

   // === Storing additionnal Geant4 MC Truth Information for the primary track only ===
   int NTrTrajPts[kMaxPrimaryPart];							 //<--Nb. of true points in the true primary trajectories
   double MidPosX[kMaxPrimaryPart][kMaxTruePrimaryPts];//<--X position of a point in the true primary trajectory
   double MidPosY[kMaxPrimaryPart][kMaxTruePrimaryPts];//<--Y position of a point in the true primary trajectory
   double MidPosZ[kMaxPrimaryPart][kMaxTruePrimaryPts];//<--Z position of a point in the true primary trajectory
   double MidPx[kMaxPrimaryPart][kMaxTruePrimaryPts];  //<- Px momentum of a point in the true primary trajectory
   double MidPy[kMaxPrimaryPart][kMaxTruePrimaryPts];  //<--Py momentum of a point in the true primary trajectory
   double MidPz[kMaxPrimaryPart][kMaxTruePrimaryPts];  //<--Pz momentum of a point in the true primary trajectory

   std::vector<int>    InteractionPoint;         //<---Geant 4 Primary Trj Point Corresponding to the Interaction
   std::vector<int>    InteractionPointType;     //<---Geant 4 Primary Interaction Type



   // === mctruth information ===
   Int_t     mcevts_truth;    //<---number of neutrino Int_teractions in the spill
   Int_t     nuPDG_truth;     //<---neutrino PDG code
   Int_t     ccnc_truth;      //<---0=CC 1=NC
   Int_t     mode_truth;      //<---0=QE/El, 1=RES, 2=DIS, 3=Coherent production
   Float_t   enu_truth;       //<---true neutrino energy
   Float_t   Q2_truth;        //<---Momentum transfer squared
   Float_t   W_truth;         //<---hadronic invariant mass
   Float_t   X_truth;
   Float_t   Y_truth;
   Int_t     hitnuc_truth;    //<---hit nucleon
   Int_t     target_truth;    //<---hit nucleus
   Float_t   nuvtxx_truth;    //<---neutrino vertex x
   Float_t   nuvtxy_truth;    //<---neutrino vertex y
   Float_t   nuvtxz_truth;    //<---neutrino vertex z
   Float_t   nu_dcosx_truth;  //<---neutrino dcos x
   Float_t   nu_dcosy_truth;  //<---neutrino dcos y
   Float_t   nu_dcosz_truth;  //<---neutrino dcos z
   Float_t   lep_mom_truth;   //<---lepton momentum
   Float_t   lep_dcosx_truth; //<---lepton dcos x
   Float_t   lep_dcosy_truth; //<---lepton dcos y
   Float_t   lep_dcosz_truth; //<---lepton dcos z
   Float_t   t0_truth;        //<--- t0
 //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 //++++++++++++++++++++adding the class of the nu interaction ++++++++++++++
   Int_t  class_truth;      //<++++ 0=numu inclusive , 1=numu_CC_QE_1track , 2= nummu_CC_QE_1proton-track, 3=numu_CC_2tracks_pi+/- ,4=numu_CC_3tracks_2pi+/- , 5=numu_CC_2tracks_pi0 ,6=numu_CC_3tracks_pi+/-&pi0 ,7=numu_CC_other ,8=nue_CC ,9=numu_CC_wrong_song ,10=NC_inclusive ,11=NC_1pi+/- ,12=NC_1pi0 ,13=NC_other


class WireDump;

class WireDump : public art::EDAnalyzer {
public:
  explicit WireDump(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WireDump(WireDump const &) = delete;
  WireDump(WireDump &&) = delete;
  WireDump & operator = (WireDump const &) = delete;
  WireDump & operator = (WireDump &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;


private:

   // === Function used to reset all the variables  ===
   void ResetVars();
  //
   // === Storing information into TTree ====
   TTree* fTree;
   
   //=== Storing Run Information ===
   int run;				//<---Run Number
   int subrun;				//<---SubRun Number
   int event;				//<---Event Number
   
   // === Storing the ADC/Time Tick info for each wire ===
   int nSamples;			//<---Number of Time Ticks
   
   int nWires;				//<---Number of Wires
   
   float WireADC[kMaxWires][4500]; //<---Wire Number and ADC
   
   
   // ### Input Module Lables ###
   std::string	fDetSimModuleLabel;///< name of module that produced the digits
   std::string fGenieGenModuleLabel;  ////Module labels to get data products

};


WireDump::WireDump(fhicl::ParameterSet const & pset) : EDAnalyzer(pset)  // , // More initializers here.
{ 
this->reconfigure(pset);

}


void WireDump::reconfigure(fhicl::ParameterSet const & pset)
{
  fDetSimModuleLabel 	= pset.get< std::string >("DetSimModuleLabel");
  fGenieGenModuleLabel  = pset.get< std::string >("GenieGenModuleLabel");
}

// #################################################
// ################### Event Loop ##################
// #################################################
void WireDump::analyze(art::Event const & evt)
{
   
   // #############################################
   // ### Reset variables before we get started ###
   // #############################################
   ResetVars();
   
   
    // === BackTracker service ===
    art::ServiceHandle<cheat::BackTracker> bt;
    const sim::ParticleList& plist = bt->ParticleList();
    
   // === Run Number ===
   run = evt.run();
   // === Sub-Run Number ===
   subrun = evt.subRun();
   // === Event Number ===
   event = evt.id().event();
  
   // ####################################
   // ### Grabbing the RawDigit Handle ###
   // ####################################
   art::Handle< std::vector<raw::RawDigit> > rdHandle;
   evt.getByLabel(fDetSimModuleLabel,rdHandle);
   
   // ### Filling a RawDigit Vector ###
   art::PtrVector<raw::RawDigit> rdvec;
   
   // ### Looping over the RawDigit Handle ###
   for(unsigned int i = 0; i < rdHandle->size(); ++i)
      {
      art::Ptr<raw::RawDigit> r(rdHandle,i);
      rdvec.push_back(r);
      }//<---End i loop
     
     
    // ### Filling the number of samples for all wires ###
    nSamples = rdvec[0]->Samples(); 
    nWires = rdvec.size();
    // #######################################
    // ### Looping over all the raw digits ###
    // #######################################
    for(unsigned int rd = 0; rd < rdvec.size(); ++rd)
       {
       //std::cout<<"Wire = "<<rd<<std::endl;
       // ### Loop over all the time ticks for this wire ###
       for(unsigned int ntick = 0; ntick < rdvec[rd]->Samples(); ntick++)
          {
	  //if(rdvec[rd]->ADC(ntick) > 2){std::cout<<rdvec[rd]->ADC(ntick)<<std::endl;}
	  WireADC[rd][ntick] = rdvec[rd]->ADC(ntick);
	  //std::cout<<"ADC = "<<WireADC[rd][ntick]<<std::endl;
	  
	  }//<---End ntick loop
       
       
       
       }//<---End rd loop
    
    // ######################################
    // ### Making a vector of MCParticles ###
    // ######################################
    std::vector<const simb::MCParticle* > geant_part;
    
    // ### Looping over all the Geant4 particles from the BackTracker ###
    for(size_t p = 0; p < plist.size(); ++p)
    {
        // ### Filling the vector with MC Particles ###
        geant_part.push_back(plist.Particle(p));
    }
    
    // ### Setting a string for primary ###
    std::string pri("primary");
    
    // ### Setting a string for PionMinusInelastic ###
    std::string PionMinusInelastic("pi-Inelastic");
    
    // ### Setting a string for NeutronInelastic ###
    std::string NeutronInelastic("neutronInelastic");
    
    // ### Setting a string for hadElastic ###
    std::string hadElastic("hadElastic");
    
    // ### Setting a string for nCapture ###
    std::string nCapture("nCapture");
    
    // This may not be called hBertiniCaptureAtRest ?
    // ### Setting a string for CHIPSNuclearCaptureAtRest ###
    std::string CHIPSNuclearCaptureAtRest("CHIPSNuclearCaptureAtRest");
    
    // ### Setting a string for Decay ###
    std::string Decay("Decay");
    
    // ### Setting a string for KaonZeroLInelastic ###
    std::string KaonZeroLInelastic("KaonZeroLInelastic");
    
    // ### Setting a string for CoulombScat ###
    std::string CoulombScat("CoulombScat");
    
    // ### Setting a string for muMinusCaptureAtRest ###
    std::string muMinusCaptureAtRest("muMinusCaptureAtRest");
    
    // ### Setting a string for ProtonInelastic ###
    std::string ProtonInelastic("protonInelastic");
    
    // ### Setting a string for Kaon Inelastic ###
    std::string KaonPlusInelastic("kaon+Inelastic");
    
    // ### Setting a string for BertiniCaptureAtRest
    std::string hBertiniCaptureAtRest("hBertiniCaptureAtRest");

    
    int primary=0;
    int geant_particle=0;
    
    // ############################################################
    // ### Determine the number of primary particles from geant ###
    // ############################################################
    for( unsigned int i = 0; i < geant_part.size(); ++i )
    {
        geant_particle++;
        // ### Counting the number of primary particles ###
        if(geant_part[i]->Process()==pri)
        { primary++;}
    }//<---End i loop
    
    
    // ### Saving the number of primary particles ###
    no_primaries=primary;
    // ### Saving the number of Geant4 particles ###
    geant_list_size=geant_particle;
    
    int iPrim = 0;
    // ### Looping over all the Geant4 particles ###
    for( unsigned int i = 0; i < geant_part.size(); ++i )
    {
        // ### If this particle is primary, set = 1 ###
        if(geant_part[i]->Process()==pri)
        {process_primary[i]=1;}
        // ### If this particle is not-primary, set = 0 ###
        else
        {process_primary[i]=0;}
        
        // ### Saving the particles mother TrackID ###
        Mother[i]=geant_part[i]->Mother();
        
        // ### Saving the particles TrackID ###
        TrackId[i]=geant_part[i]->TrackId();
        
        // ### Saving the PDG Code ###
        pdg[i]=geant_part[i]->PdgCode();
        
        // ### Saving the particles Energy ###
        Eng[i]=geant_part[i]->E();
        
        // ### Saving the Px, Py, Pz info ###
        Px[i]=geant_part[i]->Px();
        Py[i]=geant_part[i]->Py();
        Pz[i]=geant_part[i]->Pz();
        
        // ### Saving the Start and End Point for this particle ###
        StartPointx[i]=geant_part[i]->Vx();
        StartPointy[i]=geant_part[i]->Vy();
        StartPointz[i]=geant_part[i]->Vz();
        EndPointx[i]=geant_part[i]->EndPosition()[0];
        EndPointy[i]=geant_part[i]->EndPosition()[1];
        EndPointz[i]=geant_part[i]->EndPosition()[2];
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++ Computing the Track_length +++++++++++++++++++
        xi[i] = StartPointx[i] - EndPointx[i];
        yi[i] = StartPointy[i] - EndPointy[i];
        zi[i] = StartPointz[i] - EndPointz[i];
        x_i[i] = 500-StartPointx[i];
        z_i[i] = 1000-StartPointz[i];
        y_i[i] = 200-abs(StartPointy[i]);
        
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //+++++++++++++++++++++++++++computing the track_length using the functions above++++++++++++++++++++++++++++++++++++++++++++
        if (yi[i] != 0.0 and xi[i] != 0.0 and zi[i] != 0.0) {
            //when endx point is out of the boundry
            if (abs(EndPointy[i]) < 200 and EndPointz[i] < 1000 and EndPointx[i] > 500) {
                y_i[i] = EndPointy[i]+(yi[i]/xi[i])*(x_i[i]-EndPointx[i]);
                z_i[i] = EndPointz[i]+(zi[i]/xi[i])*(x_i[i]-EndPointx[i]);
            }
            //when endz point is out of the boundry
            else if (abs(EndPointy[i]) < 200 and EndPointz[i] > 1000 and EndPointx[i] < 500 ) {
                y_i[i] = EndPointy[i] +(yi[i]/zi[i])*(z_i[i]-EndPointz[i]);
                x_i[i] = EndPointx[i]+(xi[i]/zi[i])*(z_i[i]-EndPointz[i]);
            }
            //when endy point is out of the boundry
            else if (abs(EndPointy[i]) > 200 and EndPointz[i] < 1000 and EndPointx[i] < 500) {
                x_i[i] = EndPointx[i]+(xi[i]/yi[i])*(y_i[i]-EndPointy[i]);
                z_i[i] = EndPointz[i]+(zi[i]/yi[i])*(y_i[i]-EndPointy[i]);
            }
            //when end points is in the boundry
            else if (abs(EndPointy[i]) < 200 and EndPointz[i] < 1000 and EndPointx[i] < 500) {
                y_i[i] = EndPointy[i];
                z_i[i] = EndPointz[i];
                x_i[i] = EndPointx[i];
            }
            //z_i[i] = zPosition(StartPointx[i],StartPointy[i],StartPointz[i],EndPointx[i],EndPointy[i],EndPointz[i],xi[i],yi[i],zi[i]);
            //x_i[i] = xPosition(StartPointx[i],StartPointy[i],StartPointz[i],EndPointx[i],EndPointy[i],EndPointz[i],xi[i],yi[i],zi[i]);
            Track_length[i] = sqrt(pow((StartPointx[i]-x_i[i]),2)+pow((StartPointy[i]-y_i[i]),2)+pow((StartPointz[i]-z_i[i]),2));
        }
        //Wheb
        else if(yi[i] == 0.0){
            Track_length[i] = sqrt(pow((StartPointx[i]-EndPointx[i]),2)+pow((StartPointz[i]-EndPointz[i]),2));
        }
        else if(xi[i] == 0.0){
            Track_length[i] = sqrt(pow((StartPointy[i]-EndPointy[i]),2)+pow((StartPointz[i]-EndPointz[i]),2));
        }
        else if(zi[i] == 0.0){
            Track_length[i] = sqrt(pow((StartPointy[i]-EndPointy[i]),2)+pow((StartPointx[i]-EndPointx[i]),2));
        }
        else{
            Track_length[i] = -9999999.999999;
        }
        
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        
        
        // ### Saving the processes for this particle ###
        //std::cout<<"finding proc"<<std::endl;
        //G4Process.push_back( geant_part[i]->Process() );
        //G4FinalProcess.push_back( geant_part[i]->EndProcess() );
        
        // ### Saving the Start direction cosines for this particle ###
        Startdcosx[i] = geant_part[i]->Momentum(0).Px() / geant_part[i]->Momentum(0).P();
        Startdcosy[i] = geant_part[i]->Momentum(0).Py() / geant_part[i]->Momentum(0).P();
        Startdcosz[i] = geant_part[i]->Momentum(0).Pz() / geant_part[i]->Momentum(0).P();
        
        // ### Saving the number of Daughters for this particle ###
        NumberDaughters[i]=geant_part[i]->NumberDaughters();
        
        ///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        ///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        // ### Save intermediary information for the primary track
        if(geant_part[i]->Process()==pri){
            NTrTrajPts[i]=geant_part[i]->NumberTrajectoryPoints();
            simb::MCTrajectory truetraj = geant_part[i]->Trajectory();
            
            int iPrimPt = 0;
            for(auto itTraj = truetraj.begin(); itTraj != truetraj.end(); ++itTraj)
            {
                MidPosX[iPrim][iPrimPt] = truetraj.X(iPrimPt);
                MidPosY[iPrim][iPrimPt] = truetraj.Y(iPrimPt);
                MidPosZ[iPrim][iPrimPt] = truetraj.Z(iPrimPt);
                MidPx[iPrim][iPrimPt] = truetraj.Px(iPrimPt);
                MidPy[iPrim][iPrimPt] = truetraj.Py(iPrimPt);
                MidPz[iPrim][iPrimPt] = truetraj.Pz(iPrimPt);
                iPrimPt++;
            }//<--End loop on true trajectory points
            
            
            // Yet an other scheme for interaction type
            
            // ### Recording the process as a integer ###
            // 0 = NoInteractionNodaughters, thought going
            // 1 = PionMinusInelastic
            // 2 = NeutronInelastic
            // 3 = hadElastic
            // 4 = nCapture
            // 5 = CHIPSNuclearCaptureAtRest
            // 6 = Decay
            // 7 = KaonZeroLInelastic
            // 8 = CoulombScat
            // 9 = muMinusCaptureAtRest
            //10 = ProtonInelastic
            //11 = Kaon+Inelastic
            //12 = Kaon-Inelastic
            //13 = protonInelastic
            //14 = Pi+Inelastic
            
            
            auto thisTracjectoryProcessMap =  truetraj.TrajectoryProcesses();
            
            // Ok, if the size of the map is 0, all the action might happen at the end of the track
            // So we check the daugthers:
            //    - Case 1. There are daugthers:
            //               * The interesting point is the last one
            //               * The interaction type is the one that created the first daugther (this might be improved)
            //    - Case 2. There are NO daugthers:
            //              * We assign the interaction type to be 0: nothing happens, thought going particle
            //              * The interesting point is the last one (might not be in the TPC)
            if (!thisTracjectoryProcessMap.size())
            {
                int interestingPoint = (int) (NTrTrajPts[i] - 1);
                InteractionPoint.push_back(interestingPoint);
                
                if (NumberDaughters[i])
                {
                    auto thePrimaryDaughterID = geant_part[i]-> Daughter(0);
                    for( unsigned int iD = 0; iD < geant_part.size(); ++iD )
                    {
                        if (geant_part[iD]->TrackId() == thePrimaryDaughterID)
                        {
                            if(geant_part[iD]->Process() == PionMinusInelastic)
                            {InteractionPointType.push_back(1);}
                            
                            if(geant_part[iD]->Process() == NeutronInelastic)
                            {InteractionPointType.push_back(2);}
                            
                            if(geant_part[iD]->Process() == hadElastic)
                            {InteractionPointType.push_back(3);}
                            
                            if(geant_part[iD]->Process() == nCapture)
                            {InteractionPointType.push_back(4);}
                            
                            if(geant_part[iD]->Process() == CHIPSNuclearCaptureAtRest)
                            {InteractionPointType.push_back(5);}
                            
                            if(geant_part[iD]->Process() == Decay)
                            {InteractionPointType.push_back(6);}
                            
                            if(geant_part[iD]->Process() == KaonZeroLInelastic)
                            {InteractionPointType.push_back(7);}
                            
                            if(geant_part[iD]->Process() == CoulombScat)
                            {InteractionPointType.push_back(8);}
                            
                            if(geant_part[iD]->Process() == muMinusCaptureAtRest)
                            {InteractionPointType.push_back(9);}
                            
                            if(geant_part[iD]->Process() == ProtonInelastic)
                            {InteractionPointType.push_back(10);}
                            
                            if(geant_part[iD]->Process() == KaonPlusInelastic)
                            {InteractionPointType.push_back(11);}
                            
                            if(geant_part[iD]->Process() == hBertiniCaptureAtRest)
                            {InteractionPointType.push_back(12);}
                        }
                    }
                }else
                {
                    InteractionPointType.push_back(0);
                }
            }else
            {
                // The map is not zero: somthing interesting might happen in the middle of the track!!
                for(auto const& couple: thisTracjectoryProcessMap) 
                {
                    int interestingPoint = (int) couple.first;
                    InteractionPoint.push_back(interestingPoint);         	   
                    if ((truetraj.KeyToProcess(couple.second)).find("hadElastic")!= std::string::npos) InteractionPointType.push_back(3);           
                    if ((truetraj.KeyToProcess(couple.second)).find("pi-Inelastic")    != std::string::npos) InteractionPointType.push_back(1);           
                    if ((truetraj.KeyToProcess(couple.second)).find("pi+Inelastic")    != std::string::npos) InteractionPointType.push_back(14);           
                    if ((truetraj.KeyToProcess(couple.second)).find("kaon-Inelastic")  != std::string::npos) InteractionPointType.push_back(12);           
                    if ((truetraj.KeyToProcess(couple.second)).find("kaon+Inelastic")  != std::string::npos) InteractionPointType.push_back(11);           
                    if ((truetraj.KeyToProcess(couple.second)).find("protonInelastic") != std::string::npos) InteractionPointType.push_back(13);           
                    if ((truetraj.KeyToProcess(couple.second)).find("neutronInelastic")!= std::string::npos) InteractionPointType.push_back(2);           
                    
                }
            }	
            
            iPrim++;
        }//<--End if primary
            } //geant particles
    // ######################################
    // ### Making a vector of MCtruth     ###
    // ######################################
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
    art::fill_ptr_vector(mclist, mctruthListHandle);
    
    
    mcevts_truth=mclist.size();
    if (mcevts_truth){
        art::Ptr<simb::MCTruth> mctruth = mclist[0];
        if (mctruth->Origin() == simb::kBeamNeutrino){
            nuPDG_truth  = mctruth->GetNeutrino().Nu().PdgCode();
            ccnc_truth   = mctruth->GetNeutrino().CCNC();
            mode_truth   = mctruth->GetNeutrino().Mode();            //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
            Q2_truth     = mctruth->GetNeutrino().QSqr();            //Momentum transfer squared
            W_truth      = mctruth->GetNeutrino().W();
            X_truth      = mctruth->GetNeutrino().X();
            Y_truth      = mctruth->GetNeutrino().Y();
            hitnuc_truth = mctruth->GetNeutrino().HitNuc();
            target_truth = mctruth->GetNeutrino().Target();
            enu_truth    = mctruth->GetNeutrino().Nu().E();          //true neutrino energy
            nuvtxx_truth = mctruth->GetNeutrino().Nu().Vx();
            nuvtxy_truth = mctruth->GetNeutrino().Nu().Vy();
            nuvtxz_truth = mctruth->GetNeutrino().Nu().Vz();
            
            // === neutrino dcos x, y, and z ====
            if (mctruth->GetNeutrino().Nu().P()){
                nu_dcosx_truth = mctruth->GetNeutrino().Nu().Px()/mctruth->GetNeutrino().Nu().P();
                nu_dcosy_truth = mctruth->GetNeutrino().Nu().Py()/mctruth->GetNeutrino().Nu().P();
                nu_dcosz_truth = mctruth->GetNeutrino().Nu().Pz()/mctruth->GetNeutrino().Nu().P();
            }
            // === lepton dcos x, y, and z ===
            lep_mom_truth = mctruth->GetNeutrino().Lepton().P();   //lepton momentum
            if (mctruth->GetNeutrino().Lepton().P()){
                lep_dcosx_truth = mctruth->GetNeutrino().Lepton().Px()/mctruth->GetNeutrino().Lepton().P();
                lep_dcosy_truth = mctruth->GetNeutrino().Lepton().Py()/mctruth->GetNeutrino().Lepton().P();
                lep_dcosz_truth = mctruth->GetNeutrino().Lepton().Pz()/mctruth->GetNeutrino().Lepton().P();
            }
            
            if (mctruth->NParticles()){
                simb::MCParticle particle = mctruth->GetParticle(0);
                t0_truth = particle.T();
            }
            
        }   //if (mctruth->Origin() == simb::kBeamNeutrino)
    
    }    //if (mcevts_truth)
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++nuetrino interaction selection++++++++++++++++++++++++++++++++
    int N_QE_number = 0;  //number of nuetron on QE interaction
    int P_QE_number = 0;  //number of proton on QE interaction
    int Primary_track_number =0; // initializing the number of primary track
    int Pi_Plus_non_QE_number = 0; //number of Pion_plus on non_QE interaction
    int Pi_Minus_non_QE_number = 0; //number of Pion_minus on non_QE interaction
    int Pi_Plus_Minus_non_QE_number = 0; //number of Pion_plus/Pion_minus on non_QE interaction
    int Pi_Zero_non_QE_number = 0; //number of Pion_Zero on non_QE interaction
    int N_non_QE_number = 0; //number of nuetron on non_QE interaction
    int P_non_QE_number = 0; //number of proton on non_QE interaction
    int Pi_Plus_Minus_NC_number = 0; //number of Pion_plus/Pion_minus on NC interaction
    int Pi_Zero_NC_number = 0; //number of Pion_Zero on NC interaction
    //////////numu///////////////////
    if (pdg[0] == 13 && ccnc_truth ==0) {
        if (mode_truth == 0) {
            for (int PrimaryI = 1; PrimaryI < geant_list_size; PrimaryI++) {
                if (process_primary[PrimaryI]== 1) {
                    Primary_track_number = TrackId[PrimaryI];
                    if (pdg[PrimaryI]== 2112) {
                        N_QE_number = N_QE_number+1;
                    } // nuetron counts
                    if (pdg[PrimaryI]== 2212) {
                        P_QE_number = P_QE_number+1;
                    } // proton counts
                } // process primary
            } // loop over all the track of the current event
            if (P_QE_number == 1) {
                class_truth = 2;
            } // proton number_QE ==1
            else if (N_QE_number == Primary_track_number-1) {
                class_truth =1;
            } // all the other particles are
	    else {
		 class_truth =7;
		}
        } // QE
        if (mode_truth != 0) {
            for (int PrimaryI = 1; PrimaryI < geant_list_size; PrimaryI++) {
                if (process_primary[PrimaryI]== 1) {
                    if (pdg[PrimaryI]== 211) {
                        Pi_Plus_non_QE_number = Pi_Plus_non_QE_number+1;
                    } // Pion+ counts
                    if (pdg[PrimaryI]== -211) {
                        Pi_Minus_non_QE_number = Pi_Minus_non_QE_number+1;
                    } //Pion- counts
                    if (abs(pdg[PrimaryI])== 211) {
                        Pi_Plus_Minus_non_QE_number = Pi_Plus_Minus_non_QE_number+1;
                    } // Pion+/Pion- counts
                    if (pdg[PrimaryI]== 111) {
                        Pi_Zero_non_QE_number = Pi_Zero_non_QE_number+1;
                    } // Pion0 counts
                    if (pdg[PrimaryI]== 2112) {
                        N_non_QE_number = N_non_QE_number+1;
                    } // nuetron counts
                    if (pdg[PrimaryI]== 2212) {
                        P_non_QE_number = P_non_QE_number+1;
                    } // Prpton counts
                } // process primary
            }//loop over all the track of the current event
            if (Pi_Plus_Minus_non_QE_number == 1 && Pi_Zero_non_QE_number == 1) {
                class_truth = 6;
            } // 2tracks_1mu+1Pi+/Pi-+1Pi0
            else if(Pi_Plus_Minus_non_QE_number == 2){
                class_truth = 4;
            } // 3tracks_1mu+2Pi(Pi+/Pi-)
            else if(Pi_Zero_non_QE_number ==1){
                class_truth = 5;
            }//2tracks_1mu+1Pi0
            else if(Pi_Plus_Minus_non_QE_number == 1){
                class_truth = 3;
            } //2tracks_1mu+1Pi+/Pi-
            else{
                class_truth = 7;
            }// others
        }// non_QE
    }// numu pdg[0] ==13 && CC interaction
    //////////nue///////////////////
    if (pdg[0] == 11 && ccnc_truth== 0) {
        class_truth =8;
    } //nue pdg[0] == 11 && CC (Wrong sign numu CC!!!)
    if ( ccnc_truth ==1) {
        for (int PrimaryI = 1; PrimaryI < geant_list_size; PrimaryI++) {
            if (process_primary[PrimaryI]== 1) {
                if (abs(pdg[PrimaryI])== 211) {
                    Pi_Plus_Minus_NC_number = Pi_Plus_Minus_NC_number+1;
                } //Pion+/Pion- counts
                if (pdg[PrimaryI]== 111) {
                    Pi_Zero_NC_number = Pi_Zero_NC_number+1;
                } //Pion0 counts
            } //Process Primary
        } //loop over all tracks in the current event
        if (Pi_Plus_Minus_NC_number == 1) {
            class_truth = 11;
        } //1Pi+/Pi-
        else if (Pi_Zero_NC_number == 1) {
            class_truth = 12;
        } //1Pi0
        else{
            class_truth = 13;
        } // other
    } //either numu-nue && NC interaction
    
            
            
  fTree->Fill();

}

void WireDump::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("anatree","analysis tree");
    fTree->Branch("run",&run,"run/I");
    fTree->Branch("subrun",&subrun,"subrun/I");
    fTree->Branch("event",&event,"event/I");
    fTree->Branch("nSamples", &nSamples, "nSamples/I");
    fTree->Branch("nWires", &nWires, "nWires/I");
    fTree->Branch("WireADC", WireADC, "WireADC[nWires][4500]/F");
    
    // === Geant particle info ==
    fTree->Branch("no_primaries",&no_primaries,"no_primaries/I");
    fTree->Branch("geant_list_size",&geant_list_size,"geant_list_size/I");
    fTree->Branch("pdg",pdg,"pdg[geant_list_size]/I");
    fTree->Branch("Eng",Eng,"Eng[geant_list_size]/F");
    fTree->Branch("Px",Px,"Px[geant_list_size]/F");
    fTree->Branch("Py",Py,"Py[geant_list_size]/F");
    fTree->Branch("Pz",Pz,"Pz[geant_list_size]/F");
    fTree->Branch("StartPointx",StartPointx,"StartPointx[geant_list_size]/F");
    fTree->Branch("StartPointy",StartPointy,"StartPointy[geant_list_size]/F");
    fTree->Branch("StartPointz",StartPointz,"StartPointz[geant_list_size]/F");
    fTree->Branch("EndPointx",EndPointx,"EndPointx[geant_list_size]/F");
    fTree->Branch("EndPointy",EndPointy,"EndPointy[geant_list_size]/F");
    fTree->Branch("EndPointz",EndPointz,"EndPointz[geant_list_size]/F");
    fTree->Branch("Startdcosx",Startdcosx,"Startdcosx[geant_list_size]/F");
    fTree->Branch("Startdcosy",Startdcosy,"Startdcosy[geant_list_size]/F");
    fTree->Branch("Startdcosz",Startdcosz,"Startdcosz[geant_list_size]/F");
    fTree->Branch("NumberDaughters",NumberDaughters,"NumberDaughters[geant_list_size]/I");
    fTree->Branch("Mother",Mother,"Mother[geant_list_size]/I");
    fTree->Branch("TrackId",TrackId,"TrackId[geant_list_size]/I");
    fTree->Branch("process_primary",process_primary,"process_primary[geant_list_size]/I");
    //+++++++++++++++++++Track_lenght+++++++++++++++++++++++++++++++++++
    fTree->Branch("Track_length",Track_length,"Track_length[geant_list_size]/F");
    //fTree->Branch("G4Process",&G4Process);//,"G4Process[geant_list_size]");
    //fTree->Branch("G4FinalProcess",&G4FinalProcess);//,"G4FinalProcess[geant_list_size]");
    fTree->Branch("NTrTrajPts",NTrTrajPts,"NTrTrajPts[no_primaries]/I");
    fTree->Branch("MidPosX",MidPosX,"MidPosX[no_primaries][5000]/D");
    fTree->Branch("MidPosY",MidPosY,"MidPosY[no_primaries][5000]/D");
    fTree->Branch("MidPosZ",MidPosZ,"MidPosZ[no_primaries][5000]/D");
    fTree->Branch("MidPx",MidPx,"MidPx[no_primaries][5000]/D");
    fTree->Branch("MidPy",MidPy,"MidPy[no_primaries][5000]/D");
    fTree->Branch("MidPz",MidPz,"MidPz[no_primaries][5000]/D");
    fTree->Branch("InteractionPoint"         ,&InteractionPoint         );
    fTree->Branch("InteractionPointType"     ,&InteractionPointType     );
    

    
    
    
    // === Genie nu info ==
    fTree->Branch("mcevts_truth",&mcevts_truth,"mcevts_truth/I");
    fTree->Branch("nuPDG_truth",&nuPDG_truth,"nuPDG_truth/I");
    fTree->Branch("ccnc_truth",&ccnc_truth,"ccnc_truth/I");
    fTree->Branch("mode_truth",&mode_truth,"mode_truth/I");
    fTree->Branch("enu_truth",&enu_truth,"enu_truth/F");
    fTree->Branch("Q2_truth",&Q2_truth,"Q2_truth/F");
    fTree->Branch("W_truth",&W_truth,"W_truth/F");
    fTree->Branch("X_truth",&X_truth,"X_truth/F");
    fTree->Branch("Y_truth",&Y_truth,"Y_truth/F");
    fTree->Branch("hitnuc_truth",&hitnuc_truth,"hitnuc_truth/I");
    fTree->Branch("target_truth",&target_truth,"target_truth/I");
    fTree->Branch("nuvtxx_truth",&nuvtxx_truth,"nuvtxx_truth/F");
    fTree->Branch("nuvtxy_truth",&nuvtxy_truth,"nuvtxy_truth/F");
    fTree->Branch("nuvtxz_truth",&nuvtxz_truth,"nuvtxz_truth/F");
    fTree->Branch("nu_dcosx_truth",&nu_dcosx_truth,"nu_dcosx_truth/F");
    fTree->Branch("nu_dcosy_truth",&nu_dcosy_truth,"nu_dcosy_truth/F");
    fTree->Branch("nu_dcosz_truth",&nu_dcosz_truth,"nu_dcosz_truth/F");
    fTree->Branch("lep_mom_truth",&lep_mom_truth,"lep_mom_truth/F");
    fTree->Branch("lep_dcosx_truth",&lep_dcosx_truth,"lep_dcosx_truth/F");
    fTree->Branch("lep_dcosy_truth",&lep_dcosy_truth,"lep_dcosy_truth/F");
    fTree->Branch("lep_dcosz_truth",&lep_dcosz_truth,"lep_dcosz_truth/F");
    fTree->Branch("t0_truth",&t0_truth,"t0_truth/F");
    //++++++++++++++++++++++++++++++++++++++++++++++
    fTree->Branch("class_truth",&class_truth,"class_truth/I");
    
 
}


void WireDump::endJob()
{
  // Implementation of optional member function here.
}

void WireDump::ResetVars()
{
    InteractionPoint.clear();
    InteractionPointType.clear();
    
    run = -99999;
    subrun = -99999;
    event = -99999;
    
    mcevts_truth = -9999;
    nuPDG_truth = -9999;
    ccnc_truth = -9999;
    mode_truth = -9999;
    enu_truth = -9999;
    Q2_truth = -9999;
    W_truth = -9999;
    X_truth = -9999;
    Y_truth = -9999;
    hitnuc_truth = -9999;
    target_truth = -9999;
    nuvtxx_truth = -9999;
    nuvtxy_truth = -9999;
    nuvtxz_truth = -9999;
    nu_dcosx_truth = -9999;
    nu_dcosy_truth = -9999;
    nu_dcosz_truth = -9999;
    lep_mom_truth = -9999;
    lep_dcosx_truth = -9999;
    lep_dcosy_truth = -9999;
    lep_dcosz_truth = -9999;
    t0_truth = -9999;
    class_truth = -99999;
  
    for(int j = 0; j < kMaxWires; j++){
        for(int i = 0; i < kMaxTicks; i++){
            WireADC[j][i] = 0;
        }//<---End i loop
    }//<---End j loop
    
    no_primaries = -99999;
    geant_list_size=-9999;
    for (int i = 0; i<kMaxPrimaries; ++i)
    {
        pdg[i] = -99999;
        Eng[i] = -99999;
        Px[i] = -99999;
        Py[i] = -99999;
        Pz[i] = -99999;
        StartPointx[i] = -99999;
        StartPointy[i] = -99999;
        StartPointz[i] = -99999;
        EndPointx[i] = -99999;
        EndPointy[i] = -99999;
        EndPointz[i] = -99999;
        Startdcosx[i]= -99999;
        Startdcosy[i]= -99999;
        Startdcosz[i]= -99999;
        NumberDaughters[i] = -99999;
        Mother[i] = -99999;
        TrackId[i] = -99999;
        process_primary[i] = -99999;
        xi[i] = -9999999.99999999;
        yi[i] = -9999999.99999999;
        zi[i] = -9999999.99999999;
        x_i[i] = -9999999.99999999;
        y_i[i] = -9999999.99999999;
        z_i[i] = -9999999.99999999;
        Track_length[i] = -9999999.9999999;
        
    } //<---End kMaxPrimaries
    
    for(int i = 0; i<kMaxPrimaryPart; i++){
        NTrTrajPts[i] = -99999;
        for(int j = 0; j<kMaxTruePrimaryPts; j++){
            MidPosX[i][j] = -99999;
            MidPosY[i][j] = -99999;
            MidPosZ[i][j] = -99999;
            MidPx[i][j] = -99999;
            MidPy[i][j] = -99999;
            MidPz[i][j] = -99999;
        }
    }

}



DEFINE_ART_MODULE(WireDump)
