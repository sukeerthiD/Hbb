#include "Hbb.h"
#include <iostream>
#include <vector>
#include <string>
#include <TTree.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include "TROOT.h"
#include <IMPL/LCCollectionVec.h>
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/ITupleFactory.h>
#include <AIDA/ITuple.h>
//#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA




Hbb aHbb ;


Hbb::Hbb() : Processor("Hbb") {

    _description = "For now, selects H to bb events" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE,
            "MCParticleCollection" , //name
            "Name of the MCParticle collection"  , //description
            _colMCP , //my collection name
            std::string("MCParticle") //take this string from input slcio
    );

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
            "FitOutColection" , //name
            " Fit Output collection from ZH5CFit"  , //description
            _FitCol , //my collection name
            std::string("FitReco") //take this string from input slcio
    );

    registerInputCollection( LCIO::MCPARTICLE,
            "SemiLeptonicDecays" , //name
            "Semi-Leptonic Decays Collection"  , //description
            _SLDCol , //my collection name
            std::string("") //take this string from input slcio
    );

    registerProcessorParameter("outputFilename",
           "name of output file",
           _outfile,
            std::string("")
    );
            }



void Hbb::init() {

    streamlog_out(DEBUG) << "   init called  " << std::endl ;
    printParameters() ;  // usually a good idea to
    _nRun = 0 ;
    _nEvt = 0 ;
    nH=0;
    nb=0;
    nsld =0;
    nsldb =0;
    nsldc =0;
    _fout = new TFile(_outfile.c_str(),"recreate");

   HbbTree = new TTree("HbbTree","HbbTree");
   HbbTree->Branch("Hbbflag",&Hbbflag,"Hbbflag/I") ;
   HbbTree->Branch("hist_E_allH",&hist_E_allH,"hist_E_allH/F") ;
   HbbTree->Branch("hist_M_allH",&hist_M_allH,"hist_M_allH/F") ;
   HbbTree->Branch("hist_E_Hbb",&hist_E_Hbb,"hist_E_Hbb/F") ;
   HbbTree->Branch("hist_E_MCISR",&hist_E_MCISR,"hist_E_MCISR/F") ;
   HbbTree->Branch("hist_E_recISR",&hist_E_recISR,"hist_E_recISR/F") ;
   HbbTree->Branch("bestchisqu",&bestchisqu,"bestchisqu/F") ;
   HbbTree->Branch("Hmass",&Hmass,"Hmass/F") ;
   HbbTree->Branch("Zmass",&Zmass,"Zmass/F") ;
   HbbTree->Branch("bestprob",&bestprob,"bestprob/F") ;
   HbbTree->Branch("nSLDB",&nSLDB,"nSLDB/I") ;
   HbbTree->Branch("nSLDC",&nSLDC,"nSLDC/I") ;
   HbbTree->Branch("nSLDBC",&nSLDBC,"nSLDBC/I") ;
 }


void Hbb::processRunHeader( LCRunHeader* run) {
    _nRun++ ;
}

void Hbb::setZero()
  {
  hist_E_allH=0.0;
  hist_M_allH=0.0;
  hist_E_Hbb=0.0;
  hist_E_MCISR=0.0;
  hist_E_recISR=0.0;
  bestchisqu=0.0;
  bestprob =0.0;
  nSLDB =0.0;
  nSLDBC =0.0;
  nSLDC =0.0;

}

void Hbb::processEvent( LCEvent * evt ) {

  setZero();
  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber()
  << "   in run:  " << evt->getRunNumber() << std::endl ;
  // this gets called for every event
  // usually the working horse ...


  LCCollection* col{};
  col = evt->getCollection( _colMCP ) ;

  if( col == 0 ){
      streamlog_out(WARNING)  << "empty collection******************" << std::endl ;
      }
  else { //collection


      int nMCP = col->getNumberOfElements()  ;
      for(int i=0; i< nMCP ; i++){ //loop through all elements
          //nFlagHbb=0;
          MCParticle* p = dynamic_cast<MCParticle*>( col->getElementAt( i ) ) ;

          if (p->getPDG()==25) { // if the MCParticle p is Higgs
              if (i!=4) { //if Higgs is not generator element
                hist_E_allH=p->getEnergy();
                hist_M_allH=p->getMass();

                streamlog_out(DEBUG)  << "Higgs mass" << hist_M_allH << std::endl ;
                nH+=1;
                MCParticleVec Hd = p->getDaughters() ; //get daughters of the Higgs
                for (int j=0; j<Hd.size(); j++){ //loop through all daughters
                    MCParticle* Hda = Hd[j];
                    if (Hda->getPDG()==5 || Hda->getPDG()==-5) { //if Higgs decays to bb
                        nFlagHbb=1;
                        streamlog_out(DEBUG)  << "Higgs decays to a b quark" << nFlagHbb<< std::endl ;
                        nb+=1;
                        hist_E_Hbb=p->getEnergy();

                    }
                    else { //if Higgs doesnt decays to bb
                        nFlagHbb=0;
                        streamlog_out(DEBUG)  << "Higgs does not decay to a b quark" << nFlagHbb<< std::endl ;
                    }//end if Higgs does not decay to bb
                }//end-loop through all daughters
              }//end if Higgs not generator
          }// if the MCParticle p is Higgs


          if (p->getPDG()==22){ // if the MCParticle p is a photon
              //  if (p->getVertex()[0] == 0. && p->getVertex()[1]==0.) { //if ISR photon
              if (i==5 || i==6) { //if ISR photon
                streamlog_out(DEBUG)  << "******ISR***************at element: " <<i << std::endl ;
                streamlog_out(DEBUG) << " Vertex   = " << p->getVertex()[0] << " " << p->getVertex()[1] << " " << p->getVertex()[2] << std::endl;
                hist_E_MCISR+=p->getEnergy();
              } //endif ISR photon
          }// if the MCParticle p is a photon
      }//end-loop through all elements
  }//collection

  LCCollection* fitcol = NULL;
    try{
        fitcol = evt->getCollection( _FitCol );
        streamlog_out(WARNING) << _FitCol << " collection  available" << std::endl;
    }
    catch( lcio::DataNotAvailableException e )
        {
        streamlog_out(WARNING) << _FitCol << " collection not available****" << std::endl;
        fitcol = NULL;
    }
  int nfitrec = fitcol->getNumberOfElements()  ;

  for (int i=0; i<nfitrec; i++){//loop through all fit rec particles
      ReconstructedParticle* r = dynamic_cast<ReconstructedParticle*>(fitcol->getElementAt(i));
          if (r->getType()==22){//if it is Photon
          hist_E_recISR=r->getEnergy();
          streamlog_out(DEBUG) <<r->getType()  << " : IS THE ISR PHOTON?" << std::endl;
          }// end of if it is photon
          if (r->getType()==25){//if it is Higgs
          Hmass=r->getMass();
          streamlog_out(DEBUG) <<r->getMass()  << " :  HIGGS? Mass" << std::endl;
          }// end of if it is Higgs
          if (r->getType()==23){//if it is Z
          Zmass=r->getMass();
          streamlog_out(DEBUG) <<r->getMass()  << " :  Z? Mass" << std::endl;
          }// end of if it is Z
  }//end of for loop (fit rec particles)
  bestchisq=fitcol->getParameters().getFloatVal("bestchisq");
  streamlog_out(DEBUG) << bestchisq<< " bestchisquare" << std::endl;
  bestchisqu=fitcol->getParameters().getFloatVal("bestchisq");
  bestprob=fitcol->getParameters().getFloatVal("best_prob");
  Hbbflag=nFlagHbb;

  LCCollection* sldcol = NULL;
    try{
        sldcol = evt->getCollection( _SLDCol );
        streamlog_out(WARNING) << _SLDCol << " collection available*********" << std::endl;
    }
    catch( lcio::DataNotAvailableException e )
        {
        streamlog_out(WARNING) << _SLDCol << " collection not available****" << std::endl;
        sldcol = NULL;
    }
  nSLDB=sldcol->getParameters().getIntVal("nBSLD");
  nSLDC=sldcol->getParameters().getIntVal("nCSLD");
  nSLDBC=sldcol->getParameters().getIntVal("nSLD");
  //if(nFlagHbb==0){
      nsldb += nSLDB;
      nsldc += nSLDC;
      nsld += nSLDBC;
  //}
	streamlog_out(MESSAGE) << "Number of Semi-Leptonic decay of B-Hadron: " << nSLDB << std::endl;
	streamlog_out(MESSAGE) << "Number of Semi-Leptonic decay of C-Hadron: " << nSLDC << std::endl;
	streamlog_out(MESSAGE) << "Total Number of Semi-Leptonic decays: " << nSLDBC << std::endl;


  HbbTree->Fill();
  _nEvt ++ ;
  streamlog_out(DEBUG) << nFlagHbb  << " : is there a bb event" << std::endl;
} //end of the event



void Hbb::check( LCEvent * evt ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void Hbb::end(){
  streamlog_out(MESSAGE)  << "number of Higgs: " <<nH <<  ", in which number of H->bb: " <<nb/2 << std::endl ;
  streamlog_out(MESSAGE)  << "in Hbb, total number of sld: " << nsldb+nsldc <<  ", in which number of sld of B: " <<nsldb<<  ", and number of sld of C: " <<nsldc << std::endl ;
  _fout->Write(0);
  _fout->Close();
}
