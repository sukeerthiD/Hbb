#include "Hbb.h"
#include <iostream>
#include <vector>
#include <string>
#include <TTree.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
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
            "MCParticleCollection" ,
            "Name of the MCParticle collection"  ,
            _colMCP , //my collection name
            std::string("MCParticle") //take this string from input slcio
    );

    registerOutputCollection( LCIO::MCPARTICLE,
		        "HbbCol" ,
		         "Higgs to bb event collection"  ,
		         _HbbCol,
            std::string("HbbCol") ) ;
}



void Hbb::init() {

    streamlog_out(DEBUG) << "   init called  " << std::endl ;
    printParameters() ;  // usually a good idea to
    _nRun = 0 ;
    _nEvt = 0 ;
    nH=0;
    nb=0;

}


void Hbb::processRunHeader( LCRunHeader* run) {
    _nRun++ ;
}



void Hbb::processEvent( LCEvent * evt ) {


        streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber()
            << "   in run:  " << evt->getRunNumber() << std::endl ;
    // this gets called for every event
    // usually the working horse ...

 LCCollectionVec *HbbCol = new LCCollectionVec( LCIO::MCPARTICLE );
 HbbCol->setSubset( true );

    LCCollection* col{};
    col = evt->getCollection( _colMCP ) ;

      if( col == 0 ){
        streamlog_out(DEBUG)  << "empty collection******************" << std::endl ;
      }
    else { //collection


    int nMCP = col->getNumberOfElements()  ;
    for(int i=0; i< nMCP ; i++){ //loop through all elements

      MCParticle* p = dynamic_cast<MCParticle*>( col->getElementAt( i ) ) ;

      if (p->getPDG() ==25) { // if the MCParticle p is Higgs

        nH+=1;
        MCParticleVec Hd = p->getDaughters() ; //get daughters of the Higgs
        for (int j=0; j<Hd.size(); j++){ //loop through all daughters
           MCParticle* Hda = Hd[j];

           if (Hda->getPDG()==5) { //if Higgs decays to bb
             streamlog_out(DEBUG)  << "Higgs decays to a b quark" << std::endl ;
             nb+=1;
             HbbCol->addElement( p );
           }//end if Higgs decays to bb

        }//end-loop through all daughters
      }// if the MCParticle p is Higgs

    }//end-loop through all elements
  }//collection

    _nEvt ++ ;
 evt->addCollection( HbbCol, _HbbCol.c_str() );
} //end of the event



void Hbb::check( LCEvent * evt ) {
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void Hbb::end(){
    streamlog_out(DEBUG)  << "number of Higgs: " <<nH/2 <<  ", in which number of H->bb: " <<nb << std::endl ;
}
