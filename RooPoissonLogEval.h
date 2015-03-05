 /***************************************************************************** 
  * Project: RooFit                                                           * 
  *                                                                           * 
  * Simple Poisson PDF
  * author: Kyle Cranmer <cranmer@cern.ch>
  *                                                                           * 
  *****************************************************************************/ 

#ifndef ROOPOISSONLOGEVAL
#define ROOPOISSONLOGEVAL

#include "RooPoisson.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooPoissonLogEval : public RooPoisson {
public:
  RooPoissonLogEval() { _noRounding = kFALSE ; } ;
  RooPoissonLogEval(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _mean, Bool_t noRounding=kFALSE);
  RooPoissonLogEval(const RooPoissonLogEval& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooPoissonLogEval(*this,newname); }
  inline virtual ~RooPoissonLogEval() { }

  virtual Double_t getValV(const RooArgSet* set=0) const ;
  virtual Double_t getLogVal(const RooArgSet* set=0) const ;

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const {return 1;}
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const {if(code == 1) return 1.;else{assert(0);return 0;}}
  //
  //Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const;
  //void generateEvent(Int_t code);

protected:

  //RooRealProxy x ;
  //RooRealProxy mean ;
  //Bool_t  _noRounding ;
  
  Double_t evaluate() const ;
  Double_t evaluate(Double_t k) const;

  mutable Double_t _logValue ;
  

private:

  ClassDef(RooPoissonLogEval,2) // A Poisson PDF
};
 
#endif
