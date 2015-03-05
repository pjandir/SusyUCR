 /***************************************************************************** 
  * Project: RooFit                                                           * 
  *                                                                           * 
  * Simple Poisson PDF
  * author: Kyle Cranmer <cranmer@cern.ch>
  *                                                                           * 
  *****************************************************************************/ 

//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Poisson pdf
// END_HTML
//

#include <iostream> 

#include "RooPoissonLogEval.h" 
#include "RooAbsCategory.h" 

#include "RooRandom.h"
#include "RooMath.h"
#include "TMath.h"

#include <Math/SpecFuncMathCore.h>
#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFuncMathCore.h>

ClassImp(RooPoissonLogEval) 



//_____________________________________________________________________________
RooPoissonLogEval::RooPoissonLogEval(const char *name, const char *title, 
		       RooAbsReal& _x,
		       RooAbsReal& _mean,
		       Bool_t noRounding) :
   RooPoisson(name,title, _x, _mean, noRounding),_logValue(0)
{ 
  // Constructor
} 



//_____________________________________________________________________________
 RooPoissonLogEval::RooPoissonLogEval(const RooPoissonLogEval& other, const char* name) :  
   RooPoisson(other,name),_logValue(0)
   //x("x",this,other.x),
   //mean("mean",this,other.mean),
   //_noRounding(other._noRounding)
{ 
   // Copy constructor
} 




//_____________________________________________________________________________
Double_t RooPoissonLogEval::evaluate() const 
{ 
  // Implementation in terms of the TMath Poisson function
  //cout << "Trying stuff" << endl;
  //cout << "_noRounding" << _noRounding << endl;
  //cout << "x" << x << endl;
  //cout << "mean" << mean << endl;
  //cout << "floor(x)" << floor(x) << endl;
  Double_t k = _noRounding ? x : floor(x);  
  //return TMath::Poisson(k,mean) ;
  if (k<0)
    {
      _logValue = log(0.0);
      return 0;
    }
  else if (k == 0.0)
    {
      _logValue = -mean ;
      return 1./exp(mean);
    }
  else {
    Double_t lnpoisson = k*log(mean)-mean-::ROOT::Math::lgamma(k+1.);
    //cout << "the log value is " << lnpoisson << endl;
    _logValue = lnpoisson;
    return exp(lnpoisson);
  }
  
} 

Double_t RooPoissonLogEval::getValV(const RooArgSet* nset) const
{
  // Return current value, normalizated by integrating over
  // the observables in 'nset'. If 'nset' is 0, the unnormalized value. 
  // is returned. All elements of 'nset' must be lvalues
  //
  // Unnormalized values are not cached
  // Doing so would be complicated as _norm->getVal() could
  // spoil the cache and interfere with returning the cached
  // return value. Since unnormalized calls are typically
  // done in integration calls, there is no performance hit.

  // Fast-track processing of clean-cache objects
  //   if (_operMode==AClean) {
  //     cout << "RooAbsPdf::getValV(" << this << "," << GetName() << ") CLEAN  value = " << _value << endl ;
  //     return _value ;
  //   }

  // Special handling of case without normalization set (used in numeric integration of pdfs)
  if (!nset) {
    RooArgSet* tmp = _normSet ;
    _normSet = 0 ;
    Double_t val = evaluate() ;
    _normSet = tmp ;
    Bool_t error = traceEvalPdf(val) ;

    if (error) {
//       raiseEvalError() ;
      return 0 ;
    }
    return val ;
  }


  // Process change in last data set used
  Bool_t nsetChanged(kFALSE) ;
  if (nset!=_normSet || _norm==0) {
    nsetChanged = syncNormalization(nset) ;
  }

  // Return value of object. Calculated if dirty, otherwise cached value is returned.
  if (isValueDirty() || nsetChanged || _norm->isValueDirty()) {

    // Evaluate numerator
    Double_t rawVal = evaluate() ;
    Bool_t error = traceEvalPdf(rawVal) ; // Error checking and printing

    // Evaluate denominator
    Double_t normVal(_norm->getVal()) ;
    
    if (normVal<=0.) {
      error=kTRUE ;
      logEvalError("p.d.f normalization integral is zero or negative") ;  
    }

    // Raise global error flag if problems occur
    if (error) {
//       raiseEvalError() ;
      _value = 0 ;
      _logValue = log(_value);
    } else {
      _value = rawVal / normVal ;
      _logValue = _logValue - log(normVal);
//       cout << "RooAbsPdf::getValV(" << GetName() << ") writing _value = " << _value << endl ;
    }

    clearValueAndShapeDirty() ; //setValueDirty(kFALSE) ;
  } 

  return _value ;
}

Double_t RooPoissonLogEval::getLogVal(const RooArgSet* nset) const 
{
  Double_t prob = getVal(nset) ;
  //evaluate();
  //cout << "returning log value is " << _logValue << endl;
  return _logValue;
  //if(prob < 0) {
  //
  //  logEvalError("getLogVal() top-level p.d.f evaluates to a negative number") ;
  //
  //  return 0;
  //}
  //if(prob == 0) {
  //
  //  logEvalError("getLogVal() top-level p.d.f evaluates to zero") ;
  //
  //  return log((double)0);
  //}
  //
  //if (TMath::IsNaN(prob)) {
  //  logEvalError("getLogVal() top-level p.d.f evaluates to NaN") ;
  //
  //  return log((double)0);
  //  
  //}
  //return log(prob);
}
