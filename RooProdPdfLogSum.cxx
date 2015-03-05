/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id: RooProdPdfLogSum.cxx,v 1.2 2012/11/20 16:42:45 kreis Exp $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// RooProdPdf is an efficient implementation of a product of PDFs of the form 
//
//  PDF_1 * PDF_2 * ... * PDF_N
//
// PDFs may share observables. If that is the case any irreducable subset
// of PDFS that share observables will be normalized with explicit numeric
// integration as any built-in normalization will no longer be valid.
//
// Alternatively, products using conditional PDFs can be defined, e.g.
//
//    F(x|y) * G(y)
//
// meaning a pdf F(x) _given_ y and a PDF G(y). In this contruction F is only
// normalized w.r.t x and G is normalized w.r.t y. The product in this construction
// is properly normalized.
//
// If exactly one of the component PDFs supports extended likelihood fits, the
// product will also be usable in extended mode, returning the number of expected
// events from the extendable component PDF. The extendable component does not
// have to appear in any specific place in the list.
// 
// END_HTML
//

#include "RooFit.h"
#include "Riostream.h"
#include "TClass.h"

#include "TIterator.h"
#include "RooProdPdfLogSum.h"
#include "RooRealProxy.h"
#include "RooProdGenContext.h"
#include "RooGenProdProj.h"
#include "RooProduct.h"
#include "RooNameReg.h"
#include "RooMsgService.h"
#include "RooFormulaVar.h"
#include "RooRealVar.h"
#include "RooAddition.h"
#include "RooGlobalFunc.h"
#include "RooConstVar.h"
#include "RooWorkspace.h"
#include "RooRangeBoolean.h"
#include "RooCustomizer.h"
#include "RooRealIntegral.h"

#include <string.h>
#include <sstream>
#include "TSystem.h"

ClassImp(RooProdPdfLogSum)
;



//_____________________________________________________________________________
RooProdPdfLogSum::RooProdPdfLogSum() : RooProdPdf(),_logValue(0) {}

//_____________________________________________________________________________
RooProdPdfLogSum::RooProdPdfLogSum(const char *name, const char *title, Double_t cutOff) : 
  RooProdPdf(name, title, cutOff),_logValue(0) {} 

//_____________________________________________________________________________
RooProdPdfLogSum::RooProdPdfLogSum(const char *name, const char *title,
		       RooAbsPdf& pdf1, RooAbsPdf& pdf2, Double_t cutOff) : 
  RooProdPdf(name, title, pdf1, pdf2, cutOff),_logValue(0) {}

//_____________________________________________________________________________
RooProdPdfLogSum::RooProdPdfLogSum(const char* name, const char* title, const RooArgList& inPdfList, Double_t cutOff) :
  RooProdPdf(name, title, inPdfList, cutOff),_logValue(0) {}
//_____________________________________________________________________________
RooProdPdfLogSum::RooProdPdfLogSum(const char* name, const char* title, const RooArgSet& fullPdfSet,
		       const RooCmdArg& arg1, const RooCmdArg& arg2,
		       const RooCmdArg& arg3, const RooCmdArg& arg4,
		       const RooCmdArg& arg5, const RooCmdArg& arg6,
		       const RooCmdArg& arg7, const RooCmdArg& arg8) :
  RooProdPdf(name, title, fullPdfSet, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8),_logValue(0) {}

//_____________________________________________________________________________
RooProdPdfLogSum::RooProdPdfLogSum(const char* name, const char* title,
		       const RooCmdArg& arg1, const RooCmdArg& arg2,
		       const RooCmdArg& arg3, const RooCmdArg& arg4,
		       const RooCmdArg& arg5, const RooCmdArg& arg6,
		       const RooCmdArg& arg7, const RooCmdArg& arg8) :
  RooProdPdf(name, title, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8),_logValue(0) {}

//_____________________________________________________________________________
RooProdPdfLogSum::RooProdPdfLogSum(const char* name, const char* title, const RooArgSet& fullPdfSet, const RooLinkedList& cmdArgList) :
  RooProdPdf(name, title, fullPdfSet, cmdArgList),_logValue(0) {}

//_____________________________________________________________________________
RooProdPdfLogSum::RooProdPdfLogSum(const RooProdPdfLogSum& other, const char* name) :
  RooProdPdf(other,name),_logValue(0) {}

//_____________________________________________________________________________
RooProdPdfLogSum::~RooProdPdfLogSum() 
{
}

//_____________________________________________________________________________
Double_t RooProdPdfLogSum::getLogVal(const RooArgSet* nset) const 
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

//_____________________________________________________________________________
Double_t RooProdPdfLogSum::calculate(const RooArgList* partIntList, const RooLinkedList* normSetList) const
{
  // Calculate running product of pdfs terms, using the supplied
  // normalization set in 'normSetList' for each component

  RooAbsReal* partInt ;
  RooArgSet* normSet ;
  Double_t value(1.0),logValue(0.0) ;
  Int_t n = partIntList->getSize() ;

  Int_t i ;
  for (i=0 ; i<n ; i++) {
    partInt = ((RooAbsReal*)partIntList->at(i)) ;
    normSet = ((RooArgSet*)normSetList->At(i)) ;    
    Double_t piLogVal = dynamic_cast<RooAbsPdf*>(partInt)->getLogVal(normSet->getSize()>0 ? normSet : 0) ;
    //cout << "partInt(" << partInt->GetName() << ") = " << piVal << " normSet = " << normSet << " " << (normSet->getSize()>0 ? *normSet : RooArgSet()) << endl ;
    logValue += piLogVal ;
    //if (value<_cutOff) {
    //  break ;
    //}
  }
  _logValue = logValue;
  value = exp(logValue);
  return value ;
}



//_____________________________________________________________________________
Double_t RooProdPdfLogSum::calculate(const RooProdPdfLogSum::CacheElem& cache, Bool_t /*verbose*/) const 
{
  // Calculate running product of pdfs terms, using the supplied
  // normalization set in 'normSetList' for each component

  //cout << "RooProdPdfLogSum::calculate from cache" << endl ;

  Double_t value(1.0) ,logValue(0.0);

  if (cache._isRearranged) {

    if (dologD(Eval)) {
      cxcoutD(Eval) << "RooProdPdfLogSum::calculate(" << GetName() << ") rearranged product calculation" 
                    << " calculate: num = " << cache._rearrangedNum->GetName() << " = " << cache._rearrangedNum->getVal() << endl ;
//       cache._rearrangedNum->printComponentTree("",0,5) ;    
      cxcoutD(Eval) << "calculate: den = " << cache._rearrangedDen->GetName() << " = " << cache._rearrangedDen->getVal() << endl ;      
//       cache._rearrangedDen->printComponentTree("",0,5) ;    
    }

    value = cache._rearrangedNum->getVal() / cache._rearrangedDen->getVal() ;
    
  } else {
    
    cxcoutD(Eval) << "RooProdPdfLogSum::calculate(" << GetName() << ") regular product chain calculation" << endl ;
    
    RooAbsReal* partInt ;
    RooArgSet* normSet ;
    Int_t n = cache._partList.getSize() ;

    //cache._partList.Print("v");
    //cache._normList.Print("v");
    Int_t i ;
    //cout << "in the calculate loop" << endl;
    for (i=0 ; i<n ; i++) {
      //cout << "checking partInt" << endl;
      partInt = ((RooAbsReal*)cache._partList.at(i)) ;
      //cout << "checking normSet" << endl;
      normSet = ((RooArgSet*)cache._normList.At(i)) ;    
      Double_t piLogVal = dynamic_cast<RooAbsPdf*>(partInt)->getLogVal(normSet->getSize()>0 ? normSet : 0) ;
      //cout << "partInt " << partInt->GetName() << " piLogVal is " << piLogVal << endl;
      //cout << "partInt " << partInt->GetName() << " is of type " << partInt->IsA()->GetName() << endl ;
      if (dynamic_cast<RooAbsPdf*>(partInt)) {
	cxcoutD(Eval) << "product term " << partInt->GetName() << " normalized over " << (normSet?*normSet:RooArgSet())  
		      << " = " << partInt->getVal() << " / " << ((RooAbsPdf*)partInt)->getNorm(normSet) << " = " << exp(piLogVal) << endl ;
      } else {
	cxcoutD(Eval) << "product term " << partInt->GetName() << " normalized over " << (normSet?*normSet:RooArgSet()) << " = " << exp(piLogVal) << endl ;
      }
      logValue += piLogVal ;
      //if (value<_cutOff) {
      //	break ;
      //}
    }
    _logValue = logValue;
    value = exp(logValue);
  }

  cxcoutD(Eval) << "return value = " << value << endl ;
  return value ;
}


//_____________________________________________________________________________

Double_t RooProdPdfLogSum::evaluate() const 
{
  // Calculate current value of object

  Int_t code ;
  CacheElem* cache = (CacheElem*) _cacheMgr.getObj(_curNormSet,0,&code) ;
  
  // If cache doesn't have our configuration, recalculate here
  if (!cache) {
    RooArgList *plist(0) ;
    RooLinkedList *nlist(0) ;
    getPartIntList(_curNormSet,0,plist,nlist,code) ;
    cache = (CacheElem*) _cacheMgr.getObj(_curNormSet,0,&code) ;
  }

  return calculate(*cache) ;
}

//_____________________________________________________________________________

