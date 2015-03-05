/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooProdPdfLogSum.h,v 1.1 2012/10/15 08:32:20 kreis Exp $
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
#ifndef ROO_PROD_PDF_LOG_SUM
#define ROO_PROD_PDF_LOG_SUM

#include "Riosfwd.h"
#include "RooAbsPdf.h"
#include "RooProdPdf.h"
#include "RooListProxy.h"
#include "RooLinkedList.h"
#include "RooAICRegistry.h"
#include "RooCacheManager.h"
#include "RooObjCacheManager.h"
#include "RooCmdArg.h"
#include <vector>

typedef RooArgList* pRooArgList ;
typedef RooLinkedList* pRooLinkedList ;

class RooProdPdfLogSum : public RooProdPdf {
public:
  RooProdPdfLogSum() ;
  RooProdPdfLogSum(const char *name, const char *title, Double_t cutOff=0);
  RooProdPdfLogSum(const char *name, const char *title,
	    RooAbsPdf& pdf1, RooAbsPdf& pdf2, Double_t cutOff=0) ;
  RooProdPdfLogSum(const char* name, const char* title, const RooArgList& pdfList, Double_t cutOff=0) ;
  RooProdPdfLogSum(const char* name, const char* title, const RooArgSet& fullPdfSet, const RooLinkedList& cmdArgList) ;

  RooProdPdfLogSum(const char* name, const char* title, const RooArgSet& fullPdfSet,
   	     const RooCmdArg& arg1            , const RooCmdArg& arg2=RooCmdArg(),
             const RooCmdArg& arg3=RooCmdArg(), const RooCmdArg& arg4=RooCmdArg(),
             const RooCmdArg& arg5=RooCmdArg(), const RooCmdArg& arg6=RooCmdArg(),
             const RooCmdArg& arg7=RooCmdArg(), const RooCmdArg& arg8=RooCmdArg()) ;

  RooProdPdfLogSum(const char* name, const char* title, 
             const RooCmdArg& arg1,             const RooCmdArg& arg2=RooCmdArg(),
             const RooCmdArg& arg3=RooCmdArg(), const RooCmdArg& arg4=RooCmdArg(),
             const RooCmdArg& arg5=RooCmdArg(), const RooCmdArg& arg6=RooCmdArg(),
             const RooCmdArg& arg7=RooCmdArg(), const RooCmdArg& arg8=RooCmdArg()) ;

  RooProdPdfLogSum(const RooProdPdfLogSum& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooProdPdfLogSum(*this,newname) ; }
  virtual ~RooProdPdfLogSum() ;

  virtual Double_t getLogVal(const RooArgSet* set=0) const ;  
protected:

  //// The cache object
  //class CacheElem : public RooAbsCacheElement {
  //public:
  //  CacheElem() : _isRearranged(kFALSE), _rearrangedNum(0), _rearrangedDen(0) {} 
  //  virtual ~CacheElem() ;
  //  // Payload
  //  RooArgList _partList ;
  //  RooArgList _numList ;
  //  RooArgList _denList ;
  //  RooArgList _ownedList ;
  //  RooLinkedList _normList ;    
  //  Bool_t _isRearranged ;
  //  RooAbsReal* _rearrangedNum ;
  //  RooAbsReal* _rearrangedDen ;
  //  // Cache management functions
  //  virtual RooArgList containedArgs(Action) ;
  //  virtual void printCompactTreeHook(std::ostream&, const char *, Int_t, Int_t) ;
  //private:
  //  CacheElem(const CacheElem&) ;
  //} ;
  //mutable RooObjCacheManager _cacheMgr ; // The cache manager

  Double_t calculate(const RooProdPdf::CacheElem& cache, Bool_t verbose=kFALSE) const ;
  Double_t calculate(const RooArgList* partIntList, const RooLinkedList* normSetList) const ;
  mutable Double_t _logValue;
  Double_t evaluate() const;

private:

  ClassDef(RooProdPdfLogSum,4) // PDF representing a product of PDFs
};


#endif
