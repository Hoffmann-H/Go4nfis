// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__nfis
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "./TnfisUnpackEvent.h"
#include "./nfisGlobals.h"
#include "./TnfisParam.h"
#include "./TnfisAnalysisFC.h"
#include "./TnfisUnpackProc.h"
#include "./nfisHistograms.h"
#include "./TnfisAnalysis.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TnfisUnpackEvent(void *p = 0);
   static void *newArray_TnfisUnpackEvent(Long_t size, void *p);
   static void delete_TnfisUnpackEvent(void *p);
   static void deleteArray_TnfisUnpackEvent(void *p);
   static void destruct_TnfisUnpackEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TnfisUnpackEvent*)
   {
      ::TnfisUnpackEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TnfisUnpackEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TnfisUnpackEvent", ::TnfisUnpackEvent::Class_Version(), "TnfisUnpackEvent.h", 23,
                  typeid(::TnfisUnpackEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TnfisUnpackEvent::Dictionary, isa_proxy, 4,
                  sizeof(::TnfisUnpackEvent) );
      instance.SetNew(&new_TnfisUnpackEvent);
      instance.SetNewArray(&newArray_TnfisUnpackEvent);
      instance.SetDelete(&delete_TnfisUnpackEvent);
      instance.SetDeleteArray(&deleteArray_TnfisUnpackEvent);
      instance.SetDestructor(&destruct_TnfisUnpackEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TnfisUnpackEvent*)
   {
      return GenerateInitInstanceLocal((::TnfisUnpackEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TnfisUnpackEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TnfisUnpFC(void *p = 0);
   static void *newArray_TnfisUnpFC(Long_t size, void *p);
   static void delete_TnfisUnpFC(void *p);
   static void deleteArray_TnfisUnpFC(void *p);
   static void destruct_TnfisUnpFC(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TnfisUnpFC*)
   {
      ::TnfisUnpFC *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TnfisUnpFC >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TnfisUnpFC", ::TnfisUnpFC::Class_Version(), "TnfisUnpackEvent.h", 39,
                  typeid(::TnfisUnpFC), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TnfisUnpFC::Dictionary, isa_proxy, 4,
                  sizeof(::TnfisUnpFC) );
      instance.SetNew(&new_TnfisUnpFC);
      instance.SetNewArray(&newArray_TnfisUnpFC);
      instance.SetDelete(&delete_TnfisUnpFC);
      instance.SetDeleteArray(&deleteArray_TnfisUnpFC);
      instance.SetDestructor(&destruct_TnfisUnpFC);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TnfisUnpFC*)
   {
      return GenerateInitInstanceLocal((::TnfisUnpFC*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TnfisUnpFC*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TnfisUnpPreAmp(void *p = 0);
   static void *newArray_TnfisUnpPreAmp(Long_t size, void *p);
   static void delete_TnfisUnpPreAmp(void *p);
   static void deleteArray_TnfisUnpPreAmp(void *p);
   static void destruct_TnfisUnpPreAmp(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TnfisUnpPreAmp*)
   {
      ::TnfisUnpPreAmp *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TnfisUnpPreAmp >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TnfisUnpPreAmp", ::TnfisUnpPreAmp::Class_Version(), "TnfisUnpackEvent.h", 60,
                  typeid(::TnfisUnpPreAmp), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TnfisUnpPreAmp::Dictionary, isa_proxy, 4,
                  sizeof(::TnfisUnpPreAmp) );
      instance.SetNew(&new_TnfisUnpPreAmp);
      instance.SetNewArray(&newArray_TnfisUnpPreAmp);
      instance.SetDelete(&delete_TnfisUnpPreAmp);
      instance.SetDeleteArray(&deleteArray_TnfisUnpPreAmp);
      instance.SetDestructor(&destruct_TnfisUnpPreAmp);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TnfisUnpPreAmp*)
   {
      return GenerateInitInstanceLocal((::TnfisUnpPreAmp*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TnfisUnpPreAmp*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TnfisAbsorber(void *p = 0);
   static void *newArray_TnfisAbsorber(Long_t size, void *p);
   static void delete_TnfisAbsorber(void *p);
   static void deleteArray_TnfisAbsorber(void *p);
   static void destruct_TnfisAbsorber(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TnfisAbsorber*)
   {
      ::TnfisAbsorber *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TnfisAbsorber >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TnfisAbsorber", ::TnfisAbsorber::Class_Version(), "TnfisUnpackEvent.h", 102,
                  typeid(::TnfisAbsorber), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TnfisAbsorber::Dictionary, isa_proxy, 4,
                  sizeof(::TnfisAbsorber) );
      instance.SetNew(&new_TnfisAbsorber);
      instance.SetNewArray(&newArray_TnfisAbsorber);
      instance.SetDelete(&delete_TnfisAbsorber);
      instance.SetDeleteArray(&deleteArray_TnfisAbsorber);
      instance.SetDestructor(&destruct_TnfisAbsorber);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TnfisAbsorber*)
   {
      return GenerateInitInstanceLocal((::TnfisAbsorber*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TnfisAbsorber*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TnfisVetoLength(void *p = 0);
   static void *newArray_TnfisVetoLength(Long_t size, void *p);
   static void delete_TnfisVetoLength(void *p);
   static void deleteArray_TnfisVetoLength(void *p);
   static void destruct_TnfisVetoLength(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TnfisVetoLength*)
   {
      ::TnfisVetoLength *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TnfisVetoLength >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TnfisVetoLength", ::TnfisVetoLength::Class_Version(), "TnfisUnpackEvent.h", 124,
                  typeid(::TnfisVetoLength), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TnfisVetoLength::Dictionary, isa_proxy, 4,
                  sizeof(::TnfisVetoLength) );
      instance.SetNew(&new_TnfisVetoLength);
      instance.SetNewArray(&newArray_TnfisVetoLength);
      instance.SetDelete(&delete_TnfisVetoLength);
      instance.SetDeleteArray(&deleteArray_TnfisVetoLength);
      instance.SetDestructor(&destruct_TnfisVetoLength);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TnfisVetoLength*)
   {
      return GenerateInitInstanceLocal((::TnfisVetoLength*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TnfisVetoLength*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TnfisTrigAccClk(void *p = 0);
   static void *newArray_TnfisTrigAccClk(Long_t size, void *p);
   static void delete_TnfisTrigAccClk(void *p);
   static void deleteArray_TnfisTrigAccClk(void *p);
   static void destruct_TnfisTrigAccClk(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TnfisTrigAccClk*)
   {
      ::TnfisTrigAccClk *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TnfisTrigAccClk >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TnfisTrigAccClk", ::TnfisTrigAccClk::Class_Version(), "TnfisUnpackEvent.h", 147,
                  typeid(::TnfisTrigAccClk), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TnfisTrigAccClk::Dictionary, isa_proxy, 4,
                  sizeof(::TnfisTrigAccClk) );
      instance.SetNew(&new_TnfisTrigAccClk);
      instance.SetNewArray(&newArray_TnfisTrigAccClk);
      instance.SetDelete(&delete_TnfisTrigAccClk);
      instance.SetDeleteArray(&deleteArray_TnfisTrigAccClk);
      instance.SetDestructor(&destruct_TnfisTrigAccClk);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TnfisTrigAccClk*)
   {
      return GenerateInitInstanceLocal((::TnfisTrigAccClk*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TnfisTrigAccClk*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TnfisEvtCounter(void *p = 0);
   static void *newArray_TnfisEvtCounter(Long_t size, void *p);
   static void delete_TnfisEvtCounter(void *p);
   static void deleteArray_TnfisEvtCounter(void *p);
   static void destruct_TnfisEvtCounter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TnfisEvtCounter*)
   {
      ::TnfisEvtCounter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TnfisEvtCounter >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TnfisEvtCounter", ::TnfisEvtCounter::Class_Version(), "TnfisUnpackEvent.h", 182,
                  typeid(::TnfisEvtCounter), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TnfisEvtCounter::Dictionary, isa_proxy, 4,
                  sizeof(::TnfisEvtCounter) );
      instance.SetNew(&new_TnfisEvtCounter);
      instance.SetNewArray(&newArray_TnfisEvtCounter);
      instance.SetDelete(&delete_TnfisEvtCounter);
      instance.SetDeleteArray(&deleteArray_TnfisEvtCounter);
      instance.SetDestructor(&destruct_TnfisEvtCounter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TnfisEvtCounter*)
   {
      return GenerateInitInstanceLocal((::TnfisEvtCounter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TnfisEvtCounter*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TnfisCounter(void *p = 0);
   static void *newArray_TnfisCounter(Long_t size, void *p);
   static void delete_TnfisCounter(void *p);
   static void deleteArray_TnfisCounter(void *p);
   static void destruct_TnfisCounter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TnfisCounter*)
   {
      ::TnfisCounter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TnfisCounter >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TnfisCounter", ::TnfisCounter::Class_Version(), "TnfisUnpackEvent.h", 224,
                  typeid(::TnfisCounter), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TnfisCounter::Dictionary, isa_proxy, 4,
                  sizeof(::TnfisCounter) );
      instance.SetNew(&new_TnfisCounter);
      instance.SetNewArray(&newArray_TnfisCounter);
      instance.SetDelete(&delete_TnfisCounter);
      instance.SetDeleteArray(&deleteArray_TnfisCounter);
      instance.SetDestructor(&destruct_TnfisCounter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TnfisCounter*)
   {
      return GenerateInitInstanceLocal((::TnfisCounter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TnfisCounter*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TnfisTimer(void *p = 0);
   static void *newArray_TnfisTimer(Long_t size, void *p);
   static void delete_TnfisTimer(void *p);
   static void deleteArray_TnfisTimer(void *p);
   static void destruct_TnfisTimer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TnfisTimer*)
   {
      ::TnfisTimer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TnfisTimer >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TnfisTimer", ::TnfisTimer::Class_Version(), "TnfisUnpackEvent.h", 259,
                  typeid(::TnfisTimer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TnfisTimer::Dictionary, isa_proxy, 4,
                  sizeof(::TnfisTimer) );
      instance.SetNew(&new_TnfisTimer);
      instance.SetNewArray(&newArray_TnfisTimer);
      instance.SetDelete(&delete_TnfisTimer);
      instance.SetDeleteArray(&deleteArray_TnfisTimer);
      instance.SetDestructor(&destruct_TnfisTimer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TnfisTimer*)
   {
      return GenerateInitInstanceLocal((::TnfisTimer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TnfisTimer*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TnfisParamGlobal(void *p = 0);
   static void *newArray_TnfisParamGlobal(Long_t size, void *p);
   static void delete_TnfisParamGlobal(void *p);
   static void deleteArray_TnfisParamGlobal(void *p);
   static void destruct_TnfisParamGlobal(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TnfisParamGlobal*)
   {
      ::TnfisParamGlobal *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TnfisParamGlobal >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TnfisParamGlobal", ::TnfisParamGlobal::Class_Version(), "TnfisParam.h", 16,
                  typeid(::TnfisParamGlobal), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TnfisParamGlobal::Dictionary, isa_proxy, 4,
                  sizeof(::TnfisParamGlobal) );
      instance.SetNew(&new_TnfisParamGlobal);
      instance.SetNewArray(&newArray_TnfisParamGlobal);
      instance.SetDelete(&delete_TnfisParamGlobal);
      instance.SetDeleteArray(&deleteArray_TnfisParamGlobal);
      instance.SetDestructor(&destruct_TnfisParamGlobal);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TnfisParamGlobal*)
   {
      return GenerateInitInstanceLocal((::TnfisParamGlobal*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TnfisParamGlobal*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TnfisParamQDC(void *p = 0);
   static void *newArray_TnfisParamQDC(Long_t size, void *p);
   static void delete_TnfisParamQDC(void *p);
   static void deleteArray_TnfisParamQDC(void *p);
   static void destruct_TnfisParamQDC(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TnfisParamQDC*)
   {
      ::TnfisParamQDC *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TnfisParamQDC >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TnfisParamQDC", ::TnfisParamQDC::Class_Version(), "TnfisParam.h", 57,
                  typeid(::TnfisParamQDC), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TnfisParamQDC::Dictionary, isa_proxy, 4,
                  sizeof(::TnfisParamQDC) );
      instance.SetNew(&new_TnfisParamQDC);
      instance.SetNewArray(&newArray_TnfisParamQDC);
      instance.SetDelete(&delete_TnfisParamQDC);
      instance.SetDeleteArray(&deleteArray_TnfisParamQDC);
      instance.SetDestructor(&destruct_TnfisParamQDC);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TnfisParamQDC*)
   {
      return GenerateInitInstanceLocal((::TnfisParamQDC*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TnfisParamQDC*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TnfisParamToF(void *p = 0);
   static void *newArray_TnfisParamToF(Long_t size, void *p);
   static void delete_TnfisParamToF(void *p);
   static void deleteArray_TnfisParamToF(void *p);
   static void destruct_TnfisParamToF(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TnfisParamToF*)
   {
      ::TnfisParamToF *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TnfisParamToF >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TnfisParamToF", ::TnfisParamToF::Class_Version(), "TnfisParam.h", 78,
                  typeid(::TnfisParamToF), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TnfisParamToF::Dictionary, isa_proxy, 4,
                  sizeof(::TnfisParamToF) );
      instance.SetNew(&new_TnfisParamToF);
      instance.SetNewArray(&newArray_TnfisParamToF);
      instance.SetDelete(&delete_TnfisParamToF);
      instance.SetDeleteArray(&deleteArray_TnfisParamToF);
      instance.SetDestructor(&destruct_TnfisParamToF);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TnfisParamToF*)
   {
      return GenerateInitInstanceLocal((::TnfisParamToF*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TnfisParamToF*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_nfisHistograms(void *p = 0);
   static void *newArray_nfisHistograms(Long_t size, void *p);
   static void delete_nfisHistograms(void *p);
   static void deleteArray_nfisHistograms(void *p);
   static void destruct_nfisHistograms(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::nfisHistograms*)
   {
      ::nfisHistograms *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::nfisHistograms >(0);
      static ::ROOT::TGenericClassInfo 
         instance("nfisHistograms", ::nfisHistograms::Class_Version(), "nfisHistograms.h", 17,
                  typeid(::nfisHistograms), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::nfisHistograms::Dictionary, isa_proxy, 4,
                  sizeof(::nfisHistograms) );
      instance.SetNew(&new_nfisHistograms);
      instance.SetNewArray(&newArray_nfisHistograms);
      instance.SetDelete(&delete_nfisHistograms);
      instance.SetDeleteArray(&deleteArray_nfisHistograms);
      instance.SetDestructor(&destruct_nfisHistograms);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::nfisHistograms*)
   {
      return GenerateInitInstanceLocal((::nfisHistograms*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::nfisHistograms*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TnfisAnalysis(void *p = 0);
   static void *newArray_TnfisAnalysis(Long_t size, void *p);
   static void delete_TnfisAnalysis(void *p);
   static void deleteArray_TnfisAnalysis(void *p);
   static void destruct_TnfisAnalysis(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TnfisAnalysis*)
   {
      ::TnfisAnalysis *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TnfisAnalysis >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TnfisAnalysis", ::TnfisAnalysis::Class_Version(), "TnfisAnalysis.h", 23,
                  typeid(::TnfisAnalysis), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TnfisAnalysis::Dictionary, isa_proxy, 4,
                  sizeof(::TnfisAnalysis) );
      instance.SetNew(&new_TnfisAnalysis);
      instance.SetNewArray(&newArray_TnfisAnalysis);
      instance.SetDelete(&delete_TnfisAnalysis);
      instance.SetDeleteArray(&deleteArray_TnfisAnalysis);
      instance.SetDestructor(&destruct_TnfisAnalysis);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TnfisAnalysis*)
   {
      return GenerateInitInstanceLocal((::TnfisAnalysis*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TnfisAnalysis*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TnfisMakeUnp(void *p = 0);
   static void *newArray_TnfisMakeUnp(Long_t size, void *p);
   static void delete_TnfisMakeUnp(void *p);
   static void deleteArray_TnfisMakeUnp(void *p);
   static void destruct_TnfisMakeUnp(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TnfisMakeUnp*)
   {
      ::TnfisMakeUnp *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TnfisMakeUnp >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TnfisMakeUnp", ::TnfisMakeUnp::Class_Version(), "TnfisUnpackProc.h", 29,
                  typeid(::TnfisMakeUnp), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TnfisMakeUnp::Dictionary, isa_proxy, 4,
                  sizeof(::TnfisMakeUnp) );
      instance.SetNew(&new_TnfisMakeUnp);
      instance.SetNewArray(&newArray_TnfisMakeUnp);
      instance.SetDelete(&delete_TnfisMakeUnp);
      instance.SetDeleteArray(&deleteArray_TnfisMakeUnp);
      instance.SetDestructor(&destruct_TnfisMakeUnp);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TnfisMakeUnp*)
   {
      return GenerateInitInstanceLocal((::TnfisMakeUnp*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TnfisMakeUnp*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TnfisFCAnalysis(void *p = 0);
   static void *newArray_TnfisFCAnalysis(Long_t size, void *p);
   static void delete_TnfisFCAnalysis(void *p);
   static void deleteArray_TnfisFCAnalysis(void *p);
   static void destruct_TnfisFCAnalysis(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TnfisFCAnalysis*)
   {
      ::TnfisFCAnalysis *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TnfisFCAnalysis >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TnfisFCAnalysis", ::TnfisFCAnalysis::Class_Version(), "TnfisAnalysisFC.h", 23,
                  typeid(::TnfisFCAnalysis), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TnfisFCAnalysis::Dictionary, isa_proxy, 4,
                  sizeof(::TnfisFCAnalysis) );
      instance.SetNew(&new_TnfisFCAnalysis);
      instance.SetNewArray(&newArray_TnfisFCAnalysis);
      instance.SetDelete(&delete_TnfisFCAnalysis);
      instance.SetDeleteArray(&deleteArray_TnfisFCAnalysis);
      instance.SetDestructor(&destruct_TnfisFCAnalysis);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TnfisFCAnalysis*)
   {
      return GenerateInitInstanceLocal((::TnfisFCAnalysis*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TnfisFCAnalysis*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TnfisAnaFCEvent(void *p = 0);
   static void *newArray_TnfisAnaFCEvent(Long_t size, void *p);
   static void delete_TnfisAnaFCEvent(void *p);
   static void deleteArray_TnfisAnaFCEvent(void *p);
   static void destruct_TnfisAnaFCEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TnfisAnaFCEvent*)
   {
      ::TnfisAnaFCEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TnfisAnaFCEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TnfisAnaFCEvent", ::TnfisAnaFCEvent::Class_Version(), "TnfisAnalysisFC.h", 70,
                  typeid(::TnfisAnaFCEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TnfisAnaFCEvent::Dictionary, isa_proxy, 4,
                  sizeof(::TnfisAnaFCEvent) );
      instance.SetNew(&new_TnfisAnaFCEvent);
      instance.SetNewArray(&newArray_TnfisAnaFCEvent);
      instance.SetDelete(&delete_TnfisAnaFCEvent);
      instance.SetDeleteArray(&deleteArray_TnfisAnaFCEvent);
      instance.SetDestructor(&destruct_TnfisAnaFCEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TnfisAnaFCEvent*)
   {
      return GenerateInitInstanceLocal((::TnfisAnaFCEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TnfisAnaFCEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TnfisAnaFC(void *p = 0);
   static void *newArray_TnfisAnaFC(Long_t size, void *p);
   static void delete_TnfisAnaFC(void *p);
   static void deleteArray_TnfisAnaFC(void *p);
   static void destruct_TnfisAnaFC(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TnfisAnaFC*)
   {
      ::TnfisAnaFC *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TnfisAnaFC >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TnfisAnaFC", ::TnfisAnaFC::Class_Version(), "TnfisAnalysisFC.h", 91,
                  typeid(::TnfisAnaFC), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TnfisAnaFC::Dictionary, isa_proxy, 4,
                  sizeof(::TnfisAnaFC) );
      instance.SetNew(&new_TnfisAnaFC);
      instance.SetNewArray(&newArray_TnfisAnaFC);
      instance.SetDelete(&delete_TnfisAnaFC);
      instance.SetDeleteArray(&deleteArray_TnfisAnaFC);
      instance.SetDestructor(&destruct_TnfisAnaFC);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TnfisAnaFC*)
   {
      return GenerateInitInstanceLocal((::TnfisAnaFC*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TnfisAnaFC*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TnfisAnaPreAmp(void *p = 0);
   static void *newArray_TnfisAnaPreAmp(Long_t size, void *p);
   static void delete_TnfisAnaPreAmp(void *p);
   static void deleteArray_TnfisAnaPreAmp(void *p);
   static void destruct_TnfisAnaPreAmp(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TnfisAnaPreAmp*)
   {
      ::TnfisAnaPreAmp *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TnfisAnaPreAmp >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TnfisAnaPreAmp", ::TnfisAnaPreAmp::Class_Version(), "TnfisAnalysisFC.h", 111,
                  typeid(::TnfisAnaPreAmp), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TnfisAnaPreAmp::Dictionary, isa_proxy, 4,
                  sizeof(::TnfisAnaPreAmp) );
      instance.SetNew(&new_TnfisAnaPreAmp);
      instance.SetNewArray(&newArray_TnfisAnaPreAmp);
      instance.SetDelete(&delete_TnfisAnaPreAmp);
      instance.SetDeleteArray(&deleteArray_TnfisAnaPreAmp);
      instance.SetDestructor(&destruct_TnfisAnaPreAmp);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TnfisAnaPreAmp*)
   {
      return GenerateInitInstanceLocal((::TnfisAnaPreAmp*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TnfisAnaPreAmp*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TnfisUnpackEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TnfisUnpackEvent::Class_Name()
{
   return "TnfisUnpackEvent";
}

//______________________________________________________________________________
const char *TnfisUnpackEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisUnpackEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TnfisUnpackEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisUnpackEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TnfisUnpackEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisUnpackEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TnfisUnpackEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisUnpackEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TnfisUnpFC::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TnfisUnpFC::Class_Name()
{
   return "TnfisUnpFC";
}

//______________________________________________________________________________
const char *TnfisUnpFC::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisUnpFC*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TnfisUnpFC::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisUnpFC*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TnfisUnpFC::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisUnpFC*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TnfisUnpFC::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisUnpFC*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TnfisUnpPreAmp::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TnfisUnpPreAmp::Class_Name()
{
   return "TnfisUnpPreAmp";
}

//______________________________________________________________________________
const char *TnfisUnpPreAmp::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisUnpPreAmp*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TnfisUnpPreAmp::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisUnpPreAmp*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TnfisUnpPreAmp::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisUnpPreAmp*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TnfisUnpPreAmp::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisUnpPreAmp*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TnfisAbsorber::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TnfisAbsorber::Class_Name()
{
   return "TnfisAbsorber";
}

//______________________________________________________________________________
const char *TnfisAbsorber::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisAbsorber*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TnfisAbsorber::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisAbsorber*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TnfisAbsorber::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisAbsorber*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TnfisAbsorber::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisAbsorber*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TnfisVetoLength::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TnfisVetoLength::Class_Name()
{
   return "TnfisVetoLength";
}

//______________________________________________________________________________
const char *TnfisVetoLength::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisVetoLength*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TnfisVetoLength::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisVetoLength*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TnfisVetoLength::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisVetoLength*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TnfisVetoLength::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisVetoLength*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TnfisTrigAccClk::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TnfisTrigAccClk::Class_Name()
{
   return "TnfisTrigAccClk";
}

//______________________________________________________________________________
const char *TnfisTrigAccClk::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisTrigAccClk*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TnfisTrigAccClk::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisTrigAccClk*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TnfisTrigAccClk::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisTrigAccClk*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TnfisTrigAccClk::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisTrigAccClk*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TnfisEvtCounter::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TnfisEvtCounter::Class_Name()
{
   return "TnfisEvtCounter";
}

//______________________________________________________________________________
const char *TnfisEvtCounter::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisEvtCounter*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TnfisEvtCounter::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisEvtCounter*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TnfisEvtCounter::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisEvtCounter*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TnfisEvtCounter::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisEvtCounter*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TnfisCounter::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TnfisCounter::Class_Name()
{
   return "TnfisCounter";
}

//______________________________________________________________________________
const char *TnfisCounter::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisCounter*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TnfisCounter::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisCounter*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TnfisCounter::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisCounter*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TnfisCounter::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisCounter*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TnfisTimer::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TnfisTimer::Class_Name()
{
   return "TnfisTimer";
}

//______________________________________________________________________________
const char *TnfisTimer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisTimer*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TnfisTimer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisTimer*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TnfisTimer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisTimer*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TnfisTimer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisTimer*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TnfisParamGlobal::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TnfisParamGlobal::Class_Name()
{
   return "TnfisParamGlobal";
}

//______________________________________________________________________________
const char *TnfisParamGlobal::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisParamGlobal*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TnfisParamGlobal::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisParamGlobal*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TnfisParamGlobal::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisParamGlobal*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TnfisParamGlobal::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisParamGlobal*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TnfisParamQDC::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TnfisParamQDC::Class_Name()
{
   return "TnfisParamQDC";
}

//______________________________________________________________________________
const char *TnfisParamQDC::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisParamQDC*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TnfisParamQDC::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisParamQDC*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TnfisParamQDC::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisParamQDC*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TnfisParamQDC::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisParamQDC*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TnfisParamToF::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TnfisParamToF::Class_Name()
{
   return "TnfisParamToF";
}

//______________________________________________________________________________
const char *TnfisParamToF::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisParamToF*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TnfisParamToF::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisParamToF*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TnfisParamToF::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisParamToF*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TnfisParamToF::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisParamToF*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr nfisHistograms::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *nfisHistograms::Class_Name()
{
   return "nfisHistograms";
}

//______________________________________________________________________________
const char *nfisHistograms::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::nfisHistograms*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int nfisHistograms::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::nfisHistograms*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *nfisHistograms::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::nfisHistograms*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *nfisHistograms::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::nfisHistograms*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TnfisAnalysis::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TnfisAnalysis::Class_Name()
{
   return "TnfisAnalysis";
}

//______________________________________________________________________________
const char *TnfisAnalysis::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisAnalysis*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TnfisAnalysis::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisAnalysis*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TnfisAnalysis::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisAnalysis*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TnfisAnalysis::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisAnalysis*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TnfisMakeUnp::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TnfisMakeUnp::Class_Name()
{
   return "TnfisMakeUnp";
}

//______________________________________________________________________________
const char *TnfisMakeUnp::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisMakeUnp*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TnfisMakeUnp::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisMakeUnp*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TnfisMakeUnp::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisMakeUnp*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TnfisMakeUnp::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisMakeUnp*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TnfisFCAnalysis::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TnfisFCAnalysis::Class_Name()
{
   return "TnfisFCAnalysis";
}

//______________________________________________________________________________
const char *TnfisFCAnalysis::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisFCAnalysis*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TnfisFCAnalysis::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisFCAnalysis*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TnfisFCAnalysis::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisFCAnalysis*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TnfisFCAnalysis::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisFCAnalysis*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TnfisAnaFCEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TnfisAnaFCEvent::Class_Name()
{
   return "TnfisAnaFCEvent";
}

//______________________________________________________________________________
const char *TnfisAnaFCEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisAnaFCEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TnfisAnaFCEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisAnaFCEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TnfisAnaFCEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisAnaFCEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TnfisAnaFCEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisAnaFCEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TnfisAnaFC::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TnfisAnaFC::Class_Name()
{
   return "TnfisAnaFC";
}

//______________________________________________________________________________
const char *TnfisAnaFC::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisAnaFC*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TnfisAnaFC::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisAnaFC*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TnfisAnaFC::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisAnaFC*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TnfisAnaFC::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisAnaFC*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TnfisAnaPreAmp::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TnfisAnaPreAmp::Class_Name()
{
   return "TnfisAnaPreAmp";
}

//______________________________________________________________________________
const char *TnfisAnaPreAmp::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisAnaPreAmp*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TnfisAnaPreAmp::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TnfisAnaPreAmp*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TnfisAnaPreAmp::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisAnaPreAmp*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TnfisAnaPreAmp::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TnfisAnaPreAmp*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TnfisUnpackEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class TnfisUnpackEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TnfisUnpackEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(TnfisUnpackEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TnfisUnpackEvent(void *p) {
      return  p ? new(p) ::TnfisUnpackEvent : new ::TnfisUnpackEvent;
   }
   static void *newArray_TnfisUnpackEvent(Long_t nElements, void *p) {
      return p ? new(p) ::TnfisUnpackEvent[nElements] : new ::TnfisUnpackEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_TnfisUnpackEvent(void *p) {
      delete ((::TnfisUnpackEvent*)p);
   }
   static void deleteArray_TnfisUnpackEvent(void *p) {
      delete [] ((::TnfisUnpackEvent*)p);
   }
   static void destruct_TnfisUnpackEvent(void *p) {
      typedef ::TnfisUnpackEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TnfisUnpackEvent

//______________________________________________________________________________
void TnfisUnpFC::Streamer(TBuffer &R__b)
{
   // Stream an object of class TnfisUnpFC.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TnfisUnpFC::Class(),this);
   } else {
      R__b.WriteClassBuffer(TnfisUnpFC::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TnfisUnpFC(void *p) {
      return  p ? new(p) ::TnfisUnpFC : new ::TnfisUnpFC;
   }
   static void *newArray_TnfisUnpFC(Long_t nElements, void *p) {
      return p ? new(p) ::TnfisUnpFC[nElements] : new ::TnfisUnpFC[nElements];
   }
   // Wrapper around operator delete
   static void delete_TnfisUnpFC(void *p) {
      delete ((::TnfisUnpFC*)p);
   }
   static void deleteArray_TnfisUnpFC(void *p) {
      delete [] ((::TnfisUnpFC*)p);
   }
   static void destruct_TnfisUnpFC(void *p) {
      typedef ::TnfisUnpFC current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TnfisUnpFC

//______________________________________________________________________________
void TnfisUnpPreAmp::Streamer(TBuffer &R__b)
{
   // Stream an object of class TnfisUnpPreAmp.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TnfisUnpPreAmp::Class(),this);
   } else {
      R__b.WriteClassBuffer(TnfisUnpPreAmp::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TnfisUnpPreAmp(void *p) {
      return  p ? new(p) ::TnfisUnpPreAmp : new ::TnfisUnpPreAmp;
   }
   static void *newArray_TnfisUnpPreAmp(Long_t nElements, void *p) {
      return p ? new(p) ::TnfisUnpPreAmp[nElements] : new ::TnfisUnpPreAmp[nElements];
   }
   // Wrapper around operator delete
   static void delete_TnfisUnpPreAmp(void *p) {
      delete ((::TnfisUnpPreAmp*)p);
   }
   static void deleteArray_TnfisUnpPreAmp(void *p) {
      delete [] ((::TnfisUnpPreAmp*)p);
   }
   static void destruct_TnfisUnpPreAmp(void *p) {
      typedef ::TnfisUnpPreAmp current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TnfisUnpPreAmp

//______________________________________________________________________________
void TnfisAbsorber::Streamer(TBuffer &R__b)
{
   // Stream an object of class TnfisAbsorber.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TnfisAbsorber::Class(),this);
   } else {
      R__b.WriteClassBuffer(TnfisAbsorber::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TnfisAbsorber(void *p) {
      return  p ? new(p) ::TnfisAbsorber : new ::TnfisAbsorber;
   }
   static void *newArray_TnfisAbsorber(Long_t nElements, void *p) {
      return p ? new(p) ::TnfisAbsorber[nElements] : new ::TnfisAbsorber[nElements];
   }
   // Wrapper around operator delete
   static void delete_TnfisAbsorber(void *p) {
      delete ((::TnfisAbsorber*)p);
   }
   static void deleteArray_TnfisAbsorber(void *p) {
      delete [] ((::TnfisAbsorber*)p);
   }
   static void destruct_TnfisAbsorber(void *p) {
      typedef ::TnfisAbsorber current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TnfisAbsorber

//______________________________________________________________________________
void TnfisVetoLength::Streamer(TBuffer &R__b)
{
   // Stream an object of class TnfisVetoLength.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TnfisVetoLength::Class(),this);
   } else {
      R__b.WriteClassBuffer(TnfisVetoLength::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TnfisVetoLength(void *p) {
      return  p ? new(p) ::TnfisVetoLength : new ::TnfisVetoLength;
   }
   static void *newArray_TnfisVetoLength(Long_t nElements, void *p) {
      return p ? new(p) ::TnfisVetoLength[nElements] : new ::TnfisVetoLength[nElements];
   }
   // Wrapper around operator delete
   static void delete_TnfisVetoLength(void *p) {
      delete ((::TnfisVetoLength*)p);
   }
   static void deleteArray_TnfisVetoLength(void *p) {
      delete [] ((::TnfisVetoLength*)p);
   }
   static void destruct_TnfisVetoLength(void *p) {
      typedef ::TnfisVetoLength current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TnfisVetoLength

//______________________________________________________________________________
void TnfisTrigAccClk::Streamer(TBuffer &R__b)
{
   // Stream an object of class TnfisTrigAccClk.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TnfisTrigAccClk::Class(),this);
   } else {
      R__b.WriteClassBuffer(TnfisTrigAccClk::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TnfisTrigAccClk(void *p) {
      return  p ? new(p) ::TnfisTrigAccClk : new ::TnfisTrigAccClk;
   }
   static void *newArray_TnfisTrigAccClk(Long_t nElements, void *p) {
      return p ? new(p) ::TnfisTrigAccClk[nElements] : new ::TnfisTrigAccClk[nElements];
   }
   // Wrapper around operator delete
   static void delete_TnfisTrigAccClk(void *p) {
      delete ((::TnfisTrigAccClk*)p);
   }
   static void deleteArray_TnfisTrigAccClk(void *p) {
      delete [] ((::TnfisTrigAccClk*)p);
   }
   static void destruct_TnfisTrigAccClk(void *p) {
      typedef ::TnfisTrigAccClk current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TnfisTrigAccClk

//______________________________________________________________________________
void TnfisEvtCounter::Streamer(TBuffer &R__b)
{
   // Stream an object of class TnfisEvtCounter.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TnfisEvtCounter::Class(),this);
   } else {
      R__b.WriteClassBuffer(TnfisEvtCounter::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TnfisEvtCounter(void *p) {
      return  p ? new(p) ::TnfisEvtCounter : new ::TnfisEvtCounter;
   }
   static void *newArray_TnfisEvtCounter(Long_t nElements, void *p) {
      return p ? new(p) ::TnfisEvtCounter[nElements] : new ::TnfisEvtCounter[nElements];
   }
   // Wrapper around operator delete
   static void delete_TnfisEvtCounter(void *p) {
      delete ((::TnfisEvtCounter*)p);
   }
   static void deleteArray_TnfisEvtCounter(void *p) {
      delete [] ((::TnfisEvtCounter*)p);
   }
   static void destruct_TnfisEvtCounter(void *p) {
      typedef ::TnfisEvtCounter current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TnfisEvtCounter

//______________________________________________________________________________
void TnfisCounter::Streamer(TBuffer &R__b)
{
   // Stream an object of class TnfisCounter.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TnfisCounter::Class(),this);
   } else {
      R__b.WriteClassBuffer(TnfisCounter::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TnfisCounter(void *p) {
      return  p ? new(p) ::TnfisCounter : new ::TnfisCounter;
   }
   static void *newArray_TnfisCounter(Long_t nElements, void *p) {
      return p ? new(p) ::TnfisCounter[nElements] : new ::TnfisCounter[nElements];
   }
   // Wrapper around operator delete
   static void delete_TnfisCounter(void *p) {
      delete ((::TnfisCounter*)p);
   }
   static void deleteArray_TnfisCounter(void *p) {
      delete [] ((::TnfisCounter*)p);
   }
   static void destruct_TnfisCounter(void *p) {
      typedef ::TnfisCounter current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TnfisCounter

//______________________________________________________________________________
void TnfisTimer::Streamer(TBuffer &R__b)
{
   // Stream an object of class TnfisTimer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TnfisTimer::Class(),this);
   } else {
      R__b.WriteClassBuffer(TnfisTimer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TnfisTimer(void *p) {
      return  p ? new(p) ::TnfisTimer : new ::TnfisTimer;
   }
   static void *newArray_TnfisTimer(Long_t nElements, void *p) {
      return p ? new(p) ::TnfisTimer[nElements] : new ::TnfisTimer[nElements];
   }
   // Wrapper around operator delete
   static void delete_TnfisTimer(void *p) {
      delete ((::TnfisTimer*)p);
   }
   static void deleteArray_TnfisTimer(void *p) {
      delete [] ((::TnfisTimer*)p);
   }
   static void destruct_TnfisTimer(void *p) {
      typedef ::TnfisTimer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TnfisTimer

//______________________________________________________________________________
void TnfisParamGlobal::Streamer(TBuffer &R__b)
{
   // Stream an object of class TnfisParamGlobal.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TnfisParamGlobal::Class(),this);
   } else {
      R__b.WriteClassBuffer(TnfisParamGlobal::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TnfisParamGlobal(void *p) {
      return  p ? new(p) ::TnfisParamGlobal : new ::TnfisParamGlobal;
   }
   static void *newArray_TnfisParamGlobal(Long_t nElements, void *p) {
      return p ? new(p) ::TnfisParamGlobal[nElements] : new ::TnfisParamGlobal[nElements];
   }
   // Wrapper around operator delete
   static void delete_TnfisParamGlobal(void *p) {
      delete ((::TnfisParamGlobal*)p);
   }
   static void deleteArray_TnfisParamGlobal(void *p) {
      delete [] ((::TnfisParamGlobal*)p);
   }
   static void destruct_TnfisParamGlobal(void *p) {
      typedef ::TnfisParamGlobal current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TnfisParamGlobal

//______________________________________________________________________________
void TnfisParamQDC::Streamer(TBuffer &R__b)
{
   // Stream an object of class TnfisParamQDC.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TnfisParamQDC::Class(),this);
   } else {
      R__b.WriteClassBuffer(TnfisParamQDC::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TnfisParamQDC(void *p) {
      return  p ? new(p) ::TnfisParamQDC : new ::TnfisParamQDC;
   }
   static void *newArray_TnfisParamQDC(Long_t nElements, void *p) {
      return p ? new(p) ::TnfisParamQDC[nElements] : new ::TnfisParamQDC[nElements];
   }
   // Wrapper around operator delete
   static void delete_TnfisParamQDC(void *p) {
      delete ((::TnfisParamQDC*)p);
   }
   static void deleteArray_TnfisParamQDC(void *p) {
      delete [] ((::TnfisParamQDC*)p);
   }
   static void destruct_TnfisParamQDC(void *p) {
      typedef ::TnfisParamQDC current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TnfisParamQDC

//______________________________________________________________________________
void TnfisParamToF::Streamer(TBuffer &R__b)
{
   // Stream an object of class TnfisParamToF.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TnfisParamToF::Class(),this);
   } else {
      R__b.WriteClassBuffer(TnfisParamToF::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TnfisParamToF(void *p) {
      return  p ? new(p) ::TnfisParamToF : new ::TnfisParamToF;
   }
   static void *newArray_TnfisParamToF(Long_t nElements, void *p) {
      return p ? new(p) ::TnfisParamToF[nElements] : new ::TnfisParamToF[nElements];
   }
   // Wrapper around operator delete
   static void delete_TnfisParamToF(void *p) {
      delete ((::TnfisParamToF*)p);
   }
   static void deleteArray_TnfisParamToF(void *p) {
      delete [] ((::TnfisParamToF*)p);
   }
   static void destruct_TnfisParamToF(void *p) {
      typedef ::TnfisParamToF current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TnfisParamToF

//______________________________________________________________________________
void nfisHistograms::Streamer(TBuffer &R__b)
{
   // Stream an object of class nfisHistograms.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(nfisHistograms::Class(),this);
   } else {
      R__b.WriteClassBuffer(nfisHistograms::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_nfisHistograms(void *p) {
      return  p ? new(p) ::nfisHistograms : new ::nfisHistograms;
   }
   static void *newArray_nfisHistograms(Long_t nElements, void *p) {
      return p ? new(p) ::nfisHistograms[nElements] : new ::nfisHistograms[nElements];
   }
   // Wrapper around operator delete
   static void delete_nfisHistograms(void *p) {
      delete ((::nfisHistograms*)p);
   }
   static void deleteArray_nfisHistograms(void *p) {
      delete [] ((::nfisHistograms*)p);
   }
   static void destruct_nfisHistograms(void *p) {
      typedef ::nfisHistograms current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::nfisHistograms

//______________________________________________________________________________
void TnfisAnalysis::Streamer(TBuffer &R__b)
{
   // Stream an object of class TnfisAnalysis.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TnfisAnalysis::Class(),this);
   } else {
      R__b.WriteClassBuffer(TnfisAnalysis::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TnfisAnalysis(void *p) {
      return  p ? new(p) ::TnfisAnalysis : new ::TnfisAnalysis;
   }
   static void *newArray_TnfisAnalysis(Long_t nElements, void *p) {
      return p ? new(p) ::TnfisAnalysis[nElements] : new ::TnfisAnalysis[nElements];
   }
   // Wrapper around operator delete
   static void delete_TnfisAnalysis(void *p) {
      delete ((::TnfisAnalysis*)p);
   }
   static void deleteArray_TnfisAnalysis(void *p) {
      delete [] ((::TnfisAnalysis*)p);
   }
   static void destruct_TnfisAnalysis(void *p) {
      typedef ::TnfisAnalysis current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TnfisAnalysis

//______________________________________________________________________________
void TnfisMakeUnp::Streamer(TBuffer &R__b)
{
   // Stream an object of class TnfisMakeUnp.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TnfisMakeUnp::Class(),this);
   } else {
      R__b.WriteClassBuffer(TnfisMakeUnp::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TnfisMakeUnp(void *p) {
      return  p ? new(p) ::TnfisMakeUnp : new ::TnfisMakeUnp;
   }
   static void *newArray_TnfisMakeUnp(Long_t nElements, void *p) {
      return p ? new(p) ::TnfisMakeUnp[nElements] : new ::TnfisMakeUnp[nElements];
   }
   // Wrapper around operator delete
   static void delete_TnfisMakeUnp(void *p) {
      delete ((::TnfisMakeUnp*)p);
   }
   static void deleteArray_TnfisMakeUnp(void *p) {
      delete [] ((::TnfisMakeUnp*)p);
   }
   static void destruct_TnfisMakeUnp(void *p) {
      typedef ::TnfisMakeUnp current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TnfisMakeUnp

//______________________________________________________________________________
void TnfisFCAnalysis::Streamer(TBuffer &R__b)
{
   // Stream an object of class TnfisFCAnalysis.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TnfisFCAnalysis::Class(),this);
   } else {
      R__b.WriteClassBuffer(TnfisFCAnalysis::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TnfisFCAnalysis(void *p) {
      return  p ? new(p) ::TnfisFCAnalysis : new ::TnfisFCAnalysis;
   }
   static void *newArray_TnfisFCAnalysis(Long_t nElements, void *p) {
      return p ? new(p) ::TnfisFCAnalysis[nElements] : new ::TnfisFCAnalysis[nElements];
   }
   // Wrapper around operator delete
   static void delete_TnfisFCAnalysis(void *p) {
      delete ((::TnfisFCAnalysis*)p);
   }
   static void deleteArray_TnfisFCAnalysis(void *p) {
      delete [] ((::TnfisFCAnalysis*)p);
   }
   static void destruct_TnfisFCAnalysis(void *p) {
      typedef ::TnfisFCAnalysis current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TnfisFCAnalysis

//______________________________________________________________________________
void TnfisAnaFCEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class TnfisAnaFCEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TnfisAnaFCEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(TnfisAnaFCEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TnfisAnaFCEvent(void *p) {
      return  p ? new(p) ::TnfisAnaFCEvent : new ::TnfisAnaFCEvent;
   }
   static void *newArray_TnfisAnaFCEvent(Long_t nElements, void *p) {
      return p ? new(p) ::TnfisAnaFCEvent[nElements] : new ::TnfisAnaFCEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_TnfisAnaFCEvent(void *p) {
      delete ((::TnfisAnaFCEvent*)p);
   }
   static void deleteArray_TnfisAnaFCEvent(void *p) {
      delete [] ((::TnfisAnaFCEvent*)p);
   }
   static void destruct_TnfisAnaFCEvent(void *p) {
      typedef ::TnfisAnaFCEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TnfisAnaFCEvent

//______________________________________________________________________________
void TnfisAnaFC::Streamer(TBuffer &R__b)
{
   // Stream an object of class TnfisAnaFC.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TnfisAnaFC::Class(),this);
   } else {
      R__b.WriteClassBuffer(TnfisAnaFC::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TnfisAnaFC(void *p) {
      return  p ? new(p) ::TnfisAnaFC : new ::TnfisAnaFC;
   }
   static void *newArray_TnfisAnaFC(Long_t nElements, void *p) {
      return p ? new(p) ::TnfisAnaFC[nElements] : new ::TnfisAnaFC[nElements];
   }
   // Wrapper around operator delete
   static void delete_TnfisAnaFC(void *p) {
      delete ((::TnfisAnaFC*)p);
   }
   static void deleteArray_TnfisAnaFC(void *p) {
      delete [] ((::TnfisAnaFC*)p);
   }
   static void destruct_TnfisAnaFC(void *p) {
      typedef ::TnfisAnaFC current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TnfisAnaFC

//______________________________________________________________________________
void TnfisAnaPreAmp::Streamer(TBuffer &R__b)
{
   // Stream an object of class TnfisAnaPreAmp.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TnfisAnaPreAmp::Class(),this);
   } else {
      R__b.WriteClassBuffer(TnfisAnaPreAmp::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TnfisAnaPreAmp(void *p) {
      return  p ? new(p) ::TnfisAnaPreAmp : new ::TnfisAnaPreAmp;
   }
   static void *newArray_TnfisAnaPreAmp(Long_t nElements, void *p) {
      return p ? new(p) ::TnfisAnaPreAmp[nElements] : new ::TnfisAnaPreAmp[nElements];
   }
   // Wrapper around operator delete
   static void delete_TnfisAnaPreAmp(void *p) {
      delete ((::TnfisAnaPreAmp*)p);
   }
   static void deleteArray_TnfisAnaPreAmp(void *p) {
      delete [] ((::TnfisAnaPreAmp*)p);
   }
   static void destruct_TnfisAnaPreAmp(void *p) {
      typedef ::TnfisAnaPreAmp current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TnfisAnaPreAmp

namespace ROOT {
   static TClass *vectorlEstringgR_Dictionary();
   static void vectorlEstringgR_TClassManip(TClass*);
   static void *new_vectorlEstringgR(void *p = 0);
   static void *newArray_vectorlEstringgR(Long_t size, void *p);
   static void delete_vectorlEstringgR(void *p);
   static void deleteArray_vectorlEstringgR(void *p);
   static void destruct_vectorlEstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<string>*)
   {
      vector<string> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<string>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<string>", -2, "vector", 214,
                  typeid(vector<string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEstringgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<string>) );
      instance.SetNew(&new_vectorlEstringgR);
      instance.SetNewArray(&newArray_vectorlEstringgR);
      instance.SetDelete(&delete_vectorlEstringgR);
      instance.SetDeleteArray(&deleteArray_vectorlEstringgR);
      instance.SetDestructor(&destruct_vectorlEstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<string> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<string>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<string>*)0x0)->GetClass();
      vectorlEstringgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEstringgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string> : new vector<string>;
   }
   static void *newArray_vectorlEstringgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string>[nElements] : new vector<string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEstringgR(void *p) {
      delete ((vector<string>*)p);
   }
   static void deleteArray_vectorlEstringgR(void *p) {
      delete [] ((vector<string>*)p);
   }
   static void destruct_vectorlEstringgR(void *p) {
      typedef vector<string> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<string>

namespace ROOT {
   static TClass *vectorlElonggR_Dictionary();
   static void vectorlElonggR_TClassManip(TClass*);
   static void *new_vectorlElonggR(void *p = 0);
   static void *newArray_vectorlElonggR(Long_t size, void *p);
   static void delete_vectorlElonggR(void *p);
   static void deleteArray_vectorlElonggR(void *p);
   static void destruct_vectorlElonggR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<long>*)
   {
      vector<long> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<long>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<long>", -2, "vector", 214,
                  typeid(vector<long>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlElonggR_Dictionary, isa_proxy, 4,
                  sizeof(vector<long>) );
      instance.SetNew(&new_vectorlElonggR);
      instance.SetNewArray(&newArray_vectorlElonggR);
      instance.SetDelete(&delete_vectorlElonggR);
      instance.SetDeleteArray(&deleteArray_vectorlElonggR);
      instance.SetDestructor(&destruct_vectorlElonggR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<long> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<long>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlElonggR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<long>*)0x0)->GetClass();
      vectorlElonggR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlElonggR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlElonggR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<long> : new vector<long>;
   }
   static void *newArray_vectorlElonggR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<long>[nElements] : new vector<long>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlElonggR(void *p) {
      delete ((vector<long>*)p);
   }
   static void deleteArray_vectorlElonggR(void *p) {
      delete [] ((vector<long>*)p);
   }
   static void destruct_vectorlElonggR(void *p) {
      typedef vector<long> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<long>

namespace {
  void TriggerDictionaryInitialization_libGo4UserAnalysis_Impl() {
    static const char* headers[] = {
"./TnfisUnpackEvent.h",
"./nfisGlobals.h",
"./TnfisParam.h",
"./TnfisAnalysisFC.h",
"./TnfisUnpackProc.h",
"./nfisHistograms.h",
"./TnfisAnalysis.h",
0
    };
    static const char* includePaths[] = {
"/opt/Go4/go4-5.3.2/include",
"/opt/Go4/go4-5.3.2",
"/opt/root/root-6.18.04/include",
"/gpfs/home/hoffma93/Programme/Go4nfis/offline/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libGo4UserAnalysis dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$./TnfisUnpackEvent.h")))  TnfisUnpackEvent;
class __attribute__((annotate("$clingAutoload$./TnfisUnpackEvent.h")))  TnfisUnpFC;
class __attribute__((annotate("$clingAutoload$./TnfisUnpackEvent.h")))  TnfisUnpPreAmp;
class __attribute__((annotate("$clingAutoload$./TnfisUnpackEvent.h")))  TnfisAbsorber;
class __attribute__((annotate("$clingAutoload$./TnfisUnpackEvent.h")))  TnfisVetoLength;
class __attribute__((annotate("$clingAutoload$./TnfisUnpackEvent.h")))  TnfisTrigAccClk;
class __attribute__((annotate("$clingAutoload$./TnfisUnpackEvent.h")))  TnfisEvtCounter;
class __attribute__((annotate("$clingAutoload$./TnfisUnpackEvent.h")))  TnfisCounter;
class __attribute__((annotate("$clingAutoload$./TnfisUnpackEvent.h")))  TnfisTimer;
class __attribute__((annotate("$clingAutoload$./TnfisParam.h")))  TnfisParamGlobal;
class __attribute__((annotate("$clingAutoload$./TnfisParam.h")))  TnfisParamQDC;
class __attribute__((annotate("$clingAutoload$./TnfisParam.h")))  TnfisParamToF;
class __attribute__((annotate("$clingAutoload$nfisHistograms.h")))  __attribute__((annotate("$clingAutoload$./TnfisAnalysisFC.h")))  nfisHistograms;
class __attribute__((annotate("$clingAutoload$TnfisAnalysis.h")))  __attribute__((annotate("$clingAutoload$./TnfisAnalysisFC.h")))  TnfisAnalysis;
class __attribute__((annotate("$clingAutoload$TnfisUnpackProc.h")))  __attribute__((annotate("$clingAutoload$./TnfisAnalysisFC.h")))  TnfisMakeUnp;
class __attribute__((annotate("$clingAutoload$./TnfisAnalysisFC.h")))  TnfisFCAnalysis;
class __attribute__((annotate("$clingAutoload$./TnfisAnalysisFC.h")))  TnfisAnaFCEvent;
class __attribute__((annotate("$clingAutoload$./TnfisAnalysisFC.h")))  TnfisAnaFC;
class __attribute__((annotate("$clingAutoload$./TnfisAnalysisFC.h")))  TnfisAnaPreAmp;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libGo4UserAnalysis dictionary payload"

#ifndef Linux
  #define Linux 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "./TnfisUnpackEvent.h"
#include "./nfisGlobals.h"
#include "./TnfisParam.h"
#include "./TnfisAnalysisFC.h"
#include "./TnfisUnpackProc.h"
#include "./nfisHistograms.h"
#include "./TnfisAnalysis.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TnfisAbsorber", payloadCode, "@",
"TnfisAnaFC", payloadCode, "@",
"TnfisAnaFCEvent", payloadCode, "@",
"TnfisAnaPreAmp", payloadCode, "@",
"TnfisAnalysis", payloadCode, "@",
"TnfisCounter", payloadCode, "@",
"TnfisEvtCounter", payloadCode, "@",
"TnfisFCAnalysis", payloadCode, "@",
"TnfisMakeUnp", payloadCode, "@",
"TnfisParamGlobal", payloadCode, "@",
"TnfisParamQDC", payloadCode, "@",
"TnfisParamToF", payloadCode, "@",
"TnfisTimer", payloadCode, "@",
"TnfisTrigAccClk", payloadCode, "@",
"TnfisUnpFC", payloadCode, "@",
"TnfisUnpPreAmp", payloadCode, "@",
"TnfisUnpackEvent", payloadCode, "@",
"TnfisVetoLength", payloadCode, "@",
"nfisHistograms", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libGo4UserAnalysis",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libGo4UserAnalysis_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libGo4UserAnalysis_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libGo4UserAnalysis() {
  TriggerDictionaryInitialization_libGo4UserAnalysis_Impl();
}
