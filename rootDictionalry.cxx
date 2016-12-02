// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME objdIrootDictionalry

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "include/eclCrystalDB.h"
#include "include/fileFuncs.h"

// Header files passed via #pragma extra_include

namespace fileFuncs {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *fileFuncs_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("fileFuncs", 0 /*version*/, "include/fileFuncs.h", 7,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &fileFuncs_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *fileFuncs_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace ROOT {
   static TClass *eclCrystalDB_Dictionary();
   static void eclCrystalDB_TClassManip(TClass*);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::eclCrystalDB*)
   {
      ::eclCrystalDB *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::eclCrystalDB));
      static ::ROOT::TGenericClassInfo 
         instance("eclCrystalDB", "include/eclCrystalDB.h", 9,
                  typeid(::eclCrystalDB), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &eclCrystalDB_Dictionary, isa_proxy, 0,
                  sizeof(::eclCrystalDB) );
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::eclCrystalDB*)
   {
      return GenerateInitInstanceLocal((::eclCrystalDB*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::eclCrystalDB*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *eclCrystalDB_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::eclCrystalDB*)0x0)->GetClass();
      eclCrystalDB_TClassManip(theClass);
   return theClass;
   }

   static void eclCrystalDB_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
} // end of namespace ROOT for class ::eclCrystalDB

namespace {
  void TriggerDictionaryInitialization_rootDictionalry_Impl() {
    static const char* headers[] = {
"include/eclCrystalDB.h",
"include/fileFuncs.h",
0
    };
    static const char* includePaths[] = {
"/home/hershen/basf2/externals/v01-02-02/Linux_x86_64/opt/root/include",
"/media/sf_PhD/Root/commonRootFunctions/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "rootDictionalry dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$include/eclCrystalDB.h")))  eclCrystalDB;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "rootDictionalry dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "include/eclCrystalDB.h"
#include "include/fileFuncs.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"eclCrystalDB", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("rootDictionalry",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_rootDictionalry_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_rootDictionalry_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_rootDictionalry() {
  TriggerDictionaryInitialization_rootDictionalry_Impl();
}
