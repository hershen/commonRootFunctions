CC=g++
DEBUG=-g -Wall -Wextra -Wpedantic
CXX11=-std=c++1y
CPPCHECK=cppcheck --enable=warning,style,performance --inconclusive --language=c++ -q --std=c++11

CCT=$(CC) $(DEBUG) $(CXX11)

SRC_DIR=src
OBJ_DIR=obj

COMMON_ROOT_INC_DIR=/home/hershen/PhD/Root/commonRootFunctions/include
COMMON_ROOT_MATHFUNCS=/home/hershen/PhD/Root/commonRootFunctions/obj/mathFuncs.o
COMMON_ROOT_FILEFUNCS=/home/hershen/PhD/Root/commonRootFunctions/obj/fileFuncs.o
COMMON_ROOT_VT=/home/hershen/PhD/Root/commonRootFunctions/obj/MVectorTemplate.o


COMMON_ROOT_OBJ=$(COMMON_ROOT_FILEFUNCS) $(COMMON_ROOT_VT) $(COMMON_ROOT_MATHFUNCS)


ROOT_HEADERS=-I`root-config --incdir`
ROOT_LIBS=`root-config --glibs` -lMathMore
ROOTS=$(ROOT_HEADERS) $(ROOT_LIBS)
#
MY_ROOT_FUNCS=-I$(COMMON_ROOT_INC_DIR) $(COMMON_ROOT_OBJ)

#make the text red - used to distinguish cppcheck from g++
TEXT_RED=tput setaf 1;
TEXT_RESET=tput sgr0;

all: testMVectorTemplate.exe
	
clean:
	@rm testMVectorTemplate.exe

testMVectorTemplate.exe: testMVectorTemplate.cxx $(COMMON_ROOT_VT)
	@$(TEXT_RED) $(CPPCHECK) testMVectorTemplate.cxx 
	@$(TEXT_RESET)
	$(CCT) -o testMVectorTemplate.exe testMVectorTemplate.cxx $(MY_ROOT_FUNCS) $(ROOTS)
