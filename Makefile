CC=g++
CLANG_TIDY=clang-tidy
DEBUG=-Wall -Wextra -Wpedantic -Wunused-variable
CXX11=-std=c++1y
OPTIM=-O3
CATCH_TESTS=catch2Compile/tests-main.o

CCT=$(CC) $(OPTIM) $(DEBUG) $(CXX11)

SRC_DIR=src
OBJ_DIR=obj
TEST_DIR=tests
INC_DIR=include
SHARED_DIR=lib
TB_DIR=testbeam
EXEC_DIR=executables

#rootana
ifdef ROOTANASYS
ROOTANAINC = -I$(ROOTANASYS)/include
ROOTANALIBS = $(ROOTANASYS)/lib/librootana.a
else
ROOTANAINC = -I../include
ROOTANALIBS = ../lib/librootana.a
endif

EXTERNALS_INCLUDE = -I$(BELLE2_EXTERNALS_DIR)/include

INC=-I$(INC_DIR) $(ROOTANAINC)

ROOT_HEADERS=-I`root-config --incdir`

#make the text red - used to distinguish cppcheck from g++
TEXT_RED=tput setaf 1;
TEXT_RESET=tput sgr0;

OBJ_FILES=$(patsubst $(SRC_DIR)/%.cxx,$(OBJ_DIR)/%.o,$(wildcard $(SRC_DIR)/*.cxx))
TEST_EXEC_FILES=$(patsubst $(TEST_DIR)/%.cxx,$(TEST_DIR)/$(EXEC_DIR)/%.exe, $(shell find $(TEST_DIR)/ -type f -name '*.cxx'))
TB_OBJ_FILES=$(patsubst $(SRC_DIR)/$(TB_DIR)/%.cxx,$(OBJ_DIR)/$(TB_DIR)/%.o,$(wildcard $(SRC_DIR)/$(TB_DIR)/*.cxx))

#Inspired by https://stackoverflow.com/questions/2394609/makefile-header-dependencies
DEP = $(TEST_EXEC_FILES:%.exe=%.d)

all: $(OBJ_FILES) $(TB_OBJ_FILES) $(SHARED_DIR)/libbasf2Tools.so $(SHARED_DIR)/libtestBeam.so

clean:
	-@rm $(OBJ_DIR)/* || true
	-@rm $(OBJ_DIR)/$(TB_DIR)/* || true
	-@rm $(SHARED_DIR)/* || true
	-@rm -r $(TEST_DIR)/$(EXEC_DIR)/* || true


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cxx $(INC_DIR)/%.h
	$(CCT) -c $(INC) $(ROOT_HEADERS) $(EXTERNALS_INCLUDE) -fPIC -o $@ $<
	clang-check $< -- $(INC) $(ROOT_HEADERS) $(CXX11)

# $(OBJ_DIR)/rootDictionalry.o:
# 	-rootcint -f $(patsubst %.o,%.cxx,$@) $(INC_DIR)/eclCrystalDB.h $(INC_DIR)/fileFuncs.h
# 	mv $(patsubst %.o,%.cxx_tmp,$@) $(patsubst %.o,%.cxx,$(@F))   #Hack because rootcint doesnt change file name
# 	mv $(patsubst %.o,%_rdict.pcm,$@) $(SHARED_DIR)
# 	$(info var is $@)
# 	$(CCT) -c -fPIC -o $@ $(patsubst %.o,%.cxx,$(@F))
# 	rm $(patsubst %.o,%.cxx,$(@F))                                #remove dictionart.cxx - not needed for anything

$(SHARED_DIR)/libbasf2Tools.so: $(OBJ_DIR)/myRootStyle.o $(OBJ_DIR)/fileFuncs.o $(OBJ_DIR)/eclCrystalDB.o $(OBJ_DIR)/fileFuncs.o $(OBJ_DIR)/mathFuncs.o $(OBJ_DIR)/histFuncs.o $(OBJ_DIR)/stringFuncs.o $(OBJ_DIR)/generalFuncs.o $(OBJ_DIR)/fftFuncs.o $(OBJ_DIR)/Windows.o $(OBJ_DIR)/MVectorTemplate.o $(OBJ_DIR)/ParameterComparisonPlot.o $(OBJ_DIR)/filterFuncs.o
	$(CCT) -shared -o $@ $^ `root-config --glibs`

$(SHARED_DIR)/libtestBeam.so: $(OBJ_DIR)/$(TB_DIR)/testbeamFuncs.o $(OBJ_DIR)/$(TB_DIR)/RunDB.o $(OBJ_DIR)/$(TB_DIR)/TOFtiming.o $(OBJ_DIR)/$(TB_DIR)/EventLoopBase.o $(OBJ_DIR)/$(TB_DIR)/MidasLoop.o $(OBJ_DIR)/$(TB_DIR)/Waveform.o
	$(CCT)  -shared -o $@ $^ -lbasf2Tools -L$(ROOTANASYS)/lib -lrootana `root-config --glibs` -lXMLParser -lXMLIO #XML libraries are for rootana


##################################
#catch2
##################################
catch2: $(CATCH_TESTS)

$(CATCH_TESTS):
	g++ catch2Compile/tests-main.cpp -c -o $@

##################################
#Tests
##################################
tests: $(TEST_EXEC_FILES)

# Include all .d files
-include $(DEP)

$(TEST_DIR)/$(EXEC_DIR)/%.exe: $(TEST_DIR)/%.cxx
	-$(shell mkdir -p $(dir $@))
	$(CCT) -MMD $(CATCH_TESTS) $(INC) $(ROOT_HEADERS) $(EXTERNALS_INCLUDE) -L$(SHARED_DIR) -lbasf2Tools `root-config --glibs` -o $@ $<

print-%:
	@echo '$*=$($*)'
