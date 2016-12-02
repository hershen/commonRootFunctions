CC=g++
DEBUG=-Wall -Wextra -Wpedantic -Wunused-variable
CXX11=-std=c++1y
OPTIM=-O3


CCT=$(CC) $(OPTIM) $(DEBUG) $(CXX11)

SRC_DIR=src
OBJ_DIR=obj
INC_DIR=include
SHARED_DIR=lib

INC=-I$(INC_DIR)


ROOT_HEADERS=-I`root-config --incdir`
# ROOTlibs=`root-config --glibs` -lMathMore
# ROOTS=$(ROOTheaders) $(ROOTlibs)

#make the text red - used to distinguish cppcheck from g++
TEXT_RED=tput setaf 1;
TEXT_RESET=tput sgr0;

OBJ_FILES=$(patsubst $(SRC_DIR)/%.cxx,$(OBJ_DIR)/%.o,$(wildcard $(SRC_DIR)/*.cxx))
# SHARED_LIBS=$(SHARED_DIR)/eclCrystalDB.so $(SHARED_DIR)/fileFuncs.so
# sensitiveFiles=$(OBJ_FILES) $(SHARED_LIBS)
# $(info var is $(OBJ_FILES))
#all: $(OBJ_DIR)/fileFuncs.o $(OBJ_DIR)/filterFuncs.o $(OBJ_DIR)/mathFuncs.o $(OBJ_DIR)/btFuncs.o $(OBJ_DIR)/MVectorTemplate.o
all: $(OBJ_FILES)  $(SHARED_DIR)/libbasf2Tools.so #$(OBJ_DIR)/rootDictionalry.o

clean:
	-@rm $(OBJ_DIR)/*
	-@rm $(SHARED_DIR)/*
	
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cxx $(INC_DIR)/%.h
	$(CCT) -c $(INC) $(ROOT_HEADERS) -fPIC -o $@ $<
	
# $(OBJ_DIR)/eclCrystalDB.o: $(SRC_DIR)/eclCrystalDB.cxx $(INC_DIR)/eclCrystalDB.h
# 	@$(TEXT_RED) $(CPPCHECK) $<
# 	@$(TEXT_RESET)
# 	$(CCT) -c $(INC) -fPIC -o $@ $<
# 	
# $(OBJ_DIR)/fileFuncs.o: $(SRC_DIR)/fileFuncs.cxx $(INC_DIR)/fileFuncs.h
# 	@$(TEXT_RED) $(CPPCHECK) $<
# 	@$(TEXT_RESET)
# 	$(CCT) -c $(INC) -fPIC -o $@ $<

# $(OBJ_DIR)/rootDictionalry.o: 
# 	-rootcint -f $(patsubst %.o,%.cxx,$@) $(INC_DIR)/eclCrystalDB.h $(INC_DIR)/fileFuncs.h 
# 	mv $(patsubst %.o,%.cxx_tmp,$@) $(patsubst %.o,%.cxx,$(@F))   #Hack because rootcint doesnt change file name
# 	mv $(patsubst %.o,%_rdict.pcm,$@) $(SHARED_DIR)
# 	$(info var is $@)
# 	$(CCT) -c -fPIC -o $@ $(patsubst %.o,%.cxx,$(@F))
# 	rm $(patsubst %.o,%.cxx,$(@F))                                #remove dictionart.cxx - not needed for anything
	
$(SHARED_DIR)/libbasf2Tools.so: $(OBJ_DIR)/eclCrystalDB.o $(OBJ_DIR)/fileFuncs.o $(OBJ_DIR)/histFuncs.o $(OBJ_DIR)/stringFuncs.o #$(OBJ_DIR)/rootDictionalry.o
	$(info var is $^)
	$(CCT) -shared -o $@ $^ `root-config --glibs`
	
	
