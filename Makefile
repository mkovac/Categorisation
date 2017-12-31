#makefile

CXXFLAGS = -g -I. -m64 $(shell root-config --cflags) -I include
LDFLAGS = $(shell root-config --libs) -lm -lGenVector
CXX = g++

EXTLIBS = ./ext/cConstants_cc.so ./ext/FinalStates_cc.so ./ext/bitops_cc.so ./ext/Discriminants_cc.so

VPATH = ./src/ ./include/

SRCPP_CATEGORISATION = run_categorisation.cpp\
                       Categorisation.cpp\
                       Histograms.cpp\
                       Variables.cpp\
                       Tree.cpp\
                       Counters.cpp\
                       Utilities.cpp\
                       ROC.cpp\
                       CMS_lumi.cpp

INCLUDES = Categorisation.h\
           Histograms.h\
           Variables.h\
           Tree.h\
           Counters.h\
           Utilities.h\
           ROC.h\
           CMS_lumi.h

    
OBJCPP_CATEGORISATION = $(patsubst %.cpp,obj/%.o,$(SRCPP_CATEGORISATION))


all: run_categorisation
categorisation: run_categorisation


obj/%.o: %.cpp $(INCLUDES)
	@echo ">> compiling $*"
	@mkdir -p obj/
	@$(CXX) -c $< ${CXXFLAGS} -o $@
   

run_categorisation: $(OBJCPP_CATEGORISATION)
	@echo ">> linking..."
	@$(CXX) $^ $(EXTLIBS) ${LDFLAGS} ${CXXFLAGS}  -o $@


clean:
	@echo ">> Cleaning objects and executable..."
	@rm -f obj/*.o
	@rm -f run_categorisation


uninstall:
	@echo ">> Uninstalling..."
	@rm -f obj/*.o
	@rm -f ext/*.so ext/*.d ext/*.pcm
	@rm -f run_categorisation
