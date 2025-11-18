CXX = g++
ROOTFLAGS = $(shell root-config --cflags --libs)
JADELIBS = -L$(JADEROOT)/lib
JADEINC = -I$(JADEROOT)/include
JADEFLAGS = $(JADEINC) $(JADELIBS) -lJADE
OPTFLAGS = -O0
DEBUGFLAGS = -g

TARGETS = preProcessingForMichel update_michel_pair
CXXFLAGS = -Wall $(OPTFLAGS) $(DEBUGFLAGS) $(ROOTFLAGS) $(JADEFLAGS) -lMinuit

all: $(TARGETS)

%: %.cc
	$(CXX) -o $@ $^ $(CXXFLAGS)

clean:
	rm -f $(TARGETS) *.o
