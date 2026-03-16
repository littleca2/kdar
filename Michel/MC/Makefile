CXX = g++
ROOTFLAGS = $(shell root-config --cflags --libs)
JADELIBS = -L$(JADEROOT)/lib
JADEINC = -I$(JADEROOT)/include
OPTFLAGS = -O0
DEBUGFLAGS = -g

TARGETS = mu_weight_tree
CXXFLAGS = -Wall $(OPTFLAGS) $(DEBUGFLAGS) $(ROOTFLAGS) -lMinuit

all: $(TARGETS)

%: %.cc
	$(CXX) -o $@ $^ $(CXXFLAGS)

clean:
	rm -f $(TARGETS) *.o
