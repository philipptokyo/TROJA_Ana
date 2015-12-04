
.PHONY: clean
.PHONY: all

CXX = g++

ifdef ROOTSYS
CFLAGS = -c -g -Wall

LDFLAGS  := `root-config --libs`
CPPFLAGS := `root-config --cflags`
ROOTLIBS  = -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lMinuit -lSpectrum

else
all:
	@echo "ROOT not found!"
endif




all: analysis

analysis: *.o
	$(CXX) $^ -o $@ $(LDFLAGS)

%.o: %.cc
	$(CXX) $(CFLAGS) $(CPPFLAGS) $^

	
clean:
	rm -f *.o *.so analysis
