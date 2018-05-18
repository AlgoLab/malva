CFLAGS	= -Wno-char-subscripts -DNDEBUG -Wall -O3 -std=c++11 -I. -I./sdsl-lite/build/include -I./htslib/htslib -I./KMC -fopenmp
CXXFLAGS= ${CFLAGS}
LIBS = -L./sdsl-lite/build/lib -L./sdsl-lite/build/external/libdivsufsort/lib -L./htslib/ -lhts -lz -lsdsl -ldivsufsort -ldivsufsort64

.PHONY: all

all: malva

malva: main.o MurmurHash3.o \
						./KMC/kmc_api/kmc_file.o ./KMC/kmc_api/kmer_api.o ./KMC/kmc_api/mmer.o
	@echo "* Linking malva"
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(LDFLAGS)

%.o: %.cpp
	@echo '* Compiling $<'
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -rf *.o
	rm -rf KMC/kmc_api/*.o
