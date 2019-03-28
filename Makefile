CXXFLAGS=-Wno-char-subscripts -DNDEBUG -Wall -O3 -std=c++11 -I. -I./sdsl-lite/build/include -I./htslib/htslib -I./KMC -fopenmp
LIBS=-L./sdsl-lite/build/lib -L./sdsl-lite/build/external/libdivsufsort/lib -L./htslib/ -lhts -lz -lsdsl -ldivsufsort -ldivsufsort64

.PHONY: all

all: malva-geno

malva-geno: main.o MurmurHash3.o \
						./KMC/kmc_api/kmc_file.o ./KMC/kmc_api/kmer_api.o ./KMC/kmc_api/mmer.o
	@echo "* Linking $@"
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp
	@echo '* Compiling $<'
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -rf *.o
	rm -rf KMC/kmc_api/*.o
