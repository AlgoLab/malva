
CXXFLAGS=-Wno-char-subscripts -DNDEBUG -Wall -O3 -std=c++11
LIBS=-L./lib -L./htslib/ -lhts

.PHONY: all

all: malva-test

malva-test: main.o malva_test.o 
	@echo "* Linking $@"
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp
	@echo '* Compiling $<'
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -rf *.o