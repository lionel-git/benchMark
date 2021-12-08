SOURCEFILES = benchMark.cpp fft_example.cpp

all : benchMark benchMark_clang

benchMark: $(SOURCEFILES)
	g++ -O3 -march=native $(SOURCEFILES) -lpthread -o benchMark

benchMark_clang: $(SOURCEFILES)
	clang++ -O3 $(SOURCEFILES) -lpthread -o benchMark_clang

clean:
	rm -f benchMark benchMark_clang
