SOURCEFILES = benchMark.cpp fft_example.cpp invert_matrix.cpp polynom2.cpp
DEPFILES = $(SOURCEFILES) features_check.h polynom2.h

all : benchMark benchMark_clang

benchMark: $(DEPFILES)
	g++ -O3 -march=native $(SOURCEFILES) -lpthread -o benchMark

benchMark_clang: $(DEPFILES)
	clang++ -O3 $(SOURCEFILES) -lpthread -o benchMark_clang

clean:
	rm -f benchMark benchMark_clang
