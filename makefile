SOURCEFILES = benchMark.cpp fft_example.cpp invert_matrix.cpp polynom2.cpp
DEPFILES = $(SOURCEFILES) features_check.h polynom2.h

all : benchMark benchMark_clang

benchMark: $(DEPFILES)
	g++ -O3 -march=native -Wall -Wpedantic $(SOURCEFILES) -lpthread -o benchMark

benchMark_clang: $(DEPFILES)
	clang++ -O3 -Weverything -Wno-c++98-compat -Wno-missing-prototypes -Wno-c++98-compat-pedantic $(SOURCEFILES) -lpthread -o benchMark_clang

#-Weverything
#-Wc++98-compat

clean:
	rm -f benchMark benchMark_clang
