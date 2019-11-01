
all : benchMark benchMark_clang

benchMark: benchMark.cpp
	g++ -O3 benchMark.cpp -lpthread -o benchMark

benchMark_clang: benchMark.cpp
	clang++ -O3 benchMark.cpp -lpthread -o benchMark_clang

clean:
	rm -f benchMark benchMark_clang
