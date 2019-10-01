
benchMark: benchMark.cpp
	g++ -O3 benchMark.cpp -lpthread -o benchMark

clean:
	rm benchMark
