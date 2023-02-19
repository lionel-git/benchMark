SOURCEFILES = benchMark.cpp fft_example.cpp invert_matrix.cpp polynom2.cpp
DEPFILES = $(SOURCEFILES) features_check.h polynom2.h

MACHINE = $(shell /usr/bin/uname -m)
ifeq ($(MACHINE), aarch64)
	MARCH=armv8-a
	CLANG_TARGET=benchMark_clang
else ifeq ($(MACHINE), riscv64)
	MARCH=rv64gc
else
	MARCH=native
	CLANG_TARGET=benchMark_clang
endif

all: benchMark $(CLANG_TARGET)

benchMark: $(DEPFILES)
	g++ -O3 -march=$(MARCH) -Wall -Wpedantic $(SOURCEFILES) -lpthread -o benchMark

benchMark_clang: $(DEPFILES)
	clang++ -O3 -march=$(MARCH) -Weverything -Wno-c++98-compat -Wno-missing-prototypes -Wno-c++98-compat-pedantic $(SOURCEFILES) -lpthread -o benchMark_clang

clean:
	rm -f benchMark benchMark_clang
