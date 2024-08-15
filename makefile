SOURCEFILES = benchMark.cpp fft_example.cpp invert_matrix.cpp polynom2.cpp
DEPFILES = $(SOURCEFILES) features_check.h polynom2.h

# default values
MARCH=native
CLANG_TARGET=benchMark_clang
ICX_TARGET=benchMark_icx

# Need sourcing /opt/intel/oneapi/setvars.sh for icx to be found
CHECKICX=$(shell /usr/bin/which icx > /dev/null 2>&1 && echo "FOUND" )
#$(info checkicx = x$(CHECKICX)x)
ifneq ("$(CHECKICX)", "FOUND")
  ICX_TARGET=
endif
#$(info icxtarget = $(ICX_TARGET))

MACHINE = $(shell /usr/bin/uname -m)
ifeq ($(MACHINE), aarch64)
  MARCH=armv8-a
else ifeq ($(MACHINE), riscv64)
  MARCH=rv64gc
  CLANG_TARGET=
else ifneq ($(MACHINE), x86_64)
  $(error unsupported architecture: '$(MACHINE)')
endif

all: benchMark $(CLANG_TARGET) $(ICX_TARGET)

benchMark: $(DEPFILES)
	g++ -O3 -march=$(MARCH) -Wall -Wpedantic $(SOURCEFILES) -lpthread -o benchMark

benchMark_clang: $(DEPFILES)
	clang++ -O3 -march=$(MARCH) -Weverything -Wno-c++98-compat -Wno-missing-prototypes -Wno-c++98-compat-pedantic $(SOURCEFILES) -lpthread -o benchMark_clang

benchMark_icx: $(DEPFILES)
	icx -lstdc++ -O3 -march=$(MARCH) -Wall -Wpedantic $(SOURCEFILES) -lpthread -o benchMark_icx

clean:
	rm -f benchMark benchMark_clang benchMark_icx
