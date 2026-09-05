CXX ?= c++
FASTCPD_ROOT ?= ../..
R_ABSEIL_INCLUDE := $(shell Rscript -e 'cat(system.file("include", package = "abseil"))' 2>/dev/null)
R_ARMADILLO_INCLUDE := $(shell Rscript -e 'cat(system.file("include", package = "RcppArmadillo"))' 2>/dev/null)
ABSL_INCLUDE ?= $(R_ABSEIL_INCLUDE)
ARMADILLO_INCLUDE ?= $(R_ARMADILLO_INCLUDE)
FASTCPD_CPP := $(FASTCPD_ROOT)/cpp/fastcpd.cc

CPPFLAGS += -DNO_RCPP -DARMA_DONT_USE_WRAPPER -I$(FASTCPD_ROOT)/include -I$(FASTCPD_ROOT)/src
ifneq ($(strip $(ABSL_INCLUDE)),)
CPPFLAGS += -I$(ABSL_INCLUDE)
endif
ifneq ($(strip $(ARMADILLO_INCLUDE)),)
CPPFLAGS += -I$(ARMADILLO_INCLUDE)
endif
CXXFLAGS += -std=c++17 -O2
LDLIBS += -llapack -lblas

.PHONY: all clean

all: mean_change interface_contract shared_contract

mean_change: mean_change.cc fastcpd.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< fastcpd.o -o $@ $(LDLIBS)

interface_contract: interface_contract.cc fastcpd.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< fastcpd.o -o $@ $(LDLIBS)

shared_contract: shared_contract.cc fastcpd.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< fastcpd.o -o $@ $(LDLIBS)

fastcpd.o: $(FASTCPD_CPP) $(FASTCPD_ROOT)/include/fastcpd/fastcpd.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(FASTCPD_CPP) -o $@

clean:
	rm -f mean_change interface_contract shared_contract fastcpd.o
