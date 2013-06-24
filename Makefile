
LIB_VERSION = $(shell grep '[[:digit:]].[[:digit:]].[[:digit:]]' VERSION)

MODULES = plink_binary.pm
EXECUTABLES = bed_to_tped plink_binary_to_tab tab_to_plink_binary snp_af_sample_cr_bed pairwise_concordance_bed
LIBS = libplinkbin.so libplinkbin.a
TARGETS = $(EXECUTABLES) $(LIBS)
INCLUDES = utilities.h exceptions.h individual.h plink_binary.h snp.h
CXXTEST_ROOT ?= /usr/local/lib/cxxtest

INSTALL_ROOT = ...
INSTALL_INC = $(INSTALL_ROOT)/include
INSTALL_LIB = $(INSTALL_ROOT)/lib
INSTALL_BIN = $(INSTALL_ROOT)/bin

CXX = g++
CXXFLAGS = -O3 -Wall -fPIC
AR = ar
LIBPATH = -L./
LDFLAGS = $(LIBPATH) -lplinkbin

.PHONY: test clean

%.o : %.cpp
	$(CXX) -c $(CXXFLAGS) -o $@ $<

all: $(TARGETS) $(MODULES)

plink_binary_to_tab: plink_binary_to_tab.o libplinkbin.so
	$(CXX) $< $(LDFLAGS) -o $@
tab_to_plink_binary: tab_to_plink_binary.o libplinkbin.so
	$(CXX) $< $(LDFLAGS) -o $@
bed_to_tped: bed_to_tped.o libplinkbin.so
	$(CXX) $< $(LDFLAGS) -o $@
snp_af_sample_cr_bed: snp_af_sample_cr_bed.o libplinkbin.so
	$(CXX) $< $(LDFLAGS) -o $@

pairwise_concordance_bed: pairwise_concordance_bed.o
	$(CXX) $< $(LDFLAGS) -o $@

plink_binary.pm: plink_binary.i
	swig -c++ -perl plink_binary.i
	$(CXX) $(CXXFLAGS) -c plink_binary.cpp plink_binary_wrap.cxx `perl -MExtUtils::Embed -e ccopts`
	$(CXX) $(CXXFLAGS) -shared `perl -MExtUtils::Embed -e ldopts` utilities.o plink_binary.o plink_binary_wrap.o -o plink_binary.so

libplinkbin.so: utilities.o plink_binary.o
	$(CXX) -shared utilities.o plink_binary.o -o $@

libplinkbin.a: utilities.o plink_binary.o
	$(AR) rcs $@ $^

runner.cpp: test_plink_binary.h
	$(CXXTEST_ROOT)/bin/cxxtestgen -o $@ --error-printer $^

runner: runner.cpp libplinkbin.so
	$(CXX) -Wall -g -I$(CXXTEST_ROOT) $< $(LDFLAGS) -o $@

test: runner
	LD_LIBRARY_PATH=. ./runner

clean:
	rm -f *.o *.a *.so *.cxx $(TARGETS) plink_binary.pm

install: all
	cp $(INCLUDES) $(INSTALL_INC)
	cp $(LIBS) $(INSTALL_LIB)
	cp $(MODULES:pm=so) $(INSTALL_LIB)
	cp $(MODULES) $(INSTALL_LIB)
	cp $(EXECUTABLES) $(INSTALL_BIN)
