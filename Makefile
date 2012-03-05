
LIB_NAME = plinkbin
LIB_VERSION = $(shell grep '[[:digit:]].[[:digit:]].[[:digit:]]' VERSION)
LIB_BASE = $(LIB_NAME)
LIB = lib$(LIB_BASE).so

LIBS = $(LIB)
MODULES = plink_binary.pm
OBJS = bed_to_tped plink_binary_to_tab tab_to_plink_binary
TARGETS = $(OBJS) $(LIBS)
INCLUDES = utilities.h exceptions.h individual.h plink_binary.h snp.h
CXXTEST_ROOT := /usr/local/lib/cxxtest

PERL_CORE = $(shell perl -e 'print join ":", map{ <$$_/*/CORE> } @INC')

INSTALL_ROOT = ...
INSTALL_INC = $(INSTALL_ROOT)/include
INSTALL_LIB = $(INSTALL_ROOT)/lib
INSTALL_BIN = $(INSTALL_ROOT)/bin

CC = g++
CFLAGS = -O3 -Wall -fPIC
LIBPATH = -L./
LDFLAGS = -l$(LIB_BASE)

.PHONY: test

%.o : %.cpp
	$(CC) -c $(CFLAGS) -o $@ $<

all: $(TARGETS) $(MODULES)

plink_binary_to_tab: plink_binary_to_tab.o $(LIB)
	$(CC) $(LDFLAGS) $(LIBPATH) $^ -o $@
tab_to_plink_binary: tab_to_plink_binary.o $(LIB)
	$(CC) $(LDFLAGS) $(LIBPATH) $^ -o $@
bed_to_tped: bed_to_tped.o $(LIB)
	$(CC) $(LDFLAGS) $(LIBPATH) $^ -o $@

plink_binary.pm: plink_binary.i
	swig -c++ -perl plink_binary.i
	$(CC) $(CFLAGS) -c -I$(PERL_CORE) plink_binary.cpp plink_binary_wrap.cxx `perl -MExtUtils::Embed -e ccopts`
	$(CC) $(CFLAGS) -shared -L$(PERL_CORE) utilities.o plink_binary.o plink_binary_wrap.o -o plink_binary.so -lperl

$(LIB): utilities.o plink_binary.o
	$(CC) -shared utilities.o plink_binary.o -o $@

runner.cpp: test_plink_binary.h
	$(CXXTEST_ROOT)/bin/cxxtestgen -o $@ --error-printer $^

runner: runner.cpp $(LIB)
	$(CC) -Wall -g -I$(CXXTEST_ROOT) $(LDFLAGS) $(LIBPATH) $^ -o $@

test: runner
	LD_LIBRARY_PATH=. ./runner

clean:
	rm -f *.o *.so *.cxx $(TARGETS) plink_binary.pm

install:
	cp $(INCLUDES) $(INSTALL_INC)
	cp $(LIBS) $(INSTALL_LIB)
	cp $(MODULES:pm=so) $(INSTALL_LIB)
	cp $(MODULES) $(INSTALL_LIB)
	cp $(OBJS) $(INSTALL_BIN)
