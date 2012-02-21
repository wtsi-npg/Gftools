CFLAGS = -O3 -Wall -fPIC
LIBS = libplinkutils.so libplinkbin.so
MODULES = plink_binary.pm plink_utils.pm
OBJS = bed_to_tped plink_binary_to_tab tab_to_plink_binary
TARGETS = $(OBJS) $(LIBS)
INCLUDES = exceptions.h individual.h plink_binary.h plink_utils.h snp.h

LDFLAGS =
INSTALL_ROOT = ...
INSTALL_INC = $(INSTALL_ROOT)/include
INSTALL_LIB = $(INSTALL_ROOT)/lib
INSTALL_BIN = $(INSTALL_ROOT)/bin
CC = g++

%.o : %.cpp
	$(CC) -c $(CFLAGS) -o $@ $<

all: $(TARGETS) $(MODULES)

plink_binary_to_tab: plink_binary_to_tab.o plink_binary.o
	g++ plink_binary_to_tab.o plink_binary.o -o plink_binary_to_tab
tab_to_plink_binary: tab_to_plink_binary.o plink_binary.o
	g++ tab_to_plink_binary.o plink_binary.o -o tab_to_plink_binary
bed_to_tped: bed_to_tped.o plink_binary.o
	g++ bed_to_tped.o plink_binary.o -o bed_to_tped

plink_binary.pm: plink_binary.i
	swig -c++ -perl plink_binary.i
	$(CC) $(CFLAGS) -c -I/software/perl-5.8.8/lib/5.8.8/x86_64-linux-thread-multi/CORE/ plink_binary.cpp plink_binary_wrap.cxx `perl -MExtUtils::Embed -e ccopts`
	$(CC) $(CFLAGS) -shared plink_binary.o plink_binary_wrap.o -o plink_binary.so

plink_utils.pm: plink_utils.i
	swig -c++ -perl plink_utils.i
	$(CC) $(CFLAGS) -c -I/software/perl-5.8.8/lib/5.8.8/x86_64-linux-thread-multi/CORE/ plink_utils.cpp plink_utils_wrap.cxx `perl -MExtUtils::Embed -e ccopts`
	$(CC) $(CFLAGS) -shared plink_utils.o plink_utils_wrap.o -o plink_utils.so

libplinkbin.so: plink_binary.o
	$(CC) -shared plink_binary.o -o libplinkbin.so
	ar rcs libplinkbin.a plink_binary.o

libplinkutils.so: plink_utils.o
	$(CC) -shared plink_utils.o -o libplinkutils.so
	ar rcs libplinkutils.a plink_utils.o

clean:
	rm *.o *.a *.cxx $(TARGETS) plink_binary.pm plink_utils.pm

install:
	cp $(INCLUDES) $(INSTALL_INC)
	cp $(LIBS) $(INSTALL_LIB)
	cp $(LIBS:so=a) $(INSTALL_LIB)
	cp $(MODULES:pm=so) $(INSTALL_LIB)
	cp $(MODULES) $(INSTALL_LIB)
	cp $(OBJS) $(INSTALL_BIN)
