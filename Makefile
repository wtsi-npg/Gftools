CFLAGS = -O3 -Wall -fPIC
LIBS = libplinkbin.so
MODULES = plink_binary.pm
OBJS = bed_to_tped plink_binary_to_tab tab_to_plink_binary
TARGETS = $(OBJS) $(LIBS)
INCLUDES = exceptions.h individual.h plink_binary.h snp.h

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

libplinkbin.so: plink_binary.o
	$(CC) -shared plink_binary.o -o libplinkbin.so
	ar rcs libplinkbin.a plink_binary.o

clean:
	rm *.o *.a *.cxx $(TARGETS) plink_binary.pm

install:
	cp $(INCLUDES) $(INSTALL_INC)
	cp $(LIBS) $(INSTALL_LIB)
	cp $(LIBS:so=a) $(INSTALL_LIB)
	cp $(MODULES:pm=so) $(INSTALL_LIB)
	cp $(MODULES) $(INSTALL_LIB)
	cp $(OBJS) $(INSTALL_BIN)
