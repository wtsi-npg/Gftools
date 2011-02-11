CFLAGS = -O3 -Wall -fPIC
LIBS = libplinkbin.so libplinkbin.a plink_binary.so
TARGETS = bed_to_tped plink_binary_to_tab tab_to_plink_binary libplinkbin.a plink_binary.so
CC = g++

%.o : %.cpp
	$(CC) -c $(CFLAGS) -o $@ $<

all: $(TARGETS)

plink_binary_to_tab: plink_binary_to_tab.o plink_binary.o
tab_to_plink_binary: tab_to_plink_binary.o plink_binary.o
bed_to_tped: bed_to_tped.o plink_binary.o

plink_binary.so: plink_binary.o plink_binary.i
	swig -c++ -perl plink_binary.i
	$(CC) $(CFLAGS) -c -I/software/perl-5.8.8/lib/5.8.8/x86_64-linux-thread-multi/CORE/ plink_binary.cpp plink_binary_wrap.cxx `perl -MExtUtils::Embed -e ccopts`
	$(CC) $(CFLAGS) -shared plink_binary.o plink_binary_wrap.o -o plink_binary.so

libplinkbin.a: plink_binary.o
	$(CC) -shared plink_binary.o -o libplinkbin.so
	ar rcs libplinkbin.a plink_binary.o

clean:
	rm *.o $(TARGETS) plink_binary.pm plink_binary_wrap.cxx libplinkbin.so
