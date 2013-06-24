#include "plink_binary.h"
#include <iostream>
#include <vector>
#include <cstdio>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>


/*
 * Construct a binary ped file from genotype data in a matrix format
 */

using namespace std;

void line_split(char *s, vector<string> &toks);
char *getline_unlimited(FILE *f, int ignore_cr, int strip_nl);

int main(int argc, char *argv[])
{
    char c;
    string chromosome = "0";
    char missing = 'N';

    while ((c = getopt(argc, argv, "c:m:")) != -1) {
        switch (c) {
            case 'c':
                chromosome = optarg;
                break;
            case 'm':
                missing = *optarg;
                break;
        }
    }

    if (argc <= optind) {
        cout << "Usage: " << argv[0] << " [ options ] PLINK_BINARY" << endl;
        cout << "Options: -c  chromosome number" << endl;
        cout << "         -m  missing genotype character (default " << missing << ")" << endl;
        return 1;
    }

    char *buffer;
    plink_binary *pb = new plink_binary();
    buffer = getline_unlimited(stdin, 1, 1);
    vector<string> data;
    line_split(buffer, data);
    for (unsigned int i = 0; i < data.size(); i++) {
        gftools::individual ind;
        ind.name = data[i];
        pb->individuals.push_back(ind);
    }

    pb->open(argv[optind], 1);
    pb->missing_genotype = missing;
    while ((buffer = getline_unlimited(stdin, 1, 1))) {
        line_split(buffer, data);
        gftools::snp snp;
        snp.name = data[0];
        snp.chromosome = chromosome;
        vector<string> genotypes;
        for (unsigned int i = 1; i < data.size(); i++)
            genotypes.push_back(data[i].substr(0, 2));
        pb->write_snp(snp, genotypes);
    }
    pb->close();
    delete pb;
}

void line_split(char *s, vector<string> &toks)
{
    toks.resize(0);
    char *p = s;
    if (p[strlen(p) - 1] == '\n') p[strlen(p) - 1] = '\0';

    p = strtok(p, "\t");
    while (p != NULL) {
        toks.push_back(string(p));
        p = strtok(NULL, "\t");
    }
    free(p);
}

// Read an arbitrarily long line
char *getline_unlimited(FILE *f, int ignore_cr, int strip_nl)
{
    size_t len = 0;
    static char *buf = NULL;
    static size_t buflen = 0;
    int c;

    while (! feof(f) && (c = fgetc(f)) != '\n' && c != EOF) {
        if (ignore_cr && c == '\r') continue;
        if (len >= buflen) {
            buflen += BUFSIZ;
            if ((buf = (char *)realloc(buf, buflen)) == NULL) {
                fprintf(stderr, "Memory alloc problem reading\n");
                exit(5);
            }
            if (len == 0) buf[0] = 0;
        }
        buf[len++] = c;
    }
    if (buf) {
        if (strip_nl && buf[len] == '\n')
            buf[len] = 0;
        else
            buf[len++] = 0;
    }
    if (feof(f)) return NULL;

    return buf;
}
