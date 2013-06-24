/*
 * Calculate basic SNP stats (CR/AF) and sample CR/het
 *
 * Usage: snp_af_sample_cr [ options ] PLINK_BINARY
*/

#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <iomanip>
#include <getopt.h>
#include "plink_binary.h"

using namespace std;

void allele_counts(vector<int> genotypes, int &a_count, int &b_count, int &nn_count);
void usage(char *progname);
bool sort_by_cr(struct sample s1, struct sample s2);

// to store results so that they can be sorted on call rate
struct sample
{
    string name;
    float cr;
    float het_aut;
    float het_x;
};

int main (int argc, char *argv[])
{
    const char* const short_options = "r:s:va:m:";
    const struct option long_options[] = {
        { "snp", 1, NULL, 'r' },
        { "sample", 1, NULL, 's' },
        { "verbose", 1, NULL, 'v' },
        { "min_snp_cr", 1, NULL, 'm' },
        { NULL, 0, NULL, 0 }
    };
    int opt;

    string snp_file("snp_cr_af.txt");
    string sample_file("sample_cr_het.txt");
    float min_snp_cr = 0.95;
    bool verbose;

    do {
        opt = getopt_long(argc, argv, short_options, long_options, NULL);
        switch(opt) {
            case 'v':
                verbose = 1;
                break;
            case 's':
                sample_file = optarg;
                break;
            case 'm':
                // min_snp_cr = optarg;
                break;
            case 'r':
                snp_file = optarg;
                break;
        }
    } while (opt != -1);

    if (optind >= argc) {
        usage(argv[0]);
        cout << "Usage: " << argv[0] << " [options] BED_FILE" << endl;
        exit(0);
    }

    ofstream out_snp(snp_file.c_str());
    ofstream out_sample(sample_file.c_str());

    out_snp << "#SNP" << "\t" << "CR" << "\t" << "major_allele" << "\t" << "major_allele_freq" << "\t" << "minor_allele" << "\t" << "minor_allele_freq" << endl;
    out_sample << "#Sample" << "\t" << "CR" << "\t" << "autosomal_het" << "\t" << "x_het" << endl;

    plink_binary *pb = new plink_binary(argv[optind]);
    vector<gftools::individual> samples = pb->individuals;
    vector<int> genotypes;
    int total_snps = 0, good_snps = 0;

    vector<int> sample_x_het, sample_x_total;
    vector<int> sample_aut_het, sample_aut_total, sample_other_total;
    sample_x_het.resize(pb->individuals.size());
    sample_aut_het.resize(pb->individuals.size());
    sample_x_total.resize(pb->individuals.size());
    sample_aut_total.resize(pb->individuals.size());
    sample_other_total.resize(pb->individuals.size());

    for (unsigned int snp = 0; snp < pb->snps.size(); snp++) {
        pb->read_snp(snp, genotypes);
        total_snps++;
        int na, nn, nb;
        allele_counts(genotypes, na, nb, nn);
        float snp_cr = (float)(na + nb) / (2 * nn + na + nb);
        out_snp << fixed << pb->snps[snp].name << "\t" << setprecision(4) << snp_cr;
        if (na + nb == 0) {
            // zero CR
            out_snp << "\t.\t.\t.\t." << endl;
            continue;
        }

        // if A is the minor allele - exchange A/B below (*)
        string major = pb->snps[snp].allele_a;
        string minor = pb->snps[snp].allele_b;
        if (major == "0") major = ".";
        if (minor == "0") minor = ".";

        float a_freq = (float)na / (na + nb);
        // *
        if (a_freq < 0.5) {
            a_freq = 1 - a_freq;
            string tmp = major;
            major = minor;
            minor = tmp;
        }
        if (na == 0 || nb == 0)
            out_snp << "\t" << major << "\t" << 1 << "\t" << minor << "\t" << 0 << endl;
        else
            out_snp << "\t" << major << "\t" << setprecision(4) << a_freq << "\t" << minor << "\t" << setprecision(4) << 1 - a_freq << endl;

        if (snp_cr < min_snp_cr)
            continue;
        good_snps++;

        bool x_snp = pb->snps[snp].chromosome == "X" ||
	                 pb->snps[snp].chromosome == "23";
        bool other_snp = pb->snps[snp].chromosome == "XY" ||
                         pb->snps[snp].chromosome == "MT" ||
                         pb->snps[snp].chromosome == "Y" ||
                         pb->snps[snp].chromosome == "24" ||
                         pb->snps[snp].chromosome == "25" ||
                         pb->snps[snp].chromosome == "26";

        for (unsigned int ind = 0; ind < genotypes.size(); ind++) {
            if (genotypes[ind] == 0) continue;
	    if (other_snp) {
                sample_other_total[ind]++;
	    } else if (x_snp) {
                sample_x_total[ind]++;
                if (genotypes[ind] == 2) sample_x_het[ind]++;
            } else {
                sample_aut_total[ind]++;
                if (genotypes[ind] == 2) sample_aut_het[ind]++;
            }
        }
    }
    out_snp.close();

    out_sample << "#stats from " << good_snps << "/" << total_snps << " SNPs with CR >= " << setprecision(2) << 100 * min_snp_cr << "%" << endl;
    out_sample << fixed;

    vector <struct sample> results;
    for (unsigned int ind = 0; ind < samples.size(); ind++) {
        struct sample result = {
            samples[ind].name,
            (float)(sample_aut_total[ind] + sample_x_total[ind] + sample_other_total[ind]) / good_snps,
            (float)sample_aut_het[ind] / sample_aut_total[ind], 
            (float)sample_x_het[ind] / sample_x_total[ind]
        };
        results.push_back(result);
    }
    sort(results.begin(), results.end(), sort_by_cr);
    vector <struct sample>::iterator it;
    for (it = results.begin(); it != results.end(); it++) {
        out_sample << it->name << setprecision(6) << "\t" << it->cr << "\t" << setprecision(4) << it->het_aut << "\t" << setprecision(4) << it->het_x << endl;
    }
}

// return NN, allele counts
void allele_counts(vector<int> genotypes, int &a_count, int &b_count, int &nn_count)
{
    nn_count = a_count = b_count = 0;

    for (unsigned int i = 0; i < genotypes.size(); i++) {
        switch (genotypes[i]) {
            case 0:
                nn_count++;
                break;
            case 1:
                a_count += 2;
                break;
            case 2:
                a_count += 1;
                b_count += 1;
                break;
            case 3:
                b_count += 2;
                break;
        }
    }
}

// sort high to low
bool sort_by_cr(struct sample s1, struct sample s2)
{
    return s1.cr > s2.cr;
}

void usage(char *progname)
{
    cout << "Usage: " << progname << " [options] BED_FILE" << endl;
    cout << "Options: -snp         output snp_file" << endl;
    cout << "         -sample      output sample_file" << endl;
    cout << "         -verbose     verbose" << endl;
    cout << "         -min_snp_cr  min snp call rate" << endl;
}

