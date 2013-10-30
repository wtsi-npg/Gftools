/*
 * Calculate pairwise genotype concordance using a subset of SNPs
 *
 * Usage: snp_af_sample_cr [ options ] PLINK_BINARY
 * NB: stores intermediate results in a short int, hence limited to 65535 SNPs
*/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <set>
#include <getopt.h>
#include "plink_binary.h"

using namespace std;

void usage(char *progname);

int main (int argc, char *argv[])
{
    const char* const short_options = "d:n:r:f:m:";
    const struct option long_options[] = {
        { "snp", 1, NULL, 'n' },
        { "full", 1, NULL, 'f' },
        { "summary", 1, NULL, 'm' },
        { "duplicate", 1, NULL, 'd' },
        { NULL, 0, NULL, 0 }
    };
    int opt;

    string snp_file;
    string full_file("duplicate_full.txt");
    string summary_file("duplicate_summary.txt");
    float dup_threshold = 0.98;

    do {
        opt = getopt_long(argc, argv, short_options, long_options, NULL);
        switch(opt) {
            case 'f':
                full_file = optarg;
                break;
            case 'm':
                summary_file = optarg;
                break;
            case 'n':
                snp_file = optarg;
                break;
            case 'd':
                // duplicate = optarg;
                break;
                break;
        }
    } while (opt != -1);

    if (optind >= argc) {
        usage(argv[0]);
        exit(0);
    }

    ofstream out_full(full_file.c_str());
    ofstream out_summary(summary_file.c_str());

    out_full << "#(D)uplicate/(R)elated\tSample_1\tSample_2\tConcordance (#SNPs)\n";
    out_summary << "#(D)uplicate/(R)elated\tSample_1\tSample_2\tConcordance (#SNPs)\n";

    plink_binary *pb = new plink_binary(argv[optind]);
    vector<gftools::individual> samples = pb->individuals;
    vector<int> genotypes;
    set <string> snps_to_check;

    ifstream snps;
    snps.open(snp_file.c_str());
    string s;
    while (getline(snps, s)) {
        if (snps_to_check.size() >= 65535) { // unsigned short
            cout << "Maximum SNP count reached" << endl;
            exit (1);
        }
        snps_to_check.insert(s);
    }
    snps.close();

    int n_samples = samples.size();
    int n_pairs = n_samples * (n_samples - 1) / 2;

    vector <unsigned short> checked_for_pair, matched_for_pair;
    checked_for_pair.resize(n_pairs);
    matched_for_pair.resize(n_pairs);
    for (int i = 0; i < n_pairs; i++) {
        checked_for_pair[i] = matched_for_pair[i] = 0;
    }

    for (unsigned int snp = 0; snp < pb->snps.size(); snp++) {
        if (snps_to_check.size() == 0)
            // terminate if all requested SNPs analysed
            break;
        // this SNP requested?
        if (snps_to_check.find(pb->snps[snp].name) != snps_to_check.end())
            snps_to_check.erase(pb->snps[snp].name);
        else
            continue;
        pb->read_snp(snp, genotypes);

        int pair = 0;
        for (int ind_1 = 0; ind_1 < n_samples; ind_1++) {
            for (int ind_2 = 1 + ind_1; ind_2 < n_samples; ind_2++) {
                if (genotypes[ind_1] && genotypes[ind_2]) {
                    checked_for_pair[pair]++;
                    if (genotypes[ind_1] == genotypes[ind_2]){
                      matched_for_pair[pair]++;
                    }
                }

                pair++;
            }
        }
    }

    int pair = 0;
    for (int ind_1 = 0; ind_1 < n_samples; ind_1++) {
        for (int ind_2 = 1 + ind_1; ind_2 < n_samples; ind_2++) {
            float match = (float)matched_for_pair[pair] / checked_for_pair[pair];
            stringstream match_str;
            match_str << fixed;
            match_str << pb->individuals[ind_1].name << "\t" << pb->individuals[ind_2].name << "\t";
            match_str << setprecision(4) << match << "\t (" << checked_for_pair[pair] << ")" << endl;
            string match_type;
            if (match > dup_threshold) {
                out_full << "D\t" << match_str.str();
                out_summary << "D\t" << match_str.str();
            } else if (match > dup_threshold) {
                out_full << "R\t" << match_str.str();
            } else {
                out_full << "0\t" << match_str.str();
            }

            pair++;
        }
    }

    out_full.close();
    out_summary.close();
}

void usage(char *progname)
{
    cout << "Usage: " << progname << " [options] BED_FILE" << endl;
}

