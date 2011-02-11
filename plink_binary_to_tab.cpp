#include "plink_binary.h"
#include <iostream>

/*
 * Extract genotypes and output in a simple matrix
 */

using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " DATA_FILE" << endl;
        return 1;
    }

    plink_binary *pb;
    try {
        pb = new plink_binary(argv[1]);
    } catch (exception &e) {
        cout << "Error opening: " << e.what() << endl;
        return 1;
    }
    pb->missing_genotype = '0';

    for (unsigned int i = 0; i < pb->individuals.size(); i++)
        cout << "\t" << pb->individuals[i].name;
    cout << endl;

    vector<string> genotypes;
    gftools::snp snp;
    while (pb->next_snp(snp, genotypes)) {
        cout << snp.name;
        for (unsigned int i = 0; i < genotypes.size(); i++)
            cout << "\t" << genotypes[i];
        cout << endl;
    }
}
