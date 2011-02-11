#include "plink_binary.h"
#include <fstream>
#include <vector>

using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 3) return 1;

    plink_binary *pb = new plink_binary(argv[1]);
    pb->missing_genotype = '0';

    string data(argv[2]);
    string fn(data + ".tfam");
    ofstream tfam(fn.c_str());
    for (unsigned int i = 0; i < pb->individuals.size(); i++) {
        tfam << pb->individuals[i].family << " ";
        tfam << pb->individuals[i].name << " ";
        tfam << pb->individuals[i].father << " ";
        tfam << pb->individuals[i].mother << " ";
        tfam << pb->individuals[i].sex << " ";
        tfam << pb->individuals[i].phenotype << endl;
    }
    tfam.close();

    fn = data + ".tped";
    ofstream tped(fn.c_str());

    gftools::snp snp;
    vector<string> gt;

    while (pb->next_snp(snp, gt)) {
        tped << snp.chromosome << " ";
        tped << snp.name << " ";
        tped << snp.genetic_position << " ";
        tped << snp.physical_position;
        for (unsigned int i = 0; i < gt.size(); i++) {
            tped << " " << gt[i][0] << " " << gt[i][1];
        }
        tped << endl;
    }

    tped.close();
    pb->close();
    delete pb;
}
