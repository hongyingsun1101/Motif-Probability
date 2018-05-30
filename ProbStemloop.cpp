//Code to calculate the probability of exterior loop divided by partition function Q
#include "../src/pfunction.h"
#include "ProbScan.h"
#include "../src/structure.h"
#include <algorithm>
#include <functional>
#include <numeric>
#include <iomanip>
#include <assert.h>
#include <cmath>
using std::vector;
using std::string;
//modified by Hongying based on Probscan.

const static int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},{0,1,0,1,0,0},{0,0,0,0,0,0}};//array representing legal base pairs
const static int max_internal_loop=30;

ProbScan::ProbScan(std::string sequence, bool isRNA):RNA(sequence.c_str(),isRNA)
{
		PartitionFunction();
}

ProbScan::ProbScan(const char filename[], bool from_sequence_file, bool isRNA):RNA(filename,from_sequence_file?FILE_SEQ:FILE_PFS,isRNA)
{
	if (from_sequence_file)
		PartitionFunction();//calculate the partition function if it hasn't been
												//done yet
}

double ProbScan::probability_of_hairpin(int i,int j)
{
	return (double) (v->f(j,i+GetSequenceLength()) //V'(i,j)
				 * erg3(i,j,GetStructure(),pfdata,0)) //K for hairpin
				 / (w5[GetSequenceLength()]*pfdata->scaling*pfdata->scaling); //Q
//divide by scaling^2 so the closing nucs aren't double counted
}

vector<hairpin_t> ProbScan::probability_of_all_hairpins(int min_size, int max_size,double threshold)
{
	vector<hairpin_t> hairpins;
	structure* st = GetStructure();
//#pragma omp parallel for
//parallelization causes datarace because of push_back method
	for(int i=1;i<GetSequenceLength()-min_size-1;i++){//search over all 0<i<j<n
		for(int j=i+min_size+1;j<std::max(i+max_size,GetSequenceLength());j++){
			if(inc[st->numseq[i]][st->numseq[j]]){//if i and j can pair
				//get probability
				double probability = probability_of_hairpin(i,j);
				if (probability>threshold){ //add to the list if p>threshold
					hairpins.push_back(hairpin(probability,i,j));
				}
			}
		}
	}
	//hairpins now contains every hairpin where p>threshold
	return hairpins;
}


//functions for dealing with structures
hairpin_t hairpin(double p,int i, int j)
{
	hairpin_t h;
	h.probability = p;
	h.i = i;
	h.j = j;
	return h;
}





void show_hairpins(vector<hairpin_t> hairpins)//print hairpin output
{
	cout <<"--hairpins--"<<endl;
	cout << "prob i j" <<endl;
	for(vector<hairpin_t>::const_reverse_iterator it=hairpins.rbegin();it!=hairpins.rend();++it)
		cout << std::fixed<<std::setprecision(3)<<it->probability << " " << it->i << " " << it->j <<endl;
	cout<< "--hairpins end--"<<endl <<endl;
}
