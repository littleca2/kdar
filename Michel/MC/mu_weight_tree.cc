// This is made to take the mu+ and mu- MC root files (with pair_tree)
// and create weight trees to be used later to more accurately reflect mu+/mu- ratio
//
// So ultimately, the output will be same same as each original file, but with the inclusion
// of a weight tree. When the mu+ & mu- files are merged with hadd, it will also merge the 
// weight tree which can be used during analysis by declaring it as a friend of pair_tree

#include <TTree.h>
#include <TFile.h>

#include <iostream>
#include <sstream>
#include <string>

int main(int argc, char** argv) {
	if(argc != 3) {
		std::cout << "Must provide input filename for mu plus MC and mu minus MC" << std::endl;
		return -1;
	}

	std::string plus_filename = argv[1];
	std::string minus_filename = argv[2];

	std::string output_filename_plus = "friend_" + plus_filename;
	std::string output_filename_minus = "friend_" + minus_filename;

	TFile* f_plus = new TFile(plus_filename.c_str());
	TTree* p_tree_plus = (TTree*)f_plus->Get("pair_tree");
	TFile* f_minus = new TFile(minus_filename.c_str());
	TTree* p_tree_minus = (TTree*)f_minus->Get("pair_tree");

	double mu_weight;

	TFile* ff_plus = new TFile(output_filename_plus.c_str(), "recreate");
	TTree* w_tree_plus = new TTree("weight_tree", "Mu weights");
	w_tree_plus->Branch("mu_weight", &mu_weight, "mu_weight/D");

	TFile* ff_minus = new TFile(output_filename_minus.c_str(), "recreate");
	TTree* w_tree_minus = new TTree("weight_tree", "Mu weights");
	w_tree_minus->Branch("mu_weight", &mu_weight, "mu_weight/D");

	// Calculate the weights
	int nentries_plus = p_tree_plus->GetEntries();
	int nentries_minus = p_tree_minus->GetEntries();
	int total = nentries_plus + nentries_minus;

	double plus_weight = nentries_plus/(double)total;
	double minus_weight = nentries_minus/(double)total;

	// Ratio mu+/mu- = 1.21 from PRD 74.082006
	plus_weight *= 1.21;
	minus_weight *= 1.0;

	for (int i=0; i<nentries_plus; i++) {
		mu_weight = plus_weight;
		w_tree_plus->Fill();
	}
	w_tree_plus->Write();
	ff_plus->Write();
	ff_plus->Close();

	for (int i=0; i<nentries_minus; i++) {
		mu_weight = minus_weight;
		w_tree_minus->Fill();
	}
	w_tree_minus->Write();
	ff_minus->Write();
	ff_minus->Close();
}

