# This code takes the combined MC tree files from "build_trees.py" and puts the info into a form that
# the spatial analysis code (get_spatial_correction.py) can use

from rat import RAT
from rat import ROOT

from array import array

# JADE
ROOT.gSystem.Load('libJADE.so')
ROOT.gInterpreter.ProcessLine('using namespace JADE;')
print("Loaded the Libraries")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(prog="Convert MC tree files into a format to be used by the spatial analysis correction code, get_spatial_correction.py", usage='%(prog)s [options]')
    parser.add_argument("plusname", type=str, help="Filename of combined mu+ MC tree files")
    parser.add_argument("minusname", type=str, help="Filename of combined mu- MC tree files")
    parser.add_argument("outfile", type=str, help="Name of the output file")
    args = parser.parse_args()
    muplus_name = args.plusname
    muminus_name = args.minusname
    outName = args.outfile

    outFile = ROOT.TFile(outName, "RECREATE")

    energy_tree = ROOT.TNtuple("energy_tree", "energy_tree", "E_true:flux:x:y:z:PDG")

    # Read in the tree files
    treeFile_plus = ROOT.TFile.Open(muplus_name, "READ")
    treeFile_minus = ROOT.TFile.Open(muminus_name, "READ")

    secondaryTree_plus = treeFile_plus.Get("secondary_truth_tree")
    secondaryTree_minus = treeFile_minus.Get("secondary_truth_tree")

    pairTree_plus = treeFile_plus.Get("pair_tree")
    pairTree_minus = treeFile_minus.Get("pair_tree")

    numPlusEvents = secondaryTree_plus.GetEntries()
    numMinusEvents = secondaryTree_minus.GetEntries()

    for i in range(0, numPlusEvents):
        secondaryTree_plus.GetEntry(i)
        energy_true = secondaryTree_plus.KE
        flux_sum = sum(pairTree_plus.flux_arr_b)
        xpos = secondaryTree_plus.x
        ypos = secondaryTree_plus.y
        zpos = secondaryTree_plus.z
 
        muonpdg = 13
        energy_tree.Fill(energy_true, flux_sum, xpos, ypos, zpos, muonpdg)

    for i in range(0, numMinusEvents):
        secondaryTree_minus.GetEntry(i)
        energy_true = secondaryTree_minus.KE
        flux_sum = sum(pairTree_minus.flux_arr_b)
        xpos = secondaryTree_minus.x
        ypos = secondaryTree_minus.y
        zpos = secondaryTree_minus.z
 
        muonpdg = -13
        energy_tree.Fill(energy_true, flux_sum, xpos, ypos, zpos, muonpdg)
    outFile.cd()
    energy_tree.Write()
    outFile.Close()
    treeFile_plus.Close()
    treeFile_minus.Close()

