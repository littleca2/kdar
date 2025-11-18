# This code takes the combined MC tree files from "build_trees.py" and puts the info into a form that
# the spatial analysis code (ana_sys.py) can use

from rat import RAT
from rat import ROOT

# JADE
ROOT.gSystem.Load('libJADE.so')
ROOT.gInterpreter.ProcessLine('using namespace JADE;')
print("Loaded the Libraries")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(prog="Convert MC tree files into a format to be used by the spatial analysis correction code, ana_sys.py", usage='%(prog)s [options]')
    parser.add_argument("plusname", type=str, help="Filename of combined mu+ MC tree files")
    parser.add_argument("minusname", type=str, help="Filename of combined mu- MC tree files")
    parser.add_argument("outfile", type=str, help="Name of the output file")
    args = parser.parse_args()
    muplus_name = args.plusname
    muminus_name = args.minusname
    outName = args.outfile

    outFile = ROOT.TFile(outName, "RECREATE")

    energy_tree = ROOT.TNtuple("energy_tree", "energy_tree", "E:x:y:z:PDG")

    # Read in the tree files
    treeFile_plus = ROOT.TFile.Open(muplus_name, "READ")
    treeFile_minus = ROOT.TFile.Open(muminus_name, "READ")

    secondaryTree_plus = treeFile_plus.Get("secondary_truth_tree")
    secondaryTree_minus = treeFile_minus.Get("secondary_truth_tree")

    numPlusEvents = secondaryTree_plus.GetEntries()
    numMinusEvents = secondaryTree_minus.GetEntries()

    for i in range(0, numPlusEvents):
        secondaryTree_plus.GetEntry(i)
        energy = secondaryTree_plus.KE
        xpos = secondaryTree_plus.x
        ypos = secondaryTree_plus.y
        zpos = secondaryTree_plus.z
 
        muonpdg = 13
        energy_tree.Fill(energy, xpos, ypos, zpos, muonpdg)

    for i in range(0, numMinusEvents):
        secondaryTree_minus.GetEntry(i)
        energy = secondaryTree_minus.KE
        xpos = secondaryTree_minus.x
        ypos = secondaryTree_minus.y
        zpos = secondaryTree_minus.z
 
        muonpdg = -13
        energy_tree.Fill(energy, xpos, ypos, zpos, muonpdg)
    outFile.cd()
    energy_tree.Write()
    outFile.Close()
    treeFile_plus.Close()
    treeFile_minus.Close()

