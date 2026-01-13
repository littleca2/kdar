import ROOT
import numpy as np
# This file simply takes an input "analysis_tree" file and scrapes the michel
# events out of it, then puts the flux,x,y,z values for those events in the output
# file. This doesn't perform "selection" per-se however.
# Tbh, I don't really understand how Alex demarcates events in the "pair_tree"
# so how exactly Michels are being identified is a bit fuzzy to me


total_flux_hist = ROOT.TH1D("total_flux", "total_flux", 100, 0, 45000)
back_flux_hist = ROOT.TH1D("back_flux", "back_flux", 100, 0, 45000)
michel_flux_hist = ROOT.TH1D("michel_flux", "michel_flux", 100, 0, 45000)
total_tree = ROOT.TNtuple("total_tree", "total_tree", "f:edep:x:y:z:mc_x:mc_y:mc_z")
michel_tree = ROOT.TNtuple("michel_tree", "michel_tree", "f:edep:x:y:z:mc_x:mc_y:mc_z")
back_tree = ROOT.TNtuple("back_tree", "back_tree", "f:edep:x:y:z:mc_x:mc_y:mc_z")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Extract the michel event data from an input analysis tree")
    parser.add_argument("input", type=str, help="input ROOT filename, can include wildcard to concatenate files")
    parser.add_argument("output", type=str, help="output ROOT filename")
    args = parser.parse_args()
    input_filename = args.input
    output_filename = args.output

    pair_tree = ROOT.TChain("pair_tree")
    pair_tree.Add(input_filename)

    for i in range(pair_tree.GetEntries()):
        pair_tree.GetEntry(i)
        E_d = pair_tree.E_d_b
        E_dep = pair_tree.E_dep_b
        t_z = pair_tree.z_d_b
        # reco info
        trig_arr = pair_tree.Trig_arr_b
        sub_arr = pair_tree.SubE_arr_b
        flux_arr = pair_tree.flux_arr_b
        TTT_arr = pair_tree.TTT_arr_b
        x_arr = pair_tree.x_reco_arr_b
        y_arr = pair_tree.y_reco_arr_b
        z_arr = pair_tree.z_reco_arr_b

        mc_x = pair_tree.x_d_b
        mc_y = pair_tree.y_d_b
        mc_z = pair_tree.z_d_b

        for f, x, y, z, trig in zip(flux_arr, x_arr, y_arr, z_arr, trig_arr):
            if trig >= 0:
                total_flux_hist.Fill(f)
                total_tree.Fill(f, E_dep, x, y, z, mc_x, mc_y, mc_z)
        if E_d == -1:
            # Not a michel
            for f, x, y, z, trig in zip(flux_arr, x_arr, y_arr, z_arr, trig_arr):
                if trig >= 0:
                    back_flux_hist.Fill(f)
                    back_tree.Fill(f, E_dep, x, y, z, mc_x, mc_y, mc_z)
        else:
            for f, x, y, z, TTT, subtrig in zip(flux_arr, x_arr, y_arr, z_arr, TTT_arr, sub_arr):
                if(f<=0):
                    continue
                if TTT > 0:
                    back_flux_hist.Fill(f)
                    back_tree.Fill(f, E_dep, x, y, z, mc_x, mc_y, mc_z)
                elif subtrig == 0 and t_z < 1500:
                    # A michel!
                    michel_flux_hist.Fill(f)
                    michel_tree.Fill(f, E_dep, x, y, z, mc_x, mc_y, mc_z)


    C1 = ROOT.TCanvas("Michel Info")

    total_flux_hist.GetXaxis().SetTitle("Flux")
    total_flux_hist.SetTitle("")
    total_flux_hist.SetLineColor(1)
    total_flux_hist.Draw("HIST")
    michel_flux_hist.SetLineColor(4)
    michel_flux_hist.Draw("HIST Same")
    back_flux_hist.SetLineColor(2)
    back_flux_hist.Draw("HIST Same")

    leg = ROOT.TLegend(0.6, 0.5, 0.9, 0.7)
    leg.AddEntry(total_flux_hist, "Total Flux")
    leg.AddEntry(michel_flux_hist, "Michel Spectrum")
    leg.AddEntry(back_flux_hist, "Background")
    leg.Draw("Same")
    C1.Update()

    outFile = ROOT.TFile(output_filename, "RECREATE")
    outFile.cd()

    total_tree.Write()
    michel_tree.Write()
    back_tree.Write()
    outFile.Close()
