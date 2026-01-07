import numpy as np
import ROOT
import os
import csv
from scipy.optimize import curve_fit
import json
import sys
import pandas as pd
sys.path.insert(0, '/home/littleca/kdar/Michel/flux2mev_correction/')
import fitting as fit

JSON_FILE = "/home/littleca/kdar/correction_values.json"
OUTPUT_PATH = "/home/littleca/kdar/Michel/spatial_correction/"

def delta_sigma(E, R_data, R_MC):
    E_ep = fit.MICHEL_E_ENDPOINT	# MeV
    end_point_smear = np.sqrt(R_data**2 - R_MC**2)
    return np.sqrt(E_ep/E)*end_point_smear*E

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Get the MC simulation energy corrections using the flux to MeV correction values from data.")
    parser.add_argument("data_input", type=str, help="input ROOT filename (updateMichelPair output; michel_pair_combined.root)")
    parser.add_argument("mc_input", type=str, help="input ROOT filename")
    parser.add_argument("version", type=str, help="title of JSON version for this group of data")
    parser.add_argument("--output", type=str, help="output ROOT filename")
    args = parser.parse_args()

    data_filename = args.data_input
    mc_filename = args.mc_input
    versionID = args.version

    if args.output:
        output_filename = args.output
    else :
        outDir = OUTPUT_PATH+str(versionID)
        if not os.path.isdir(outDir) :
            os.makedirs(outDir)
        output_filename = outDir+"/spatialCorr_"+str(versionID)+".root"
    
    ROOT.TH1.AddDirectory(0)

    data_file = ROOT.TFile(data_filename, "READ")
    bin_r = 12
    bin_z = 10

    R_max = 1600 # mm
    Rs_max = (R_max**2) # mm^2
    z_max = 1250 # mm
    R_FV = 1400
    Z_FV_LOW = -1000
    Z_FV_HIGH = 1000


    # Energy histogram binning
    min_E = 20
    max_E = 70
    bwidth = 0.25 # bin width
    n_bin = int((max_E - min_E)/bwidth)


    # Read in correction values from json file
    with open(JSON_FILE, "r") as f:
        data_vals = json.load(f)
    print("JSON correction values loaded")

    data_json_vals = data_vals[versionID]
    df_json = pd.DataFrame(data_json_vals)
    n_periods = int(df_json.max(axis=1)["RunPeriod"]+1)

    print("Found %i run periods" % (n_periods))

    # We read in the run data for runs that are in the JSON file
    data = []
    prev_run = 0
    no_period = 0
    for evt in data_file["event_tree"]:
        if evt.run != prev_run:
            prev_run = evt.run
            try:
                period = data_json_vals["run_"+str(int(evt.run))]["RunPeriod"]
                no_period = 0
            except KeyError:
                no_period = 1
                if evt.run != None:
                    print("Could not find run period information for run %i; Ignoring run." % (int(evt.run)))
                    continue
                else:
                    continue
        if not no_period:
            data.append(np.array([evt.x, evt.y, evt.z, int(period), evt.flux, 0, evt.run, evt.subrun]))

    data = np.array(data)

    data[:, :3] *= 1e3	# Convert data coordinates to millimeters

    for i, value_dict in enumerate(data_json_vals.values()):
        np.divide(data[:, 4], value_dict["Flux2MeV"], out=data[:,5], where=data[:,6]==value_dict["Run"])
    # List to hold event energies for each run period
    inner_E_data = [[] for i in range(n_periods)]
    # Lists that hold histograms to be plotted later
    location_data_hist = []
    E_data_hist = []
    for p in range(n_periods):
        location_data_hist.append(fit.PositionBasedHistograms("Data", "", 12, 0, 2.56e6, 10, -1.25e3, 1.25e3, n_bin, min_E, max_E))
        E_data_hist.append(ROOT.TH1D("energy_inner_data", "", n_bin, min_E, max_E))

    for i, (x,y,z,p) in enumerate(data[:, :4]):
        if np.sqrt(x**2 + y**2) < R_FV and Z_FV_LOW < z < Z_FV_HIGH :
            E = data[i, 5]
            inner_E_data[int(p)].append(E)

            location_data_hist[int(p)].Fill(x**2+y**2, z, E)
            E_data_hist[int(p)].Fill(E)

    # Fit the MeV converted data for all the runs
    set_hist = {"n_bins": n_bin, "min": min_E, "max": max_E}

    scale_data = [[] for i in range(n_periods)]
    scale_data_err = [[] for i in range(n_periods)]
    res_data = [[] for i in range(n_periods)]
    res_data_err = [[] for i in range(n_periods)]
    data_fit_graphs = []
    for p in range(n_periods):
        data_fit_vals, data_fit_errors, data_fit_graph, data_fit_chi2 = fit.do_edep_fit(np.array(inner_E_data[p]), set_hist, escale=1.0)

        data_fit_graphs.append(data_fit_graph)

        scale_data[p] = data_fit_vals[1]
        scale_data_err[p] = np.sqrt(data_fit_errors[1][1])
        res_data[p] = data_fit_vals[2]
        res_data_err[p] = np.sqrt(data_fit_errors[2][2])

        print("\nRun Period %i Inner Fit E Data Scale %0.5f ± %0.5f" % (p, scale_data[p], scale_data_err[p]))
        print("Run Period %i Inner Fit E Data Resolution %0.5f ± %0.5f" % (p, res_data[p], res_data_err[p]))

    scale_data = np.array(scale_data)
    scale_data_err = np.array(scale_data_err)
    res_data = np.array(res_data)
    res_data_err = np.array(res_data_err)

    # Read in MC data
    MC_file = ROOT.TFile(mc_filename, "READ")
    MC_tree = MC_file["energy_tree"]
    x_MC, y_MC, z_MC, E_MC, ptype_MC = zip(*[[evt.x, evt.y, evt.z, evt.E, evt.PDG] for evt in MC_tree])
    MC_file.Close()
    print("MC data loaded")

    # Calculate the MC weights for the mu-plus and mu-minus proportion
    plus_count = len([x for x in ptype_MC if x==13])
    minus_count = len([x for x in ptype_MC if x==-13])
    plus_weight = float(plus_count)/(plus_count + minus_count)
    minus_weight = float(minus_count)/(plus_count + minus_count)
    # Ratio mu+/mu- = 1.21 from PRD 74.082006
    plus_weight *= 1.21
    minus_weight *= 1.

    ##### Perform Fit of the Inner Energy Hist for MC Data #####
    inner_MC = [[E, (x, y, z), ptype] for E,x,y,z,ptype in zip(E_MC, x_MC, y_MC, z_MC, ptype_MC) if np.sqrt(x**2 + y**2) < R_FV and Z_FV_LOW < z < Z_FV_HIGH and E > 20.0]

    inner_pos_MC = [pos for E, pos, ptype in inner_MC]
    # We need to apply the mu+/- weights to the list of event energies so the fit matches with the hist
    # So we go ahead and bin the data according to our set_hist parameters
    # and then "re-create" the now weighted data by representing it with a list of the 
    # bin center values for each histogram entry.
    inner_MC_arr = np.array([[E, ptype] for E, pos, ptype in inner_MC])
    E_MC_bins = np.linspace(min_E, max_E, n_bin+1)
    E_MC_bin_centers = np.linspace((min_E + bwidth/2), (max_E - bwidth/2), n_bin)
    E_MC_nphist = np.histogram(inner_MC_arr[:,0], bins=E_MC_bins, weights=inner_MC_arr[:,1])[0]
    inner_E_MC_weighted_centers = []
    for binN, yval in enumerate(E_MC_nphist):
        for count in range(int(yval)):
            inner_E_MC_weighted_centers.append(E_MC_bin_centers[binN])

    MC_fit_vals, MC_fit_errors, MC_fit_graph, MC_fit_chi2 = fit.do_edep_fit(inner_E_MC_weighted_centers, set_hist, escale=1.0)

    scale_MC = MC_fit_vals[1]
    scale_MC_err = np.sqrt(MC_fit_errors[1][1])
    res_MC = MC_fit_vals[2]
    res_MC_err = np.sqrt(MC_fit_errors[2][2])

    print("\nInner Fit E MC Scale %0.5f ± %0.5f" % (scale_MC, scale_MC_err))
    print("Inner Fit E MC Resolution %0.5f ± %0.5f" % (res_MC, res_MC_err))

    # Fill the MC events position based histograms
    location_mc_hist = fit.PositionBasedHistograms("MC", "", 12, 0, 2.56e6, 10, -1.25e3, 1.25e3, n_bin, min_E, max_E)
    E_mc_hist = ROOT.TH1D("energy_inner_MC", "MC_Energy", n_bin, min_E, max_E)

    [location_mc_hist.Fill(x**2+y**2, z, E) for E, (x, y, z), ptype in inner_MC]
    [E_mc_hist.Fill(E, plus_weight if ptype==13 else minus_weight) for E, pos, ptype in inner_MC]

    # Compare data and MC energy scales
    scale_ratio = []
    scale_ratio_err = []
    for period in range(n_periods):
        scale_ratio.append(scale_MC/scale_data[period])
        scale_ratio_err.append(scale_ratio[period]*np.sqrt( (scale_MC_err/scale_MC)**2 + (scale_data_err[period]/scale_data[period])**2 ))

        print("Run Period %i Scale ratio: %0.5f ± %0.5f" % (period, scale_ratio[period], scale_ratio_err[period]))

    corr_h2 = []
    for period in range(n_periods):
        location_mc_hist.Fit()
        location_data_hist[period].Fit()

        loc_data_hist_nbinsX = location_data_hist[period].scale_h2.GetNbinsX()
        loc_data_hist_nbinsY = location_data_hist[period].scale_h2.GetNbinsY()

        data_scale = np.zeros((loc_data_hist_nbinsX, loc_data_hist_nbinsY))
        mc_scale = np.zeros((loc_data_hist_nbinsX, loc_data_hist_nbinsY))
        corr = np.zeros((loc_data_hist_nbinsX, loc_data_hist_nbinsY))

        corr_h2.append(ROOT.TH2D("mc_energy_scale_correction", "", 12, 0, 2.56e3, 10, -1.25e3, 1.25e3))
        corr_avg = []

        for i in range(loc_data_hist_nbinsX):
            for j in range(loc_data_hist_nbinsY):
                d = location_data_hist[period].scale_h2.GetBinContent(i+1, j+1)
                mc = location_mc_hist.scale_h2.GetBinContent(i+1, j+1)
                data_scale[i, j] = d
                mc_scale[i, j] = mc
                corr[i, j] = d/mc if mc !=0 else 1.0
                corr_h2[period].SetBinContent(i+1, j+1, corr[i, j])

                if mc != 0:
                    corr_avg.append(d/mc)

        corr_avg = np.array(corr_avg)

        print("Average scale factor: "+str(np.average(corr_avg)))

    # Create smeared MC data and fit it
    #bin_centers = (np.linspace(min_E, max_E, n_bin) + (max_E-min_E)/(2*n_bin))[:-1]
    #sigma_smear = delta_sigma(E, res_data, res_MC)
    res_data_avg = np.average(res_data)
    sigma_smear = np.sqrt(res_data_avg**2 - res_MC**2)
    smeared_MC = fit.apply_smearing(E_MC_nphist, E_MC_bin_centers, sigma_smear, 0)

    E_mc_hist_s = ROOT.TH1D("smeared_MC", "smeared_MC", n_bin, min_E, max_E)
    [E_mc_hist_s.Fill(E) for E in smeared_MC]

    MC_fit_vals_s, MC_fit_err_s, MC_fit_graph_s, MC_fit_chi2_s = fit.do_edep_fit(smeared_MC, set_hist, escale=1.0)

    scale_MC_s = MC_fit_vals_s[1]
    scale_MC_err_s = np.sqrt(MC_fit_err_s[1][1])
    res_MC_s = MC_fit_vals_s[2]
    res_MC_err_s = np.sqrt(MC_fit_err_s[2][2])

    # Write to outfile
    outFile = ROOT.TFile(output_filename, "RECREATE")
    outFile.cd()

    ntuple_out = ROOT.TNtuple("energy_tree", "energy_tree", "E:x:y:z")
    [ntuple_out.Fill(E,x,y,z) for E, (x, y, z) in zip(smeared_MC, inner_pos_MC)]
    ntuple_out.Write()

    c6 = []
    for period in range(n_periods):
        c6.append(ROOT.TCanvas("MC Energy Scale Correction wrt Run Period %i" % (period)))
        corr_h2[period].Draw("TEXT COLZ")
        c6[period].Write()

    c0 = ROOT.TCanvas("Inner Michel Reco E: MC Data")
    E_mc_hist.SetTitle("MC")
    E_mc_hist.GetXaxis().SetTitle("Reconstructed Energy [MeV]")
    E_mc_hist.GetYaxis().SetTitle("Events/%0.2f MeV" % bwidth)
    E_mc_hist.Draw("HIST")
    E_mc_hist.SetStats(0)
    MC_fit_graph.SetLineColor(2)
    MC_fit_graph.SetLineWidth(2)
    MC_fit_graph.Draw("sameL")
    leg0 = ROOT.TLegend(0.60, 0.65, 0.88, 0.86)
    leg0.AddEntry(E_mc_hist, "MC Data", "l")
    leg0.AddEntry(MC_fit_graph, "Michel Fit", "l")
    leg0.AddEntry("", "Scale: %0.5f #pm %0.5f" % (scale_MC, scale_MC_err), "")
    leg0.AddEntry("", "Ep Res: %0.3f%% #pm %0.3f" % (res_MC, res_MC_err), "")
    leg0.Draw("same")
    c0.Write()

    c1 = []
    leg1 = []
    for period in range(n_periods):
        c1.append(ROOT.TCanvas("Inner Michel Reco E: Real Data Period %i" % (period)))
        E_data_hist[period].SetTitle("Data Period %i" % (period))
        E_data_hist[period].GetXaxis().SetTitle("Reconstructed Energy [MeV]")
        E_data_hist[period].GetYaxis().SetTitle("Events/%0.2f MeV" % (bwidth))
        E_data_hist[period].SetStats(0)
        E_data_hist[period].SetLineColor(ROOT.kBlack)
        E_data_hist[period].SetLineWidth(2)
        E_data_hist[period].Draw("HISTE")
        data_fit_graphs[period].SetLineColor(2)
        data_fit_graphs[period].SetLineWidth(2)
        data_fit_graphs[period].Draw("Same")
        leg1.append(ROOT.TLegend(0.60, 0.65, 0.88, 0.86))
        leg1[period].AddEntry(E_data_hist[period], "Real Data", "l")
        leg1[period].AddEntry(data_fit_graphs[period], "Michel Fit", "l")
        leg1[period].AddEntry("", "Scale: %0.5f #pm %0.5f" % (scale_data[period], scale_data_err[period]), "")
        leg1[period].AddEntry("", "Ep Res: %0.3f%% #pm %0.3f" % (res_data[period]*100, res_data_err[period]*100), "")
        leg1[period].Draw("Same")
        c1[period].Write()

    c2 = ROOT.TCanvas("Inner Michel Reco E: Smeared MC")
    E_mc_hist_s.SetTitle("Smeared MC")
    E_mc_hist_s.GetXaxis().SetTitle("Reconstructed Energy [MeV]")
    E_mc_hist_s.GetYaxis().SetTitle("Events/%0.2f MeV" % (bwidth))
    E_mc_hist_s.SetStats(0)
    E_mc_hist_s.SetLineColor(ROOT.kBlue)
    E_mc_hist_s.SetLineWidth(2)
    E_mc_hist_s.Draw("HISTE")
    MC_fit_graph_s.SetLineColor(2)
    MC_fit_graph_s.SetLineWidth(2)
    MC_fit_graph_s.Draw("Same")
    leg2 = ROOT.TLegend(0.60, 0.65, 0.88, 0.86)
    leg2.AddEntry(E_mc_hist_s, "Smeared MC", "l")
    leg2.AddEntry(MC_fit_graph_s, "Michel Fit", "l")
    leg2.AddEntry("", "Scale: %0.5f #pm %0.5f" % (scale_MC_s, scale_MC_err_s), "")
    leg2.AddEntry("", "Ep Res: %0.5f%% #pm %0.5f" % (res_MC_s, res_MC_err_s), "")
    leg2.Draw("Same")
    c2.Write()
    
    c3 = []
    norm_MC = E_mc_hist.Clone("norm_MC")
    norm_MC.GetXaxis().SetTitle("Reconstructed Energy [MeV]")
    norm_MC.GetYaxis().SetTitle("Normalized")
    norm_MC.SetLineColor(ROOT.kBlue)
    norm_MC.Scale(1.0/norm_MC.Integral())

    norm_smear = E_mc_hist_s.Clone("norm_smear")
    norm_smear.SetLineColor(ROOT.kGreen+1)
    norm_smear.Scale(1.0/norm_smear.Integral())

    for period in range(n_periods):
        c3.append(ROOT.TCanvas("Normalized Comparison (Run Period %i)" % (period)))

        norm_MC.Draw("HIST")
        norm_smear.Draw("HIST Same")

        norm_data = E_data_hist[period].Clone("norm_data")
        norm_data.SetLineColor(1)
        norm_data.Scale(1.0/norm_data.Integral())
        norm_data.Draw("HIST Same")

        leg3 = ROOT.TLegend(0.12, 0.65, 0.4, 0.86)
        leg3.AddEntry(norm_MC, "MC")
        leg3.AddEntry(norm_data, "Data Run Period %i " % (period))
        leg3.AddEntry(norm_smear, "Smeared MC")
        leg3.Draw("Same")
        c3[period].Write()

    c4_0 = []
    c4_1 = []
    for period in range(n_periods):
        c4_0.append(ROOT.TCanvas("Data Fit Scale Run Period %i" % (period)))
        c4_1.append(ROOT.TCanvas("Data Fit Resolution Run Period %i" % (period)))
        location_data_hist[period].DrawFit(c4_0[period], c4_1[period])
        c4_0[period].Update()
        c4_1[period].Update()
        c4_0[period].Draw()
        c4_1[period].Draw()
    
    c5_0 = ROOT.TCanvas("MC Fit Scale")
    c5_1 = ROOT.TCanvas("MC Fit Resolution")
    location_mc_hist.DrawFit(c5_0, c5_1)
    c5_0.Update()
    c5_1.Update()
    c5_0.Draw()
    c5_1.Draw()
