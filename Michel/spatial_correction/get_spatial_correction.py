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
from ROOT import gStyle

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

    # Get events only in the detectors inner volume
    for i, (x,y,z,p) in enumerate(data[:, :4]):
        if np.sqrt(x**2 + y**2) < R_FV and Z_FV_LOW < z < Z_FV_HIGH :
            E = data[i, 5]
            if E > 20.0 : 
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

    # Get scale and resolution for each area divided in the detector

    # Read in MC data
    MC_file = ROOT.TFile(mc_filename, "READ")
    MC_tree = MC_file["michel_tree"]
    flux_MC, Edep_MC, x_MC, y_MC, z_MC, w_MC = zip(*[[evt.f, evt.edep, evt.x, evt.y, evt.z, evt.mu_weight] for evt in MC_tree])
    MC_file.Close()
    print("MC data loaded")

    # Calculate energy in MeV from flux
    flux_min = 15e3
    flux_max = 60e3

    inner_flux_MC = []
    inner_flux_w_MC = []
    for f,x,y,z,w in zip(flux_MC, x_MC, y_MC, z_MC, w_MC):
        if np.sqrt(x**2 + y**2) < R_FV and Z_FV_LOW < z < Z_FV_HIGH:
            inner_flux_MC.append(f)
            inner_flux_w_MC.append(w)
    # Get the optimal numbe rof flux bins using the same methos as in get_data_params.y
    inner_flux_MC_slice = [i for i in inner_flux_MC if i < flux_max and i > flux_min]
    n_bins_flux = len(np.histogram_bin_edges(inner_flux_MC_slice, bins='rice')) - 1
    set_MC_hist = {"n_bins": n_bins_flux, "min": flux_min, "max": flux_max}

    plot_MC_flux = ROOT.TH1D("flux_inner_MC", "MC_Flux", n_bins_flux, flux_min, flux_max)
    for f, w in zip(inner_flux_MC, inner_flux_w_MC):
        plot_MC_flux.Fill(f,w)

    MC_flux_fit_vals, MC_flux_errors, MC_flux_fit_graph, MC_flux_fit_chi2 = fit.do_edep_fit(inner_flux_MC, set_MC_hist, weights=inner_flux_w_MC)

    MC_conv = MC_flux_fit_vals[1]
    MC_conv_err = np.sqrt(MC_flux_errors[1][1])
    MC_res_scale = MC_flux_fit_vals[2]
    MC_res_scale_err = np.sqrt(MC_flux_errors[2][2])

    print("\nMC flux to MeV conversion %0.2f ± %0.2f" % (MC_conv, MC_conv_err))

    # Like with the real data, convert MC flux values to MeV using the conv we just found
    E_MC = [ flux/float(MC_conv) for flux in flux_MC ]

    # Get the events in the detectors inner volume
    location_mc_hist = fit.PositionBasedHistograms("MC", "", 12, 0, 2.56e6, 10, -1.25e3, 1.25e3, n_bin, min_E, max_E)
    E_mc_hist = ROOT.TH1D("energy_inner_MC", "MC_Energy", n_bin, min_E, max_E)

    inner_E_MC = []
    inner_pos_MC = []
    inner_w_MC = []
    for E,x,y,z,w in zip(E_MC, x_MC, y_MC, z_MC, w_MC):
        if np.sqrt(x**2 + y**2) < R_FV and Z_FV_LOW < z < Z_FV_HIGH and E > 20.0:
            inner_E_MC.append(E)
            inner_w_MC.append(w)
            inner_pos_MC.append((x,y,z))

            location_mc_hist.Fill(x**2+y**2, z, E, w)
            E_mc_hist.Fill(E, w)

    # Fit the MC inner energy
    MC_fit_vals, MC_fit_errors, MC_fit_graph, MC_fit_chi2 = fit.do_edep_fit(inner_E_MC, set_hist, weights=inner_w_MC, escale=1.0)

    scale_MC = MC_fit_vals[1]
    scale_MC_err = np.sqrt(MC_fit_errors[1][1])
    res_MC = MC_fit_vals[2]
    res_MC_err = np.sqrt(MC_fit_errors[2][2])

    print("\nInner Fit E MC Scale %0.5f ± %0.5f" % (scale_MC, scale_MC_err))
    print("Inner Fit E MC Resolution %0.5f ± %0.5f" % (res_MC, res_MC_err))

    # Create smeared MC data and fit it
    #sigma_smear = delta_sigma(E, res_data, res_MC)
    res_data_avg = np.average(res_data)
    sigma_smear = np.sqrt(res_data_avg**2 - res_MC**2)

    # Bin the data
    #bin_centers = (np.linspace(min_E, max_E, n_bin) + (max_E-min_E)/(2*n_bin))[:-1]
    MC_bins = np.linspace(min_E, max_E, n_bin+1)
    bin_centers = (MC_bins[:-1] + MC_bins[1:])/2
    MC_bin_idx = np.digitize(inner_E_MC, MC_bins)
    binned_MC = np.zeros(len(MC_bins)-1)
    for i, bin_id in enumerate(MC_bin_idx):
        if bin_id == 0 or bin_id == len(MC_bins):
            continue
        binned_MC[bin_id-1] += inner_w_MC[i]

    smeared_MC = fit.apply_smearing(binned_MC, bin_centers, sigma_smear, 0)

    E_mc_g_s = ROOT.TGraph()
    for i, (x,y) in enumerate(zip(bin_centers, smeared_MC)):
        E_mc_g_s.SetPoint(i, x, y)

    E_mc_hist_s = ROOT.TH1D("energy_inner_MC_smeared", "Smeared_MC_Energy", n_bin, min_E, max_E)
    smeared_MC_data = []
    for i, center in enumerate(bin_centers):
        for yval in range(int(smeared_MC[i])):
            E_mc_hist_s.Fill(center)
            smeared_MC_data.append(center)

    MC_fit_vals_s, MC_fit_err_s, MC_fit_graph_s, MC_fit_chi2_s = fit.do_edep_fit(smeared_MC_data, set_hist, escale=1.0)

    scale_MC_s = MC_fit_vals_s[1]
    scale_MC_err_s = np.sqrt(MC_fit_err_s[1][1])
    res_MC_s = MC_fit_vals_s[2]
    res_MC_err_s = np.sqrt(MC_fit_err_s[2][2])

    print("\nInner Fit E Smeared MC Scale %0.5f ± %0.5f" % (scale_MC_s, scale_MC_err_s))
    print("Inner Fit E Smeared MC Resolution %0.5f ± %0.5f" % (res_MC_s, res_MC_err_s))
   
    # Compare data and MC energy scales
    scale_ratio = []
    scale_ratio_err = []
    for period in range(n_periods):
        scale_ratio.append(scale_MC/scale_data[period])
        scale_ratio_err.append(scale_ratio[period]*np.sqrt( (scale_MC_err/scale_MC)**2 + (scale_data_err[period]/scale_data[period])**2 ))

        print("Run Period %i Scale ratio: %0.5f ± %0.5f" % (period, scale_ratio[period], scale_ratio_err[period]))

    corr_h2 = []
    diff_h2 = []
    h_scale = []
    for period in range(n_periods):
        location_mc_hist.Fit()
        location_data_hist[period].Fit()

        loc_data_hist_nbinsX = location_data_hist[period].scale_h2.GetNbinsX()
        loc_data_hist_nbinsY = location_data_hist[period].scale_h2.GetNbinsY()

        data_scale = np.zeros((loc_data_hist_nbinsX, loc_data_hist_nbinsY))
        mc_scale = np.zeros((loc_data_hist_nbinsX, loc_data_hist_nbinsY))
        corr = np.zeros((loc_data_hist_nbinsX, loc_data_hist_nbinsY))
        diff = np.zeros((loc_data_hist_nbinsX, loc_data_hist_nbinsY))

        corr_h2.append(ROOT.TH2D("mc_energy_scale_correction", "", 12, 0, 2.56e3, 10, -1.25e3, 1.25e3))
        diff_h2.append(ROOT.TH2D("data_mc_energy_scale_difference", "", 12, 0, 2.56e3, 10, -1.25e3, 1.25e3))
        h_scale.append(ROOT.TH1D("scale_distribution_for_period", "", 50, 1.5, 0.5))
        corr_avg = []

        for i in range(loc_data_hist_nbinsX):
            for j in range(loc_data_hist_nbinsY):
                d = location_data_hist[period].scale_h2.GetBinContent(i+1, j+1)
                mc = location_mc_hist.scale_h2.GetBinContent(i+1, j+1)
                data_scale[i, j] = d
                mc_scale[i, j] = mc
                diff[i, j] = d-mc
                diff_h2[period].SetBinContent(i+1, j+1, diff[i, j])
                corr[i, j] = d/mc if mc !=0 else 1.0
                corr_h2[period].SetBinContent(i+1, j+1, corr[i, j])
                if d > 0.5 :
                    h_scale[period].Fill(d)

                if mc != 0:
                    corr_avg.append(d/mc)

        corr_avg = np.array(corr_avg)

        print("Average scale factor: "+str(np.average(corr_avg)))

    # Compare the scales from each year to 2021
    corr_21_22 = corr_h2[1].Clone()
    corr_21_22.Divide(corr_h2[0])

    corr_21_23 = corr_h2[2].Clone()
    corr_21_23.Divide(corr_h2[0])

    # Apply the spatial correction to the MC data
    pixel_edges_r = np.linspace(0, 2.56e6, 13)
    pixel_edges_z = np.linspace(-1.25e3, 1.25e3, 11)

    rz_inner_pos_MC = [np.array([(x**2 + y**2), z]) for x,y,z in inner_pos_MC]
    rz_inner_pos_MC = np.array(rz_inner_pos_MC)

    x_bin_idx = np.digitize(rz_inner_pos_MC[:,0], pixel_edges_r)
    y_bin_idx = np.digitize(rz_inner_pos_MC[:,1], pixel_edges_z)

    location_mc_corr_hist = []

    for period in range(n_periods):
        location_mc_corr_hist.append(fit.PositionBasedHistograms("MC", "", 12, 0, 2.56e6, 10, -1.25e3, 1.25e3, n_bin, min_E, max_E))
        for i, E in enumerate(inner_E_MC):
            corr_bin = corr_h2[period].GetBin(int(x_bin_idx[i]), int(y_bin_idx[i]))
            corr_val = corr_h2[period].GetBinContent(corr_bin)
            E_corr = E*(1/corr_val)
            x, y, z = inner_pos_MC[i]
            w = inner_w_MC[i]

            location_mc_corr_hist[period].Fill((x**2 + y**2), z, E_corr, w)

    # Get the corrected MC plots
    diff_corr_scale_h2 = []
    diff_corr_res_h2 = []
    for period in range(n_periods):
        location_mc_corr_hist[period].Fit()

        loc_data_hist_nbinsX = location_data_hist[period].scale_h2.GetNbinsX()
        loc_data_hist_nbinsY = location_data_hist[period].scale_h2.GetNbinsY()

        diff_scale = np.zeros((loc_data_hist_nbinsX, loc_data_hist_nbinsY))
        diff_res = np.zeros((loc_data_hist_nbinsX, loc_data_hist_nbinsY))

        diff_corr_scale_h2.append(ROOT.TH2D("corrected_data_mc_energy_scale_difference", "", 12, 0, 2.56e3, 10, -1.25e3, 1.25e3))
        diff_corr_res_h2.append(ROOT.TH2D("corrected_data_mc_energy_resolution_difference", "", 12, 0, 2.56e3, 10, -1.25e3, 1.25e3))

        for i in range(loc_data_hist_nbinsX):
            for j in range(loc_data_hist_nbinsY):
                d_scale = location_data_hist[period].scale_h2.GetBinContent(i+1, j+1)
                mc_scale = location_mc_corr_hist[period].scale_h2.GetBinContent(i+1, j+1)
                diff_scale[i, j] = d_scale-mc_scale
                diff_corr_scale_h2[period].SetBinContent(i+1, j+1, diff_scale[i, j])

                d_res = location_data_hist[period].resolution_h2.GetBinContent(i+1, j+1)
                mc_res = location_mc_corr_hist[period].resolution_h2.GetBinContent(i+1, j+1)
                diff_res[i, j] = d_res-mc_res
                diff_corr_res_h2[period].SetBinContent(i+1, j+1, diff_res[i, j])

   # Write to outfile
    outFile = ROOT.TFile(output_filename, "RECREATE")
    outFile.cd()

    for period in range(n_periods):
        c5_4 = ROOT.TCanvas("Period %i Data Fit Scale" % period)
        c5_5 = ROOT.TCanvas("Period %i Data Fit Resolution" % period)
        location_data_hist[period].DrawFit(c5_4, c5_5)
        c5_4.Update()
        c5_5.Update()
        c5_4.Draw()
        c5_5.Draw()
        c5_4.Write()
        c5_5.Write()

    c5_0 = ROOT.TCanvas("MC Fit Scale")
    c5_1 = ROOT.TCanvas("MC Fit Resolution")
    location_mc_hist.DrawFit(c5_0, c5_1)
    c5_0.Update()
    c5_1.Update()
    c5_0.Draw()
    c5_1.Draw()
    c5_0.Write()
    c5_1.Write()

    for period in range(n_periods):
        c5_6 = ROOT.TCanvas("Period %i Corrected MC Fit Scale" % period)
        c5_7 = ROOT.TCanvas("Period %i Corrected MC Fit Resolution" % period)
        location_mc_corr_hist[period].DrawFit(c5_6, c5_7)
        c5_6.Update()
        c5_7.Update()
        c5_6.Draw()
        c5_7.Draw()
        c5_6.Write()
        c5_7.Write()

    # Compare data and smeared MC energy scales
    scale_ratio_s = []
    scale_ratio_err_s = []
    for period in range(n_periods):
        scale_ratio_s.append(scale_MC_s/scale_data[period])
        scale_ratio_err_s.append(scale_ratio_s[period]*np.sqrt( (scale_MC_err_s/scale_MC_s)**2 + (scale_data_err[period]/scale_data[period])**2 ))

        print("Run Period %i Scale ratio: %0.5f ± %0.5f" % (period, scale_ratio_s[period], scale_ratio_err_s[period]))

    corr_h2_s = []
    for period in range(n_periods):
        location_mc_hist.SmearFit(sigma_smear)

        loc_data_hist_nbinsX = location_data_hist[period].scale_h2.GetNbinsX()
        loc_data_hist_nbinsY = location_data_hist[period].scale_h2.GetNbinsY()

        data_scale = np.zeros((loc_data_hist_nbinsX, loc_data_hist_nbinsY))
        mc_scale = np.zeros((loc_data_hist_nbinsX, loc_data_hist_nbinsY))
        corr = np.zeros((loc_data_hist_nbinsX, loc_data_hist_nbinsY))

        corr_h2_s.append(ROOT.TH2D("mc_energy_scale_correction", "", 12, 0, 2.56e3, 10, -1.25e3, 1.25e3))
        corr_avg_s = []

        for i in range(loc_data_hist_nbinsX):
            for j in range(loc_data_hist_nbinsY):
                d = location_data_hist[period].scale_h2.GetBinContent(i+1, j+1)
                mc = location_mc_hist.scale_h2.GetBinContent(i+1, j+1)
                data_scale[i, j] = d
                mc_scale[i, j] = mc
                corr[i, j] = d/mc if mc !=0 else 1.0
                corr_h2_s[period].SetBinContent(i+1, j+1, corr[i, j])

                if mc != 0:
                    corr_avg_s.append(d/mc)
                # TODO Save the correction values to json file

        corr_avg_s = np.array(corr_avg_s)

        print("Average scale factor for smeared MC: "+str(np.average(corr_avg_s)))

    outFile.cd()

    c5_2 = ROOT.TCanvas("Smeared MC Fit Scale")
    c5_3 = ROOT.TCanvas("Smeared MC Fit Resolution")
    location_mc_hist.DrawFit(c5_2, c5_3)
    c5_2.Update()
    c5_3.Update()
    c5_2.Draw()
    c5_3.Draw()
    c5_2.Write()
    c5_3.Write()

    # MC flux to MeV
    c8 = ROOT.TCanvas("MC Flux to MeV")
    plot_MC_flux.GetXaxis().SetTitle("Flux")
    plot_MC_flux.GetYaxis().SetTitle("Counts")
    plot_MC_flux.GetYaxis().SetRangeUser(0, 1.2*plot_MC_flux.GetMaximum())
    plot_MC_flux.SetStats(0)
    plot_MC_flux.SetLineColor(ROOT.kBlack)
    plot_MC_flux.SetLineWidth(2)
    MC_flux_fit_graph.SetLineColor(2)
    MC_flux_fit_graph.SetLineWidth(2)
    MC_flux_fit_graph.Draw("al")
    plot_MC_flux.Draw("SAMEE")
    c8.Write()
  

    ntuple_out = ROOT.TNtuple("energy_tree", "energy_tree", "E:x:y:z")
    [ntuple_out.Fill(E,x,y,z) for E, (x, y, z) in zip(smeared_MC, inner_pos_MC)]
    ntuple_out.Write()

    c6 = []
    for period in range(n_periods):
        c6.append(ROOT.TCanvas("MC Energy Scale Correction wrt Run Period %i" % (period)))
        corr_h2[period].SetStats(0)
        corr_h2[period].Draw("TEXT COLZ")
        c6[period].Write()

    c6_1 = []
    for period in range(n_periods):
        c6_1.append(ROOT.TCanvas("Data - MC Energy Scale Difference wrt Run Period %i" % (period)))
        diff_h2[period].SetStats(0)
        diff_h2[period].Draw("TEXT COLZ")
        c6_1[period].Write()

    c6_2 = []
    for period in range(n_periods):
        c6_2.append(ROOT.TCanvas("Smeared MC Energy Scale Correction wrt Run Period %i" % (period)))
        corr_h2_s[period].SetStats(0)
        corr_h2_s[period].Draw("TEXT COLZ")
        c6_2[period].Write()

    c7 = []
    for period in range(n_periods):
        c7.append(ROOT.TCanvas("Data - Corrected MC Energy Resolution Difference wrt Run Period %i" % (period)))
        diff_corr_res_h2[period].SetStats(0)
        diff_corr_res_h2[period].Draw("TEXT COLZ")
        c7[period].Write()

    c7_1 = []
    for period in range(n_periods):
        c7_1.append(ROOT.TCanvas("Data - Corrected MC Energy Scale Difference wrt Run Period %i" % (period)))
        diff_corr_scale_h2[period].SetStats(0)
        diff_corr_scale_h2[period].Draw("TEXT COLZ")
        c7_1[period].Write()

    c6_3 = ROOT.TCanvas("2022/2021 Energy Scale Comparison")
    corr_21_22.SetStats(0)
    corr_21_22.Draw("TEXT COLZ")
    c6_3.Write()

    c6_4 = ROOT.TCanvas("2023/2021 Energy Scale Comparison")
    corr_21_23.SetStats(0)
    corr_21_23.Draw("TEXT COLZ")
    c6_4.Write()

    c6_5 = []
    for period in range(n_periods):
        c6_5.append(ROOT.TCanvas("1D Energy Scale Hist Period %i" % (period)))
        h_scale[period].SetStats(0)
        h_scale[period].Draw("TEXT COLZ")
        c6_5[period].Write()


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
    leg0.AddEntry("", "Ep Res: %0.3f%% #pm %0.3f" % (res_MC*100, res_MC_err*100), "")
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
    E_mc_g_s.SetTitle("Smeared MC")
    E_mc_g_s.GetXaxis().SetTitle("Reconstructed Energy [MeV]")
    E_mc_g_s.GetYaxis().SetTitle("Events/%0.2f MeV" % (bwidth))
    E_mc_g_s.SetStats(0)
    E_mc_g_s.SetLineColor(ROOT.kBlue)
    E_mc_g_s.SetLineWidth(2)
    E_mc_g_s.Draw("AP")
    MC_fit_graph_s.SetLineColor(2)
    MC_fit_graph_s.SetLineWidth(2)
    MC_fit_graph_s.Draw("Same")
    leg2 = ROOT.TLegend(0.60, 0.65, 0.88, 0.86)
    leg2.AddEntry(E_mc_hist_s, "Smeared MC", "l")
    #leg2.AddEntry(E_mc_g_s, "Smeared MC", "l")
    leg2.AddEntry(MC_fit_graph_s, "Michel Fit", "l")
    leg2.AddEntry("", "Scale: %0.5f #pm %0.5f" % (scale_MC_s, scale_MC_err_s), "")
    leg2.AddEntry("", "Ep Res: %0.5f%% #pm %0.5f" % (res_MC_s*100, res_MC_err_s*100), "")
    leg2.Draw("Same")
    c2.Write()
    
    c3 = []
    norm_MC = E_mc_hist.Clone("norm_MC")
    norm_MC.GetXaxis().SetTitle("Reconstructed Energy [MeV]")
    norm_MC.GetYaxis().SetTitle("Normalized")
    norm_MC.SetLineColor(ROOT.kBlue)
    norm_MC.Scale(1.0/norm_MC.Integral())

    norm_smear = E_mc_hist_s.Clone("norm_smear")
    #norm_smear = E_mc_g_s.Clone("norm_smear")
    norm_smear.SetLineColor(ROOT.kGreen+1)
    norm_smear.Scale(1.0/norm_smear.Integral())

    for period in range(n_periods):
        c3.append(ROOT.TCanvas("Normalized Comparison (Run Period %i)" % (period)))

        norm_MC.Draw("HIST")
        norm_smear.Draw("HIST Same")
        #norm_smear.Draw("Same")

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
    
#    c5_0 = ROOT.TCanvas("MC Fit Scale")
#    c5_1 = ROOT.TCanvas("MC Fit Resolution")
#    location_mc_hist.DrawFit(c5_0, c5_1)
#    c5_0.Update()
#    c5_1.Update()
#    c5_0.Draw()
#    c5_1.Draw()
