# -*- coding: utf-8 -*-


import numpy as np
import ROOT
import os
import csv
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import fitting as fit
from array import array
import json
#debug
import time
import sys
import pdb

JSON_FILE = "/home/littleca/kdar/correction_values.json"
	
w_pixel = 2 # bin width for indiv pixels
pixel_bins = int(70/w_pixel)

def Sigma(E, R_ep, constant_term=0.02):
    E_ep = fit.MICHEL_E_ENDPOINT # MeV
    p = constant_term
    x = E/E_ep
    sigma = E_ep*x*np.sqrt( ((R_ep**2 - p**2)/x) + p**2 )
    return sigma

def delta_sigma(E, R_data, R_MC):
    E_ep = fit.MICHEL_E_ENDPOINT # MeV
    end_point_smear = np.sqrt(R_data**2 - R_MC**2)
    return np.sqrt(E_ep/E)*end_point_smear*E

def chandong_correction(flux):
    p0 = 1.00039e+00
    p1 = 4.28201e-05
    p2 = -4.27272e+00
    return p0 - np.exp(-p1*flux+p2)

MASKED_PMTS_RUNS = np.array([1609,1785,1788,1790,1807,1809,1810,1811,1812,1813,1814,1815,1817,1828,1837,1843,1850,1878,1879,1909,1910,1912,1914,1916,1920,1922,1923,1924,1932,1935,1936,1940,1942,1944,1949,1950,1952])

def apply_chandong_correction(run, flux):
    mask = np.in1d(run, MASKED_PMTS_RUNS)
    factors = chandong_correction(flux[mask])
    flux[mask] = flux[mask]/factors


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("data_input", type=str, help="input ROOT filename (updateMichelPair output; michel_pair_combined.root)")
    parser.add_argument("mc_input", type=str, help="input ROOT filename")
    parser.add_argument("output", type=str, help="output ROOT filename")
    parser.add_argument("version", type=str, help="title of JSON version for this group of data")
    args = parser.parse_args()

    data_filename = args.data_input
    mc_filename = args.mc_input
    output_filename = args.output
    versionID = args.version
    
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

    location_data_hist = fit.PositionBasedHistograms("Data", "", 12, 0, 2.56e6, 10, -1.25e3, 1.25e3, n_bin, min_E, max_E)


	# Read in correction values from json file
    with open(JSON_FILE, "r") as f:
        data_vals = json.load(f)

    print("JSON correction values loaded")
    #pdb.set_trace()
    # Create a new dictionary with the Flux2MeV values for each run in the specified version
    f2mev_dict = {value["Run"]: value["Flux2MeV"] for key, value in data_vals[versionID].items()}

    # We read in the run data for runs that are in the JSON file
    # If a run is not in the JSON file, it means there were no events that passed the Michel cuts
    data_dtype = np.dtype([('pos', np.float, (3,)), ('flux', np.float), ('energy', np.float), ('run', np.int), ('subrun', np.int)])
    data = np.array([((evt.x, evt.y, evt.z), evt.flux, 0, evt.run, evt.subrun) for evt in data_file.event_tree if evt.run in f2mev_dict], dtype=data_dtype)

    # Create an array of f2mev values that match with the order of events in "data"
    f2mev = np.zeros(data.size)

    prev_data_run = 9999
    count_run = 0
    count_subrun = 0
    for i, data_run in enumerate(data['run']) :
        if data_run == prev_data_run :
            count_subrun += 1
        elif i > 0 :
            count_run += 1

        f2mev[count_run+count_subrun] = f2mev_dict[data_run]
        prev_data_run = data_run

    # Check that all the runs have a flux to MeV value
    if 0.0 in f2mev :
        indicies = np.where( f2mev == 0.0 )[0]
        print("Run(s) that have no flux to MeV conversion recorded: ")
        print(indicies)

    apply_chandong_correction(data['run'], data['flux'])

    data['energy'] = data['flux'] / f2mev
    data = data[data['energy'] > 20.0] # Apply lower energy threshold
    #x_data, y_data, z_data, E_data = [[evt.x, evt.y, evt.z, evt.flux] for evt in tree]
    data_file.Close()

    # convert data coordinates from meters to milli-meters
    #x_data, y_data, z_data = [[x*1000.0 for x in arr] for arr in [x_data, y_data, z_data]]
    data["pos"] *= 1e3 # convert data coordinates to milli-meters

    MC_file = ROOT.TFile(mc_filename, "READ")
    MC_tree = MC_file.energy_tree
    x_MC, y_MC, z_MC, E_MC, ptype_MC = zip(*[[evt.x, evt.y, evt.z, evt.E, evt.PDG] for evt in MC_tree])
    MC_file.Close()
    print("MC data loaded")

    # Calcuate the MC weights for the mu-plus and mu-minus proportion
    plus_count = len([x for x in ptype_MC if x==13])
    minus_count = len([x for x in ptype_MC if x==-13])
    plus_weight = float(plus_count)/(plus_count + minus_count)
    minus_weight = float(minus_count)/(plus_count + minus_count)
    # Ratio mu+/mu- = 1.21 from PRD 74,082006
    plus_weight *= 1.21
    minus_weight *= 1.
    print("Mu-minus weight = %f" % minus_weight)
    print("Mu-plus weight = %f" % plus_weight)

    # Let's create an array of the reco michel spectra
    name_fmt = "energy_data_%i_%i"

    for E, pos in data[['energy', 'pos']]:
        x,y,z = pos
        Rs_reco = x**2+y**2
        #print("DEBUG - x=%f y=%f z=%f r=%f" % (x,y,z,Rs_reco))
        if Rs_reco < R_FV**2 and Z_FV_LOW < z < Z_FV_HIGH:
            location_data_hist.Fill(Rs_reco, z, E)
        #else:
           # print("DEBUG - could not fill location_data_hist")


    inner_E_data = [E for E,(x,y,z) in data[["energy", "pos"]]\
                      if np.sqrt(x**2 + y**2) < R_FV and Z_FV_LOW < z < Z_FV_HIGH]
    energy_inner_data_hist = ROOT.TH1D("energy_inner_data", "", n_bin, min_E, max_E)
    energy_inner_data_hist.Sumw2()
    [energy_inner_data_hist.Fill(v) for v in inner_E_data]

    ########## Perform Fit of the Inner Energy Hist for Real Data ##########
    fit_vals_data, fit_err_data, inner_E_fit_data_graph = fit.do_michel_fit(energy_inner_data_hist)
    scale_data, R_data, ampl_data = fit_vals_data
    scale_data_err, res_data_err, ampl_data_err = fit_err_data
    print("\nInner_Fit_E_Data Scale %0.5f ± %0.5f" % (scale_data, scale_data_err))
    print("Inner_Fit_E_Data Res %0.5f ± %0.5f" % (R_data, res_data_err))
    print("Inner_Fit_E_Data Ampl %0.5f ± %0.5f\n" % (ampl_data, ampl_data_err))

    ########## Perform Fit of the Inner Energy Hist for MC Data ##########
    location_mc_hist = fit.PositionBasedHistograms("MC", "", 12, 0, 2.56e6, 10, -1.25e3, 1.25e3, n_bin, min_E, max_E)
    hist_MC = ROOT.TH1D("energy_inner_MC", "MC Energy", n_bin, min_E, max_E)
    [hist_MC.Fill(E, plus_weight if ptype==1 else minus_weight) for E,x,y,z,ptype in zip(E_MC, x_MC, y_MC, z_MC, ptype_MC)\
                     if np.sqrt(x**2 + y**2) < R_FV and Z_FV_LOW < z < Z_FV_HIGH and E > 20.0]

    [location_mc_hist.Fill(x**2+y**2, z, E, plus_weight if ptype==1 else minus_weight) for E,x,y,z,ptype in zip(E_MC, x_MC, y_MC, z_MC, ptype_MC)\
                     if np.sqrt(x**2 + y**2) < R_FV and Z_FV_LOW < z < Z_FV_HIGH and E > 20.0]
    fit_vals_mc, fit_err_mc, inner_E_fit_MC_graph = fit.do_michel_fit(hist_MC)
    scale_MC, R_MC, ampl_MC = fit_vals_mc
    scale_MC_err, res_MC_err, ampl_MC_err = fit_err_mc
    print("\nInner_Fit_E_MC Scale %0.5f ± %0.5f" % (scale_MC, scale_MC_err))
    print("Inner_Fit_E_MC Res %0.5f ± %0.5f" % (R_MC, res_MC_err))
    print("Inner_Fit_E_MC Ampl %0.5f ± %0.5f\n" % (ampl_MC, ampl_MC_err))


    # Compare MC and data energy scale
    scale_ratio = scale_MC/scale_data
    scale_ratio_err = scale_ratio*np.sqrt( (scale_MC_err/scale_MC)**2 + (scale_data_err/scale_data)**2 )
    print("#################### CENTRAL VALUE SCALING/SMEARING ####################")
    print("R_data: %0.5f ± %0.5f" % (R_data, res_data_err))
    print("R_MC: %0.5f ± %0.5f" % (R_MC, res_MC_err))
    print("Scale_ratio: %0.5f ± %0.5f" % (scale_ratio, scale_ratio_err))
    print("########################################################################")

    location_mc_hist.Fit()
    location_data_hist.Fit()
    data_scale = np.zeros((location_data_hist.scale_h2.GetNbinsX(), location_data_hist.scale_h2.GetNbinsY()))
    mc_scale = np.zeros((location_data_hist.scale_h2.GetNbinsX(), location_data_hist.scale_h2.GetNbinsY()))
    corr = np.zeros((location_data_hist.scale_h2.GetNbinsX(), location_data_hist.scale_h2.GetNbinsY()))
    corr_h2 = ROOT.TH2D("mc_energy_scale_correction", "", 12, 0, 2.56e3, 10, -1.25e3, 1.25e3)
    corr_avg = []
    for i in range(location_data_hist.scale_h2.GetNbinsX()):
        for j in range(location_data_hist.scale_h2.GetNbinsY()):
            d = location_data_hist.scale_h2.GetBinContent(i+1, j+1)
            mc = location_mc_hist.scale_h2.GetBinContent(i+1, j+1)
            data_scale[i,j] = d
            mc_scale[i,j] = mc
            corr[i, j] = d/mc if mc !=0 else 1.0
            corr_h2.SetBinContent(i+1, j+1, corr[i,j])

            if mc != 0:
                corr_avg.append(d/mc)

    corr_avg = np.array(corr_avg)
    print("Average scale factor: "+np.average(corr_avg))

    # Now create an array of the Smeared MC michel spectra
    # Total FV energy histogram
    smeared_MC = ROOT.TH1D("smeared_MC", "smeared_MC", n_bin, min_E, max_E)
    smeared_MC.Sumw2()

    smeared_data_out = []
    for E,x,y,z,ptype in zip(E_MC, x_MC, y_MC, z_MC, ptype_MC):
        R = np.sqrt(x**2 + y**2)
        Rs_reco = R**2
        #E *= scale_ratio
        weight = plus_weight if ptype==1 else minus_weight
        sigma_smear = delta_sigma(E, R_data, R_MC)
        N_TRIALS = 1
        E *= corr_h2.GetBinContent(corr_h2.FindBin(Rs_reco, z))
        smeared_energies = np.random.normal(loc=E, scale=sigma_smear, size=N_TRIALS) # do the smearing
        for E in smeared_energies:
            smeared_data_out.append([E,x,y,z,weight/N_TRIALS])
            if R < R_FV and Z_FV_LOW < z < Z_FV_HIGH and E >= 20:
                smeared_MC.Fill(E, weight/N_TRIALS)

    # Fit the MC Smeared Data
    fit_vals_s, fit_err_s, smeared_fit_graph = fit.do_michel_fit(smeared_MC)
    scale_smeared, res_smeared, _ = fit_vals_s
    scale_smeared_err, res_smeared_err, _ = fit_err_s


    ########## MAKE NICE PLOTS AND WRITE DATA TO A FILE ##################
    outFile = ROOT.TFile(output_filename, "RECREATE")
    outFile.cd()

    ntuple_out = ROOT.TNtuple("energy_tree", "energy_tree", "E:x:y:z:weight")
    [ntuple_out.Fill(E,x,y,z,w) for E,x,y,z,w in smeared_data_out]
    ntuple_out.Write()

    # Put Draw setting stuff
    #corr_h2.Draw("TEXT COLZ")
    corr_h2.Write()

    #location_data_hist.Draw("TEXT COLZ")
    #location_data_hist.Write()

    #location_mc_hist.Draw("TEXT COLZ")
    #location_mc_hist.Write()

    c0 = ROOT.TCanvas()
    hist_MC.Draw("HIST")
    hist_MC.SetStats(0)
    inner_E_fit_MC_graph.SetLineColor(2)
    inner_E_fit_MC_graph.SetLineWidth(2)
    inner_E_fit_MC_graph.Draw("sameL")
    leg0 = ROOT.TLegend(0.60, 0.65, 0.88, 0.86)
    leg0.AddEntry(hist_MC, "MC Data", "l")
    leg0.AddEntry(inner_E_fit_MC_graph, "Michel Fit", "l")
    leg0.AddEntry("", "Scale: %0.5f #pm %0.5f" %(scale_MC, scale_smeared_err), "")
    leg0.AddEntry("", "Ep Res: %0.3f%% #pm %0.3f" % (R_MC*100, res_MC_err*100), "")
    leg0.Draw("same")

    
    C1 = ROOT.TCanvas("Inner Reco Energy: Real Data")
    energy_inner_data_hist.SetTitle("Data")
    energy_inner_data_hist.GetXaxis().SetTitle("Reconstructed Energy [MeV]")
    energy_inner_data_hist.GetYaxis().SetTitle("Events/%0.2f MeV" % bwidth)
    energy_inner_data_hist.GetYaxis().SetTitleOffset(1.3)
    energy_inner_data_hist.SetTitle("")
    energy_inner_data_hist.SetStats(0)
    energy_inner_data_hist.SetLineColor(ROOT.kBlack)
    energy_inner_data_hist.SetLineWidth(2)
    energy_inner_data_hist.Draw("HISTE")
    inner_E_fit_data_graph.SetLineColor(2)
    inner_E_fit_data_graph.SetLineWidth(2)
    inner_E_fit_data_graph.Draw("Same")
    leg1 = ROOT.TLegend(0.60, 0.65, 0.88, 0.86)
    leg1.AddEntry(energy_inner_data_hist, "Real Data", "l")
    leg1.AddEntry(inner_E_fit_data_graph, "Michel Fit", "l")
    leg1.AddEntry("", "Scale: %0.5f #pm %0.5f" %(scale_data, scale_data_err), "")
    leg1.AddEntry("", "Ep Res: %0.3f%% #pm %0.3f" % (R_data*100, res_data_err*100), "")
    leg1.Draw("Same")
    C1.Write()

    C2 = ROOT.TCanvas("Smeared MC")
    smeared_MC.GetXaxis().SetTitle("Reconstructed Energy [MeV]")
    smeared_MC.GetYaxis().SetTitle("Events/"+str(bwidth)+" MeV")
    smeared_MC.GetYaxis().SetTitleOffset(1.3)
    smeared_MC.SetTitle("")
    smeared_MC.SetStats(0)
    smeared_MC.SetLineColor(ROOT.kBlue)
    smeared_MC.SetLineWidth(2)
    smeared_MC.Draw("HISTE")
    smeared_fit_graph.SetLineColor(2)
    smeared_fit_graph.SetLineWidth(2)
    smeared_fit_graph.Draw("Same")
    leg2 = ROOT.TLegend(0.60, 0.65, 0.88, 0.86)
    leg2.AddEntry(smeared_MC, "Smeared MC", "l")
    leg2.AddEntry(smeared_fit_graph, "Michel Fit", "l")
    leg2.AddEntry("", "Scale: %0.5f #pm %0.5f" % (scale_smeared, scale_smeared_err), "")
    leg2.AddEntry("", "Ep Res: %0.5f%% #pm %0.5f" % (res_smeared*100, res_smeared_err*100), "")
    leg2.Draw("Same")
    C2.Write()

    C5 = ROOT.TCanvas("Normalized Comparison")
    norm_MC = hist_MC.Clone("norm_MC")
    norm_MC.SetLineColor(ROOT.kBlue)
    norm_data = energy_inner_data_hist.Clone("norm_data")
    norm_smear = smeared_MC.Clone("norm_smear")
    norm_MC.Scale(1.0/norm_MC.Integral())
    norm_data.Scale(1.0/norm_data.Integral())
    norm_smear.Scale(1.0/norm_smear.Integral())
    norm_MC.GetYaxis().SetTitle("Normalized")
    norm_MC.GetYaxis().SetTitleOffset(1.4)
    norm_MC.Draw("HIST")
    norm_data.SetLineColor(1)
    norm_data.Draw("HIST Same")
    norm_smear.SetLineColor(ROOT.kGreen+1)
    norm_smear.Draw("HIST Same")

    leg = ROOT.TLegend(0.12, 0.65, 0.4, 0.86)
    leg.AddEntry(norm_MC, "MC")
    leg.AddEntry(norm_data, "Data")
    leg.AddEntry(norm_smear, "Smeared MC")
    leg.Draw("Same")
    C5.Write()

    c99  = ROOT.TCanvas()
    c9  = ROOT.TCanvas()
    location_mc_hist.DrawFit(c9, c99)
    c9.Update()
    c99.Update()
    c9.Write()
    c99.Write()
    c88, c8 = ROOT.TCanvas(), ROOT.TCanvas()
    location_data_hist.DrawFit(c8, c88)
    c8.Update()
    c88.Update()
    c8.Write()
    c88.Write()

