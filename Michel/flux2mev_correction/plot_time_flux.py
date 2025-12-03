"""
	Originally from McCaffrey:
	/home/marzece/KDAR_Analysis/MichelAnalysis/f2mev_correction/plot_time_flux.py

	Takes the michel fit values for the flux to MeV conversion from correction_values.json
	and plots the the flux to MeV over time. Data taking periods are automatically assigned
	based on the dates associated with the runs and a user defined length of time between periods.

	The input conversion values from correction_values.json are compared to those from HyoungKu Jeon, Eric Marzec, and Taku Dodo (all 2021)
	as well as a single user defined version from the json file.

"""

# coding: utf-8
import ROOT
import numpy as np
from datetime import datetime, timedelta
from collections import defaultdict
import os
import time
import csv
import json
import math
import pdb
import re
import itertools
import pandas as pd

JSONFILE = "/home/littleca/kdar/correction_values.json"
OUTFILE_PATH = "/home/littleca/kdar/Michel/flux2mev_correction/output_fluxCorr/"
WORK_PATH= "/home/littleca/kdar/Michel/flux2mev_correction/"

FIRST_RUN = 1597
ALPHA = 0.35
TREND_WIDTH = 2.0*24*3600
DOWNTIME_BTWN_PERIODS = 5.26e6    # Seconds

# From https://github.com/JSNS2/AnalysisTools/blob/main/EnergyCorrection/20240312/TimeCorrection_2021.txt
dodo_correction = np.loadtxt("/home/littleca/kdar/Michel/flux2mev_correction/Dodo_correction.txt")

def gaussian(x,mu=0.0,sigma=1.0):
    return np.exp(-0.5*((x-mu)/sigma)**2)/(sigma*np.sqrt(2*np.pi))

def get_trend_line_simple(points):
    """ Input argument points should be a (N, 3) array where the columns are
        x, y, and error values. """
    sigma = TREND_WIDTH
    xax = np.linspace(np.min(points[:,0]), np.max(points[:, 0]), 1000)
    peak_val = gaussian(0,0,sigma=sigma)
    results = np.zeros((len(xax), 3))
    results[:, 0] = xax

    point_err_weight = 1.0/points[:, 2]**2
    # If there are large periods of time between runs, the time_weight will be nonsense
    # So we will calculate the trend line for each seperate run period
    time_weights = np.array([gaussian(points[:,0], mu=x0, sigma=sigma) for x0 in xax])

    mean_vals = np.array([np.average(points[:, 1], weights=w*point_err_weight) for w in time_weights])
    results[:, 1] = mean_vals

    variance = np.array([np.average((points[:,1]-mean)**2,weights=w*point_err_weight)  for w, mean in zip(time_weights,mean_vals)])
    std_dev = np.sqrt(variance)

    weighted_counts = time_weights/peak_val
    counts = np.sum(weighted_counts, axis=1)
    mean_err = std_dev/np.sqrt(counts)
    results[:, 2] = mean_err
    return results

def get_trend_line(points):
    sigma = trend_width
    peak_val = gaussian(0,0,sigma=sigma)
    xax = np.linspace(points[0,0], points[-1, 0], 1000)
    result = np.zeros((len(xax), 3))
    result[:, 0] = xax
    point_vals = points[:, 1]
    point_err = points[:, 2]
    point_err_weight = 1.0/point_err**2

    weights = np.array([gaussian(points[:, 0], mu=x0, sigma=sigma) for x0 in xax])
    weighted_counts = weights/peak_val

    #mean_vals = [point_vals[ww>1e-20]*ww[ww>1e-20]*(1.0/2*(point_err[ww>1e-20]*ww[ww>1e-20])**2) / np.sum(ww[ww>1e-20]/2*(point_err[ww>1e-20]*ww[ww>1e-20])**2) for ww in weights]
    mean_vals = []
    width_vals = []
    for ww in weights:
        ww = ww/np.sum(ww)
        mask = ww > 1e-80
        ww = ww[mask]
        vals = point_vals[mask]
        errs = point_err[mask]
        mu = np.sum(vals*ww/(2*(errs)**2))/(np.sum(ww/(2*(errs)**2)))
        width = np.sqrt(1.0/np.sum(ww/errs**2))
        width_vals.append(width)
        #mu = np.sum(vals*ww)
        mean_vals.append(mu)


    #mean_vals2 = np.array([np.average(points[:, 1], weights=w*point_err_weight) for w in weights])
    result[:, 1] = np.array(mean_vals)

    #result[:, 1] = np.array([np.sum(points[:, 1]*point_err**2*w)/(np.sum(point_err**2)*np.sum(w))      for w in weights])

    #error_sum = 1.0/np.sqrt(np.array([np.sum((1.0/point_err**2)*w/np.sum(w)) for w in weights]))
    error_sum = np.array(width_vals)

    variance = np.array([np.average((points[:,1]-mean)**2,weights=w*point_err_weight)  for w, mean in zip(weights,result[:,1])])
    std_dev = np.sqrt(variance)
    counts = np.sum(weighted_counts, axis=1)
    mean_err = std_dev/np.sqrt(counts)

    err = np.sqrt(std_dev**2 + 0*error_sum**2)
    result[:, 2] = err
    return result

def get_run_dates():
    vals = list(csv.reader(open(WORK_PATH+"run_dates.txt")))
    run_times = defaultdict(lambda : (None, None))
    for i, (run, start_date, start_time, end_date, end_time) in enumerate(vals):
        run = int(run)
        if start_date:
            start_date = datetime.strptime(start_date.strip(), "%Y/%m/%d")
            if start_time:
                start_time = datetime.strptime(start_time.strip(), "%H:%M:%S")
                start_date = datetime(start_date.year, start_date.month, start_date.day,
                                      start_time.hour, start_date.minute, start_time.second)
        else:
            start_date = None

        if end_date:
            end_date = datetime.strptime(end_date.strip(), "%Y/%m/%d")
            if end_time:
                end_time = datetime.strptime(end_time.strip(), "%H:%M:%S")
                end_date = datetime(end_date.year, end_date.month, end_date.day,
                                      end_time.hour, end_date.minute, end_time.second)
        else:
            end_date = None
        if(start_date):
            run_times[run] = (start_date, end_date)
    return run_times

def get_subrun_times():
    subrun_data = [[int(v) for v in row] for row in csv.reader(open("RUN_TIMESTAMP_DATA.txt"))]
    subrun_dict = defaultdict(dict)
    for rn, srn, ts in subrun_data:
        subrun_dict[rn].update({srn:ts*8e-9})
    return subrun_dict

def get_time(run, subrun):
    get_time.run_dates = get_run_dates()
    get_time.subrun_times = get_subrun_times()

    run_date = get_time.run_dates[run]
    if run_date is None or run_date[0] is None:
        return None
    try:
        subrun_time = get_time.subrun_times[run][subrun]
    except KeyError as e:
        print(run ,subrun)
        available_srns = get_time.subrun_times[run].keys()
        if not available_srns:
            return run_date[0]
        subrun = max(available_srns)
        subrun_time = get_time.subrun_times[run][subrun]
    return run_date[0] + timedelta(seconds=subrun_time)

def time_to_rn(dt):
    run_dates = get_run_dates()
    subrun_times = get_subrun_times()
    result = []
    for rn, run_start_time, in sorted(run_dates.items(), key=lambda x: x[0]):
        if not subrun_times[rn]:
            continue
        if not run_start_time or not run_start_time[0]:
            continue
        run_start_time = run_start_time[0]
        max_srn = max(subrun_times.keys())
        for subrun, subrun_dt in sorted(subrun_times[rn].items(), key=lambda x: x[0]):
            this_time = run_start_time + timedelta(seconds=subrun_dt)
            result.append((rn, subrun, max_srn, this_time))
    return result


def get_hyoungkus_correction(run_times):
    path = "/home/littleca/kdar/Michel/flux2mev_correction/HyoungKu_correction/"
    hk_vals = []
    for filename in os.listdir(path):
        if re.match(r'Run_(\d+)_timeCorr.dat', filename):
            fin = open(os.path.join(path, filename), 'r')
            lines = fin.readlines()
            val_list = [float(x) for x in lines[0].strip().split('\t')]
            run = val_list[0]

            timestamp = (time.mktime(run_times[run][0].timetuple()) - ROOT.gStyle.GetTimeOffset())
            # Run #, timestamp, flux to MeV energy scale, scale error
            hk_vals.append([run, timestamp, val_list[2], val_list[3]])
    hk_vals = np.array(hk_vals)
    return hk_vals

def get_json_correction(versionID, run_times):
    # Read in data from json file
    with open(JSONFILE, "r") as f:
        data_vals = json.load(f)

    vals = np.zeros((len(data_vals[versionID]), 4))
    for i, (entry_type, value_dict) in enumerate(data_vals[versionID].items()):
        if not isinstance(value_dict["Run"],int):
            continue
        if run_times[value_dict["Run"]][0] is None:
            continue
        timestamp = (time.mktime(run_times[value_dict["Run"]][0].timetuple()) - ROOT.gStyle.GetTimeOffset())
        # Val columns are run #, timestamp, energy scale, scale error
        vals[i,:] = value_dict["Run"], timestamp, value_dict["Flux2MeV"], value_dict["Flux2MeVErr"]

    # Trim the vals array so there are no empty rows (due to not having a date)
    mask = vals == 0
    countVals = (~mask).sum(axis=1)
    zeroIndex = np.where(countVals==0)[0]
    vals = np.delete(vals, zeroIndex, 0)

    vals_per_period = get_period_values(vals)

    # Find the nominal flux ot MeV (the weighted average value) and it's error for all run periods and for each individual period
    # The weight is 1/flux to mev error
    # We will normalize our results later with respect to this average
    all_vals = list(itertools.chain.from_iterable(vals_per_period.values()))
    all_vals = np.array(all_vals)

    all_f2mev = np.average(all_vals[:,2], weights=list(map(lambda x: 1.0/x, all_vals[:,3])))
    all_f2mev_err = np.average(all_vals[:,3], weights=list(map(lambda x: 1.0/x, all_vals[:,3])))

    nominal_f2mev_dict = defaultdict(lambda: [])
    nominal_f2mev_dict['all'] = np.array([all_f2mev, all_f2mev_err])

    for p_i in vals_per_period:
        n_f2mev = np.average(vals_per_period[p_i][:,2], weights=list(map(lambda x: 1.0/x, vals_per_period[p_i][:,3])))
        n_f2mev_err = np.average(vals_per_period[p_i][:,3], weights=list(map(lambda x: 1.0/x, vals_per_period[p_i][:,3])))

        nominal_f2mev_dict[p_i] = np.array([n_f2mev, n_f2mev_err])

    return vals_per_period, nominal_f2mev_dict

def get_period_values(vals):
    # Group the run values into periods so we know when there is a long time between runs
    # We set what the minimum time between run periods is in seconds

    vals_t_s = vals[vals[:,1].argsort()]
    prev_start = -1*DOWNTIME_BTWN_PERIODS
    period = -1
    vals_per_period = defaultdict(lambda : [])
    for i, (r, start, v, err) in enumerate(vals_t_s):
        if (start - prev_start) >= DOWNTIME_BTWN_PERIODS:
            period += 1
        prev_start = start

        vals_per_period[period].append(np.array([r, start, v, err]))

    period = period+1
    print("Found %i run periods." % (period))
    period_ends = []
    for i in range(period):
        vals_per_period[i] = np.array(vals_per_period[i])
        period_ends.append(vals_per_period[i][:,1].max())
    period_ends.sort()
    print([time.strftime('%Y-%m-%d', time.gmtime(i + ROOT.gStyle.GetTimeOffset())) for i in period_ends])

    return vals_per_period

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Plot flux to MeV corrections over time given from the corrctions JSON file and compare to other correction values, included (Dodo's corrections, HyoungKu's corrections) or optionally input.")
    parser.add_argument("version", type=str, help="version ID from correction values' top level key")
    parser.add_argument("--michelEndpoint", type=float, help="Michel endpoint value [MeV] for this versions data")
    parser.add_argument("--compare", type=str, help="version ID to compare input version to")
    parser.add_argument("--compareMichelEndpoint", type=float, help="Michel endpoint value [MeV] for the the comparison version data")
    args = parser.parse_args()

    versionID = args.version
    MICHEL_ENDPOINT = args.michelEndpoint
    compareID = args.compare
    compare_MICHEL_ENDPOINT = args.compareMichelEndpoint

    # Set the timing offset to the first run to normalize the rest of the times
    # We've hardcoded the first run to apply for all the data sets
    run_times = get_run_dates()
    ROOT.gStyle.SetTimeOffset(time.mktime(run_times[FIRST_RUN][0].timetuple()))
    first_time = (time.mktime(run_times[FIRST_RUN][0].timetuple()) - ROOT.gStyle.GetTimeOffset())
    last_time = (time.mktime(run_times[max(run_times.keys())][0].timetuple()) - ROOT.gStyle.GetTimeOffset())

    # Histograms for comparing the current version's flux to MeV with our other data sets
    eric_comp_h = ROOT.TH1D("Fraction of Eric's Data","Fraction of Eric's Data",500,-1,2)
    dodo_comp_h = ROOT.TH1D("Fraction of Dodo's Data","Fractio of Dodo's Data",500,-1,2)
    hk_comp_h = ROOT.TH1D("Fraction of HK's Data","Fraction of HK'd Data",500,-1,2)
    comp_h = ROOT.TH1D("Fraction of Version %s's Data" % (compareID),"Fraction of Version %s's Data" % (compareID),500,-1,2)

    # Histograms for comparing the current version's error with our other data sets
    eric_comp_err_h = ROOT.TH1D("Difference with Eric's Data","Fraction of Eric's Data",500,-100,100)
    hk_comp_err_h = ROOT.TH1D("Difference with HK's Data","Fraction of HK'd Data",500,-100,100)
    comp_err_h = ROOT.TH1D("Difference with Version %s's Data" % (compareID),"Fraction of Version %s's Data" % (compareID),500,-100,100)

    # Get numerical comparison between the different corrections
    # HyoungKu, Dodo, and Eric's data uses a michel endpoint = 52.8 MeV
    endpoint_scale = 52.8/MICHEL_ENDPOINT
    if args.compare:
        compare_endpoint_scale = compare_MICHEL_ENDPOINT/MICHEL_ENDPOINT

    # Fill TGraph data for various data sets
    gabs = []						# Array of TGraphErrors of current version data
    gmine_indiv = []					# Array of TGraphErros of current version data (normalized wrt average flux to MeV per period)
    gmine_all = []					# Array of TGraphErros of current version data (normalized wrt average flux to MeV for all periods)
    ghk = ROOT.TGraphErrors()				# HyoungKu's data (normalized wrt average flux to MeV)
    gdodo = ROOT.TGraph()				# Dodo's data (normalized wrt first run's flux to MeV)
    geric = ROOT.TGraphErrors()				# Eric's data (normalized wrt average flux to MeV)

    # ===== Current Version Data =====
    vals_per_period, nominal_flux2mev_dict = get_json_correction(versionID, run_times)

    for p_i, p_v in vals_per_period.items():
        gmine_indiv.append(ROOT.TGraphErrors())
        gmine_all.append(ROOT.TGraphErrors())
        gabs.append(ROOT.TGraphErrors())
        for i, (r, start, v, err) in enumerate(p_v):
            gmine_indiv[p_i].SetPoint(i, start, v/nominal_flux2mev_dict[p_i][0])
            gmine_indiv[p_i].SetPointError(i,0, err/nominal_flux2mev_dict[p_i][0])

            gmine_all[p_i].SetPoint(i, start, v/nominal_flux2mev_dict['all'][0])
            gmine_all[p_i].SetPointError(i,0, err/nominal_flux2mev_dict['all'][0])

            gabs[p_i].SetPoint(i, start, v)
            gabs[p_i].SetPointError(i,0, err)

    # ===== HyoungKu's Data =====
    hk_vals = get_hyoungkus_correction(run_times)

    # To normalize HK's data with respect to the average flux to mev value
    # Change hk_AV_FLUX2MEV = 1 if we want to normalize with respect to HKs first run
    hk_AVG_FLUX2MEV = np.average(hk_vals[:,2], weights=list(map(lambda x: 1.0/x, hk_vals[:,3])))
    #hk_AVG_FLUX2MEV_ERR = np.average(hk_vals[:,3], weights=list(map(lambda x: 1.0/x, hk_vals[:,3])))

    for i, (r, start, v, err) in enumerate(hk_vals):
        ghk.SetPoint(i, start, v/hk_AVG_FLUX2MEV)
        ghk.SetPointError(i, 0, err/hk_AVG_FLUX2MEV)

        # Make sure we're getting the correct run data from vals_per_period
        run_idx = np.where(vals_per_period[0][:,0] == r)[0]
        if len(run_idx) == 1 :
            hk_comp_h.Fill((vals_per_period[0][run_idx[0], 2]/nominal_flux2mev_dict[0][0])/(v/hk_AVG_FLUX2MEV)*endpoint_scale)
            hk_comp_err_h.Fill(vals_per_period[0][run_idx[0], 3] - err)

    # ===== Dodo's Data =====
    # Note: For comparisons sake, all of the data is normalized like the json/current data. Dodo's corrections have already been normalized, most likely with respect to the flux to MeV from the first run.
    # So comparisons against Dodo's data may not be super accurate.
    # Clean up the data from Dodo at the start of the code
    dodo_correction = [(get_time(r, srn), v, r) for r, srn, _, _, v in dodo_correction]
    dodo_correction_withRun = {time.mktime(t.timetuple()) - ROOT.gStyle.GetTimeOffset(): (v, r) for t, v, r in dodo_correction if t is not None}	# So I can have the run # to use with comparison histogram
    dodo_correction_forPlots = {time.mktime(t.timetuple()) - ROOT.gStyle.GetTimeOffset(): v for t, v, r in dodo_correction if t is not None}

    for i, (d, vr_tup) in enumerate(sorted(dodo_correction_withRun.items(), key=lambda x: x[0])):
        gdodo.SetPoint(i, d, vr_tup[0])

		# Find corresponding runs flux to MeV value for current version
        # Make sure we're getting the correct run data from vals_per_period
        run_idx = np.where(vals_per_period[0][:,0] == vr_tup[1])[0]
        if len(run_idx) == 1 :
            dodo_comp_h.Fill((vals_per_period[0][run_idx[0], 2]/nominal_flux2mev_dict[0][0])/vr_tup[0]*endpoint_scale)
 
    # ===== Eric's Data =====
    # Read in Eric's correction values to compare against
    eric_data = []
    with open("/home/littleca/kdar/Michel/flux2mev_correction/Eric_correction/temp.txt", 'r') as efile:
        csvFile = csv.reader(efile)
        for line in csvFile:
            eric_data.append([float(i) for i in line])

    for i, run_data in enumerate(eric_data):
        eric_timestamp = (time.mktime(run_times[run_data[0]][0].timetuple()) - ROOT.gStyle.GetTimeOffset())
        eric_data[i].insert(1, eric_timestamp)
        # Run	Time	Flux2MeV	Err

    # Get the average flux to MeV value 
    eric_data_list = list(zip(*eric_data))
    eric_AVG_FLUX2MEV = np.average(list(eric_data_list[2]), weights=list(map(lambda x: 1.0/x, list(eric_data_list[3]))))
    #eric_AVG_FLUX2MEV_ERR = np.average(list(eric_data_list[3]), weights=list(map(lambda x: 1.0/x, list(eric_data_list[3]))))

    for i, (r, start, v, err) in enumerate(eric_data):
        geric.SetPoint(i, start, v/eric_AVG_FLUX2MEV)
        geric.SetPointError(i,0, err/eric_AVG_FLUX2MEV)

		# Find corresponding runs flux to MeV value for current version
        # Make sure we're getting the correct run data from vals_per_period
        run_idx = np.where(vals_per_period[0][:,0] == r)[0]
        if len(run_idx) == 1 :
            eric_frac = vals_per_period[0][run_idx[0], 2]/v*endpoint_scale
            eric_comp_h.Fill(eric_frac)
            eric_comp_err_h.Fill(vals_per_period[0][run_idx[0], 3] - err)

    # ===== Comparison Data =====
    # Read in corrections for optional JSON version to compare against
    if args.compare :
        print("Comparing to : " + compareID)
        comp_vals, comp_nominal_flux2mev_dict = get_json_correction(compareID, run_times)

        gcomp_indiv = []
        for p_i, p_v in comp_vals.items():
            gcomp_indiv.append(ROOT.TGraphErrors())
            for i, (r, start, v, err) in enumerate(p_v):
                gcomp_indiv[p_i].SetPoint(i, start, v/comp_nominal_flux2mev_dict[p_i][0])
                gcomp_indiv[p_i].SetPointError(i,0, err/comp_nominal_flux2mev_dict[p_i][0])

		    # Find corresponding runs flux to MeV value for current version
                # Make sure we're getting the correct run data from vals_per_period
                run_idx = np.where(vals_per_period[0][:,0] == r)[0]
                if len(run_idx) == 1 :
                    comp_h.Fill(vals_per_period[0][run_idx[0], 2]/v*compare_endpoint_scale)
                    comp_err_h.Fill(vals_per_period[0][run_idx[0], 3] - err)

    OUTFILE = OUTFILE_PATH+"/"+str(versionID)+"/time_flux_plots_"+str(versionID)+".root" 
    outFile = ROOT.TFile(OUTFILE, "RECREATE")
 
    eric_comp_h.Write()
    dodo_comp_h.Write()
    hk_comp_h.Write()
    comp_h.Write()

    eric_comp_err_h.Write()
    dodo_comp_err_h.Write()
    hk_comp_err_h.Write()
    comp_err_h.Write()

    # Get the trend line of current version correction data for each period
    gtrend_indiv = []
    gtrend_all = []
    gtrend_abs = []
    for p_i, p_v in vals_per_period.items():
        gtrend_indiv.append(ROOT.TGraphErrors())
        gtrend_all.append(ROOT.TGraphErrors())
        gtrend_abs.append(ROOT.TGraphErrors())

        trend = get_trend_line_simple(p_v[:, 1:4])
        for i, (x, y, err) in enumerate(trend):
            gtrend_abs[p_i].SetPoint(i,x,y)
            gtrend_abs[p_i].SetPointError(i,0,err)
            gtrend_indiv[p_i].SetPoint(i,x,y/nominal_flux2mev_dict[p_i][0])
            gtrend_indiv[p_i].SetPointError(i,0,err/nominal_flux2mev_dict[p_i][0])
            gtrend_all[p_i].SetPoint(i,x,y/nominal_flux2mev_dict['all'][0])
            gtrend_all[p_i].SetPointError(i,0,err/nominal_flux2mev_dict['all'][0])
        gtrend_abs[p_i].SetLineColor(ROOT.kBlue)
        gtrend_abs[p_i].SetFillColorAlpha(ROOT.kBlue, ALPHA)
        gtrend_abs[p_i].SetFillStyle(3001)
        gtrend_indiv[p_i].SetLineColor(ROOT.kBlue)
        gtrend_indiv[p_i].SetFillColorAlpha(ROOT.kBlue, ALPHA)
        gtrend_indiv[p_i].SetFillStyle(3001)
        gtrend_all[p_i].SetLineColor(ROOT.kViolet)
        gtrend_all[p_i].SetFillColorAlpha(ROOT.kViolet, ALPHA)
        gtrend_all[p_i].SetFillStyle(3001)

    # Plot the current version corrections flux to MeV and its trend line
    c1 = ROOT.TCanvas("%s F2MeV Absolute" % (versionID))
    mg0 = ROOT.TMultiGraph()
    for p_i in vals_per_period:
        mg0.Add(gabs[p_i])
        mg0.Add(gtrend_abs[p_i], "lineE3")
        mg0.Add(gtrend_all[p_i], "lineE3")
    mg0.Draw("APE*")
    mg0.GetXaxis().SetTimeDisplay(1)
    mg0.GetYaxis().SetTitle("Michel Flux to MeV")
    mg0.GetXaxis().SetTimeFormat("%Y/%m/%d")
    leg0 = ROOT.TLegend(0.1, 0.66, 0.45, 0.9)
    leg0.AddEntry(gabs[0], "%s" % (versionID), "PE")
    leg0.Draw('same')
    #c1.Update()
    c1.Write()

    # Plot the normalized conversion with HyoungKu's results
    gnominal_indiv = []
    gnominal_all = []
    for p_i in vals_per_period:
        gnominal_indiv.append(ROOT.TGraphErrors())
        gnominal_indiv[p_i].SetPoint(0, first_time, 1.0)
        gnominal_indiv[p_i].SetPoint(1, last_time, 1.0)
        gnominal_indiv[p_i].SetPointError(0, 0, nominal_flux2mev_dict[p_i][1]/nominal_flux2mev_dict[p_i][0])
        gnominal_indiv[p_i].SetPointError(1, 0, nominal_flux2mev_dict[p_i][1]/nominal_flux2mev_dict[p_i][0])
        gnominal_indiv[p_i].SetFillColorAlpha(p_i+5, ALPHA)
        gnominal_indiv[p_i].SetFillStyle(3001)
        gnominal_indiv[p_i].SetLineColor(p_i+5)

        gnominal_all.append(ROOT.TGraphErrors())
        gnominal_all[p_i].SetPoint(0, first_time, 1.0)
        gnominal_all[p_i].SetPoint(1, last_time, 1.0)
        gnominal_all[p_i].SetPointError(0, 0, nominal_flux2mev_dict['all'][1]/nominal_flux2mev_dict['all'][0])
        gnominal_all[p_i].SetPointError(1, 0, nominal_flux2mev_dict['all'][1]/nominal_flux2mev_dict['all'][0])
        gnominal_all[p_i].SetFillColorAlpha(ROOT.kGray, ALPHA)
        gnominal_all[p_i].SetFillStyle(3001)
        gnominal_all[p_i].SetLineColor(ROOT.kGray+2)
    ghk.SetMarkerColor(2)
    ghk.SetLineColor(2)
        # Plot just the first run period which is shared between results
    c2 = ROOT.TCanvas("HyoungKu Compare 2021")
    mg = ROOT.TMultiGraph()
    mg.Add(ghk)
    mg.Add(gnominal_indiv[0], "lineE3")
    mg.Add(gmine_indiv[0])
    mg.Add(gtrend_indiv[0], "lineE3")
    mg.Draw("APE*")
    mg.GetXaxis().SetTimeDisplay(1)
    mg.GetXaxis().SetTimeFormat("%Y/%m/%d")
    mg.GetYaxis().SetTitle("Relative Flux per MeV")
    leg2 = ROOT.TLegend(0.1, 0.66, 0.45, 0.9)
    leg2.AddEntry(gmine_indiv[0], "New measurement, %s" % (versionID), "pe")
    leg2.AddEntry(ghk, "Hyoungku's measurement", "pe")
    leg2.AddEntry(gnominal_indiv[0], "Flux-to-MeV (2021) = %0.3f #pm %0.3f" % (nominal_flux2mev_dict[0][0], nominal_flux2mev_dict[0][1]),  "f")
    leg2.Draw('same')
    c2.Update()
    c2.Write()
        # Plot all the run periods
    c2_1 = ROOT.TCanvas("HyoungKu Compare All Runs (Nominal Flux-to-MeV Per Run)")
    mg_1 = ROOT.TMultiGraph()
    for p_i in vals_per_period:
        mg_1.Add(gmine_indiv[p_i])
        mg_1.Add(gtrend_indiv[p_i], "lineE3")
        mg_1.Add(gnominal_indiv[p_i], "lineE3")
    mg_1.Draw("APE*")
    mg_1.GetXaxis().SetTimeDisplay(1)
    mg_1.GetXaxis().SetTimeFormat("%Y/%m/%d")
    mg_1.GetYaxis().SetTitle("Relative Flux per MeV")
    leg2_1 = ROOT.TLegend(0.1, 0.66, 0.45, 0.9)
    leg2_1.AddEntry(gmine_indiv[0], "New measurement, %s" % (versionID), "pe")
    leg2_1.AddEntry(ghk, "Hyoungku's measurement", "pe")
    for p_i in vals_per_period:
        leg2_1.AddEntry(gnominal_indiv[p_i], "Flux-to-MeV (Period %s) = %0.3f #pm %0.3f" % (p_i, nominal_flux2mev_dict[p_i][0], nominal_flux2mev_dict[p_i][1]),  "f")
    leg2_1.Draw('same')
    c2_1.Update()
    c2_1.Write()

    c2_2 = ROOT.TCanvas("HyoungKu Compare All Runs (Nominal Flux-to-MeV For All Runs)")
    mg_2 = ROOT.TMultiGraph()
    for p_i in vals_per_period:
        mg_2.Add(gmine_all[p_i])
        mg_2.Add(gtrend_all[p_i], "lineE3")
        mg_2.Add(gnominal_all[p_i], "lineE3")
    mg_2.Draw("APE*")
    mg_2.GetXaxis().SetTimeDisplay(1)
    mg_2.GetXaxis().SetTimeFormat("%Y/%m/%d")
    mg_2.GetYaxis().SetTitle("Relative Flux per MeV")
    leg2_2 = ROOT.TLegend(0.1, 0.66, 0.45, 0.9)
    leg2_2.AddEntry(gmine_all[0], "New measurement, %s" % (versionID), "pe")
    leg2_2.AddEntry(ghk, "Hyoungku's measurement", "pe")
    leg2_2.AddEntry(gnominal_all[0], "Flux-to-MeV (All runs) = %0.3f #pm %0.3f" % (nominal_flux2mev_dict['all'][0], nominal_flux2mev_dict['all'][1]),  "f")
    leg2_2.Draw('same')
    c2_2.Update()
    c2_2.Write()

    # Plot the normalized conversion with Dodo's results
    c3 = ROOT.TCanvas("Dodo Compare 2021")
    gdodo.SetMarkerColor(ROOT.kGreen)
    gdodo.SetLineColor(ROOT.kGreen)
        # Plot just the first run period which is shared between results
    mg3 = ROOT.TMultiGraph()
    mg3.Add(gdodo, "P*")
    mg3.Add(gmine_indiv[0].Clone(), "P*")
    #mg3.Add(gtrend_indiv[0].Clone(), "lineE3")
    mg3.Draw('APE*')
    mg3.GetXaxis().SetTimeDisplay(1)
    mg3.GetXaxis().SetTimeFormat("%Y/%m/%d")
    mg3.GetYaxis().SetTitle("Relative Flux per MeV")
    leg3 = ROOT.TLegend(0.1, 0.66, 0.45, 0.9)
    leg3.AddEntry(gdodo, "Dodo's measurement (normalized by first run)", "pe")
    leg3.AddEntry(gmine_indiv[0], "New measurement, %s" % (versionID), "pe")
    leg3.Draw('same')
    c3.Update()
    c3.Write()
        # Plot all the run periods

    dodo_vals = np.zeros((len(dodo_correction_forPlots), 3))
    dodo_vals[:, :2] = np.array(sorted(dodo_correction_forPlots.items(), key=lambda x: x[0]))
    dodo_vals[:, 2] = 0.001
    dodo_trend = get_trend_line_simple(dodo_vals)
    g_dodo_trend = ROOT.TGraphErrors()
    for i, (x, y, err), in enumerate(dodo_trend):
        g_dodo_trend.SetPoint(i,x,y)
        g_dodo_trend.SetPointError(i,0,err)
    g_dodo_trend.SetLineColor(ROOT.kGreen)
    g_dodo_trend.SetFillColorAlpha(ROOT.kGreen, ALPHA)
    g_dodo_trend.SetFillStyle(3001)
    
    c4 = ROOT.TCanvas("Dodo Compare All Runs (Nominal Flux-to-MeV Per Run)")
    mg4 = ROOT.TMultiGraph()
    for p_i in vals_per_period:
        mg4.Add(gmine_indiv[p_i].Clone(), "P*")
        mg4.Add(gtrend_indiv[p_i].Clone(), "lineE3")
        mg4.Add(gnominal_indiv[p_i].Clone(), "lineE3")
    mg4.Add(g_dodo_trend, "lineE3")
    mg4.Add(gdodo.Clone(), "P*")

    mg4.Draw("A")
    mg4.GetYaxis().SetTitle("Relative Flux per MeV")
    mg4.GetXaxis().SetTimeDisplay(1);
    mg4.GetXaxis().SetTimeFormat("%Y/%m/%d")
    leg4 = ROOT.TLegend(0.1, 0.66, 0.45, 0.9)
    leg4.AddEntry(gdodo, "Dodo's measurement (normalized by first run)", "pe")
    leg4.AddEntry(gmine_indiv[0], "New measurement, %s" % (versionID), "pe")
    for p_i in vals_per_period:
        leg4.AddEntry(gnominal_indiv[p_i], "Flux-to-MeV (Period %s) = %0.3f #pm %0.3f" % (p_i, nominal_flux2mev_dict[p_i][0], nominal_flux2mev_dict[p_i][1]),  "f")
    leg4.Draw('same')
    c4.Update()
    c4.Write()

#    c5 = ROOT.TCanvas()
#    temp_vals = vals[:, 1:4]
#    temp_vals[:, 1:3] /= NOMINAL_FLUX2MEV
#    combined_dset = np.concatenate((dodo_vals, temp_vals))
#    combined_dset = combined_dset[np.argsort(combined_dset[:, 0])]
#    #combined_trend = get_trend_line_simple(combined_dset)
#    g_combined_trend = ROOT.TGraphErrors()
#    for i, (x, y, err), in enumerate(combined_trend):
#        g_combined_trend.SetPoint(i,x,y)
#        g_combined_trend.SetPointError(i,0,err)
#    g_combined_trend.SetLineColor(ROOT.kRed)
#    g_combined_trend.SetFillColorAlpha(ROOT.kRed, ALPHA)
#    g_combined_trend.SetFillStyle(3001)
#    mg5 = ROOT.TMultiGraph()
#    mg5.Add(g_dodo_trend.Clone(), "lineE3")
#    mg5.Add(gtrend_indiv[0].Clone(), "lineE3")
#    mg5.Add(g_combined_trend, "lineE3")
#    mg5.Draw("A")
#    mg5.GetYaxis().SetTitle("Relative Flux per MeV")
#    mg5.GetXaxis().SetTimeDisplay(1);
#    mg5.GetXaxis().SetTimeFormat("%m/%d")
#    leg5 = ROOT.TLegend(0.14, 0.66, 0.44, 0.88)
#    leg5.AddEntry(g_dodo_trend, "Dodo's nGd measurement (normalized by first run)", "LF")
#    leg5.AddEntry(gtrend_indiv[0], "New Michel measurement", "LF")
#    leg5.AddEntry(g_combined_trend, "Combined measurement", "LF")
#    leg5.Draw()
#    c5.Update()
#    c5.Write()
#
#    c6 = ROOT.TCanvas()
#    mg6 = ROOT.TMultiGraph()
#    gt1 = ROOT.TGraphErrors()
#    gt2 = ROOT.TGraphErrors()
#    gt3 = ROOT.TGraphErrors()
#    gp_dodo = ROOT.TGraphErrors()
#    gp_me = ROOT.TGraphErrors()
#    for i, (d, v) in enumerate(sorted(dodo_correction_forPlots.items(), key=lambda x: x[0])):
#        try:
#            corr = combined_trend[combined_trend[:, 0] > d][0, 1]
#        except IndexError:
#            corr = combined_trend[-1, 1]
#        gp_dodo.SetPoint(i, d, v-corr)
#        gp_dodo.SetPointError(i, 0, 0.001)
#    count = 0
#    for i, (d, v, err) in enumerate(vals[:, 1:4]):
#        if (err > 0.001):
#            continue
#        try:
#            corr = combined_trend[combined_trend[:, 0] > d][0, 1]
#        except IndexError:
#            corr = combined_trend[-1, 1]
#        gp_me.SetPoint(count, d, v-corr)
#        gp_me.SetPointError(count, 0, err)
#        count +=1
#    for i, ((x,y, err), (x0, y0, _)) in enumerate(zip(combined_trend, combined_trend)):
#        gt1.SetPoint(i, x, y-y0)
#        gt1.SetPointError(i, 0, err)
#    for i, ((x,y, err), (x0, y0, _)) in enumerate(zip(dodo_trend, combined_trend)):
#        gt2.SetPoint(i, x, y-y0)
#        #gt2.SetPointError(i, 0, err)
#    for i, ((x, y, err), (x0, y0, _)) in enumerate(zip(trend, combined_trend)):
#        gt3.SetPoint(i, x, y/NOMINAL_FLUX2MEV-y0)
#        #gt3.SetPointError(i, 0, err/NOMINAL_FLUX2MEV)
#    for g, color in [(gt1, ROOT.kRed), (gt2, ROOT.kGreen), (gt3, ROOT.kBlue)]:
#        g.SetLineColor(color)
#        g.SetLineWidth(2)
#        g.SetFillColorAlpha(color, 0.2)
#        g.SetFillStyle(3001)
#    gp_dodo.SetLineColor(ROOT.kGreen)
#    gp_me.SetLineColor(ROOT.kBlue)
#    gp_dodo.SetMarkerColor(ROOT.kGreen)
#    gp_me.SetMarkerColor(ROOT.kBlue)
#    mg6.Add(gt1, "line E3")
#    mg6.Add(gt2, "line E3")
#    mg6.Add(gt3, "line E3")
#    mg6.Add(gp_dodo, "*E")
#    mg6.Add(gp_me, "*E")
#    mg6.Draw("A")
#    mg6.GetYaxis().SetTitle("Relative Flux per Mev Difference")
#    mg6.GetYaxis().SetTitleOffset(1.25)
#    mg6.GetXaxis().SetTimeDisplay(1);
#    mg6.GetXaxis().SetTimeFormat("%Y/%m/%d")
#    c6.Update()
#    c6.Write()

    # Compare to Eric's data
        # Plot just the first run period which is shared between results
    eric_trend = get_trend_line_simple(eric_data[:, 1:4])
    geric_trend = ROOT.TGraphErrors()
    for i, (x, y, err), in enumerate(eric_trend):
        geric_trend.SetPoint(i,x,y/eric_AVG_FLUX2MEV)
        geric_trend.SetPointError(i,0,err/eric_AVG_FLUX2MEV)
    geric_trend.SetLineColor(ROOT.kGreen)
    geric_trend.SetFillColorAlpha(ROOT.kGreen, ALPHA)
    geric_trend.SetFillStyle(3001)

    c7 = ROOT.TCanvas("Eric Compare 2021")
    geric.SetMarkerColor(ROOT.kGreen)
    geric.SetLineColor(ROOT.kGreen)
    mg7 = ROOT.TMultiGraph()
    mg7.Add(geric)
    mg7.Add(geric_trend, "lineE3")
    #mg7.Add(gnominal_indiv[0].Clone(), "lineE3")
    mg7.Add(gmine_indiv[0].Clone())
    mg7.Add(gtrend_indiv[0].Clone(), "lineE3")
    mg7.Draw("APE*")
    mg7.GetXaxis().SetTimeDisplay(1)
    mg7.GetXaxis().SetTimeFormat("%Y/%m/%d")
    mg7.GetYaxis().SetTitle("Relative Flux per MeV")
    leg7 = ROOT.TLegend(0.1, 0.66, 0.45, 0.9)
    leg7.AddEntry(gmine_indiv[0], "New measurement, %s" % (versionID), "pe")
    leg7.AddEntry(geric, "Eric's measurement", "pe")
    #leg7.AddEntry(geric_trend, "Eric's Flux-to-MeV")
    #leg7.AddEntry(gnominal_indiv[0], "Flux-to-MeV (2021) = %0.3f #pm %0.3f" % (nominal_flux2mev_dict[0][0], nominal_flux2mev_dict[0][1]),  "f")
    leg7.Draw('same')
    c7.Update()
    c7.Write()
        # Plot all the run periods
    c7_1 = ROOT.TCanvas("Eric Compare All Runs (Nominal Flux-to-MeV Per Run)")
    mg7_1 = ROOT.TMultiGraph()
    for p_i in vals_per_period:
        mg7_1.Add(gmine_indiv[p_i].Clone())
        mg7_1.Add(gtrend_indiv[p_i].Clone(), "lineE3")
        mg7_1.Add(gnominal_indiv[p_i].Clone(), "lineE3")
    mg7_1.Add(geric_trend.Clone(), "lineE3")
    mg7_1.Add(geric.Clone())
    mg7_1.Draw("APE*")
    mg7_1.GetXaxis().SetTimeDisplay(1)
    mg7_1.GetXaxis().SetTimeFormat("%Y/%m/%d")
    mg7_1.GetYaxis().SetTitle("Relative Flux per MeV")
    leg7_1 = ROOT.TLegend(0.1, 0.66, 0.45, 0.9)
    leg7_1.AddEntry(gmine_indiv[0], "New measurement, %s" % (versionID), "pe")
    leg7_1.AddEntry(geric, "Eric's measurement", "pe")
    leg7_1.AddEntry(geric_trend, "Eric's Flux-to-MeV")
    for p_i in vals_per_period:
        leg7_1.AddEntry(gnominal_indiv[p_i], "Flux-to-MeV (Period %s) = %0.3f #pm %0.3f" % (p_i, nominal_flux2mev_dict[p_i][0], nominal_flux2mev_dict[p_i][1]),  "f")
    leg7_1.Draw('same')
    c7_1.Update()
    c7_1.Write()

    # Compare to optional version
    if compareID :
        c8 = ROOT.TCanvas()
        mg8 = ROOT.TMultiGraph()
        for p_i in vals_per_period:
            gcomp_indiv[p_i].SetMarkerColor(ROOT.kRed)
            gcomp_indiv[p_i].SetLineColor(ROOT.kRed)
            mg8.Add(gmine_indiv[p_i].Clone(), "P*")
            mg8.Add(gcomp_indiv[p_i], "P*")
            #mg8.Add(gnominal_indiv[p_i].Clone(), "lineE3")
        mg8.Draw('APE*')
        mg8.GetXaxis().SetTimeDisplay(1)
        mg8.GetXaxis().SetTimeFormat("%m/%d")
        mg8.GetYaxis().SetTitle("Relative Flux per MeV")
        leg8 = ROOT.TLegend(0.14, 0.66, 0.44, 0.88)
        leg8.AddEntry(gcomp_indiv[0], "Version %s" % (compareID), "PE")
        leg8.AddEntry(gmine_indiv[0], "Version %s" % (versionID), "PE")
        #for p_i in vals_per_period:
        #    leg8.AddEntry(gnominal_indiv[p_i], "Flux-to-MeV (Period %s) = %0.3f #pm %0.3f" % (p_i, nominal_flux2mev_dict[p_i][0], nominal_flux2mev_dict[p_i][1]),  "f")
        leg8.Draw()
        c8.Update()
        c8.Write()

        gcomp_trend = []
        for p_i, p_v in comp_vals.items():
            gcomp_trend.append(ROOT.TGraphErrors())

            trend = get_trend_line_simple(p_v[:, 1:4])
            for i, (x, y, err) in enumerate(trend):
                gcomp_trend[p_i].SetPoint(i,x,y/comp_nominal_flux2mev_dict[p_i][0])
                gcomp_trend[p_i].SetPointError(i,0,err/comp_nominal_flux2mev_dict[p_i][0])
            gcomp_trend[p_i].SetLineColor(ROOT.kRed)
            gcomp_trend[p_i].SetFillColorAlpha(ROOT.kRed, ALPHA)
            gcomp_trend[p_i].SetFillStyle(3001)


    # Plot the trendlines only for current, Eric's, and compareison version if applicable
    c9 = ROOT.TCanvas("Trends")
    mg9 = ROOT.TMultiGraph()
    for p_i in vals_per_period:
        mg9.Add(gtrend_indiv[p_i].Clone(), "lineE3")
        if compareID:
            mg9.Add(gcomp_trend[p_i].Clone(), "lineE3")
    mg9.Add(geric_trend.Clone(), "lineE3")
    mg9.Draw('A')
    mg9.GetXaxis().SetTimeDisplay(1)
    mg9.GetXaxis().SetTimeFormat("%m/%d")
    mg9.GetYaxis().SetTitle("Relative Flux per MeV")
    leg9 = ROOT.TLegend(0.1, 0.66, 0.45, 0.9)
    leg9.AddEntry(gtrend_indiv[0], "New measurement trend, %s" % (versionID), "LF")
    leg9.AddEntry(geric_trend, "Eric's measurement trend", "LF")
    if compareID:
        leg9.AddEntry(gcomp_trend[0], "Comparison version trend, %s" % (compareID), "LF")
    leg9.Draw()
    c9.Update()
    c9.Write()

    outFile.Close()


