"""
Taken from McCaffrey:
/home/marzece/KDAR_Analysis/MichelAnalysis/f2mev_correction/plot_time_flux.py

description tbd

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

JSONFILE = "/home/littleca/kdar/correction_values.json"
OUTFILE_PATH = "/home/littleca/kdar/cleanerKDAR/Michel/flux2mev_correction/output_fluxCorr/"
WORK_PATH= "/home/littleca/kdar/cleanerKDAR/Michel/flux2mev_correction/"

FIRST_RUN = 1597
ALPHA = 0.35
TREND_WIDTH = 2.0*24*3600
DOWNTIME_BTWN_PERIODS = 5.26e6    # Seconds

# From https://github.com/JSNS2/AnalysisTools/blob/main/EnergyCorrection/20240312/TimeCorrection_2021.txt
dodo_correction = [[1596,0,1600,2531,1], [1600,2532,1609,2817,1.00182], [1609,2818,1625,1791,1.0058], [1625,1792,1625,6791,1.00486], [1625,6792,1626,915,1.00604], [1626,916,1634,810,1.00505], [1634,811,1672,637,1.00625], [1672,638,1678,544,1.00466], [1678,545,1678,5544,1.00406], [1678,5545,1684,1607,1.00375], [1684,1608,1684,6607,1.00273], [1684,6608,1689,9,1.00395], [1689,10,1689,5009,1.00644], [1689,5010,1689,10009,1.00338], [1689,10010,1692,1719,1.00691], [1692,1720,1717,1996,1.00915], [1717,1997,1718,1104,1.01197], [1718,1105,1718,6104,1.01135], [1718,6105,1722,1281,1.01209], [1722,1282,1725,330,1.01121], [1725,331,1726,1637,1.01254], [1726,1638,1729,3434,1.01154], [1729,3435,1749,3244,1.0146], [1749,3245,1769,681,1.0108], [1769,682,1785,1724,1.00911], [1785,1725,1787,4457,1.01003], [1787,4458,1788,4147,1.01104], [1788,4148,1790,3500,1.01153], [1790,3501,1814,825,1.01248], [1814,826,1815,4482,1.00921], [1815,4483,1818,1196,1.00951], [1818,1197,1826,1543,1.01118], [1826,1544,1828,2560,1.01114], [1828,2561,1830,2419,1.0121], [1830,2420,1831,1358,1.01458], [1831,1359,1837,993,1.01308], [1837,994,1839,2722,1.01577], [1839,2723,1839,7722,1.01571], [1839,7723,1843,4380,1.02102], [1843,4381,1850,3640,1.02025], [1850,3641,1850,8640,1.02026], [1850,8641,1850,13640,1.02502], [1850,13641,1864,357,1.02644], [1864,358,1867,239,1.02869], [1867,240,1874,302,1.02494], [1874,303,1878,2860,1.02688], [1878,2861,1878,7860,1.0261], [1878,7861,1879,1355,1.02415], [1879,1356,1909,2369,1.02553], [1909,2370,1909,7369,1.02369], [1909,7370,1910,4198,1.02243], [1910,4199,1914,1908,1.02243], [1914,1909,1916,3446,1.02531], [1916,3447,1916,8446,1.02344], [1916,8447,1917,2033,1.02223], [1917,2034,1922,151,1.02324], [1922,152,1923,3166,1.02255], [1923,3167,1930,54,1.02148], [1930,55,1936,1808,1.02211], [1936,1809,1940,3594,1.01947], [1940,3595,1944,1517,1.01794], [1944,1518,1944,6517,1.01917], [1944,6518,1950,193,1.01657], [1950,194,1952,1348,1.0176], [1952,1349,1953,220,1.01794]]


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

    # debug
    for x0 in xax:
        w = gaussian(points[:,0], mu=x0, sigma=sigma)
    #    if sum(w) < 1e-50:
    #        print(sum(w), x0)

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
    path = "/home/marzece/KDAR_Analysis/MichelAnalysis/f2mev_correction/FluxMeVCorr_time_blessed"
    hk_vals = []
    for filename in os.listdir(path):
        if re.match("Run_\d+_timeCorr.dat", filename):
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

    # Find the average value for flux to MeV and it's error for all run periods and for each individual period
    # We will normalize our results later with respect to this average
    nominal_f2mev_dict = defaultdict(lambda: [])

    # TODO remove these after we see that  using vals_per_period works
    #flux2mev_vals = [value_dict["Flux2MeV"] for entry_type, value_dict in data_vals[versionID].items() if isinstance(value_dict["Run"], int)]
    #flux2mev_err_vals = [value_dict["Flux2MeVErr"] for entry_type, value_dict in data_vals[versionID].items() if isinstance(value_dict["Run"],int)]
    # debug
    #print(vals_per_period.values())
    all_vals = vals_per_period.values()[0]

    all_f2mev = np.average(all_vals[:,2], weights=list(map(lambda x: 1.0/x, all_vals[:,2])))
    all_f2mev_err = np.average(all_vals[:,3], weights=list(map(lambda x: 1.0/x, all_vals[:,3])))

    nominal_f2mev_dict['all'] = np.array([all_f2mev, all_f2mev_err])

    for p_i in vals_per_period:
        n_f2mev = np.average(vals_per_period[p_i][:,2], weights=list(map(lambda x: 1.0/x, vals_per_period[p_i][:,2])))
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

    # Histograms for comparing the current version's flux 2 MeV with our other data sets
    eric_comp_h = ROOT.TH1D("Fraction of Eric's Data","Fraction of Eric's Data",500,-1,2)
    dodo_comp_h = ROOT.TH1D("Fraction of Dodo's Data","Fractio of Dodo's Data",500,-1,2)
    hk_comp_h = ROOT.TH1D("Fraction of HK's Data","Fraction of HK'd Data",500,-1,2)
    comp_h = ROOT.TH1D("Fraction of Version %s's Data" % (compareID),"Fraction of Version %s's Data" % (compareID),500,-1,2)

    # Get numerical comparison between the different corrections
    # HyoungKu, Dodo, and Eric's data uses a michel endpoint = 52.8 MeV
    endpoint_scale = 52.8/MICHEL_ENDPOINT
    if args.compare:
        compare_endpoint_scale = compare_MICHEL_ENDPOINT/MICHEL_ENDPOINT

    # Fill TGraph data for various data sets
    gabs = []						# Array of TGraphErrors of current version data
    gmine_indiv = []				# Array of TGraphErros of current version data (normalized wrt average flux to MeV per period)
    gmine_all = []					# Array of TGraphErros of current version data (normalized wrt average flux to MeV for all periods)
    ghk = ROOT.TGraphErrors()		# HyoungKu's data (normalized wrt average flux to MeV)
    gdodo = ROOT.TGraph()			# Dodo's data (normalized wrt first run's flux to MeV)
    geric = ROOT.TGraphErrors()		# Eric's data (normalized wrt average flux to MeV)

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

    # Let's get the POT for each period we're analyzing

    # ===== HyoungKu's Data =====
    hk_vals = get_hyoungkus_correction(run_times)

    # To normalize HK's data with respect to the average flux to mev value
    # Change hk_AV_FLUX2MEV = 1 if we want to normalize with respect to HKs first run
    hk_AVG_FLUX2MEV = np.average(hk_vals[:,2], weights=list(map(lambda x: 1.0/x, hk_vals[:,3])))

    for i, (r, start, v, err) in enumerate(hk_vals):
        ghk.SetPoint(i, start, v/hk_AVG_FLUX2MEV)
        ghk.SetPointError(i, 0, err/hk_AVG_FLUX2MEV)

        # Make sure we're getting the correct run data from vals_per_period
        run_idx = np.where(vals_per_period[0][:,0] == r)[0]
        if len(run_idx) == 1 :
            hk_comp_h.Fill((vals_per_period[0][run_idx[0], 2]/nominal_flux2mev_dict[0][0])/(v/hk_AVG_FLUX2MEV)*endpoint_scale)

    # ===== Dodo's Data =====
    # TODO Double check I can't normalize wrt average like the others. I know it's already normalized but can I, like, do it again?
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

        #vals_v_l = [val_list[2] for val_list in vals if val_list[0]==vr_tup[1]]
        #if len(vals_v_l) == 1:
        #    vals_v = vals_v_l[0]
        #    dodo_comp_h.Fill((vals_v/NOMINAL_FLUX2MEV)/vr_tup[0]*endpoint_scale)
        # TODO Were these error catching things necessary?
        #if len(current_vals_v) > 1 : 
        #    print("Error: multiple runs with the same number found wile filling Eric's comparison histogram")
        #elif len(current_vals_v) == 1 :
        #    dodo_comp_h.Fill((current_vals_v[0]/NOMINAL_FLUX2MEV)/v*endpoint_scale)
        #else :
        #    print("Nothing to compare Dodo's run: %i to. Current flux2mev was found to be " % (vr_tup[1]), current_vals_v) 
 
    # ===== Eric's Data =====
    # Read in Eric's correction values to compare against
    eric_data = []
    with open("/home/marzece/KDAR_Analysis/MichelAnalysis/f2mev_correction/temp.txt", 'r') as efile:
        csvFile = csv.reader(efile)
        for line in csvFile:
            eric_data.append([float(i) for i in line])

    for i, run_data in enumerate(eric_data):
        eric_timestamp = (time.mktime(run_times[run_data[0]][0].timetuple()) - ROOT.gStyle.GetTimeOffset())
        eric_data[i].insert(1, eric_timestamp)
        # Run	Time	Flux2MeV	Err

    # Get the average flux to MeV value 
    eric_data_list = zip(*eric_data)
    eric_AVG_FLUX2MEV = np.average(list(eric_data_list[2]), weights=list(map(lambda x: 1.0/x, list(eric_data_list[3]))))

    # debug
    eric_frac_dict = defaultdict(lambda : [])

    for i, (r, start, v, err) in enumerate(eric_data):
        geric.SetPoint(i, start, v/eric_AVG_FLUX2MEV)
        geric.SetPointError(i,0, err/eric_AVG_FLUX2MEV)

		# Find corresponding runs flux to MeV value for current version
        # Make sure we're getting the correct run data from vals_per_period
        run_idx = np.where(vals_per_period[0][:,0] == r)[0]
        if len(run_idx) == 1 :
            eric_frac = vals_per_period[0][run_idx[0], 2]/v*endpoint_scale
            eric_comp_h.Fill(eric_frac)
            # debug
            eric_frac_dict[r].append([eric_frac, v-vals_per_period[0][run_idx[0], 2]])

    # debug
    # Write out differences between each runs flux to mev for current and Erics version
    with open(OUTFILE_PATH+"/"+str(versionID)+"_EricCompareRuns.txt", "a") as fe:
        for r, v in eric_frac_dict.items() :
            fe.write(r+"\t"+v[0]+"\t"+v[1])   

    # ===== Comparison Data =====
    # Read in corrections for optional JSON version to compare against
    if args.compare :
        print("Comparing to : " + compareID)
        gcomp = ROOT.TGraphErrors()
        comp_vals, comp_NOMINAL_FLUX2MEV, comp_NOMINAL_FLUX2MEV_ERR = get_json_correction(compareID, run_times)

        for i, (r, start, v, err) in enumerate(comp_vals):
            gcomp.SetPoint(count, start, v/comp_NOMINAL_FLUX2MEV)
            gcomp.SetPointError(count,0, err/comp_NOMINAL_FLUX2MEV)

		    # Find corresponding runs flux to MeV value for current version
            # Make sure we're getting the correct run data from vals_per_period
            run_idx = np.where(vals_per_period[0][:,0] == r)[0]
            if len(run_idx) == 1 :
                comp_h.Fill(vals_per_period[0][run_idx[0], 2]/v*compare_endpoint_scale)
            #vals_v_l = [val_list[2] for val_list in vals if val_list[0]==r]
            #if len(vals_v_l) == 1:
            #    vals_v = vals_v_l[0]
            #    comp_h.Fill(vals_v/v*compare_endpoint_scale)
            # TODO Were these error catching things necessary?
            #if len(current_vals_v) > 1 : 
            #    print("Error: multiple runs with the same number found wile filling Eric's comparison histogram")
            #elif len(current_vals_v) == 1 :
            #    comp_h.Fill(current_vals_v[0]/v*compare_endpoint_scale)
            #else :
            #    print("Nothing to compare %s's run: %i to. Current flux2mev was found to be " % (compareID, r), current_vals_v) 


    # debug? Checking the difference between points
    # Specifically between the ericBinDiv02 and Eric's data
#    gComp = ROOT.TGraph()
#    hComp = ROOT.TH1D("ericBinDiv02 - Eric's data","ericBinDiv02 - Eric's data",1000,-2,2)
#    for i, (r, start, v, err) in enumerate(vals):
#        formatDate = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(start))
#        #print(r," Start: ",start," ",formatDate," Flux2MevScale: ",v," Error: ", err)
#        gComp.SetPoint(i, start, (v/NOMINAL_FLUX2MEV)-(v_eric_data[i]/eric_AVG_FLUX2MEV))
#        hComp.Fill((v/NOMINAL_FLUX2MEV)-(v_eric_data[i]/eric_AVG_FLUX2MEV))

    OUTFILE = OUTFILE_PATH+"/"+str(versionID)+"/time_flux_plots_"+str(versionID)+".root" 
    outFile = ROOT.TFile(OUTFILE, "RECREATE")
 
    eric_comp_h.Write()
    dodo_comp_h.Write()
    hk_comp_h.Write()
    comp_h.Write()
    #debug?
 #   hComp.Write()

    # debug? Plot the difference between points
 #   c10 = ROOT.TCanvas()
 #   gComp.GetXaxis().SetTimeDisplay(1)
 #   gComp.GetXaxis().SetTimeFormat("%m/%d")
 #   gComp.GetYaxis().SetTitle("New ericBinDiv02 - Eric's data")
 #   gComp.Draw("AP*")
  #  c10.Update()
  #  c10.Write()

    # TODO Add the distinction between normalizing wrt all period f2mev vs individual periods
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
    # TODO plot each run individually, too
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
    mg3.Add(gmine[0].Clone(), "P*")
    #mg3.Add(gtrend.Clone(), "lineE3")
    mg3.Draw('APE*')
    mg3.GetXaxis().SetTimeDisplay(1)
    mg3.GetXaxis().SetTimeFormat("%Y/%m/%d")
    mg3.GetYaxis().SetTitle("Relative Flux per MeV")
    leg3 = ROOT.TLegend(0.1, 0.66, 0.45, 0.9)
    leg3.AddEntry(gdodo, "Dodo's measurement (normalized by first run)", "pe")
    leg3.AddEntry(gmine, "New measurement, %s" % (versionID), "pe")
    leg3.Draw('same')
    c3.Update()
    c3.Write()
        # Plot all the run periods

    dodo_vals = np.zeros((len(dodo_correction_forPlots), 3))
    dodo_vals[:, :2] = np.array(sorted(dodo_correction_forPlots.items(), key=lambda x: x[0]))
    dodo_vals[:, 2] = 0.001
    #dodo_trend = get_trend_line_simple(dodo_vals)
    g_dodo_trend = ROOT.TGraphErrors()
    for i, (x, y, err), in enumerate(dodo_trend):
        g_dodo_trend.SetPoint(i,x,y)
        g_dodo_trend.SetPointError(i,0,err)
    g_dodo_trend.SetLineColor(ROOT.kGreen)
    g_dodo_trend.SetFillColorAlpha(ROOT.kGreen, ALPHA)
    g_dodo_trend.SetFillStyle(3001)
    c4 = ROOT.TCanvas()
    mg4 = ROOT.TMultiGraph()
    mg4.Add(g_dodo_trend, "lineE3")
    mg4.Add(gtrend.Clone(), "lineE3")

    mg4.Add(gdodo.Clone(), "P*")
    mg4.Add(gmine.Clone(), "P*")

    mg4.Draw("A")
    mg4.GetYaxis().SetTitle("Relative Flux per MeV")
    mg4.GetXaxis().SetTimeDisplay(1);
    mg4.GetXaxis().SetTimeFormat("%Y/%m/%d")
    leg4 = ROOT.TLegend(0.1, 0.66, 0.45, 0.9)
    leg4.AddEntry(gdodo, "Dodo's measurement (normalized by first run)", "pe")
    leg4.AddEntry(gmine, "New measurement, %s" % (versionID), "pe")
    leg2.AddEntry(gtrend, "Flux-to-MeV = %0.3f #pm %0.3f" % (NOMINAL_FLUX2MEV, NOMINAL_FLUX2MEV_ERR),  "f")
    leg4.Draw('same')
    c4.Update()
    c4.Write()

    c5 = ROOT.TCanvas()
    temp_vals = vals[:, 1:4]
    temp_vals[:, 1:3] /= NOMINAL_FLUX2MEV
    combined_dset = np.concatenate((dodo_vals, temp_vals))
    combined_dset = combined_dset[np.argsort(combined_dset[:, 0])]
    #combined_trend = get_trend_line_simple(combined_dset)
    g_combined_trend = ROOT.TGraphErrors()
    for i, (x, y, err), in enumerate(combined_trend):
        g_combined_trend.SetPoint(i,x,y)
        g_combined_trend.SetPointError(i,0,err)
    g_combined_trend.SetLineColor(ROOT.kRed)
    g_combined_trend.SetFillColorAlpha(ROOT.kRed, ALPHA)
    g_combined_trend.SetFillStyle(3001)
    mg5 = ROOT.TMultiGraph()
    mg5.Add(g_dodo_trend.Clone(), "lineE3")
    mg5.Add(gtrend.Clone(), "lineE3")
    mg5.Add(g_combined_trend, "lineE3")
    mg5.Draw("A")
    mg5.GetYaxis().SetTitle("Relative Flux per MeV")
    mg5.GetXaxis().SetTimeDisplay(1);
    mg5.GetXaxis().SetTimeFormat("%m/%d")
    leg5 = ROOT.TLegend(0.14, 0.66, 0.44, 0.88)
    leg5.AddEntry(g_dodo_trend, "Dodo's nGd measurement (normalized by first run)", "LF")
    leg5.AddEntry(gtrend, "New Michel measurement", "LF")
    leg5.AddEntry(g_combined_trend, "Combined measurement", "LF")
    leg5.Draw()
    c5.Update()
    c5.Write()

    c6 = ROOT.TCanvas()
    mg6 = ROOT.TMultiGraph()
    gt1 = ROOT.TGraphErrors()
    gt2 = ROOT.TGraphErrors()
    gt3 = ROOT.TGraphErrors()
    gp_dodo = ROOT.TGraphErrors()
    gp_me = ROOT.TGraphErrors()
    for i, (d, v) in enumerate(sorted(dodo_correction_forPlots.items(), key=lambda x: x[0])):
        try:
            corr = combined_trend[combined_trend[:, 0] > d][0, 1]
        except IndexError:
            corr = combined_trend[-1, 1]
        gp_dodo.SetPoint(i, d, v-corr)
        gp_dodo.SetPointError(i, 0, 0.001)
    count = 0
    for i, (d, v, err) in enumerate(vals[:, 1:4]):
        if (err > 0.001):
            continue
        try:
            corr = combined_trend[combined_trend[:, 0] > d][0, 1]
        except IndexError:
            corr = combined_trend[-1, 1]
        gp_me.SetPoint(count, d, v-corr)
        gp_me.SetPointError(count, 0, err)
        count +=1
    for i, ((x,y, err), (x0, y0, _)) in enumerate(zip(combined_trend, combined_trend)):
        gt1.SetPoint(i, x, y-y0)
        gt1.SetPointError(i, 0, err)
    for i, ((x,y, err), (x0, y0, _)) in enumerate(zip(dodo_trend, combined_trend)):
        gt2.SetPoint(i, x, y-y0)
        #gt2.SetPointError(i, 0, err)
    for i, ((x, y, err), (x0, y0, _)) in enumerate(zip(trend, combined_trend)):
        gt3.SetPoint(i, x, y/NOMINAL_FLUX2MEV-y0)
        #gt3.SetPointError(i, 0, err/NOMINAL_FLUX2MEV)
    for g, color in [(gt1, ROOT.kRed), (gt2, ROOT.kGreen), (gt3, ROOT.kBlue)]:
        g.SetLineColor(color)
        g.SetLineWidth(2)
        g.SetFillColorAlpha(color, 0.2)
        g.SetFillStyle(3001)
    gp_dodo.SetLineColor(ROOT.kGreen)
    gp_me.SetLineColor(ROOT.kBlue)
    gp_dodo.SetMarkerColor(ROOT.kGreen)
    gp_me.SetMarkerColor(ROOT.kBlue)
    mg6.Add(gt1, "line E3")
    mg6.Add(gt2, "line E3")
    mg6.Add(gt3, "line E3")
    mg6.Add(gp_dodo, "*E")
    mg6.Add(gp_me, "*E")
    mg6.Draw("A")
    mg6.GetYaxis().SetTitle("Relative Flux per Mev Difference")
    mg6.GetYaxis().SetTitleOffset(1.25)
    mg6.GetXaxis().SetTimeDisplay(1);
    mg6.GetXaxis().SetTimeFormat("%Y/%m/%d")
    c6.Update()
    c6.Write()

    c7 = ROOT.TCanvas()
    geric.SetMarkerColor(ROOT.kGreen)
    mg3 = ROOT.TMultiGraph()
    mg3.Add(geric, "P*")
    mg3.Add(gmine.Clone(), "P*")
    mg3.Draw('APE*')
    mg3.GetXaxis().SetTimeDisplay(1)
    mg3.GetXaxis().SetTimeFormat("%Y/%m/%d")
    mg3.GetYaxis().SetTitle("Relative Flux per MeV")
    leg7 = ROOT.TLegend(0.14, 0.66, 0.44, 0.88)
    leg7.AddEntry(geric, "Eric's measurement (Endpoint=52.8 MeV, Bin*~0.2)", "PE")
    leg7.AddEntry(gmine, "New Michel measurement, %s" % (versionID), "PE")
    leg7.Draw()
    c7.Update()
    c7.Write()

    c9 = ROOT.TCanvas()
    eric_trend = get_trend_line_simple(eric_data[:, 1:4])
    g_eric_trend = ROOT.TGraphErrors()
    for i, (x, y, err), in enumerate(eric_trend):
        g_eric_trend.SetPoint(i,x,y)
        g_eric_trend.SetPointError(i,0,err)
    g_eric_trend.SetLineColor(ROOT.kGreen)
    g_eric_trend.SetFillColorAlpha(ROOT.kGreen, ALPHA)
    g_eric_trend.SetFillStyle(3001)
    mg9 = ROOT.TMultiGraph()
    mg9.Add(g_dodo_trend.Clone(), "lineE3")
    mg9.Add(gtrend.Clone(), "lineE3")
    mg9.Draw("A")
    mg9.GetYaxis().SetTitle("Relative Flux per MeV")
    mg9.GetXaxis().SetTimeDisplay(1);
    mg9.GetXaxis().SetTimeFormat("%Y/%m/%d")
    leg9 = ROOT.TLegend(0.14, 0.66, 0.44, 0.88)
    leg9.AddEntry(g_eric_trend, "Eric's measurement", "LF")
    leg9.AddEntry(gtrend, "New Michel measurement", "LF")
    leg9.Draw()
    c5.Update()
    c5.Write()

    if compareID :
        
        c8 = ROOT.TCanvas()
        gcomp.SetMarkerColor(ROOT.kGreen)
        mg4 = ROOT.TMultiGraph()
        mg4.Add(gcomp, "P*")
        mg4.Add(gmine.Clone(), "P*")
        mg4.Draw('APE*')
        mg4.GetXaxis().SetTimeDisplay(1)
        mg4.GetXaxis().SetTimeFormat("%m/%d")
        mg4.GetYaxis().SetTitle("Relative Flux per MeV")
        leg8 = ROOT.TLegend(0.14, 0.66, 0.44, 0.88)
        leg8.AddEntry(gcomp, "Version %s" % (compareID), "PE")
        leg8.AddEntry(gmine, "Version %s" % (versionID), "PE")
        leg8.Draw()
        c8.Update()
        c8.Write()



    # Prepare for writing values to output file
    def time_to_trend_val(t):
        if(t >= combined_trend[-1, 0]):
            idx = -1
        else:
            idx = np.where(combined_trend[:, 0] > t)[0][0]
        return combined_trend[idx,1], combined_trend[idx, 2]
    zz = time_to_rn(None)
    info = [(rn, srn, dt) for rn, srn, max_srn, dt in zz if srn in [0, 5000, 10000]]
    info2 = [(rn, srn, time.mktime(dt.timetuple()) - ROOT.gStyle.GetTimeOffset()) for rn, srn, dt in info]
    info3 = []
#    for rn, srn, tv in info2:
        #print(rn,srn,tv)
#        info3.append((rn,srn, time_to_trend_val(tv)))
#    info4 = [[rn, srn, NOMINAL_FLUX2MEV*v, verr*NOMINAL_FLUX2MEV] for rn, srn, v, verr in info3]

