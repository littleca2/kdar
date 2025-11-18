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

# Endpoint value for data, 0 = 52.8 MeV		1 = 53.3 MeV
endpoint_ver = 1
comp_endpoint_ver = 1

JSONFILE = "/home/littleca/kdar/correction_values.json"
OUTFILE_PATH = "/home/littleca/kdar/cleanerKDAR/Michel/flux2mev_correction/output_fluxCorr/"
WORK_PATH= "/home/littleca/kdar/cleanerKDAR/Michel/flux2mev_correction/"

ALPHA = 0.35
TREND_WIDTH = 2.0*24*3600

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

    # I dont think this is the correct way to fix there being large gaps in dates
    time_weights_list = []
    for x0 in xax:
        weight_list = gaussian(points[:, 0], mu=x0, sigma=sigma)
        if sum(weight_list) < 1e-100:
            #i = 2
            #while sum(weight_list)==0:
            #    weight_list = gaussian(points[:, 0], mu=x0, sigma=sigma*i)
            #    i+=1
            weight_list = time_weights_list[-1]
        time_weights_list.append(weight_list)

    time_weights = np.array(time_weights_list)
    print(time_weights.shape)

    # debug
    for x0 in xax:
        s = gaussian(points[:,0], mu=x0, sigma=sigma)
        if sum(s)==0:
            formatDate = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(start))
            #print("times ", points[:0])
            #print(formatDate)
            #print("mu ", x0)

    mean_vals = np.array([np.average(points[:, 1], weights=w*point_err_weight) for w in time_weights])
    results[:, 1] = mean_vals

    variance = np.array([np.average((points[:,1]-mean)**2,weights=w*point_err_weight)  for w, mean in zip(time_weights,mean_vals)])
    std_dev = np.sqrt(variance)

    weighted_counts = time_weights/peak_val
    counts = np.sum(weighted_counts, axis=1)
    mean_err = std_dev/np.sqrt(counts)
    #debug
    for i, err in enumerate(mean_err):
        if err > 10e100:
            print("Mean Error: ", err)
            print("Time weight: ", time_weights[i])
            print("index: ", i)

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
    #import ipdb;ipdb.set_trace()

    err = np.sqrt(std_dev**2 + 0*error_sum**2)
    result[:, 2] = err
    return result

def get_run_dates():
    # debug
    #pdb.set_trace()

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


def get_hyoungkus_correction(run):
    path = "/home/marzece/KDAR_Analysis/MichelAnalysis/f2mev_correction/FluxMeVCorr_time_blessed"
    fn = "Run_%i_timeCorr.dat" % run
    f_join = os.path.join(path, fn)

    if os.path.isfile(f_join):
        fin = open(os.path.join(path, fn))
        lines = fin.readlines()
        return [float(x) for x in lines[0].strip().split('\t')][2:]
    else :
        print("COULD NOT OPEN ", f_join)
        return [None,None]

def get_json_correction(versionID):
    # Read in data from json file
    with open(JSONFILE, "r") as f:
        data_vals = json.load(f)

    # Get the name of the "sum"key for this version
    # There should be only be one value for sum
#    sumlist = [key for key, value in data_vals[versionID].items() if 'sum' in key]
#    if len(sumlist) != 1 :
#        print("ERROR: Multiple sum entries (", len(sumlist),") present for version ", versionID)
#        exit()

    # Finding the flux2mev of the earliest run to normalize against
    firstRun = 9999
    for runID in list(data_vals[versionID]) :
        run = data_vals[versionID][runID]["Run"]
        if run is None:
            continue
        elif run == 1597:		# run < firstRun
            firstRun = run
            firstRunID = runID
#    NOMINAL_FLUX2MEV = data_vals[versionID][runID]["Flux2MeV"]			# Normalize all flux to MeV vlaues so that the first run is "1", currently hardcoded to bypass point with large error
#    NOMINAL_FLUX2MEV_ERR = data_vals[versionID][runID]["Flux2MeVErr"]

 #   print("First run: ", firstRun)

    # We will normalize the flux to MeV values based on their average value
    flux2mev_vals = [value_dict["Flux2MeV"] for entry_type, value_dict in data_vals[versionID].items() if isinstance(value_dict["Run"], int)]
    flux2mev_err_vals = [value_dict["Flux2MeVErr"] for entry_type, value_dict in data_vals[versionID].items() if isinstance(value_dict["Run"],int)]
#    NOMINAL_FLUX2MEV = sum(flux2mev_vals)/float(len(flux2mev_vals))
#    NOMINAL_FLUX2MEV_ERR = sum(flux2mev_err_vals)/float(len(flux2mev_err_vals))
    NOMINAL_FLUX2MEV = np.average(flux2mev_vals, weights=list(map(lambda x: 1.0/x, flux2mev_err_vals)))
    NOMINAL_FLUX2MEV_ERR = np.average(flux2mev_err_vals, weights=list(map(lambda x: 1.0/x, flux2mev_err_vals)))

    run_times = get_run_dates()

    # Set the timing offset to the first run to normalize the rest of the times
    ROOT.gStyle.SetTimeOffset(time.mktime(run_times[firstRun][0].timetuple()))

    vals = np.zeros((len(data_vals[versionID]), 4))
    for i, (entry_type, value_dict) in enumerate(data_vals[versionID].items()):
        if not isinstance(value_dict["Run"],int):
            continue
        if run_times[value_dict["Run"]][0] is None:
            continue
        timestamp = (time.mktime(run_times[value_dict["Run"]][0].timetuple()) - ROOT.gStyle.GetTimeOffset())
        vals[i,:] = value_dict["Run"], timestamp, value_dict["Flux2MeV"], value_dict["Flux2MeVErr"]
        #debug
        #formatDatedebug = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(timestamp))
        #print("Run ", value_dict["Run"], " Time: ", formatDatedebug)

    # Trim the vals array so there are no empty rows (due to not having a date)
    mask = vals == 0
    countVals = (~mask).sum(axis=1)
    zeroIndex = np.where(countVals==0)[0]
    # Val columns are run #, timestamp, energy scale, scale_error
    vals = np.delete(vals, zeroIndex, 0)

    return vals, NOMINAL_FLUX2MEV, NOMINAL_FLUX2MEV_ERR

if __name__ == "__main__":
    #TODO make it so that i can have the optional arguement of compareID
    import argparse
    parser = argparse.ArgumentParser(prog="Plot flux to MeV corrections over time given from the corrctions JSON file and compare to other correction values, included (Dodo's corrections, HyoungKu's corrections) or optionally input.", usage='%(prog)s [options]')
    parser.add_argument("version", nargs='?', type=str, help="version ID from correction values' top level key")
    parser.add_argument("compareVersion", type=str, help="version ID to compare input version to")
    args = parser.parse_args()
    versionID = args.version
    compareID = args.compareVersion
    compareID = str(compareID)

    vals, NOMINAL_FLUX2MEV, NOMINAL_FLUX2MEV_ERR = get_json_correction(versionID)

    # debug
    print("Number of value entries, ",vals.shape)
    print("Max & Min Flux2MeV, ", vals[:,2].max(), " ", vals[:,2].min())

    # Histograms for comparing the current version's flux 2 MeV with our other data sets
    eric_comp_h = ROOT.TH1D("Fraction of Eric's Data","Fraction of Eric's Data",500,-1,2)
    dodo_comp_h = ROOT.TH1D("Fraction of Dodo's Data","Fractio of Dodo's Data",500,-1,2)
    hk_comp_h = ROOT.TH1D("Fraction of HK's Data","Fraction of HK'd Data",500,-1,2)
    comp_h = ROOT.TH1D("Fraction of Version %s's Data" % (compareID),"Fraction of Version %s's Data" % (compareID),500,-1,2)

    # Get numerical comparison between the different corrections
    # TODO there has to be a better way to do this without putting in another args variable
    if endpoint_ver == 1:	# endpoint = 53.3 MeV
        older_endpoint_scale = 0.99062

        if comp_endpoint_ver == 0 :
            comp_endpoint_scale = older_endpoint_scale
        else :
            comp_endpoint_scale = 1

    elif endpoint_ver == 0:
        older_endpoint_scale = 1

        if comp_endpoint_ver == 0 :
            comp_endpoint_scale = 1
        else :
            comp_endpoint_scale = 0.99062

    # Fill TGraph data for current version and HyoungKu's data
    gmine = ROOT.TGraphErrors()
    gabs = ROOT.TGraphErrors()
    ghk = ROOT.TGraphErrors()
    gdodo = ROOT.TGraph()
    geric = ROOT.TGraphErrors()
    gdebug = ROOT.TGraphErrors()

    count = 0
    first_time = None
    last_time = None

    # To normalize HK's data with respect to the average flux to mev value
    v_hk_list = []
    err_hk_list = []
    for i, (r, start, v, err) in enumerate(vals):
        v_hk, err_hk = get_hyoungkus_correction(r)
        if v_hk is None:
            continue
        v_hk_list.append(v_hk)
        err_hk_list.append(err_hk)
        
    # Change hk_AV_FLUX2MEV = 1 if we want to normalize with respect to HKs first run
    #hk_AVG_FLUX2MEV = sum(v_hk_list)/float(len(v_hk_list))
    hk_AVG_FLUX2MEV = np.average(v_hk_list, weights=list(map(lambda x: 1.0/x, err_hk_list)))

    for i, (r, start, v, err) in enumerate(vals):
        # debug
        formatDate = time.strftime('%Y-%m-%d', time.gmtime(start))
        if r > 2000 :
            gdebug.SetPoint(count, start, v)
            gdebug.SetPointError(count, 0, err)
            print(r," Start: ",start," ",formatDate)#," Flux2MevScale: ",v," Error: ", err)

        gmine.SetPoint(count, start, v/NOMINAL_FLUX2MEV)
        gmine.SetPointError(count,0, err/NOMINAL_FLUX2MEV)

        gabs.SetPoint(count, start, v)
        gabs.SetPointError(count,0, err)

        r = int(r)
        v_hk, err_hk = get_hyoungkus_correction(r)
        if v_hk is None:
            continue

        ghk.SetPoint(count, start, v_hk/hk_AVG_FLUX2MEV)
        ghk.SetPointError(count, 0, err_hk/hk_AVG_FLUX2MEV)

        hk_comp_h.Fill((v/NOMINAL_FLUX2MEV)/(v_hk/hk_AVG_FLUX2MEV)*older_endpoint_scale)

        count+=1

        # Find the time of the first and last run
        if first_time is None :
            if start is not None :
                first_time = start
            else :
                continue
        else :
            if start is not None :
                min(first_time, start)
            else :
                continue

        if last_time is None :
            if start is not None :
                last_time = start
            else :
                continue
        else :
            if start is not None :
                max(last_time, start)
            else :
                continue

    # Read in corrections for optional JSON version to compare against
    if compareID :
        print("Comparing to : "+compareID)
        gcomp = ROOT.TGraphErrors()
        comp_vals, comp_NOMINAL_FLUX2MEV, comp_NOMINAL_FLUX2MEV_ERR = get_json_correction(compareID)
        count = 0
        for i, (r, start, v, err) in enumerate(comp_vals):
            # debug
            #formatDate = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(start))

            gcomp.SetPoint(count, start, v/comp_NOMINAL_FLUX2MEV)
            gcomp.SetPointError(count,0, err/comp_NOMINAL_FLUX2MEV)

		    # Find corresponding runs flux to MeV value for current version
            current_vals_v = [val_tuple[2] for val_tuple in vals if val_tuple[0]==r]
            if len(current_vals_v) > 1 : 
                print("Error: multiple runs with the same number found wile filling Eric's comparison histogram")
            elif len(current_vals_v) == 1 :
                comp_h.Fill(current_vals_v[0]/v*comp_endpoint_scale)
            else :
                print("Nothing to compare %s's run: %i to. Current flux2mev was found to be " % (compareID, r), current_vals_v) 

    # Clean up the data from Dodo at the start of the code
    dodo_correction = [(get_time(r, srn), v, r) for r, srn, _, _, v in dodo_correction]
    dodo_correction_withRun = {time.mktime(t.timetuple()) - ROOT.gStyle.GetTimeOffset(): (v, r) for t, v, r in dodo_correction if t is not None}	# So I can have the run # to use with comparison histogram
    dodo_correction_forPlots = {time.mktime(t.timetuple()) - ROOT.gStyle.GetTimeOffset(): v for t, v, r in dodo_correction if t is not None}

    for i, (d, vr_tup) in enumerate(sorted(dodo_correction_withRun.items(), key=lambda x: x[0])):
        gdodo.SetPoint(i, d, vr_tup[0])

		# Find corresponding runs flux to MeV value for current version
        current_vals_v = [val_tuple[2] for val_tuple in vals if val_tuple[0]==vr_tup[1]]
        if len(current_vals_v) > 1 : 
            print("Error: multiple runs with the same number found wile filling Eric's comparison histogram")
        elif len(current_vals_v) == 1 :
            dodo_comp_h.Fill((current_vals_v[0]/NOMINAL_FLUX2MEV)/v*older_endpoint_scale)
        else :
            print("Nothing to compare Dodo's run: %i to. Current flux2mev was found to be " % (vr_tup[1]), current_vals_v) 
 
    # Read in Eric's correction values to compare against
    eric_data = []
    with open("/home/marzece/KDAR_Analysis/MichelAnalysis/f2mev_correction/temp.txt", 'r') as efile:
        csvFile = csv.reader(efile)
        for line in csvFile:
            eric_data.append([float(i) for i in line])

    run_times = get_run_dates()

    for i, run_data in enumerate(eric_data):
        eric_timestamp = (time.mktime(run_times[run_data[0]][0].timetuple()) - ROOT.gStyle.GetTimeOffset())
        eric_data[i].append(eric_timestamp)
        # Run	Flux2MeV	Err		Time

    # Get the average flux 2 mev value 
    eric_data_list = zip(*eric_data)
    v_eric_data = list(eric_data_list[1])
    err_eric_data = list(eric_data_list[2])
#    eric_AVG_FLUX2MEV = sum(v_eric_list)/float(len(v_eric_list))
    eric_AVG_FLUX2MEV = np.average(v_eric_data, weights=list(map(lambda x: 1.0/x, err_eric_data)))

    for i, (r, v, err, start) in enumerate(eric_data):
        formatDate = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(start))
        geric.SetPoint(i, start, v/eric_AVG_FLUX2MEV)
        geric.SetPointError(i,0, err/eric_AVG_FLUX2MEV)

		# Find corresponding runs flux to MeV value for current version
        current_vals_v = [val_tuple[2] for val_tuple in vals if val_tuple[0]==r]
        if len(current_vals_v) > 1 : 
            print("Error: multiple runs with the same number found wile filling Eric's comparison histogram")
        elif len(current_vals_v) == 1:
            eric_comp_h.Fill(current_vals_v[0]/v*older_endpoint_scale)
        else :
            print("Nothing to compare Eric's run: %i to. Current flux2mev was found to be " % (vr_tup[1]), current_vals_v) 

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

    # Get the trend line of current version correction data
    gtrend = ROOT.TGraphErrors()
    gtrend2 = ROOT.TGraphErrors()

    trend = get_trend_line_simple(vals[:, 1:4])
    for i, (x, y, err) in enumerate(trend):
        gtrend.SetPoint(i,x,y)
        gtrend.SetPointError(i,0,err)
        gtrend2.SetPoint(i,x,y/NOMINAL_FLUX2MEV)
        gtrend2.SetPointError(i,0,err/NOMINAL_FLUX2MEV)
        #gtrend2.SetPoint(i,x,y)
        #gtrend2.SetPointError(i,0,err)
    gtrend.SetLineColor(ROOT.kBlue)
    gtrend.SetFillColorAlpha(ROOT.kBlue, ALPHA)
    gtrend.SetFillStyle(3001)
    gtrend2.SetLineColor(ROOT.kBlue)
    gtrend2.SetFillColorAlpha(ROOT.kBlue, ALPHA)
    gtrend2.SetFillStyle(3001)

    #debug
    c555 = ROOT.TCanvas("debug2000")
    gdebug.Draw("AP*")
    c555.Write()
    #c666 = ROOT.TCanvas("debug2")
    #mgdebug = ROOT.TMultigraph()
    #gdebug.SetMarkerColor(2)
    #mgdebug.Add(gdebug)
    #mgdebug.Add(gabs)
    #mgdebug.Draw("APE*")
    #c666.Update()
    #c666.write()

    # Plot the current version corrections flux to MeV and its trend line
    c1 = ROOT.TCanvas()
    mg0 = ROOT.TMultiGraph()
    mg0.Add(gabs)
    mg0.Add(gdebug)
    mg0.Add(gtrend, "lineE3")
    mg0.Draw("APE*")
    mg0.GetXaxis().SetTimeDisplay(1)
    mg0.GetYaxis().SetTitle("Michel Flux to MeV")
    mg0.GetXaxis().SetTimeFormat("%Y/%m/%d")
    leg0 = ROOT.TLegend(0.1, 0.66, 0.45, 0.9)
    leg0.AddEntry(gabs, "%s" % (versionID), "PE")
    leg0.Draw('same')
    c1.Update()
    c1.Write()

    # 
    gnominal = ROOT.TGraphErrors()
    gnominal.SetPoint(0, first_time, 1.0)
    gnominal.SetPoint(1, last_time, 1.0)
    gnominal.SetPointError(0, 0, NOMINAL_FLUX2MEV_ERR/NOMINAL_FLUX2MEV)
    gnominal.SetPointError(1, 0, NOMINAL_FLUX2MEV_ERR/NOMINAL_FLUX2MEV)
    gnominal.SetFillColorAlpha(ROOT.kGray, ALPHA)
    gnominal.SetLineColor(ROOT.kGray)
    ghk.SetMarkerColor(2)
    ghk.SetLineColor(2)
    c2 = ROOT.TCanvas()
    mg = ROOT.TMultiGraph()
    mg.Add(gmine)
    mg.Add(ghk)
    mg.Add(gnominal, "E3")
    mg.Add(gtrend2, "lineE3")
    mg.Draw("APE*")
    mg.GetXaxis().SetTimeDisplay(1)
    mg.GetXaxis().SetTimeFormat("%Y/%m/%d")
    mg.GetYaxis().SetTitle("Relative Flux per MeV")
    leg2 = ROOT.TLegend(0.1, 0.66, 0.45, 0.9)
    leg2.AddEntry(gmine, "New measurement, %s" % (versionID), "pe")
    leg2.AddEntry(ghk, "Hyoungku's measurement", "pe")
    leg2.AddEntry(gnominal, "Flux-to-MeV = %0.3f #pm %0.3f" % (NOMINAL_FLUX2MEV, NOMINAL_FLUX2MEV_ERR),  "f")
    leg2.Draw('same')
    c2.Update()
    c2.Write()

    c3 = ROOT.TCanvas()
    gdodo.SetMarkerColor(ROOT.kGreen)
    gdodo.SetLineColor(ROOT.kGreen)
    mg3 = ROOT.TMultiGraph()
    mg3.Add(gdodo, "P*")
    mg3.Add(gmine.Clone(), "P*")
    #mg3.Add(gtrend2.Clone(), "lineE3")
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
    mg4.Add(gtrend2.Clone(), "lineE3")

    mg4.Add(gdodo.Clone(), "P*")
    mg4.Add(gmine.Clone(), "P*")

    mg4.Draw("A")
    mg4.GetYaxis().SetTitle("Relative Flux per MeV")
    mg4.GetXaxis().SetTimeDisplay(1);
    mg4.GetXaxis().SetTimeFormat("%Y/%m/%d")
    leg4 = ROOT.TLegend(0.1, 0.66, 0.45, 0.9)
    leg4.AddEntry(gdodo, "Dodo's measurement (normalized by first run)", "pe")
    leg4.AddEntry(gmine, "New measurement, %s" % (versionID), "pe")
    leg2.AddEntry(gtrend2, "Flux-to-MeV = %0.3f #pm %0.3f" % (NOMINAL_FLUX2MEV, NOMINAL_FLUX2MEV_ERR),  "f")
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
    mg5.Add(gtrend2.Clone(), "lineE3")
    mg5.Add(g_combined_trend, "lineE3")
    mg5.Draw("A")
    mg5.GetYaxis().SetTitle("Relative Flux per MeV")
    mg5.GetXaxis().SetTimeDisplay(1);
    mg5.GetXaxis().SetTimeFormat("%m/%d")
    leg5 = ROOT.TLegend(0.14, 0.66, 0.44, 0.88)
    leg5.AddEntry(g_dodo_trend, "Dodo's nGd measurement (normalized by first run)", "LF")
    leg5.AddEntry(gtrend2, "New Michel measurement", "LF")
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
    mg9.Add(gtrend2.Clone(), "lineE3")
    mg9.Draw("A")
    mg9.GetYaxis().SetTitle("Relative Flux per MeV")
    mg9.GetXaxis().SetTimeDisplay(1);
    mg9.GetXaxis().SetTimeFormat("%Y/%m/%d")
    leg9 = ROOT.TLegend(0.14, 0.66, 0.44, 0.88)
    leg9.AddEntry(g_eric_trend, "Eric's measurement", "LF")
    leg9.AddEntry(gtrend2, "New Michel measurement", "LF")
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

