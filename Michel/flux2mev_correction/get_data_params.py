"""
	Originally from Alex on McCaffrey at /home/marzece/KDAR_Analysis/MichelAnalysis/AlexMichel/
	
	Takes the parsed michel data from updateMichelPair and returns an addition to the json file
	correction_values.json as well as individual ROOT files with the fits.

	On running, the user must declare if it will find the fit for
	0: all the data in michel_pair_combined.root
	1: Each run
	2: Each subrun (not implemented)
"""

import ROOT
import numpy as np
from scipy.optimize import curve_fit
import fitting as fit
from array import array
from scipy.linalg import cholesky
import os
import csv
from datetime import datetime
import calendar
import time
from collections import defaultdict
import json

OUTPUT_PATH = "/home/littleca/kdar/Michel/flux2mev_correction/output_fluxCorr/"
JSON_NAME="/home/littleca/kdar/correction_values.json"

MICHEL_ENDPOINT_ENERGY = 53.3 # MeV
# Fiducial volume and energy ROI selection values
FV_R_MAX = 1400.0 # milli-meters
FV_Z_MAX = 1000.0 # milli-meters
FV_Z_MIN =-1000.0 # milli-meters
AV_R = 1600.0
AV_Z_MAX = 1250.0 
AV_Z_MIN =-1250.0
PMT_Z_MAX = 1470 
PMT_Z_ = -1470
R_PMT = 1850.0
MICHEL_EMIN_CUT = 20 # MeV
MICHEL_EMAX_CUT = 60 # MeV
MULTISIM_NTRIALS = int(1e5) # number of parameterizations
#mc_michel = fit.MCMichelSpectrum()
delta_t_fitf = ROOT.TF1("", "(x>[0])*([1]*exp(-(x-[0])/[2]) + [3])", 0, 10)


def fit_delta_t(h):
    delta_t_fitf.SetParameter(0, 2.0)
    delta_t_fitf.SetParameter(1, h.GetMaximum())
    delta_t_fitf.SetParameter(2, 2.2)
    delta_t_fitf.SetParameter(3, h.GetBinContent(h.GetNbinsX()))
    h.Fit(delta_t_fitf, "S")
    return delta_t_fitf.Clone()

def get_run_dates():
    vals = list(csv.reader(open("/home/littleca/kdar/Michel/flux2mev_correction/run_dates.txt")))
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

def get_hyoungkus_correction(run):
    path = "/home/littleca/kdar/Michel/flux2mev_correction/HyoungKu_correction/"
    fn = "Run_%i_timeCorr.dat" % run
    if os.path.isfile(os.path.join(path, fn)) :
        fin = open(os.path.join(path, fn))
        lines = fin.readlines()
        return [float(x) for x in lines[0].strip().split('\t')][2:]
    else :
        return None

class RunBasedHistogram():
    def __init__(self, name, title, nz, zlow, zhigh):
        def new_th1d():
            return ROOT.TH1D("", "", nz, zlow, zhigh)
        self.run_hists = defaultdict(new_th1d)
        self.graph = ROOT.TGraphErrors()
        self.run2date = get_run_dates()
        self.fit_vals = {}
    def Fill(self, run, zval):
        self.run_hists[run].Fill(zval)
    def DoFit(self):
        for r,h in self.run_hists.items():
            if h.Integral() < 1e3:
                continue
            E_vals, N_vals, err = fit.GetData(h)
            guess = [0.04, h.GetMaximum(), MICHEL_ENDPOINT_ENERGY*FLUX2MEV]
            fv, fv_err = curve_fit(fit.Michel_Conv, E_vals, N_vals, p0=guess, sigma=err)
            f2mev = fv[2]/MICHEL_ENDPOINT_ENERGY
            f2mev_err = np.sqrt(fv_err[2][2])/MICHEL_ENDPOINT_ENERGY
            self.fit_vals[r] =  (f2mev, f2mev_err)

    def Draw(self, scale_to_t0=True):
        t0 = None
        tf = None
        first_val = None
        self.hk_graph = ROOT.TGraphErrors()
        count = 0
        for i, (r, fv) in enumerate(sorted(self.fit_vals.items(), key=lambda x: x[0])):
            date = self.run2date[r]
            if(date[0] is None):
                continue
            date = time.mktime(date[0].timetuple())-ROOT.gStyle.GetTimeOffset()
            t0 = min(t0, date) if t0 else date
            tf = max(tf, date) if tf else date
            #inverse_val = 1.0/fv[0][0]
            #fractional_error = fv[1][0]/fv[0][0]
            #value = inverse_val

            value = fv[0]
            fractional_error = fv[1]/value


            first_val = value if (first_val is None) else first_val
            if(fractional_error*value < 5.0): 
                plot_value = value/first_val;
                self.graph.SetPoint(count, date, plot_value)
                self.graph.SetPointError(count, 0, plot_value*fractional_error)
                count+=1

        for i, r in enumerate(sorted(self.fit_vals.keys())):
            hk_val = get_hyoungkus_correction(r)
            if hk_val is None:
                continue
            date = self.run2date[r]
            if(date[0] is None):
                continue
            date = time.mktime(date[0].timetuple())-ROOT.gStyle.GetTimeOffset()
            self.hk_graph.SetPoint(i, date, hk_val[0])
            self.hk_graph.SetPointError(i, 0, hk_val[1])

        self.hk_graph.GetXaxis().SetTimeDisplay(1)
        self.hk_graph.GetXaxis().SetTimeFormat("%m/%d")
        self.hk_graph.GetYaxis().SetTitle("Relative Flux-To-MeV")
        self.hk_graph.SetMarkerColor(2)
        self.hk_graph.SetLineColor(2)
        self.graph.GetXaxis().SetTimeDisplay(1)
        self.graph.GetXaxis().SetTimeFormat("%m/%d")
        self.graph.GetYaxis().SetTitle("Relative Flux-To-MeV")
        self.tline = ROOT.TLine(t0, 1, tf, 1)
        self.tline.SetLineColor(2)
        self.graph.Draw('APE')
        self.hk_graph.Draw('samePE')
        self.tline.Draw('same')


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str, help="input ROOT filename (update_michel_pair.cc output; michel_pair_{run#}.root if run type = 1 or michel_pair_combined.root if run type = 0)")
    parser.add_argument("runType", type=int, help="0: Sum, 1: Run, 2: Subrun")
    parser.add_argument("version", type=str, help="Name to be used as version ID")
    args = parser.parse_args()
    inputName = args.input
    RUN_TYPE = args.runType
    versionID = args.version

    fid_max = 70
    fid_min = 0
    NBINS_R = 7
    NBINS_Z = 15
    Z_MIN = -1.5e3
    Z_MAX = 2.25e3
    R_MAX = 1.85e3**2
    #flux_min = 10e3 # FOR MDAQ 
    #flux_max = 60e3 # FOR MDAQ 
    flux_min = 15e3 # FOR JADE V3
    flux_max = 60e3 # FOR JADE V3
    #flux_min = 12e3 # FOR JADE V0
    #flux_max = 48e3 # FOR JADE V0

    DELT_MIN =  2.0 # micro-seconds
    KICKER_DEAD_TIME = 5000 # ns

    pos_hist = ROOT.TH2D("", ";#rho^{2}/#rho_{PMT}^{2};Z [m]", 100, 0, 1.1, 100, -1.5, 1.5)
    
    #FLUX2MEV = 2*392.1 # MDAQ ONLY
    FLUX2MEV = 960 #956.0776 # From earlier michel calibration
    if (RUN_TYPE==1):
        # Get the Flux2MeV from the json file
        with open(JSON_NAME, "r") as f:
           data_vals = json.load(f)

        # Check that this version has a sum correction
        # If so, use it for the FLUX2MEV value
        # If not, prompt the use to confirm the use of the default value
        sumlist = []
        if versionID in data_vals.keys():
            # Get the name of the "sum"key for this version
            # There should be only be one value for sum
            sumlist = [key for key, value in data_vals[versionID].items() if 'sum' in key]
            if len(sumlist) > 1 :
                print("ERROR: Multiple sum entries (%i) present for version %s" %(len(sumlist),  versionID))
                exit()

        if len(sumlist) == 1 :
            FLUX2MEV = data_vals[versionID][sumlist[0]]["Flux2MeV"]

        else :
            while True :
                cont_prompt = input("Version \"%s\" does not have a sum correction value, continue using the default flux to MeV of %i? [Y/n]" % (versionID, FLUX2MEV))
                valid = {"" : True, "Yes" : True, "Y" : True, "y" : True, "yes" : True, "No" : False, "N" : False, "n" : False, "no" : False}
                if cont_prompt in valid.keys() :
                    if valid[cont_prompt] :
                        # Check if we have a directory to save the output to
                        newDir = OUTPUT_PATH+str(versionID)
                        if os.path.isdir(newDir) :
                            print("Found output directory %s" % (newDir))
                            break
                        else :
                            os.makedirs(newDir) 
                            print("Created output directory %s" % (newDir))
                            break
                    else :
                        exit()
                else :
                    print("Please respond with 'yes' or 'no' or press enter to continue")

    location_flux_hist = fit.PositionBasedHistograms("", "", NBINS_R, 0, R_MAX, NBINS_Z, Z_MIN, Z_MAX, 100, flux_min, flux_max)
    location_deltat_hist = fit.PositionBasedHistograms("", "", NBINS_R, 0, R_MAX, NBINS_Z, Z_MIN, Z_MAX, 100, 2, 10)
    delt_hist = ROOT.TH1D("", "", 200, 0, 10)


    # Grab the Data to analyze
    print("Opening %s" % (inputName))
    data_file = ROOT.TFile.Open(inputName, "READ")
    energy_tree = data_file["event_tree"]
    n_entries = energy_tree.GetEntries()
    data = np.zeros(n_entries, dtype=[("x", float), ("y", float), ("z", float),
                                      ("flux", float), ("nsat", int), ("delt", float), ("delvtx", float),
                                      ('run', int), ("subrun", int), ('entry', int)])
    flux_list = []
    for i in range(n_entries):
        energy_tree.GetEntry(i)
        #if(energy_tree.time_since_kicker < KICKER_DEAD_TIME):
        #    continue
        if(energy_tree.delt < DELT_MIN):
            continue
        # Event positions are stored as meters, convert the to mm
        x = energy_tree.x*1000.
        y = energy_tree.y*1000.
        z = energy_tree.z*1000.
        R = np.sqrt(x**2 + y**2)
        del_vtx = energy_tree.delvtx

        data[i]['x'] =  x
        data[i]['y'] =  y
        data[i]['z'] =  z
        data[i]['flux'] =  energy_tree.flux
        data[i]['nsat'] =  energy_tree.nsat
        data[i]['delt'] =  energy_tree.delt
        data[i]['delvtx'] =  del_vtx

        data[i]['run'] =  energy_tree.run
        data[i]['subrun'] =  energy_tree.subrun
        data[i]['entry'] =  energy_tree.index
        location_deltat_hist.Fill(R**2, z, energy_tree.delt)

        if z < FV_Z_MAX and z > FV_Z_MIN and R < FV_R_MAX:
            flux_list.append(energy_tree.flux)

    flux_list = np.array(flux_list)
    flux_list_slice = np.array([i for i in flux_list if i < flux_max and i > flux_min])

    # Compute the number of bins
    n_bins_flux = len(np.histogram_bin_edges(flux_list_slice, bins='rice')) - 1
    #n_bins_flux = 100
    set_hist = {"n_bins" : n_bins_flux, "min" : flux_min, "max" : flux_max}

    flux_hist = ROOT.TH1D("flux_hist", "Michel Flux;Flux;Counts", n_bins_flux, flux_min, flux_max)

    for flux in flux_list:
        flux_hist.Fill(flux)

    # Fit the data for y-axis scaling (unused), flux to MeV energy scaling, and energy resolution scaling factor
    flux_fit_vals, flux_errors, flux_fit_graph = fit.do_edep_fit(flux_list, set_hist)

    conv = flux_fit_vals[1]
    conv_err = np.sqrt(flux_errors[1][1])
    res_scale = flux_fit_vals[2]
    res_scale_err = np.sqrt(flux_errors[2][2])

    print("File = %s, Conversion = %f +- %f" % (inputName, conv, conv_err))
    flux_fit_graph.SetLineColor(2)
    flux_fit_graph.SetLineWidth(2)

    time_e_hist = RunBasedHistogram("", "", 100, flux_min/conv, flux_max/conv)

    fid_list = []
    # Go through events again, convert their energy to MeV and apply michel cuts
    for i in range(n_entries):
        energy_tree.GetEntry(i)
        #if(energy_tree.time_since_kicker < KICKER_DEAD_TIME):
        #    continue
        if(energy_tree.delt < DELT_MIN):
            continue
        x = energy_tree.x*1000.
        y = energy_tree.y*1000.
        z = energy_tree.z*1000.
        R = np.sqrt(x**2 + y**2)
        E = energy_tree.flux/conv
        #rho_sqrd = R**2/R_PMT**2
        rho_sqrd = R**2
        location_flux_hist.Fill(rho_sqrd, z, energy_tree.flux)

        if(E <MICHEL_ENDPOINT_ENERGY and E>MICHEL_EMIN_CUT):
            pos_hist.Fill(rho_sqrd/R_PMT**2, z/1000.)
        if E<MICHEL_EMAX_CUT and E > MICHEL_EMIN_CUT and z < FV_Z_MAX and z > FV_Z_MIN and R < FV_R_MAX:
            time_flux_hist.Fill(energy_tree.run, energy_tree.flux)
            delt_hist.Fill(energy_tree.delt)
            fid_list.append(E)
     
    fid_list = np.array(fid_list)
    fid_list_slice = np.array([i for i in fid_list if i < fid_max and i > fid_min])
    
    # Compute the number of bins
    n_bins_fid = len(np.histogram_bin_edges(fid_list_slice, bins='rice')) - 1
    #n_bins_fid = 100
    fid_hist = ROOT.TH1D("fid_hist", "fid_hist", n_bins_fid, fid_min, fid_max)

    for E in fid_list:
        fid_hist.Fill(E)

    # Write fit values to json file
    output_run=0
    if (RUN_TYPE==0) :
        # Looking at a combined number of runs
        output_type = "sum"
        output_run = None
        output_run_name = output_type
        output_subrun = None
        # Make output directory for that version ID
        newDir = OUTPUT_PATH+str(versionID)
        os.makedirs(newDir) 
        outFile = ROOT.TFile(newDir+"/fluxCorr_combined.root", "RECREATE")
    elif (RUN_TYPE==1) :
        # Looking on per-run basis
        output_type = "run"
        output_run = int(data[0]['run'])
        output_run_name = output_type+"_"+str(output_run)
        output_subrun = None
        outFile = ROOT.TFile(OUTPUT_PATH+str(versionID)+"/fluxCorr_"+str(output_run)+".root", "RECREATE")
    elif (RUN_TYPE==3) :
        # Looking on per-subrun basis
        output_type = "subrun"
        output_run = int(data[0]['run'])
        output_run_name = output_type+"_"+str(output_run)
        output_subrun = int(dat[0]['subrun'])
        outFile = ROOT.TFile(OUTPUT_PATH+str(versionID)+"/fluxCorr_"+str(output_run)+"_"+str(output_subrun)+".root", "RECREATE")
    
    out_vals = {output_run_name : {"Run" : output_run, "Subrun" : output_subrun, "StartingFlux2MeV": FLUX2MEV, "Flux2MeV" : conv, "Flux2MeVErr" : conv_err, "Flux2MeVScale" : None, "Flux2MeVScaleErr" : None}}
    version_update = {str(versionID) : out_vals}

    if (os.path.isfile(JSON_NAME)) :
        with open(os.path.abspath(JSON_NAME), "r") as jfile:
            json_data = json.load(jfile)
        if (str(versionID) in json_data) :
            json_data[str(versionID)].update(out_vals)
        else :
            json_data.update(version_update)

        with open(os.path.abspath(JSON_NAME), "w") as jfile:
            json.dump(json_data, jfile, indent=4)
    else:
        with open(os.path.abspath(JSON_NAME), "w") as jfile:
            json.dump(version_update, jfile, indent=4)

    outFile.cd()
    c0 = ROOT.TCanvas("MC Conv")
    flux_hist.GetXaxis().SetTitle("Flux")
    flux_hist.GetYaxis().SetTitle("Counts")
    flux_hist.GetYaxis().SetRangeUser(0, 1.2*flux_hist.GetMaximum())
    flux_hist.SetTitle("")
    flux_hist.SetStats(0)
    flux_hist.SetLineColor(ROOT.kBlack)
    flux_hist.SetLineWidth(2)
    flux_fit_graph.Draw("al")
    flux_hist.Draw("SAMEE")
    latex = ROOT.TLatex()
    latex.DrawLatex(1000, 0.8*flux_hist.GetMaximum(), "Flux/MeV: %0.3f #pm %0.3f" % (conv, conv_err))
    latex.DrawLatex(5, 0.8*fid_hist.GetMaximum(), "R_{ep} Scale = %0.4f%% #pm %0.4f" % (res_scale*100, res_scale_err*100))
    c0.Write()

    C = ROOT.TCanvas("Conv Michel E")
    fid_hist.GetXaxis().SetTitle("Reconstructed Energy [MeV]")
    fid_hist.GetYaxis().SetTitle("Counts/%0.1fMeV" % fid_hist.GetBinWidth(1))
    fid_hist.GetYaxis().SetRangeUser(0, 1.5*fid_hist.GetMaximum())
    fid_hist.SetTitle("")
    fid_hist.SetStats(0)
    fid_hist.SetLineColor(ROOT.kBlack)
    fid_hist.SetLineWidth(2)
    fid_hist.Draw("SAMEE")

    leg = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    leg.AddEntry(fid_hist, "Sym Fid: ( #mu^{+} + #mu^{-} )")
   # leg.Draw("Same")

    latex = ROOT.TLatex()
    latex.DrawLatex(5, 0.9*fid_hist.GetMaximum(), "Flux/MeV = %0.4f #pm %0.4f" % (conv, conv_err))
    latex.DrawLatex(5, 0.8*fid_hist.GetMaximum(), "R_{ep} Scale = %0.4f%% #pm %0.4f" % (res_scale*100, res_scale_err*100))
    C.Write()

    c2 = ROOT.TCanvas("MC Spatial")
    c2.SetLogz()
    pos_hist.SetStats(0)
    pos_hist.Draw("colz")
    # Draw the fiducial volume
    ll1 = ROOT.TLine(0, FV_Z_MAX/1e3, FV_R_MAX**2/R_PMT**2, FV_Z_MAX/1e3) 
    ll2 = ROOT.TLine(FV_R_MAX**2/R_PMT**2, FV_Z_MAX/1e3, FV_R_MAX**2/R_PMT**2, FV_Z_MIN/1e3) 
    ll3 = ROOT.TLine(FV_R_MAX**2/R_PMT**2, FV_Z_MIN/1e3, 0, FV_Z_MIN/1e3) 
    lines = [ll1, ll2, ll3]
    [ll.SetLineWidth(2) for ll in lines]
    [ll.SetLineColor(2) for ll in lines]
    [ll.Draw('same') for ll in lines]

    # Draw the acrylic vessel
    av_ll1 = ROOT.TLine(0, AV_Z_MAX/1e3, AV_R**2/R_PMT**2, AV_Z_MAX/1e3) 
    av_ll2 = ROOT.TLine(AV_R**2/R_PMT**2, AV_Z_MAX/1e3, AV_R**2/R_PMT**2, AV_Z_MIN/1e3) 
    av_ll3 = ROOT.TLine(AV_R**2/R_PMT**2, AV_Z_MIN/1e3, 0, AV_Z_MIN/1e3) 
    lines = [av_ll1, av_ll2, av_ll3]
    [ll.SetLineWidth(2) for ll in lines]
    [ll.SetLineColor(ROOT.kBlue) for ll in lines]
    [ll.Draw("same") for ll in lines]
    
    c2.Write()

    c3 = ROOT.TCanvas("DeltatHist")
    delt_hist.Draw()
    c3.Write()

    outFile.Close()

#    c4 = ROOT.TCanvas()
#    c4.SetRightMargin(0.0)
#    c4.SetLeftMargin(0.0)
#    c4.SetTopMargin(0.0)
#    c4.SetBottomMargin(0.0)
#    location_flux_hist.Draw(c4, fitf=lambda h: fit.do_michel_fit(h)[2])
#
#    ROOT.gStyle.SetPaintTextFormat("0.4g")
#    location_flux_hist.Fit()
#    c5 = ROOT.TCanvas()
#    c6 = ROOT.TCanvas()
#    location_flux_hist.DrawFit(c5, c6)

#    c9 = ROOT.TCanvas("TimeFluxHist")
#    time_e_hist.DoFit()
#    time_e_hist.Draw()
#    c9.Write()

#    c10 = ROOT.TCanvas()
#    location_deltat_hist.Draw()
#    c10.Write()

    # Create Multisim parameter sets for Systematics
    # Decompose the cov matrix for sampling
    #cky = cholesky(fit_errors, lower=True)
    #mc_smear_param_tree = ROOT.TNtuple("data_param_tree", "data_param_tree", "scale:res:ampl")
    #gaus = np.random.normal(loc=0, scale=1, size=(MULTISIM_NTRIALS, 3))
    #for random_vec in gaus:
        # TODO there's probably a more numpy way of calculating this
    #    p = fit_vals + np.dot(cky, random_vec)
    #    mc_smear_param_tree.Fill(p[0], p[1], p[2])
    #mc_smear_param_tree.Write()
    #outFile.Close()

