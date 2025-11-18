import ROOT
import numpy as np
import fitting as fit
import os
import csv

# The point of this script is to
def get_time_correction_data():
    data_path = "/home/marzece/KDAR_Analysis/KDARSelection/dat"
    files = [x for x in os.listdir(data_path) if x[-4:] == ".dat"]
    result = {}
    for fn in files:
        with open(os.path.join(data_path, fn)) as fin:
            data = list(csv.reader(fin, delimiter="\t"))
        data = {(int(rn), int(srn)): (float(v), float(err)) for rn, srn, v, err in data}
        result.update(data)
    return result

def apply_flux_to_mev_correction(run, subrun):
    NOMINAL_DATA_FLUX2MEV = 942.183
    try:
        val = apply_flux_to_mev_correction.data_table[run, subrun][0]
    except AttributeError:
        apply_flux_to_mev_correction.data_table = get_time_correction_data()
        val = apply_flux_to_mev_correction.data_table[run, subrun][0]
    except KeyError:
        val = apply_flux_to_mev_correction.data_table[run, 0][0]
    return NOMINAL_DATA_FLUX2MEV*val


if __name__ == "__main__":
    ROOT.TH1.AddDirectory(0)
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("michel_input", type=str, help="input ROOT filename for MC Michel data")
    parser.add_argument("data_input", type=str, help="input ROOT filename for (real) Michel data")
    parser.add_argument("output", type=str, help="output ROOT filename")
    args = parser.parse_args()
    mc_input_file = args.michel_input
    data_input_file = args.data_input

    n_bins = 150
    canvii = set()

    fid_fit_energy_hist = ROOT.TH1D("fid_energy_fit", "MC Michel Energy;Energy [MeV];Counts", n_bins, 20, 60)
    fid_fit_energy_hist.Sumw2()

    mc_location_flux_hist = fit.PositionBasedHistograms("MC", "MC", 12, 0, 2.56e6, 10, -1.25e3, 1.25e3, 100, 20, 60)
    data_location_flux_hist = fit.PositionBasedHistograms("Data", "Data", 12, 0, 2.56e6, 10, -1.25e3, 1.25e3, 100, 20, 60)

    michel_data = ROOT.TFile(mc_input_file, "READ")

    tree = michel_data.energy_tree
    michel_data_dtype = np.dtype([('pos', np.float, (3,)), ('energy', np.float), ('weight', np.float)])
    michel_event_data = np.zeros(tree.GetEntries(), dtype=michel_data_dtype)

    for i, event in enumerate(tree):
        E = event.E
        R = np.sqrt(event.x**2 + event.y**2)
        Rs_reco = R**2
        michel_event_data[i] = ((event.x, event.y, event.z), E, event.weight)
        mc_location_flux_hist.Fill(Rs_reco, event.z, E, event.weight)

    data_file = ROOT.TFile.Open(data_input_file, "READ")
    real_data_tree = data_file.event_tree
    n_entries = real_data_tree.GetEntries()
    data = np.zeros(n_entries, dtype=[("x", np.float), ("y", np.float), ("z", np.float), ("energy", np.float)])

    for i in range(n_entries):
        real_data_tree.GetEntry(i)
        # Event positions are stored as meters, convert the to mm
        x = real_data_tree.x*1000.
        y = real_data_tree.y*1000.
        z = real_data_tree.z*1000.

        flux =  real_data_tree.flux
        run =  real_data_tree.run
        subrun=  real_data_tree.subrun
        f2mev = apply_flux_to_mev_correction(run, subrun)
        energy = flux/f2mev

        data[i]['x'] =  x
        data[i]['y'] =  y
        data[i]['z'] =  z
        data[i]['energy'] = energy
        data_location_flux_hist.Fill(x**2+y**2, z, energy)

    outFile = ROOT.TFile(args.output, "RECREATE")
    outFile.cd()


    c = ROOT.TCanvas()
    canvii.add(c)
    c.SetWindowSize(c.GetWw()*2, c.GetWh()*2)
    c.SetRightMargin(0.0)
    c.SetLeftMargin(0.0)
    c.SetTopMargin(0.0)
    c.SetBottomMargin(0.0)
    mc_location_flux_hist.Draw(c, fitf=lambda h: fit.do_michel_fit(h)[2])

    ROOT.gStyle.SetPaintTextFormat("0.4g")
    mc_location_flux_hist.Fit()
    c5 = ROOT.TCanvas()
    c6 = ROOT.TCanvas()
    canvii.add(c5)
    canvii.add(c6)
    mc_location_flux_hist.DrawFit(c5, c6)

    mc_location_flux_hist.scale_h2.SetTitle("MC, Energy Scale")
    mc_location_flux_hist.resolution_h2.SetTitle("MC, Energy Resolution")
    mc_location_flux_hist.scale_h2.Write()
    mc_location_flux_hist.resolution_h2.Write()

    c = ROOT.TCanvas()
    canvii.add(c)
    c.SetWindowSize(c.GetWw()*2, c.GetWh()*2)
    c.SetRightMargin(0.0)
    c.SetLeftMargin(0.0)
    c.SetTopMargin(0.0)
    c.SetBottomMargin(0.0)
    data_location_flux_hist.Draw(c, fitf=lambda h: fit.do_michel_fit(h)[2])

    ROOT.gStyle.SetPaintTextFormat("0.4g")
    data_location_flux_hist.Fit()
    c5 = ROOT.TCanvas()
    c6 = ROOT.TCanvas()
    canvii.add(c5)
    canvii.add(c6)
    data_location_flux_hist.DrawFit(c5, c6)
    data_location_flux_hist.scale_h2.Write()
    data_location_flux_hist.scale_h2.SetTitle("Data, Energy Scale")
    data_location_flux_hist.resolution_h2.SetTitle("Data, Energy Resolution")
    data_location_flux_hist.resolution_h2.Write()
    outFile.Close()


    me_vs_mc = np.array([[data_location_flux_hist.scale_h2.GetBinContent(i+1, j+1) - mc_location_flux_hist.scale_h2.GetBinContent(i+1, j+1) for i in range(12)] for j in range(10)])
    me_vs_mc_err = np.array([[data_location_flux_hist.scale_h2.GetBinError(i+1, j+1)**2 + mc_location_flux_hist.scale_h2.GetBinError(i+1, j+1)**2 for i in range(12)] for j in range(10)])
    h2_comp = ROOT.TH2D("", "", 12, 0 ,2.56e6, 10, -1.25e3, 1.25e3)
    fv_mean_diff = np.mean(me_vs_mc[1:-1, :-3])
    for j in range(10):
        for i in range(12):
            h2_comp.SetBinContent(i+1, j+1, (me_vs_mc[j,i] - fv_mean_diff)*100)
            h2_comp.SetBinError(i+1, j+1, me_vs_mc_err[j,i]*100)

    c99 = ROOT.TCanvas()
    h2_comp.SetStats(0)
    ROOT.gStyle.SetPaintTextFormat("0.3g")
    c99.SetWindowSize( c99.GetWw()*2, c99.GetWh()*2)
    ll1 = ROOT.TLine(0, 1e3, 1400**2, 1e3)
    ll2 = ROOT.TLine(1400**2, 1e3, 1400**2, -1e3)
    ll3 = ROOT.TLine(1400**2, -1e3, 0, -1e3)
    [ll.SetLineWidth(2) for ll in [ll1, ll2, ll3]]
    [ll.SetLineColor(2) for ll in [ll1, ll2, ll3]]
    h2_comp.GetXaxis().SetTitle("R^{2} [mm]")
    h2_comp.GetYaxis().SetTitle("Z [mm]")
    h2_comp.Draw("colztextE")
    [ll.Draw('same') for ll in [ll1, ll2, ll3]]
    h2_comp.GetXaxis().SetTitle("R^{2} [mm]")
    h2_comp.GetYaxis().SetTitle("Z [mm]")

    c98 = ROOT.TCanvas()
    me_vs_mc_res = np.array([[data_location_flux_hist.resolution_h2.GetBinContent(i+1, j+1) - mc_location_flux_hist.resolution_h2.GetBinContent(i+1, j+1) for i in range(12)] for j in range(10)])
    me_vs_mc_res_err = np.array([[data_location_flux_hist.resolution_h2.GetBinError(i+1, j+1)**2 + mc_location_flux_hist.resolution_h2.GetBinError(i+1, j+1)**2 for i in range(12)] for j in range(10)])
    h2_comp_res = ROOT.TH2D("", "", 12, 0 ,2.56e6, 10, -1.25e3, 1.25e3)
    for j in range(10):
        for i in range(12):
            h2_comp_res.SetBinContent(i+1, j+1, me_vs_mc_res[j,i])
            h2_comp_res.SetBinError(i+1, j+1, me_vs_mc_res_err[j,i])
    h2_comp_res.Draw('colztexte')
    h2_comp_res.SetStats(0)
    h2_comp_res.GetXaxis().SetTitle("R^{2} [mm]")
    h2_comp_res.GetYaxis().SetTitle("Z [mm]")
    [ll.Draw('same') for ll in [ll1, ll2, ll3]]

    escale_rms_fv = np.sqrt(np.mean(me_vs_mc[1:-1, :-3]**2))
    eres_rms_fv = np.sqrt(np.mean(me_vs_mc_res[1:-1, :-3]**2))
    print("FV RMS Discrepancy - Energy Scale: %f%%" % (escale_rms_fv*100.0,))
    print("FV RMS Discrepancy - Energy Resolution: %f%%" % eres_rms_fv)
