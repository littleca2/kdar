

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




