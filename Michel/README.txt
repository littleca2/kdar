
This describes the workflow for getting michel electron information from JSNS2 data.
Using the michel electron's energy spectrum, we calibrate against a monte carlo simulation to find the flux to Mev ratio.

Begin with pre-processed run data (e.g. /mnt/braxton/KDAR_JADE_v3/preprod/)

1. preProcessingForMichel
./preProcessingForMichel input_filename(pre-processed JADE data; preprod_r00{run#}.root) output filename
 * This takes pre-processed data and parses through it to create prompt and delayed michel event candidates.
 *
 * The output of this is used in:
 * update_michel_pair.cc
 * (/home/marzece/KDAR_Analysis/MichelAnalysis/HyoungkuMichel/update_michel_pair.cc)

=======================================================
2. updateMichelPair
./updateMichelPair input_filename(preProcessingForMichel output; michel_pre_{run#}.root) output filename

 * This takes the output of preProcessingForMichel.cc and finds michel events.
 * It outputs these events as well as graphs their energy spectrums.

 * NOTE: The user needs to combine the output files using something like ROOT's $(hadd)
 * to create a file named "michel_pair_combined.root"

=======================================================
(./flux2mev_corrections)
3. get_data_params.py
python get_data_params.py input_filename(updateMichelPair output; michle_pair_{run#}.root) runType(0: Sum, 1: Run, 2:Subrun) version(Name to be used to discern this set of data)

 * This takes the outputs of updateMichelPair.cc and creates a best fit of the flux data using the
 * muon MC simulation data found in ~/kdar/Michel/MC/MC_mu_[minus/plus]_FV_edep_vals.nptxt.
 * It updates the ~/kdar/correction_values.json file which contains all of the FluxtoMev (+other) information for different variations of analysis (versions)
 * It also saves individual fit information in a directory in /output_fluxCorr named after the user input version.

 * The runType variable decides if it will calculate
 *	0: The michel_pair_combined.root data to find a starting FluxtoMeV value for subsequent runs or subruns which is saved in ~/kdar/correction_values.json
 *	1: For each run in updateMichelPair's output, it will find the FluxtoMeV values (+other) to append to ~/kdar/correction_values.json
 *	2: For each subrun in updateMichelPair's output, it will find the FluxtoMeV values (+other) to append to ~/kdar/correction_values.json. (NOTE: This probably isn't totally implemented, yet)

 * It uses the definitions in fitting.py which find the flux to MeV conversion by:
 *	Reading in the flux data
 *	Taking the muon MC energy deposited data and attempting to manipulate it to match the flux data using the variables:
 *		- Y-axis scaling
 *		- Energy scale (Scaling factor for flux bin centers -> MeV bin centers)
 *		- Energy resolution scaling ("a" in the energy resolution model in docDB 3009-v5 eq. 10)
 *	Then smearing the MC data using a gaussian function with the calculated energy resolution as the sigma.
 *	Multiplying the MC datas X-axis [MeV] by the energy scale gives the converted MeV->flux bin centers for the MC data.
 *	Using interpolation, we try and re-create the flux data with the MC X and Y values we've found.
 *	optimize.curve_fit uses this function to find our best fit variables.

 * The energy resoluton scaling is used as an analog for the actual energy resolution.
 * We ignore the constant term since it usually ~ 0.

=======================================================
(./flux2mev_corrections)
4. plot_time_flux.py
python plot_time_flux.py version --michelEndpoint Eep (For version) --compare compareID (version name for data you'd like to compare "version" against) --compareMichelEndpoint Eepc (Michel endpoint for comparison version)

 * This takes the values in ~/kdar/correction_values.json and plots the Flux to MeV conversions (among other things) over time. It discerns individual run periods by using a user defined time limit between runs.
 * It compares to ~2021 energy correction values from HyoungKu, Dodo, and Eric. 
 * NOTE: Dodo's normalization/errors may not be totally accurate with respect to how it calculates the rest of the data.
 * The output is saved in the same place as get_data_params.py (/output_fluxCorr/${version})


