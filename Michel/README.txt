
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


2. updateMichelPair
./updateMichelPair input_filename(preProcessingForMichel output; michel_pre_{run#}.root) output filename

 * This takes the output of preProcessingForMichel.cc and finds michel events.
 * It outputs these events as well as graphs their energy spectrums.


(./flux2mev_corrections)
3. get_data_params.py



4. plot_time_flux.py



