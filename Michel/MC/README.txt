To generate the MC and reconfigure it for use in the spatial analysis,
there are a few steps one needs to do.

If you need to generate the MC dataset from RAT:

michel_rat.sh
Calls RAT to use a macro in "mac/" to generate inital MC files to be processed in jade.
Currently is set-up to create a mu+ and mu- set of data


michel_rat_proc.sh
Processes MC data with JADE.


michel_rat_reco.sh
Completes JADE reconstruction with the processed MC data.

----
michel_build_trees.sh
Runs "build_trees.py" for the defined file indexes. This reformats all of the MC files generated 
for a given run so far into a root file that contains event information for michel pairs.

After completing this for every MC file, use 'hadd' to combine the outputs for each kind of mu, so that mu+ and mu- each have their own file that contains all the data.

mu_weight_tree.cc
Takes the combined MC root files from michel_build_trees.sh (with pair_tree)
and creates weight trees to be used later to more accurately reflect mu+/mu- ratio.

So ultimately, the output will be same same as each original file, but with the inclusion
of a weight tree.**EXPLAIN WHEN THIS HAPPENS**  When the mu+ & mu- files are merged with hadd, it will also merge the 
weight tree which can be used during analysis by declaring it as a friend of pair_tree


extract_mc_michel_info.py
Scrapes the data in the combined pair_tree file and its accompanying "friend" file and formats 
the event info into a michel_tree that will be used in the spatial analysis.


 
-----------------
This folder also contains the MC info for the energy deposited fit.

