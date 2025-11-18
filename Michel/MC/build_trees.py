# June 20 20223, Eric M, Updated to add saturation info from the process file
# to the output data file

import numpy as np
from array import array
from rat import RAT
from rat import ROOT

### JADE ###
ROOT.gSystem.Load('libJADE.so')
ROOT.gInterpreter.ProcessLine('using namespace JADE;')
print("Loaded the Libraries")


def get_pdg_code(mcparticle):
    # Use abs() here b/c I want to use this code both for mu+ and mu- events
    # TODO, should think of some more elegant solution...
    return abs(mcparticle.GetPDGCode())

# This is meant to be streamlined code for producing tree objects for higher level analysis
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("output", help="Output ROOT filename")
    parser.add_argument("rat", help="Input RAT MC file")
    parser.add_argument("proc", help="Input JADE process file")
    parser.add_argument("reco", help="Input JADE reconstruction file")
    parser.add_argument("index", type=int, help="file index, for debugging purposes only")
    args = parser.parse_args()
    outName = args.output
    rat_name = args.rat
    proc_name = args.proc
    reco_name = args.reco
    file_index = args.index

    # Put the rat indices of bad events in this tree, hopefully it rarely happens
    event_status_tree = ROOT.TNtuple("bad_events", "bad_events", "F_idx:Rat_idx")

    primary_truth_tree = ROOT.TNtuple("primary_truth_tree", "primary_truth_tree", "F_idx:Rat_idx:PDG:KE:t:x:y:z")
    secondary_truth_tree = ROOT.TNtuple("secondary_truth_tree", "secondary_truth_tree", "F_idx:Rat_idx:PDG:KE:t:x:y:z")

    # This is the important tree for the KDAR Analysis
    outFile = ROOT.TFile(outName, "RECREATE")
    pair_tree = ROOT.TTree("pair_tree", "pair_tree")

    # Rat info/branches
    F_idx = array('i', [0])
    R_idx = array('i', [0])
    E_numu = array('f', [0])
    E_vis = array('f', [0])
    E_dep = array('f', [0])
    x_p = array('f', [0])
    y_p = array('f', [0])
    z_p = array('f', [0])
    t_p = array('f', [0])

    E_d = array('f', [0])
    E_d_dep = array('f', [0])
    x_d = array('f', [0])
    y_d = array('f', [0])
    z_d = array('f', [0])
    t_d = array('f', [0])

    pair_tree.Branch("F_idx_b", F_idx, "F_idx_b/I")
    pair_tree.Branch("R_idx_b", R_idx, "R_idx_b/I")
    pair_tree.Branch("E_numu_b", E_numu, "E_numu_b/F")
    pair_tree.Branch("E_dep_b", E_dep, "E_dep_b/F")
    # Prompt Info
    pair_tree.Branch("E_vis_b", E_vis, "E_vis_b/F")
    pair_tree.Branch("x_p_b", x_p, "x_p_b/F")
    pair_tree.Branch("y_p_b", y_p, "y_p_b/F")
    pair_tree.Branch("z_p_b", z_p, "z_p_b/F")
    pair_tree.Branch("t_p_b", t_p, "t_p_b/F")
    # Delayed Info
    pair_tree.Branch("E_d_b", E_d, "E_d_b/F")
    pair_tree.Branch("E_d_dep_b", E_d_dep, "E_d_dep_b/F")
    pair_tree.Branch("x_d_b", x_d, "x_d_b/F")
    pair_tree.Branch("y_d_b", y_d, "y_d_b/F")
    pair_tree.Branch("z_d_b", z_d, "z_d_b/F")
    pair_tree.Branch("t_d_b", t_d, "t_d_b/F")


    # Jade Info/Branches
    max_L = 20 # Maximum Number of SubEvents per RAT event --> Would be absolutely shocked if this could be exceeded

    Trig_arr = array('i', max_L*[-1])
    SubE_arr = array('i', max_L*[-1])
    flux_arr = array('f', max_L*[-1.])
    TTT_arr = array('f', max_L*[-1.])
    t_reco_arr = array('f', max_L*[-1.])
    x_reco_arr = array('f', max_L*[-1.])
    y_reco_arr = array('f', max_L*[-1.])
    z_reco_arr = array('f', max_L*[-1.])
#    reco_min_arr = array('f', max_L*[-1.])
    saturation_arr = array('i', max_L*96*[-1])
    temp_str = "PMTSaturation[%i]/I" % (96*max_L,)
    pair_tree.Branch("Trig_arr_b", Trig_arr, "Trig_arr_b["+str(max_L)+"]/I")
    pair_tree.Branch("SubE_arr_b", SubE_arr, "SubE_arr_b["+str(max_L)+"]/I")
    pair_tree.Branch("flux_arr_b", flux_arr, "flux_arr_b["+str(max_L)+"]/F")
    pair_tree.Branch("TTT_arr_b", TTT_arr, "TTT_arr_b["+str(max_L)+"]/F")
    pair_tree.Branch("PMTSaturation", saturation_arr, temp_str)
    pair_tree.Branch("t_reco_arr_b", t_reco_arr, "t_reco_arr_b["+str(max_L)+"]/F")
    pair_tree.Branch("x_reco_arr_b", x_reco_arr, "x_reco_arr_b["+str(max_L)+"]/F")
    pair_tree.Branch("y_reco_arr_b", y_reco_arr, "y_reco_arr_b["+str(max_L)+"]/F")
    pair_tree.Branch("z_reco_arr_b", z_reco_arr, "z_reco_arr_b["+str(max_L)+"]/F")
#    pair_tree.Branch("reco_min_arr_b", reco_min_arr, "reco_min_arr_b["+str(max_L)+"]/F")


    ###################### Load the Rat Stuff ##################
    runTree = ROOT.TChain("runT")
    runTree.Add(rat_name)
    RAT.DS.RunStore.SetReadTree(runTree)
    run = RAT.DS.Run(RAT.DS.RunStore.Get().GetRun(1))
    if (run == 0):
        print("Run not Found!!!!!!")
    ds = RAT.DSReader(rat_name)

    # Particle PDG Codes
    numu_code = 14
    proton_code = 2212
    neutron_code = 2112
    muon_code = 13 # negative
    electron_code = 11

    # Read in the Files to be analyzed
    procFile = ROOT.TFile.Open(proc_name, "READ")
    recoFile = ROOT.TFile.Open(reco_name, "READ")
    triggerTree = procFile.Get("triggerTree")
    recoTriggerTree = recoFile.Get("recoTriggerTree")
    numRatEvents = ds.GetTotal()
    numTriggers = triggerTree.GetEntries()
    print("\nNumber of Rat Events: ", numRatEvents)
    print("Number of Proc Triggers:", numTriggers, "Reco Triggers:", recoTriggerTree.GetEntries())

    iTrig = 0
    for RatEvent in range(0, numRatEvents):
        pair = False
        root = ds.GetEvent(RatEvent)
        michel_energy = 0
        # Monte Carlo Information for this event
        mc = root.GetMC()
        summary = mc.GetMCSummary()
        veto_edep = summary.GetEnergyLossByVolume("GCVeto")
        Edep = summary.GetTotalScintEdep() - veto_edep
        # Number of Triggered Rat Events
        EVCount = root.GetEVCount()

        numParticles = mc.GetMCParticleCount()
        numParents = mc.GetMCParentCount()
        numSecondaries = mc.GetMCSecondaryCount()

        trig_list = [] # Let's grab all the trigger indices for this RAT Event
        for num in range(EVCount):
            trig_list.append(iTrig)
            iTrig += 1
        if(numSecondaries == 0):
            continue

        count = 0
        F_idx[0] = file_index
        R_idx[0] = RatEvent

        for index in trig_list:
            triggerTree.GetEntry(index)
            recoTriggerTree.GetEntry(index)

            numSubEvents = triggerTree.trigger.GetEventCountID()
            numRecoEvents = recoTriggerTree.recoTrigger.GetRecoEventCount()
            print("In Reco Loop:")
            print("NumSubEvents:", numSubEvents)
            if numSubEvents != numRecoEvents:
                print("PROBLEM->Sub Events and Reco Events do not agree !!!!!!!!!")

            for i in range(0, numSubEvents):
                event = triggerTree.trigger.GetEventID(i)
                recoEvent = recoTriggerTree.recoTrigger.GetRecoEvent(i)
                eventTime = (event.GetTime()*2.0)/1000.0 # microsec
                TTT = triggerTree.trigger.GetTriggerBoardTimeTag() # clock ticks?
                TTT = TTT*8e-3 # microsec
                totalTimeOffset = (TTT + eventTime) # microsec
                x_reco = recoEvent.GetRecoVertex().X()*1000 # mm
                y_reco = recoEvent.GetRecoVertex().Y()*1000 # mm
                z_reco = recoEvent.GetRecoVertex().Z()*1000 # mm
                flux = recoEvent.GetRecoFlux()
                for ipmt in range(96):
                    saturation_arr[96*count +ipmt] = 1 if event.GetPMTSaturatedFlag(ipmt) else 0
                Trig_arr[count] = index
                SubE_arr[count] = i
                flux_arr[count] = flux
                TTT_arr[count] = TTT
                t_reco_arr[count] = totalTimeOffset
                x_reco_arr[count] = x_reco
                y_reco_arr[count] = y_reco
                z_reco_arr[count] = z_reco
#                reco_min_arr[count] = sum(recoEvent.GetRecoMinPMT(ipmt) for ipmt in range(96))
                count += 1

        if numRatEvents % 100 == 0:
            print("RatEvent:", RatEvent, "Triggered ->", EVCount)

        try:
            code = get_pdg_code(mc.GetMCParticle(0))
        except:
            print("Bad Event !!??")
            event_status_tree.Fill(RatEvent)
            continue

        # Find the MC Truth Muon
        prompt_muon = 0
        for i in range(0, numParticles):
            code = get_pdg_code(mc.GetMCParticle(i))
            if abs(code) == muon_code:
                prompt_muon = mc.GetMCParticle(i)
                break

        print("Made it: %i" % RatEvent)
        muon_time = prompt_muon.GetTime()/1000.0 # us
        prompt_x = prompt_muon.GetPosition().X() # mm
        prompt_y = prompt_muon.GetPosition().Y() # mm
        prompt_z = prompt_muon.GetPosition().Z() # mm
        R_prompt = np.sqrt(prompt_x**2 + prompt_y**2)
        prompt_energy = 0

        michel_energy = 0
        michel_time = 0
        delayed_x = 0
        delayed_y = 0
        delayed_z = 0

        for i in range(numParticles):
            part = mc.GetMCParticle(i)
            code = get_pdg_code(part)
            primary_truth_tree.Fill(file_index, RatEvent, code, part.GetKE(),
                                    part.GetTime()/1000.0, part.GetPosition().X(),
                                    part.GetPosition().Y(), part.GetPosition().Z())
            if code in [muon_code, proton_code]:
                prompt_energy += part.GetKE()

        E_vis[0] = prompt_energy
        E_dep[0] = Edep
        x_p[0] = prompt_x
        y_p[0] = prompt_y
        z_p[0] = prompt_z
        t_p[0] = muon_time

        numu_energy = 0
        for i in range(numParents):
            if get_pdg_code(mc.GetMCParent(i)) == numu_code:
                numu_energy = mc.GetMCParent(i).GetKE() # MeV
        E_numu[0] = numu_energy

        for p in range(numSecondaries):
            part = mc.GetMCSecondary(p)
            code = get_pdg_code(part)
            secondary_truth_tree.Fill(file_index, RatEvent, code, part.GetKE(),
                                      part.GetTime()/1000.0, part.GetPosition().X(),
                                      part.GetPosition().Y(), part.GetPosition().Z())

            if get_pdg_code(mc.GetMCSecondary(p)) == electron_code:
                # Found a delayed event that may have triggered !
                michel = mc.GetMCSecondary(p)
                michel_energy = mc.GetMCSecondary(p).GetKE()
                michel_time = mc.GetMCSecondary(p).GetTime()/1000.0
                del_t = michel_time - muon_time
                delayed_x = michel.GetPosition().X() # mm
                delayed_y = michel.GetPosition().Y() # mm
                delayed_z = michel.GetPosition().Z() # mm

                mcs = mc.GetMCSummary()
                michel_dep_energy = mcs.GetTotalScintEdep()
                #if (abs(prompt_z) <= 1250): # Don't allow Chimney prompt events. This is a quality control criteria
                pair = True
                E_d[0] = michel_energy
                E_d_dep[0] = michel_dep_energy
                x_d[0] = delayed_x
                y_d[0] = delayed_y
                z_d[0] = delayed_z
                t_d[0] = michel_time

        if not pair:
            # No Michel elextron, may have triggered though?
            # -> Muon Nuclear Capture Event ~ 3% of the time
            E_d[0] = -1
        pair_tree.Fill()

        # Reset Storage Variables for the next event
        for num in range(max_L):
            Trig_arr[num] = -1
            SubE_arr[num] = -1
            flux_arr[num] = -1.
            TTT_arr[num] = -1.
            t_reco_arr[num] = -1.
            x_reco_arr[num] = -1.
            y_reco_arr[num] = -1.
            z_reco_arr[num] = -1.
            for ipmt in range(96):
                saturation_arr[num*96 + ipmt] = -1

    # Write Trees to an ouput file
    outFile.cd()
    event_status_tree.Write()
    pair_tree.Write()
    primary_truth_tree.Write()
    secondary_truth_tree.Write()
    outFile.Close()
    procFile.Close()
    recoFile.Close()
