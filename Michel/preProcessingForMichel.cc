/*
 * Taken from McCaffrey:
 * /home/marzece/KDAR_Analysis/MichelAnalysis/HyoungkuMichel/preProcessingForMichel.cc
 * 
 * This takes pre-processed data and parses through it to create prompt and delayed michel event candidates.
 *
 * The output of this is used in:
 * update_michel_pair.cc
 * (/home/marzece/KDAR_Analysis/MichelAnalysis/HyoungkuMichel/update_michel_pair.cc)
 *
 * NEED TO UPDATE
 * TODO
 * - Use energy deposited
 */
#include <iostream>
#include <sstream>
#include <fstream>
#include <signal.h>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

using namespace std;
const double TTT_to_ns = 8.0; // From CAEN digitizer documentation
#define JADE_V0 0

sig_atomic_t sigint_recieved = 0;
void signal_handler(int sig) {
    const char msg[] = "SIGINT recieved, will be exiting main loop...\n";
    write(0, msg, sizeof(msg));
    if(sigint_recieved) {
        _exit(0);
    }
    sigint_recieved = 1;
}

int main(int argc, char ** argv){
    if(argc != 3) {
        cout << "Must provide 2 arguments, input filename (pre-processed JADE data; preprod_r00{run#}.root) & output filename" << endl;
        return -1;
    }
    const char* input_file_str = argv[1];
    const char* output_file_str = argv[2];

    struct sigaction sa;
    sigset_t sigset;
    sigemptyset(&sigset);
    sa.sa_handler = signal_handler;
    sa.sa_flags = 0;
    sa.sa_mask = sigset;
    sigaction(SIGINT, &sa, NULL);

    int TrigNum;
    vector<double>* totQ_TV = new vector<double>;
    vector<double>* totQ_BV = new vector<double>;
    vector<double>* QmaxQtot = new vector<double>;
    vector<double>* QTot = new vector<double>;
    vector<double>* recoFlux = new vector<double>;
    vector<double>* recoX = new vector<double>;
    vector<double>* recoY = new vector<double>;
    vector<double>* recoZ = new vector<double>;
    vector<int>* nsat_pmts = new vector<int>;
    vector<double>* inner_time = new vector<double>;

    int runNum = 0;
    int subRunNum = 0;
    int kicker_trigger = 0;
    double event_time = 0;

    unsigned int TTT;
    double time_since_kicker = -1;
    double last_kicker_time = -1;

    TChain * track = new TChain("InfoTree");
    track->Add(input_file_str);
    track->LoadTree(0);
    if(!(track->GetEntries() > 0)) {
        cout << "Input file(s) '"<< input_file_str << "' does not seem to have any data\n";
        return -1;
    }
    track->SetBranchAddress("t_runNum", &runNum);
    track->SetBranchAddress("t_subRunNum", &subRunNum);
    track->SetBranchAddress("t_trigID", &TrigNum);
    track->SetBranchAddress("TriggerBoardTimeTag", &TTT);
    track->SetBranchAddress("TV_TotalQ", &totQ_TV);
    track->SetBranchAddress("BV_TotalQ", &totQ_BV);
    track->SetBranchAddress("inner_QmaxQtot", &QmaxQtot);
    track->SetBranchAddress("inner_TotalQ", &QTot);
    track->SetBranchAddress("inner_RecoFlux", &recoFlux);
    track->SetBranchAddress("inner_Reco_Vertex_X", &recoX);
    track->SetBranchAddress("inner_Reco_Vertex_Y", &recoY);
    track->SetBranchAddress("inner_Reco_Vertex_Z", &recoZ);
    track->SetBranchAddress("SatPMTs", &nsat_pmts);
    track->SetBranchAddress("kicker_trigger", &kicker_trigger);
    track->SetBranchAddress("inner_Time", &inner_time);


    int tTrig = track->GetEntries();

    cout << " Total Origin Trigger entry : " << tTrig << endl;


    // Size each input vector according to how many entries are available.
    // Some entries might be empty (i.e. a CAEN readout occurred but no events
    // were made from the waveforms), so those will be skipped.
    vector<int> V_runNum(tTrig);
    vector<int> V_subRunNum(tTrig);
    vector<int> V_TrigID(tTrig);
    vector<double> V_TTT(tTrig);
    vector<int> V_KickerTrigger(tTrig);
    vector<double> V_EventT(tTrig);
    vector<double> V_RecoFlux(tTrig);
    vector<double> V_tot_TVQ(tTrig);
    vector<double> V_tot_BVQ(tTrig);
    vector<double> V_QmaxQtot(tTrig);
    vector<double> V_QTot(tTrig);
    vector<double> V_RecoVtxX(tTrig);
    vector<double> V_RecoVtxY(tTrig);
    vector<double> V_RecoVtxZ(tTrig);
    vector<double> V_NSatPMTS(tTrig);


    int trig_count = 0;
    // Go through all events and stuff them into arrays/vectors
    for(int iTrig = 0; iTrig < tTrig; iTrig++) {
        track->GetEntry( iTrig );

        if(recoFlux->empty()) {
            continue;
        }
        if((iTrig+1) % 500000 == 0 ) {
            cout << (iTrig+1) << " events was converted as a vector." << endl;
        }

        V_runNum[trig_count] =  runNum;
        V_subRunNum[trig_count] =  subRunNum;
        V_TrigID[trig_count] =  TrigNum;
        V_TTT[trig_count] = (double)TTT;
        V_RecoFlux[trig_count] = recoFlux->at(0);
        V_tot_TVQ[trig_count] = totQ_TV->at(0);
        V_tot_BVQ[trig_count] = totQ_BV->at(0);
        V_QmaxQtot[trig_count] = QmaxQtot->at(0);
        V_QTot[trig_count] = QTot->at(0);
        V_RecoVtxX[trig_count] = recoX->at(0);
        V_RecoVtxY[trig_count] = recoY->at(0);
        V_RecoVtxZ[trig_count] = recoZ->at(0);
        V_NSatPMTS[trig_count] = nsat_pmts->at(0);
        V_KickerTrigger[trig_count] = kicker_trigger;
        V_EventT[trig_count] = inner_time->at(0);
        trig_count += 1;
        if(sigint_recieved) {
            exit(1);
        }
    }
    // Resize all the data vectors to match the actual amount of data found
    V_runNum.resize(trig_count);
    V_subRunNum.resize(trig_count);
    V_TrigID.resize(trig_count);
    V_TTT.resize(trig_count);
    V_KickerTrigger.resize(trig_count);
    V_EventT.resize(trig_count);
    V_RecoFlux.resize(trig_count);
    V_tot_TVQ.resize(trig_count);
    V_tot_BVQ.resize(trig_count);
    V_QmaxQtot.resize(trig_count);
    V_QTot.resize(trig_count);
    V_RecoVtxX.resize(trig_count);
    V_RecoVtxY.resize(trig_count);
    V_RecoVtxZ.resize(trig_count);
    V_NSatPMTS.resize(trig_count);

    double FluxMeV = 0.;
    /*
    if( runNum >= 1400 ) { FluxMeV = 555.0; }
    if( runNum == 1500 ) { FluxMeV = 590.1; }
    if( runNum > 1500 ) { FluxMeV = 590.1; }
    if( runNum >= 1551 ) { FluxMeV = 5200/8.; }
    if( runNum >= 1600 ) { FluxMeV = 5200/8.; }
    */
#if JADE_V0
    FluxMeV = 659.96 // From JADE V0 Michel analysis..... (central value averaged over entire dataset)
#else
    FluxMeV = 942.183; // From Michel analysis by EDM (DocDB-2811)
#endif

    cout << endl;
    cout << " Current run number : " << runNum << endl;
    cout << " Flux to MeV conv F : " << FluxMeV << " flux/MeV " << endl;

    int count = 0;
    int prompt_count = 0;
    int delayed_count = 0;

    double TotalLiveTime = 0.;
    double now_TTT = 0.;
    double prev_TTT = 0.;
    double temp_TTT = 0.;

    int s1_index = 0;
    int s1_index_prev = 0;
    int s1_runNum_prev = 0;
    int s1_subRunNum_prev = 0;
    int s2_index = 0;
    int s1_runNum = 0;
    int s1_subRunNum = 0;
    int s2_runNum = 0;
    int s2_subRunNum = 0;
    double s1_Flux = 0.;
    double s1_Energy = 0.;
    double s1_FluxPrev = 0.;
    double s1_QTot = 0;
    double s1_EnergyPrev = 0.;
    double s2_Flux = 0.;
    double s2_Energy = 0.;
    double s2_QTot = 0;
    double s1_T = 0.;
    double s1_Tprev = 0.;
    double s2_T = 0.;
    double s1_TimeSinceKicker = 0.;
    double s2_TimeSinceKicker = 0.;

    double s1_TV = 0.;
    double s1_BV = 0.;
    double s1_TVprev = 0.;
    double s1_BVprev = 0.;
    double s2_TV = 0.;
    double s2_BV = 0.;

    double s1_x = 0.;
    double s1_y = 0.;
    double s1_z = 0.;
    double s1_nsat = 0.;
    double s1prev_x = 0.;
    double s1prev_y = 0.;
    double s1prev_z = 0.;
    int s1prev_nsat = 0.;
    double s2_x = 0.;
    double s2_y = 0.;
    double s2_z = 0.;
    int s2_nsat = 0;

    cout << " Total # of Triggers : " << tTrig << endl;

    // Event type flags
    int Flag_mu = false;
    bool Flag_prompt = false;
    bool Flag_delayed = false;


    // Cut value
    const double MUON_TOP_VETO_CHARGE_MINIMUM = 100; // PE
    const double MUON_BOTTOM_VETO_CHARGE_MINIMUM = 100; // PE

    const double total_t_width = 100; // 100 us
    const double PROMPT_MINIMUM_ENERGY = 0; // MeV
    const double PROMPT_MAXIMUM_ENERGY = 800; // MeV
    const double DELAYED_MINIMUM_ENERGY = 0; // MeV
    const double DELAYED_MAXIMUM_ENERGY = 100; // MeV

    int TrigID;
    int pairing_s1Index;
    int pairing_s2Index;
    int pairing_s1RunNum;
    int pairing_s2RunNum;
    int pairing_s1SubRunNum;
    int pairing_s2SubRunNum;
    double pairing_s1Energy;
    double pairing_s1Flux;
    double pairing_s1QTot;
    double pairing_s1T;
    double pairing_s2Energy;
    double pairing_s2Flux;
    double pairing_s2T;
    double pairing_s2QTot;

    double pairing_s1TV;
    double pairing_s1BV;
    double pairing_s2TV;
    double pairing_s2BV;

    double pairing_s1x;
    double pairing_s1y;
    double pairing_s1z;
    double pairing_s2x;
    double pairing_s2y;
    double pairing_s2z;
    double pairing_s1_nsat;
    double pairing_s2_nsat;
    double pairing_s1_timesincekicker;
    double pairing_s2_timesincekicker;


    TFile * outputFile = new TFile(output_file_str, "RECREATE");
    TTree * outputTree = new TTree("InfoTree","InfoTree");
    outputTree->Branch("TrigID", &TrigID, "TrigID/I");
    outputTree->Branch("pairing_s1Index", &pairing_s1Index, "pairing_s1Index/I");
    outputTree->Branch("pairing_s1RunNum", &pairing_s1RunNum, "pairing_s1RunNum/I");
    outputTree->Branch("pairing_s1SubRunNum", &pairing_s1SubRunNum, "pairing_s1SubRunNum/I");
    outputTree->Branch("pairing_s1Energy", &pairing_s1Energy, "pairing_s1Energy/D");
    outputTree->Branch("pairing_s1Flux", &pairing_s1Flux, "pairing_s1Flux/D");
    outputTree->Branch("pairing_s1QTot", &pairing_s1QTot, "pairing_s1QTot/D");
    outputTree->Branch("pairing_s1T", &pairing_s1T, "pairing_s1T/D");
    outputTree->Branch("pairing_s2Index", &pairing_s2Index, "pairing_s2Index/I");
    outputTree->Branch("pairing_s2RunNum", &pairing_s2RunNum, "pairing_s2RunNum/I");
    outputTree->Branch("pairing_s2SubRunNum", &pairing_s2SubRunNum, "pairing_s2SubRunNum/I");
    outputTree->Branch("pairing_s2Energy", &pairing_s2Energy, "pairing_s2Energy/D");
    outputTree->Branch("pairing_s2QTot", &pairing_s2QTot, "pairing_s2QTot/D");
    outputTree->Branch("pairing_s2Flux", &pairing_s2Flux, "pairing_s2Flux/D");
    outputTree->Branch("pairing_s2T", &pairing_s2T, "pairing_s2T/D");
    outputTree->Branch("pairing_s1TV", &pairing_s1TV, "pairing_s1TV/D");
    outputTree->Branch("pairing_s1BV", &pairing_s1BV, "pairing_s1BV/D");
    outputTree->Branch("pairing_s2TV", &pairing_s2TV, "pairing_s2TV/D");
    outputTree->Branch("pairing_s2BV", &pairing_s2BV, "pairing_s2BV/D");

    outputTree->Branch("pairing_s1x", &pairing_s1x, "pairing_s1x/D");
    outputTree->Branch("pairing_s1y", &pairing_s1y, "pairing_s1y/D");
    outputTree->Branch("pairing_s1z", &pairing_s1z, "pairing_s1z/D");
    outputTree->Branch("pairing_s2x", &pairing_s2x, "pairing_s2x/D");
    outputTree->Branch("pairing_s2y", &pairing_s2y, "pairing_s2y/D");
    outputTree->Branch("pairing_s2z", &pairing_s2z, "pairing_s2z/D");

    outputTree->Branch("pairing_s1_nsat", &pairing_s1_nsat, "pairing_s1_nsat/I");
    outputTree->Branch("pairing_s2_nsat", &pairing_s2_nsat, "pairing_s2_nsat/I");
    outputTree->Branch("pairing_s1_timesincekicker", &pairing_s1_timesincekicker, "pairing_s1_timesincekicker/D");
    outputTree->Branch("pairing_s2_timesincekicker", &pairing_s2_timesincekicker, "pairing_s2_timesincekicker/D");


    outputTree->Branch("EffLiveTime", &TotalLiveTime, "TotalLiveTime/D");


    TH1D * evtRate = new TH1D("totalEvt", "", 500, 0, 100);
    TH1D * muonRate_s = new TH1D("muonEvt", "", 500, 0, 100);
    TH1D * muonRate = new TH1D("muonRate", "", 500, 0, 1000);

    int cosmu_count = 0;
    int cosmu_countTV = 0;
    int cosmu_countBV = 0;
    int cosmu_countTVBV = 0;

    for(size_t iTrig = 0; iTrig < V_TrigID.size(); iTrig++) {
        if ( (iTrig+1) % 20000 == 0 ) {
            cout << (iTrig+1) << " done. " << endl;
        }

        now_TTT = V_TTT[iTrig];
        kicker_trigger = V_KickerTrigger[iTrig];
        event_time = V_EventT[iTrig];

        if( now_TTT < prev_TTT) { // Indicates a TTT roll-over
            // temp_TTT is the "global" timestamp (since the start of the run)
            // but it only updates when roll-overs occurs
            temp_TTT += prev_TTT;
        }
        TotalLiveTime = (now_TTT + temp_TTT)*TTT_to_ns;
        prev_TTT = now_TTT;

        if(V_RecoFlux[iTrig] <= 0) {
            continue;
        }
        if(V_QmaxQtot[iTrig] > 1) {
            continue;
        }

        if(kicker_trigger) {
            time_since_kicker = 0;
            last_kicker_time = TotalLiveTime + event_time*2;
        }
        else if(last_kicker_time > 0) {
            time_since_kicker = TotalLiveTime + event_time*2 - last_kicker_time;
        }

        double RecoEnergy = V_RecoFlux[iTrig]/FluxMeV;
        evtRate->Fill(RecoEnergy);

        // Mark this event as a muon if it's top/bottom veto charge is above some threshold.
        Flag_mu = V_tot_TVQ[iTrig] >= MUON_TOP_VETO_CHARGE_MINIMUM ||
                  V_tot_BVQ[iTrig] >= MUON_BOTTOM_VETO_CHARGE_MINIMUM;

        if(Flag_mu) {
            cosmu_count++;
            muonRate_s->Fill(RecoEnergy);
            muonRate->Fill(RecoEnergy);
            if(V_tot_TVQ[iTrig] >= MUON_TOP_VETO_CHARGE_MINIMUM) {
                cosmu_countTV++;
            }
            if(V_tot_BVQ[iTrig] >= MUON_BOTTOM_VETO_CHARGE_MINIMUM) {
                cosmu_countBV++;
            }
            if(V_tot_TVQ[iTrig] >= MUON_TOP_VETO_CHARGE_MINIMUM && V_tot_BVQ[iTrig] >= MUON_BOTTOM_VETO_CHARGE_MINIMUM) {
                cosmu_countTVBV++;
            }
        }

        count++;
        if(Flag_mu && RecoEnergy > PROMPT_MINIMUM_ENERGY && RecoEnergy < PROMPT_MAXIMUM_ENERGY) {
            prompt_count++;
            s1_index = V_TrigID[iTrig];
            s1_runNum = V_runNum[iTrig];
            s1_subRunNum = V_subRunNum[iTrig];

            s1_Flux = V_RecoFlux[iTrig];
            s1_Energy = RecoEnergy;
            s1_T = TotalLiveTime + event_time*2;
            s1_TimeSinceKicker = time_since_kicker;
            s1_QTot = V_QTot[iTrig];

            s1_TV = V_tot_TVQ[iTrig];
            s1_BV = V_tot_BVQ[iTrig];

            s1_x = V_RecoVtxX[iTrig];
            s1_y = V_RecoVtxY[iTrig];
            s1_z = V_RecoVtxZ[iTrig];
            s1_nsat = V_NSatPMTS[iTrig];

            Flag_prompt = true;
        }
        if(!Flag_mu && RecoEnergy > DELAYED_MINIMUM_ENERGY && RecoEnergy < DELAYED_MAXIMUM_ENERGY) { delayed_count++;
            s2_index = V_TrigID[iTrig];
            s2_runNum = V_runNum[iTrig];
            s2_subRunNum = V_subRunNum[iTrig];
            s2_Flux = V_RecoFlux[iTrig];
            s2_Energy = RecoEnergy;
            s2_QTot = V_QTot[iTrig];
            s2_T = TotalLiveTime + event_time*2;
            s2_TimeSinceKicker = time_since_kicker;

            s2_TV = V_tot_TVQ[iTrig];
            s2_BV = V_tot_BVQ[iTrig];

            s2_x = V_RecoVtxX[iTrig];
            s2_y = V_RecoVtxY[iTrig];
            s2_z = V_RecoVtxZ[iTrig];
            s2_nsat = V_NSatPMTS[iTrig];

            Flag_delayed = true;
        }
        double dT = (s2_T - s1_Tprev)/1000.; // Find delta-t and convert from nano to micro-seconds
        if( dT > 0 && dT < total_t_width && Flag_prompt && Flag_delayed){
            // A valid pair of events that passes all cuts!
            TrigID = iTrig;
            pairing_s1Index = s1_index_prev;
            pairing_s1RunNum = s1_runNum_prev;
            pairing_s1SubRunNum = s1_subRunNum_prev;
            pairing_s1Energy = s1_EnergyPrev;
            pairing_s1Flux = s1_FluxPrev;
            pairing_s1QTot = s1_QTot;
            pairing_s1T = s1_Tprev;
            pairing_s2Index = s2_index;
            pairing_s2RunNum = s2_runNum;
            pairing_s2SubRunNum = s2_subRunNum;
            pairing_s2Flux = s2_Flux;
            pairing_s2Energy = s2_Energy;
            pairing_s2QTot = s2_QTot;
            pairing_s2T = s2_T;

            pairing_s1TV = s1_TVprev;
            pairing_s1BV = s1_BVprev;
            pairing_s2TV = s2_TV;
            pairing_s2BV = s2_BV;

            pairing_s1x = s1prev_x;
            pairing_s1y = s1prev_y;
            pairing_s1z = s1prev_z;
            pairing_s2x = s2_x;
            pairing_s2y = s2_y;
            pairing_s2z = s2_z;

            pairing_s1_nsat = s1prev_nsat;
            pairing_s2_nsat = s2_nsat;

            pairing_s1_timesincekicker = s1_TimeSinceKicker;
            pairing_s2_timesincekicker = s2_TimeSinceKicker;

            Flag_prompt  = false;
            Flag_delayed  = false;
            outputTree->Fill();
        }

        s1_index_prev = s1_index;
        s1_runNum_prev = s1_runNum;
        s1_subRunNum_prev = s1_subRunNum;
        s1_FluxPrev = s1_Flux;
        s1_EnergyPrev = s1_Energy;
        s1_Tprev = s1_T;
        s1prev_x = s1_x;
        s1prev_y = s1_y;
        s1prev_z = s1_z;
        s1prev_nsat = s1_nsat;

        s1_TVprev = s1_TV;
        s1_BVprev = s1_BV;

        if(sigint_recieved) {
            break;
        }
    }

    TH1D * cosMuRateSave = new TH1D("cosMuCount","", 1, 0, 1);
    cosMuRateSave->SetBinContent(0, cosmu_count );

    cout << " Total Events : " << count << endl;
    cout << " Total Live Time : " << TotalLiveTime*1e-9 << " s " << endl;
    cout << " Single rate : " << count/(TotalLiveTime*1e-9) << endl << endl;

    double effLiveT = TotalLiveTime*1e-9;
    double cosmu_rate = cosmu_count/effLiveT;
    double cosmu_rate_TVOnly = cosmu_countTV/effLiveT;
    double cosmu_rate_BVOnly = cosmu_countBV/effLiveT;
    double cosmu_rate_Both = cosmu_countTVBV/effLiveT;
    cout << " Cosmic muon rate : " << cosmu_rate << " Hz +- " << sqrt(cosmu_rate/effLiveT) << endl;
    cout << " Cosmic muon rate [TV] : " << cosmu_rate_TVOnly << " Hz " << endl;
    cout << " Cosmic muon rate [BV] : " << cosmu_rate_BVOnly << " Hz " << endl;
    cout << " Cosmic muon rate [Both] : " << cosmu_rate_Both << " Hz " << endl;

    cout << endl << " Prompt Candidate : " << prompt_count << endl;
    cout << " Prompt Candidate-rate: " << prompt_count/(TotalLiveTime*1e-9) << endl;
    cout << " Delayed Candidate : " << delayed_count << endl;
    cout << " Delayed Candidate-rate : " << delayed_count/(TotalLiveTime*1e-9) << endl;

    TH1D * oriTot = (TH1D*)evtRate->Clone();
    oriTot->SetName("OriTotal");
    TH1D * oriMu = (TH1D*)muonRate_s->Clone();
    oriMu->SetName("OriMu");

    evtRate->Scale(1./effLiveT/1000.);
    muonRate_s->Scale(1./effLiveT/1000.);

    cosMuRateSave->Write();
    evtRate->Write();
    muonRate_s->Write();
    muonRate->Write();
    oriTot->Write();
    oriMu->Write();
    outputTree->Write();
    outputFile->Write();
    outputFile->Close();
}
