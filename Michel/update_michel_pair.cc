/*
 * Taken from:
 * /home/marzece/KDAR_Analysis/MichelAnalysis/HyoungkuMichel/update_michel_pair.cc
 * 
 * This takes the output of preProcessingForMichel.cc and finds michel events.
 * It outputs these events as well as graphs their energy spectrums.
 * 
 * NEED TO UPDATE
 * TODO
 * - use energy deposited
 */

#include <iostream>
#include <sstream>

#include <TStyle.h>
#include <TChain.h>
#include <TVector3.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TNtuple.h>

using std::cout;
using std::endl;
using std::string;

double fcn_2exp(double *x, double *p){
    return exp(p[0]-x[0]/p[1]) + exp(p[2]-x[0]/p[3]);
}

void style_histogram(TH1D* h, int line_color=-1, int marker_style=-1, double ytitle_offset=-1, double ylabel_offset=-1, double xtitle_offset=-1, double xlabel_offset=-1) {
    if(line_color >= 0) {
        h->SetLineColor(line_color);
        h->SetMarkerColor(line_color);
    }
    if(marker_style >= 0){
        h->SetMarkerStyle(marker_style);
    }
    if(ytitle_offset >= 0) {
        h->GetYaxis()->SetTitleOffset(ytitle_offset);
    }
    if(ylabel_offset >= 0) {
        h->GetYaxis()->SetLabelOffset(ylabel_offset);
    }
    if(xtitle_offset >= 0) {
        h->GetXaxis()->SetTitleOffset(xtitle_offset);
    }
    if(xlabel_offset >= 0) {
        h->GetXaxis()->SetLabelOffset(xlabel_offset);
    }
}

const double MICHEL_DELTA_T_MAX = 10; // 10us
const double MICHEL_DELTA_T_MIN = 0;
const double ACCIDENTAL_DELTA_T_MIN = MICHEL_DELTA_T_MAX; // Once the michel window ends the "accidental" window starts
const double ACCIDENTAL_DELTA_T_MAX = 20; // 20us
const double DELTA_VERTEX_THRESHOLD = 1.300; // meters

bool CheckSinglePos(const TVector3 start){
    const double TankHeight = 1250.;
    const double TankRadius = 1580.;
    const double sX = start.X();
    const double sY = start.Y();
    const double sZ = start.Z();
    const double sR = sqrt(sX*sX + sY*sY);
    return  sR < TankRadius && abs(sZ) < TankHeight;
}

int main(int argc, char** argv) {
    if(argc != 3) {
        cout << "Must provide input filename (preProcessingForMichel output; michel_pre_{run#}.root) & output filename as arguments" << endl;
        return -1;
    }
    string input_filename = argv[1];
    string output_filename = argv[2];
    int TrigID;
    double pairing_s1Energy;
    double pairing_s1Flux;
    double pairing_s1T;
    double pairing_s2Energy;
    double pairing_s2Flux;
    double pairing_s2T;
    double pairing_s1_timesincekicker;
    int pairing_s1Index, pairing_s2Index;
    int pairing_s1RunNum, pairing_s2RunNum;
    int pairing_s1SubRunNum ,pairing_s2SubRunNum;

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
    int pairing_s1_nsat;
    int pairing_s2_nsat;
    double pairing_s1_qtot;
    double pairing_s2_qtot;
    TNtuple* michel_output = new TNtuple("event_tree", "Michel Data", "x:y:z:flux:E:delvtx:delt:time_since_kicker:nsat:qtot:run:subrun:index");

    double EffLiveTime = 0.;
    TChain * track = new TChain("InfoTree");
    track->Add(input_filename.c_str());
    track->LoadTree(0);
    track->SetBranchAddress("TrigID", &TrigID);
    track->SetBranchAddress("pairing_s1Index", &pairing_s1Index);
    track->SetBranchAddress("pairing_s2Index", &pairing_s2Index);
    track->SetBranchAddress("pairing_s1RunNum", &pairing_s1RunNum);
    track->SetBranchAddress("pairing_s2RunNum", &pairing_s2RunNum);
    track->SetBranchAddress("pairing_s1SubRunNum", &pairing_s1SubRunNum);
    track->SetBranchAddress("pairing_s2SubRunNum", &pairing_s2SubRunNum);
    track->SetBranchAddress("pairing_s1Energy", &pairing_s1Energy);
    track->SetBranchAddress("pairing_s1Flux", &pairing_s1Flux);
    track->SetBranchAddress("pairing_s1T", &pairing_s1T);
    track->SetBranchAddress("pairing_s2Energy", &pairing_s2Energy);
    track->SetBranchAddress("pairing_s2Flux", &pairing_s2Flux);
    track->SetBranchAddress("pairing_s2T", &pairing_s2T);
    track->SetBranchAddress("pairing_s1TV", &pairing_s1TV);
    track->SetBranchAddress("pairing_s1BV", &pairing_s1BV);
    track->SetBranchAddress("pairing_s2TV", &pairing_s2TV);
    track->SetBranchAddress("pairing_s2BV", &pairing_s2BV);
    track->SetBranchAddress("pairing_s1_timesincekicker", &pairing_s1_timesincekicker);

    track->SetBranchAddress("pairing_s1x", &pairing_s1x);
    track->SetBranchAddress("pairing_s1y", &pairing_s1y);
    track->SetBranchAddress("pairing_s1z", &pairing_s1z);
    track->SetBranchAddress("pairing_s2x", &pairing_s2x);
    track->SetBranchAddress("pairing_s2y", &pairing_s2y);
    track->SetBranchAddress("pairing_s2z", &pairing_s2z);
    track->SetBranchAddress("pairing_s1_nsat", &pairing_s1_nsat);
    track->SetBranchAddress("pairing_s2_nsat", &pairing_s2_nsat);
    track->SetBranchAddress("pairing_s1QTot", &pairing_s1_qtot);
    track->SetBranchAddress("pairing_s2QTot", &pairing_s2_qtot);
    track->SetBranchAddress("EffLiveTime", &EffLiveTime);

    TH1D * deltaT = new TH1D("deltaT",";#Delta t [us];Entries",100,0,30);

    const int EbinDiv = 1;
    // Prompt event energy histogram range & binning
    const double pEmax = 1000;
    const double pEmin = 0;
    const double pEbin = (pEmax-pEmin)*EbinDiv;
    TH1D* promptE = new TH1D("prompt_energy","Michel prompt Energy;Reco Energy[MeV];Rate/kHz/0.2MeV",pEbin,pEmin,pEmax);
    TH1D* promptE_acc = new TH1D("prompt_accidentals","Prompt Accidentals Energy;Energy[MeV];Entries",pEbin,pEmin,pEmax);
    TH1D* promptE_sub = new TH1D("prompt_energy_subtracted",";Energy[MeV];Entries",pEbin,pEmin,pEmax);

    // Prompt event energy histogram range & binning
    const double dEmax = 100;
    const double dEmin = 0;
    const double dEbin = (dEmax-dEmin)*EbinDiv;
    TH1D* delayedE = new TH1D("delayed_energy",";Reco Energy[MeV];Rate/kHz/0.2MeV",dEbin,dEmin,dEmax);
    TH1D* delayedE_acc = new TH1D("delayed_accidentals",";Energy[MeV];Entries",dEbin,dEmin,dEmax);
    TH1D* delayedE_sub = new TH1D("delayed_energy_subtracted",";Energy[MeV];Entries",dEbin,dEmin,dEmax);

    TH2D* h_energyComp = new TH2D("h_energyComp",";Charge /pe;Reco E /MeV", 100, 0, 2000, 50, 0, 50);
    TH2D* h_posR2Z = new TH2D("h_posR2Z", ";R^{2} /m^{2};Z /m", 100, 0, 4, 100, -2, 2);
    TH2D* h_posXY = new TH2D("h_posXY",";X /m;Y /m", 100, -2, 2, 100, -2, 2);
    TH1D* h_numHitPMT = new TH1D("h_numHitPMT",";# of Hit PMTs;", 100, 0, 100);

    int pair_count = 0;

    cout << " Total size : " << track->GetEntries() << endl;
    for(int i = 0; i < track->GetEntries(); i++){
        track->GetEntry(i);

        if( (i+1) % 50000 == 0 ) {
            cout << (i+1) << " done. " << endl;
        }

        double prompt_T =  pairing_s1T;
        double promptEnergy =  pairing_s1Energy;
        double prompt_x =  pairing_s1x;
        double prompt_y =  pairing_s1y;
        double prompt_z =  pairing_s1z;

        double delay_T = pairing_s2T;
        double delayEnergy = pairing_s2Energy;
        double delayFlux = pairing_s2Flux;
        double delay_x = pairing_s2x;
        double delay_y = pairing_s2y;
        double delay_z = pairing_s2z;
        int delayed_nsat = pairing_s2_nsat;
        double delayed_QTot = pairing_s2_qtot;
        int delayed_run_num = pairing_s2RunNum;
        int delayed_subrun_num = pairing_s2SubRunNum;
        int delayed_trigID = pairing_s2Index;

        double delT = (delay_T - prompt_T)/1000.; // Convert from nano-seconds to micro-seconds

        TVector3 pos_p( prompt_x, prompt_y, prompt_z );
        TVector3 pos_d( delay_x ,delay_y, delay_z);
        TVector3 diffVTX = pos_d - pos_p;

        if( delayEnergy > 20. && delayEnergy < 60.){
            deltaT->Fill( delT );
        }

        if( delT < MICHEL_DELTA_T_MAX && delT > MICHEL_DELTA_T_MIN  && diffVTX.Mag() < DELTA_VERTEX_THRESHOLD){
        //if( delT < MICHEL_DELTA_T_MAX && delT > MICHEL_DELTA_T_MIN  && diffVTX.Mag() < DELTA_VERTEX_THRESHOLD && pairing_s1_timesincekicker < 100e3){
            pair_count++;
            promptE->Fill(promptEnergy);
            delayedE->Fill(delayEnergy);
            const double delay_r = delay_x*delay_x + delay_y*delay_y;
            h_posR2Z->Fill(delay_r, delay_z);
            h_posXY->Fill(delay_x, delay_y);
            michel_output->Fill(delay_x, delay_y, delay_z, delayFlux, delayEnergy,
                                diffVTX.Mag(), delT, pairing_s1_timesincekicker, delayed_nsat, delayed_QTot,
                                delayed_run_num, delayed_subrun_num, delayed_trigID);
        }
        if( delT > ACCIDENTAL_DELTA_T_MIN && delT < ACCIDENTAL_DELTA_T_MAX){
            promptE_acc->Fill(promptEnergy);
            delayedE_acc->Fill(delayEnergy);
        }
    }

    double liveT = EffLiveTime*1e-9; // Scale from nano-seconds to seconds
    double scale_factor = 1/liveT;

    TH1D* TT = new TH1D("TT", "Live Time", 1, 0, 1);
    TT->SetBinContent(0, liveT);

    delayedE->Sumw2();
    delayedE_acc->Sumw2();
    TH1D* oriSub = (TH1D*)delayedE->Clone();
    oriSub->Add(delayedE_acc, -1);
    oriSub->SetName("NoScaleDelayedSub");

    promptE->Scale(scale_factor/1000.); // kHz
    promptE_acc->Scale(scale_factor/1000.); // kHz

    delayedE->Scale(scale_factor/1000.); // kHz
    delayedE_acc->Scale(scale_factor/1000.); // kHz

    TCanvas* can2 = new TCanvas("can2","can2",800,600);
    can2->SetGrid();
    can2->SetLogy();
    can2->cd();
    deltaT->Draw("E");
    deltaT->SetStats(0);
    deltaT->GetXaxis()->SetRangeUser(0,30);
    deltaT->SetLineColor(1);
    deltaT->SetLineWidth(2);
    TF1* fcn1 = new TF1("fcn1", "exp([0]-x/[1])",  0, 50);
    fcn1->SetParameter(0,5);
    fcn1->SetParameter(1,22);
    fcn1->SetLineColor(2);
    fcn1->SetLineStyle(9);
    deltaT->Fit(fcn1,"","R",2,8);
    TF1* fcn2 = new TF1("fcn2", "exp([0]-x/[1])",  0, 50);
    fcn2->SetParameter(0,1);
    fcn2->SetParameter(1,20);
    fcn2->SetLineColor(kGreen+1);
    fcn2->SetLineStyle(9);
    deltaT->Fit(fcn2,"","R",20,50);

    TF1* fcn_double = new TF1("fcn_double", fcn_2exp, 0, 50, 4);
    fcn_double->FixParameter(0, fcn1->GetParameter(0));
    fcn_double->FixParameter(1, fcn1->GetParameter(1));
    fcn_double->FixParameter(2, fcn2->GetParameter(0));
    fcn_double->FixParameter(3, fcn2->GetParameter(1));
    deltaT->Fit(fcn_double,"","R",0.5,30);
    fcn1->Draw("same");
    fcn2->Draw("same");
    fcn_double->Draw("same");
    fcn_double->SetLineColor(4);
    fcn_double->SetLineWidth(2);

    double AreaCorr = fcn1->Integral(1,10)/ fcn1->Integral(2,20);
    cout << " Area correction Factor : " << AreaCorr << endl;

    TCanvas * can4 = new TCanvas("can4","can4",1200,600);
    can4->Divide(2,1);

    can4->cd(1)->SetLogy();
    can4->cd(1)->SetGrid();
    promptE->Draw("hist");
    promptE->SetStats(0);
    promptE->GetYaxis()->SetTitleOffset(1.3);
    promptE->GetYaxis()->SetLabelSize(0.03);
    promptE->SetMarkerStyle(8);
    promptE->SetMarkerColor(1);
    promptE->SetMarkerSize(1);

    promptE_acc->Draw("histsame");
    promptE_acc->SetLineColor(2);
    promptE_acc->SetMarkerColor(2);
    promptE_acc->SetMarkerSize(1);
    promptE_acc->SetMarkerStyle(8);

    for(int i = 0; i < pEbin; i++){
        double a1 = promptE->GetBinContent(i);
        double a2 = promptE_acc->GetBinContent(i);
        promptE_sub->SetBinContent(i, a1-a2);
    }
    promptE_sub->Draw("histsame");
    promptE_sub->SetLineColor(4);
    promptE_sub->SetLineWidth(1);

    can4->cd(2)->SetLogy();
    can4->cd(2)->SetGrid();
    delayedE->Draw("PE+");
    delayedE->SetStats(0);
    delayedE->GetYaxis()->SetTitleOffset(1.8);
    delayedE->GetYaxis()->SetTitleSize(0.03);
    delayedE->GetYaxis()->SetLabelSize(0.025);
    delayedE->GetXaxis()->SetRangeUser(0,70);
    delayedE->SetMarkerStyle(20);
    delayedE->SetMarkerColor(1);
    delayedE->SetMarkerSize(1);

    delayedE_acc->Draw("samePE+");
    delayedE_acc->SetLineColor(2);
    delayedE_acc->SetMarkerColor(2);
    delayedE_acc->SetMarkerSize(1);
    delayedE_acc->SetMarkerStyle(20);

    // delayed_sub is empty, set it equal to the delayed spectrum
    // minus the accidental spectrum.
    delayedE_sub->Add(delayedE);
    delayedE_sub->Add( delayedE_acc, -1);

    delayedE_sub->Draw("samePE+");
    delayedE_sub->SetLineColor(4);
    delayedE_sub->SetLineWidth(2);
    delayedE_sub->SetMarkerColor(4);
    delayedE_sub->SetMarkerSize(1);
    delayedE_sub->SetMarkerStyle(20);

    // Integration regions to determine the prompt & delayed rate
    const double PROMPT_HIGH = 800; // MeV
    const double PROMPT_LOW = 10; // MeV
    const double DELAYED_HIGH = 60; // MeV
    const double DELAYED_LOW = 20; // MeV
    if(promptE->IsBinOverflow(PROMPT_HIGH)) {
        cout << "Warning! "<<PROMPT_HIGH <<" MeV is histogram overflow bin, rate estimates may not be accurate" << endl;
    }
    double RateP = promptE->Integral(promptE->FindBin(PROMPT_LOW), promptE->FindBin(PROMPT_HIGH)) * 1000.; // Hz
    double RateD = delayedE->Integral(delayedE->FindBin(DELAYED_LOW), delayedE->FindBin(DELAYED_HIGH)) * 1000.; // Hz
    double targetRateP = promptE_sub->Integral(promptE_sub->FindBin(PROMPT_LOW), promptE_sub->FindBin(PROMPT_HIGH)) * 1000.; // Hz
    double targetRateD = delayedE_sub->Integral(delayedE_sub->FindBin(DELAYED_LOW), delayedE_sub->FindBin(DELAYED_HIGH)) * 1000.; // Hz

    cout << "Total Effective Live Time : " << liveT << endl;
    cout << "Before accidental background subtraction" << endl;
    cout << "Prompt Rate[10 to 800 MeV] : " << RateP << " +- " << sqrt(RateP*liveT)/liveT <<  endl;
    cout << "Delayed Rate[20 to 60 MeV] : " << RateD << " +- " << sqrt(RateD*liveT)/liveT <<  endl;
    cout << endl << "After subtraction " << endl;
    cout << "Prompt Rate[10 to 800 MeV] : " << targetRateP << " +- " << sqrt(targetRateP*liveT)/liveT << endl;
    cout << "Delayed Rate[20 to 60 MeV] : " << targetRateD << " +- " << sqrt(targetRateD*liveT)/liveT << endl;

    TCanvas * tt = new TCanvas("tt","tt",1200,1000);
    tt->Divide(2,2);
    tt->cd(1)->SetGrid();
    tt->cd(1)->SetLogz();
    h_energyComp->Draw("colz");
    h_energyComp->SetStats(0);
    gStyle->SetPalette(55);

    tt->cd(2)->SetGrid();
    tt->cd(2)->SetLogy();
    h_numHitPMT->Draw();
    h_numHitPMT->SetStats(0);
    h_numHitPMT->SetLineColor(1);
    h_numHitPMT->SetLineWidth(2);

    tt->cd(3)->SetGrid();
    tt->cd(3)->SetLogz();
    h_posR2Z->Draw("colz");
    h_posR2Z->SetStats(0);
    gStyle->SetPalette(55);

    tt->cd(4)->SetGrid();
    tt->cd(4)->SetLogz();
    h_posXY->Draw("colz");
    h_posXY->SetStats(0);
    gStyle->SetPalette(55);

    TFile* output = new TFile(output_filename.c_str(), "RECREATE");
    michel_output->Write();
    deltaT->Write();
    promptE->Write();
    promptE_acc->Write();
    promptE_sub->Write();
    delayedE->Write();
    delayedE_acc->Write();
    delayedE_sub->Write();
    TT->Write();
    oriSub->Write();
    output->Close();
    return 0;
}
