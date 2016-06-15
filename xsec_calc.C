#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TBranch.h"
#include "TF1.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TStyle.h"

using namespace std;

void xsec_calc()
{

  TFile *myfile = TFile::Open("NuMIFlux.root");
  if (myfile->IsOpen() == true){cout << "File opened" << endl;}

  TH1D * h1 = (TH1D*)myfile->Get("nueFluxHisto"); 

  const int mybins = 6;
  double edges [mybins+1] = {0.0, 0.4, 0.6, 0.8, 1.0, 1.4, 2.5};

  const int my_phi_bins = 8;
  double phi_edges [my_phi_bins+1] = {-180.0, -100.0, -60.0, -20.0, 0.0, 20.0, 60.0, 100.0, 180.0};

  const int my_theta_bins = 8;
  double theta_edges [my_theta_bins+1] = {-180.0/2., -100.0/2., -60.0/2., -20.0/2., 0.0, 20.0/2., 60.0/2., 100.0/2., 180.0/2.};

  TH1D *h_flux_vbin = new TH1D("flux_vbin", "NuMI nue flux, variable bins", mybins, edges);
  TH1D *h_nu_E = new TH1D("nu_E_vbin", "Truth Electron Neutrino Energy", mybins, edges);
  TH1D *h_xsec = new TH1D("nu_E_xsec", "Electron Neutrino Cross Section", mybins, edges);
  TH1D *h_e_eng = new TH1D("electron_E", "Electron Energy", mybins, edges);
  TH1D *h_phi = new TH1D("electron_phi", "Electron Phi", my_phi_bins, phi_edges);
  TH1D *h_theta = new TH1D("electron_theta", "Electron Theta", my_theta_bins, theta_edges);
  TH1D *h_xsec_e_eng = new TH1D("xsec_e_eng", "Differential Electron Neutrino Cross Section / electron energy", mybins, edges);
  TH1D *h_xsec_phi = new TH1D("xsec_phi", "Differential Electron Neutrino Cross Section / electron phi", my_phi_bins, phi_edges);
  TH1D *h_xsec_theta = new TH1D("xsec_theta", "Differential Electron Neutrino Cross Section / electron theta", my_theta_bins, theta_edges);

  THStack *error_stack_E = new THStack("error_stack_E", "Diff. Cross Section Error; Electron Energy [GeV]; d#sigma/dE Uncertainty [cm^2/GeV]");
  TH1D *error1_E = new TH1D("error1_E", "error1_E", mybins, edges);
  TH1D *error2_E = new TH1D("error2_E", "error2_E", mybins, edges);
  TH1D *error3_E = new TH1D("error3_E", "error3_E", mybins, edges);

  THStack *error_stack_theta = new THStack("error_stack_theta", "Diff. Cross Section Error; Electron #theta; d#sigma/d#theta Uncertainty [cm^2]");
  TH1D *error1_theta = new TH1D("error1_theta", "error1_theta", my_theta_bins, theta_edges);
  TH1D *error2_theta = new TH1D("error2_theta", "error2_theta", my_theta_bins, theta_edges);
  TH1D *error3_theta = new TH1D("error3_theta", "error3_theta", my_theta_bins, theta_edges);

  THStack *error_stack_phi = new THStack("error_stack_phi", "Diff. Cross Section Error; Electron #phi; d#sigma/d#phi Uncertainty [cm^2]");
  TH1D *error1_phi = new TH1D("error1_phi", "error1_phi", my_phi_bins, phi_edges);
  TH1D *error2_phi = new TH1D("error2_phi", "error2_phi", my_phi_bins, phi_edges);
  TH1D *error3_phi = new TH1D("error3_phi", "error3_phi", my_phi_bins, phi_edges);

  int nbins = h1->GetNbinsX();
  cout << "Number of bins: " << nbins << endl;

  TFile *myfile2 = TFile::Open("ana_hist.root");
  if (myfile2->IsOpen() == true){cout << "File opened" << endl;}

  TTree *tree = (TTree*)myfile2->Get("analysistree/anatree");
  TTree *tree2 = (TTree*)myfile2->Get("analysistree/pottree");


  /// 5822 nu events for 3.8322x10^19 POT w/ 225 nu_e, 167 nu_e CC
  /// NuMIFlux.root results for 6.6x10^20 POT/cm^2

  double measured;
  double flux;
  double xsec_val = 0;
  double total_xsec_val;
  double total_flux = 0;
  double neutrino_flux = 0;
  double total_measured;
  double total_pot = 0;
  double nu_e_counter = 0;
  //double ar_nucleons = 40; /// Naturally occuring is 39.948
  double ar_nucleons = 4.76*pow(10,31); //This is #nucleons

  float enu_truth[2];
  int nuPDG_truth[2];
  int mcevts_truth;
  int ccnc_truth[2];
  float lep_mom_truth[2];
  float lep_dcosx_truth[2];
  float lep_dcosy_truth[2];
  float lep_dcosz_truth[2];

  double pot;

  double mass_e = 0.000511; /// GeV

  tree->SetBranchAddress("enu_truth", enu_truth);
  tree->SetBranchAddress("nuPDG_truth", nuPDG_truth);
  tree->SetBranchAddress("mcevts_truth", &mcevts_truth);
  tree->SetBranchAddress("ccnc_truth", ccnc_truth);
  tree->SetBranchAddress("lep_mom_truth", lep_mom_truth);
  tree->SetBranchAddress("lep_dcosx_truth", lep_dcosx_truth);
  tree->SetBranchAddress("lep_dcosy_truth", lep_dcosy_truth);
  tree->SetBranchAddress("lep_dcosz_truth", lep_dcosz_truth);
  tree2->SetBranchAddress("pot", &pot);

  int nentries = tree->GetEntries();
  cout << "Number of entries: " << nentries << endl;

//  total_measured = 225./(3.86375*pow(10,19)); /// total # nu_e
  total_measured = 167./(3.86375*pow(10,19)); /// total # nu_e CC

  /// Fill histogram of neutrino energy with same binning as flux
  for(int i = 0; i < tree2->GetEntries(); i++){tree2->GetEntry(i); total_pot += pot;}

  for(int i = 0; i < nentries; i++)
  {
    tree->GetEntry(i);
    //total_pot += pot;
    for(int j = 0; j < mcevts_truth; j++)
    {
      if(nuPDG_truth[j] == 12 && ccnc_truth[j] == 0)
      {
	h_nu_E->Fill(enu_truth[j]);
	nu_e_counter++;

	double phi = atan2(lep_dcosx_truth[j], lep_dcosz_truth[j]);
	double theta = asin(lep_dcosy_truth[j]);
	double energy = sqrt((lep_mom_truth[j]*lep_mom_truth[j])+(mass_e*mass_e));

	h_e_eng->Fill(energy);
	h_phi->Fill(phi*180./3.141592654);
	h_theta->Fill(theta*180./3.141592654);

      }
    }
  }


  for(int i = 1; i < nbins; i++)
  {

    /// Flux for this energy bin [cm^2]^-1
    flux = h1->GetBinContent(i)/(6.6*pow(10,20));
    neutrino_flux += flux;
    h_flux_vbin->Fill(h1->GetBinLowEdge(i), flux);
  }
  for(int i = 1; i < mybins; i++)
  {
    /// Scale for proper POT - nue events
    measured = h_nu_E->GetBinContent(i)/(3.8322*pow(10,19));
    double v_flux = h_flux_vbin->GetBinContent(i);


    if(v_flux != 0){xsec_val = total_measured/(v_flux*ar_nucleons);}//*(edges[i+1]-edges[i])
    if(v_flux == 0){xsec_val = 0;}

    //total_flux += v_flux;

    h_xsec->Fill(h_flux_vbin->GetBinCenter(i), xsec_val);

  }/// End loop over flux histogram

  /// Differential Cross Section - electron energy
  for(int i = 1; i < h_e_eng->GetNbinsX()+1; i++)
  {
    /// For differential cross section - consider bin width too
    double diff_xsec_val = h_e_eng->GetBinContent(i)/(3.8322*pow(10,19)*neutrino_flux*ar_nucleons*(edges[i]-edges[i-1]));

    /// Stat error + flux error + atom error
    error1_E->SetBinContent(i, 0.02 * diff_xsec_val);
    error2_E->SetBinContent(i, 0.2 * diff_xsec_val);
    error3_E->SetBinContent(i, 1/sqrt(h_e_eng->GetBinContent(i))*diff_xsec_val);

    double binerror = diff_xsec_val*sqrt(pow(1/sqrt(h_e_eng->GetBinContent(i)),2)+0.04+0.004);

    h_xsec_e_eng->Fill(h_e_eng->GetBinCenter(i), diff_xsec_val);
    h_xsec_e_eng->SetBinError(i,binerror);
  }


  /// Differential Cross Section - electron phi
  for(int i = 1; i < h_phi->GetNbinsX()+1; i++)
  {
    /// For differential cross section - consider bin width too
    double diff_xsec_val = h_phi->GetBinContent(i)/(3.8322*pow(10,19)*neutrino_flux*ar_nucleons*(phi_edges[i]-phi_edges[i-1]));

    /// Stat error + flux error + atom error
    error1_phi->SetBinContent(i, 0.02 * diff_xsec_val);
    error2_phi->SetBinContent(i, 0.2 * diff_xsec_val);
    error3_phi->SetBinContent(i, 1/sqrt(h_e_eng->GetBinContent(i))*diff_xsec_val);

    double binerror = diff_xsec_val*sqrt(pow(1/sqrt(h_phi->GetBinContent(i)),2)+0.04+0.004);
    h_xsec_phi->Fill(h_phi->GetBinCenter(i), diff_xsec_val);
    h_xsec_phi->SetBinError(i,binerror);
  }

  /// Differential Cross Section - electron theta
  for(int i = 1; i < h_theta->GetNbinsX()+1; i++)
  {
    /// For differential cross section - consider bin width too
    double diff_xsec_val = h_theta->GetBinContent(i)/(3.8322*pow(10,19)*neutrino_flux*ar_nucleons*(theta_edges[i]-theta_edges[i-1]));

    /// Stat error + flux error + atom error
    error1_theta->SetBinContent(i, 0.02 * diff_xsec_val);
    error2_theta->SetBinContent(i, 0.2 * diff_xsec_val);
    error3_theta->SetBinContent(i, 1/sqrt(h_e_eng->GetBinContent(i))*diff_xsec_val);

    double binerror = diff_xsec_val*sqrt(pow(1/sqrt(h_theta->GetBinContent(i)),2)+0.04+0.004);
    h_xsec_theta->Fill(h_theta->GetBinCenter(i), diff_xsec_val);
    h_xsec_theta->SetBinError(i,binerror);
  }

  //total_xsec_val = total_measured/total_flux;
  cout << "Num nu_e: " << nu_e_counter << endl;
  cout << "Total POT: " << total_pot << endl;
  cout << "Total Measured: " << total_measured << ", Total Flux: " << neutrino_flux << endl;

  total_xsec_val = total_measured/neutrino_flux;
  cout << "Total nue CC xsec: " << total_xsec_val/ar_nucleons << endl;

  TCanvas *can1 = new TCanvas();
  can1->cd();
  h1->Draw();
  can1->Print("nueFluxHisto.pdf");

  TCanvas *can2 = new TCanvas();
  can2->cd();
  h_nu_E->Draw();
  can2->Print("nueEnergyHisto.pdf");

  TCanvas *can3 = new TCanvas();
  can3->cd();	
  h_xsec->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
  h_xsec->Draw();
  can3->Print("nueXsecHisto.pdf");

  TCanvas *can4 = new TCanvas();
  can4->cd();
  h_xsec_e_eng->GetXaxis()->SetTitle("Electron Energy [GeV]");
  h_xsec_e_eng->GetYaxis()->SetTitle("d#sigma/dE [cm^2/GeV]");
  h_xsec_e_eng->Draw("E1");
  can4->Print("nueXsec_e_eng.pdf");

  TCanvas *can5 = new TCanvas();
  can5->cd();
  h_xsec_phi->GetXaxis()->SetTitle("Electron #phi");
  h_xsec_phi->GetYaxis()->SetTitle("d#sigma/d#phi [cm^2]");
  h_xsec_phi->Draw("E1");
  can5->Print("nueXsec_phi.pdf");

  TCanvas *can6 = new TCanvas();
  can6->cd();
  h_xsec_theta->GetXaxis()->SetTitle("Electron #theta");
  h_xsec_theta->GetYaxis()->SetTitle("d#sigma/d#theta [cm^2]");
  h_xsec_theta->Draw("E1");
  can6->Print("nueXsec_theta.pdf");

  TCanvas *can7 = new TCanvas();
  can7->cd();
  error1_E->SetFillColor(9);
  error_stack_E->Add(error1_E);
  error2_E->SetFillColor(8);
  error_stack_E->Add(error2_E);
  error3_E->SetFillColor(46);
  error_stack_E->Add(error3_E);
  error_stack_E->Draw();
  TLegend *leg1 = new TLegend(0.65, 0.65, 0.85, 0.85);
  leg1->SetHeader("Error Sources");
  leg1->AddEntry(error1_E, "Atomic", "f");
  leg1->AddEntry(error2_E, "Flux", "f");
  leg1->AddEntry(error3_E, "Stat.", "f");
  leg1->Draw();
  can7->Print("nueXsec_e_eng_error.pdf");

  TCanvas *can8 = new TCanvas();
  can8->cd();
  error1_theta->SetFillColor(9);
  error_stack_theta->Add(error1_theta);
  error2_theta->SetFillColor(8);
  error_stack_theta->Add(error2_theta);
  error3_theta->SetFillColor(46);
  error_stack_theta->Add(error3_theta);
  error_stack_theta->Draw();
  TLegend *leg2 = new TLegend(0.65, 0.65, 0.85, 0.85);
  leg2->SetHeader("Error Sources");
  leg2->AddEntry(error1_theta, "Atomic", "f");
  leg2->AddEntry(error2_theta, "Flux", "f");
  leg2->AddEntry(error3_theta, "Stat.", "f");
  leg2->Draw();
  can8->Print("nueXsec_e_theta_error.pdf");

  TCanvas *can9 = new TCanvas();
  can9->cd();
  error1_phi->SetFillColor(9);
  error_stack_phi->Add(error1_phi);
  error2_phi->SetFillColor(8);
  error_stack_phi->Add(error2_phi);
  error3_phi->SetFillColor(46);
  error_stack_phi->Add(error3_phi);
  error_stack_phi->Draw();
  TLegend *leg3 = new TLegend(0.65, 0.65, 0.85, 0.85);
  leg3->SetHeader("Error Sources");
  leg3->AddEntry(error1_phi, "Atomic", "f");
  leg3->AddEntry(error2_phi, "Flux", "f");
  leg3->AddEntry(error3_phi, "Stat.", "f");
  leg3->Draw();
  can9->Print("nueXsec_e_phi_error.pdf");

}

int main(){xsec_calc();}
