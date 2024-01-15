#include <iostream>
#include <fstream>

Double_t gen_mc_prim(Double_t *x, Double_t *params) {
    // Function parameters
    Double_t alpha = params[0];
    Double_t s = params[1];
    Double_t ml = params[2];
    Double_t costheta = x[0];

    // Function calculations
    Double_t frac1 = TMath::Power(alpha, 2) / (4 * s);
    Double_t frac2 = 1 - (4 * TMath::Power(ml, 2)) / s;
    Double_t frac3 = 1 + (4 * TMath::Power(ml, 2)) / s;

    Double_t result = frac1 * TMath::Sqrt(frac2) * (frac2 + frac3 * TMath::Power(costheta, 2));

    return result;
}

Double_t gen_data_prim(Double_t *x, Double_t *params) {
    // Function parameters
    Double_t alpha = params[0];
    Double_t s = params[1];
    Double_t m_P = params[2];
    Double_t form_factor = params[3];
    Double_t costheta = x[0];

    // Function calculations
    Double_t sigma_P = TMath::Sqrt(1 - 4 * TMath::Power(m_P, 2)) / s;
    Double_t result = TMath::Power(alpha, 2) / (8 * s) * TMath::Power(sigma_P, 3) * TMath::Power(form_factor, 2) * (1 - TMath::Power(costheta, 2));

    return result;
}

void tester(){
    gROOT->SetBatch(kTRUE); // No popup window

    int nEvents_data = 10000;
    int nEvents_mc = 100000;

    int nTFractionFitter_repetitions = 1;

    Double_t P0 = 0.3;  // signal probability

	TH1F *data;     // full data histogram                              
    TH1F *real_bkg; // real background histogram
    TH1F *real_signal;  // real signal histogram

    TH1F *signal;   // signal mc
    TH1F *bkg;      // background mc
                               
    data = new TH1F("data", "Data angle distribution", 100, 0, 1);
    real_bkg = new TH1F("real_bkg", "real background angle distribution", 100, 0, 1);
    real_signal = new TH1F("real_signal", "real signal angle distribution", 100, 0, 1);

    bkg = new TH1F("signal", "generated background angle distribution", 100, 0, 1);
    signal = new TH1F("bkg", "generated signal angle distribution", 100, 0, 1);

    auto bkg_function = new TF1("gen_mc_prim", gen_mc_prim, 0., 1, 3);
    auto signal_function = new TF1("gen_data_prim", gen_data_prim, 0., 1, 4);

    bkg_function->SetParameters(0.2, 1, 0.105);   // in GeV
    signal_function->SetParameters(0.2, 1, 0.139, 1);


    // Calculating TfrctionFitter result multiple times
    std::ofstream outputFile("output_new.txt", std::ios::app);

    if (!outputFile.is_open()) {
        std::cerr << "Error opening the file!" << std::endl;
        exit(0);
    }
    outputFile << "BKG PROBABILITY;BKG ERROR;SIGNAL PROBABILITY;SIGNAL ERROR;CORR;MAN_CORR" << std::endl;


    for (int j = 0; j < nTFractionFitter_repetitions; j++) {

    // Data generation
    Double_t x, p;
    for (int i = 0; i < nEvents_data; ++i) {
        p = gRandom->Uniform();
		if( p<P0 ) { 
			x = signal_function->GetRandom();
			real_signal->Fill(x);
		}
		else { 
			x = bkg_function->GetRandom();
			real_bkg->Fill(x);
		}
        data->Fill(x);
    }

    // generate MC samples
	bkg->SetXTitle("bkg generated");
	bkg->SetLineColor(2);
	bkg->SetMarkerColor(2);
	bkg->SetMarkerStyle(24);
	bkg->SetMarkerSize(.7);
	for( Int_t i=0; i<nEvents_mc; i++) {
		bkg->Fill(bkg_function->GetRandom() ); 
	}

	signal->SetXTitle("signal generated");
	signal->SetLineColor(2);
	signal->SetMarkerColor(2);
	signal->SetMarkerStyle(24);
	signal->SetMarkerSize(.7);
	for( Int_t i=0; i<nEvents_mc; i++) {
		signal->Fill(signal_function->GetRandom() ); 
	}

    TCanvas *canvas = new TCanvas("canvas", "Data and MC Plots", 1200, 800);
    canvas->Divide(3, 1); // Divide canvas into three pads (two plots side by side)

    // Draw real mc1 histogram 
    canvas->cd(1);
    data->SetLineColor(1); 
    data->SetXTitle("cos theta");
    data->Draw();

    canvas->cd(2);
    real_signal->SetLineColor(2);
    real_signal->SetXTitle("cos theta"); 
    real_signal->Draw();

    canvas->cd(3);
    real_bkg->SetLineColor(3); 
    real_bkg->SetXTitle("cos theta");
    real_bkg->Draw();

    // Draw histograms for data
    TCanvas *canvas1 = new TCanvas("canvas1", "MC Plots", 1200, 800);
    canvas1->Divide(2, 1); // Divide canvas into three pads (two plots side by side)

    // Draw MC signal
    canvas1->cd(1);
    bkg->SetLineColor(2); 
    bkg->Draw();

    // Draw MC bkg
    canvas1->cd(2);
    signal->SetLineColor(3); 
    signal->Draw();

	// FractionFitter
	TObjArray *mc = new TObjArray(2);        // MC histograms are put in this array
	mc->Add(bkg);
	mc->Add(signal);
	TFractionFitter* fitter = new TFractionFitter(data, mc); // initialise
	fitter->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
	fitter->Constrain(0,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
	//fit->SetRangeX(1,15);                    // use only the first 15 bins in the fit
	TFitResultPtr fit_result = fitter->Fit();               // perform the fit
	cout << "fit status: " << fit_result << endl;
    Double_t p0, p1, err0, err1, corr_coeff;

    TH1F* result = (TH1F*) fitter->GetPlot();
    TCanvas *c1 = new TCanvas("c1", "Result Canvas", 1800, 800);
    c1->Divide(3, 1);
    c1->cd(1);
    result->Draw("C");
    c1->cd(2);
    data->Draw();
    // data->SetFillColor(kBlue);    // Blue color for data histogram
    // result->SetFillColor(kRed); 


    // TH1F* r = (TH1F*)data->Clone("r");  // Clone the data histogram
    // r->Add(result, -1); 

    auto rp1 = new TRatioPlot(data, result, "diffsig");
    c1->cd(3);
    rp1->Draw("B");
    // c1->Update();
    c1->SaveAs("res.png");

    // Added covariance matrix 
    TMatrixDSym corr = fit_result->GetCorrelationMatrix();
    TMatrixDSym cov = fit_result->GetCovarianceMatrix();

    cout << "first row: " << cov(0, 0) << ' ' << cov(0, 1) << endl;
    cout << "second row: " << cov(1, 0) << ' ' << cov(1, 1) << endl;

    double man_correlation = cov(1, 0) / sqrt(cov(0, 0) * cov(1, 1));
    cout << man_correlation << endl;

    corr_coeff = corr(0, 1);

    fitter->GetResult(0, p0, err0);
    fitter->GetResult(1, p1, err1);

    outputFile << p0 << ";" << err0 << ";";
    outputFile << p1 << ";" << err1 << ";";
    outputFile << corr_coeff << ";" << man_correlation << endl;

    // canvas1->SaveAs("mc_contribution.png");
    // canvas->SaveAs("data_contribution.png");
    }
    outputFile.close();
    std::cout << "Data appended to file successfully!" << std::endl;

    // delete data;
    // delete bkg;
    // delete signal;
    // delete real_bkg;
    // delete real_signal;


    // gApplication->Terminate();
}