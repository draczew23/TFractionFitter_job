float gen_mc(float alpha, float s, float ml, float costheta) {
    float frac1 = pow(alpha, 2) / (4*s);
    float frac2 = 1 - (4*pow(ml, 2)) / s;
    float frac3 = 1 + (4*pow(ml, 2)) / s;

    float result;
    result = frac1 * sqrt(frac2) * (frac2 + frac3 * pow(costheta, 2));

    return result;
}

float gen_data(float alpha, float s, float m_P, float form_factor, float costheta) {
    float result;
    float sigma_P = sqrt(1 - 4 * pow(m_P, 2)) / s;
    result = pow(alpha, 2) / (8*s) * pow(sigma_P, 3) * pow(form_factor, 2) * (1 - pow(costheta, 2));

    return result;
}

// ml - masa leptonu -> mion
// s -> 1 GeV
// alpha -> stala struktury subtelnej

// 1 + cos^2 theta z grubsza to nasz rozklad mc

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


void fiiter_macro(){

    gROOT->SetBatch(kTRUE); // No popup window 

    int nEvents = 10000;
    int nData = 1000;

    Double_t P0 = 0.5;
	Double_t P1 = 1 - P0;
    
    // contribution 0
	TF1 *f0 = new TF1("f0", "[0]*(1-cos(x))/TMath::Pi()", 0., TMath::Pi());
	f0->SetParameter(0,1.);
	f0->SetLineColor(2);
	Double_t int0 = f0->Integral( 0., TMath::Pi());

	// contribution 1
	TF1 *f1 = new TF1("f1", "[0]*(1-cos(x)*cos(x))*2./TMath::Pi()", 0., TMath::Pi());
	f1->SetParameter(0,1.);
	f1->SetLineColor(3);
	Double_t int1 = f1->Integral( 0., TMath::Pi());

    // Set up a random number generator
    TRandom randomGenerator;

	TH1F *data;                              //data histogram
	TH1F *true_mc0;
    TH1F *true_mc1;                               // first MC histogram
    TH1F *mc0;                               // first MC histogram
	TH1F *mc1;                               // second MC histogram

    data = new TH1F("data0", "Data angle distribution", 100, 0, TMath::Pi());
    true_mc0 = new TH1F("data1", "mc0 angle distribution", 100, 0, TMath::Pi());
    true_mc1 = new TH1F("data1", "mc0 angle distribution", 100, 0, TMath::Pi());

    // Data generation
    // Generate and fill the histogram with random numbers following a uniform distirbution
    Double_t p, x;
    for (int i = 0; i < nEvents; ++i) {
        p = gRandom->Uniform();
		if( p<P0 ) { 
			x = f0->GetRandom();
			true_mc0->Fill(x);
		}
		else { 
			x = f1->GetRandom();
			true_mc1->Fill(x);
		}
        data->Fill(x);
    }

	// generate MC samples
	mc0 = new TH1F("mc0", "MC sample 0 angle distribution", 100, 0, TMath::Pi());
	mc0->SetXTitle("x");
	mc0->SetLineColor(2);
	mc0->SetMarkerColor(2);
	mc0->SetMarkerStyle(24);
	mc0->SetMarkerSize(.7);
	for( Int_t i=0; i<nEvents; i++) {
		mc0->Fill( f0->GetRandom() ); 
	}

	mc1 = new TH1F("mc1", "MC sample 1 angle distribution", 100, 0, TMath::Pi());
	mc1->SetXTitle("x");
	mc1->SetLineColor(3);
	mc1->SetMarkerColor(3);
	mc1->SetMarkerStyle(24);
	mc1->SetMarkerSize(.7);
	for( Int_t i=0; i<nEvents; i++) {
		mc1->Fill( f1->GetRandom() ); 
	}

    // Draw histograms for data
    TCanvas *canvas = new TCanvas("canvas", "Data and MC Plots", 1200, 800);
    canvas->Divide(3, 1); // Divide canvas into three pads (two plots side by side)

    // Draw real mc1 histogram 
    canvas->cd(1);
    true_mc0->SetLineColor(2); 
    true_mc0->Draw();
    f0->Draw("same"); 

    // Draw real mc0 contribution
    canvas->cd(2);
    true_mc1->SetLineColor(3); 
    true_mc1->Draw();
    f1->Draw("same"); 
    
    // Draw data histogram
    canvas->cd(3);
    data->SetLineColor(4); // Set line color for data1
    data->Draw();

    // Save the canvas as an image file (optional)
    canvas->SaveAs("plots.png");

    // Draw histograms for data
    TCanvas *canvas1 = new TCanvas("canvas1", "Data and MC Plots", 1200, 800);
    canvas1->Divide(2, 1); // Divide canvas into three pads (two plots side by side)

    // Draw real mc1 histogram 
    canvas1->cd(1);
    mc0->SetLineColor(2); 
    mc0->Draw();
    f0->Draw("same"); 

    // Draw real mc0 contribution
    canvas1->cd(2);
    mc1->SetLineColor(3); 
    mc1->Draw();
    f1->Draw("same"); 

    canvas1->SaveAs("plots1.png");

	// FractionFitter
	TObjArray *mc = new TObjArray(2);        // MC histograms are put in this array
	mc->Add(mc0);
	mc->Add(mc1);
	TFractionFitter* fit = new TFractionFitter(data, mc); // initialise
	fit->Constrain(0,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
	fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
	//fit->SetRangeX(1,15);                    // use only the first 15 bins in the fit
	Int_t status = fit->Fit();               // perform the fit
	cout << "fit status: " << status << endl;

	// Display
	gStyle->SetOptStat(0);
	TCanvas c("c", "FractionFitter example", 700, 700);
	c.Divide(2,2);
	
	c.cd(1);
	f0->DrawClone();
	f0->GetHistogram()->SetTitle("Original MC distributions");
	f1->DrawClone("same");
	
	c.cd(2);
	data->SetTitle("Data distribution with true contributions");
	data->DrawClone("EP");
	true_mc0->Draw("same");
	true_mc1->Draw("same");


	c.cd(3);
	mc0->SetTitle("MC generated samples with fit predictions");
	mc0->Draw("PE");
	mc1->Draw("PEsame");

	TH1F *mcp0, *mcp1;
	if (status == 0) {                       // check on fit status
		mcp0 = (TH1F*)fit->GetMCPrediction(0);
		mcp0->SetLineColor(2);
		mcp0->Draw("same");
		mcp1 = (TH1F*)fit->GetMCPrediction(1);
		mcp1->SetLineColor(3);
		mcp1->Draw("same");
	}

    c.cd(4);
	Double_t p0, p1, errP0, errP1;
	TLatex l;
	l.SetTextSize(.035);
	Char_t texte[200];

if (status == 0) {                       // check on fit status
		TH1F* result = (TH1F*) fit->GetPlot();
		fit->GetResult( 0, p0, errP0);
		printf(" Parameter %d: true %.3f, estim. %.3f +/- %.3f\n", 0, P0, p0, errP0);
		fit->GetResult( 1, p1, errP1);
		printf(" Parameter %d: true %.3f, estim. %.3f +/- %.3f\n", 1, P1, p1, errP1);
		data->SetTitle("Data distribution with fitted contributions");
		data->DrawClone("Ep");
		result->Draw("same");
		f0->SetParameter(0,nEvents*p0/int0*data->GetBinWidth(1));
		f0->SetLineStyle(2);
		f0->DrawClone("same");
		f1->SetParameter(0,nEvents*p1/int1*data->GetBinWidth(1));
		f1->SetLineStyle(2);
		f1->DrawClone("same");
		sprintf( texte, "%d: true %.2f, estimated %.2f +/- %.2f\n", 0, P0, p0, errP0);
		l.DrawTextNDC( .45, .30, texte);
		sprintf( texte, "%d: true %.2f, estimated %.2f +/- %.2f\n", 1, P1, p1, errP1);
		l.DrawTextNDC( .45, .25, texte);
	}

    c.SaveAs("plots2.png");

    // Clean up
    delete canvas;
    delete data;
    delete true_mc0;
    delete true_mc1;
    delete f0;
    delete f1;

    // If you are running in batch mode, you can close the ROOT application
   gApplication->Terminate();
}
