import ROOT

# Assume you have data and MC histograms defined
data = ROOT.TH1F("data", "Data Histogram", 100, 0, 100)
mc0 = ROOT.TH1F("mc0", "MC Component 0", 100, 0, 100)
mc1 = ROOT.TH1F("mc1", "MC Component 1", 100, 0, 100)
mc2 = ROOT.TH1F("mc2", "MC Component 2", 100, 0, 100)

# Fill histograms with example data (replace this with your actual data)
data.FillRandom("gaus", 1000)
mc0.FillRandom("gaus", 500)
mc1.FillRandom("gaus", 300)
mc2.FillRandom("gaus", 200)

# Create a TObjArray for MC histograms
mc = ROOT.TObjArray(3)
mc.Add(mc0)
mc.Add(mc1)
mc.Add(mc2)

# Initialize TFractionFitter
fit = ROOT.TFractionFitter(data, mc)

# Constrain fraction 1 to be between 0 and 1 (optional)
fit.Constrain(1, 0.0, 1.0)

# Set fit range (optional)
fit.SetRangeX(1, 15)

# Perform the fit
status = fit.Fit()

# Check fit status
if status == 0:
    # Get the fitted result
    result = fit.GetPlot()

    # Draw the original data and fitted result
    data.Draw("Ep")
    result.Draw("same")

    # Display the plot
    ROOT.gPad.Modified()
    ROOT.gPad.Update()

    # Print fit fractions
    fractions = [fit.GetFitter().GetParameter(i) for i in range(mc.GetEntries())]
    print("Fit Fractions:", fractions)

    # Keep the plot window open
    ROOT.gApplication.Run()
else:
    print("Fit failed with status:", status)
