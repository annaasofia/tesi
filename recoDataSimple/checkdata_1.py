import ROOT

ROOT.ROOT.EnableImplicitMT() 

file = 8650
filename = "recoDataSimple_" + str(file) + "_xtalMerging.root"
df = ROOT.RDataFrame("simpleEvent", filename)

ROOT.gStyle.SetOptStat(111111)
ROOT.gStyle.SetPalette(ROOT.kBird) 

# filtering
df_clean = df.Filter("SingleTrack == 1", "Single Tracks").Filter("Tracks.chi2_x < 5.0 && Tracks.chi2_y < 5.0", "Good Chi2")

# beam profile
h_beam = df_clean.Histo2D(("h_beam", "Beam Profile on Detector; d0_x; d0_y", 100, -20, 30, 100, -20, 15), "Tracks.d0_x", "Tracks.d0_y")
h_beam_angle = df_clean.Histo2D(("h_beam_angle", "Beam Divergence Profile; #theta_{in,x}; #theta_{in,y}", 500, -0.002, 0.002, 500, -0.002, 0.002), "Tracks.thetaIn_x", "Tracks.thetaIn_y")

# deflection angles
# Calcoliamo delta_theta usando i rami Tracks
df_physics = df_clean.Define("delta_theta_x", "Tracks.thetaOut_x - Tracks.thetaIn_x").Define("delta_theta_y", "Tracks.thetaOut_y - Tracks.thetaIn_y")
    
h_deflection_x = df_physics.Histo1D(("h_deflection_x", "Deflection Angle #Delta#theta_{x} = #theta_{out} - #theta_{in}; #Delta#theta_{x}; Events", 1000, -0.01, 0.01), "delta_theta_x")
h_deflection_y = df_physics.Histo1D(("h_deflection_y", "Deflection Angle #Delta#theta_{y} = #theta_{out} - #theta_{in}; #Delta#theta_{y}; Events", 1000, -0.01, 0.01), "delta_theta_y")
h_deflection_xy = df_physics.Histo2D(("h_deflection_xy", "Deflections x and y; #Delta#theta_{x}; #Delta#theta_{y}", 500, -0.0005, 0.0005, 500, -0.0005, 0.0005), "delta_theta_x", "delta_theta_y")
    
# angular scan
h_scan_x = df_physics.Histo2D(("h_scan_x", "#Delta#theta_{x} vs #theta_{In, x}; #theta_{In, x}; #Delta#theta_{x}", 500, -0.001, 0.001, 500, -0.002, 0.008), "Tracks.thetaIn_x", "delta_theta_x")
h_scan_y = df_physics.Histo2D(("h_scan_y", "#Delta#theta_{y} vs #theta_{In, y}; #theta_{In, y}; #Delta#theta_{y}", 500, -0.002, 0.002, 500, -0.002, 0.002), "Tracks.thetaIn_y", "delta_theta_y")


c1 = ROOT.TCanvas("c1", "analysis", 1600, 900)
c1.Divide(3, 2) 

c1.cd(1)
h_beam.Draw("COLZ") 

c1.cd(2)
ROOT.gPad.SetLogy()
h_deflection_x.Draw()

c1.cd(3)
ROOT.gPad.SetLogy()
h_deflection_y.Draw()

c1.cd(4)
h_beam_angle.Draw("COLZ")

c1.cd(5)
ROOT.gPad.SetLogz()
h_scan_x.Draw("SURF2")

c1.cd(6)
ROOT.gPad.SetLogz()
h_scan_y.Draw("SURF2")

#c1.cd(6)
#h_deflection_xy.Draw("COLZ")

c1.Update()

