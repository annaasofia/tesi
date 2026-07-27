import ROOT
import math
import sys
import numpy as np
from array import array

ROOT.ROOT.EnableImplicitMT() 
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kBird)

def get_run_parameters(file_id):
    if file_id in [8430, 8431]:
        parameters = {
            "theta_L": 12.992,
            "deflection_peak": 6010.0,
            "width": 12.8,
            "height": 2,
            "delta_x": 0.5 * 74 * 6.01 * pow(10, -3),  # 0.222 mm
            "max_value": 8000
        }
    elif file_id in [8650]:
        parameters = {
            "theta_L": 14.294,
            "deflection_peak": 6080.0,
            "width": 12.8,
            "height": 2,
            "delta_x": 0.5 * 74 * 6.01 * pow(10, -3),  # 0.222 mm
            "max_value": 8000
        }
    elif file_id in [8655, 8656]:
        parameters = {
            "theta_L": 12.992,
            "deflection_peak": 6130.0,
            "width": 12.8,
            "height": 2,
            "delta_x": 0.5 * 74 * 6.13 * pow(10, -3),  # 0.222 mm
            "max_value": 8000
        }
    else:
        print(f"Run {file_id} not found...")
        sys.exit(1)
    return parameters

def df_convert_to_urad_define_deltatheta(df):
    return df.Define("thetaIn_x", "Tracks.thetaIn_x * 1e6").Define("Deltatheta_x", "(Tracks.thetaOut_x - Tracks.thetaIn_x) * 1e6")

def filter_message(number, count_before, count_after):
    return f"Filter {number}: {((count_before - count_after) / count_before * 100):.2f}% events out of {count_before} got discarded."

def preliminary_cut_on_deltatheta(df, deflection_peak):
    # Applies a preliminary cut on Deltatheta_x to select channeled particles.
    h_defl_before = df.Histo1D(("h_defl_before", "Angular Deflection; #Delta#theta_{x} [#murad]; Counts", 500, -2000, deflection_peak + 1000), "Deltatheta_x") 
    h_defl_val = h_defl_before.GetValue()

    pre_fit = ROOT.TF1("pre_fit", "gaus", deflection_peak - 500, deflection_peak + 500)
    h_defl_val.Fit(pre_fit, "RQ0") # R=Range, Q=Quiet, 0=NoDraw
    pre_mean = pre_fit.GetParameter(1); pre_sigma = pre_fit.GetParameter(2)
    preliminary_cut = pre_mean - 3.0 * pre_sigma
    df_cut = df.Filter(f"Deltatheta_x > {preliminary_cut}", "Preliminary cut on Deltatheta_x to select channeled particles")
    return df_cut, preliminary_cut

def compute_spatial_cut_bounds(df_phys, parameters):
    # histograms pre spatial cut
    h_d0_xy = df_phys.Histo2D(("h_d0_xy", "Incoming d0_x vs d0_y of all particles; d0_x [mm]; d0_y [mm]", 200, -10, 20, 200, -20, 20), "Tracks.d0_x", "Tracks.d0_y")
    h_d0_Out_xy = df_phys.Histo2D(("h_d0_Out_xy", "Outgoing d0_x vs d0_y of all particles; d0_x [mm]; d0_y [mm]", 200, -10, 20, 200, -20, 20), "Tracks.d0Out_x", "Tracks.d0Out_y")
    
    # c0 = ROOT.TCanvas("c0", "d0_x and d0_y of all particles", 1600, 900)
    # c0.Divide(2, 2)
    # c0.cd(1); h_d0_xy.Draw("COLZ")
    # x1, y1, x2, y2 = -2.0, -6.0, 2.0, 8.0
    # box1 = ROOT.TBox(x1, y1, x2, y2); box1.SetLineColor(ROOT.kMagenta); box1.SetLineWidth(2); box1.SetFillStyle(0); box1.Draw("SAME")
    # c0.cd(2); h_d0_Out_xy.Draw("COLZ")
    # box2 = ROOT.TBox(x1, y1, x2, y2); box2.SetLineColor(ROOT.kMagenta); box2.SetLineWidth(2); box2.SetFillStyle(0); box2.Draw("SAME")

    # i choose a cut on Deltatheta_x to select only channeled particles
    df_cut, _ = preliminary_cut_on_deltatheta(df_phys, parameters["deflection_peak"])

    # i look at the beam profile of d0 of those who channeled, and find the region (width x length)
    h_d0_xy_ch = df_cut.Histo2D(("h_d0_xy_ch", "Incoming beam - channeled particles; d0_x [mm]; d0_y [mm]", 500, -2, 3, 200, -6, 8), "Tracks.d0_x", "Tracks.d0_y")
    h_d0_Out_xy_ch = df_cut.Histo2D(("h_d0_Out_xy_ch", "Outgoing beam - channeled particles; d0Out_x [mm]; d0Out_y [mm]", 500, -2, 3, 200, -6, 8), "Tracks.d0Out_x", "Tracks.d0Out_y")

    h_d0_x_ch = df_cut.Histo1D(("h_d0_x_ch", "d0_x of channeled particles; d0_x [mm]; Counts", 500, -2, 3), "Tracks.d0_x")
    h_d0_y_ch = df_cut.Histo1D(("h_d0_y_ch", "d0_y of channeled particles; d0_y [mm]; Counts", 200, -6, 8), "Tracks.d0_y")
    h_d0_Out_x_ch = df_cut.Histo1D(("h_d0_Out_x_ch", "d0Out_x of channeled particles; d0Out_x [mm]; Counts", 500, -2, 3), "Tracks.d0Out_x")
    h_d0_Out_y_ch = df_cut.Histo1D(("h_d0_Out_y_ch", "d0Out_y of channeled particles; d0Out_y [mm]; Counts", 200, -6, 8), "Tracks.d0Out_y")

    y_c_in = h_d0_y_ch.GetMean(); y_c_out = h_d0_Out_y_ch.GetMean(); y_c = (y_c_in + y_c_out) / 2.0 # crystal centre (mm)
    y_min, y_max = y_c - parameters["width"] / 2.0, y_c + parameters["width"] / 2.0

    # sliding window method - maximize the integral in order to find x edges
    h_x = h_d0_Out_x_ch.GetValue()
    max_particles = -1
    best_x_min = 0
    for i in range(1, h_x.GetNbinsX() + 1):
        current_x_start = h_x.GetBinLowEdge(i)
        current_x_end = current_x_start + parameters["height"]
        
        bin_start = i
        bin_end = h_x.FindBin(current_x_end)
        
        particles_in_window = h_x.Integral(bin_start, bin_end)
        
        if particles_in_window > max_particles:
            max_particles = particles_in_window
            best_x_min = current_x_start

    best_x_max = best_x_min + parameters["height"]
    # print(f"\tCut in x found: [{best_x_min:.2f} mm, {best_x_max:.2f} mm]")

    x_min, x_max = best_x_min - parameters["delta_x"], best_x_max - parameters['delta_x']
    print(f"\tFinal cut in x (after shift of {parameters['delta_x']:.3f} mm): x = [{x_min:.4f} mm, {x_max:.4f} mm], y = [{y_min:.4f} mm, {y_max:.4f} mm]")

    # c0.cd(3); h_d0_xy_ch.Draw("COLZ")
    # box3 = ROOT.TBox(x_min, y_min, x_max, y_max); box3.SetLineColor(ROOT.kRed); box3.SetLineWidth(2); box3.SetFillStyle(0); box3.Draw("SAME")
    # c0.cd(4); h_d0_Out_xy_ch.Draw("COLZ")
    # box4 = ROOT.TBox(x_min+parameters["delta_x"], y_min, x_max+parameters["delta_x"], y_max); box4.SetLineColor(ROOT.kRed); box4.SetLineWidth(2); box4.SetFillStyle(0); box4.Draw("SAME")
    # box5 = ROOT.TBox(x_min, y_min, x_max, y_max); box5.SetLineColor(ROOT.kRed); box5.SetLineWidth(1); box5.SetLineStyle(2); box5.SetFillStyle(0); box5.Draw("SAME")
    # c0.Update()

    # c0b = ROOT.TCanvas("c0b", "d0_x and d0_y of channeled particles", 1200, 900)
    # h_d0_xy.Draw("COLZ"); h_d0_xy.SetTitle("Incoming beam of all particles; d0_x [mm]; d0_y [mm]")
    # box_selected = ROOT.TBox(x_min, y_min, x_max, y_max); box_selected.SetLineColor(ROOT.kRed); box_selected.SetLineWidth(2); box_selected.SetFillStyle(0); box_selected.Draw("SAME")
    # c0b.Update()

    # c1 = ROOT.TCanvas("c1", "d0_x and d0_y of channeled particles", 1000, 900)
    # pad_center = ROOT.TPad("pad_center", "pad_center", 0, 0, 0.65, 0.65)
    # pad_center.Draw()
    # pad_top = ROOT.TPad("pad_top", "pad_top", 0.0, 0.60, 0.65, 1.0)
    # pad_top.Draw()
    # pad_right = ROOT.TPad("pad_right", "pad_right", 0.60, 0.0, 1.0, 0.65)
    # pad_right.Draw()
    # h_d0_Out_xy_ch.SetTitle(""); pad_center.cd(); h_d0_Out_xy_ch.Draw("COL"); box4.Draw("SAME")
    # h_d0_Out_x_ch.SetTitle(""); h_d0_Out_x_ch.GetXaxis().SetTitle(""); pad_top.cd(); h_d0_Out_x_ch.SetFillColor(ROOT.kAzure-3); h_d0_Out_x_ch.Draw("BAR X+")
    # l1 = ROOT.TLine(x_min+parameters["delta_x"], 0, x_min+parameters["delta_x"], h_d0_Out_x_ch.GetMaximum()); l1.SetLineColor(ROOT.kRed); l1.SetLineStyle(2); l1.SetLineWidth(2); l1.Draw("SAME")
    # l2 = ROOT.TLine(x_max+parameters["delta_x"], 0, x_max+parameters["delta_x"], h_d0_Out_x_ch.GetMaximum()); l2.SetLineColor(ROOT.kRed); l2.SetLineStyle(2); l2.SetLineWidth(2); l2.Draw("SAME")
    # h_d0_Out_y_ch.SetTitle(""); h_d0_Out_y_ch.GetXaxis().SetTitle(""); pad_right.cd(); h_d0_Out_y_ch.SetFillColor(ROOT.kAzure-3); h_d0_Out_y_ch.Draw("HBAR Y+")
    # l3 = ROOT.TLine(0, y_min, h_d0_Out_y_ch.GetMaximum(), y_min); l3.SetLineColor(ROOT.kRed); l3.SetLineStyle(2); l3.SetLineWidth(2); l3.Draw("SAME")
    # l4 = ROOT.TLine(0, y_max, h_d0_Out_y_ch.GetMaximum(), y_max); l4.SetLineColor(ROOT.kRed); l4.SetLineStyle(2); l4.SetLineWidth(2); l4.Draw("SAME")
    # c1.Update()

    # c2 = ROOT.TCanvas("c2", "d0_x and d0_y of channeled particles", 1000, 900)
    # pad_center = ROOT.TPad("pad_center", "pad_center", 0, 0, 0.65, 0.65)
    # pad_center.Draw()
    # pad_top = ROOT.TPad("pad_top", "pad_top", 0.0, 0.60, 0.65, 1.0)
    # pad_top.Draw()
    # pad_right = ROOT.TPad("pad_right", "pad_right", 0.60, 0.0, 1.0, 0.65)
    # pad_right.Draw()
    # h_d0_xy_ch.SetTitle(""); pad_center.cd(); h_d0_xy_ch.Draw("COL"); box3.Draw("SAME")
    # h_d0_x_ch.SetTitle(""); h_d0_x_ch.GetXaxis().SetTitle(""); pad_top.cd(); h_d0_x_ch.SetFillColor(ROOT.kAzure-3); h_d0_x_ch.Draw("BAR X+"); 
    # l5 = ROOT.TLine(x_min, 0, x_min, h_d0_x_ch.GetMaximum()); l5.SetLineColor(ROOT.kRed); l5.SetLineStyle(2); l5.SetLineWidth(2); l5.Draw("SAME")
    # l6 = ROOT.TLine(x_max, 0, x_max, h_d0_x_ch.GetMaximum()); l6.SetLineColor(ROOT.kRed); l6.SetLineStyle(2); l6.SetLineWidth(2); l6.Draw("SAME")
    # h_d0_y_ch.SetTitle(""); h_d0_y_ch.GetXaxis().SetTitle(""); pad_right.cd(); h_d0_y_ch.SetFillColor(ROOT.kAzure-3); h_d0_y_ch.Draw("HBAR Y+"); l3.Draw("SAME"); l4.Draw("SAME")
    # c2.Update()

    return x_min, x_max, y_min, y_max

def torsion_map(df, parameters, x_min, x_max, y_min, y_max, x_cut_margin, y_cut_margin, nx_slices, ny_slices, restricted=False, linear=False, chosen_model="full_quadratic"):

    # Create 3D histograms (one for all and one for channeled only) to avoid looping over RDataFrame
    _, preliminary_cut = preliminary_cut_on_deltatheta(df, parameters["deflection_peak"])
    h3_all = df.Histo3D(("h3_all", "", 1000, -100, 100, nx_slices, x_min, x_max, ny_slices, y_min, y_max), "thetaIn_x", "Tracks.d0_x", "Tracks.d0_y").GetValue()
    h3_chan = df.Filter(f"Deltatheta_x > {preliminary_cut}").Histo3D(("h3_chan", "", 1000, -100, 100, nx_slices, x_min, x_max, ny_slices, y_min, y_max), "thetaIn_x", "Tracks.d0_x", "Tracks.d0_y").GetValue()

    # To store the final map
    h2_torsion_map = ROOT.TH2D("h2_torsion_map", "2D Torsion Map (bins 0.1 x 0.1 mm); x at crystal surface [mm]; y at crystal surface [mm]; angle shift #theta_{0} [#murad]", nx_slices, x_min, x_max, ny_slices, y_min, y_max)
    h2_eff_map = ROOT.TH2D("h2_eff_map", "2D Efficiency Map; x [mm]; y [mm]; Local Efficiency [%]", nx_slices, x_min, x_max, ny_slices, y_min, y_max)

    h2_torsion_map.Sumw2()
    h2_eff_map.Sumw2()

    scan_min, scan_max, step = -80, 80, 1.0  # Scan range and step for theta_in

    # Extract local theta_0 for each (x, y) grid square
    for ix in range(1, nx_slices + 1):
        for iy in range(1, ny_slices + 1):
            
            # Project 3D histo onto 1D theta_in axis for this specific (x,y) bin
            h1_all = h3_all.ProjectionX(f"all_proj_{ix}_{iy}", ix, ix, iy, iy)
            h1_chan = h3_chan.ProjectionX(f"chan_proj_{ix}_{iy}", ix, ix, iy, iy)
            
            best_theta_0_raw = 0.0
            max_eff_raw = -1.0
            
            theta_vals, eff_vals = [], []
            current_theta = scan_min
            
            # Angular scan for this specific square
            while current_theta <= scan_max:
                bin_min = h1_all.FindBin(current_theta - parameters["theta_L"] / 2.0)
                bin_max = h1_all.FindBin(current_theta + parameters["theta_L"] / 2.0)
                
                n_tot_test = h1_all.Integral(bin_min, bin_max)
                n_ch_test  = h1_chan.Integral(bin_min, bin_max)
                
                eff_test = 0.0
                if n_tot_test > 100: # minimum statistics
                    eff_test = n_ch_test / n_tot_test * 100.0
                    if eff_test > max_eff_raw:
                        max_eff_raw = eff_test
                        best_theta_0_raw = current_theta
                        
                theta_vals.append(current_theta)
                eff_vals.append(eff_test)
                current_theta += step
                
            # Fit the local efficiency peak
            if max_eff_raw > 0:
                gr_eff_slice = ROOT.TGraph(len(theta_vals), array('d', theta_vals), array('d', eff_vals))
                gaus_eff_slice = ROOT.TF1(f"gaus_eff_{ix}_{iy}", "gaus", best_theta_0_raw - 5.0, best_theta_0_raw + 20.0)
                gaus_eff_slice.SetParameters(max_eff_raw, best_theta_0_raw, 5.0)
                gaus_eff_slice.SetParLimits(1, -80.0, 80.0) # Limit the mean to be within the scan range
                
                gr_eff_slice.Fit(gaus_eff_slice, "RQ0")
                
                local_theta_0 = gaus_eff_slice.GetParameter(1)
                local_theta_0_err = gaus_eff_slice.GetParError(1)
                local_max_eff = gaus_eff_slice.GetParameter(0)
                local_max_eff_err = gaus_eff_slice.GetParError(0)
                
                # Fill the 2D map
                h2_torsion_map.SetBinContent(ix, iy, local_theta_0)
                h2_torsion_map.SetBinError(ix, iy, local_theta_0_err)
                h2_eff_map.SetBinContent(ix, iy, local_max_eff)
                h2_eff_map.SetBinError(ix, iy, local_max_eff_err)

    c_torsion_2d = ROOT.TCanvas("c_torsion_2d", "2D Torsion Map", 900, 700)
    c_torsion_2d.SetRightMargin(0.15) # room for colorbar
    h2_torsion_map.SetStats(0)
    h2_torsion_map.GetZaxis().SetRangeUser(-50, 20) # Adjust based on expected torsion range
    h2_torsion_map.Draw("surf2")
    c_torsion_2d.Update()
    
    # Fit the Torsion Map with a xy function
    current_x_min = (x_min + x_cut_margin) if restricted else x_min
    current_x_max = (x_max - x_cut_margin) if restricted else x_max
    current_y_min = (y_min + y_cut_margin) if restricted else y_min
    current_y_max = (y_max - y_cut_margin) if restricted else y_max

    fit_models = {
        "linear": "[0] + [1]*x + [2]*y",
        "pure_parabolic_y": "[0] + [1]*y + [2]*y*y",
        "parabolic_y": "[0] + [1]*x + [2]*y + [3]*y*y",
        "full_quadratic": "[0] + [1]*x + [2]*y + [3]*x*x + [4]*y*y + [5]*x*y"
    }
    torsion_fit_2d = ROOT.TF2("torsion_fit_2d", fit_models[chosen_model], current_x_min, current_x_max, current_y_min, current_y_max)
    fit_result = h2_torsion_map.Fit(torsion_fit_2d, "RQ0S")

    # chi2 = fit_result.Chi2()
    # ndf = fit_result.Ndf()        # Number of Degrees of Freedom
    # p_value = fit_result.Prob()   # Fit probability
    # print(f"Chi2 / NDF: {chi2 / ndf:.2f}")
    # print(f"p-value: {p_value:.4f}")

    match chosen_model:
        case "linear":
            theta_0_baseline = torsion_fit_2d.GetParameter(0)
            tau_x = torsion_fit_2d.GetParameter(1)
            tau_y = torsion_fit_2d.GetParameter(2)
            theta_0_baseline_err = torsion_fit_2d.GetParError(0)
            tau_x_err = torsion_fit_2d.GetParError(1)
            tau_y_err = torsion_fit_2d.GetParError(2)

            
            fit_params = (theta_0_baseline, tau_x, tau_y)
            fit_errors = (theta_0_baseline_err, tau_x_err, tau_y_err)
            surface_expr = f"({theta_0_baseline} + {tau_x}*Tracks.d0_x + {tau_y}*Tracks.d0_y)"

        case "parabolic_y":
            theta_0_baseline = torsion_fit_2d.GetParameter(0)
            tau_x = torsion_fit_2d.GetParameter(1)
            tau_y = torsion_fit_2d.GetParameter(2)
            c_yy = torsion_fit_2d.GetParameter(3)
            theta_0_baseline_err = torsion_fit_2d.GetParError(0)
            tau_x_err = torsion_fit_2d.GetParError(1)
            tau_y_err = torsion_fit_2d.GetParError(2)
            c_yy_err = torsion_fit_2d.GetParError(3)
            
            fit_params = (theta_0_baseline, tau_x, tau_y, c_yy)
            fit_errors = (theta_0_baseline_err, tau_x_err, tau_y_err, c_yy_err)
            surface_expr = f"({theta_0_baseline} + {tau_x}*Tracks.d0_x + {tau_y}*Tracks.d0_y + {c_yy}*Tracks.d0_y*Tracks.d0_y)"

        case "pure_parabolic_y":
            theta_0_baseline = torsion_fit_2d.GetParameter(0)
            tau_y = torsion_fit_2d.GetParameter(1)
            c_yy = torsion_fit_2d.GetParameter(2)
            theta_0_baseline_err = torsion_fit_2d.GetParError(0)
            tau_y_err = torsion_fit_2d.GetParError(1)
            c_yy_err = torsion_fit_2d.GetParError(2)
            
            fit_params = (theta_0_baseline, tau_y, c_yy)
            fit_errors = (theta_0_baseline_err, tau_y_err, c_yy_err)
            surface_expr = f"({theta_0_baseline} + {tau_y}*Tracks.d0_y + {c_yy}*Tracks.d0_y*Tracks.d0_y)"

        case "full_quadratic":
            theta_0_baseline = torsion_fit_2d.GetParameter(0)
            tau_x = torsion_fit_2d.GetParameter(1)
            tau_y = torsion_fit_2d.GetParameter(2)
            c_xx = torsion_fit_2d.GetParameter(3)
            c_yy = torsion_fit_2d.GetParameter(4)
            c_xy = torsion_fit_2d.GetParameter(5)
            theta_0_baseline_err = torsion_fit_2d.GetParError(0)
            tau_x_err = torsion_fit_2d.GetParError(1)
            tau_y_err = torsion_fit_2d.GetParError(2)
            c_xx_err = torsion_fit_2d.GetParError(3) 
            c_yy_err = torsion_fit_2d.GetParError(4)
            c_xy_err = torsion_fit_2d.GetParError(5)

            
            fit_params = (theta_0_baseline, tau_x, tau_y, c_xx, c_yy, c_xy)
            fit_errors = (theta_0_baseline_err, tau_x_err, tau_y_err, c_xx_err, c_yy_err, c_xy_err)
            surface_expr = f"({theta_0_baseline} + {tau_x}*Tracks.d0_x + {tau_y}*Tracks.d0_y + {c_xx}*Tracks.d0_x*Tracks.d0_x + {c_yy}*Tracks.d0_y*Tracks.d0_y + {c_xy}*Tracks.d0_x*Tracks.d0_y)"
            
        case _:
            raise ValueError(f"Failed: model '{chosen_model}' not supported.")

    rdf_surface_expr = f"({theta_0_baseline} + {tau_y}*Tracks.d0_y)" if linear else surface_expr

    print(f"\tBaseline Theta_0: {theta_0_baseline:.2f} urad")
    print(f"\tLinear tau_x: {tau_x:.2f} urad/mm")
    print(f"\tLinear tau_y: {tau_y:.2f} urad/mm")
    # print(f"\tc_xx: {c_xx:.2f}")
    # print(f"\tc_yy: {c_yy:.2f}")
    # print(f"\tc_xy: {c_xy:.2f}")

    # torsion_plot_2d = ROOT.TF2("torsion_plot_2d", plot_formula, x_min, x_max, y_min, y_max)
    torsion_plot_2d = torsion_fit_2d.Clone("torsion_plot_2d")
    torsion_plot_2d.SetRange(current_x_min, current_y_min, current_x_max, current_y_max)
    torsion_plot_2d.SetTitle("Continuous 2D Torsion Map; x at crystal surface [mm]; y at crystal surface [mm]; angle shift #theta_{0} [#murad]")
    c_torsion_smooth = ROOT.TCanvas("c_torsion_smooth", "Continuous 2D Torsion Map", 900, 700)
    c_torsion_smooth.SetRightMargin(0.15)
    torsion_plot_2d.Draw("surf2") 
    c_torsion_smooth.Update()

    ROOT.SetOwnership(h2_torsion_map, False)
    ROOT.SetOwnership(torsion_plot_2d, False)
    ROOT.SetOwnership(c_torsion_2d, False)
    ROOT.SetOwnership(c_torsion_smooth, False)

    h2_residuals = ROOT.TH2D("h2_residuals", f"Torsion Residuals (Data - {chosen_model}); x [mm]; y [mm]; #Delta#theta_{0} [#murad]", nx_slices, x_min, x_max, ny_slices, y_min, y_max)
    h2_pulls = ROOT.TH2D("h2_pulls", f"Torsion Pulls (Data - {chosen_model})/#sigma; x [mm]; y [mm]; Pull", nx_slices, x_min, x_max, ny_slices, y_min, y_max)
    for ix in range(1, nx_slices + 1):
        for iy in range(1, ny_slices + 1):
            data_val = h2_torsion_map.GetBinContent(ix, iy)
            if data_val != 0: # Evita di sottrarre in bin vuoti
                x_center = h2_torsion_map.GetXaxis().GetBinCenter(ix)
                y_center = h2_torsion_map.GetYaxis().GetBinCenter(iy)
                fit_val = torsion_fit_2d.Eval(x_center, y_center)
                residual = data_val - fit_val
                h2_residuals.SetBinContent(ix, iy, residual)
                data_err = h2_torsion_map.GetBinError(ix, iy)
                if data_err > 0:
                    h2_pulls.SetBinContent(ix, iy, residual / data_err)


    c_eff_map = ROOT.TCanvas("c_eff_map", "2D Efficiency Map", 900, 700)
    c_eff_map.SetRightMargin(0.15)
    h2_eff_map.SetStats(0)
    h2_eff_map.Draw("COLZ")
    c_eff_map.Update()

    c_residuals = ROOT.TCanvas("c_residuals", "2D Residuals Map", 900, 700)
    c_residuals.SetRightMargin(0.15)
    h2_residuals.SetStats(0)
    h2_residuals.Draw("COLZ")
    c_residuals.Update()

    c_pulls = ROOT.TCanvas("c_pulls", "2D Pulls Map", 900, 700)
    c_pulls.SetRightMargin(0.15)
    h2_pulls.SetStats(0)
    h2_pulls.Draw("COLZ")
    c_pulls.Update()
    
    # Ownership dei nuovi grafici
    ROOT.SetOwnership(h2_eff_map, False)
    ROOT.SetOwnership(c_eff_map, False)
    ROOT.SetOwnership(h2_residuals, False)
    ROOT.SetOwnership(c_residuals, False)
    ROOT.SetOwnership(h2_pulls, False)
    ROOT.SetOwnership(c_pulls, False)

    return fit_params, fit_errors, h2_torsion_map, h2_eff_map, rdf_surface_expr

def channeling_efficiency(df, parameters, best_theta_0):
    N_tot = df.Count().GetValue()

    h_defl_cut = df.Histo1D(("h_defl_cut", "Angular Deflection cut at #pm #theta_{L}/2; #Delta#theta_{x} [#murad]; No. particles", 5000, -2000, parameters["max_value"]), "Deltatheta_x")
    h_cut_value = h_defl_cut.GetValue().Clone("h_cut_value_cloned")

    # GAUSSIAN FIT
    eff_ch = 0.0
    eff_err = 0.0
    N_ch = 0.0
    fit_mean = 0.0
    fit_sigma = 0.0
    
    if N_tot > 0:
        # fitting only right side of the peak (cleanest one)
        fit_min = parameters["deflection_peak"] - 15
        fit_max = parameters["deflection_peak"] + 100
        # fit_max = parameters["max_value"]
        gaus_fit = ROOT.TF1("gaus_fit", "gaus", fit_min, fit_max)
        gaus_fit.SetLineColor(ROOT.kRed)
        
        h_cut_value.Fit(gaus_fit, "RQ0")

        fit_mean = gaus_fit.GetParameter(1)
        fit_mean_error = gaus_fit.GetParError(1)
        fit_sigma = gaus_fit.GetParameter(2)
        fit_sigma_error = gaus_fit.GetParError(2)

        fit_parameters = [fit_mean, fit_mean_error, fit_sigma, fit_sigma_error]
        
        # COMPUTE N CH FROM GAUSSIAN INTEGRAL
        bin_width = h_cut_value.GetBinWidth(1)
        # method 1: integral of the gaussian fit
        # N_ch = gaus_fit.Integral(5950, parameters["max_value"]) / bin_width
        # method 2: integral of the histogram
        bin_min = h_cut_value.FindBin(fit_mean - 3.0 * fit_sigma)
        # bin_max = h_cut_value.FindBin(fit_mean + 4.0 * fit_sigma)
        bin_max = h_cut_value.FindBin(parameters["max_value"])
        N_ch = h_cut_value.Integral(bin_min, bin_max)
        
        # EFFICIENCY WITH ITS ERROR (binomial distribution)
        eff_ch = (N_ch / N_tot) * 100.0
        eff_err = math.sqrt(eff_ch/100.0 * (1.0 - eff_ch/100.0) / N_tot) * 100.0


    c5 = ROOT.TCanvas("c5", "Channeling Efficiency Fit", 1400, 900)
    h_cut_value.SetFillColorAlpha(ROOT.kOrange, 0.6)
    h_cut_value.SetLineColor(ROOT.kOrange)
    h_cut_value.Draw("HIST")
    c5.SetLogy()
    c5.Update()


    ROOT.SetOwnership(c5, False)
    ROOT.SetOwnership(h_cut_value, False)

    c6 = ROOT.TCanvas("c6", "Channeling Efficiency Fit", 1400, 900)
    h_cut_value.SetFillColorAlpha(ROOT.kOrange, 0.6)
    h_cut_value.SetLineColor(ROOT.kOrange)
    h_cut_value.GetXaxis().SetRangeUser(5900, 6150)
    # h_cut_value.GetXaxis().SetRangeUser(6000, 6200)
    h_cut_value.Draw("HIST")
    if N_tot > 0:
        gaus_fit.Draw("SAME")
    legend = ROOT.TPaveText(0.70, 0.75, 0.85, 0.85, "NDC")
    legend.SetBorderSize(1)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetTextAlign(12)
    legend.AddText(f"#epsilon_{{ch}} = {eff_ch:.1f} #pm {eff_err:.1f} %")
    legend.AddText(f"Fit mean: {fit_mean:.1f} #murad")
    legend.Draw()
    c6.Update()

    ROOT.SetOwnership(c6, False)
    ROOT.SetOwnership(h_cut_value, False)
    ROOT.SetOwnership(gaus_fit, False)
    ROOT.SetOwnership(legend, False)

    return eff_ch, eff_err, fit_parameters

def plot_deltatheta_theta(df, parameters, best_theta_0):
    h_deltatheta = df.Histo1D(("h_deltatheta", "Angular Deflection; #Delta#theta_{x} [#murad]; No. particles", 5000, -2000, parameters["max_value"]), "Deltatheta_x")
    h_thetax = df.Histo1D(("h_deltatheta", "Incident Angle; Incident angle #theta_{In, x} [#murad]; No. particles", 500, -150, 150), "thetaIn_x")
    h_deflection = df.Histo2D(("h_deflection", "Angular deflection of the particles as a function of the incident angle; Incident angle #theta_{In, x} [#murad]; Deflection #Delta#theta_{x} [#murad]", 500, -150, 150, 500, -2000, parameters["max_value"]), "thetaIn_x", "Deltatheta_x")
    c_deflection = ROOT.TCanvas("c_deflection","Angular Deflection vs Incident Angle", 800, 600)

    # line1 = ROOT.TLine(best_theta_0 - parameters["theta_L"]/2.0, -2000, best_theta_0 - parameters["theta_L"]/2.0, parameters["max_value"])
    # line2 = ROOT.TLine(best_theta_0 + parameters["theta_L"]/2.0, -2000, best_theta_0 + parameters["theta_L"]/2.0, parameters["max_value"])
    # line1.SetLineColor(ROOT.kRed); line1.SetLineStyle(2); line1.SetLineWidth(2); line1.Draw("SAME")
    # line2.SetLineColor(ROOT.kRed); line2.SetLineStyle(2); line2.SetLineWidth(2); line2.Draw("SAME")

    h_deflection.Draw("COLZ")
    c_deflection.Update()

    ROOT.SetOwnership(c_deflection, False)
    ROOT.SetOwnership(h_deflection, False)

def plot_mean_impact_angle(df, y_min, y_max):
    # TProfile (mean impact angle theta in fuction of impact position y) -> should be linear
    h_prof = df.Profile1D(("h_prof", "Mean Impact Angle vs Y; Impact position d0_y [mm]; Mean Incident Angle #theta_{in, x} [#murad]", 50, y_min, y_max), "Tracks.d0_y", "thetaIn_x")
    
    c_prof = ROOT.TCanvas("c_prof", "Mean Impact Angle", 800, 600)
    h_prof.SetLineColor(ROOT.kBlue)
    h_prof.SetMarkerStyle(20)
    h_prof.SetMarkerColor(ROOT.kBlue)
    h_prof.Draw("PE")
    
    c_prof.Update()
    
    ROOT.SetOwnership(c_prof, False)
    ROOT.SetOwnership(h_prof, False)


def plot_global_efficiency_curve(df, parameters, fit_params, rdf_surface_expr):

    corr_expr = f"thetaIn_x - ({rdf_surface_expr})"
    
    df_corr = df.Define("thetaIn_x_corr", corr_expr)
    df_cut, _ = preliminary_cut_on_deltatheta(df_corr, parameters["deflection_peak"])

    h_all = df_corr.Histo1D(("h_all_corr", "", 1000, -150, 150), "thetaIn_x_corr").GetValue()
    h_chan = df_cut.Histo1D(("h_chan_corr", "", 1000, -150, 150), "thetaIn_x_corr").GetValue()

    scan_min = - 80
    scan_max = + 80
    step = 2.0

    theta_vals, eff_vals, err_vals = [], [], []

    current_theta = scan_min
    while current_theta <= scan_max:
        bin_min = h_all.FindBin(current_theta - parameters["theta_L"] / 2.0)
        bin_max = h_all.FindBin(current_theta + parameters["theta_L"] / 2.0)

        n_tot = h_all.Integral(bin_min, bin_max)
        n_ch  = h_chan.Integral(bin_min, bin_max)

        if n_tot > 50: # Evitiamo rumore statistico sulle code del fascio
            eff = n_ch / n_tot * 100.0
            err = math.sqrt(eff/100.0 * (1.0 - eff/100.0) / n_tot) * 100.0
            
            theta_vals.append(current_theta)
            eff_vals.append(eff)
            err_vals.append(err)

        current_theta += step

    n_points = len(theta_vals)
    gr_eff = ROOT.TGraphErrors(n_points, array('d', theta_vals), array('d', eff_vals), array('d', [0]*n_points), array('d', err_vals))
    
    gr_eff.SetTitle("Global Angular Acceptance (Torsion Corrected); Torsion-Corrected Incoming Angle [#murad]; Channeling Efficiency [%]")
    gr_eff.SetMarkerStyle(20); gr_eff.SetMarkerColor(ROOT.kBlue+1); gr_eff.SetLineColor(ROOT.kBlue+1)

    c_eff_curve = ROOT.TCanvas("c_eff_curve", "Global Efficiency Curve", 800, 600)
    gr_eff.GetXaxis().SetRangeUser(-80,80)
    gr_eff.Draw("AP")

    gaus_eff = ROOT.TF1("gaus_eff", "gaus", -20, +20)
    gaus_eff.SetLineColor(ROOT.kRed)
    gr_eff.Fit(gaus_eff, "RQ0")
    gaus_eff.Draw("SAME")

    max_eff = gaus_eff.GetParameter(0)
    max_eff_err = gaus_eff.GetParError(0)
    mean = gaus_eff.GetParameter(1)
    mean_err = gaus_eff.GetParError(1)
    sigma = gaus_eff.GetParameter(2)
    compatibility = abs(mean)/pow(sigma,0.5)
    
    leg = ROOT.TPaveText(0.53, 0.78, 0.88, 0.88, "NDC")
    leg.SetFillColor(ROOT.kWhite); leg.SetBorderSize(1)
    leg.AddText(f"Max Global Efficiency = ({max_eff:.1f} +/- {max_eff_err:.1f}) %")
    leg.AddText(f"Mean = {mean:.3f} +/- {mean_err:.3f} (comp = {compatibility:.2f})")
    leg.Draw("SAME")

    c_eff_curve.Update()

    ROOT.SetOwnership(gr_eff, False)
    ROOT.SetOwnership(c_eff_curve, False)
    ROOT.SetOwnership(leg, False)

    return max_eff, max_eff_err


def scan_y_margins(df_phys, parameters, x_min, x_max, y_min, y_max, h2_torsion_map, rdf_surface_expr):
    print("\n" + "="*50)
    print("Running Sliding Window Scan for Tau_y and Efficiency...")
    
    window_width = 2 
    step = 0.25
    
    y_centers = []
    tau_ys = []
    tau_y_errs = []
    efficiencies = []
    eff_errs = []
    
    current_y_min = y_min
    
    while current_y_min + window_width <= y_max:
        current_y_max = current_y_min + window_width
        y_center = (current_y_min + current_y_max) / 2.0
        
        fit_func = ROOT.TF2(f"fit_{y_center:.2f}", "[0] + [1]*x + [2]*y", x_min, x_max, current_y_min, current_y_max)
        h2_torsion_map.Fit(fit_func, "RQ0")
        local_theta_0 = fit_func.GetParameter(0)
        local_tau_x = fit_func.GetParameter(1)
        local_tau_y = fit_func.GetParameter(2)
        local_tau_y_err = fit_func.GetParError(2)

        local_fit_params = (local_theta_0, local_tau_x, local_tau_y, 0.0, 0.0, 0.0)
        
        # 2. Estraiamo l'efficienza LOCALE
        df_window = filter2_spatial_cut(df_phys, x_min, x_max, current_y_min, current_y_max)
        # Applichiamo il Lindhard cut dinamico usando i parametri appena trovati
        df_chan = filter3_Lindhard_cut(df_window, parameters, local_fit_params, rdf_surface_expr)
        # Calcoliamo l'efficienza locale riutilizzando la funzione
        eff, err, _ = channeling_efficiency(df_chan, parameters, local_theta_0)
        
        if eff > 0:
            y_centers.append(y_center)
            tau_ys.append(local_tau_y)
            tau_y_errs.append(local_tau_y_err)
            efficiencies.append(eff)
            eff_errs.append(err)
            
        current_y_min += step

    n_pts = len(y_centers)
    arr_y_centers = array('d', y_centers)
    
    # Grafico Tau_y vs Posizione
    gr_tau = ROOT.TGraphErrors(n_pts, arr_y_centers, array('d', tau_ys), array('d', [0]*n_pts), array('d', tau_y_errs))
    gr_tau.SetTitle("Local Torsion vs Y Cut Position; Center of Y-Cut [mm]; Local #tau_{y} [#murad/mm]")
    gr_tau.SetMarkerStyle(20); gr_tau.SetMarkerColor(ROOT.kRed)
    
    c_scan_tau = ROOT.TCanvas("c_scan_tau", "Torsion Scan", 800, 600)
    gr_tau.Draw("APL") # A=Axis, P=Points, L=Line
    c_scan_tau.Update()
    
    # Grafico Efficienza vs Posizione
    gr_eff = ROOT.TGraphErrors(n_pts, arr_y_centers, array('d', efficiencies), array('d', [0]*n_pts), array('d', eff_errs))
    gr_eff.SetTitle("Local Efficiency vs Y Cut Position; Center of Y-Cut [mm]; Local Channeling Efficiency [%]")
    gr_eff.SetMarkerStyle(20); gr_eff.SetMarkerColor(ROOT.kBlue+1)
    
    c_scan_eff = ROOT.TCanvas("c_scan_eff", "Efficiency Scan", 800, 600)
    gr_eff.Draw("APL")
    c_scan_eff.Update()

    ROOT.SetOwnership(c_scan_tau, False); ROOT.SetOwnership(gr_tau, False)
    ROOT.SetOwnership(c_scan_eff, False); ROOT.SetOwnership(gr_eff, False)

def filter1_initial(df):
    df_filtered = df.Filter("SingleTrack == 1")
    return df_filtered

def filter2_spatial_cut(df, x_min, x_max, y_min, y_max, x_cut_margin=0, y_cut_margin=0, restricted=True):
    x_min_new = x_min + x_cut_margin if restricted else x_min
    x_max_new = x_max - x_cut_margin if restricted else x_max
    y_min_new = y_min + y_cut_margin if restricted else y_min
    y_max_new = y_max - y_cut_margin if restricted else y_max
    spatial_cut = f"Tracks.d0_x > {x_min_new} && Tracks.d0_x < {x_max_new} && Tracks.d0_y > {y_min_new} && Tracks.d0_y < {y_max_new}"
    df_filtered = df.Filter(spatial_cut, "Spatial Cut (Crystal Area)")
    return df_filtered

def filter3_Lindhard_cut(df, parameters, fit_params, rdf_surface_expr):
    Lindhard_cut = f"abs(thetaIn_x - {rdf_surface_expr}) <= {parameters['theta_L'] / 2.0}"
    df_filtered = df.Filter(Lindhard_cut, f"Torsion-corrected Lindhard cut")
    return df_filtered

# ==========================================================================================

def main():

    # file = input("File number: ")
    file = 8650
    files = ["recoDataSimple_8430_xtalMerging.root", "recoDataSimple_8431_xtalMerging.root"]
    parameters = get_run_parameters(file)
    filename = "recoDataSimple_" + str(file) + "_xtalMerging.root"

    df = ROOT.RDataFrame("simpleEvent", filename)
    # df = ROOT.RDataFrame("simpleEvent", files)
    print(f"Analyzing {filename} ...")
    print("="*50)

    # CONVERT angles from rad to urad and define Deltatheta_x
    df = df_convert_to_urad_define_deltatheta(df)
    count_0 = df.Count()

    # FILTER 1: select single tracks
    df_phys = filter1_initial(df)
    count_1 = df_phys.Count()

    # FILTER 2: spatial cut (d0_x and d0_y within the crystal area)
    x_cut_margin = 0.2
    y_cut_margin = 0.8
    x_min, x_max, y_min, y_max = compute_spatial_cut_bounds(df_phys, parameters)
    df_phys = filter2_spatial_cut(df_phys, x_min, x_max, y_min, y_max, x_cut_margin, y_cut_margin, restricted=False)
    count_2 = df_phys.Count()

    # plot_mean_impact_angle(df_phys, y_min, y_max)

    # FILTER 3: 2D torsion mapping and dynamic Lindhard cut
    fit_params, fit_errors, h2_torsion_map, h2_eff_map, rdf_surface_expr = torsion_map(df_phys, parameters, x_min, x_max, y_min, y_max, 
                                                                                       x_cut_margin, y_cut_margin, nx_slices=10, ny_slices=40, 
                                                                                       restricted=True, linear=False, chosen_model="parabolic_y")

    theta_0_baseline, tau_x, tau_y, c_yy = fit_params
    theta_0_baseline_err, tau_x_err, tau_y_err, _ = fit_errors 


    # not needed anymore because the torsion fit is not linear but parabolic
    # scan_y_margins(df_phys, parameters, x_min, x_max, y_min, y_max, h2_torsion_map, rdf_surface_expr)

    max_eff_global, max_eff_global_err = plot_global_efficiency_curve(df_phys, parameters, fit_params, rdf_surface_expr)

    df_phys = filter3_Lindhard_cut(df_phys, parameters, fit_params, rdf_surface_expr)
    count_3 = df_phys.Count()

    # CHANNELING EFFICIENCY
    eff_ch, eff_err_stat, fit_efficiency = channeling_efficiency(df_phys, parameters, best_theta_0=fit_params[0])

    surface_expr_up = rdf_surface_expr.replace(f"({theta_0_baseline}", f"({theta_0_baseline + theta_0_baseline_err}", 1)
    surface_expr_down = rdf_surface_expr.replace(f"({theta_0_baseline}", f"({theta_0_baseline - theta_0_baseline_err}", 1)

    df_up = filter3_Lindhard_cut(df_phys, parameters, fit_params, surface_expr_up)
    df_down = filter3_Lindhard_cut(df_phys, parameters, fit_params, surface_expr_down)
    eff_up, _, _ = channeling_efficiency(df_up, parameters, best_theta_0=theta_0_baseline + theta_0_baseline_err)
    eff_down, _, _ = channeling_efficiency(df_down, parameters, best_theta_0=theta_0_baseline - theta_0_baseline_err)
 
    eff_err_syst = abs(eff_up - eff_down) / 2.0
    eff_err_total = math.sqrt(eff_err_stat**2 + eff_err_syst**2)

    # print("="*50)
    # print(filter_message(1, count_0.GetValue(), count_1.GetValue()))
    # print(filter_message(2, count_1.GetValue(), count_2.GetValue()))
    # print(filter_message(3, count_2.GetValue(), count_3.GetValue()))
    # print("="*50)

    print(filter_message("1+2+3", count_0.GetValue(), count_3.GetValue()))
    print(f"Computed channeling efficiency = ({eff_ch:.1f} +/- {eff_err_stat:.1f} [stat] +/- {eff_err_syst:.1f} [syst, torsion map]) %")
    print(f"Total error = +/- {eff_err_total:.1f} %")
    print(f"Channeling peak = ({fit_efficiency[0]:.1f} +/- {fit_efficiency[1]:.1f}) urad , sigma = ({fit_efficiency[2]:.1f} +/- {fit_efficiency[3]:.1f}) urad")
    print(f"Torsion tau_x = {tau_x:.2f} +/- {tau_x_err:.2f} urad/mm")
    print(f"Torsion tau_y = {tau_y:.2f} +/- {tau_y_err:.2f} urad/mm")


if __name__ == "__main__":
    main()




