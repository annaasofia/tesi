import ROOT
import sys
import math
from array import array

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

def preliminary_cut_on_deltatheta(df, deflection_peak, sigma_mult=3.0, tag="nominal"):
    h_defl_before = df.Histo1D((f"h_defl_before_{tag}", "Angular Deflection; #Delta#theta_{x} [#murad]; Counts", 500, -2000, deflection_peak + 1000), "Deltatheta_x") 
    h_defl_val = h_defl_before.GetValue()

    pre_fit = ROOT.TF1(f"pre_fit_{tag}", "gaus", deflection_peak - 500, deflection_peak + 500)
    h_defl_val.Fit(pre_fit, "RQ0")
    pre_mean = pre_fit.GetParameter(1); pre_sigma = pre_fit.GetParameter(2)
    preliminary_cut = pre_mean - sigma_mult * pre_sigma
    df_cut = df.Filter(f"Deltatheta_x > {preliminary_cut}", "Preliminary cut on Deltatheta_x to select channeled particles")
    return df_cut, preliminary_cut

def compute_spatial_cut_bounds(df_phys, parameters, sigma_mult=3.0, tag="nominal"):
    h_d0_xy = df_phys.Histo2D((f"h_d0_xy_{tag}", "Incoming d0_x vs d0_y of all particles; d0_x [mm]; d0_y [mm]", 1700, -10, 20, 1700, -20, 20), "Tracks.d0_x", "Tracks.d0_y")
    h_d0_Out_xy = df_phys.Histo2D((f"h_d0_Out_xy_{tag}", "Outgoing d0_x vs d0_y of all particles; d0_x [mm]; d0_y [mm]", 1700, -10, 20, 1700, -20, 20), "Tracks.d0Out_x", "Tracks.d0Out_y")
    
    df_cut, _ = preliminary_cut_on_deltatheta(df_phys, parameters["deflection_peak"], sigma_mult=sigma_mult, tag=tag)

    h_d0_xy_ch = df_cut.Histo2D((f"h_d0_xy_ch_{tag}", "Incoming beam - channeled particles; d0_x [mm]; d0_y [mm]", 5000, -2, 3, 1700, -8, 9), "Tracks.d0_x", "Tracks.d0_y")
    h_d0_Out_xy_ch = df_cut.Histo2D((f"h_d0_Out_xy_ch_{tag}", "Outgoing beam - channeled particles; d0Out_x [mm]; d0Out_y [mm]", 5000, -2, 3, 1700, -8, 9), "Tracks.d0Out_x", "Tracks.d0Out_y")

    h_d0_x_ch = df_cut.Histo1D((f"h_d0_x_ch_{tag}", "d0_x of channeled particles; d0_x [mm]; Counts", 5000, -2, 3), "Tracks.d0_x")
    h_d0_y_ch = df_cut.Histo1D((f"h_d0_y_ch_{tag}", "d0_y of channeled particles; d0_y [mm]; Counts", 1700, -8, 9), "Tracks.d0_y")
    h_d0_Out_x_ch = df_cut.Histo1D((f"h_d0_Out_x_ch_{tag}", "d0Out_x of channeled particles; d0Out_x [mm]; Counts", 5000, -2, 3), "Tracks.d0Out_x")
    h_d0_Out_y_ch = df_cut.Histo1D((f"h_d0_Out_y_ch_{tag}", "d0Out_y of channeled particles; d0Out_y [mm]; Counts", 1700, -8, 9), "Tracks.d0Out_y")

    f_gaus_in = ROOT.TF1(f"f_gaus_in_{tag}", "gaus", -8, 9)
    h_d0_y_ch.Fit(f_gaus_in, "RQ0")
    y_c_in, y_c_in_err = f_gaus_in.GetParameter(1), f_gaus_in.GetParError(1)
    f_gaus_out = ROOT.TF1(f"f_gaus_out_{tag}", "gaus", -8, 9)
    h_d0_Out_y_ch.Fit(f_gaus_out, "RQ0")
    y_c_out, y_c_out_err = f_gaus_out.GetParameter(1), f_gaus_out.GetParError(1)

    h_y_in = h_d0_y_ch.GetValue()
    max_particles = -1
    best_y_min_in = 0
    for i in range(1, h_y_in.GetNbinsX() + 1):
        current_y_start = h_y_in.GetBinLowEdge(i)
        current_y_end = current_y_start + parameters["width"]
        bin_start = i
        bin_end = h_y_in.FindBin(current_y_end)
        particles_in_window = h_y_in.Integral(bin_start, bin_end)
        if particles_in_window > max_particles:
            max_particles = particles_in_window
            best_y_min_in = current_y_start

    y_min_in = best_y_min_in
    y_max_in = best_y_min_in + parameters["width"]

    h_y_out = h_d0_Out_y_ch.GetValue()
    max_particles = -1
    best_y_min_out = 0
    for i in range(1, h_y_out.GetNbinsX() + 1):
        current_y_start = h_y_out.GetBinLowEdge(i)
        current_y_end = current_y_start + parameters["width"]
        bin_start = i
        bin_end = h_y_out.FindBin(current_y_end)
        particles_in_window = h_y_out.Integral(bin_start, bin_end)
        if particles_in_window > max_particles:
            max_particles = particles_in_window
            best_y_min_out = current_y_start

    y_min_out = best_y_min_out
    y_max_out = best_y_min_out + parameters["width"]

    y_min, y_max = (y_min_in + y_min_out)/2.0, (y_max_in + y_max_out)/2.0

    h_xy = h_d0_Out_xy_ch.GetValue() if hasattr(h_d0_Out_xy_ch, 'GetValue') else h_d0_Out_xy_ch
    y_bin_min = h_xy.GetYaxis().FindBin(y_min)
    y_bin_max = h_xy.GetYaxis().FindBin(y_max)
    h_x_restricted = h_xy.ProjectionX(f"h_x_restricted_{tag}", y_bin_min, y_bin_max)

    max_particles = -1
    best_x_min = 0
    for i in range(1, h_x_restricted.GetNbinsX() + 1):
        current_x_start = h_x_restricted.GetBinLowEdge(i)
        current_x_end = current_x_start + parameters["height"]
        bin_start = i
        bin_end = h_x_restricted.FindBin(current_x_end)
        particles_in_window = h_x_restricted.Integral(bin_start, bin_end)
        if particles_in_window > max_particles:
            max_particles = particles_in_window
            best_x_min = current_x_start

    best_x_max = best_x_min + parameters["height"]

    x_min, x_max = best_x_min - parameters["delta_x"], best_x_max - parameters['delta_x']
    print(f"Final cut in x (after shift of {parameters['delta_x']:.3f} mm): x = [{x_min:.4f} mm, {x_max:.4f} mm], y = [{y_min:.4f} mm, {y_max:.4f} mm]")

    return x_min, x_max, y_min, y_max

def compute_torsion_grid_bounds(df_phys, parameters, sigma_mult=2.5, prelim_sigma_mult=3.0, tag="nominal"):
    df_cut, _ = preliminary_cut_on_deltatheta(df_phys, parameters["deflection_peak"], sigma_mult=prelim_sigma_mult, tag=tag)

    h_d0_y_ch = df_cut.Histo1D((f"h_d0_y_ch_grid_{tag}", "d0_y of channeled particles", 1700, -8, 9), "Tracks.d0_y").GetValue()

    f_gaus = ROOT.TF1(f"f_gaus_grid_{tag}", "gaus", -8, 9)
    h_d0_y_ch.Fit(f_gaus, "RQ0")
    mu, sigma = f_gaus.GetParameter(1), f_gaus.GetParameter(2)

    y_min_fit = mu - sigma_mult * sigma
    y_max_fit = mu + sigma_mult * sigma

    probs = array('d', [0.005, 0.995])
    quant = array('d', [0., 0.])
    h_d0_y_ch.GetQuantiles(2, quant, probs)
    y_min_q, y_max_q = quant

    return y_min_fit, y_max_fit

def torsion_map(df, parameters, x_min, x_max, y_min, y_max, x_cut_margin, y_cut_margin, nx_slices, ny_slices, restricted=False, linear=False, chosen_model="full_quadratic", sigma_mult=3.0, tag="nominal"):

    _, preliminary_cut = preliminary_cut_on_deltatheta(df, parameters["deflection_peak"], sigma_mult=sigma_mult, tag=tag)
    h3_all = df.Histo3D((f"h3_all_{tag}", "", 1000, -100, 100, nx_slices, x_min, x_max, ny_slices, y_min, y_max), "thetaIn_x", "Tracks.d0_x", "Tracks.d0_y").GetValue()
    h3_chan = df.Filter(f"Deltatheta_x > {preliminary_cut}").Histo3D((f"h3_chan_{tag}", "", 1000, -100, 100, nx_slices, x_min, x_max, ny_slices, y_min, y_max), "thetaIn_x", "Tracks.d0_x", "Tracks.d0_y").GetValue()

    h2_torsion_map = ROOT.TH2D(f"h2_torsion_map_{tag}", "2D Torsion Map (bins 0.1 x 0.1 mm); x at crystal surface [mm]; y at crystal surface [mm]; angle shift #theta_{0} [#murad]", nx_slices, x_min, x_max, ny_slices, y_min, y_max)
    h2_eff_map = ROOT.TH2D(f"h2_eff_map_{tag}", "2D Efficiency Map; x [mm]; y [mm]; Local Efficiency [%]", nx_slices, x_min, x_max, ny_slices, y_min, y_max)

    h2_torsion_map.Sumw2()
    h2_eff_map.Sumw2()

    scan_min, scan_max, step = -80, 80, 1.0

    for ix in range(1, nx_slices + 1):
        for iy in range(1, ny_slices + 1):
            
            h1_all = h3_all.ProjectionX(f"all_proj_{tag}_{ix}_{iy}", ix, ix, iy, iy)
            h1_chan = h3_chan.ProjectionX(f"chan_proj_{tag}_{ix}_{iy}", ix, ix, iy, iy)
            
            best_theta_0_raw = 0.0
            max_eff_raw = -1.0
            
            theta_vals, eff_vals = [], []
            current_theta = scan_min
            
            while current_theta <= scan_max:
                bin_min = h1_all.FindBin(current_theta - parameters["theta_L"] / 2.0)
                bin_max = h1_all.FindBin(current_theta + parameters["theta_L"] / 2.0)
                
                n_tot_test = h1_all.Integral(bin_min, bin_max)
                n_ch_test  = h1_chan.Integral(bin_min, bin_max)
                
                eff_test = 0.0
                if n_tot_test > 200:
                    eff_test = n_ch_test / n_tot_test * 100.0
                    if eff_test > max_eff_raw:
                        max_eff_raw = eff_test
                        best_theta_0_raw = current_theta
                        
                theta_vals.append(current_theta)
                eff_vals.append(eff_test)
                current_theta += step
                
            if max_eff_raw > 0:
                gr_eff_slice = ROOT.TGraph(len(theta_vals), array('d', theta_vals), array('d', eff_vals))
                gaus_eff_slice = ROOT.TF1(f"gaus_eff_{tag}_{ix}_{iy}", "gaus", best_theta_0_raw - 5.0, best_theta_0_raw + 20.0)
                gaus_eff_slice.SetParameters(max_eff_raw, best_theta_0_raw, 5.0)
                gaus_eff_slice.SetParLimits(1, -80.0, 80.0)
                
                gr_eff_slice.Fit(gaus_eff_slice, "RQ0")
                
                local_theta_0 = gaus_eff_slice.GetParameter(1)
                local_theta_0_err = gaus_eff_slice.GetParError(1)
                local_max_eff = gaus_eff_slice.GetParameter(0)
                local_max_eff_err = gaus_eff_slice.GetParError(0)
                
                h2_torsion_map.SetBinContent(ix, iy, local_theta_0)
                h2_torsion_map.SetBinError(ix, iy, local_theta_0_err)
                h2_eff_map.SetBinContent(ix, iy, local_max_eff)
                h2_eff_map.SetBinError(ix, iy, local_max_eff_err)

    # c_torsion_2d = ROOT.TCanvas(f"c_torsion_2d_{tag}", "2D Torsion Map", 900, 700)
    # c_torsion_2d.SetRightMargin(0.15)
    # h2_torsion_map.SetStats(0)
    # h2_torsion_map.Draw("colz")
    # c_torsion_2d.Update()
    
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
    torsion_fit_2d = ROOT.TF2(f"torsion_fit_2d_{tag}", fit_models[chosen_model], current_x_min, current_x_max, current_y_min, current_y_max)
    fit_result = h2_torsion_map.Fit(torsion_fit_2d, "RQ0S")

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

    # torsion_plot_2d = torsion_fit_2d.Clone(f"torsion_plot_2d_{tag}")
    # torsion_plot_2d.SetRange(current_x_min, current_y_min, current_x_max, current_y_max)
    # torsion_plot_2d.SetTitle("Continuous 2D Torsion Map; x at crystal surface [mm]; y at crystal surface [mm]; angle shift #theta_{0} [#murad]")
    # c_torsion_smooth = ROOT.TCanvas(f"c_torsion_smooth_{tag}", "Continuous 2D Torsion Map", 900, 700)
    # c_torsion_smooth.SetRightMargin(0.15)
    # torsion_plot_2d.Draw("surf2") 
    # c_torsion_smooth.Update()

    # ROOT.SetOwnership(h2_torsion_map, False)
    # ROOT.SetOwnership(torsion_plot_2d, False)
    # ROOT.SetOwnership(c_torsion_2d, False)
    # ROOT.SetOwnership(c_torsion_smooth, False)

    h2_residuals = ROOT.TH2D(f"h2_residuals_{tag}", f"Torsion Residuals (Data - {chosen_model}); x [mm]; y [mm]; #Delta#theta_{0} [#murad]", nx_slices, x_min, x_max, ny_slices, y_min, y_max)
    h2_pulls = ROOT.TH2D(f"h2_pulls_{tag}", f"Torsion Pulls (Data - {chosen_model})/#sigma; x [mm]; y [mm]; Pull", nx_slices, x_min, x_max, ny_slices, y_min, y_max)
    for ix in range(1, nx_slices + 1):
        for iy in range(1, ny_slices + 1):
            data_val = h2_torsion_map.GetBinContent(ix, iy)
            if data_val != 0:
                x_center = h2_torsion_map.GetXaxis().GetBinCenter(ix)
                y_center = h2_torsion_map.GetYaxis().GetBinCenter(iy)
                fit_val = torsion_fit_2d.Eval(x_center, y_center)
                residual = data_val - fit_val
                h2_residuals.SetBinContent(ix, iy, residual)
                data_err = h2_torsion_map.GetBinError(ix, iy)
                if data_err > 0:
                    h2_pulls.SetBinContent(ix, iy, residual / data_err)

    # c_eff_map = ROOT.TCanvas(f"c_eff_map_{tag}", "2D Efficiency Map", 900, 700)
    # c_eff_map.SetRightMargin(0.15)
    # h2_eff_map.SetStats(0)
    # h2_eff_map.Draw("COLZ")
    # c_eff_map.Update()

    # c_residuals = ROOT.TCanvas(f"c_residuals_{tag}", "2D Residuals Map", 900, 700)
    # c_residuals.SetRightMargin(0.15)
    # h2_residuals.SetStats(0)
    # h2_residuals.Draw("COLZ")
    # c_residuals.Update()

    # c_pulls = ROOT.TCanvas(f"c_pulls_{tag}", "2D Pulls Map", 900, 700)
    # c_pulls.SetRightMargin(0.15)
    # h2_pulls.SetStats(0)
    # h2_pulls.Draw("COLZ")
    # c_pulls.Update()
    
    # ROOT.SetOwnership(h2_eff_map, False)
    # ROOT.SetOwnership(c_eff_map, False)
    # ROOT.SetOwnership(h2_residuals, False)
    # ROOT.SetOwnership(c_residuals, False)
    # ROOT.SetOwnership(h2_pulls, False)
    # ROOT.SetOwnership(c_pulls, False)

    return fit_params, fit_errors, h2_torsion_map, h2_eff_map, rdf_surface_expr

def channeling_efficiency(df, parameters, best_theta_0, n_sigma_low=3.0, tag="nominal"):
    N_tot = df.Count().GetValue()

    gaus_fit = None
    fit_parameters = [0.0, 0.0, 0.0, 0.0]

    if N_tot == 0:
        print("WARNING: channeling_efficiency called on an empty dataframe, returning zeros.")
        return 0.0, 0.0, fit_parameters, None

    h_defl_cut = df.Histo1D((f"h_defl_cut_{tag}", "Angular Deflection cut at #pm #theta_{L}/2; #Delta#theta_{x} [#murad]; No. particles", 5000, -2000, parameters["max_value"]), "Deltatheta_x")
    h_cut_value = h_defl_cut.GetValue().Clone(f"h_cut_value_cloned_{tag}")

    eff_ch = 0.0
    eff_err = 0.0
    N_ch = 0.0
    fit_mean = 0.0
    fit_sigma = 0.0
    
    if N_tot > 0:
        fit_min = parameters["deflection_peak"] - 20
        fit_max = parameters["deflection_peak"] + 100
        gaus_fit = ROOT.TF1(f"gaus_fit_{tag}", "gaus", fit_min, fit_max)
        gaus_fit.SetLineColor(ROOT.kRed)
        
        h_cut_value.Fit(gaus_fit, "RQ0")

        fit_mean = gaus_fit.GetParameter(1)
        fit_mean_error = gaus_fit.GetParError(1)
        fit_sigma = gaus_fit.GetParameter(2)
        fit_sigma_error = gaus_fit.GetParError(2)

        fit_parameters = [fit_mean, fit_mean_error, fit_sigma, fit_sigma_error]
        
        bin_min = h_cut_value.FindBin(fit_mean - n_sigma_low * fit_sigma)
        bin_max = h_cut_value.FindBin(parameters["max_value"])
        N_ch = h_cut_value.Integral(bin_min, bin_max)
        
        eff_ch = (N_ch / N_tot) * 100.0
        eff_err = math.sqrt(eff_ch/100.0 * (1.0 - eff_ch/100.0) / N_tot) * 100.0

    # c5 = ROOT.TCanvas(f"c5_{tag}", "Channeling Efficiency Fit", 1400, 900)
    # h_cut_value.SetFillColorAlpha(ROOT.kOrange, 0.6)
    # h_cut_value.SetLineColor(ROOT.kOrange)
    # h_cut_value.Draw("HIST")
    # c5.SetLogy()
    # c5.Update()

    # ROOT.SetOwnership(c5, False)
    # ROOT.SetOwnership(h_cut_value, False)

    # c6 = ROOT.TCanvas(f"c6_{tag}", "Channeling Efficiency Fit", 1400, 900)
    # h_cut_value.SetFillColorAlpha(ROOT.kOrange, 0.6)
    # h_cut_value.SetLineColor(ROOT.kOrange)
    # h_cut_value.GetXaxis().SetRangeUser(5900, 6150)
    # h_cut_value.Draw("HIST")
    # if N_tot > 0:
    #     gaus_fit.Draw("SAME")
    # legend = ROOT.TPaveText(0.70, 0.75, 0.85, 0.85, "NDC")
    # legend.SetBorderSize(1)
    # legend.SetFillColor(ROOT.kWhite)
    # legend.SetTextAlign(12)
    # legend.AddText(f"#epsilon_{{ch}} = {eff_ch:.1f} #pm {eff_err:.1f} %")
    # legend.AddText(f"Fit mean: {fit_mean:.1f} #murad")
    # legend.Draw()
    # c6.Update()

    # ROOT.SetOwnership(c6, False)
    # ROOT.SetOwnership(h_cut_value, False)
    # if gaus_fit is not None:
    #     ROOT.SetOwnership(gaus_fit, False)
    # ROOT.SetOwnership(legend, False)

    return eff_ch, eff_err, fit_parameters, h_cut_value

def plot_global_efficiency_curve(df, parameters, fit_params, rdf_surface_expr, tag="global_curve"):

    corr_expr = f"thetaIn_x - ({rdf_surface_expr})"
    
    df_corr = df.Define("thetaIn_x_corr", corr_expr)
    df_cut, _ = preliminary_cut_on_deltatheta(df_corr, parameters["deflection_peak"], tag=tag)

    h_all = df_corr.Histo1D((f"h_all_corr_{tag}", "", 1000, -150, 150), "thetaIn_x_corr").GetValue()
    h_chan = df_cut.Histo1D((f"h_chan_corr_{tag}", "", 1000, -150, 150), "thetaIn_x_corr").GetValue()

    scan_min = -80
    scan_max = +80
    step = 2.0

    theta_vals, eff_vals, err_vals = [], [], []

    current_theta = scan_min
    while current_theta <= scan_max:
        bin_min = h_all.FindBin(current_theta - parameters["theta_L"] / 2.0)
        bin_max = h_all.FindBin(current_theta + parameters["theta_L"] / 2.0)

        n_tot = h_all.Integral(bin_min, bin_max)
        n_ch  = h_chan.Integral(bin_min, bin_max)

        if n_tot > 50:
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

    c_eff_curve = ROOT.TCanvas(f"c_eff_curve_{tag}", "Global Efficiency Curve", 800, 600)
    gr_eff.GetXaxis().SetRangeUser(-80,80)
    gr_eff.Draw("AP")

    gaus_eff = ROOT.TF1(f"gaus_eff_{tag}", "gaus", -20, +20)
    gaus_eff.SetLineColor(ROOT.kRed)
    gr_eff.Fit(gaus_eff, "RQ0")
    gaus_eff.Draw("SAME")

    max_eff = gaus_eff.GetParameter(0)
    max_eff_err = gaus_eff.GetParError(0)
    mean = gaus_eff.GetParameter(1)
    mean_err = gaus_eff.GetParError(1)
    sigma = gaus_eff.GetParameter(2)
    compatibility = abs(mean)/mean_err
    
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

def filter3_Lindhard_cut(df, parameters, fit_params, rdf_surface_expr, halfwidth_shift=0.0):
    halfwidth = parameters['theta_L'] / 2.0 + halfwidth_shift
    Lindhard_cut = f"abs(thetaIn_x - {rdf_surface_expr}) <= {halfwidth}"
    df_filtered = df.Filter(Lindhard_cut, f"Torsion-corrected Lindhard cut")
    return df_filtered

def evaluate_channeling_efficiency_variant(df_source, parameters, fit_params, x_bounds, y_bounds, surface_expr, halfwidth_shift=0.0, tag="variant"):
    df = filter2_spatial_cut(df_source, x_bounds[0], x_bounds[1], y_bounds[0], y_bounds[1], restricted=False)
    df = filter3_Lindhard_cut(df, parameters, fit_params, surface_expr, halfwidth_shift=halfwidth_shift)
    return channeling_efficiency(df, parameters, best_theta_0=fit_params[0], tag=tag)

def shift_and_rerun_systematic(name, df_source, parameters, fit_params, nominal_kwargs, up_overrides, down_overrides):
    kwargs_up = {**nominal_kwargs, **up_overrides}
    kwargs_down = {**nominal_kwargs, **down_overrides}

    eff_up, _, _, _ = evaluate_channeling_efficiency_variant(df_source, parameters, fit_params, tag=f"{name}_up".replace(" ", "_"), **kwargs_up)
    eff_down, _, _, _ = evaluate_channeling_efficiency_variant(df_source, parameters, fit_params, tag=f"{name}_down".replace(" ", "_"), **kwargs_down)

    syst = abs(eff_up - eff_down) / 2.0
    print(f"\tsyst ({name}): up={eff_up:.3f}%  down={eff_down:.3f}%  -> +/-{syst:.3f}%")
    return syst, eff_up, eff_down

def compute_nsigma_systematic(df_phys, parameters, fit_params, n_sigma_up=2.5, n_sigma_down=3.5):
    eff_up, _, _, _ = channeling_efficiency(df_phys, parameters, best_theta_0=fit_params[0], n_sigma_low=n_sigma_up, tag="nsigma_up")
    eff_down, _, _, _ = channeling_efficiency(df_phys, parameters, best_theta_0=fit_params[0], n_sigma_low=n_sigma_down, tag="nsigma_down")
    syst = abs(eff_up - eff_down) / 2.0
    print(f"\tsyst (n_sigma window): 2.5sigma={eff_up:.3f}%  3.5sigma={eff_down:.3f}%  -> +/-{syst:.3f}%")
    return syst

def compute_efficiency_systematics(df_singletrack, df_phys, parameters, fit_params, fit_errors, x_min, x_max, y_min, y_max, rdf_surface_expr):
    nominal_kwargs = dict(x_bounds=(x_min, x_max), y_bounds=(y_min, y_max), surface_expr=rdf_surface_expr, halfwidth_shift=0.0)

    systematics = {}

    # --- 1. Torsion baseline (theta_0) fit uncertainty ---
    surface_expr_up = rdf_surface_expr.replace(f"({fit_params[0]}", f"({fit_params[0] + fit_errors[0]}", 1)
    surface_expr_down = rdf_surface_expr.replace(f"({fit_params[0]}", f"({fit_params[0] - fit_errors[0]}", 1)
    systematics["theta0_shift"], _, _ = shift_and_rerun_systematic(
        "theta0 shift", df_singletrack, parameters, fit_params, nominal_kwargs,
        up_overrides=dict(surface_expr=surface_expr_up),
        down_overrides=dict(surface_expr=surface_expr_down))

    # --- 2. Angular (theta_in) resolution on the Lindhard window ---
    mean_theta_res = df_phys.Mean("Tracks.thetaInErr_x").GetValue() * 1e6  # urad
    shift = min(mean_theta_res, parameters['theta_L'] / 4.0)
    print(f"Mean theta resolution = {mean_theta_res:.3f} urad (applied shift = {shift:.3f} urad)")
    systematics["theta_resolution"], _, _ = shift_and_rerun_systematic(
        "theta resolution", df_singletrack, parameters, fit_params, nominal_kwargs,
        up_overrides=dict(halfwidth_shift=+shift),
        down_overrides=dict(halfwidth_shift=-shift))

    # --- 3. Spatial box boundary (d0 resolution), x and y treated separately ---
    mean_d0err_x = df_phys.Mean("Tracks.d0Err_x").GetValue()
    mean_d0err_y = df_phys.Mean("Tracks.d0Err_y").GetValue()

    systematics["spatial_box_x"], _, _ = shift_and_rerun_systematic(
        "spatial box x", df_singletrack, parameters, fit_params, nominal_kwargs,
        up_overrides=dict(x_bounds=(x_min - mean_d0err_x, x_max + mean_d0err_x)),
        down_overrides=dict(x_bounds=(x_min + mean_d0err_x, x_max - mean_d0err_x)))

    systematics["spatial_box_y"], _, _ = shift_and_rerun_systematic(
        "spatial box y", df_singletrack, parameters, fit_params, nominal_kwargs,
        up_overrides=dict(y_bounds=(y_min - mean_d0err_y, y_max + mean_d0err_y)),
        down_overrides=dict(y_bounds=(y_min + mean_d0err_y, y_max - mean_d0err_y)))

    systematics["n_sigma_window"] = compute_nsigma_systematic(df_phys, parameters, fit_params)

    return systematics


def compute_margin_systematic(df_phys, parameters, x_min, x_max, y_min_grid, y_max_grid,
                               nx_slices, ny_slices, eff_nominal,
                               margin_up=(0.35, 0.8), margin_down=(0.1, 0.25)):
    effs = []
    for i, (margin_x, margin_y) in enumerate((margin_up, margin_down)):
        variant_tag = f"margin_{i}"
        fit_params_v, fit_errors_v, _, _, surface_expr_v = torsion_map(
            df_phys, parameters, x_min, x_max, y_min_grid, y_max_grid,
            margin_x, margin_y, nx_slices=nx_slices, ny_slices=ny_slices,
            restricted=True, linear=False, chosen_model="parabolic_y", tag=variant_tag)

        df_v = filter3_Lindhard_cut(df_phys, parameters, fit_params_v, surface_expr_v)
        eff_v, _, _, _ = channeling_efficiency(df_v, parameters, best_theta_0=fit_params_v[0], tag=variant_tag)
        effs.append(eff_v)
        print(f"\tmargins ({margin_x}, {margin_y}) -> eff = {eff_v:.3f}%")

    syst_margin = abs(effs[0] - effs[1]) / 2.0
    print(f"\tsyst (torsion-fit margins): +/-{syst_margin:.3f}%  (nominal = {eff_nominal:.3f}%)")
    return syst_margin


def compute_torsion_variant_systematic(name, df_phys, parameters, x_min, x_max, y_min_grid, y_max_grid,
                                        x_cut_margin, y_cut_margin, variant_a_kwargs, variant_b_kwargs):
    effs = []
    base_tag = name.replace(" ", "_")
    for i, kwargs in enumerate((variant_a_kwargs, variant_b_kwargs)):
        variant_tag = f"{base_tag}_{i}"
        params = dict(nx_slices=10, ny_slices=65, restricted=True, linear=False, chosen_model="parabolic_y")
        params.update(kwargs)

        fit_params_v, fit_errors_v, _, _, surface_expr_v = torsion_map(
            df_phys, parameters, x_min, x_max, y_min_grid, y_max_grid,
            x_cut_margin, y_cut_margin, tag=variant_tag, **params)

        df_v = filter3_Lindhard_cut(df_phys, parameters, fit_params_v, surface_expr_v)
        eff_v, _, _, _ = channeling_efficiency(df_v, parameters, best_theta_0=fit_params_v[0], tag=variant_tag)
        effs.append(eff_v)

    syst = abs(effs[0] - effs[1]) / 2.0
    print(f"\tsyst ({name}): {effs[0]:.3f}% vs {effs[1]:.3f}%  -> +/-{syst:.3f}%")
    return syst


def compute_preliminary_cut_systematic(df_phys, parameters, nx_slices, ny_slices,
                                        x_cut_margin, y_cut_margin,
                                        sigma_mult_up=2.5, sigma_mult_down=3.5):
    effs = []
    for i, sigma_mult in enumerate((sigma_mult_up, sigma_mult_down)):
        variant_tag = f"prelim_{i}"
        x_min_v, x_max_v, y_min_v, y_max_v = compute_spatial_cut_bounds(
            df_phys, parameters, sigma_mult=sigma_mult, tag=variant_tag)

        y_min_grid_v, y_max_grid_v = compute_torsion_grid_bounds(
            df_phys, parameters, sigma_mult=2.0, prelim_sigma_mult=sigma_mult, tag=variant_tag)

        df_v = filter2_spatial_cut(df_phys, x_min_v, x_max_v, y_min_v, y_max_v, restricted=False)

        fit_params_v, fit_errors_v, _, _, surface_expr_v = torsion_map(
            df_v, parameters, x_min_v, x_max_v, y_min_grid_v, y_max_grid_v,
            x_cut_margin, y_cut_margin, nx_slices=nx_slices, ny_slices=ny_slices,
            restricted=True, linear=False, chosen_model="parabolic_y", sigma_mult=sigma_mult, tag=variant_tag)

        df_v = filter3_Lindhard_cut(df_v, parameters, fit_params_v, surface_expr_v)
        eff_v, _, _, _ = channeling_efficiency(df_v, parameters, best_theta_0=fit_params_v[0], tag=variant_tag)
        effs.append(eff_v)
        print(f"\tpreliminary cut sigma_mult={sigma_mult} -> eff = {eff_v:.3f}%")

    syst = abs(effs[0] - effs[1]) / 2.0
    print(f"\tsyst (preliminary cut threshold): +/-{syst:.3f}%")
    return syst

def compute_block_stability_check(df_phys, parameters, fit_params, n_blocks=10):

    # based on event number, split the dataset into n_blocks contiguous blocks, compute channeling efficiency in each block
    # and compare the observed spread (RMS) to the average binomial statistical error
    # if rms_observed >> mean_stat_err, there is time-correlation / overdispersion not captured by the simple binomial formula (e.g. beam drift, slow alignment changes during the run)
    min_evt = df_phys.Min("Event.evtnum").GetValue()
    max_evt = df_phys.Max("Event.evtnum").GetValue()
    block_width = (max_evt - min_evt) / n_blocks

    effs, errs = [], []
    print(f"\n\tBlock stability check ({n_blocks} blocks, evtnum {min_evt:.0f}-{max_evt:.0f}):")
    for i in range(n_blocks):
        lo = min_evt + i * block_width
        hi = min_evt + (i + 1) * block_width
        df_block = df_phys.Filter(f"Event.evtnum >= {lo} && Event.evtnum < {hi}")
        eff, err, _, _ = channeling_efficiency(df_block, parameters, best_theta_0=fit_params[0], tag=f"block_{i}")
        if eff > 0:
            effs.append(eff)
            errs.append(err)
            print(f"\t\tblock {i}: eff = {eff:.3f} +/- {err:.3f} %")

    if len(effs) < 2:
        print("\tWARNING: not enough non-empty blocks for stability check.")
        return 0.0, 0.0

    mean_eff = sum(effs) / len(effs)
    rms_observed = math.sqrt(sum((e - mean_eff)**2 for e in effs) / (len(effs) - 1))
    mean_stat_err = sum(errs) / len(errs)

    ratio = rms_observed / mean_stat_err if mean_stat_err > 0 else float('nan')
    print(f"\tRMS across blocks     = {rms_observed:.3f} %")
    print(f"\tMean binomial error   = {mean_stat_err:.3f} %")
    print(f"\tRatio (RMS/stat_err)  = {ratio:.2f}  ({'consistent with binomial' if ratio < 1.5 else 'possible overdispersion / time correlation'})")

    return rms_observed, mean_stat_err

def compute_fit_uncertainty_systematic(df_phys, parameters, fit_efficiency, tag="fitunc"):

    # error on n ch propagated from the uncertainty on the fit parameters used to define the integration window for n_ch
    fit_mean, fit_mean_err, fit_sigma, fit_sigma_err = fit_efficiency

    N_tot = df_phys.Count().GetValue()
    if N_tot == 0:
        return 0.0

    h_defl_cut = df_phys.Histo1D((f"h_defl_cut_{tag}", "Angular Deflection; #Delta#theta_{x} [#murad]; No. particles", 5000, -2000, parameters["max_value"]), "Deltatheta_x")
    h_cut_value = h_defl_cut.GetValue().Clone(f"h_cut_value_cloned_{tag}")

    def eff_from_bounds(mean, sigma, n_sigma_low=3.0):
        bin_min = h_cut_value.FindBin(mean - n_sigma_low * sigma)
        bin_max = h_cut_value.FindBin(parameters["max_value"])
        N_ch = h_cut_value.Integral(bin_min, bin_max)
        return (N_ch / N_tot) * 100.0

    eff_up = eff_from_bounds(fit_mean + fit_mean_err, fit_sigma + fit_sigma_err)
    eff_down = eff_from_bounds(fit_mean - fit_mean_err, fit_sigma - fit_sigma_err)

    syst = abs(eff_up - eff_down) / 2.0
    print(f"\tstat (fit uncertainty on N_ch window): up={eff_up:.3f}%  down={eff_down:.3f}%  -> +/-{syst:.3f}%")
    return syst
