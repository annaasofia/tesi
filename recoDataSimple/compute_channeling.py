import ROOT
import math
import sys
import os
import numpy as np
from array import array
import matplotlib.pyplot as plt
import channeling_functions as chfun
import plotting_utils as pu

ROOT.ROOT.EnableImplicitMT() 
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kBird)

def main():

    file = 8430
    parameters = chfun.get_run_parameters(file)
    filename = "data/recoDataSimple_" + str(file) + "_xtalMerging.root"

    # global PLOT_DIR
    PLOT_DIR = f"./plots/plots_{file}_py"
    os.makedirs(PLOT_DIR, exist_ok=True)
    chfun.PLOT_DIR = PLOT_DIR

    df = ROOT.RDataFrame("simpleEvent", filename)
    print(f"Analyzing {filename} ...")
    print("=" * 50)

    df = chfun.df_convert_to_urad_define_deltatheta(df)
    count_0 = df.Count()

    # ===================== FILTERS =====================
    # FILTER 1: single tracks
    df_phys = chfun.filter1_initial(df)
    count_1 = df_phys.Count()
    df_singletrack = df_phys  # kept for systematic variants (pre-spatial-cut)

    # FILTER 2: nominal spatial cut
    x_cut_margin = 0.2
    y_cut_margin = 0.5
    x_min, x_max, y_min, y_max = chfun.compute_spatial_cut_bounds(df_phys, parameters, tag="nominal")
    df_phys = chfun.filter2_spatial_cut(df_phys, x_min, x_max, y_min, y_max,
                                   x_cut_margin, y_cut_margin, restricted=False)
    count_2 = df_phys.Count()

    # FILTER 3: torsion map + dynamic Lindhard cut (nominal)
    y_min_grid, y_max_grid = chfun.compute_torsion_grid_bounds(df_phys, parameters, sigma_mult=2.0, tag="nominal")
    nx_slices, ny_slices = 10, 65
    fit_params, fit_errors, h2_torsion_map, h2_eff_map, rdf_surface_expr = chfun.torsion_map(
        df_phys, parameters, x_min, x_max, y_min_grid, y_max_grid,
        x_cut_margin, y_cut_margin, nx_slices=nx_slices, ny_slices=ny_slices,
        restricted=True, linear=False, chosen_model="parabolic_y", tag="nominal")

    _, _ = chfun.plot_global_efficiency_curve(df_phys, parameters, fit_params, rdf_surface_expr, tag="nominal")
    chfun.plot_impact_vs_deltatheta(df_phys, parameters, axis="x", tag="nominal", xlim=(x_min, x_max))
    chfun.plot_impact_vs_deltatheta(df_phys, parameters, axis="y", tag="nominal", xlim=(y_min, y_max))
    chfun.plot_deflection_map(df_phys, parameters, x_min, x_max, y_min, y_max, nx_slices=10, ny_slices=65, tag="nominal")

    df_phys = chfun.filter3_Lindhard_cut(df_phys, parameters, fit_params, rdf_surface_expr)
    count_3 = df_phys.Count()

    # NOMINAL CHANNELING EFFICIENCY
    eff_ch, eff_err_stat_binomial, fit_efficiency, histo_final = chfun.channeling_efficiency(df_phys, parameters, best_theta_0=fit_params[0], tag="nominal")

    # ===================== STATISTICAL ERROR BUDGET =====================
    print("\n" + "=" * 50)
    print("Statistical error breakdown:")

    eff_err_stat_fitunc = chfun.compute_fit_uncertainty_systematic(df_phys, parameters, fit_efficiency)

    eff_err_stat = math.sqrt(eff_err_stat_binomial**2 + eff_err_stat_fitunc**2)
    print(f"\tbinomial              = {eff_err_stat_binomial:.3f} %")
    print(f"\tfit uncertainty       = {eff_err_stat_fitunc:.3f} %")
    print(f"\tcombined statistical  = {eff_err_stat:.3f} %")

    # rms_blocks, mean_stat_err_blocks = chfun.compute_block_stability_check(df_phys, parameters, fit_params, n_blocks=10)

    # ===================== SAVE FINAL HISTO =====================

    # out_filename = f"final_histo_run_{file}.root"
    # out_file = ROOT.TFile(out_filename, "RECREATE")
    # histo_final.Write(f"h_defl_run_{file}")
    # out_file.Close()

    # ===================== SYSTEMATICS =====================
    # print("\n" + "=" * 50)
    # print("Systematic error breakdown:")
    # systematics = chfun.compute_efficiency_systematics(
    #     df_singletrack, df_phys, parameters, fit_params, fit_errors,
    #     x_min, x_max, y_min, y_max, rdf_surface_expr)

    # systematics["fit_model"] = compute_torsion_variant_systematic(
    #     "fit model choice", df_phys, parameters, x_min, x_max, y_min_grid, y_max_grid,
    #     x_cut_margin, y_cut_margin,
    #     variant_a_kwargs=dict(chosen_model="parabolic_y"),
    #     variant_b_kwargs=dict(chosen_model="full_quadratic"))

    # systematics["grid_binning"] = compute_torsion_variant_systematic(
    #     "grid bin choice", df_phys, parameters, x_min, x_max, y_min_grid, y_max_grid,
    #     x_cut_margin, y_cut_margin,
    #     variant_a_kwargs=dict(nx_slices=5, ny_slices=40),
    #     variant_b_kwargs=dict(nx_slices=15, ny_slices=90))

    # systematics["torsion_margins"] = compute_margin_systematic(
    #      df_phys, parameters, x_min, x_max, y_min_grid, y_max_grid,
    #      nx_slices, ny_slices, eff_nominal=eff_ch)

    # systematics["preliminary_cut"] = compute_preliminary_cut_systematic(
    #     df_phys, parameters, nx_slices, ny_slices, x_cut_margin, y_cut_margin)

    eff_err_syst, eff_err_total = 0.0, 0.0
    # eff_err_syst = math.sqrt(sum(v**2 for v in systematics.values()))
    # eff_err_total = math.sqrt(eff_err_stat**2 + eff_err_syst**2)

    # ===================== REPORT =====================
    print("=" * 50)
    print(chfun.filter_message("1+2+3", count_0.GetValue(), count_3.GetValue()))
    print("=" * 50)
    print(f"Computed channeling efficiency = ({eff_ch:.1f} +/- {eff_err_stat:.1f} [stat] +/- {eff_err_syst:.3f} [syst]) %")
    print(f"Total error = +/- {eff_err_total:.1f} %")
    # for name, val in systematics.items():
    #     print(f"\tsyst ({name}) = {val:.3f} %")
    print(f"Channeling peak = ({fit_efficiency[0]:.1f} +/- {fit_efficiency[1]:.1f}) urad, sigma = ({fit_efficiency[2]:.1f} +/- {fit_efficiency[3]:.1f}) urad")
    print(f"Torsion tau_x = {fit_params[1]:.2f} +/- {fit_errors[1]:.2f} urad/mm")
    print(f"Torsion tau_y = {fit_params[2]:.2f} +/- {fit_errors[2]:.2f} urad/mm")


if __name__ == "__main__":
    main()