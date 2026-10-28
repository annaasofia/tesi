1. [Preliminary analysis](#preliminary-analysis)
    - ✅ [Calibrate the detectors and alignment check](#calibrate-the-detectors-and-alignment-check) 
    - ✅ [To do preliminary](#to-do-preliminary)
    - ✅ [Crystal edges found for each run](#crystal-edges-found-for-each-run)
2. [Analysis](#analysis)
    - ✅ [Measure the channeling efficiency $\epsilon_{ch}=\frac{N_{ch}}{N_{tot}}\cdot 100$](#measure-the-channeling-efficiency)
    - [Measure the torsion correction](#measure-the-torsion-correction)
    - [Mapping the torsion](#mapping-the-torsion)
    - [To do about torsion](#to-do-about-torsion)
3. [Error analysis](#error-analysis)
    - [Taking into account d and theta errors](#taking-into-account-d-and-theta-errors)
    - [Computing efficiency errors](#computing-efficiency-errors)
4. [Multiple coulomb scattering](#multiple-coulomb-scattering)
5. [Data results](#data-results)
6. [Simulation](#simulation)

# Preliminary Analysis
## Calibrate the detectors and alignment check
Calibrate the detectors with [check_alignment_detectors.py](../recoDataSimple/preliminary_analysis/check_alignment_detectors.py):  
    1) compute histograms in $x$ and $y$ showing $d0_{out} - d0_{in}$: if the offset is smaller than $\sim 1–2\,\sigma$ of its own fit error, treat it as statistically consistent with zero;  
    2) fit linearly $d0_{in}$ vs $d0_{out}$, which slope should be $\sim 1$, in this way we also know there is no rotation and no scale mismatch between the two detectors planes (confirmation of alignment);  
    3) also check fit quality skewness (asymmetry of a data distribution around its mean) to ensure no asymmetric tail or contaminating subpopulation.  

Results:  
|  | **x** |  |  |  |  |  | **y** |  |  |  |  |  |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
|  | **offset** | **sigma** | **X2/ndf** | **skewness** | **slope** | **intercept** | **offset** | **sigma** | **X2/ndf** | **skewness** | **slope** | **intercept** |
| **8430** | 0.0028 +/- 0.0001 mm | 0.0911 mm | 1.521 | -0.0074 | 0.9954 | 0.0065 mm | 0.0046 +/- 0.0001 mm | 0.0915 mm | 1.527 | -0.0058 | 0.9987 | 0.0056 mm |
| **8431** | 0.0021 +/- 0.0001 mm | 0.0912 mm | 1.548 | -0.0014 | 0.9957 | 0.0051 mm | -0.0049 +/- 0.0001 mm | 0.0918 mm | 1.574 | 0.0160 | 0.9988 | -0.0043 mm |
| **8650** | 0.0092 +/- 0.0001 mm | 0.0989 mm | 1.145 | -0.0173 | 0.9966 | 0.0102 mm | 0.0040 +/- 0.0001 mm | 0.0997 mm | 1.172 | -0.0110 | 0.9989 | 0.0044 mm |
| **8655** | 0.0109 +/- 0.0001 mm | 0.0987 mm | 1.335 | -0.0202 | 0.9965 | 0.0111 mm | 0.0148 +/- 0.0001 mm | 0.0993 mm | 1.422 | -0.0403 | 0.9985 | 0.0151 mm |
| **8656** | 0.0080 +/- 0.0001 mm | 0.0988 mm | 1.348 | -0.0132 | 0.9965 | 0.0078 mm | 0.0109 +/- 0.0001 mm | 0.0996 mm | 1.368 | -0.0296 | 0.9987 | 0.0110 mm |

## To do preliminary
- explore 8650:  
    $\to$ compared in [compare_histo.py](../recoDataSimple/preliminary_analysis/compare_histo.py)
    - Luigi and Melanie said there were some problems during the data taking and with one magnet: at energy 150GeV/c and by beamline hardware constraint also the current was lowered (otherwise particles would not follow the intended trajectory)
    - wider critical angle
    - channeling peak sigma is wider, as lower momentum particles get kicked around more (more multiple coulomb scattering, as mcs scales as $1/p$)
    - lower efficiency, as dechanneling probability increases at lower energy (trajectory more easily perturbed out of the channel)
    - why channeling peak shifted by $~70\,\mu\text{rad}$? the tracker calibration itself is not the source of the anomaly (see [offset results](#calibrate-the-detectors-and-alignment-check))

- look into the spikes: select bins with spikes vs the ones that have not: 
    - the periodic spikes seen in the downstream impact position distribution of channeled particles are an instrumental artifact (purely geometric/instrumental), not a beam or crystal effect: each tracking plane measures a hit only to the precision of one readout strip, so the position from any single plane is quantized in steps of the strip pitch $p$. A track's slope is reconstructed from the difference between hits on two planes separated by baseline $L$, so the slope itself is quantized in steps of $p/L$. When this slope is extrapolated a distance $D$ back to the crystal exit (or any reference plane), the reconstructed position inherits a quantization step of order $Δx ≈ p × D / L$ (a short baseline $L$ amplifies this step, a long baseline suppresses it). In our setup, the upstream arm spans $≈10\text{ m}$ while the downstream arm spans only $≈0.5\text{ m}$: a factor of $~20$ difference in $L$. For a comparable extrapolation distance $D$, this means the downstream reconstructed position is quantized in steps roughly $20×$ coarser than upstream. That coarse grid becomes visible as periodic spikes in `d0Out_x`, while the upstream arm's much finer slope resolution keeps its position distribution effectively continuous.
- look for each bin of efficiency mapping how different are the plots deltatheta vs theta for different efficiency

## Crystal edges found for each run  

Along $y$:
- gaussian distribution both incoming and outgoing
- consistency check: mean and edges of incoming/outgoing are compatible
- they are statistically identical since $\Delta y < \sigma_\mu$, compatibility of 0.3
- doing also here a sliding window that mazimize the integral

Along $x$:
- incoming/outgoing cannot be averaged due to $\Delta x = 0.22\,\text{mm}$
- only using outgoing distribution because has cleaner edges
- i need a quantitative reason of "more cleaner" like the width of the transition or the slope

$\Delta x$:  
Knowing that the detectors are aligned, the shift we see has only physics reasons.  
We can compute precisely the $x$ shift for channeled particles:  
- curvature radius $R=\frac{L}{\theta}\sim 12$ m  
- the lateral displacement at the end of an arc is $\Delta x = R (1 - \cos\theta)$
- considering $\cos\theta=1-\theta^2/2$ (small angles)
- we obtain $\Delta x = \frac{L\cdot\theta}{2}$

| **run** | **$x_{min}$** |**$x_{max}$**|**$y_{min}$** |**$y_{max}$** |
|-|-|-|-|-|
|**8430**|[-1.1044 ± 0.1055 mm | 0.8965 ± 0.0933 mm]|[-5.0250 ± 0.1055 mm | 7.7750 ± 0.0933 mm]|
|**8431**|[-1.1044 ± 0.1055 mm | 0.8956 ± 0.0933 mm]|[-5.1000 ± 0.1055 mm | 7.7000 ± 0.0933 mm]|
|**8650**|[-0.2744 ± 0.1155 mm | 1.7256 ± 0.1117 mm]|[-5.7450 ± 0.1156 mm | 7.0550 ± 0.1117 mm]|
|**8655**|[0.5536 ± 0.1155 mm | 2.5536 ± 0.1117 mm]|[-6.9150 ± 0.1156 mm | 5.8850 ± 0.1117 mm]|
|**8656**|[0.5546 ± 0.1155 mm | 2.5546 ± 0.1117 mm]|[-6.8750 ± 0.1156 mm | 5.9250 ± 0.1117 mm]|

So also we can correct the Lindhard angle estimation:
$$\theta_L=\bigg(1-\frac{\rho_c}{\rho}\bigg)\sqrt{\frac{2U_0}{E}}$$
where:
- $\rho_c=\frac{E}{U'(x_c)}$ is the critical radius
- $E = p\beta c = 180\, \text{GeV}$ or $150\, \text{GeV}$
- $U'(x_c)=5.7\text{ GeV/cm}$ for the (110) plane in silicon
- $U_0=16\text{ eV}$ for silicon
- $\rho=\frac{L}{\theta_b}$

so for the $L=74\text{ mm}$ long crystals, with a $\theta_b=6\text{ mrad}$:  
$\big(1-\frac{180\cdot 10\cdot 6 \cdot 10^{-3}}{74\cdot5.7}\big)\sqrt{\frac{2\cdot 16}{180\cdot 10^9}}=12.992\,\mu\text{rad}$ or  
$\big(1-\frac{150\cdot 10\cdot 6 \cdot 10^{-3}}{74\cdot5.7}\big)\sqrt{\frac{2\cdot 16}{150\cdot 10^9}}=14.294\,\mu\text{rad}$

# Analysis
## Measure the channeling efficiency $\epsilon_{ch}=\frac{N_{ch}}{N_{tot}}\cdot 100$

- center the angle and choosing $\theta_0$ for torsion correction $\tau_y$, which is the variation of the crystalline plane orientation along the direction ($y$) perpendicular to the bending plane ($x$), and must be accounted for (the different orientations of the crystallographic planes along the crystal front surface influence the angular range in which particles can be channelled, i.e., the cut at $\pm 1/2 \,\theta_L$ may not be centred on zero, but a correction parameter $\theta_0 = \theta_0 (y)$)
    - in [compute_edges1.py](../recoDataSimple/compute_edges1.py) we treat the crystal as if it was ideal, and we find a single global $\theta_{optimal}$ from the peak of the angular scan: we find the function of $\theta_0(x)$ and find the peak - the one that maximize the efficiency.
    - what the articles does is that it takes the beam “spot” (`d0_x` vs `d0_y`) and divides it into a checkerboard grid. For each square on the grid, it performs an angular scan (like my previous h_scan) and finds the specific peak angle $\theta_0$ for that square. It then reapplies a correction to each particle by subtracting the local $\theta_0$
- $N_{tot}:$ filter the data set, choosing the events that produced a single track, have a good $\chi^2$ value (?), have actually entered the crystal, and have an angle within $\pm\frac{1}{2}\theta_L$ (actually we choose them $|\theta_{in,x}-\theta_0|\leq\frac{1}{2}\theta_L$):
    - to filter spatially those who entered the crystal (of dimension width $\times$ length), i look at `d0` variables $x$ and $y$ (`d0` and `d0Out`), total and from those events that will present some channeling (so $\Delta\theta_x>\sim 5850 \mu$ rad, number found as $\mu-4\sigma$ by preliminary fit) and from these last ones filtered `d0Out`, which form an almost perfect rectangular, i infere the position of the crystal with respect to the beam through a sliding window method (which maximize the integral).
    - to filter the particles that have entered with the right angle, we first need to find $\theta_0$ in order to apply the filter $|\theta_{in,x}-\theta_0|\leq\frac{1}{2}\theta_L$:
- identify $N_{ch}$: on the $\Delta\theta_x$ filtered histogram perform a gaussian fit on the channeling peak, only from the right side of the peak, then do the integral around the peak
- compute $\epsilon_{ch}=\frac{N_{ch}}{N_{tot}}\cdot 100$ amd its error - we choose a binomial distribution

## Measure the torsion correction
Checkerboard grid search (1 mm $\times$ 1 mm) of torsion correction (see [week 2 analysis](#measure-the-channeling-efficiency)) $\to$ [compute_channeling.py](../recoDataSimple/compute_channeling.py):  
- for every small square (the $xy$ beam distribution is divided into) we find at which $\theta_{x,in}$ the efficiency from the cut $\frac{1}{2}\theta_L$ is maximized
- do a 2D torsion map ($x$, $y$, and angle shift). these values are then fitted ussing a function along $x$ and $y$ in the entire crystal surface to extract the continuous distribution of the map ($z=p_0+p_1 x + p_2 y + p_3 y^2$ where $p_0$ is $\theta_{baseline}$, $p_1$ is $\tau_x$, and $p_2$ is $\tau_y$). we will obtain an average value for the torsion on the $y$ direction, while the torsion along the $x$ direction is negligible ($\tau_x\sim 0$).  
- cut will be $|\theta_{in,x}-(\theta_{0,center}+\tau_y\cdot d0_y)|\leq\frac{1}{2}\theta_L$ (this means that you are accepting only those particles whose entry angle ($\theta_{in}$) is at most half the Lindhard angle ($0.5 \cdot \theta_L$) away from the optimal local angle of the crystal plane at that height $y$ (which is calculated as $\theta_{off} + \tau_y \cdot y$). This is the most elegant and correct way to apply the cut, taking into account both the goniometer offset ($\theta_{off}$) and the twist ($\tau_y$) simultaneously.)  
- computing the final $\epsilon_{ch}$, after applying an angular shift to the incoming direction of the particles, accordingly to the value of the torsion map at their impact position on the crystal surface  
    $\to$ we obtain a global efficiency curve


## Mapping the torsion
I choose a fit parabolic function (i want to use the simplest model that can accurately describe the data, but not so simple that it hides the actual physics):  
$\theta(x,y) = \theta_{baseline} + \tau_x \cdot x + \tau_y \cdot y + p_3\cdot y^2$  
- the `full_quadratic` function with 6 parameters can overfit the data and its noise, while a `pure_parabolic_y` is automatically imposing $\tau_x=0$ and it is too risky  
$\to$ first, i verify with the `full_quadratic` function that the curvature along $x$ is negligible, and then for the final fit i use the `parabolic_y` function to have a more robust fit  
- plot of how uniform is $\tau_y$ across $(x,y)$: plot the residuals from the fit  
- plot of how uniform is $\epsilon_{ch}$ across $(x,y)$  
- understand the meaning of the global efficiency curve plot: how aligned is the crystal after the torsion correction, the efficiency peak should be centered at zero. The curve is obtained from all the data, through a scan along $\theta_{in}$.

Final results:  
|_parabolic model_|     | $\tau_y$        | $\tau_x$         | $\epsilon_{ch}$                     | angle               | sigma             | events discarded |
|---|-------------------|-----------------|------------------|-------------------------------------|---------------------|-------------------|------------------|
| 1 | 8430              | 6.25(4) urad/mm | -0.28(4) urad/mm | (19.6 ± 0.1 [stat] ± 0.0 [syst])%   | (6008.4 ± 0.2) urad | (21.0 ± 0.2) urad | 98%              |
| 1 | 8431              | 5.71(4) urad/mm |  0.28(4) urad/mm | (19.1 ± 0.1 [stat] ± 0.0 [syst])%   | (6008.4 ± 0.2) urad | (20.4 ± 0.2) urad | 98%              |
| 1 | 8430 & 8431       | 7.35(3) urad/mm | -0.70(4) urad/mm | (19.6 ± 0.0 [stat] ± 0.0 [syst])%   | (6008.2 ± 0.1) urad | (21.0 ± 0.1) urad | 98%              |
| 1 | 8650              | 3.89(1) urad/mm | -0.19(4) urad/mm | (17.1 ± 0.1 [stat] ± 0.0 [syst])%   | (6081.7 ± 0.1) urad | (26.0 ± 0.2) urad | 96%              |
| 2 | 8655              | 3.67(2) urad/mm |  0.83(5) urad/mm | (18.0 ± 0.1 [stat] ± 0.0 [syst])%   | (6122.6 ± 0.3) urad | (32.2 ± 0.2) urad | 97%              |
| 2 | 8656              | 3.80(1) urad/mm |  0.87(4) urad/mm | (18.0 ± 0.1 [stat] ± 0.0 [syst])%   | (6122.0 ± 0.3) urad | (32.4 ± 0.2) urad | 97%              |

## To do about torsion

- do not consider the edges in the analysis (cannot measure a local torsion angle where there's no beam so only considering the beam-illuminated window):
    - need a different consideration for edges for `filter2_spatial_cut` and for the fit: crop to where the data lives
    - torsion-map grid range: i can consider along y the area $\mu\pm 3 \sigma$ from the beam or the quantiles (does not assume gaussian tails), should be a separate, narrower range that tracks where the beam actually has statistics
    - if i restrict also the window for the spatial cut i am losing statistics, but i am also considering the data (maybe not from the fit but) still for the efficiency calculation


# Error Analysis
- Week 12: la sezione “Digits” (cifre significative), dentro “To do about MCS”.

## Taking into account d and theta errors

- consider angle errors:
    - $\theta_{in,x}$ and $\theta_{in,y}$ are $\sim8.89\,\mu\text{rad}$ for 8430/8431 and $\sim9.67\,\mu\text{rad}$ for 8650/8655/8656
    - $\theta_{out,x}$ and $\theta_{out,y}$ are $\sim0\,\mu\text{rad}$
    - considering as it was an error and also downstream angle have that error? so $\Delta\theta$ gets an error of $\sim12.57\,\mu\text{rad}$ or $\sim13.67\,\mu\text{rad}$
    - no need of taking them into account in the edge finding

- consider `d0` errors:
    - $d_{in,x}$ and $d_{in,y}$ have $\sim0.105473\text{ mm}$ error for 8430/8431 and $\sim0.115547\text{ mm}$ for 8650/8655/8656
    - $d_{out,x}$ and $d_{out,y}$ have $\sim0.093297\text{ mm}$ error for 8430/8431 and $\sim0.111709\text{ mm}$ for 8650/8655/8656

in [`compute_edges1.py`](../recoDataSimple/compute_edges1.py):
- adding in quadrature: `d0` errors (tracker resolution on the coordinates) + uncertainty due to binning choice
- no error from fitting 

in [`compute_edges2.py`](../recoDataSimple/compute_edges2.py):
- resolution subtraction: since my mcs width is a convolution of physics and tracker resolution (the resolution adds in quadrature to a gaussian width, not linearly):
    - same order of 
- the transition parameter (4) is where d0 resolution shows up for edge finding - it's already fit from data but i can do a cross-check (if the fitted transition is noticeably larger than the mean resolution, that's telling that there's real edge roughness beyond pure detector blur):
    - fitted transition x is $\sim 0.092$ and resolution is $0.105$ so the edge is what we expect from the tracker resolution (good sign)
- systematic smearing: quantify how much my cuts are affected by finite d0 and theta resolution

in [`compute_channeling.py`](../recoDataSimple/compute_channeling.py):
- error on lindhard cut is computed doing a cut within $\theta_L/2\,\pm ...$ instead of just $\theta_L/2$ in order to compute an additional sys error for theta

## Computing efficiency errors

- torsion error (by the fit):
    - now every square bin has local $\theta_0\pm\sigma$ and local $\epsilon_{ch}\pm\sigma$
    - now `h2_torsion_map` carries real per-bin errors (from the local Gaussian fit above) instead of ROOT's default sqrt(content) fallback. This makes the fit below a proper inverse-variance-weighted chi2 fit: bins near the edges of the acceptance get a large local_theta_0_err and are automatically down-weighted, instead of contributing to the surface fit with the same weight as clean core bins.
    - `eff_err_stat` is the binomial (statistical) error on eff_ch, as before. We add a systematic term from the torsion-map baseline uncertainty theta_0_baseline_err: re-run the Lindhard cut + efficiency with theta_0 shifted by +/-1 sigma and take half the resulting spread in eff_ch. This is deliberately a re-evaluation rather than an analytic slope propagation, because the sensitivity of the acceptance is NOT uniform: it is large where the local efficiency amplitude is still near its ~25% core value, and vanishes towards the edges where the local efficiency amplitude drops to ~0%. Re-running the actual cut captures that shape automatically instead of assuming one fixed slope.

- efficiency error for single bins, and then combined (std dev which takes into account the spread between 15 to 25%):  
this would be wrong, because it would be before torsion correction and lindhard critical angle cut, and also the bins have different statistics, it is correct instead to recompute the efficiency over the full dataset  
    - we can treat `h2_eff_map` as validation to check spatial uniformity across the crystal surface, and if done after is to diagnose if the ploynomial is correct or is missing a term or if a region of the crystal is channeling differently (edges/miscuts)
    - also `plot_global_efficiency_curve` can be used as diagnostic: if the peak efficiency is centered at $\theta=0$ after correction, the torsion fit is correct (a mean of $\sim 0.56 \pm 0.02$ it is not compatible with zero but is considered negligible compared to the width of the gaussian distribution); it is checking that the correction surface correctly recenters the whole angular distribution

- the fit uncertainty on the gaussian fits gives <0.001%, probably due to the high statistics   
(i shifted mean and sigma by their fit error and tried to redo the efficiency computation, errors so small that nothing changes)

- also checked for time correlation via block splitting: divide the run in 10 blocks and redo the calculation of the efficiency, if ratio >1 there could be non homogeneous behaviour during the run  
RMS across blocks     = 0.238 %  
Mean binomial error   = 0.200 %  
Ratio (RMS/stat_err)  = 1.19  (consistent with binomial)  


## Systematics 

*the systematic errors are propagated on the final result, not on the intermediate fit parameters*

HOW:  
1) recompute `filter3_Lindhard_cut()` 
2) recompute `channeling_efficiency()` 
3) systematic error as mean between the two computed efficiencies $|\epsilon_{up}-\epsilon_{down}|/2$  
*And to find the total error:*  
4) add the errors in quadrature, assuming indipendence.Two of these are not obviously independent: the spatial-box shift changes which particles enter the torsion-map fit, which in turn slightly changes fit_params - but i am reusing the nominal fit_params for the box variation, which effectively assumes independence by construction. *"systematics treated as uncorrelated; box and torsion-map re-fit not jointly propagated."*

MAIN SOURCES:
- *fit model choice* - $\tau_x,\,\tau_y,\,c_y$ dependence:  
`parabolic_y` vs `full_quadratic` differ in $\tau_y$, rerun the whole efficiency chain and take the eff difference (one-sided difference)
- torsion map/grid bin choice - nx_slices, ny_slices dependence:  
rerun torsion_map with a coarser (nx_slices=5) or finer x-binning and check $\tau_x$ stability (if $\tau_x$ shifts by more than its statistical error, that's a systematic worth quoting)
- n sigma integration for $N_{ch}$ - bin min and bin max dependence:  
Varying that lower bound (e.g. $2.5σ$ vs $3.5σ$) and taking half the spread in efficiency is a clean shift-and-rerun systematic, and it's cheap since it doesn't require re-running the spatial/torsion/Lindhard chain — only the final channeling_efficiency integration step.
- preliminary channeling selection cut $\mu-3\sigma$:  
This threshold (preliminary_cut_on_deltatheta) determines who counts as "channeled" for computing the torsion map and spatial box in the first place — it's used upstream of everything. Varying it (e.g. $μ − 2.5σ$ vs $μ − 3.5σ$) and checking how much $τ_x$, $τ_y$, and the final edges shift would test whether your whole downstream chain is sensitive to this somewhat arbitrary choice.
- crystal dimension box:  
shifting the x margins of $\pm$ the $d_0$ error ($0.1\text{ mm}=5%$)  
as expected, the efficiency is higher if the margins are narrower and lower if margins are broader, meaning that the statistical distribution decreases there and/or because the local angular acceptance is less precise at the edges (miscut, edge damage) *(consistent with the fact that we cut the edges when computing the efficiency - the signal is less reliable near the edges)*: this is not an artifact of the method: 0.105 mm out of 2 mm is literally the spatial resolution relative to the smallest dimension of the crystal, so it makes sense that the efficiency would be sensitive to this.
- crystal cut margins:  
defines where i compute the final efficiency, subset of the entrance surface - needa refit of the torsion map with different margins
- critical angle $\theta_L$
- bkg subtraction? NO
- shift of $\theta_0$ by $\pm$ their fit error?   
this is testing how much the efficiency is changed if we shift the baseline by $\pm 1\sigma$ of its fit error

what i did:  
*(easy)*  
1. shift of $\theta_0$ by $\pm$ their fit error  
even if small is worth mentioning it in the total budget  
0.004%  
2. shifting lindhard angle by $\theta_in$ error:  
my angular resolution (~8.9 µrad) is larger than $θ_L/2$ ($~6.5 µ\text{rad}$), so treating "shift the cut edge by one full resolution sigma" as systematic is too aggressive: a more defensible (and numerically stable) choice is to shift by a fraction of the resolution or a smaller fixed step, rather than the full σ  
0.981%
3. shifting the box margins by $d_0$ error along x *(main source)*:  
1.098%  
-  shifting the box margins by $d_0$ error along y (shift of 2.5% of y length):  
0.001%
4. n sigma integration window (changing bin min (2.5σ - 3.5σ) from which i integrate n_ch):  
0.054%

*(heavier)*
1. fit model choice (full quadratic vs parabolic y):  
0.031%
2. nx/ny slice choice:  
0.003%
3. crystal margins cut:
0.077%
4. preliminary cut threshold:  
0.063%

| %        | **total** | binomial| $\theta_0$ | $\theta_L$ | shift $x$ | shift $y$ | bin min | fit model choice | nx ny choice | crystal margins cut | preliminary cut |
|----------|-----------|---------|------------|------------|-----------|-----------|---------|------------------|--------------|---------------------|-----------------|
| **8430** | **1.5**   | 0.063   | 0.004      | 0.981      | 1.098     | 0.001     | 0.054   | 0.031            | 0.003        | 0.077               | 0.063           |
| **8431** | **1.5**   | 0.064   | 0.007      | 0.957      | 1.079     | 0.001     | 0.054   | 0.188            | 0.130        | 0.017               | 0.072           |
| **8650** | **1.3**   | 0.051   | 0.003      | 0.755      | 1.001     | 0.001     | 0.066   | 0.004            | 0.034        | 0.009               | 0.000           |
| **8655** | **1.4**   | 0.062   | 0.009      | 0.737      | 1.187     | 0.002     | 0.076   | 0.016            | 0.199        | 0.180               | 0.036           |
| **8656** | **1.4**   | 0.058   | 0.006      | 0.814      | 1.191     | 0.001     | 0.080   | 0.010            | 0.104        | 0.164               | 0.094           |

## PDG on errors

**40.2.6**  
Including additional parameters may eliminate or at least reduce the effect of systematic uncertainties, their presence will result in increased statistical uncertainties for the parameters of interest. This occurs because the estimators for the nuisance parameters and those of interest will in general be correlated.  
($\to$ correlation between $\tau_y$ and $\epsilon$ uncertainties)  
$\to$ standard method is *profile likelihood*: fixing all nuissance parameters to their optimal value (instead of their nominal value)

**40.2.2.1**  
In the large-sample limit the likelihood function itself takes on a Gaussian shape, so that the log-likelihood becomes a quadratic function.
$\to$ la matrice di covarianza dei parametri del fit si ottiene dalle derivate seconde della log-likelihood al minimo, esattamente quello che fit_result. (GetCovarianceMatrix("S"))

**40.3.3.1**  
Check $\chi^2/ndf$ of the torsion fit.  
Poor goodness-of-fit, however, does not mean that one will have large statistical errors for parameter estimates. If, for example, the error bars (or covariance matrix) used in constructing the $\chi^2$ are underestimated, then this will lead to underestimated statistical errors for the fitted parameters and an increased value of the minimized $\chi^2$.  

**40.4.2.4**  
Estimating the binomial error for the efficiency is fine since $N$ is large and $p$ not so near 0% or 100% (we have ~20%).  
Also exist the alternative Clopper-Pearson, for less statistics (es. se calcoli efficienze locali per singolo bin della torsion map con poca statistica, dove l'approssimazione gaussiana può fallire).

**40.4.2.3** and **Eq.40.80**  
If N independent measurements result in log-likelihood functions $\ln Li(θ)$, then the combined log-likelihood function is simply the sum, $\ln L(θ) = \sum \ln Li(θ)$. This can then be used to determine an approximate confidence interval or region.  
Questo è il modo statisticamente corretto di combinare stime da run diverse (es. τ_y da 8430 e 8431 separatamente) invece di fare una media pesata ad-hoc — rilevante se in futuro si vuole combinare i risultati dei diversi run in un unico numero finale per la tesi.

## Digits
The golden rule of metrology (ISO/GUM) stipulates that the decimal places of the central value must align with the position of the first significant digit of the total uncertainty associated with that estimate.  
Check the detector precision and change the uncertainty di conseguenza, not possible to have so many significance digits 0.1 um:
- trackers have resolution of $7\,\mu\text{m}=0.007\,\text{mm}$
- writing $0.0001\,\text{mm}$ is like claiming we have a sensibility 70 times smaller than it really is
- it is better to round up to $0.01\,\text{mm}=10\,\mu\text{m}$ so compared to the $7\,\mu\text{m}$ resolution, we are safe  

**$\to$ approximate everything to 2 decimal digits if in mm**

Also about the angular uncertainty, that actually came from the position uncertainty: the angular resolution $\Delta\theta$ is determined by the spatial resolution of the detector (those $7\,\mu\text{m}$) divided by the distance between the tracker planes:
$$\theta=\frac{x_2-x_1}{L}\iff\sigma_\theta=\frac{\sqrt{2}\,\sigma_x}{L}\sim5-10\,\mu\text{rad}$$
$\to$ so in this case: **approximate everything to none decimal digits if in urad**

*(In experimental physics and engineering, significant figures are not just a mathematical formality: they represent the physical limit of what your instrument is capable of measure. Writing more significant figures than necessary is equivalent to claiming to have information on a scale that your detector is physically incapable of measuring.)*



# Multiple Coulomb Scattering
- Week 6: solo i punti su compute_edges2.py in “Taking into account d and theta errors”, se ti servono.
- Week 8: tutta la sezione “Find the crystal edges: scattering method”.
- Week 12: “To do about MCS”, solo la parte “Amorphous scattering” (calcolo di σ_bsl, MCS in aria, ecc.).
- Week 13: “About MCS: adding different contributions” e “My analysis”.

The question we need to ask ourselves is: on average, how much will this particle be deflected at the end of its path? Most of the particles will be deviated of an angle near zero (center of the gaussian distribution).  
The width of this curve is the key value: it tells us how much, on average, the particle is deflected. This value is called the root-mean-square (RMS) angle (angolo quadratico medio) and is denoted by $\theta_0$, and can be computed through the Highland formula, which is derived by Rutherford scattering cross section formula.  
In the formula we have: $\theta_0\propto 1/p$, $\propto 1/\beta c$ (the higher the velocity and more massive the particle is, the less deviated it gets), and $\propto z$. Then we also have the material's characteristics: the width $x$ and the radiation length $X_0$ (depends on the 'density').  

$\to$ from [PDG: Passage of Particles Through Matter](https://pdg.lbl.gov/2023/reviews/rpp2023-rev-passage-particles-matter.pdf#section.34.3)

A charged particle traversing a medium is deflected by many small-angle scatters. Most of this deflection is due to Coulomb scattering from nuclei as described by the Rutherford cross section. (However, for hadronic projectiles, the strong interactions also contribute to multiple scattering.)  
For many small-angle scatters the net scattering and displacement distributions are Gaussian via the central limit theorem. Less frequent “hard” scatters produce non-Gaussian tails *(= particelle che prendono grandi spallate perche colpiscono nucleo in pieno)*. These Coulomb scattering distributions are well-represented by the theory of Molière (relativistic pions, kaons, and protons). If we define  

$\theta_0=\theta_{plane}^{rms}=\frac{1}{\sqrt{2}}\theta_{space}^{rms}$  

then it is sufficient for many applications to use a Gaussian approximation for the central 98% of the projected angular distribution, with an rms width given by Lynch & Dahl:  

$\theta_0=\frac{13.6\text{ MeV}}{\beta c p}z\sqrt{\frac{x}{X_0}}\big[1+0.088\log_{10}\big(\frac{x z^2}{X_0 \beta^2}\big)\big]$  
$=\frac{13.6\text{ MeV}}{\beta c p}z\sqrt{\frac{x}{X_0}}\big[1+0.038\ln\big(\frac{x z^2}{X_0 \beta^2}\big)\big]$  

Here $p$, $βc$, and $z$ are the momentum, speed, and charge number of the incident particle, and $x/X_0$ is the thickness of the scattering medium in radiation lengths. This takes into account the $p$ and $z$ dependence quite well at small $Z$, but for large $Z$ and small $x$ the $β$-dependence is not well represented. This equation describes scattering from a single material, while the usual problem involves the multiple scattering of a particle traversing many different layers and mixtures. Since it is from a fit to a Molière distribution, it is incorrect to add the individual $θ_0$ contributions in quadrature *(= non si puo' sommare e basta i contributi aria + silicio + rame ecc, ci sarebbe grave sottostima deviazione)*; the result is systematically too small. It is much more accurate to apply this equation once, after finding $x$ and $X_0$ for the combined scatterer.  

The nonprojected (space) and projected (plane) angular distributions are given approximately by  

$\frac{1}{2\pi\theta_0^2}\exp\big(-\frac{\theta^2_{space}}{2\theta_0^2}\big)d\Omega$  
$\frac{1}{\sqrt{2\pi}\theta_0}\exp\big(-\frac{\theta^2_{plane}}{2\theta_0^2}\big)d\theta_{plane}$  

where $θ$ is the deflection angle. In this approximation, $\theta^2_{space} ≈ (θ_{plane,x}^2 + θ_{plane,y}^2)$, where the $x$ and $y$ axes are orthogonal to the direction of motion, and $dΩ ≈ dθ_{plane,x} dθ_{plane,y}$. Deflections into $θ_{plane,x}$ and $θ_{plane,y}$ are independent and identically distributed. Fig. 34.10 shows these and other quantities sometimes used to describe multiple Coulomb scattering. They are  

$\psi^{rms}_{plane}=\frac{1}{\sqrt{3}}\theta_{plane}^{rms}=\frac{1}{\sqrt{3}}\theta_0$    
$y^{rms}_{plane}=\frac{1}{\sqrt{3}}x\theta_{plane}^{rms}=\frac{1}{\sqrt{3}}x\theta_0$  
$s^{rms}_{plane}=\frac{1}{4\sqrt{3}}x\theta_{plane}^{rms}=\frac{1}{4\sqrt{3}}x\theta_0$  

All the quantitative estimates in this section apply only in the limit of small $θ^{rms}_{plane}$ and in the absence of large-angle scatters. The random variables $s$, $ψ$, $y$, and $θ$ in a given plane are correlated. Obviously, $y ≈ xψ$. In addition, $y$ and $θ$ have the correlation coefficient $ρ_{yθ} = \sqrt{3}/2 ≈ 0.87$ *(1 would be perfect correlation)*.  

For Monte Carlo generation of a joint ($y_{plane}$, $θ_{plane}$) distribution, or for other calculations, it may be most convenient to work with independent Gaussian random variables ($z_1$, $z_2$) with mean zero and variance one, and then set  

$y_{plane}=z_1\, x\, \theta_0(1-\rho^2_{y\theta})^{1/2}/\sqrt{3}+z_2\,\rho_{y\theta}\, x\, \theta_0/\sqrt{3}=z_1\, x\, \theta_0/\sqrt{12}+z_2\, x\, \theta_0/2$  
$\theta_{plane}=z_2\,\theta_0$  

In this way the computer instantaneously gets an exit angle and a final position, which are correlated (no need of simulating all the small scatterings happening inside).  
Note that the second term for $y_{plane}$ equals $x\,θ_{plane}/2$ and represents the displacement that would have occurred had the deflection $θ_{plane}$ all occurred at the single point $x/2$.  

$\to$ **multiple coulomb scattering RMS of protons ($z=1$) for the $L=74\text{ mm}$ long silicon crystals ($X_0\sim9.37\text{ cm}$):**  
$\theta_0=66.5 \,\mu\text{rad}$

## Measuring multiple coulomb scattering

$\to$ performed in [compute_edges2.py](../recoDataSimple/compute_edges2.py) (this method provides indipendence from channeling efficiency or torsion)

My current method ([compute_edges1.py](../recoDataSimple/compute_edges1.py)) to find the crystal edges is isolating the channeled particles, then use their impact position distribution to find where the population is the densest and call that the crystal footprint, but this only uses the channeled fraction, so ~20% of all the particles that actually crossed the crystal. what Luigi was suggesting instead is instead of relying on channeling at all, i can use the fact that every particle that physically traverses the crystal, picks up extra multiple coulomb scatterings that a particle passing beside the crystal does not.

- in small bins $(x,y)$ across the beam spot, compute the width (gaussian sigma or rms) of the outgoing angular distribution $\Delta\theta_x$ (or $\theta_{out}$) rather than its channeled peak position
- i am considering $\Delta\theta_x$ only on a range $\pm 150\,\mu\text{rad}$, in order to exclude channeling particles, while mcs is a purely statistic process and can be described by a gaussian curve around zero.
- map the local width as a function of impact position: outside the crystal, the width sits at the baseline, inside the crystal it should jump to $\theta_0\sim66.5\,\mu\text{rad}$ so the edges are were there should be a sharp, well-defined transition - that can be fitted through a step function  

1. in reality $\sigma_{measured}$ has also detector resolution contribution (usually measured through a run without target) (added in quadrature): $\sigma_{meas}=\sqrt{\sigma_{mcs}^2+\sigma_{track,in}^2+\sigma_{track,out}^2}$
2. why outside the crystal we don't get zero but $\sim12\,\mu\text{rad}$? even without the crystal, particles would still not travel with a perfect linear track, the beam has an intrinsic angular divergence (should be of $\sim12-15\,\mu\text{rad}$)

RESULTS:  
- i obtain a width in y of $\sim 8.4\text{ mm}$ instead of $12.8\text{ mm}$: re-try considering only one slide, but anyway could be because what we're getting it the *primary clean beam* (don't worry that for the torsion analysis we are not considering the edges anyway) $\to$ infact we are detecting the clamps and the clean beam
- along x i obtain exactly $\sim2\text{ mm}$
- along x outside the crystal i have $\sigma\sim12\,\mu\text{rad}$ which corresponds to the detector resolution, and inside the crystal i have   


## Find the crystal edges: scattering method

Results from [compute_edges2.py](../recoDataSimple/compute_edges2.py):

$\theta_{MCS} = \sqrt{\sigma_{high}^2​−\sigma_{bsl}^2​} = \sqrt{(\sigma_{bsl} + \sigma_{jump})^2​−\sigma_{bsl}^2​}$  
*but the jump $(\sigma_{high} - \sigma_{bsl})$ seems more compatible with $\theta_{mcs}^{th}$*

**Considerations:**
- besides computing baseline and the high level, which difference in quadrature should give the $\theta_{mcs}$ i can also:
- consider the angle resolution (of $\Delta\theta_x$) of $\sim 8.89$ (8430/8431) or $9.67$ (8650/8655/8656)
- compare the impact position resolution/uncertainty with the width of the erf transition (how blurred and not perfectly sharp) so is a quality check: we expect that gradient not to be arbitrary, but to be explained by the tracker’s spatial resolution (d0Err_x/y, the uncertainty with which you measure the impact position). If the two numbers are comparable, it means that the “blurred” edge you see in the fit is simply a known instrumental resolution effect, not a true, physically gradual edge of the crystal (i.e., it is not edge damage or extensive mechanical damage). If, on the other hand, the fitted transition were much wider than d0Err, it would be a sign that something else is going on.
- **material of the clamps:** from Highland formula $p, \beta, z$ depend on the beam, while $x, X_0$ on the material traversed. Supposing $x$ (width in $z$) as the same of the crystal or slightly more ($\sim 80\,\text{mm}$) i can resolve for $X_0$:
    - stainless steel (Fe based) $X_0\sim 1.76\,\text{cm}$ $\to$ $\theta_{MCS} = 191.99$ (for $l=100\,\text{mm}$) 
    - titanium $X_0\sim 3.56\,\text{cm}$ $\to$ $\theta_{MCS} = 131.60$ (for $l=100\,\text{mm}$)
    - brass/copper $X_0\sim 1.4-1.5\,\text{cm}$ $\to$ $\theta_{MCS} = 209.90/213.77$ (for $l=100\,\text{mm}$)
    - aluminum $X_0\sim 8.9\,\text{cm}$ $\to$ $\theta_{MCS} = 80.44$ (for $l=100\,\text{mm}$)
    - tungsten $X_0\sim 0.35\,\text{cm}$ $\to$ $\theta_{MCS} = 455.32$ (for $l=100\,\text{mm}$)


## To do about MCS:

**Amorphous scattering**
- the slices distribution should always be centered at zero - so now i'm seeing different effects summed up
- only take amorpheous scattering to make sure i'm only selecting particles that never came close to alignment with any plane at any depth so they genuinely never channel, volume-reflect or dechannel - so we are sure they undergo pure random multiple coulomb scattering the whole way through
- are we sure we have to add in quadrature?  
check [week 13: about MCS (subsection)](#about-mcs-adding-different-contributions)
- check rms in air - not exactly zero
- along $y$ do fewer slices

$$\theta_{RMS} = \frac{13.6 \text{ MeV}}{180000 \text{ MeV/c}} \sqrt{T_{tot}} [1 + 0.038 \ln(T_{tot})]$$
and we can use $T_{tot}\approx T_{cryst}$ obtaining $\theta_{RMS} \approx 66 \,\mu\text{rad}$.

$$\theta_{MCS,air}=\frac{13.6\,\text{MeV}}{180\,\text{GeV}} \sqrt{\frac{x}{X_0}} \bigg[1+0.038\ln\Big(\frac{x}{X_0}\Big)\bigg]$$
and choosing $X_0\sim300\,\text{m}$ and $x=74\,\text{mm}$ (does not change much even if it is more) we get $\theta_{MCS,air}=0.812\,\mu\text{rad}$ so the term from the tracker angular resolution is still the dominant one.

$$\sigma_{bsl}=\sqrt{\sigma^2_{tracker,res} + \sigma^2_{MCS, air}}$$
and using the known $\sigma_{tracker}=12.57-13.68\,\mu\text{rad}$ we obtain:  
$\sigma_{bsl}=12.6-13.7\,\mu\text{rad}$ to be compared with the one coming from the fit.   

NB: $\sigma_{div,x}=25.4\,\mu\text{rad}$ and $\sigma_{div,y}=39.2\,\mu\text{rad}$ (which come from $\sigma_{\theta}=\sqrt{\sigma_{divergence}^2+\sigma_{tracker,res}^2}$, where $\sigma_{\theta}$ is from the data - width of the gaussian fit) do not enter the calculations because the angular divergence of the beam is common for $\theta_{in}$ and $\theta_{out}$ for a linear track and so in the difference $\Delta\theta$ it gets cancelled out.


And for the crystal material:
$$\sigma_{crystal}=\sqrt{\sigma^2_{MCS} + \sigma^2_{bsl} +(- \sigma^2_{74mm, air})}$$
even though the $\sigma^2_{74mm, air}$ contribution is really small (well below angular resolution):
$$\Delta\theta=\frac{\sigma^2_{74mm, air}}{2\theta_{cryst}}=\frac{0.812^2}{2\cdot70}\sim0.005\,\mu\text{rad}$$

To have the same kind of comparison inside/outside (since the processes are independent, we can extract the pure scattering from the crystal by inverting the sum of the squares), compute:  
- $\theta_{MCS,air}=\sqrt{\sigma^2_{fit,bsl} - \sigma^2_{tracker}}$ to be compared with $\theta_{MCS,air}^{th}=0.8\,\mu\text{rad}$
- $\theta_{MCS,cryst}=\sqrt{\sigma^2_{fit,high} - \sigma^2_{fit,bsl}}$ to be compared with $\theta_{MCS,cryst}^{th}=66.5\,\mu\text{rad}$

## About MCS: adding different contributions

We are looking at the distribution of the angular deflection, with mean always zero, and we are measuring its spread (of the distribution), which is the variance. In statistic, when you combine indipendent random variables, it is their variances that add together in quadrature. Indeed, what we are measuring is the convolution of two gaussians (widths get add in quuadrature, while means get simle sum).  
Since the RMS angle ($\theta$) represents the standard deviation (the square root of the variance), the basic probability math dictates that the contributions combine via the square root of the sum of their squares. This is the definition of adding in quadrature. If you were to sum the RMS contributions linearly ($\theta_{total} = \theta_1 + \theta_2$), you would mathematically be assuming perfect positive correlation. Physically, this implies that if a particle scatters $1^\circ$ to the left in the first material, it "remembers" this and preferentially scatters to the left again in the next material, actively steering the particle in a continuous curve.  
$$\sigma_{total}^2 = \sigma_1^2 + \sigma_2^2 + 2 \rho \sigma_1 \sigma_2$$
$$\rho=0\,\,\text{(indipendent)}: \sigma_{total}^2 = \sigma_1^2 + \sigma_2^2$$
$$\rho=1\,\,\text{(correlated)}: \sigma_{total}^2 = (\sigma_1 + \sigma_2)^2$$
While adding in quadrature is technically an approximation due to the Highland logarithmic correction, summing the RMS angles linearly is physically and statistically incorrect.  

Adding in quadrature assumes that the angular deflections perfectly follow a Gaussian (normal) distribution. If multiple Coulomb scattering were purely Gaussian, adding the variances in quadrature would be 100% mathematically exact. However, the scattering distribution is not purely Gaussian.  

While the vast majority of scatters are tiny "glancing" interactions that build up into a Gaussian core, particles occasionally undergo a single, hard collision with a nucleus. This is governed by Rutherford scattering, which creates heavy tails in the distribution.

Because of these large-angle Rutherford scatters, the overall shape of the distribution (described by Molière theory) changes depending on how thick the material is. To define a useful "width" (RMS) for the central 98% of this non-Gaussian distribution, physicists use the Highland formula, which introduces a logarithmic correction factor to account for the thickness:
$$\theta_0 = \frac{13.6 \text{ MeV}}{\beta c p} z \sqrt{\frac{x}{X_0}} \left[ 1 + 0.038 \ln\left(\frac{x}{X_0}\right) \right]$$
Because of that $\ln(x/X_0)$ term, the mathematical rule of quadrature breaks. If you take two layers, calculate their Highland $\theta_0$ separately, and add them in quadrature, the logarithmic terms interact incorrectly. Adding Highland angles in quadrature systematically overestimates the true scattering angle, because you end up artificially inflating the impact of the Rutherford tails.  

Adding in quadrature is still widely used in physics as a practical approximation. It holds (or is "good enough") in the following cases: first-order estimates, highly uniform materials, macroscopic vacuum gaps.

**Correct formula:** If you have a particle passing through multiple layers of different materials and you want the exact RMS scattering angle according to the Particle Data Group (PDG), you do not calculate individual angles at all.Instead, you sum the fractional radiation lengths of the materials first, and then apply the Highland formula exactly once.  
1) Calculate the total material thickness in terms of radiation lengths ($T$):$$T = \sum_{i} \frac{x_i}{X_{0,i}}$$
2) Plug the total $T$ into the Highland formula:$$\theta_{total} = \frac{13.6 \text{ MeV}}{\beta c p} z \sqrt{T} \left[ 1 + 0.038 \ln(T) \right]$$By summing the radiation lengths before doing the math, you ensure the logarithmic correction scales correctly for the total amount of material traversed, accurately reflecting how the Gaussian core and Rutherford tails behave in reality.

*For the crystal:*
$$T_{tot} = \frac{x_{aria1}}{X_{0,aria}} + \frac{x_{cristallo}}{X_{0,cristallo}} + \frac{x_{aria2}}{X_{0,aria}}$$
con $X_0^{air} \approx 30420 \text{ cm}$ and $X_0^{Si} \approx 9.36 \text{ cm}$.  
different contributions: with $74 \text{ mm}$ of crystal: $T_{cryst} = \frac{7.4}{9.36} \approx 0.79$ and assumin even ($100 \text{ cm}$) of air: $T_{air} = \frac{100}{30420} \approx 0.003$.

We expect that along $y$ we will have a cleaner measure, because channeling, volume reflection and dechanneling are all phenomenons happening along the bending direction, $x$, while along $y$ we will have pure mcs.

## My analysis:

Analyze the tracking deflections across transverse coordinates to measure MCS.  
By mapping the local angular spread $\sigma(\Delta\theta)$ as a function of the transverse track position ($x$, $y$), material transitions appear as sharp steps convoluted with the telescope's spatial resolution, modeled as error functions.

Channeling, volume reflection, and dechanneling are phenomena related to the orientation of the crystal planes, which in our crystal are cut or bent to deflect in the $x$-direction (this is why the channeling peak appears in $Δθ_x$, not $Δθ_y$). The $y$-direction is transverse to the bending plane, so it is not affected by these effects: it is sensitive only to isotropic MCS (projection of the same multiple Coulomb scattering) and to instrumental resolution.

1. generate 2d distributions of track position versus deflection angle distribution, to check if some cut/filter on the incoming angle is needed (we dont want VR or channeled particles).
2. to filter out channeled particles we simply select $\Delta\theta << \theta_b$, but unlike the short crystals there is no way to select single regions (VR, channeled, dechanneled, amorph) along $\Delta\theta_x$ and indeed the analysis using $\Delta\theta_x$ produces unreliable results (distributions not gaussian and not centered at zero and too large with respect to the theoretical expectations) - consistent with the expectation that no filtering of $θ_{in}$ can remove the VR.
3. using $\Delta\theta_y$ because is not the axis along which channeling happens so the interactions should be pure MCS (isotropic) + detector resolution. 
4. iterative gaussian slices: project each spatial bin of the crystal along the deflection axis ($\Delta\theta_y$), execute an initial gaussian fit within $\pm 2.5\times RMS$ and automatically refines the fit range to $\mu\pm2.5\sigma$ using fitted parameters (an rejecting unstable fits). 
5. finally it records the resulting profile: bin centers ($x$ or $y$) vs the fitted width $\sigma$. 
6. erf fit along $x$ (air $\to$ crystal $\to$ air) and $y$ (air $\to$ clamp $\to$ crystal $\to$ clamp $\to$ air).
7. because the spatial slice in $y$ depende on knowing the $x$-boundaries of the crystal and vice-versa, the function iterate until convergence (alternating convergence loop):
$$\text{Slice } y \in [y_{e2}^{(k)}, y_{e3}^{(k)}] \longrightarrow \text{Fit } x\text{-edges} \longrightarrow \text{Slice } x \in [x_{\text{lo}}^{(k+1)}, x_{\text{hi}}^{(k+1)}] \longrightarrow \text{Fit } y\text{-edges}$$
(checking physical consistency at every step and monitoring the maximum edge coordinate shift)  
8. from the converged fits, we can extract $\theta_{MCS, crystal}$, $\theta_{MCS, air}$, and $\theta_{MCS, clamp}$
9. check on telescope tracking performance: the fitted transition parameter with the track impact parameter uncertainty

$\to$ *“I identified that $x$ is structurally unsuitable for this purpose, I predicted that $y$ would be clean, and the agreement with the theory confirms this.”*

## MCS final results

```
for v in 8430 8431 8650 8655 8656; do python compute_edges2.py "$v"; done
```

run [mm] | fitted transition $x$ | fitted transition $y$ | mean position error $d_0$ 
-------- | --------------------- | --------------------- | -------------------------
8430     |  0.09 ± 0.01      |  0.29 ± 0.01      | 0.1055 
8431     |  0.08 ± 0.01      |  0.28 ± 0.01      | 0.1055
<!-- 8430/8431|  0.0900 ± 0.0002      |  0.2617 ± 0.0064      | 0.1055 -->
<!-- 8650     |  0.0068 ± 0.0001      |  0.1375 ± 0.0062      | 0.1055
8655     |  0.0740 ± 0.0005      |  2.6090 ± 0.0784      | 0.1155 
8656     |  0.0731 ± 0.0005      | 2.6705 ± 0.0557       | 0.1155  -->

run [urad] | resolution | **x:** $\sigma_{bsl}$ |  $\sqrt{\sigma_{high}^2​−\sigma_{bsl}^2​}$ | **y:** $\sigma_{bsl}$ | $\theta_{mcs}=\sqrt{\sigma_{high}^2​−\sigma_{bsl}^2​}$ | $\theta_{clamp}=\sqrt{\sigma_{clamp}^2​−\sigma_{bsl}^2​}$
--         | --         | --                    |   --             | --                    | --             | --             
8430       | 8.9       | 12.4 ± 0.1          |  65.1 ± 0.1    | 11.5 ± 0.1          | 62.9 ± 0.1   | 100.5 ± 0.7 
8431       | 8.9       |  12.4 ± 0.1         | 65.1 ± 0.1    | 11.5 ± 0.1          | 62.9 ± 0.1   | 101.1 ± 0.7 
<!-- 8430/8431  | 8.89       | 12.75 ± 0.00          |  78.44 ± 0.05    | 11.80 ± 0.06          | 63.09 ± 0.03   | 99.35 ± 0.49  -->
<!-- 8650       | 9.67       | 14.73 ± 0.01          |   95.02 ± 0.14   | 13.34 ± 0.17          | 73.54 ± 0.06   | 118.80 ± 1.68 
8655       | 9.67       | 14.74 ± 0.01          |  86.82 ± 0.12    | ?                     | 70.50 ± 0.09   | ?
8656       | 9.67       | 14.76 ± 0.01          |  86.66 ± 0.11    | ?                     | 70.13 ± 0.09   | ? -->

Even if the resolution is $\sim 7$ it is still correct to use decimal digits in the results if we are talking about the statistic uncertainty of the fit, since the estimate is accompanied by sufficient statistical data to justify that level of precision for the parameter.  
$\delta\sigma \approx \frac{\sigma}{\sqrt{2N}}$ if $N \approx 10^6$ eventi, for $\sigma \approx 12$ we obtain $\delta\sigma \approx 12 / \sqrt{2 \cdot 10^6} \approx 0.01$.

| **run**         | **$x$ [mm]**    | $x_{min}$        | $x_{max}$       | **$y$ [mm]**     | $y_{clamp}$       | $y_{min}$        | $y_{max}$       |
|-----------------|-----------------|------------------|-----------------|------------------|-------------------|------------------|-----------------|
| **8430**        | 2.11 ± 0.01 | -1.07 ± 0.0003 | 1.05 ± 0.01 | 14.09 ± 0.04 | -11.07 ± 0.01 | -5.43 ± 0.01 | 8.66 ± 0.03 |
| **8431**        | 2.11 ± 0.01 | -1.07 ± 0.0003 | 1.05 ± 0.01 | 14.03 ± 0.04 | -11.08 ± 0.01 | -5.41 ± 0.01 | 8.62 ± 0.03 |
<!-- | **8430 & 8431** | 2.0907 ± 0.0003 | -1.0600 ± 0.0002 | 1.0308 ± 0.0002 | 13.9346 ± 0.0195 | -11.1049 ± 0.0059 | -5.3966 ± 0.0091 | 8.5380 ± 0.0173 |
| **8650**        | 2.0971 ± 0.0001 | -0.2709 ± 0.0001 | 1.8261 ± 0.0000 | 13.7271 ± 0.0156 | -11.8176 ± 0.0094 | -6.3653 ± 0.0156 | 7.3618 ± 0.0015 |
| **8655**        | 2.1200 ± 0.0006 | 0.6031 ± 0.0005  | 2.7232 ± 0.0004 | 12.6059 ± 0.0717 | -10.3721 ± 0.0581 | -6.9081 ± 0.0570 | 5.6978 ± 0.0436 |
| **8656**        | 2.1224 ± 0.0006 | 0.6015 ± 0.0004  | 2.7239 ± 0.0004 | 13.0073 ± 0.1430 | -10.2544 ± 0.0557 | -7.1044 ± 0.0575 | 5.9029 ± 0.1309 | -->



# Data Results 
- Week 13: “MCS final results”. Le tabelle della Week 8 nel logbook nuovo non ci sono più, quindi vale solo questa.

Final results:  
|_parabolic model_|     | $\tau_y$        | $\tau_x$         | $\epsilon_{ch}$                     | angle               | sigma             | events discarded |
|---|-------------------|-----------------|------------------|-------------------------------------|---------------------|-------------------|------------------|
| 1 | 8430              | 6.25(4) urad/mm | -0.28(4) urad/mm | (19.6 ± 0.1 [stat] ± 0.0 [syst])%   | (6008.4 ± 0.2) urad | (21.0 ± 0.2) urad | 98%              |
| 1 | 8431              | 5.71(4) urad/mm |  0.28(4) urad/mm | (19.1 ± 0.1 [stat] ± 0.0 [syst])%   | (6008.4 ± 0.2) urad | (20.4 ± 0.2) urad | 98%              |
| 1 | 8430 & 8431       | 7.35(3) urad/mm | -0.70(4) urad/mm | (19.6 ± 0.0 [stat] ± 0.0 [syst])%   | (6008.2 ± 0.1) urad | (21.0 ± 0.1) urad | 98%              |
| 1 | 8650              | 3.89(1) urad/mm | -0.19(4) urad/mm | (17.1 ± 0.1 [stat] ± 0.0 [syst])%   | (6081.7 ± 0.1) urad | (26.0 ± 0.2) urad | 96%              |
| 2 | 8655              | 3.67(2) urad/mm |  0.83(5) urad/mm | (18.0 ± 0.1 [stat] ± 0.0 [syst])%   | (6122.6 ± 0.3) urad | (32.2 ± 0.2) urad | 97%              |
| 2 | 8656              | 3.80(1) urad/mm |  0.87(4) urad/mm | (18.0 ± 0.1 [stat] ± 0.0 [syst])%   | (6122.0 ± 0.3) urad | (32.4 ± 0.2) urad | 97%              |


what i did:  
*(easy)*  
1. shift of $\theta_0$ by $\pm$ their fit error  
even if small is worth mentioning it in the total budget  
0.004%  
2. shifting lindhard angle by $\theta_in$ error:  
my angular resolution (~8.9 µrad) is larger than $θ_L/2$ ($~6.5 µ\text{rad}$), so treating "shift the cut edge by one full resolution sigma" as systematic is too aggressive: a more defensible (and numerically stable) choice is to shift by a fraction of the resolution or a smaller fixed step, rather than the full σ  
0.981%
3. shifting the box margins by $d_0$ error along x *(main source)*:  
1.098%  
-  shifting the box margins by $d_0$ error along y (shift of 2.5% of y length):  
0.001%
4. n sigma integration window (changing bin min (2.5σ - 3.5σ) from which i integrate n_ch):  
0.054%

*(heavier)*
1. fit model choice (full quadratic vs parabolic y):  
0.031%
2. nx/ny slice choice:  
0.003%
3. crystal margins cut:
0.077%
4. preliminary cut threshold:  
0.063%

| %        | **total** | binomial| $\theta_0$ | $\theta_L$ | shift $x$ | shift $y$ | bin min | fit model choice | nx ny choice | crystal margins cut | preliminary cut |
|----------|-----------|---------|------------|------------|-----------|-----------|---------|------------------|--------------|---------------------|-----------------|
| **8430** | **1.5**   | 0.063   | 0.004      | 0.981      | 1.098     | 0.001     | 0.054   | 0.031            | 0.003        | 0.077               | 0.063           |
| **8431** | **1.5**   | 0.064   | 0.007      | 0.957      | 1.079     | 0.001     | 0.054   | 0.188            | 0.130        | 0.017               | 0.072           |
| **8650** | **1.3**   | 0.051   | 0.003      | 0.755      | 1.001     | 0.001     | 0.066   | 0.004            | 0.034        | 0.009               | 0.000           |
| **8655** | **1.4**   | 0.062   | 0.009      | 0.737      | 1.187     | 0.002     | 0.076   | 0.016            | 0.199        | 0.180               | 0.036           |
| **8656** | **1.4**   | 0.058   | 0.006      | 0.814      | 1.191     | 0.001     | 0.080   | 0.010            | 0.104        | 0.164               | 0.094           |

## MCS final results

```
for v in 8430 8431 8650 8655 8656; do python compute_edges2.py "$v"; done
```

run [mm] | fitted transition $x$ | fitted transition $y$ | mean position error $d_0$ 
-------- | --------------------- | --------------------- | -------------------------
8430     |  0.09 ± 0.01      |  0.29 ± 0.01      | 0.1055 
8431     |  0.08 ± 0.01      |  0.28 ± 0.01      | 0.1055
<!-- 8430/8431|  0.0900 ± 0.0002      |  0.2617 ± 0.0064      | 0.1055 -->
<!-- 8650     |  0.0068 ± 0.0001      |  0.1375 ± 0.0062      | 0.1055
8655     |  0.0740 ± 0.0005      |  2.6090 ± 0.0784      | 0.1155 
8656     |  0.0731 ± 0.0005      | 2.6705 ± 0.0557       | 0.1155  -->

run [urad] | resolution | **x:** $\sigma_{bsl}$ |  $\sqrt{\sigma_{high}^2​−\sigma_{bsl}^2​}$ | **y:** $\sigma_{bsl}$ | $\theta_{mcs}=\sqrt{\sigma_{high}^2​−\sigma_{bsl}^2​}$ | $\theta_{clamp}=\sqrt{\sigma_{clamp}^2​−\sigma_{bsl}^2​}$
--         | --         | --                    |   --             | --                    | --             | --             
8430       | 8.9       | 12.4 ± 0.1          |  65.1 ± 0.1    | 11.5 ± 0.1          | 62.9 ± 0.1   | 100.5 ± 0.7 
8431       | 8.9       |  12.4 ± 0.1         | 65.1 ± 0.1    | 11.5 ± 0.1          | 62.9 ± 0.1   | 101.1 ± 0.7 
<!-- 8430/8431  | 8.89       | 12.75 ± 0.00          |  78.44 ± 0.05    | 11.80 ± 0.06          | 63.09 ± 0.03   | 99.35 ± 0.49  -->
<!-- 8650       | 9.67       | 14.73 ± 0.01          |   95.02 ± 0.14   | 13.34 ± 0.17          | 73.54 ± 0.06   | 118.80 ± 1.68 
8655       | 9.67       | 14.74 ± 0.01          |  86.82 ± 0.12    | ?                     | 70.50 ± 0.09   | ?
8656       | 9.67       | 14.76 ± 0.01          |  86.66 ± 0.11    | ?                     | 70.13 ± 0.09   | ? -->

Even if the resolution is $\sim 7$ it is still correct to use decimal digits in the results if we are talking about the statistic uncertainty of the fit, since the estimate is accompanied by sufficient statistical data to justify that level of precision for the parameter.  
$\delta\sigma \approx \frac{\sigma}{\sqrt{2N}}$ if $N \approx 10^6$ eventi, for $\sigma \approx 12$ we obtain $\delta\sigma \approx 12 / \sqrt{2 \cdot 10^6} \approx 0.01$.

| **run**         | **$x$ [mm]**    | $x_{min}$        | $x_{max}$       | **$y$ [mm]**     | $y_{clamp}$       | $y_{min}$        | $y_{max}$       |
|-----------------|-----------------|------------------|-----------------|------------------|-------------------|------------------|-----------------|
| **8430**        | 2.11 ± 0.01 | -1.07 ± 0.0003 | 1.05 ± 0.01 | 14.09 ± 0.04 | -11.07 ± 0.01 | -5.43 ± 0.01 | 8.66 ± 0.03 |
| **8431**        | 2.11 ± 0.01 | -1.07 ± 0.0003 | 1.05 ± 0.01 | 14.03 ± 0.04 | -11.08 ± 0.01 | -5.41 ± 0.01 | 8.62 ± 0.03 |
<!-- | **8430 & 8431** | 2.0907 ± 0.0003 | -1.0600 ± 0.0002 | 1.0308 ± 0.0002 | 13.9346 ± 0.0195 | -11.1049 ± 0.0059 | -5.3966 ± 0.0091 | 8.5380 ± 0.0173 |
| **8650**        | 2.0971 ± 0.0001 | -0.2709 ± 0.0001 | 1.8261 ± 0.0000 | 13.7271 ± 0.0156 | -11.8176 ± 0.0094 | -6.3653 ± 0.0156 | 7.3618 ± 0.0015 |
| **8655**        | 2.1200 ± 0.0006 | 0.6031 ± 0.0005  | 2.7232 ± 0.0004 | 12.6059 ± 0.0717 | -10.3721 ± 0.0581 | -6.9081 ± 0.0570 | 5.6978 ± 0.0436 |
| **8656**        | 2.1224 ± 0.0006 | 0.6015 ± 0.0004  | 2.7239 ± 0.0004 | 13.0073 ± 0.1430 | -10.2544 ± 0.0557 | -7.1044 ± 0.0575 | 5.9029 ± 0.1309 | -->



# Simulation
- Cristallo lungo: Week 12, “Simulation for the long crystal” e “Submit the jobs on the cluster” (parti tecniche, valuta se metterle in appendice). Poi Week 14, “Analysis and plots for the long crystal simulations”.
- Applicazione dello smear: Week 15, “MCS in simulations”, e Week 16, “Plotting long crystal final simulation” (confronto con/senza smear).
- Xcoll, corto e lungo: Week 16, “Xcoll”. Per ora è solo un elenco di cose da fare, quindi la sezione andrà scritta quando avrai i risultati.


## BDSIM Simulation

**Geant4** is a Monte Carlo toolkit for simulating particles passing through matter. Given a geometry, a set of materials, and a physics list (which processes are switched on — multiple scattering, ionization, hadronic interactions, etc.), it steps particles through space and simulates their interactions stochastically. `g4FTFP_BERT` is a standard hadronic physics list (Fritiof model at high energy, Bertini cascade at low energy) — general-purpose, no crystal physics.  

**BDSIM** is a Geant4-based application purpose-built for accelerator beamlines: it turns an accelerator lattice description into Geant4 geometry and tracks beam through it, including realistic material interactions in collimators, targets, and (via an extension) bent-crystal channeling. Its input language, *gmad*, looks like MAD-X syntax but describes beamline elements, materials, beam parameters, and physics options.  

**BDSIM-Link** (what both scripts use) is a hybrid-tracking mode: an external code (here, xtrack/xsuite) does the fast, analytic tracking through magnets, and only hands particles to Geant4/BDSIM at the few elements where detailed matter interaction actually matters (a collimator jaw, a target, a crystal). This is much cheaper than tracking the whole ring in Geant4.

in the folder [simulation/bdsim_test](../simulation/bdsim_test/) there is a [demo.py](../simulation/bdsim_test/demo.py) file which is the script:
- tracking protons in a ring accelerator
- with 2 objects: a helium target (generate a halo) and a tungsten jaw collimator (absorb the halo)
- 100 turns, each turn particles are sent through coll1, emittance are recorded, then xtrack propagates the surving particles through the arc to coll2, then back to the start $\to$ surviving particle counts and coordinated at both elements are logged each turn
- needs the file [trackerInterface.gmad](../simulation/bdsim_test/trackerInterface.gmad) for specifics

in the folder [simulation/bdsim_crystal](../simulation/bdsim_crystal/) there is a [demo_crystal.py](../simulation/bdsim_crystal/demo_crystal.py) file which is the script:
- tracking 180 GeV protons in short/long crystal
- needs the file [trackerInterface.gmad](../simulation/bdsim_crystal/trackerInterface.gmad) for specifics

$\to$ ✅ try to run a simulation of ~20000 particles for a short crystal and plot the angular scan (deflection vs incoming angle)

## Study the beam distribution

what i can do, is study the beam distribution from a run of real data, so the impact distributions of $(x,y)$ and $(\theta_x, \theta_y)$, three out of four are perfect gaussian distributions, while the histogram of the position $x$ seems descrivable by a uniform distribution, but is indeed a normal distribution too, just on a wider range. can i simply sample their results and put them in the simulation?  
$\to$ in principel yes, but i need to be careful:  
fitting the 1d marginal distributions and then sampling them indipendently is correct only if the variables are uncorrelated with one another. If, on the other hand, there is a correlation (e.g. $x-\theta_x$, which is very common in a real beam due to the upstream optics), sampling $x$ and $\theta$ from two separate `np.random.normal()` distributions destroys the correlation, even if each marginal distribution individually appears perfect.

I check the covariance and correlation matrices not with all data but within $\pm 5 \sigma$ of the gaussian fit (98% of data), otherwise the widths of the angular distributions would be too high.

|                  | $\mu_x$      | $\sigma_x$  | $\mu_y$      | $\sigma_y$  | $\mu_{\theta,x}$ | $\sigma_{\theta,x}$ | $\mu_{\theta,y}$ | $\sigma_{\theta,y}$ |
|------------------|--------------|-------------|--------------|-------------|------------------|---------------------|------------------|---------------------|
| 8430             | -0.21 ± 0.00 | 2.32 ± 0.00 | 0.60 ± 0.00  | 2.40 ± 0.00 | -0.93 ± 0.01     | 26.97 ± 0.01        | 6.56 ± 0.01      | 41.09 ± 0.01        |
| 8431             | -0.35 ± 0.00 | 2.37 ± 0.00 | 0.66 ± 0.00  | 2.40 ± 0.00 | -1.60 ± 0.01     | 27.08 ± 0.01        | 7.95 ± 0.01      | 41.02 ± 0.01        |
| 8650             | -0.60 ± 0.00 | 2.51 ± 0.00 | 0.55 ± 0.00  | 2.62 ± 0.00 | -2.79 ± 0.01     | 29.82 ± 0.01        | 5.64 ± 0.02      | 44.30 ± 0.01        |
| 8655             | -0.42 ± 0.00 | 2.22 ± 0.00 | -0.02 ± 0.00 | 2.53 ± 0.00 | -1.73 ± 0.01     | 29.63 ± 0.01        | -1.51 ± 0.01     | 43.59 ± 0.01        |
| 8656             | -0.38 ± 0.00 | 2.23 ± 0.00 | -0.04 ± 0.00 | 2.53 ± 0.00 | -1.37 ± 0.01     | 29.66 ± 0.01        | -1.24 ± 0.01     | 43.58 ± 0.01        |
| weighted average | -0.37867394  | 2.01454317  | 0.3615005    | 2.80776263  | -1.62272038      | 30.89146025         | 3.64729974       | 45.734049           |

|**correct** *8430 Correlation matrix after cut:*  ||||
|--|--|--|--|
|[ 1.         |  0.01347898 |  0.51028308 | 0.02852061]|
|[ 0.01347898 |  1.         | -0.0796839  | 0.23902702]|
|[ 0.51028308 | -0.0796839  |  1.         | 0.05384073]|
|[ 0.02852061 |  0.23902702 |  0.05384073 | 1.        ]|

|*8430 Correlation matrix before cut:*  ||||
|--|--|--|--|
|[ 1.         | 0.00249682  | 0.56065383  | 0.04752554]  |
|[ 0.00249682 |  1.         | -0.04600791 |  0.34909413] |
|[ 0.56065383 | -0.04600791 |  1.         |  0.07679966] |
|[ 0.04752554 |  0.34909413 |  0.07679966 |  1.        ] |


|**correct** *8430 Covariance matrix after cut:*||||
|--|--|--|--|
|[ 3.82592775e+00 |  7.08532402e-02 |  2.94087313e+01 | 2.44291188e+00]|
|[ 7.08532402e-02 |  7.22218856e+00 | -6.30959887e+00 | 2.81295028e+01]|
|[ 2.94087313e+01 | -6.30959887e+00 |  8.68147418e+02 | 6.94685664e+01]|
|[ 2.44291188e+00 |  2.81295028e+01 |  6.94685664e+01 | 1.91761445e+03]|

|*8430 Covariance matrix before cut:*||||
|--|--|--|--|
|[   4.74346824 |   0.01606491 |   89.41765855 |    7.18923536]|
|[   0.01606491 |   8.72738918 |   -9.9530266  |   71.62957414]|
|[  89.41765855 |  -9.9530266  | 5362.42227158 |  390.61366706]|
|[   7.18923536 |  71.62957414 |  390.61366706 | 4824.08887831]|


- $x-y$: correlation 0.0025 $\to$ negligible
- $x-\theta_x$: correlation 0.56 $\to$ NOT negligible
- $y-\theta_y$: correlation 0.35 $\to$ NOT negligible
- $\theta_x-\theta_y$: correlation 0.076 $\to$ small - can be statistic rumor

It is correct to have positive correlation between the two, a greater angle corresponds to a greater incoming position.

## Study the beam distribution (part 2)

how to consider the correlation?
1. Paired resampling from real data (the simplest and most robust method; it does not require assuming a functional form for the correlation):  
```python
idx = np.random.choice(len(x_meas), npart, replace=True)
x_impact  = x_meas[idx]
thetax_impact = thetax[idx]
# y e thetay possono restare fit+sampling indipendente, o anch'essi paired se preferisci uniformità di metodo
y_impact  = y_meas[idx]
thetay_impact = thetay[idx]
```
By using the same idx for all four, you automatically preserve all the actual correlations (even the small ones), without having to decide on a case-by-case basis which ones to model.  

2. **actually used**: I combined all 5 covariance matrices (correlation and covariance matrix obtained from each run), into a single one by doing a weighted average - weighted by their number of entries $w=n/n_{tot}$ - and the same for the gaussian means and sigmas. then by doing a multivariate random distribution I get my distribution.  
`rng = np.random.default_rng(seed=82)`  
`samples = rng.multivariate_normal(mean=mu_avg, cov=cov_avg, size=npart)`  

## Study the beam distribution (part 3)
- it is correct to have correlation actually:  
In beam optics, position and angle at a given location are correlated whenever the beam isn't exactly at a **waist** (a point where the beta function is at a local minimum and the Twiss `alpha` parameter is zero). The phase-space ellipse (x vs. px) has a tilt described by `alpha`: $<x·p_x> = -\alpha * emittance$.  
If $\alpha ≠ 0$ at the crystal location (which is the general case — you'd need to specifically design the optics to have a waist exactly there), then x and px are correlated by construction. This isn't an artifact or noise — it's a direct consequence of how the beam was transported (focusing, drifts) upstream of the measurement point. A correlation of ~0.15 is modest but completely plausible for a real beamline where the crystal isn't sitting exactly at a waist.  
Since $\alpha ≠ 0$ is now expected, you could cross-check your measured correlation against the beamline's actual Twiss parameters at the crystal location, if you have access to the optics model (MAD-X/xsuite lattice for that beamline). The theoretical relationship is: $corr(x, p_x) = -\alpha / \sqrt{1 + \alpha²}$.  
If you know $\alpha$ at the crystal from the lattice design, you can compute the expected correlation and compare it to your measured ~0.15 as a sanity check — if they roughly agree, that's strong confirmation the correlation is real beam optics and not, say, a reconstruction/tracking artifact in your telescope analysis.

- about the spikes check if rebinning it changes, if it shows also on the "all" particles and if the remaining particles - the "not entered" - are compatible as the others

## Channeling Analysis on the simulated short crystal

(i am currently not saving the variables `delta` - momentum divergence, `zeta` - longitudinal coordinate, `ptau` - energy divergence, `pdg_id` and `weight`)

- no single track selection or detector alignment checks needed.
- selecting the survived particles (`state_out == 1`)(also `state_in == 1` needed? just confirms a slot is a real, active particle at the start) we are selecting the particles that survived without being absorbed (made it through).s
- what is the aperture? it describes the surrounding vacuum-pipe/housing aperture - the mechanical structure around the crystal. think about changing `cry1.aper1` which is 50 mm meters (was 100 meters before!! so a lot of particles were not passing though the crystal (survived with a $\Delta\theta=0$)): keeping it greater than the crystal's dimension is also keeping the `state_out` filter physically meaningful (decoupling what survived the aperture and what hit the crystal)
- *spatial cut* (selection of the ones that actually enter the crystal) on the impact position needed: $|x_{in}| < x_{size}/2$ and $|y_{in}| < y_{size}/2$ .
- *Lindhard cut:* $|\theta_{in}| < \theta_L/2$ so selecting the particles with the potential to be channeled.
- now we can compute the channeling efficiency on the "subset" (= went inside the crystal and entered within the angular acceptance)

**Results:**  
with [18 sep run](../simulation/bdsim_crystal/cry1_20260918_10.npz) of 200.000 particles and through [plot_simulation.py](../simulation/bdsim_crystal/plot_simulation.py) we obtain:  
- channeling efficiency: $\epsilon_{ch}=26190/32201=81.3\%$ if area below, or $82.4\%$ if within $\pm 3 \sigma$
- channeling peak: $\theta_b=(50.04\pm0.09)\,\mu\text{rad}$ and $\sigma=(7.12\pm0.09)\,\mu\text{rad}$

## Simulation results - short crystal

Out of 200k simulated particles for the short crystal:
- survived = all
    - what gets the particles absorbed, is the inelastic nuclear interaction, that goes like $1-e^{-L/\lambda}$, where $\lambda$ is the length for the nuclear interaction in silicon ($\sim 46.2\,\text{cm}$):
    - short crystal ($L=4\,\text{mm}$): $(1-e^{-4/462})=0.9\%$
    - long crystal ($L=74\,\text{mm}$): $(1-e^{-74/462})=15\%$
- within critical angle ~ 120  
    - 120/650 ~ 18% of survivors land within the lindhard window, and that is what i expect from my beam divergence: $\sigma(\theta_x)$ is around $30 \,\mu\text{rad}$ and $\theta_L/2\sim6.6\,\mu\text{rad}$, for a gaussian of that width, the fraction falling inside $\pm\theta_L/2$ works out about 18-19% $\to$ so this part of the pipeline is behaving exactly as physics predicts, my beam is simply more divergent than the crystal's angular acceptance.
    - for a gaussian angular distribution with width $\sigma_\theta$, the fraction of particles falling inside a symmetric window $\pm a$ around the mean is $f=erf(\frac{a}{\sigma_\theta\sqrt{2}})$. in our case $\pm a=\pm\theta_L/2$ and using $\sigma_\theta\sim27\,\mu\text{rad}$ we get $f\sim19.4%$


## Simulation for the long crystal
Through `merge_npz.py` I merge all the runs of the three condors:  
- condor: 242 runs $\to$ 2'400'400 particles
- condor2: 1927 runs $\to$ 19'270'000 particles
- condor3: 1600 runs $\to$ 16'000'000 particles

**Key points for the analysis:**
1. filter survived particles: select `state_out == 1`
2. geometrical cut: i know for sure from the simulation script that the crystal is centered around zero so i can cut like a box - $(|x_{in}|<cry_x/2, |y_{in}|<cry_y/2)$
3. critical angle selection: within $\pm 1/2 \,\theta_L$
4. channeling analysis: gaus fit around the channeled peak and integral - remember to include in $N_{tot}$ also particles with the potential to be channeled that did not survive
- multiple coulomb scattering
- volume reflection fraction: near $\Delta\theta \approx -θ_L$ $\to$ $14\%$
- dechanneling: histo of $\Delta\theta_x$ between 0 and the peak, on the attemped population
- energy lost in the crystal

**Results:**  
||survived|geom hit|$\theta_b$ [urad]|$\sigma_\theta$ [urad]|$\epsilon_{eff}$|
|--|--|--|--|--|--|
|1|90.16%|99.1%|5999 ± 1|6.7 ± 0.1|(35.1± 0.1)%|
|2|90.15%|99.2%|5999 ± 1|6.7 ± 0.1|(34.9± 0.1)%|
|3|90.15%|99.2%|5999 ± 1|6.7 ± 0.1|(34.9± 0.1)%|

## MCS in simulations

MCS in the simulation for Silicon is already taken into account - also nuclear interactions (for dechanneled and volume reflected particles - while for channeled particles their effect is given by the channeling model).  
What I have to add manually is the **tracker uncertainty contribution**.

## How to apply the smear
$\to$ doubt: isn't the effect already taken into account in the beam distribution i got from real data?  
Yes, but still the single values of the angles in the simulation does not have a error value, so the addition is necessary.  
We can subtract the $\sigma_{track}$ from the incoming beam distribution (but we would need to re-run the simulation) and then reapply the error, if we don't do it, it shouldn't be a big problem: the contribution of the effect is:
- $θ_x: \sqrt{(27² − 9²)} ≈ 25.5 µrad$ vs $27$, $\sim 5\%$ less
- $θ_y: \sqrt{(40² − 9²)} ≈ 39 µrad$ vs $40$, $\sim 2.5\%$ less

If we want to be correct we just subtract $9^2 \,\mu\text{rad}^2$ from angular diagonal of covariance matrix before the beam gets generated:  
```python
cov_true = cov_avg.copy()
cov_true[2, 2] -= (9e-6)**2     # theta_x
cov_true[3, 3] -= (9e-6)**2     # theta_y
```

The measures are independent.  
$\theta_{in} = \theta_{in} + N(0,\sigma_{in})$ with $\sigma_{in}^2=\sigma_{track}^2+\sigma^2_{air, up}$  
$\to$ so for each particle I extract a casual number from a gaussian distribution (`rng.normal(0, s_in, theta_x_in.size)`)  
$\theta_{out} = \theta_{out} + N(0,\sigma_{out})$ with $\sigma_{out}^2=\sigma_{track}^2+\sigma^2_{air, down}$  
$\Delta\theta = \theta_{out}^{new} - \theta_{in}^{new}$ so $\sigma(\Delta\theta)\approx\sqrt{\sigma_{in}^2+\sigma_{out}^2}$  


After, the Lindhard cut will be done on $\theta_{in}^{new}$, since $\sigma_{track}$ is compatible with $\theta_L$.

About the air contribution:  
the factor $\sqrt{x/X_0}$ sums exactly in quadrature (variances become linear in $x/X_0$), but the logarithmic correction does not. For two equal segments of length $x$, adding the two $θ_0$ values in quadrature instead of using a single $θ_0$ with $2x$ changes the logarithmic factor by $0.038\cdot\ln 2 ≈ 2.6\%$. 

**using:**  
- $\theta_{in} = \theta_{out} = 9\,\mu\text{rad}$
- $L_{air}=0$
- seed $=83$

## Plotting long crystal final simulation
$\to$ [cry2_merged5.npz](../simulation/bdsim_longcrystal/runs/cry2_merged5.npz) [6 GB!!]

Through [plot_simulation.py](../simulation/bdsim_longcrystal/plot_simulation.py):  
- incoming beam features
- added smear taken from data $\sigma_{tracker}=9\,\mu\text{rad}$
- compute efficiency (taking into account in the denominator also the ones that could have channeled but were absorbed/lost)
- both dataset are kept: comparison with/without smear

**RESULTS**:
||without ||| with $\sigma=9$|||
|-|-|-|-|-|-|-|
|**particles**|$\varepsilon_{ch}$|**peak**|sigma |$\varepsilon_{ch}$|**peak**|sigma|
|19740000|(34.4± 0.0)%|5999.15 ± 0.03 urad|6.75 ± 0.03|(21.8± 0.0)%|5999.18 ± 0.02 urad |11.19 ± 0.02|

## Xcoll

Simulation through XSuite, with the XColl routine, designed for the short collimation crystals.  
Here I can find an example: [https://github.com/xsuite/xcoll/blob/main/examples/everest_crystal.py](https://github.com/xsuite/xcoll/blob/main/examples/everest_crystal.py), which i can adapt for my use.

- simulate short crystal
- simulate long crystal
- my beam distribution
- adding the smear
