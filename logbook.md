# Traineeship al CERN

1. [week 1](#week-1)
2. [week 2](#week-2)  
    - [measuring the channeling efficiency](#measure-the-channeling-efficiency)
3. [week 3](#week-3)
4. [week 4](#week-4)

## WEEK 1  

books/articles to read: 
- [_Crystal channeling and its application at high-energy accelerators_](./aboutcrystals/crystal%20channelling.pdf)  by Biryukov, Chesnokov, Kotov
- [_Performance of short and long bent crystals for the TWOCRYST experiment at the Large Hadron Collider_](./aboutcrystals/s10052-025-15092-y.pdf) in The European Physics Journal (May 2025)
- [_New direction for bent crystals_](./aboutcrystals/New%20directions%20for%20bent%20crystals%20–%20CERN%20Courier.pdf) in Cern Courier by P. Hermes, S. Redaelli (March 2026)  

$\to$ summary in [books.md](./aboutcrystals/books.md)

data to look at:  
[recoDataSimple_8430_xtalMerging.root](https://cernbox.cern.ch/files/spaces/eos/user/p/pahermes/018_crystals/260701_Anna_Background_Infos/recoDataSimple_8430_xtalMerging.root):
- 1 tree (`simpleEvent`)  
- 9 branches  
    - `Time`
    - `Date`
    - `Event` (run/I, evtnum/I, nuclear/I, nuclearRaw/I)
    - `GonioPos` (x/D, y/D, z/D)
    - `MultiHits` (p_nHits[6]/I, thetaIn_x/D, thetaIn_y/D, thetaInErr_x/D, thetaInErr_y/D, d0_x/D, d0_y/D, d0Err_x/D, d0Err_y/D)
    - `Tracks` (thetaIn_x/D, thetaIn_y/D, thetaOut_x/D, thetaOut_y/D, thetaInErr_x/D, thetaInErr_y/D, thetaOutErr_x/D, thetaOutErr_y/D, d0_x/D, d0_y/D, d0Err_x/D, d0Err_y/D, d0Out_x/D, d0Out_y/D, d0OutErr_x/D, d0OutErr_y/D, chi2_x/D, chi2_y/D, chi2Out_x/D, chi2Out_y/D)
    - `MultiTracks` (d0In1_x/D, d0Out1_x/D, d0In2_x/D, d0Out2_x/D, d0In3_x/D, d0Out3_x/D, thetaIn1_x/D, thetaOut1_x/D, thetaIn2_x/D, thetaOut2_x/D, thetaIn3_x/D, thetaOut3_x/D, d0In1_y/D, d0Out1_y/D, d0In2_y/D, d0Out2_y/D, d0In3_y/D, d0Out3_y/D, thetaIn1_y/D, thetaOut1_y/D, thetaIn2_y/D, thetaOut2_y/D, thetaIn3_y/D, thetaOut3_y/D)
    - `SingleTrack` (singleTrackEvent/I)
    - `MultiHit` event (multiHitEvent/I)
- time and date are `Char_t` so characters $\to$ to see them write on root terminal `root[1] simpleEvent->Scan("Time:Date")`
- some leaves are `TLeafI` so integers, others `TLeafD` so double

$\to$ code to open them [open.py](./recoDataSimple/open.py) or [open_as_df.py](./recoDataSimple/open_as_df.py)



## WEEK 2

data:
| long crystal 1 | data   | time   | long crystal 2 | data   | time   |
|----------------|--------|--------|----------------|--------|--------|
| 8430           | 151024 | 184745 | 8655           | 281024 | 012716 |
| 8431           | 151024 | 201923 | 8656           | 281024 | 021740 |
| 8430 & 8431    | 151024 |        |                |        |        |
| 8650           | 271024 | 181615 |                |        |        |

what type of data analysis can i do?  
- isolate clean events (`singleTrackEvent == 1`) which prevent from considering events from "pile-up effect", and check $\chi^2$ (`chi2_x`, `chi2_y`)
- beam divergence profile: how is the beam arriving to the detector: `thetaIn_y` vs `thetaIn_x`, `d0_y` vs `d0_x` (?)
- $\Delta\theta_x=\theta_{out}-\theta_{in}$ of `singleTrackEvent` (where are the peaks? i expect 0 and some $\mu$ rad)
- crystal acceptance $\Delta\theta_x$ vs $\theta_{in,x}$  

from [checkdata_1.py](./recoDataSimple/checkdata_1.py):
![alt text](image.png)

what to do next:
- ✅ check through [checkdata_0.py](./recoDataSimple/checkdata_0.py) that within the same run the position of the goniometer ($x,y,z$) does not change
- ✅ check about `SingleTrack`, `MultiHit`: they are either 0 or 1 and they always match in every run
- ✅ look at the angle values from the [articles](./aboutcrystals/s10052-025-15092-y.pdf)  
    | | | TCCP | TCCPA |
    |-|-|-----|-------|
    |length|[mm]|70|70.5|
    |width|[mm]|8|22.5|
    |height|[mm]|2|2|
    |bend radius $\rho$ | [m]|10|5.3|
    |bend angle $\theta_b$ | [mrad]|7.0|13.3|
    |$\theta_L$ @ 180 GeV/c | [$\mu$ rad]|12.9|12.5|  
$\theta_L = (1-\frac{\rho_c}{\rho})\sqrt{\frac{2U_0}{E}}$ where: $\rho_c=\frac{E}{U'(x_c)}$, $U'(x_c)=5.7$ GeV/cm, $U_0=16$ eV for Silicon (110) 

### Measure the channeling efficiency $\epsilon_{ch}=\frac{N_{ch}}{N_{tot}}\cdot 100$

- center the angle and choosing $\theta_0$ for torsion correction $\tau_y$, which is the variation of the crystalline plane orientation along the direction ($y$) perpendicular to the bending plane ($x$), and must be accounted for (the different orientations of the crystallographic planes along the crystal front surface influence the angular range in which particles can be channelled, i.e., the cut at $\pm 1/2 \,\theta_L$ may not be centred on zero, but a correction parameter $\theta_0 = \theta_0 (y)$)
    - ✅ in [checkdata_2.py](./recoDataSimple/checkdata_2.py) we treat the crystal as if it was ideal, and we find a single global $\theta_{optimal}$ from the peak of the angular scan: we find the function of $\theta_0(x)$ and find the peak - the one that maximize the efficiency.
    - ❌ what the articles does is that it takes the beam “spot” (`d0_x` vs `d0_y`) and divides it into a checkerboard grid. For each square on the grid, it performs an angular scan (like my previous h_scan) and finds the specific peak angle $\theta_0$ for that square. It then reapplies a correction to each particle by subtracting the local $\theta_0$
- $N_{tot}:$ filter the data set, choosing the events that produced a single track, have a good $\chi^2$ value (?), have actually entered the crystal, and have an angle within $\pm\frac{1}{2}\theta_L$ (actually we choose them $|\theta_{in,x}-\theta_0|\leq\frac{1}{2}\theta_L$):
    - to filter spatially those who entered the crystal (of dimension width $\times$ length), i look at `d0` variables $x$ and $y$ (`d0` and `d0Out`), total and from those events that will present some channeling (so $\Delta\theta_x>\sim 5850 \mu$ rad, number found as $\mu-4\sigma$ by preliminary fit) and from these last ones filtered `d0Out`, which form an almost perfect rectangular, i infere the position of the crystal with respect to the beam through a sliding window method (which maximize the integral).
    - to filter the particles that have entered with the right angle, we first need to find $\theta_0$ in order to apply the filter $|\theta_{in,x}-\theta_0|\leq\frac{1}{2}\theta_L$:
- identify $N_{ch}$: on the $\Delta\theta_x$ filtered histogram perform a gaussian fit on the channeling peak, only from the right side of the peak, then do the integral around the peak
- compute $\epsilon_{ch}=\frac{N_{ch}}{N_{tot}}\cdot 100$ amd its error - we choose a binomial distribution

$\to$ [slides week 2](./slides/week2.pdf)

## WEEK 3

to do:
- checkerboard grid search (1 mm $\times$ 1 mm) of torsion correction (see [week 2 analysis](#measure-the-channeling-efficiency)) $\to$ [checkdata_3.py](./recoDataSimple/checkdata_3.py):
    - cut will be $|\theta_{in,x}-(\theta_{0,center}+\tau_y\cdot d0_y)|\leq\frac{1}{2}\theta_L$
    - plotting $\theta_0$ vs $y$, we should see a linear correlation, and thus find $\tau_y$ [$\mu$ rad/mm]  
    $\to$ do a 2D torsion map ($x$, $y$, and angle shift). these values are then fitted ussing a linear interpolator along $x$ and $y$ in the entire crystal surface to extract the continuous distribution of the map, we will obtain an average value for the torsion on the $y$ direction, while the torsion along the $x$ direction is negligible  
        $\to$ we expect $\tau_y(\text{TCCP})> 5\tau_y(\text{TCCPA})$
- computing the final $\epsilon_{ch}$, after applying an angular shift to the incoming direction of the particles, accordingly to the value of the torsion map at thier impact position on the crystal surface  
    $\to$ we obtain a global efficiency curve
- improve $\Delta x$ spatial shift: if $\Delta z$ is the distance between the two detectors, the shift should be $\Delta x=\Delta z \cdot\tan(\theta_{defl})\sim\Delta z\cdot\theta_{defl}$
- remove bkg?:
    - removing instrumental noise, volume reflection, multiple scattering bkg from $N_{tot}$(?) by using the particles that didn't hit the crystal as gauge, plot their $\Delta\theta$ distribution, normalize it and subtract it (so that its height matches the zero-peak of my data) from my actual data
    - removing dechanneling bkg (gradually decreasing) from $N_{ch}$: can be modelled by a crystal ball function (power law tail) and subtracted it from the $N_{ch}$ integral  
    $\to$ the article says that since the two peaks are very well separated, the dechanneled background is negligible, because it is spread over a much broader angular range
- found a systematical uncertainty in $\epsilon_{ch}$ and in $N_{ch}$
- analyse all the datasets

