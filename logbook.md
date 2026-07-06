# Traineeship al CERN

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
- look at the angle values from the [articles](./aboutcrystals/s10052-025-15092-y.pdf)
- measure the efficiency $\epsilon_{ch}$ doing a gaussian fit only from the right side of the peak

