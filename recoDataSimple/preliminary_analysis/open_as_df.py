import ROOT
import sys


try:
    df = ROOT.RDataFrame("simpleEvent", f"../data/recoDataSimple_{sys.argv[1]}_xtalMerging.root")
    count = df.Count().GetValue()
    print(f"Dataset has {count} events.")
    df_clean = df.Filter("SingleTrack == 1")
    # print(f"{df.Count().GetValue() - df_clean.Count().GetValue()} events ({(df.Count().GetValue() - df_clean.Count().GetValue())/(df.Count().GetValue())*100}%) out of {df.Count().GetValue()} got discarded.")

except Exception as e:
    raise Exception("Error: file has not opened correctly.") from e


# columns = df.GetColumnNames()
# for col in columns:
#     print(col)

# tabella delle prime tot righe
# df.Display(["GonioPos.x", "GonioPos.y", "GonioPos.z"], 10).Print()
# df.Display(["Event.run", "Event.evtnum", "Event.nuclear", "Event.nuclearRaw"], 15).Print()
# df.Display(["Event.run", "Event.evtnum", "SingleTrack", "MultiHit", "MultiHits.p_nHits"], 25).Print()
df = df.Define("thetaIn_x", "Tracks.thetaIn_x * 1e6").Define("thetaOut_x", "Tracks.thetaOut_x * 1e6")\
    .Define("thetaIn_y", "Tracks.thetaIn_y * 1e6").Define("thetaOut_y", "Tracks.thetaOut_y * 1e6")\
    .Define("thetaInErr_x", "Tracks.thetaInErr_x * 1e6")\
    .Define("thetaInErr_y", "Tracks.thetaInErr_y * 1e6")\
    .Define("thetaOutErr_x", "Tracks.thetaOutErr_x * 1e6")\
    .Define("thetaOutErr_y", "Tracks.thetaOutErr_y * 1e6")\
    .Define("Deltatheta_x", "(Tracks.thetaOut_x - Tracks.thetaIn_x) * 1e6")\
    .Define("Deltatheta_y", "(Tracks.thetaOut_y - Tracks.thetaIn_y) * 1e6")\
    .Define("DeltathetaErr_x", "sqrt(Tracks.thetaInErr_x * Tracks.thetaInErr_x + Tracks.thetaOutErr_x * Tracks.thetaOutErr_x) * 1e6")\
    .Define("DeltathetaErr_y", "sqrt(Tracks.thetaInErr_y * Tracks.thetaInErr_y + Tracks.thetaOutErr_y * Tracks.thetaOutErr_y) * 1e6")

# df.Display(["thetaIn_x", "thetaInErr_x", "thetaOut_y", "thetaOutErr_y"], 20).Print()
# df.Display(["Deltatheta_x", "DeltathetaErr_x", "Deltatheta_y", "DeltathetaErr_y"], 20).Print()
# df.Display(["Tracks.d0_x", "Tracks.d0Err_x", "Tracks.d0_y", "Tracks.d0Err_y"], 20).Print()
# df.Display(["Tracks.d0Out_x", "Tracks.d0OutErr_x", "Tracks.d0Out_y", "Tracks.d0OutErr_y"], 20).Print()


# istogramma di una variabile
#histo = df.Histo1D("Tracks.d0_y")
#histo.Draw()

# istogramma di due variabili
#histo2d = df.Histo2D(('map','map',1000,-10,-10,1000,-10,-10),"Tracks.d0_x","Tracks.d0_y")
#histo2d.Draw("COLZ")


n_all = df.Count().GetValue()
n_single = df.Filter("SingleTrack == 1").Count().GetValue()
print(f"{n_all} events, {n_all - n_single} multi-hit/multi-track ({100*(1-n_single/n_all):.2f}%)")


#df.Get("simpleEvent").Scan("Time:Date", "", "", 10)
# df.Display(["Time","Date"],10).Print()
