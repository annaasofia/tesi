import ROOT
import sys


try:
    df = ROOT.RDataFrame("simpleEvent", "recoDataSimple_8430_xtalMerging.root")
    count = df.Count().GetValue()
    print(f"Dataset has {count/10000000}*e7 events.")
    df_clean = df.Filter("SingleTrack == 1")
    print(f"{df.Count().GetValue() - df_clean.Count().GetValue()} events ({(df.Count().GetValue() - df_clean.Count().GetValue())/(df.Count().GetValue())*100}%) out of {df.Count().GetValue()} got discarded.")

except Exception as e:
    raise Exception("Error: file has not opened correctly.") from e


# columns = df.GetColumnNames()
# for col in columns:
#     print(col)

# tabella delle prime 5 righe
#df.Display(["GonioPos.x", "GonioPos.y", "GonioPos.z"], 5).Print()
#df.Display(["Event.run", "Event.evtnum", "Event.nuclear", "Event.nuclearRaw"], 15).Print()
# df.Display(["Event.run", "Event.evtnum", "SingleTrack", "MultiHit", "MultiHits.p_nHits"], 25).Print()

# istogramma di una variabile
#histo = df.Histo1D("Tracks.d0_y")
#histo.Draw()

# istogramma di due variabili
#histo2d = df.Histo2D(('map','map',1000,-10,-10,1000,-10,-10),"Tracks.d0_x","Tracks.d0_y")
#histo2d.Draw("COLZ")
