import ROOT

df = ROOT.RDataFrame("simpleEvent", "recoDataSimple_8430_xtalMerging.root")

columns = df.GetColumnNames()
for col in columns:
    print(col)

# crea una visualizzazione in stile tabella delle prime 5 righe
df.Display(["GonioPos.x", "GonioPos.y", "GonioPos.z"], 5).Print()

# grafico/istogramma di una variabile
#istogramma = df.Histo1D("variabile1")
#istogramma.Draw() # Lo disegna a schermo
