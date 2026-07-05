import ROOT
import sys

f = ROOT.TFile.Open("root://eosuser.cern.ch//eos/user/p/pahermes/018_crystals/260701_Anna_Background_Infos/recoDataSimple_8430_xtalMerging.root")
#df = ROOT.RDataFrame("simpleEvent", "root://eosuser.cern.ch//eos/user/p/pahermes/018_crystals/260701_Anna_Background_Infos/recoDataSimple_8430_xtalMerging.root")

if not f or f.IsZombie():
    print("Error: file has not opened correctly.")
    sys.exit(1)

tree = f.Get("simpleEvent")

if not tree:
    print("Error: tree not found.")
    sys.exit(2)


#columns = df.GetColumnNames()
#for col in columns:
#    print(col)

# Crea una visualizzazione in stile tabella delle prime 5 righe
#df.Display(["variabile1", "variabile2"], 5).Print()

# 3. Vuoi fare un grafico/istogramma di una variabile al volo?
#istogramma = df.Histo1D("variabile1")
#istogramma.Draw() # Lo disegna a schermo
