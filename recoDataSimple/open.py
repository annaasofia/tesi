import ROOT

file_remoto = ROOT.TFile.Open("root://eosuser.cern.ch//eos/user/p/pahermes/018_crystals/260701_Anna_Background_Infos/recoDataSimple_8430_xtalMerging.root")
# https://cernbox.cern.ch/files/spaces/eos/user/p/pahermes/018_crystals/260701_Anna_Background_Infos/recoDataSimple_8430_xtalMerging.root

#df = ROOT.RDataFrame("nome_del_tree", "root://eosuser.cern.ch//eos/user/t/tuousername/data/miofile.root")

# Esempio: fai un istogramma di una variabile
#h = df.Histo1D("mia_variabile")
#h.Draw()

if not file_remoto or file_remoto.IsZombie():
    print("Errore: impossibile accedere al file remoto.")
#else:
    # Prendi il tree e analizza
    #tree = file_remoto.Get("nome_del_tree")
    #print(f"Il tree contiene {tree.GetEntries()} eventi.")
