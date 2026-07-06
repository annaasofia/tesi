import ROOT
import sys

f = ROOT.TFile.Open("recoDataSimple_8430_xtalMerging.root")

if not f or f.IsZombie():
    print("Error: file has not opened correctly.")
    sys.exit(1)
else:
    print("File has opened correctly!")

tree = f.Get("simpleEvent")

if not tree:
    print("Error: tree not found.")
    sys.exit(2)


