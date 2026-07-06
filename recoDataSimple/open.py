import ROOT
import sys

f = ROOT.TFile.Open("recoDataSimple_8430_xtalMerging.root")

if not f or f.IsZombie():
    print("Error: file has not opened correctly.")
    sys.exit(1)

tree = f.Get("simpleEvent")

if not tree:
    print("Error: tree not found.")
    sys.exit(2)

events = tree.GetEntries()
print(f"File has opened correctly!\nTree contains {events} events.")
print("="*50)

list_branches = tree.GetListOfBranches()

for branch in list_branches:
    name_branch = branch.GetName()
    print(f"Branch: {name_branch}")

    list_leaves = branch.GetListOfLeaves()
    for leaf in list_leaves:
        name_leaf = leaf.GetName()
        type_leaf = leaf.GetTypeName()
        print(f"     Leaf: {name_leaf} - type {type_leaf}")
    
    print('-'*40)
