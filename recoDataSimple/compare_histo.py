import ROOT


files = []
h = []
runs = [8430, 8650]

for run in runs:
    filename = f"final_histo_run_{run}.root"
    current_file = ROOT.TFile(filename, "READ")
    files.append(current_file)

    histoname = f"h_defl_run_{run}"
    current_hist = current_file.Get(histoname)
    if not current_hist:
        print(f"ERROR: histogram {histoname} not found.")
        continue

    current_hist.SetDirectory(0)
    h.append(current_hist)

ROOT.gStyle.SetOptStat(0)

colors = [ROOT.kBlue, ROOT.kRed]
labels = ["Run 8430", "Run 8650"]

for i in range(len(h)):
    # normalize the area
    h[i].Rebin(5)
    if h[i].Integral() > 0:
        h[i].Scale(1.0 / h[i].Integral())

    
    h[i].SetLineColor(colors[i])
    h[i].SetFillColor(0)
    h[i].SetLineWidth(2)

canvas = ROOT.TCanvas("canvas", "Confronto Channeling", 900, 700)

h[0].SetTitle("Comparison channeling peaks; #Delta#theta_{x} [#murad]; Events")

# h[0].GetXaxis().SetRangeUser(5900, 6250)

h[0].Draw("HIST")
h[1].Draw("HIST SAME")
leg = ROOT.TLegend(0.75, 0.75, 0.88, 0.88)
leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.035)
for i in range(len(h)):
    leg.AddEntry(h[i], labels[i], "l")
leg.Draw("SAME")

canvas.Update()

ROOT.SetOwnership(canvas, False)
ROOT.SetOwnership(leg, False)
