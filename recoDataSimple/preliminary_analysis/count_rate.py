import sys
import numpy as np
import matplotlib.pyplot as plt
import ROOT

from rate_stats import spill_statistics, print_report

plt.rcParams.update({
    "font.size": 16,
    "axes.titlesize": 18,
    "axes.labelsize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,
})
 
run = sys.argv[1] if len(sys.argv) > 1 else "8430"
file_path = f"../data/recoDataSimple_{run}_xtalMerging.root"

print(f"Apertura del file {file_path}...")

# Funzione C++ compilata al volo: gestisce il TTree a basso livello saltando i limiti di RDataFrame
ROOT.gInterpreter.Declare("""
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>

struct RateResult {
    long t_start;
    long t_stop;
    long duration_sec;
    std::vector<long long> counts;
};

long parse_time_raw(const char* t) {
    if (!t) return -1;
    std::string digits = "";
    for (int i = 0; t[i] != '\\0' && digits.size() < 6; ++i) {
        if (isdigit(t[i])) digits += t[i];
    }
    if (digits.size() < 6) return -1;
    int hh = std::stoi(digits.substr(0, 2));
    int mm = std::stoi(digits.substr(2, 2));
    int ss = std::stoi(digits.substr(4, 2));
    return hh * 3600 + mm * 60 + ss;
}

RateResult analyze_tree(const char* filepath, const char* treename) {
    RateResult res{0, 0, 0, {}};
    TFile* f = TFile::Open(filepath, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Errore apertura file ROOT!" << std::endl;
        return res;
    }

    TTree* tree = (TTree*)f->Get(treename);
    if (!tree) {
        std::cerr << "Errore: TTree non trovato!" << std::endl;
        f->Close();
        return res;
    }

    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("Time", 1);

    Char_t time_buf[80];
    tree->SetBranchAddress("Time", time_buf);

    Long64_t nentries = tree->GetEntries();
    if (nentries == 0) {
        f->Close();
        return res;
    }

    tree->GetEntry(0);
    long t_start = parse_time_raw(time_buf);

    tree->GetEntry(nentries - 1);
    long t_stop = parse_time_raw(time_buf);

    if (t_stop < t_start) {
        t_stop += 24 * 3600;
    }

    long duration_sec = t_stop - t_start + 1;
    std::vector<long long> counts(duration_sec, 0);

    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        long t_cur = parse_time_raw(time_buf);
        if (t_cur < 0) continue;
        if (t_cur < t_start) t_cur += 24 * 3600;

        long rel_sec = t_cur - t_start;
        if (rel_sec >= 0 && rel_sec < duration_sec) {
            counts[rel_sec]++;
        }
    }

    res.t_start = t_start;
    res.t_stop = t_stop;
    res.duration_sec = duration_sec;
    res.counts = counts;

    f->Close();
    return res;
}
""")

res = ROOT.analyze_tree(file_path, "simpleEvent")

duration_sec = res.duration_sec
if duration_sec <= 0:
    print("Errore durante l'elaborazione del TTree.")
    sys.exit(1)

counts_per_sec = np.array(res.counts, dtype=np.int64)
dur_min = duration_sec // 60
dur_rem_sec = duration_sec % 60

t_start_day = res.t_start % 86400
t_stop_day = res.t_stop % 86400

print("=" * 50)
print(f"Inizio run: {t_start_day//3600:02d}:{(t_start_day%3600)//60:02d}:{t_start_day%60:02d} | "
      f"Fine run: {t_stop_day//3600:02d}:{(t_stop_day%3600)//60:02d}:{t_stop_day%60:02d}")
print(f"Durata totale: {duration_sec} s ({dur_min} min {dur_rem_sec} s)")
print("=" * 50)

empty_mask = counts_per_sec == 0
empty_seconds = int(np.sum(empty_mask))
pct_empty = (empty_seconds / duration_sec) * 100
print(f"Secondi senza alcun evento: {empty_seconds} su {duration_sec} ({pct_empty:.2f}%)")

# ---- NEW: spill statistics ----
print("=" * 50)
stats = spill_statistics(counts_per_sec)
print_report(stats)
print("=" * 50)

# ---- Plots ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

time_axis = np.arange(duration_sec)
ax1.step(time_axis, counts_per_sec, where="mid", color="royalblue", lw=0.8)
ax1.set_xlabel("Relative time from start [s]")
ax1.set_ylabel("Events / second")
ax1.set_title(f"Rate of the run {run}")
ax1.grid(alpha=0.3)

ax2.hist(counts_per_sec, bins=100, color="royalblue", edgecolor="royalblue", alpha=0.8)
ax2.set_yscale("log")
ax2.set_xlabel("Rate [events / s]")
ax2.set_ylabel("Number of 1-s intervals (log scale)")
ax2.set_title("Distribution of rate")
ax2.grid(alpha=0.3)

plt.tight_layout()
out_png = f"rate_analysis_{run}.png"
plt.savefig(out_png, dpi=200)
plt.savefig(f"rate_analysis_{run}.pdf")
print(f"Grafico salvato in {out_png}")

# ---- NEW: zoom on ~6 periods, so the spill structure is visible ----
period = stats["period_autocorr"] or stats["period_median_runs"] or 35
win = int(min(duration_sec, max(120, 6 * period)))
start = max(0, min(duration_sec // 2, duration_sec - win))   # window in the middle of the run
fig2, axz = plt.subplots(figsize=(7, 4.5))
axz.step(time_axis[start:start + win], counts_per_sec[start:start + win],
         where="mid", color="royalblue", lw=1.2)
axz.axhline(0.5 * stats["ref_level"], color="gray", ls="--", lw=0.8, label="spill threshold (50%)")
axz.set_xlabel("Relative time from start [s]")
axz.set_ylabel("Events / second")
axz.set_title(f"Run {run}: zoom on {win} s")
axz.legend(frameon=False)
axz.grid(alpha=0.3)
fig2.tight_layout()
fig2.savefig(f"rate_zoom_{run}.png", dpi=200)
fig2.savefig(f"rate_zoom_{run}.pdf")
print(f"Zoom salvato in rate_zoom_{run}.png")

plt.show()