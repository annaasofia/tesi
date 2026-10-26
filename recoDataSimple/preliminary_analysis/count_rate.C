#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>

using namespace std;

// Funzione ausiliaria per convertire "HHMMSS" in secondi della giornata
long timeToSeconds(const string& t_str) {
    if (t_str.length() < 6) return -1;
    int hh = stoi(t_str.substr(0, 2));
    int mm = stoi(t_str.substr(2, 2));
    int ss = stoi(t_str.substr(4, 2));
    return hh * 3600 + mm * 60 + ss;
}

void count_rate(TTree* tree = nullptr) {
    if (!tree) {
        tree = (TTree*)gDirectory->Get("simpleEvent");
        if (!tree) {
            cerr << "Errore: TTree 'simpleEvent' non trovato!" << endl;
            return;
        }
    }

    // Disattiva tutti i branch e attiva solo 'Time' per massimizzare la velocita' di I/O
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("Time", 1);

    Char_t time_buf[80];
    tree->SetBranchAddress("Time", time_buf);

    Long64_t nentries = tree->GetEntries();
    if (nentries == 0) return;

    // 1. Estrai il tempo iniziale e finale
    tree->GetEntry(0);
    string start_str(time_buf);
    long t_start = timeToSeconds(start_str);

    tree->GetEntry(nentries - 1);
    string stop_str(time_buf);
    long t_stop = timeToSeconds(stop_str);

    // Gestione passaggio della mezzanotte se t_stop risulta numericamente minore
    if (t_stop < t_start) {
        t_stop += 24 * 3600;
    }

    long duration_sec = t_stop - t_start + 1;
    long dur_min = duration_sec / 60;
    long dur_rem_sec = duration_sec % 60;

    cout << "========================================" << endl;
    cout << "Inizio run: " << start_str << " | Fine run: " << stop_str << endl;
    cout << "Durata totale: " << duration_sec << " s (" 
         << dur_min << " min " << dur_rem_sec << " s)" << endl;
    cout << "========================================" << endl;

    // Vettore dei conteggi per secondo [0, duration_sec - 1]
    vector<Long64_t> counts_per_sec(duration_sec, 0);

    cout << "Scansione degli eventi in corso..." << endl;
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        long t_cur = timeToSeconds(time_buf);
        if (t_cur < t_start) t_cur += 24 * 3600;

        long rel_sec = t_cur - t_start;
        if (rel_sec >= 0 && rel_sec < duration_sec) {
            counts_per_sec[rel_sec]++;
        }
    }

    // 2. Controllo dei secondi a zero eventi
    long empty_seconds = 0;
    vector<long> empty_indices;
    for (long s = 0; s < duration_sec; ++s) {
        if (counts_per_sec[s] == 0) {
            empty_seconds++;
            if (empty_indices.size() < 10) {
                empty_indices.push_back(s);
            }
        }
    }

    cout << "Secondi senza alcun evento: " << empty_seconds << " su " << duration_sec;
    if (duration_sec > 0) {
        cout << " (" << (double)empty_seconds / duration_sec * 100.0 << "%)" << endl;
    } else {
        cout << endl;
    }

    if (empty_seconds > 0) {
        cout << "Esempi di secondi senza eventi (secondi relativi dall'inizio): ";
        for (auto idx : empty_indices) cout << idx << "s ";
        if (empty_seconds > 10) cout << "...";
        cout << endl;
    }

    // 3. Istogrammi
    Long64_t max_rate = *max_element(counts_per_sec.begin(), counts_per_sec.end());
    TH1D* h_rate = new TH1D("h_rate", "Distribuzione Eventi/s;Eventi per secondo;Numero di secondi", 
                            100, 0, max_rate * 1.05);

    TH1D* h_time_profile = new TH1D("h_time_profile", "Rate temporale del run;Tempo relativo [s];Eventi/s", 
                                    duration_sec, 0, duration_sec);

    for (long s = 0; s < duration_sec; ++s) {
        h_rate->Fill(counts_per_sec[s]);
        h_time_profile->SetBinContent(s + 1, counts_per_sec[s]);
    }

    TCanvas* c1 = new TCanvas("c_rate", "Analisi Rate Temporale", 1200, 500);
    c1->Divide(2, 1);

    c1->cd(1);
    gPad->SetGrid();
    h_time_profile->SetLineColor(kAzure + 1);
    h_time_profile->Draw("HIST");

    c1->cd(2);
    gPad->SetGrid();
    gPad->SetLogy();
    h_rate->SetLineColor(kTeal + 2);
    h_rate->SetFillColorAlpha(kTeal + 2, 0.35);
    h_rate->Draw("HIST");

    // Ripristina la lettura su tutti i branch del TTree
    tree->SetBranchStatus("*", 1);
}