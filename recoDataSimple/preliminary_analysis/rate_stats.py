"""Spill statistics from a per-second event-count array (numpy only)."""
import numpy as np


def _run_lengths(mask):
    """Lengths of consecutive True runs and False runs, dropping the first and
    last run of the array (they are truncated by the start/stop of the run)."""
    mask = np.asarray(mask, dtype=bool)
    change = np.flatnonzero(np.diff(mask.astype(np.int8))) + 1
    edges = np.concatenate(([0], change, [len(mask)]))
    lengths = np.diff(edges)
    values = mask[edges[:-1]]
    if len(lengths) > 2:
        lengths, values = lengths[1:-1], values[1:-1]
    else:
        lengths, values = lengths[:0], values[:0]
    return lengths[values], lengths[~values]


def autocorr_period(counts, lag_min=5, lag_max=300, min_peak=0.2):
    """First significant peak of the autocorrelation (period in seconds).
    Returns (None, ac) if no peak has autocorrelation >= min_peak."""
    x = np.asarray(counts, dtype=float)
    x = x - x.mean()
    n = len(x)
    if n < 4 * lag_min or x.std() == 0:
        return None, None
    f = np.fft.rfft(x, 2 * n)
    ac = np.fft.irfft(f * np.conj(f))[:n]
    ac /= ac[0]
    hi = min(lag_max, n // 2)
    lags = np.arange(lag_min, hi)
    seg = ac[lag_min:hi]
    is_peak = (seg[1:-1] > seg[:-2]) & (seg[1:-1] >= seg[2:])
    peak_idx = np.flatnonzero(is_peak) + 1
    if len(peak_idx) == 0:
        return None, ac
    best = seg[peak_idx].max()
    if best < min_peak:  # no significant periodicity (e.g. noise or isolated gaps only)
        return None, ac
    # first, not highest: avoids multiples of the period
    good = peak_idx[(seg[peak_idx] >= 0.5 * best) & (seg[peak_idx] >= min_peak)]
    if len(good) == 0:
        return None, ac
    return int(lags[good[0]]), ac


def spill_statistics(counts, thr_frac=0.5, full_frac=0.8):
    """
    counts: events per second.
    - reference level = 95th percentile of counts (robust against the max)
    - 'in spill' second:  counts > thr_frac  * ref   (used for duty cycle / durations)
    - 'plateau' second:   counts > full_frac * ref   (full seconds only, used for the rate:
                          the first/last second of a spill is only partially filled)
    """
    c = np.asarray(counts, dtype=float)
    ref = np.percentile(c, 95)
    in_spill = c > thr_frac * ref
    plateau = c > full_frac * ref

    out = dict(
        n_sec=len(c),
        total_events=float(c.sum()),
        mean_rate_run=float(c.mean()),
        ref_level=float(ref),
        max_rate=float(c.max()),
        duty_cycle=float(in_spill.mean()),
        n_plateau=int(plateau.sum()),
    )
    if plateau.any():
        out["rate_spill_mean"] = float(c[plateau].mean())
        out["rate_spill_std"] = float(c[plateau].std(ddof=1)) if plateau.sum() > 1 else 0.0
        out["rate_spill_median"] = float(np.median(c[plateau]))
    on, off = _run_lengths(in_spill)
    out["spill_len_median"] = float(np.median(on)) if len(on) else None
    out["gap_len_median"] = float(np.median(off)) if len(off) else None
    out["n_spills"] = int(len(on))
    out["period_median_runs"] = (
        out["spill_len_median"] + out["gap_len_median"]
        if len(on) and len(off) else None)
    out["period_autocorr"], out["_ac"] = autocorr_period(c)
    return out


def print_report(s):
    print(f"Reference level (95th pct)         : {s['ref_level']:.0f} events/s")
    print(f"Max rate in a 1-s bin              : {s['max_rate']:.0f} events/s")
    print(f"Mean rate over the whole run       : {s['mean_rate_run']:.1f} events/s "
          f"({s['total_events']:.3e} events / {s['n_sec']} s)")
    print(f"Duty cycle (seconds with beam)     : {s['duty_cycle']*100:.1f} %")
    if "rate_spill_mean" in s:
        print(f"Rate in spill (plateau, {s['n_plateau']} s)   : "
              f"mean {s['rate_spill_mean']:.0f}, median {s['rate_spill_median']:.0f}, "
              f"std {s['rate_spill_std']:.0f} events/s")
    if s["spill_len_median"] is not None:
        print(f"Spill length / gap (median, runs)  : {s['spill_len_median']:.0f} s / "
              f"{s['gap_len_median']:.0f} s  ({s['n_spills']} complete spills)")
    if s["period_autocorr"] is not None:
        print(f"Period from autocorrelation        : {s['period_autocorr']} s")
    else:
        print("Period from autocorrelation        : none (no significant periodicity)")
    if s["duty_cycle"] > 0.98:
        print("NOTE: duty cycle ~100%: no spill structure visible in this run, "
              "spill/gap/period numbers are not meaningful.")
    print("Note: 1-s time resolution -> durations are quantized to +/-1 s.")