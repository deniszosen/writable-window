"""Metrics derived from the generative model.

The observable everywhere is the one an experiment actually has: does this
recovered array carry at least one exposure symbol. Alongside it:

  - the realised write-day density, which generalises the writing-opportunity
    density to designs where writing is not exclusively division-coupled;
  - the number of independent writing events behind a set of marked nuclei,
    because an edit made in a progenitor is inherited by its whole subclone
    and those nuclei are not independent observations;
  - between-culture variance of the mark rate and the sample size a contrast
    actually needs, which replaces counting marked cells;
  - landmark readers that are honest about what a counted landmark can and
    cannot bound, for landmark alphabets smaller than the number of pulses.
"""
import numpy as np
from recorder import (Params, simulate, readout, UNWRITTEN, CLOCK, FIDUCIAL,
                      FID_BASE, EXPOSURE_A, EXPOSURE_B, ERRORSYM)

Z_ALPHA_2 = 1.959963985
Z_BETA_80 = 0.8416212336


def _h(x):
    x = np.clip(np.asarray(x, float), 1e-12, 1 - 1e-12)
    return -(x * np.log2(x) + (1 - x) * np.log2(1 - x))


def mutual_information(pi, phi, prior=0.5):
    """Bits per recovered cell about a binary exposure."""
    m = prior * np.asarray(pi, float) + (1 - prior) * np.asarray(phi, float)
    return _h(m) - (prior * _h(pi) + (1 - prior) * _h(phi))


# --------------------------------------------------------------- mark rates --
def mark_rate(p: Params, seed=None):
    q = Params(**{**p.__dict__, "seed": (p.seed if seed is None else seed)})
    s = simulate(q)
    sym, got, _, _ = readout(s, q)
    if got.sum() == 0:
        return np.nan, s
    return float((sym[got] == EXPOSURE_A).any(1).mean()), s


def vehicle_mark_rate(p: Params, seed=None, paired=True):
    """The same schedule with the agent never applied.

    paired=True gives the vehicle culture the same seed as the exposed one, so
    the two share a population realisation and the difference of the two mark
    rates has far less Monte Carlo noise than it would with the vehicle as an
    independent simulation.
    """
    sd = (p.seed if seed is None else seed)
    q = Params(**{**p.__dict__, "epochs": (), "extra_exposures": (),
                  "seed": sd if paired else sd + 5000})
    s = simulate(q)
    sym, got, _, _ = readout(s, q)
    return float((sym[got] == EXPOSURE_A).any(1).mean()) if got.sum() else np.nan


def contrast(p: Params, seeds=(0, 1, 2), paired=True):
    """pi, phi and bits per recovered cell, with a standard error over seeds.

    Mutual information is computed per seed and then averaged, rather than from
    the averaged rates, because the map from rates to bits is not linear.
    """
    pis, phis, bits = [], [], []
    for sd in seeds:
        pi, _ = mark_rate(p, seed=p.seed + sd)
        phi = vehicle_mark_rate(p, seed=p.seed + sd, paired=paired)
        pis.append(pi)
        phis.append(phi)
        bits.append(float(mutual_information(pi, phi)))
    n = max(len(seeds), 1)
    sem = lambda v: float(np.std(v, ddof=1) / np.sqrt(n)) if n > 1 else 0.0
    return dict(pi=float(np.mean(pis)), phi=float(np.mean(phis)),
                pi_sem=sem(pis), phi_sem=sem(phis),
                net=float(np.mean(np.array(pis) - np.array(phis))),
                net_sem=sem(list(np.array(pis) - np.array(phis))),
                bits=float(np.mean(bits)), bits_sem=sem(bits),
                bits_from_mean_rates=float(mutual_information(np.mean(pis), np.mean(phis))),
                n_seeds=n)


# ------------------------------------------------------- writing over time --
def writable_profile(sim):
    """Per day: fraction of terminal cells whose lineage divided that day."""
    return sim["div_day"].mean(0)


def writing_opportunity_mass(sim):
    w = writable_profile(sim)
    return w / w.sum()


def write_day_density(sim):
    """Normalised density of symbols actually written, over developmental time.

    This is the quantity the closed form needs. It equals the normalised
    writing-opportunity density when writing is purely division-coupled, and it
    stays correct when it is not.
    """
    wd = sim["wday"]
    v = wd[wd >= 0]
    T = sim["params"].T
    h = np.bincount(v.astype(int), minlength=T)[:T].astype(float)
    return h / max(h.sum(), 1.0)


def per_epoch_recording(sim):
    return sim["recorded"].mean(0)


def epoch_gradient(sim):
    """Ratio of the first to the last epoch's recording probability."""
    r = per_epoch_recording(sim)
    return float(r[0] / r[-1]) if r[-1] > 0 else np.inf


def f_on(sim, p: Params, realised=True):
    """Share of a cell's writing opportunities that fall inside an exposure
    window. `realised=True` uses the days symbols were actually written, which
    generalises to division-independent writing."""
    m = write_day_density(sim) if realised else writing_opportunity_mass(sim)
    days = []
    for e in p.epochs:
        days.extend(range(int(e), int(e) + int(p.epoch_duration)))
    if not days:
        return 0.0
    d = np.asarray(days, int)
    d = d[d < len(m)]
    return float(m[d].sum())


# ------------------------------------------------- channel and error floors --
def confusion_rate(p: Params):
    """Probability that a written symbol is read as an exposure symbol through
    miscalling alone. Zero under the sentinel error model."""
    if p.symbol_error <= 0 or p.error_model == "sentinel":
        return 0.0
    a = p.alphabet().size
    return float(p.symbol_error / max(a - 1, 1))


def required_dynamic_range(sim, p: Params, alpha=0.1, realised=True):
    """Expected false writes per cell from channel leak are N_w (1 - f) s_off
    and true writes are N_w f s_on, with s_on = e/(c+e) and s_off =
    (e/D)/(c+e/D). Requiring the ratio to be at most alpha gives, for large D,

        D  >=  (1 - f) (c + e) / (c f alpha).
    """
    f = f_on(sim, p, realised=realised)
    if f <= 0:
        return f, np.inf
    return f, (1 - f) * (p.clock_share + p.exp_share) / (p.clock_share * f * alpha)


def dynamic_range_ceiling(sim, p: Params, realised=True):
    """The dynamic range beyond which channel leak is no longer what limits the
    false-mark rate, because symbol miscalling produces more false exposure
    symbols than the channel does. Setting the per-write leak e/(cD) equal to
    the per-write confusion rate gives D = e/(c eps)."""
    eps = confusion_rate(p)
    if eps <= 0:
        return np.inf
    return float(p.exp_share / (p.clock_share * eps))


def expected_marks(sim, p: Params, realised=True):
    """Closed-form mark rates including the miscalling floor. Returns
    (pi_hat, phi_hat)."""
    f = f_on(sim, p, realised=realised)
    n_w = float(np.minimum(sim["ndiv"] * p.q_div, p.L).mean()) if p.q_nondiv == 0 \
        else float(np.minimum(sim["slot"], p.L).mean())
    s_on = p.exp_share / (p.clock_share + p.exp_share)
    b = p.exp_share / max(p.dynamic_range, 1e-9)
    s_off = b / (p.clock_share + b)
    eps = confusion_rate(p)
    keep = 1.0 - p.symbol_error
    true_w = n_w * f * s_on * keep
    false_w = n_w * (1 - f) * s_off * keep
    conf_w = n_w * eps                      # any other symbol misread as exposure
    return (float(1 - np.exp(-(true_w + false_w + conf_w))),
            float(1 - np.exp(-(n_w * s_off * keep + conf_w))))


# --------------------------------------------------- independence of records --
def event_stats(sim, p: Params, sym=None, got=None, eid=None):
    """How many independent writing events lie behind the marked nuclei.

    An edit made in a progenitor is inherited by every descendant, so marked
    nuclei from one subclone are one observation, not many.
    """
    if sym is None:
        sym, got, _, eid = readout(sim, p)
    is_exp = sym == EXPOSURE_A
    rows = np.flatnonzero(got & is_exp.any(1))
    if rows.size == 0:
        return dict(n_marked=0, n_events=0, pseudoreplication=np.nan,
                    largest_clone_share=np.nan, median_clone=np.nan)
    ids = eid[rows][is_exp[rows]]
    ids = ids[ids > 0]
    uniq, counts = np.unique(ids, return_counts=True)
    return dict(n_marked=int(rows.size), n_events=int(uniq.size),
                pseudoreplication=float(rows.size / max(uniq.size, 1)),
                median_clone=float(np.median(counts)),
                mean_clone=float(counts.mean()),
                max_clone=int(counts.max()),
                largest_clone_share=float(counts.max() / counts.sum()),
                top_decile_share=float(np.sort(counts)[::-1][:max(1, uniq.size // 10)].sum()
                                       / counts.sum()))


def between_culture_variance(p: Params, seeds=range(8), paired=True):
    """Variance of a culture's population mark rate across independent
    cultures. This is the variance component that sequencing more nuclei from
    one culture cannot reduce."""
    nets, pis = [], []
    for sd in seeds:
        pi, _ = mark_rate(p, seed=p.seed + 100 * (sd + 1))
        phi = vehicle_mark_rate(p, seed=p.seed + 100 * (sd + 1), paired=paired)
        pis.append(pi)
        nets.append(pi - phi)
    return dict(pi_mean=float(np.mean(pis)), pi_var=float(np.var(pis, ddof=1)),
                pi_sd=float(np.std(pis, ddof=1)),
                net_mean=float(np.mean(nets)), net_sd=float(np.std(nets, ddof=1)),
                n_cultures=len(list(seeds)))


def design_effect(pi, var_between, n_nuclei):
    """Var(p_hat) = var_between + pi(1-pi)/n. Expressed as a multiplier on the
    binomial variance, the number of nuclei is worth less and less as it grows."""
    binom = pi * (1 - pi) / max(n_nuclei, 1)
    return float((var_between + binom) / binom) if binom > 0 else np.inf


def n_per_group(p1, p2, alpha=0.05, power=0.8, deff=1.0):
    """Nuclei per condition for a two-proportion contrast at the given power,
    inflated by the design effect."""
    if p1 == p2:
        return np.inf
    pbar = 0.5 * (p1 + p2)
    z_a = Z_ALPHA_2
    z_b = Z_BETA_80 if abs(power - 0.8) < 1e-9 else _z(power)
    num = (z_a * np.sqrt(2 * pbar * (1 - pbar))
           + z_b * np.sqrt(p1 * (1 - p1) + p2 * (1 - p2))) ** 2
    return float(num / (p1 - p2) ** 2 * deff)


def _z(pp):
    from math import sqrt, log
    # Moro/Beasley-Springer style inverse normal, adequate for power values
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    q = pp - 0.5
    if abs(q) <= 0.425:
        r = 0.180625 - q * q
        return q * (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) / \
            (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)
    r = pp if q < 0 else 1 - pp
    r = sqrt(-log(r))
    v = 2.32121276858 + r * (4.85014127135 + r * (-2.29796479134 + r * 0.00772454011)) \
        if r <= 5 else r
    return v if q > 0 else -v


def n_cultures_needed(net, net_sd, effect_fraction=0.5, alpha=0.05, power=0.8):
    """Cultures per arm to detect a change of `effect_fraction` in the net mark
    rate, given between-culture standard deviation."""
    if net_sd <= 0:
        return 1.0
    d = effect_fraction * net / net_sd
    return float(2 * ((Z_ALPHA_2 + Z_BETA_80) / d) ** 2) if d > 0 else np.inf


# ---------------------------------------------------------------- landmarks --
def _pulse_symbols(p: Params):
    return np.array([p.fid_symbol(j) for j in range(len(p.fiducial_days()))])


def brackets(sim, p: Params, sym=None, got=None, wday=None, reader="conservative",
             max_rows=20000, rng=None):
    """Bracket each recovered exposure symbol with the landmarks around it.

    reader='naive'        counts landmarks and assumes the count identifies
                          the pulse, which over-states precision whenever a
                          pulse was missed.
    reader='conservative' uses only what a counted landmark can support: the
                          earliest pulse consistent with the observed sequence
                          of landmark symbols on the left, and the latest one
                          consistent on the right. Bounds are valid by
                          construction unless the landmark channel leaks.
    """
    rng = rng or np.random.default_rng(p.seed + 31)
    if sym is None:
        sym, got, wday, _ = readout(sim, p)
    fdays = np.asarray(sorted(set(int(d) for d in p.fiducial_days())))
    psym = _pulse_symbols(p)
    is_exp = sym == EXPOSURE_A
    rows = np.flatnonzero(got & is_exp.any(1))
    if rows.size == 0 or fdays.size == 0:
        return dict(n=0, reader=reader, frac_two_sided=np.nan, acc_all=np.nan,
                    acc_two_sided=np.nan, median_width_all=np.nan,
                    median_width_two_sided=np.nan)
    if rows.size > max_rows:
        rows = rng.choice(rows, max_rows, replace=False)
    # pulse indices available for each landmark symbol
    by_sym = {int(s): np.flatnonzero(psym == s) for s in np.unique(psym)}
    n_p = psym.size

    lo_a, hi_a, ok_a, two_a, lowonly_a = [], [], [], [], []
    S = sym[rows]
    W = wday[rows]
    K = is_exp[rows].argmax(1)
    for r in range(rows.size):
        row, k = S[r], int(K[r])
        pos = np.flatnonzero(row >= FID_BASE)
        left = [int(row[i]) for i in pos if i < k]
        right = [int(row[i]) for i in pos if i > k]
        true_day = int(W[r, k])
        if reader == "naive":
            m = len(left)
            lo = fdays[min(m - 1, n_p - 1)] if m >= 1 else 0
            two = len(right) > 0 and m < n_p
            hi = fdays[min(m, n_p - 1)] if two else p.T
        else:
            # earliest chain on the left
            j = -1
            okl = True
            for s in left:
                cand = by_sym[s][by_sym[s] > j]
                if cand.size == 0:
                    okl = False
                    break
                j = int(cand[0])
            lo = fdays[j] if (okl and j >= 0) else 0
            # latest chain on the right
            jr = n_p
            okr = True
            for s in reversed(right):
                cand = by_sym[s][by_sym[s] < jr]
                if cand.size == 0:
                    okr = False
                    break
                jr = int(cand[-1])
            two = bool(right) and okr and jr < n_p
            # the latest day that pulse could have written is the last day of
            # its window, not the first
            hi = min(int(fdays[jr]) + int(p.fid_duration) - 1, p.T) if two else p.T
        lo_a.append(lo)
        hi_a.append(hi)
        two_a.append(bool(two))
        lowonly_a.append(bool(not two and lo > 0))
        ok_a.append(bool(lo <= true_day <= hi))
    lo_a = np.array(lo_a); hi_a = np.array(hi_a)
    ok_a = np.array(ok_a); two_a = np.array(two_a)
    lowonly = np.array(lowonly_a)
    width = hi_a - lo_a
    w2 = width[two_a]
    out = dict(n=int(rows.size), reader=reader,
               frac_two_sided=float(two_a.mean()),
               frac_lower_bound_only=float(lowonly.mean()),
               frac_no_landmark=float((~two_a & ~lowonly).mean()),
               acc_all=float(ok_a.mean()),
               acc_two_sided=float(ok_a[two_a].mean()) if two_a.any() else np.nan,
               median_width_all=float(np.median(width)),
               median_width_two_sided=float(np.median(w2)) if two_a.any() else np.nan,
               right_censored=float((hi_a >= p.T).mean()))
    # the distribution of two-sided bracket widths, so a reader can see the
    # spread rather than only the median
    out["width_quantiles_two_sided"] = (
        {str(q): float(np.percentile(w2, q)) for q in (10, 25, 50, 75, 90)}
        if w2.size >= 20 else None)
    return out


def epoch_assignment(sim, p: Params, sym=None, got=None, wday=None):
    """The naive counting rule, kept so the two readers can be compared."""
    return brackets(sim, p, sym=sym, got=got, wday=wday, reader="naive")


def pulse_capture(sim, p: Params, sym=None, got=None):
    """Probability that a lineage records a given landmark pulse."""
    if sym is None:
        sym, got, _, _ = readout(sim, p)
    psym = _pulse_symbols(p)
    if p.fid_alphabet == 0:
        return float(np.mean([(sym == s).any(1).mean() for s in psym]))
    # with a shared or cyclic alphabet a cell cannot say which pulse it caught;
    # report the per-pulse rate from the ground truth day of each landmark
    return float(np.mean([(sym == s).any(1).mean() for s in np.unique(psym)]))


# ------------------------------------------------------------------- order --
def order_accuracy(p: Params, day0, delta, seeds=(0,)):
    accs, ns = [], []
    for sd in seeds:
        q = Params(**{**p.__dict__, "epochs": (),
                      "extra_exposures": ((day0, 1, EXPOSURE_A), (day0 + delta, 1, EXPOSURE_B)),
                      "seed": p.seed + 900 + sd})
        s = simulate(q)
        sym, got, _, _ = readout(s, q)
        a, b = sym == EXPOSURE_A, sym == EXPOSURE_B
        both = got & a.any(1) & b.any(1)
        ns.append(int(both.sum()))
        if both.sum() >= 20:
            accs.append(float((a[both].argmax(1) < b[both].argmax(1)).mean()))
    return (float(np.mean(accs)) if accs else np.nan), int(np.sum(ns))


def clock_variability(sim, p: Params):
    sym = sim["sym"]
    rng = np.random.default_rng(3)
    if sym.shape[0] > 6000:
        sym = sym[rng.choice(sym.shape[0], 6000, replace=False)]
    fdays = np.asarray(sorted(set(int(d) for d in p.fiducial_days())))
    rows = []
    for j in range(len(fdays) - 1):
        vals = []
        for row in sym:
            f = np.flatnonzero(row >= FID_BASE)
            if f.size > j + 1:
                vals.append(int((row[f[j] + 1:f[j + 1]] == CLOCK).sum()))
        if len(vals) >= 50:
            v = np.asarray(vals, float)
            rows.append(dict(t0=int(fdays[j]), t1=int(fdays[j + 1]), n=len(vals),
                             mean=float(v.mean()),
                             cv=float(v.std() / v.mean()) if v.mean() > 0 else np.nan))
    return rows


def usable_cells(pi, recovery, nuclei_per_condition=15500, competent=0.50):
    """A yield metric that does not survive scrutiny, kept only so the
    comparison can be run. It multiplies nuclei by the net mark rate, which is
    neither a sample size nor a count of marked cells; the power calculation
    above is what should be used instead."""
    return nuclei_per_condition * recovery * competent * pi


# ================================================= readers, power, matching ==
#   - readers richer than "does this array carry an exposure symbol", which
#     turns the information lower bound into a measured gap;
#   - a power calculation that does not pretend the between-culture variance is
#     known, since pretending otherwise is itself a pseudoreplication error;
#   - total expansion as a constraint, so a population variant can be matched
#     to the observed expansion instead of being compared at a different one.

def _entropy_bits(counts, correct=True):
    """Plug-in entropy in bits, with the Miller-Madow bias correction.

    The correction matters because the richer readers below spread the same
    number of arrays over many more bins, and the plug-in estimator is biased
    downwards in entropy and therefore upwards in mutual information.
    """
    c = np.asarray(counts, float)
    n = c.sum()
    if n <= 0:
        return 0.0
    pr = c[c > 0] / n
    h = float(-(pr * np.log2(pr)).sum())
    if correct:
        h += (pr.size - 1) / (2.0 * n * np.log(2.0))
    return h


def _discrete_mi(x_exp, x_veh, correct=True):
    """Bits per recovered cell about a binary exposure, for a discrete
    observable, with a uniform prior over the two conditions."""
    keys = np.union1d(np.unique(x_exp), np.unique(x_veh))
    ce = np.array([(x_exp == k).sum() for k in keys], float)
    cv = np.array([(x_veh == k).sum() for k in keys], float)
    pe, pv = ce / max(ce.sum(), 1), cv / max(cv.sum(), 1)
    mix = 0.5 * pe + 0.5 * pv
    h_mix = _entropy_bits(mix * (ce.sum() + cv.sum()), correct)
    return float(h_mix - 0.5 * _entropy_bits(ce, correct)
                 - 0.5 * _entropy_bits(cv, correct))


def _observables(sym, got, p: Params):
    """The three readers, for the recovered arrays of one condition.

    binary     does the array carry an exposure symbol. The one used throughout.
    count      how many exposure symbols it carries.
    joint      how many exposure symbols, and how many monomers are written at
               all. The second is what makes the first interpretable: one
               exposure symbol among two written monomers is strong evidence,
               and one among twelve is weak, because the miscall floor scales
               with the number of written monomers.
    """
    a = sym[got]
    n_exp = (a == EXPOSURE_A).sum(1)
    n_written = (a != UNWRITTEN).sum(1)
    cap = int(min(p.L, 40))
    return dict(binary=(n_exp > 0).astype(np.int64),
                count=np.minimum(n_exp, 8),
                joint=np.minimum(n_exp, 8) * 64 + np.minimum(n_written, cap))


def reader_information(p: Params, seeds=(0, 1, 2), correct=True, n_perm=3):
    """Bits per recovered cell for each reader, with a permutation floor.

    The permutation floor is the mutual information obtained after shuffling
    the condition labels across the pooled arrays. It is zero in expectation
    and its size is the estimation bias that remains after correction, so a
    reader's gain is only believable if it is large compared with the floor.
    """
    out = {k: [] for k in ("binary", "count", "joint")}
    floor = {k: [] for k in out}
    for sd in seeds:
        q = Params(**{**p.__dict__, "seed": p.seed + sd})
        se = simulate(q)
        sym_e, got_e, _, _ = readout(se, q)
        v = Params(**{**q.__dict__, "epochs": (), "extra_exposures": ()})
        sv = simulate(v)
        sym_v, got_v, _, _ = readout(sv, v)
        oe, ov = _observables(sym_e, got_e, q), _observables(sym_v, got_v, v)
        rng = np.random.default_rng(1000 + sd)
        for k in out:
            out[k].append(_discrete_mi(oe[k], ov[k], correct))
            pool = np.concatenate([oe[k], ov[k]])
            fl = []
            for _ in range(n_perm):
                idx = rng.permutation(pool.size)
                fl.append(_discrete_mi(pool[idx[:oe[k].size]],
                                       pool[idx[oe[k].size:]], correct))
            floor[k].append(float(np.mean(fl)))
    n = max(len(seeds), 1)
    sem = lambda v_: float(np.std(v_, ddof=1) / np.sqrt(n)) if n > 1 else 0.0
    return {k: dict(bits=float(np.mean(out[k])), bits_sem=sem(out[k]),
                    permutation_floor=float(np.mean(floor[k])))
            for k in out}


def expansion(sim):
    """Total fold expansion of the founder population, cycling and
    post-mitotic together. This is the quantity the population model is
    calibrated against, so a variant that changes it is not comparable to the
    central design until it is matched."""
    return float(sim["scale"] * sim["sym"].shape[0] / sim["params"].n_founder)


def match_expansion(p: Params, target, key="pp_start", lo=0.02, hi=0.60,
                    tol=0.02, iters=12, seed=None):
    """Bisect one fate probability so the variant expands like the reference.

    Added because a population variant that adds divisions also adds cells, and
    comparing its recording gradient against the central design at a different
    total expansion compares two cultures that would not both look like an
    organoid. Section 2.2 excludes a variant on exactly this ground; matching
    is the constructive form of the same argument.
    """
    sd = p.seed if seed is None else seed
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        e = expansion(simulate(Params(**{**p.__dict__, key: mid, "seed": sd})))
        if abs(np.log(e / target)) < tol:
            break
        if e > target:
            hi = mid
        else:
            lo = mid
    return float(mid), float(e)
