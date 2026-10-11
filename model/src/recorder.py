"""
Generative model of DNA-based exposure recording in a differentiating human
cortical progenitor population.

Two coupled processes:
  (1) a branching process of progenitor division and cell-cycle exit,
  (2) sequential symbol insertion into a tandem prime-editing array, which
      occurs on division and, at a separately specified rate, independently
      of division.

Each modelling choice that could reasonably go the other way is a switch, so
the alternative stays runnable and the cost of the choice is measurable:

  error_model      "confusion" lets a miscalled symbol be read as another
                   symbol in the alphabet, so clock symbols can be misread as
                   exposure symbols. "sentinel" makes a miscall produce a
                   unique detectable code that can never generate a false
                   positive, which is optimistic about the readout.
  contraction_mode "delete" excises the monomer and shifts the remainder left,
                   which is what recombination between tandem repeats does.
                   "erase" blanks a monomer in place and so preserves the
                   positional order of everything else.
  write_timing     "pre" writes before the division is resolved, so both
                   daughters inherit the symbol. "post" gives the symbol to
                   one daughter only. Section 2.7 reports what the choice costs.
  writes_per_div   "poisson" allows more than one monomer to be filled in a
                   single cycle, which is what in vivo DNA Typewriter data
                   require at high editing rates. "bernoulli" allows at most one.
  fid_alphabet     size of the landmark symbol alphabet. 0 makes every pulse
                   distinguishable, which assumes as many orthogonal inducible
                   channels as pulses. 1 is a single shared landmark symbol.
                   k cycles k symbols over the pulses.
  fid_dyn_range    leak of the landmark channel between pulses. Infinite is a
                   perfectly silent landmark channel.

Event identity is tracked through the tree: every write gets a unique id that
descendants inherit, so the number of independent recording events behind a
set of marked nuclei can be counted rather than assumed.

Nothing is fitted to a single dataset. Every parameter is exposed.
"""
from dataclasses import dataclass
import numpy as np

UNWRITTEN, CLOCK, EXPOSURE_A, EXPOSURE_B, ERRORSYM = 0, 1, 2, 3, 4
FID_BASE = 10          # landmark pulse j writes FID_BASE + (j mod alphabet)
FIDUCIAL = FID_BASE

def is_fiducial(x):
    return x >= FID_BASE


@dataclass
class Params:
    # ---- population -------------------------------------------------------
    T: int = 150
    n_founder: int = 1500
    tc_start: float = 30.0
    tc_end: float = 110.0
    pp_start: float = 0.30
    pp_end: float = 0.06
    nn_start: float = 0.06
    nn_end: float = 0.13
    fate_midpoint: float = 45.0
    fate_steep: float = 0.06
    # spread, in days, of the per-lineage fate-switch midpoint. Zero gives
    # every lineage the same schedule, which produces far less clonal
    # heterogeneity than cerebral organoids show:
    # Lindenhofer et al. 2024 find lineage sizes spanning five orders of
    # magnitude, with 80 per cent of progeny from 5 per cent of lineages.
    fate_mid_sd: float = 0.0
    # a second progenitor compartment that stays proliferative longer, as a
    # stand-in for outer radial glia and other late-cycling progenitors. A
    # fraction of founders get their fate-switch midpoint shifted later by
    # late_mid_shift days, inherited by their descendants.
    late_fraction: float = 0.0
    late_mid_shift: float = 0.0
    # outer radial glia self-renew while producing neurons, so a persistent
    # compartment divides almost entirely asymmetrically and cycles slowly.
    # Modelling it as "same fate probabilities, later switch" instead makes the
    # population expand without bound and is not what the tissue does.
    late_pp: float = 0.02        # P(symmetric proliferative) in the late compartment
    late_nn: float = 0.02        # P(symmetric neurogenic) in the late compartment
    late_tc_factor: float = 1.5  # its cell-cycle length, relative to the main one
    # Day on which the second compartment starts behaving differently. Before
    # it, flagged cells are indistinguishable from the main population. This is
    # what makes a gliogenic switch expressible: glia divide after the neurons
    # are made, not from the beginning. Zero reproduces the original behaviour.
    late_start_day: float = 0.0

    # ---- array and writing -------------------------------------------------
    L: int = 50
    q_div: float = 0.18          # writing probability per division
    q_nondiv: float = 0.0        # writing probability per day, post-mitotic cell
    # Writing per unit time in EVERY cell, cycling or not. q_nondiv above
    # applies only after cycle exit; this is the regime the Fisher-information
    # analyses of a polymerase ticker tape assume, where the clock runs on wall
    # time rather than on division. Zero in every design but the benchmark.
    q_time: float = 0.0          # writing probability per day, any cell
    clock_share: float = 1.0
    fid_share: float = 8.0
    exp_share: float = 6.0
    dynamic_range: float = 4.0

    # ---- schedule ----------------------------------------------------------
    n_fiducial: int = 8
    fid_duration: int = 2
    fiducial_distinct: bool = True
    fid_alphabet: int = 0        # 0 = all distinct, 1 = shared, k = cyclic
    epochs: tuple = (10, 30, 55, 85)
    epoch_duration: int = 1
    extra_exposures: tuple = ()

    # ---- readout -----------------------------------------------------------
    recovery: float = 0.60
    symbol_error: float = 0.02
    contraction: float = 0.01
    # amplicon recovery of a long tandem array is incomplete: Choi et al. 2022
    # recovered 12.5 +- 4.3 of 20 monomers and 8.4 +- 3.3 of 12. span_mean = 1
    # makes a recovered array recovered whole, which makes extra length free.
    span_mean: float = 1.0
    span_sd: float = 0.0

    # ---- construct and population mechanisms ------------------------------
    # Each is off at its default, so the central design is reproduced exactly
    # when none is set. Section 2.2 sweeps them.
    #
    # Editor silencing. The writing rate is multiplied by 2**(-t/q_halflife),
    # so a progressively silenced cassette writes less as the experiment runs.
    # Infinite keeps the rate constant.
    q_halflife: float = np.inf
    # Position-dependent editing efficiency along the array. The write
    # probability at monomer s is scaled by 2**(-s/q_pos_halflife). Published
    # DNA Typewriter saturation profiles show later monomers editing less
    # efficiently than earlier ones. Infinite is uniform.
    q_pos_halflife: float = np.inf
    # Death. A per-day hazard on post-mitotic cells, plus an extra hazard that
    # falls only on cells which left the cycle before `death_split`. The second
    # is early-born-biased death, whose direction is not obvious in advance.
    death_pm: float = 0.0
    death_pm_early: float = 0.0
    death_split: int = 40
    # Transit-amplifying intermediate progenitors. A fraction of the daughters
    # that would have become post-mitotic instead divide `ipc_rounds` further
    # times on a shorter cycle and then exit. This is a different object from
    # the persistent asymmetric compartment of section 2.3: it amplifies a
    # progenitor division rather than self-renewing.
    ipc_fraction: float = 0.0
    ipc_rounds: int = 0
    ipc_tc_factor: float = 0.6

    # ---- modelling switches, see the module docstring ----------------------
    error_model: str = "confusion"
    contraction_mode: str = "delete"
    write_timing: str = "pre"
    writes_per_div: str = "bernoulli"
    fid_dyn_range: float = np.inf

    # ---- bookkeeping -------------------------------------------------------
    cap: int = 40000
    seed: int = 0
    track_events: bool = True

    def fiducial_days(self):
        if self.n_fiducial <= 0:
            return np.array([], dtype=int)
        return np.unique(np.linspace(8, self.T - 8, self.n_fiducial).round().astype(int))

    def fid_symbol(self, j):
        """Symbol written by pulse j, given the landmark alphabet."""
        k = self.fid_alphabet
        if k == 0:
            return FID_BASE + j if self.fiducial_distinct else FID_BASE
        return FID_BASE + (j % max(k, 1))

    def alphabet(self):
        """Every symbol code the design can produce, for the confusion model."""
        codes = [CLOCK, EXPOSURE_A]
        codes += [int(c) for (_, _, c) in self.extra_exposures]
        codes += [self.fid_symbol(j) for j in range(len(self.fiducial_days()))]
        return np.array(sorted(set(codes)), dtype=np.int8)


def _logistic(t, lo, hi, mid, steep):
    return lo + (hi - lo) / (1.0 + np.exp(-steep * (t - mid)))


def _sel(x, i):
    """Index a per-cell parameter, or pass a scalar through unchanged."""
    return x[i] if isinstance(x, np.ndarray) else x


def simulate(p: Params):
    rng = np.random.default_rng(p.seed)
    L = p.L

    fid_days, fid_pulse = {}, {}
    fdl = p.fiducial_days()
    for j, d in enumerate(fdl):
        for k in range(p.fid_duration):
            fid_days[int(d) + k] = p.fid_symbol(j)
            fid_pulse[int(d) + k] = j
    epoch_days = {}
    for i, e in enumerate(p.epochs):
        for k in range(p.epoch_duration):
            epoch_days[int(e) + k] = (i, EXPOSURE_A)
    for (d, dur, code) in p.extra_exposures:
        for k in range(int(dur)):
            epoch_days[int(d) + k] = (len(p.epochs) + list(
                dict.fromkeys(c for _, _, c in p.extra_exposures)).index(code), int(code))

    n = p.n_founder
    sym = np.zeros((n, L), dtype=np.int8)
    wday = np.full((n, L), -1, dtype=np.int16)
    eid = np.zeros((n, L), dtype=np.int32)       # identity of the writing event
    slot = np.zeros(n, dtype=np.int16)
    ndiv = np.zeros(n, dtype=np.int16)
    cyc = np.ones(n, dtype=bool)
    founder = np.arange(n, dtype=np.int32)
    late_flag = (rng.random(n) < p.late_fraction) if p.late_fraction > 0 else None
    if p.fate_mid_sd > 0 or p.late_fraction > 0:
        fmid = np.full(n, float(p.fate_midpoint))
        if p.fate_mid_sd > 0:
            fmid = rng.normal(p.fate_midpoint, p.fate_mid_sd, size=n)
        if p.late_fraction > 0:
            late = rng.random(n) < p.late_fraction
            fmid = fmid + late * p.late_mid_shift
        fmid = np.clip(fmid, 8.0, 400.0)
    else:
        fmid = None
    exposed = np.zeros((n, len(p.epochs)), dtype=bool)
    recorded = np.zeros((n, len(p.epochs)), dtype=bool)
    # day the cell left the cycle, -1 while still cycling. Only needed for
    # early-born-biased death, but it is one array and it makes the write-day
    # density decomposable by birth cohort.
    exitd = np.full(n, -1, dtype=np.int16)
    # remaining transit-amplifying divisions. 0 is an apical progenitor.
    tac = np.zeros(n, dtype=np.int8)

    cyc_frac, pool_size, writable = [], [], []
    div_day = np.zeros((n, p.T), dtype=bool)
    nondiv_writes = np.zeros(p.T, dtype=np.int64)
    scale, next_eid = 1.0, 1

    def shares_for(t):
        """Channel activities on day t, as (codes, normalised weights)."""
        on = t in epoch_days
        exp_active = on and epoch_days[t][1] == EXPOSURE_A
        s = [(CLOCK, p.clock_share),
             (EXPOSURE_A, p.exp_share if exp_active
              else p.exp_share / max(p.dynamic_range, 1e-9))]
        if on and epoch_days[t][1] != EXPOSURE_A:
            s.append((epoch_days[t][1], p.exp_share))
        if t in fid_days:
            s.append((fid_days[t], p.fid_share))
        elif np.isfinite(p.fid_dyn_range) and fdl.size:
            # a leaky landmark channel writes its symbol between pulses too.
            # the leaked symbol is the one whose cassette is closest in time,
            # which is the pulse the reader will attribute it to.
            j = int(np.argmin(np.abs(fdl - t)))
            s.append((p.fid_symbol(j), p.fid_share / max(p.fid_dyn_range, 1e-9)))
        codes = np.array([x[0] for x in s])
        w = np.array([x[1] for x in s], dtype=float)
        return codes, w / w.sum()

    def _subset(keep):
        """Restrict every per-cell array to `keep`, a mask or an index array.

        Used by both the population cap, where it is subsampling and the scale
        factor compensates, and by death, where it is real attrition and the
        scale factor must not.
        """
        nonlocal sym, wday, eid, slot, ndiv, cyc, founder, exposed, recorded
        nonlocal div_day, fmid, late_flag, exitd, tac
        sym, wday, eid = sym[keep], wday[keep], eid[keep]
        slot, ndiv, cyc = slot[keep], ndiv[keep], cyc[keep]
        founder, exposed, recorded = founder[keep], exposed[keep], recorded[keep]
        div_day, exitd, tac = div_day[keep], exitd[keep], tac[keep]
        if fmid is not None:
            fmid = fmid[keep]
        if late_flag is not None:
            late_flag = late_flag[keep]

    for t in range(p.T):
        mid = fmid if fmid is not None else p.fate_midpoint
        tc = _logistic(t, p.tc_start, p.tc_end, mid, p.fate_steep)
        pp = _logistic(t, p.pp_start, p.pp_end, mid, p.fate_steep)
        nn = _logistic(t, p.nn_start, p.nn_end, mid, p.fate_steep)
        if late_flag is not None and t >= p.late_start_day:
            tc = np.where(late_flag, np.asarray(tc) * p.late_tc_factor, tc)
            pp = np.where(late_flag, p.late_pp, pp)
            nn = np.where(late_flag, p.late_nn, nn)
        p_div = np.minimum(1.0, 24.0 / tc)
        pn = np.maximum(0.0, 1.0 - pp - nn)

        cyc_frac.append(cyc.mean())
        pool_size.append(scale * cyc.sum())
        codes, w = shares_for(t)

        idx = np.flatnonzero(cyc)
        if idx.size and p.ipc_rounds > 0:
            # transit-amplifying progenitors run a shorter cycle than the
            # apical progenitors that produced them
            pd = np.asarray(_sel(p_div, idx), dtype=float) * np.ones(idx.size)
            p_div_ipc = np.minimum(
                1.0, 24.0 / (np.asarray(tc, dtype=float) * max(p.ipc_tc_factor, 1e-9)))
            pd = np.where(tac[idx] > 0,
                          np.asarray(_sel(p_div_ipc, idx), dtype=float) * np.ones(idx.size),
                          pd)
            dividing = idx[rng.random(idx.size) < pd]
        else:
            dividing = idx[rng.random(idx.size) < _sel(p_div, idx)] if idx.size else idx
        writable.append(scale * dividing.size)

        # ------------- resolve the division first, so that a sibling row
        # ------------- exists before any symbol is written into it
        sib = None
        if dividing.size:
            u = rng.random(dividing.size)
            pp_d, pn_d = _sel(pp, dividing), _sel(pn, dividing)
            kind = np.where(u < pp_d, 2, np.where(u < pp_d + pn_d, 1, 0))
            tac_d = tac[dividing]
            par_tac = np.zeros(dividing.size, dtype=np.int8)
            sib_tac = np.zeros(dividing.size, dtype=np.int8)
            pr = np.zeros(dividing.size, dtype=bool)   # parent stays cycling as an IPC
            sr = np.zeros(dividing.size, dtype=bool)   # sibling starts cycling as an IPC
            if p.ipc_rounds > 0:
                is_ipc = tac_d > 0
                # an IPC division is symmetric: both daughters continue only
                # while amplifying rounds remain, and exit together when they
                # are spent
                kind = np.where(is_ipc, np.where(tac_d > 1, 2, 0), kind)
                par_tac = np.where(is_ipc, tac_d - 1, 0).astype(np.int8)
                sib_tac = par_tac.copy()
                # an apical division's post-mitotic daughter becomes a
                # transit-amplifying progenitor instead, with probability
                # ipc_fraction
                ap = ~is_ipc
                pr = ap & (kind == 0) & (rng.random(dividing.size) < p.ipc_fraction)
                sr = ap & (kind <= 1) & (rng.random(dividing.size) < p.ipc_fraction)
                par_tac = np.where(pr, p.ipc_rounds, par_tac).astype(np.int8)
                sib_tac = np.where(sr, p.ipc_rounds, sib_tac).astype(np.int8)
            ndiv[dividing] += 1
            div_day[dividing, t] = True
            n_before = sym.shape[0]
            sym = np.concatenate([sym, sym[dividing]])
            wday = np.concatenate([wday, wday[dividing]])
            eid = np.concatenate([eid, eid[dividing]])
            slot = np.concatenate([slot, slot[dividing]])
            ndiv = np.concatenate([ndiv, ndiv[dividing]])
            cyc = np.concatenate([cyc, (kind == 2) | sr])
            tac = np.concatenate([tac, sib_tac])
            exitd = np.concatenate([exitd, exitd[dividing]])
            founder = np.concatenate([founder, founder[dividing]])
            if fmid is not None:
                fmid = np.concatenate([fmid, fmid[dividing]])
            if late_flag is not None:
                late_flag = np.concatenate([late_flag, late_flag[dividing]])
            exposed = np.concatenate([exposed, exposed[dividing]])
            recorded = np.concatenate([recorded, recorded[dividing]])
            div_day = np.concatenate([div_day, div_day[dividing]])
            sib = np.arange(n_before, n_before + dividing.size)
            cyc[dividing] = (kind >= 1) | pr
            tac[dividing] = par_tac

        # ------------- writing on division ---------------------------------
        # editor silencing scales the rate by day, position-dependent
        # efficiency scales it by how far along the array the next monomer is
        q_t = (p.q_div * 2.0 ** (-t / p.q_halflife)
               if np.isfinite(p.q_halflife) else p.q_div)
        if dividing.size:
            q_eff = np.full(dividing.size, q_t, dtype=float)
            if np.isfinite(p.q_pos_halflife):
                q_eff = q_eff * 2.0 ** (-slot[dividing] / p.q_pos_halflife)
            if p.writes_per_div == "poisson":
                mult = rng.poisson(q_eff, size=dividing.size)
            else:
                mult = (rng.random(dividing.size) < q_eff).astype(int)
            # for a post-replicative edit the carrier daughter is chosen once
            # per division, so successive monomers stay in the same cell
            to_par = rng.random(dividing.size) < 0.5
            for step in range(int(mult.max()) if mult.size else 0):
                act = (mult > step)
                par = dividing[act]
                sb = sib[act]
                keep_par = to_par[act]
                free = slot[par] < L
                par, sb, keep_par = par[free], sb[free], keep_par[free]
                if par.size == 0:
                    continue
                pick = codes[rng.choice(len(codes), size=par.size, p=w)]
                ids = np.arange(next_eid, next_eid + par.size, dtype=np.int32)
                next_eid += par.size
                s = slot[par]
                if p.write_timing == "pre":
                    # the edit precedes the division, so both daughters carry it
                    sym[par, s] = pick;  sym[sb, s] = pick
                    wday[par, s] = t;    wday[sb, s] = t
                    eid[par, s] = ids;   eid[sb, s] = ids
                    slot[par] = s + 1;   slot[sb] = s + 1
                else:
                    # the flap is installed on one strand, so one daughter carries it
                    rows_ = np.where(keep_par, par, sb)
                    sym[rows_, s] = pick
                    wday[rows_, s] = t
                    eid[rows_, s] = ids
                    slot[rows_] = s + 1
                if t in epoch_days:
                    ei, code = epoch_days[t]
                    if ei < recorded.shape[1]:
                        hit = pick == code
                        if p.write_timing == "pre":
                            recorded[par[hit], ei] = True
                            recorded[sb[hit], ei] = True
                        else:
                            recorded[rows_[hit], ei] = True

        if t in epoch_days and epoch_days[t][0] < exposed.shape[1]:
            exposed[:, epoch_days[t][0]] = True

        # ------------- division-independent writing -------------------------
        # q_nondiv writes only in cells that have left the cycle; q_time writes
        # in every cell. They are separate switches and may both be on.
        for _rate, _target in ((p.q_nondiv, "postmitotic"), (p.q_time, "all")):
          if _rate > 0:
            pmi = (np.flatnonzero(~cyc) if _target == "postmitotic"
                   else np.arange(cyc.size))
            if pmi.size:
                qn = (_rate * 2.0 ** (-t / p.q_halflife)
                      if np.isfinite(p.q_halflife) else _rate)
                qn = np.full(pmi.size, qn, dtype=float)
                if np.isfinite(p.q_pos_halflife):
                    qn = qn * 2.0 ** (-slot[pmi] / p.q_pos_halflife)
                wr = pmi[rng.random(pmi.size) < qn]
                wr = wr[slot[wr] < L]
                if wr.size:
                    pick = codes[rng.choice(len(codes), size=wr.size, p=w)]
                    ids = np.arange(next_eid, next_eid + wr.size, dtype=np.int32)
                    next_eid += wr.size
                    s = slot[wr]
                    sym[wr, s] = pick
                    wday[wr, s] = t
                    eid[wr, s] = ids
                    slot[wr] = s + 1
                    nondiv_writes[t] = wr.size
                    if t in epoch_days:
                        ei, code = epoch_days[t]
                        if ei < recorded.shape[1]:
                            recorded[wr[pick == code], ei] = True

        # ------------- record the day each cell left the cycle ---------------
        # a cell leaves once and never returns, so the first day it is seen
        # non-cycling is its exit day
        just_left = (~cyc) & (exitd < 0)
        if just_left.any():
            exitd[just_left] = t

        # ------------- death -------------------------------------------------
        # applied after writing, so a cell that wrote today still contributes
        # that write to its surviving relatives but not to the endpoint itself
        if p.death_pm > 0 or p.death_pm_early > 0:
            pm = ~cyc
            haz = np.where(pm, p.death_pm, 0.0)
            if p.death_pm_early > 0:
                haz = haz + np.where(pm & (exitd >= 0) & (exitd < p.death_split),
                                     p.death_pm_early, 0.0)
            die = rng.random(cyc.size) < haz
            if die.any():
                _subset(~die)

        # ------------- cap ---------------------------------------------------
        if sym.shape[0] > p.cap:
            keep = rng.choice(sym.shape[0], p.cap, replace=False)
            scale *= sym.shape[0] / p.cap
            _subset(keep)

    return dict(late=late_flag, exit_day=exitd, tac=tac,
                sym=sym, wday=wday, eid=eid, slot=slot, ndiv=ndiv, cyc=cyc,
                div_day=div_day, founder=founder, exposed=exposed, recorded=recorded,
                scale=scale, cyc_frac=np.array(cyc_frac), pool=np.array(pool_size),
                writable=np.array(writable), nondiv_writes=nondiv_writes,
                params=p, fid_days=np.array(sorted(fid_days)),
                fid_pulse_of_day=fid_pulse, n_events=next_eid - 1)


def readout(sim, p: Params, rng=None):
    """Apply recovery, symbol miscalling and array contraction.

    Returns (sym, got, wday, eid) with all four positionally aligned, because
    contraction in "delete" mode moves symbols and any analysis that looks up
    the day a monomer was written has to follow the move.
    """
    rng = rng or np.random.default_rng(p.seed + 7717)
    sym = sim["sym"].copy()
    wday = sim["wday"].copy()
    eid = sim["eid"].copy()
    n, L = sym.shape
    got = rng.random(n) < p.recovery
    written = sym != UNWRITTEN

    if p.symbol_error > 0:
        flip = written & (rng.random(sym.shape) < p.symbol_error)
        if p.error_model == "sentinel":
            sym[flip] = ERRORSYM
        else:
            # a miscall is read as some other symbol in the alphabet. this is
            # the channel through which an abundant clock symbol becomes a
            # false exposure symbol.
            alpha = p.alphabet()
            if alpha.size > 1:
                r = np.flatnonzero(flip.ravel())
                cur = sym.ravel()[r]
                draw = alpha[rng.integers(0, alpha.size - 1, size=r.size)]
                # map a self-draw onto the last element, so the read symbol
                # is always different from the true one
                same = draw == cur
                draw[same] = alpha[-1]
                still = draw == cur
                draw[still] = alpha[0]
                flat = sym.ravel().copy()
                flat[r] = draw
                sym = flat.reshape(sym.shape)

    if p.span_mean < 1.0:
        # only a prefix of the array is recovered from the amplicon, and the
        # longer the array the more of it is lost
        if p.span_sd > 0:
            m, s = p.span_mean, p.span_sd
            v = max(min(s * s, m * (1 - m) * 0.99), 1e-9)
            k = m * (1 - m) / v - 1
            span = rng.beta(max(m * k, 1e-3), max((1 - m) * k, 1e-3), size=n)
        else:
            span = np.full(n, p.span_mean)
        keepto = np.ceil(span * L).astype(int)
        beyond = np.arange(L)[None, :] >= keepto[:, None]
        sym[beyond] = UNWRITTEN
        wday[beyond] = -1
        eid[beyond] = 0
        written = sym != UNWRITTEN

    if p.contraction > 0:
        lost = written & (rng.random(sym.shape) < p.contraction)
        if p.contraction_mode == "erase":
            sym[lost] = UNWRITTEN
        else:
            # recombination excises the monomer; everything distal shifts up
            order = np.argsort(lost.astype(np.int8), axis=1, kind="stable")
            sym = np.take_along_axis(sym, order, axis=1)
            wday = np.take_along_axis(wday, order, axis=1)
            eid = np.take_along_axis(eid, order, axis=1)
            nlost = lost.sum(1)
            tail = np.arange(L)[None, :] >= (L - nlost)[:, None]
            sym[tail] = UNWRITTEN
            wday[tail] = -1
            eid[tail] = 0
    return sym, got, wday, eid
