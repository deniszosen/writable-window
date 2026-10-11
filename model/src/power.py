"""Power for a two-arm comparison of culture-level mark rates.

The obvious way to compute the smallest detectable difference is

    delta = (z_0.975 + z_0.8) * sqrt(2) * sqrt(sigma_b^2/k + pi(1-pi)/(k n))

evaluated from one culture per arm upwards. Three things are wrong with that,
and all three flatter the design, so none is used here:

1. The between-culture variance is not known, it is estimated from the same
   cultures. The test statistic is therefore a two-sample t on culture-level
   rates with 2(k-1) degrees of freedom, not a z. At two cultures per arm the
   difference is a factor of 1.8 in the detectable effect.

2. At one culture per arm there are zero degrees of freedom. The between-culture
   component is not estimable and the arm effect is perfectly confounded with
   culture identity, so there is no test and no detectable effect. Reporting
   a number there is the same error the paper diagnoses in nuclei, committed
   one level up, in cultures.

3. sigma_b is itself estimated from a handful of cultures, so it carries
   chi-squared uncertainty that is strongly asymmetric: from eight cultures the
   90 per cent interval runs from 0.71 to 1.80 times the estimate. Quoting a
   single curve implies the variance is pinned down when the upper tail is
   nearly twice the point estimate.

This module computes the detectable effect exactly, by solving the noncentral t
power equation, and returns the chi-squared band alongside the point estimate.
It needs scipy; nothing else in the pipeline does.
"""
import numpy as np
from scipy import stats, optimize


def se_diff(sd_between, pi, n_arrays, k):
    """Standard error of the difference in mean mark rate between two arms of
    k cultures each, with n_arrays recovered arrays sequenced per culture.

    Var(one culture's observed rate) = sigma_b^2 + pi(1-pi)/n. The second term
    is binomial in n because sampling arrays from a culture whose true rate is
    p_c is genuinely binomial in p_c; the clonal correlation that makes nuclei
    non-independent is what makes p_c vary between cultures, and it is therefore
    already inside sigma_b^2 rather than missing from this expression.
    """
    return float(np.sqrt(2.0 * (sd_between ** 2 + pi * (1 - pi) / n_arrays) / k))


def detectable_delta(sd_between, pi, n_arrays, k, power=0.80, alpha=0.05,
                     exact=True):
    """Smallest detectable difference in mark rate. None when k < 2.

    exact=True solves the noncentral t power equation. exact=False uses the
    (t_{1-alpha/2} + t_{power}) approximation, which is within a few per cent
    except at two cultures per arm, where it is optimistic.
    """
    df = 2 * (k - 1)
    if df < 1:
        return None
    se = se_diff(sd_between, pi, n_arrays, k)
    t_crit = stats.t.ppf(1 - alpha / 2, df)
    if not exact:
        return float((t_crit + stats.t.ppf(power, df)) * se)

    def shortfall(d):
        ncp = d / se
        upper = stats.nct.sf(t_crit, df, ncp)
        lower = stats.nct.cdf(-t_crit, df, ncp)
        # `lower` is the probability of rejecting in the wrong direction.
        # scipy (1.18) returns NaN for it at moderate and large positive ncp
        # instead of underflowing to zero, which aborts the root finder. A
        # noncentral t is (Z + ncp)/sqrt(V/df) with a positive denominator, so
        # T < 0 requires Z < -ncp and the term is bounded above by Phi(-ncp);
        # every ncp at which scipy returns NaN here puts that bound below
        # 1e-17, so zero is exact in double precision.
        if not np.isfinite(lower):
            lower = 0.0
        return upper + lower - power

    return float(optimize.brentq(shortfall, 1e-9, 50 * se))


def sigma_band(sd_between, m_cultures, conf=0.90):
    """Confidence interval on sigma_between from m cultures.

    (m-1) s^2 / sigma^2 is chi-squared with m-1 degrees of freedom, so the
    interval on sigma runs from s * sqrt((m-1)/chi2_upper) to
    s * sqrt((m-1)/chi2_lower), and is markedly asymmetric for the m of five to
    eight simulated cultures used here.
    """
    df = m_cultures - 1
    a = (1 - conf) / 2
    lo = sd_between * np.sqrt(df / stats.chi2.ppf(1 - a, df))
    hi = sd_between * np.sqrt(df / stats.chi2.ppf(a, df))
    return float(lo), float(hi)


def curve(sd_between, pi, net, n_arrays, k_max=8, m_cultures=8, power=0.80,
          alpha=0.05):
    """Detectable difference and fold change against cultures per arm.

    The fold change is the ratio of the two arms' net mark rates the design can
    resolve, net / (net - delta). It is undefined once delta reaches net, which
    is what happens at small k, and is reported as None rather than dropped.
    """
    lo_sd, hi_sd = sigma_band(sd_between, m_cultures)
    pts = []
    for k in range(1, k_max + 1):
        row = dict(cultures=k)
        for tag, sd in (("", sd_between), ("_sd_low", lo_sd), ("_sd_high", hi_sd)):
            d = detectable_delta(sd, pi, n_arrays, k, power, alpha)
            row["min_detectable_delta" + tag] = d
            row["min_detectable_fold" + tag] = (
                float(net / (net - d)) if d is not None and d < net else None)
        dz = float((stats.norm.ppf(1 - alpha / 2) + stats.norm.ppf(power))
                   * se_diff(sd_between, pi, n_arrays, k))
        row["min_detectable_delta_z"] = dz
        row["min_detectable_fold_z"] = (
            float(net / (net - dz)) if dz < net else None)
        pts.append(row)
    return dict(sd_between=float(sd_between), sd_between_ci=[lo_sd, hi_sd],
                m_cultures=int(m_cultures), arrays_per_culture=int(n_arrays),
                pi=float(pi), net=float(net), points=pts)


def cultures_needed(sd_between, pi, net, n_arrays, fold, m_cultures=8,
                    power=0.80, alpha=0.05, k_max=40):
    """Cultures per arm that resolve the given fold change, at the point
    estimate of sigma and at the pessimistic end of its interval."""
    target = net - net / fold
    out = {}
    for tag, sd in (("point", sd_between),
                    ("sigma_high", sigma_band(sd_between, m_cultures)[1])):
        out[tag] = None
        for k in range(2, k_max + 1):
            d = detectable_delta(sd, pi, n_arrays, k, power, alpha)
            if d is not None and d <= target:
                out[tag] = k
                break
    return out


def demo():
    """Self-check: every correction points the same way, and the z-based
    calculation is the optimistic end of each comparison."""
    sd, pi, net = 0.01401, 0.2204, 0.1783
    c = curve(sd, pi, net, 1000)
    p1, p2, p3 = c["points"][0], c["points"][1], c["points"][2]
    assert p1["min_detectable_fold"] is None, "one culture per arm has no test"
    assert p1["min_detectable_fold_z"] is not None, "the z form reports one"
    assert p2["min_detectable_fold"] > p2["min_detectable_fold_z"], "t beats z"
    assert p3["min_detectable_fold"] > p3["min_detectable_fold_z"]
    lo, hi = c["sd_between_ci"]
    assert lo < sd < hi and (hi - sd) > (sd - lo), "the band is asymmetric upwards"
    # the approximation is optimistic against the exact solution, and the gap
    # closes as the degrees of freedom grow
    g2 = (detectable_delta(sd, pi, 1000, 2)
          - detectable_delta(sd, pi, 1000, 2, exact=False))
    g8 = (detectable_delta(sd, pi, 1000, 8)
          - detectable_delta(sd, pi, 1000, 8, exact=False))
    assert g2 > g8 >= 0, (g2, g8)
    k = cultures_needed(sd, pi, net, 1000, 1.33)
    assert k["point"] is not None and k["sigma_high"] > k["point"], k
    print(f"ok  one culture per arm: no test (the z form reports "
          f"{p1['min_detectable_fold_z']:.2f}-fold)")
    print(f"ok  three cultures per arm: {p3['min_detectable_fold']:.2f}-fold "
          f"exact against {p3['min_detectable_fold_z']:.2f}-fold under z")
    print(f"ok  sigma_between {sd:.4f}, 90 per cent interval {lo:.4f} to {hi:.4f}")
    print(f"ok  cultures per arm for 1.33-fold: {k['point']} at the point "
          f"estimate, {k['sigma_high']} at the pessimistic end")


if __name__ == "__main__":
    demo()
