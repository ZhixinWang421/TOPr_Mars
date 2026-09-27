# TOPr_Mars

`TOPr_Mars` is a self-contained MATLAB implementation of a moving-window
total orbital power ratio (TOPr) with Monte Carlo sensitivity analysis for
Mars orbital-parameter time series.

For every moving window, the reported ratio is

\[
\mathrm{TOPr}=
\frac{\text{power in the union of the target-frequency bands}}
     {\text{power between }F_{t,\min}\text{ and }F_{t,\max}}.
\]

The Monte Carlo ensemble varies selected analysis choices (by default, the
window length and one shared orbital-band width factor) and reports pointwise
percentiles. It does not add observational noise, recenter or clip the
ensemble, or convert TOPr to `1 - TOPr`.

## Files

- `TOPr_Mars.m` — complete calculation, Monte Carlo analysis, and plotting.
- `README.md` — usage and method summary.

`TOPr_Mars.m` contains the required power-decomposition logic and does not
depend on a separate `pdan.m` file.

## Requirements

- MATLAB
- Signal Processing Toolbox (`pmtm`)

The input must be a finite, evenly sampled, two-column matrix:

1. time in kyr, strictly increasing;
2. observed value.

All period, window, and sampling-interval settings are consequently expressed
in kyr, and frequencies are in kyr\(^{-1}\).

## Quick start

```matlab
data = readmatrix('mars_series.txt');
[topr,mc] = TOPr_Mars(data,1000);

% Recommended representative curve: [time_kyr, median_topr]
writematrix(topr,'TOPr_Mars_median.csv');

% Full pointwise percentile summary
writetable(mc.topr_percentiles,'TOPr_Mars_percentiles.csv');
```

For non-interactive or batch use:

```matlab
cfg = struct('MakePlot',false,'ShowProgress',false);
[topr,mc] = TOPr_Mars(data,1000,cfg);
```

## Default configuration

| Setting | Default |
|---|---:|
| Target periods (kyr) | `[125 94.6 54.4 1260 116.4 65.8 118.4 52.6 37.2]` |
| Fixed-reference window | `1260 kyr` |
| Monte Carlo window range | `945–1575 kyr` |
| Multitaper time-bandwidth product, `NW` | `2` |
| Fixed-reference bandwidth factor | `1.2` |
| Monte Carlo bandwidth-factor range | `0.9–1.2` |
| Total-power frequency range | `0–0.1 kyr^-1` |
| Window step | `10 samples` |
| FFT length | `4096` |
| Percentiles | `[5 25 50 75 95]` |
| Random seed | `20260912` |

The bandwidth factor is shared by all nine target periods within each Monte
Carlo realization, keeping the band definitions internally coherent.
Sampling interval and `NW` remain fixed unless `SampleRangeKyr` or multiple
`NWValues` are explicitly supplied.

Example with additional sampling-interval and `NW` sensitivity:

```matlab
cfg = struct( ...
    'SampleRangeKyr',[4.5 5.5], ...
    'NWValues',[2 2.5 3], ...
    'Seed',20260912, ...
    'PlotBaseline',true);

[topr,mc] = TOPr_Mars(data,1000,cfg);
```

Any default can be overridden with a same-named field in `cfg`. The function
rejects unknown field names so that spelling mistakes do not silently change
an analysis.

## Main outputs

- `topr` — two-column representative series: time and pointwise Monte Carlo
  median TOPr, restricted to the requested common-support interval.
- `mc.topr_percentiles` — pointwise percentile table.
- `mc.topr_simulations` — every Monte Carlo realization (returned by default).
- `mc.deterministic_baseline` — fixed-parameter reference calculation.
- `mc.parameter_draws` — the sampled parameter values.
- `mc.realized_lower_cutoff_kyr_inv` and
  `mc.realized_upper_cutoff_kyr_inv` — realized target-band limits.
- `mc.config` — the resolved configuration, including defaults.

The random-number generator is local to the function, so a fixed `Seed`
reproduces the same parameter draws without changing MATLAB's global random
stream.

## Interpretation

The percentile envelopes quantify sensitivity to the method parameters that
were sampled. They are not p-values, confidence intervals for measurement
error, or age-model uncertainty intervals.

If a later analysis needs spectral uncertainty, analyze every column of
`mc.topr_simulations` and summarize the resulting spectra. The spectrum of the
pointwise median TOPr curve is generally not the same as the median spectrum
of all Monte Carlo realizations.

## Method background

The moving-window power-decomposition concept follows the approach documented
for `pdan` in:

Li, M. et al. (2016), *Geology*. https://doi.org/10.1130/G37970.1

The Monte Carlo layer borrows the general idea of propagating plausible
analysis-parameter choices through repeated calculations; its intervals have
the specific conditional interpretation described above.

