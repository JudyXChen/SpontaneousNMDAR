> Current layout: recorded data is in `results/calcium_total/spont/conditions/1_cluster_5_NMDAR/`.
> For current run commands and shared settings, see [the shared model guide](../../../../doc/shared_model.md).
> The notes below describe the historical condition.

# Calibrated spontaneous release with 5 NMDARs at -67 mV

This leaf is a matched derivative of
`model/calcium_total/spont/1_cluster_10_NMDAR`. It changes only the NMDAR
count and receptor-patch walls needed for the 5-versus-10 comparison. It
retains the same geometry, 278 AMPARs, 500-molecule glutamate packet,
calcium-handling network, and calibrated explicit NMDAR magnesium block
model.

The five NMDARs occupy PSD walls 1688 and 1696, with a combined area of
1,328.5 nm². This is 4.3% above half of the 10-NMDAR patch area (1,273.4 nm²),
so the receptor densities differ by approximately 4.2%. The glutamate release
site is centered above this patch.

## Fixed-voltage protocol

The postsynaptic membrane is clamped at **-67 mV**. There is no action
potential waveform and no NEURON-generated variable-rate table. All
voltage-dependent rates are scalar parameters in `model.bngl`:

- eight VSCC state-transition rates;
- the open-VSCC calcium-event rate;
- calibrated NMDAR magnesium block and unblock rates; and
- the open, unblocked NMDAR calcium-event rate.

The NMDAR block calculation retains the calibrated model's inherited
**1 mM extracellular magnesium** condition. AMPAR binding, opening, closing,
and desensitization are ligand-dependent and are unchanged.

All receptors, channels, calcium, buffers, pumps, and ER machinery are
initialized at time zero. The 500-molecule glutamate packet is released at
10 ms, after a resting equilibration interval. A matched comparison run uses
450,000 iterations at the 1e-7 s timestep, giving 10 ms before release and
35 ms after release.

## Rate equations

VSCC gating uses

```text
a_i(V) = alpha_i * exp(V / V_i)
b_i(V) = beta_i  * exp(-V / V_i)
```

with the alpha and beta values in s^-1 and voltage constants in mV supplied
for the parent model. Open-channel calcium entry uses

```text
k_VSCC(V) = gamma_VSCC * N_A/(2*F) * V
            * (0.393 - exp(-V/80.36)) / (1 - exp(V/80.36))
```

where the leading `V` is converted to volts. The positive exponent in the
denominator is intentional.

The calibrated NMDAR block/unblock equations are

```text
k_B(V, Mg) = 1130.6388608604038 * Mg * exp(-V/20.091025673263097)
k_U(V)     = 11462.54604245 * exp(V/81.78933361)
```

with V in mV, Mg in mM, and rates in s^-1. The NMDAR calcium-event rate uses
the inherited 4.5 pS calcium-specific conductance and +3 mV reversal
potential.

## Running

The current shared entry point uses the 45 ms protocol in parameters.json.
Run from this directory with the CellBlender Python interpreter:

```bash
"$PYTHON" -B model.py --seeds 1
```

That command writes only inside this spontaneous leaf.
