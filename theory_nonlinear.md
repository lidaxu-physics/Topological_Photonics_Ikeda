# Nonlinear Kerr-Comb Solver — Theory and Implementation

This document explains the nonlinear part of the app: how `TMM_app_nonlinear.py`
generates Kerr frequency combs in a lattice of ring resonators. The solver is a
**dispersive Ikeda map** driven by the **same z-resolved transfer-matrix propagator**
that the linear app already uses. There is no effective tight-binding Hamiltonian and
no hidden link rings — the comb is propagated through the actual geometry, one round
trip at a time, with the Kerr nonlinearity added in the ring waveguides. Beyond the
detuning sweep, a **slow-time analysis** (Section 16) evolves the settled field at a fixed
detuning to expose the temporal dynamics (stationary soliton, breather, or chaos), and the
whole run can be saved and reloaded (Section 17).

The sweep lives in `lattice_kerr_sweep`; the slow-time evolution in `slow_time_evolve`.

### Notation

$E$ intracavity field envelope; $\xi$ normalized position around a ring ($0$ to $1$ per
round trip); $\tau$ fast time within one round-trip period; $\mu$ comb-mode index;
$\delta$ round-trip detuning; $\alpha$ intrinsic loss; $\kappa$ bus coupling; $\kappa_J$
link coupling; $\gamma$ Kerr coefficient; $D_2$ second-order dispersion; $s=-1$ anomalous
dispersion sign; $N_s$ grid points per site ring; $dz_s = 1/N_s$ and $dz_\ell$ the site
and link step lengths.

---

## 1. Scope and philosophy

The linear app computes the cold-cavity drop spectrum of an IQH/AQH ring lattice with a
transfer-matrix method (TMM): it builds a sparse propagator $R(\omega)$ from the geometry
(site rings, link rings, directional couplers, synthetic-flux phases) and either solves
the steady state $(I-R)E = s$ or iterates the time map below.

The nonlinear solver **reuses that exact propagator** and adds three things: a
fast-time / comb-mode axis on every grid point so the field can carry a comb; the Kerr
nonlinearity in the site-ring waveguides; and an explicit second-order dispersion $D_2$
(the bare geometry has a flat FSR). The guiding principle: at zero nonlinearity
($\gamma = 0$) the solver must reduce, exactly, to the linear TMM (verified in Section 15).

---

## 2. The round-trip (Ikeda) map

Light in a ring circulates and returns to the coupler every round-trip time
$\tau_R = L_R/v_g$ (ring length over group velocity; both normalized to $1$, so
$\tau_R = 1$ and the free spectral range is $\Omega_R = 2\pi/\tau_R = 2\pi$). The Ikeda
map treats the field as a sequence of round trips indexed by $m$:

$$E^{m+1} = \mathcal{B}[\mathcal{P}(E^{m})]$$

where $\mathcal{P}$ propagates the field once around the ring under the nonlinear wave
equation (Section 4), and $\mathcal{B}$ groups the **boundary** operations — the directional
couplers (the scattering that builds the lattice), the detuning phase, and the CW pump. What
makes this an *Ikeda map* (rather than a mean-field / LLE model) is that the couplers act as
**discrete scattering events** — the field is scattered where it physically meets each
coupler — instead of being smeared into a continuous, round-trip-averaged evolution. This
discreteness is what faithfully captures coupled-resonator lattices, where the per-round-trip
coupling and the supermode splittings are a finite fraction of the FSR. This is the framework
used by Hashemi & Mittal for 2D ring-lattice combs.

In this implementation the round trip is **z-resolved**, and this is where "boundary applied
once per round trip" becomes literal only for a *single* ring: each ring is discretized into
$N_s$ points, and one application of the propagator $R$ advances every grid point by one step
$dz_s = 1/N_s$. The couplers sit at fixed grid points, so **the boundary is applied at each
$dz$ step, at the coupler's location** — the field is scattered when it passes that point,
not lumped at the end of the lap. A full round trip is $N_s$ applications of $R$; each coupler
fires once per round trip, but *where* it sits, interleaved with the propagation. The detuning
($e^{i\delta}$, distributed as $e^{i\delta dz_s}$ per step) and the pump are likewise applied
every step. By default the Kerr is distributed too — interleaved with the $R$
applications by Strang splitting (Section 8) — so the linear propagation, the boundary, and
the nonlinearity are all spread across the sub-steps; nothing is a single per-round-trip
kick (that is the optional `kerr_lumped` mode). The abstract $\mathcal{B}[\mathcal{P}]$
form above is the conceptual picture; the z-resolved scheme realizes it as one continuous
split-step integration.

---

## 3. The physical device

The lattice is a 2D array of **site rings** (the resonant elements) coupled through
**link rings** (auxiliary rings held anti-resonant) that mediate the coupling and imprint
the synthetic gauge phase. A CW pump enters through an **input bus** coupled to one site
ring; the comb is collected at a **drop bus** coupled to a different site ring.

```
   input bus ──►(input ring)══(link)══(ring)══ ... ══(drop ring)──► drop bus (comb out)
                     ║            ║                        ║
                  (link)       (link)        ...        (link)
                     ║            ║                        ║
                  (ring)══(link)══(ring)══ ...
```

Couplings are realized by **directional couplers**. A coupler between two waveguides with
field coupling $\kappa$ is the lossless unitary scattering matrix

$$\begin{pmatrix} a_{\mathrm{out}} \\ b_{\mathrm{out}} \end{pmatrix} = \begin{pmatrix} t & i\kappa \\ i\kappa & t \end{pmatrix} \begin{pmatrix} a_{\mathrm{in}} \\ b_{\mathrm{in}} \end{pmatrix}$$

with $t = \sqrt{1-\kappa^2}$, so a fraction $\kappa^2$ of the power crosses over and
$t^2 = 1-\kappa^2$ stays. Three coupler types appear: **bus couplers** $\kappa$ (input bus
to input ring, drop ring to drop bus) set the external loss of the two bus rings; **link
couplers** $\kappa_J$ (site ring to link ring) set the inter-ring hopping; and an
**intrinsic loss** $\alpha$ damps every segment by $e^{-\alpha dz/2}$ per step.

The link rings are detuned to anti-resonance (parameter $\beta_0\eta$) and carry a
**position-dependent Peierls phase**: the round-trip phase of a link depends on its
location, which is a synthetic magnetic flux through each plaquette (IQH) or the
brick-wall Floquet drive (AQH). That phase, baked into $R$, gives the lattice its
topological band. The comb is read only at the drop ring:

$$s_{\mathrm{drop}} = i\kappa E[\mathrm{drop}]$$

---

## 4. The intra-cavity nonlinear Schrodinger equation

Inside a site ring the slowly varying envelope $E(\xi,\tau)$ obeys

$$\frac{\partial E}{\partial \xi} = -\frac{\alpha}{2} E - i s D_2 \frac{\partial^2 E}{\partial \tau^2} + i \gamma |E|^2 E$$

The three terms are intrinsic loss, group-velocity dispersion, and the Kerr effect
(self-phase modulation). In the comb-mode basis $\partial^2/\partial\tau^2 \to -\mu^2$, so
the dispersion becomes a per-mode phase $+i s D_2 \mu^2$ (Section 9); $s=-1$ is anomalous
dispersion, required for **bright** solitons. The detuning and the inter-ring coupling are
not in this equation — they live in the round-trip boundary. Over one round trip, mode
$\mu$ acquires an amplitude factor $e^{-\alpha/2}$ and a dispersion phase
$e^{i s D_2 \mu^2}$; the boundary then applies the detuning $e^{i\delta}$, the couplers,
and the pump.

---

## 5. The linear transfer-matrix machine that is reused

The propagator is built by `build_propagator` from a template, and the discrete-time map
is

$$E^{n+1} = R(\omega) E^{n} + s$$

$R(\omega)$ is sparse over the waveguide grid points; $s$ injects the pump at the input
bus. Three facts make it directly reusable.

**(i) One application of $R$ is one $dz_s$**, so $N_s$ applications make one round trip.

**(ii) $R(\omega)$ factorizes in frequency.** Every entry has the form

$$R(\omega)_{ab} = c_{ab} e^{i \omega dz_{\mathrm{kind}}}$$

where the prefactor $c_{ab}$ holds the directional-coupler amplitude, the structural
(Peierls) phase, and the loss, and the only $\omega$-dependence is $e^{i \omega dz_s}$ on
site edges or $e^{i \omega dz_\ell}$ on link edges.

**(iii) The drop port** is

$$s_{\mathrm{drop}} = i\kappa p_s E[\mathrm{drop}]$$

with $p_s = e^{i \omega dz_s - \alpha dz_s/2}$ — the identical formula used by `solve_one`
and `time_evolve`.

---

## 6. Field representation: a comb on every grid point

The nonlinear field is a 2D array $E[g,\mu]$ over waveguide grid point $g$ and comb mode
$\mu$. The entry $E[g,\mu]$ is the **mode-domain amplitude** of comb line $\mu$ at grid
point $g$; the fast time $\tau$ and the mode $\mu$ are Fourier conjugates. The pump is CW,
so it lives entirely in mode $0$:

$$E[\mathrm{src},0] \to E[\mathrm{src},0] + p_{\mathrm{in}}\cdot \mathrm{src\_vals}$$

where $p_{\mathrm{in}}$ is the pump field amplitude. The first $n_{\mathrm{sites}} N_s$
grid points are the site rings; the rest are link rings. Both are nonlinear waveguides, so
the Kerr effect acts on all of them (Section 8).

---

## 7. The per-mode factorization of R

Because $R(\omega)$ factorizes, the whole comb propagates with two sparse multiplies per
step plus a per-mode phase. Comb mode $\mu$ sees frequency

$$\omega_\mu = \delta + s D_2 \mu^2$$

so the per-mode site and link propagation phases are

$$m_s(\mu) = e^{i \omega_\mu dz_s}$$

$$m_\ell(\mu) = e^{i \omega_\mu dz_\ell}$$

and the $\omega$-independent coupling matrices are

$$C^{\mathrm{site}}_{ab} = -c_{ab} e^{-\alpha dz_s/2}$$

$$C^{\mathrm{link}}_{ab} = -c_{ab} e^{i \beta_0\eta/N_\ell - \alpha dz_\ell/2 + i \phi_{ab}}$$

with $\phi_{ab}$ the per-link Peierls phase. The linear update is

$$E \leftarrow m_s \odot (C^{\mathrm{site}} E) + m_\ell \odot (C^{\mathrm{link}} E)$$

This is exact because, for a single mode $\mu$,

$$R(\omega_\mu) E = e^{i \omega_\mu dz_s}(C^{\mathrm{site}} E) + e^{i \omega_\mu dz_\ell}(C^{\mathrm{link}} E)$$

every site edge sharing $dz_s$ (so its phase factors out of the site sum), and likewise
for link edges. At $\mu=0$ this reproduces $R(\delta)$ entry by entry — which is precisely
why the $\gamma = 0$ (zero-nonlinearity) limit is the linear TMM (Section 15). For $\mu \ne 0$ the only
difference is the dispersion phase $s D_2 \mu^2$.

---

## 8. The Kerr nonlinearity and four-wave mixing

Under pure Kerr the intensity at each point in fast time is conserved, so the exact
solution over a step is a pointwise phase rotation:

$$E(\tau) \to E(\tau) e^{i \gamma |E(\tau)|^2 dz}$$

In the mode basis the same operation $|E|^2 E$ is a convolution — four-wave mixing:

$$|E|^2 E \longleftrightarrow \sum_{\mu_1+\mu_2-\mu_3=\mu} \hat{E}_{\mu_1} \hat{E}_{\mu_2} \hat{E}^{*}_{\mu_3}$$

Two photons at modes $\mu_1, \mu_2$ annihilate and create photons at $\mu_3$ and
$\mu = \mu_1 + \mu_2 - \mu_3$, transferring power from the pump line ($\mu=0$) into new
comb teeth and cascading outward. The process has parametric gain; once the intracavity
power exceeds a threshold set by the cavity loss, the CW background goes modulationally
unstable and a comb grows from noise. Re-applying the same self-phase-modulation phase
every round trip then builds and shapes the comb — Turing patterns, chaos, breathers, and,
in the right regime, cavity solitons.

**Split-step accuracy.** Dispersion is diagonal in frequency and Kerr is diagonal in fast
time, so the map handles them in alternation — each sub-step exact in its own basis. By
default (`kerr_lumped = False`) the Kerr is **distributed across the round trip by
second-order Strang splitting**: a half-Kerr, then $N_s$ pairs of [linear step,
full-Kerr] with the last Kerr a half-step. Consecutive half-Kerrs at round-trip boundaries
merge into full steps, so the whole sweep is one continuous split-step integration with a
splitting error $O(dz^2)$ — exactly Ekström's scheme (his Eq 5.7–5.9). Setting
`kerr_lumped = True` instead applies a single Kerr kick once per round trip ($\mathrm{kdz}=1$):
$N_s$× fewer FFTs and faster, valid when the per-round-trip nonlinear phase is small, but a
coarser (first-order) approximation.

**A normalization subtlety.** $E[g,\mu]$ holds mode amplitudes, but the Kerr phase needs
the physical time-domain field. The convention is that a DC amplitude $E_0$ corresponds to
a constant physical field $E_0$, i.e. $\mathrm{field} = \mathrm{FFT}(E)$ and
$E = \mathrm{IFFT}(\mathrm{field})$. Using IFFT first (which carries a $1/N$ factor) would
shrink the intensity the Kerr sees by $N^2$ and erase the nonlinearity. The Kerr effect is
a material property of the waveguide, so it is applied at **every** ring grid point — site
rings and link rings alike (site and link rings share the same normalized length, so a
single round-trip Kerr phase $\gamma |E|^2$ applies consistently to both).

---

## 9. Dispersion

The bare geometry has a constant FSR — the propagation phase is linear in $\omega$, so
equally spaced longitudinal modes stay equally spaced and there is no second-order
dispersion. The solver adds anomalous group-velocity dispersion explicitly, as a per-mode
round-trip phase

$$\phi_D(\mu) = s D_2 \mu^2$$

folded into $\omega_\mu$ (and hence $m_s, m_\ell$). The quadratic form is the leading
dispersion; $s=-1$ (taken whenever $D_2 \ge 0$) is the sign where dispersion balances the
Kerr self-focusing to form bright solitons. High modes acquire rapidly growing phase
($\propto \mu^2$), which sets the comb bandwidth and pushes far modes off resonance.

---

## 10. The supermode / band picture

A single ring has one resonance per FSR. Coupling $N$ site rings splits it into $N$
collective **supermodes** whose resonance frequencies spread into the topological band.
Each comb mode $\mu$ sees this same lattice, so the full mode structure is
two-dimensional: a longitudinal axis ($\mu$, the comb) and a transverse axis (the
supermode index, the band). The cold-cavity drop spectrum is exactly this band — sweeping
the probe frequency, the drop power peaks at each supermode resonance. The bulk supermodes
cluster into the bulk band(s); the topological **edge modes** appear as a ladder of
resonances that traverse the bulk gap — spectrally *inside* the gap between bulk bands, not
at the band edges. It is these gap-crossing edge resonances, spatially localized on the
lattice boundary, that make the topological pumping robust. Because
the nonlinear solver uses the real propagator, this band emerges automatically: at
$\gamma = 0$ the pump-line sweep reproduces the drop spectrum exactly (Section 15), at low
power it does so approximately, and a comb forms when the pump is brought near a supermode
resonance.

---

## 11. Pump injection and the detuning convention

The CW pump is injected each step at the input bus, with amplitude
$p_{\mathrm{in}} = \sqrt{\mathrm{pump\_power}}$ and the bus coupling $i\kappa$ already
carried in $\mathrm{src\_vals}$. The detuning is the round-trip phase, identical to the
frequency variable $\omega$ of `solve_one`:

$$\delta = 2\pi (\omega/2\pi)$$

(the dialog axis is $\omega/2\pi$ in FSR units). The solver uses $e^{+i\delta}$, so
$\delta>0$ is blue (pump above the cold resonance) and $\delta<0$ is red. As intracavity
power rises, the Kerr phase shifts the effective resonance toward negative $\delta$: the
resonance leans red and becomes bistable, and cavity solitons sit on the red-detuned upper
branch.

---

## 12. Bistability, sweep direction, and the soliton step

The Kerr-shifted resonance is multivalued (an S-shaped intracavity-power-versus-$\delta$
curve), so the cavity is **bistable** over a range of red detunings: a low-power CW branch
and a high-power branch coexist, with solitons on the upper branch. To reach them you start
blue (single-valued CW), tune down across the resonance where modulation instability fills
the cavity, and continue red until the field collapses onto one or a few solitons. The comb
power versus $\delta$ then shows a **soliton step** — a plateau that drops in discrete jumps
as solitons are lost. Hence the sweep direction is high to low $\delta$ (blue to red). The
solver always sweeps descending and carries the field forward, from
$\max(\delta_{\mathrm{start}}, \delta_{\mathrm{end}})$ down to
$\min(\delta_{\mathrm{start}}, \delta_{\mathrm{end}})$. (Verified: blue to red reaches
deeper into the red comb region with far richer combs than red to blue, confirming
branch-following.)

---

## 13. Observables (read at the drop port)

Everything plotted is read at the drop ring:

$$s_{\mathrm{drop}}(\mu) = i\kappa p_d(\mu) E[\mathrm{drop},\mu]$$

with $p_d(\mu) = e^{i \omega_\mu dz_s - \alpha dz_s/2}$. The quantity
$s_{\mathrm{drop}}$ is already the drop-port comb spectrum (mode domain), from which

$$P_{\mathrm{pump}} = |s_{\mathrm{drop}}(0)|^2$$

$$P_{\mathrm{comb}} = \sum_{\mu \ne 0} |s_{\mathrm{drop}}(\mu)|^2$$

and the comb spectrum is displayed as $|s_{\mathrm{drop}}(\mu)|^2$ in dB.
$P_{\mathrm{comb}}$ is time-averaged over the last 30 percent of round trips to smooth
breather flicker. At low power $P_{\mathrm{pump}}$ versus $\delta$ traces the cold drop
spectrum; at high power $P_{\mathrm{comb}}$ grows and the pump line tilts (the Kerr shift).

> Domain note: $s_{\mathrm{drop}}$ is a spectrum indexed by $\mu$. The pump line is
> $|s_{\mathrm{drop}}(0)|^2$ (not a time-mean) and the spectrum is
> $|s_{\mathrm{drop}}(\mu)|^2$ directly (no further FFT). Treating it as a time field would
> read the pump as $P_{\mathrm{pump}}/N^2$ and the comb as $P_{\mathrm{comb}}/N$ — a
> spurious two-orders-of-magnitude gap.

---

## 14. Modulation-instability seeding

A pure CW field ($\mu=0$ only) is a fixed point: Kerr on a flat field is a global phase and
makes no sidebands. Real combs nucleate when modulation instability amplifies noise, so the
solver injects a small noise field on the **first detuning only** ($i=0$):

$$E \to E + \sigma (\mathcal{N} + i \mathcal{N})$$

with $\sigma = \mathrm{seed\_noise}$ and $\mathcal{N}$ standard normal. Once the comb
nucleates, the carried-forward field rides its own branch — no further seeding needed.
Seeding every detuning re-nucleates the comb and washes out the hysteresis; this stands in
for the spontaneous noise that seeds combs in real devices.

---

## 15. The gamma = 0 validation

With the Kerr off, only $\mu=0$ is ever populated and the map collapses to

$$E \leftarrow R(\delta) E + s$$

the linear TMM. The pump line versus $\delta$ must then equal the exact drop spectrum from
`solve_one`. Verified numerically: on a 4x4 IQH lattice and on the single ring, the
$\gamma=0$ pump-line sweep reproduces `solve_one` with **correlation 1.0000** — identical
peak positions, heights, and linewidths. This guarantees the nonlinear solver drives the
same machine as the linear app; the band emerges from the real geometry, not from any
fitted parameter.

---

## 16. Slow-time analysis (dynamics at a fixed detuning)

The sweep gives one *settled* comb per detuning, but not how that state behaves over time.
The slow-time analysis takes the settled field at a selected detuning, **holds $\delta$
fixed**, and evolves the round-trip map for a chosen number of round trips
(`slow_time_evolve`, the same Strang machinery), recording the intracavity field each round
trip. This reveals whether the state is a stationary soliton, a breather, or chaotic.

The initial condition is the **carried-forward settled field** at that detuning, not a fresh
start — that matters, because starting from noise could settle onto a different branch (CW
instead of the soliton). The settled field of every detuning is therefore kept during the
sweep (`self._fields`, one `complex64` snapshot per detuning), so slow-time is instant.

Two views come out of the same evolution:

**The 2D map** — the drop-ring fast-time waveform $|E(\tau)|^2$ with fast time $\tau$ on the
horizontal axis and round-trip number (slow time) on the vertical. A stationary soliton is a
straight bright line; a breather oscillates; chaos is turbulent.

**The ring movie** — the fast time $\tau$ *is* the position $z$ around the ring (one
round-trip period = one lap), so the waveform maps onto the ring perimeter. Each site ring is
drawn as a rounded square by the **identical linear-app renderer** (`plot_field_distribution`
/ `_draw_ring`, via `plot_waveform_distribution`) and coloured by its own waveform: a soliton
is one bright arc, Turing rolls are several, animated over round trips with a slider and
MP4/GIF export. Mapping the *fast-time waveform* onto the perimeter (not the per-grid total
power, which is ~uniform because the pulse circulates rigidly) is what makes the structure
visible. For a lattice each ring carries its own waveform, so a super-soliton circulating the
edge appears as the bright arc hopping ring to ring.

`slow_time_evolve` returns two arrays: the drop waveform `wf[n_rt, N]` (the 2D map) and the
per-ring waveform `wfr[n_rt, n_rings, N]` (the ring movie).

---

## 17. Software and workflow

- **Threaded sweep.** The detuning sweep runs off the GUI thread in a
  `NonlinearSweepWorker(QThread)`. Each detuning emits a signal carrying the comb spectrum,
  pump line, soliton metric and settled field; a main-thread slot does the plotting and
  storage. NumPy releases the GIL during the FFTs and matmuls, so the window stays fully
  responsive (drag, click, Stop) during a run, with a progress bar and an estimated finish
  time. Stop aborts at the next detuning boundary, keeping every comb collected so far.

- **Save / load.** A run saves to one `.npz`: all dialog parameters plus the lattice **build
  parameters** (so the exact template can be rebuilt), the drop spectrum, the sweep curves,
  the per-detuning comb spectra and settled fields, and the last slow-time result. Loading
  rebuilds the template from the build parameters and replicates every plot — verified
  bit-identical to the original run.

---

## 18. Parameter mapping (dialog knobs to symbols)

| Dialog field            | Symbol         | Role |
|-------------------------|----------------|------|
| Number of modes         | $N$            | comb lines (fast-time resolution) |
| Pump power              | $\mathrm{pump\_power}$ | pump field $p_{\mathrm{in}} = \sqrt{\mathrm{pump\_power}}$ |
| Kerr $\gamma$           | $\gamma$       | Kerr coefficient (self-phase modulation) |
| Dispersion $D_2$        | $D_2$          | second-order dispersion; $s=-1$ if $D_2 \ge 0$ (anomalous) |
| Detuning start / end    | $\delta$ range | swept high to low (blue to red) |
| Detuning step           | $d\delta$      | sweep resolution |
| Round trips / detuning  | $n_{\mathrm{rt}}$ | iterations to settle at each $\delta$ |

Cavity and lattice parameters come from the main window (the same ones the linear TMM
uses): $\kappa$ (`spn_kex`), $\alpha$ (`spn_alpha`), $\beta_0\eta$ (`spn_beta0eta`),
lattice size, flux, and the template itself. The inter-ring coupling is not re-entered — it
is already in the propagator $R$.

---

## 19. Known limitations and numerical artifacts

- **Chaos vs. soliton.** Whether the comb is a Turing pattern, chaos, a breather, or a clean
  soliton depends on the operating point $(\mathrm{pump}, \delta, D_2)$ — a tuning question,
  not a solver error.

- **Kerr split-step cost.** The default distributed Strang scheme does $\approx N_s$ FFT-pairs
  per round trip; the lumped option (`kerr_lumped = True`) does one, trading second-order
  accuracy for speed on long runs or quick previews.

- **Convergence.** For very low loss the cavity is high-Q; increase $n_{\mathrm{rt}}$ if the
  low-power pump line does not fully resolve the band.

---

## 20. References

- S. D. Hashemi and S. Mittal, *Floquet topological dissipative Kerr solitons and
  incommensurate frequency combs*, Nat. Commun. (2024).
- S. D. Hashemi and S. Mittal, *Reconfigurable non-Hermitian soliton combs using
  dissipative couplings and topological windings*, Sci. Adv. (2025).
- C. J. Flower et al., *Observation of topological frequency combs*, Science (2024).
- L. Xu et al., *On-chip multi-timescale spatiotemporal optical synchronization*,
  Sci. Adv. (2025).
- A. Ekström, *Modelling of optical frequency comb generation in microresonators*
  (MSc thesis, 2022).
- K. Ikeda, *Multiple-valued stationary state and its instability of the transmitted light
  by a ring cavity system*, Opt. Commun. (1979).
- T. Herr et al., *Temporal solitons in optical microresonators*, Nat. Photonics (2014).

See `README.md` for the linear TMM theory that this solver builds on.
