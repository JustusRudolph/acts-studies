# navigator_studies

Design work on track finding for **looping tracks** — low-$p_T$ particles whose
transverse radius `r_max = |c| + R` falls inside the tracking volume, so they
curve back towards the beamline instead of leaving the detector.

The main deliverable is a design note arguing for a **segment-and-merge**
architecture: stop the CKF at the radial turning point and merge the resulting
segments afterwards, rather than building a navigator that flies through the
turn. The supporting argument is geometric — consecutive radial turning points
are separated by exactly half a revolution of the track circle, so every segment
is monotonic in `r` and sits inside the regime the stock `Navigator` already
handles correctly.

## Layout

```
navigator_studies/
├── Makefile              build everything
├── README.md             this file
├── tex/                  sources (tracked)
│   ├── loopers.tex           the design note
│   └── fig-NN-*.tex          standalone TikZ figures
└── figures/              build output (gitignored)
    ├── loopers.pdf           <- the note
    ├── fig-NN-*.pdf
    └── fig-NN-*-1.png        (only after `make png`)
```

## Build

```sh
make              # figures/loopers.pdf, with figures
make figs         # the standalone figures only
make png          # PNG rasters (DPI=150; override with DPI=300)
make clean        # wipe figures/
```

Requires `pdflatex` with `standalone`, `tikz`, `newtxtext`, `booktabs`,
`enumitem`, `caption`, `hyperref`, `geometry`. `pdftoppm` is needed for
`make png`.

`figures/` matches the `*/figures/` rule in the top-level `.gitignore`, so build
products stay untracked — consistent with `geo_studies/` and `pgun_studies/`.

## Figures

| File | Shows |
|---|---|
| `fig-01-propagator-loop` | The three-way handshake. Navigator and CKF never call each other; `Propagator.ipp:99-193` mediates, and `state.navigation.currentSurface` is the whole contract. |
| `fig-02-tangent-vs-chord` | Tangent vs. chord vs. true arc at increasing curvature-per-step, and the degeneracy at the turning point. |
| `fig-03-turning-point` | The failure, in the transverse plane. At `T` the tangent is purely azimuthal, so candidates are scored along a ray that leaves the reachable region `r <= r_max`. |
| `fig-04-chord-lookbehind` | The retrospective sweep from `TryAllNavigator.hpp:340-346`: search backwards along the chord with `nearLimit = -|AB|`, `farLimit = 0`. |
| `fig-05-ckf-branch-tree` | The CKF depth-first branch walk and the `navigator.initialize()` calls at every backtrack (`CombinatorialKalmanFilter.hpp:471-538`). |
| `fig-06-policy-lifecycle` | `INavigationPolicy` call sites across a volume traversal, and why `isValid` defaulting to `true` means Gen3 never re-resolves inside a volume. |

## Proposed implementation order

1. **Turning-point `branchStopper`** — `StopAndKeep` on the sign flip of the
   transverse radial direction cosine. ~60 lines, no propagator or navigator
   changes.
2. **Directed continuation** — analytic extrapolation across the turn feeding a
   second CKF pass. Removes the need for inbound-capable seeding.
3. **Segment matching for orphans** — helix invariants `(c_x, c_y, R, tanλ, q)`
   for pre-selection, combined refit for confirmation.
4. **Acceptance widening** — only with χ² pull distributions in hand.
5. **A custom navigator** — only for whatever residual remains after 1–4.

## Key facts the figures and note encode

- `Surface::intersect` is **straight-line only**. There is no helix–surface
  intersection anywhere in ACTS Core.
- `Navigator::nextTarget(state, position, direction)` and Gen3's
  `NavigationArguments` carry no `q/p`, charge or field, so the navigation layer
  is curvature-blind by interface.
- Gen3 candidates are resolved once per volume and consumed by a monotonically
  increasing index (`Navigator.cpp:574-586`); `isValid` defaulting to `true`
  (`INavigationPolicy.hpp:265-271`) means no concrete policy in Core ever forces
  a refresh.
- Loop protection caps the path at `loopFraction * 2*pi*p/B` = half a helix turn
  by default (`PropagatorOptions.hpp:46-48`) — which is *exactly* one segment,
  so segmentation makes it a backstop rather than an antagonist.
- `maxSteps = 1000` is a budget for the entire branch-tree walk, and exhausting
  it discards every track found so far.
- Dynamic track columns are copied on branch split (`TrackProxy.hpp:551`), so
  per-branch state in a `branchStopper` is safe.

## Provenance

Line references are against ACTS commit `664c82a5` (19 August 2026), checked out
at `../../ACTS`. Appendix B of the note carries the full reference table; if the
checkout moves, that table and the figure annotations are what need refreshing.
