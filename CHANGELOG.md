# Changelog

## 5.1.0 — 2026-09-10

- Fixed partial-year NAO parsing; the current snapshot includes August 2026.
- Explicitly selected and refreshed ERA5. Incomplete monthly variables are withheld; incomplete regional years are excluded from extreme-value analysis.
- Removed unsupported baseline substitutes; corrected tied ranks, calendar-aware smoothing and lag matching, selected-test q-values and wet-day threshold sensitivity.
- Separated rare-event units and displayed bootstrap intervals. Added source dates, a glossary, three reading paths and expandable technical sections, retaining the existing layout.
- Labelled the 1,000-borehole layer as a bounded preview. Source publication delays remain visible.
- Updated analysis downloads and checksums. Offline annual indices, trends and deterministic bootstrap results exactly matched the published snapshot. Import regression tests and browser interactions passed.
- The monthly update timer remains active. No new Zenodo archive was published.

[Website](https://sebastienspiess.ch/commissions/bristol-climate/) · [Source](https://github.com/sebastiengeorgesspiess-source/bristol-climate-ground-conditions)


## Earlier releases

The linked Zenodo 5.0.0 archive is historical; this repository starts with the corrected 5.1.0 source.
