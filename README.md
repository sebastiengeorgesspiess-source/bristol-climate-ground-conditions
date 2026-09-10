# Bristol Climate & Ground Conditions — 5.1.0

Reviewed 10 September 2026. Regional observations, local stations, ERA5 model estimates and geological indexes remain separate.

## Reproduce annual extremes

Download the daily JSON and `bristol_extremes_reproducible.py`. Install numpy and scipy, then run:

    python bristol_extremes_reproducible.py bristol-climate-daily.json

The snapshot date fixes the complete-year cutoff. The script reproduces regional annual indices, exploratory trend diagnostics and deterministic stationary GEV bootstrap estimates. Each variable requires a complete year. Intervals are conditional on that model, not forecasts or a complete uncertainty assessment.

## Monthly comparisons

The interactive dashboard uses calendar-aware matching, 1991–2020 baselines with at least 24 valid years per calendar month, tied ranks, and a selected-test BH q-value within 15 smoothing/lag variants. Inference uses standardised values and an approximate lag-1 effective sample size. It does not correct every user-selected comparison or establish causality.

Use the dashboard CSV to export the selected transformed values and settings. `bristol_climate_analysis.py` is only a raw overlapping-month export helper, not an independent reproduction of interactive inference. Its notebook states this scope.

## Data updates and limits

The monthly server job runs on day 6 at 04:35 UTC with up to 30 minutes of scheduling delay. A successful retrieval does not imply every source is current. Observation dates are displayed separately. ERA5 is explicitly selected and refreshed; partial monthly variables are null. NOAA NAO partial-year rows are supported. Local station context requires at least 300 days, is not homogenised and is not used for regional inference. Groundwater daily means use available qualified readings.

The borehole layer is a labelled 1,000-record preview, not the full BGS catalogue. Source publication delays and access-gated datasets remain visible. No new Zenodo release is implied: the linked 5.0.0 archive is historical.

## Changes

- Fixed NAO partial-year parsing and explicitly selected ERA5.
- Excluded incomplete regional years and incomplete monthly variables.
- Removed invented baseline fallbacks; corrected tied ranks, calendar gaps and selected q-values.
- Made wet-day sensitivity thresholds effective.
- Separated return-level units and displayed bootstrap intervals.
- Added source freshness, glossary and short reading paths; retained the layout and charts.

Raw input is retained. Earlier code and snapshots were backed up before deployment.

## Map configuration

The live CARTO browser credential is intentionally removed from this public archive. Replace `YOUR_CARTO_BROWSER_KEY` only in your deployed copy with your own appropriately restricted browser key, or configure a compatible tile provider. Do not commit credentials. This redaction does not alter the live site.

## Repository scope

This is a source archive, not a standalone deployment of the multi-project production server. `templates/` and `static/` preserve the live dashboard code; shared website navigation, artwork, vendor assets, server configuration and upstream data are not bundled. The offline analysis scripts can be run separately with the public daily JSON. The importer uses public upstream services and the existing production data directory; adapt that path for a separate installation. Do not run it against a production directory without a backup.

Live dashboard: https://sebastienspiess.ch/commissions/bristol-climate/
Daily JSON: https://sebastienspiess.ch/commissions/download/bristol-daily

Tests: `python test_quality.py` and `node test_math.cjs`. The NAO regression fixture is dated September 2026.
