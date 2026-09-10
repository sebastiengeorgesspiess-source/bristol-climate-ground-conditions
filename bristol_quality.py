"""Public metadata and validation for Bristol 5.1.0."""
from datetime import datetime, timezone

def complete_row(row):
    if 'value' in row:
        return row['value'] is not None
    values=[v for k,v in row.items() if k not in ('date','dateTime','valid_days','expected_days','n')]
    return bool(values) and all(v is not None for v in values)

def finish_document(doc, groundwater, storm_rain, flow):
    doc['schema_version'] = '5.1.0'
    doc['updated_at'] = datetime.now(timezone.utc).isoformat()
    doc['refresh'] = {'schedule': 'Monthly, on day 6, 04:35–05:05 UTC',
                      'meaning': 'Successful retrieval is not the observation date. Sources publish at different speeds.',
                      'failure_policy': 'Keep the last complete dashboard snapshot if a required source or calculation fails.'}
    coverage = {}
    for key, rows in doc['series'].items():
        dates = sorted(r.get('date') or r.get('dateTime', '')[:10] for r in rows)
        valid = sorted(r.get('date') or r.get('dateTime', '')[:10] for r in rows
                       if complete_row(r))
        coverage[key] = {'first': dates[0] if dates else None, 'last_received': dates[-1] if dates else None,
                         'last_value': valid[-1] if valid else None, 'records': len(rows),
                         'incomplete_periods': sum(not complete_row(r) for r in rows)}
    for key, packet in [('felton_groundwater_daily_mAOD', groundwater),
                        ('clifton_rain_15m_extremes_mm', storm_rain),
                        ('frenchay_flow_15m_extremes_m3s', flow)]:
        dates = sorted(x.get('date') or x.get('dateTime', '')[:10] for x in packet[0])
        coverage[key]['last_received'] = dates[-1] if dates else None
        coverage[key]['retrieval_limit_reached'] = len(packet[0]) >= 100000
    doc['coverage'] = coverage
    for key, collection in doc['geology'].items():
        if isinstance(collection, dict) and 'features' in collection:
            collection['returned'] = len(collection['features'])
            total = collection.get('numberMatched')
            collection['truncated'] = total is not None and int(total) > collection['returned']
    for source in doc['sources']:
        if source['publisher'] == 'ECMWF Copernicus via Open-Meteo':
            source['dataset'] = 'ERA5 historical reanalysis for Bristol grid point · model explicitly selected'
            source['caveat'] = 'Consistent ERA5 selection, not a weather-model blend. Grid estimates are not station measurements. Incomplete monthly variables are withheld.'
    doc['variables']['bristol_reanalysis']['label'] = 'ERA5 Bristol grid cell'
    doc['variables']['bristol_reanalysis']['missing'] = 'Each monthly variable requires every calendar day; otherwise null. Daily counts are retained.'
    for artifact in doc['artifacts']:
        if artifact['id'] == 'bristol_reanalysis':
            artifact['version'] = 'Explicit ERA5 daily API responses; refreshed this run'
    doc['methods'].update({
        'aggregation': 'Complete calendar months only for local monthly totals and averages; incomplete periods retain counts and a null value.',
        'annual_completeness': 'Regional extremes: completed years only, with every calendar day for the relevant variable. Local station context uses at least 300 days and is not homogenised. Raw daily input remains available.',
        'seasonality': 'Comparisons use 1991–2020 calendar-month baselines only where each month has at least 24 valid years. No substitute baseline is invented.',
        'inference': 'Exploratory Pearson and tie-corrected Spearman associations. Fisher intervals and p-values use a lag-1 effective-sample approximation, not a general dependence correction.',
        'multiple_testing': 'The selected comparison has its own BH-adjusted q-value within 15 predeclared smoothing/lag variants. This is not a correction for all user-selected series or date ranges.',
        'extremes': 'Stationary GEV fits use complete years. Bootstrap intervals are model-conditional, not forecasts; climate trends may violate stationarity.',
        'sensitivity': 'Wet-day threshold changes the reference wet-day percentile. RX1day trend uses the same full-year analysis period and is independent of that threshold.',
    })
    doc['limitations'].extend([
        'Borehole map is a labelled bounded preview, not the complete BGS catalogue.',
        'Annual trend tests and change points are exploratory; their p-values do not account for all serial dependence or multiple searches.',
        'The latest timestamp of one source must not be read as freshness of every series.',
    ])
    # Structural regression gates: fail before replacing the public snapshot.
    assert doc['series']['nao_index'], 'NAO empty'
    for row in doc['extremes']['swep_annual']:
        assert row['year'] < datetime.now(timezone.utc).year, 'Partial year in inference'
    assert all(r.get('valid_days') for r in doc['series']['bristol_reanalysis']), 'Old unvalidated reanalysis cache'
