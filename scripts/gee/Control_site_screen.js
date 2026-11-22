// ============================================================================
// Control-site screening for the NO2 / ventilation article.  NO NO2 DATA IS USED HERE.
// A grid of candidate points (0.25 deg) over northern Algeria is described with static layers only:
//   built-up share (ESA WorldCover 2021), population (WorldPop 2020), night-light (VIIRS 2019-2025),
//   elevation and its standard deviation (SRTM), and distance to the six study sites.
// The control site is then picked from this table by the fixed rule in select_control.py,
// so the choice cannot depend on the NO2 that the site shows.
// Run: paste in the Code Editor, press Run, then press RUN on the single task "Control_candidates".
// ============================================================================
var SITES = [[36.75, 3.06], [35.70, -0.63], [35.85, -0.31], [36.37, 6.61], [36.19, 5.41], [31.68, 6.07]];
var siteGeoms = SITES.map(function (s) { return ee.Geometry.Point([s[1], s[0]]); });
var pop = ee.ImageCollection('WorldPop/GP/100m/pop')
            .filter(ee.Filter.eq('country', 'DZA')).filter(ee.Filter.eq('year', 2020)).mosaic().select('population');
var lights = ee.ImageCollection('NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG')
               .filterDate('2019-01-01', '2026-01-01').select('avg_rad').mean();
var elev = ee.Image('USGS/SRTMGL1_003').select('elevation');
var built = ee.ImageCollection('ESA/WorldCover/v200').first().select('Map').eq(50).rename('built');
var water = ee.ImageCollection('ESA/WorldCover/v200').first().select('Map').eq(80).rename('water');

var pts = [];
for (var lat = 33.50; lat <= 37.001; lat += 0.25) {
  for (var lon = -1.50; lon <= 8.501; lon += 0.25) { pts.push(ee.Feature(ee.Geometry.Point([lon, lat]), {lat: lat, lon: lon})); }
}
var cand = ee.FeatureCollection(pts);

function red(img, buf, reducer, scale) {
  return function (f) {
    return img.reduceRegion({reducer: reducer, geometry: f.geometry().buffer(buf), scale: scale, maxPixels: 1e9, bestEffort: true});
  };
}
var out = cand.map(function (f) {
  var g = f.geometry();
  var d10 = g.buffer(10000), d50 = g.buffer(50000);
  var b10 = built.reduceRegion({reducer: ee.Reducer.mean(), geometry: d10, scale: 100, maxPixels: 1e9, bestEffort: true}).get('built');
  var b50 = built.reduceRegion({reducer: ee.Reducer.mean(), geometry: d50, scale: 500, maxPixels: 1e9, bestEffort: true}).get('built');
  var w10 = water.reduceRegion({reducer: ee.Reducer.mean(), geometry: d10, scale: 100, maxPixels: 1e9, bestEffort: true}).get('water');
  var p10 = pop.reduceRegion({reducer: ee.Reducer.sum(), geometry: d10, scale: 100, maxPixels: 1e9, bestEffort: true}).get('population');
  var l10 = lights.reduceRegion({reducer: ee.Reducer.mean(), geometry: d10, scale: 500, maxPixels: 1e9, bestEffort: true}).get('avg_rad');
  var l50 = lights.reduceRegion({reducer: ee.Reducer.mean(), geometry: d50, scale: 500, maxPixels: 1e9, bestEffort: true}).get('avg_rad');
  var e10 = elev.reduceRegion({reducer: ee.Reducer.mean().combine(ee.Reducer.stdDev(), '', true), geometry: d10, scale: 200, maxPixels: 1e9, bestEffort: true});
  var dmin = ee.List(siteGeoms.map(function (s) { return g.distance(s, 100); })).reduce(ee.Reducer.min());
  return f.set({built_10km: b10, built_50km: b50, water_10km: w10, pop_10km: p10, nightlight_10km: l10, nightlight_50km: l50,
                elev_mean_10km: e10.get('elevation_mean'), elev_sd_10km: e10.get('elevation_stdDev'),
                dist_nearest_site_km: ee.Number(dmin).divide(1000)});
});
print('first candidates', out.limit(5));
Export.table.toDrive({collection: out, description: 'Control_candidates', folder: 'NO2_GEE', fileFormat: 'CSV'});
