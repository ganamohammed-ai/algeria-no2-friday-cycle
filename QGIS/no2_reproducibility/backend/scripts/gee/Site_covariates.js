// ============================================================================
// Static site descriptors for the NO2 / ventilation article (one CSV, six rows):
//   population (WorldPop 2020) within 10 km and 50 km, mean night-light radiance 2019-2025,
//   mean elevation and its standard deviation (terrain complexity), built-up share (ESA WorldCover 2021).
// Run: paste in the Code Editor, press Run, then press RUN on the single task "Site_covariates".
// ============================================================================
var SITES = {
  algiers:        [36.75,  3.06],
  oran:           [35.70, -0.63],
  arzew:          [35.85, -0.31],
  constantine:    [36.37,  6.61],
  setif:          [36.19,  5.41],
  hassi_messaoud: [31.68,  6.07]
};
var pop = ee.ImageCollection('WorldPop/GP/100m/pop')
            .filter(ee.Filter.eq('country', 'DZA')).filter(ee.Filter.eq('year', 2020)).mosaic().select('population');
var lights = ee.ImageCollection('NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG')
               .filterDate('2019-01-01', '2026-01-01').select('avg_rad').mean();
var elev = ee.Image('USGS/SRTMGL1_003').select('elevation');
var built = ee.ImageCollection('ESA/WorldCover/v200').first().select('Map').eq(50).rename('built');

function stats(img, geom, reducer, scale) {
  return img.reduceRegion({reducer: reducer, geometry: geom, scale: scale, maxPixels: 1e10, bestEffort: true});
}
var feats = Object.keys(SITES).map(function (name) {
  var pt = ee.Geometry.Point([SITES[name][1], SITES[name][0]]);
  var d10 = pt.buffer(10000), d50 = pt.buffer(50000);
  var p10 = stats(pop, d10, ee.Reducer.sum(), 100).get('population');
  var p50 = stats(pop, d50, ee.Reducer.sum(), 100).get('population');
  var nl10 = stats(lights, d10, ee.Reducer.mean(), 500).get('avg_rad');
  var nl50 = stats(lights, d50, ee.Reducer.mean(), 500).get('avg_rad');
  var e10 = stats(elev, d10, ee.Reducer.mean().combine(ee.Reducer.stdDev(), '', true), 90);
  var e50 = stats(elev, d50, ee.Reducer.mean().combine(ee.Reducer.stdDev(), '', true), 90);
  var b10 = stats(built, d10, ee.Reducer.mean(), 100).get('built');
  var b50 = stats(built, d50, ee.Reducer.mean(), 100).get('built');
  return ee.Feature(null, {
    site: name, lat: SITES[name][0], lon: SITES[name][1],
    pop_10km: p10, pop_50km: p50,
    nightlight_10km: nl10, nightlight_50km: nl50,
    elev_mean_10km: e10.get('elevation_mean'), elev_sd_10km: e10.get('elevation_stdDev'),
    elev_mean_50km: e50.get('elevation_mean'), elev_sd_50km: e50.get('elevation_stdDev'),
    built_share_10km: b10, built_share_50km: b50
  });
});
var fc = ee.FeatureCollection(feats);
print('covariates', fc);
Export.table.toDrive({collection: fc, description: 'Site_covariates', folder: 'NO2_GEE', fileFormat: 'CSV'});
