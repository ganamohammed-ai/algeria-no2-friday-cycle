// ============================================================================
// Daily TROPOMI NO2 stacks for the two CONTROL sites (same format as the NO2_daily_<site>_<year>.tif files).
// One band per day, 0.01 deg grid, box +-1.5 deg latitude and +-1.5/cos(lat) longitude, umol/m2, -9999 = no data,
// pixels with cloud_fraction >= 0.3 masked, several orbits of a day averaged.
// Sites were chosen by select_control.py --amended from static layers only (no NO2 was used).
// HOW TO RUN: paste in the Code Editor, Run, then press RUN on each of the 28 tasks (NO2_daily_<site>_<year>).
// Copy the 14 .tif files from Drive folder NO2_GEE_CTLX into  article NO2\control\  and tell Claude.
// ============================================================================
var SITES = {ctl_x01: [35.75, 2.25], ctl_x02: [35.75, 3.5], ctl_x03: [36.0, 8.0], ctl_x04: [35.25, 7.0]};   // [lat, lon]
var YEARS = [2019, 2020, 2021, 2022, 2023, 2024, 2025];
var FIRST_DAY = '2019-08-07';
var NO2 = 'tropospheric_NO2_column_number_density';
var col = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2').select([NO2, 'cloud_fraction']);

function dayImage(d) {
  d = ee.Date(d);
  var c = col.filterDate(d, d.advance(1, 'day')).map(function (im) {
    return im.select(NO2).updateMask(im.select('cloud_fraction').lt(0.3));
  });
  var im = ee.Image(ee.Algorithms.If(c.size().gt(0),
      c.mean().multiply(1e6).unmask(-9999), ee.Image.constant(-9999)));
  return im.rename('no2').toFloat();
}

Object.keys(SITES).forEach(function (name) {
  var la = SITES[name][0], lo = SITES[name][1];
  var hwLat = 1.5, hwLon = Math.round(1.5 / Math.cos(la * Math.PI / 180) * 100) / 100;
  var lon0 = Math.round((lo - hwLon) * 100) / 100, lat0 = Math.round((la + hwLat) * 100) / 100;
  var lon1 = Math.round((lo + hwLon) * 100) / 100, lat1 = Math.round((la - hwLat) * 100) / 100;
  var rect = ee.Geometry.Rectangle([lon0, lat1, lon1, lat0], 'EPSG:4326', false);
  YEARS.forEach(function (y) {
    var start = ee.Date(y === 2019 ? FIRST_DAY : y + '-01-01');
    var end = ee.Date((y + 1) + '-01-01');
    var n = end.difference(start, 'day');
    var days = ee.List.sequence(0, n.subtract(1)).map(function (i) { return start.advance(i, 'day'); });
    var stack = ee.ImageCollection.fromImages(days.map(dayImage)).toBands();
    Export.image.toDrive({
      image: stack, description: 'NO2_daily_' + name + '_' + y, folder: 'NO2_GEE_CTLX', fileNamePrefix: 'NO2_daily_' + name + '_' + y,
      region: rect, crs: 'EPSG:4326', crsTransform: [0.01, 0, lon0, 0, -0.01, lat0], maxPixels: 1e10, fileFormat: 'GeoTIFF'
    });
  });
  print(name, 'box lon', lon0, lon1, 'lat', lat1, lat0);
});
