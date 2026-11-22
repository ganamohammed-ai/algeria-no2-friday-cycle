// ============================================================================
// Daily UV Aerosol Index (dust / smoke flag) and cloud fraction at the six sites, 2019-08-07 .. 2025-12-31.
// Same product and cloud rule as the NO2 exports (COPERNICUS/S5P/OFFL/L3_NO2 carries the AAI band).
// Output: one CSV per year, one row per site and day:
//   uvai_inner10km, uvai_box (about 167 km radius), cloud_fraction_box, n_valid_uvai_pixels
// Run: paste in the Code Editor, press Run, then press RUN on the 7 tasks (S5P_AAI2_2019 ... 2025).
// Files go to the Drive folder DRIVE_FOLDER. Put them next to the other files in "article NO2".
// ============================================================================
var SITES = {
  algiers:        [36.75,  3.06],
  oran:           [35.70, -0.63],
  arzew:          [35.85, -0.31],
  constantine:    [36.37,  6.61],
  setif:          [36.19,  5.41],
  hassi_messaoud: [31.68,  6.07]
};
var YEARS = [2019, 2020, 2021, 2022, 2023, 2024, 2025];
var FIRST_DAY = '2019-08-07', LAST_DAY = '2025-12-31';
var HALF_DEG = 1.5, INNER_KM = 10, SCALE = 1113.2;
var DRIVE_FOLDER = 'NO2_GEE';
var COL = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2');

function yearRange(y) {
  var a = (y + '-01-01' < FIRST_DAY) ? FIRST_DAY : (y + '-01-01');
  var b = (y + '-12-31' > LAST_DAY) ? LAST_DAY : (y + '-12-31');
  var n = Math.round(((new Date(b + 'T00:00:00Z')) - (new Date(a + 'T00:00:00Z'))) / 86400000) + 1;
  return {start: a, n: n};
}

function siteTable(name, y) {
  var la = SITES[name][0], lo = SITES[name][1];
  var pt = ee.Geometry.Point([lo, la]);
  var inner = pt.buffer(INNER_KM * 1000);
  var halfLon = HALF_DEG / Math.cos(la * Math.PI / 180);
  var box = ee.Geometry.Rectangle([lo - halfLon, la - HALF_DEG, lo + halfLon, la + HALF_DEG]);
  var r = yearRange(y);
  return ee.List.sequence(0, r.n - 1).map(function (i) {
    var d = ee.Date(r.start).advance(i, 'day');
    var day = COL.filterBounds(box).filterDate(d, d.advance(1, 'day'))
                 .select(['absorbing_aerosol_index', 'cloud_fraction']);
    // a day without any image gives a fully masked image, so the numbers come out empty instead of failing
    var empty = ee.Image.constant([0, 0]).rename(['absorbing_aerosol_index', 'cloud_fraction'])
                  .updateMask(ee.Image.constant(0));
    var img = ee.Image(ee.Algorithms.If(day.size().gt(0), day.mean(), empty));
    var bx = ee.Dictionary(img.reduceRegion({reducer: ee.Reducer.mean().combine(ee.Reducer.count(), '', true),
                               geometry: box, scale: SCALE, maxPixels: 1e9}));
    var inn = ee.Dictionary(img.select('absorbing_aerosol_index').reduceRegion({reducer: ee.Reducer.mean(),
                               geometry: inner, scale: SCALE, maxPixels: 1e9}));
    return ee.Feature(null, {
      site: name, date: d.format('YYYY-MM-dd'),
      uvai_inner10km: inn.get('absorbing_aerosol_index', null),
      uvai_box: bx.get('absorbing_aerosol_index_mean', null),
      cloud_fraction_box: bx.get('cloud_fraction_mean', null),
      n_valid_uvai_pixels: bx.get('absorbing_aerosol_index_count', 0)
    });
  });
}

YEARS.forEach(function (y) {
  var lists = Object.keys(SITES).map(function (n) { return siteTable(n, y); });
  var all = ee.List(lists[0]);
  for (var k = 1; k < lists.length; k++) { all = all.cat(lists[k]); }
  Export.table.toDrive({
    collection: ee.FeatureCollection(all),
    description: 'S5P_AAI2_' + y,
    folder: DRIVE_FOLDER,
    fileFormat: 'CSV',
    selectors: ['site', 'date', 'uvai_inner10km', 'uvai_box', 'cloud_fraction_box', 'n_valid_uvai_pixels']
  });
});
print('7 tasks created: open the Tasks tab and press RUN on each.');
