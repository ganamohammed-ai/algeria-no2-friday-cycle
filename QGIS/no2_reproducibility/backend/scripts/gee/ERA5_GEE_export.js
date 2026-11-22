// ============================================================================
// ERA5 hourly fields for the NO2 / ventilation study, exported from Earth Engine.
// Replaces the Copernicus CDS download (which is queued for days).
//
// Datasets (both 0.25 deg, hourly, 1940 -> present):
//   ECMWF/ERA5/HOURLY                   single levels : 10 m wind, boundary layer height (BLH),
//                                        2 m temperature, surface pressure, solar radiation
//   ECMWF/ERA5/HOURLY_PRESSURE_LEVELS   u and v wind at the pressure levels (for the layer wind)
//
// Only the 12:00 and 13:00 UTC images are used (the satellite passes at about 12:30 UTC).
// For every site the script samples a 3 x 3 set of points (offsets -0.25, 0, +0.25 degrees), so that
// the four ERA5 cells around the site are always included; the python script then interpolates
// bilinearly, exactly as it does for the CDS files. The cell centre of each sample comes with the
// data (columns latitude / longitude).
//
// HOW TO RUN
//   1. Paste into code.earthengine.google.com and press Run.
//   2. Open the Tasks tab: 7 tasks appear (ERA5X_SL_2019 ... ERA5X_SL_2025). Press RUN on each.
//   3. The CSV files appear in your Google Drive folder DRIVE_FOLDER. Copy the 7 files next to the
//      older ERA5_SL_*.csv files and tell Claude.
// ============================================================================

var SITES = {
  algiers:        [36.75,  3.06],   // [lat, lon]
  oran:           [35.70, -0.63],
  arzew:          [35.85, -0.31],
  constantine:    [36.37,  6.61],
  setif:          [36.19,  5.41],
  hassi_messaoud: [31.68,  6.07]
};
var YEARS = [2019, 2020, 2021, 2022, 2023, 2024, 2025];
var DO_SL = true;
// The pressure-level collection has no images for 2019-2025 in Earth Engine (checked), so it is off.
var DO_PL = false;
var DRIVE_FOLDER = 'ERA5_GEE';

var SL = ee.ImageCollection('ECMWF/ERA5/HOURLY');
var PL = ee.ImageCollection('ECMWF/ERA5/HOURLY_PRESSURE_LEVELS');
var SL_BANDS = ['u_component_of_wind_10m', 'v_component_of_wind_10m', 'boundary_layer_height',
                'temperature_2m', 'surface_pressure', 'surface_solar_radiation_downwards',
                // wind higher up (the single-level collection also carries these):
                'u_component_of_wind_100m', 'v_component_of_wind_100m',
                'u_component_of_wind_850hPa', 'v_component_of_wind_850hPa',
                'u_component_of_wind_500hPa', 'v_component_of_wind_500hPa'];
var LEVELS = [1000, 975, 950, 925, 900, 850, 500];
var PL_BANDS = [];
LEVELS.forEach(function (p) {
  PL_BANDS.push('u_component_of_wind_' + p + 'hPa');
  PL_BANDS.push('v_component_of_wind_' + p + 'hPa');
});
var SCALE = 27830;
var OFFS = [-0.25, 0, 0.25];

// sample points: 9 per site
var feats = [];
Object.keys(SITES).forEach(function (name) {
  var la = SITES[name][0], lo = SITES[name][1];
  OFFS.forEach(function (dy) {
    OFFS.forEach(function (dx) {
      feats.push(ee.Feature(ee.Geometry.Point([lo + dx, la + dy]), {site: name, slat: la, slon: lo}));
    });
  });
});
var PTS = ee.FeatureCollection(feats);

function extract(col, bands, year) {
  var c = col.filterDate(year + '-01-01', (year + 1) + '-01-01')
             .filter(ee.Filter.calendarRange(12, 13, 'hour'))
             .select(bands);
  var fc = c.map(function (img) {
    var t = img.date().format('YYYY-MM-dd HH');
    var im = img.addBands(ee.Image.pixelLonLat());
    return im.reduceRegions({collection: PTS, reducer: ee.Reducer.first(), scale: SCALE})
             .map(function (f) { return f.set('time_utc', t).setGeometry(null); });
  }).flatten();
  return fc;
}

var COLS_HEAD = ['site', 'slat', 'slon', 'time_utc', 'latitude', 'longitude'];

YEARS.forEach(function (y) {
  if (DO_SL) {
    Export.table.toDrive({
      collection: extract(SL, SL_BANDS, y),
      description: 'ERA5X_SL_' + y,
      folder: DRIVE_FOLDER,
      fileFormat: 'CSV',
      selectors: COLS_HEAD.concat(SL_BANDS)
    });
  }
  if (DO_PL) {
    Export.table.toDrive({
      collection: extract(PL, PL_BANDS, y),
      description: 'ERA5_PL_' + y,
      folder: DRIVE_FOLDER,
      fileFormat: 'CSV',
      selectors: COLS_HEAD.concat(PL_BANDS)
    });
  }
});

// quick look: one image, to check that the band names exist (see the Console)
print('SL first image bands', SL.first().bandNames().size());
print('test point values 2022-01-15 12:00 (algiers)',
      SL.filterDate('2022-01-15T12:00:00', '2022-01-15T13:00:00').first()
        .select(SL_BANDS).reduceRegion(ee.Reducer.first(), ee.Geometry.Point([3.06, 36.75]), SCALE));
