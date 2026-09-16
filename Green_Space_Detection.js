// ================================================================
// DELHI GREEN SPACE DETECTION
// Sentinel-2 + Sentinel-1 + Landsat 8/9
// Output:
//   1 = Green Space
//   0 = Non-Green Space
// ================================================================


// ---------------------------------------------------------------
// 1. DELHI BOUNDARY
// ---------------------------------------------------------------

var delhi = ee.FeatureCollection('FAO/GAUL/2015/level1')
  .filter(ee.Filter.eq('ADM1_NAME', 'Delhi'));

var roi = delhi.geometry();

Map.centerObject(roi, 10);

Map.addLayer(
  delhi,
  {color: 'red'},
  'Delhi Boundary',
  false
);


// ---------------------------------------------------------------
// 2. ANALYSIS PERIOD
// ---------------------------------------------------------------

// You can change these dates.
// For better vegetation detection, use a full growing/annual period.

var startDate = '2025-01-01';
var endDate   = '2025-12-31';


// ---------------------------------------------------------------
// 3. SENTINEL-2 SURFACE REFLECTANCE
// ---------------------------------------------------------------

function maskS2(image) {

  // QA60 cloud mask
  var qa = image.select('QA60');

  var cloudBitMask  = 1 << 10;
  var cirrusBitMask = 1 << 11;

  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
    .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image
    .updateMask(mask)
    .divide(10000)
    .copyProperties(image, ['system:time_start']);
}


var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(roi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', 10))
  .map(maskS2);

print('Sentinel-2 images:', s2.size());


// ---------------------------------------------------------------
// 4. SENTINEL-2 MEDIAN COMPOSITE
// ---------------------------------------------------------------

var s2Composite = s2.median().clip(roi);


// ---------------------------------------------------------------
// 5. SENTINEL-2 VEGETATION INDICES
// ---------------------------------------------------------------

// NDVI
var ndviS2 = s2Composite
  .normalizedDifference(['B8', 'B4'])
  .rename('NDVI_S2');

// EVI
var eviS2 = s2Composite.expression(
  '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))',
  {
    'NIR': s2Composite.select('B8'),
    'RED': s2Composite.select('B4'),
    'BLUE': s2Composite.select('B2')
  }
).rename('EVI_S2');

// NDMI
var ndmiS2 = s2Composite
  .normalizedDifference(['B8', 'B11'])
  .rename('NDMI_S2');

// NDBI
var ndbiS2 = s2Composite
  .normalizedDifference(['B11', 'B8'])
  .rename('NDBI_S2');

// MNDWI
var mndwiS2 = s2Composite
  .normalizedDifference(['B3', 'B11'])
  .rename('MNDWI_S2');


// ---------------------------------------------------------------
// 6. SENTINEL-1 SAR
// ---------------------------------------------------------------

var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(roi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(
    ee.Filter.listContains(
      'transmitterReceiverPolarisation',
      'VV'
    )
  )
  .filter(
    ee.Filter.listContains(
      'transmitterReceiverPolarisation',
      'VH'
    )
  )
  .select(['VV', 'VH']);

print('Sentinel-1 images:', s1.size());


// Median SAR composite
var s1Composite = s1.median().clip(roi);


// VV and VH
var VV = s1Composite.select('VV').rename('VV');
var VH = s1Composite.select('VH').rename('VH');


// VH/VV ratio
var vhvv = VH
  .subtract(VV)
  .rename('VH_VV');


// ---------------------------------------------------------------
// 7. LANDSAT 8 + LANDSAT 9
// ---------------------------------------------------------------

// Landsat 8
var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2');

// Landsat 9
var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2');


// Cloud mask + scaling
function maskLandsat(image) {

  var qa = image.select('QA_PIXEL');

  // Cloud shadow
  var cloudShadow = 1 << 4;

  // Snow
  var snow = 1 << 5;

  // Cloud
  var cloud = 1 << 3;

  // Cirrus
  var cirrus = 1 << 2;

  var mask = qa.bitwiseAnd(cloudShadow).eq(0)
    .and(qa.bitwiseAnd(snow).eq(0))
    .and(qa.bitwiseAnd(cloud).eq(0))
    .and(qa.bitwiseAnd(cirrus).eq(0));

  // Surface reflectance scaling
  var optical = image
    .select([
      'SR_B2',
      'SR_B3',
      'SR_B4',
      'SR_B5',
      'SR_B6',
      'SR_B7'
    ])
    .multiply(0.0000275)
    .add(-0.2);

  return image
    .addBands(optical, null, true)
    .updateMask(mask)
    .copyProperties(image, ['system:time_start']);
}


l8 = l8
  .filterBounds(roi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lte('CLOUD_COVER', 10))
  .map(maskLandsat);


l9 = l9
  .filterBounds(roi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lte('CLOUD_COVER', 10))
  .map(maskLandsat);


// Merge Landsat 8 + 9
var landsat = l8.merge(l9);

print('Landsat images:', landsat.size());


// ---------------------------------------------------------------
// 8. LANDSAT COMPOSITE
// ---------------------------------------------------------------

var landsatComposite = landsat
  .median()
  .clip(roi);


// ---------------------------------------------------------------
// 9. LANDSAT VEGETATION INDICES
// ---------------------------------------------------------------

// Landsat NDVI
var ndviLandsat = landsatComposite
  .normalizedDifference(['SR_B5', 'SR_B4'])
  .rename('NDVI_Landsat');


// Landsat EVI
var eviLandsat = landsatComposite.expression(
  '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))',
  {
    'NIR': landsatComposite.select('SR_B5'),
    'RED': landsatComposite.select('SR_B4'),
    'BLUE': landsatComposite.select('SR_B2')
  }
).rename('EVI_Landsat');


// ---------------------------------------------------------------
// 10. COMBINE OPTICAL + RADAR FEATURES
// ---------------------------------------------------------------

// Reproject Landsat to Sentinel-2 10 m grid
var ndviLandsat10m = ndviLandsat
  .resample('bilinear')
  .reproject({
    crs: ndviS2.projection(),
    scale: 10
  });

var eviLandsat10m = eviLandsat
  .resample('bilinear')
  .reproject({
    crs: ndviS2.projection(),
    scale: 10
  });


// Feature stack
var features = ee.Image.cat([
  ndviS2,
  eviS2,
  ndmiS2,
  ndbiS2,
  mndwiS2,

  ndviLandsat10m,
  eviLandsat10m,

  VV,
  VH,
  vhvv
]).clip(roi);


print('Feature stack:', features);


// ---------------------------------------------------------------
// 11. VEGETATION DETECTION RULE
// ---------------------------------------------------------------

// Main vegetation conditions
//
// NDVI > 0.30
// EVI  > 0.15
//
// Additional conditions remove:
//   water
//   strongly built-up surfaces

var greenSpace = ndviS2.gt(0.30)
  .and(eviS2.gt(0.15))
  .and(ndviLandsat10m.gt(0.25))
  .and(ndbiS2.lt(0.25))
  .and(mndwiS2.lt(0.20))
  .rename('Green_Space');


// Convert to 0/1 raster
greenSpace = greenSpace
  .selfMask()
  .uint8();


// ---------------------------------------------------------------
// 12. REMOVE SMALL NOISY PATCHES
// ---------------------------------------------------------------

var connectedPixels = greenSpace
  .connectedPixelCount(100, true);

var greenSpaceClean = greenSpace
  .updateMask(connectedPixels.gte(9))
  .rename('Green_Space_Clean');


// ---------------------------------------------------------------
// 13. DENSE VEGETATION / TREE-LIKE AREAS
// ---------------------------------------------------------------

var denseVegetation = ndviS2.gt(0.50)
  .and(eviS2.gt(0.30))
  .and(ndviLandsat10m.gt(0.40))
  .and(ndmiS2.gt(0.05))
  .and(ndbiS2.lt(0.20))
  .rename('Dense_Vegetation');

denseVegetation = denseVegetation
  .selfMask()
  .uint8();


// ---------------------------------------------------------------
// 14. MODERATE / LOW VEGETATION
// ---------------------------------------------------------------

var lowModerateVegetation = ndviS2.gt(0.30)
  .and(ndviS2.lte(0.50))
  .and(eviS2.gt(0.15))
  .and(ndbiS2.lt(0.25))
  .and(mndwiS2.lt(0.20))
  .rename('Moderate_Vegetation');


// ---------------------------------------------------------------
// 15. VISUALIZATION
// ---------------------------------------------------------------

// Sentinel-2 RGB
Map.addLayer(
  s2Composite,
  {
    bands: ['B4', 'B3', 'B2'],
    min: 0,
    max: 0.3
  },
  'Sentinel-2 RGB'
);


// Sentinel-2 NDVI
Map.addLayer(
  ndviS2,
  {
    min: 0,
    max: 0.8,
    palette: [
      'brown',
      'yellow',
      'lightgreen',
      'green',
      'darkgreen'
    ]
  },
  'Sentinel-2 NDVI',
  false
);


// Sentinel-1 VH
Map.addLayer(
  VH,
  {
    min: -25,
    max: 0
  },
  'Sentinel-1 VH',
  false
);


// Green space
Map.addLayer(
  greenSpaceClean,
  {
    palette: ['00A000']
  },
  'GREEN SPACE - FINAL'
);


// Dense vegetation
Map.addLayer(
  denseVegetation,
  {
    palette: ['006400']
  },
  'DENSE VEGETATION / TREES'
);


// Moderate vegetation
Map.addLayer(
  lowModerateVegetation.selfMask(),
  {
    palette: ['7CFC00']
  },
  'MODERATE VEGETATION',
  false
);


// ---------------------------------------------------------------
// 16. GREEN SPACE AREA
// ---------------------------------------------------------------

var greenArea = greenSpaceClean
  .multiply(ee.Image.pixelArea())
  .reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: roi,
    scale: 10,
    maxPixels: 1e13
  });

print(
  'Green Space Area (m²):',
  greenArea
);

print(
  'Green Space Area (hectares):',
  ee.Number(greenArea.get('Green_Space_Clean'))
    .divide(10000)
);

print(
  'Green Space Area (km²):',
  ee.Number(greenArea.get('Green_Space_Clean'))
    .divide(1e6)
);


// ---------------------------------------------------------------
// 17. EXPORT GREEN SPACE RASTER
// ---------------------------------------------------------------

Export.image.toDrive({
  image: greenSpaceClean,
  description: 'Delhi_Green_Space_10m',
  folder: 'GEE_Delhi_GreenSpace',
  fileNamePrefix: 'Delhi_Green_Space_10m',
  region: roi,
  scale: 10,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});


// ---------------------------------------------------------------
// 18. EXPORT DENSE VEGETATION
// ---------------------------------------------------------------

Export.image.toDrive({
  image: denseVegetation,
  description: 'Delhi_Dense_Vegetation_10m',
  folder: 'GEE_Delhi_GreenSpace',
  fileNamePrefix: 'Delhi_Dense_Vegetation_10m',
  region: roi,
  scale: 10,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});


// ---------------------------------------------------------------
// 19. EXPORT NDVI
// ---------------------------------------------------------------

Export.image.toDrive({
  image: ndviS2,
  description: 'Delhi_NDVI_Sentinel2',
  folder: 'GEE_Delhi_GreenSpace',
  fileNamePrefix: 'Delhi_NDVI_Sentinel2',
  region: roi,
  scale: 10,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
