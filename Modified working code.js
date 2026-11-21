// JavaScript code to be implemented in Google Earth Engine(c) developed by H.A. Orengo
// to accompany the paper:

// Orengo, H.A. and Petrie, C.A. 2017. Large-scale, multi-temporal remote sensing of
// palaeo-river networks: a case study from northwest India and its implications for
// the Indus Civilisation. Remote Sensing, 9(7): 735; doi: 10.3390/rs9070735

//                ------------- ooo -------------

// TO EXECUTE THESE SCRIPTS PASTE THIS CODE INTO GOOGLE EARTH ENGINE CODE EDITOR AND PRESS RUN

//                ------------- ooo -------------

// The different algorithms are written here following the order adopted in the journal article:
// 1. NDVSI (Normalised Difference Vegetation Seasonality Index), line 102
// 2. SMTVI (Seasonal Multi-Temporal Vegetation Indices),         line 152
// 3. Wet Months Multi-Temporal Tasselled Cap Transformation,     line 236
// 4. Dry Months Multi-Temporal Tasselled Cap Transformation,     line 295
// 5. Wet Months Multi-Temporal Principal Component Analysis,     line 346
// 6. Dry Months Multi-Temporal Principal Component Analysis,     line 457

// For more information on how these scripts work and how to apply them refer to the text of the 
// article, which is freely available at: http://www.mdpi.com/2072-4292/9/7/735
// Suggestions and code improvementents are welcome! Please, contact Hector A. Orengo

// NOTES, READ CAREFULLY!
// 1. This is a modified version of the original script. GEE made some changes to the Landsat 5
//    collection that rendered previous code unusable. The changes can be summarised in:
//          a.  Use of the atmospherically corrected surface reflectance instead of TOA.
//          b.  Inclussion of a bit cloud mask instead of the Fmask algortihm used in the previous
//              TOA dataset.
//    This update uses Landsat Collection 2. SMTVI retains the original TOA basis.
// 1. Visualisation parameters for each layer can be adjusted in the Layers dialogue in the map
//    area of Goggle Earth Engine. This might be convenient when visualising at different scales.
// 2. To apply these analyses in other areas simply change the 'geometry' variable and the
//    'Map.setCenter' coordinates. To analyse the current map extent, use
//    var geometry = ee.Geometry.Rectangle(Map.getBounds()); in place of the polygon below.
//    The user can also draw his/her own polygon using the tool provided in the top left of Google
//    Earth Engine's map window. This will create a polygon named 'geometry' by default. The user
//    can then delete the variable called 'geometry' in lines 63-71 and the analysis will
//    be performed in the newly defined area. Make sure the central coordinates defined by
//    'Map.setCenter' (line 75) are also changed so the map window will zoom to these
//    coordinates when executing the scripts. 
// 3. The execution of the algorithms can take some time (particularly the PCAs). To shorten their
//    visualisation, select only those you require in the layers panel of the Map window.
// 4. The images resulting from these scripts can be transferred to your Google Drive running the
//    analysis from the 'task' tab on the top right of the screen.
// 5. All these algorithms can be run separately, however, there are variables employed by several
//    of those that only appear the first time that they are required and are not repeated in
//    subsequent algorithms. The reader is advised to check carefully which variables are employed
//    by each algorithm and to what part of the code they refer before attempting to run a single
//    algorithm separately from the rest of the code written here.


//                ------------- ooo -------------


// In the section below a set of variables defining (1) the collection of all Landsat 5 images available
// and filtering the clouds, (2) the Area of Interest (AOI) and (3) the central coordinate for visualisation
// will be defined for the use of the different algorithms

// Define your AOI as a polygon with a set of WGS84 coordinates
var geometry = ee.Geometry.Polygon(
    [[[77.59,28.8],
      [77.59,30.29],
      [77.06,30.56],
      [76.83,30.74],
      [76.55,31.09],
      [74.45,31.09],
      [73.84,30.21],
      [73.84,28.8]]]);

// Define a central point in your study area (as X,Y WGS84 decimal degrees) and a scale,
// just for visualisation purposes
Map.setCenter(75.564, 29.874, 8);

// Set all the available Landsat 5 atmospherically corrected surface reflectance images as a variable
var L5_cloud = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(geometry);

// Collection 2 QA_PIXEL flags for fill (0), high-confidence cloud (3) and cloud shadow (4).
function clearL5Pixels(image) {
  return image.select('QA_PIXEL').bitwiseAnd((1 << 0) | (1 << 3) | (1 << 4)).eq(0);
}

var cloudMaskL457 = function(image) {
  var optical = image.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'],
                             ['B1', 'B2', 'B3', 'B4', 'B5', 'B7']);
  // Remove edge pixels that don't occur in all bands
  var mask2 = optical.mask().reduce(ee.Reducer.min());
  // Convert to surface reflectance, then retain the original reflectance x 10000 units for TCT.
  return optical.multiply(0.0000275).subtract(0.2).multiply(10000).toFloat()
      .updateMask(clearL5Pixels(image)).updateMask(mask2)
      .set('system:time_start', image.get('system:time_start'));
};

var L5 = L5_cloud.map(cloudMaskL457);


//                ------------- ooo -------------


// NDVSI (Normalised Difference Vegetation Seasonality Index)

// Load a Landsat 5 collection (1984-2012) and select the bands of interest: B3 (red) and B4 (NIR)
var L5_NDVSI = L5.select(['B3', 'B4']);

// Select the L5 images corresponding to the wet months and average their values
var wet1 = L5_NDVSI.filter(ee.Filter.dayOfYear(1,120))
    .mean();
var wet2 = L5_NDVSI.filter(ee.Filter.dayOfYear(181,240))
    .mean();
var wet = (wet1.add(wet2)).divide(2);

// Develop a simple vegetation ratio for the wet months
var dvi_wet = wet.expression(
    'b("B4") / (b("B4") + b("B3"))');

// Select the L5 images corresponding to the dry months and average their values
var dry1 = L5_NDVSI.filter(ee.Filter.dayOfYear(121,180))
    .mean();
var dry2 = L5_NDVSI.filter(ee.Filter.dayOfYear(241,360))
    .mean();
var dry = (dry1.add(dry2)).divide(2);

// Develop a simple vegetation ratio for the dry months
var dvi_dry = dry.expression(
    'b("B4") / (b("B4") + b("B3"))');

// Develop the NDVSI and clip it using the AOI
var ndvsi = ((dvi_dry.subtract(dvi_wet)))
    .divide((dvi_dry.add(dvi_wet)))
    .clip(geometry);

var ndvsi_palette =
    '011301, 011d01, 012e01, 023b01, 004c00, 056201, 207401, 3e8601, 529400,' +
    '66a000, 74a901, 99b718, fcd163, f1b555, df923d, ce7e45, 9d350e';

Map.addLayer(ndvsi, {min: -0.0705, max: 0.00408, palette: ndvsi_palette}, 'NDVSI');

Export.image.toDrive({
  image: ndvsi,
  description: 'NDVSI',
  scale: 30,
  maxPixels: 1e9,
  region: geometry
});


//                ------------- ooo -------------


// SMTVI (Seasonal Multi-Temporal Vegetation Indices)

// Rebuild the Landsat 5 eight-day EVI composites from Tier 1 daytime scenes over the AOI.
var L5_raw = ee.ImageCollection('LANDSAT/LT05/C02/T1')
    .filterBounds(geometry).filter(ee.Filter.gt('SUN_ELEVATION', 0));

print('Landsat 5 scenes intersecting the AOI (before pixel masking)', ee.Dictionary({
  surface_reflectance: L5_cloud.size(), daytime_level_1_for_EVI: L5_raw.size()
}));

function sceneEVI(image) {
  var toa = ee.Algorithms.Landsat.TOA(image);
  var nir = toa.select('B4');
  var red = toa.select('B3');
  var blue = toa.select('B1');
  var denominator = nir.add(red.multiply(6)).subtract(blue.multiply(7.5)).add(1);
  var evi = nir.subtract(red).multiply(2.5).divide(denominator).rename('EVI');
  return evi.updateMask(clearL5Pixels(image)).updateMask(denominator.neq(0))
      .updateMask(evi.gte(-1).and(evi.lte(1))).toFloat()
      .set('system:time_start', image.get('system:time_start'));
}

// Start a new eight-day sequence on 1 January each year, including the overlapping final window.
var windowStarts = ee.List.sequence(1984, 2012).map(function(year) {
  var january1 = ee.Date.fromYMD(year, 1, 1);
  return ee.List.sequence(0, 360, 8).map(function(offset) {
    return january1.advance(offset, 'day').millis();
  });
}).flatten();

// A fully masked EVI band gives empty windows zero weight in the two-month means.
var emptyEVI = ee.Image.constant(0).rename('EVI').toFloat().selfMask()
    .set('system:time_start', 0);
var eviL5 = ee.ImageCollection.fromImages(windowStarts.map(function(millis) {
  var start = ee.Date(ee.Number(millis));
  var scenes = L5_raw.filterDate(start, start.advance(8, 'day')).map(sceneEVI);
  // The most recent valid observation at each pixel has priority in the original composites.
  return ee.ImageCollection([emptyEVI]).merge(scenes).sort('system:time_start').mosaic()
      .set('system:time_start', start.millis());
}));

// Filter them by two-moths periods and extract the average values 
var eviL5_JanFeb = eviL5.filter(ee.Filter.dayOfYear(1,60))
    .mean();
var eviL5_MarApr = eviL5.filter(ee.Filter.dayOfYear(61,120))
    .mean();
var eviL5_MayJun = eviL5.filter(ee.Filter.dayOfYear(121,180))
    .mean();
var eviL5_JulAug = eviL5.filter(ee.Filter.dayOfYear(181,240))
    .mean();
var eviL5_SepOct = eviL5.filter(ee.Filter.dayOfYear(241,300))
    .mean();
var eviL5_NovDec = eviL5.filter(ee.Filter.dayOfYear(301,360))
    .mean();

// Create a composite image with all the two-month average EVI composites and clip the area of analysis
// according to your AOI (the AOI is defined by the geometry variable, which can be created/modified using
// the map window below)
var composite = ee.Image([eviL5_JanFeb, eviL5_MarApr, eviL5_MayJun, eviL5_JulAug, eviL5_SepOct, eviL5_NovDec])
    .clip(geometry);

// Rename the composite bands so they can be easier to interpret
var SMTVI = composite.select(
    ['EVI', 'EVI_1', 'EVI_2', 'EVI_3', 'EVI_4', 'EVI_5'], // old names
    ['EVI_JanFeb', 'EVI_MarApr', 'EVI_MayJun', 'EVI_JulAug', 'EVI_SepOct', 'EVI_NovDec'] // new names
);

// Add the newly created layer to the map window. In this case we have created a RGB composite joining
// the bi-month EVI average values from July-August (R), March-April (G) and January-February (B) 
Map.addLayer(SMTVI, {bands: ["EVI_JulAug","EVI_MarApr","EVI_JanFeb"], min: 0.12340243216502006, max: 0.5330676407686823}, 'SMTVI');

// Export the SMTVI image (with all bi-month EVI averages) to your Google drive
Export.image.toDrive({
  image: SMTVI,
  description: 'SMTVI',
  scale: 30,
  maxPixels: 1e9,
  region: geometry
});


//                ------------- ooo -------------


// Wet Months Multi-Temporal Tasselled Cap Transformation

// Define an Array of Tasselled Cap coefficients.
var coefficients = ee.Array([
  [0.3037, 0.2793, 0.4743, 0.5585, 0.5082, 0.1863],
  [-0.2848, -0.2435, -0.5436, 0.7243, 0.0840, -0.1800],
  [0.1509, 0.1973, 0.3279, 0.3406, -0.7112, -0.4572],
  [-0.8242, 0.0849, 0.4392, -0.0580, 0.2012, -0.2768],
  [-0.3280, 0.0549, 0.1075, 0.1855, -0.4357, 0.8085],
  [0.1084, -0.9022, 0.4120, 0.0573, -0.0251, 0.0238]
]);

// Load a Landsat 5 collection (1984-2012) and select the bands of interest
var L5_bands = L5.select(['B1', 'B2', 'B3', 'B4', 'B5', 'B7']);

// Select the L5 images corresponding to the wet months, average their values and clip to AOI 
var L5_bands_Wet1 = L5_bands.filter(ee.Filter.dayOfYear(1,120))
    .mean();
var L5_bands_Wet2 = L5_bands.filter(ee.Filter.dayOfYear(181,240))
    .mean();
var L5_bands_Wet = (L5_bands_Wet1.add(L5_bands_Wet2)).divide(2)
    .clip(geometry);

// Make an Array Image, with a 1-D Array per pixel.
var arrayImage1D_Wet = L5_bands_Wet.toArray();

// Make an Array Image with a 2-D Array per pixel, 6x1.
var arrayImage2D_Wet = arrayImage1D_Wet.toArray(1);

// Do a matrix multiplication: 6x6 times 6x1.
var TCT_LS5_Wet = ee.Image(coefficients)
  .matrixMultiply(arrayImage2D_Wet)
  // Get rid of the extra dimensions.
  .arrayProject([0])
  .arrayFlatten(
    [['brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth']]);

// Display the first three bands of the result and the input imagery.
var vizTCTwet = {
  bands: ['brightness', 'greenness', 'wetness'],
  min: [2915.8, -35.067, -1376.2], max: [6010, 1612.2, 442.56]
};

// Add the resulting layer
Map.addLayer(TCT_LS5_Wet, vizTCTwet, 'TCT_LS5_Wet');

// Export the Multitemporal TCT of the wet months to your Google Drive
Export.image.toDrive({
  image: TCT_LS5_Wet,
  description: 'TCT_L5_Wet',
  scale: 30,
  maxPixels: 1e9,
  region: geometry
});


//                ------------- ooo -------------


// Dry Months Multi-Temporal Tasselled Cap Transformation

// NOTES!
// 1. The Array of Tasselled Cap coefficients has already been defined in the previous section
// 2. The bands of interest from the Landsat 5 collection have already been selected in the
//    previous section

// Select the L5 images corresponding to the dry months, average their values and clip to AOI
var L5_bands_Dry1 = L5_bands.filter(ee.Filter.dayOfYear(121,180))
    .mean();
var L5_bands_Dry2 = L5_bands.filter(ee.Filter.dayOfYear(241,360))
    .mean();
var L5_bands_Dry = (L5_bands_Dry1.add(L5_bands_Dry2)).divide(2)
    .clip(geometry);

// Make an Array Image, with a 1-D Array per pixel.
var arrayImage1D_Dry = L5_bands_Dry.toArray();

// Make an Array Image with a 2-D Array per pixel, 6x1.
var arrayImage2D_Dry = arrayImage1D_Dry.toArray(1);

// Do a matrix multiplication: 6x6 times 6x1.
var TCT_LS5_Dry = ee.Image(coefficients)
  .matrixMultiply(arrayImage2D_Dry)
  // Get rid of the extra dimensions.
  .arrayProject([0])
  .arrayFlatten(
    [['brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth']]);

// Display the first three bands of the result and the input imagery.
var vizTCTdry = {
  bands: ['brightness', 'greenness', 'wetness'],
  min: [3259.5, -253.58, -1692.1], max: [6306.6, 1596.27, -342.66]
};

// Add the resulting layer
Map.addLayer(TCT_LS5_Dry, vizTCTdry, 'TCT_LS5_Dry');

// Export the Multitemporal TCT of the dry months to your Google Drive
Export.image.toDrive({
  image: TCT_LS5_Dry,
  description: 'TCT_L5_Dry',
  scale: 30,
  maxPixels: 1e9,
  region: geometry
});


//                ------------- ooo -------------


// Wet Months Multi-Temporal Principal Component Analysis

// NOTES!
// 1. The bands of interest from the Landsat 5 collection have already been selected in the
//    previous section
// 2. The filter of seasonal images has already been done in previous sections


// Get some information about the input to be used later.
// Use a Landsat 30 m grid for both PCAs. A multi-scene mean can have a one-degree default grid.
var scale = 30;
var pcaProjection = ee.Image(L5.first()).select('B1').projection();
var bandNames = L5_bands_Wet.bandNames();

// Mean center the data to enable a faster covariance reducer
// and an SD stretch of the principal components.
var meanDict = L5_bands_Wet.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: geometry,
    crs: pcaProjection,
    scale: scale,
    tileScale: 4,
    maxPixels: 1e9
});
var means = ee.Image.constant(meanDict.values(bandNames));
var centered = L5_bands_Wet.subtract(means);

// This helper function returns a list of new band names.
var getNewBandNames = function(prefix) {
  var seq = ee.List.sequence(1, bandNames.length());
  return seq.map(function(b) {
    return ee.String(prefix).cat(ee.Number(b).format('%d'));
  });
};

// This function accepts mean centered imagery, a scale and
// a region in which to perform the analysis.  It returns the
// Principal Components (PC) in the region as a new image.
var getPrincipalComponents = function(centered, scale, geometry) {
  // Collapse the bands of the image into a 1D array per pixel.
  var arrays = centered.toArray();

  // Compute the covariance of the bands within the region.
  var covar = arrays.reduceRegion({
    reducer: ee.Reducer.centeredCovariance(),
    geometry: geometry,
    crs: pcaProjection,
    scale: scale,
    tileScale: 4,
    maxPixels: 1e9
  });

  // Get the 'array' covariance result and cast to an array.
  // This represents the band-to-band covariance within the region.
  var covarArray = ee.Array(covar.get('array'));

  // Perform an eigen analysis and slice apart the values and vectors.
  var eigens = covarArray.eigen();

  // This is a P-length vector of Eigenvalues.
  var eigenValues = eigens.slice(1, 0, 1);
  // This is a PxP matrix with eigenvectors in rows.
  var eigenVectors = eigens.slice(1, 1);

  // Convert the array image to 2D arrays for matrix computations.
  var arrayImage = arrays.toArray(1);

  // Left multiply the image array by the matrix of eigenvectors.
  var principalComponents = ee.Image(eigenVectors).matrixMultiply(arrayImage);

  // Turn the square roots of the Eigenvalues into a P-band image.
  var sdImage = ee.Image(eigenValues.sqrt())
    .arrayProject([0]).arrayFlatten([getNewBandNames('sd')]);

  // Turn the PCs into a P-band image, normalized by SD.
  return principalComponents
    // Throw out an an unneeded dimension, [[]] -> [].
    .arrayProject([0])
    // Make the one band array image a multi-band image, [] -> image.
    .arrayFlatten([getNewBandNames('pc')])
    // Normalize the PCs by their SDs.
    .divide(sdImage);
};

// Get the PCs at the specified scale and in the specified region
var pcaL5wet = getPrincipalComponents(centered, scale, geometry);

//Define the visualisation parameters
var vizPCAwet = {
  bands: ['pc2', 'pc4', 'pc1'],
  min: [-3.2339], max: [3.2845]
};

// Display a composite with PCs 2,4,1 for visualisation
// Please change the visualization paramenters to combine different bands
// and/or reproduce the band combinations of the paper's figures
// NOTE! The exported image incorporates all Principal Components
Map.addLayer(pcaL5wet, vizPCAwet, 'PCA_composite_(2,4,1)_Wet');

Export.image.toDrive({
  image: pcaL5wet,
  description: 'PCA_L5_Wet',
  scale: 30,
  maxPixels: 1e9,
  region: geometry
});


//                ------------- ooo -------------


// Dry Months Multi-Temporal Principal Component Analysis

// NOTES!
// 1. The bands of interest from the Landsat 5 collection have already been selected in the
//    previous section
// 2. The filter of seasonal images has already been done in previous sections
// 3. region, scale, and bandNames variables have already been created in the previous section


// Mean center the data to enable a faster covariance reducer
// and an SD stretch of the principal components.
var meanDict_Dry = L5_bands_Dry.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: geometry,
    crs: pcaProjection,
    scale: scale,
    tileScale: 4,
    maxPixels: 1e9
});
var means_Dry = ee.Image.constant(meanDict_Dry.values(bandNames));
var centered_Dry = L5_bands_Dry.subtract(means_Dry);

// This helper function returns a list of new band names.
var getNewBandNames_Dry = function(prefix) {
  var seq_Dry = ee.List.sequence(1, bandNames.length());
  return seq_Dry.map(function(b) {
    return ee.String(prefix).cat(ee.Number(b).format('%d'));
  });
};

// This function accepts mean centered imagery, a scale and
// a region in which to perform the analysis.  It returns the
// Principal Components (PC) in the region as a new image.
var getPrincipalComponents_Dry = function(centered_Dry, scale, geometry) {
  // Collapse the bands of the image into a 1D array per pixel.
  var arrays_Dry = centered_Dry.toArray();

  // Compute the covariance of the bands within the region.
  var covar_Dry = arrays_Dry.reduceRegion({
    reducer: ee.Reducer.centeredCovariance(),
    geometry: geometry,
    crs: pcaProjection,
    scale: scale,
    tileScale: 4,
    maxPixels: 1e9
  });

  // Get the 'array' covariance result and cast to an array.
  // This represents the band-to-band covariance within the region.
  var covarArray_Dry = ee.Array(covar_Dry.get('array'));

  // Perform an eigen analysis and slice apart the values and vectors.
  var eigens_Dry = covarArray_Dry.eigen();

  // This is a P-length vector of Eigenvalues.
  var eigenValues_Dry = eigens_Dry.slice(1, 0, 1);
  // This is a PxP matrix with eigenvectors in rows.
  var eigenVectors_Dry = eigens_Dry.slice(1, 1);

  // Convert the array image to 2D arrays for matrix computations.
  var arrayImage_Dry = arrays_Dry.toArray(1);

  // Left multiply the image array by the matrix of eigenvectors.
  var principalComponents_Dry = ee.Image(eigenVectors_Dry).matrixMultiply(arrayImage_Dry);

  // Turn the square roots of the Eigenvalues into a P-band image.
  var sdImage_Dry = ee.Image(eigenValues_Dry.sqrt())
    .arrayProject([0]).arrayFlatten([getNewBandNames_Dry('sd')]);

  // Turn the PCs into a P-band image, normalized by SD.
  return principalComponents_Dry
    // Throw out an an unneeded dimension, [[]] -> [].
    .arrayProject([0])
    // Make the one band array image a multi-band image, [] -> image.
    .arrayFlatten([getNewBandNames_Dry('pc')])
    // Normalize the PCs by their SDs.
    .divide(sdImage_Dry);
};

// Get the PCs at the specified scale and in the specified region
var pcaL5dry = getPrincipalComponents_Dry(centered_Dry, scale, geometry);

//Define the visualisation parameters
var vizPCAdry = {
  bands: ['pc2', 'pc4', 'pc1'],
  min: [-2.2689], max: [2.3412]
};

// Display a composite with PCs 2,4,1 for visualisation
// Please change the visualization paramenters to combine different bands
// and/or reproduce the band combinations of the paper's figures
// NOTE! The exported image incorporates all Principal Components
Map.addLayer(pcaL5dry, vizPCAdry, 'PCA_composite_(2,4,1)_Dry');

Export.image.toDrive({
  image: pcaL5dry,
  description: 'PCA_L5_Dry',
  scale: 30,
  maxPixels: 1e9,
  region: geometry
});
