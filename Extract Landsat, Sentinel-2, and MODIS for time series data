//set StudyArea
var StudyPoint = ee.Geometry.Point(lat,lon);
//-----------------------------------------------------------------MaskCloud--------------------------------------------------------

//Mask Landsat 5,7,8 Cloud

// Mask cloud and shadow using SR CFMASK and TOA FMask in Landsat5 SR
function landsat5SRCloudShadowCFmask_FMask(image) // For USGS version SR
{
	var timeStart = image.get('system:time_start'); // Data availability (time): Jan 1, 1984 - May 5, 2012
  var toaImageList = ee.ImageCollection('LANDSAT/LT5_L1T_TOA_FMASK') /*USGS Landsat 5 TOA Reflectance (Orthorectified) with Fmask*/
                       .filterMetadata('system:time_start','equals',timeStart)
                       .toList(5);
  // If TOA image existed
  var fmaskTOA = ee.Image(ee.Algorithms.If(toaImageList.size().gt(0),
                          ee.Image(toaImageList.get(0)).select('fmask'),
                          ee.Image(1) /* if TOA not exist , creat an all-1 image */
                          ));
  var cfmaskSR = image.select('cfmask');
  var cfmaskSR_conf = image.select('cfmask_conf');
  var mask = cfmaskSR.eq(0)/* 0=clear, set all 0 as 1 while non-0 pixels as 0  */
            .or(cfmaskSR.eq(1)) /* 1=water, then add all 1 pixels to the above image */
            .and(cfmaskSR_conf.lte(1)) /*0=none, 1=cloud confidence ≤ 12.5% , select1 pixels in above image  with cloud mask confidence*/
            .and(fmaskTOA.eq(0).or(fmaskTOA.eq(1))); 
  return image.updateMask(mask);
}


// Mask cloud and shadow using SR CFMASK and TOA FMask in Landsat7 SR
function landsat7SRCloudShadowCFmask_FMask(image) // For USGS version SR
{
	var timeStart = image.get('system:time_start');
  var toaImageList = ee.ImageCollection('LANDSAT/LE7_L1T_TOA_FMASK') /*USGS Landsat 7 TOA Reflectance (Orthorectified) with Fmask*/
                       .filterMetadata('system:time_start','equals',timeStart)
                       .toList(5);
  // If TOA image existed
  var fmaskTOA = ee.Image(ee.Algorithms.If(toaImageList.size().gt(0),
                          ee.Image(toaImageList.get(0)).select('fmask'),
                          ee.Image(1) /* if TOA not exist , creat an all-1 image */
                          ));
  var cfmaskSR = image.select('cfmask');
  var cfmaskSR_conf = image.select('cfmask_conf');
  var mask = cfmaskSR.eq(0)/* 0=clear, set all 0 as 1 while non-0 pixels as 0  */
            .or(cfmaskSR.eq(1)) /* 1=water, then add all 1 pixels to the above image */
            .and(cfmaskSR_conf.lte(1)) /*0=none, 1=cloud confidence ≤ 12.5% , select1 pixels in above image  with cloud mask confidence*/
            .and(fmaskTOA.eq(0).or(fmaskTOA.eq(1))); 
  return image.updateMask(mask);
}


// Mask cloud and shadow using SR CFMASK and TOA FMask in Landsat8 SR
function landsat8SRCloudShadowCFmask_FMask(image) // For USGS version SR
{ 
	var timeStart = image.get('system:time_start');
  var toaImageList = ee.ImageCollection('LANDSAT/LC8_L1T_TOA_FMASK')
                   .filterMetadata('system:time_start','equals',timeStart)
                   .toList(5);
  // If TOA image existed
  var fmaskTOA = ee.Image(ee.Algorithms.If(toaImageList.size().gt(0),
                          ee.Image(toaImageList.get(0)).select('fmask'),
                          ee.Image(1)));
  var cfmaskSR = image.select('cfmask');
  var cfmaskSR_conf = image.select('cfmask_conf');
  var mask = cfmaskSR.eq(0).or(cfmaskSR.eq(1)).and(cfmaskSR_conf.lte(1))
                     .and(fmaskTOA.eq(0).or(fmaskTOA.eq(1)));
  return image.updateMask(mask);
}

////////Mask Modis Cloud

// Filtering bad observations
// See Sur_refl_state_500m in https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod09a1
// Note MODIS products usually come with Quality Assurance layer with values coded/packed in bits.
// Create a function to extract bits from a position i to j
var getQABits = function(image, start, end, newName) {
  // Compute the bits we need to extract
  var pattern = 0;
  for (var i = start; i <= end; i++) {
    pattern += Math.pow(2,i);
  }
  // Return a single band image of the extracted QA bits, giving the band a new name.
  return image.select([0], [newName])
              .bitwiseAnd(pattern)
              .rightShift(start);
};

// Filter based on Sur_refl_state_500m band
var filterBadObs = function(image){
  // [bit 0-1] "MOD35 cloud" =  0 ("clear")
  // Get the QA band
    var cloudQA = getQABits(image.select('StateQA'),0,1,'cloud');

    // Mask the input image with the above QA bands
    var maskedImage = ee.Image(0).where(cloudQA.eq(0),1);
    
    // Return only good pixels
    return image.mask(maskedImage);  
};
 
// Function to mask clouds using the Sentinel-2 QA band.
function maskS2clouds(image) {
  var qa = image.select('QA60');
  
  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = Math.pow(2, 10);
  var cirrusBitMask = Math.pow(2, 11);
  
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(
             qa.bitwiseAnd(cirrusBitMask).eq(0));

  // Return the masked and scaled data.
  return image.updateMask(mask);
}
//----------------------------------------------------------GetDataCollection-----------------------------------------------

//Load collections of surface reflectance for Landsat 5,7 and 8 in the area of interest and cloud masked  
var L5coll = ee.ImageCollection('LANDSAT/LT5_SR')
.filterDate('2000-01-01', '2012-05-05')//Time series(2000.1.1-2012.5.5)
.filterBounds(StudyPoint)
.map(landsat5SRCloudShadowCFmask_FMask)
.sort('system:time_start', true);

var L7coll = ee.ImageCollection('LANDSAT/LE7_SR')
.filterDate('2000-01-01', '2017-08-01')//Time series(2000.1.1-2017.8.1)
.filterBounds(StudyPoint)
.map(landsat7SRCloudShadowCFmask_FMask)
.sort('system:time_start', true);

var L8coll = ee.ImageCollection('LANDSAT/LC8_SR')
.filterDate('2013-04-11', '2017-08-01')//Time series(2013.4.11-2017.8.1)
.filterBounds(StudyPoint)
.map(landsat8SRCloudShadowCFmask_FMask)
.sort('system:time_start', true);
//load modis
var ModisCollection = ee.ImageCollection('MODIS/MOD09A1')
    .filterDate(ee.Date('2000-02-18'), ee.Date('2017-08-01'))//Time series(2000.2.18-2017.8.1)
    // .select(MODIS_BANDS, MODIS_NAMES)
    .filterBounds(StudyPoint)
    .sort('system:time_start', true)
    .map(filterBadObs);
//load sentinel2
var composite = ee.ImageCollection('COPERNICUS/S2')
                  .filterDate('2014-01-01', '2017-08-01')//Time series(2014.1.1-2017.8.1)
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
                  .filterBounds(StudyPoint)
                  .select(['B2','B3','B4','B8','B11','QA60'],['blue','green','red','nir','swir','QA60'])
                  .map(maskS2clouds);
//-------------------------------------------------------------ComputeVI&AddBands--------------------------------------------------------
// calculate LE5 NDVI 
// multiply 0.0001 is needed because the SR data is scaled by 10000
var L5_ndvi = L5coll.map(
  function(image) {
    var ndvi = image.expression(
      '(nir - red) / (nir + red)',
        {
          red: image.select('B3').multiply(0.0001),    // 620-670nm, RED
          nir: image.select('B4').multiply(0.0001),    // 841-876nm, NIR
          blue: image.select('B1').multiply(0.0001)    // 459-479nm, BLUE
        });
        
    var lswi=image.normalizedDifference(['B4','B5']);//lswi=(nir-swir)/(nir+swir)
    var evi=image.expression(
                  '2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 10000)',
             {
              red: image.select('B3'),    // 620-670nm, RED
              nir: image.select('B4'),    // 841-876nm, NIR
              blue: image.select('B1')    // 459-479nm, BLUE
             }
            );        
    // // Rename that band to something appropriate
    // return ndvi.select([0], ['ndvi']).set('system:time_start', image.get('system:time_start'));
    //creating new image collections by andding bands of ndvi,evi and lswi
    return image.select(["B1","B2","B3","B4","B5","B7",'cfmask','cfmask_conf'])
              .addBands([ndvi,evi,lswi])
              .select([0,1,2,3,4,5,6,7,8,9,10],['blue','green','red','nir','swir1','swir2','cfmask','cfmask_conf','ndvi','evi','lswi'])
              .set('system:time_start', image.get('system:time_start'));
  }
);

// calculate LE7 EVI 
// multiply 0.0001 is needed because the SR data is scaled by 10000
var L7_ndvi = L7coll.map(
  function(image) {
    var ndvi = image.expression(
      '(nir - red) / (nir + red)',
        {
          red: image.select('B3').multiply(0.0001),    // 620-670nm, RED
          nir: image.select('B4').multiply(0.0001),    // 841-876nm, NIR
          blue: image.select('B1').multiply(0.0001)    // 459-479nm, BLUE
        });
    var lswi=image.normalizedDifference(['B4','B5']);//lswi=(nir-swir)/(nir+swir)
    var evi=image.expression(
                  '2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 10000)',
             {
              red: image.select('B3'),    // 620-670nm, RED
              nir: image.select('B4'),    // 841-876nm, NIR
              blue: image.select('B1')    // 459-479nm, BLUE
             }
            ); 

    // // Rename that band to something appropriate
    // return ndvi.select([0], ['ndvi']).set('system:time_start', image.get('system:time_start'));
    //creating new image collections by andding bands of ndvi,evi and lswi
    return image.select(["B1","B2","B3","B4","B5","B7",'cfmask','cfmask_conf'])
              .addBands([ndvi,evi,lswi])
              .select([0,1,2,3,4,5,6,7,8,9,10],['blue','green','red','nir','swir1','swir2','cfmask','cfmask_conf','ndvi','evi','lswi'])    
              .set('system:time_start', image.get('system:time_start'));
  }
);

// calculate LC8 EVI - note the bands are different for LC8
// multiply 0.0001 is needed because the SR data is scaled by 10000
var L8_ndvi = L8coll.map(
  function(image) {
    var ndvi = image.expression(
      '(nir - red) / (nir + red)',
        {
          red: image.select('B4').multiply(0.0001),    // 620-670nm, RED
          nir: image.select('B5').multiply(0.0001),    // 841-876nm, NIR
          blue: image.select('B2').multiply(0.0001)    // 459-479nm, BLUE
        });
    var lswi=image.normalizedDifference(['B5','B6']);//lswi=(nir-swir)/(nir+swir)
    var evi=image.expression(
                  '2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 10000)',
             {
              red: image.select('B4'),    // 620-670nm, RED
              nir: image.select('B5'),    // 841-876nm, NIR
              blue: image.select('B2')    // 459-479nm, BLUE
             }
            ); 

    // // Rename that band to something appropriate
    // return ndvi.select([0], ['ndvi']).set('system:time_start', image.get('system:time_start'));
    //creating new image collections by andding bands of ndvi,evi and lswi
    return image.select(["B2","B3","B4","B5","B6","B7",'cfmask','cfmask_conf'])
              .addBands([ndvi,evi,lswi])
              .select([0,1,2,3,4,5,6,7,8,9,10],['blue','green','red','nir','swir1','swir2','cfmask','cfmask_conf','ndvi','evi','lswi'])    
              .set('system:time_start', image.get('system:time_start'));    
  }
);

// calculate MODIS MOD09A1 NDVI,LSWI and EVI
var ModisVi = ModisCollection.map(
  function(image){
    var ndvi = image.expression(
      '(nir - red) / (nir + red)',
      {
        red:image.select('sur_refl_b01').multiply(0.0001),//
        nir:image.select('sur_refl_b02').multiply(0.0001),//
      });
    var lswi=image.normalizedDifference(['sur_refl_b02','sur_refl_b06']);//lswi=(nir-swir)/(nir+swir)
    var evi=image.expression(
      '2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 10000)',
      {
        red: image.select('sur_refl_b01'),    // 620-670nm, RED
        nir: image.select('sur_refl_b02'),    // 841-876nm, NIR
        blue: image.select('sur_refl_b03')    // 459-479nm, BLUE
      }
      ); 
     return image.select(["sur_refl_b01","sur_refl_b02","sur_refl_b03","sur_refl_b04","sur_refl_b05","sur_refl_b06","sur_refl_b07",'QA','StateQA'])
              .addBands([ndvi,evi,lswi])
              .select([0,1,2,3,4,5,6,7,8,9,10,11],['red','nir1','blue','green','nir2','swir1','swir2','QA','StateQA','ndvi','evi','lswi'])
              .set('system:time_start', image.get('system:time_start'));
  }
  );
// calculate Sentinel2 NDVI,LSWI and EVI  
var S2Vi = composite.map(
  function(image){
    var ndvi = image.expression(
      '(nir - red) / (nir + red)',
      {
        red:image.select('red'),//
        nir:image.select('nir')//
      }).rename('ndvi');
    var lswi=image.normalizedDifference(['nir','swir']).rename('lswi');//lswi=(nir-swir)/(nir+swir)
    var evi=image.expression(
      '2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 10000)',
      {
        red: image.select('red'),    // 620-670nm, RED
        nir: image.select('nir'),    // 841-876nm, NIR
        blue: image.select('blue')    // 459-479nm, BLUE
      }
      ).rename('evi'); 
     return image.addBands([ndvi,evi,lswi]).set('system:time_start', image.get('system:time_start'))
}
  );

//---------------------------------------------------------Sample------------------------------------------------------------

// Add a time band for each EVI collection
var L5time = L5_ndvi.map(function(image) {
  var time = ee.Image.constant(image.date().millis()).rename("time");
// Get a collection containing values for all of the grid points where there are values.
  return image.addBands(time).sampleRegions({
    collection: StudyPoint,
    properties: ["name"],
  });
}).flatten();

var L7time = L7_ndvi.map(function(image) {
  var time = ee.Image.constant(image.date().millis()).rename("time");
// Get a collection containing values for all of the grid points where there are values.
  return image.addBands(time).sampleRegions({
    collection: StudyPoint,
    properties: ["name"],
  });
}).flatten();
var L8time = L8_ndvi.map(function(image) {
  var time = ee.Image.constant(image.date().millis()).rename("time");
// Get a collection containing values for all of the grid points where there are values.
  return image.addBands(time).sampleRegions({
    collection: StudyPoint,
    properties: ["name"],
  });
}).flatten();

//Sample MODIS
var ModTime = ModisVi.map(function(image) {
// Get a collection containing values for all of the grid points where there are values.
  return image.sampleRegions({
    collection: StudyPoint,
    properties: ["name"],
  });
}).flatten();

//sample sentinel2
var S2Time = S2Vi.map(function(image) {
// Get a collection containing values for all of the grid points where there are values.
  return image.sampleRegions({
    collection: StudyPoint,
    properties: ["name"],
  });
}).flatten();

//Merge over a single collection
// var Final_collection = ee.FeatureCollection(L5time.merge(L7time)).merge(L8time)
// var Final2 = ee.FeatureCollection(ModTime.merge(S2Time));
// var final = ee.FeatureCollection(Final_collection.merge(Final2));
//------------------------------------------------------Export-------------------------------------------------------
// //Export all the collections
// Export.table.toDrive({
//   collection: final,
//   description: 'AllVIs_Peru', 
//   fileFormat: 'CSV'
// });

//Export the results for L5
Export.table.toDrive({
  collection: L5time,
  description: 'Landsat5VIs_Peru', 
  fileFormat: 'CSV'
});

//Export the results for L7
Export.table.toDrive({
  collection: L7time,
  description: 'Landsat7VIs_Peru', 
  fileFormat: 'CSV'
});

//Export the results for L8
Export.table.toDrive({
  collection: L8time,
  description: 'Landsat8VIs_Peru', 
  fileFormat: 'CSV'
});

//Export MODIS
Export.table.toDrive({
  collection: ModTime,
  description: 'AllModisVIs_Peru', 
  fileFormat: 'CSV'
});

//Export sentinel
Export.table.toDrive({
  collection: S2Time,
  description: 'AllSentinel2VIs_Peru', 
  fileFormat: 'CSV'
});
