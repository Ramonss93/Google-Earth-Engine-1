//Map.setCenter(118.8, 32)
var start_date = '2015-01-01';
var end_date = '2017-12-31';
// var StudyArea = ee.FeatureCollection(YZbound);
var StudyArea = ee.FeatureCollection(Nanjing);
//Map.addLayer(Test_Yz)

//-------------------------------------------------------------------CloudMask------------------------------------------------------------------------
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

//-------------------------------------------------------------------------------------
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

function WaterMask(image)
{

  return image.select('MNDWI').gt(image.select('evi'))
    .or(image.select('MNDWI').gt(image.select('ndvi')))
    .and(image.select('evi').lt(0.1));

}
//-------------------------------------------------------------------Function&addBands-------------------------------------------------------------------

// calculate LE5 NDVI 
// multiply 0.0001 is needed because the SR data is scaled by 10000
  function L5_VI(image) {
    var ndvi = image.expression(
      '(nir - red) / (nir + red)',
        {
          red: image.select('B3').multiply(0.0001),    // 620-670nm, RED
          nir: image.select('B4').multiply(0.0001),    // 841-876nm, NIR
          blue: image.select('B1').multiply(0.0001)    // 459-479nm, BLUE
        }).rename('ndvi').set('system:time_start',image.get('system:time_start'));
        
    var lswi=image.normalizedDifference(['B4','B5'])
          .rename('lswi').set('system:time_start',image.get('system:time_start'));//lswi=(nir-swir)/(nir+swir)
    var evi=image.expression(
                  '2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 10000)',
             {
              red: image.select('B3').multiply(0.0001),    // 620-670nm, RED
              nir: image.select('B4').multiply(0.0001),    // 841-876nm, NIR
              blue: image.select('B1').multiply(0.0001)    // 459-479nm, BLUE
             }
            ).rename('evi').set('system:time_start',image.get('system:time_start'));
    var mndwi = image.expression('(green - mir) / (green + mir)',
           {
            green: image.select('B2').multiply(0.0001),
            mir: image.select('B5').multiply(0.0001)
            }).rename('MNDWI').set('system:time_start',image.get('system:time_start'));
    // // Rename that band to something appropriate
    // return ndvi.select([0], ['ndvi']).set('system:time_start', image.get('system:time_start'));
    //creating new image collections by andding bands of ndvi,evi and lswi
    return image.select(["B1","B2","B3","B4","B5","B7",'cfmask','cfmask_conf'])
              .addBands([ndvi,evi,lswi,mndwi])
              .set('system:time_start', image.get('system:time_start'));
  }

// calculate LE7 EVI 
// multiply 0.0001 is needed because the SR data is scaled by 10000
  function L7_VI(image) {
    var ndvi = image.expression(
      '(nir - red) / (nir + red)',
        {
          red: image.select('B3').multiply(0.0001),    // 620-670nm, RED
          nir: image.select('B4').multiply(0.0001),    // 841-876nm, NIR
          blue: image.select('B1').multiply(0.0001)    // 459-479nm, BLUE
        }).rename('ndvi').set('system:time_start',image.get('system:time_start'));
        
    var lswi=image.normalizedDifference(['B4','B5'])
          .rename('lswi').set('system:time_start',image.get('system:time_start'));//lswi=(nir-swir)/(nir+swir)
    var evi=image.expression(
                  '2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 10000)',
             {
              red: image.select('B3').multiply(0.0001),    // 620-670nm, RED
              nir: image.select('B4').multiply(0.0001),    // 841-876nm, NIR
              blue: image.select('B1').multiply(0.0001)    // 459-479nm, BLUE
             }
            ).rename('evi').set('system:time_start',image.get('system:time_start'));
    var mndwi = image.expression('(green - mir) / (green + mir)',
           {
            green: image.select('B2').multiply(0.0001),
            mir: image.select('B5').multiply(0.0001)
            }).rename('MNDWI').set('system:time_start',image.get('system:time_start'));
    // // Rename that band to something appropriate
    // return ndvi.select([0], ['ndvi']).set('system:time_start', image.get('system:time_start'));
    //creating new image collections by andding bands of ndvi,evi and lswi
    return image.select(["B1","B2","B3","B4","B5","B7",'cfmask','cfmask_conf'])
              .addBands([ndvi,evi,lswi,mndwi])
              .set('system:time_start', image.get('system:time_start'));
  }

// calculate LC8 EVI - note the bands are different for LC8
// multiply 0.0001 is needed because the SR data is scaled by 10000
  function L8_VI(image) {
    var ndvi = image.expression(
      '(nir - red) / (nir + red)',
        {
          red: image.select('B4').multiply(0.0001),    // 620-670nm, RED
          nir: image.select('B5').multiply(0.0001),    // 841-876nm, NIR
          blue: image.select('B2').multiply(0.0001)    // 459-479nm, BLUE
        }).rename('ndvi').set('system:time_start',image.get('system:time_start'));
    var lswi=image.normalizedDifference(['B5','B6']).rename('evi')
        .rename('lswi').set('system:time_start',image.get('system:time_start'));//lswi=(nir-swir)/(nir+swir)
    var evi=image.expression(
                  '2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 10000)',
             {
              red: image.select('B4').multiply(0.0001),    // 620-670nm, RED
              nir: image.select('B5').multiply(0.0001),    // 841-876nm, NIR
              blue: image.select('B2').multiply(0.0001)    // 459-479nm, BLUE
             }
            ).rename('evi').set('system:time_start',image.get('system:time_start')); 
    var mndwi = image.expression('(green - mir) / (green + mir)',
           {
            green: image.select('B3').multiply(0.0001),
            mir: image.select('B6').multiply(0.0001)
            }).rename('MNDWI').set('system:time_start',image.get('system:time_start'));
    // // Rename that band to something appropriate
    // return ndvi.select([0], ['ndvi']).set('system:time_start', image.get('system:time_start'));
    //creating new image collections by andding bands of ndvi,evi and lswi
    return image.select(["B2","B3","B4","B5","B6","B7",'cfmask','cfmask_conf'])
              .addBands([ndvi,evi,lswi,mndwi])
              .set('system:time_start', image.get('system:time_start'));    
  }

//-------------------------------------------------------------------Load collections-----------------------------------------------------------------


//Load collections of surface reflectance for Landsat 5,7 and 8 in the area of interest and cloud masked  
//L5,from 2000-01-01 to 2012-05-05
var L5coll = ee.ImageCollection('LANDSAT/LT5_SR')
.filterDate(start_date, end_date)
.filterBounds(StudyArea)
.sort('system:time_start', true)
.map(landsat5SRCloudShadowCFmask_FMask)
.map(L5_VI);
print('Landsat5')
print(L5coll.size());

//L7,FROM 2000-01-01 to 2016-12-31
var L7coll = ee.ImageCollection('LANDSAT/LE7_SR')
.filterDate(start_date, end_date)
.filterBounds(StudyArea)
.sort('system:time_start', true)
.map(landsat7SRCloudShadowCFmask_FMask)
.map(L7_VI);

print('Landsat7')
print(L7coll.size());

//L8, FROM 2013-04-11
var L8coll = ee.ImageCollection('LANDSAT/LC8_SR')
.filterDate(start_date, end_date)
.filterBounds(StudyArea)
.sort('system:time_start', true)
.map(landsat8SRCloudShadowCFmask_FMask)
.map(L8_VI);

print('Landsat8')
print(L8coll.size());

//combine two collection
var combined1 = ee.ImageCollection(L5coll.merge(L7coll));
var combined2 = ee.ImageCollection(combined1.merge(L8coll));
print('combined Landsat')
print(combined2.size());
//-------------------------------------------------------------------Frequency-----------------------------------------------------------------

//compute innudation freq
var imgFreqWater = combined2.map(WaterMask).sum().multiply(100)
     .divide(combined2.select('MNDWI').count()).rename('imgFreqWater').clip(StudyArea);

//-------------------------------------------------------------------Legend-----------------------------------------------------------------
// set palette
var vis = {min: 0,
          max: 100,
          palette: ['FFFFFF', '87CEFA', '1E90FF', '4169E1', '0000FF',
  '0000CD']};
// visualize map
Map.centerObject(StudyArea,8);
Map.addLayer(imgFreqWater,vis,'imgFreqWater');
//-------------------------------------------------------------------Export-----------------------------------------------------------------
// var valid = combined2.select('MNDWI').count();
// Create a task that you can launch from the Tasks tab.
Export.image.toDrive({
  image: imgFreqWater,
  description: 'Freq_3YGap_2016',
  scale: 30,
  maxPixels: 1e12 
});
