//Getting the multiply bands about VI and greenness,inundation and cannopy.
// Load a Landsat 8 collection for a single path-row.
//variable
var start_date = '2014-01-01';
var end_date = '2016-07-23';
var point_x = 109.07;
var point_y = 21.455;
var studyArea = Hainan1;
var VIthreshold = 0.3;
var DEM = ee.Image('USGS/SRTMGL1_003');
var slope = ee.Terrain.slope(DEM);

//function: mask cloud
function landsat8SRCloudShadowCFmask_FMask(image)
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

// Mask cloud and shadow using SR CFMASK and TOA FMask
function landsat7SRCloudShadowCFmask_FMask(image) // For USGS version SR
{
	var timeStart = image.get('system:time_start');
  var toaImageList = ee.ImageCollection('LANDSAT/LE7_L1T_TOA_FMASK')
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


function ND_VI(image,b1,b2,bName)
{
  var VI = image.normalizedDifference([b1,b2]).rename(bName);
  return VI.updateMask(VI.gt(-1).and(VI.lt(1)));
}

function funEVI(image,B1,B2,B3)
{
 
	var VI = image.expression('2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)',
    {
      blue: image.select(B1).multiply(0.0001),   
      red:  image.select(B2).multiply(0.0001),   
      nir:  image.select(B3).multiply(0.0001)
    }).rename('EVI');
  return VI.updateMask(VI.gt(-1).and(VI.lt(1)));

}


//function: greenness
function greenness01Msk(image)
{
  return image.select('EVI').gt(0.2).and(image.select('LSWI').gt(0.0));
}

//function innudation
function innuadation01Msk(image)
{
  return image.select('LSWI').gte(image.select('EVI'))
    .or(image.select('LSWI').gte(image.select('NDVI')));
}

//function mangrove
function mangrove01Msk(image)
{
  return image.select('NDVI').gt(0.3)
    .and(image.select('LSWI').gt(0.3));
}
//function: mndwi>0
function mndwi01Msk(image)
{
  return image.select('mNDWI').gt(0.0);
}



//function lswi < 0
function lswi01Msk(image)
{
  return image.select('LSWI').lt(0.0);
}

// Add NDVI, LSWI,EVI
// function addLandsat8VIs(img)
// {
//   var NDVI = ND_VI(img,'B5','B4','NDVI');
//   var EVI = funEVI(img,'B2','B4','B5');
//   var LSWI = ND_VI(img,'B5','B6','LSWI');
//   var mNDWI = ND_VI(img,'B3','B6','mNDWI');
//   return img.addBands(NDVI).addBands(EVI).addBands(LSWI).addBands(mNDWI);
// }

function addLandsat7VIs(img)
{
  var NDVI = ND_VI(img,'B4','B3','NDVI');
  var EVI = funEVI(img,'B1','B3','B4');
  var LSWI = ND_VI(img,'B4','B5','LSWI');
  var mNDWI = ND_VI(img,'B2','B5','mNDWI');
  return img.addBands(NDVI).addBands(EVI).addBands(LSWI).addBands(mNDWI);
}

//get Landsat8 collction
var collection1 = ee.ImageCollection('LANDSAT/LC8_SR')
    .filterBounds(studyArea)
    .filterDate(start_date, end_date)
    .map(landsat8SRCloudShadowCFmask_FMask)
    .select( // Rename LC8 bands same as L5/7
      ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'cfmask','cfmask_conf']
      ,['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'cfmask','cfmask_conf'])
      .map(addLandsat7VIs);
print('Landsat8')
print(collection1.size());

//get Landsat7 collction
var collection2 = ee.ImageCollection('LANDSAT/LE7_SR')
    .filterBounds(studyArea)
    .filterDate(start_date, end_date)
    .map(landsat7SRCloudShadowCFmask_FMask)
    .select( 
      ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'cfmask','cfmask_conf'])
     .map(addLandsat7VIs);
print('Landsat7')
print(collection2.size());


//Combine two collection
var  colVIs = ee.ImageCollection(collection2.merge(collection1));
// var colVIs = collection.map(addLandsatVIs); 

// Map.addLayer(ee.Image(colVIs),{},'firstImg');
var CntValid = colVIs.select('NDVI').count().rename('CntValid');
//compute greenness freq
var CntEvgreen = colVIs.map(greenness01Msk).sum().rename('CntEvgreen');
Map.addLayer(ee.Image(CntEvgreen),{min:0,max:100},'CntEvgreen',false); 

//compute innudation freq
var CntWetland = colVIs.map(innuadation01Msk).sum().rename('CntWetland');
Map.addLayer(ee.Image(CntWetland),{min:0,max:100},'CntWetland',false);      
//compute cannopy of mangrove freq
var CntMangrove = colVIs.map(mangrove01Msk).sum().rename('CntMangrove');
Map.addLayer(ee.Image(CntMangrove),{min:0,max:100},'CntMangrove',true); 
// compute mndwi freq
var CntMndwi = colVIs.map(mndwi01Msk).sum().rename('CntMndwi');

//compute built-up or barren land Frequency :imgFreqBarSoRD
var CntBarSoRD = colVIs.map(lswi01Msk).sum().rename('CntBarSoRD');
//Map.addLayer(ee.Image(colVIs.select('NDVI').count()),{min:0,max:100},'Landsat_Count'); 
//annual mean ndvi
var imgNDVIavg = colVIs.select('NDVI').reduce(ee.Reducer.mean()).rename('imgNDVIavg');

//annual mean lswi
var imgLSWIavg = colVIs.select('LSWI').reduce(ee.Reducer.mean()).rename('imgLSWIavg');

var mvcImage = CntValid.addBands(CntEvgreen).addBands(CntWetland)
                       .addBands(CntMangrove).addBands(CntMndwi)
                       .addBands(CntBarSoRD).addBands(imgNDVIavg)
                       .addBands(imgLSWIavg);
var ExportImg = mvcImage.clip(studyArea);

Map.addLayer(ExportImg,{color:'00ff00'},'ExportImg',false);
// 
Export.image.toAsset({
  image: ExportImg,
  description: 'mvc_Hainan2_2015',
  scale: 30,
  maxPixels:1e12
});
