//-----------------------------------------------------------step1------------------------------------------------------
//Getting the multiply bands about VI and greenness,inundation and cannopy using all landsat8 images 
//during 2014-01-01 to 2016-07-23.

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


//-----------------------------------------------------------step2------------------------------------------------------
//Loading Sentinel1 images and computing permanent water frequency.
//variable
var start_date = '2014-01-01';
var end_date = '2016-07-23';
var studyArea = FJ3;
var VIthreshold = 0.3;
var DEM = ee.Image('USGS/SRTMGL1_003');
var slope = ee.Terrain.slope(DEM);


//function VH<-19
function s1_water(image)
{
  return image.select('VH').lt(-19);
}

// S1_GRD_VH
var sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  // Filter to get images with VV and VH dual polarization.
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .select('VH')
  // // Filter to get images collected in interferometric wide swath mode.
  // .filter(ee.Filter.eq('instrumentMode', 'IW'))
  // Select only high quality images
  .filter(ee.Filter.eq('resolution', 'H'))
  // Filter to extract the correct area
  .filterBounds(studyArea)
  // Filter to extract the correct period
  .filterDate(start_date, end_date)
  // Filter to extract the correct period
  .sort('system:time_start', true);
//compute sentinel water
print('Sentinel')
print(sentinel1.size());
var waterFreq = sentinel1.map(s1_water).sum().multiply(100.0)
  .divide(sentinel1.count()).rename('waterFreq').clip(studyArea)
Map.addLayer(waterFreq,{color:'00ff00'},'ExportImg',true);
// // 
Export.image.toAsset({
  image: waterFreq,
  description: 'S1WaterFreqFJ3',
  scale: 10,
  maxPixels:1e13
});

//-----------------------------------------------------------step3------------------------------------------------------
//loading all collections got from previous steps and computing and exporting the mangrove vectors.
var VIthreshold = 0.3;
var DEM = ee.Image('USGS/SRTMGL1_003');
var slope = ee.Terrain.slope(DEM);

//load collection got from previous steps
var mvcImage = ee.ImageCollection([image]).mosaic();
Map.addLayer(mvcImage,{min:0,max:100},'mvdImage',false);

var WaterFreq = ee.ImageCollection([image2]).mosaic();
// Map.addLayer(WaterFreq,{"palette":["ff0000"]},'WaterFreq',false);


var waterS = WaterFreq.gt(80);// Water from Sentinel-1A
Map.addLayer(waterS.updateMask(waterS),{"palette":["ff0000"]},'waterS',false);


var imgNDVIavg = mvcImage.select('imgNDVIavg');//Anual mean ndvi
var imgLSWIavg = mvcImage.select('imgLSWIavg');//Anual mean lswi
var imgFreqmNDWIWt = mvcImage.select('CntMndwi').multiply(100).divide(mvcImage.select('CntValid')); //mndwi frequency
var imgFreqBarSoRD = mvcImage.select('CntBarSoRD').multiply(100).divide(mvcImage.select('CntValid')); //bare land et al frequency
var imgFreqWetland = mvcImage.select('CntWetland').multiply(100).divide(mvcImage.select('CntValid'));//inundation(wetland) frequency
var imgFreqMangrov = mvcImage.select('CntMangrove').multiply(100).divide(mvcImage.select('CntValid'));//mangrove(canopy) frequency
var imgFreqEvgreen = mvcImage.select('CntEvgreen').multiply(100).divide(mvcImage.select('CntValid'));//evergreen frequency
var permanentWater = waterS.eq(1).and(imgFreqmNDWIWt.gt(10)); //permanent water
Map.addLayer(imgFreqmNDWIWt,{"palette":["ff0000"]},'imgFreqmNDWIWt',false);
Map.addLayer(imgFreqBarSoRD,{"palette":["ff0000"]},'imgFreqBarSoRD',false);
Map.addLayer(imgFreqWetland,{"palette":["ff0000"]},'imgFreqWetland',false);
Map.addLayer(imgFreqMangrov,{"palette":["ff0000"]},'imgFreqMangrov',false);
Map.addLayer(imgFreqEvgreen,{"palette":["ff0000"]},'imgFreqEvgreen',false);


// Maskout water and bare soil/rock/dead vegetation
var possibleMangroveRegion = permanentWater.eq(0).and(imgFreqBarSoRD.lt(50));  
// Map.addLayer(possibleMangroveRegion.updateMask(possibleMangroveRegion),{"palette":["ff0000"]},'permanentWater',true);

//greenness
var thresholdEg = imgNDVIavg.multiply(211.47).subtract(33.455);
    thresholdEg = thresholdEg.where(thresholdEg.gt(80),80);
    thresholdEg = thresholdEg.where(thresholdEg.lt(20),20);

//inundation  
var thresholdWD = imgNDVIavg.multiply(-155.81).add(114.01);
    thresholdWD = thresholdWD.where(thresholdWD.gt(70),70);
    thresholdWD = thresholdWD.where(thresholdWD.lt(0),0);

//cannopy of mangrove  
var thresholdMg = imgNDVIavg.multiply(193.85).subtract(43.784);
    thresholdMg = thresholdMg.where(thresholdMg.gt(70),70);
    thresholdMg = thresholdMg.where(thresholdMg.lt(10),10);

//Evergreen Region
var mgEvergreenRegion = imgFreqEvgreen.gt(thresholdEg).and(imgNDVIavg.gte(VIthreshold).or(imgLSWIavg.gte(VIthreshold)));
  
//Wetland region
var mgWetlandRegion = imgFreqWetland.gte(thresholdWD).and(imgNDVIavg.gte(VIthreshold).or(imgLSWIavg.gte(VIthreshold)));

// Mangrove region
var mgMangroveRegion = imgFreqMangrov.gt(thresholdMg).and(imgNDVIavg.gte(VIthreshold).or(imgLSWIavg.gte(VIthreshold)));

var Mangrove = possibleMangroveRegion.eq(1)
              .and(mgEvergreenRegion.eq(1))
              .and(mgWetlandRegion.eq(1))
              .and(mgMangroveRegion.eq(1))
              .and(DEM.lt(10)).and(slope.lt(10));

Map.addLayer(Mangrove.updateMask(Mangrove),{"palette":["ff0000"]},'Mangrove-Classified',false);

// Remove small pathes.
var mangroveConCur = Mangrove.updateMask(Mangrove).connectedPixelCount(4,false);
mangroveConCur = mangroveConCur.updateMask(mangroveConCur.gte(4)); // Remove small patches
Map.addLayer(mangroveConCur.updateMask(mangroveConCur),{"palette":["00ff00"]},'mangroveConCur');

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
var BufferPieces = ee.FeatureCollection(FJ3_COL).toList(100);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
print(BufferPieces);
var permanentWater = permanentWater.updateMask(permanentWater);
var mangroveVecList = ee.List([0]);
var first = mangroveVecList;
mangroveVecList = ee.List(BufferPieces.iterate(filterMangroveByList,first));
mangroveVecList = mangroveVecList.remove(mangroveVecList.get(0));
var mangrove = ee.FeatureCollection(mangroveVecList.flatten());
Export.table.toDrive({collection:mangrove,description:'Mangrove2015-Fj03',fileFormat:'KML'}); 

// Maping over list
function filterMangroveByList(feature,list)
{
  var curRegion = ee.Feature(feature).geometry();
  var mangroveConVecCur = mangroveConCur.reduceToVectors({
          geometry:curRegion,
          crs:mangroveConCur.projection(),
          scale:30,
          eightConnected:false,
          labelProperty:'Mangrove',
          maxPixels:1e10,
          reducer:ee.Reducer.countEvery()
        });
        
  // Get the sea vector boundary

  // Map.addLayer(permanentWater,{"palette":["ff0000"]},'permanentWater',true);

  var permanentWaterCur = permanentWater.clip(curRegion);
  permanentWaterCur = permanentWaterCur.focal_median(10,'square','meters').toInt();
  permanentWaterCur = permanentWaterCur.updateMask(permanentWaterCur);
  var permanentWaterConCur = permanentWaterCur.connectedPixelCount(512,false);
  permanentWaterConCur = permanentWaterConCur.updateMask(permanentWaterConCur.gte(512));

  // Reduce to vector
  var permanentWaterVecCur = permanentWaterConCur.reduceToVectors({
        geometry:curRegion,
        crs:permanentWaterConCur.projection(),
        scale:10,
        eightConnected:false,
        labelProperty:'Water',
        maxPixels:1e12,
        reducer:ee.Reducer.countEvery()
        });
  permanentWaterVecCur = permanentWaterVecCur.filter(ee.Filter.gt('count',10000));

  var seaBuff100Cur = permanentWaterVecCur.geometry().buffer(100);
  var mangroveVecCur = mangroveConVecCur.filterBounds(seaBuff100Cur);
  return ee.List(list).add(mangroveVecCur.toList(500));
}


