## If you make any changes, document them here with dates #######################################################
## First version   --Zhenhua --20171002


##
################################ Only make changes to the below section ########################################################################
#input region boundary
UpLat = 40
LowLat = 35
LeftLon = -100
RightLon = -95
#input year
ini_year = 2016
end_year = 2016
#blocksize
yBlockSize = 5
xBlockSize = 5
################################ Only make changes to the above section ########################################################################
import ee
ee.Initialize()
# remove bad observations
def filterBadObs(img):
  bandminMax = img.select([0,1,2,3,4]).reduce(ee.Reducer.minMax())
  band0Mask = ee.Image(0).where(bandminMax.select([0]).eq(bandminMax.select([1])),1)
  azimuth = ee.Number(img.get('solar_azimuth_angle'))
  zenith = ee.Number(img.get('solar_zenith_angle'))
  dem = ee.Image('USGS/SRTMGL1_003')
  shadow_original = ee.Terrain.hillShadow(dem,azimuth,zenith,100) 
  shadow = shadow_original.focal_min(1,"square","pixels",1)
  cfmask = img.select(['cfmask'])
  cfmask01_shadow_band0 = ee.Image(0).where(cfmask.lt (2).bitwiseAnd(shadow.eq(1)).bitwiseAnd(band0Mask.eq(0)),1)
  result = img.mask(cfmask01_shadow_band0)
  return result
# vegetation calculation
def getVIs_Shadow(img):
  ndvi = img.expression('(nir - red)/(nir + red)',
          {
            'red':img.select(['red']),
            'nir':img.select(['nir'])
          })
  mndwi = img.expression('(green - swir1)/(green + swir1)',
          {
            'green':img.select(['green']),
            'swir1':img.select(['swir1'])
          })  
  evi = img.expression('2.5*(nir - red)/(nir + (6.0 * red) - (7.5 * blue) + 10000.0)',
          {
            'red': img.select(['red']),
            'nir': img.select(['nir']),
            'blue': img.select(['blue'])
          })
  return img.select([]).addBands([ndvi,evi,mndwi]).select([0,1,2],['ndvi','evi','mndwi']).copyProperties(img,img.propertyNames())
#Water detection
def maskPosNegWater(img):
  mndwi = img.select(['mndwi'])
  ndvi = img.select(['ndvi'])
  evi = img.select(['evi'])
  maskPos = ee.Image(0).where(evi.lt(0.1).bitwiseAnd(mndwi.gt(evi).bitwiseOr(mndwi.gt(ndvi))),1)
  maskNeg = ee.Image(0).where(evi.gte(0.1).bitwiseOr(mndwi.lte(evi).bitwiseAnd(mndwi.lte(ndvi))),1)
  return img.select([]).addBands([maskPos,maskNeg]).select([0,1],['Pos','Neg']).copyProperties(img,img.propertyNames())

############## main function #######################################
def calculate_and_output(LeftLoni,LowLati,RightLoni,UpLati,yeari):
  # get the time 
  startDate = str(yeari)+'-01-01'   # change this line
  endDate = str(yeari+1)+'-01-01'      # change this line
  region = ee.Geometry.Rectangle(LeftLoni,LowLati,RightLoni,UpLati);
  region_Bound = str(region.bounds(0.001,'EPSG:4326').getInfo()['coordinates'][0])
  # set output name
  freqName = 'Freq_usgs578_cfmask01_TShadow_band0_'+str(yeari)+'_'+str(LeftLoni).zfill(4)+'_'+str(LowLati).zfill(3)
  qualName = 'Qual_usgs578_cfmask01_TShadow_band0_'+str(yeari)+'_'+str(LeftLoni).zfill(4)+'_'+str(LowLati).zfill(3)
  # Get image collection
  ImageCollectionL5 = ee.ImageCollection('LANDSAT/LT5_SR').filterBounds(region).filterDate(startDate,endDate).sort('system:time_start').select(['B1','B2','B3','B4','B5','B7','cfmask','cfmask_conf'],['blue','green','red','nir','swir1','swir2','cfmask','cfmask_conf'])
  ImageCollectionL7 = ee.ImageCollection('LANDSAT/LE7_SR').filterBounds(region).filterDate(startDate,endDate).sort('system:time_start').select(['B1','B2','B3','B4','B5','B7','cfmask','cfmask_conf'],['blue','green','red','nir','swir1','swir2','cfmask','cfmask_conf'])
  ImageCollectionL8 = ee.ImageCollection('LANDSAT/LC8_SR').filterBounds(region).filterDate(startDate,endDate).sort('system:time_start').select(['B2','B3','B4','B5','B6','B7','cfmask','cfmask_conf'],['blue','green','red','nir','swir1','swir2','cfmask','cfmask_conf'])
  ImageCollection = ee.ImageCollection(ImageCollectionL5.merge(ImageCollectionL7).merge(ImageCollectionL8))
  print('ImageCollection: ',ImageCollection)
  # filter bad observations
  GoodObsers = ImageCollection.map(filterBadObs) 
  print('GoodObsers: ',GoodObsers)
  # calculate vegetation indices
  VIs = GoodObsers.map(getVIs_Shadow)
  print('VIs: ',VIs)
  # water detection
  WaterPosNeg = VIs.map(maskPosNegWater)
  print('WaterPosNeg: ',WaterPosNeg)
  # water frequency and quality band (number of good observations)
  total_Pos = WaterPosNeg.select(['Pos']).sum()
  total_Neg = WaterPosNeg.select(['Neg']).sum()
  flooding_freq = ee.Image(0).expression('((Pos)/(Pos + Neg))*254+1',
                        {
                            'Pos':total_Pos.select([0]),
                            'Neg':total_Neg.select([0])
                        })
  goodGuals = ee.Image(0).expression('(Pos + Neg)',
                        {
                            'Pos':total_Pos.select([0]),
                            'Neg':total_Neg.select([0])
                        })
  freqTask = ee.batch.Export.image.toDrive(image=flooding_freq.toUint8(), crs='EPSG:4326', description=freqName,region=region_Bound, maxPixels=100000000000,scale=30,folder='Water')
  qualTask = ee.batch.Export.image.toDrive(image=goodGuals.toInt16(), crs='EPSG:4326', description=qualName,region=region_Bound,maxPixels=100000000000,scale=30,folder='Water')
  freqTask.start()
  qualTask.start()
################################################################################################################################################
######### use the main function  ################ don't chang the section below
for LowLati in range(LowLat,UpLat,yBlockSize):
  if LowLati + yBlockSize < UpLat:
    UpLati = LowLati + yBlockSize
  else:
    UpLati = UpLat
  for LeftLoni in range(LeftLon,RightLon,xBlockSize):
    if LeftLoni + xBlockSize < RightLon:
      RightLoni = LeftLoni+xBlockSize
    else:
      RightLoni = RightLon
    for yeari in range(ini_year,end_year+1):
      calculate_and_output(LeftLoni,LowLati,RightLoni,UpLati,yeari)
