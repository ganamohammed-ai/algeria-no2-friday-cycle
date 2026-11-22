"""Python equivalents of all 13 supplied Earth Engine export scripts.

Default is an offline task-plan preview. --submit explicitly starts Drive exports
using the user's Earth Engine credentials/project. No authentication is embedded.
WorldCover v200 (2021) is preserved as used in the supplied scripts, despite the
2020 reference in the manuscript. Do not silently substitute v100.
"""
from pathlib import Path
import argparse,json,math,datetime as dt

ROOT=Path(__file__).resolve().parents[2]
CONFIG=ROOT/'metadata/earthengine_tasks.json'
SL_BANDS=['u_component_of_wind_10m','v_component_of_wind_10m','boundary_layer_height','temperature_2m','surface_pressure','surface_solar_radiation_downwards','u_component_of_wind_100m','v_component_of_wind_100m','u_component_of_wind_850hPa','v_component_of_wind_850hPa','u_component_of_wind_500hPa','v_component_of_wind_500hPa']
def jsround(x):return math.floor(x+.5)
def bounds(lat,lon):
    half=jsround(1.5/math.cos(math.radians(lat))*100)/100
    return [jsround((lon-half)*100)/100,jsround((lat-1.5)*100)/100,jsround((lon+half)*100)/100,jsround((lat+1.5)*100)/100]
def plan(config):
    tasks=[];kind=config['kind']
    if kind in ['site_covariates','control_screen']:
        return [dict(task='Site_covariates' if kind=='site_covariates' else 'Control_candidates',kind=kind,worldcover='ESA/WorldCover/v200',worldcover_year=2021)]
    for year in config['years']:
        if kind=='no2':
            for site,(lat,lon) in config['sites'].items():
                if config.get('retry') and year not in config['retry'].get(site,[]):continue
                start='2019-08-07' if year==2019 else f'{year}-01-01';end='2026-11-15' if year==2026 else f'{year+1}-01-01'
                tasks.append(dict(kind=kind,task=f'NO2_daily_{site}_{year}',site=site,lat=lat,lon=lon,year=year,start=start,end_exclusive=end,bounds=bounds(lat,lon),collection='COPERNICUS/S5P/OFFL/L3_NO2',cloud_fraction_lt=.3,column_multiplier=1e6,nodata=-9999,grid_degrees=.01,folder=config['folder']))
        elif kind=='era5':
            tasks.append(dict(kind=kind,task=f'ERA5X_SL_{year}',year=year,start=f'{year}-01-01',end_exclusive=f'{year+1}-01-01',collection='ECMWF/ERA5/HOURLY',hours_utc=[12,13],sites=config['sites'],folder=config['folder'],bands=SL_BANDS,scale_m=27830))
        elif kind=='dust':tasks.append(dict(kind=kind,task=f'S5P_AAI2_{year}',year=year,start='2019-08-07' if year==2019 else f'{year}-01-01',end_exclusive=f'{year+1}-01-01',folder=config['folder'],sites=config['sites']))
    return tasks

def no2_task(ee,t):
    collection=ee.ImageCollection(t['collection']).select(['tropospheric_NO2_column_number_density','cloud_fraction'])
    start,end=ee.Date(t['start']),ee.Date(t['end_exclusive'])
    def day(i):
        d=start.advance(i,'day')
        c=collection.filterDate(d,d.advance(1,'day')).map(lambda im:im.select('tropospheric_NO2_column_number_density').updateMask(im.select('cloud_fraction').lt(.3)))
        return ee.Image(ee.Algorithms.If(c.size().gt(0),c.mean().multiply(1e6).unmask(-9999),ee.Image.constant(-9999))).rename('no2').toFloat()
    days=ee.List.sequence(0,end.difference(start,'day').subtract(1));stack=ee.ImageCollection.fromImages(days.map(day)).toBands();x0,y0,x1,y1=t['bounds']
    return ee.batch.Export.image.toDrive(image=stack,description=t['task'],folder=t['folder'],fileNamePrefix=t['task'],region=ee.Geometry.Rectangle([x0,y0,x1,y1],'EPSG:4326',False),crs='EPSG:4326',crsTransform=[.01,0,x0,0,-.01,y1],maxPixels=1e10,fileFormat='GeoTIFF')

def era5_task(ee,t):
    features=[ee.Feature(ee.Geometry.Point([lon+dx,lat+dy]),dict(site=s,slat=lat,slon=lon)) for s,(lat,lon) in t['sites'].items() for dy in [-.25,0,.25] for dx in [-.25,0,.25]]
    pts=ee.FeatureCollection(features)
    col=ee.ImageCollection(t['collection']).filterDate(t['start'],t['end_exclusive']).filter(ee.Filter.calendarRange(12,13,'hour')).select(SL_BANDS)
    def sample(img):
        stamp=img.date().format('YYYY-MM-dd HH')
        return img.addBands(ee.Image.pixelLonLat()).reduceRegions(collection=pts,reducer=ee.Reducer.first(),scale=27830).map(lambda f:f.set('time_utc',stamp).setGeometry(None))
    fc=col.map(sample).flatten()
    return ee.batch.Export.table.toDrive(collection=fc,description=t['task'],folder=t['folder'],fileFormat='CSV',selectors=['site','slat','slon','time_utc','latitude','longitude']+SL_BANDS)

def dust_task(ee,t):
    col=ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2');start=ee.Date(t['start']);end=ee.Date(t['end_exclusive']);lists=[]
    for site,(lat,lon) in t['sites'].items():
        half=1.5/math.cos(math.radians(lat));box=ee.Geometry.Rectangle([lon-half,lat-1.5,lon+half,lat+1.5]);inner=ee.Geometry.Point([lon,lat]).buffer(10000)
        # List.map constructs the graph immediately for this site's closure.
        def feature(i):
            d=start.advance(i,'day');day=col.filterBounds(box).filterDate(d,d.advance(1,'day')).select(['absorbing_aerosol_index','cloud_fraction'])
            empty=ee.Image.constant([0,0]).rename(['absorbing_aerosol_index','cloud_fraction']).updateMask(ee.Image.constant(0))
            im=ee.Image(ee.Algorithms.If(day.size().gt(0),day.mean(),empty))
            bx=ee.Dictionary(im.reduceRegion(reducer=ee.Reducer.mean().combine(ee.Reducer.count(),'',True),geometry=box,scale=1113.2,maxPixels=1e9))
            inn=ee.Dictionary(im.select('absorbing_aerosol_index').reduceRegion(reducer=ee.Reducer.mean(),geometry=inner,scale=1113.2,maxPixels=1e9))
            return ee.Feature(None,dict(site=site,date=d.format('YYYY-MM-dd'),uvai_inner10km=inn.get('absorbing_aerosol_index',None),uvai_box=bx.get('absorbing_aerosol_index_mean',None),cloud_fraction_box=bx.get('cloud_fraction_mean',None),n_valid_uvai_pixels=bx.get('absorbing_aerosol_index_count',0)))
        lists.append(ee.List.sequence(0,end.difference(start,'day').subtract(1)).map(feature))
    fc=ee.FeatureCollection(ee.List(lists).flatten())
    return ee.batch.Export.table.toDrive(collection=fc,description=t['task'],folder=t['folder'],fileFormat='CSV',selectors=['site','date','uvai_inner10km','uvai_box','cloud_fraction_box','n_valid_uvai_pixels'])

def static_task(ee,config):
    pop=ee.ImageCollection('WorldPop/GP/100m/pop').filter(ee.Filter.eq('country','DZA')).filter(ee.Filter.eq('year',2020)).mosaic().select('population')
    lights=ee.ImageCollection('NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG').filterDate('2019-01-01','2026-01-01').select('avg_rad').mean()
    elev=ee.Image('USGS/SRTMGL1_003').select('elevation');wc=ee.ImageCollection('ESA/WorldCover/v200').first().select('Map');built=wc.eq(50).rename('built');water=wc.eq(80).rename('water')
    kind=config['kind'];out=[]
    if kind=='site_covariates':points=[(s,lat,lon) for s,(lat,lon) in config['sites'].items()]
    else:points=[(None,33.5+i*.25,-1.5+j*.25) for i in range(15) for j in range(41)]
    geoms=[ee.Geometry.Point([lon,lat]) for lat,lon in config['study_coords']]
    for site,lat,lon in points:
        g=ee.Geometry.Point([lon,lat]);d10=g.buffer(10000);d50=g.buffer(50000)
        def red(im,geom,reducer,scale):return im.reduceRegion(reducer=reducer,geometry=geom,scale=scale,maxPixels=1e10 if site else 1e9,bestEffort=True)
        em=ee.Reducer.mean().combine(ee.Reducer.stdDev(),'',True)
        if site:
            e10=red(elev,d10,em,90);e50=red(elev,d50,em,90)
            props=dict(site=site,lat=lat,lon=lon,pop_10km=red(pop,d10,ee.Reducer.sum(),100).get('population'),pop_50km=red(pop,d50,ee.Reducer.sum(),100).get('population'),nightlight_10km=red(lights,d10,ee.Reducer.mean(),500).get('avg_rad'),nightlight_50km=red(lights,d50,ee.Reducer.mean(),500).get('avg_rad'),elev_mean_10km=e10.get('elevation_mean'),elev_sd_10km=e10.get('elevation_stdDev'),elev_mean_50km=e50.get('elevation_mean'),elev_sd_50km=e50.get('elevation_stdDev'),built_share_10km=red(built,d10,ee.Reducer.mean(),100).get('built'),built_share_50km=red(built,d50,ee.Reducer.mean(),100).get('built'))
        else:
            e10=red(elev,d10,em,200)
            dist=ee.List([g.distance(s,100) for s in geoms]).reduce(ee.Reducer.min())
            props=dict(lat=lat,lon=lon,built_10km=red(built,d10,ee.Reducer.mean(),100).get('built'),built_50km=red(built,d50,ee.Reducer.mean(),500).get('built'),water_10km=red(water,d10,ee.Reducer.mean(),100).get('water'),pop_10km=red(pop,d10,ee.Reducer.sum(),100).get('population'),nightlight_10km=red(lights,d10,ee.Reducer.mean(),500).get('avg_rad'),nightlight_50km=red(lights,d50,ee.Reducer.mean(),500).get('avg_rad'),elev_mean_10km=e10.get('elevation_mean'),elev_sd_10km=e10.get('elevation_stdDev'),dist_nearest_site_km=ee.Number(dist).divide(1000))
        out.append(ee.Feature(None if site else g,props))
    name='Site_covariates' if kind=='site_covariates' else 'Control_candidates'
    return ee.batch.Export.table.toDrive(collection=ee.FeatureCollection(out),description=name,folder='NO2_GEE',fileFormat='CSV')

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--task',required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--project',default='');ap.add_argument('--submit',action='store_true');a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    configs=json.loads(CONFIG.read_text(encoding='utf-8'));config=configs[a.task];tasks=plan(config);manifest=dict(original_script=a.task,config=config,tasks=tasks,submitted=False,scope='Python translation; authenticated remote exports were not exercised in publication validation.')
    if a.submit:
        if not a.project:raise SystemExit('An Earth Engine Cloud project is required for submission.')
        import ee
        ee.Initialize(project=a.project)
        live=[]
        if config['kind'] in ['site_covariates','control_screen']:live=[static_task(ee,config)]
        else:
            for t in tasks:live.append({'no2':no2_task,'era5':era5_task,'dust':dust_task}[config['kind']](ee,t))
        for task in live:task.start()
        manifest.update(submitted=True,task_ids=[t.id for t in live])
    (a.out/(Path(a.task).stem+'_tasks.json')).write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8');print('Export tasks:',len(tasks),'submitted:',manifest['submitted'])
if __name__=='__main__':main()
