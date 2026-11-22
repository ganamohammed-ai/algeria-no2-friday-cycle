"""Native QGIS Processing entry points with an isolated, versioned scientific backend."""
from pathlib import Path
import csv,json,os,sys,subprocess,time,shutil
from qgis.core import (QgsProcessingProvider,QgsProcessingAlgorithm,QgsProcessingException,QgsProcessingParameterFile,QgsProcessingParameterString,QgsProcessingParameterBoolean,QgsProcessingParameterNumber,QgsProcessingParameterEnum,QgsProcessingParameterFolderDestination,QgsProcessingOutputString,QgsProcessingParameterFeatureSink,QgsProcessing,QgsFields,QgsField,QgsFeature,QgsGeometry,QgsPointXY,QgsCoordinateReferenceSystem,QgsWkbTypes,QgsSettings)
from qgis.PyQt.QtCore import QVariant

HERE=Path(__file__).resolve().parent
BACKEND=HERE/'backend'
CATALOG=json.loads((HERE/'catalog.json').read_text(encoding='utf-8'))

def backend_default():
    saved=QgsSettings().value('no2_reproducibility/python','')
    if saved and Path(saved).is_file():return saved
    for p in [BACKEND/'.venv/Scripts/python.exe',BACKEND/'.venv/bin/python']:
        if p.is_file():return str(p)
    for parent in HERE.parents:
        for p in [parent/'Local_QGIS_Runtime/Scripts/python.exe',parent/'Local_QGIS_Runtime/bin/python']:
            if p.is_file():return str(p)
    return ''

def execute(cmd,work,feedback,label):
    work=Path(work);work.mkdir(parents=True,exist_ok=True);log=work/(label+'.log')
    env=os.environ.copy()
    for k in ['PYTHONHOME','PYTHONPATH','PYTHONEXECUTABLE','PYTHON_EXECUTABLE','PYTHONUSERBASE','PYTHONSTARTUP','__PYVENV_LAUNCHER__','QGIS_PREFIX_PATH','QT_PLUGIN_PATH','QT_QPA_PLATFORM_PLUGIN_PATH','GDAL_DRIVER_PATH']:env.pop(k,None)
    env.update(PYTHONIOENCODING='utf-8',PYTHONDONTWRITEBYTECODE='1',MPLBACKEND='Agg',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MPLCONFIGDIR=str(work/'matplotlib_cache'))
    flags=getattr(subprocess,'CREATE_NO_WINDOW',0) if os.name=='nt' else 0
    with log.open('w',encoding='utf-8') as f:
        proc=subprocess.Popen([str(x) for x in cmd],cwd=work,env=env,stdout=f,stderr=subprocess.STDOUT,creationflags=flags)
        offset=0
        while proc.poll() is None:
            if feedback.isCanceled():
                proc.terminate();proc.wait();raise QgsProcessingException('Cancelled. Outputs already written are retained in '+str(work))
            time.sleep(.15)
            with log.open(encoding='utf-8',errors='replace') as r:
                r.seek(offset);new=r.read();offset=r.tell()
            if new:feedback.pushInfo(new.rstrip())
        if proc.returncode:
            tail=log.read_text(encoding='utf-8',errors='replace')[-6000:]
            raise QgsProcessingException('Backend failed. Log: '+str(log)+'\n'+tail)
    feedback.pushInfo('Completed: '+str(work));return log

class BaseAlgorithm(QgsProcessingAlgorithm):
    def common(self,data=True):
        self.addParameter(QgsProcessingParameterFile('PYTHON','Scientific Python executable (blank: use prepared backend)',defaultValue=backend_default(),optional=True))
        if data:self.addParameter(QgsProcessingParameterFile('DATA_ROOT','Processed data root (contains outputs and data)',behavior=QgsProcessingParameterFile.Folder,optional=True))
        self.addParameter(QgsProcessingParameterFolderDestination('OUTPUT_DIR','Independent output folder'))
        self.addOutput(QgsProcessingOutputString('LOG','Execution log'))
    def engine(self,p,c):
        v=self.parameterAsFile(p,'PYTHON',c) or backend_default()
        if not v or not Path(v).is_file():raise QgsProcessingException('Select a Python executable from the prepared scientific environment.')
        return v

class ReplayAlgorithm(BaseAlgorithm):
    def name(self):return 'reproduce_study'
    def displayName(self):return '01 Reproduce the study and verify reference results'
    def group(self):return 'Study workflows'
    def groupId(self):return 'workflows'
    def createInstance(self):return ReplayAlgorithm()
    def shortHelpString(self):return 'Runs the original estimators from archived daily tables. Quick, extended and full modes are described in the bundled README. Full ML/injection is expensive and was not fully rerun for package validation. The output directory must be separate from the input tree.'
    def initAlgorithm(self,config=None):
        self.common();self.addParameter(QgsProcessingParameterEnum('MODE','Replay scope',['quick','extended','full'],defaultValue=1))
    def processAlgorithm(self,p,c,f):
        data=self.parameterAsFile(p,'DATA_ROOT',c);out=self.parameterAsString(p,'OUTPUT_DIR',c)
        if not data or not (Path(data)/'outputs').is_dir():raise QgsProcessingException('DATA_ROOT must contain the extracted outputs folder.')
        mode=['quick','extended','full'][self.parameterAsEnum(p,'MODE',c)]
        log=execute([self.engine(p,c),'-u',BACKEND/'reproduce.py','--mode',mode,'--data-root',data,'--out',out],out,f,'qgis_reproduction')
        return {'OUTPUT_DIR':out,'LOG':str(log)}

class SetupAlgorithm(QgsProcessingAlgorithm):
    def name(self):return 'prepare_backend'
    def displayName(self):return '00 Prepare isolated scientific Python environment'
    def group(self):return 'Study workflows'
    def groupId(self):return 'workflows'
    def createInstance(self):return SetupAlgorithm()
    def initAlgorithm(self,config=None):
        self.addParameter(QgsProcessingParameterFile('BASE_PYTHON','Python 3.12 with working SSL (standalone installation or existing environment)',defaultValue=backend_default()))
        self.addParameter(QgsProcessingParameterFolderDestination('ENVIRONMENT','New environment folder (not the QGIS installation)'))
        self.addParameter(QgsProcessingParameterBoolean('EARTH_ENGINE','Also install Earth Engine API for optional remote exports',defaultValue=True))
        self.addOutput(QgsProcessingOutputString('PYTHON','Prepared Python executable'))
    def shortHelpString(self):return 'Creates a separate virtual environment and installs the tested requirements. Does not modify QGIS Python or the original project. Choose a standalone Python 3.12 or an existing scientific environment with SSL; some QGIS Python distributions cannot support a standalone SSL-enabled environment. Internet is needed. Earth Engine authentication is a separate user action.'
    def processAlgorithm(self,p,c,f):
        base=self.parameterAsFile(p,'BASE_PYTHON',c);dest=Path(self.parameterAsString(p,'ENVIRONMENT',c)).resolve()
        if not base:raise QgsProcessingException('Choose Python 3.12.')
        if dest==Path(sys.prefix).resolve() or (dest.exists() and any(dest.iterdir())):raise QgsProcessingException('Choose a new empty environment folder.')
        execute([base,'-c','import sys,ssl; assert sys.version_info[:2] == (3,12), "Python 3.12 is required"; print(ssl.OPENSSL_VERSION)'],dest.parent,f,'python_ssl_preflight')
        execute([base,'-m','venv',dest],dest.parent,f,'create_no2_environment')
        py=dest/('Scripts/python.exe' if os.name=='nt' else 'bin/python')
        requirement=BACKEND/('requirements-earthengine.txt' if self.parameterAsBool(p,'EARTH_ENGINE',c) else 'requirements-tested.txt')
        execute([py,'-m','pip','install','-r',requirement],dest.parent,f,'install_no2_requirements')
        QgsSettings().setValue('no2_reproducibility/python',str(py))
        return {'ENVIRONMENT':str(dest),'PYTHON':str(py)}

class SitesAlgorithm(QgsProcessingAlgorithm):
    def name(self):return 'sites_layer'
    def displayName(self):return '02 Create the 19 study and reference sites as a QGIS layer'
    def group(self):return 'Study workflows'
    def groupId(self):return 'workflows'
    def createInstance(self):return SitesAlgorithm()
    def initAlgorithm(self,config=None):self.addParameter(QgsProcessingParameterFeatureSink('OUTPUT','Sites',type=QgsProcessing.TypeVectorPoint))
    def processAlgorithm(self,p,c,f):
        fields=QgsFields()
        for n in ['site','group','selection_basis']:fields.append(QgsField(n,QVariant.String))
        sink,dest=self.parameterAsSink(p,'OUTPUT',c,fields,QgsWkbTypes.Point,QgsCoordinateReferenceSystem('EPSG:4326'))
        with (BACKEND/'metadata/sites.csv').open(encoding='utf-8-sig',newline='') as h:
            for r in csv.DictReader(h):
                feature=QgsFeature(fields);feature.setAttributes([r['site'],r['group'],r['selection_basis']]);feature.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(float(r['longitude_deg_east']),float(r['latitude_deg_north']))));sink.addFeature(feature)
        return {'OUTPUT':dest}

def key(flag):return 'ARG_'+flag.strip('-').upper().replace('-','_')
def automatic(script,flag,data):
    if not data:return ''
    d=Path(data);o=d/'outputs'
    defaults={'--table':o/'no2_ventilation_daily_v2.csv','--controls':o/'control/control_no2_ventilation_daily.csv','reference':o/'friday_cycle.csv','--reference':o/'friday_cycle.csv','--results':o,'--effects':o/('friday_cycle.csv' if 'friday' in script else 'lockdown_effect.csv'),'--l2':o/('l2_v240_daily_q75all.csv' if script=='lockdown_l2' else 'l2_ak_daily.csv'),'--l3':str(o/'no2_ventilation_daily_v2.csv')+';'+str(o/'control/control_no2_ventilation_daily.csv'),'--rotation':o/'wind_rotation_results.csv','--input-dir':o/'wind_rotation','--aai-dir':d/'data/04_dust_aai','CANDIDATES_POS':d/'data/05_site_descriptors/Control_candidates.csv','TABLES_POS':str(o/'no2_ventilation_daily_v2.csv')+';'+str(o/'control/control_no2_ventilation_daily.csv')}
    if flag=='table':return str(o/('study_2026/no2_ventilation_daily.csv' if 'holdout' in script else 'cities/no2_ventilation_daily.csv'))
    if flag=='--tif-dir':return str(d/'data/01_tropomi_L3_daily_geotiff')
    if flag=='--era5-dir':return str(d/'data/03_era5')
    if flag=='--input-root':return str(d)
    if flag=='--nc':return str(next((d/'data/06_edgar').glob('*.nc'),''))
    if flag=='--clips':return str(d/'data/02_tropomi_L2_orbit_clips')
    if flag=='--progress':return str(d/'data/02_tropomi_L2_orbit_clips/progress.csv')
    if script=='run_extra_detect' and flag=='--table':return str(o/'control_extra/no2_ventilation_daily.csv')
    if script=='fig_elasticity_controls':
        if flag=='--study':return str(o/'regression_results_v2_3measures.csv')
        if flag=='--controls':return str(o/'control/control_regression_results.csv')
    return str(defaults.get(flag,''))

class ScriptAlgorithm(BaseAlgorithm):
    def __init__(self,spec):super().__init__();self.spec=spec
    def name(self):return self.spec['id'].lower()
    def displayName(self):return self.spec['label']
    def group(self):return self.spec['group']
    def groupId(self):return self.spec['group'].lower().replace(' ','_')
    def createInstance(self):return ScriptAlgorithm(self.spec)
    def shortHelpString(self):
        if self.spec.get('gee'):return 'Python translation of '+self.spec['original']+'. Default generates an offline export-task manifest. Submission requires the Earth Engine API, an authenticated user and their Cloud project. Exports are sent to Drive; no remote jobs were started during validation. Static exports preserve WorldCover v200 (2021).'
        return 'Native Processing parameters for '+self.spec['script']+'. Scientific calculations reuse the versioned original backend. Blank file parameters are inferred from DATA_ROOT when possible. For multiple files/folders use semicolon-separated paths. All writes occur in the chosen output folder.'
    def initAlgorithm(self,config=None):
        self.common(data=not self.spec.get('gee'))
        if self.spec.get('gee'):
            self.addParameter(QgsProcessingParameterBoolean('SUBMIT','Start remote Earth Engine export tasks (otherwise preview only)',defaultValue=False))
            self.addParameter(QgsProcessingParameterString('PROJECT','Authenticated Earth Engine Cloud project',optional=True));return
        for a in self.spec['arguments']:
            typ=a['type'];flag=a['flag'];k=key(flag);label=flag.lstrip('-').replace('_POS','').replace('-',' ')
            if typ=='output':continue
            if typ=='bool':param=QgsProcessingParameterBoolean(k,label,defaultValue=bool(a.get('default',False)))
            elif typ=='enum':
                vals=a['choices'];default=vals.index(a['default']) if a.get('default') in vals else 0
                param=QgsProcessingParameterEnum(k,label,[str(v) for v in vals],defaultValue=default)
            elif typ in ['int','float']:
                param=QgsProcessingParameterNumber(k,label,type=QgsProcessingParameterNumber.Integer if typ=='int' else QgsProcessingParameterNumber.Double,defaultValue=a.get('default'),optional=not a.get('required',False))
            elif typ in ['file','folder']:
                param=QgsProcessingParameterFile(k,label+' (blank: auto from data root)',behavior=QgsProcessingParameterFile.Folder if typ=='folder' else QgsProcessingParameterFile.File,optional=True)
            else:
                default='' if typ=='list' else a.get('default','')
                param=QgsProcessingParameterString(k,label+(' (semicolon-separated values)' if typ=='list' else ''),defaultValue='' if default is None else str(default),optional=True)
            self.addParameter(param)
    def processAlgorithm(self,p,c,f):
        out=Path(self.parameterAsString(p,'OUTPUT_DIR',c));out.mkdir(parents=True,exist_ok=True)
        cmd=[self.engine(p,c),'-u',BACKEND/self.spec['script']]
        if self.spec.get('gee'):
            cmd+=['--task',self.spec['original'],'--out',out]
            if self.parameterAsBool(p,'SUBMIT',c):cmd+=['--submit','--project',self.parameterAsString(p,'PROJECT',c)]
        else:
            data=self.parameterAsFile(p,'DATA_ROOT',c)
            if data and self.spec['id'] in ['analyze_no2_ventilation','sensitivity_analyses']:
                # These legacy CLIs read their inputs from the output folder.
                # Stage read-only copies there instead of writing into DATA_ROOT.
                for original,target in [('no2_ventilation_daily_v2.csv','no2_ventilation_daily.csv'),('no2_daily_table.csv','no2_daily_table.csv')]:
                    source=Path(data)/'outputs'/original
                    if source.is_file() and not (out/target).exists():shutil.copy2(source,out/target)
            if data and self.spec['id']=='wind_rotation':
                stage=next(a for a in self.spec['arguments'] if a['flag']=='--stage')
                if stage['choices'][self.parameterAsEnum(p,key('--stage'),c)]=='fit':
                    source=Path(data)/'outputs/wind_rotation'
                    if source.is_dir():shutil.copytree(source,out,dirs_exist_ok=True)
            for a in self.spec['arguments']:
                typ,flag=a['type'],a['flag'];k=key(flag)
                if typ=='output':value=str(out/(self.spec['id']+'.png')) if self.spec['id'].startswith('fig_') else str(out)
                elif typ=='bool':
                    if self.parameterAsBool(p,k,c):cmd.append(flag)
                    continue
                elif typ=='enum':value=str(a['choices'][self.parameterAsEnum(p,k,c)])
                elif typ in ['int','float']:
                    if p.get(k,a.get('default')) in [None,'']:
                        if a.get('required'):raise QgsProcessingException('Missing parameter '+flag)
                        continue
                    value=str(self.parameterAsInt(p,k,c) if typ=='int' else self.parameterAsDouble(p,k,c))
                else:
                    value=self.parameterAsString(p,k,c).strip() or automatic(self.spec['id'],flag,data)
                    if not value:
                        if a.get('required'):raise QgsProcessingException('Select '+flag+' or supply a matching DATA_ROOT.')
                        continue
                if flag.startswith('--'):cmd.append(flag)
                if typ=='list':cmd += [v.strip() for v in value.split(';') if v.strip()]
                else:cmd.append(value)
        log=execute(cmd,out,f,self.name())
        manifest=out/'qgis_outputs.json';manifest.write_text(json.dumps({'algorithm':self.id(),'files':[str(x) for x in out.rglob('*') if x.is_file()]},indent=2),encoding='utf-8')
        return {'OUTPUT_DIR':str(out),'LOG':str(log)}

class NO2Provider(QgsProcessingProvider):
    def id(self):return 'no2algeria'
    def name(self):return 'NO2 Algeria Reproducibility'
    def loadAlgorithms(self):
        self.addAlgorithm(SetupAlgorithm());self.addAlgorithm(ReplayAlgorithm());self.addAlgorithm(SitesAlgorithm())
        for spec in CATALOG:self.addAlgorithm(ScriptAlgorithm(spec))
