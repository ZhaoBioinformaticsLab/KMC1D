import tempfile, sqlite3, traceback, re, os, sys, subprocess
from flask import Flask, Response, render_template, make_response, request, g, session, redirect, url_for, render_template_string, abort, flash, send_file
from settings import *
from biotools import ApprunClient, UniqSessionID2YearMonth, sizeof_fmt
from werkzeug import secure_filename
from os.path import isfile, join
from webapp import app
from collections import OrderedDict
from mailportal import Gmail
from upload import merge_slice
import urllib

@app.before_request
def before_request():
    g.apprun = ApprunClient(APP_ROOT, APPRUN_HOST, APPRUN_PORT, APPRUN_PIN)
    g.app_name = APP_NAME
    g.app_desc = APP_DESC
    g.app_root = APP_ROOT
    g.app_templates = APP_TEMPLATES
    
@app.route('/home')
def home():
    params={}
    params['entryname']='Home'
    return render_template('home.html', **params)
    
@app.route('/download')
def download():
    params={}
    params['entryname']='Download'
    return render_template('download.html', **params)

@app.route('/help')
def help():
    params={}
    params['entryname']='Help'
    return render_template('help.html', **params)

@app.route('/email')
def email():
    params={}
    params['entryname']='Contact us'
    if request.args.get('dosubmit'):
        gmail=Gmail()
        toaddrs  = EMAILRECPT
        email=request.args.get('email','').strip()
        if email:
            replytos = [email]
        else:
            replytos = []
        subject=request.args.get('subject','').strip()
        content=request.args.get('content','').strip()
        gmail=Gmail()
        try:
            gmail.sendmail(toaddrs, replytos, [], subject, content)
        except:
            params['failedtosend']=True
        else:
            params['done']=True
    return render_template('email.html', **params)

@app.route('/info')
def info():
    params={}
    params['entryname']='Info'
    return render_template('info.html', **params)


def downloadfile(fullfilename, remotefilename):
    def generate():
        with open(fullfilename) as f:
            for line in f:
                yield line
    headers = {"Content-Disposition": "attachment; filename=%s" % remotefilename}
    return Response(generate(), mimetype='text/plain', headers=headers)


@app.route('/', methods=['GET', 'POST'])
@app.route('/analysis', methods=['GET', 'POST'])
def analysis():
    params={'HELPTIP':HELPTIP}
    try:
        simplepage=elementValue('simplepage',None)
        if simplepage==None and 'simplepage' not in session:
            session['simplepage']=False
        if simplepage!=None:
            if simplepage in ['0','1']:
                session['simplepage']=bool(int(simplepage))
            else:
                raise Exception('simplepage paramter can only be 0 or 1')
                        
        #calculate switch html5/legacy browser
        switchurl=url_for('analysis', simplepage=(0 if session['simplepage'] else 1))
        switchtip ='HTML5' if session['simplepage'] else "Legacy"
        if session['simplepage']:
            tooltip="You are on traditional uploading page compatible with legacy browser. Click here switch to modern HTML5 based page which capable of uploading large file."
        else:
            tooltip="You are on HTML5 uploading page capable of uploading large file. Click here if you are using legacy browser(e.g. IE<9) or having compatibility issue with HTML5"
        params['switchstr']=' [<a data-placement=\'bottom\' data-toggle=\'tooltip\'  title=\'%s\' href=\'%s\'>%s browser</a>]' % (tooltip, switchurl, switchtip)

        params['entryname']='Analysis'
        if request.method != "POST":
            return analysis_input(params)
        cmdlist, numq, sessionid = makecmd()
        print(cmdlist)
        cmdstr=" ".join(cmdlist)
        print cmdstr
        #sessionrelpath=os.path.join( UniqSessionID2YearMonth(sessionid),  sessionid] )
        cleancmd=""   #"%s/cleansession.sh %s" % (APP_BIN, sessionrelpath)   ###clean session folder except session.props for life span
        params['sessionid']=g.apprun.run(cmdstr, sessionid, 3600*24*APRRUN_LIFESPAN, 9, APP_NAME+'_q%s' % numq, cleancmd)
        params['message']='Waiting for response of backend queue management......'
        params['waittime']=2000   ##in milli seconds
        params['newurl']=url_for('result',sessionid=sessionid)
        return render_template("show_progress.html", **params)
    except Exception as e:
        traceback.print_exc(file=sys.stdout)
        params['message']=str(e)
        return render_template("prompt.html", **params)

def analysis_input(params):
 
    params['GENO_UPLIMIT']=kMC1D_LIMITS['GENO_UPLIMIT']   
    params['geno_uplimit_str']=sizeof_fmt(kMC1D_LIMITS['GENO_UPLIMIT'])
    params['MIN_ACCESSION_NUM']=kMC1D_LIMITS['MIN_ACCESSION_NUM']
    params['MAX_ACCESSION_NUM']=kMC1D_LIMITS['MAX_ACCESSION_NUM']
    params['MIN_MARKER_NUM']=kMC1D_LIMITS['MIN_MARKER_NUM']
    params['MAX_MARKER_NUM']=kMC1D_LIMITS['MAX_MARKER_NUM']

    return render_template('analysis.html', **params)

@app.route('/result', methods=['GET', 'POST'])
def result():
    params={'HELPTIP':HELPTIP}
    params['entryname']='Analysis result'
    params['sessionid']=request.args.get('sessionid')
    if not params['sessionid']:
        params['message']='Please give sessionid'
        return render_template("prompt.html", **params)
    #params['sessionpath']=os.path.join(APP_SESSIONBASE, UniqSessionID2YearMonth(params['sessionid']),  params['sessionid']) 
    try:  #test if it is valid apprun session
        status = g.apprun.getStatu(params['sessionid'])           ## if sessioid is invalid, will throw exception
        tmpcmd=g.apprun.getSessionProps(params['sessionid'])[7]
        if not re.search(MAIN_PIPELINE,tmpcmd):
            params['message']='This is not a %s analysis session' % APP_NAME
            return render_template("prompt.html", **params)
    except Exception as e:  # it is not an effective apprun session
        params['message']='Failed to find the session!'
        return render_template("prompt.html", **params)
    else:  #as an effective apptun session
        if status in [ApprunClient.RUNNING, ApprunClient.PENDING]:
            params['waittime']=5000   ##in milli-seconds
            params['newurl']=url_for(request.endpoint,sessionid=params['sessionid'])
            if status == ApprunClient.RUNNING:
                progressline=None
                for x in g.apprun.getFullerr(params['sessionid']).split("\n"):
                    if not x:
                        continue
                    if x.startswith('Stage'):
                        params['message']=x
                        progressline=None
                    if x.startswith('Finished:') or x.startswith('Iteration'):
                        progressline=x
                if progressline:
                    m=re.search('Finished:\s+([0-9.]+)\s+%;\s*running/total\s*nodes:\s*([0-9]+)/([0-9]+)\s*$', progressline)   ##Finished:    99 %; running/total nodes: 53/662
                    if m:
                        params['message']=params['message']+"   "+  "{}% finished; {} of total {} nodes.".format(m.group(1),m.group(2),m.group(3))
                    else:
                        params['message']=params['message']+"   "+progressline
            else:
                params['message']='%i job(s) ahead of you in queue' % g.apprun.getnumofsessionsaheadof(params['sessionid'])
            return render_template("show_progress.html", **params)
        elif status == ApprunClient.FAILED:
            fullerr=g.apprun.getFullerr(params['sessionid'])
            for x in fullerr.split('\n'):
                if not x:
                    continue
                if re.search(u'(?i)^ERROR',x):
                    params['message']=x
            return render_template("prompt.html", **params)
        elif status == ApprunClient.FINISHED:
            return display_result(params)

def display_result(params):
    filename=None    

    properties=dict()
    for line in g.apprun.getFullout(params['sessionid']).split("\n"):
        m=re.search(u'^(\S+)=(\S+)',line)
        if m:
            properties[m.group(1)]=m.group(2)

    if request.args.get('contentype','')=='csv':
        result=elementValue("result",None,type=str)
        if result:
            fn="%s_%s.csv" % (result, params['sessionid'])
        else:
            params['message']='Please specify result'
            return render_template("prompt.html", **params)
        if result not in properties:
            params['message']='Need a valid result'
            return render_template("prompt.html", **params)
        return send_file(properties[result], as_attachment=True, attachment_filename=fn)
    else:
#       url=url_for('result', sessionid=params['sessionid'], contentype='csv', result='KinshipMatrix1D')
#       params['kinshipmatrixurl'] = urllib.pathname2url(url)
        return render_template("result_display.html", **params)

def elementValue(keyname,defvalue=None,type=str):
    if request.method=='GET':
        return request.args.get(keyname,defvalue,type)
    ret= defvalue
    try:
        ret=type(request.form[keyname])
    except:
        pass
    return ret

def getFileLine(filename):
    try:
        cc=sum(1 for line in open(filename))
    except:
        return None
    else:
        return cc


def makecmd(): 
    #python main_pipeline.py  --sessionid test --phenofile demo_data3/pheHybrid1500.csv --genofile demo_data3/gen.csv  --phenocolname KGW --randomcolname Hybrid 
    #--varcovarfile demo_data3/gmat.csv --fixedeffect1 Location-c --fixedeffect2 YD-n
    cmd=[MAIN_PIPELINE]

    genofile, uploadsize=treatupload2('geno', kMC1D_LIMITS['GENO_UPLIMIT'], 'genotype')
   
    if not genofile :
        raise Exception('Please check if the genotypic file is in a correct form.')
    if genofile:
        cmd.append("--genofile")
        cmd.append(genofile)
    
    '''
    #pre-check genotype file format
    
    if genofile:
        if os.path.isdir(genofile):
            genofile=os.path.join(genofile,'part.1')
        with open(genofile) as f:
            accession_columns =f.readline().strip().split(',')
            accessionnum      =len(accession_columns)
            
            if(accessionnum<kMC1D_LIMITS['MIN_ACCESSION_NUM']):
				accession_columns =f.readline().strip().split('\t')
				accessionnum      =len(accession_columns)
            
            #print accessionnames,accessionnum,MIN_ACCESSION_NUM
            if(accessionnum<kMC1D_LIMITS['MIN_ACCESSION_NUM']):
                raise Exception("Only %s accession columns in your genotype file. You need %s accessions at least" % (accessionnum, kMC1D_LIMITS['MIN_ACCESSION_NUM']))
            if(accessionnum>kMC1D_LIMITS['MAX_ACCESSION_NUM']):
                raise Exception("At most %s columns (accessions) are allowed in your genotype csv file. You have submitted %s accessions" % (kMC1D_LIMITS['MAX_ACCESSION_NUM'], accessionnum))
        markernum=getFileLine(genofile)
        if(markernum<kMC1D_LIMITS['MIN_MARKER_NUM']):
            raise Exception("Your genotypic file must have at least %s markers." % kMC1D_LIMITS['MIN_MARKER_NUM'])
        if(markernum>kMC1D_LIMITS['MAX_MARKER_NUM']):
            raise Exception("Your genotypic file must have at maximum %s markers." % kMC1D_LIMITS['MAX_MARKER_NUM'])
            
    cmd.append("--accessionnum");
    cmd.append(str(accessionnum));
    '''
        
    argumentname='Marker_Block'
    tmpv=elementValue(argumentname,10000,type=int)
    if tmpv>=10 :
        cmd.append("--markerblock")
        cmd.append(str(tmpv))
    else:
        raise Exception('Invalid marker block size. The markers are divided into blocks for GPU paralleing computing')
        
    cmd.append("--runmode");
    cmd.append("0");  # 0,1,2 for GPU Parallel,  CPU Serial, and Read Through 
            
    if "X-Real-IP" in request.headers:
        remoteaddr=request.headers["X-Real-IP"]
    else:
        remoteaddr=request.remote_addr
    cmd.append("--remoteaddr")
    cmd.append(remoteaddr)

    sessionid=g.apprun.getuniqsessionid()
    cmd.append("--sessionid")
    cmd.append(sessionid)

    workload= uploadsize
    workloadbase = 1000*1000
    numq=1
    if (workload >= 10*workloadbase) :
        numq=2
    if (workload >= 200*workloadbase) :
        numq=3
    if(workload >= 20000*workloadbase) :
        numq=4

    return (cmd,numq,sessionid)    

def treatupload2(element_name, uplimit, tip="content"):  #HTML5 resumable compatible
    if request.method != 'POST':
        return None
    
    ##upload by text area in web page
    con=elementValue(element_name+"_content",'').strip()
    if con!="":
        #if len(con) == 0:
        #    raise Exception('Uploaded %s is empty' % tip)
        if len(con) > uplimit:
            raise Exception('Error: Uploaded %s is more than %s ' % (tip, uplimit))
        f=tempfile.NamedTemporaryFile(dir=UPLOAD_FOLDER,delete=False,prefix="pythontmp_")
        f.write(con)
        fullfilename=f.name
        f.close()
        os.chmod(fullfilename,0777)
        return fullfilename, len(con)

    ###for resumable upload
    uploaded_fn=elementValue(element_name+"_uploaded",'').strip()   ###pattern of uploaded_fn:  1524609577762-119093422
    if uploaded_fn!="":
        m=re.search('^([0-9]+)-([0-9]+)$',uploaded_fn)  ##check format of identifier
        if not m:
            raise Exception('Invalid resumerable upload identifier')
        sizeonclientsize=int(m.group(2))
        if sizeonclientsize > uplimit:
            raise Exception('Error: Uploaded %s is more than %s ' % (tip, uplimit))
        if sizeonclientsize > 200*1000*1000:
            fullfilename=merge_slice(uploaded_fn, nomerge=True)   ###nomerge=True will only return slices folder name
            subprocess.call(["chmod", "-Rf", '0777', fullfilename])
        else:
            fullfilename=merge_slice(uploaded_fn)   ###nomerge=True will only return slices folder name
            os.chmod(fullfilename,0777)
        return fullfilename, sizeonclientsize
    
    ###for traditional upload
    myfile = request.files.get(element_name)
    if myfile:
        f=tempfile.NamedTemporaryFile(dir=UPLOAD_FOLDER,prefix="pythontmp_")
        fullfilename=f.name
        f.close()
        myfile.save(fullfilename)
        if os.path.getsize(fullfilename) > uplimit:
            os.remove(fullfilename)
            raise Exception('Uploaded %s is more than %s ' % (tip, uplimit))
        os.chmod(fullfilename,0777)
        return fullfilename, os.path.getsize(fullfilename)

    return None 


def treatupload(file_element_name, content_element_name, uplimit, tip="content"):
    if request.method == 'POST':
        myfile = request.files[file_element_name]
        if myfile:
            f=tempfile.NamedTemporaryFile(dir=UPLOAD_FOLDER,prefix="pythontmp_")
            fullfilename=f.name
            f.close()
            myfile.save(fullfilename)
            if os.path.getsize(fullfilename) > uplimit:
                os.remove(fullfilename)
                raise Exception('Uploaded %s is more than %s ' % (tip, uplimit))
            return fullfilename
        else:
            con=elementValue(content_element_name,'')
            if len(con) == 0:
                raise Exception('Uploaded %s is empty' % tip)
            if len(con) > uplimit:
                raise Exception('Uploaded %s is more than %s ' % (tip, uplimit))
            f=tempfile.NamedTemporaryFile(dir=UPLOAD_FOLDER,prefix="pythontmp_")
            f.write(con)
            fullfilename=f.name
            f.close()
            return fullfilename
            
