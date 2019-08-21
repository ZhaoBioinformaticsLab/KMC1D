import os,sys,re, shutil
from flask import request, url_for, render_template, abort, session
from settings import *
from werkzeug import secure_filename
from settings import UPLOAD_FOLDER
from webapp import app

if not os.path.exists(UPLOAD_FOLDER):	
    os.makedirs(UPLOAD_FOLDER, 0777)    

@app.route('/simpleupload', methods=['GET','POST'])
def simpleupload():
    if request.method == 'POST':
        file = request.files['file']
        if file:
            if not os.path.isdir(temp_dir):
                os.makedirs(temp_dir, 0777)
            filename = secure_filename(file.filename)
            fullfilename=os.path.join(UPLOAD_FOLDER, filename)
            file.save(fullfilename)
            return fullfilename

###HTML5 resmable upload
@app.route('/upload', methods=['GET'])
def upload_query():
    #resumable upload function
    identifier = request.args.get('resumableIdentifier', type=str)
    filename = request.args.get('resumableFilename', type=str)
    chunk_number = request.args.get('resumableChunkNumber', type=int)

    if not identifier or not filename or not chunk_number:
        # Parameters are missing or invalid
        abort(500, 'Parameter error')

    # chunk folder path based on the parameters
    temp_dir = "{}/{}.dir".format(UPLOAD_FOLDER, identifier)
    # chunk path based on the parameters
    chunk_file = "{}/part.{}".format(temp_dir, chunk_number)

    app.logger.debug('Getting chunk: %s', chunk_file)

    if os.path.isfile(chunk_file):
        #Let resumable.js know this chunk already exists
        return 'Found chunk'     ### 'OK' is data body, status code is 200
    else:
        #Let resumable.js know this chunk does not exists and needs to be uploaded
        return 'No chunk', 204

def argvalue(keyname,defvalue=None,type=str):
    ret = request.args.get(keyname,None,type)
    if ret!=None:
        return ret
    ret = defvalue
    try:
        ret=type(request.form[keyname])
    except:
        pass
    return ret

###HTML5 resmable upload
@app.route('/upload', methods=['POST'])
def upload_create():
    #resumable upload function, for method="octet" only
    identifier = request.args.get('resumableIdentifier', type=str)
    filename = request.args.get('resumableFilename', type=str)
    chunk_number = request.args.get('resumableChunkNumber', type=int)
    chunk_size = request.args.get('resumableChunkSize', type=int)
    current_chunk_size = request.args.get('resumableCurrentChunkSize', type=int)
    total_size = request.args.get('resumableTotalSize', type=int)
    total_chunks = request.args.get('resumableTotalChunks', type=int)
    # or not chunk_size or not current_chunk_size or not total_size:
    if not identifier or not filename or not chunk_number or not chunk_size or not current_chunk_size or not total_size or not total_chunks:
        # Parameters are missing or invalid
        abort(400, 'Parameter error')

    # chunk folder path based on the parameters
    temp_dir = "{}/{}.dir".format(UPLOAD_FOLDER, identifier)
    # chunk path based on the parameters
    chunk_file = "{}/part.{}".format(temp_dir, chunk_number)

    app.logger.debug('Creating chunk: %s', chunk_file)

    try:
        # Create directory of not exists
        if not os.path.isdir(temp_dir):
            os.makedirs(temp_dir, 0777)
        
        with open(chunk_file,'wb') as f:
            f.write(request.data)
 
    except Exception as e:
        abort(500, 'Failed to save chunk')

    return 'OK'


def merge_slice(identifier, nomerge=False):
    m=re.search('^([0-9]+)-([0-9]+)$',identifier)  ##check format of identifier
    if not m:
        raise Exception('Invalid identifier')
    sizeonclientside=int(m.group(2))
    
    target_file_name = "{}/{}".format(UPLOAD_FOLDER, identifier)
    if os.path.isfile(target_file_name):
        if os.path.getsize(target_file_name)==sizeonclientside:
            return target_file_name 
        else:
            raise Exception('The size of existing %s is not equal to the one on client side' % identifier)

    temp_dir = "{}/{}.dir".format(UPLOAD_FOLDER, identifier)     ##temp folder for slice files
    if not os.path.exists(temp_dir):
        raise Exception('Couldn\'t find slice temp folder for resumable identifier {}'.format(identifier))
    if nomerge:
        return temp_dir

    with open(target_file_name, "wb") as target_file:
        chunk_number=0
        while True:
            chunk_number+=1
            stored_chunk_file_name = "{}/part.{}".format(temp_dir, str(chunk_number))
            if not os.path.isfile(stored_chunk_file_name):
                break
            stored_chunk_file = open(stored_chunk_file_name, 'rb')
            # Write chunk into target file
            target_file.write( stored_chunk_file.read() )
            stored_chunk_file.close()

    if os.path.getsize(target_file_name)!=sizeonclientside:
        raise Exception('The size of uploaded file is not equal to the one on client side')
    shutil.rmtree(temp_dir)
    #print target_file_name
    #sys.stdout.flush()
    return target_file_name            
