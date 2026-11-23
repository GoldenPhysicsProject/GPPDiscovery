"""Download Platt's zero tables from the shared Drive folder into ZETA_FOLDER/data and platt_zeros.py next to this
script's workdir.  Uses the service account in $GDRIVE_SA_JSON when set, otherwise gdown (folder shared by link).
Usage: python fetch_zeros.py ZETA_FOLDER FOLDER_ID PLATT_PY_OUT"""
import os
import sys
import json
import shutil
from pathlib import Path

zeta = Path(sys.argv[1]); fid = sys.argv[2]; platt_out = Path(sys.argv[3])
data = zeta / 'data'
data.mkdir(parents=True, exist_ok=True)
expect = int(os.environ.get('EXPECT_DAT', '0'))   # number of zeros_*.dat files the folder is known to hold (strict check)
if True:
    sa = os.environ.get('GDRIVE_SA_JSON')
    if sa:
        from google.oauth2 import service_account
        from googleapiclient.discovery import build
        from googleapiclient.http import MediaIoBaseDownload
        creds = service_account.Credentials.from_service_account_info(
            json.loads(sa), scopes=['https://www.googleapis.com/auth/drive.readonly'])
        svc = build('drive', 'v3', credentials=creds, cache_discovery=False)
        files, tok = [], None
        while True:
            r = svc.files().list(q=f"'{fid}' in parents and trashed=false", fields='nextPageToken,files(id,name,size)',
                                 pageToken=tok, pageSize=200).execute()
            files += r['files']; tok = r.get('nextPageToken')
            if not tok: break
        dats = [f for f in files if f['name'].endswith('.dat')]
        print('service account lists', len(dats), 'zeros_*.dat files:', sorted(f['name'] for f in dats), flush=True)
        if expect and len(dats) != expect:
            sys.exit(f'INCOMPLETE LISTING: expected {expect} .dat files, the service account sees {len(dats)} '
                     '(Drive may still be propagating the folder permission to its children; re-run in a few minutes)')
        for f in files:
            if not (f['name'].endswith('.dat') or f['name'] == 'platt_zeros.py'):
                continue
            dest = platt_out if f['name'] == 'platt_zeros.py' else data / f['name']
            if dest.exists() and dest.stat().st_size == int(f['size']):
                continue
            with open(dest, 'wb') as fh:
                dl = MediaIoBaseDownload(fh, svc.files().get_media(fileId=f['id']), chunksize=64 << 20)
                done = False
                while not done:
                    _, done = dl.next_chunk()
            print('fetched', f['name'], dest.stat().st_size, flush=True)
    else:
        import gdown
        tmp = zeta / 'dl'
        tmp.mkdir(exist_ok=True)
        gdown.download_folder(id=fid, output=str(tmp), quiet=True, use_cookies=False)
        for p in tmp.rglob('*'):
            if p.is_file():
                if p.name == 'platt_zeros.py':
                    shutil.move(str(p), platt_out)
                elif p.suffix == '.dat':
                    shutil.move(str(p), data / p.name)
        shutil.rmtree(tmp, ignore_errors=True)
n_local = len(list(data.glob('zeros_*.dat')))
print('data files:', n_local, 'platt_zeros.py:', platt_out.exists())
if expect and n_local != expect:
    sys.exit(f'INCOMPLETE DOWNLOAD: {n_local} of {expect} .dat files present')
