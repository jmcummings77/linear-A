#!/usr/bin/env python3
"""Render saved benchmark snapshots for descriptive comparison, without measuring."""
import argparse
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
from publication import PublicSanitizer
ROOT=Path(__file__).resolve().parents[1]
HERE=Path(__file__).resolve().parent/'comparison'

def select_public(data):
    if not isinstance(data,dict) or data.get("schema_version") != 1 or not isinstance(data.get("results"),list) or not isinstance(data.get("implementations"),list):
        raise ValueError("Expected a schema-version-1 benchmark report")
    sanitizer=PublicSanitizer(root=ROOT)
    # Explicit schema: no environment, profile payload, runner command or diagnostics.
    report={key:data.get(key) for key in ('schema_version','revision','source_sha256','dirty','created_at','suite','seed','machine','methodology','implementations','results')}
    report=sanitizer.report(report)
    report['machine']={k:report.get('machine',{}).get(k) for k in ('os','release','architecture','logical_cpus')}
    report['methodology']={k:report.get('methodology',{}).get(k,'') for k in ('numeric_type','matrix_layout','timing','warmup','sampling','sample_unit','workload_version')}
    report['implementations']=[{k:i.get(k) for k in ('id','name','status','toolchain','build_commands')} for i in report.get('implementations',[])]
    report['results']=[{**{k:r.get(k) for k in ('implementation','operation','size','status')},'samples':[{k:s.get(k) for k in ('elapsed_ns','iterations','ns_per_op')} for s in r.get('samples',[])]} for r in report.get('results',[])]
    return report

def git_snapshot(ref,path):
    relative=PurePosixPath(path)
    if relative.is_absolute() or '..' in relative.parts or '\\' in path or not path.endswith('.json'):raise ValueError('Use a repository-relative JSON path')
    revision=subprocess.check_output(['git','rev-parse','--verify','--end-of-options',ref+'^{commit}'],cwd=ROOT,text=True).strip()
    content=subprocess.check_output(['git','show',revision+':'+path],cwd=ROOT,text=True)
    return json.loads(content),revision

def render(snapshots,destination):
    if not 2<=len(snapshots)<=30:raise ValueError('Need 2–30 snapshots')
    data=[{'label':PublicSanitizer(root=ROOT).text(label),'data':select_public(report)} for label,report in snapshots]
    encoded=json.dumps(data,allow_nan=False).replace('<','\\u003c').replace('&','\\u0026')
    # Embed the pure model and UI together, with no runtime network dependencies.
    model=(HERE/'model.mjs').read_text().replace('export ','')
    replacements={'__SNAPSHOTS__':encoded,'__SCRIPT__':model+'\n'+(HERE/'app.mjs').read_text()}
    html=re.sub('|'.join(replacements),lambda m:replacements[m[0]],(HERE/'template.html').read_text())
    destination.parent.mkdir(parents=True,exist_ok=True);destination.write_text(html)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('files',type=Path,nargs='*')
    parser.add_argument('--from-git',nargs=2,metavar=('BASE_REF','CANDIDATE_REF'))
    parser.add_argument('--path',default='benchmarks/reports/latest/results.json')
    parser.add_argument('--catalog',action='store_true')
    parser.add_argument('--output',type=Path,default=HERE/'index.html')
    args=parser.parse_args();snapshots=[]
    if sum((bool(args.files),bool(args.from_git),args.catalog))!=1:parser.error('Choose files, --from-git, or --catalog')
    if args.from_git:
        for ref in args.from_git:
            data,containing=git_snapshot(ref,args.path);snapshots.append(('Saved in '+containing[:12],data))
    else:
        paths=sorted((ROOT/'benchmarks/reports').glob('*/results.json')) if args.catalog else args.files
        for path in paths:
            data=json.loads(path.read_text())
            if args.catalog and not data.get('results'):continue
            snapshots.append(((path.parent.name+' / '+path.name) if args.catalog else 'Snapshot '+str(len(snapshots)+1),data))
    render(snapshots,args.output)
    print('Rendered %d snapshots; no benchmarks executed'%len(snapshots))
if __name__=='__main__':main()
