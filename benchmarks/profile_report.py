#!/usr/bin/env python3
"""Add separately captured stack profiles to a saved report, preserving its timings."""
import argparse
import datetime
import json
from pathlib import Path
import subprocess

from profiles import collect_profiles
from publication import PublicSanitizer
from report import render
from run import ROOT, build_one, fingerprint, verify
from vector_reference import vector_fixtures
from reference import fixtures


def add_profiles(path, operation, size, no_build=False):
    path=Path(path)
    data=json.loads(path.read_text())
    expected_size={'cross':3,'rotation2d':2,'rotation3d':3}.get(operation)
    if operation not in ('multiply','cross','rotation2d','rotation3d') or type(size) is not int or size<1 or (expected_size and size!=expected_size):
        raise ValueError('invalid profile workload or size')
    if not any(r.get('operation')==operation and r.get('size')==size for r in data.get('results',[])):
        raise ValueError('choose a workload present in the saved report')
    implementations=[]
    for item in data['implementations']:
        if item['status']!='passed':continue
        print('Preparing profile: '+item['id'],flush=True)
        implementation=build_one(item['id'],no_build)
        cases=vector_fixtures() if expected_size else fixtures()
        if not verify(implementation,cases):raise RuntimeError('profile correctness checks failed: '+item['id'])
        implementations.append(implementation)
    if not implementations:raise ValueError('no verified implementations to profile')
    provenance={
        'created_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'dirty':bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True)),
        'source_sha256':fingerprint(),
    }
    profiles=collect_profiles(implementations,path.parent,data['seed'],operation,size)
    toolchains={i['id']:i['toolchain'] for i in implementations}
    for profile in profiles:
        profile.update(provenance,toolchain=toolchains[profile['implementation']])
    data['profiles']=PublicSanitizer().report({'profiles':profiles})['profiles']
    # Only replace profile data: the original run's timing/provenance fields stay intact.
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    render(data,path.with_name('index.html'))
    print('%d/%d profiles available'%(sum(p['status']=='available' for p in profiles),len(profiles)),flush=True)
    return profiles


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report',type=Path,help='saved results.json to enrich')
    parser.add_argument('--operation',choices=['multiply','cross','rotation2d','rotation3d'],required=True)
    parser.add_argument('--size',type=int,required=True)
    parser.add_argument('--no-build',action='store_true')
    args=parser.parse_args()
    result=add_profiles(args.report,args.operation,args.size,args.no_build)
    return 0 if all(p['status']=='available' for p in result) else 1

if __name__=='__main__':raise SystemExit(main())
