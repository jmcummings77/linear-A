import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'benchmarks'))
import profile_report
import profiles

class ProfileReportTests(unittest.TestCase):
    def test_enrichment_preserves_every_recorded_timing_and_original_provenance(self):
        data={'created_at':'old time','revision':'old revision','dirty':False,'source_sha256':'old digest','seed':17,
              'implementations':[{'id':'c','status':'passed','toolchain':'old compiler'}],
              'results':[{'operation':'cross','size':3,'median_ns':123,'samples':[{'elapsed_ns':246,'iterations':2}]}], 'profiles':[]}
        item={'id':'c','toolchain':'new compiler'}
        samples=[{'implementation':'c','status':'available','stacks':[['cross']],'weights':[3],'chronological':False,'workload':{'operation':'cross','size':3,'iterations':10,'seed':17}}]
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'results.json';path.write_text(json.dumps(data))
            with patch.object(profile_report,'build_one',return_value=item),patch.object(profile_report,'verify',return_value=True),patch.object(profile_report,'collect_profiles',return_value=copy.deepcopy(samples)) as collect,patch.object(profile_report,'fingerprint',return_value='new digest'),patch.object(profile_report.subprocess,'check_output',side_effect=['new revision\n','dirty']),patch.object(profile_report,'render'):
                profile_report.add_profiles(path,'cross',3,True)
            updated=json.loads(path.read_text());new=updated.pop('profiles');old=copy.deepcopy(data);old.pop('profiles')
            self.assertEqual(updated,old)
            self.assertEqual(new[0]['revision'],'new revision')
            self.assertEqual(new[0]['source_sha256'],'new digest')
            self.assertEqual(new[0]['toolchain'],'new compiler')
            self.assertTrue(new[0]['dirty'])
            self.assertEqual(collect.call_args.args,([item],path.parent,17,'cross',3))

    def test_incompatible_workload_does_not_change_the_report(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'results.json';original='{"results": [{"operation":"cross","size":3}]}'
            path.write_text(original)
            for operation,size in [('cross',4),('multiply',48),('missing',3)]:
                with self.assertRaises(ValueError):profile_report.add_profiles(path,operation,size)
                self.assertEqual(path.read_text(),original)

    def test_profile_calibration_uses_requested_vector_workload(self):
        with patch('run.measure',return_value={'ns_per_op':1000}) as measure:
            count,checksum=profiles._calibrate({'id':'c'},17,'cross',3)
        self.assertGreater(count,0)
        self.assertTrue(all(call.args[1:3]==('cross',3) for call in measure.call_args_list))
        from vector_reference import vector_checksum
        self.assertEqual(checksum,vector_checksum('cross',3,17))
