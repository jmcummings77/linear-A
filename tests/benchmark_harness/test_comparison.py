import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'benchmarks'))
import compare

class ComparisonPublicationTests(unittest.TestCase):
    def sample(self):
        return json.loads((ROOT/'benchmarks/reports/vectors/results.json').read_text())
    def test_renderer_preserves_samples_and_drops_private_payloads(self):
        data=self.sample();data['environment']={'API_KEY':'secret-value'};data['profiles']=[{'private':'secret-value'}]
        data['implementations'][0]['runner']=['/Users/name/private/runner'];data['results'][0]['error']='secret-value'
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'index.html';compare.render([('one',data),('two',data)],p);html=p.read_text()
        encoded=html.split('<script id="snapshots" type="application/json">')[1].split('</script>')[0]
        saved=json.loads(encoded)[0]['data'];self.assertNotIn('secret-value',html);self.assertNotIn('/Users/',html)
        self.assertEqual(saved['results'][0]['samples'][0]['elapsed_ns'],data['results'][0]['samples'][0]['elapsed_ns'])
        self.assertEqual(saved['results'][0]['samples'][0]['ns_per_op'],data['results'][0]['samples'][0]['ns_per_op'])
    def test_embedded_labels_cannot_terminate_script(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'index.html';compare.render([('</script>',self.sample()),('two',self.sample())],p)
            encoded=p.read_text().split('<script id="snapshots" type="application/json">')[1].split('</script>')[0]
            self.assertEqual(json.loads(encoded)[0]['label'],'</script>')
    def test_git_refs_resolve_to_commits_before_reading_files_without_checkout(self):
        sha='a'*40
        with patch.object(compare.subprocess,'check_output',side_effect=[sha+'\n',json.dumps(self.sample())]) as command:
            data,resolved=compare.git_snapshot('HEAD~2','benchmarks/reports/latest/results.json')
        self.assertEqual(resolved,sha);self.assertIn('--end-of-options',command.call_args_list[0].args[0])
        self.assertEqual(command.call_args_list[1].args[0],['git','show',sha+':benchmarks/reports/latest/results.json'])
        for path in ['/tmp/private.json','../private.json','a\\b.json','README.md']:
            with self.assertRaises(ValueError):compare.git_snapshot('main',path)
    def test_wrong_schema_and_missing_report_fields_fail_explicitly(self):
        for data in ({},[],{'schema_version':2,'results':[],'implementations':[]}):
            with self.assertRaises(ValueError):compare.select_public(data)
