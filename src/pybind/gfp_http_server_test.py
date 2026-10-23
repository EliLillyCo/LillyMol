import os
import sys
import unittest


def _add_pybind_runfiles_to_path() -> None:
  candidates = ['pybind']
  for envvar in ('TEST_SRCDIR', 'RUNFILES_DIR'):
    root = os.environ.get(envvar)
    if root:
      candidates.extend([
          os.path.join(root, '_main/pybind'),
          os.path.join(root, 'pybind'),
      ])
  for candidate in candidates:
    if os.path.isdir(candidate) and candidate not in sys.path:
      sys.path.insert(0, candidate)


_add_pybind_runfiles_to_path()

from Utilities.GFP_Tools import nn_request_pb2

try:
  from fastapi.testclient import TestClient
  from Utilities.GFP_Tools.gfp_http_server import create_app
  _IMPORT_ERROR = None
except ModuleNotFoundError as err:
  TestClient = None
  create_app = None
  _IMPORT_ERROR = err


def _testdata_path() -> str:
  candidates = ['pybind/testdata/rand10.gfp']
  test_srcdir = os.environ.get('TEST_SRCDIR')
  if test_srcdir:
    candidates.extend([
        os.path.join(test_srcdir, '_main/pybind/testdata/rand10.gfp'),
        os.path.join(test_srcdir, 'pybind/testdata/rand10.gfp'),
    ])
  runfiles_dir = os.environ.get('RUNFILES_DIR')
  if runfiles_dir:
    candidates.extend([
        os.path.join(runfiles_dir, '_main/pybind/testdata/rand10.gfp'),
        os.path.join(runfiles_dir, 'pybind/testdata/rand10.gfp'),
    ])

  for candidate in candidates:
    if os.path.exists(candidate):
      return candidate
  raise FileNotFoundError('Cannot find pybind/testdata/rand10.gfp')


class GfpHttpServerTest(unittest.TestCase):

  @classmethod
  def setUpClass(cls):
    super().setUpClass()
    if _IMPORT_ERROR is not None:
      raise unittest.SkipTest(
          f'gfp_http_server_test requires optional HTTP server dependencies: {_IMPORT_ERROR}')
    cls.app = create_app(_testdata_path())
    cls.client = TestClient(cls.app)
    cls.query = {
        'smiles': 'O(C1=CC=C2C(=C1)C=CC=C2)C(C(=O)O)CC',
        'id': 'query',
        'nbrs': 3,
    }

  def test_healthz(self):
    response = self.client.get('/healthz')
    self.assertEqual(response.status_code, 200)
    self.assertEqual(response.json()['status'], 'ok')
    self.assertEqual(response.json()['pool_size'], 10)

  def test_search_json(self):
    response = self.client.post('/search.json', json=self.query)
    self.assertEqual(response.status_code, 200)
    payload = response.json()
    self.assertEqual(payload['status'], 'OK')
    self.assertEqual(payload['result']['name'], 'query')
    self.assertEqual(len(payload['result']['nbr']), 3)
    self.assertEqual(payload['result']['nbr'][0]['id'], 'CHEMBL2141072')
    self.assertAlmostEqual(payload['result']['nbr'][0]['dist'], 0.0)

  def test_search_binary_proto(self):
    request = nn_request_pb2.NnRequest(**self.query)
    response = self.client.post(
        '/search',
        content=request.SerializeToString(),
        headers={'content-type': 'application/x-protobuf'},
    )
    self.assertEqual(response.status_code, 200)

    reply = nn_request_pb2.NnReply()
    reply.ParseFromString(response.content)
    self.assertEqual(reply.status, nn_request_pb2.OK)
    self.assertEqual(len(reply.result.nbr), 3)
    self.assertEqual(reply.result.nbr[0].id, 'CHEMBL2141072')
    self.assertAlmostEqual(reply.result.nbr[0].dist, 0.0)

  def test_batch_json(self):
    response = self.client.post('/batch.json', json={
        'request': [
            self.query,
            {'smiles': 'not a smiles', 'id': 'bad', 'nbrs': 1},
        ],
    })
    self.assertEqual(response.status_code, 200)
    payload = response.json()
    self.assertEqual([reply['status'] for reply in payload['reply']], ['OK', 'BAD_SMILES'])

  def test_batch_binary_proto(self):
    request = nn_request_pb2.NnBatchRequest()
    request.request.add(**self.query)
    request.request.add(smiles='not a smiles', id='bad', nbrs=1)

    response = self.client.post(
        '/batch',
        content=request.SerializeToString(),
        headers={'content-type': 'application/x-protobuf'},
    )
    self.assertEqual(response.status_code, 200)

    reply = nn_request_pb2.NnBatchReply()
    reply.ParseFromString(response.content)
    self.assertEqual(len(reply.reply), 2)
    self.assertEqual(reply.reply[0].status, nn_request_pb2.OK)
    self.assertEqual(reply.reply[1].status, nn_request_pb2.BAD_SMILES)

  def test_invalid_json(self):
    response = self.client.post('/search.json', json={'unknown_field': 1})
    self.assertEqual(response.status_code, 400)


if __name__ == '__main__':
  unittest.main()
