import contextlib
import io
import os
import sys
import tempfile
import unittest


def _add_runfiles_to_path() -> None:
  candidates = ['.']
  for envvar in ('TEST_SRCDIR', 'RUNFILES_DIR'):
    root = os.environ.get(envvar)
    if root:
      candidates.extend([
          os.path.join(root, '_main'),
          root,
      ])
  for candidate in candidates:
    if os.path.isdir(candidate) and candidate not in sys.path:
      sys.path.insert(0, candidate)


_add_runfiles_to_path()

from Utilities.GFP_Tools import gfp_client_app
from Utilities.GFP_Tools import nn_request_pb2


def _temp_smiles_file(contents: str) -> str:
  handle = tempfile.NamedTemporaryFile(mode='w', suffix='.smi', delete=False)
  with handle:
    handle.write(contents)
  return handle.name


class GfpHttpClientTest(unittest.TestCase):

  def test_read_smiles_file(self):
    fname = _temp_smiles_file('CCO ethanol\n\nCCN ethylamine extra\n')
    self.addCleanup(lambda: os.path.exists(fname) and os.unlink(fname))

    self.assertEqual(
        list(gfp_client_app._read_smiles_file(fname)),
        [('CCO', 'ethanol'), ('CCN', 'ethylamine')],
    )

  def test_read_smiles_file_rejects_missing_identifier(self):
    fname = _temp_smiles_file('CCO\n')
    self.addCleanup(lambda: os.path.exists(fname) and os.unlink(fname))

    with self.assertRaisesRegex(ValueError, 'expected SMILES and identifier'):
      list(gfp_client_app._read_smiles_file(fname))

  def test_read_batches(self):
    fname = _temp_smiles_file('CCO ethanol\nCCN ethylamine\nC methane\n')
    self.addCleanup(lambda: os.path.exists(fname) and os.unlink(fname))

    self.assertEqual(
        list(gfp_client_app._read_batches(fname, 2)),
        [[('CCO', 'ethanol'), ('CCN', 'ethylamine')], [('C', 'methane')]],
    )

  def test_build_batch_request(self):
    args = gfp_client_app._build_parser().parse_args(['--nbrs', '4'])
    request = gfp_client_app._build_batch_request(
        args, [('CCO', 'ethanol'), ('CCN', 'ethylamine')])

    self.assertEqual(len(request.request), 2)
    self.assertEqual(request.request[0].smiles, 'CCO')
    self.assertEqual(request.request[0].id, 'ethanol')
    self.assertEqual(request.request[0].nbrs, 4)
    self.assertEqual(request.request[1].smiles, 'CCN')

  def test_json_alias(self):
    args = gfp_client_app._build_parser().parse_args(['--output-format', 'tdt', '--json'])
    self.assertEqual(gfp_client_app._output_format(args), 'json')

  def test_tdt_output_with_neighbour(self):
    request = nn_request_pb2.NnRequest(smiles='CCO', id='ethanol', nbrs=2)
    reply = nn_request_pb2.NnReply(status=nn_request_pb2.OK)
    reply.result.smiles = 'CCO'
    reply.result.name = 'ethanol'
    nbr = reply.result.nbr.add()
    nbr.smi = 'CCN'
    nbr.id = 'ethylamine'
    nbr.dist = 0.25

    stdout = io.StringIO()
    with contextlib.redirect_stdout(stdout):
      gfp_client_app._write_reply(reply, 'tdt', request)

    self.assertEqual(
        stdout.getvalue(),
        '$SMI<CCO>\nPCN<ethanol>\n$SMI<CCN>\nPCN<ethylamine>\nDIST<0.25>\n|\n',
    )

  def test_tdt_zero_neighbour_output(self):
    request = nn_request_pb2.NnRequest(smiles='CCO', id='ethanol', nbrs=2)
    reply = nn_request_pb2.NnReply(status=nn_request_pb2.OK)

    stdout = io.StringIO()
    with contextlib.redirect_stdout(stdout):
      gfp_client_app._write_reply(reply, 'tdt', request)

    self.assertEqual(stdout.getvalue(), '$SMI<CCO>\nPCN<ethanol>\n|\n')

  def test_tdt_error_warns_and_suppresses_output(self):
    request = nn_request_pb2.NnRequest(smiles='bad', id='badmol', nbrs=1)
    reply = nn_request_pb2.NnReply(status=nn_request_pb2.BAD_SMILES, message='bad smiles')

    stdout = io.StringIO()
    stderr = io.StringIO()
    with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
      gfp_client_app._write_reply(reply, 'tdt', request)

    self.assertEqual(stdout.getvalue(), '')
    self.assertEqual(stderr.getvalue(), 'Warning: badmol: bad smiles\n')

  def test_main_uses_batch_endpoint(self):
    fname = _temp_smiles_file('CCO ethanol\nCCN ethylamine\nC methane\n')
    self.addCleanup(lambda: os.path.exists(fname) and os.unlink(fname))

    seen = []
    written = []
    original_post_batch = gfp_client_app._post_batch_proto
    original_write_reply = gfp_client_app._write_reply

    def fake_post_batch(url, request, timeout):
      seen.append((url, [(r.smiles, r.id, r.nbrs) for r in request.request], timeout))
      reply = nn_request_pb2.NnBatchReply()
      for _ in request.request:
        reply.reply.add(status=nn_request_pb2.OK)
      return reply

    def fake_write_reply(reply, output_format, request=None):
      written.append((output_format, request.id if request is not None else None))

    gfp_client_app._post_batch_proto = fake_post_batch
    gfp_client_app._write_reply = fake_write_reply
    try:
      rc = gfp_client_app.main([
          '--url', 'http://example.test',
          '--input', fname,
          '--batch-size', '2',
          '--nbrs', '4',
          '--timeout', '7',
          '--output-format', 'tdt',
      ])
    finally:
      gfp_client_app._post_batch_proto = original_post_batch
      gfp_client_app._write_reply = original_write_reply

    self.assertEqual(rc, 0)
    self.assertEqual(seen, [
        ('http://example.test/batch', [('CCO', 'ethanol', 4), ('CCN', 'ethylamine', 4)], 7.0),
        ('http://example.test/batch', [('C', 'methane', 4)], 7.0),
    ])
    self.assertEqual(written, [('tdt', 'ethanol'), ('tdt', 'ethylamine'), ('tdt', 'methane')])


if __name__ == '__main__':
  unittest.main()
