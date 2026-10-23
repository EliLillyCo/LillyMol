#!/usr/bin/env python3
"""Command-line client for the GFP HTTP near-neighbour server."""

from __future__ import annotations

import argparse
import contextlib
import sys
import urllib.error
import urllib.parse
import urllib.request

from google.protobuf.json_format import MessageToJson
from google.protobuf.message import DecodeError

from Utilities.GFP_Tools import nn_request_pb2

_PROTOBUF_MEDIA_TYPE = "application/x-protobuf"


def _build_parser() -> argparse.ArgumentParser:
  parser = argparse.ArgumentParser(
      description="Query a LillyMol GFP HTTP near-neighbour server",
  )
  parser.add_argument(
      "--url",
      default="http://127.0.0.1:8000",
      help="Base server URL [default: http://127.0.0.1:8000]",
  )
  parser.add_argument(
      "--smiles",
      default="C1=NC(=CC2=C1C=CC=C2)NC(=O)NCC",
      help="Query SMILES",
  )
  parser.add_argument(
      "--id",
      default="CHEMBL4076111",
      help="Query identifier",
  )
  parser.add_argument(
      "--input",
      help="Input SMILES file. Column 1 is SMILES, column 2 is identifier. Use '-' for stdin",
  )
  parser.add_argument(
      "--batch-size",
      type=int,
      default=1,
      help="Number of input records to send per /batch request [default: 1]",
  )
  group = parser.add_mutually_exclusive_group(required=False)
  group.add_argument(
      "--nbrs",
      type=int,
      default=1,
      help="Number of neighbours to find [default: 1]",
  )
  group.add_argument(
      "--distance",
      type=float,
      help="Return all neighbours within this distance",
  )
  parser.add_argument(
      "--niter",
      type=int,
      default=1,
      help="Number of repeated lookups, useful for simple timing [default: 1]",
  )
  parser.add_argument(
      "--output-format",
      choices=["textproto", "json", "tdt"],
      default="textproto",
      help="Output format [default: textproto]",
  )
  parser.add_argument(
      "--json",
      action="store_true",
      help="Alias for --output-format json",
  )
  parser.add_argument(
      "--timeout",
      type=float,
      default=60.0,
      help="HTTP timeout in seconds [default: 60]",
  )
  return parser


def _search_url(base_url: str) -> str:
  return urllib.parse.urljoin(base_url.rstrip('/') + '/', 'search')


def _batch_url(base_url: str) -> str:
  return urllib.parse.urljoin(base_url.rstrip('/') + '/', 'batch')


def _build_request(args: argparse.Namespace, smiles: str | None = None,
                   identifier: str | None = None) -> nn_request_pb2.NnRequest:
  request = nn_request_pb2.NnRequest()
  request.smiles = smiles if smiles is not None else args.smiles
  request.id = identifier if identifier is not None else args.id

  if args.distance is not None:
    request.distance = args.distance
  else:
    request.nbrs = args.nbrs

  return request


def _post_proto(url: str, request: nn_request_pb2.NnRequest, timeout: float) -> nn_request_pb2.NnReply:
  http_request = urllib.request.Request(
      url,
      data=request.SerializeToString(),
      headers={
          'Content-Type': _PROTOBUF_MEDIA_TYPE,
          'Accept': _PROTOBUF_MEDIA_TYPE,
      },
      method='POST',
  )

  try:
    with urllib.request.urlopen(http_request, timeout=timeout) as response:
      payload = response.read()
  except urllib.error.HTTPError as err:
    detail = err.read().decode('utf-8', errors='replace')
    raise RuntimeError(f'HTTP {err.code}: {detail}') from err
  except urllib.error.URLError as err:
    raise RuntimeError(f'Cannot contact server {url}: {err.reason}') from err

  reply = nn_request_pb2.NnReply()
  try:
    reply.ParseFromString(payload)
  except DecodeError as err:
    raise RuntimeError('Server returned invalid NnReply protobuf') from err

  return reply


def _post_batch_proto(url: str, request: nn_request_pb2.NnBatchRequest,
                      timeout: float) -> nn_request_pb2.NnBatchReply:
  http_request = urllib.request.Request(
      url,
      data=request.SerializeToString(),
      headers={
          'Content-Type': _PROTOBUF_MEDIA_TYPE,
          'Accept': _PROTOBUF_MEDIA_TYPE,
      },
      method='POST',
  )

  try:
    with urllib.request.urlopen(http_request, timeout=timeout) as response:
      payload = response.read()
  except urllib.error.HTTPError as err:
    detail = err.read().decode('utf-8', errors='replace')
    raise RuntimeError(f'HTTP {err.code}: {detail}') from err
  except urllib.error.URLError as err:
    raise RuntimeError(f'Cannot contact server {url}: {err.reason}') from err

  reply = nn_request_pb2.NnBatchReply()
  try:
    reply.ParseFromString(payload)
  except DecodeError as err:
    raise RuntimeError('Server returned invalid NnBatchReply protobuf') from err

  return reply


def _output_format(args: argparse.Namespace) -> str:
  if args.json:
    return 'json'
  return args.output_format


def _status_name(status: int) -> str:
  return nn_request_pb2.Status.Name(status)


def _write_tdt(reply: nn_request_pb2.NnReply,
               request: nn_request_pb2.NnRequest | None) -> None:
  if reply.status != nn_request_pb2.OK:
    message = reply.message or _status_name(reply.status)
    if request is not None and request.id:
      print(f'Warning: {request.id}: {message}', file=sys.stderr)
    else:
      print(f'Warning: {message}', file=sys.stderr)
    return

  result = reply.result
  smiles = result.smiles
  identifier = result.name
  if request is not None:
    if not smiles:
      smiles = request.smiles
    if not identifier:
      identifier = request.id

  print(f'$SMI<{smiles}>')
  print(f'PCN<{identifier}>')
  for nbr in result.nbr:
    if nbr.smi:
      print(f'$SMI<{nbr.smi}>')
    print(f'PCN<{nbr.id}>')
    print(f'DIST<{nbr.dist:g}>')
  print('|')


def _write_reply(reply: nn_request_pb2.NnReply, output_format: str,
                 request: nn_request_pb2.NnRequest | None = None) -> None:
  if output_format == 'json':
    print(MessageToJson(
        reply,
        always_print_fields_with_no_presence=True,
        preserving_proto_field_name=True,
    ))
  elif output_format == 'tdt':
    _write_tdt(reply, request)
  else:
    print(reply)


def _read_smiles_file(fname: str):
  if fname == '-':
    input_stream = contextlib.nullcontext(sys.stdin)
  else:
    input_stream = open(fname, 'r', encoding='utf-8')

  with input_stream as reader:
    for line_number, line in enumerate(reader, 1):
      line = line.strip()
      if not line:
        continue

      fields = line.split()
      if len(fields) < 2:
        raise ValueError(f'{fname}:{line_number}: expected SMILES and identifier')

      yield fields[0], fields[1]


def _read_batches(fname: str, batch_size: int):
  batch = []
  for smiles, identifier in _read_smiles_file(fname):
    batch.append((smiles, identifier))
    if len(batch) == batch_size:
      yield batch
      batch = []

  if batch:
    yield batch


def _build_batch_request(args: argparse.Namespace,
                         records: list[tuple[str, str]]) -> nn_request_pb2.NnBatchRequest:
  batch = nn_request_pb2.NnBatchRequest()
  for smiles, identifier in records:
    batch.request.append(_build_request(args, smiles, identifier))
  return batch


def main(argv: list[str] | None = None) -> int:
  parser = _build_parser()
  args = parser.parse_args(argv)

  if args.niter < 1:
    parser.error('--niter must be positive')
  if args.nbrs is not None and args.nbrs < 1:
    parser.error('--nbrs must be positive')
  if args.distance is not None and args.distance < 0.0:
    parser.error('--distance must be non-negative')
  if args.batch_size < 1:
    parser.error('--batch-size must be positive')

  output_format = _output_format(args)

  try:
    if args.input:
      if args.niter != 1:
        parser.error('--niter is only supported for single-molecule queries')
      if args.batch_size == 1:
        url = _search_url(args.url)
        for smiles, identifier in _read_smiles_file(args.input):
          request = _build_request(args, smiles, identifier)
          reply = _post_proto(url, request, args.timeout)
          _write_reply(reply, output_format, request)
      else:
        url = _batch_url(args.url)
        for records in _read_batches(args.input, args.batch_size):
          request = _build_batch_request(args, records)
          reply = _post_batch_proto(url, request, args.timeout)
          if len(reply.reply) != len(records):
            raise RuntimeError(
                f'Server returned {len(reply.reply)} replies for {len(records)} requests')
          for item, original_request in zip(reply.reply, request.request):
            _write_reply(item, output_format, original_request)
    else:
      url = _search_url(args.url)
      request = _build_request(args)
      for _ in range(args.niter):
        reply = _post_proto(url, request, args.timeout)
        _write_reply(reply, output_format, request)
  except (RuntimeError, ValueError) as err:
    print(err, file=sys.stderr)
    return 1

  return 0


if __name__ == '__main__':
  sys.exit(main())
