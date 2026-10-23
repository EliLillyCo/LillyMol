#!/usr/bin/env python3
"""HTTP service for GFP near-neighbour searching.

The service wraps the `lillymol_gfp_server.GfpNearNeighbourServer` extension class.
It supports both binary protobuf endpoints and protobuf-JSON endpoints.

Typical use:

  PYTHONPATH=bazel-bin/nanobind:bazel-bin \
    ~/.venv/bin/python Utilities/GFP_Tools/gfp_http_server.py \
      --gfp rand10.gfp --port 8000

Binary protobuf endpoints:
  POST /search  NnRequest -> NnReply
  POST /batch   NnBatchRequest -> NnBatchReply

JSON endpoints use protobuf JSON field names:
  POST /search.json
  POST /batch.json
"""

from __future__ import annotations

import argparse
from typing import Any

from fastapi import FastAPI, HTTPException, Request, Response
from starlette.concurrency import run_in_threadpool
from google.protobuf.json_format import MessageToDict, ParseDict, ParseError
from google.protobuf.message import DecodeError
import uvicorn

import lillymol_gfp_server
from Utilities.GFP_Tools import nn_request_pb2

_PROTOBUF_MEDIA_TYPE = "application/x-protobuf"


def _message_to_dict(message: Any) -> dict[str, Any]:
  return MessageToDict(
      message,
      always_print_fields_with_no_presence=True,
      preserving_proto_field_name=True,
  )


def create_app(
    gfp_file: str,
    pool_size_hint: int = 0,
    store_smiles: bool = True,
) -> FastAPI:
  """Create a FastAPI app backed by an in-memory GFP search pool."""

  engine = lillymol_gfp_server.GfpNearNeighbourServer(
      gfp_file,
      pool_size_hint=pool_size_hint,
      store_smiles=store_smiles,
  )

  app = FastAPI(
      title="LillyMol GFP near-neighbour server",
      version="0.1",
  )
  app.state.engine = engine
  app.state.gfp_file = gfp_file

  @app.get("/healthz")
  def healthz() -> dict[str, Any]:
    return {
        "status": "ok",
        "pool_size": app.state.engine.pool_size(),
        "gfp_file": app.state.gfp_file,
    }

  @app.post("/search")
  async def search_binary(request: Request) -> Response:
    body = await request.body()
    reply = await run_in_threadpool(app.state.engine.search_proto, body)
    return Response(content=reply, media_type=_PROTOBUF_MEDIA_TYPE)

  @app.post("/batch")
  async def batch_binary(request: Request) -> Response:
    body = await request.body()
    reply = await run_in_threadpool(app.state.engine.search_batch_proto, body)
    return Response(content=reply, media_type=_PROTOBUF_MEDIA_TYPE)

  @app.post("/search.json")
  async def search_json(request: Request) -> dict[str, Any]:
    try:
      payload = await request.json()
      proto = ParseDict(payload, nn_request_pb2.NnRequest())
    except (DecodeError, ParseError, ValueError, TypeError) as err:
      raise HTTPException(status_code=400, detail=str(err)) from err

    serialized_reply = await run_in_threadpool(
        app.state.engine.search_proto, proto.SerializeToString())
    reply = nn_request_pb2.NnReply()
    reply.ParseFromString(serialized_reply)
    return _message_to_dict(reply)

  @app.post("/batch.json")
  async def batch_json(request: Request) -> dict[str, Any]:
    try:
      payload = await request.json()
      proto = ParseDict(payload, nn_request_pb2.NnBatchRequest())
    except (DecodeError, ParseError, ValueError, TypeError) as err:
      raise HTTPException(status_code=400, detail=str(err)) from err

    serialized_reply = await run_in_threadpool(
        app.state.engine.search_batch_proto, proto.SerializeToString())
    reply = nn_request_pb2.NnBatchReply()
    reply.ParseFromString(serialized_reply)
    return _message_to_dict(reply)

  return app


def _build_parser() -> argparse.ArgumentParser:
  parser = argparse.ArgumentParser(
      description="Serve LillyMol GFP near-neighbour searches over HTTP",
  )
  parser.add_argument(
      "--gfp",
      required=True,
      help="Input GFP file loaded into memory at startup",
  )
  parser.add_argument(
      "--host",
      default="127.0.0.1",
      help="Host/interface to bind [default: 127.0.0.1]",
  )
  parser.add_argument(
      "--port",
      type=int,
      default=8000,
      help="Port to listen on [default: 8000]",
  )
  parser.add_argument(
      "--pool-size-hint",
      type=int,
      default=0,
      help="Optional fingerprint count hint; if omitted the file is scanned",
  )
  parser.add_argument(
      "--no-store-smiles",
      action="store_true",
      help="Do not retain pool SMILES in memory or return neighbour smiles",
  )
  return parser


def main() -> None:
  parser = _build_parser()
  args = parser.parse_args()

  app = create_app(
      args.gfp,
      pool_size_hint=args.pool_size_hint,
      store_smiles=not args.no_store_smiles,
  )
  uvicorn.run(
      app,
      host=args.host,
      port=args.port,
  )


if __name__ == "__main__":
  main()
