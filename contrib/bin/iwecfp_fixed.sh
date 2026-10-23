#!/usr/bin/env bash
# Invoke iwecfp for descriptor generation
exec iwecfp -R 3 -P C -Y desc=1024 "$@"
