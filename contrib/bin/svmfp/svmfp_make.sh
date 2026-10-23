#!/usr/bin/env bash
ruby_script="${0%%.sh}.rb"
exec ruby ${ruby_script} "$@"
