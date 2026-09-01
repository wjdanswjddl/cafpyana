#!/usr/bin/env bash
# Wrapper — canonical script lives on PRL data volume.
exec bash /exp/sbnd/data/users/munjung/PRL_data/status "$@"
