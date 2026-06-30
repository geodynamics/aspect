#!/bin/bash

if [ "$1" = "screen-output" ]; then
  grep -v "btl_tcp_component_create_listen] bind() failed: Operation not permitted"
else
  cat
fi
