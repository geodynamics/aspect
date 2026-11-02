#!/bin/sh
# The plugin variant in the loading message is not a numerical result.
"$(dirname "$0")/../cmake/default" "$@" | sed -E 's/\.(debug|release)\.so>/.so>/g'
