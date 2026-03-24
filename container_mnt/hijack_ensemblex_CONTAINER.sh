#!/bin/bash

# Loop through arguments until we find the one that matches CONTAINER1
while [[ $# -gt 0 ]]; do
    arg="$1"
    shift

    if [[ $arg == $CONTAINER1 ]]; then
        break
    fi

done

# Now execute
real_cmd="$@"
eval "$real_cmd"
exit $?
