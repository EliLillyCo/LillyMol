#!/bin/bash

# $Id$

if [[ -d "../ruby" ]]
then # place here all your code
    # define executable
    program="../ruby/dopattern.rb"
    # check that executable exists
    if [ ! -s "$program" ]
        then
        echo "Cannot access executable '$program'" >&2
        exit 1
    fi

    # this will execute and exit, returning the exit status of the command
    exec $program "$@"

else
    echo $(basename $0)": Required libraries are not found under ruby folder" >&2 && exit 1
fi
