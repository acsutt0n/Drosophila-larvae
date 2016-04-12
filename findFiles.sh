#!/bin/bash

# usage: ./findFiles dir ext outfilename

# This will write all filenames and paths (recursively) from a directory
# with a given extension and write them to a new txt file

if [ "$#" -ne 3 ]; then
  echo "Need 3 arguments: dir(1) ext(2) outfilename(#)"
else
  find $1 -type f -name "*$2" > $3
  echo "$3 written"
fi


