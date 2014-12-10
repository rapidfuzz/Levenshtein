#!/bin/sh

if test -f genextdoc.py
then
    python genextdoc.py Levenshtein
else
    echo "Could not found genextdoc.py; run gendoc.sh in the docs directory"
fi
