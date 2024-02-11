#!/bin/sh
curdir="${0%/*}"

generate_cython()
{
  python -m cython -I "$curdir" --cplus "$curdir"/"$1".pyx -o "$curdir"/"$1".cxx || exit 1
  echo "Generated $curdir/$1.cxx"
}

generate_cython levenshtein_cpp
