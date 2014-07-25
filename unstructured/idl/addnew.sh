#!/bin/sh

for i in `find . -maxdepth 1 -mtime -1`; do
    svn add $i
done
