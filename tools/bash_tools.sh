#!/bin/bash

tail -n 1 *.out
find . -name \*.out -type f -delete
find . -type f -name "*.glo" -print0 | xargs -0 sed -i 's/given/g/g'
find . -name '*.csv' | xargs wc -l # number of lines
find . -type f -name "101.csv" -exec awk 'FNR == 1 && /Imimic/ {print FILENAME; nextfile}' {} + # find file with specific string
rsync -avm --include='*.glo' -f 'hide,! */' continuous/* lcontinuous/
