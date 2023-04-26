#!/bin/bash

tail -n 1 *.out
find . -name \*.out -type f -delete
find . -type f -name "*.glo" -print0 | xargs -0 sed -i 's/given/g/g'
find . -name '*.csv' | xargs wc -l # number of lines
rsync -avm --include='*.sh' -f 'hide,! */' continuous/* lcontinuous/
