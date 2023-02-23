#!/bin/bash

source_folder="$HOME/code/gnr/results/" # source
destination_folder="$STORE/gnr/results/" # destination

rsync --archive --verbose --progress --exclude='*.out' $source_folder $destination_folder 

