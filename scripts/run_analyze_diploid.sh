#!/bin/sh

echo analyze diploid maps

for d in diploid*; do
    cd $d
    pwd
    ../analyze_diploid_maps.sh > analyze_diploid_maps.out
    cd ..
done

echo done
