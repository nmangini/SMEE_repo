#!/bin/bash

for cat in Q HR RO3
do
    for theta in 0 90
    do
        ./bhextractor_pca.py ${cat} ${theta}
    done
done
