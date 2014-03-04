#!/usr/bin/bash
#
# BHextractor HTML page
#
# Generates an html page to contain beta posterior plots for a given injection

# --- Define some environment variables to contain appropriately named
# directories etc...

injection_name="Q9" # Injection catalogue and waveform number
catalogue_name="Q" # Reconstruction catalogue
filename=${injection_name}-${catalogue_name}
logB=
SNR=
numactive=
nits=

# HTML string:
for seed in {1..20}
do

    echo """
    <html>
    
    <h1>Injection: ${injection_name}, Catalogue: ${catalogue_name}</h1>
    
    <h2>Bayes Factor</h2>
    logB = ${logB}
    
    <h2>Beta Posteriors</h2>
    <table>
    <tr>
    <td><b>Beta 1</b></td>
    </tr>
    <tr>
    <td><img src="${filename}-seed${seed}_postbeta.png"></td>
    </tr>
    </table>
    """ > ${filename}.html

   done

   # Call python plotting script to make figures:
    python /data/nmangini/SMEE_repo/Utilities/plot_betas.py ${filename}-seed

