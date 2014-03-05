#!/usr/bin/bash
#
# BHextractor HTML page
#
# Generates an html page to contain beta posterior plots for a given injection

# --- Define some environment variables to contain appropriately named
# directories etc...

injection_name="RO311" # Injection catalogue and waveform number
catalogue_name="RO3" # Reconstruction catalogue
filename=${injection_name}-${catalogue_name}
current_directory=pwd
logB=
SNR=10
numactive=50
nits=10

# HTML string:
echo """
<html>

<h1>Injection: ${injection_name}, Catalogue: ${catalogue_name}</h1>

<h1>Input Information</h1>
<ul>
<li>SNR = ${SNR}</li>
<li>Live Points = ${numactive}</li>
<li>MCMC Samples = ${nits}</li>
</ul>
""" > ${filename}.html

for seed in {1..20}
do

    echo """
    <html>

    <h2>Run #${seed}: Analysis Information</h2>
    <ul>
    <li>logB = ${logB}</li>
    <li>Seed = ${seed}</li>
    </ul>
    
    <h2>Beta Posterior</h2>
    <table>
    <tr>
    <td><img src="${filename}-seed${seed}_postbetas.png"></td>
    </tr>
    </table>
    """ >> ${filename}.html

   done

   # Call python plotting script to make figures:
    python /data/nmangini/SMEE_repo/Utilities/plot_betas.py ${filename}-seed

