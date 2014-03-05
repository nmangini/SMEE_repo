#!/usr/bin/bash
#
# BHextractor HTML page
#
# Generates an html page to contain beta posterior plots for a given injection

# --- Define some environment variables to contain appropriately named
# directories etc...

injection_name="HR13" # Injection catalogue and waveform number
injection="HR"
filename=${injection_name}-resultspage
current_directory=$(pwd)
logB=
SNR=10
numactive=50
nits=10

# HTML string:
echo """
<html>

<h1>Injection: ${injection_name}</h1>

<h2>Input Information</h2>
<ul>
<li>SNR = ${SNR}</li>
<li>Live Points = ${numactive}</li>
<li>MCMC Samples = ${nits}</li>
</ul>

<h2>logB vs. Seed #</h2>
<table>
<tr>
<td><img src="${current_directory}/${injection}-bayes.png"></td>
</tr>
</table>
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
    
    <h2>Beta Posteriors</h2>
    <table>
    <tr>
    <th><b>Q Reconstruction</b></th>
    <th><b>HR Reconstruction</b></th>
    <th><b>RO3 Reconstruction</b></th>
    </tr>
    <tr>
    <td><img src="${current_directory}/${injection_name}-Q-seed${seed}_postbetas.png"></td>
    <td><img src="${current_directory}/${injection_name}-HR-seed${seed}_postbetas.png"></td>
    <td><img src="${current_directory}/${injection_name}-RO3-seed${seed}_postbetas.png"></td>
    </tr>
    </table>
    """ >> ${filename}.html

   done

   # Call python plotting script to make figures:
#    python /data/nmangini/SMEE_repo/Utilities/plot_betas.py ${filename}-seed

