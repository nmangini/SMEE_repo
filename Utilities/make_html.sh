#!/usr/bin/bash
#
# BHextractor HTML page
#
# Generates an html page to contain beta posterior plots for a given injection

# --- Define some environment variables to contain appropriately named
# directories etc...

injection_name="RO311" # Injection catalogue and waveform number
injection="RO3"
filename=${injection_name}-resultspage
image_directory="/data/nmangini/SMEE_repo/SMEE/SMEE_BBH/snr-10-increasedprior"
logB=
SNR=10
numactive=50
nits=10

# HTML string:
echo """
<html>

<h1>Injection: ${injection_name}</h1>

<h2>Input Information</h2>
<table border="1">
<tr align="center">
<td></td><b>Seed</b><td><b>MCMC Samples</b></td><td><b>Live Points</b></td><td><b>SNR</b></td><td><b>logB</b></td>
</tr>
<tr align="center">
<td>1</td><td>${nits}</td><td>${numactive}</td><td>${SNR}</td><td>${logB}</td>
</tr>
<tr align="center">
<td>2</td><td>${nits}</td><td>${numactive}</td><td>${SNR}</td><td>${logB}</td>
</tr>
<tr align="center">
<td>3</td><td>${nits}</td><td>${numactive}</td><td>${SNR}</td><td>${logB}</td>
</tr>
<tr align="center">
<td>4</td><td>${nits}</td><td>${numactive}</td><td>${SNR}</td><td>${logB}</td>
</tr>
<tr align="center">
<td>5</td><td>${nits}</td><td>${numactive}</td><td>${SNR}</td><td>${logB}</td>
</tr>
<tr align="center">
<td>6</td><td>${nits}</td><td>${numactive}</td><td>${SNR}</td><td>${logB}</td>
</tr>
<tr align="center">
<td>7</td><td>${nits}</td><td>${numactive}</td><td>${SNR}</td><td>${logB}</td>
</tr>
<tr align="center">
<td>8</td><td>${nits}</td><td>${numactive}</td><td>${SNR}</td><td>${logB}</td>
</tr>
<tr align="center">
<td>9</td><td>${nits}</td><td>${numactive}</td><td>${SNR}</td><td>${logB}</td>
</tr>
<tr align="center">
<td>10</td><td>${nits}</td><td>${numactive}</td><td>${SNR}</td><td>${logB}</td>
</tr>
<tr align="center">
<td>11</td><td>${nits}</td><td>${numactive}</td><td>${SNR}</td><td>${logB}</td>
</tr>
<tr align="center">
<td>12</td><td>${nits}</td><td>${numactive}</td><td>${SNR}</td><td>${logB}</td>
</tr>
<tr align="center">
<td>13</td><td>${nits}</td><td>${numactive}</td><td>${SNR}</td><td>${logB}</td>
</tr>
<tr align="center">
<td>14</td><td>${nits}</td><td>${numactive}</td><td>${SNR}</td><td>${logB}</td>
</tr>
<tr align="center">
<td>15</td><td>${nits}</td><td>${numactive}</td><td>${SNR}</td><td>${logB}</td>
</tr>
<tr align="center">
<td>16</td><td>${nits}</td><td>${numactive}</td><td>${SNR}</td><td>${logB}</td>
</tr>
<tr align="center">
<td>17</td><td>${nits}</td><td>${numactive}</td><td>${SNR}</td><td>${logB}</td>
</tr>
<tr align="center">
<td>18</td><td>${nits}</td><td>${numactive}</td><td>${SNR}</td><td>${logB}</td>
</tr>
<tr align="center">
<td>19</td><td>${nits}</td><td>${numactive}</td><td>${SNR}</td><td>${logB}</td>
</tr>
<tr align="center">
<td>20</td><td>${nits}</td><td>${numactive}</td><td>${SNR}</td><td>${logB}</td>
</tr>
</table>

<h2>logB vs. Seed #</h2>
<table border="1">
<tr>
<td><img src="${image_directory}/${injection}-bayes.png"></td>
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
    <table border="1">
    <tr>
    <th><b>Q Reconstruction</b></th>
    <th><b>HR Reconstruction</b></th>
    <th><b>RO3 Reconstruction</b></th>
    </tr>
    <tr>
    <td><img src="${image_directory}/${injection_name}-Q-seed${seed}_postbetas.png"></td>
    <td><img src="${image_directory}/${injection_name}-HR-seed${seed}_postbetas.png"></td>
    <td><img src="${image_directory}/${injection_name}-RO3-seed${seed}_postbetas.png"></td>
    </tr>
    </table>
    """ >> ${filename}.html

   done

   # Call python plotting script to make figures:
#    python /data/nmangini/SMEE_repo/Utilities/plot_betas.py ${filename}-seed

