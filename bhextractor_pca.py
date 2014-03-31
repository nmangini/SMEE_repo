#!/usr/bin/env python
# Copyright (C) 2014-2015 James Clark <james.clark@ligo.org>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""
"""

from __future__ import division

from os import environ,listdir
from os.path import isdir, isfile, join
import sys 

__author__ = "James Clark <james.clark@ligo.org>"

import numpy as np
import scipy.signal as signal
import scipy.io as sio

import lal
import lalsimulation as lalsim
import pmns_utils

#import matplotlib
#from matplotlib import pyplot as pl
#from matplotlib.mlab import PCA

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INPUT 

catalogue_name=sys.argv[1]
theta=float(sys.argv[2])
phi=0.0

# Quick 'n' dirty sanity check (should use option parser really)
if catalogue_name not in ['Q','HR','RO3']:
    print >> sys.stderr, \
            "Error, catalogue '%s' not recognised "\
            "(must be Q, HR or RO3)"%catalogue_name

Mtot=250.0
Dist=1 # Gpc
fs=16384
targetlen=1.0*fs # length of waveforms in catalogue
NRD_sampls=500 # number of samples to retain after peak

# Dictionary with the maximum waveform data lengths
maxlengths={'Q':10959, 'HR':31238, 'RO3':16493}

# Dictionary of NR sample rates. XXX Should probably get this from the data,
# really
NR_deltaT={'Q':0.155, 'HR':0.08, 'RO3':2./15}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Construct Waveforms

# --- Identify waveforms in catalogue
catalogue_path=environ['SMEEBBH_PREFIX']+'/Waveforms/'+catalogue_name+'-series'
waveforms = [ f for f in listdir(catalogue_path) if isdir(join(catalogue_path,f)) ]

# --- Load waveforms

# Preallocate a matrix for the catalogue
wflen=maxlengths[catalogue_name]
catalogue=np.zeros(shape=(wflen,len(waveforms)))
times=np.zeros(shape=(wflen,len(waveforms)))

# Get waveforms in this catalogue
for w,waveform in enumerate(waveforms):
    print 'Building %s waveform'%waveform

    # Get files (modes) in this waveform dir
    waveform_dir = catalogue_path+'/'+waveform
    modes = [ f for f in listdir(waveform_dir) if isfile(join(waveform_dir,f)) ]

    # --- Sum modes
    # See line 95 of NRWaveInject.c
    hplus=np.zeros(wflen)
    hcross=np.zeros(wflen)
    for mode in modes:

        # Load mode data
        mode_data = np.loadtxt(waveform_dir+'/'+mode)

        # Extract l,m from filename
        spharm_degree = int(mode.split('_')[2].replace('l',''))
        spharm_order  = int(mode.split('_')[3].replace('m',''))

        # Compute spin -2 weighted spherical harmonic
        sYlm = lal.SpinWeightedSphericalHarmonic(theta, phi, 
                -2, spharm_degree, spharm_order)

        # Orient Waveforms
        hplus[0:len(mode_data)]  += mode_data[:,1]*np.real(sYlm) + mode_data[:,2]*np.imag(sYlm)
        hcross[0:len(mode_data)] += mode_data[:,2]*np.real(sYlm) - mode_data[:,1]*np.imag(sYlm)

    # XXX: for now we only consider optimally oriented signals so we'll only add
    # the plus polarisation to the catalogue
    catalogue[:,w] = hplus
    

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Align Peaks

print 'aligning peak times'

# Time axis
time=mode_data[:,0]

# Find peak indices
peak_indices=np.argmax(abs(catalogue),axis=0)

# Align all waveform peaks to the latest-time peak
align_to_idx=max(peak_indices)
#align_to_idx=np.floor(wflen/2)

for w in range(len(waveforms)):

    # Temp array to hold aligned waveform
    tmp = np.zeros(len(catalogue[:,w]))
    wf  = np.copy(catalogue[:,w])

    # Get the lengths of current waveform data to the left/right of the peak of
    # this waveform
    llen = len(wf[:peak_indices[w]])
    rlen = len(tmp[align_to_idx:])

    # populate left side of peak
    tmp[align_to_idx-llen:align_to_idx] = wf[:peak_indices[w]]
    # populate right side of peak
    tmp[align_to_idx:] = wf[peak_indices[w]:peak_indices[w]+rlen]

    # Re-populate catalogue with the aligned waveform
    catalogue[:,w] = np.copy(tmp)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Truncate (and taper) Inspiral
# 
# Let's retain 1.0 seconds of a 250 Msun inspiral, walking backwards from the
# peak time
#
# This is also a convenient place to apply the mass scaling and resample the PCs
# to 16384 Hz.  XXX: We probably want to do this elsewhere (e.g., in the
# evidence calculation itself) in future.

print 'truncating and tapering'

inspiral_time = 0.5
inspiral_len = np.ceil((inspiral_time / (Mtot*lal.LAL_MTSUN_SI)) /
        np.diff(time)[0])
turn_on_idx = align_to_idx - inspiral_len
catalogue[0:turn_on_idx] = 0.0

# ---------- MASS DEPENDENT CALCULATIONS -----------
# Determine current sample rate.  XXX: This assumes a mass scale
NR_deltaT = Mtot * lal.LAL_MTSUN_SI * NR_deltaT[catalogue_name]
Mscale = Mtot * lal.LAL_MRSUN_SI / (Dist * 1e9 * lal.LAL_PC_SI)

# Resampled catalogue
#resamp_catalogue = np.zeros(shape=(wflen*NR_deltaT*fs,len(waveforms)))
resamp_catalogue = np.zeros(shape=(wflen*NR_deltaT*fs,len(waveforms)))

# Resample and taper
for w in range(len(waveforms)):

        # --- Get rid of late-time numerical noise
        #zero_tail_idx = \
        #        np.argwhere(abs(catalogue[:,w])>0.01*max(abs(catalogue[:,w])))[-1]
        zero_tail_idx = np.argmax(abs(catalogue[:,w]))+NRD_sampls
        catalogue[zero_tail_idx:,w] = 0.0

        # --- Tapering is done on lal TimeSeries objects
        hoft = lal.CreateREAL8TimeSeries('hoft', lal.LIGOTimeGPS(), 0.0,
                NR_deltaT, lal.lalStrainUnit, maxlengths[catalogue_name])
        hoft.data.data = catalogue[:,w]
        lalsim.SimInspiralREAL8WaveTaper(hoft.data,
                lalsim.LAL_SIM_INSPIRAL_TAPER_STARTEND)

        # --- Resampling (Note that LAL only handles downsampling by a factor of 2)
        # XXX apply mass scaling here.
        resampled_wf = Mscale*signal.resample(hoft.data.data, wflen*NR_deltaT * fs)

        # --- Compute SNR and rescale to SNR=1
        hoft = lal.CreateREAL8TimeSeries('hoft', lal.LIGOTimeGPS(), 0.0,
                1/float(fs), lal.lalStrainUnit, len(resampled_wf))
        hoft.data.data = resampled_wf
        rho, hrss, _, _, _ = pmns_utils.optimal_snr(hoft,freqmin=10,freqmax=100)
        resampled_wf *= 1.0/rho
        #resampled_wf *= 1.0/hrss

        # Populate catalogue with scaled, resampled, tapered data. 
        resamp_catalogue[:,w] = resampled_wf

# Eliminate memory-intensive and unncessary data
#del catalogue, hoft, hplus, hcross

# Finally, remove leading/trailing zeros
# We'll just remove everything that's <0.1% of the peak amplitude.  One of the
# earlier steps, probably tapering, means things are formally non-zero
lead_nonzero_idx = np.zeros(len(waveforms))
for w in range(len(waveforms)):
    lead_nonzero_idx[w] = \
            np.argwhere(abs(resamp_catalogue[:,w])>0.01*max(abs(resamp_catalogue[:,w])))[0]
resamp_catalogue = resamp_catalogue[min(lead_nonzero_idx):,:]

# XXX This doesn't work so well for the end
trail_nonzero_idx = np.zeros(len(waveforms))
for w in range(len(waveforms)):
    trail_nonzero_idx[w] = \
            np.argwhere(abs(resamp_catalogue[:,w])>0.001*max(abs(resamp_catalogue[:,w])))[-1]

resamp_catalogue = resamp_catalogue[:min(trail_nonzero_idx),:]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create standardised catalogue 
# --- all waveforms have 1s duration with peak aligned to 0.75s

waveform_catalogue = np.zeros(shape=(targetlen,len(waveforms)))
peak_indices=np.argmax(resamp_catalogue,axis=0)
align_to_idx=np.floor(0.75*targetlen)
for w in range(len(waveforms)):

    # current waveform
    wf=resamp_catalogue[:,w]

    # Get the lengths of current waveform data to the left/right of the peak of
    # this waveform
    llen = len(wf[:peak_indices[w]])
    rlen = len(wf[peak_indices[w]:])

    # populate left side of peak
    waveform_catalogue[align_to_idx-llen:align_to_idx,w] = wf[:peak_indices[w]]
    # populate right side of peak
    waveform_catalogue[align_to_idx:align_to_idx+rlen,w] = wf[peak_indices[w]:peak_indices[w]+rlen]

del resamp_catalogue

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PCA
#
# This should be a straight copy of the matlab PCA scripts.  It looks like there
# are some differences in eig() between numpy and matlab but the results should
# come out the same.
#
# H = catalogue matrix (columns=waveforms)
# 
print 'Performing PCA'

# --- Make catalogue a matrix for matrix arithmetic
H = np.matrix(waveform_catalogue)

# --- 1) Compute catalogue covariance
C = H.T * H / np.shape(H)[0]

# --- 2) Compute eigenvectors (V) and eigenvalues (S) of C
S, V = np.linalg.eig(C)

# --- 3) Sort eigenvectors in descending order
idx = np.argsort(S)[::-1]
V = V[:,idx]

# --- 4) Compute the eigenvectors of the real covariance matrix U and normalise
U = H*V 
U /= np.linalg.norm(U,axis=0)

# Convert PCs back to an numpy.array
U = np.array(U)

# --- 5) Calculate the beta values by projecting A onto U
# Don't think we really need this for now
#betas = np.zeros(shape=(np.shape(catalogue)))
#recWf = np.zeros(shape=(np.shape(catalogue),len(waveforms)))

#for i in range(len(waveforms)):
#    for j in range(len(waveforms)):
#        betas[i,:] = C.T[:,i] * U


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save results

PCA_outname=catalogue_name + '_PCs_' + 'theta-%.0f'%theta
PC_dict={'PCs_final':U}
sio.savemat(PCA_outname, PC_dict)

catalogue_outname=catalogue_name + '_catalogue_' + 'theta-%.0f'%theta
waveform_dict={'MDC_final':waveform_catalogue}
sio.savemat(catalogue_outname, waveform_dict)

print "PCA complete"
print "PC matrix written to %s.mat"%PCA_outname
print "Catalogue matrix written to %s.mat"%catalogue_outname










