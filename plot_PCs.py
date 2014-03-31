#!/usr/bin/env python

import sys
import numpy as np
import scipy.io as sio
from scipy import signal
import HHT
from matplotlib import pyplot as pl

inputfile=sys.argv[1]
outname=inputfile.replace('.mat','')
data=sio.loadmat(inputfile)

fs=16384.

ZDHP=np.loadtxt('/Users/jclark/Projects/SMEEBBH/SMEE_repo/SRDs/ZERO_DET_high_P.txt')
noise_freqs=ZDHP[:,0]
Sf=ZDHP[:,1]**2

# --- HHT
fig,ax=pl.subplots(np.shape(data['PCs_final'])[1],ncols=2,figsize=(8,20),sharex='col')
for r,row in enumerate(ax):

    # Take HHT
    time=np.arange(0,len(data['PCs_final'][:,r])/fs,1/fs)
    h=HHT.HHT(time,data['PCs_final'][:,r],N=0)

    row[0].plot(time,h[0]/max(h[0]))
    row[0].set_ylim(0,1.1)

    row[1].plot(time,h[1])
    row[1].set_ylim(0,100)

    if r==0:
        row[0].set_title('Instaneous Amplitude [arb units]')
        row[1].set_title('Instaneous Frequency [Hz]')

row[0].set_xlabel('Time [s]')
row[1].set_xlabel('Time [s]')

pl.suptitle('%s Hilbert Huang Transform'%outname)

fig.subplots_adjust(hspace=0)

pl.savefig('%s_princcomps_HHT.png'%outname)

# --- Time series
fig,ax=pl.subplots(np.shape(data['PCs_final'])[1],figsize=(8,12),sharex='col')
for r,row in enumerate(ax):
    row.plot(data['PCs_final'][:,r]/max(data['PCs_final'][:,r]))
    row.set_ylim(-1.1,1.1)
    row.set_yticklabels('')
pl.xlabel('Time sample')
pl.suptitle('%s [arb y-units]'%outname)
pl.subplots_adjust(hspace=0.2)
pl.savefig('%s_princcomps_TD.png'%outname)

# --- PSDs
fig,ax=pl.subplots(np.shape(data['PCs_final'])[1],figsize=(8,12),sharex='col')
for r,row in enumerate(ax):
    freq, Pxx_den = signal.periodogram(data['PCs_final'][:,r], fs)
    row.plot(freq,Pxx_den/max(Pxx_den))
    row.set_ylim(0,1.1)
    row.set_xlim(5,50)
    row.axvline(10,color='k',label=r'$f_{\mathrm{low}}$')
pl.xlabel('Frequency [Hz]')
pl.suptitle('%s [arb y-units]'%outname)
pl.subplots_adjust(hspace=0.2)

pl.savefig('%s_princcomps_PSD.png'%outname)

# --- Whitened PSDs

# Interpolate the noise PSD to data freqs
Sf_interp=np.interp(freq,noise_freqs,Sf)

fig,ax=pl.subplots(np.shape(data['PCs_final'])[1],figsize=(8,12),sharex='col')
for r,row in enumerate(ax):
    freq, Pxx_den = signal.periodogram(data['PCs_final'][:,r], fs)
    row.plot(freq,Pxx_den/Sf_interp)
    row.set_yticklabels('')
    row.set_xlim(5,100)
    row.axvline(10,color='k',label=r'$f_{\mathrm{low}}$')
pl.xlabel('Frequency [Hz]')
pl.suptitle('Whitened %s [arb y-units]'%outname)
pl.subplots_adjust(hspace=0.2)

pl.savefig('%s_princcomps_whitenedPSD.png'%outname)

# --- HHT

# --- Complex Spectra
fig,ax=pl.subplots(np.shape(data['PCs_final'])[1],ncols=2,figsize=(8,12),sharex='col')
for r,row in enumerate(ax):

    # Take FFT
    Fspec=np.fft.fft(data['PCs_final'][:,r])
    freq=np.fft.fftfreq(len(data['PCs_final'][:,r]),1./fs)

    # Take +ve freq
    Fspec=Fspec[freq>=0]
    freq=freq[freq>=0]

    row[0].plot(freq,Fspec/max(Fspec))
    row[0].set_ylim(-1.1,1.1)
    row[0].set_xlim(0,100)
    row[0].axvline(10,color='k',label=r'$f_{\mathrm{low}}$')
    row[1].plot(freq,Fspec/max(Fspec))
    row[1].set_ylim(-1.1,1.1)
    row[1].set_xlim(0,100)
    row[0].set_yticklabels('')
    row[1].set_yticklabels('')
    row[1].axvline(10,color='k',label=r'$f_{\mathrm{low}}$')

    if r==0:
        row[0].set_title('Re[FFT]')
        row[1].set_title('Im[FFT]')

row[0].set_xlabel('Frequency [Hz]')
row[1].set_xlabel('Frequency [Hz]')
pl.suptitle('%s [arb y-units]'%outname)

fig.subplots_adjust(hspace=0)

pl.savefig('%s_princcomps_FD.png'%outname)

# --- Whitened Complex Spectra
fig,ax=pl.subplots(np.shape(data['PCs_final'])[1],ncols=2,figsize=(8,12),sharex='col')
for r,row in enumerate(ax):

    # Take FFT
    Fspec=np.fft.fft(data['PCs_final'][:,r])
    freq=np.fft.fftfreq(len(data['PCs_final'][:,r]),1./fs)

    # Take +ve freq
    Fspec=Fspec[freq>=0]
    freq=freq[freq>=0]

    # Interpolate the noise PSD to data freqs
    Sf_interp=np.interp(freq,noise_freqs,Sf)

    # whitened complex spectrum
    FspecWhite=Fspec/np.sqrt(Sf_interp)

    row[0].plot(freq,np.real(FspecWhite)/max(np.real(FspecWhite)))
    row[0].set_ylim(-1.1,1.1)
    row[0].set_xlim(0,100)
    row[0].axvline(10,color='k',label=r'$f_{\mathrm{low}}$')
    row[1].plot(freq,np.imag(FspecWhite)/max(np.imag(FspecWhite)))
    row[1].set_ylim(-1.1,1.1)
    row[1].set_xlim(0,100)
    row[0].set_yticklabels('')
    row[1].set_yticklabels('')
    row[1].axvline(10,color='k',label=r'$f_{\mathrm{low}}$')

    if r==0:
        row[0].set_title('Re[FFT]')
        row[1].set_title('Im[FFT]')

row[0].set_xlabel('Frequency [Hz]')
row[1].set_xlabel('Frequency [Hz]')
pl.suptitle('%s [arb y-units]'%outname)

fig.subplots_adjust(hspace=0)

pl.savefig('%s_princcomps_whitenedFD.png'%outname)

