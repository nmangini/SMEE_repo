#!/usr/bin/env python
# plot_betas.py

#import scipy.io as sio 
#sio.savemat('pmns_catalogue.mat',{'waveform_catalogue':waveform_catalogue_td})

import sys
from glob import glob
import numpy as np
import scipy.io as sio 
from matplotlib import pyplot as pl

def beta_priors(idx):
    """
    Define a list so we can easily draw the prior limits on the beta
    posterior plots (this might be in the output file anyway but this is fairly
    handy)
    """
#    beta_lims = [
#            [-6,7],
#            [-7,7],
#            [-7,8],
#            [-9,3],
#            [-7,6],
#            [-5,9],
#            [-5,8],
#            [-9,5],
#            [-5,8],
#            [-5,8],
#            [-8,3]
#            ]

    beta_lims = [
            [-50,50],
            [-50,50],
            [-50,50],
            [-50,50],
            [-50,50],
            [-50,50],
            [-50,50],
            [-50,50],
            [-50,50],
            [-50,50],
            [-50,50]
            ]

    return beta_lims[idx]

#
# Read file list
#
globpattern=sys.argv[1]

matfiles=glob(globpattern+'*/*mat')
if len(matfiles)==0:
    print >> sys.stderr, \
            "error, no mat files found matching %s*/*mat"%globpattern

for m,matfile in enumerate(matfiles):

    # strip the mat extension so we can conveniently name the figure
    #savename=matfile.replace('.mat','_postbetas.png')
    savename= globpattern + str(m+1) + '_postbetas.png'

    #
    # Load file
    #
    data=sio.loadmat(matfile)

    postbetas=data['postbetas']
    nbetas=np.shape(postbetas)[1]

    fig, axes = pl.subplots(nrows=nbetas,ncols=2,figsize=(10,10))
#    fig.tight_layout()
    for r,row in enumerate(axes):

        # Histogram posterior samples
        row[0].hist(data['postbetas'][:,r])
        # Plot the samples
        row[1].plot(data['postbetas'][:,r],'.')

        # draw the prior range
        prior_range = beta_priors(r)
        row[0].axvline(prior_range[0],color='r')
        row[0].axvline(prior_range[1],color='r')
        row[0].set_xlim(1.1*prior_range[0],1.1*prior_range[1])

        row[1].axhline(prior_range[0],color='r')
        row[1].axhline(prior_range[1],color='r')
        row[1].set_ylim(1.1*prior_range[0],1.1*prior_range[1])

    pl.savefig(savename)



