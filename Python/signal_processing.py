#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 21:24:17 2018

@author: Zhaoyi.Shen
"""
import sys
sys.path.append('/home/z1s/py/lib/')
from lanczos_filter import lanczos_filter
import numpy as np
import scipy as sp
from scipy.signal import butter, lfilter, filtfilt

def lfca(x, cutoff, truncation, scale, **kwargs):
    if x.ndim!=2:
        return
    if 'covtot' in kwargs.keys():
        covtot = kwargs['covtot']
    else:
        covtot = np.cov(x,rowvar=False)
    (n,p) = x.shape
    if covtot.shape!=(p,p):
        return
    
    # center data
    x = x - np.nanmean(x,0)[np.newaxis,...]
    xs = x * np.transpose(scale)
    
    # eigendecomposition of covariance matrix
    #scale.shape = (1,p)
    covtot = np.transpose(scale)*covtot*scale
    pcvec, evl, rest = peigs(covtot, min(n-1,p))
    trcovtot = np.trace(covtot)
    #scale.shape = (p)
    
    # percent of total sample variation accounted for by ead EOF
    pvar = evl/trcovtot*100
    # principal component time series
    pcs = np.dot(xs, pcvec)
    # return EOFs in original scaling as patterns (row vectors)
    eof = np.transpose(pcvec)/np.transpose(scale)
                      
    # truncation of EOFs
    ntr = truncation
    
    # whitening transformation
    f = np.sqrt(np.squeeze(evl)[0:ntr])
    # get transformation matrices
    s = np.dot(pcvec[:,0:ntr], np.diag(1./f))
    sadj = np.dot(np.diag(f), np.transpose(pcvec[:,0:ntr]))
    
    # filter data matrix
    b,a = butter(5,1./cutoff,btype='low')
    t = np.arange(1,n+1)
    #t.shape = (1,n)
    #t = np.transpose(t)
    x_f = xs.copy()
    for i in range(xs.shape[1]):
        p = np.polyfit(t,xs[:,i],1)
        tmp = xs[t-1,i]-p[0]*t-p[1]
        tmp1 = np.concatenate((np.flipud(tmp),tmp,np.flipud(tmp)))
        #tmp_filt = filtfilt(b,a,tmp)
        tmp_filt = lanczos_filter(tmp1,1,1./cutoff)[0]
        x_f[:,i] = tmp_filt[np.size(tmp_filt)/3:2*np.size(tmp_filt)/3]+p[0]*t+p[1]
        #x_f[:,i] = tmp_filt+p[0]*t+p[1]
        
    # whiten variables
    y = np.dot(x_f, s)
    # slow covariance matrix of whitened variables
    gamma = np.cov(y,rowvar=False)
    # SVD of slow variance matrix
    dummy, r, v = csvd(gamma)
    
    # weight vectors and patterns
    weights = scale * np.dot(s, v)
    lfps = np.dot(np.transpose(v), sadj)/np.transpose(scale)
                 
    # choose signs of patterns, weights, eofs, and pcs
    #scale.shape = (1,p)
    for j in range(lfps.shape[0]):
        if np.dot(lfps[j,:][np.newaxis,...], scale)<0:
            lfps[j,:] = -lfps[j,:]
            weights[:,j] = -weights[:,j]
    for j in range(eof.shape[0]):
        if np.dot(eof[j,:][np.newaxis,...], scale)<0:
            eof[j,:] = -eof[j,:]
            pcs[:,j] = -pcs[:,j]
    #scale.shape = (p)
            
    # low-frequency components
    xs = xs/np.transpose(scale)
    lfcs = np.dot(xs, weights)
    
    # slow covariance of untruncated state space
    cov_slow = np.cov(x_f,rowvar=False)
    trcovslow = np.trace(cov_slow)
    w = weights/scale
    p = lfps*np.transpose(scale)
    
    pw_diag = np.diag(np.dot(p,w))
    slow_var = np.diag(np.dot(np.dot(p,cov_slow),w))/pw_diag
    tot_var = np.diag(np.dot(np.dot(p,covtot),w))/pw_diag
    
    pcvec_diag = np.diag(np.dot(np.transpose(pcvec),pcvec))
    slow_var_eofs = np.diag(np.dot(np.dot(np.transpose(pcvec),cov_slow),pcvec))/pcvec_diag
    tot_var_eofs = np.diag(np.dot(np.dot(np.transpose(pcvec),covtot),pcvec))/pcvec_diag
                          
    # slow variance and total variance in each LFC
    pvar_slow = slow_var/trcovslow*100
    pvar_lfc = tot_var/trcovtot*100
    r_eofs = slow_var_eofs/tot_var_eofs
    pvar_slow_eofs = slow_var_eofs/trcovslow*100
    
    return lfcs, lfps, weights, r, pvar, pcs, eof, ntr, pvar_slow, pvar_lfc, r_eofs, pvar_slow_eofs

def peigs(a, rmax):
    
    (m,n) = a.shape
    if rmax>min(m,n):
        rmax = min(m,n)
    
    if rmax<min(m,n)/10.:
        (d,v) = sp.sparse.linalg.eigs(a, rmax)
    else:
        (d,v) = np.linalg.eig(a)
    
    if d.size>max(d.shape):
        d = np.diag(d)
    
    # ensure that eigenvalues are monotonically decreasing    
    i = np.argsort(-d)
    d = -np.sort(-d)
    v = v[:,i]
    # estimate number of positive eigenvalues of a
    d_min = max(d)*max(m,n)*np.spacing(1)
    r = np.sum(d>d_min)
    # discard eigenpairs with eigenvalues that are close to or less than zero
    d = d[:r]
    v = v[:,:r]
    d = d[:]
    
    return v, d, r

def csvd(a):
    
    (m,n) = a.shape
    if m>=n:
        (u,s,v) = np.linalg.svd(a,0)
        v = np.transpose(v)
    else:
        (v,s,u) = np.linalg.svd(a.transpose(),0)
        u = np.transpose(u)
        
    return u, s, v
    
    