#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 13:37:16 2018

@author: Zhaoyi.Shen
"""
import numpy as np

def lanczos_filter(x, dt, cf):
    
    nf = 1./(2*dt)
    m = 100
    
    cf = cf/nf
    coef = np.squeeze(lanczos_filter_coef(cf,m)[:,0])
    window, ff = spectral_window(coef, np.size(x))
    ff = ff*nf
    xmean = np.nanmean(x)
    x[np.where(np.isnan(x))] = xmean
    y, cx = spectral_filtering(x, window)
    
    return y, coef, window, cx, ff

def lanczos_filter_coef(cf, m):
    
    hkcs = lowpass_cosine_filter_coef(cf, m)
    sigma = np.concatenate(([1],np.sin(np.pi*np.arange(1,m+1)/m)/(np.pi*np.arange(1,m+1)/m)))
    hkb = hkcs*sigma
    hka = -hkb
    hka[0] = hka[0]+1
    coef = np.concatenate((hkb[:,np.newaxis], hka[:,np.newaxis]),axis=1)
    
    return coef
    
def lowpass_cosine_filter_coef(cf, m):
    
    coef = cf*np.concatenate(([1],np.sin(np.pi*np.arange(1,m+1)*cf)/(np.pi*np.arange(1,m+1)*cf)))
    return coef

def spectral_window(coef, n):
    
    ff = np.arange(0,1+np.spacing(1),2./n)
    window = np.zeros(np.size(ff))
    for i in range(np.size(ff)):
        window[i] = coef[0] + 2*np.sum(coef[1:]*np.cos(np.arange(1,np.size(coef))*np.pi*ff[i]))
    
    return window, ff

def spectral_filtering(x, window):
    
    nx = np.size(x)
    cx = np.fft.fft(x)[:np.int(np.floor(nx/2.)+1)]
    cxh = np.zeros(nx, dtype=np.complex_)
    cxh[:np.size(cx)] = cx*window
    cxh[np.size(cx):nx] = np.conj(cxh[np.arange(nx-np.size(cx),0,-1)])
    y = np.real(np.fft.ifft(cxh))
    
    return y, cx