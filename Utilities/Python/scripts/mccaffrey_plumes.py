# -*- coding: utf-8 -*-
"""
Created on Fri May 17 08:02:55 2024

@author: jhodges
"""

import numpy as np
import pandas as pd
import scipy.interpolate
import os, sys, inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

from fdsplotlib import plot_to_fig

def getCase(resolution, Q):
    if Q == 22: qnamespace = '22kW'
    elif Q == 45: qnamespace = '45kW'
    if resolution == 0.050: rnamespace = '50mm'
    elif resolution == 0.025: rnamespace = '25mm'
    elif resolution == 0.012: rnamespace = '12mm'
    elif resolution == 0.006: rnamespace = '06mm'
    expname = '../../../../exp/McCaffrey_Plume/McCaffrey_%s.csv'%(qnamespace)
    resultDir = '../../../../out/McCaffrey_Plume/'
    chid = 'McCaffrey_%s_%s'%(qnamespace,rnamespace)
    return resultDir, chid, expname

def calculateV(Q, z, g=9.81):
    try:
        z[0]
        z = np.array(z)
    except:
        z = np.array([z])
    # Drysdale 2011 Table 4.2
    heightFactor = z/(Q**(2/5))
    u = np.zeros_like(z)
    for i, z0 in enumerate(z):
        if heightFactor[i] <0.08: #z <= H*0.6:
            # Continuous flaming region
            k = 6.8 # m(1/2)/s
            eta = 0.5
        elif heightFactor[i] <= 0.2: #z <= H*1.2:
            # Intermittent flaming region
            k = 1.9 # m/kW(1/5)s
            eta = 0
        else:
            # Plume above flame
            k = 1.1 # m(4/3)/kW(1/3)s
            eta = -1/3
        u[i] = k*Q**(1/5)*(z[i]/Q**(2/5))**eta
    return u

def getExperimentalData(Q, quantity):
    if quantity == 'Temperature':
        _, _, expname = getCase(0.05, Q)
        d = pd.read_csv(expname)
        X = d['X']
        Z = d['Z']
        T = d['T']
        x, z = np.meshgrid(np.linspace(0,0.2,100), np.linspace(0,1.5,100))
        points = np.array([X.values, Z.values]).T
        data = scipy.interpolate.griddata(points, T, (x, z), method='linear')
        levels = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200]
    elif quantity == 'Velocity':
        alpha = 0.16
        z = np.linspace(0.0,1.5,101)
        v0 = calculateV(Q, z, g=9.81)
        
        z_over_Q2p5 = z/(Q**0.4)
        
        x1e = 0.46*z_over_Q2p5+0.013*(Q**0.5)
        
        x = np.linspace(0, 0.2,101)
        
        data = np.zeros((z.shape[0], x.shape[0]))
        for i in range(0, z.shape[0]):
            x_over_z = x/z[i] if z[i] > 0 else 0
            if z_over_Q2p5[i] >= 0.2:
                v_over_v0 = np.exp(-(((5/6)/alpha*x_over_z)**2))
            elif z_over_Q2p5[i] >= 0.08:
                v_over_v0 = np.exp(-((x/x1e[i])**(3/2)))
            else:
                factor = 2.5-z_over_Q2p5[i]/0.08
                v_over_v0 = np.exp(-((x/x1e[i])**(factor)))
            data[i, :] = v_over_v0 * v0[i]
        levels = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0]
        x, z = np.meshgrid(x, z)
    return x, z, data, levels

def getModelData(Q, resolution, quantity, x, z):
    resultDir, chid, expname = getCase(resolution, Q)
    
    d_mod = pd.read_csv(resultDir + os.sep+chid+"_line.csv", header=1)
    columns = d_mod.columns
    columns = [x.replace('velocity-z','W-0p000-z').replace('velocity','W-0p000').replace('temperature','T-0p000') for x in columns]
    d_mod.columns = columns
    locs = [i for i in d_mod.columns if ('W' in i and ('z') in i)]
    
    zms = []
    xms = []
    Tms = []
    Wms = []
    for loc in locs:
        zm = d_mod[loc].values
        xm = np.zeros_like(zm) + float(loc.split('-')[1].replace('p','.'))
        Tm = d_mod[loc.replace('W','T').replace('-z','')].values
        Wm = d_mod[loc.replace('-z','')].values
        
        xms.extend(xm)
        zms.extend(zm)
        Tms.extend(Tm)
        Wms.extend(Wm)
    
    xms = np.array(xms)
    zms = np.array(zms)
    Tms = np.array(Tms)
    pointsm = np.array([xms, zms]).T
    
    if quantity == 'Temperature':
        model_data = scipy.interpolate.griddata(pointsm, Tms, (x, z), method='linear')
    elif quantity == 'Velocity':
        model_data = scipy.interpolate.griddata(pointsm, Wms, (x, z), method='linear')
    with open(resultDir+os.sep+chid+"_git.txt","r") as f: version_string = (f.read()).strip()
    return model_data, version_string

if __name__ == "__main__":
    resolutions = [0.050, 0.025, 0.012] #, 0.006]
    
    pltdir = '../../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/McCaffrey_Plume/'
    isDir = os.path.isdir(pltdir)
    if not isDir:
        os.mkdir(pltdir)
    for Q in [45]:
        for quantity in ['Velocity','Temperature']:
            x, z, data, levels = getExperimentalData(Q, quantity)
            
            f = None
            for i, resolution in enumerate(resolutions):
                model_data, version_string = getModelData(Q, resolution, quantity, x, z)
                f = plot_to_fig(x, z, contour_data=data, data_label='Exp',
                                plot_type='line_contour', levels=levels, line_color='k',
                                x_label='Burner Radius (m)', y_label='Height (m)', 
                                x_min=0.0, x_max=0.2, x_nticks=11, figure_handle=f)
                f = plot_to_fig(x, z, contour_data=model_data, data_label='FDS', 
                                show_legend=True, legend_location=4,
                                plot_type='line_contour', levels=levels, line_color='r',
                                x_label='Burner Radius (m)', y_label='Height (m)', 
                                revision_label=version_string,
                                figure_handle=f)
                f.savefig(pltdir + 'McCaffrey_Radial_Plume_%s_%s_kW_%0.0fmm.pdf'%(quantity,Q, resolution*1e3), backend='pdf')
                f = None
            