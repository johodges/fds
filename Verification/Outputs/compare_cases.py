import matplotlib.pyplot as plt
import numpy as np
import pyfdstools as fds
import os

if __name__ == "__main__":
    chid = "duct_flow_both"
    quantity = "U-VELOCITY"
    working_dir = os.path.dirname(__file__) + os.sep
    
    time, dt = 30, None
    tol=1e-10
    qnty_mn, qnty_mx, cbarnumticks, levels = -2, 2, 11, 100
    axis, value = 3, 4.6
    
    # make directories
    for directory in ['SMV', 'VTK']:
        try:
            os.mkdir(directory)
        except:
            pass
    
    # read smv slice
    data, unit = fds.query2dAxisValue(working_dir, chid, quantity, axis, value, time=time, dt=dt)
    clabel = "%s (%s)"%(quantity, unit)
    fig1, ax1 = plt.subplots(1, 1, figsize=(8, 6))
    fig1, ax1 = fds.plotSlice(data['x'], data['z'], data['datas'][:, :, -1], axis, 
                            levels=levels, qnty_mn=qnty_mn, qnty_mx=qnty_mx, cbarnumticks=cbarnumticks, fig=fig1, ax=ax1, clabel=clabel)
    fig1.savefig(working_dir + os.sep + "SMV" + os.sep + chid+'.png', dpi=300)
    
    # read vtkhdf slice
    data_vtk = fds.query2dAxisValue_vtkhdf(working_dir, chid, quantity, axis, value, time=time, dt=dt)
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    fds.plotSlice(data_vtk['x'], data_vtk['z'], data_vtk['datas'], axis, 
                  levels=levels, qnty_mn=qnty_mn, qnty_mx=qnty_mx, cbarnumticks=cbarnumticks, fig=fig, ax=ax, clabel=clabel)
    fig.savefig(working_dir + os.sep + "VTK" + os.sep + chid+'.png', dpi=300)
    