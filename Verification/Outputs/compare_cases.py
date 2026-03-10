import matplotlib.pyplot as plt
import numpy as np
import pyfdstools as fds
import os, subprocess

def safe_run(script_path):
    try:
        result = subprocess.run(
            ["bash", script_path, '-1', 'SMV', '-2', 'VTK'],
            capture_output=True,
            text=True
        )
        
        print(result.stdout)
    except Exception as exc:
        print(f"Error in {script_path}: {exc}")

def getCases():
    cases = {
        'duct_flow_both_uvelocity': {'chid': 'duct_flow_both', 'quantity': 'U-VELOCITY', 'dtype': 'slcf', 'axis': 3, 'value': 4.6, 'qnty_mn': -2, 'qnty_mx': 2, 'cbarnumticks': 11, 'levels': 100, 'time': 30, 'dt': None},
        'duct_flow_both_vvelocity': {'chid': 'duct_flow_both', 'quantity': 'V-VELOCITY', 'dtype': 'slcf', 'axis': 3, 'value': 4.6, 'qnty_mn': -2, 'qnty_mx': 2, 'cbarnumticks': 11, 'levels': 100, 'time': 30, 'dt': None},
            }
    return cases

def getCase(case, cases):
    chid = cases[case]['chid']
    dtype = cases[case]['dtype']
    axis = cases[case]['axis']
    value = cases[case]['value']
    quantity = cases[case]['quantity']
    qnty_mn = cases[case]['qnty_mn']
    qnty_mx = cases[case]['qnty_mx']
    cbarnumticks = cases[case]['cbarnumticks']
    levels = cases[case]['levels']
    time = cases[case]['time']
    dt = cases[case]['dt']
    name = case
    return name, chid, dtype, axis, value, quantity, qnty_mn, qnty_mx, cbarnumticks, levels, time, dt

def compare_case(working_dir, name, chid, dtype, axis, value, quantity, qnty_mn, qnty_mx, cbarnumticks, levels, time, dt):
    
    # read data
    if dtype == 'slcf':
        data, unit = fds.query2dAxisValue(working_dir, chid, quantity, axis, value, time=time, dt=dt)
        clabel = "%s (%s)"%(quantity, unit)
        print(working_dir, chid, quantity, axis, value)
        data_vtk = fds.query2dAxisValue_vtkhdf(working_dir, chid, quantity, axis, value, time=time, dt=dt)
    
    # generate smv figure
    fig1, ax1 = plt.subplots(1, 1, figsize=(8, 6))
    fig1, ax1 = fds.plotSlice(data['x'], data['z'], data['datas'][:, :, -1], axis, 
                            levels=levels, qnty_mn=qnty_mn, qnty_mx=qnty_mx, cbarnumticks=cbarnumticks, fig=fig1, ax=ax1, clabel=clabel)
    fig1.savefig(working_dir + os.sep + "SMV" + os.sep + name+'.png', dpi=300)
    
    # generate vtkhdf figure
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    fds.plotSlice(data_vtk['x'], data_vtk['z'], data_vtk['datas'], axis, 
                  levels=levels, qnty_mn=qnty_mn, qnty_mx=qnty_mx, cbarnumticks=cbarnumticks, fig=fig, ax=ax, clabel=clabel)
    fig.savefig(working_dir + os.sep + "VTK" + os.sep + name+'.png', dpi=300)
    

if __name__ == "__main__":
    working_dir = os.path.dirname(__file__) + os.sep
    # make directories
    for directory in ['SMV', 'VTK','Summary']:
        try:
            os.mkdir(directory)
        except:
            pass
    
    # Generate images for each case
    cases = getCases()
    for case in list(cases.keys()):
        name, chid, dtype, axis, value, quantity, qnty_mn, qnty_mx, cbarnumticks, levels, time, dt = getCase(case, cases)
        compare_case(working_dir, name, chid, dtype, axis, value, quantity, qnty_mn, qnty_mx, cbarnumticks, levels, time, dt)
        
    # Run image comparison
    compare_images_script = os.path.join(working_dir,'..','..','..','bot','Firebot','compare_images.sh')
    safe_run(compare_images_script)
    
    