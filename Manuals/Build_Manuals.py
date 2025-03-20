import os, subprocess, sys, shutil

def make_one_guide(firemodels, guide):
    manuals_dir = os.path.join(firemodels,'fds','Manuals')
    manual_dir = os.path.join(firemodels,'fds','Manuals',guide)
    command = 'python ' + os.path.join(firemodels,'fds','Manuals','scripts','make_guide.py') + ' --file %s --clean'%(guide)
    print("Building %s at: %s"%(guide.replace('_',' '), manual_dir))
    p = subprocess.Popen(command, cwd=manual_dir, stdout=sys.stdout, stderr=sys.stderr, shell=True, close_fds=True, env=os.environ)
    txt = p.communicate()
    
    manual_file = os.path.join(manual_dir, guide + '.pdf')
    manuals_file = os.path.join(manuals_dir, guide + '.pdf')
    shutil.copy2(manual_file, manuals_file)
    
firemodels = os.path.join(os.path.dirname(__file__),'..','..','..')

make_one_guide(firemodels, 'FDS_Config_Management_Plan')
make_one_guide(firemodels, 'FDS_User_Guide')
make_one_guide(firemodels, 'FDS_Technical_Reference_Guide')
make_one_guide(firemodels, 'FDS_Verification_Guide')
make_one_guide(firemodels, 'FDS_Validation_Guide')
