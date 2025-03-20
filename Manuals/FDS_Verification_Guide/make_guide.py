import os, subprocess, sys, shutil

guide = 'FDS_Verification_Guide'
firemodels = os.path.join(os.path.dirname(__file__),'..','..','..')
manual_dir = os.path.join(firemodels,'fds','Manuals',guide)
command = 'python ' + os.path.join(firemodels,'fds','Manuals','scripts','make_guide.py') + ' --file %s --clean'%(guide)

userguide_aux = os.path.join(firemodels, 'fds','Manuals','FDS_User_Guide', 'FDS_User_Guide.aux')
verguide_aux = os.path.join(firemodels, 'fds','Manuals','FDS_Verification_Guide', 'FDS_User_Guide.aux')

if os.path.exists(userguide_aux) is False:
    print("***warning: $AUXUSER does not exist. Build the FDS User's")
    print("            Guide before building the FDS Verification Guide")
else:
    shutil.copy2(userguide_aux, verguide_aux)

p = subprocess.Popen(command, cwd=manual_dir, stdout=sys.stdout, stderr=sys.stderr, shell=True, close_fds=True, env=os.environ)
txt = p.communicate()