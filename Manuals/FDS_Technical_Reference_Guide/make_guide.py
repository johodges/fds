import os, subprocess, sys

guide = 'FDS_Technical_Reference_Guide'
firemodels = os.path.join(os.path.dirname(__file__),'..','..','..')
manual_dir = os.path.join(firemodels,'fds','Manuals',guide)
command = 'python ' + os.path.join(firemodels,'fds','Manuals','scripts','make_guide.py') + ' --file %s --clean'%(guide)

p = subprocess.Popen(command, cwd=manual_dir, stdout=sys.stdout, stderr=sys.stderr, shell=True, close_fds=True, env=os.environ)
txt = p.communicate()