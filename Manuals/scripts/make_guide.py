import os, platform, sys, argparse, subprocess, re, glob

def run_command(command, working_dir, output_file, mode):
    working_dir = os.path.abspath(working_dir)
    if 'Windows' in platform.platform():
        #working_dir = working_dir.replace('\\','/')
        working_dir = working_dir[0].upper() + working_dir[1:]
        
    p = subprocess.Popen(command, cwd=working_dir, stdout=subprocess.PIPE, shell=True, close_fds=True, env=os.environ)
    txt = p.communicate()[0].decode('utf-8')
    txt = txt.replace('\r\n','\n')
    if 'Windows' in platform.platform():
        print(command)
        print(working_dir)
        print(txt)
        
    """Run a shell command and redirect output to a file."""
    if output_file != '':
        with open(output_file, mode) as f:
            f.write(txt)
        
    return txt

def check_errors_in_file(file_path, patterns, ignore_patterns=[]):
    """Check for specific error patterns in a log file."""
    errors_found = []
    with open(file_path, "r") as f:
        for line in f:
            if any(re.search(pattern, line) for pattern in patterns) and not any(re.search(ignore, line) for ignore in ignore_patterns):
                errors_found.append(line.strip())
    return errors_found

def get_manuals_datafile(firemodels, file):
    scripts_dir = os.path.join(firemodels, 'fds', 'Manuals', 'scripts')
    check_manuals_datafile = ''
    if file == 'FDS_User_Guide':
        check_manuals_datafile = '--datafile ' + os.path.join(scripts_dir,'files_to_check_usr.txt')
    elif file == 'FDS_Verification_Guide':
        check_manuals_datafile = '--datafile ' + os.path.join(scripts_dir,'files_to_check_ver.txt')
    elif file == 'FDS_Validation_Guide':
        check_manuals_datafile = '--datafile ' + os.path.join(scripts_dir,'files_to_check_val.txt')
    elif file == 'FDS_Technical_Reference_Guide':
        check_manuals_datafile = '--datafile ' + os.path.join(scripts_dir,'files_to_check_tech.txt')
    elif file == 'FDS_Config_Management_Plan':
        check_manuals_datafile = '--datafile ' + os.path.join(scripts_dir,'files_to_check_cfg.txt')
    return check_manuals_datafile

def delete_files_by_extension(working_dir, ext):
    files = glob.glob(os.path.join(working_dir,'*'+ext))
    for file in files:
        os.remove(file)
    
if __name__ == "__main__":

    args = sys.argv
    parser = argparse.ArgumentParser(prog='make_guide',
                                     description='builds FDS manual')
    parser.add_argument('call')
    parser.add_argument('--file', help='filename to build', required=True)
    parser.add_argument('--dir', help='directory which contains file', default='')
    parser.add_argument('--clean', help='delete latex temporary files at start', action='store_true')
    
    cmdargs = parser.parse_args(args)    
    
    # Initialize
    firemodels = os.path.join(os.path.dirname(__file__),'..','..','..')
    clean_build = True
    
    # Add LaTeX search path
    texinputs = os.path.abspath(os.path.join(firemodels, 'fds', 'Manuals', 'LaTeX_Style_Files')) + os.sep
    texinputs = '.:..%sLaTeX_Style_Files:'%(os.sep)
    os.environ["TEXINPUTS"] = texinputs
    
    # Set directory information
    tex_file = cmdargs.file
    if tex_file[-4:] == '.tex': tex_file = tex_file[-4:]
    output_log = tex_file + '.err'
    if cmdargs.dir != '':
        manual_dir = cmdargs.dir
    else:
        manual_dir = os.path.join(firemodels,'fds','Manuals',tex_file)
    os.chdir(manual_dir)
    
    # Clean if requested
    if cmdargs.clean:
        delete_files_by_extension(manual_dir, '.aux')
        delete_files_by_extension(manual_dir, '.lof')
        delete_files_by_extension(manual_dir, '.log')
        delete_files_by_extension(manual_dir, '.lot')
        delete_files_by_extension(manual_dir, '.out')
        delete_files_by_extension(manual_dir, '.pdf')
        delete_files_by_extension(manual_dir, '.toc')
        delete_files_by_extension(manual_dir, '.err')
        delete_files_by_extension(manual_dir, '.bbl')
        delete_files_by_extension(manual_dir, '.blg')
        delete_files_by_extension(manual_dir, '.brf')
    
    # Get Git revision
    git_revision = subprocess.getoutput("git describe --abbrev=7 --long --dirty")
    gitfile = os.path.join(firemodels, 'fds', 'Manuals', 'Bibliography','gitrevision.tex')
    with open(gitfile, "w") as f:
        f.write(f"\\newcommand{{\\gitrevision}}{{{git_revision}}}\n")
    
    # Run LaTeX build process
    mode = 'w'
    working_dir = os.path.join(firemodels, 'fds', 'Manuals', tex_file)
    for i in range(0, 4):
        print("pass %d"%(i+1))
        run_command("pdflatex -interaction nonstopmode %s"%(tex_file), working_dir, output_log, mode)
        if i == 0:
            run_command("bibtex %s"%(tex_file), working_dir, output_log, mode)
    
    # Check if the guide exists
    if not os.path.exists(f"{tex_file}.pdf"):
        clean_build = False
        print("***error: the %s failed to build!"%(tex_file.replace('_',' ')))
    
    # Scan and report any errors
    latex_error_patterns = [
        "Too many", "Undefined control sequence", "Error:", "Fatal error",
        "! LaTeX Error:", "Paragraph ended before", "Missing \\ inserted", "Misplaced"
    ]
    ignore_patterns = ["xpdf supports version 1.5"]
    errors = check_errors_in_file(output_log, latex_error_patterns, ignore_patterns)
    
    if errors:
        print("LaTeX errors detected:")
        print("\n".join(errors))
        clean_build = False
    
    # Check for warnings
    warning_patterns = ["undefined", "multiply defined", "multiply-defined"]
    warnings = check_errors_in_file(output_log, warning_patterns)
    
    if warnings:
        print("LaTeX warnings detected:")
        print("\n".join(warnings))
        clean_build = False
    
    # Run additional manual checks
    manual_check_log = "%s_py.err"%(tex_file)
    check_manuals_datafile = get_manuals_datafile(firemodels, tex_file)
    check_manuals_file = os.path.join(firemodels, 'fds', 'Manuals', 'scripts', 'check_manuals.py')
    
    command = "python %s %s --outname %s --suppressconsole"%(check_manuals_file, check_manuals_datafile, manual_check_log)
    run_command(command, working_dir, '', 'w')
    
    manual_errors = check_errors_in_file(manual_check_log, ["ERROR:", "WARNING:"])
    
    if manual_errors:
        print("Other errors, warnings, or misspellings identified:")
        print("\n".join(manual_errors))
        clean_build = False
    
    # Final status message
    if clean_build:
        print("%s built successfully!"%(tex_file.replace('_',' ')))
