#!$FIREMODELS/fds/.github/fds_python_env/bin/python3
# McDermott
# 2 April 2024

import subprocess

# Scripts to run prior to dataplot

print("ignition_delay...");   subprocess.run(["python3","./scripts/ignition_delay.py"])

# Dataplot and scatplot options

# Statistics output options

# Run dataplot and scatplot scripts

# Special cases

print("verification scripts completed successfully!")
