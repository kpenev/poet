"""
Generates the MESA stellar evolution tracks used by the POET package --

2013 Nov 12 -- Written by Brian Jackson (decaelus@gmail.com)

"""

import sys, os, shutil as shu, update_inlist as ul, numpy as np

def main():

  delta_mass = 0.05
  min_mass = 0.5
  max_mass = 1.2

  stellar_masses = np.arange(min_mass, max_mass + delta_mass, delta_mass)

  initial_filename = 'inlist_project'

  for i in range(0, len(stellar_masses)):
    inlist_filename = 'inlist_project_initial-mass' + str(stellar_masses[i])

    if(not os.path.isfile(inlist_filename)):
      shu.copyfile(initial_filename, inlist_filename)

    command_string = 'python update_inlist.py -m initial_mass ' + \
      str(stellar_masses[i]) + ' ' + inlist_filename
#   print command_string
    os.system(command_string)

    command_string = "./mk"
    os.system(command_string)

    log_file = "LOGS_initial-mass" + str(stellar_masses[i])

    #Check whether model has already been run
    if(not os.path.isdir(log_file)):
      print "Running MESA with initial_mass " + str(stellar_masses[i])
      command_string = "./rn"
      os.system(command_string)

      print "Move tracks"
      command_string = "cp -r LOGS " + log_file
      os.system(command_string)

if __name__=="__main__":
    main()

