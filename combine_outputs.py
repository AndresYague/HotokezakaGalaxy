# Python packages
import glob

# Path results
path = './'

# Read list of times
with open(path+'tArray.in', 'r') as ff:
    ff.readline() # Skip header
    line = ff.readline() # Actual values
    
    line = (line.replace('D','E')).split()
    tArray = line

# Read list of tau
with open(path+'tauList.in', 'r') as ff:
    ff.readline() # Skip header
    line = ff.readline() # Actual values
    
    # Create one file for each tau
    line = line.replace('D','e')
    tauList = line.split()
    for tau in tauList:
        with open(path + "tau_{}.txt".format(tau), "w") as fwrite:
            for tt in tArray:
                fwrite.write("{} ".format(tt))

# Get output files
files = sorted(glob.glob("%s/Output*" %path))
    
# Read evolution of radioactive abundances
for the_file in files:
    with open(the_file, 'r') as ff:
        not_header = False
        i_tau = -1
        for line in ff:
            if '#' in line:
                i_tau += 1
                if i_tau > 0:
                    fwrite.close()
                not_header = True
                fwrite = open(path + "tau_{}.txt".format(tauList[i_tau]), "a")
            elif not_header:
                line = " ".join(line.split()[1:])
                fwrite.write("\n" + line)
        
        fwrite.close()
