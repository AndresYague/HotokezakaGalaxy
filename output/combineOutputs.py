# Python packages
import glob

# Path results
pathInput = '../input/'
pathOutput = './'

# Read list of times
with open(pathInput+'tArray.in', 'r') as ff:
    ff.readline() # Skip header
    line = ff.readline() # Actual values

    line = (line.replace('D','E')).split()
    tArray = line

# Read list of tau
with open(pathInput+'tauList.in', 'r') as ff:
    ff.readline() # Skip header
    line = ff.readline() # Actual values

    # Create one file for each tau
    line = line.replace('D','e')
    tauList = line.split()
    for tau in tauList:
        with open(pathOutput + "tau_{}.out".format(tau), "w") as fwrite:
            for tt in tArray:
                fwrite.write("{} ".format(tt))

# Get output files
files = sorted(glob.glob("%s/Output*" %pathOutput))

# Read evolution of radioactive abundances
for the_file in files:
    with open(the_file, 'r') as ff:
        while True:
            line = ff.readline()
            if "#" in line:
                for tau in tauList:
                    fwrite = open(pathOutput + "tau_{}.out".format(tau), "a")
                    line = " ".join(ff.readline().split()[1:])
                    fwrite.write("\n" + line)
                    fwrite.close()

            elif len(line) == 0:
                break
