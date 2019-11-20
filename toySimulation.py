import numpy as np
import matplotlib.pyplot as plt
import os

def rateFunction(r0, rAvg, totTime):
    '''
    Define the change of rate with time
    -r0 is the present rate in events/Myr
    -rAvg is the average rate in events/Myr
    -totTime is the total simulation time in Myr before present
    '''
    
    # Suppose a linear function
    rate = lambda t: rAvg + (rAvg - r0)*(1 - 2*t/totTime)
    
    return rate

class SimulationObj():
    '''
    The instantiation values for this class should be:

    -alpha: the mixing lenght parameter
    -vt: the turbulent velocity (in km/s)
    -hscale: the height scale (in kpc)
    -circumf: the total circle length considered (in kpc)
    -width: the circle width (in kpc)
    -time: the total time of the simulation (in Myr)
    '''

    def __init__(self, alpha, vt, hscale, circumf, width, time):
        '''Initialization function'''
        
        # Calculate diffusion coefficient in kpc^2/Myr
        self.diff = alpha*vt/7*hscale/0.2*1e-3
        
        # Keep everything else the same
        self.hscale = hscale
        self.circumf = circumf
        self.width = width
        self.timeMyr = time
    
    def runSimulation(self, tau, rateFunc, sampleDt):
        '''Generate the events according to the parameters
        introduced by the user and the rate function'''
        
        # Get the simulation timepoints
        simulTime = np.arange(0, self.timeMyr + sampleDt, sampleDt)
        lenSimulTime = len(simulTime)
        simulOutpt = np.zeros(lenSimulTime)
        
        # Get the random events
        
        # For each Myr, produce R events randomly distributed in
        # said Myr:
        tt = 0; times = []; distances = []
        while tt <= self.timeMyr:
            nEvents = int(np.round(rateFunc(tt)*sampleDt))
            
            # Get the times
            times = np.append(times, np.sort(sampleDt*np.random.random(nEvents)) + tt)
            
            # Get positions
            xx = self.circumf*np.random.random(nEvents) - 0.5*self.circumf
            yy = self.width*np.random.random(nEvents) - 0.5*self.width
            zz = np.random.laplace(0, self.hscale, nEvents)
            
            # Calculate the distances from origin
            distances = np.append(distances, np.sqrt(xx**2 + yy**2 + zz**2))
            
            tt += sampleDt
        
        jj = 0
        # Update all values
        for ii in range(len(times)):
            if times[ii] > simulTime[jj]:
                jj += 1
                if jj > lenSimulTime - 1:
                    break
            
            dt = simulTime[jj:] - times[ii]
            r2 = self.diff*dt
            
            dilution = np.exp(-distances[ii]**2/(4*r2) - dt/tau)
            dilution /= np.minimum((4*np.pi*r2)**1.5, 8*np.pi*self.hscale*r2)
            
            simulOutpt[jj:] += dilution
        
        return simulTime, simulOutpt
    
def main():
    '''Toy simulation from Hotokezaka2015'''
    
    # Check that input file exists
    inFile = "inputToy.in"
    inExample = "inputToy.example"
    if not os.path.isfile(inFile):
        s = 'File "{}" not found, '.format(inFile)
        s += 'please use "{}" as a template'.format(inExample)
        print(s)
        return 1
    
    # Read input file
    inputArgs = {}
    with open("inputToy.in", "r") as fread:
        for line in fread:
            lnlst = line.split()
            
            # Ignore comments and empty lines
            if len(lnlst) == 0 or lnlst[0] == "#":
                continue
            
            # Fill all the arguments
            val = float(lnlst[0]); key = lnlst[1]
            inputArgs[key] = val
    
    # Check at least correct number of arguments
    if len(inputArgs) < 9:
        s = 'File "{}" has an incorrect number of entries, '.format(inFile)
        s += 'please use "{}" as a template'.format(inExample)
        print(s)
        return 1
    
    
    r0 = inputArgs["r0"]       # Present-day rate of events per Myr
    rate = rateFunction(r0, 5*r0, inputArgs["time"])
    
    # Set tau (in Myr)
    tau = inputArgs["tau"]
    
    # Set number of runs and sample Dt (in Myr)
    sampleDt = inputArgs["sampleDt"]
    
    # Initialize and run simulation
    simul = SimulationObj(alpha = inputArgs["alpha"],
                          vt = inputArgs["vt"],
                          hscale = inputArgs["hscale"],
                          circumf = inputArgs["circumf"],
                          width = inputArgs["width"],
                          time = inputArgs["time"])
    
    simulTime, simulOutpt = simul.runSimulation(tau, rate, sampleDt)
    
    # Plot results
    plt.plot(simulTime, simulOutpt)
    plt.xlim([2000, inputArgs["time"]])
    plt.ylim([1e0, 1e3])
    plt.yscale("log")
    
    plt.xlabel("Time (Myr)")
    plt.ylabel("Mass (Normalized)")
    plt.show()

if __name__ == "__main__":
    main()
