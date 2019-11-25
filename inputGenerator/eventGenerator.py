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

class EventsObj():
    '''
    The instantiation values for this class should be:

    -alpha: the mixing lenght parameter
    -vt: the turbulent velocity (in km/s)
    -hscale: the height scale (in kpc)
    -circumf: the total circle length considered (in kpc)
    -width: the circle width (in kpc)
    -time: the total time of the simulation (in Myr)
    '''

    def __init__(self, hscale, circumf, width, time):
        '''Initialization function'''
        
        # Keep everything else the same
        self.hscale = hscale
        self.circumf = circumf
        self.width = width
        self.simulationTime = time
    
    def getEvents(self, rateFunc, sampleDt, nRuns):
        '''Generate the events according to the parameters
        introduced by the user and the rate function'''
        
        # Temp file for times:
        timesFile = "tempTime.in"
        fwriteTime = open(timesFile, "w")
        
        # Temp file for distances:
        distFile = "tempDist.in"
        fwriteDist = open(distFile, "w")
        
        # Divide the nRuns in chunks of 10
        tempLen = 0
        tempRuns = 10; totRuns = 0
        while True:
            if totRuns + tempRuns > nRuns:
                tempRuns = nRuns - totRuns
            if tempRuns == 0:
                break
            
            # Get the random events
            size = self.__getRuns(fwriteTime, fwriteDist, rateFunc, sampleDt,
                                  tempRuns)
            
            # Update tempLen and totRuns
            tempLen += size
            totRuns += tempRuns
        
        fwriteTime.close()
        fwriteDist.close()
        
        return size, timesFile, distFile
    
    def __getRuns(self, fwriteTime, fwriteDist, rateFunc, sampleDt, nRuns):
        '''Inner function to generate the files'''
        
        # Get the random events
        
        # For each Myr, produce R events randomly distributed in
        # said Myr:
        tt = 0
        times = [[] for x in range(nRuns)]
        distances = [[] for x in range(nRuns)]
        tempSampleDt = sampleDt
        while tt < self.simulationTime:
            nEvents = int(np.round(rateFunc(tt)*tempSampleDt))
            
            # Check that the number of events changes meaningfully
            tempSampleDt = sampleDt
            while tt + tempSampleDt <= self.simulationTime:
                nEvents = int(np.round(rateFunc(tt)*tempSampleDt))
                futureNEvents = int(np.round(rateFunc(tt +
                                                tempSampleDt)*tempSampleDt))
                if nEvents != futureNEvents:
                    break
                tempSampleDt *= 10
            
            if tt + tempSampleDt > self.simulationTime:
                tempSampleDt = self.simulationTime - tt
                nEvents = int(np.round(rateFunc(tt)*tempSampleDt))
            
            # Get the times
            times = np.append(times,
                    np.sort(tempSampleDt*np.random.random((nRuns, nEvents)) +
                             tt), axis = 1)
            
            # Get positions
            xx = self.circumf*np.random.random((nRuns, nEvents)) - 0.5*self.circumf
            yy = self.width*np.random.random((nRuns, nEvents)) - 0.5*self.width
            zz = np.random.laplace(0, self.hscale, (nRuns, nEvents))
            
            # Calculate the distances from origin
            distances = np.append(distances, np.sqrt(xx**2 + yy**2 + zz**2),
                    axis = 1)
            
            tt += tempSampleDt
        
        # Write to temporal files
        for timeArr in times:
            fwriteTime.write(" ".join([str(x) for x in timeArr]))
            fwriteTime.write("\n")
        
        for distArr in distances:
            fwriteDist.write(" ".join([str(x) for x in distArr]))
            fwriteDist.write("\n")
        
        return len(times[0])
    
def main():
    '''Monte Carlo event generator'''
    
    # Check that input file exists
    inFile = "inputEventGenerator.in"
    inExample = "inputEventGenerator.example"
    if not os.path.isfile(inFile):
        s = 'File "{}" not found, '.format(inFile)
        s += 'please use "{}" as a template'.format(inExample)
        print(s)
        return 1
    
    # Read input file
    inputArgs = {}
    with open(inFile, "r") as fread:
        for line in fread:
            lnlst = line.split()
            
            # Ignore comments and empty lines
            if len(lnlst) == 0 or lnlst[0] == "#":
                continue
            
            # Fill all the arguments
            val = float(lnlst[0]); key = lnlst[1]
            inputArgs[key] = val
    
    # Check at least correct number of arguments
    if len(inputArgs) < 8:
        s = 'File "{}" has an incorrect number of entries, '.format(inFile)
        s += 'please use "{}" as a template'.format(inExample)
        print(s)
        return 1
    
    # Calculate the solar circumference
    circumf = inputArgs["rSun"]*2*np.pi
    
    # Calculate the simulation fractional mass by approximating:
    # dM = dr*rsun*exp(-rsun/rd)/(rd*rd)
    rScaleSun = inputArgs["rSun"]/inputArgs["rd"]
    
    fraction = np.exp(-rScaleSun)*rScaleSun
    fraction *= inputArgs["width"]/inputArgs["rd"]
    
    # Get the rate
    r0 = inputArgs["r0"]*fraction
    rate = rateFunction(r0, 5*r0, inputArgs["time"])
    
    # Set number of runs and sample Dt for the rate function (in Myr)
    nRuns = int(inputArgs["nRuns"])
    sampleDt = inputArgs["sampleDt"]
    
    # Initialize simulation
    events = EventsObj(hscale = inputArgs["hscale"],
                       circumf = circumf,
                       width = inputArgs["width"],
                       time = inputArgs["time"])
    sizeRun, timesFile, distFile = events.getEvents(rate, sampleDt, nRuns)
    
    fileName = "../input/Run"
    fileName += "".join(["_{}_{:.2f}".format(key, inputArgs[key])
                        for key in inputArgs])
    fileName += ".in"
    with open(fileName, "w") as fwrite:
        fwrite.write("{} {} {}\n".format(nRuns, sizeRun, inputArgs["hscale"]))
        
        # Put all times
        with open(timesFile, "r") as fread:
            for line in fread:
                # Write times
                fwrite.write(line)
        
        # Now put all distances
        with open(distFile, "r") as fread:
            for line in fread:
                
                # Write times
                fwrite.write(line)
        
        # Write the production factors
        fwrite.write("{} {}\n{}".format(1, 1, 1))
    
    # Remove the temporal files
    os.remove(timesFile)
    os.remove(distFile)
    
    # Create a suitable tArray.in for this simulation
    step = 10
    tArray = np.arange(0, inputArgs["time"] + step, step)
    with open(fileName + "tArray.in", "w") as fwrite:
        fwrite.write("{}\n".format(len(tArray)))
        fwrite.write(" ".join([str(x) for x in tArray]))

if __name__ == "__main__":
    main()
