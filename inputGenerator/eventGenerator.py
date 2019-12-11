import numpy as np
import os

def cucciatiFunctionStep(r0, totTime, lst):
    '''
    Define the change of rate with time according to Cucciati et al. 2012
    -r0 is the present rate in events/Myr
    -totTime is the total simulation time in Myr before present
    '''
    
    # Hubble constant (km/(s Mpc))
    H0 = 70
    
    # Inverse of transformed constant in 1/Myr
    H0m1 = 1/(H0*3.24e-7*np.pi)
    
    # Total life of the universe in Myr
    lifeUniverse = 14000
    
    # Function from time to redshift
    zz = lambda t: np.sqrt(2*H0m1/(lifeUniverse + (t - totTime)) - 1) - 1
    
    def rate(t):
        zt = zz(t)
        for elem in lst:
            if elem[0] < zt and zt <= elem[1]:
                return 10**elem[2]*r0/(10**lst[0][2])
    
    # Give the rate normalized to current rate (r0)
    return rate

def cucciatiFunctionPolyfit(r0, totTime, lst):
    '''
    Define the change of rate with time according to Cucciati et al. 2012
    -r0 is the present rate in events/Myr
    -totTime is the total simulation time in Myr before present
    '''
    
    # Hubble constant (km/(s Mpc))
    H0 = 70
    
    # Inverse of transformed constant in 1/Myr
    H0m1 = 1/(H0*3.24e-7*np.pi)
    
    # Total life of the universe in Myr
    lifeUniverse = 14000

    # List with the average bin in redshift
    xx = [(x[0] + x[1])*0.5 for x in lst]
    
    # List with the rate values
    yy = [x[2] for x in lst]
    
    # Fit with parabola
    fitCoef = np.polyfit(xx, yy, 3)
    pp = np.poly1d(fitCoef)
    
    # Function from time to redshift
    zz = lambda t: np.sqrt(2*H0m1/(lifeUniverse + (t - totTime)) - 1) - 1
    
    # Give the rate normalized to current rate (r0)
    return lambda t: 10**pp(zz(t))*r0/10**(-1.65)

def cucciatiFunctionInterpol(r0, totTime, lst0):
    '''
    Define the change of rate with time according to Cucciati et al. 2012
    -r0 is the present rate in events/Myr
    -totTime is the total simulation time in Myr before present
    '''
    
    # Hubble constant (km/(s Mpc))
    H0 = 70
    
    # Inverse of transformed constant in 1/Myr
    H0m1 = 1/(H0*3.24e-7*np.pi)
    
    # Total life of the universe in Myr
    lifeUniverse = 14000
    
    lst = [[lst0[0][0], lst0[0][2]]]
    for elem in lst0:
        lst.append([(elem[0] + elem[1])*0.5, elem[2]])
    lst.append([lst0[-1][1], lst0[-1][2]])
    
    # Function from time to redshift
    zz = lambda t: np.sqrt(2*H0m1/(lifeUniverse + (t - totTime)) - 1) - 1
    
    def rate(t):
        zt = zz(t)
        for ii in range(len(lst) - 1):
            if lst[ii][0] <= zt and zt < lst[ii + 1][0]:
                x0 = lst[ii][0]; x1 = lst[ii + 1][0]
                y0 = lst[ii][1]; y1 = lst[ii + 1][1]
                
                mm = (y1 - y0)/(x1 - x0)
                yy = mm*zt + y1 - mm*x1
                return 10**yy*r0/(10**lst[0][1])
    
    # Give the rate normalized to current rate (r0)
    return rate

def wandermanFunction(r0, totTime):
    '''
    Define the change of rate with time according to Wanderman et al. 2015
    -r0 is the present rate in events/Myr
    -totTime is the total simulation time in Myr before present
    '''
    
    # Hubble constant (km/(s Mpc))
    H0 = 60
    
    # Inverse of transformed constant in 1/Myr
    H0m1 = 1/(H0*3.24e-7*np.pi)
    
    # Total life of the universe in Myr
    lifeUniverse = 14000
    
    # Function from time to redshift
    zz = lambda t: np.sqrt(2*H0m1/(lifeUniverse + (t - totTime)) - 1) - 1
    
    # Give the rate normalized to current rate (r0)
    def rate(t):
        normConst = r0*np.exp(0.9/0.39)
        if zz(t) > 0.9:
            return np.exp(-(zz(t) - 0.9)/0.26)*normConst
        else:
            return np.exp((zz(t) - 0.9)/0.39)*normConst
    
    return rate

def hopiknsFunction(r0, totTime):
    '''
    Define the change of rate with time according to Hopkins et al. 2006
    -r0 is the present rate in events/Myr
    -totTime is the total simulation time in Myr before present
    '''
    
    # Parameters for function
    a=0.017; b=0.13
    c=3.3; d=5.3; h=0.7
    
    # Hubble constant (km/(s Mpc))
    H0 = 70
    
    # Inverse of transformed constant in 1/Myr
    H0m1 = 1/(H0*3.24e-7*np.pi)
    
    # Total life of the universe in Myr
    lifeUniverse = 14000
    
    # Function from time to redshift
    z = lambda t: np.sqrt(2*H0m1/(lifeUniverse + (t - totTime)) - 1) - 1
    
    # Rate normalized to present (r0)
    rate = lambda t: (a + b*z(t))*h/(1 + (z(t)/c)**d) * (r0/a*h)
    
    return rate

def linearFunction(r0, rAvg, totTime):
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
    
    def getEvents(self, rateFunc, nRuns):
        '''Generate the events according to the parameters
        introduced by the user and the rate function'''
        
        # Temp file for times:
        timesFile = "tempTime.in"
        fwriteTime = open(timesFile, "w")
        
        # Temp file for distances:
        distFile = "tempDist.in"
        fwriteDist = open(distFile, "w")
        
        # First integrate the function to know how many events in total
        maxDistrib = 0; sampleDt = 10
        prevNum = None; totNum = 0
        while prevNum != totNum:
            halfSample = sampleDt*0.5
            
            # Advance prevNum
            prevNum = totNum
            
            # Do a trapezoidal quadrature with sampleDt
            totNum = 0; tt = 0
            fun0 = rateFunc(0); fun1 = rateFunc(sampleDt)
            while True:
                # Integrate step
                totNum += (fun0 + fun1)*halfSample
                
                # Store maximum of distribution
                maxDistrib = max(maxDistrib, fun0, fun1)
                
                # Check if in final steps
                if tt < self.simulationTime:
                    tt += min(sampleDt, self.simulationTime - tt)
                    fun0 = fun1
                    fun1 = rateFunc(tt + sampleDt)
                else:
                    break
            
            # Advance integration refinement
            totNum = int(np.round(totNum))
            sampleDt = halfSample
        
        # Divide the nRuns in chunks of 10
        tempLen = 0
        tempRuns = 10; totRuns = 0
        while True:
            if totRuns + tempRuns > nRuns:
                tempRuns = nRuns - totRuns
            if tempRuns == 0:
                break
            
            # Get the random events
            self.__getRuns(fwriteTime, fwriteDist, rateFunc, totNum, maxDistrib,
                           tempRuns)
            
            # Update totRuns
            totRuns += tempRuns
        
        fwriteTime.close()
        fwriteDist.close()
        
        return totNum, timesFile, distFile
    
    def __getRuns(self, fwriteTime, fwriteDist, rateFunc, totNum, maxDistrib,
                  nRuns):
        '''Inner function to generate the files'''
        
        # Get the random events
        
        # For each Myr, produce R events randomly distributed in
        # said Myr:
        times = [[] for x in range(nRuns)]
        distances = [[] for x in range(nRuns)]
        for ii in range(nRuns):
            tempNum = totNum
            while tempNum > 0:
                xx = np.random.random()*self.simulationTime
                yy = np.random.random()*maxDistrib
                
                # If a number is selected
                if yy <= rateFunc(xx):
                    tempNum -= 1
                    
                    # Store time
                    times[ii].append(xx)
                    
                    # Get position
                    xx_pos = self.circumf*np.random.random() - 0.5*self.circumf
                    yy_pos = self.width*np.random.random() - 0.5*self.width
                    zz_pos = np.random.laplace(0, self.hscale)
                    
                    # Calculate the distances from origin
                    distances[ii].append(np.sqrt(xx_pos**2 + yy_pos**2 +
                                                 zz_pos**2))
            
            # Sort the time array
            times[ii].sort()
        
        # Write to temporal files
        for timeArr in times:
            fwriteTime.write(" ".join([str(x) for x in timeArr]))
            fwriteTime.write("\n")
        
        for distArr in distances:
            fwriteDist.write(" ".join([str(x) for x in distArr]))
            fwriteDist.write("\n")
    
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
            try:
                val = float(lnlst[0])
            except ValueError:
                val = lnlst[0]
            except:
                raise
            
            key = lnlst[1]
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
    if inputArgs["rateFunct"] == "linear":
        rate = linearFunction(r0, 5*r0, inputArgs["time"])
        
    elif inputArgs["rateFunct"] == "constant":
        rate = linearFunction(r0, r0, inputArgs["time"])
        
    elif inputArgs["rateFunct"] == "hopkins":
        rate = hopiknsFunction(r0, inputArgs["time"])
        
    elif inputArgs["rateFunct"] == "wanderman":
        rate = wandermanFunction(r0, inputArgs["time"])
        
    elif "cucciati" in inputArgs["rateFunct"]:
        
        # Tabulated values of Cucciati rate function
        lst = [[0.0, 0.2, -1.65],
               [0.2, 0.4, -1.44],
               [0.4, 0.6, -1.34],
               [0.6, 0.8, -1.15],
               [0.8, 1.0, -0.90],
               [1.0, 1.2, -0.85],
               [1.2, 1.7, -0.85],
               [1.7, 2.5, -0.62],
               [2.5, 3.5, -0.86],
               [3.5, 4.5, -1.37]]
        
        if "Interpol" in inputArgs["rateFunct"]:
            rate = cucciatiFunctionInterpol(r0, inputArgs["time"], lst)
        elif "Polyfit" in inputArgs["rateFunct"]:
            rate = cucciatiFunctionPolyfit(r0, inputArgs["time"], lst)
        elif "Step" in inputArgs["rateFunct"]:
            rate = cucciatiFunctionStep(r0, inputArgs["time"], lst)
        else:
            print("Unknown rate function")
            return 1
       
    else:
        print("Unknown rate function")
        return 1
    
    # Set number of runs
    nRuns = int(inputArgs["nRuns"])
    
    # Initialize simulation
    events = EventsObj(hscale = inputArgs["hscale"],
                       circumf = circumf,
                       width = inputArgs["width"],
                       time = inputArgs["time"])
    sizeRun, timesFile, distFile = events.getEvents(rate, nRuns)
    
    # Create filename
    fileName = "../input/Run"
    
    # Sort keys so they always go in same order
    sortedKeys = [key for key in inputArgs]; sortedKeys.sort()
    try:
        fileName += "".join(["_{}_{:.2f}".format(key, inputArgs[key])
                            for key in sortedKeys])
    except ValueError:
        fileName += "".join(["_{}_{}".format(key, inputArgs[key])
                            for key in sortedKeys])
    except:
        raise
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
    step = inputArgs["sampleDt"]
    tArray = np.arange(0, inputArgs["time"] + step, step)
    with open(fileName + "tArray.in", "w") as fwrite:
        fwrite.write("{}\n".format(len(tArray)))
        fwrite.write(" ".join([str(x) for x in tArray]))

if __name__ == "__main__":
    main()
