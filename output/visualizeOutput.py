import sys, random
import matplotlib.pyplot as plt
import numpy as np

def main():
    if len(sys.argv) < 2:
        print("Use: python3 {} <outputFile> [outputFile2]".format(sys.argv[0]))
        s = "When the second output file is included, the ratio between the"
        s += " two solutions is shown."
        print(s)
        return 1
    
    # Get input file
    inputFile = sys.argv[1]
    if len(sys.argv) > 2:
        inputFile2 = sys.argv[2]
        isSynchronous = input("Synchronous division? (y/n): ")
        isSynchronous = True if isSynchronous == "y" else False
    else:
        inputFile2 = None
        isSynchronous = None
    
    # Read the runs
    with open(inputFile, "r") as fread:
        tArray = [float(x) for x in fread.readline().split()]
        
        if inputFile2 is not None:
            fread2 = open(inputFile2, "r")
            if isSynchronous:
                fread2.readline() # Remove the tArray from this one
        
        # Read each run
        results = []
        for line in fread:
            if inputFile2 is not None:
                line2 = fread2.readline().split()
                results.append([float(x)/(1e-200 + float(y)) for x, y in
                                    zip(line.split(), line2)])
            else:
                results.append([float(x) for x in line.split()])
        
        if inputFile2 is not None:
            fread2.close()
            # If not synchronous, remove first result (divided by tArray!)
            if not isSynchronous:
                results.pop(1)
    
    # Take one random example for the plot
    oneExample = results[np.random.choice(len(results))]
    
    # Transpose the results for calculations
    results = np.transpose(results)
    
    # Calculate statistics
    median = [np.median(x) for x in results]
    sig1p = [np.percentile(x, 50 + 34) for x in results]
    sig1n = [np.percentile(x, 50 - 34) for x in results]
    sig2p = [np.percentile(x, 50 + 47.5) for x in results]
    sig2n = [np.percentile(x, 50 - 47.5) for x in results]
    maxNum = [max(x) for x in results]
    minNum = [min(x) for x in results]
    
    # Choose RGB colors for:
    # 1 sigma
    innerColor = (0.1, 0.1, 0.1)
    
    # 2 Sigma
    outerColor = (0.4, 0.4, 0.4)
    
    # Extremes
    gray = (0.8, 0.8, 0.8)
    
    # Get minimum time for plotting
    lowTime = 2000
    for ii in range(len(tArray)):
        if tArray[ii] >= lowTime:
            break
    
    # Plot outside in
    plt.fill_between(tArray[ii:], minNum[ii:], maxNum[ii:], color = gray,
                     label = "Full")
    plt.fill_between(tArray[ii:], sig2n[ii:], sig2p[ii:], color = outerColor,
                     label = "95%")
    plt.fill_between(tArray[ii:], sig1n[ii:], sig1p[ii:], color = innerColor,
                     label = "68%")
    plt.plot(tArray[ii:], median[ii:], "y-", label = "Median")
    plt.plot(tArray[ii:], oneExample[ii:], "b-", label = "One run")
    
    # Plot log only if there is more than one order of magnitude
    if min(sig2n[ii:]) < 0.1*max(sig2p[ii:]):
        plt.yscale("log")
    
    plt.xlabel("time [Myr]")
    plt.ylabel("Mass (arbitrary)")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
