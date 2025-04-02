import argparse
import random


def generateTimePoints(horizon, nbmin, nbmax):
    """
    Place a random number of time points between @nbmin and @nbmax
    in the time frame [0, @horizon] 
    """
    nb = random.randint(nbmin, nbmax)
    timepoints = []
    for i in range(nb):
        timepoints.append(random.randint(0, horizon))
    return sorted(set(timepoints))


def generate01Template(nbpoints):
    """
    Generate randomly 0/1 template that determinates     
    """
    points01 = [] 
    for i in range(nbpoints):
        if(i > 0 and points01[-1] != 0):
            rand = random.uniform(0, 1)
            if(rand <= 0.7):
                points01.append(0)
            else:
                points01.append(1)
        else:
            points01.append(1)
    return points01


def generateCandidateValues(nb, min, max):
    """
    Generate @nb candidate values between @min and @max 
    """
    values = []
    for i in range(nb):
        values.append(random.randint(min, max))
    return values


def generateValues(template, vals):
    """
    Generate @nbtimepoints random event values in candidate vals that can be 0 or a value between @min and @max
    """
    values = []
    for i in range(len(template)):
        val = random.choice(vals)
        values.append(template[i] * val)
    return values


def generateBuffer(horizon, minTimepoints, maxTimepoints, nbValues, minValue, maxValue):
    """
    Generate the buffer events, capacities and initial state
    """
    timepoints = generateTimePoints(horizon, minTimepoints, maxTimepoints)
    template = generate01Template(len(timepoints))
    candidates = generateCandidateValues(nbValues, minValue, maxValue)
    values = generateValues(template, candidates) 
    return {
        "timepoints": timepoints,
        "values": values
    }



def computeSumData(buffer, horizon):
    """
    Compute the total production for each buffer
    """
    sum = 0
    cur = buffer["values"][0]
    for i in range(1, len(buffer["timepoints"])):
        sum += cur * (buffer["timepoints"][i] - buffer["timepoints"][i-1])
        cur = buffer["values"][i]
    sum += cur * (horizon - buffer["timepoints"][-1])
    return sum


def computeDownlinkWindows(buffers, nbWindows, horizon, percentVisibility):
    """
    Compute the total data production and generate downlink windows to downlink it 
    """
    totalProduction = 0
    for buffer in buffers:
        totalProduction += computeSumData(buffer, horizon)
    values = []
    for i in range(nbWindows):
        target = totalProduction/nbWindows 
        values.append(round(random.uniform(0.8*target, 1.2*target)))
    visibilitySize = round(horizon / nbWindows * percentVisibility)
    nonvisibilitySize = round(horizon / nbWindows * (1 - percentVisibility))
    windows = []
    currentTime = nonvisibilitySize
    for i in range(nbWindows):
        windows.append(
            {
                "start": currentTime,
                "end": currentTime + visibilitySize,
                "value": round(values[i] / visibilitySize)
            }
        )
        currentTime += visibilitySize + nonvisibilitySize
    return windows


def computeCapacity(buffer, nbWindows, horizon):
    """
    Compute the capacity of a given buffer. 
    We target sumProduction * 2 / nbWindow => 
    The buffer cannot hold data for more than 2 downlink window more or less
    """
    target = computeSumData(buffer, horizon) * 3 / nbWindows
    return round(random.uniform(target * 0.7, target * 1.3)) 



def writeInstance(path, nbBuffers, nbWindows):

    horizon = 10000
    minmaxvalueslb = 10
    minmaxvaluesub = 100 
    timepointslb = 10
    timepointsub = 100
    sizeWindows = 0.85

    buffers = []
    capacities = []
    initialMemState = []

    for i in range(nbBuffers):
        nbValues = random.randint(1, 10)
        minmaxvalues = sorted([random.randint(minmaxvalueslb, minmaxvaluesub) for _ in range(2)])
        buffers.append(generateBuffer(horizon, timepointslb, timepointsub, nbValues, minmaxvalues[0], minmaxvalues[1]))
    
    windows = computeDownlinkWindows(buffers, nbWindows, horizon, sizeWindows)

    for i in range(nbBuffers):
        capacities.append(computeCapacity(buffers[i], nbWindows, horizon))
        initialMemState.append(round(random.uniform(0, 0.2 * capacities[i])))

    # print(buffers)
    # print(windows)
    # print(capacities)
    # print(initialMemState)

    with open(path, "w") as f:
        f.write(str(nbBuffers)+ " buffers\n")
        for i in range(nbBuffers):
            f.write(str(capacities[i]) + " " + str(initialMemState[i]) + "\n")
        f.write(str(nbWindows) + " downlinks\n")
        for i in range(nbWindows):
            f.write(str(i) + " " + str(windows[i]["start"]) + " " + str(windows[i]["end"]) + " " + str(windows[i]["value"]) + "\n")
        for i in range(nbBuffers):
            f.write(str(len(buffers[i]["timepoints"])) + " events for " + str(i) + "\n")
            for e in range(len(buffers[i]["timepoints"])):
                f.write(str(buffers[i]["timepoints"][e]) + " " + str(buffers[i]["values"][e]) + "\n")


def main(args):
    """
    Generate @numinstances for a given number of @buffers and @windows
    Write the instance file as txt in @folder. 
    """

    buffers = [4]
    windows = [4]

    for b in buffers:
            for w in windows:
                for i in range(args.numinstances):
                    writeInstance(args.folder+"/"+str(b)+"-"+str(w)+"-"+str(i)+".txt", b, w)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate instance files for oMDP.")
    parser.add_argument('-n', '--numinstances', help="number of instances to generate",     type=int)
    parser.add_argument("-b", "--buffers",      help="number of buffers for each instance", type=int)
    parser.add_argument("-w", "--windows",      help="number of downlink windows",          type=int)
    parser.add_argument("-f", "--folder",       help="path to store the instances",         type=str)
    args = parser.parse_args()
    main(args)
