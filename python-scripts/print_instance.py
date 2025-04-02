import matplotlib.pyplot as plt
import sys
import matplotlib.ticker as ticker

# lines = []
# downlinks = []
# capacity = []
# events = []

# # Read instance file
# with open("instances/" + sys.argv[1], 'r') as f:
#     lines = f.readlines()

# indicator = "buffers"
# event_index = -1

# for line in lines:
#     if "events" in line:
#         indicator = "events"
#         events.append([])
#         event_index += 1

#     if "downlinks" in line:
#         indicator = "downlink"

#     if indicator == "buffers":
#         if not("buffers" in line):
#             capacity.append(float(line.split()[0]))

#     if indicator == "downlink":
#         if not("downlink" in line):
#             items = line.split()
#             downlinks.append([float(items[1]), float(items[2]), float(items[3])])
    
#     if indicator == "events":
#         if not("events" in line):
#             items = line.split()
#             events[event_index].append([float(items[0]), float(items[1])])


# # Read simulation from logs
# with open("logs/" + sys.argv[1], 'r') as f:
#     data = [list(map(float, line.strip().split())) for line in f]

# # Extract time and values for each item
# times = [d[0] for d in data]
# num_items = len(data[0]) - 1
# inflexion_points = [[] for _ in range(num_items)]
# for d in data:
#     for i in range(num_items):
#         inflexion_points[i].append(d[i+1]*capacity[i]/100)
# max_peaks = []
# for index,i in enumerate(inflexion_points):
#     max_peaks.append([times[i.index(max(i, key=lambda x: x/capacity[index]*100))], max(i, key=lambda x: x/capacity[index]*100)])

# max_peak = [times[max_peaks.index(max(max_peaks, key=lambda x: x[1]))], max(max_peaks, key=lambda x: x[1])[1]]

# # Read solution
# solution = []
# with open("solutions/" + sys.argv[1], 'r') as f:
#     lines = f.readlines()

# for line in lines:
#     items = line.split()
#     priority_window = []
#     for i in items:
#         priority_window.append(int(i))
#     solution.append(priority_window)

# print(downlinks)
# print(capacity)
# print(events)

# num_items = len(events)
# horizon = max([e[-1][0] for e in events])

########
# PLOT #
########

# num_items = 3
# downlinks = [[1, 6, 9], [8, 11, 8]]
# capacity = [10, 10, 10]
# events = [
#     [[0, 5], [1, 2], [4, 4], [6, 0], [7, 2], [9, 0], [10, 1], [11, 0]],
#     [[1, 3], [3, 0], [4, 3], [8, 2], [12, 0]],
#     [[1, 2], [2, 0], [3, 8], [5, 0], [7, 2], [8, 4], [10, 3], [11, 0]]
# ]
# transfer = [
#     [0, 5, 3, 0, 0, 0, 1, 2, 2, 4, 6, 6, 6, 7, 7, 7],
#     [0, 0, 0, 0, 0, 0, 0, 0, 3, 6, 8, 8.5, 7, 4, 6, 6],
#     [0, 0, 0, 0, 0, 1, 6, 3, 3, 5, 1, 0, 0, 0, 0, 0]
# ]
# times = [0, 1, 2, 2+3/4, 3, 4, 5, 6, 7, 8, 9, 9 + 1/4, 10, 11, 12, 13]
# max_peak = [1, 8.5]

# transfer = [
#     [0, 5, 3, 0, 0, 0, 1, 1.5, 1, 1, 3, 5, 5, 6, 6, 6],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 6, 4, 2, 4, 6, 6],
#     [0, 0, 0, 0, 0, 1, 6, 4, 4, 4, 6, 6, 6, 1, 1, 1]
# ]
# times = [0, 1, 2, 2+3/4, 3, 4, 5, 5.5, 6, 7, 8, 9, 10, 11, 12, 13]
# max_peak = [1, 8.5]

# transfer = [
#     [0, 5, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 3, 6, 4, 3, 0, 0, 0, 2, 2],
#     [0, 0, 0, 0, 0, 2, 0.66, 0, 0, 0, 0, 0, 0]
# ]
# times = [0, 1, 1+7/9, 6, 7, 8, 8+2/6, 8+2/4, 8 + 6/6, 9, 11, 12, 13]
# max_peak = [1, 6]

# transfer = [
#     [0, 5, 3, 0, 0, 0, 1, 2, 2, 4, 6, 6, 7, 7, 7],
#     [0, 0, 0, 0, 0, 0, 0, 0, 3, 6, 4, 2, 4, 6, 6],
#     [0, 0, 0, 0, 0, 1, 6, 3, 3, 5, 5, 5, 0, 0, 0]
# ]
# times = [0, 1, 2, 2+3/4, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
# max_peak = [1, 8.5]

# stops = [6, 6, 5.5]

# transfer = [
#     [0, 5, 3, 0, 0, 0, 1, 2, 2, 4],
#     [0, 0, 0, 0, 0, 0, 0, 0, 3, 6],
#     [0, 0, 0, 0, 0, 1, 6, 3, 3, 5]
# ]
# times = [0, 1, 2, 2+3/4, 3, 4, 5, 6, 7, 8]
# max_peak = [1, 6]

# horizon = 13
# priorities = [
#     [1, 1, 1],
#     [2, 1, 1]
# ]

# num_items = 3
# downlinks = [[1, 6, 9], [8, 11, 8]]
# capacity = [10, 10, 10]
# events = [
#     [[0, 5], [1, 2], [6, 0], [7, 2], [9, 0], [10, 3], [11, 0]],
#     [[1, 3], [3, 0], [8, 2], [12, 0]],
#     [[1, 2], [2, 0], [3, 3], [5, 0], [7, 2], [11, 0]]
# ]
# transfer = [
#     [0, 5, 7, 9, 10, 10, 10, 10, 10, 10, 10, 10],
#     [0, 0, 3, 6, 6,   6,  6,  6,  6,  8, 10, 10],
#     [0, 0, 2, 2, 3.5, 5,  8,  8, 10, 10, 10, 10]
# ]
# overflows = [3.5, 10, 8]
# times = [0, 1, 2, 3, 3.5, 4, 5, 7, 8, 9, 10, 13]
# max_peak = [1, 8.5]
# horizon = 13
# priorities = [
#     [1, 1, 1],
#     [3, 2, 1]
# ]

# num_items = 3
# downlinks = [[1, 6, 9]]
# capacity = [10, 10, 10]
# events = [
#     [[0, 5], [1, 2], [4, 4], [6, 0]],
#     [[1, 3], [3, 0], [4, 3], [6, 0]],
#     [[1, 2], [2, 0], [3, 8], [5, 0]]
# ]
# transfer = [
#     [0, 5, 0, 0],
#     [3, 3, 3+3*(7/9), 0, 0],
#     [2, 2, 0, 0, 1, 4, 0, 0]
# ]
# # times = [0, 1, 1+7/9, 6]
# # times = [0, 1, 1+7/9, 1+7/9+(3+3*(7/9))/9, 6]
# times = [0, 1, 1+4/4.5, 3, 4, 5, 5 + 4/5, 6]
# horizon = 6
# max_peak = [2, 4]

# num_items = 2
# downlinks = [[1, 3, 3], [4, 8, 4]]
# capacity = [10, 10]
# events = [
#     [[0, 1], [1, 2], [3, 0], [4, 6], [5, 0], [6, 4], [8, 0]],
#     [[1, 3], [2, 0], [3, 2], [4, 1], [6, 0]]
# ]

# transfer = [
#     [0, 1, 0, 0, 0, 4, 2, 4, 4],
#     [0, 0, 3, 2, 4, 3, 2, 0, 0]
# ]

# times = [0, 1, 2, 3, 4, 5, 6, 7, 8]
# max_peak = [0, 4]
# horizon = 8

# priorities = [
#     [1,2],
#     [1,1]
# ]

# num_items = 3
# downlinks = [[1, 7, 3]]
# capacity = [10, 10, 10]
# events = [
#     [[0, 5], [1, 3], [3, 0]],
#     [[0, 5], [1, 0], [3, 3], [5, 0]],
#     [[0, 5], [1, 0],         [5, 3], [7, 0]]
# ]

# transfer = [
#     [0, 5, 8, 5, 5, 5],
#     [0, 5, 2, 5, 5, 5],
#     [0, 5, 5, 5, 5, 5]
# ]

# transfer = [
#     [0, 5, 5, 5, 5, 5],
#     [0, 5, 5, 5, 5, 5],
#     [0, 5, 5, 5, 5, 5]
# ]
# stops = [3, 5, 7]   

# times = [0, 1, 3, 5, 7, 8]
# # max_peak = [0, 5]
# horizon = 8

# priorities = [
#     [1,2,3]
# ]

# num_items = 2
# downlinks = [[1, 5, 2]]
# capacity = [10, 10]
# events = [
#     [[0, 2],         [3, 0]],
#     [[0, 2], [1, 0],         [4, 2], [5, 0]]
# ]

# transfer = [
#     [0, 2, 3, 4, 2, 0, 0],
#     [0, 2, 1, 0, 0, 2, 2]
# ]

# times = [0, 1, 2, 3, 4, 5, 6]
# horizon = 6

# priorities = [
#     [1,1]
# ]
# stops = [5, 3]   


##########################
# Bandwidth loss example #
##########################

# num_items = 2
# downlinks = [[0, 4, 1]]
# capacity = [2.2, 2.2]
# events = [
#     [[0, 2], [1, 0],                   [3, 0.5], [4, 0]],                 
#     [                 [1.5, 1], [2, 0]]
# ]

# transfer = [
#     [0, 1, 0.5, 0.5, 0.5, 1],
#     [1, 1, 1,   1,   0,   0]
# ]

# times = [0, 1, 1.5, 2, 3, 4]
# horizon = 4

# priorities = [
#     [1,2]
# ]
# stops = [1.5, 100] 



##########################
# Single Window NON OPTI #
##########################

# SW1
num_items = 3
downlinks = [[0, 6, 1]]
capacity = [2.2, 2.2, 2.2]
real_capacity = ["1", "7/6", "5/12"]
events = [ 
    [[1, 1.34], [2, 0], [4, 2], [5,0]],                 
    [[0, 1],    [2, 0]],
    [[0, 1/2],  [2, 0], [5, 17/12], [6, 0]]
]
eventstxt = [ 
    [[1, "4/3+e"], [2, 0], [4,2], [5,0]],                 
    [[0, 1],    [2, 0]],
    [[0, "1/2"],  [2, 0], [5, "17/12"], [6, 0]]
]
# overflowMarker = [0, 2, 1.01]

# transfer = [
#     [0, 0, 1.01],
#     [0, 1/2, 7/6],
#     [0, 0, 1/6]
# ]

# times = [0,1,2]
horizon = 6

# priorities = [
#     [3,3,3]
# ]
# stops = [100, 100, 100] 


# # SW2
# num_items = 3
# downlinks = [[0, 6, 1]]
# capacity = [2.2, 2.2, 2.2]
# real_capacity = ["1", "7/6", "5/12"]
# events = [ 
#     [[1, 1.34], [2, 0], [4, 2], [5,0]],                 
#     [[0, 1],    [2, 0]],
#     [[0, 1/2],  [2, 0], [5, 17/12], [6, 0]]
# ]
# eventstxt = [ 
#     [[1, "4/3+e"], [2, 0], [4,2], [5,0]],                 
#     [[0, 1],    [2, 0]],
#     [[0, "1/2"],  [2, 0], [5, "17/12"], [6, 0]]
# ]
# overflowMarker = [2, 1+1/3, 5/12]

# transfer = [
#     [0, 0, 0, 1/12],
#     [0, 0, 1/6, 3/6],
#     [0, 5/12, 4/12, 5/12]
# ]

# times = [0,5/6,1,1+1/3]
# horizon = 6

# priorities = [
#     [2,2,3]
# ]
# stops = [100, 5/6, 100] 


##########################
# Illustrative Example   #
##########################

num_items = 3
downlinks = [[1, 4, 2], [6, 9, 3]]
capacity = [4, 4, 4]
real_capacity = [4, 4, 4]
events = [ 
    [[0, 2], [2, 0], [4, 2], [5,0], [6, 3], [7, 0]],                 
    [[1, 1], [4, 0], [7, 3], [8,0]], 
    [[1, 1], [4, 0], [8, 3], [9,0]]
]

eventstxt = [ 
    [[0, 2], [2, 0], [4, 2], [5,0], [6, 3], [7, 0]],                 
    [[1, 1], [5, 0], [7, 3], [8,0]], 
    [[1, 1], [5, 0], [8, 3], [9,0]]
]

transfer = [
    [0, 2, 2, 0, 0, 2, 2, 2, 2, 2, 2],
    [0, 0, 1, 2, 2, 2, 2, 2, 2, 2, 2],
    [0, 0, 1, 2, 2, 2, 2, 2, 2, 2, 2],
]

times = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
horizon = 10

priorities = [
    [1,2,2],
    [1,2,3]
]
stops = [7, 8, 100] 


##########################
# Illustrative Example   #
##########################

# num_items = 3
# downlinks = [[1, 7, 3]]
# capacity = [6, 6, 6]
# real_capacity = ["10", "10", "10"]
# events = [ 
#     [[0, 5], [1, 0], [1, 3], [3, 0]],                 
#     [[0, 5], [1, 0], [3, 3], [5, 0]],
#     [[0, 5], [1, 0], [5, 3], [7, 0]]
# ]
# eventstxt = [ 
#     [[0, 5], [1, 0], [1, 3], [3, 0]],                 
#     [[0, 5], [1, 0], [3, 3], [5, 0]],
#     [[0, 5], [1, 0], [5, 3], [7, 0]]
# ]
# # overflowMarker = [0, 2, 1.01]

# transfer = [
#     [0, 5, 5],
#     [0, 5, 5],
#     [0, 5, 5]
# ]

# priorities = [[1,2,3]]

# stops = [3, 5, 100] 
# times = [0,1,7]
# horizon = 7


fontSize = 18

# for step in range(len(transfer[1]), len(transfer[1])+1):
for step in range(1):

    # Specify figure dimensions in inches
    fig_size = (16, 6)  # Adjust dimensions as needed

    # Create plot with subplots for each item
    fig, axes = plt.subplots(nrows=num_items, ncols=1, sharex=True, sharey=True, figsize=fig_size)

    # Set y-axis tick interval and force them to integers
    for ax in axes:
        ax.yaxis.set_major_locator(ticker.MultipleLocator(base=5.0))  # set the interval to 10
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1.0))  # set the interval to 1
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{int(x)}'))  # force to integers
        ax.tick_params(axis='both', which='major', labelsize=fontSize-4)

    # Plot each item in a separate subplot
    for i in range(num_items):
        # if i == max_peak[0]:
        #     axes[i].axhline(y=max_peak[1], linestyle='dashed', color='red')
            # axes[i].text(0.2, max_peak[1], f"rmax = {int(max_peak[1]*100/capacity[i])}%", color="red", fontsize=14, verticalalignment='bottom')
            # axes[i].text(1.1, max_peak[1], f"rmax = {int(max_peak[1]*100/capacity[i])}%", color="red", fontsize=14, verticalalignment='bottom')
        axes[i].plot(times, transfer[i], color="red", linewidth=2)
        axes[i].plot(stops[i], 0, marker="*", markersize=14, clip_on=False, alpha=1)
        axes[i].set_ylabel(r"$C_{}={}$".format(i+1, real_capacity[i]), fontsize=fontSize, labelpad=30)
        axes[i].yaxis.label.set_rotation(0)
        axes[i].set_ylim([0, capacity[i]])
        axes[i].set_xlim([0, horizon])

    colors = ["green", "gold", "blue", "red", "pink", "silver"]

    for i in range(len(events)):
        j=0
        while j < len(events[i]):
            start = events[i][j][0]
            end = events[i][j + 1][0] if j + 1 < len(events[i]) else start + 1
            height = events[i][j][1]
            label = str(eventstxt[i][j][1])

            axes[i].fill_between([events[i][j][0], events[i][j+1][0]], 0, events[i][j][1], color=colors[i], edgecolor='white')
            axes[i].text((start + end) / 2, height / 2, label, ha='center', va='center', fontsize=fontSize, color='white')
            if events[i][j+1][1] == 0:
                j+=2
            else:
                j+=1

    # Add shaded regions and labels to all subplots
    for index,d in enumerate(downlinks):
        start, end, value = d[0], d[1], d[2]
        for i in range(num_items):
            axes[i].fill_betweenx([0, axes[i].get_ylim()[1]], start, end, alpha=0.2, color='green')
           #  axes[i].text(start + (end-start)/2, axes[0].get_ylim()[1]*0.6, r'$p_{{{}}}={}$ '.format(i+1, priorities[index][i]), ha='center', fontsize=fontSize)
            axes[i].text(start + (end-start)/2, axes[0].get_ylim()[1]*0.85, r'$p_{{{}, {}}}={}$ '.format(i+1, index+1, priorities[index][i]), ha='center', fontsize=14)    
        axes[0].text(start + (end-start)/2, axes[0].get_ylim()[1]*1.1, r'$\delta_{}={}$'.format(index+1, round(d[2])), ha='center', fontsize=fontSize+2)
        # axes[0].text(start + (end-start)/2, axes[0].get_ylim()[1]*1.05, r'$\delta_{}=0$'.format(index+1), ha='center', color="red", fontsize=14)

    # Add new domains for simulation optimistic
    # axes[0].text(start + (end-start)/2, axes[0].get_ylim()[1]*0.85, r'$p_{{i,1}}=3$ ', ha='center', color="red", fontsize=14) 
    # axes[1].text(start + (end-start)/2, axes[0].get_ylim()[1]*0.85, r'$p_{{i,2}}=2$ ', ha='center', color="red", fontsize=14) 
    # axes[2].text(start + (end-start)/2, axes[0].get_ylim()[1]*0.85, r'$p_{{i,3}}=1$ ', ha='center', color="red", fontsize=14) 

    # axes[0].text(0.5, axes[0].get_ylim()[1]*0.85, r'$mem_{{i,1}}=0$ ', ha='center', color="red", fontsize=14) 
    # axes[1].text(0.5, axes[0].get_ylim()[1]*0.85, r'$mem_{{i,2}}=1$ ', ha='center', color="red", fontsize=14) 
    # axes[2].text(0.5, axes[0].get_ylim()[1]*0.85, r'$mem_{{i,3}}=2$ ', ha='center', color="red", fontsize=14) 

    # axes[0].text(start + (end-start)/2, axes[0].get_ylim()[1]*0.85, r'$p_{{i,1}}=[1]$ ', ha='center', color="black", fontsize=13) 
    # axes[1].text(start + (end-start)/2, axes[0].get_ylim()[1]*0.85, r'$p_{{i,2}}=[2,3]$ ', ha='center', color="black", fontsize=13) 
    # axes[2].text(start + (end-start)/2, axes[0].get_ylim()[1]*0.85, r'$p_{{i,3}}=[1,2,3]$ ', ha='center', color="black", fontsize=13) 

    # axes[0].text(0.5, axes[0].get_ylim()[1]*0.85, r'$mem_{{i,1}}=[0, 10]$ ', ha='center', color="black", fontsize=13) 
    # axes[1].text(0.5, axes[0].get_ylim()[1]*0.85, r'$mem_{{i,2}}=[3, 10]$ ', ha='center', color="black", fontsize=13) 
    # axes[2].text(0.5, axes[0].get_ylim()[1]*0.85, r'$mem_{{i,3}}=[2, 10]$ ', ha='center', color="black", fontsize=13) 
    
    # axes[0].text(end -.5, axes[0].get_ylim()[1]*0.85, r'$mem_{{i+1,1}}=[0, 10]$', ha='center', color="black", fontsize=13) 
    # axes[1].text(end -.5, axes[0].get_ylim()[1]*0.85, r'$mem_{{i+1,2}}=[0, 10]$', ha='center', color="black", fontsize=13) 
    # axes[2].text(end -.5, axes[0].get_ylim()[1]*0.85, r'$mem_{{i+1,3}}=[0, 10]$', ha='center', color="black", fontsize=13) 

    # axes[0].text(start + (end-start)/2, axes[0].get_ylim()[1]*0.85, r'$p_{{i,1}}=1$ ', ha='center', color="red", fontsize=14) 
    # axes[1].text(start + (end-start)/2, axes[0].get_ylim()[1]*0.85, r'$p_{{i,2}}=3$ ', ha='center', color="red", fontsize=14) 
    # axes[2].text(start + (end-start)/2, axes[0].get_ylim()[1]*0.85, r'$p_{{i,3}}=1$ ', ha='center', color="red", fontsize=14) 

    # axes[0].text(0.5, axes[0].get_ylim()[1]*0.85, r'$mem_{{i,1}}=0$ ', ha='center', color="red", fontsize=14) 
    # axes[1].text(0.5, axes[0].get_ylim()[1]*0.85, r'$mem_{{i,2}}=3$ ', ha='center', color="red", fontsize=14) 
    # axes[2].text(0.5, axes[0].get_ylim()[1]*0.85, r'$mem_{{i,3}}=2$ ', ha='center', color="red", fontsize=14) 
    
    # axes[2].text(end -.5, axes[2].get_ylim()[1]*0.85, r'$mem_{{i+1,3}}=[0, 10]$', ha='center', color="red", fontsize=14) 
    # axes[2].text(end -.5, axes[2].get_ylim()[1]*0.55, r'$r_{{i,3}}=[40\%, 100\%]$', ha='center', color="red", fontsize=14)

    # Set x-axis label for last subplot
    axes[num_items-1].set_xlabel('Time', fontsize=fontSize)

    # plt.show()
    # plt.savefig("sw1.pdf", format="pdf")
    # plt.savefig("sw2.pdf", format="pdf")
    # plt.savefig("sw3.pdf", format="pdf")
    # plt.savefig("sw.pdf", format="pdf")
    # plt.savefig("swnonopti.pdf", format="pdf")
    # plt.savefig("swopti.pdf", format="pdf")
    plt.savefig("plots/illustrative.pdf", format="pdf")