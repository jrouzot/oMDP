import os
import re
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict


def getStats(folder):
    for root, _, files in os.walk(folder):
        nb_files = 0
        branches = 0
        failures = 0
        optimal = 0
        time = 0
        for file in files:
            filepath = os.path.join(folder, file)
            nb_files += 1
            with open(filepath, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if line.startswith("branches"):
                        branches += int(line.strip().split("=")[1])
                    if line.startswith("failures"):
                        failures += int(line.strip().split("=")[1])
                    if line.startswith("real_time"):
                        time += int(line.strip().split("=")[1])
                    if line.startswith("status"):
                        if line.strip().split("=")[1] == "OPTIMAL":
                            optimal += 1

        print(optimal)
        print(branches / nb_files / 1000000) 
        print(failures / nb_files / 1000000)
        print(time / nb_files / 1000 )


def makeResultDict(method, folders):
    """
    Make a dict from the raw result files
    """

    # Data storage
    data = defaultdict(lambda: {"total": 0, "optimal": 0, "obj_values": [], "lb_values": []})

    # Pattern matching
    pattern = re.compile(r"(buffers|windows|obj|real_obj|status|original_lb)=([\w\d\.]+)")

    # Process files
    for folder in folders:
        for root, _, files in os.walk(folder):
            for file in files:

                filepath = os.path.join(folder, file)
            
                if not os.path.isfile(filepath):
                    continue  # Skip non-files

                with open(filepath, "r") as f:
                    content = f.read()
            
                # Extract key-value pairs
                matches = dict(pattern.findall(content))
            
                # Convert relevant values
                buffers = int(matches.get("buffers", -1))
                windows = int(matches.get("windows", -1))
                if method.startswith("ours"):
                    windows -= 1
                    obj = float(matches.get("real_obj", np.nan))
                else:
                    obj = float(matches.get("obj", np.nan))
                status = matches.get("status", "")
                original_lb = float(matches.get("original_lb", np.nan))     # Extract original_lb

                # Skip invalid entries
                if buffers == -1 or windows == -1 or np.isnan(obj) or np.isnan(original_lb):
                    continue

                # Update statistics
                key = (buffers, windows)
                data[key]["total"] += 1
                if method.startswith("ours"):
                    data[key]["obj_values"].append(obj / 100)
                    data[key]["lb_values"].append(original_lb / 1000)
                    if status == "OPTIMAL":
                        data[key]["optimal"] += 1
                else:
                    data[key]["obj_values"].append(obj)
                    data[key]["lb_values"].append(original_lb)
                    if obj <= original_lb + 0.0001:
                        data[key]["optimal"] += 1

    return data



def makeTable(data):
    """
    Make a proportion optimal table and a mean gap table
    """

    # Prepare data for tables
    table_data = []
    gap_data = []

    for (buffers, windows), stats in data.items():
        proportion_optimal = stats["optimal"] / stats["total"] if stats["total"] > 0 else 0
        avg_gap = np.mean([(obj - lb) / lb for obj, lb in zip(stats["obj_values"], stats["lb_values"]) if lb > 0])
        
        table_data.append([buffers, windows, proportion_optimal])
        gap_data.append([buffers, windows, avg_gap])

    # Convert to Pandas DataFrame for better display
    df_optimal = pd.DataFrame(table_data, columns=["Buffers", "Windows", "Proportion Optimal"])
    df_gap = pd.DataFrame(gap_data, columns=["Buffers", "Windows", "Average Optimality Gap"])

    # Print the tables
    print("\nProportion of Optimal Solutions:")
    print(df_optimal.pivot(index="Windows", columns="Buffers", values="Proportion Optimal"))

    print("\nAverage Optimality Gap:")
    print(df_gap.pivot(index="Windows", columns="Buffers", values="Average Optimality Gap"))



def makeOptimalityBarPlot(data1, data2, data3):
    """
    Make a plot bar according to the data sets
    """

    # Prepare data for tables
    table_data_1 = []
    table_data_2 = []
    # table_data_3 = []

    for (buffers, windows), stats in data1.items():
        proportion_optimal = stats["optimal"] / stats["total"] if stats["total"] > 0 else 0
        table_data_1.append([buffers, windows, proportion_optimal])

    for (buffers, windows), stats in data2.items():
        proportion_optimal = stats["optimal"] / stats["total"] if stats["total"] > 0 else 0
        table_data_2.append([buffers, windows, proportion_optimal])

    # for (buffers, windows), stats in data3.items():
    #     proportion_optimal = stats["optimal"] / stats["total"] if stats["total"] > 0 else 0
    #     table_data_3.append([buffers, windows, proportion_optimal])


    # Convert to Pandas DataFrame for better display
    df_optimal_1 = pd.DataFrame(table_data_1, columns=["Buffers", "Windows", "Proportion Optimal"])
    df_optimal_2 = pd.DataFrame(table_data_2, columns=["Buffers", "Windows", "Proportion Optimal"])
    # df_optimal_3 = pd.DataFrame(table_data_3, columns=["Buffers", "Windows", "Proportion Optimal"])

    # Get sorted unique values
    buffer_values = sorted(set(df_optimal_1["Buffers"]))
    window_values_1 = sorted(set(df_optimal_1["Windows"]))
    window_values_2 = sorted(set(df_optimal_2["Windows"]))
    # window_values_3 = sorted(set(df_optimal_3["Windows"]))

    # Bar width
    bar_width = 0.15  
    x_indexes = np.arange(len(window_values_1))

    # Plot results
    plt.figure(figsize=(10, 6))
    # colors = ["red", "blue", "green", "gold"]  # Assign colors for buffers

    for i, buffers in enumerate(buffer_values):
        y_values = df_optimal_1[df_optimal_1["Buffers"] == buffers].set_index("Windows").reindex(window_values_1)["Proportion Optimal"].fillna(0)
        plt.bar(x_indexes - (len(buffer_values)/2 - i)*bar_width, y_values, width=bar_width, label=f"{buffers} Buffers")

    for i, buffers in enumerate(buffer_values):
        y_values = df_optimal_2[df_optimal_2["Buffers"] == buffers].set_index("Windows").reindex(window_values_2)["Proportion Optimal"].fillna(0)
        plt.bar(x_indexes - (len(buffer_values)/2 - i)*bar_width, y_values, width=bar_width, color="black", alpha=0.3)

    # for i, buffers in enumerate(buffer_values):
    #     y_values = df_optimal_3[df_optimal_3["Buffers"] == buffers].set_index("Windows").reindex(window_values_3)["Proportion Optimal"].fillna(0)
    #     plt.bar(x_indexes - (len(buffer_values)/2 - i)*bar_width, y_values, width=bar_width, color="black", alpha=0.3)


    plt.rcParams.update({'font.size': 12})
    plt.xlabel("Number of Downlink Windows", fontsize=13)
    plt.ylabel("Proportion of Optimal Solutions", fontsize=13)
    plt.xticks(x_indexes - bar_width/2, window_values_1, fontsize=11)  # Center labels
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], fontsize=11)
    plt.legend()
    plt.grid(axis="y", linestyle="--", alpha=0.7)

    # Save the plot
    # plt.show()
    # plt.savefig("../plots/reduced.pdf", format="pdf")
    # plt.savefig("../plots/new-instances/bar.pdf", format="pdf")
    plt.savefig("../plots/new-propagator/bar.pdf", format="pdf")
    # plt.close()


def cactus(folders, methods):
    """
    Make cactus plot for random, firstfail and downlink count methods for the given instance
    """
    results = {}
    instances = ["MTP011", "MTP012", "MTP013", "MTP014"]
    for ins in instances:
        results[ins] = {}
    for i,folder in enumerate(folders):
        for root, _, files in os.walk(folder):
            for file in files:
                instance = "MTP01"+file.split(".")[0]
                with open(folder+"/"+file, 'r') as f:
                    results[instance][methods[i]] = {}
                    results[instance][methods[i]]["obj"] = []
                    results[instance][methods[i]]["time"] = []
                    for line in f.readlines():
                        if line.startswith("obj"):
                            results[instance][methods[i]]["obj"].append(int(line.split("=")[1]) / 10)
                        if line.startswith("time"):
                            results[instance][methods[i]]["time"].append(int(line.split("=")[1]) / 60000)
                        if line.startswith("peak-ub"):
                            results[instance][methods[i]]["ub"] = int(line.split("=")[1]) / 10
    

    plt.rcParams.update({'font.size': 14})

    markers = ['o', 's', '^']

    for instance in instances:
        for i,method in enumerate(methods):
            x = results[instance][method]["time"]
            y = results[instance][method]["obj"]
            # Insert starting point*
            x.insert(0, 0)
            y.insert(0, results[instance][method]["ub"])
            plt.plot(x, y, marker=markers[i], markersize=8, linewidth=2.2, label=method)
        plt.grid(True)
        plt.title(instance)
        plt.xlabel("time (min)")
        plt.ylabel("rmax (%)")
        plt.legend()
        plt.subplots_adjust(bottom=0.15)
        plt.savefig("../plots/"+instance+".pdf", format="pdf")
        plt.close()

        # plt.legend()
        # plt.grid(axis="y", linestyle="--", alpha=0.7)
        # plt.subplots_adjust(bottom=0.1) # or whatever

        # # Save the plot
        # # plt.show()
        # plt.savefig("../plots/new-propagator/small-"+str(buffers)+"buffers.pdf", format="pdf")
        # # plt.savefig("../plots/old-instances/"+str(buffers)+"buffers.pdf", format="pdf")
        # plt.close()


# def makeGapBarPlot(data1, data2):

#     gap_data_1 = []
#     gap_data_2 = []

#     for (buffers, windows), stats in data1.items():
#         avg_gap = np.mean([(obj - lb) / lb for obj, lb in zip(stats["obj_values"], stats["lb_values"]) if lb > 0])
#         gap_data_1.append([buffers, windows, avg_gap])
    
#     for (buffers, windows), stats in data2.items():
#         avg_gap = np.mean([(obj - lb) / lb for obj, lb in zip(stats["obj_values"], stats["lb_values"]) if lb > 0])
#         gap_data_2.append([buffers, windows, avg_gap])

#     df_gap_1 = pd.DataFrame(gap_data_1, columns=["Buffers", "Windows", "Average Optimality Gap"])
#     df_gap_2 = pd.DataFrame(gap_data_2, columns=["Buffers", "Windows", "Average Optimality Gap"])


def makeObjectiveValuePlot(datas, names):
    """
    Make a line plot for the objective value for data1, data2 and data3, data4
    """

    objs = []

    for i,data in enumerate(datas):
        objs.append([])
        for (buffers, windows), stats in data.items():
            avg_obj = np.mean(stats["obj_values"])
            objs[i].append([buffers, windows, avg_obj])

    dfs = []
    for obj in objs:
        dfs.append(pd.DataFrame(obj, columns=["Buffers", "Windows", "Objective Value"]))

    print(dfs)

    buffer_values = sorted(set(dfs[0]["Buffers"]))
    window_values = sorted(set(dfs[0]["Windows"]))
    x_indexes = np.arange(len(window_values))

    plt.rcParams.update({'font.size': 12})

    markers = ['o', 's', '^']

    for i, buffers in enumerate(buffer_values):
        # Plot results
        plt.figure(figsize=(10, 6))
        for i,df in enumerate(dfs):
            y_values = df[df["Buffers"] == buffers].set_index("Windows").reindex(window_values)["Objective Value"].fillna(0)
            plt.plot(x_indexes, y_values, marker=markers[i], label=f"{buffers} Buffers, {names[i]}")

        plt.xlabel("Number of Downlink Windows", fontsize=12)
        plt.ylabel("Mean Objective Value", fontsize=12)
        # plt.title("Mean Objective Value by Windows")
        plt.xticks(x_indexes, [w for w in window_values], fontsize=11)  # Center labels
        plt.yticks(fontsize=11)
        plt.legend()
        plt.grid(axis="y", linestyle="--", alpha=0.7)
        plt.subplots_adjust(bottom=0.1) # or whatever

        # Save the plot
        # plt.show()
        plt.savefig("../plots/large-"+str(buffers)+"buffers.pdf", format="pdf")
        # plt.savefig("../plots/old-instances/"+str(buffers)+"buffers.pdf", format="pdf")
        plt.close()


# data1 = makeResultDict("ours", ["../results/Ours/small-1h"])
# data2 = makeResultDict("repair", ["../results/RepairDescent/small-1h"])
# data3 = makeResultDict("iterative", ["../results/IterativeLeveling/small-1h"])

data1 = makeResultDict("ours", ["../results/Ours/large-1h"])
data2 = makeResultDict("repair", ["../results/RepairDescent/large-1h"])
data3 = makeResultDict("iterative", ["../results/IterativeLeveling/large-1h"])

makeObjectiveValuePlot([data1, data2, data3], ["CP", "Repair", "Iterative"])

cactus(["../results/Ours/mtp", "../results/Ours/mtp-firstfail", "../results/Ours/mtp-random"], ["downlink count", "min dom" ,"random"])