from Rebuild_Linear_Regression_big import rebuild_dataset_big
from Rebuild_Linear_Regression_big import rebuild_linear_regression_big

from Rebuild_Linear_Regression_small import rebuild_dataset_small
from Rebuild_Linear_Regression_small import rebuild_linear_regression_small

from Simple_Linear_Regression_big import simple_dataset_big
from Simple_Linear_Regression_big import simple_linear_regression_big

from Simple_Linear_Regression_small import simple_dataset_small
from Simple_Linear_Regression_small import simple_linear_regression_small


import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
from scipy.stats import pearsonr as pearson


def make_scatter(measurements, predictions, title):
    """makes a scatter plot"""
    # braucht messungswerte, gefittete werte und vorhergesagte werte; geht nur für einzelne Vorhersagen; z.B. aus Kreuzvalidierung
    fig, ax = plt.subplots(figsize=(10,6))
    fig.canvas.manager.set_window_title(title)
    plt.scatter(measurements, predictions, s=0.3, color="dodgerblue")
    plt.xlabel("measurements")
    plt.ylabel("predictions")

    m, b = np.polyfit(measurements, predictions, 1)

    angle = np.arctan(m)
    angle_degree = np.degrees(angle)

    plt.plot(measurements, m*np.array(measurements)+b, color="midnightblue", label="linear regression")
    plt.plot(measurements, measurements, "k:", alpha=0.15, label=rf"expected line (45°)")

    plt.xlim(0,0.8)
    #plt.ylim(-0.5,2)
    plt.legend(loc="upper left")
    plt.show()


def make_boxplot(data, title):
    """makes a box plot for pearson
 correlations"""
    
    fig, ax = plt.subplots(figsize=(10,6))
    fig.canvas.manager.set_window_title(title)
    fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

    bp = ax.boxplot(data, notch=False, sym="+", vert=True, whis=1.5)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+')
    # Add a horitzontal grid
    ax.yaxis.grid(True, linestyle="-", which="major", color="lightgrey", alpha=0.5)
    ax.set(
        axisbelow=True,
        title=title,
        xlabel="Distribution",
        ylabel="Pearson value",
        )
    # fill boxes with desired color
    box_colors = ['dodgerblue', 'midnightblue']
    num_boxes = len(data)
    medians = np.empty(num_boxes)
    for i in range(num_boxes):
        box = bp["boxes"][i]
        box_x = []
        box_y = []
        for j in range(5):
            box_x.append(box.get_xdata()[j])
            box_y.append(box.get_ydata()[j])
        box_coords = np.column_stack([box_x, box_y])
        # Alternate color
        ax.add_patch(Polygon(box_coords, facecolor=box_colors[i % 2]))

        # median lines back over what we just filled
        med = bp["medians"][i]
        median_x = []
        median_y = []
        for j in range(2):
            median_x.append(med.get_xdata()[j])
            median_y.append(med.get_ydata()[j])
            ax.plot(median_x, median_y, "k")
        medians[i] = median_y[0]
        # samples average
        ax.plot(np.average(med.get_xdata()), np.average(data[i]), color="w", marker="*", markeredgecolor="k")
    # set axes range and axes labels
    ax.set_xlim(0.5, num_boxes + 0.5)
    top = max([max(val) for val in data]) + 0.02
    bottom = min([min(val) for val in data]) - 0.02
    ax.set_ylim(bottom, top)
    ax.set_xticklabels(["simple small", "simple big", "rebuild small", "rebuild big"],
                       rotation=45, fontsize=8)

    # add upper x-axis tick labels with the sample medians
    pos = np.arange(num_boxes) + 1
    upper_labels = [str(round(s, 2)) for s in medians]
    weights = ["bold", "semibold"]
    for tick, label in zip(range(num_boxes), ax.get_xticklabels()):
                           k = tick % 2
                           ax.text(pos[tick], .95, upper_labels[tick],
                                    transform=ax.get_xaxis_transform(),
                                    horizontalalignment="center", size="x-small",
                                    weight=weights[k], color=box_colors[k])
    # add basic legend
    fig.text(0.80, 0.08, "75 training points", backgroundcolor=box_colors[0], color="black", weight="roman", size="x-small")
    fig.text(0.80, 0.045, "3065 training points", backgroundcolor=box_colors[1], color="white", weight="roman", size="x-small")
    fig.text(0.80, 0.013, "*", color="white", backgroundcolor="silver", weight="bold", size="x-small")
    fig.text(0.815, 0.013, " Average Value", color="black", weight="roman", size="x-small")
    
    plt.show()
        
    # all pearson correlations, nach der Kreuzvalidierung und alle zu vergleichenden Werte auf einmal
    return



def make_histogram(helix_percentage, sheet_percentage, coil_percentage, title):
    fig = plt.figure(figsize=(10,6))
    fig.canvas.manager.set_window_title(title)
    all_percentage = np.array(helix_percentage) + np.array(sheet_percentage) + np.array(coil_percentage)
    print("Average: ",sum(all_percentage)/len(all_percentage))
    print("Median: ", np.median(all_percentage))

    
    y, x, _ = plt.hist(all_percentage, 100, density=True, facecolor="midnightblue")
    
    plt.axvline(all_percentage.mean(), color='white', linestyle='dashed', linewidth=1)
    plt.annotate(f"Mean: {round(all_percentage.mean(),2)}", (all_percentage.mean()+0.1, y.max()-0.1))

    plt.axvline(np.median(all_percentage), color='red', linestyle='dashed', linewidth=1)
    plt.annotate(f"Median: {round(np.median(all_percentage),2)}", (np.median(all_percentage)+0.1, y.max()-0.2), color="red")
    
    plt.show()
    

    
def make_barplots(measurements, predictions, title):
    """makes a barplot"""
    fig, ax = plt.subplots(figsize=(10,6))
    fig.canvas.manager.set_window_title(title)
    bar_width = 0.3
    r1 = np.arange(3)
    r2 = [x + bar_width for x in r1]

    sum_all_measurements = sum([sum(measurement) for measurement in measurements])
    sum_single_measurements = [(sum(measurement)/sum_all_measurements)*100 for measurement in measurements]
    sum_all_predictions = sum([sum(prediction) for prediction in predictions])
    sum_single_predictions = [(sum(prediction)/sum_all_predictions)*100 for prediction in predictions]

    print(sum_single_measurements)
    print(sum_single_predictions)

    
    plt.bar(r1, sum_single_measurements, width = bar_width, color = "midnightblue", edgecolor = 'black', label="measurements")
 
    plt.bar(r2, sum_single_predictions, width = bar_width, color = "dodgerblue", edgecolor = 'black', label="predictions")

    plt.xticks([r + (bar_width/2) for r in range(3)], ["Helix", "Sheet", "Coil"])
    plt.ylabel("Secondary Structure in %")
    plt.legend(loc="upper left")
    plt.show()
    
    





def linear_regression_results(n, function_dataset, function_linear_regression):
    helix_percentage, sheet_percentage, coil_percentage, matrix = function_dataset()
    
    h_pearson_ = []
    s_pearson_ = []
    c_pearson_ = []

    h_rmse_ = []
    s_rmse_ = []
    c_rmse_ = []

    h_predictions_ = []
    s_predictions_ = []
    c_predictions_ = []

    h_measurements_ = []
    s_measurements_ = []
    c_measurements_ = []

    h_r2_ = []
    s_r2_ = []
    c_r2_ = []
    
    for i in range(n):
        h_measurements, h_predictions, h_pearson, h_rmse, h_r2 = function_linear_regression(helix_percentage, matrix)
        s_measurements, s_predictions, s_pearson, s_rmse, s_r2 = function_linear_regression(sheet_percentage, matrix)
        c_measurements, c_predictions, c_pearson, c_rmse, c_r2 = function_linear_regression(coil_percentage, matrix)

        
        h_predictions_ += list(h_predictions)
        s_predictions_ += list(s_predictions)
        c_predictions_ += list(c_predictions)

        h_measurements_ += list(h_measurements)
        s_measurements_ += list(s_measurements)
        c_measurements_ += list(c_measurements)
        
        h_pearson_ += list(h_pearson)
        s_pearson_ += list(s_pearson)
        c_pearson_ += list(c_pearson)

        h_rmse_ += list(h_rmse)
        s_rmse_ += list(s_rmse)
        c_rmse_ += list(c_rmse)

        h_r2_ += list(h_r2)
        s_r2_ += list(s_r2)
        c_r2_ += list(c_r2)

    pearson_combi = [h_pearson_, s_pearson_, c_pearson_]
    rmse_combi = [h_rmse_, s_rmse_, c_rmse_]
    predictions_combi = [h_predictions_, s_predictions_, c_predictions_]
    measurements_combi = [h_measurements_, s_measurements_, c_measurements_]
    r2_combi = [h_r2_, s_r2_, c_r2_]
    
    return pearson_combi, rmse_combi, predictions_combi, measurements_combi, r2_combi




if __name__=="__main__":
    n = 10
    
    pearson_rebuild_big, rmse_rebuild_big, predictions_rebuild_big, measurements_rebuild_big, r2_rebuild_big = linear_regression_results(n, rebuild_dataset_big, rebuild_linear_regression_big)

    pearson_rebuild_small, rmse_rebuild_small, predictions_rebuild_small, measurements_rebuild_small, r2_rebuild_small = linear_regression_results(n, rebuild_dataset_small, rebuild_linear_regression_small)

    pearson_simple_big, rmse_simple_big, predictions_simple_big, measurements_simple_big, r2_simple_big = linear_regression_results(n, simple_dataset_big, simple_linear_regression_big)

    pearson_simple_small, rmse_simple_small, predictions_simple_small, measurements_simple_small, r2_simple_small = linear_regression_results(n, simple_dataset_small, simple_linear_regression_small)
        
    print("Helix R^2: ", np.mean(r2_rebuild_big[0]), np.mean(r2_rebuild_small[0]), np.mean(r2_simple_big[0]), np.mean(r2_simple_small[0]))
    print("Sheet R^2: ", np.mean(r2_rebuild_big[1]), np.mean(r2_rebuild_small[1]), np.mean(r2_simple_big[1]), np.mean(r2_simple_small[1]))
    print("Coil R^2: ", np.mean(r2_rebuild_big[0]), np.mean(r2_rebuild_small[0]), np.mean(r2_simple_big[0]), np.mean(r2_simple_small[0]))

    
    plt.style.use("seaborn")
    #-------------------barplots-----------------------#
    make_barplots(measurements_rebuild_big, predictions_rebuild_big, "rebuild big")
    make_barplots(measurements_rebuild_small, predictions_rebuild_small, "rebuild small")
    make_barplots(measurements_simple_big, predictions_simple_big, "simple big")
    make_barplots(measurements_simple_small, predictions_simple_small, "simple small")
   

    #-------------------histogram plots-----------------------#
    # rebuild big
    make_histogram(predictions_rebuild_big[0], predictions_rebuild_big[1], predictions_rebuild_big[2], "rebuild big")
    # rebuild small
    make_histogram(predictions_rebuild_small[0], predictions_rebuild_small[1], predictions_rebuild_small[2], "rebuild small")
    # simple big
    make_histogram(predictions_simple_big[0], predictions_simple_big[1], predictions_simple_big[2], "simple big")
    # simple small
    make_histogram(predictions_simple_small[0], predictions_simple_small[1], predictions_simple_small[2], "simple small")

    
    #-------------------scatter plots-----------------------#
    # helix scatter plots
    make_scatter(measurements_rebuild_big[0], predictions_rebuild_big[0], "rebuild big helix")
    make_scatter(measurements_rebuild_small[0], predictions_rebuild_small[0], "rebuild small helix")
    make_scatter(measurements_simple_big[0], predictions_simple_big[0], "simple big helix")
    make_scatter(measurements_simple_small[0], predictions_simple_small[0], "simple small helix")

    # sheet scatter plots
    make_scatter(measurements_rebuild_big[1], predictions_rebuild_big[1], "rebuild big sheet")
    make_scatter(measurements_rebuild_small[1], predictions_rebuild_small[1], "rebuild small sheet")
    make_scatter(measurements_simple_big[1], predictions_simple_big[1], "simple big sheet")
    make_scatter(measurements_simple_small[1], predictions_simple_small[1], "simple small sheet")

    # coil scatter plots
    make_scatter(measurements_rebuild_big[2], predictions_rebuild_big[2], "rebuild big coil")
    make_scatter(measurements_rebuild_small[2], predictions_rebuild_small[2], "rebuild small coil")
    make_scatter(measurements_simple_big[2], predictions_simple_big[2], "simple big coil")
    make_scatter(measurements_simple_small[2], predictions_simple_small[2], "simple small coil")

    #-------------------pearson boxplots-----------------------#
    # helix boxplots pearson 
    make_boxplot([pearson_simple_small[0], pearson_simple_big[0], pearson_rebuild_small[0], pearson_rebuild_big[0]], "Comparison of Pearson Values Distribution for Linear Models Helix Prediction")
    # sheet boxplots pearson
    make_boxplot([pearson_simple_small[1], pearson_simple_big[1], pearson_rebuild_small[1], pearson_rebuild_big[1]], "Comparison of Pearson Values Distribution for Linear Models Sheet Prediction")
    # coil boxplots pearson
    make_boxplot([pearson_simple_small[2], pearson_simple_big[2], pearson_rebuild_small[2], pearson_rebuild_big[2]], "Comparison of Pearson Values Distribution for Linear Models Coil Prediction")


    #-------------------RMSE boxplots-----------------------#
    # helix boxplots rmse 
    make_boxplot([rmse_simple_small[0], rmse_simple_big[0], rmse_rebuild_small[0], rmse_rebuild_big[0]], "Comparison of RMSE Values Distribution for Linear Models Helix Prediction")
    # sheet boxplots rmse
    make_boxplot([rmse_simple_small[1], rmse_simple_big[1], rmse_rebuild_small[1], rmse_rebuild_big[1]], "Comparison of RMSE Values Distribution for Linear Models Sheet Prediction") 
    # coil boxplots rmse
    make_boxplot([rmse_simple_small[2], rmse_simple_big[2], rmse_rebuild_small[2], rmse_rebuild_big[2]], "Comparison of RMSE Values Distribution for Linear Models Coil Prediction")
                  


   
