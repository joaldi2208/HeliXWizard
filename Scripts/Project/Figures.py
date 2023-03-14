import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
import seaborn as sns
from sklearn.metrics import r2_score

def make_one_scatter(measurements, predictions_, title):
    """makes a scatter plot"""
    fig, ax = plt.subplots(figsize=(10,6))
    fig.canvas.manager.set_window_title(title)
    
    plt.gca().set_xlim(0,1)
    plt.gca().set_ylim(0,1)
    plt.grid()

    line = [x/10 for x in range(10)]
    plt.plot(line, line)
    names = ["rebuild big", "simple big", "rebuild small", "simple small"]
    for name, predictions in zip(names, predictions_):
        sns.regplot(measurements, predictions, scatter=False, label=name)

    plt.xlabel("measurements")
    plt.ylabel("predictions")

    plt.legend(loc="upper left")
    plt.savefig("../../Drafts/Project-Thesis/img/" + "Pearson_scatter_" + title)
    plt.show()
    plt.close()


    
def make_scatter(measurements, predictions, title):
    """makes a scatter plot"""
    # braucht messungswerte, gefittete werte und vorhergesagte werte; geht nur für einzelne Vorhersagen; z.B. aus Kreuzvalidierung
    fig, ax = plt.subplots(figsize=(10,6))
    fig.canvas.manager.set_window_title(title)
    plt.scatter(measurements, predictions, s=8, color="dodgerblue")
    plt.xlabel("measurements", fontsize=20)
    plt.ylabel("predictions", fontsize=20)

    m, b = np.polyfit(measurements, predictions, 1)

    print("r2_score: ", title, ": ",r2_score(measurements, predictions))
    plt.plot(measurements, m*np.array(measurements)+b, color="midnightblue")

    
    #plt.xlim(0,0.8)
    #plt.ylim(-0.5,2)
    #plt.legend(loc="upper left")
    plt.savefig("../../Drafts/Project-Thesis/img/" + "Pearson_scatter_" + title)
    #plt.show()
    plt.close()
    


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

    plt.savefig("../../Drafts/Project-Thesis/img/" + "Boxplot_" + title)
    #plt.show()
    plt.close()
    
        
    # all pearson correlations, nach der Kreuzvalidierung und alle zu vergleichenden Werte auf einmal
    return



def make_histogram(helix_percentage, sheet_percentage, coil_percentage, title):
    fig = plt.figure(figsize=(10,6))
    fig.canvas.manager.set_window_title(title)
    all_percentage = np.array(helix_percentage) + np.array(sheet_percentage) + np.array(coil_percentage)
    #print("Average: ",sum(all_percentage)/len(all_percentage))
    #print("Median: ", np.median(all_percentage))

    #for i in range(len(helix_percentage)):
    #    print(helix_percentage[i],",", sheet_percentage[i],",", coil_percentage[i])
    
    y, x, _ = plt.hist(all_percentage, 100, density=True, facecolor="midnightblue")
    #----sns.kdeplot(all_percentage, shade=True, color="grey", bw_method="silverman")
    
    #plt.axvline(all_percentage.mean(), color='white', linestyle='dashed', linewidth=1)
    #plt.annotate(f"Mean: {round(all_percentage.mean(),2)}", (all_percentage.mean()+0.1, y.max()-0.1))

    #plt.axvline(np.median(all_percentage), color='red', linestyle='dashed', linewidth=1)
    #plt.annotate(f"Median: {round(np.median(all_percentage),2)}", (np.median(all_percentage)+0.1, y.max()-0.2), color="red")

    plt.xlim(-1,3)
    #######plt.savefig("../../Drafts/Project-Thesis/img/" + "Sum_percentage_" + title)
    plt.show()
    #plt.close()
    
    

    
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

    print(title)
    print(sum_single_measurements)
    print(sum_single_predictions)

    
    plt.bar(r1, sum_single_measurements, width = bar_width, color = "midnightblue", edgecolor = 'black', label="measurements")
 
    plt.bar(r2, sum_single_predictions, width = bar_width, color = "dodgerblue", edgecolor = 'black', label="predictions")

    plt.xticks([r + (bar_width/2) for r in range(3)], ["Helix", "Sheet", "Coil"])
    plt.ylabel("Secondary Structure in %")
    plt.legend(loc="upper left")
    #plt.savefig("../../Drafts/Project-Thesis/img/" + "Measure_vs_Predict_" + title)
    plt.show()
    #plt.close()
    


def make_histoplots(measurements, predictions, title):

    names = ["helix", "sheet", "coil"]
    
    
    for j, (sec_struc_measure, sec_struc_preds) in enumerate(zip(measurements, predictions)):
        pred_error_ = []
        sec_struc_measure = list(sec_struc_measure) * 1000
        for i, (measurement, pred) in enumerate(zip(sec_struc_measure, sec_struc_preds)):
            pred_error = measurement - pred
            pred_error_.append(pred_error)

        fig = plt.figure(figsize=(10,6))
        fig.canvas.manager.set_window_title(title + names[j])

        plt.xlim(-0.75,0.75)
        plt.hist(pred_error_, 100, density=True)
        #---sns.kdeplot(pred_error_, color="grey", shade=True, bw_method="silverman")
        ####plt.savefig("../../Drafts/Project-Thesis/img/" + "Error_Predict_" + str(j) + title)
        plt.show()
        plt.close()
        


def make_bar_std_plots(pearson_values, title):
    fig = plt.figure(figsize=(10,6))
    fig.canvas.manager.set_window_title(title)

    mean_ = []
    std_ = []
    for values in pearson_values:
        mean_.append(np.mean(values))
        std_.append(np.std(values))

    for i in range(4):
        plt.bar(f"model {i}", mean_[i])
        plt.errorbar(f"model {i}", mean_[i], yerr=std_[i], fmt="o", color="grey", capsize=100) 

    plt.savefig("../../Drafts/Project-Thesis/img/" + "Bar_predict_metric_" + title)
    #plt.show()
    plt.close()
    

    
from run_LR import load_dataset
from run_LR import build_simple_models, build_rebuild_models

from LR_data_preprocessing import seperate_secondary_structures

if __name__=="__main__":

    # test data
    X_test, y_test = load_dataset(488)
    
    ## big models
    X_train, y_train = load_dataset(2353)
    helix, sheet, coil = seperate_secondary_structures(y_test)
    # simple big model
    res3 = build_simple_models(X_train, X_test, y_train, y_test, modelsize="big")
    # rebuild big model
    res4 = build_rebuild_models(X_train, X_test, y_train, y_test, modelsize="big")

    #print()
    
    ## small models
    X_train, y_train = load_dataset(401)
    # simple small model
    res1 = build_simple_models(X_train, X_test, y_train, y_test, modelsize="small")
    # rebuild simple model
    res2 = build_rebuild_models(X_train, X_test, y_train, y_test, modelsize="small")

    
    
    
    plt.style.use("default")
    #----------------make barplots std pearson---------------------#
    #make_bar_std_plots([res1[0][1], res3[0][1], res2[0][1], res4[0][1]], "Comparison of Pearson Values Distribution for Linear Models Helix Prediction")
    #make_bar_std_plots([res1[1][1], res3[1][1], res2[1][1], res4[1][1]], "Comparison of Pearson Values Distribution for Linear Models Sheet Prediction")
    #make_bar_std_plots([res1[2][1], res3[2][1], res2[2][1], res4[2][1]], "Comparison of Pearson Values Distribution for Linear Models Coil Prediction")
    #----------------make barplots std rmse---------------------#
    #make_bar_std_plots([res1[0][2], res3[0][2], res2[0][2], res4[0][2]], "Comparison of RMSE Values Distribution for Linear Models Helix Prediction")
    #make_bar_std_plots([res1[1][2], res3[1][2], res2[1][2], res4[1][2]], "Comparison of RMSE Values Distribution for Linear Models Sheet Prediction")
    #make_bar_std_plots([res1[2][2], res3[2][2], res2[2][2], res4[2][2]], "Comparison of RMSE Values Distribution for Linear Models Coil Prediction")
   
    #--------------histoplots error distribution----------------#
    make_histoplots([np.array(helix).ravel(), np.array(sheet).ravel(), np.array(coil).ravel()],
                  [res4[0][0], res4[1][0], res4[2][0]],
                  "rebuild big")

    make_histoplots([np.array(helix).ravel(), np.array(sheet).ravel(), np.array(coil).ravel()],
                  [res3[0][0], res3[1][0], res3[2][0]],
                  "simple big")

    make_histoplots([np.array(helix).ravel(), np.array(sheet).ravel(), np.array(coil).ravel()],
                  [res2[0][0], res2[1][0], res2[2][0]],
                  "rebuild small")
    make_histoplots([np.array(helix).ravel(), np.array(sheet).ravel(), np.array(coil).ravel()],
                  [res1[0][0], res1[1][0], res1[2][0]],
                  "simple small")

    #--------------scatter plots------------------#

    #make_one_scatter(np.array(helix).ravel(), [res4[0][4],res3[0][4],res2[0][4],res1[0][4]], "all_in_one_helix")
    #make_one_scatter(np.array(sheet).ravel(), [res4[1][4],res3[1][4],res2[1][4],res1[1][4]], "all_in_one_sheet")
    #make_one_scatter(np.array(coil).ravel(), [res4[2][4],res3[2][4],res2[2][4],res1[2][4]], "all_in_one_coil")
    #print(len(helix))
    #make_scatter(np.array(helix).ravel(), res4[0][4], "rebuild big helix")
    #make_scatter(np.array(sheet).ravel(), res4[1][4], "rebuild big sheet")
    #make_scatter(np.array(coil).ravel(), res4[2][4], "rebuild big coil")

    #make_scatter(np.array(helix).ravel(), res3[0][4], "simple big helix")
    #make_scatter(np.array(sheet).ravel(), res3[1][4], "simple big sheet")
    #make_scatter(np.array(coil).ravel(), res3[2][4], "simple big coil")

    #make_scatter(np.array(helix).ravel(), res2[0][4], "rebuild small helix")
    #make_scatter(np.array(sheet).ravel(), res2[1][4], "rebuild small sheet")
    #make_scatter(np.array(coil).ravel(), res2[2][4], "rebuild small coil")

    #make_scatter(np.array(helix).ravel(), res1[0][4], "simple small helix")
    #make_scatter(np.array(sheet).ravel(), res1[1][4], "simple small sheet")
    #make_scatter(np.array(coil).ravel(), res1[2][4], "simple small coil")

    
    #-------------------barplots-----------------------#
    
    #make_barplots([np.array(helix).ravel(), np.array(sheet).ravel(), np.array(coil).ravel()],
    #              [res4[0][0], res4[1][0], res4[2][0]],
    #              "rebuild big")

    #make_barplots([np.array(helix).ravel(), np.array(sheet).ravel(), np.array(coil).ravel()],
    #              [res3[0][0], res3[1][0], res3[2][0]],
    #              "simple big")

    #make_barplots([np.array(helix).ravel(), np.array(sheet).ravel(), np.array(coil).ravel()],
    #              [res2[0][0], res2[1][0], res2[2][0]],
    #              "rebuild small")

    #make_barplots([np.array(helix).ravel(), np.array(sheet).ravel(), np.array(coil).ravel()],
    #              [res1[0][0], res1[1][0], res1[2][0]],
    #              "simple small")

    #-------------------histogram plots-----------------------#
    #for i in y_test:
    #    print(i[0], ",", i[1], ",", i[2])
    #print()
    make_histogram(res4[0][0], res4[1][0], res4[2][0], "rebuild big")
    #print()
    make_histogram(res3[0][0], res3[1][0], res3[2][0], "simple big")
    #print()
    make_histogram(res2[0][0], res2[1][0], res2[2][0], "rebuild small")
    #print()
    make_histogram(res1[0][0], res1[1][0], res1[2][0], "simple small")
    
    

     #-------------------pearson boxplots-----------------------#
    #make_boxplot([res1[0][1], res3[0][1], res2[0][1], res4[0][1]], "Comparison of Pearson Values Distribution for Linear Models Helix Prediction")
    #make_boxplot([res1[1][1], res3[1][1], res2[1][1], res4[1][1]], "Comparison of Pearson Values Distribution for Linear Models Sheet Prediction")
    #make_boxplot([res1[2][1], res3[2][1], res2[2][1], res4[2][1]], "Comparison of Pearson Values Distribution for Linear Models Coil Prediction")

     #-------------------rmse boxplots-----------------------#
    #make_boxplot([res1[0][2], res3[0][2], res2[0][2], res4[0][2]], "Comparison of RMSE Values Distribution for Linear Models Helix Prediction")
    #make_boxplot([res1[1][2], res3[1][2], res2[1][2], res4[1][2]], "Comparison of RMSE Values Distribution for Linear Models Sheet Prediction")
    #make_boxplot([res1[2][2], res3[2][2], res2[2][2], res4[2][2]], "Comparison of RMSE Values Distribution for Linear Models Coil Prediction")

