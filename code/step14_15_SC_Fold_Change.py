import numpy as np
import math
import pandas as pd
import seaborn as sns
from tqdm import tqdm
import matplotlib.pyplot as plt
from statannot import add_stat_annotation

def SC_Fold_Change(firstDIFN, firstScOFN, secondDIFN, secondScOFN, foldChangeOFN):
    # first
    firstDF = pd.read_csv(firstDIFN)
    firstDF['Density'] = np.round(firstDF['Density'], 9)
    firstDF['mockDensity'] = firstDF['Density']
    firstDF.loc[firstDF['mockDensity'] == 0, 'mockDensity'] = 0.00001
    firstDF['sc'] = np.round(np.log2(firstDF['mockDensity']), 9)
    firstDF.set_index('Gene ID', inplace = True)
    # print(firstDF)


    # second
    secondDF = pd.read_csv(secondDIFN)
    secondDF['Density'] = np.round(secondDF['Density'], 9)
    secondDF['mockDensity'] = secondDF['Density']
    secondDF.loc[secondDF['mockDensity'] == 0, 'mockDensity'] = 0.00001
    secondDF['sc'] = np.round(np.log2(secondDF['mockDensity']), 9) 
    secondDF.set_index('Gene ID', inplace = True)
    # print(secondDF)
    

    foldChangePlot = pd.DataFrame()
    foldChangePlot['Gene name'] = firstDF['Gene name']
    alpha = 0.0001748
    foldChangePlot['fold'] = secondDF['sc'] - firstDF['sc']
    # foldChangePlot['fold'] = np.log2((secondDF['mockDensity'] + alpha) / (firstDF['mockDensity'] + alpha))
    FDF = pd.read_csv(firstDIFN)
    FDF['Density'] = np.round(FDF['Density'], 9)
    SDF = pd.read_csv(secondDIFN)
    SDF['Density'] = np.round(SDF['Density'], 9)
    FDF.loc[FDF['Density'] == 0, 'Density'] = None
    SDF.loc[SDF['Density'] == 0, 'Density'] = None
    FDF.loc[FDF['Density'] != 0, 'fold_without0'] = np.log2(SDF['Density'] / FDF['Density'])
    FDF.set_index('Gene name', inplace = True)
    # print(firstDF)
    foldChangePlot.reset_index(inplace = True)
    foldChangePlot.set_index('Gene name', inplace = True)
    foldChangePlot.drop("Gene ID", axis=1, inplace = True) 
    foldChangePlot['fold_without0'] = FDF['fold_without0']
    foldChangePlot['fold_without0'].fillna(value='NULL', inplace=True)
    foldChangePlot.to_csv(foldChangeOFN)
    # print(foldChangePlot)

    firstDF.drop("mockDensity", axis=1, inplace = True) 
    firstDF.to_csv(firstScOFN)
    secondDF.drop("mockDensity", axis=1, inplace = True)
    secondDF.to_csv(secondScOFN)

    threePlot = []
    threePlot.append(firstDF)
    threePlot.append(secondDF)
    threePlot.append(foldChangePlot)
    return(threePlot)

def preN2sc(plotN, mergeIndexName, Plot_1, plotN_1, Plot_2, plotN_2):
    extration1 = pd.DataFrame()
    extration1[mergeIndexName] = Plot_1[mergeIndexName]
    extration1[plotN_1] = Plot_1['sc']
    extration2 = pd.DataFrame()
    extration2[mergeIndexName] = Plot_2[mergeIndexName]
    extration2[plotN_2] = Plot_2['sc']
    extration = pd.merge(extration1, plotN, on=[mergeIndexName])
    extration = pd.merge(extration, extration2, on=[mergeIndexName])
    extration.set_index(mergeIndexName, inplace = True)
    return(extration)

def preN2box(plotN, foldChange, mergeIndexName):
    extration = pd.merge(foldChange, plotN, on=[mergeIndexName])
    extration.set_index(mergeIndexName, inplace = True)
    return(extration)

def preComparison(target, yLabel, col_1, col_2):
    plot = pd.DataFrame()
    plot['%s'%yLabel] = target['%s'%col_1]
    plot['%s'%col_2] = target['%s'%col_2]
    return(plot)

def scatter_Plot(scattePlot, plotN_1, plotN_2, axes, plotWhere):
    scatterPlotPic = scattePlot.plot.scatter(x = plotN_1,y = plotN_2, ax=axes[plotWhere], s = 1, c = 'red', figsize=(15, 10))

    lineX = list(scatterPlotPic.get_xlim())
    lineY = list(scatterPlotPic.get_ylim())
    scatterPlotPic.plot(lineX, lineY, color = 'black')
    scatterPlotPic.plot([lineX[0] + 1, lineX[1] + 1], [lineY[0], lineY[1]], linestyle = '--', color = 'black')
    scatterPlotPic.plot([lineX[0], lineX[1]], [lineY[0] + 1, lineY[1] + 1], linestyle = '--', color = 'black')

def foldChange_Plot(targetGeneName, targetGeneTranscript, n_1, n_2, plotN_1, plotN_2, IFN1_cat, IFN2_cat, mergeIndexName, foldChangeN):

    scPlot_1 = pd.read_csv(n_1) #head
    scPlot_2 = pd.read_csv(n_2)
    plotN = pd.read_csv(targetGeneTranscript) # head

    scattePlot = preN2sc(plotN, mergeIndexName, scPlot_1.copy(), plotN_1, scPlot_2.copy(), plotN_2)

    # print(scattePlot)

    fig, axes = plt.subplots(nrows=1, ncols=2)
    fig.suptitle('%s v.s. %s %s target\n#FULL GENE=%s'%(IFN1_cat, IFN2_cat, targetGeneName, len(plotN)))

    plotWhere = 0
    scatter_Plot(scattePlot, plotN_1, plotN_2, axes, plotWhere)


    foldChange = pd.read_csv(foldChangeN)
    foldChange.set_index(mergeIndexName, inplace = True)
    FC_box = preN2box(plotN, foldChange, mergeIndexName)
    # print(foldChange)
    boxplot = FC_box.boxplot(column=['fold', 'fold_without0'], ax=axes[1], showfliers = False, figsize=(15, 10))
    plt.savefig('%s_FoldChange.png'%targetGeneName, bbox_inches='tight')
    plt.close()

def comparison_Plot(targetGeneName_1, targetGeneName_2, targetGeneTranscript_1, targetGeneTranscript_2, n_1, n_2, plotN_1, plotN_2, IFN1_cat, IFN2_cat, mergeIndexName, foldChangeN):
    
    scPlot_1_1 = pd.read_csv(n_1) 
    scPlot_2_1 = pd.read_csv(n_2)
    plotN_first = pd.read_csv(targetGeneTranscript_1)
    scattePlot = preN2sc(plotN_first, mergeIndexName, scPlot_1_1.copy(), plotN_1, scPlot_2_1.copy(), plotN_2)

    fig, axes = plt.subplots(nrows=1, ncols=2)
    scatterPlotPic_1 = scattePlot.plot.scatter(x = plotN_1,y = plotN_2, ax=axes[0], s = 1, c = 'red', title='%s target'%targetGeneName_1, figsize=(15, 10))
    lineX = list(scatterPlotPic_1.get_xlim())
    lineY = list(scatterPlotPic_1.get_ylim())
    scatterPlotPic_1.plot(lineX, lineY, color = 'black')
    scatterPlotPic_1.plot([lineX[0] + 1, lineX[1] + 1], [lineY[0], lineY[1]], linestyle = '--', color = 'black')
    scatterPlotPic_1.plot([lineX[0], lineX[1]], [lineY[0] + 1, lineY[1] + 1], linestyle = '--', color = 'black')

    scPlot_1_2 = pd.read_csv(n_1) #head
    scPlot_2_2 = pd.read_csv(n_2)
    plotN_second = pd.read_csv(targetGeneTranscript_2) # head
    scattePlot = preN2sc(plotN_second, mergeIndexName, scPlot_1_2.copy(), plotN_1, scPlot_2_2.copy(), plotN_2)
    scatterPlotPic_2 = scattePlot.plot.scatter(x = plotN_1,y = plotN_2, ax=axes[1], s = 1, c = 'red', title='%s target'%targetGeneName_2, figsize=(15, 10))
    lineX = list(scatterPlotPic_2.get_xlim())
    lineY = list(scatterPlotPic_2.get_ylim())
    scatterPlotPic_2.plot(lineX, lineY, color = 'black')
    scatterPlotPic_2.plot([lineX[0] + 1, lineX[1] + 1], [lineY[0], lineY[1]], linestyle = '--', color = 'black')
    scatterPlotPic_2.plot([lineX[0], lineX[1]], [lineY[0] + 1, lineY[1] + 1], linestyle = '--', color = 'black')

    fig.suptitle('%s v.s. %s, %s target v.s. %s target\n#%s target=%s\n#%s target=%s'%(IFN1_cat, IFN2_cat, targetGeneName_1, targetGeneName_2, targetGeneName_1, len(plotN_first), len(plotN_second), targetGeneName_2))
    plt.savefig('ScatterPlot_Comparison.png')
    plt.close()

    foldChange = pd.read_csv(foldChangeN)
    foldChange.set_index(mergeIndexName, inplace = True)
    plotN_first = pd.read_csv(targetGeneTranscript_1)
    plotN_second = pd.read_csv(targetGeneTranscript_2)
    FC_1 = preN2box(plotN_first, foldChange, mergeIndexName)
    FC_2 = preN2box(plotN_second, foldChange, mergeIndexName)
    FC_1['DataSet'] = '%s target'%targetGeneName_1
    FC_2['DataSet'] = '%s target'%targetGeneName_2

    col_1 = 'fold'
    col_2 = 'DataSet'
    yLabel = 'log2(%s + alpha / %s + alpha)'%(IFN2_cat, IFN2_cat)
    with0_Plot_1 = preComparison(FC_1, yLabel, col_1, col_2)
    with0_Plot_2 = preComparison(FC_2, yLabel, col_1, col_2)
    with0_Plot_all = with0_Plot_1.copy()
    with0_Plot_all = with0_Plot_all.append(with0_Plot_2)
    x = "DataSet"
    y = yLabel
    N_1 = '%s target'%targetGeneName_1
    N_2 = '%s target'%targetGeneName_2
    sns.set(rc={'figure.figsize':(10, 8)})
    order = [N_1, N_2]
    ax = sns.boxplot(data = with0_Plot_all, x = x, y = y, order = order, showfliers = False)

    test_results = add_stat_annotation(ax, data = with0_Plot_all, x = x, y = y, order = order,
    box_pairs=[(N_1, N_2)],
    test='Mann-Whitney', text_format = 'full',loc = 'inside', verbose = 2)
    plt.savefig('With0_Comparison.png')
    plt.title('%s v.s. %s\n %s target v.s. %s target\n#%s target=%s\n#%s target=%s\n#alpha=%s'%(IFN1_cat, IFN2_cat, targetGeneName_1, targetGeneName_2, targetGeneName_1, len(plotN_first), targetGeneName_2, len(plotN_second), '5.2e-05'))
    plt.close()

    foldChange = pd.read_csv(foldChangeN)
    foldChange.set_index(mergeIndexName, inplace = True)
    plotN_first = pd.read_csv(targetGeneTranscript_1)
    plotN_second = pd.read_csv(targetGeneTranscript_2)
    FC_1 = preN2box(plotN_first, foldChange, mergeIndexName)
    FC_2 = preN2box(plotN_second, foldChange, mergeIndexName)
    FC_1['DataSet'] = '%s target'%targetGeneName_1
    FC_2['DataSet'] = '%s target'%targetGeneName_2

    col_1 = 'fold_without0'
    col_2 = 'DataSet'
    yLabel = 'log2(%s + alpha / %s + alpha)'%(IFN2_cat, IFN2_cat)
    with0_Plot_1 = preComparison(FC_1, yLabel, col_1, col_2)
    with0_Plot_2 = preComparison(FC_2, yLabel, col_1, col_2)
    with0_Plot_all = with0_Plot_1.copy()
    with0_Plot_all = with0_Plot_all.append(with0_Plot_2)
    x = "DataSet"
    y = yLabel
    N_1 = '%s target'%targetGeneName_1
    N_2 = '%s target'%targetGeneName_2
    sns.set(rc={'figure.figsize':(10, 8)})
    order = [N_1, N_2]
    ax = sns.boxplot(data = with0_Plot_all, x = x, y = y, order = order, showfliers = False)

    test_results = add_stat_annotation(ax, data = with0_Plot_all, x = x, y = y, order = order,
    box_pairs=[(N_1, N_2)],
    test='Mann-Whitney', text_format = 'full',loc = 'inside', verbose = 2)
    plt.savefig('Without0_Comparison.png')
    plt.title('%s v.s. %s\n %s target v.s. %s target\n#%s target=%s\n#%s target=%s\n#alpha=%s'%(IFN1_cat, IFN2_cat, targetGeneName_1, targetGeneName_2, targetGeneName_1, len(plotN_first), targetGeneName_2, len(plotN_second), '5.2e-05'))
    plt.close()