import math
import pandas as pd
import numpy as np
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
from statannot import add_stat_annotation

# FN = File Name OFN = Output File Name
def Density_Metagene(sourceFN, mRNALibFN, metaGeneOFN, densityOFN):
    normF = pd.read_csv(sourceFN)
    # print(normF)

    # d = 0
    # print(normF.loc[normF["ref_id"] == "R05D3.3.1"]["even_read_count"].sum())

    targetGene = pd.read_csv(mRNALibFN)
    targetGene.drop("Type", axis=1, inplace = True) 
    targetGene.drop("sequence", axis=1, inplace = True)
    # print(targetGene)


    # Densityplot data
    totalCount = []
    density = []
    len5 = []
    lenCDS = []
    len3 = []
    count5 = []
    countCDS = []
    count3 = []
    binRC = [[]] * len(targetGene["Gene name"])
    binSize = []
    binCol = []
    timesCount = 1
    while timesCount < 101:
        binCol.append('B%s'%timesCount)
        timesCount += 1
    metaGenePlot = pd.DataFrame()
    # metaGenePlot.columns = binCol
    binCol.append('Bin_Size')
    metaGenePlot['Gene name'] = targetGene['Gene name']
    metaGenePlot.set_index('Gene name', inplace = True)
    metaGenePlot = metaGenePlot.reindex(columns = binCol)
    # for i in binCol:
    #     metaGenePlot[i] = 0.0
    # metaGenePlot['Bin_Size'] = 0.0

    # Progress Visualization
    progess = tqdm(total = len(targetGene['Gene name']))
    wbreak = 0
    # Gene by Gene
    for GeneName in targetGene["Gene name"]:
        # total Read Count
        nowTotalRC = normF.loc[normF["ref_id"] == GeneName]["even_read_count"].sum()
        # totalCount.append(float(np.round(nowTotalRC, 9)))
        totalCount.append(nowTotalRC)
        # Info about Gene
        nowTargetGeneInfo = targetGene.loc[targetGene["Gene name"] == GeneName]
        # Gene total Length
        nowGeneSeqLength = int(nowTargetGeneInfo["sequence_length"])
        # density
        # density.append(float(np.round(nowTotalRC / nowGeneSeqLength, 9)))
        density.append(nowTotalRC / nowGeneSeqLength)
        # Bin Length
        nowBinSize = nowGeneSeqLength / 100
        # Gene CDS Start & End
        nowGeneCDSStrat = int(nowTargetGeneInfo["CDS start"])
        nowGeneCDSEnd = int(nowTargetGeneInfo["CDS end"])
        # Gene Len5 & LenCDS & Len3
        len5.append(nowGeneCDSStrat - 1)
        lenCDS.append(nowGeneCDSEnd - nowGeneCDSStrat + 1)
        len3.append(nowGeneSeqLength - nowGeneCDSEnd)
    
        # input in gene
        nowGene = normF.loc[normF["ref_id"] == GeneName]
        # Zerolize
        nowCount5 = 0
        nowCountCDS = 0
        nowCount3 = 0
        thisBinRC = [0.0] * 100

        nowGene.reset_index(inplace = True)

        for inputID in nowGene["index"]:
            # print(inputID)
            nowInputGene = nowGene.loc[nowGene["index"] == inputID]
            # print(nowInputGene)
            thisSeqStart = int(nowInputGene["init_pos"])
            thisSeqEnd = int(nowInputGene["end_pos"])
            thisSeqLength = thisSeqEnd - thisSeqStart + 1
            # print(thisSeqLength)
            thisGeneEvenRC = float(nowInputGene["even_read_count"])

            # check in length
            thisLen5Part = 0
            if thisSeqStart < nowGeneCDSStrat:
                if thisSeqEnd <= nowGeneCDSStrat:
                    thisLen5Part = thisSeqLength
                else:
                    thisLen5Part = nowGeneCDSStrat - thisSeqStart
            else:
                thisLen5Part = 0

            thisLen3Part = 0
            if thisSeqEnd > nowGeneCDSEnd:
                if thisSeqStart >= nowGeneCDSEnd:
                    thisLen3Part = thisSeqLength
                else:
                    thisLen3Part = thisSeqEnd - nowGeneCDSEnd
            else:
                thisLen3Part = 0

            # Length & Count
            thisLenCDS = thisSeqLength - thisLen5Part - thisLen3Part
            nowCount5 += thisGeneEvenRC * thisLen5Part / thisSeqLength
            nowCountCDS += thisGeneEvenRC * thisLenCDS / thisSeqLength
            nowCount3 += thisGeneEvenRC * thisLen3Part / thisSeqLength

            # bin check
            timesCount = 0
            center = int(((thisSeqEnd - thisSeqStart) / 2) + thisSeqStart)
            binNumber = math.ceil(center / nowBinSize)
            thisBinRC[binNumber - 1] += thisGeneEvenRC

        # count5.append(float(np.round(nowCount5, 9)))
        # countCDS.append(float(np.round(nowCountCDS, 9)))
        # count3.append(float(np.round(nowCount3, 9)))
        count5.append(nowCount5)
        countCDS.append(nowCountCDS)
        count3.append(nowCount3)
        timesCount = 0
        # binSize.append(nowBinSize)
        thisBinRC.append(nowBinSize)
        metaGenePlot.loc[GeneName] = thisBinRC
        progess.update(1)

        # break point check
        # if wbreak > 100:
        #     break
        # else:
        #     wbreak += 1

    metaGenePlot.to_csv(metaGeneOFN)

    len5PD = pd.DataFrame(len5, columns = ['Len 5'])
    targetGene["Len 5"] = len5PD

    lenCDSPD = pd.DataFrame(lenCDS, columns = ['Len CDS'])
    targetGene["Len CDS"] = lenCDSPD

    len3PD = pd.DataFrame(len3, columns = ['Len 3'])
    targetGene["Len 3"] = len3PD

    count5PD = pd.DataFrame(count5, columns = ['Count 5'])
    targetGene["Count 5"] = count5PD

    countCDSPD = pd.DataFrame(countCDS, columns = ['Count CDS'])
    targetGene["Count CDS"] = countCDSPD

    count3PD = pd.DataFrame(count3, columns = ['Count 3'])
    targetGene["Count 3"] = count3PD

    saveTargetGene = targetGene
    normColOrder = ['sequence_length', 'Gene name', 'CDS start', 'CDS end', 'Len 5', 'Len CDS', 'Len 3', 'Count 5', 'Count CDS', 'Count 3', 'Total Counts', 'Density']

    totalCountPD = pd.DataFrame(totalCount, columns = ['Total Counts'])
    targetGene["Total Counts"] = totalCountPD

    densityPD = pd.DataFrame(density, columns = ['Density'])
    targetGene["Density"] = densityPD

    targetGene.set_index("Gene ID", inplace = True)
    targetGene = targetGene.reindex(columns = normColOrder)
    targetGene.to_csv(densityOFN)

    twoFile = []
    twoFile.append(targetGene)
    twoFile.append(metaGenePlot)
    return(twoFile)

def preDensity(targetGroupName, densityPlot, mergeIndexName, groupName):
    targetGroupName = pd.merge(densityPlot, targetGroupName, on=[mergeIndexName])
    targetGroupName.set_index(['Gene ID'], inplace = True)
    targetGroupName['DataSet'] = groupName
    targetGroupName['log10(read counts/mRNA length * M)'] = np.log10((targetGroupName['Total Counts'] / targetGroupName['sequence_length'] * 10**6))
    return(targetGroupName)

def pre4Density(targetGroupName, densityPlot, mergeIndexName, groupName, targetGeneName):
    targetGroupName = pd.merge(densityPlot, targetGroupName, on=[mergeIndexName])
    targetGroupName.set_index(['Gene ID'], inplace = True)
    targetGroupName['DataSet'] = groupName + '\n' + targetGeneName
    targetGroupName['log10(read counts/mRNA length * M)'] = np.log10((targetGroupName['Total Counts'] / targetGroupName['sequence_length'] * 10**6))
    return(targetGroupName)

def plot2Density(targetGeneTranscript, targetGeneName, IFN1_SD, IFN1_category, IFN2_SD, IFN2_category):
    tGT = pd.read_csv(targetGeneTranscript)
    # DP = Density Plot
    first_DP = pd.read_csv(IFN1_SD)
    second_DP = pd.read_csv(IFN2_SD)
    tGTD_first = preDensity(tGT, first_DP.copy(), 'Gene name', IFN1_category)
    tGTD_second = preDensity(tGT, second_DP.copy(), 'Gene name', IFN2_category)
    # print(tGTD_first)
    # print(tGTD_second)
    bothDP = tGTD_first.copy()
    bothDP = bothDP.append(tGTD_second)
    # print(bothDP)

    x = "DataSet"
    y = "log10(read counts/mRNA length * M)"
    # sns.set(rc={'figure.figsize':(5, 5)})
    order = [IFN1_category, IFN2_category]
    ax = sns.boxplot(data = bothDP, x = x, y = y, order = order, showfliers = False)

    test_results = add_stat_annotation(ax, data = bothDP, x = x, y = y, order = order,
    box_pairs=[(IFN1_category, IFN2_category)],
    test='Mann-Whitney', text_format = 'full',loc = 'inside', verbose = 2)

    plt.title('%s target N = %s'%(targetGeneName, tGT.shape[0]))
    plt.savefig('%s_densityPlot.png'%targetGeneName)
    plt.close()

def plot4Density(targetGeneTranscript_1, targetGeneTranscript_2, targetGeneName_1, targetGeneName_2, IFN1_SD_1, IFN1_SD_2, IFN1_category, IFN2_SD_1, IFN2_SD_2, IFN2_category):
    tGT_1 = pd.read_csv(targetGeneTranscript_1)
    # DP = Density Plot
    first_DP_1 = pd.read_csv(IFN1_SD_1)
    second_DP_1 = pd.read_csv(IFN2_SD_1)
    tGTD_first_1 = pre4Density(tGT_1, first_DP_1.copy(), 'Gene name', IFN1_category, targetGeneName_1)
    tGTD_second_1 = pre4Density(tGT_1, second_DP_1.copy(), 'Gene name', IFN2_category, targetGeneName_1)

    tGT_2 = pd.read_csv(targetGeneTranscript_2)
    # DP = Density Plot
    first_DP_2 = pd.read_csv(IFN1_SD_2)
    second_DP_2 = pd.read_csv(IFN2_SD_2)
    tGTD_first_2 = pre4Density(tGT_2, first_DP_2.copy(), 'Gene name', IFN1_category, targetGeneName_2)
    tGTD_second_2 = pre4Density(tGT_2, second_DP_2.copy(), 'Gene name', IFN2_category, targetGeneName_2)

    allDP = tGTD_first_1.copy()
    allDP = allDP.append(tGTD_second_1)
    allDP = allDP.append(tGTD_first_2)
    allDP = allDP.append(tGTD_second_2)
    print(allDP)
    # print(bothDP)

    x = "DataSet"
    y = "log10(read counts/mRNA length * M)"
    sns.set(rc={'figure.figsize':(8, 10)})
    N_1 =  IFN1_category + '\n' + targetGeneName_1
    N_2 =  IFN2_category + '\n' + targetGeneName_1
    N_3 =  IFN1_category + '\n' + targetGeneName_2
    N_4 =  IFN2_category + '\n' + targetGeneName_2
    order = [N_1, N_2, N_3, N_4]
    ax = sns.boxplot(data = allDP, x = x, y = y, order = order, showfliers = False)

    test_results = add_stat_annotation(ax, data = allDP, x = x, y = y, order = order,
    box_pairs=[(N_1, N_2), (N_1, N_3), (N_3, N_4), (N_4, N_2)],
    test='Mann-Whitney', text_format = 'full',loc = 'inside', verbose = 2)

    plt.title('%s target v.s. %s target\n%s v.s. %s\n%s target N = %s\n%s target N = %s'%(targetGeneName_1, targetGeneName_2, IFN1_category, IFN2_category, targetGeneName_1, tGT_1.shape[0], targetGeneName_2, tGT_2.shape[0]))
    plt.savefig('%s_target_VS_%s_target.png'%(targetGeneName_1, targetGeneName_2))
    # plt.show()
    plt.close()

def preMetaGene(targetGroupName, metagenePlot, tCIndex, mergeIndexName, groupName):
    targetGroupName = pd.merge(metagenePlot, targetGroupName, on=[mergeIndexName])
    targetGroupName = pd.merge(tCIndex, targetGroupName, on=[mergeIndexName])
    targetGroupName.set_index([mergeIndexName], inplace = True)
    return(targetGroupName)

def preN2C(densityPlot):
    extration = pd.DataFrame()
    extration['Gene name'] = densityPlot['Gene name']
    extration['Total Counts'] = densityPlot['Total Counts']
    return(extration)

def drawSTE(STE, tarColor, xLab, yLab, subPlot, xScale, axes):
    STE = pd.DataFrame(STE, columns=['-STE', '+STE'])
    # print(perBinSTE_1)

    STE_P = STE.plot(xlabel = xLab, ylabel = yLab, legend=False, linestyle = '--', c = tarColor, linewidth = 0.5,
    ax=axes[subPlot], xlim = (0, xScale - 1))
    STE_P.fill_between(x = range(xScale),y1 = STE['-STE'].values.tolist(), y2 = STE['+STE'].values.tolist(), color = tarColor, alpha = 0.3)
    return(STE_P)

def plotMetaGene(targetGeneTranscript, targetGeneName, IFN1_Meta, IFN1_Den, IFN2_Meta, IFN2_Den, IFN1_category, IFN2_category):
    
    metagenePlot1 = pd.read_csv(IFN1_Meta)
    metagenePlot2 = pd.read_csv(IFN2_Meta)
    densityPlot1 = pd.read_csv(IFN1_Den)
    densityPlot2 = pd.read_csv(IFN2_Den)
    tCIndex1 = preN2C(densityPlot1)
    tCIndex1.set_index(['Gene name'], inplace = True)
    tCIndex2 = preN2C(densityPlot2)
    tCIndex2.set_index(['Gene name'], inplace = True)


    tMG_Trans = pd.read_csv(targetGeneTranscript)
    # print(wago)

    csr = pd.read_csv(targetGeneTranscript)
    # print(csr)

    # tMG = target MetaGenePlot
    tMG = preMetaGene(tMG_Trans, metagenePlot1.copy(), tCIndex1, 'Gene name', 'WAGO-1 target')
    csrMetaGene = preMetaGene(csr, metagenePlot2.copy(), tCIndex2, 'Gene name', 'CSR-1 target')
    # print(tMG)
    # print(csrMetaGene)

    tMGN = tMG.copy()
    csrMetaGeneN = csrMetaGene.copy()
    binOrder = 1
    while binOrder < 101:
        tMGN['B%s'%binOrder] = tMGN['B%s'%binOrder]/tMGN['Total Counts']
        csrMetaGeneN['B%s'%binOrder] = csrMetaGeneN['B%s'%binOrder]/csrMetaGeneN['Total Counts']
        binOrder += 1
    # print(tMG)
    # print(tMGN)
    # print(csrMetaGene)
    # print(csrMetaGeneN)

    perBinRatio = [[]] * 100
    perBinMV = [[]] * 100
    perBinWholeMV = [[]] * 100
    perBinSTE_1 = [[]] * 100
    perBinSTE_2 = [[]] * 100
    perBinSTEN_1 = [[]] * 100
    perBinSTEN_2 = [[]] * 100
    perBinMVN = [[]] * 100
    perBinSTEN = [[]] * 100
    binOrder = 1
    csrMetaGeneDF = len(csrMetaGene)
    tMGDF = len(tMG)
    progess = tqdm(total = 100, desc = targetGeneName + ' MetaGene')
    for binOrder in range(1, 101):
        # Calculate Ratio
        csrMetaGeneCount = ((csrMetaGene['B%s'%binOrder] != 0.0).sum()) / csrMetaGeneDF
        tMGCount = ((tMG['B%s'%binOrder] != 0.0).sum()) / tMGDF
        pair = []
        pair.append(csrMetaGeneCount)
        pair.append(tMGCount)
        perBinRatio[binOrder - 1] = pair
        # Calculate MeanValue & STE
        csrMV = csrMetaGene['B%s'%binOrder].sum()/len(csrMetaGene)
        csrSTE = np.std(csrMetaGene['B%s'%binOrder], ddof=1) / np.sqrt(np.size(csrMetaGene['B%s'%binOrder]))
        wholePair = []
        value = csrMV - csrSTE
        wholePair.append(value)
        # wholePair.append(csrMV)
        wholePair.append(csrMV + csrSTE)
        perBinSTE_1[binOrder - 1] = wholePair
        wagoMV = tMG['B%s'%binOrder].sum()/len(tMG)
        wagoSTE = np.std(tMG['B%s'%binOrder], ddof=1) / np.sqrt(np.size(tMG['B%s'%binOrder]))
        wholePair = []
        wholePair.append((wagoMV - wagoSTE))
        # wholePair.append(wagoMV)
        wholePair.append((wagoMV + wagoSTE))
        perBinSTE_2[binOrder - 1] = wholePair
        MVPair = []
        MVPair.append(csrMV)
        MVPair.append(wagoMV)
        perBinMV[binOrder - 1] = MVPair
        # Normalize
        csrMVN = csrMetaGeneN['B%s'%binOrder].sum()/len(csrMetaGeneN)
        csrSTEN = np.std(csrMetaGeneN['B%s'%binOrder], ddof=1) / np.sqrt(np.size(csrMetaGeneN['B%s'%binOrder]))
        # print(csrSTEN)
        csrNPair = []
        csrNPair.append(csrMVN - csrSTEN)
        # csrNPair.append(csrMVN)
        csrNPair.append(csrMVN + csrSTEN)

        wagoMVN = tMGN['B%s'%binOrder].sum()/len(tMGN)
        wagoSTEN = np.std(tMGN['B%s'%binOrder], ddof=1) / np.sqrt(np.size(tMGN['B%s'%binOrder]))
        wagoNPair = []
        wagoNPair.append(wagoMVN - wagoSTEN)
        # wagoNPair.append(wagoMVN)
        wagoNPair.append(wagoMVN + wagoSTEN)
        MVNPair = []
        MVNPair.append(csrMVN)
        MVNPair.append(wagoMVN)
        perBinMVN[binOrder - 1] = MVNPair
        perBinSTEN_1[binOrder - 1] = csrNPair
        perBinSTEN_2[binOrder - 1] = wagoNPair
        progess.update(1)

    fig, axes = plt.subplots(nrows=3, ncols=1)

    # Ratio
    perBinRatio = pd.DataFrame(perBinRatio, columns=['%s mRNA site ratio (site/position)'%IFN1_category, '%s mRNA site ratio (site/position)'%IFN2_category])
    # print(perBinRatio)

    # plt.legend(bbox_to_anchor =(1.65, 1.15))
    ratioChart = perBinRatio.plot(title='%s target N = %s'%(targetGeneName, tMG_Trans.shape[0]),
    xlabel='region(%)', ylabel='ratio mRNA', figsize=(3, 10), ax=axes[0], xlim = (0, 99))
    plt.xlim(0, 99)
    ratioChart.legend(bbox_to_anchor =(1, 1))
    # plt.savefig('RatioChart.png', bbox_inches='tight')



    # Read counts

    perBinMV = pd.DataFrame(perBinMV, columns=['CSR-1 target AVG (+/-STE)', 'WAGO-1 target AVG (+/-STE)'])
    # print(perBinMV)

    ste1_1 = perBinMV.plot(xlabel='region(%)', ylabel='read counts', linewidth = 1, ax=axes[1], xlim = (0, 99))

    xLab = 'region(%)'
    yLab = 'read counts'
    subPlot = 1
    xScale = 100
    tarColor = 'royalblue'
    rCC_STE_1 = drawSTE(perBinSTE_1, tarColor, xLab, yLab, subPlot, xScale, axes)
    tarColor = 'darkorange'
    rCC_STE_2 = drawSTE(perBinSTE_2, tarColor, xLab, yLab, subPlot, xScale, axes)

    # axes[1].get_legend().remove()
    axes[1].legend(bbox_to_anchor =(1, 1))
    # plt.savefig('perBinMV.png', bbox_inches='tight')


    # Read counts Distribution
    perBinMVN = pd.DataFrame(perBinMVN, columns=['CSR-1 target AVG (+/-STE)', 'WAGO-1 target AVG (+/-STE)'])
    # print(perBinMVN)

    rCountNChart = perBinMVN.plot(xlabel='region(%)', ylabel='read counts Distribution', linewidth = 1, ax=axes[2], xlim = (0, 99))

    xLab = 'region(%)'
    yLab = 'read counts Distribution'
    subPlot = 2
    xScale = 100
    tarColor = 'royalblue'
    rCCN_STE_1 = drawSTE(perBinSTEN_1, tarColor, xLab, yLab, subPlot, xScale, axes)
    tarColor = 'darkorange'
    rCCN_STE_2 = drawSTE(perBinSTEN_2, tarColor, xLab, yLab, subPlot, xScale, axes)

    plt.legend(bbox_to_anchor =(1, 1))
    plt.savefig('%s_MetaGene.png'%targetGeneName, bbox_inches='tight')
    plt.close()