from http.client import REQUESTED_RANGE_NOT_SATISFIABLE
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt

# RA = Require Amount
def Head_Tail(sourceFN, mRNALibFN, headRA, tailRA, headPlotOFN, tailPlotOFN):
    normF = pd.read_csv(sourceFN)
    # print(normF)

    # d = 0
    # print(normF.loc[normF["ref_id"] == "R05D3.3.1"]["even_read_count"].sum())

    targetGene = pd.read_csv(mRNALibFN)
    targetGene.drop("Type", axis=1, inplace = True) 
    targetGene.drop("sequence", axis=1, inplace = True)
    # print(targetGene)


    # Head or Tail data
    totalCount = []
    density = []
    # Require Amount Col
    rACol = []
    for timesCount in range(0, headRA):
        rACol.append('%s'%timesCount)
    # print(rACol)

    headPlot = pd.DataFrame()
    headPlot['Gene name'] = targetGene['Gene name']
    headPlot.set_index('Gene name', inplace = True)
    headPlot = headPlot.reindex(columns = rACol)
    # print(headPlot)

    rACol = []
    for timesCount in range(0, tailRA):
        rACol.append('%s'%timesCount)

    tailPlot = pd.DataFrame()
    tailPlot['Gene name'] = targetGene['Gene name']
    tailPlot.set_index('Gene name', inplace = True)
    tailPlot = tailPlot.reindex(columns = rACol)
    # print(TailPlot)

    # Progress Visualization
    progess = tqdm(total = len(targetGene['Gene name']))
    wbreak = 0
    # Gene by Gene
    for GeneName in targetGene["Gene name"]:
        # input in gene
        nowGene = normF.loc[normF["ref_id"] == GeneName]
        # Info about Gene
        nowTargetGeneInfo = targetGene.loc[targetGene["Gene name"] == GeneName]
        # Gene total Length
        nowGeneSeqLength = int(nowTargetGeneInfo["sequence_length"])
        # Zerolize
        thisHeadRC = [0.0] * headRA
        thisTailRC = [0.0] * tailRA

        nowGene.reset_index(inplace = True)

        # each mRNA
        for inputID in nowGene["index"]:
            # print(inputID)
            nowInputGene = nowGene.loc[nowGene["index"] == inputID]
            # print(nowInputGene)
            thisSeqStart = int(nowInputGene["init_pos"])
            thisSeqEnd = int(nowInputGene["end_pos"])
            thisSeqLength = thisSeqEnd - thisSeqStart + 1
            # print(thisSeqLength)
            thisGeneEvenRC = float(np.round(nowInputGene["even_read_count"], 3))

            # CT = Coodinate Transform
            thisSeqStartHCT = thisSeqStart - 1
            thisSeqEndHCT = thisSeqEnd - 1
            regionEnd = headRA - 1
            regionStart = 0

            startSeq = 0
            closeSeq = 0
            overlap = 0
            if thisSeqStartHCT <= regionEnd:
                overlap = 1
                startSeq = thisSeqStartHCT
                if thisSeqEndHCT < regionEnd:
                    closeSeq = thisSeqEndHCT
                else:
                    closeSeq = regionEnd
            else:
                overlap = 0
            
            if overlap:
                for timesCount in range(startSeq + 1, closeSeq + 1):
                    thisHeadRC[timesCount] += thisGeneEvenRC
            
            # tail
            thisSeqStartHCT = thisSeqStart - nowGeneSeqLength
            thisSeqEndHCT = thisSeqEnd - nowGeneSeqLength
            regionEnd = 1 - headRA
            regionStart = 0

            startSeq = 0
            closeSeq = 0
            overlap = 0
            if thisSeqEndHCT >= regionEnd:
                overlap = 1
                startSeq = abs(thisSeqEndHCT)
                if thisSeqStartHCT < regionEnd:
                    closeSeq = abs(regionEnd)
                else:
                    closeSeq = abs(thisSeqStartHCT)
            else:
                overlap = 0
            
            if overlap:
                for timesCount in range(startSeq + 1, closeSeq + 1):
                    thisTailRC[timesCount] += thisGeneEvenRC

        headPlot.loc[GeneName] = thisHeadRC
        tailPlot.loc[GeneName] = thisTailRC
        progess.update(1)

        # break point check
        # if wbreak > 50:
        #     break
        # else:
        #     wbreak += 1

    headPlot.reset_index(inplace=True)
    headPlot.drop("Gene name", axis=1, inplace = True) 
    headPlot.set_index('0', inplace = True)
    headPlot.to_csv(headPlotOFN)

    tailPlot.reset_index(inplace=True)
    tailPlot.drop("Gene name", axis=1, inplace = True) 
    tailPlot.set_index('0', inplace = True)
    tailPlot.to_csv(tailPlotOFN)

    HTPlot = []
    HTPlot.append(headPlot)
    HTPlot.append(tailPlot)
    return(HTPlot)

def preGene(targetGroupName, metagenePlot, tCIndex, mergeIndexName):
    tCIndex.reset_index(inplace=True)
    metagenePlot['Gene name'] = tCIndex['Gene name']
    metagenePlot['Total Counts'] = tCIndex['Total Counts']
    targetGroupName = pd.merge(metagenePlot, targetGroupName, on=[mergeIndexName])
    targetGroupName.set_index([mergeIndexName], inplace = True)
    return(targetGroupName)

def preN2C(densityPlot):
    extration = pd.DataFrame()
    extration['Gene name'] = densityPlot['Gene name']
    extration['Total Counts'] = densityPlot['Total Counts']
    return(extration)

def drawSTE(STE, tarColor, xLab, yLab, subPlotX, subPlotY, xScale, xRangeStart, xRangeEnd, axes):
    STE = pd.DataFrame(STE, columns=['-STE', '+STE'])
    # print(perBinSTE_1)

    STE_P = STE.plot(xlabel = xLab, ylabel = yLab, legend=False, linestyle = '--', c = tarColor, linewidth = 0.5,
    ax=axes[subPlotX, subPlotY], xlim = (xRangeStart, xRangeEnd))
    STE_P.fill_between(x = range(xScale),y1 = STE['-STE'].values.tolist(), y2 = STE['+STE'].values.tolist(), color = tarColor, alpha = 0.3)
    return(STE_P)

def draw4Plot(HD_1, HD_2, targetGeneName, IFN1_category, IFN2_category, axes, posN, xRangeStart, xRangeEnd, fileN):
    HD_N_1 = HD_1.copy()
    HD_N_2 = HD_2.copy()
    binOrder = 0
    # print(HD_N_1.sum(axis = 1))
    HD_N_1['accCounts'] = HD_N_1.sum(axis = 1)
    HD_N_2['accCounts'] = HD_N_2.sum(axis = 1)
    HD_N_factor_1 = HD_N_1['accCounts'].sum()
    HD_N_factor_2 = HD_N_2['accCounts'].sum()

    progess = tqdm(total = 100, desc = 'File Preprocess')
    for binOrder in range(0, 100):
        HD_N_1['%s'%binOrder] = HD_N_1['%s'%binOrder]/HD_N_factor_1
        HD_N_2['%s'%binOrder] = HD_N_2['%s'%binOrder]/HD_N_factor_2
        progess.update(1)

    perBinRatio = [[]] * 100
    perBinMV = [[]] * 100
    perBinSTE_1 = [[]] * 100
    perBinSTE_2 = [[]] * 100
    perBinSTEN_1 = [[]] * 100
    perBinSTEN_2 = [[]] * 100
    perBinMVN = [[]] * 100
    sitePosition = [[]] * 100
    progess = tqdm(total = 100, desc = fileN)
    for binOrder in range(0, 100):
        # Calculate Ratio
        # position
        position_1 = HD_1.loc[HD_1['%s'%binOrder] != 'N']['%s'%binOrder].count()
        # print(position_1)
        position_2 = HD_2.loc[HD_2['%s'%binOrder] != 'N']['%s'%binOrder].count()
        
        # site
        Site_1 = HD_1.loc[HD_1['%s'%binOrder] > 0]['%s'%binOrder].count()
        Site_2 = HD_2.loc[HD_2['%s'%binOrder] > 0]['%s'%binOrder].count()

        pair = []
        pair.append(position_1)
        pair.append(Site_1)
        pair.append(position_2)
        pair.append(Site_2)
        sitePosition[binOrder] = pair

        R_1 = Site_1 / position_1
        R_2 = Site_2 / position_2
        pair = []
        pair.append(R_1)
        pair.append(R_2)
        perBinRatio[binOrder] = pair
        # Calculate MeanValue & STE
        MV_1 = HD_1['%s'%binOrder].sum()/len(HD_1)
        STE_1 = np.std(HD_1['%s'%binOrder], ddof=1) / np.sqrt(np.size(HD_1['%s'%binOrder]))
        wholePair = []
        wholePair.append(MV_1 - STE_1)
        # wholePair.append(MV_1)
        wholePair.append(MV_1 + STE_1)
        perBinSTE_1[binOrder] = wholePair
        MV_2 = HD_2['%s'%binOrder].sum()/len(HD_2)
        STE_2 = np.std(HD_2['%s'%binOrder], ddof=1) / np.sqrt(np.size(HD_2['%s'%binOrder]))
        wholePair = []
        wholePair.append((MV_2 - STE_2))
        # wholePair.append(MV_2)
        wholePair.append((MV_2 + STE_2))
        perBinSTE_2[binOrder] = wholePair
        MVPair = []
        MVPair.append(MV_1)
        MVPair.append(MV_2)
        perBinMV[binOrder] = MVPair
        # Normalize
        HD_N_M_1 = HD_N_1['%s'%binOrder].sum()/len(HD_N_1)
        STE_1N = np.std(HD_N_1['%s'%binOrder], ddof=1) / np.sqrt(np.size(HD_N_1['%s'%binOrder]))
        Pair_N_1 = []
        Pair_N_1.append(HD_N_M_1 - STE_1N)
        Pair_N_1.append(HD_N_M_1 + STE_1N)

        HD_N_M_2 = HD_N_2['%s'%binOrder].sum()/len(HD_N_2)
        STE_2N = np.std(HD_N_2['%s'%binOrder], ddof=1) / np.sqrt(np.size(HD_N_2['%s'%binOrder]))
        Pair_N_2 = []
        Pair_N_2.append(HD_N_M_2 - STE_2N)
        Pair_N_2.append(HD_N_M_2 + STE_2N)
        MVNPair = []
        MVNPair.append(HD_N_M_1)
        MVNPair.append(HD_N_M_2)
        perBinMVN[binOrder] = MVNPair
        perBinSTEN_1[binOrder] = Pair_N_1
        perBinSTEN_2[binOrder] = Pair_N_2
        progess.update(1)

    # site position
    sitePosition = pd.DataFrame(sitePosition, columns=['%s mRNA have position'%IFN2_category, '%s mRNA have position'%IFN1_category, '%s mRNA have site'%IFN2_category, '%s mRNA have site'%IFN1_category])
    # print(perBinRatio)

    sipoChart = sitePosition.plot(title='%s target_HEAD100'%targetGeneName, legend = False,
    xlabel='region(%)', ylabel='# of mRNA', figsize=(18, 14), ax=axes[0, posN], xlim = (xRangeStart, xRangeEnd))
    plt.xlim(0, 99)

    # Ratio
    perBinRatio = pd.DataFrame(perBinRatio, columns=['%s mRNA site ratio (site/position)'%IFN1_category, '%s mRNA site ratio (site/position)'%IFN2_category])

    ratioChart = perBinRatio.plot(legend = False,
    xlabel='region(%)', ylabel='ratio mRNA', figsize=(18, 14), ax=axes[1, posN], xlim = (xRangeStart, xRangeEnd))
    plt.xlim(0, 99)



    # Read counts

    perBinMV = pd.DataFrame(perBinMV, columns=['%s target AVG (+/-STE)'%IFN1_category, '%s target AVG (+/-STE)'%IFN2_category])
    # print(perBinMV)

    ste1_1 = perBinMV.plot(xlabel='region(%)', ylabel='read counts', linewidth = 1, ax=axes[2, posN], xlim = (xRangeStart, xRangeEnd), legend = False)

    xLab = 'region(%)'
    yLab = 'read counts'
    subPlotX = 2
    subPlotY = posN
    xScale = 100
    tarColor = 'royalblue'
    rCC_STE_1 = drawSTE(perBinSTE_1, tarColor, xLab, yLab, subPlotX, subPlotY, xScale, xRangeStart, xRangeEnd, axes)
    tarColor = 'darkorange'
    rCC_STE_2 = drawSTE(perBinSTE_2, tarColor, xLab, yLab, subPlotX, subPlotY, xScale, xRangeStart, xRangeEnd, axes)

    # Read counts Distribution
    perBinMVN = pd.DataFrame(perBinMVN, columns=['%s target AVG (+/-STE)'%IFN1_category, '%s target AVG (+/-STE)'%IFN2_category])
    # print(perBinMVN)

    rCountNChart = perBinMVN.plot(xlabel='region(%)', ylabel='read counts Distribution', linewidth = 1, ax=axes[3, posN], xlim = (xRangeStart, xRangeEnd), legend = False)

    xLab = 'region(%)'
    yLab = 'read counts Distribution'
    subPlotX = 3
    subPlotY = posN
    xScale = 100
    tarColor = 'royalblue'
    rCCN_STE_1 = drawSTE(perBinSTEN_1, tarColor, xLab, yLab, subPlotX, subPlotY, xScale, xRangeStart, xRangeEnd, axes)
    tarColor = 'darkorange'
    rCCN_STE_2 = drawSTE(perBinSTEN_2, tarColor, xLab, yLab, subPlotX, subPlotY, xScale, xRangeStart, xRangeEnd, axes)

    allPlot = []
    allPlot.append(sipoChart)
    allPlot.append(ratioChart)
    allPlot.append(ste1_1)
    allPlot.append(rCountNChart)
    return(allPlot)


def plotHeadTail(targetGeneTranscript, targetGeneName, IFN1_Head, IFN1_Tail, IFN1_Den, IFN2_Head, IFN2_Tail, IFN2_Den, IFN1_category, IFN2_category):
    
    headPlot1 = pd.read_csv(IFN1_Head)
    tailPlot1 = pd.read_csv(IFN1_Tail)
    headPlot2 = pd.read_csv(IFN2_Head)
    tailPlot2 = pd.read_csv(IFN2_Tail)
    densityPlot1 = pd.read_csv(IFN1_Den)
    densityPlot2 = pd.read_csv(IFN2_Den)
    tCIndex1 = preN2C(densityPlot1)
    tCIndex1.set_index(['Gene name'], inplace = True)
    tCIndex2 = preN2C(densityPlot2)
    tCIndex2.set_index(['Gene name'], inplace = True)


    GeneTranscript = pd.read_csv(targetGeneTranscript)

    # tMG = target MetaGenePlot
    HD_1 = preGene(GeneTranscript, headPlot1.copy(), tCIndex1, 'Gene name')
    HD_2 = preGene(GeneTranscript, headPlot2.copy(), tCIndex2, 'Gene name')

    fig, axes = plt.subplots(nrows=4, ncols=2)

    posN = 0
    xRangeStart = 0
    xRangeEnd = 99
    fileN = targetGeneName + ' HEAD'
    draw4Plot(HD_1, HD_2, targetGeneName, IFN1_category, IFN2_category, axes, posN, xRangeStart, xRangeEnd, fileN)

# ####
    HD_1 = preGene(GeneTranscript, tailPlot1.copy(), tCIndex1, 'Gene name')
    HD_2 = preGene(GeneTranscript, tailPlot2.copy(), tCIndex2, 'Gene name')

    posN = 1
    xRangeStart = 99
    xRangeEnd = 0
    fileN = targetGeneName + ' TAIL'
    allPolt = draw4Plot(HD_1, HD_2, targetGeneName, IFN1_category, IFN2_category, axes, posN, xRangeStart, xRangeEnd, fileN)

    progess = tqdm(total = 4, desc = 'Adjust Plot and Save Plot')
    for index in range(0, 4):
        allPolt[index].legend(bbox_to_anchor =(1, 1))
        progess.update(1)
# ####
    plt.legend(bbox_to_anchor =(1, 1))
    plt.savefig('%s_Head_Tail.png'%targetGeneName, bbox_inches='tight')
    plt.close()