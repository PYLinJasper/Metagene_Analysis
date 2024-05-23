import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt

def checkNull(codon, RA, RC, GeneLen):
    # place check forward
    if codon <= RA:
        leftBound = 1 - codon
        for timesCount in range(0, RA - codon + 1):
            RC[timesCount] = None
    else:
        leftBound = -RA
        
    # place check backward
    if codon + RA > GeneLen:
        rightBound = GeneLen - codon
        for timesCount in range(GeneLen - codon + RA + 1, RA * 2 + 1):
            RC[timesCount] = None
    else:
        rightBound = RA
    
    # bound -> Coodinate Transform
    bound = []
    bound.append(leftBound)
    bound.append(rightBound)
    return(bound)

def countRC(SS, SE, ERC, CS, bound, RA, RC):
    # SS = Seq Start, SE = Seq End, ERC = Even Read Count, CS = Codon Start/Stop, CT = Coodinate Transform
    SSHCT = SS - CS
    SEHCT = SE - CS - 1

    startSeq = 0
    closeSeq = 0
    overlap = 0
    regionStart = bound[0]
    regionEnd = bound[1]
    if SEHCT >= regionStart:
        if SSHCT <= regionEnd:
            overlap = 1
            if SSHCT >= regionStart:
                startSeq = SSHCT
            else:
                startSeq = regionStart
                    
            if SEHCT >= regionEnd:
                closeSeq = regionEnd
            else:
                closeSeq = SEHCT
        else:
            overlap = 0
    else:
        overlap = 0
            
    startSeq += RA
    closeSeq += RA
    if overlap:
        for timesCount in range(startSeq, closeSeq + 1):
            RC[timesCount] += ERC

# RA = Require Amount
def Start_Stop(sourceFN, mRNALibFN, startRA, stopRA, startPlotOFN, stopPlotOFN):
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
    for timesCount in range(-startRA, startRA + 1):
        rACol.append('%s'%timesCount)
    # print(rACol)

    startPlot = pd.DataFrame()
    startPlot['Gene name'] = targetGene['Gene name']
    startPlot.set_index('Gene name', inplace = True)
    startPlot = startPlot.reindex(columns = rACol)
    # print(startPlot)

    rACol = []
    for timesCount in range(-stopRA, stopRA + 1):
        rACol.append('%s'%timesCount)

    stopPlot = pd.DataFrame()
    stopPlot['Gene name'] = targetGene['Gene name']
    stopPlot.set_index('Gene name', inplace = True)
    stopPlot = stopPlot.reindex(columns = rACol)
    # print(stopPlot)

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
        # codon info
        nowGeneCodonStart = int(nowTargetGeneInfo["CDS start"])
        nowGeneCodonStop = int(nowTargetGeneInfo["CDS end"])
        # Zerolize
        thisStartRC = [0] * (startRA * 2 + 1)
        thisStopRC = [0] * (stopRA * 2 + 1)

        nowGene.reset_index(inplace = True)

        boundStart = checkNull(nowGeneCodonStart, startRA, thisStartRC, nowGeneSeqLength)
        boundStop = checkNull(nowGeneCodonStop, stopRA, thisStopRC, nowGeneSeqLength)

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

            countRC(thisSeqStart, thisSeqEnd, thisGeneEvenRC, nowGeneCodonStart, boundStart, startRA, thisStartRC)
            countRC(thisSeqStart, thisSeqEnd, thisGeneEvenRC, nowGeneCodonStop, boundStop, stopRA, thisStopRC)

        startPlot.loc[GeneName] = thisStartRC
        stopPlot.loc[GeneName] = thisStopRC
        progess.update(1)

        # break point check
        # if wbreak > 130:
        #     break
        # else:
        #     wbreak += 1

    startPlot.reset_index(inplace=True)
    startPlot.drop("Gene name", axis=1, inplace = True) 
    startPlot.set_index('-100', inplace = True)
    startPlot.to_csv(startPlotOFN)
    # print(startPlot)

    stopPlot.reset_index(inplace=True)
    stopPlot.drop("Gene name", axis=1, inplace = True) 
    stopPlot.set_index('-100', inplace = True)
    stopPlot.to_csv(stopPlotOFN)
    # print(stopPlot)

    HTPlot = []
    HTPlot.append(startPlot)
    HTPlot.append(stopPlot)
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
    for binOrder in range(-100, 101):
        HD_1.loc[HD_1['%s'%binOrder] == 'N', '%s'%binOrder] = -1
        HD_2.loc[HD_2['%s'%binOrder] == 'N', '%s'%binOrder] = -1

    HD_1 = HD_1.astype(float)
    HD_2 = HD_2.astype(float)
    for binOrder in range(-100, 101):
        HD_1.loc[HD_1['%s'%binOrder] == -1, '%s'%binOrder] = None
        HD_2.loc[HD_2['%s'%binOrder] == -1, '%s'%binOrder] = None

    HD_N_1 = HD_1.copy()
    HD_N_2 = HD_2.copy()
    HD_N_1['accCounts'] = HD_N_1.sum(axis = 1)
    HD_N_2['accCounts'] = HD_N_2.sum(axis = 1)
    HD_N_factor_1 = HD_N_1['accCounts'].sum()
    HD_N_factor_2 = HD_N_2['accCounts'].sum()

    # print(HD_N_1)

    progess = tqdm(total = 201, desc = 'File Preprocess')
    for binOrder in range(-100, 101):
        HD_N_1['%s'%binOrder] = HD_N_1['%s'%binOrder]/HD_N_factor_1
        HD_N_2['%s'%binOrder] = HD_N_2['%s'%binOrder]/HD_N_factor_2
        progess.update(1)

    # print(HD_N_1)
    perBinRatio = [[]] * 201
    perBinMV = [[]] * 201
    perBinSTE_1 = [[]] * 201
    perBinSTE_2 = [[]] * 201
    perBinSTEN_1 = [[]] * 201
    perBinSTEN_2 = [[]] * 201
    perBinMVN = [[]] * 201
    sitePosition = [[]] * 201
    progess = tqdm(total = 201, desc = fileN)
    for binOrder in range(0, 201):
        # Calculate Ratio
        realOrder = binOrder - 100
        # position
        position_1 = HD_1.loc[HD_1['%s'%realOrder] != None]['%s'%realOrder].count()
        # print(position_1)
        position_2 = HD_2.loc[HD_2['%s'%realOrder] != None]['%s'%realOrder].count()
        
        # site
        Site_1 = HD_1.loc[HD_1['%s'%realOrder] > 0]['%s'%realOrder].count()
        Site_2 = HD_2.loc[HD_2['%s'%realOrder] > 0]['%s'%realOrder].count()

        # position_1 = 0
        # position_2 = 0
        # Site_1 = 0
        # Site_2 = 0

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
        MV_1 = HD_1['%s'%realOrder].sum()/len(HD_1)
        STE_1 = np.std(HD_1['%s'%realOrder], ddof=1) / np.sqrt(np.size(HD_1['%s'%realOrder]))
        wholePair = []
        wholePair.append(MV_1 - STE_1)
        # wholePair.append(MV_1)
        wholePair.append(MV_1 + STE_1)
        perBinSTE_1[binOrder] = wholePair
        MV_2 = HD_2['%s'%realOrder].sum()/len(HD_2)
        STE_2 = np.std(HD_2['%s'%realOrder], ddof=1) / np.sqrt(np.size(HD_2['%s'%realOrder]))
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
        HD_N_M_1 = HD_N_1['%s'%realOrder].sum()/len(HD_N_1)
        STE_1N = np.std(HD_N_1['%s'%realOrder], ddof=1) / np.sqrt(np.size(HD_N_1['%s'%realOrder]))
        Pair_N_1 = []
        Pair_N_1.append(HD_N_M_1 - STE_1N)
        Pair_N_1.append(HD_N_M_1 + STE_1N)

        HD_N_M_2 = HD_N_2['%s'%realOrder].sum()/len(HD_N_2)
        STE_2N = np.std(HD_N_2['%s'%realOrder], ddof=1) / np.sqrt(np.size(HD_N_2['%s'%realOrder]))
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
    # print(sitePosition)

    sipoChart = sitePosition.plot(title='%s target_HEAD100'%targetGeneName, legend = False,
    xlabel='region(%)', ylabel='# of mRNA', figsize=(18, 14), ax=axes[0, posN], xlim = (xRangeStart, xRangeEnd))
    # plt.xlim(-100, 100)
    min = sitePosition.min().values.min()
    max = sitePosition.max().values.max()
    sipoChart.vlines(100, min, max, linestyle = '--', color = 'black')

    # Ratio
    perBinRatio = pd.DataFrame(perBinRatio, columns=['%s mRNA site ratio (site/position)'%IFN1_category, '%s mRNA site ratio (site/position)'%IFN2_category])

    ratioChart = perBinRatio.plot(legend = False,
    xlabel='region(%)', ylabel='ratio mRNA', figsize=(18, 14), ax=axes[1, posN], xlim = (xRangeStart, xRangeEnd))
    # plt.xlim(-100, 100)
    min = perBinRatio.min().values.min()
    max = perBinRatio.max().values.max()
    ratioChart.vlines(100, min, max, linestyle = '--', color = 'black')

    # Read counts

    perBinMV = pd.DataFrame(perBinMV, columns=['%s target AVG (+/-STE)'%IFN1_category, '%s target AVG (+/-STE)'%IFN2_category])
    # print(perBinMV)

    ste1_1 = perBinMV.plot(xlabel='region(%)', ylabel='read counts', linewidth = 1, ax=axes[2, posN], xlim = (xRangeStart, xRangeEnd), legend = False)

    xLab = 'region(%)'
    yLab = 'read counts'
    subPlotX = 2
    subPlotY = posN
    xScale = 201
    tarColor = 'royalblue'
    rCC_STE_1 = drawSTE(perBinSTE_1, tarColor, xLab, yLab, subPlotX, subPlotY, xScale, xRangeStart, xRangeEnd, axes)
    tarColor = 'darkorange'
    rCC_STE_2 = drawSTE(perBinSTE_2, tarColor, xLab, yLab, subPlotX, subPlotY, xScale, xRangeStart, xRangeEnd, axes)
    min = perBinMV.min().values.min()
    max = perBinMV.max().values.max()
    ste1_1.vlines(100, min, max, linestyle = '--', color = 'black')

    # Read counts Distribution
    perBinMVN = pd.DataFrame(perBinMVN, columns=['%s target AVG (+/-STE)'%IFN1_category, '%s target AVG (+/-STE)'%IFN2_category])
    # print(perBinMVN.max().values[0])

    rCountNChart = perBinMVN.plot(xlabel='region(%)', ylabel='read counts Distribution', linewidth = 1, ax=axes[3, posN], xlim = (xRangeStart, xRangeEnd), legend = False)

    xLab = 'region(%)'
    yLab = 'read counts Distribution'
    subPlotX = 3
    subPlotY = posN
    xScale = 201
    tarColor = 'royalblue'
    rCCN_STE_1 = drawSTE(perBinSTEN_1, tarColor, xLab, yLab, subPlotX, subPlotY, xScale, xRangeStart, xRangeEnd, axes)
    tarColor = 'darkorange'
    rCCN_STE_2 = drawSTE(perBinSTEN_2, tarColor, xLab, yLab, subPlotX, subPlotY, xScale, xRangeStart, xRangeEnd, axes)
    min = perBinMVN.min().values.min()
    max = perBinMVN.max().values.max()
    rCountNChart.vlines(100, min, max, linestyle = '--', color = 'black')


    allPlot = []
    allPlot.append(sipoChart)
    allPlot.append(ratioChart)
    allPlot.append(ste1_1)
    allPlot.append(rCountNChart)
    return(allPlot)


def plotStartStop(targetGeneTranscript, targetGeneName, IFN1_Start, IFN1_Stop, IFN1_Den, IFN2_Start, IFN2_Stop, IFN2_Den, IFN1_category, IFN2_category):
    startPlot1 = pd.read_csv(IFN1_Start)
    stopPlot1 = pd.read_csv(IFN1_Stop)
    startPlot2 = pd.read_csv(IFN2_Start)
    stopPlot2 = pd.read_csv(IFN2_Stop)
    densityPlot1 = pd.read_csv(IFN1_Den)
    densityPlot2 = pd.read_csv(IFN2_Den)
    tCIndex1 = preN2C(densityPlot1)
    tCIndex1.set_index(['Gene name'], inplace = True)
    tCIndex2 = preN2C(densityPlot2)
    tCIndex2.set_index(['Gene name'], inplace = True)


    GeneTranscript = pd.read_csv(targetGeneTranscript)

    # tMG = target MetaGenePlot
    HD_1 = preGene(GeneTranscript, startPlot1.copy(), tCIndex1, 'Gene name')
    HD_2 = preGene(GeneTranscript, startPlot2.copy(), tCIndex2, 'Gene name')

    fig, axes = plt.subplots(nrows=4, ncols=2)

    posN = 0
    xRangeStart = 0
    xRangeEnd = 201
    fileN = targetGeneName + ' START'
    draw4Plot(HD_1, HD_2, targetGeneName, IFN1_category, IFN2_category, axes, posN, xRangeStart, xRangeEnd, fileN)

# ####
    HD_1 = preGene(GeneTranscript, stopPlot1.copy(), tCIndex1, 'Gene name')
    HD_2 = preGene(GeneTranscript, stopPlot2.copy(), tCIndex2, 'Gene name')

    posN = 1
    xRangeStart = 0
    xRangeEnd = 201
    fileN = targetGeneName + ' STOP'
    allPolt = draw4Plot(HD_1, HD_2, targetGeneName, IFN1_category, IFN2_category, axes, posN, xRangeStart, xRangeEnd, fileN)

    progess = tqdm(total = 4, desc = 'Adjust Plot and Save Plot')
    for index in range(0, 4):
        allPolt[index].legend(bbox_to_anchor =(1, 1))
        progess.update(1)
# ####
    plt.legend(bbox_to_anchor =(1, 1))
    plt.savefig('%s_Start_Stop.png'%targetGeneName, bbox_inches='tight')
    plt.close()
