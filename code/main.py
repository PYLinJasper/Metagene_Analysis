import step6_7_Adjust_Even
import step8_9_Density_Metagene
import step10_11_Head_Tail
import step12_13_Start_Stop
import step14_15_SC_Fold_Change

# HW1 Section 1 produce right raw data

# NFile = pd.read_csv('mapping_WT__HRDEIP_norm.csv')
# target = NFile.drop('evenly_rc', 1)
# target = target.set_index('input_id', 1)
# target.to_csv('mapping_WT__HRDEIP.csv')
# print(target)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#### step6_7_Adjust_Even
# targetFN = 'final_result_WT1_HRDEIP.csv'
# scalar = 2.9683
# outputFN = 'final_result_WT1_HRDEIP_Adjust_Even.csv'
# outputFN_PD = step6_7_Adjust_Even.Adjust_Even(targetFN, scalar, outputFN)
# print(outputFN_PD)

# targetFN = 'final_result_CSR_KO1_HRDEIP.csv'
# scalar = 6.70947
# outputFN = 'final_result_CSR_KO1_HRDEIP_Adjust_Even.csv'
# outputFN_PD = step6_7_Adjust_Even.Adjust_Even(targetFN, scalar, outputFN)
# print(outputFN_PD)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#### step8_9_desity_MetaGene_fileProduce

# sourceFN = 'final_result_WT1_HRDEIP_Adjust_Even.csv'
# mRNALibFN = 'mRNA_WS275.csv'
# metaGeneOFN = 'WT1_HRDEIP_all_mRNA_tool2_100.csv'
# densityOFN = 'WT1_HRDEIP_all_mRNA_tool1.csv'
# targetGene1 = step8_9_Density_Metagene.Density_Metagene(sourceFN, mRNALibFN, metaGeneOFN, densityOFN)
# print('Density')
# print(targetGene1[0])
# print('MetaGene')
# print(targetGene1[1])

# sourceFN = 'final_result_CSR_KO1_HRDEIP_Adjust_Even.csv'
# metaGeneOFN = 'CSR_KO1_HRDEIP_all_mRNA_tool2_100.csv'
# densityOFN = 'CSR_KO1_HRDEIP_all_mRNA_tool1.csv'
# targetGene2 = step8_9_Density_Metagene.Density_Metagene(sourceFN, mRNALibFN, metaGeneOFN, densityOFN)
# print('Density')
# print(targetGene2[0])
# print('MetaGene')
# print(targetGene2[1])

#### + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

#### Plot Density
targetGeneTranscript = 'CSR-1_target_transcripts.csv'
targetGeneName = 'CSR-1'
IFN1_SD = 'WT1_HRDEIP_all_mRNA_tool1.csv'
IFN1_category = 'WT_HRDEIP'
IFN2_SD = 'CSR_KO1_HRDEIP_all_mRNA_tool1.csv'
IFN2_category = 'CSR_KO_HRDEIP'
step8_9_Density_Metagene.plot2Density(targetGeneTranscript, targetGeneName, IFN1_SD, IFN1_category, IFN2_SD, IFN2_category)

targetGeneTranscript = 'WAGO-1_target_transcripts.csv'
targetGeneName = 'WAGO-1'
step8_9_Density_Metagene.plot2Density(targetGeneTranscript, targetGeneName, IFN1_SD, IFN1_category, IFN2_SD, IFN2_category)

#### + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

#### Plot MetaGene

targetGeneTranscriptN = 'CSR-1_target_transcripts.csv'
targetGeneN = 'CSR-1'
IFN1_Meta = 'WT1_HRDEIP_all_mRNA_tool2_100.csv'
IFN1_Den = 'WT1_HRDEIP_all_mRNA_tool1.csv'
IFN2_Meta = 'CSR_KO1_HRDEIP_all_mRNA_tool2_100.csv'
IFN2_Den = 'CSR_KO1_HRDEIP_all_mRNA_tool1.csv'
IFN2_cat = 'CSR_KO_HRDEIP'
IFN1_cat = 'WT_HRDEIP'
step8_9_Density_Metagene.plotMetaGene(targetGeneTranscriptN, targetGeneN, IFN1_Meta, IFN1_Den, IFN2_Meta, IFN2_Den, IFN1_cat, IFN2_cat)

targetGeneTranscriptN = 'WAGO-1_target_transcripts.csv'
targetGeneN = 'WAGO-1'
step8_9_Density_Metagene.plotMetaGene(targetGeneTranscriptN, targetGeneN, IFN1_Meta, IFN1_Den, IFN2_Meta, IFN2_Den, IFN1_cat, IFN2_cat)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#### step10_11_Head_Tail_fileProduce
#### RA = Require Amount
# sourceFN = 'final_result_WT1_HRDEIP_Adjust_Even.csv'
# mRNALibFN = 'mRNA_WS275.csv'
# headRA = 100
# tailRA = 100
# headPlotOFN = 'WT1_HRDEIP_all_mRNA_tool3_HEAD_100.csv'
# tailPlotOFN = 'WT1_HRDEIP_all_mRNA_tool3_TAIL_100.csv'
# targetGene3 = step10_11_Head_Tail.Head_Tail(sourceFN, mRNALibFN, headRA, tailRA, headPlotOFN, tailPlotOFN)
# print('Head')
# print(targetGene3[0])
# print('Tail')
# print(targetGene3[1])

# sourceFN = 'final_result_CSR_KO1_HRDEIP_Adjust_Even.csv'
# headRA = 100
# tailRA = 100
# headPlotOFN = 'CSR_KO1_HRDEIP_all_mRNA_tool3_HEAD_100.csv'
# tailPlotOFN = 'CSR_KO1_HRDEIP_all_mRNA_tool3_TAIL_100.csv'
# targetGene4 = step10_11_Head_Tail.Head_Tail(sourceFN, mRNALibFN, headRA, tailRA, headPlotOFN, tailPlotOFN)
# print('Head')
# print(targetGene4[0])
# print('Tail')
# print(targetGene4[1])

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

#### Plot Head Tail # #
targetGeneTranscriptN = 'CSR-1_target_transcripts.csv'
targetGeneN = 'CSR-1'
IFN1_Head = 'WT1_HRDEIP_all_mRNA_tool3_HEAD_100.csv'
IFN1_Tail = 'WT1_HRDEIP_all_mRNA_tool3_TAIL_100.csv'
IFN1_Den = 'WT1_HRDEIP_all_mRNA_tool1.csv'
IFN2_Head = 'CSR_KO1_HRDEIP_all_mRNA_tool3_HEAD_100.csv'
IFN2_Tail = 'CSR_KO1_HRDEIP_all_mRNA_tool3_TAIL_100.csv'
IFN2_Den = 'CSR_KO1_HRDEIP_all_mRNA_tool1.csv'
IFN2_cat = 'CSR_KO_HRDEIP'
IFN1_cat = 'WT_HRDEIP'
step10_11_Head_Tail.plotHeadTail(targetGeneTranscriptN, targetGeneN, IFN1_Head, IFN1_Tail, IFN1_Den, IFN2_Head, IFN2_Tail, IFN2_Den, IFN1_cat, IFN2_cat)

targetGeneTranscriptN = 'WAGO-1_target_transcripts.csv'
targetGeneN = 'WAGO-1'
step10_11_Head_Tail.plotHeadTail(targetGeneTranscriptN, targetGeneN, IFN1_Head, IFN1_Tail, IFN1_Den, IFN2_Head, IFN2_Tail, IFN2_Den, IFN1_cat, IFN2_cat)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#### step12_13_Start_Stop_fileProduce
#### RA = Require Amount
# sourceFN = 'final_result_WT1_HRDEIP_Adjust_Even.csv'
# mRNALibFN = 'mRNA_WS275.csv'
# startRA = 100
# stopRA = 100
# startPlotOFN = 'WT1_HRDEIP_all_mRNA_tool3_START_100.csv'
# stopPlotOFN = 'WT1_HRDEIP_all_mRNA_tool3_STOP_100.csv'
# targetGene5 = step12_13_Start_Stop.Start_Stop(sourceFN, mRNALibFN, startRA, stopRA, startPlotOFN, stopPlotOFN)
# print('Start')
# print(targetGene5[0])
# print('Stop')
# print(targetGene5[1])

# sourceFN = 'final_result_CSR_KO1_HRDEIP_Adjust_Even.csv'
# startRA = 100
# stopRA = 100
# startPlotOFN = 'CSR_KO1_HRDEIP_all_mRNA_tool3_START_100.csv'
# stopPlotOFN = 'CSR_KO1_HRDEIP_all_mRNA_tool3_STOP_100.csv'
# targetGene6 = step12_13_Start_Stop.Start_Stop(sourceFN, mRNALibFN, startRA, stopRA, startPlotOFN, stopPlotOFN)
# print('Start')
# print(targetGene6[0])
# print('Stop')
# print(targetGene6[1])

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

## Plot start stop
targetGeneTranscriptN = 'CSR-1_target_transcripts.csv'
targetGeneN = 'CSR-1'
IFN1_Head = 'WT1_HRDEIP_all_mRNA_tool3_START_100.csv'
IFN1_Tail = 'WT1_HRDEIP_all_mRNA_tool3_STOP_100.csv'
IFN1_Den = 'WT1_HRDEIP_all_mRNA_tool1.csv'
IFN2_Head = 'CSR_KO1_HRDEIP_all_mRNA_tool3_START_100.csv'
IFN2_Tail = 'CSR_KO1_HRDEIP_all_mRNA_tool3_STOP_100.csv'
IFN2_Den = 'CSR_KO1_HRDEIP_all_mRNA_tool1.csv'
IFN2_cat = 'CSR_KO_HRDEIP'
IFN1_cat = 'WT_HRDEIP'
step12_13_Start_Stop.plotStartStop(targetGeneTranscriptN, targetGeneN, IFN1_Head, IFN1_Tail, IFN1_Den, IFN2_Head, IFN2_Tail, IFN2_Den, IFN1_cat, IFN2_cat)

targetGeneTranscriptN = 'WAGO-1_target_transcripts.csv'
targetGeneN = 'WAGO-1'
step12_13_Start_Stop.plotStartStop(targetGeneTranscriptN, targetGeneN, IFN1_Head, IFN1_Tail, IFN1_Den, IFN2_Head, IFN2_Tail, IFN2_Den, IFN1_cat, IFN2_cat)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#### fold_change_FileProduce
# firstDIFN = 'WT1_HRDEIP_all_mRNA_tool1.csv'
# firstScOFN = 'sc_WT.csv'
# secondDIFN = 'CSR_KO1_HRDEIP_all_mRNA_tool1.csv'
# secondScOFN = 'sc_CSR_KO.csv'
# foldChangeOFN = 'fold_change.csv'
# # log2(secondSc / firstSc)
# targetGene7 = step14_15_SC_Fold_Change.SC_Fold_Change(firstDIFN, firstScOFN, secondDIFN, secondScOFN, foldChangeOFN)
# print(targetGene7[0])
# print(targetGene7[1])
# print(targetGene7[2])

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

#### Plot Fold Change

targetGeneTranscript = 'CSR-1_target_transcripts.csv'
targetGeneName = 'CSR-1'
n_1 = 'sc_WT.csv'
n_2 = 'sc_CSR_KO.csv'
plotN_1 = 'log2(read count) WT_HRDEIP'
plotN_2 = 'log2(read count) CSR_CO_HRDEIP'
IFN1_cat = 'WT_HRDEIP'
IFN2_cat = 'CSR_KO_HRDEIP'
mergeIndexName = 'Gene name'
foldChangeN = 'fold_change.csv'

step14_15_SC_Fold_Change.foldChange_Plot(targetGeneName, targetGeneTranscript, n_1, n_2, plotN_1, plotN_2, IFN1_cat, IFN2_cat, mergeIndexName, foldChangeN)

targetGeneTranscript = 'WAGO-1_target_transcripts.csv'
targetGeneName = 'WAGO-1'

step14_15_SC_Fold_Change.foldChange_Plot(targetGeneName, targetGeneTranscript, n_1, n_2, plotN_1, plotN_2, IFN1_cat, IFN2_cat, mergeIndexName, foldChangeN)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Plot 4 Density Comparison
targetGeneTranscript_1 = 'CSR-1_target_transcripts.csv'
targetGeneName_1 = 'CSR-1'
IFN1_SD_1 = 'WT1_HRDEIP_all_mRNA_tool1.csv'
IFN1_category = 'WT_HRDEIP'
IFN2_SD_1 = 'CSR_KO1_HRDEIP_all_mRNA_tool1.csv'
IFN2_category = 'CSR_KO_HRDEIP'

targetGeneTranscript_2 = 'WAGO-1_target_transcripts.csv'
targetGeneName_2 = 'WAGO-1'
IFN1_SD_2 = 'WT1_HRDEIP_all_mRNA_tool1.csv'
IFN2_SD_2 = 'CSR_KO1_HRDEIP_all_mRNA_tool1.csv'
step8_9_Density_Metagene.plot4Density(targetGeneTranscript_1, targetGeneTranscript_2, targetGeneName_1, targetGeneName_2, IFN1_SD_1, IFN1_SD_2, IFN1_category, IFN2_SD_1, IFN2_SD_2, IFN2_category)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Plot scatter boxplot comparison

targetGeneTranscript_1 = 'CSR-1_target_transcripts.csv'
targetGeneName_1 = 'CSR-1'
n_1 = 'sc_WT.csv'
n_2 = 'sc_CSR_KO.csv'
plotN_1 = 'log2(read count) WT_HRDEIP'
plotN_2 = 'log2(read count) CSR_CO_HRDEIP'
IFN1_cat = 'WT_HRDEIP'
IFN2_cat = 'CSR_KO_HRDEIP'
mergeIndexName = 'Gene name'
foldChangeN = 'fold_change.csv'
targetGeneTranscript_2 = 'WAGO-1_target_transcripts.csv'
targetGeneName_2 = 'WAGO-1'


step14_15_SC_Fold_Change.comparison_Plot(targetGeneName_1, targetGeneName_2, targetGeneTranscript_1, targetGeneTranscript_2, n_1, n_2, plotN_1, plotN_2, IFN1_cat, IFN2_cat, mergeIndexName, foldChangeN)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #