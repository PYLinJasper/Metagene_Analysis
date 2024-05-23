from cProfile import label
import pandas as pd
import numpy as np
from tqdm import tqdm
import time

def Adjust_Even(targetFN, scalar, outputFN):
    originF = pd.read_csv(targetFN)
    # print(originF)
    # Adjust
    # originF['read_count'] = originF['read_count']*scalar

    # inputID accumulate occur times
    count = originF['input_id'].value_counts()
    count = pd.DataFrame(count, columns=['input_id'])
    count.reset_index(inplace = True)
    count.rename(columns={'index': 'input_id', 'input_id': 'occur_count'}, inplace=True)
    # print(count)

    # merge Count -> prepare to even
    mergeIndexName = 'input_id'
    originF = pd.merge(originF, count, on=[mergeIndexName])
    # print(originF)

    # even
    normF = originF.copy()
    normF['even_read_count'] = normF["read_count"] / normF['occur_count']
    normF['even_read_count'] = normF['even_read_count']*scalar
    # print(normF)

    normF = normF.set_index('input_id')
    # print(normF)
    normF.to_csv(outputFN)

    return(normF)