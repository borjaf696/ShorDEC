#Borja :)
import sys, re
import altair as alt
import pandas as pd
import numpy as np

LOWERBOUND = 0
UPPEROUND = 490

def __histogramToDf(histogram):
    x = np.arange(len(histogram))
    return pd.DataFrame({'x':x,'histogram':histogram})

def __buildReverseCumulative(histogram):
    size, cumm = len(histogram), 0
    cummHisto = [0]*size
    for i in range(size, 0, -1):
        value = histogram[i-1]
        cummHisto[i-1] = value + cumm
        cumm += value
    return cummHisto

def __buildHistogram(pathHisto, upper = UPPEROUND, lower = LOWERBOUND):
    histogram, numLine = [],0
    with open(pathHisto, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if re.search('^[0-9]+',line) and ',' not in line:
                if numLine >= lower and numLine <= upper:
                    histogram.append(int(line))
                    numLine += 1
                else:
                    numLine+=1
    return histogram

def __meadianFilter(histogram):
    '''
    Create a median filter in the histogram
    :param histogram:
    :return:
    '''
    window = 12
    for i, value in enumerate(histogram[window:(len(histogram)-window-1)]):
        if i <= window:
            continue
        vect = np.asarray(histogram[i-window:i+window])
        histogram[i] = np.median(vect)

def __studyHistogram(histogram):
    window, check= 1, False
    for i,val in enumerate(histogram):
        if i == 0:
            continue
        if histogram[i] < histogram[i-1]:
            check = True
        if check:
            vectTend = [False]*window
            selection = min(window,len(histogram)-i)
            for j in range(min(window, len(histogram)-i-1)):
                vectTend[j] = (histogram[i+j+1] > histogram[i])
            sumTend = 0
            for b in vectTend:
                if b:
                    sumTend += 1
                    if sumTend > (selection / 2):
                        return i+LOWERBOUND
    return 0

def __exportChart(df,nameFile, type = 'line'):
    chart = alt.Chart(df)
    if type == 'line':
        chart = chart.mark_line().encode(
            x='x',
            y='histogram')
    if type == 'bar':
        chart = chart.mark_bar().encode(
            x='x',
            y='histogram'
        )
    chart.save(nameFile+'.html')

if __name__=='__main__':
    print('Stats script')

    cmd = sys.argv[1]
    if cmd == 'histogram':
        output_dir = 'stats/'
        pathHisto = sys.argv[2]
        histogram = __buildHistogram(pathHisto)
        print('Processing histogram from: ',pathHisto)
        df = __histogramToDf(histogram)
        __exportChart(df, output_dir+'histogram','bar')
        print('Process reverse cumulative histogram: ', pathHisto)
        __exportChart(__histogramToDf(__buildReverseCumulative(histogram)),output_dir+'cumRevHistogram')
    if cmd == 'analyse':
        pathHisto = sys.argv[2]
        output_dir = 'stats/'
        histogram = __buildHistogram(pathHisto, 1000)
        __meadianFilter(histogram)
        __exportChart(__histogramToDf(histogram),output_dir+'smooth_histogram','bar')
        print('Threshold: ', __studyHistogram(histogram))

