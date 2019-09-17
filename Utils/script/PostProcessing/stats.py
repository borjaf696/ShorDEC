#Borja :)
import sys, re
import altair as alt
import pandas as pd
import numpy as np

LOWERBOUND = 20
UPPEROUND = 300

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

def __buildHistogram(pathHisto):
    histogram, numLine = [],0
    with open(pathHisto, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if re.search('^[0-9]+',line) and ',' not in line:
                if numLine > LOWERBOUND and numLine < UPPEROUND:
                    histogram.append(int(line))
                    numLine += 1
                else:
                    numLine+=1
    return histogram

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

