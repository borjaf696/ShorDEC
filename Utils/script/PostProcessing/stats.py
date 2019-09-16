#Borja :)
import sys, re
import altair as alt
import pandas as pd
import numpy as np

LOWERBOUND = 10
UPPEROUND = 200

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
    x = np.arange(len(histogram))
    return pd.DataFrame({'x':x,'histogram':histogram})

def __exportChart(df,nameFile, type = 'line'):
    chart = alt.Chart(df)
    if type == 'line':
        chart = chart.mark_line().encode(
            x='x',
            y='histogram')

    chart.save(nameFile+'.html')

if __name__=='__main__':
    print('Stats script')

    cmd = sys.argv[1]
    if cmd == 'histogram':
        pathHisto = sys.argv[2]
        print('Processing histogram from: ',pathHisto)
        df = __buildHistogram(pathHisto)
        __exportChart(df, 'histogram')
