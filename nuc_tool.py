import numpy as np
import math, re, random, csv
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from pylab import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

def SavePlot(y, filename, xtitle = '', ytitle = '',  title = ''):
    plot = Figure(figsize=(12, 3))
    ax =plot.add_subplot(111)
#    plot.grid(True)
    ax.set_title(title)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
#    ax.axis(ymax=1, ymin =-3)
    plot.subplots_adjust(left=0.1, bottom=0.2)

    x = range(len(y))
#    ax.plot(x, nu.log2(Y))
    ax.plot(x, y)
    FigureCanvasAgg(ax.get_figure()).print_figure(filename, dpi=120)
    return

def Save2DArray(d, filename, names =[] ):
    f = open(filename, mode='wb')
    csvwriter = csv.writer(f, dialect='excel-tab')
    for i in transposed2(d):
        csvwriter.writerow(i)
    f.close
    return

def DisplayPlot(y):
    x = range(len(y))
    ax = plt.subplot(111)
    plt.subplots_adjust(bottom=0.2)
    plt.plot(x,y)
    plt.show()
    return

def transposed2(lists, defval=0):
   if not lists: return []
   return map(lambda *row: [elem or defval for elem in row], *lists)

def CleanSeq(dna):
    dna = dna.upper()                 # uppercase
    dna = dna.replace('U','T')        # use T instead of U (in case we get RNA)
    dna = re.sub(r'[^GATC]','', dna)  # remove every character (including whitespace) that is not G, A, T, or C
    return dna

def cshift(l, offset):
    offset %= len(l)
    return np.concatenate((l[-offset:], l[:-offset]))

def SaveSequence(dna, filename):
    f = open(filename, "w")
    print >>f, dna
    f.close()
    return

def base2index(base):
    i=0
    if base=='A':
        i = 0
    if base=='C':
        i = 1
    if base=='G':
        i = 2
    if base=='T':
        i = 3
    return i

def getweight(w,p,b):
    x = np.arange (w)

    AA = 0.25 + b * np.sin(2 * math.pi * x / p)
    AC = 0.25 - b * np.sin(2 * math.pi * x / p) / 3
    AG = AC
    AT = AC

    CA = x*0 + 0.25
    CC = CA
    CG = CA
    CT = CA

    GA = 0.25 + b * np.sin(2 * math.pi * x / p) / 3
    GC = 0.25 - b * np.sin(2 * math.pi * x / p)
    GG = GA
    GT = GA

    TA = 0.25 + b * np.sin(2 * math.pi * x / p)
    TC = 0.25 - b * np.sin(2 * math.pi * x / p)
    TG = TC
    TT = TA
    return [[AA, AC, AG, AT],[CA, CC, CG, CT],[GA, GC, GG, GT],[TA, TC, TG, TT]]

def calcE(seq, w, B, p):
    prob_array = getweight(w, p ,B)
    p_f = []
    p_r = []
    for i in range( len(seq)-w):
        p_s_f = 1.
        p_s_r = 1.
        for s in range(w):
            p_s_f = p_s_f * prob_array[base2index(seq[i+s-1])][base2index(seq[i+s])][s]
            p_s_r = p_s_r * prob_array[3-base2index(seq[i+w-s])][3-base2index(seq[i+w-s-1])][s]
        p_f = np.append(p_f, p_s_f)
        p_r = np.append(p_r, p_s_r)
    p_f = p_f * 4.**w
    p_r = p_r * 4.**w
    p_r = cshift(p_r,-1)
    E = (p_r * np.log(p_r) + p_f * np.log(p_f)) / (p_r + p_f)
    return E

def smooth(x,window_len):
    s=np.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]
    w=np.ones(window_len, 'd')
    y=np.convolve(w / w.sum(), s, mode='valid')
    return y[len(x[window_len-1:0:-1]):len(x[window_len-1:0:-1])+len(x)+1]

def vanderlick(Energy, mu):
    E_out = Energy - mu
    footprint = 147

    forward = np.zeros(len(Energy))
    for i in range( len(Energy) ):
        tmp = sum(forward[ max( i - footprint , 0):i])
        forward[i] = np.exp(E_out[i] - tmp)

    backward = np.zeros(len(Energy))
    r_forward = forward[::-1]
    for i in range( len(Energy) ):
        backward[i] = 1 - sum(r_forward[max(i - footprint , 0):i]* backward[max(i - footprint,0):i])

    P = forward * backward[::-1]
    return P

def CreateDNA(dnalength):
    dna601 = 'ACAGGATGTATATATCTGACACGTGCCTGGAGACTAGGGAGTAATCCCCTTGGCGGTTAAAACGCGGGGGACAGCGCGTACGTGCGTTTAAGCGGTGCTAGAGCTGTCTACGACCAATTGAGCGGCCTCGGCACCGGGATTCTCCAG'
    dna = ''.join(random.choice('ACGT') for x in range((dnalength-len(dna601))/2)) + dna601 +''.join(random.choice('ACGT') for x in range((dnalength-len(dna601))/2))
    dna = CleanSeq(dna)
    return dna

def CalcNucPositions(dna, w =146, mu = -1.5, B=0.2, period=10.1, penalty = None):
    dna = 'A'*73 +dna + 'A'*73
    E_n = calcE(dna, w, B, period)

    E = smooth(E_n,10)
    if penalty is not None:
        E +=penalty
    P = vanderlick(E, mu)
    # P = nu.concatenate(  (nu.zeros(math.ceil(w/2.)), P,  nu.zeros(w/2) ) )
    N = np.convolve(P, np.ones(146), mode ='same')
    return N

if __name__ == '__main__':
    dna = CreateDNA(2000)
    res = CalcNucPositions(dna)
    DisplayPlot(res[2])
    # SavePlot(res[3], 'E:\\tmp\\numpytest.jpg' ,'position (bp)', 'P', 'Nucleosome occupancy')