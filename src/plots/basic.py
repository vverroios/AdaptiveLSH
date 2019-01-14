'''
Created on Jun 21, 2016

@author: vasilis verroios
'''

import sys

import matplotlib.pyplot as plt

from matplotlib.pyplot import xticks, savefig, close, legend, ylabel, xlabel

from numbers import Number

from lsh.probabilities import cliqueRndJaccard

PLOTFOLDER = './plots/'

def plot_pss_jaccard(cliqueSizes):
    sims = [i/100.0 for i in range(100)]
    for q in cliqueSizes:
        plt.plot(sims, 
                 [cliqueRndJaccard(sim, q)
                  for sim in sims])
        
    plt.show()
    
def plot_AndOr_scheme_jaccard(cliqueSizes, z, b):
    sims = [i/100.0 for i in range(100)]
    if isinstance(z, Number):
        z = [z]*len(cliqueSizes)
    for i in range(len(cliqueSizes)):
        pss = [cliqueRndJaccard(sim, 
                                cliqueSizes[i])
               for sim in sims]
        psb = [1.0 - (1.0-(p**z[i]))**(b/z[i])
               for p in pss]
        plt.plot(sims, psb)
            
    plt.show()

def plot_pss_cosine():
    
    plt.plot()
    
def lsh(s,x,y): 
    return 1.0-(1.0-((1.0-s)**x))**y

def plot_paper_blockExample():
    
    MSIZE = 15
    
    xs = [x/100.0 for x in range(101)]
    # set the locations and labels of the xticks
    xticks( [(x)/180.0 for x in [15, 180]], 
            ('15/180', '1.0'), 
            fontsize=20)
    
    plt.plot(xs,[lsh(s,15,140) for s in xs], 
             color="gray",
             linewidth=4.0, linestyle="-.", #marker="s",
             markersize=MSIZE, label="w=15, z=140")
    
    plt.plot(xs,[lsh(s,30,70) for s in xs],
             color="green",
             linewidth=2.0, linestyle="-", #marker="s",
             markersize=MSIZE, label="w=30, z=70")
    
    plt.plot(xs,[lsh(s,60,35) for s in xs],
             color="blue",
             linewidth=2.0, linestyle="--", #marker="s",
             markersize=MSIZE, label="w=60, z=35")
    
    LOCATION = 'upper right'
    lgd = legend(loc=LOCATION, prop={'size':22})
    font = {'family' : 'normal', 'size'   : 26}
    plt.rc('font', **font)
        
    ylabel("Probability of hashing\n to the same bucket", fontsize=26)
#     xlabel("Cosine distance (degrees)", fontsize=26)
    xlabel("Normalized Angle", fontsize=26)
#     ylabel("Probability of the two records\n hashing to the same bucket", fontsize=26)
#     xlabel("Cosine distance between two records (degrees)", fontsize=26)
    
    savefig(PLOTFOLDER + 'blockExampleTest.eps', dpi=300, 
            bbox_extra_artists=(lgd,), bbox_inches='tight')    
    close()    
#     plt.show()

    
def main(argv=None):
    plot_paper_blockExample()
    
if __name__ == '__main__':
    sys.exit(main(sys.argv))

