import numpy as np
import math
# NCD function from eq. 1 in Bitarello, B. D. et al. Signatures of long-term balancing selection in human genomes. Genome Biol. Evol. 10, 939–955 (2018).
def ncd(MAF, TF):
    ncd = math.sqrt(sum(TF-MAF)**2/len(MAF)) 
    return(ncd)
