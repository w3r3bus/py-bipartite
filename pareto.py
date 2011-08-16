# ===================================================================
# Pareto distribution pseudo-random number generation
# (Python's random version does not support the scale parameter)
#
# (C) 2011, Bruno Quoitin (bruno.quoitin@umons.ac.be)
# ===================================================================

import random


class Pareto:

    def __init__(self, alpha, beta=1.0):
        self.__alpha= alpha
        self.__beta= beta

    def get_alpha(self):
        return self.__alpha

    def get_beta(self):
        return self.__beta

    def pdf(self, x):
        return (1.0 * self.__alpha) \
               * pow(self.__beta, self.__alpha) \
               / pow(x, self.__alpha+1)

    def cdf(self, x):
        return 1-pow(1.0 * self.__beta / x, self.__alpha)

    def random(self):
        while True:
            u= random.random()
            if (u > 0) and (u < 1):
                break
        return 1.0 * self.__beta / pow(1-u, 1.0 / self.__alpha)

    def expectation(self):
        if self.__alpha == 1:
            return float('nan')
        return (1.0 * self.__alpha * self.__beta) \
               / (self.__alpha - 1)

    def variance(self):
        return (1.0 * pow(self.__beta, 2) * self.__alpha) \
               / (pow(self.__alpha - 1, 2) * (self.__alpha - 2))

    # Compute maximum-likelihood estimator for alpha.
    # It is assumed that the beta parameter is known.
    @classmethod
    def mle(cls, x, beta=1.0):
        n= len(x)
        y= 0
        # note assumption that MLE for beta = 1
        for i in x:
            y+= (math.log(i)-math.log(beta))
        return (1.0*n)/y


# -----[ Test application ]------------------------------------------
# For each exponent 'a', will generate 'num_samples' samples.
# Will compute and print frequency distribution for each set of
# samples.
# -------------------------------------------------------------------
def test_pareto_sample(range_a, b, num_samples, rounding="ROUND"):
    freq= {}
    count= 0
    total= num_samples*len(range_a)
    for a in range_a:
        p= Pareto(a+1, b)
        print "# exponent=%.20f" % (p.get_alpha())
        print "#   expectation=%.20f" % (p.expectation())
        mean= 0
        for i in range(num_samples):
            if count % 1000 == 0:
                print >> sys.stderr, "\r%.2f %%" % (100.0*count/total),
            if rounding == "TRUNC":
                r= int(p.random())
            elif rounding == "ROUND":
                r= round(p.random())
            else:
                raise "Error: unsupported rounding"
            mean+= r
            if r not in freq:
                freq[r]= {}
                freq[r][a]= 1
            else:
                if a not in freq[r]:
                    freq[r][a]= 1
                else:
                    freq[r][a]+= 1
            count+= 1
        print "#   mean=%.20f" % (mean/num_samples)
        print >> sys.stderr, "\r%.2f %%" % (100.0*count/total)
                    
    for i in sorted(freq.keys()):
        print "%.20f" % (i),
        for a in range_a:
            if a in freq[i]:
                print "\t%.20f" % ((1.0*freq[i][a])/num_samples),
            else:
                print "\t0",
        print
    return
    

def main():
    num_samples= 1000000
    b= 1
    #test_pareto_sample(range(3), b, num_samples)
    test_pareto_sample(range(3), b, num_samples, rounding="TRUNC")
    

if __name__ == "__main__":
    import sys
    main()
