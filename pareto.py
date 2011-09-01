# ===================================================================
# Pareto distribution pseudo-random number generation
# (Python's random version does not support the scale parameter)
#
# (C) 2011, Bruno Quoitin (bruno.quoitin@umons.ac.be)
# ===================================================================

import random
import math

class Pareto:

    def __init__(self, alpha, beta=1.0):
        self.__alpha= float(alpha)
        self.__beta= float(beta)

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
            if (u >= 0) and (u < 1):
                break
        return 1.0 * self.__beta / pow(1.0-u, 1.0 / self.__alpha)

    def expectation(self):
        if self.__alpha == 1:
            return float('nan')
        return (1.0 * self.__alpha * self.__beta) \
               / (self.__alpha - 1)

    def variance(self):
        if (self.__alpha == 1) or (self.__alpha == 2):
            return float('nan')
        return (1.0 * pow(self.__beta, 2) * self.__alpha) \
               / (pow(self.__alpha - 1, 2) * (self.__alpha - 2))

    # Compute maximum-likelihood estimator for alpha.
    # It is assumed that the beta parameter is known.
    # Notes:
    #   - if there is a 0-valued sample, NaN is returned
    #   - if there is a sample with value less than beta, it is ignored
    @classmethod
    def mle(cls, x, beta=1.0):
        n= len(x)
        y= 0
        for i in x:
            if (i == 0):
                return float('nan');
            if (i < beta):
                n-= 1
                continue
            y+= (math.log(i)-math.log(beta))
        return (1.0*n)/y


# -----[ Test application ]------------------------------------------
# For each exponent 'a', will generate 'num_samples' samples.
# Will compute and print frequency distribution for each set of
# samples.
# -------------------------------------------------------------------
def test_pareto_sample(range_a, range_b, num_samples, rounding="ROUND"):
    freq= {}
    count= 0
    total= num_samples*len(range_a)*len(range_b)
    for a in range_a:
        for b in range_b:
            data= []
            p= Pareto(a, b)
            print "# exponent=%.20f" % (p.get_alpha())
            print "# beta=%.20f" % (p.get_beta())
            print "#   expectation=%.20f" % (p.expectation())
            print "#   variance=%.20f" % (p.variance())
            mean= 0
            for i in range(num_samples):
                if count % 1000 == 0:
                    print >> sys.stderr, "\r%.2f %%" % (100.0*count/total),
                r= p.random()
                mean+= r
                data.append(r)
                if rounding == "TRUNC":
                    r= int(r)
                elif rounding == "ROUND":
                    r= round(r)
                else:
                    raise "Error: unsupported rounding"
                id= "%f;%f" % (a, b)
                if r not in freq:
                    freq[r]= {}
                    freq[r][id]= 1
                else:
                    if id not in freq[r]:
                        freq[r][id]= 1
                    else:
                        freq[r][id]+= 1
                count+= 1
            print "#   mean=%.20f" % (mean/num_samples)
            print >> sys.stderr, "\r%.2f %%" % (100.0*count/total)
            print "#   number of samples generated=%d" % (len(data))
            # Compute MLE for generated data
            print "#   maximum likelihood estimate=%.20f" % (Pareto.mle(data, b))
            
                    
    for i in sorted(freq.keys()):
        print "%.20f" % (i),
        for a in range_a:
            for b in range_b:
                id= "%f;%f" % (a, b)
                if id in freq[i]:
                    print "\t%.20f" % ((1.0*freq[i][id])/num_samples),
                else:
                    print "\t0",
        print
    return
    

def main():
    #num_samples= 1000000
    num_samples=1000
    #test_pareto_sample(range(3), range(2), num_samples)
    test_pareto_sample([1.5, 2.0, 2.5], [1.0, 2.0], num_samples, rounding="ROUND")
    

if __name__ == "__main__":
    import sys
    main()
