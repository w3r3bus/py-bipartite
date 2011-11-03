# ===================================================================
# Pareto distribution pseudo-random number generation
# (Python's random version does not support the xmin parameter)
#
# (C) 2011, Bruno Quoitin (bruno.quoitin@umons.ac.be)
# ===================================================================

import random
import math
import util

class Pareto2:

    def __init__(self, alpha, xmin=1.0):
        self.__alpha= float(alpha)
        self.__xmin= float(xmin)

    def get_alpha(self):
        return self.__alpha

    def get_xmin(self):
        return self.__xmin

    def pdf(self, x):
        if x < self.__xmin:
            return 0
        return (self.__alpha-1.0)/self.__xmin * \
               pow(x / self.__xmin, -self.__alpha)

    def cdf(self, x):
        if x < self.__xmin:
            return 0
        return pow(x / self.__xmin, -self.__alpha)

    def random(self):
        u= random.uniform(1,0)
        return self.__xmin * pow(u, -1.0 / (self.__alpha - 1.0))

    def expectation(self):
        return float('nan')

    def variance(self):
        return float('nan')

    # Compute maximum-likelihood estimator for alpha.
    # It is assumed that the beta parameter is known.
    # Notes:
    #   - if there is a sample with value less than beta, it is ignored
    #   - to avoid numerical cancellation, samples are first sorted
    #     from smallest to largest (smallest value will have the least
    #     log difference with beta and must be summed first).
    @classmethod
    def mle(cls, x, beta=1.0):
        n= len(x)
        y= 0.0
        for i in sorted(x):
            if (i < beta):
                n-= 1
                continue
            y+= (math.log(i)-math.log(beta))
        return 1.0 + (1.0*n)/y


def main(argv):
    if len(argv) < 2:
        sys.stderr.write("not enough arguments\n");
        sys.exit(-1);
    mode= argv[1]

    if mode == "mle":
        print "Reading data from stdin..."
        samples= util.read_samples(sys.stdin)
        print "MLE = %f" % (Pareto2.mle(samples))

    else:

        if len(argv) < 4:
            sys.stderr.write("not enough arguments\n")
            sys.exit(-1)

        alpha= float(argv[2])
        beta= float(argv[3])

        p= Pareto2(alpha, beta)

        if mode == "exp":
            print "%f" % (p.expectation())

        else:

            if len(argv) < 5:
                sys.stderr.write("not enough arguments\n")
                sys.exit(-1)

            size= int(argv[4])

            if mode == "pdf":
                for x in range(size):
                    print "%d\t%.20f" % (x, p.pdf(x))

            elif mode == "cdf":
                for x in range(size):
                    print "%d\t%.20f" % (x, p.cdf(x))

            elif mode == "random":
                freq= {}
                samples= []
                for i in range(size):
                    if (i % 1000) == 0:
                        print >> sys.stderr, "\rGenerating samples %d/%d" % \
                              (i+1, size),
                    r= p.random()
                    samples.append(r)
                print >> sys.stderr, "\rGenerating samples %d/%d" % \
                      (i+1, size)
                mean= 0
                for r in sorted(samples):
                    mean+= r
                mean/= size
                print "# mean = %f" % (mean)
                for r in samples:
                    print "%d" % (r)
            
            else:
                sys.stderr.write("unknown mode \"%s\"\n" % (mode))

if __name__ == "__main__":
    import sys
    main(sys.argv)
