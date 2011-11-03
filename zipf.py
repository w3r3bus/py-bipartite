#!/usr/bin/env python
import math
import random
import sys
import util

__nth_hn_cache= {}

# -----[ nth_hn ]----------------------------------------------------
# Returns the Nth generalized harmonic number Hn,s
# A cache is used to hold already computed harmonic numbers
def nth_hn(alpha, n, xmin=1):
    n= int(n)
    # Lookup in cache
    hn_id= "%d_%f_%d" % (n, alpha, xmin)
    if not(hn_id in __nth_hn_cache):
        # Summation starts from the smallest components (i.e.
        # those with the largest exponent) to minimize
        # floating-point cancellation
        hn= 0.0
        while (n >= xmin):
            hn= hn + 1.0/pow(n, alpha)
            n-= 1
            # Store in cache
        __nth_hn_cache[hn_id]= hn
    else:
        hn= __nth_hn_cache[hn_id]
    return hn


class Zipf:

    def __init__(self, alpha, N, xmin=1):
        self.__alpha= float(alpha)
        self.__N= int(N)
        self.__xmin= int(xmin)
        return

    # -----[ pmf ]---------------------------------------------------
    # Zipf probability mass function: returns the frequency of the
    # k^th element
    def pmf(self, k):
        k= int(k)
        if (k < self.__xmin) or (k > self.__N):
            return 0
        return 1.0/(pow(k, self.__alpha) * \
                    nth_hn(self.__alpha, self.__N, self.__xmin))

    def pdf(self, k):
        return self.pmf(k)

    # -----[ cdf ]---------------------------------------------------
    # Zipf cumulative distribution function: returns the cumulative
    # frequency of the elements < k
    def cdf(self, k):
        k= int(k)
        # Summation starts from the smallest components (i.e. the
        # frequency of highest ranks) to minimize floating-point
        # cancellation
        cdf= 0
        while (k > 0):
            cdf+= self.pmf(k)
            k-= 1
        return cdf


    # -----[ random ]------------------------------------------------
    # Generate values of a random variable having a Zipf probability mass
    # function. The inverse transform method is used.
    def random(self):
        u= random.uniform(0, 1)

        # Perform inverse transform
        k= self.__xmin
        sum= 0.0
        while k <= self.__N:
            p= self.pmf(k)
            sum= sum + p
            if sum >= u:
                return k
            k+= 1
        return None

    def expectation(self):
        return nth_hn(self.__alpha-1.0, self.__N, self.__xmin) / \
               nth_hn(self.__alpha, self.__N, self.__xmin)

    #def variance(self):
    #    return float('NaN')

    @classmethod
    def mle(cls, samples, xmin=1.0):
        # Discard samples lower than xmin
        samples= [s for s in samples if s >= xmin]
        y= sum(math.log(s/(xmin-(1.0/2.0))) for s in sorted(samples))
        return 1.0 + len(samples) * (1.0/y)


def main(mode, argv):

    if mode == "mle":

        if len(argv) > 0:
            xmin= int(argv[0])
        else:
            xmin= 1
        
        samples= util.read_samples(sys.stdin)
        print "xmin = %d" % (xmin)
        print "MLE = %f" % (Zipf.mle(samples, xmin))
        return 0

    else:

        if len(argv) < 3:
            sys.stderr.write("not enough arguments\n")
            return -1

        alpha= float(argv[0])
        N= int(argv[1])
        xmin= int(argv[2])

        z= Zipf(alpha, N, xmin)

        if mode == "exp":
            print "%f" % (z.expectation())
            return 0

        elif mode == "pdf":
            for i in range(xmin, N+1):
                print "%d\t%.20f" % (i, z.pmf(i))
            return 0

        elif mode == "cdf":
            for i in range(xmin, N+1):
                print "%d\t%.20f" % (i, z.cdf(i))
            return 0

        elif mode == "random":

            if len(argv) < 4:
                sys.stderr.write("not enough arguments\n")
                return -1

            size= int(argv[3])
            freq= {}
            samples= []
            for i in range(size):
                if (i % 1000) == 0:
                    print >> sys.stderr, "\rGenerating samples %d/%d" % \
                          (i+1, size),
                r= z.random()
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
            #print >> sys.stderr, "\rfinished :-)"
        
        else:
            sys.stderr.write("invalid mode \"%s\"\n" % (mode));
        
    return


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write("not enough arguments\n");
        sys.exit(-1);

    mode= sys.argv[1]
    
    main(mode, sys.argv[2:])
