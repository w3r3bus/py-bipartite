#!/usr/bin/env python2.7
import sys


def read_samples(f):
    samples= []
    for line in sys.stdin:
        line.strip("\r\n")
        if (line[0] == '#'): continue
        sample= float(line)
        samples.append(sample)
    return samples


def freq(samples):
    freq_dict= {}
    for r in samples:
        if r not in freq_dict:
            freq_dict[r]= 1
        else:
            freq_dict[r]+= 1
    return freq_dict


def summary(samples):
        mean= 0
        min= None
        max= None
        samples= sorted(samples)
        for r in samples:
            mean+= r
            if (max == None) or (r > max):
                max= r
            if (min == None) or (r < min):
                min= r
        mean/= len(samples)
        size= len(samples)
        if size % 2 == 0:
            median= (samples[size/2-1]+samples[size/2])/2.0
        else:
            median= samples[size/2]
        return {'mean'  : mean,
                'min'   : min,
                'max'   : max,
                'median': median,
                'size'  : size}


def main(mode, argv):
    if mode == "freq":
        samples= read_samples(sys.stdin)
        if len(argv) > 0:
            xmin= int(argv[0])
            samples= [s for s in samples if s >= xmin]
        freq_dict= freq(samples)
        total=0
        for k in sorted(freq_dict.keys()):
            total+= freq_dict[k]
            print "%d\t%.20f\t%d\t%.20f\t%d" % \
                  (k, (1.0*freq_dict[k])/len(samples), freq_dict[k], \
                   (1.0*total)/len(samples), total)

    elif mode == "summary":
        samples= read_samples(sys.stdin)
        if len(argv) > 0:
            xmin= int(argv[0])
            samples= [s for s in samples if s >= xmin]
        summary_dict= summary(samples)
        print "#",
        for k in sorted(summary_dict):
            print "%s\t" %(k),
        print
        for k in sorted(summary_dict):
            print "%s\t" % (summary_dict[k]),
        print

    else:
        sys.stderr.write("unknown mode \"%s\"\n" % (mode))

    return


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.stderr.write("not enough arguments\n")
        sys.exit(-1)
        
    mode= sys.argv[1]
    main(mode, sys.argv[2:])
