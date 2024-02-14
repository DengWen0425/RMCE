#! /usr/bin/env python3

from tgraph import TGraph
import sys, time
import optparse
import resource
import csv

def main():

    usage = "Usage: ./temporal-cliques.py [options]"
    fmt = optparse.IndentedHelpFormatter(max_help_position=50, width=100)
    parser = optparse.OptionParser(usage=usage, formatter=fmt)
    parser.add_option("-f", "--file", type="string", dest="filename", default="",
                      help="Input file")
    parser.add_option("-d", "--delta", type="int", dest="delta", default=0,
                      help="Delta, default is 0.")
    parser.add_option("-p", "--pivoting", type="int", dest="pivoting_method", default=0,
                      help="Set the algorithm used for pivoting by setting pivoting_method to one of the following. "
                      "0: Do not use pivoting. "
                      "1: Pick a single arbitrary pivot. "
                      "2: Pick a single pivot maximizing the number of removed elements. "
                      "3: Pick an arbitrary maximal set of pivots. "
                      "4: Pick a set of pivots greedily according to the maximum number of removed elements. "
                      "5: Pick a set of pivots maximizing the number of removed elements.")

    (options, args) = parser.parse_args()
    delta = options.delta
    
    graph = TGraph(delta, options.pivoting_method)
    if len(options.filename) == 0:
        print("error: no input file specified")
    sys.stdout.write("Reading file: %s\n" % options.filename)
    graph.from_file(options.filename)

    # Start execution
    sys.stdout.write("Starting execution\n")
    t1 = time.time()
    graph.run()
    t2 = time.time()
    sys.stdout.write("Accumulated time to calculate Cliques: %s seconds \n" % (t2 - t1))
    sys.stdout.write("# maximal Delta-Cliques: %d \n" % (graph.clique_count))
if __name__ == "__main__":
    sys.exit(main())
