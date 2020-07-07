import os, sys
import shapefile
import numpy 
import scipy


def process(shapein, csvout):
    sf = shapefile.Reader(shapein)
    shapes = sf.shapes()

    with open(csvout,'w') as f:
        for i in range(0,len(shapes)):
            f.write("NaN"+' '+"NaN"+'\n')
            for point in shapes[i].points:
                f.write(str(point[0])+' '+str(point[1])+'\n')

        f.write("NaN"+' '+"NaN"+'\n')
    f.close()
    print("Done! Check for file: ",csvout)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        prog='shape2coast.py', usage='%(prog)s shapein csvout')
    #(i.e python3 ../../ana/pitt_ana.py '*.txt' 1 6 )
    ## main arguments
	
    parser.add_argument('shapein', type=str,help='input file name + path')
    parser.add_argument('csvout', type=str,help='output file name + path')
    args = parser.parse_args()
	

	### PRINT IT ALL
    print('input name : %s' % (args.shapein))
    print('output name : %s' % (args.csvout))
	
    process(args.shapein, args.csvout)
