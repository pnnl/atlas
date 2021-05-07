#!/usr/bin/env python
import os
import sys
import re



def main(jgi_file):
    # parsing input
    header = {}
    col2keep = ['contigName', 'contigLen', 'totalAvgDepth']
    with open(jgi_file) as inF:
        for i,line in enumerate(inF):
            line = line.rstrip().split('\t')
            if i == 0:
                header = {x:ii for ii,x in enumerate(line)}
                col2keep += [x for x in line if x.endswith('.bam')]
                print('\t'.join(col2keep))
                continue
            elif line[0] == '':
                continue
            # contig ID
            contig = line[header['contigName']]
            # collect per-sample info
            out = []
            for col in col2keep:
                out.append(line[header[col]])
            print('\t'.join(out))

if __name__ == '__main__':

    if 'snakemake' in globals():

        with open(snakemake.log[0],'w') as log:
            sys.stderr = log

            with open(snakemake.output[0],'w') as outf:
                sys.stdout= outf

                main(snakemake.input[0])


    else:


        import argparse
        import logging
        logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
        class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                              argparse.RawDescriptionHelpFormatter):
            pass

        desc = 'Converting jgi_summarize_bam_contig_depths output to format used by VAMB'
        epi = """DESCRIPTION:
        Output format: contigName<tab>contigLen<tab>totalAvgDepth<tab>SAMPLE1.sort.bam<tab>Sample2.sort.bam<tab>...
        Output written to STDOUT
        """
        parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                         formatter_class=CustomFormatter)
        argparse.ArgumentDefaultsHelpFormatter
        parser.add_argument('jgi_file', metavar='jgi_file', type=str,
                            help='jgi_summarize_bam_contig_depths output table')
        parser.add_argument('--version', action='version', version='0.0.1')


        args = parser.parse_args()
        main(args.jgi_file)
