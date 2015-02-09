#!/usr/local/sci/bin/python2.7
"""
Script to run CPT converter from command line.

"""
import argparse
from cpt_converter_module import cpt_converter

parser = argparse.ArgumentParser(description='Convert netcdf, grib or pp into'\
                                 ' CPT format. Output is a text file.')
parser.add_argument('filename', help='The full file path to the data to '\
                    'be converted.')
parser.add_argument('savename', help='The full file path to where the '\
                    'converted data is saved.')
parser.add_argument('variable', default=None, nargs='?', help='Optional. '\
                    'Convert only the given meteorological variable.')
parser.add_argument('-s', action='store_true', default=False, dest='simple',
                    help='Optional. Use simple formatting.')

args = parser.parse_args()
savepath = cpt_converter(args.filename, args.savename, args.variable, args.simple)
print 'Converted data saved to: %s' % savepath
