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
parser.add_argument('config_file', default=None, nargs='?', help='Optional. '\
                    'A .txt configuration file can be used to overwrite or '\
                    'add cpt tags.')
parser.add_argument('-s', action='store_true', default=False, dest='simple',
                    help='Optional. Use simple formatting.')

args = parser.parse_args()

filename = args.filename
savename = args.savename
if args.variable == None:
    variable    = None
    config_file = None
else:
    if args.variable[-4:] == ".txt":
        config_file = args.variable
        if args.config_file is not None:
            # Arguments have been provided backwards.
            if args.config_file[-4:] == ".txt":
                raise UserWarning("Too many files given.")
            else:
                variable = args.config_file
        else:
            variable = None
    else:
        variable = args.variable
        if args.config_file is not None:
            if args.config_file[-4:] != ".txt":
                raise UserWarning("Too many arguments given.")
        config_file = args.config_file



savepath = cpt_converter(filename, savename, variable, config_file, args.simple)
print 'Converted data saved to: %s' % savepath
