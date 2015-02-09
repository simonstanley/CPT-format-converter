#!/usr/local/sci/bin/python2.7
"""
Simple script to merge NetCDF, pp or grib files into a single file.

"""
import iris
import argparse
import os

def check_filepath(filepath):
    """
    Check if file exists, if so, ask the user if they wish to overwrite it.

    Args:

    * filepath: string
        Full path to the file to check

    """
    if os.path.isfile(filepath):
        response = raw_input("This file already exists, would you like to "\
                             "overwrite?")
        if response.lower() in ['y', 'yes']:
            os.remove(filepath)
        elif response.lower() in ['n', 'no']:
            raise UserWarning('Program stopped.')
        else:
            print "Please respond yes or no"
            check_filepath(filepath)

def main(filenames, savename, variable=None):
    check_filepath(savename)
    data = iris.load(filenames, variable)
    iris.save(data, savename)
    print "Done"


parser = argparse.ArgumentParser(description='Merge netcdf, grib or pp into '\
                                             'a single file')
parser.add_argument('filenames', help='The full file paths to the data to '\
                    'be converted.', nargs='+')
parser.add_argument('savename', help='The full file path to where the '\
                    'merged file is saved. Use the extention to determine '\
                    'output format. E.g. /path/to/example_file.nc will save '\
                    'the merged file as a NetCDF file. .pp and .grib are also'\
                    ' valid formats.')

args = parser.parse_args()
main(args.filenames, args.savename)
