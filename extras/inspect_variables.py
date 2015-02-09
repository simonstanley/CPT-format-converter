#!/usr/local/sci/bin/python2.7
"""
Simple script to print variables with a NetCDF, pp or grib file.

"""
import iris
import argparse

def main(filename):
    data = iris.load(filename)
    for data_cube in data:
        print data_cube.name()

parser = argparse.ArgumentParser(description='List all the variable names '\
                                 'within a NetCDF, pp or grib file.')
parser.add_argument('filename', help='The full file paths to the file to '\
                    'be inspected.')

args = parser.parse_args()
main(args.filename)
