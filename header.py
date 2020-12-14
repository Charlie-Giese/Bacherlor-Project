from astropy.io import fits
import argparse
import os

parser = argparse.ArgumentParser(description='Plot an image in fits format.')
parser.add_argument("-i", "--inputfile", help="Path of input fits file relative to current working directory")
args = parser.parse_args()
inputfile=args.inputfile


hdu = fits.open(os.getcwd()+inputfile)
hrd=hdu[0].header
print(hrd)
