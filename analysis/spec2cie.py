#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt
import colour

def read_spectra(filename):
    data = np.loadtxt(filename)
    wavelengths = data[:, 0]
    values = data[:, 1]
    return dict(zip(wavelengths, values))

def calculate_xyz(spectrum):
    cmfs = colour.MSDS_CMFS['CIE 1931 2 Degree Standard Observer']
    illuminant = colour.SDS_ILLUMINANTS['D65']
    
    sd = colour.SpectralDistribution(spectrum, name='Sample')
    sd = sd.copy().align(cmfs.shape)  # Align to cmf wavelengths

    XYZ = colour.sd_to_XYZ(sd, cmfs=cmfs, illuminant=illuminant)
    xy = colour.XYZ_to_xy(XYZ)
    return XYZ, xy

def plot_chromaticity_diagram(xy):
    colour.plotting.plot_chromaticity_diagram_CIE1931(show=False)
    plt.plot(xy[0], xy[1], 'o', color='red', markersize=10)
    plt.text(xy[0] + 0.01, xy[1], f"({xy[0]:.3f}, {xy[1]:.3f})", fontsize=12)
    plt.title("CIE 1931 Chromaticity Diagram")
    plt.tight_layout()
    plt.show()

def main():
    if len(sys.argv) != 2:
        print("Usage: python command.py file.spectra")
        sys.exit(1)

    filename = sys.argv[1]
    spectrum = read_spectra(filename)
    XYZ, xy = calculate_xyz(spectrum)

    print("CIE 1931 XYZ: ", XYZ)
    print("Chromaticity xy: ", xy)
    plot_chromaticity_diagram(xy)

if __name__ == "__main__":
    main()

