#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""
Read a grid written by g_thickness and display it as an image.
"""

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FormatStrFormatter
__version__ = '0.1'
__author__ = "Jonathan Barnoud"

color_YM = [
    (0.0, (0.0,0.3,0.9)), #blue
    # (0.2, (0.7,0.8,1.0)),
    (0.5, (1.0,1.0,1.0)), #white
    (0.7, (1.0,1.0,0.0)), #yellow
    (1.0, (0.9,0.1,0.1))] #red
cmap_YM = LinearSegmentedColormap.from_list("own_cmap", color_YM)

def _check_fraction(value):
    """
    Check if a fraction is in the [0; 1[ interval for use in argparse.
    """
    value = float(value)
    if 0 <= value < 1:
        return value
    raise argparse.ArgumentTypeError("Expect a value beween 0 included "
                                     "and 1 excluded.")


def _check_file(value):
    """
    Check if a file for use in argparse.
    """
    if os.path.exists(value):
        return value
    else:
        raise argparse.ArgumentTypeError("Do not find {}".format(value))


def handle_user():
    """
    Get and process user command line inputs.
    """
    parser = argparse.ArgumentParser(
        description='Draw a picture from grid data.')
    parser.add_argument('data_file', metavar='DATA', type=_check_file,
                        help='The data file describing the grid.')
    parser.add_argument('--sampling', '-s', metavar='SAMPING',
                        type=_check_file, default=None,
                        help='The file describing the grid sampling.')
    parser.add_argument('output_file', metavar='OUTPUT', type=str,
                        help='The output image.')
    parser.add_argument('--format', '-f', type=str, default=None,
                        choices=["png", "pdf", "ps", "eps", "svg"],
                        help=('The output image file format, default is '
                              'by extension.'))
    parser.add_argument('--min_fraction', '-m', type=_check_fraction,
                        default=0.5,
                        help=('The minimum fraction of the maximum '
                              'sampling for a cell to be drawn.'))
    parser.add_argument('--nlevels', '-n', type=int, default=None,
                        help=('Number of color levels, default is i'
                              'automatic.'))
    parser.add_argument('--discrete', '-d', action="store_true", default=False,
                        help=('Draw discrete colors instead of an '
                              'interpolated image.'))
    parser.add_argument('--center', '-c', action="store_true", default=False,
                        help=('Put the origin at the center of the grid.'))
    parser.add_argument('--zlim', '-z', type=float, nargs=2,
                        metavar=("ZMIN", "ZMAX"), default=(None, None),
                        help="Value limits for the colorscale")
    parser.add_argument('--xlim', '-x', type=float, nargs=2,
                        metavar=("XMIN", "XMAX"), default=(None, None),
                        help="Picture limits in X")
    parser.add_argument('--ylim', '-y', type=float, nargs=2,
                        metavar=("YMIN", "YMAX"), default=(None, None),
                        help="Picture limits in Y")
    parser.add_argument('--font_size', type=float, default=12, metavar="SIZE",
                        help="Font size, default is 12")
    args = parser.parse_args()
    return args


def read_data(infile_fn):
    """
    Read a grid file and build both a numpy grid and a dictionary of metadata.

    The function gets the path to the file to read. The metadata keys are:

    * xwidth, ywidth: the average dimension of the trajectory box on the axis
      represented by x and y on the plot;
    *  xlabel, ylabel: the axis labels;
    * legend: the colorbar legend.
    """
    template = {"xwidth": float, "ywidth": float,
                "xlabel": str, "ylabel": str,
                "legend": str}
    metadata = {"legend": "thickness (nm)"}
    data = []
    with open(infile_fn) as infile:
        for line in infile:
            splitted = line[:-1].split()
            if line[0] == "@":
                key = splitted[0][1:]
                value = template[key](" ".join(splitted[1:]))
                metadata[key] = value
            elif line[0] == "&":
                break
            elif line[0] != "#":
                data.append(splitted)
    return np.array(data, dtype=float), metadata


def draw(data, metadata,
         nlevels=None, image=True, center=False,
         zlim=(None, None), xlim=(None, None), ylim=(None, None),
         colormap=cmap_YM, colorbar_format="%.2f"):
    """
    Draw the grid with the color levels and the contours.
    """
    if center:
        extent = [-1 * metadata["xwidth"] / 2, metadata["xwidth"] / 2,
                  -1 * metadata["ywidth"] / 2, metadata["ywidth"] / 2]
    else:
        extent = [0, metadata["xwidth"], 0, metadata["ywidth"]]
    if image:
        fill = plt.imshow(data.transpose(), extent=extent,
                          interpolation="bicubic",cmap=colormap,
                          vmin=zlim[0], vmax=zlim[1])
    elif nlevels is None:
        fill = plt.contourf(data.transpose(),  extent=extent, origin="image", cmap=colormap)
    else:
        fill = plt.contourf(data.transpose(), nlevels - 1, extent=extent,
                            origin="image",cmap=colormap)

    if nlevels is None:
        plt.contour(data.transpose(), colors='k', extent=extent,
                    origin="image")
    else:
        plt.contour(data.transpose(), nlevels - 1, colors='k', extent=extent,
                    origin="image")
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xlabel(metadata["xlabel"])
    plt.ylabel(metadata["ylabel"])

   ##############
    # Hide all elements except the figure
    plt.xticks([])
    plt.yticks([])
    plt.xlabel('')
    plt.ylabel('')
    #plt.gca().set_frame_on(False)  # Hide the frame
    
    plt.draw()
    ############    
    
#    colorbar = plt.colorbar(fill)
#    colorbar.set_label(metadata["legend"])
#        # Format the colorbar ticks
#    colorbar.formatter = FormatStrFormatter(colorbar_format)
#    colorbar.update_ticks()  # Update colorbar with formatted ticks
    
    


def cache_low_sampling(bad, metadata, center=False):
    """
    Draw white cells  on the plot.

    The bad array has to contain 0 in the cell to draw and numpy.NaN in the
    others.
    """
    if center:
        extent = [-1 * metadata["xwidth"] / 2, metadata["xwidth"] / 2,
                  -1 * metadata["ywidth"] / 2, metadata["ywidth"] / 2]
    else:
        extent = [0, metadata["xwidth"], 0, metadata["ywidth"]]
    plt.imshow(bad, extent=extent, cmap=cm.gray_r, interpolation="nearest")


def filter_sampling(data, sampling, minimum_fraction):
    """
    Remove the cells where the sampling is bellog a given fraction of the
    maximum sampling. The removed cells contain numpy.NaN.

    Returns the filtered array and a cache array with 0 in the cells of low
    sampling and numpy.NaN in the others.
    """
    filtered = np.array(data)
    low_sampling = (sampling < minimum_fraction * sampling.max())
    filtered[low_sampling] = np.NaN
    bad = np.zeros(data.shape)
    bad *= np.NaN
    bad[low_sampling] = 0
    return filtered, bad


def plot_grid(data, sampling, metadata, minimum_fraction=0.3,
              nlevels=None, image=True, center=False,
              zlim=(None, None), xlim=(None, None), ylim=(None, None)):
    """
    Draw the data and cache the cells with low sampling if asked. To ignore
    sampling, set the "sampling" argument to None.
    """
    if not sampling is None:
        filtered, bad = filter_sampling(data, sampling, minimum_fraction)
    else:
        filtered = data
    draw(filtered, metadata,
         nlevels=nlevels, image=image, center=center,
         zlim=zlim, xlim=xlim, ylim=ylim)
    if (image and not sampling is None
            and (bad.shape[0] * bad.shape[1] - np.sum(np.isnan(bad)) != 0)):
        cache_low_sampling(bad, metadata, center=center)


def save(output, file_format=None):
    """
    Save the picture. If the format is not specified then matplotlib will
    guess it from extension.
    """
    if file_format is None:
        file_format = os.path.splitext(output)[-1][1:]
    plt.savefig(output, format=file_format, bbox_inches="tight")


def main():
    """
    Run the program.
    """
    # Read user input
    args = handle_user()
    # Read files
    if not args.sampling is None:
        sampling, metadata = read_data(args.sampling)
    else:
        sampling = None
    data, metadata = read_data(args.data_file)
    # Do the plot
    rc('font', **{"size": args.font_size})
    plot_grid(data, sampling, metadata, args.min_fraction,
              args.nlevels, not args.discrete, args.center,
              args.zlim, args.xlim, args.ylim)
    # Save the figure
    save(args.output_file, args.format)


if __name__ == "__main__":
    main()
