from matplotlib.colors import BoundaryNorm
from matplotlib.colors import LightSource

from preprocessing.ASCtoSTLdict import DEM_file
from preprocessing.ASCtoSTLdict import xc
from preprocessing.ASCtoSTLdict import yc

import re
import vtk
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from linecache import getline
import pandas as pd
from scipy.interpolate import griddata
import netCDF4

sys.path.insert(0, './preprocessing')

cell_size_map = 25.0

map_size = 500.0


# Print iterations progress
def printProgressBar(iteration,
                     total,
                     prefix='',
                     suffix='',
                     decimals=1,
                     bar_length=20):
    """
    Call in a loop to create terminal progress bar

    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : number of decimals in percent complete (Int)
        bar_length  - Optional  : character length of bar (Int)
    """
    str_format = "{0:." + str(decimals) + "f}"
    percents = str_format.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = '█' * filled_length + '-' * (bar_length - filled_length)

    sys.stdout.write('\r%s |%s| %s%s %s' %
                     (prefix, bar, percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()


def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


def readerVTK(filename, dispersed_phase):

    reader = vtk.vtkDataSetReader()

    reader.SetFileName(filename)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.ReadAllTensorsOn()
    reader.ReadAllFieldsOn()
    reader.Update()

    data = reader.GetOutput()

    npoints = data.GetNumberOfPoints()

    d = data.GetPointData()

    phaseName = "alpha." + dispersed_phase

    alpha_array = d.GetArray(phaseName)

    position = np.zeros((npoints, 3))
    alpha_particles = np.zeros(npoints)

    for n in range(npoints):

        position[n, :] = data.GetPoint(n)
        alpha_particles[n] = np.array(alpha_array.GetTuple(n))

    return position, alpha_particles


def readPhaseProperties():

    print('Reading file: phaseProperties')

    word1 = 'phases'
    word2 = 'continuousPhase'

    with open('./constant/phaseProperties', 'r') as fp:
        # read all lines using readline()
        lines = fp.readlines()
        for row in lines:
            # check if string present on a current line

            line = row.lstrip()
            find1 = row.find(word1)
            if find1 != -1:

                if line.split()[0] == word1:
                    print(word1 + ' exists in file')
                    matches = re.findall(r'\((.*?)\)', line)
                    phases = matches[0].split()
                    print(word1 + ' = ', phases)

            find2 = row.find(word2)
            if find2 != -1:
                print(word2 + ' exists in file')
                continuousPhase = (row.split()[-1].replace(';', '')).strip()
                print(word2 + ' = ', continuousPhase)

    dispersedPhases = [phase for phase in phases if phase != continuousPhase]

    return dispersedPhases


def readControlDict():

    print('Reading control file: controldict')

    word1Found = False
    word2Found = False

    with open('./system/controlDict', 'r') as fp:
        # read all lines using readline()
        lines = fp.readlines()
        for row in lines:
            # check if string present on a current line

            word1 = 'writeControl'
            find1 = row.find(word1)
            if find1 != -1:
                print(word1 + ' exists in file')
                print('line Number:', lines.index(row))
                writeControl = (row.split()[-1].replace(';', '')).strip()
                print(word1 + ' = ', writeControl)
                word1Found = True

            word2 = 'writeInterval'
            find2 = row.find(word2)
            if find2 != -1:
                print(word2 + ' exists in file')
                print('line Number:', lines.index(row))
                writeInterval = float((row.split()[-1].replace(';',
                                                               '')).strip())
                print(word2 + ' = ', writeInterval)
                word2Found = True

    dt = 1.0

    if not word1Found:

        print(word1, ' not found')

    elif writeControl == 'adjustableRunTime':

        if not word2Found:

            print(word2, ' not found')

        else:

            dt = writeInterval

    return dt


def readASC(DEM_file):

    print('Reading DEM file: ' + DEM_file)
    # Parse the topography header
    hdr = [getline(DEM_file, i) for i in range(1, 7)]
    values = [float(h.split()[-1].strip()) for h in hdr]
    cols, rows, lx, ly, cell, nd = values
    cols = int(cols)
    rows = int(rows)

    xs_DEM = lx + 0.5 * cell + np.linspace(0, (cols - 1) * cell, cols)
    ys_DEM = ly + 0.5 * cell + np.linspace(0, (rows - 1) * cell, rows)

    extent = lx, lx + cols * cell, ly, ly + rows * cell

    # Load the topography into a numpy array
    DEM = pd.read_table(DEM_file,
                        delim_whitespace=True,
                        header=None,
                        skiprows=6).astype(float).values
    DEM = np.flipud(DEM)
    DEM[DEM == nd] = 0.0

    xinit = xs_DEM
    yinit = ys_DEM

    xmin = np.amin(xinit)
    xmax = np.amax(xinit)

    print('xmin,xmax', xmin, xmax)

    ymin = np.amin(yinit)
    ymax = np.amax(yinit)

    print('ymin,ymax', ymin, ymax)

    Xinit, Yinit = np.meshgrid(xinit, yinit)
    Zinit = DEM

    return Xinit, Yinit, Zinit, cell, extent


def main():

    dt = readControlDict()

    print('dt = ', dt)

    dispersedPhases = readPhaseProperties()

    print('dispersedPhases = ', dispersedPhases)

    current_dir = os.getcwd()
    current_dir_name = current_dir.split('/')[-1]

    check_dir = current_dir + '/VTK'

    VTKexists = os.path.exists(check_dir)

    foamCommand = "foamToVTK -fields '()' -noInternal " + \
        "-noFaceZones -excludePatches '(atm top terrain_in terrain_out)'"

    print(foamCommand)

    if VTKexists:

        print("")
        print(check_dir + " found")

    else:

        os.system(foamCommand)

    print('DEM_file', DEM_file)
    filename = './preprocessing/' + DEM_file
    Xinit, Yinit, Zinit, cell, extent = readASC(filename)

    print(extent)

    ls = LightSource(azdeg=45, altdeg=45)

    working_dir = current_dir + '/VTK/terrain_out'
    files = os.listdir(working_dir)

    file_idx = []
    for filename in files:

        file_idx.append(float(filename.split('_')[-1].rsplit('.', 1)[0]))

    sorted_files = [files[i] for i in np.argsort(file_idx)]

    full_filename = working_dir + '/' + sorted_files[-1]

    position, alpha_max = readerVTK(full_filename, dispersedPhases[0])

    position[:, 0] += xc
    position[:, 1] += yc

    points = position[:, 0:2]

    x_map_min = np.amin(position[:, 0])
    x_map_max = np.amax(position[:, 0])
    y_map_min = np.amin(position[:, 1])
    y_map_max = np.amax(position[:, 1])

    x_map = np.arange(x_map_min, x_map_max, step=cell_size_map)
    y_map = np.arange(y_map_min, y_map_max, step=cell_size_map)

    extent_map = x_map_min - 0.5 * cell_size_map, x_map_max + 0.5 * \
        cell_size_map, y_map_min - 0.5 * cell_size_map, \
        y_map_max + 0.5 * cell_size_map

    print('extent_map', extent_map)

    xx, yy = np.meshgrid(x_map, y_map)

    ncfilename = current_dir_name + '_solidMap.nc'

    ncfile = netCDF4.Dataset(ncfilename, mode='w', format='NETCDF4_CLASSIC')

    ncfile.createDimension('x', xx.shape[1])
    ncfile.createDimension('y', xx.shape[0])
    ncfile.createDimension('time', None)

    ncfile.title = current_dir_name + ' output'

    ncfile.Conventions = "CF-1.0"
    ncfile.subtitle = "My model data subtitle"
    ncfile.anything = "write anything"

    x = ncfile.createVariable('x', np.float64, ('x', ))
    x.long_name = 'x dim'
    x.units = 'meters'

    y = ncfile.createVariable('y', np.float64, ('y', ))
    y.long_name = 'y dim'
    y.units = 'meters'

    time = ncfile.createVariable('time', np.float64, ('time', ))
    time.long_name = 'Time'
    time.units = 'seconds'

    # note: unlimited dimension is leftmost
    alphaP = ncfile.createVariable('alphaP', np.float64, ('time', 'y', 'x'))
    alphaP.standard_name = 'alfa particles'
    alphaP.units = 'volume fraction'

    x[:] = xx[0, :]
    y[:] = yy[:, 0]

    zz = np.zeros_like(xx)

    for i, filename in enumerate(sorted_files[:]):

        printProgressBar(i, len(sorted_files) - 1)

        time[i] = dt * i
        full_filename = working_dir + '/' + filename
        # print(filename)

        for dispersed_phase in dispersedPhases:

            alpha_particles = readerVTK(full_filename,dispersed_phase)[1]
            zz += griddata(points, alpha_particles, (xx, yy), method='nearest')

        alphaP[i, :, :] = zz

        alpha_max = np.maximum(alpha_max, alpha_particles)

    print(ncfile)

    ncfile.close()
    print('Dataset is closed!')

    print(np.amin(alpha_max))
    print(np.amax(alpha_max))

    values = alpha_max[:]

    zz = griddata(points, values, (xx, yy), method='nearest')

    zz[zz < 1.e-6] = np.nan

    # print(np.amin(zz))
    # print(np.amax(zz))

    xmin = np.amin(Xinit) - 0.5 * cell
    xmax = np.amax(Xinit) + 0.5 * cell

    ymin = np.amin(Yinit) - 0.5 * cell
    ymax = np.amax(Yinit) + 0.5 * cell

    fig, ax = plt.subplots()

    ax.imshow(ls.hillshade(np.flipud(Zinit), vert_exag=1.0, dx=cell, dy=cell),
              cmap='gray',
              extent=extent)

    ax.set_aspect('equal', 'box')

    levels = np.linspace(-6.0, 1.0, 11)
    label_str = 'Probabilty [0;1]'

    cmap = plt.get_cmap('viridis')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    p1 = ax.imshow(np.flipud(np.log10(zz)),
                   cmap=cmap,
                   norm=norm,
                   interpolation='nearest',
                   extent=extent_map,
                   alpha=0.5)

    ax.set_xlim(xc - map_size, xc + map_size)
    ax.set_ylim(yc - map_size, yc + map_size)

    ax.tick_params(axis='both', which='major', labelsize=8)

    clb = plt.colorbar(p1)

    label_str = 'Log10 of Volume fraction'
    clb.set_label(label_str, labelpad=5, y=0.5, rotation=90)

    title = 'Solid particles'
    plt.title(title)

    png_file = current_dir_name + '_solidMap.png'

    plt.savefig(png_file, dpi=200)

    nd = -9999

    zz[np.isnan(zz)] = nd

    asc_file = current_dir_name + '_solidMap.asc'

    header = "ncols     %s\n" % xx.shape[1]
    header += "nrows    %s\n" % xx.shape[0]
    header += "xllcorner " + str(x_map_min - 0.5 * cell_size_map) + "\n"
    header += "yllcorner " + str(y_map_min - 0.5 * cell_size_map) + "\n"
    header += "cellsize " + str(cell_size_map) + "\n"
    header += "NODATA_value " + str(nd)

    np.savetxt(asc_file,
               np.flipud(zz),
               header=header,
               fmt='%1.5f',
               comments='')

    plt.show()
    plt.close(fig)


if __name__ == '__main__':

    main()
