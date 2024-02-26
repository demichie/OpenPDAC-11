# from matplotlib.colors import BoundaryNorm
from matplotlib.colors import LightSource

from preprocessing.ASCtoSTLdict import xc
from preprocessing.ASCtoSTLdict import yc
from preprocessing.ASCtoSTLdict import DEM_file

import vtk
import os
import sys
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from linecache import getline
import pandas as pd

sys.path.insert(0, './preprocessing')

toll = 1.0
step_dens = 50.0

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


def readerVTK(filename):

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

    U_array = d.GetArray("U")
    d_array = d.GetArray("d")
    rho_array = d.GetArray("rho")

    origId_array = d.GetArray("origId")

    position = np.zeros((npoints, 3))
    U = np.zeros((npoints, 3))
    d = np.zeros(npoints)
    rho = np.zeros(npoints)
    origId = np.zeros(npoints).astype(int)

    for n in range(npoints):

        position[n, :] = data.GetPoint(n)
        U[n, :] = U_array.GetTuple(n)
        d[n] = np.array(d_array.GetTuple(n))
        rho[n] = np.array(rho_array.GetTuple(n))
        origId[n] = np.array(origId_array.GetTuple(n))

    return origId, d, U, position, rho


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

    print('dt= ', dt)

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

    ls = LightSource(azdeg=45, altdeg=45)

    working_dir = current_dir + '/VTK/lagrangian/cloud'
    files = os.listdir(working_dir)
    n_files = len(files)
    file_idx = []
    for filename in files:
        # print(filename)
        file_idx.append(float(filename.split('_')[1].rsplit('.', 1)[0]))

    sorted_files = [files[i] for i in np.argsort(file_idx)]
    sorted_times = [dt * i for i in range(len(file_idx))]

    # print('Times read',sorted_times)

    full_filename = working_dir + '/' + sorted_files[-1]

    print('file', full_filename)

    # get Ids and diameters of ballistic from last output
    origIdLast, d = readerVTK(full_filename)[0:2]

    diams = np.unique(d)

    nballistics = d.shape[0]
    n_times = n_files
    print('nballistics', nballistics)

    # matr is an array with one row for each timestep:
    # column 0,1,2: ballistic coordiantes
    # column 3,4,5: ballistic velocity components
    # column 6: ballistic velocity magnitude
    matr = np.zeros((n_times, 7, nballistics))

    for i, filename in enumerate(sorted_files[:]):

        printProgressBar(i, len(sorted_files) - 1)

        full_filename = working_dir + '/' + filename
        origId, d, U, position, rho = readerVTK(full_filename)

        # sort Ids to match those from last output
        sorter = np.argsort(origId)
        subset = sorter[np.searchsorted(origId, origIdLast, sorter=sorter)]

        for k in range(nballistics):
            matr[i, 0:3, k] = position[subset[k], :]
            matr[i, 3:6, k] = U[subset[k], :]
            matr[i, 6, k] = LA.norm(matr[i, 3:6, k])

    # Determine impact time for all clasts as the first time index
    # such that velocity (in norm) falls below toll tolerance and,
    # in the previous time step, velocity along z-axis is negative.
    # Otherwise, set impact time to be zero

    t_impact = np.zeros((nballistics, 1)).astype(int)
    time_impact = np.zeros((nballistics, 1)).astype(float)

    for ib in range(nballistics):
        for it in range(3, n_times):
            if (matr[it, -1, ib] < toll and matr[it - 1, -2, ib] < 0):
                t_impact[ib] = int(it)
                time_impact[ib] = sorted_times[it]
                break

    # Calculate mean and maximum velocities along particle trajectory
    # velocities is an array with these values:
    # column 1: ballistic index
    # column 2: diameter
    # column 3: mean velocity
    # column 4: max velocity

    velocities = np.zeros((nballistics, 4))
    for k in range(nballistics):
        velocities[k, 0] = int(k)
        velocities[k, 1] = d[k]
        velocity = matr[:int(t_impact[k]) + 1, -1, k]
        mean_velocity = np.mean(velocity)
        velocities[k, 2] = mean_velocity
        max_velocity = np.amax(velocity)
        velocities[k, 3] = max_velocity

    C = ['index', 'diameter [m]', 'mean vel [m/s]', 'max vel [m/s]']
    df = pd.DataFrame(velocities)
    df[0] = df[0].astype(int)
    df.to_csv("velocities.csv", header=C, index=False)

    # Determine the mass of each block
    r = d / 2
    V = 4 / 3 * (np.pi) * (r**3)
    m = rho * V
    K = np.zeros((nballistics, 1))

    mat1 = np.zeros((nballistics, 9))

    for s in range(nballistics):
        mat1[s, 0] = s
        if t_impact[s] != 0:

            K[s] = 0.5 * m[s] * (matr[t_impact[s] - 1, -1, s]**2)

            mat1[s, 1] = d[s]
            mat1[s, 2] = rho[s]
            mat1[s, 3] = time_impact[s]
            mat1[s, 4] = matr[t_impact[s], 0, s]
            mat1[s, 5] = matr[t_impact[s], 1, s]
            mat1[s, 6] = matr[t_impact[s], 2, s]
            mat1[s, 7] = matr[t_impact[s] - 1, -1, s]
            mat1[s, 8] = K[s]

    C = [
        'index', 'diameter [m]', 'density [kg/m3]', 'impact time [s]', 'x [m]',
        'y [m]', 'z [m]', 'impact velocity [m/s]', 'landing energy [J]'
    ]

    df = pd.DataFrame(mat1)
    df[0] = df[0].astype(int)
    df.to_csv("impacts.csv", header=C, index=False, float_format='%.3f')

    x = np.array(position[:, 0])
    y = np.array(position[:, 1])
    diam = np.array(d)

    xPmin = np.amin(x) + xc
    xPmax = np.amax(x) + xc
    yPmin = np.amin(y) + yc
    yPmax = np.amax(y) + yc
    DeltaxP = xPmax - xPmin
    DeltayP = yPmax - yPmin

    xmin = np.amin(Xinit) - 0.5 * cell
    xmax = np.amax(Xinit) + 0.5 * cell

    ymin = np.amin(Yinit) - 0.5 * cell
    ymax = np.amax(Yinit) + 0.5 * cell

    step_dens = 50.0
    x_density = np.arange(xmin, xmax, step=step_dens)
    y_density = np.arange(ymin, ymax, step=step_dens)

    x_dens_min = np.amin(x_density)
    x_dens_max = np.amax(x_density)
    y_dens_min = np.amin(y_density)
    y_dens_max = np.amax(y_density)

    extent_density = x_dens_min, x_dens_max, y_dens_min, y_dens_max

    print('extent_density', extent_density)

    xx, yy = np.meshgrid(x_density, y_density)

    count_ballistic_class = np.zeros(
        (len(diams) + 1, xx.shape[0], xx.shape[1]))
    zz = np.zeros_like(xx)

    for xi, yi, di in zip(x, y, diam):

        i = (xi + xc - x_dens_min) / step_dens
        j = (yi + yc - y_dens_min) / step_dens

        i = int(np.ceil(i))
        j = int(np.ceil(j))

        count_ballistic_class[-1, j, i] += 1
        k = int(np.argwhere(di == diams)[0])
        count_ballistic_class[k, j, i] += 1

    for i in range(len(diams) + 1):

        fig, ax = plt.subplots()

        ax.imshow(ls.hillshade(np.flipud(Zinit),
                               vert_exag=1.0,
                               dx=cell,
                               dy=cell),
                  cmap='gray',
                  extent=extent)

        ax.set_aspect('equal', 'box')

        zz[:, :] = np.squeeze(count_ballistic_class[i, :, :])
        sum_zz = np.nansum(zz)
        print('sum_zz', sum_zz)
        zz = zz / np.sum(zz) * 100.0
        zz = np.log10(zz)
        zz_max = np.amax(zz)

        zz_linspace = np.linspace(0, zz_max, num=11)
        ticks = []
        for val in zz_linspace:
            ticks.append(str(val))

        # im_ratio = (ymax - ymin) / (xmax - xmin)

        label_str = 'Probabilty [0;1]'

        cmap = plt.get_cmap('terrain_r')

        # min_arr = np.amin(zz)
        # max_arr = np.amin(zz)
        # levels = np.linspace(min_arr, max_arr, 11)
        # norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        # ax.scatter(x+xc,y+yc,c=z,s=3,alpha=0.1,edgecolors='none')

        p1 = ax.imshow(np.flipud(zz),
                       cmap=cmap,
                       interpolation='nearest',
                       extent=extent_density,
                       alpha=0.65)

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

        if 'map_size' in globals():

            ax.set_xlim(xc - map_size, xc + map_size)
            ax.set_ylim(yc - map_size, yc + map_size)

        else:

            ax.set_xlim(xPmin - 1.0 * DeltaxP, xPmax + 1.0 * DeltaxP)
            ax.set_ylim(yPmin - 1.0 * DeltayP, yPmax + 1.0 * DeltayP)

        clb = plt.colorbar(p1)

        label_str = 'Log of % ballistic'
        clb.set_label(label_str, labelpad=-40, y=1.05, rotation=0)

        if i < len(diams):

            string = '_d' + str(i) + '_'
            title = 'Diameter = ' + f'{diams[i]:.2e}' + 'm'

        else:

            string = '_tot_'
            title = 'All sizes'

        plt.title(title)

        png_file = current_dir_name + string + 'ballistic.png'

        plt.savefig(png_file, dpi=200)
        plt.show()
        plt.close(fig)

        nd = -9999

        zz[zz == -np.inf] = nd

        asc_file = current_dir_name + string + 'ballistic.asc'

        header = "ncols     %s\n" % xx.shape[1]
        header += "nrows    %s\n" % xx.shape[0]
        header += "xllcorner " + str(x_dens_min) + "\n"
        header += "yllcorner " + str(y_dens_min) + "\n"
        header += "cellsize " + str(step_dens) + "\n"
        header += "NODATA_value " + str(nd)

        np.savetxt(asc_file,
                   np.flipud(zz),
                   header=header,
                   fmt='%1.5f',
                   comments='')


if __name__ == '__main__':

    main()
