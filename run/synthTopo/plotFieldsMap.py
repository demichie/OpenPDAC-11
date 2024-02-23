from matplotlib.colors import BoundaryNorm
from matplotlib.colors import LightSource
try:
    from preprocessing.ASCtoSTLdict import domain_size_x
    from preprocessing.ASCtoSTLdict import domain_size_y
except:
    print('') 
    
from preprocessing.ASCtoSTLdict import xc
from preprocessing.ASCtoSTLdict import yc
from preprocessing.ASCtoSTLdict import DEM_file


import vtk
import os
import sys
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.util import numpy_support
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from linecache import getline
import pandas as pd
from scipy.interpolate import griddata

sys.path.insert(0, './preprocessing')

toll = 1.0
step_dens = 50.0


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
    # print('npoints',npoints)
    point = data.GetPoint(0)
    d = data.GetPointData()
    
    position_array = data.GetPoint
    alpha_array = d.GetArray("alpha.particles")  
    
    position = np.zeros((npoints,3))
    alpha_particles = np.zeros(npoints)
 
    for n in range(npoints):

        position[n,:] = data.GetPoint(n)
        alpha_particles[n] = np.array(alpha_array.GetTuple(n))
         
    return position, alpha_particles


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

    print(extent)

    ls = LightSource(azdeg=45, altdeg=45)

    working_dir = current_dir + '/VTK/terrain_out'
    files = os.listdir(working_dir)
    n_files = len(files)
    file_idx = []
    for filename in files:
        
        file_idx.append(float(filename.split('_')[-1].rsplit('.', 1)[0]))

    sorted_files = [files[i] for i in np.argsort(file_idx)]
    sorted_times = [dt * i for i in range(len(file_idx))]

    # print('Times read',sorted_times)

    full_filename = working_dir + '/' + sorted_files[-1]

    position, alpha_max = readerVTK(full_filename)

    n_times = n_files

    for i, filename in enumerate(sorted_files[:]):

        full_filename = working_dir + '/' + filename
        # print(filename)
        position, alpha_particles = readerVTK(full_filename)
        alpha_max = np.maximum(alpha_max,alpha_particles)


    print(np.amin(alpha_particles))
    print(np.amax(alpha_particles))

    step_dens = 20.0

    x_dens_min = np.amin(position[:,0])
    x_dens_max = np.amax(position[:,0])
    y_dens_min = np.amin(position[:,1])
    y_dens_max = np.amax(position[:,1])

    x_density = np.arange(x_dens_min, x_dens_max, step=step_dens)
    y_density = np.arange(y_dens_min, y_dens_max, step=step_dens)

    extent_density = x_dens_min, x_dens_max, y_dens_min, y_dens_max

    print('extent_density', extent_density)

    xx, yy = np.meshgrid(x_density, y_density)

    
    points = position[:,0:2]
    values = alpha_particles[:]
    zz = griddata(points, values, (xx, yy), method='linear')

    zz[zz<1.e-6] = np.nan

    # print(np.amin(zz))
    # print(np.amax(zz))
    
    xmin = np.amin(Xinit) - 0.5 * cell
    xmax = np.amax(Xinit) + 0.5 * cell

    ymin = np.amin(Yinit) - 0.5 * cell
    ymax = np.amax(Yinit) + 0.5 * cell

    fig, ax = plt.subplots()

    im = ax.imshow(ls.hillshade(np.flipud(Zinit),
                                vert_exag=1.0,
                                dx=cell,
                                dy=cell),
                   cmap='gray',
                   extent=extent)

    ax.set_aspect('equal', 'box')
    
    im_ratio = (ymax - ymin) / (xmax - xmin)

    levels = np.linspace(-6.0, 1.0, 11)
    label_str = 'Probabilty [0;1]'

    cmap = plt.get_cmap('viridis')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    # ax.scatter(x+xc,y+yc,c=z,s=3,alpha=0.1,edgecolors='none')

    p1 = ax.imshow(np.flipud(np.log10(zz)),
                   cmap=cmap,
                   interpolation='nearest',
                   extent=extent_density,
                   alpha=0.65)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    clb = plt.colorbar(p1)

    label_str = 'Log10 of Volume fraction'
    clb.set_label(label_str, labelpad=5, y=0.5, rotation=90)

    title = 'Solid particles'
    plt.title(title)

    png_file = current_dir_name + '_solid.png'

    plt.savefig(png_file, dpi=200)
    plt.close(fig)

    nd = -9999

    zz[zz == np.nan] = nd

    asc_file = current_dir_name + '_solid.asc'

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

    # plt.show()


if __name__ == '__main__':

    main()
