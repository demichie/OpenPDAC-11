# clear previous stl surfaces and output of previous runs
rm -rf constant/triSurface/*
foamCleanCase

# create the STL surface from the ASC file
cd preprocessing
unzip synthDEM.zip
python3 ASCtoSTL.py
cd ..

# check that the surfaces for initial conditions are closed
surfaceCheck constant/triSurface/surface_crater_closed.stl
surfaceCheck constant/triSurface/surface_conduit_closed.stl
surfaceCheck constant/triSurface/surface_total_closed.stl

# create the initial uniform grid 
cp ./system/controlDict.init ./system/controlDict
blockMesh 
checkMesh -allTopology -allGeometry

# refinement of the grid
snappyHexMesh -overwrite
checkMesh -allTopology -allGeometry

# assign common names to boundary faces 
changeDictionary

# create zones based on stl files of closed surfaces (crater, conduit,...)
topoSet -dict topoSetDict-conduit

# 
cp ./system/fvSolution.init ./system/fvSolution
cp ./constant/cloudProperties.init ./constant/cloudProperties

rm -rf 0
cp -r org.0 0

#FOR PARALLEL RUN:
#sbatch MPIJob_init.script
#squeue

#FOR SCALAR RUN:
foamRun
setFields

cp ./system/controlDict.run system/controlDict
cp ./system/fvSolution.run system/fvSolution
cp ./constant/cloudProperties.run ./constant/cloudProperties

#FOR PARALLEL RUN:
#sbatch MPIJob_run.script
#squeue

#FOR PARALLEL RUN ON PC:
decomposePar
mpirun -np xx foamRun -parallel
reconstructPar -newTimes -fields '(p U.gas alpha.particles)' -lagrangianFields '(origId U d rho)'
foamToVTK -fields '()' -noInternal -noFaceZones -excludePatches '(atm top terrain_in terrain_out)'

#FOR SCALAR RUN:
foamRun


python plotBallistics.py

rm -rf VTK
foamToVTK -fields "(alpha.particles)"

