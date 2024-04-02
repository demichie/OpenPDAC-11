path = '../'
DEM_file = 'synthDEM.asc'

subsample = 2  # subsampling factor of the original DEM
dist0 = 25.0  # semi-width of the fissure/radius of the cylinder
dist_flat = 15.0
enne = 1.0  # shape parameter (1 linear, 2 spherical)
depth = 10.0  # depth of the fissure
nlevels = 4  # levels of refinements of the subsampled grid
RBF_interpolation = True # smooth the topography inside the crater
top_smooth_flag = True # smooth the top of crater volume

lagrangian_layer_depth = 10.0

xc = -140.0  # x of first point of fissure/x-center of cylinder (UTM)
yc = 30.0  # y of first point of fissure/y-center of cylinder (UTM)

xb = 0.0  # horizontal x-translation of the base of the fissure/conduit
yb = 0.0  # horizontal y-translation of the base of the fissure/conduit

# FOR CYLINDRICAL FISSURE
points = [(xc, yc), (xc, yc)]

# FOR LINEAR FISSURE
# points = [(xc, yc), (xc+20, yc+10), (xc+40, yc+15)]

conduit_radius = 5.0
conduit_length = 20.0
conduit_buffer = 3.0
conduit_shift_x = 0.0
conduit_shift_y = 0.0

xbc = 1.0
ybc = 5.0


saveDicts_flag = True
z_atm = 2000
offset_mesh = 50.0
delta_mesh = 100.0

# Relative (with respect to xc,yc) coordinates of blockMesh 
xmin = -500.0
xmax = 500.0
ymin = -500.0
ymax = 500.0

# create a surface at fixed elevation used to sample field values
z_sample = 10.0

# create probe points at elevation dz with respect to topography
xProbes = [20]
yProbes = [20]
dzProbes = [1]
