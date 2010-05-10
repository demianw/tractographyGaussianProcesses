def getMostRepresentativeFiber( fiberBundle, subsample=20, useOnlyTopPercentile=0 ):
  from numpy import eye, array, asarray
  from myPackage import SlicerFibersProcessing as sfp
  from . import fiberGPAnalysis as fga

  if hasattr( fiberBundle, '__getitem__' ):
    fibersToProcess = fiberBundle
  elif hasattr( fiberBundle, 'GetPolyData' ):
    tr = sfp.tractography()
    tr.from_vtkPolyData( fiberBundle.GetPolyData() )
    tr.subsampleFibers( subsample )
    fibersToProcess = tr.getFibersToProcess()
  else:
    raise ValueError('Fiber bundle type is not recognized')


  if useOnlyTopPercentile>0:
    if ( useOnlyTopPercentile>1.):
      raise ValueError('Percentile should be greater or equal than 0 and smaller than 1')

    fiberLengths = array([ fiberLength(f) for f in fibersToProcess ])
    fibersToProcess = asarray( fibersToProcess )[ fiberLengths > (fiberLengths.max()*(1.-useOnlyTopPercentile))  ]
  
  print "Using ",len(fibersToProcess)," fibers"

  gps = [ fga.fiberGP( f ) for  f in fibersToProcess ]

  ips = eye( len(gps), dtype='float32' )

  for i in xrange(len(gps)):
    for j in xrange(0,i+1):
      ips[i,j] = fga.innerProduct_thinPlate_R3( gps[i], gps[j] )
      ips[j,i] = ips[i,j]
 
  mrFiberId = ips.sum(1).argmax()

  if hasattr( fiberBundle, '__getitem__' ):
    return fibersToProcess[ mrFiberId ]
  else:
    return fibersToProcess[ mrFiberId ]


def slicerMRMLVolumeNode2Image_ras2ijk( volumeNode ):
  from Slicer import slicer
  from numpy import eye, ndindex

  vtk_ras2ijk = slicer.vtkMatrix4x4()
  volumeNode.GetRASToIJKMatrix( vtk_ras2ijk )

  image = volumeNode.GetImageData().ToArray().T

  ras2ijk = eye(4)
  for c in ndindex( ras2ijk.shape ):
    ras2ijk[c] = vtk_ras2ijk.GetElement( c[0], c[1] )

  return image, ras2ijk


def straightenImage( fiber, ijkImage, ras2ijk, radiusL1=50, resolution=1, exactTimes = False, timeSamples=40 ):
  from myPackage import splineTools as spt
  from numpy import linspace, dot, ones,hstack,transpose
  from scipy import interpolate as interp


  fiber_ijk = dot( ras2ijk, hstack((fiber, ones((len(fiber),1)))).T).T[:,:3]

  spline = interp.splprep( fiber_ijk.T,task=-1,t=linspace(0,1,len(fiber_ijk)) )[0]

  if exactTimes:
    mrFiber_length = spt.arcLength( spline,0,1)
    times = spt.findTimes( linspace(0, mrFiber_length, timeSamples ),spline )
  else:
    times = linspace( 0, 1, timeSamples )

  points = transpose(interp.splev( times, spline ))

  

  tangents = spt.tangents( spline, times )
  normalPlanes = spt.normalPlanes( tangents )

  return straightenImageFromPointsAndPlanes( points, normalPlanes, ijkImage, radiusL1, resolution )


def straightenImageFromPointsAndPlanes( points, normalPlanes, image,radiusL1=11,resolution=1 ):
  from scipy import ndimage
  from numpy import zeros, linspace, isscalar, asarray, diag, sqrt
  from myPackage import splineTools as spt


  if isscalar(radiusL1):
    radiusL1 = (radiusL1,radiusL1)
  if isscalar(resolution):
    resolution = (resolution, resolution)

  voxels = asarray(radiusL1)*resolution


  image_around_fiber = zeros((len(points),voxels[0],voxels[1]))

  for i in xrange(len(points)):
          planeCoords_points = spt.planeCoords(
              normalPlanes[i],  
              linspace(-radiusL1[0],radiusL1[0],voxels[0]),
              linspace(-radiusL1[1],radiusL1[1],voxels[1]) )+points[i]

          image_plane = ndimage.map_coordinates( image, planeCoords_points.T )
          image_around_fiber[i][:] = image_plane.reshape( voxels[0], voxels[1] )

  ijkss2ijk = diag( [1./resolution[0], 1./resolution[0], 1./fiberMeanStep(points) ,1 ] )
  planeCoords_points_0 = spt.planeCoords(
              normalPlanes[0],  
              linspace(-radiusL1[0],radiusL1[0],voxels[0]),
              linspace(-radiusL1[1],radiusL1[1],voxels[1]) )+points[0]

  ijkss2ijk[:3,-1] = planeCoords_points_0[0][::-1]
  

  return image_around_fiber.T, ijkss2ijk

import numpy as _np
fiberLength = lambda points: _np.sqrt( ( (points[1:]-points[:-1])**2).sum(1) ).sum()
fiberMeanStep = lambda points: _np.sqrt( ( (points[1:]-points[:-1])**2).sum(1) ).mean() 
