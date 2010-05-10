from numpy import pi
try:
  from scipy import linalg
except ImportError:
  from numpy import linalg
import numpy

from . import bounded_thinPlate
from pymc import gp as _GP

def unionBoundingBox(b1,b2):
  return numpy.vstack((numpy.vstack((b1[0],b2[0])).min(0),numpy.vstack((b1[1],b2[1])).max(0)))


def intersectionBoundingBox(b1,b2):
  one_point_b1_inside_b2 = len(((b1>=b2[0])*(b1<=b2[1])).prod(1).nonzero()[0])>0
  one_point_b2_inside_b1 = len(((b2>=b1[0])*(b2<=b1[1])).prod(1).nonzero()[0])>0
  if one_point_b1_inside_b2 or one_point_b2_inside_b1:
    return numpy.vstack((numpy.vstack((b1[0],b2[0])).max(0),numpy.vstack((b1[1],b2[1])).min(0)))
  else:
    return None


fiberDelta = lambda f: numpy.sqrt(((f[1:]-f[:-1])**2).sum(1))

def _constant(x, val = 0):
    """docstring for parabolic_fun"""
    return numpy.zeros(x.shape[:-1],dtype=float) + val


cov_thinPlate3D = _GP.covariance_function_bundle( 'thinplate3d_mat', 'tractographyGP.bounded_thinPlate', {'r':'r'} )
cov_thinPlate2D = _GP.covariance_function_bundle( 'thinplate2d_mat', 'tractographyGP.bounded_thinPlate', {'r':'r'} )

class fiberGP:
  
  _k = 0
  _precomputedVariance = None
  _precomputedMean = None
  _scaleFactor = 1
  _mean = None
  _cov = None
  _alpha = None
  _inv_cov = None

  def __init__(self,fiber,  resolution=.5, rFactor=1,epsilonFactor=1,scaleFactor=1):

    if fiber.shape[1]==2:
      self._c = cov_thinPlate2D.euclidean
#      self._c = 
    elif fiber.shape[1]==3:
      self._c = cov_thinPlate2D.euclidean
    else:
      raise 

    self._scaleFactor = scaleFactor
    self._bounds = numpy.vstack( ( fiber.min(0), fiber.max(0) ) )
    self._adjustedBounds = self._bounds.copy()
    self._resolution = resolution
    self._fiber = fiber
    self._trainingPoints = numpy.ones( len(self._fiber) )

    delta = fiberDelta( fiber ).max()
    self._epsilon = epsilonFactor*delta
    self.setK( rFactor*delta )

  def getResolution(self):

    return self._resolution

  def setResolution(self,resolution):
    self._resolution = resolution
    self.precomputeVarianceField()
    

  def getK(self):
    return self._k

  def setK( self, k ):
    if k!=self._k:
      deltaK = (k-self._k)
      self._k = k
      self._adjustedBounds[0]-=deltaK
      self._adjustedBounds[1]+=deltaK
      
      self._cov = _GP.Covariance( self._c, amp=1, r= self._k  ) 
      self._cov_matrix = self._cov( self._fiber, self._fiber )
      self._cov_inv_matrix = linalg.inv( self._cov_matrix )
      self._alpha = numpy.asmatrix(numpy.dot( numpy.ones(len(self._fiber)), self._cov_inv_matrix )).T

      self._mean = _GP.Mean( _constant, val = 0 )
      self._gp = _GP.observe( self._mean, self._cov, obs_mesh = self._fiber, obs_vals=numpy.ones( len(self._fiber)), obs_V = numpy.zeros(len(self._fiber)  )+self._epsilon )

      if self._precomputedVariance!=None:
        self.precomputeVarianceField()
      if self._precomputedMean!=None:
        self.precomputeMeanField()

  def getBoundingBox( self ):
    return self._bounds

  def getAdjustedBoundingBox( self ):
    return self._adjustedBounds

  def getGaussianProcess( self ):
    return self._gp

  def intersectionBoundingBox(self, boundingBox ):
    b1 = self._adjustedBounds
    b2 = boundingBox
    return intersectionBoundingBox(b1,b2)


  def precomputeVarianceField( self ):

    bb = self.getAdjustedBoundingBox()

    if self._fiber.shape[1] == 2:

      self._precomputedVarianceGrid = numpy.mgrid[\
          bb[0,0] :bb[1,0] :self._resolution,
          bb[0,1] :bb[1,1] :self._resolution\
              ]
      self._precomputedVariance = self._cov( self._precomputedVarianceGrid.T )

    elif self._fiber.shape[1] == 3:
   
      self._precomputedVarianceGrid = numpy.mgrid[\
          bb[0,0] :bb[1,0] :self._resolution,
          bb[0,1] :bb[1,1] :self._resolution,\
          bb[0,2] :bb[1,2] :self._resolution\
              ]
      self._precomputedVariance = self._cov( self._precomputedVarianceGrid.T )

    else:
      raise

  def getPrecomputedVarianceField(self):
    if self._precomputedVariance==None:
      self.precomputeVarianceField()
    return self.getAdjustedBoundingBox(),self._precomputedVariance

  def precomputeMeanField( self ):

    bb = self.getAdjustedBoundingBox()

    if self._fiber.shape[1] == 2:

      self._precomputedMeanGrid = numpy.mgrid[\
          bb[0,0] :bb[1,0] :self._resolution,
          bb[0,1] :bb[1,1] :self._resolution\
              ]
      self._precomputedMean = self._mean( self._precomputedMeanGrid.T )

    elif self._fiber.shape[1] == 3:
    
      self._precomputedMeanGrid = numpy.mgrid[\
          bb[0,0] :bb[1,0] :self._resolution,
          bb[0,1] :bb[1,1] :self._resolution,\
          bb[0,2] :bb[1,2] :self._resolution\
              ]
      self._precomputedMean = self._mean( self._precomputedMeanGrid.T )

    else:
      raise

  def getPrecomputedMeanField(self):
    if self._precomputedMean==None:
      self.precomputeMeanField()
    return self.getAdjustedBoundingBox(),self._precomputedMean

  def getVarianceFieldForBoundingBox( self, boundingBox, resolution=None):
    if resolution==None:
      resolution = self._resolution


    if self._fiber.shape[1] == 2:

      precomputedVarianceGrid = numpy.mgrid[\
          boundingBox[0,0] :boundingBox[1,0] :resolution,
          boundingBox[0,1] :boundingBox[1,1] :resolution\
              ]
      return self._cov( self._precomputedVarianceGrid.T )

    elif self._fiber.shape[1] == 3:
    
      precomputedVarianceGrid = numpy.mgrid[\
          boundingBox[0,0] :boundingBox[1,0] :resolution,
          boundingBox[0,1] :boundingBox[1,1] :resolution,\
          boundingBox[0,2] :boundingBox[1,2] :resolution\
              ]
      return self._cov( precomputedVarianceGrid.T )

      assaise

  def getMeanFieldForBoundingBox( self, boundingBox, resolution=None):
    if resolution==None:
      resolution = self._resolution

    if self._fiber.shape[1] == 2:

      precomputedMeanGrid = numpy.mgrid[\
          boundingBox[0,0]:boundingBox[1,0]:resolution,
          boundingBox[0,1]:boundingBox[1,1]:resolution\
              ]
      return self._cov( self._precomputedMeanGrid.T )

    elif self._fiber.shape[1] == 3:
    
      precomputedMeanGrid = numpy.mgrid[\
          boundingBox[0,0] :boundingBox[1,0]:resolution,
          boundingBox[0,1] :boundingBox[1,1]:resolution,\
          boundingBox[0,2] :boundingBox[1,2]:resolution\
              ]
      return  self._mean( precomputedMeanGrid.T )

    else:
      raise

  def getMean(self, point ):
    if self._fiber.shape[1]==3:
      c3 = self._c3(self._k)
      kstar = numpy.apply_along_axis( c3, 1, self._fiber, point )
      gradient = numpy.dot( kstar.T, self._gp._alpha )
      
      return gradient
    else:
      raise NotImplemented


  def __getstate__(self):
        from types import MethodType, FunctionType
  	return dict([ i for i in self.__dict__.items() if (type(i[1])!=MethodType) and (type(i[1])!=FunctionType) ])

#  def __setstate__(self,d):
#    self.__dict__ = d
#
#    if self._fiber.shape[1]==2:
#      self._c = self._c2
#    elif self._fiber.shape[1]==3:
#      self._c = self._c3
#    else:
#      raise 
#
#    self._gp.setCovarianceFunction( self._c(self._k))




   


nD = lambda D,v:numpy.sqrt( numpy.dot(v,numpy.dot(D,v)))
km = lambda D,x0: lambda x,y: exp(- nD(inv(D)  ,x-y) ) if allclose(x,x0) or allclose(y,x0) else ( 1e-10 if allclose(x,y) else 0 )
km2 = lambda D,x0,x1: lambda x,y: exp(- nD(inv(D) if allclose(x,x0) else inv(D) if allclose(x,x1) else eye(len(D)),(x-y)))

sDm = lambda Ds,Fs,x:(inv(Ds[ Fs==x ]) if len(Ds[Fs==x])>1 else eye(len(Ds[0])))
kmIx = lambda Ds,Fs: lambda x,y: exp(- nD( sDm(Ds,Fs,x)  ,(x-y)))
kmIx2 = lambda Ds,Fs: lambda x,y: (kmIx(Ds,Fs)(x,y)+kmIx(Ds,Fs)(y,x))/2.

def submatrixGrids( b1, b2, resolution=2. ):
  b1 = numpy.floor(b1)
  b2 = numpy.floor(b2)
  isec = intersectionBoundingBox(b1,b2)
  if isec == None:
    return None,None
  g1 = numpy.floor(numpy.vstack((   isec[0]-b1[0],  isec[1]-b1[0] ))/resolution)
  g2 = numpy.floor(numpy.vstack((   isec[0]-b2[0],  isec[1]-b2[0] ))/resolution)
  mindist =numpy.vstack( ( g1[1]-g1[0],g2[1]-g2[0] ) ).min(0)
  g1 = numpy.vstack(  (g1[0],g1[0]+mindist))
  g2 = numpy.vstack( (g2[0],g2[0]+mindist))
  g1_b = map( lambda s: slice(s[0],s[1]),g1.T)
  g2_b = map( lambda s: slice(s[0],s[1]),g2.T)
  return g1_b, g2_b

def unionCMap(p1, p2, resolution = 2):
  unionBB = unionBoundingBox( p1[0], p2[0] )
  dimensions = numpy.ceil((unionBB[1]-unionBB[0])*1./resolution)
  unionBox = numpy.zeros( dimensions )
  unionBox[:]+= p1[1].max()+p2[1].max()

  base1 = numpy.floor((p1[0][0]-unionBB[0])*1./resolution)
  end1 = base1+p1[1].shape

  unionBox[ base1[0]:end1[0], base1[1]:end1[1], base1[2]:end1[2] ] += p1[1]-p1[1].max()

  base2 = numpy.floor((p2[0][0]-unionBB[0])*1./resolution)
  end2 = base2+p2[1].shape

  unionBox[ base2[0]:end2[0], base2[1]:end2[1], base2[2]:end2[2] ] += p2[1]-p2[1].max()

  return unionBB,unionBox/4.

def unionCMapNN( lcmaps, resolution = 2):

  unionBB = reduce( lambda x,y:  unionBoundingBox( x,  y[0] ), lcmaps[1:], lcmaps[0][0] )
  dimensions = numpy.ceil((unionBB[1]-unionBB[0])*1./resolution)

  unionBox = numpy.zeros( dimensions )
  maxima = numpy.array(map( lambda p: p[1].max(), lcmaps ))
  unionBox[:]+= maxima.sum()

  for i in xrange(len(lcmaps)):
    p = lcmaps[i]
    base = numpy.floor((p[0][0]-unionBB[0])*1./resolution)
    end = base+p[1].shape
    unionBox[ base[0]:end[0], base[1]:end[1], base[2]:end[2] ] += p[1]
    unionBox[ base[0]:end[0], base[1]:end[1], base[2]:end[2] ] -= maxima[i]


  return unionBB,unionBox



def PMap(p1, p2, resolution = 2):
  unionBB = unionBoundingBox( p1[0], p2[0] )
  dimensions = numpy.ceil((unionBB[1]-unionBB[0])*1./resolution)
  Box = numpy.zeros( dimensions )

  base1 = numpy.floor((p1[0][0]-unionBB[0])*1./resolution)
  end1 = base1+p1[1].shape

  Box[ base1[0]:end1[0], base1[1]:end1[1], base1[2]:end1[2] ] += p1[1]

  base2 = numpy.floor((p2[0][0]-unionBB[0])*1./resolution)
  end2 = base2+p2[1].shape

  unionBox[ base2[0]:end2[0], base2[1]:end2[1], base2[2]:end2[2] ] += p2[1]

  return unionBB,unionBox



def jointPMap( p1, p2, resolution=2 ):
  g1,g2 = submatrixGrids( p1[0],p2[0], resolution=resolution )
  if g1==None:
    return numpy.array( [0] )
  return p1[1][g1]*p2[1][g2]

def unionPMapNN( lcmaps, resolution = 2):

  unionBB = reduce( lambda x,y:  unionBoundingBox( x,  y[0] ), lcmaps[1:], lcmaps[0][0] )
  dimensions = numpy.ceil((unionBB[1]-unionBB[0])*1./resolution)

  unionBox = numpy.zeros( dimensions )

  for i in xrange(len(lcmaps)):
    p = lcmaps[i]
    base = numpy.floor((p[0][0]-unionBB[0])*1./resolution)
    end = base+p[1].shape
    unionBox[ base[0]:end[0], base[1]:end[1], base[2]:end[2] ] += p[1]

  return unionBB,unionBox

def unionPMap( lcmaps, resolution = 2):

  unionBB = reduce( lambda x,y:  unionBoundingBox( x,  y[0] ), lcmaps[1:], lcmaps[0][0] )
  dimensions = numpy.ceil((unionBB[1]-unionBB[0])*1./resolution)

  unionBox = numpy.zeros( dimensions )

  for i in xrange(len(lcmaps)):
    p = lcmaps[i]
    base = numpy.floor((p[0][0]-unionBB[0])*1./resolution)
    end = base+p[1].shape
    unionBox[ base[0]:end[0], base[1]:end[1], base[2]:end[2] ] += p[1]

  unionBox/=len(lcmaps)

  return unionBB,unionBox


def innerProductN( p1, p2, resolution=2):
  jm = jointPMap( p1, p2 , resolution=resolution)
  return (2 * jm.sum())/((p1[1]**2).sum()+(p2[1]**2).sum())

def innerProductG( p1, p2 , resolution=2):
  jm = jointPMap( p1, p2 ,resolution=resolution)
  return (jm.sum())

def innerProductR( p1, p2, resolution=2):
  jm = jointPMap( p1, p2, resolution=resolution )
  return (jm.sum())/numpy.sqrt(((p1[1]**2).sum())*((p2[1]**2).sum()))

fibersPMap = lambda h,m: map( lambda f: (f[0],1./(2*numpy.sqrt(pi*(h**2+f[1]**2)))), m)

areaList = lambda m: map( lambda bb: area(bb[0]), m )
area = lambda bb: (bb[1]-bb[0]).prod()

probsDN =  lambda min,max:lambda h, v : (v[0], numpy.sqrt(  ((h**2+max)/(h**2+v[1]) )-1.) /numpy.sqrt(  ((h**2+max)/(h**2+min) )-1.) )

normProb = lambda m: (m-m.min())/(m.max()-m.min())
probsNN =  lambda h, v : (v[0],1./(2*numpy.sqrt(numpy.pi*(h**2+v[1]))))

probImageNN = lambda h,im : (1./(2*numpy.sqrt(numpy.pi*(h**2+im))))
probImage = lambda h,im : normProb(1./(2*numpy.sqrt(numpy.pi*(h**2+im))))


probs =  lambda h, v : (v[0], normProbImage(h,v[1]) )
probsMap = lambda h,m: map( lambda v: (v[0],normProb(1./(2*numpy.sqrt(numpy.pi*(h**2+v[1]))))),m)

def getInnerProducts( fibersProbMap, innerPs = None, fileName=None ):

  if innerPs == None:
    innerPs = numpy.zeros( (len(fibersProbMap),len(fibersProbMap)) )

  areas = areaList( fibersProbMap )
  sortedAreas = numpy.argsort( areas )

  for i0 in xrange(len(sortedAreas)):
    i = sortedAreas[i0]
    for j0 in xrange(i0,len(sortedAreas)):
      j = sortedAreas[j0]
      innerPs[i,j] = innerProductG( fibersProbMap[i], fibersProbMap[j] )
      innerPs[j,i] = innerPs[i,j]

    print '%.2f'%(i0*1./len(sortedAreas))
    if i0%numpy.round( len(sortedAreas)*.1) and fileName!=None:
      numpy.save(fileName,innerPs)

    numpy.save(fileName,innerPs)
  return innerPs

riemannNorm = lambda m: (lambda d: d*m*d.T)( (numpy.diag(m)**(-.5)).reshape((len(m),1)))


sumNorm = lambda m: (lambda d: m*2./(d+d.T))( diag(m).reshape((len(m),1)))

squareDistanceMatrix = lambda m: (lambda d: d+(-2*m)+d.T)( numpy.diag(m).reshape((len(m),1)))


def innerProduct_thinPlate_R3( fiberGP1, fiberGP2 ):
  from numpy import pi,zeros,dot,vectorize,sqrt,newaxis
  from . import bounded_thinPlate as bt

  assert( fiberGP1._fiber.shape[-1]==3 and  fiberGP2._fiber.shape[-1]==3 )

  fiber1 = fiberGP1._fiber
  fiber2 = fiberGP2._fiber

  R = fiberGP1._k
  Q = fiberGP2._k

  distmatrix = sqrt(((fiber2[newaxis,...,:]-fiber1[...,newaxis,:])**2).sum(-1))

  if len(fiber1)==len(fiber2):
    bt.innerproduct_thinplate3d_normalized( distmatrix, R, Q, symm=True)
  else:
    bt.innerproduct_thinplate3d_normalized( distmatrix, R, Q, symm=False)

  return float(dot(dot( fiberGP1._alpha.T,distmatrix),fiberGP2._alpha))


