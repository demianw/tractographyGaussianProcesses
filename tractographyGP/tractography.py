import numpy as _numpy
try:
  from scipy import interpolate as _interpolate
except:
  pass

class tractography:

  _originalFibers = []
  _originalLines = []
  _originalData = {}
  _fibers = []
  _fibersMap = []
  _subsampledFibers = []
  _quantitiyOfPointsPerfiber = None
  _interpolated = False

  def __init__(self):
    pass

  def from_vtkPolyData(self, vtkPolyData, ratio=1):
    #Fibers and Lines are the same thing
    self._originalFibers = []
    self._fiberData = {}
    self._originalLines = []
    self._originalData = {}
    self._fibers = []
    self._fibersMap = []
    self._subsampledFibers = []
    self._quantitiyOfPointsPerfiber = None
    self._interpolated = False
    self._fiberKeepRatio = ratio

    #Bug fix for the inhability to convert vtkId to a numeric type _numpy array
    linesUnsignedInt = _slicer.vtkUnsignedIntArray()
    linesUnsignedInt.DeepCopy( vtkPolyData.GetLines().GetData() )
    lines = linesUnsignedInt.ToArray().squeeze()
    points = vtkPolyData.GetPoints().GetData().ToArray()

    actualLineIndex = 0
    numberOfFibers = vtkPolyData.GetLines().GetNumberOfCells()
    for l in xrange( numberOfFibers ):
      #print '%.2f%%'%(l*1./numberOfFibers * 100.)
      self._fibers.append( points[ lines[actualLineIndex+1: actualLineIndex+lines[actualLineIndex]+1] ] )
      self._originalLines.append( _numpy.array(lines[actualLineIndex+1: actualLineIndex+lines[actualLineIndex]+1],copy=True)  )
      actualLineIndex += lines[actualLineIndex]+1

    for i in xrange( vtkPolyData.GetPointData().GetNumberOfArrays() ):
      array = vtkPolyData.GetPointData().GetArray(i)
      array_data = array.ToArray()
      self._originalData[ array.GetName() ] = array_data
      self._fiberData[ array.GetName() ] =  map( lambda f: array_data[ f ], self._originalLines )

    if (vtkPolyData.GetPointData().GetScalars()!=[]):
      self._originalData['vtkScalars']=vtkPolyData.GetPointData().GetScalars().ToArray()
    if (vtkPolyData.GetPointData().GetTensors()!=[]):
      self._originalData['vtkTensors']=vtkPolyData.GetPointData().GetTensors().ToArray()
    if (vtkPolyData.GetPointData().GetVectors()!=[]):
      self._originalData['vtkVectors']=vtkPolyData.GetPointData().GetTensors().ToArray()

    self._originalFibers = list(self._fibers)
    if self._fiberKeepRatio!=1:
      self._fibers = self._originalFibers[::_numpy.round(len(self._originalFibers)*self._fiberKeepRatio)]
    


  def from_dictionary(self, d, ratio=1 ):
    #Fibers and Lines are the same thing
    self._originalFibers = []
    self._fiberData = {}
    self._originalLines = []
    self._originalData = {}
    self._fibers = []
    self._fibersMap = []
    self._subsampledFibers = []
    self._quantitiyOfPointsPerfiber = None
    self._interpolated = False
    self._fiberKeepRatio = ratio

    #Bug fix for the inhability to convert vtkId to a numeric type _numpy array
    lines = _numpy.asarray(d['lines']).squeeze()
    points = d['points']

    actualLineIndex = 0
    numberOfFibers = d['numberOfLines']
    for l in xrange( numberOfFibers ):
      #print '%.2f%%'%(l*1./numberOfFibers * 100.)
      self._fibers.append( points[ lines[actualLineIndex+1: actualLineIndex+lines[actualLineIndex]+1] ] )
      self._originalLines.append( _numpy.array(lines[actualLineIndex+1: actualLineIndex+lines[actualLineIndex]+1],copy=True)  )
      actualLineIndex += lines[actualLineIndex]+1

    for k in d['pointData']:
      array_data = d['pointData'][k]
      self._originalData[ k ] = array_data
      self._fiberData[ k ] =  map( lambda f: array_data[ f ], self._originalLines )

#    if (vtkPolyData.GetPointData().GetScalars()!=[]):
#      self._originalData['vtkScalars']=vtkPolyData.GetPointData().GetScalars().ToArray()
#    if (vtkPolyData.GetPointData().GetTensors()!=[]):
#      self._originalData['vtkTensors']=vtkPolyData.GetPointData().GetTensors().ToArray()
#    if (vtkPolyData.GetPointData().GetVectors()!=[]):
#      self._originalData['vtkVectors']=vtkPolyData.GetPointData().GetTensors().ToArray()

    self._originalFibers = list(self._fibers)
    if self._fiberKeepRatio!=1:
      self._fibers = self._originalFibers[::_numpy.round(len(self._originalFibers)*self._fiberKeepRatio)]
    
  def to_vtkPolyData(self,vtkPolyData, selectedFibers=None):

    fibers = self.getOriginalFibers()
    if selectedFibers!=None:
      fibers = [ fibers[i] for i in selectedFibers]

    numberOfPoints = reduce( lambda x,y: x+y.shape[0], fibers,0 )
    numberOfCells = len(fibers)
    numberOfCellIndexes = numberOfPoints+numberOfCells

    linesUnsignedInt = _slicer.vtkUnsignedIntArray()
    linesUnsignedInt.SetNumberOfComponents(1)
    linesUnsignedInt.SetNumberOfTuples(numberOfCellIndexes)
    lines = linesUnsignedInt.ToArray().squeeze()
    #print lines.shape

    
    points_vtk = vtkPolyData.GetPoints()
    if points_vtk==[]:
      point_vtk = _slicer.vtkPoints()
      vtkPolyData.SetPoints(points_vtk)
#      points_vtk.Delete()
    points_vtk.SetNumberOfPoints(numberOfPoints)
    points_vtk.SetDataTypeToFloat()
    points = points_vtk.GetData().ToArray()

    actualLineIndex = 0
    actualPointIndex = 0
    for i in xrange(numberOfCells):
      lines[actualLineIndex] = fibers[i].shape[0]
      lines[actualLineIndex+1: actualLineIndex+1+lines[actualLineIndex]] = _numpy.arange( lines[actualLineIndex] )+actualPointIndex
      points[ lines[ actualLineIndex+1: actualLineIndex+1+lines[actualLineIndex] ] ] = fibers[i]
      actualPointIndex += lines[actualLineIndex]
      actualLineIndex = actualLineIndex + lines[actualLineIndex]+1


    vtkPolyData.GetLines().GetData().DeepCopy(linesUnsignedInt)

  def unsubsampleFibers(self):
    self._subsampledFibers = []
    
  def unfilterFibers(self):
    self._fibersMap = []

  def subsampleFibers(self, quantitiyOfPointsPerfiber):
    self._quantitiyOfPointsPerfiber = quantitiyOfPointsPerfiber
    self._subsampledFibers = []
    self._subsampledLines= []
    self._subsampledData = {}
    
    for k in self._fiberData:
      self._subsampledData[ k ] = []

    for i in xrange(len(self._fibers)):
      f = self._fibers[i]

      s = _numpy.linspace( 0, f.shape[0]-1,min(f.shape[0],self._quantitiyOfPointsPerfiber) ).round().astype(int)
      self._subsampledFibers.append( f[s,:] )
      self._subsampledLines.append(s)
      for k in self._fiberData:
        self._subsampledData[ k ].append( self._fiberData[k][i][s] )


    self._interpolated = False


  def subsampleInterpolatedFibers(self, quantitiyOfPointsPerfiber):
    if '_interpolate' not in dir():
      raise NotImplementedError('scipy module could not be imported, this operation is not implemented')

    self._quantitiyOfPointsPerfiber = quantitiyOfPointsPerfiber
    self._subsampledFibers = []
    s = _numpy.linspace(0,1,quantitiyOfPointsPerfiber)
    i=0.
    for f in self._fibers:
      i+=1
      tck,u = _interpolate.splprep( f.T )
      self._subsampledFibers.append( _numpy.transpose(_interpolate.spleval(s,tck)))

    self._subsampledFibers = _numpy.array( self._subsampledFibers ) 
    self._interpolated = True

  def filterFibers(self, minimumNumberOfSamples ):

    if len(self._originalFibers)==0:
      self._originalFibers = self._fibers

    self._fibers = filter( lambda f: f.shape[0]>= minimumNumberOfSamples, self._originalFibers )
    self._fibersMap = filter( lambda i: self._originalFibers[i].shape[0]>= minimumNumberOfSamples, range(len(self._originalFibers)) )

    
    if self._quantitiyOfPointsPerfiber!=None:
      if self._interpolated:
        self.subsampleInterpolatedFibers( self._quantitiyOfPointsPerfiber )
      else:
        self.subsampleFibers( self._quantitiyOfPointsPerfiber )

  def areFibersFiltered(self):
    return self._fibersMap!=[]

  def areFibersSubsampled(self):
    return self._subsampledFibers!=[]

  def areSubsampledFibersInterpolated(self):
    return self._interpolated

  def getOriginalFibersData(self):
    return self._fiberData

  def getOriginalFibers(self):
    return self._originalFibers

  def getOriginalLines(self):
    return self._originalLines

  def getOriginalData(self):
    return self._originalData

  def getFilteredFibersMap(self):
    return self._fibersMap

  def getFibersToProcess(self):
    if self._subsampledFibers!=[]:
      return self._subsampledFibers
    elif self._fibers!=[]:
      return self._fibers

  def getFibersDataToProcess(self):
    if self._subsampledData!=[]:
      return self._subsampledData
    elif self._fiberData!=[]:
      return self._fiberData


