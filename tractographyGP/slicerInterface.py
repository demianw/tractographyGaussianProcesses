import numpy as _numpy
try:
  import scipy as _scipy
  from scipy import optimize as _optimize
except ImportError:
  pass

try:
  import Slicer as _Slicer
  from Slicer import slicer as _slicer
except ImportError:
  pass

from tractography import *

class VisualizeDendrogram:
# decorators
  def refreshDisplay(meth):
      def new(self, *args, **kws):
        self.preRefreshDisplayProperties()
        res = meth(self, *args, **kws)
        self.postRefreshDisplayProperties()
        return res
      return new

  def refreshDisplayGenerator(meth):
      def new(self, *args, **kws):
        gen = meth(self, *args, **kws)
        while True:
          self.preRefreshDisplayProperties()
          n = gen.next()
          self.postRefreshDisplayProperties()
          yield n
      return new

  def __init__(self, dendrogram, slicer, fiberBundleNode, modelPrefixName='fiber_', tubePrefixName='tube_', lineNumberArrayName='lineNumber',colorTable=None):
    self._dendrogram = dendrogram
    self._slicer = slicer
    self._modelPrefixName = modelPrefixName
    self._tubePrefixName = tubePrefixName
    self._lineNumberArrayName = lineNumberArrayName
    self._fiberBundleNode = fiberBundleNode
    self._displayedModels = {}
    self._displayedModelsTube = {}
    self._tubeDisplay = False
    self._colors = (1.,0.,0.)
    self._vn = None
    self._colorByTable = (colorTable!=None)
    self._colorTable = colorTable

  def setColorByTable(self, value, list=''):
    if value:
      if list == '':
        self._colorTable = list
      else:
        if self._colorTable==None:
          raise ValueError('the class must have a colorTable')
      self._colorByTable = True
    else:
      self._colorByTable = False

  def setColorTable(self, list ):
    self._colorTable = list

  def setTubeDisplay( self, value ):
    self._tubeDisplay = value

  def tubeDisplay( self ):
    modelsToAdd = [ (k,i) for (k,i) in self._displayedModels.items() if self._tubePrefixName+k not in self._displayedModels.keys() ]

    for k,i in modelsToAdd:
      newTube = lineModelToTubeModel( self._slicer, i, self._tubePrefixName )
      self._displayedModelsTube[newTube[0]]= newTube[1] 

  @refreshDisplay
  def displayThresholdedByLevel(self, lmin, lmax ):
    self._displayedModels = modelsForSelectedTractographyFibers(\
        self._slicer,\
        dg_threshold_level_bc( self._dendrogram, lmin, lmax ),\
        self._dendrogram,  self._fiberBundleNode,\
        modelPrefix=self._modelPrefixName, lineNumberArrayName=self._lineNumberArrayName, models=self._displayedModels )

  @refreshDisplay
  def displayThresholdedByEnergy(self, lmin, lmax ):
    self._displayedModels = modelsForSelectedTractographyFibers(\
        self._slicer,\
        dg_threshold_energy_bc( self._dendrogram, lmin, lmax ),\
        self._dendrogram,  self._fiberBundleNode,\
        modelPrefix=self._modelPrefixName, lineNumberArrayName=self._lineNumberArrayName, models=self._displayedModels )

  @refreshDisplay
  def displayDendrogramElement(self, index ):
    self._displayedModels = modelsForSelectedTractographyFibers(\
        self._slicer,\
        [index],\
        self._dendrogram,  self._fiberBundleNode,\
        modelPrefix=self._modelPrefixName, lineNumberArrayName=self._lineNumberArrayName, models=self._displayedModels )

  
  @refreshDisplayGenerator
  def displayThresholdedByLevelGenerator(self,lmin,lmax):

    for cluster in dg_threshold_level_bc( self._dendrogram, lmin, lmax ):
      self._displayedModels = modelsForSelectedTractographyFibers(\
        self._slicer,\
        [cluster],\
        self._dendrogram,  self._fiberBundleNode,\
        modelPrefix=self._modelPrefixName, lineNumberArrayName=self._lineNumberArrayName, models=self._displayedModels )

      yield self._displayedModels
  
  @refreshDisplayGenerator
  def displayThresholdedByEnergyGenerator(self,emin,emax):

    for cluster in dg_threshold_energy_bc( self._dendrogram, emin, emax ):
      self._displayedModels = modelsForSelectedTractographyFibers(\
        self._slicer,\
        [cluster],\
        self._dendrogram,  self._fiberBundleNode,\
        modelPrefix=self._modelPrefixName, lineNumberArrayName=self._lineNumberArrayName, models=self._displayedModels )

      yield self._displayedModels


  @refreshDisplayGenerator
  def displayElementHistoryGenerator(self,element, backward=True, forward = True, reverse=False):

    elementSet = set(self._dendrogram[ element ][1] )

    forwardHistory = set([])
    if forward:
      forwardHistory =set([ i for i in xrange( len(self._dendrogram) ) if i>=element and len( elementSet.intersection(self._dendrogram[i][1]) )>0 ] )

    backwardHistory = set([])
    if backward:
      backwardHistory = set([ i for i in xrange( len(self._dendrogram) ) if i<=element and len( elementSet.intersection(self._dendrogram[i][1]) )>0 ] )

    elementHistory = backwardHistory.union( forwardHistory )

    for l in sorted(set([self._dendrogram[i][0] for i in elementHistory]),reverse=reverse): 
      clusters = [ i for i in elementHistory if self._dendrogram[i][0]==l ]
      l2=l-1
      while l2>=0:
        shownFibers = reduce( lambda x,y:x.union(self._dendrogram[y][1]),clusters,set()) 
        clusters+=[ i for i in elementHistory if self._dendrogram[i][0]==l2 and not shownFibers.issuperset( self._dendrogram[i][1] ) ]
        l2 = l2-1
      self._displayedModels = modelsForSelectedTractographyFibers(\
        self._slicer,\
        clusters,\
        self._dendrogram,  self._fiberBundleNode,\
        modelPrefix=self._modelPrefixName, lineNumberArrayName=self._lineNumberArrayName, models=self._displayedModels )

      yield self._displayedModels

  @refreshDisplay
  def displayClusters(self,clusters):

    self._displayedModels = modelsForSelectedTractographyFibers(\
      self._slicer,\
      clusters,\
      self._dendrogram,  self._fiberBundleNode,\
      modelPrefix=self._modelPrefixName, lineNumberArrayName=self._lineNumberArrayName, models=self._displayedModels )

    return self._displayedModels

  @refreshDisplayGenerator
  def displayClustersGenerator(self,clusters):

    for cluster in clusters:
      self._displayedModels = modelsForSelectedTractographyFibers(\
        self._slicer,\
        [cluster],\
        self._dendrogram,  self._fiberBundleNode,\
        modelPrefix=self._modelPrefixName, lineNumberArrayName=self._lineNumberArrayName, models=self._displayedModels )

      yield self._displayedModels

  def randomColorModels(self):
    from numpy import random
    random.seed(0)
    for m in [ self._displayedModels[k] for k in sorted( self._displayedModels.keys() ) ]:
      if m.GetDisplayNode()!=[]:
        m.GetDisplayNode().SetColor( 1-random.rand(), 1-random.rand(), 1-random.rand() ) 


  def colorModelsByTable(self):
    i=0
    for m in [ self._displayedModels[k] for k in sorted( self._displayedModels.keys() ) ]:
      if m.GetDisplayNode()!=[]:
        m.GetDisplayNode().SetColor( self._colorTable[i,0], self._colorTable[i,1], self._colorTable[i,2]) 
        i=(i+1)%self._colorTable.shape[0]


  def colorModels(self,colors):
    for m in self._displayedModels.values():
      if m.GetDisplayNode()!=[]:
        m.GetDisplayNode().SetColor( colors[0], colors[1], colors[2] )

  def setDefaultColor(self, colors ):
    self._colors = colors

  def removeLineModels(self):
    for m in self._displayedModels.values():
      self._slicer.MRMLScene.RemoveNode( m )

    self._displayedModels = {}

  def removeTubeModels(self):
    for m in self._displayedModelsTube.values():
      self._slicer.MRMLScene.RemoveNode( m )

    self._displayedModelsTube = {}


  def removeModels(self):
    self.removeLineModels()
    self.removeTubeModels()

  def preRefreshDisplayProperties(self):
    if self._tubeDisplay:
      self.removeTubeModels()

  def postRefreshDisplayProperties(self):
    if self._colorByTable:
      self.colorModelsByTable()
    elif len( self._displayedModels )>1:
      self.randomColorModels()
    else:
      self.colorModels( self._colors)

    if self._tubeDisplay:
      self.tubeDisplay()

  def __del__(self):
    self.removeModels()


class DendrogramVolumeLabelAnalysis:

  def __init__(self, dendrogram, slicer, fiberBundleNode, scalarVolumeLabelNode ):
    from numpy import asmatrix, asarray, zeros, ndindex,swapaxes
    self._dendrogram = dendrogram
    self._slicer = slicer
    self._fiberBundleNode = fiberBundleNode
    self._scalarVolumeLabelNode = scalarVolumeLabelNode

    self._tractography = tractography()
    self._tractography.from_vtkPolyData( fiberBundleNode.GetPolyData() )


    mat = self._slicer.vtkMatrix4x4()
    self._matRASToIJK = asmatrix(zeros((4,4)))
    self._scalarVolumeLabelNode.GetRASToIJKMatrix( mat )
    for c in ndindex( self._matRASToIJK.shape ):
       self._matRASToIJK[c] = mat.GetElement(c[0],c[1])
    self._labelsImage = asarray(swapaxes( self._scalarVolumeLabelNode.GetImageData().ToArray(), 0,2 ))


  def getDendrogramNodeLabelFrequencies( self, cluster ):
    from numpy import hstack, vstack, round, dot, asarray, bincount, ones

    bundle = vstack([ self._tractography.getOriginalFibers()[i] for i in self._dendrogram[cluster][1] ] )
    fibs_all_homogeneous = hstack( (bundle, ones((len(bundle),1))) )
    fibs_all_IJK = round(dot( self._matRASToIJK, fibs_all_homogeneous.T ).T[:,:-1]).astype(int)

    fibs_all_IJK_labels = self._labelsImage[ tuple(fibs_all_IJK.T) ].squeeze()

    counts = bincount( fibs_all_IJK_labels )
    freqs = counts*1./counts.sum()

    return freqs

  def getFiberCounts( self ):
    from numpy import hstack, vstack,floor, round, dot, asarray, bincount, ones,zeros,histogram
    from scipy import sparse

    originalFibers = self._tractography.getOriginalFibers()

    labelCrossings = sparse.lil_matrix( (len(originalFibers), self._labelsImage.max()+1) )

    dimensions = self._labelsImage.shape

    for i in xrange( len(originalFibers) ):
      fiber = originalFibers[i]
      fibs_homogeneous = hstack( (fiber, ones((len(fiber),1))) )
      fibs_IJK = round(dot( self._matRASToIJK, fibs_homogeneous.T ).T[:,:-1]).astype(int)
      #clip points outside
      fibs_IJK[:,0][ fibs_IJK[:,0]>= dimensions[0] ] = dimensions[0]-1
      fibs_IJK[:,1][ fibs_IJK[:,1]>= dimensions[1] ] = dimensions[1]-1   
      fibs_IJK[:,2][ fibs_IJK[:,2]>= dimensions[2] ] = dimensions[2]-1

      fibs_IJK_labels = self._labelsImage[ tuple(fibs_IJK.T) ].squeeze()

      counts = bincount( fibs_IJK_labels )

      if len(counts)>1:
        labelCrossings[i,:len(counts)] = counts
      else:
         labelCrossings[i,0] = counts[0]

    return labelCrossings.tocsr()

  def getLabelCrossingsDictionaryHieThresholdedByLevel(self, min, max ):
    return dict( [(t, self.getDendrogramNodeLabelFrequencies(t)) for t in dg_threshold_level_bc(self._dendrogram,min,max) ] )

  def getTractsCrossingSectionsQuerier(self,min=10):
    return TractsCrossingSectionsQuerier( self, self.getFiberCounts(), min=min )


class TractsCrossingSectionsQuerier:
    def __init__(self, dendrogramVolumeLabelAnalysis, labelCrossings, min=10):
      self._dendrogramVolumeLabelAnalysis = dendrogramVolumeLabelAnalysis
      self._dendrogram = self._dendrogramVolumeLabelAnalysis._dendrogram
      self._labelCrossings = labelCrossings
      self._min = min

    def getSingleLevelDendrogramCrossings( self, level ):
        from numpy import asarray

        clusters =[ (t, self._dendrogram[t][1]) for t in xrange(len(self._dendrogram)) if self._dendrogram[t][0] == level  ]
        return dict( [(t[0], (lambda x:x*1./x.sum())(asarray(self._labelCrossings[t[1],:].sum(0)).squeeze() )) for t in clusters ] )

    def getTractOptimizingSectionCrossing(self,sections,threshold=.01, start=None):
      from numpy import asarray,argmax,max
      if start==None:
        start = self._min

      maxHeight = max([ v[0] for v in self._dendrogram ])

      lastValue = 0
      selBundle = None
      while start<=maxHeight:
        singleLevelDendrogramCrossings = self.getSingleLevelDendrogramCrossings( start ) 
        bundlesCrossing = asarray([ (i[0], i[1][sections].prod() ) for i in singleLevelDendrogramCrossings.items() if all( i[1][sections]>threshold )  ])
        if len(bundlesCrossing)>0:
          maxPos = argmax( bundlesCrossing[:,1] )
          maxItem = bundlesCrossing[maxPos,:]
          if maxItem[1]>lastValue:
            lastValue = maxItem[1]
            selBundle = int(maxItem[0])
        start=start+1

      return selBundle, lastValue

    def getTractsCrossingSections(self, sections, threshold=0 ):
      from numpy import max,all

      return  [ (i[0], i[1][sections]) for i in self._labelCrossingsDictionary.items() if len(i[1])>=max(sections)+1 and all( i[1][sections]>threshold )  ]



dg_threshold_level = lambda dg,i,j: filter( lambda d:  ( dg[d][0]>=i and dg[d][0]<=j ) , xrange(len(dg)))
dg_threshold_energy = lambda dg,i,j: filter( lambda d: len(dg[d])>=3 and ( dg[d][3]>i and dg[d][3]<j ) , xrange(len(dg)))
dg_biggest_clusters = lambda dg,indexes: filter( lambda y: all( [ len(set(dg[y][1]).intersection( dg[i][1]))==0 for i in indexes if dg[i][-1]<dg[y][-1] ] ),  indexes  ) 
dg_smallest_clusters = lambda dg,indexes: filter( lambda y: all( [ len(set(dg[y][1]).intersection( dg[i][1]))==0 for i in indexes if dg[i][-1]>dg[y][-1] ] ),  indexes  ) 


dg_threshold_energy_bc = lambda dg,i,j: dg_biggest_clusters( dg, dg_threshold_energy( dg,i,j))
dg_threshold_energy_sc = lambda dg,i,j: dg_smallest_clusters( dg, dg_threshold_energy( dg,i,j))

dg_threshold_level_bc = lambda dg,i,j: dg_biggest_clusters( dg, dg_threshold_level( dg,i,j))
dg_threshold_level_sc = lambda dg,i,j: dg_smallest_clusters( dg, dg_threshold_level( dg,i,j))

def modelsForSelectedTractographyFibers( slicer, fiberClusters, hie,  fiberBundleNode,  modelPrefix='fiber_', lineNumberArrayName='lineNumber', models=dict(), tube=False, numberModelByClusterNumber = False ):
  from Slicer import Plugin
  from random import random
  import numpy
  import sys

  names = [ (modelPrefix+'%0.4d')%i for i in (xrange(len(fiberClusters)) if not numberModelByClusterNumber else fiberClusters ) ]

  for m in ( n for n in models.keys() if n not in names ):
    slicer.MRMLScene.RemoveNode( models[m] )
    del models[ m ]

  for i in xrange(len(fiberClusters)):
    name = names[i]
  
    if name in models.keys():
      continue 
    
    cluster = hie[ fiberClusters[i] ][1]

    modelNode = slicer.vtkMRMLModelNode()
    modelDisplayNode = slicer.vtkMRMLModelDisplayNode()

    pd = slicer.vtkPolyData()
    pd.SetPoints( fiberBundleNode.GetPolyData().GetPoints() )
    pd.SetLines( slicer.vtkCellArray() )

    if tube:
      tubeFilter = slicer.vtkTubeFilter()
      tubeFilter.SetInput(pd)
      tubeFilter.Update()
      modelNode.SetAndObservePolyData( tubeFilter.GetOutput() )
    else:
      modelNode.SetAndObservePolyData( pd )

    if tube:
      print tubeFilter

    modelNode.SetName(name)

    extractSelectedTractographyFibers( slicer, cluster, fiberBundleNode, modelNode )

    slicer.MRMLScene.AddNode( modelNode )
    slicer.MRMLScene.AddNode( modelDisplayNode )
    modelNode.SetAndObserveDisplayNodeID( modelDisplayNode.GetID() )

    models[name]=modelNode

  return models

def extractSelectedTractographyFibers( slicer, linesToExtract, fiberBundleNode,  modelNode ):
  from Slicer import Plugin
  from random import random
  import numpy
  import sys

  model_pd = modelNode.GetPolyData()
  bundle_pd = fiberBundleNode.GetPolyData()
  assert( model_pd!=[] and bundle_pd!=[] )
  assert( bundle_pd.GetNumberOfCells() == bundle_pd.GetNumberOfLines() )

  if model_pd.GetPoints()!= bundle_pd.GetPoints():
    model_pd.SetPoints( bundle_pd.GetPoints() )
#  model_pd.SetLines(  slicer.vtkCellArray() )
  model_lines = model_pd.GetLines()

  model_lines.Reset()
  for lineNo in linesToExtract:
      line = bundle_pd.GetCell( lineNo )
      model_lines.InsertNextCell( line )

  model_lines.Modified()
  model_pd.Modified()
  modelNode.Modified()

  return modelNode

def lineModelToTubeModel(slicer, model, modelNamePrefix='tube_'):
  tubeFilter = slicer.vtkTubeFilter()
  tubeFilter.SetInput( model.GetPolyData() )
  modelNode = slicer.vtkMRMLModelNode()
  modelDisplayNode = slicer.vtkMRMLModelDisplayNode()
  modelNode.SetAndObservePolyData( tubeFilter.GetOutput() )
  modelNode.SetName( modelNamePrefix+model.GetName() )
  slicer.MRMLScene.AddNode( modelNode )
  slicer.MRMLScene.AddNode( modelDisplayNode )
  modelNode.SetAndObserveDisplayNodeID( modelDisplayNode.GetID() )

  color = model.GetDisplayNode().GetColor()
  modelDisplayNode.SetColor( color[0], color[1], color[2] )
  
  return (modelNamePrefix+model.GetName(),modelNode)




def fiberGPFieldToMRMLScalarVN( fib, vn ):
  from Slicer import slicer
  matrix = slicer.vtkMatrix4x4()
  matrix.SetElement(0,0,0)
  matrix.SetElement(1,1,0)
  matrix.SetElement(2,2,0)

  matrix.SetElement(0,2,2)
  matrix.SetElement(1,1,2)
  matrix.SetElement(2,0,2)

  matrix.SetElement(0,3,fib[0][0,0])
  matrix.SetElement(1,3,fib[0][0,1])
  matrix.SetElement(2,3,fib[0][0,2])
  matrix.Invert()

  vn.SetRASToIJKMatrix(matrix)

  id = vn.GetImageData()

  fb = fib[1] 
  shape = fb.shape
  id.SetDimensions( shape[-1], shape[-2], shape[-3] )
  id.SetNumberOfScalarComponents(1)
  id.AllocateScalars()

  id.ToArray()[:] = fb

  id.Modified()

  vn.Modified()



