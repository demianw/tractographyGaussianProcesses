class DendrogramVolumeLabelAnalysis:

  def __init__(self, tractography,  atlasImage, atlasRASToIJK, dendrogram):
    from numpy import asmatrix, asarray, zeros, ndindex,swapaxes
    self._dendrogram = dendrogram
    self._tractography = tractography

    self._matRASToIJK = atlasRASToIJK
    self._labelsImage = atlasImage


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




#Dendrogram (dg) thresholding functions
dg_threshold_level = lambda dg,i,j: filter( lambda d:  ( dg[d][0]>=i and dg[d][0]<=j ) , xrange(len(dg)))
dg_threshold_energy = lambda dg,i,j: filter( lambda d: len(dg[d])>=3 and ( dg[d][3]>i and dg[d][3]<j ) , xrange(len(dg)))
dg_biggest_clusters = lambda dg,indexes: filter( lambda y: all( [ len(set(dg[y][1]).intersection( dg[i][1]))==0 for i in indexes if dg[i][-1]<dg[y][-1] ] ),  indexes  ) 
dg_smallest_clusters = lambda dg,indexes: filter( lambda y: all( [ len(set(dg[y][1]).intersection( dg[i][1]))==0 for i in indexes if dg[i][-1]>dg[y][-1] ] ),  indexes  ) 


dg_threshold_energy_bc = lambda dg,i,j: dg_biggest_clusters( dg, dg_threshold_energy( dg,i,j))
dg_threshold_energy_sc = lambda dg,i,j: dg_smallest_clusters( dg, dg_threshold_energy( dg,i,j))

dg_threshold_level_bc = lambda dg,i,j: dg_biggest_clusters( dg, dg_threshold_level( dg,i,j))
dg_threshold_level_sc = lambda dg,i,j: dg_smallest_clusters( dg, dg_threshold_level( dg,i,j))


