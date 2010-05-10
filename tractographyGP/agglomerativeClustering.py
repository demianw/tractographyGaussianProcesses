def innerProductHierarchicalClustering( innerProductMatrix, checkPointFile=None, checkPointCount=0, loadFromCheckpoint=False, normalized=False):
  import sys 
  import os
  import numpy
  from heapq import heappush, heappop


  innerProductMatrix = innerProductMatrix.astype('float32')

  try:
    from scipy import sparse 
    if  sparse.issparse(innerProductMatrix):
      innerProductMatrix = innerProductMatrix.todense()
  except ImportError:
    pass
  except:
    raise

  innerProductHeap = []

  cycles = 0
  if not (loadFromCheckpoint and os.path.exists(checkPointFile)):


	  print 'Initializing'
         
          N = innerProductMatrix.shape[0]

          norms = numpy.zeros( N )
          norms[:N] = numpy.sqrt(innerProductMatrix.diagonal())

	  active = set(range(N))
	  dendogram = map( lambda i:(0,(i,),),range(N))


          for (i,j) in numpy.transpose(innerProductMatrix.nonzero()).squeeze():
            if (i>j):
              heappush( innerProductHeap, (-innerProductMatrix[i,j]/(1. if not normalized else (norms[i]*norms[j]) ), (i,j)))
          print 'Loaded InnerProducts from matrix'
          sys.stdout.flush()


	  if checkPointFile!=None:
		print '%s checkpointing %d'%(checkPointFile,cycles)
                if  not normalized:
                  cPickle.dump( (active,dendogram,innerProductHeap,innerProductMatrix),open( checkPointFile,'w'))
                else:
                  cPickle.dump( (active,dendogram,innerProductHeap,innerProductMatrix,norms),open( checkPointFile,'w'))

  else:
	print 'Loading from checkpoint'
        if not normalized:
          (active,dendogram,innerProductHeap,innerProductMatrix)=cPickle.load( open( checkPointFile))
        else:
          (active,dendogram,innerProductHeap,innerProductMatrix,norms)=cPickle.load( open( checkPointFile))


  print 'Starting the clustering'

  newIndex=0
  nextFreeElment = CircularBuffer( innerProductMatrix.shape[0] )
  elementMap =  numpy.arange( innerProductMatrix.shape[0], dtype='int' )
  newLine_0 = numpy.zeros( innerProductMatrix.shape[0] )
  newLine_1 = numpy.zeros( innerProductMatrix.shape[0] )
  newLine =  numpy.zeros( innerProductMatrix.shape[0] )

  while len(active)>1:
    
    while len(innerProductHeap)>0:
      nextToJoin = heappop( innerProductHeap )
      nextToJoin_0 = nextToJoin[1][0]
      nextToJoin_1 = nextToJoin[1][1]
      if ( nextToJoin_0 in active) and (nextToJoin_1  in active ):
        break

    if len(innerProductHeap)==0:
      break

    active.remove( nextToJoin_0 )
    active.remove( nextToJoin_1 )

    dn_0 = dendogram[ nextToJoin_0 ]
    dn_1 = dendogram[ nextToJoin_1 ]
    clusteredElements = dn_0[1]+dn_1[1]
    n_0 = len( dn_0[1] )
    n_1 = len( dn_1[1] )

    newIndex = len(dendogram) 
    print "%6d: Joining %6d and %6d #elements %6d, ip:%f"%(newIndex,nextToJoin_0,nextToJoin_1,len(clusteredElements),-nextToJoin[0])

    nextToJoin_0_index = ( elementMap== nextToJoin_0 ).nonzero()[0][0]
    nextToJoin_1_index = ( elementMap== nextToJoin_1 ).nonzero()[0][0]

    nextFreeElment.append( nextToJoin_0_index )
    nextFreeElment.append( nextToJoin_1_index ) 

    newIndex_index = nextFreeElment.pop()
    elementMap[ newIndex_index ] = newIndex

    newLine_0[:] = n_0*innerProductMatrix[ nextToJoin_0_index,:]
    newLine_1[:] = n_1*innerProductMatrix[ nextToJoin_1_index,:]
    newLine[:] = (newLine_0+newLine_1)

    newLine[newIndex_index] = (n_0*innerProductMatrix[ nextToJoin_0_index,nextToJoin_0_index]+n_1*innerProductMatrix[ nextToJoin_1_index,nextToJoin_1_index ]+(n_0+n_1)*innerProductMatrix[ nextToJoin_0_index,nextToJoin_1_index ])

    elementMap[ newIndex_index ] = newIndex

    innerProductMatrix[ newIndex_index,: ] = newLine
    innerProductMatrix[ newIndex_index,: ]/=len(clusteredElements)
    innerProductMatrix[ :, newIndex_index ]= innerProductMatrix[ newIndex_index , :].T

    norms[newIndex_index] = numpy.sqrt(innerProductMatrix[newIndex_index,newIndex_index])

    candidates =  elementMap[ newLine.nonzero() ]
    for i in active.intersection( candidates ):
      pos = ( elementMap== i ).nonzero()[0][0]
      assert( innerProductMatrix[newIndex_index, pos ]!=0)
      innerP = innerProductMatrix[newIndex_index,pos]/(1. if not normalized else  (norms[pos]*norms[newIndex_index]) )
      heappush( innerProductHeap, (-innerP, ( i , newIndex )))

    active.add( newIndex )
    dendogram.append( (max(dn_0[0],dn_1[0])+1, clusteredElements,(nextToJoin_0,nextToJoin_1),-nextToJoin[0])  )

    if  checkPointCount>0 and checkPointFile!=None and cycles%checkPointCount==0:
        print '%s checkpointing %d'%(checkPointFile,cycles); sys.stdout.flush()
        if not normalized:
          cPickle.dump( (active,dendogram,innerProductHeap,innerProductMatrix),open( checkPointFile,'w'))
        else:
          cPickle.dump( (active,dendogram,innerProductHeap,innerProductMatrix,norms),open( checkPointFile,'w'))
    cycles+=1


  return dendogram,innerProductMatrix,elementMap


class CircularBuffer:
    
    _array = None
    _act = 0
    _size = 0

    def __init__(self,n):
        import numpy
        self._array = numpy.zeros(n,dtype='int')
        self._act = 0
        self._size = 0

    def pop(self):
        if self._size==0:
            raise ValueErrorException()
        else:
            self._act=(self._act+1)%len(self._array)
            self._size-=1
            return self._array[self._act-1]
    def append(self, n ):
        if self._size==len(self._array):
            raise MemoryError()
        else:
            self._array[ (self._act+self._size)%len(self._array)]=n
            self._size+=1

    def __len__(self):
        return self._size


