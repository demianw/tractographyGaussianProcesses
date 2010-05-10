import readVtkPolyData_ext as _rvtk

_header = """# vtk DataFile Version 3.0
vtk output
ASCII
"""
_polyDataType = "DATASET POLYDATA\n"
_pointsHeader = "POINTS %d float\n"
_linesHeader = "LINES %d %d\n"
_pointDataHeader = "POINT_DATA %d\n"
_pointDataAttributeHeader = "%s %s float"

def readVtkPolyData(fileName, xmlFormat=False ):
  if xmlFormat:
    return _rvtk.c_pdrXML(fileName)
  else:
    return _rvtk.c_pdr(fileName)

def writeLinesToVtkPolyData( fileName, lines, pointData={} ):
  f = open(fileName,'w')
  f.write(_header)
  f.write(_polyDataType)

  numberOfPoints = sum( [ len(l) for l in lines ] )

  f.write(_pointsHeader%numberOfPoints)
  for line in lines:
    for point in line:
      f.write(str(point).strip()[1:-1]+'\n')

  numberOfLines = len(lines)
  f.write(_linesHeader%(numberOfLines,numberOfLines+numberOfPoints))
  pointsForLineSaved = 0
  for line in lines:
    f.write("%d %s \n"%(len(line),reduce( lambda x,y:x+' %d'%(y+pointsForLineSaved), xrange(len(line)),'' )))
    pointsForLineSaved+=len(line)

  if pointData:
    f.write(_pointDataHeader%numberOfPoints)
    
    print pointData.keys()
    for key,data in pointData.items():
      if key.startswith('vtk'):
        title = key[3:].upper()
        name = key[3:].lower()
      else:
        title = 'SCALARS'
        name = key


      
      if title == 'SCALARS':
        numberOfComponents = len(data[0][0])
        f.write(_pointDataAttributeHeader%(title,name)+' %d\n'%numberOfComponents)
        f.write('LOOKUP_TABLE default\n')
      else:
        f.write(_pointDataAttributeHeader%(title,name)+' %d\n')

      for line in data:
        for attribute in line:
          f.write(str(attribute).strip()[1:-1]+'\n')

  
  f.flush()
  f.close()

