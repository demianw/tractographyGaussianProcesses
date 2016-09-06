import tractographyGP as tgp
from tract_querier import tractography
import numpy

# Load the tractography file
tr = tractography.tractography_from_file('MNI_tq_01248_af.left.vtk')

# Subsample files to 20 points per tract
tr.subsample_tracts(20)


# Generate the GP representation for each tract
gps = [tgp.fiberGPAnalysis.fiberGP(f) for f in tr.tracts()]


# Generate the inner product matrix
ips = numpy.zeros((len(gps), len(gps)), dtype='float32')
for i in xrange(len(gps)):
    for j in xrange(0, i+1):
        ips[i, j] = (
            tgp.fiberGPAnalysis.innerProduct_thinPlate_R3(gps[i], gps[j])
        )
        ips[j, i] = ips[i, j]
    print(
        "Inner product matrix: %3.0f %% done (%d of %d)" %
        (((i+1)**2.)/(len(gps)**2.)*100, i, len(gps))
    )


# Produce the clustering tree.
# The format of `hie` is a list of tuples
# (
#    node height in the tree,
#    tracts joined on the left,
#    tracts joined on the right,
#    energy of the merge
# )
hie = tgp.agglomerativeClustering.innerProductHierarchicalClustering(ips)[0]
