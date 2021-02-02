from __future__ import print_function
from astropy.coordinates import SkyCoord
from dustmaps.bayestar import BayestarQuery

coords = SkyCoord(180., 0., unit='deg', frame='galactic')

# Note that below, we could use version='bayestar2017' to get the newer
# version of the map. Note, however, that the reddening units are not
# identical in the two versions of the map. See Green et al. (2018) for
# an explanation of the units.
bayestar = BayestarQuery(max_samples=2, version='bayestar2015')

ebv = bayestar(coords, mode='random_sample')

print(ebv)
