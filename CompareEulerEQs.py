from sqlalchemy import create_engine, func, inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.schema import Table,MetaData
from sqlalchemy import Column, Integer, Float, ForeignKey
from geoalchemy2 import Geometry
from geoalchemy2.functions import GenericFunction
from geoalchemy2.elements import WKTElement, WKBElement
from sqlalchemy.orm import sessionmaker, relationship, backref, aliased
from math import sqrt, atan2, pi, log10, log, sin, cos, radians
from Scientific.Geometry import Vector
from scipy import spatial
from sklearn import neighbors
import numpy as np
import sys

#Testing things


# sqlalchemy vodoo
Base = declarative_base()

# Hooking things up to the database system
db = 'postgresql://frank:f00bar@localhost:5433/frank'
engine = create_engine('%s'%db, echo=False)
Session = sessionmaker(bind=engine)
session = Session()
connect = engine.connect()
    
meta = MetaData()


# This is a black magic function, that hooks up an existing database table, but that still allows
# for python object access to the database data. 
# We will hook up the Euler solution points
class ADKBGAEuler(Base):
	__table__ = Table('ADK_BGA_Euler_Solutions', meta, autoload=True, autoload_with=engine)

# We will hook up the earthquake hypocenters
class ADKMergedEQs(Base):
    __table__ = Table('adk_merged_eqs', meta, autoload=True, autoload_with=engine)

# A function that converts latitude and longitudes (in degrees)
# for 2 different points into Great Circle distances in kilometers.
def gc_dist(lat1,lon1,lat2,lon2):
    # cribbed from <http://code.activestate.com/recipes/
    # 576779-calculating-distance-between-two-geographic-points/>
    # Radius of a sphere with the equivalent volume to the Earth
    R = 6371.0
    lat1 = radians(lat1)
    lon1 = radians(lon1)
    lat2 = radians(lat2)
    lon2 = radians(lon2)
    
    dlon = (lon2 - lon1)
    dlat = (lat2 - lat1)
    a = (sin(dlat/2.))**2 + cos(lat1) * cos(lat2) * (sin(dlon/2.))**2
    c = 2. * atan2(sqrt(a), sqrt(1.-a))
    return R * c



# Utility function: how many degrees away is something km apart on the surface of the Earth
def kmToDegrees(km):
    # 6371 is again the radius of the Earth
    return 360. * km / (6371.*2.*pi)



# Pulling in euler points
euler_query = session.query(ADKBGAEuler).filter(ADKBGAEuler.depth <= 7500.)

# Turning euler points into numpy array
euler_pt_coords = np.array([[e.xeuler,e.yeuler,e.depth] for e in euler_query])

# Creating scikit-learn KDTree to speed up earthquake-euler point comparison
#euler_kd = neighbors.KDTree(euler_pt_coords,leaf_size=100)

eq_query = session.query(ADKMergedEQs,
                         func.ST_Transform(ADKMergedEQs.geom,32618).ST_X(),
                         func.ST_Transform(ADKMergedEQs.geom,32618).ST_Y() )

km_10_degs = kmToDegrees(10.)


r = 10000.


min_dist_to_nodes = []
eq_depths = []
	
for p,p_lon,p_lat in eq_query:
#.filter(ADKMergedEQs._Depth_km_ <= 7.5):
    
    # depth must be in meters!
    #eq_pt = [p_lon,p_lat,1000.*p._Depth_km_]
    
    # New scikit_learn.neighbors implementation of the query
    #wq,dq = euler_kd.query_radius(eq_pt,r=r,return_distance = True,sort_results=True)
    
    # Displays earthquakes outside the range
    #if wq[0].shape[0] == 0:
    #    print "No Euler points within %f meters."%r
    #    continue
    #min_dist_to_nodes += [dq[0][0]]
    
    if type(p._Depth_km_) != float:
    	print p._Depth_km_
    	continue
    
    eq_depths += [p._Depth_km_]
    
    #sys.stdout.flush()
    #print 'NEW EARTHQUAKE'
    
print "Done"



