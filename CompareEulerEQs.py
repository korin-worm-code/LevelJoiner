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

euler_points = 'adk_bga_euler_new'

earthquakes = 'merged_ta_neic_eqs'

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
class Eulers(Base):
	__table__ = Table(euler_points, meta, autoload=True, autoload_with=engine)

# We will hook up the earthquake hypocenters
class EQs(Base):
    __table__ = Table(earthquakes, meta, autoload=True, autoload_with=engine)

# Pulling in euler points
euler_query = session.query(Eulers).filter(Eulers.depth <= 15000.)

# Turning euler points into numpy array
euler_pt_coords = np.array([[e.x_euler,e.y_euler,e.depth] for e in euler_query])

# Creating scikit-learn KDTree to speed up earthquake-euler point comparison
euler_kd = neighbors.KDTree(euler_pt_coords,leaf_size=100)

eq_query = session.query(EQs,
                         EQs.geom.ST_X(),
                         EQs.geom.ST_Y() )



r = 10000.


min_dist_to_nodes = []
eq_depths = []
	
for p,p_lon,p_lat in eq_query.filter(EQs._DepthMeters_ <= 15000., EQs._DepthMeters_ != 0.):
    
    if type(p._DepthMeters_) != float:
    	print p._DepthMeters_
       	continue
    
    # depth must be in meters!
    eq_pt = [p_lon,p_lat,p._DepthMeters_]
    
    # New scikit_learn.neighbors implementation of the query
    wq,dq = euler_kd.query_radius(eq_pt,r=r,return_distance = True,sort_results=True)
    
    # Displays earthquakes outside the range
    if wq[0].shape[0] == 0:
        print "No Euler points within %f meters."%r
        continue
    
    min_dist_to_nodes += [dq[0][0]]
    

    
    if p._DepthMeters_ == 0.:
    	print "Zero depth!"
    
    eq_depths += [p._DepthMeters_]
    
    sys.stdout.flush()
    #print 'NEW EARTHQUAKE'
    
print "Done"



