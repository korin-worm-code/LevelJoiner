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

from WormDBStuff.WormDBStuff import WormDBStuffFactory

#Testing things

basename = 'ADKMergedBGA2500'

euler_points = 'adk_bga_euler_new'


WormPoint, WormLevelPoints, WormLevel, tablenames = WormDBStuffFactory(basename,to_max_grad = True)


# sqlalchemy vodoo
Base = declarative_base()



# Hooking things up to the database system
db = 'postgresql://frank:f00bar@localhost:5433/frank'
engine = create_engine('%s'%db, echo=False)
Session = sessionmaker(bind=engine)
session = Session()
connect = engine.connect()

if not engine.dialect.has_table(connect, points_name):
    raise AttributeError('The Points table is missing.')
if not engine.dialect.has_table(connect, levels_name):
    raise AttributeError('The Levels table is missing.')
if not engine.dialect.has_table(connect, levels_points_name):
    raise AttributeError('The Levels_Points table is missing.')
    
meta = MetaData()


# This is a black magic function, that hooks up an existing database table, but that still allows
# for python object access to the database data. 
# We will hook up the earthquake hypocenters (not valid anymore)
class Eulers(Base):
    __table__ = Table(euler_points, meta, autoload=True, autoload_with=engine)



# This is an example of the sqlalchemy way to encapsulate a SQL query.
# This particular query builds a database "join" (perhaps not exactly due to the sqlalchemy innards)
# where all entities returned will be the edge "end point" and "edge" data structures that match.
# This is actually the head end of more restrictive filterings of the database tables
point_query = session.query(WormPoint,WormLevelPoints).filter(WormPoint.worm_point_id == WormLevelPoints.point_id)

# Pull all worm data structures from the database; 
# returns both WormPoint and WormLevelPoints as a tuple(?) for each item
all_worm_points = point_query.all()

# It's actually simpler to dig the relevant bits out from the data structures returned by the database now
# than trying to deal with the headache of getting all of the indexing correct everywhere else.
# Think of it as a "once and only once" for getting the bloody indexing right...

# Build an array of 3-coords for each worm point to feed into the kd-tree for indexing
worm_pt_coords = np.array([[w[0].x,w[0].y,w[0].z] for w in all_worm_points])


# Creating an array out of the worm levels
worm_sgmt_levels = np.array([w[1].worm_level_id for w in all_worm_points])
# Creating an array out of the worm segments
worm_sgmt_ids = np.array([w[1].worm_seg_id for w in all_worm_points])
# Creating an array out of the sequential worm pieces
worm_sgmt_seq_num = np.array([w[1].seg_sequence_num for w in all_worm_points])


# We are building a numpy record array so that we can sort them with auxiliary sorting order.
worm_rec = np.rec.fromarrays([worm_sgmt_levels, worm_sgmt_ids, worm_sgmt_seq_num])


# Now create the ndarray of the results from the query. 
# N.B. Both the end point and the edge are contained in each element.
all_worm_data = np.array(all_worm_points,dtype=[('worm_point',WormPoint),('worm_level_points',WormLevelPoints)])


# Creating SciPy KDTree to speed up earthquake-worm point comparison
#worm_kd = spatial.KDTree(worm_pt_coords,leafsize=50)
# Updating to be runable with mag data
worm_kd = neighbors.KDTree(worm_pt_coords,leaf_size=100)

# Pulling in the Euler points from the database
euler_query = session.query(Eulers)



# This is the distance we are searching within, in meters
r = 10000.

# Let's build something for some quick stats...

min_dist_to_nodes = []
#far_eq = []



for p in euler_query.filter(Eulers.depth <= 15000.):
	# We are no longer working with earthquakes, so we don't need to sort them by magnitude
	#.filter(ADKMergedEQs._Depth_km_ == 0.).order_by(ADKMergedEQs._Magnitude_):
    #print p._latitude_, p._longitude_, p._depth_km_, p._magnitude_
    
    # depth must be in meters!
    euler_pt = [p.x_euler,p.y_euler,p.depth]
    
    # SciPy KDTrees
    #dq,wq = worm_kd.query(euler_pt,k=20,distance_upper_bound=r)
    wq,dq = worm_kd.query_radius(euler_pt,r=r,return_distance = True,sort_results=True)
    
    
    # New return style
    if wq[0].shape[0] == 0:
    #    print "No Worms within %f meters."%r
        continue
    
    # Distance to the closest worm point
    min_dist_to_nodes += [dq[0][0]]
    

    sys.stdout.flush()
    


#session.commit()







