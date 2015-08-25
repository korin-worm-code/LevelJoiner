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

basename = 'ravat_ADK_PSG1250'

WormPoint, WormLevelPoints, WormLevel, tablenames = WormDBStuffFactory(basename)

earthquakes = 'merged_ta_neic_eqs'

euler_points = 'adk_psg_euler_new'

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

# We will hook up the Euler solution points
class Euler(Base):
	__table__ = Table(euler_points, meta, autoload=True, autoload_with=engine)

# We will hook up the earthquake hypocenters
class EQs(Base):
    __table__ = Table(earthquakes, meta, autoload=True, autoload_with=engine)






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


# Trying the new scikit-learn implementation of KDTree
worm_kd = neighbors.KDTree(worm_pt_coords,leaf_size=100)


# Pulling in euler points
euler_query = session.query(Euler).filter(Euler.depth <= 15000.)

# Turning euler points into numpy array
euler_pt_coords = np.array([[e.x_euler,e.y_euler,e.depth] for e in euler_query])

# Creating scikit-learn KDTree to speed up earthquake-euler point comparison
euler_kd = neighbors.KDTree(euler_pt_coords,leaf_size=100)


eq_query = session.query(EQs,
                         EQs.geom.ST_X(),
                         EQs.geom.ST_Y() )




# THE MAIN OUTER LOOP
# We are looping over everyything in point_query, with extra restrictions, ordering, and limits...


r = 10000.

# Let's build something for some quick stats...

min_dist_to_nodes = []

distance_to_worms = []
distance_to_eulers = []

for p,p_lon,p_lat in eq_query.filter(EQs._DepthMeters_ <= 15000,EQs._DepthMeters_ != 0.,EQs._DepthMeters != 1000.,EQs._DepthsMeters_ != 5000.):
    
    # depth must be in meters!
    eq_pt = [p_lon,p_lat,p._DepthMeters_]
    
    # New scikit_learn.neighbors implementation of the query for worms
    ww,dw = worm_kd.query_radius(eq_pt,r=r,return_distance = True,sort_results=True)
	
	# New scikit_learn.neighbors implementation of the query for euler
    we,de = euler_kd.query_radius(eq_pt,r=r,return_distance = True,sort_results=True)
	
	# Displays earthquakes outside the range of worms
    if ww[0].shape[0] == 0:
    #    print "No Worms within %f meters."%r
        continue
        
    # Displays earthquakes outside the range of eulers
    if we[0].shape[0] == 0:
    #    print "No Euler points within %f meters."%r
        continue
    
    
    min_dist_to_nodes += [[p.id,dw[0][0],de[0][0]]]
    
    
    sys.stdout.flush()

	# Finding a way to compare the accuracy of worm and Euler points for specific earthquakes
	
    # Option 1: write a new column to the database, and then use a GIS tool to compare locations where distance_from_worm is greater
    # than distance_from_euler or vice-versa
    #p.distance_to_worms = dw[0][0]
    #p.distance_to_eulers = de[0][0]
    distance_to_worms += [dw[0][0]]
    distance_to_eulers += [de[0][0]]


#session.commit()
print "Done"
    

       

