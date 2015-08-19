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

# This is the base of all PostGIS table names for this project
# With a little luck, all of this "by hand" construction of tablenames
# will get fixed in the worming code shortly, but for now, let's keep on doing this.
basename = 'ADKMergedBGA2500_to_max_grad'
#basename = 'ADKMergedBGA2500'
#layer_name = basename
layer_name = 'ADKMergedBGA2500'
points_name = basename + '_points'
#levels_name = basename + '_levels'
levels_name = 'ADKMergedBGA2500' + '_levels'
levels_points_name = basename + '_levels_points'

earthquakes = 'adk_merged_eqs_distance_from_worms_and_eulers'

euler_points = 'ADK_BGA_Euler_Solutions'

# This code is an example of wrapping a PostGIS function that is not already wrapped via geoalchemy2
class ST_Collect(GenericFunction):
    name = 'ST_Collect'
    type = Geometry

# Originally Copied from WriteWormsToPostGIS module from BSDWormer; now modified

# sqlalchemy vodoo
Base = declarative_base()

# This is a class from the "declarative base" 'Object Relational Mapper' (ORM) of sqlalchemy
# It's job is to map between the database table and Python objects
# The structure essentially mimics an alread existing table, or declares a new table
# if we create it here in this code. 
# In this instance, it already exists.
#FIXME Change name to WormPoints, and correct in rest of script. Or maybe think about naming in general...
class WormPoint(Base):
    __tablename__ = points_name
    # Primary Key. Boring.
    worm_point_id = Column(Integer, primary_key=True, index=True)
    # An id from the worming code. I don't remember if it was unique, so I didn't use it as a PK.
    vtk_id = Column(Integer,index=True)
    # Coordinates of the point in some 'native' CRS
    x = Column(Float)
    y = Column(Float)
    z = Column(Float)
    # The scalar value of the magnitude of the horizontal gradient
    grad = Column(Float)
    # The height of upward continuation from which the grad and coordinates were drawn.
    height = Column(Float)
    # A PostGIS point geometry, in the native CRS
    pt = Column(Geometry('POINT'),index=True)
    # Database magic that links entries in this table with entries in another table
    level = relationship('WormLevel', secondary=levels_points_name)
    # A duplicate of pt in WGS84 coordinates; converted by PostGIS at write-time
    wgs84_pt = Column(Geometry('POINT'),index=True)

    
class WormLevel(Base):
    __tablename__ = levels_name
    # A PK, there are only ~10 entries in this table, so it's tiny, so no index.
    worm_level_id = Column(Integer, primary_key=True)
    # The actual level (prob in meters, but potentially varies...)
    level = Column(Float)
    # Database magic that links entries in this table with entries in another table
    point = relationship('WormPoint', secondary=levels_points_name)
    
class WormLevelPoints(Base):
    __tablename__ = levels_points_name
    # This table has a "composite primary key" composed of the first 2 ForeignKey entries and the internal primary key
    # This is the level_id in the external table
    worm_level_id = Column(Integer, ForeignKey(levels_name + '.worm_level_id'), primary_key=True)
    # This is the point id of the END point of a line segment.
    point_id = Column(Integer, ForeignKey(points_name + '.worm_point_id'), primary_key=True)
    # In addition to participating in a composite primary key, this field is 
    # a unique-within-a-level index for worm segments. 
    worm_seg_id = Column(Integer,primary_key=True,index=True)
    # Database magic that links entries in this table with entries in another table
    worm_level = relationship(WormLevel, backref=backref("worm_point_assoc"))
    # Database magic that links entries in this table with entries in another table
    worm_point = relationship(WormPoint, backref=backref("worm_level_assoc"))
    # This is an index number internal to each worm segment, numbering the edges
    # FIXME (maybe) This terminology needs to be cleaned up.
    seg_sequence_num = Column(Integer)
    # This holds the PostGIS geometry structure for a single edge, in some native CRS.
    line_segmt = Column(Geometry('LINESTRING'),index=True)
    # This scalar gradient value is derived from the average of the point grads on either end of the edge
    # Currently, the upstream code is doing that for the LOG(value), so this is in fact now
    # sqrt(grad(pt1) * grad(pt2))
    line_grad = Column(Float)
    # The azimuth of the edge in degrees East of North.
    azimuth = Column(Float)
    # This is the point ID in the points table of the starting point of an edge
    # FIXME (maybe) this could and probably should be an actual relation into the points table, for ease of retrieval.
    start_point_id = Column(Integer)
    # This is a duplicate of line_segmt but explicitly stored in wgs84.
    wgs84_line_segmt = Column(Geometry('LINESTRING'),index=True)



# Hooking things up to the database system
db = 'postgresql://frank:f00bar@localhost:5433/frank'
engine = create_engine('%s'%db, echo=False)
Session = sessionmaker(bind=engine)
session = Session()
connect = engine.connect()

if not engine.dialect.has_table(connect, layer_name):
    raise AttributeError('The Layer table is missing.')
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
euler_query = session.query(Euler).filter(Euler.depth <= 7500.)

# Turning euler points into numpy array
euler_pt_coords = np.array([[e.xeuler,e.yeuler,e.depth] for e in euler_query])

# Creating scikit-learn KDTree to speed up earthquake-euler point comparison
euler_kd = neighbors.KDTree(euler_pt_coords,leaf_size=100)


eq_query = session.query(EQs,
                         func.ST_Transform(EQs.geom,32618).ST_X(),
                         func.ST_Transform(EQs.geom,32618).ST_Y() )

# This is a "north unit vector" 
North = Vector(x=1., y=0., z=0.)
km_10_degs = kmToDegrees(10.)

# sqlalchemy voodoo, these keep aliases of tables for constructing "subqueries" 
wlp = aliased(WormLevelPoints)
wp = aliased(WormPoint)

# THE MAIN OUTER LOOP
# We are looping over everyything in point_query, with extra restrictions, ordering, and limits...


r = 10000.

# Let's build something for some quick stats...

min_dist_to_nodes = []


for p,p_lon,p_lat in eq_query.filter(EQs._Depth_km_ <= 7.5):
    
    # depth must be in meters!
    eq_pt = [p_lon,p_lat,1000.*p._Depth_km_]
    
    # New scikit_learn.neighbors implementation of the query for worms
    ww,dw = worm_kd.query_radius(eq_pt,r=r,return_distance = True,sort_results=True)
	
	# New scikit_learn.neighbors implementation of the query for euler
    we,de = euler_kd.query_radius(eq_pt,r=r,return_distance = True,sort_results=True)
	
	# Displays earthquakes outside the range of worms
    if ww[0].shape[0] == 0:
        print "No Worms within %f meters."%r
        continue
        
    # Displays earthquakes outside the range of eulers
    if we[0].shape[0] == 0:
        print "No Euler points within %f meters."%r
        continue
    
    
    min_dist_to_nodes += [[p.id,dw[0][0],de[0][0]]]
    
    
    
    sys.stdout.flush()

	# Finding a way to compare the accuracy of worm and Euler points for specific earthquakes
	
    # Option 1: write a new column to the database, and then use a GIS tool to compare locations where distance_from_worm is greater
    # than distance_from_euler or vice-versa
    p.distance_from_worms = dw[0][0]
    p.distance_from_eulers = de[0][0]


session.commit()
print "Done"
    

       

