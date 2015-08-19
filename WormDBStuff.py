from sqlalchemy import Column, Integer, Float, ForeignKey
from geoalchemy2 import Geometry

from sqlalchemy.orm import sessionmaker, relationship, backref, aliased

# This is the base of all PostGIS table names for this project
# With a little luck, all of this "by hand" construction of tablenames
# will get fixed in the worming code shortly, but for now, let's keep on doing this.

# We are building classes here that are to be "mixed in" with ORM declarative_base classes.
# The google group message that suggested this design is found at:
# <https://groups.google.com/d/msg/sqlalchemy/kjBpgm3CBSE/_6hS95YTKQAJ>
#
# These classes are responsible for declaring the DB columns "once and only once".
# The classes that use these as mixins elsewhere are responsible for dealing with the
# database specific infrastructure, including naming the tables.


# This is a class to be "mixed in" with  the "declarative base" 'Object Relational Mapper' 
# (ORM) of sqlalchemy
# It's job is to map between the database table and Python objects
# The structure essentially mimics an already existing table, or declares a new table
# if we create it here in this code. 
# In this instance, it already exists.
#FIXME Change name to WormPoints, and correct in rest of script. Or maybe think about naming in general...
class WormPointBase(object):
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
    # A duplicate of pt in WGS84 coordinates; converted by PostGIS at write-time
    wgs84_pt = Column(Geometry('POINT'),index=True)

    
class WormLevelBase(object):
    # A PK, there are only ~10 entries in this table, so it's tiny, so no index.
    worm_level_id = Column(Integer, primary_key=True)
    # The actual level (prob in meters, but potentially varies...)
    level = Column(Float)
    # Database magic that links entries in this table with entries in another table
    
class WormLevelPointsBase(object):
    # __tablename__ = levels_points_name
    # This table has a "composite primary key" composed of the first 2 ForeignKey entries and the internal primary key
    # This is the level_id in the external table
    #worm_level_id = Column(Integer, ForeignKey(levels_name + '.worm_level_id'), primary_key=True)
    # This is the point id of the END point of a line segment.
    #point_id = Column(Integer, ForeignKey(points_name + '.worm_point_id'), primary_key=True)
    # In addition to participating in a composite primary key, this field is 
    # a unique-within-a-level index for worm segments. 
    worm_seg_id = Column(Integer,primary_key=True,index=True)
    # Database magic that links entries in this table with entries in another table
    #worm_level = relationship(WormLevel, backref=backref("worm_point_assoc"))
    # Database magic that links entries in this table with entries in another table
    #worm_point = relationship(WormPoint, backref=backref("worm_level_assoc"))
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




