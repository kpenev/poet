#!/usr/bin/python3 -u

"""Define the data model for stellar_evolution_manager sqlalchemy based."""

from sqlalchemy import\
    Column,\
    Integer,\
    String,\
    Float,\
    ForeignKey,\
    create_engine
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

DataModelBase = declarative_base()

interpolator_nodes_association_table = 

class SerializedInterpolator(DataModelBase) :
    """The current set of serialized stellar evolution interpolators."""

    __tablename__ = 'interpolators'

    id = Column(
        Integer,
        primary_key = True,
        doc = 'A uniquie identifier for each serialized interpolator'
    )
    filename = Column(
        String,
        nullable = False,
        doc = 'The name of the file containing the serialized interpolator.'
    )
    nodes           = relationship("Nodes")
    smoothing       = relationship("Smoothing")
    tracks          = relationship("Track")

class Nodes(DataModelBase) :
    """The number of interpolation nodes used by each interpolator."""
    __tablename__ = 'nodes'

    interpolator_id = Column(
        Integer,
        ForeignKey('interpolators.id',
                   onupdate = 'CASCADE',
                   ondelete = 'RESTRICT'),
        doc = 'The interpolator this set of nodes applies to.'
    )
    quantity_id = Column(
        Integer,
        ForeignKey('quantities.id',
                   onupdate = 'CASCADE',
                   ondelete = 'RESTRICT'),
        primary_key = True,
        doc = 'The ID of the quantity to which this set of nodes applies.'
    )
    value = Column(
        name = 'nodes',
        type_ = Integer,
        nullable = False,
        doc = 'The number of nodes used when constructing the interpolator '
        'for the given quantity.'
    )
    quantity = relationship("Quantity")

    def __repr__(self) :
        return repr(self.quantity) + ':' + repr(self.value)

class Smoothing(DataModelBase) :
    """The smoothing arguments used by each interpolator."""
    __tablename__ = 'smoothing'

    interpolator_id = Column(
        Integer,
        ForeignKey('interpolators.id',
                   onupdate = 'CASCADE',
                   ondelete = 'RESTRICT'),
        primary_key = True,
        doc = 'The interpolator this set of smoothing arguments applies to.'
    )
    quantity_id = Column(
        Integer,
        ForeignKey('quantities.id',
                   onupdate = 'CASCADE',
                   ondelete = 'RESTRICT'),
        primary_key = True,
        doc = 'The ID of the quantity to which this set of smoothing '
        'arguments applies.'
    )
    value = Column(
        name = 'smoothing',
        type_ = Float,
        nullable = False,
        doc = 'The smoothing argument used when constructing the '
        'interpolator for the given quantity.'
    )
    quantity = relationship("Quantity")

class Track(DataModelBase) :
    """The grid of tracks used by each interpolator."""

    __tablename__ = 'track_grids'

    grid_id = Column(
        Integer,
        primary_key = True,
        doc = 'A unique identifier for each grid of tracks.'
    )
    interpolator_id = Column(
        Integer,
        ForeignKey('interpolators.id',
                   onupdate = 'CASCADE',
                   ondelete = 'RESTRICT'),
        doc = 'The interpolator based on this grid of tracks.'
    )
    mass = Column(
        Float,
        primary_key = True,
        doc = 'The stellar mass of the track contained in the given file.',
    )
    metallicity = Column(
        Float,
        primary_key = True,
        doc = 'The metallicity ([Fe/H])of the track contained in the given '
        'file.',
    )
    checksum = Column(
        String,
        nullable = False,
    )

class Quantity(DataModelBase) :
    """The set of quantities tracked by stellar evolution."""

    __tablename__ = 'quantities'

    id = Column(
        Integer,
        primary_key = True,
        doc = 'A unique identifier for each quantity.'
    )
    name = Column(
        String,
        nullable = False,
        doc = 'The name used for this quantity by the stellar_evolution '
        'module.'
    )

    def __repr__(self) :
        return self.name + '(%d)' % self.id
