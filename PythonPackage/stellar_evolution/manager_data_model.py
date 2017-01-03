#!/usr/bin/python3 -u

"""Define the data model for stellar_evolution_manager sqlalchemy based."""

from sqlalchemy import\
    Column,\
    Integer,\
    String,\
    Boolean,\
    Numeric,\
    ForeignKey,\
    create_engine,\
    Table
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base
from uuid import uuid4 as get_uuid

DataModelBase = declarative_base()

interpolator_tracks_table = Table(
    'interpolator_tracks',
    DataModelBase.metadata,
    Column('interpolator_id',
           Integer,
           ForeignKey('interpolators.id',
                      onupdate = 'CASCADE',
                      ondelete = 'RESTRICT')),
    Column('track_id',
           Integer,
           ForeignKey('tracks.id',
                      onupdate = 'CASCADE',
                      ondelete = 'RESTRICT'))
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

    def __str__(self) :
        return self.name + '(%d)' % self.id

class ModelSuite(DataModelBase) :
    """The various stellar evolution suites (e.g. MESA, YREC, ...)."""

    __tablename__ = 'model_suites'

    id = Column(
        Integer,
        primary_key = True,
        doc = 'A unique identifier for the suite.'
    )
    name = Column(
        String,
        nullable = False,
        doc = 'The name of the suite.'
    )

    def __str__(self) :
        return self.name + '_%d' % self.id

class Track(DataModelBase) :
    """The available stellar evolution tracks."""

    __tablename__ = 'tracks'

    id = Column(
        Integer,
        primary_key = True,
        doc = 'A unique identifier for each track.'
    )
    filename = Column(
        String,
        nullable = False,
        unique = True,
        index = True,
        doc = 'The absolute name of the file containing the track.'
    )
    mass = Column(
        Numeric(5, 3),
        primary_key = True,
        doc = 'The stellar mass of the track contained in the given file.',
    )
    metallicity = Column(
        Numeric(5, 3),
        primary_key = True,
        doc = 'The metallicity ([Fe/H])of the track contained in the given '
        'file.',
    )
    model_suite_id = Column(
        Integer,
        ForeignKey('model_suites.id',
                   onupdate = 'CASCADE',
                   ondelete = 'RESTRICT'),
        doc = 'The stellar evolution suite used to calculate the tracks.'
    )
    checksum = Column(
        String,
        nullable = False,
    )
    suite = relationship('ModelSuite')

    def __str__(self) :
        """Human readable representation."""

        return (repr(self.suite)
                +
                '(M = %g Msun, [Fe/H] = %g) [%s]' % (self.mass,
                                                     self.metallicity,
                                                     self.checksum))

class InterpolationParameters(DataModelBase) :
    """ The nodes and smoothing arguments for a given interpolator/quantity."""

    __tablename__ = 'interpolation_parameters'

    interpolator_id = Column(
        Integer,
        ForeignKey('interpolators.id',
                   onupdate = 'CASCADE',
                   ondelete = 'RESTRICT'),
        primary_key = True,
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
    nodes = Column(
        name = 'nodes',
        type_ = Integer,
        nullable = False,
        doc = 'The number of nodes used when constructing the interpolator '
        'for the given quantity.'
    )
    smoothing = Column(
        name = 'smoothing',
        type_ = Numeric(5, 3),
        nullable = True,
        doc = 'The smoothing argument used when constructing the '
        'interpolator for the given quantity.'
    )
    vs_log_age = Column(
        name = 'vs_log_age',
        type_ = Boolean,
        nullable = False,
        doc = 'Is the given quantity interpolation vs log(age) instead of '
        'age?'
    )
    log_quantity = Column(
        name = 'log_quantity',
        type_ = Boolean,
        nullable = False,
        doc = 'Is the log of the given quantity interpolated instead of '
        'the quantity?'
    )
    quantity = relationship('Quantity')
    interpolator = relationship('SerializedInterpolator',
                                back_populates = 'parameters')

    def __str__(self) :
        return (repr(self.quantity)
                +
                '(%d nodes, %g smoothing)' % (self.nodes, self.smoothing))

class SerializedInterpolator(DataModelBase) :
    """The current set of serialized stellar evolution interpolators."""

    __tablename__ = 'interpolators'

    id = Column(
        Integer,
        primary_key = True,
        doc = 'A uniquie identifier for each serialized interpolator'
    )
    name = Column(
        String,
        nullable = False,
        unique = True,
        index = True,
        doc = 'A unique human-readable name assigned to this interpolator.'
    )
    filename = Column(
        String,
        nullable = False,
        doc = 'The name of the file containing the serialized interpolator.'
    )
    checksum = Column(
        String,
        nullable = False,
        doc = 'A checksum of the serialized interpolator file.'
    )
    parameters = relationship('InterpolationParameters',
                              back_populates = 'interpolator')
    tracks = relationship("Track",
                          secondary = interpolator_tracks_table)

    def __str__(self) :
        return ('%s_%d(%s), Parameters: %s, Tracks: %s'
                %
                (self.name,
                 self.id,
                 self.filename,
                 '; '.join([repr(p) for p in self.parameters]),
                 '; '.join([repr(t) for t in self.tracks])))
