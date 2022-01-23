"""Define the data model for the p_ms coefficient expansions."""

from sqlalchemy import\
    Column,\
    Integer,\
    Float,\
    TIMESTAMP,\
    ForeignKey

from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

DataModelBase = declarative_base()

class Interpolations(DataModelBase):
    """The available coefficient expansions."""
    
    __tablename__ = 'interpolations'
    
    id = Column(
        Integer,
        primary_key=True,
        doc='A unique identifier for each coefficient.'
    )
    
    m = Column(
        Integer,
        nullable=False,
        doc='The mth coefficient; can be -2, 0, or 2.'
    )
    s = Column(
        Integer,
        nullable=False,
        server_default='',
        doc='The sth coefficient; all integer values 0 and up.'
    )
    min_interp_e = Column(
        Float,
        nullable=False,
        doc='The smallest value of e which the process was extended to.'
    )
    number_of_steps = Column(
        Integer,
        nullable=False,
        doc='The number of steps between the last directly calculated zero and 1 (inclusive).'
    )
    max_checked_e = Column(
        Float,
        nullable=False,
        doc='The highest value of e for which we compared directly calculated and interpolated values.'
    )
    interp_accuracy = Column(
        Float,
        nullable=False,
        doc='The maximum difference between calculated and interpolated, as well as the point below which a calculated value was considered zero.'
    )
    timestamp = Column(
        TIMESTAMP,
        nullable=False,
        doc='When was this record last changed.'
    )

    interpData = relationship('InterpolationData', back_populates='interps')

class InterpolationData(DataModelBase):
    """The Chebyeshev coefficients for each expansion."""

    __tablename__ = 'interpolation_data'

    p_id = Column(
        Integer,
        ForeignKey('interpolations.id',
                   onupdate='CASCADE',
                   ondelete='RESTRICT'),
        primary_key=True,
        doc='The expansion this coefficient corresponds to.'
    )
    step_number = Column(
        Integer,
        primary_key=True,
        doc='The position between the final zero and one for this value.'
    )
    y_value = Column(
        Float,
        primary_key=False,
        doc='The value of p_ms at this position.'
    )
    timestamp = Column(
        TIMESTAMP,
        nullable=False,
        doc='When was this record last changed.'
    )
    
    interps = relationship('Interpolations', back_populates='interpData')
