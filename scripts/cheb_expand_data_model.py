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

class MAndSToAccuracy(DataModelBase):
	"""The available coefficient expansions."""
	
	__tablename__ = 'm_and_s_to_accuracy'
	
	id = Column(
		Integer,
		primary_key=True,
		doc='A unique identifier for each combination of coefficient and accuracy.'
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
	accuracy = Column(
        Float,
        nullable=False,
        doc='The accuracy the coefficient was expanded to.'
	)
	timestamp = Column(
        TIMESTAMP,
        nullable=False,
        doc='When was this record last changed.'
	)

	chebCoeffs = relationship('ChebExpansionCoeffs', back_populates='MAndSToAccuracys')

class ChebExpansionCoeffs(DataModelBase):
    """The Chebyeshev coefficients for each expansion."""

    __tablename__ = 'cheb_expansion_coeffs'

    p_id = Column(
        Integer,
        ForeignKey('m_and_s_to_accuracy.id',
                   onupdate='CASCADE',
                   ondelete='RESTRICT'),
        primary_key=True,
        doc='The expansion this coefficient corresponds to.'
    )
    place_in_expansion = Column(
        Integer,
        primary_key=True,
        doc='The position in the overall expansion the coefficient should be placed.'
    )
    coefficient_value = Column(
        Float,
        primary_key=False,
        doc='The value of the coefficient.'
    )
    timestamp = Column(
        TIMESTAMP,
        nullable=False,
        doc='When was this record last changed.'
    )
	
    MAndSToAccuracys = relationship('MAndSToAccuracy', back_populates='chebCoeffs')
