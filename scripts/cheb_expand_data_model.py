"""Define the data model for class records."""
# I'm starting from the phys2325 model, forgive any irregularities please


# sqlite browser     does it come by default or do I install it separately?
#database_interface.py 
#initialize_database.py 
#poet: python, stellar_ev, manager.py (also manager data model

from sqlalchemy import\
    Column,\
    Integer,\
    String,\
    Float,\
    Date,\
    TIMESTAMP,\
    ForeignKey,\
    Index,\
    ForeignKeyConstraint

from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

#pylint false positive: this is actually a class name
#pylint: disable=invalid-name
DataModelBase = declarative_base()
#pylint: enable=invalid-name

#The standard use of SQLAlchemy ORM requires classes with no public methods.
#pylint: disable=too-few-public-methods

class PmsVAcc(DataModelBase):
	"""The available coefficient expansions."""
	
	__tablename__ = 'pVacc'
	
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
	accur = Column(
        Float,
        nullable=False,
        doc='The accuracy the coefficient was expanded to.'
	)
	timestamp = Column(
        TIMESTAMP,
        nullable=False,
        doc='When was this record last changed.'
	)

	chebCoeffs = relationship('ChebCo', back_populates='pmsvaccs')

class ChebCo(DataModelBase):
    """The Chebyeshev coefficients for each expansion."""

    __tablename__ = 'cheb'

    p_id = Column(
        Integer,
        ForeignKey('pVacc.id',
                   onupdate='CASCADE',
                   ondelete='RESTRICT'),
        primary_key=True,
        doc='The expansion this coefficient corresponds to.'
    )
    chebOrder = Column(
        Integer,
        primary_key=True,
        doc='The position in the overall expansion the coefficient should be placed.'
    )
    value = Column(
        Float,
        primary_key=False,
        doc='The value of the coefficient.'
    )
    timestamp = Column(
        TIMESTAMP,
        nullable=False,
        doc='When was this record last changed.'
    )
	
    pmsvaccs = relationship('PmsVAcc', back_populates='chebCoeffs')

#pylint: enable=too-few-public-methods
