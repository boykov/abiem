#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

"""

__author__ = "Evgeny Boykov. <artscan@list.ru>"

import unittest
from ellipsoid import ellipsoid
from test_default_arg_sqlite import *
from numpy import *

set_printoptions(precision=16)
Base = declarative_base()

class PhiWITH(Base):
    __tablename__ = 'phi'

    id = Column(Integer, primary_key=True)
    points_id = Column(Integer, ForeignKey('points.id'))
    test_phi = Column(Numeric(36, 12))

class EllipsoidWITH(Base):
    __tablename__ = 'ellipsoid'

    id = Column(Integer, primary_key=True)
    axe1 = Column(Numeric(36, 16))
    axe2 = Column(Numeric(36, 16))
    axe3 = Column(Numeric(36, 16))
    axes = Column(SqliteArray)

class PointsWITH(Base):
    __tablename__ = 'points'

    id = Column(types.Integer, primary_key=True)
    surface_id = Column(Integer, ForeignKey('ellipsoid.id'))
    numpoints = Column(types.Integer, unique=True)
    timestamp = Column(types.DateTime, default = datetime.datetime.now)
    numnodes = Column(types.Integer, unique=True)
    node_coordinates = Column(SqliteArray)
    normal_coordinates = Column(SqliteArray)
    hval = Column(Numeric(36, 16))
    nstroke_coordinates = Column(SqliteArray)

class IntegWITH(Base):
    __tablename__ = 'integ'
    id = Column(types.Integer, primary_key = True)
    points_id = Column(types.Integer, ForeignKey('points.id'))
    dim_quad = Column(types.Integer)
    intphi_over = Column(SqliteArray)

class SingularWITH(Base):
    __tablename__ = 'singular'
    id = Column(types.Integer, primary_key = True)
    points_id = Column(types.Integer, ForeignKey('points.id'))
    dim_quad = Column(types.Integer)
    k_wave = Column(Numeric(36,16))
    fsingular3 = Column(SqliteArray)

def create_session(Base):
    engine = create_engine('sqlite:///ddd.db',
                           echo=False, 
                           encoding='utf-8',  
                           pool_recycle=3600,
                           connect_args={'detect_types':
                                         sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES})
    Base.metadata.create_all(engine)
    return sessionmaker(bind=engine, autoflush=True, autocommit=False)()
