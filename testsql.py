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

class GaussWITH(Base):
    __tablename__ = 'gauss'
    id = Column(types.Integer, primary_key = True)
    integ_id = Column(types.Integer, ForeignKey('integ.id'))
    k_wave = Column(Numeric(36,16))
    gauss1 = Column(SqliteArray)
    gauss3 = Column(SqliteArray)

class AhmedWITH(Base):
    __tablename__ = 'ahmed'
    id = Column(types.Integer, primary_key = True)
    integ_id = Column(types.Integer, ForeignKey('integ.id'))
    singular_id = Column(types.Integer, ForeignKey('singular.id'))
    name_vectorb = Column(types.String(80))
    name_matrixa = Column(types.String(80))
    eps_matgen = Column(Numeric(36, 16))
    eps_aggl = Column(Numeric(36, 16))
    eps_gmres = Column(Numeric(36, 16))
    eta = Column(Numeric(36, 16))
    bmin  = Column(types.Integer)
    rankmax = Column(types.Integer)
    q_ahmed = Column(SqliteArray)

class TestWITH(Base):
    __tablename__ = 'testBIE'
    id = Column(types.Integer, primary_key = True)
    name_matrixa = Column(types.String(80))
    name_approximateu = Column(types.String(80))
    k_wave = Column(Numeric(36,16))
    numpoints = Column(types.Integer)
    slae_tol = Column(Numeric(36, 12))
    slae_places = Column(types.Integer)
    orderquad = Column(types.Integer)
    integ_places = Column(types.Integer)
    timestamp = Column(types.DateTime, default = datetime.datetime.now)


def create_session(Base):
    engine = create_engine('sqlite:///ddd.db',
                           echo=False, 
                           encoding='utf-8',  
                           pool_recycle=3600,
                           connect_args={'detect_types':
                                         sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES})
    Base.metadata.create_all(engine)
    return sessionmaker(bind=engine, autoflush=True, autocommit=False)()
