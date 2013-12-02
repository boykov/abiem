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

    def setup(self):
        self.test_phi = 182

class EllipsoidWITH(Base):
    __tablename__ = 'ellipsoid'

    id = Column(Integer, primary_key=True)
    axe1 = Column(Numeric(36, 16))
    axe2 = Column(Numeric(36, 16))
    axe3 = Column(Numeric(36, 16))
    axes = Column(SqliteArray)

    def setup(self):
        axes = zeros((3))
        axes[:] = [float(self.axe1),float(self.axe2),float(self.axe3)]
        self.axes = axes

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

    def setup(self, axes):
        self.e = ellipsoid(axes, int(D(self.numpoints)))
        self.hval = self.e.get_h()
        self.numnodes = self.e.points.shape[0]
        self.node_coordinates = self.e.points[:,:]
        self.normal_coordinates = self.e.normalvectors[:,:]

def create_session(Base):
    engine = create_engine('sqlite:///ddd.db',
                           echo=False, 
                           encoding='utf-8',  
                           pool_recycle=3600,
                           connect_args={'detect_types':
                                         sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES})
    Base.metadata.create_all(engine)
    return sessionmaker(bind=engine, autoflush=True, autocommit=False)()

class Test(unittest.TestCase):
        
    def test_types(self):
        session = create_session(Base)

        try:
            ell  = session.query(EllipsoidWITH).filter_by(axe1=0.75).first()
            if not ell:
                ell = EllipsoidWITH(axe1=0.75, axe2 = 1.0, axe3 = 0.5)
                ell.setup()
                session.add(ell)
                session.commit()

            pnts = session.query(PointsWITH).filter_by(numpoints=200).first()
            if not pnts:
                pnts = PointsWITH(numpoints=200, surface_id = ell.id)
                pnts.setup(ell.axes)
                session.add(pnts)
                session.commit()

            pnts2 = session.query(PointsWITH).filter_by(numpoints='200').first()
            self.assertEquals(pnts, pnts2)

            phi = session.query(PhiWITH).filter_by(points_id=pnts.id).first()
            if not phi:
                phi = PhiWITH(points_id = pnts.id)
                phi.setup()
                session.add(phi)
                session.commit()

            phi2 = session.query(PhiWITH).filter_by(test_phi='182').first()
            self.assertEquals(phi, phi2)
        finally:
            session.close()
            
        
