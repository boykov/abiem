import ConfigParser, os
from setuptools import setup

config = ConfigParser.ConfigParser()
config.readfp(open('defaults.cfg'))

setup(name='bie',
      version = 0.1
    )

config.get('config', 'petsc_dir')

