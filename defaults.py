import ConfigParser, os, sys

config = ConfigParser.ConfigParser()
config.readfp(open('defaults.cfg'))

print config.get('config', sys.argv[1])

