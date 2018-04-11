import ConfigParser, os, sys

hostname = os.popen("hostname -s").read().rstrip()

config = ConfigParser.ConfigParser()
config.readfp(open('../defaults.cfg'))

print config.get(hostname, sys.argv[1])

