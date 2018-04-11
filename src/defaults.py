import ConfigParser, os, sys

def get_option(name):
    hostname = os.popen("hostname -s").read().rstrip()
    if hostname[:7]=="jupiter": hostname = "jupiter"

    config = ConfigParser.ConfigParser()
    config.readfp(open('../defaults.cfg'))
    return config.get(hostname, name)

if __name__ == "__main__":
    print get_option(sys.argv[1])
