import sys

fn = "/home/mta/repos/alphaCore/data/label.airport-US2010.txt"
rfn = sys.argv[1]

acodes = []
counter = 0
with open(fn, 'rb') as f:
    for l in iter(f.readline, ''):
        l = l.strip().strip('"')
        acodes.append(l)

with open(rfn, 'rb') as f:
    f.readline()
    for l in iter(f.readline, ''):
        l = l.strip()
        data = l.split(',')
        airport = data[2].strip('"')
        print "{}\t{}".format(data[3].strip('"'), acodes[int(airport) - 1])
