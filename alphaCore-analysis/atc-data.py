# this script is used to form the ATC data from SEES lab and open flights
# for reading into alphaCore

import sys

fn_cities = "/home/mta/repos/alphaCore/data/ATC/na-cities.dat"
fn_airports = "/home/mta/repos/alphaCore/data/ATC/USairport_2010_codes.txt"
fn_data = "/home/mta/repos/alphaCore/data/ATC/US_airport_2010.txt"

# read in domestic airports
UScities = dict()
airports = [None]

with open(fn_cities, 'rb' ) as f:
    for l in iter(f.readline, ''):
        l = l.strip()
        data = l.split('|')
        UScities[data[0]] = data[1]

with open(fn_airports, 'rb') as f:
    for l in iter(f.readline, ''):
        l = l.strip()
        data = l.split(" ")
        airports.append(data[1].strip('"'))

with open(fn_data, 'rb' ) as f:
    for l in iter(f.readline, ''):
        l = l.strip()
        data = l.split(' ')
        acodef = airports[int(data[0])]
        acodet = airports[int(data[1])]

        try:
            UScities[acodef]
            UScities[acodet]
        except KeyError:
            continue

        print l

