import gzip
import sys


if __name__ == '__main__':
    node_id_old = dict()
    records = 0
    outfn = sys.argv[2]
    with open(sys.argv[1], 'rb' ) as f:
        for l in iter(f.readline, ''):
            data = l.strip().split(' ')
            nodes = [str(data[0]), str(data[1])]
            for g in nodes:
                try:
                    node_id_old[g]
                except KeyError:
                    node_id_old[g] = True
            records += 1
            if (records % 100000 == 0):
                print "parsed records: {}".format(records)
    new_id = 1
    with open(outfn, 'wb') as f:
        for gid in sorted(node_id_old.keys()):
            f.write("{}\t{}\n".format(gid, new_id))
            new_id += 1
 
