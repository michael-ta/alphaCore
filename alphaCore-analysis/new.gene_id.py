import gzip
import sys


if __name__ == '__main__':
    gene_id_old = dict()
    records = 0
    outfn = sys.argv[2]
    with gzip.open(sys.argv[1], 'rb' ) as f:
        for l in iter(f.readline, ''):
            data = l.strip().split('\t')
            genes = [int(data[0]), int(data[1])]
            for g in genes:
                try:
                    gene_id_old[g]
                except KeyError:
                    gene_id_old[g] = True
            records += 1
            if (records % 100000 == 0):
                print "parsed records: {}".format(records)
    new_id = 1
    with open(outfn, 'wb') as f:
        for gid in sorted(gene_id_old.keys()):
            f.write("{}\t{}\n".format(gid, new_id))
            new_id += 1
 
