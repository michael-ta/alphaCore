import sys
import gzip

if __name__ == '__main__':
    # read gene ids
    node_id_map = dict()
    with open(sys.argv[1], 'rb') as f:
        for l in iter(f.readline, ''):
            data = l.strip().split('\t')
            node_id_map[data[0]] = data[1]

    with open(sys.argv[2], 'rb') as f:
        for l in iter(f.readline, ''):
            data = l.strip().split(' ')
            print "{} {} {}".format(
                node_id_map[data[0]],
                node_id_map[data[1]],
                data[2])

