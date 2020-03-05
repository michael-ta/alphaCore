import sys
import numpy as np

txn = dict()

if __name__ == '__main__':
    with open(sys.argv[1], 'r') as f:
        for l in iter(f.readline, ''):
            data = l.strip().split(' ')
            # skip values that are > 10^30 due to hacking attempts
            if (int(data[2]) >= 10**30):
                continue
            try:
                txn[data[0]]
            except KeyError:
                txn[data[0]] = dict()

            try:
                txn[data[0]][data[1]]
                txn[data[0]][data[1]].append(int(data[2]))
            except KeyError:
                txn[data[0]][data[1]] = [int(data[2])]


    for key1 in txn.keys():
        for key2 in txn[key1].keys():
            print("{} {} {}".format(key1, 
                                    key2, sum(txn[key1][key2]))) 
            
