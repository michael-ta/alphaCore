import sys
import numpy as np


def calculate_time_diff(txs, mode="mean"):
    if len(txs) <= 1:
        return 0

    diff = []
    sorted_txs = sorted(txs)
    for i in range(len(txs) - 1):
        diff.append(sorted_txs[i + 1] - sorted_txs[i])

    if mode == "mean":
        return 1/(1 + np.mean(np.array(diff)))
    else:
        return 1/(1 + np.median(np.array(diff)))

if __name__ == '__main__':
    nodes_tx_out = dict()
    nodes_tx_in = dict()
    with open(sys.argv[1], 'r') as f:
        for l in iter(f.readline, ''):
            data = l.strip().split(' ')
            # ignore hacking attempts
            if int(data[3]) > 10**30:
                continue
            # get time of transaction out
            try:
                nodes_tx_out[data[0]]
                nodes_tx_out[data[0]].append(int(data[2]))
            except KeyError:
                nodes_tx_out[data[0]] = [int(data[2])]
            try:
                nodes_tx_out[data[1]]
            except KeyError:
                nodes_tx_out[data[1]] = []

           # get time of transaction in
            try:
                nodes_tx_in[data[1]]
                nodes_tx_in[data[1]].append(int(data[2]))
            except KeyError:
                nodes_tx_in[data[1]] = [int(data[2])]
            try:
                nodes_tx_in[data[0]]
            except KeyError:
                nodes_tx_in[data[0]] = []


    print ("node mean_tx_out median_tx_out mean_tx_in " + \
           "mean_tx_int mean_tx_all median_tx_all mean_tx_sep_all median_tx_sep_all")
    with open(sys.argv[2], 'r') as f:
        for l in iter(f.readline, ''):
            data = l.strip().split('\t')
            mean_time_diff = calculate_time_diff(nodes_tx_out[data[0]]) * float(len(nodes_tx_out[data[0]])) + \
                             calculate_time_diff(nodes_tx_in[data[0]]) * float(len(nodes_tx_in[data[0]]))
            median_time_diff = [calculate_time_diff(nodes_tx_out[data[0]], mode="median"), 
                                calculate_time_diff(nodes_tx_in[data[0]], mode="mdeian")]

            print ("{} {} {} {} {} {} {} {} {}".format(data[1], 
                                     calculate_time_diff(nodes_tx_out[data[0]]),
                                     calculate_time_diff(nodes_tx_out[data[0]], mode="median"),
                                     calculate_time_diff(nodes_tx_in[data[0]]),
                                     calculate_time_diff(nodes_tx_in[data[0]], mode="median"),
                                     calculate_time_diff(nodes_tx_out[data[0]] + nodes_tx_in[data[0]]),
                                     calculate_time_diff(nodes_tx_out[data[0]] + nodes_tx_in[data[0]],
                                                         mode="median"),
                                     mean_time_diff / len(nodes_tx_out[data[0]] + nodes_tx_in[data[0]]),
                                     np.median(np.array(median_time_diff)) ))

