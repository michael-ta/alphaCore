# get the publication record for the authors in the citation dataset

import sys
import collections

# author list
afn = "/home/mta/repos/alphaCore/data/authorList.txt"
# network-citation
ncfn = "/home/mta/repos/alphaCore/data/networkcitation.txt"
# author paper byadj
apfn = "/home/mta/repos/alphaCore-old/data/Citation/authorPaperBiadj.txt"

authors = []
citations = collections.defaultdict(int)
selfcitations = collections.defaultdict(int)
papers = collections.defaultdict(int)

with open(afn, 'rb') as f:
    for l in iter(f.readline, ''):
        l = l.strip().replace('"', '')
        authors.append(l)
        
with open(ncfn, 'rb') as f:
    for l in iter(f.readline, ''):
        l = l.strip()
        data = l.split('\t')
        fidx = int(data[0]) - 1
        aidx = int(data[1]) - 1
        c = int(data[-1])
        if (fidx == aidx):
            selfcitations[authors[aidx]] += c
        else:
            citations[authors[aidx]] += c

with open(apfn, 'rb') as f:
    aidx = 0
    for l in iter(f.readline, ''):
        l = l.strip()
        data = map(int, l.split('\t'))
        papers[authors[aidx]] = sum(data)
        aidx += 1

for k in authors:
    print "\t".join([k, str(papers[k]), 
                     str(citations[k]), str(selfcitations[k])])
        
