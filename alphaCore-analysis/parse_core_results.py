import sys

fn = sys.argv[1]
authorfn = "/home/mta/repos/alphaCore/data/label.citation-statistics.txt"
pubfn = "/home/mta/repos/scripts/alphaCore-analysis/authorPubRec.txt"

authors = []
pubinfo = dict()

with open(authorfn, 'rb') as f:
    for l in iter(f.readline, ''):
        l=l.strip().replace('"', '')
        authors.append(l)

with open(pubfn, 'rb') as f:
    for l in iter(f.readline, ''):
        l=l.strip()
        aname, papers, citedcount, selfcitecount = l.split("\t")
        pubinfo[aname] = {"papers": int(papers),
                          "cited": int(citedcount),
                          "selfcited": int(selfcitecount)}

group=[]
level=None

with open(fn, 'rb') as f:
    # skip first line
    f.readline()    
    for l in iter(f.readline, ''):
        l=l.strip().replace('"', '')
        data = l.split(',')
        # change this to use between k & alpha core
        aidx = int(data[2]) - 1
        if data[-1] == level:
            group.append(authors[aidx])
        else:
            # print out here in alphabetical order
            for a in sorted(group):
                print "\t".join([str(level), a, str(pubinfo[a]["papers"]), 
                                    str(pubinfo[a]["cited"]),
                                    str(pubinfo[a]["selfcited"])])

            group = [authors[aidx]]
            level = data[-1]
    for a in sorted(group):
        print "\t".join([str(level), a, str(pubinfo[a]["papers"]), 
                             str(pubinfo[a]["cited"]),
                             str(pubinfo[a]["selfcited"])])
