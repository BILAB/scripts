#!/usr/bin/env python3

import sys
from collections import defaultdict

if len(sys.argv) != 3:
    print("""usage:
    python3 pymol_mapalign.py [hmmsearch_filt.map] [mapalign.txt]
    """)
    sys.exit()

hmm = sys.argv[1]
mapalign = sys.argv[2]

aligndict = defaultdict(int)

with open(mapalign) as f:
    l = f.readlines()[-1]
    for sec in l.split():
        if ':' in sec:
            secsplit = sec.split(':')

            aligndict[ int(secsplit[1]) ] = int(secsplit[0])

with open(hmm) as f:
    l = f.readlines()[1:]
    for s in l:
        sec = s.split()

        sec[1] = int(aligndict[ int(sec[1]) ]) + 1
        sec[2] = int(aligndict[ int(sec[2]) ]) + 1

        if sec[1] * sec[2] == 0:
            sys.stderr.write('No associate given by mapalign against {}'.format(s))
            #sys.exit(1)
            continue
        print(" ".join(map(str,sec)))
