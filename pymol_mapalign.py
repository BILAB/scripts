import pymol
import sys
from collections import defaultdict

if len(sys.argv) != 4:
    print("""usage:
    python3 pymol_mapalign.py [pdbfile] [hmmsearch_filt.map] [mapalign.txt]
    
    only for python 3.6 or newer.
    """)
    sys.exit()

pymol.finish_launching()

pdb = sys.argv[1]
model = pdb.split('.')[0]
pymol.cmd.load(pdb)

pymol.stored.resdict = {}
pymol.cmd.iterate("chain A", "pymol.stored.resdict[ int(resi) ] = resn")
resnlist = pymol.stored.resdict.keys()

hmm = sys.argv[2]
mapalign = sys.argv[3]

#defaultdict(int) は val のない key を投げると 0 を返す dict
aligndict = defaultdict(int)

with open(mapalign) as f:
    l = f.readlines()[-1]
    for sec in l.split():
        if ':' in sec:
            secsplit = sec.split(':')
            # 5dir テスト用
            # aligndict[ int(secsplit[1]) ] = int(secsplit[0]) - 6
            
            # pdb file のナンバーと query_filt のナンバーが一致していることを期待しています
            # クエリの都合上ずれるならずれた分だけ定数値を増減させてください
            aligndict[ int(secsplit[1]) ] = int(secsplit[0]) + 1

with open(hmm) as f:
    l = f.readlines()[1:]
    for s in l:
        sec = s.split()
        
        a = aligndict[ int(sec[1]) ]
        b = aligndict[ int(sec[2]) ]
        
        if a * b == 0: continue
        if not a in resnlist or not b in resnlist: continue

        cara = carb = "CB"
        if pymol.stored.resdict[a] == "GLY":
            cara = "CA"
        if pymol.stored.resdict[b] == "GLY":
            carb = "CA"

        pymol.cmd.distance(f"(/{model}//A/{a}/{cara})",(f"(/{model}//A/{b}/{carb})"))

print("\n===finished===")