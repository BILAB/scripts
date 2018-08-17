'''
(c) 2018 satot & Yoshitaka Moriwaki, BILAB

License: BSD-2-Clause
'''

from __future__ import absolute_import
from __future__ import with_statement

from pymol import cmd, stored, CmdException

def contactmap(selection='', map='hmmsearch_filt.map', align='mapalignout.txt', offset=0):
    '''
DESCRIPTION

    visualize the results from map_align algorithm in PyMOL.
    https://github.com/sokrypton/map_align

ARGUMENTS

    selection = string: atom selection {default: all}

    map = string: name of contact map file {default: 'hmmsearch_filt.map'}

    align = string: name of map_align result txt file {default: 'mapalignout.txt'}

    offset = int: slide the sequence number to fit its original number manually {default: 0}

EXAMPLE

    contactmap 5dir, yourcontact.map, mapalignout.txt, offset=-7

    '''
    from collections import defaultdict

    if selection == '':
        print("ERROR: choose the model to display the contact map onto")
        raise CmdException
    model = selection

    stored.resdict = {}
    cmd.iterate(model, "stored.resdict[ int(resi) ] = resn")
    resnlist = stored.resdict.keys()

    aligndict = defaultdict(int)

    with open(align) as f:
        l = f.readlines()[-1]
        for sec in l.split():
            if ':' in sec:
                secsplit = sec.split(':')
                aligndict[ int(secsplit[1]) ] = int(secsplit[0]) + 1 + int(offset)

    with open(map) as f:
        l = f.readlines()[1:]
        for s in l:
            sec = s.split()

            a = aligndict[ int(sec[1]) ]
            b = aligndict[ int(sec[2]) ]

            if a * b == 0: continue
            if not a in resnlist or not b in resnlist: continue

            cara = carb = "CB"
            if stored.resdict[a] == "GLY":
                cara = "CA"
            if stored.resdict[b] == "GLY":
                carb = "CA"

            cmd.distance(f"{a}_{b}",f"(/{model}//A/{a}/{cara})",(f"(/{model}//A/{b}/{carb})"))

    cmd.set('dash_color', 'cyan')
    cmd.set('dash_radius', '0.25')
    print("===finished===")

cmd.extend('contactmap', contactmap)

# tab-completion of arguments
cmd.auto_arg[0].update({
    'contactmap'       : cmd.auto_arg[0]['zoom'],
})
