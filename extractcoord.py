#!/usr/bin/env python3

import re
import subprocess

path = '/Users/YoshitakaM/Desktop/4elr_HCO3.log'
pdb = '/Users/YoshitakaM/Desktop/4elr_HCO3.pdb'

# logファイル内のNAtomsの数とpdbのHETATMレコードまたはATOM  レコードの数が一致するかを確認する
cmd = "grep 'NAtoms=' "+path+" | cut -d ' ' -f 5 | tail -n 1"
Natoms = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).stdout.readlines()
cmd2 = "grep -e '^ATOM  ' -e 'HETATM' "+pdb+" | wc -l"
pdbatoms = subprocess.Popen(cmd2, stdout=subprocess.PIPE, shell=True).stdout.readlines()

if not int(Natoms[0]) == int(pdbatoms[0]):
    print("The number of atoms in the Gaussian log file does not match the atom number in PDB file.")
    exit()

# 座標書き出し処理
readso = 0
readbar = 0
logcount = 1
atomlist = []
with open(path) as f:
    for line in f:
        # Standard orientationをサーチ
        if "Standard orientation" in line:
            readso = 1    #Standard orientation内部にいるフラグ
            readbar = 0   #----カウントをリセット
        if (readso == 1 and ("--------" in line)):
            readbar += 1
            # "Standard orientation"が来ており、かつ"-----"の行が３回来るとリセット
            if (readso == 1 and readbar == 3):
                readso = 0
                readbar = 0
                logcount += 1
                continue
        if ((readso == 1) and not (readbar % 3 == 0)):
            if (re.match('^\s+[0-9]',line)):
                line_dict = line.split()
                atomlist.append(['{:>8.3f}'.format(float(line_dict[3])),
                                 '{:>8.3f}'.format(float(line_dict[4])),
                                 '{:>8.3f}'.format(float(line_dict[5]))])
                # print('{:>8.3f}'.format(float(coordx)))

j = 0

for i in range(int(logcount)-1):
    print("MODEL " + str(int(i) + 1))
    with open(pdb) as p:
        for line in p:
            if re.match('^ATOM  |^HETATM', line):
                line_edit = line.strip()
                print(line_edit[:30] + atomlist[j][0] + atomlist[j][1] + atomlist[j][2] + line_edit[54:])
                j += 1
    print("ENDMDL")
