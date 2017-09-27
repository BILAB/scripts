#!/bin/sh

inp=single2

###
### g16の結果のlogファイルをinputとして、QM領域の部分電荷を自動更新してくれるスクリプト
### 先にg16rootとAmberToolsの環境変数を設定する。

export g16root=/work/gh43/g67002/apps
. $g16root/g16/bsd/g16.profile
export GAUSS_SCRDIR=${CACHE_TMPDIR}

# taopackage
export PATH="/home/moriwaki/apps/taopackage/bin:$PATH"

# AMBERHOMEの確認。PATHが通ってなければ異常終了する。oniomrespはAMBERのプログラムを一部利用するため。
export AMBERHOME="/home/moriwaki/apps/amber16"
test -f ${AMBERHOME}/amber.sh && source ${AMBERHOME}/amber.sh
test -f ${AMBERHOME}/amber.sh || { echo "Environmental variable AMBERHOME is not set. Exit." ; exit 1 ; }

# $PBS_O_WORKDIRが定義されていればそのディレクトリに移動する（qsub時用）。
# 定義されていなければcurrent directoryで実行する
# PJM_O_WORKDIRはOakforest-packs用
# test ${PJM_O_WORKDIR} && cd ${PJM_O_WORKDIR}
# PBS_O_WORKDIRはReedbushやnezu用
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR

## main process

# backup
cp -rp ${inp}.log ${inp}.log.bak

# Gaussian 09, Gaussian 16の仕様変更で、logファイル内の'Input orientation'をStandard orientationに変更してやらないとoniomrespが認識してくれなかった気がする（要検証）
sed -i -e "s/Input orientation/Standard orientation/g" ${inp}.log

# QM領域とそのキャップした原子をGaussianのONIOMのログファイルから抽出し，Gaussianのインプットファイルに保存する
# oniomresp -m1 -g ${inp}.log -o ${inp}_4ESP.gjf を実行すると同時に
# capping atomの数字を変数capnumに記憶させる。
capnum=`oniomresp -m1 -g ${inp}.log -o ${inp}_4ESP.gjf | grep "capping atoms added" | sed -e "s/capping atoms added.//g" | sed -e "s/There are//g"`

# sedコマンドを使って.gjf内の微妙な設定変更を行うと便利
sed -i -e "s/.gjf//g" ${inp}_4ESP.gjf
sed -i -e "s/2000MB/40GB/g" ${inp}_4ESP.gjf
sed -i -e "s#nprocshared=2#nprocshared=64#g" ${inp}_4ESP.gjf

mv ${inp}.log.bak ${inp}.log

# 上の操作で${inp}_4ESP.gjfにはQM領域＋キャッピング原子のみが含まれ、
# ルートセクションにはシングルポイント計算の命令が書かれてあるはず
# これをg16で計算し、出力させる
g16 < ${inp}_4ESP.gjf > ${inp}_4ESP.log ||  { echo "g16 failed." ; exit 1 ; }

# g16の計算が終わったら、続けてoniomresp -m2を実行する
# キャッピング原子の原子数の指定を必ず行う。
oniomresp -m2 -\# ${capnum} -g ${inp}_4ESP.log -o ${inp}_4ESP_RESP.in

espgen -i ${inp}_4ESP.log -o ${inp}_4ESP.esp

resp -i ${inp}_4ESP_RESP.in -o ${inp}_4ESP_RESP.out -p ${inp}_4ESP_RESP.pch -t ${inp}_4ESP_RESP.qout -e ${inp}_4ESP.esp

oniomlog -oi -t ${inp}.gjf -fo ${inp}_oldchg.gjf -i ${inp}.log
#numm2=`expr $num2 + 1`

oniomresp -m3 -g ${inp}_oldchg.gjf -qin ${inp}_4ESP_RESP.qout -o ${inp}new.gjf -c ${inp}_chgcomp.txt

rm *.chk

