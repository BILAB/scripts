概要
====

- PaCS-MDの xcrypt & Gromacs-2016-x による実装です。
- [PaCS-MD]() と [xcrypt](https://bitbucket.org/tasuku/xcrypt) についての詳細はオリジナルをご覧ください。


ファイルの説明
==============
- `pacs_run.xcr`
    - `xcrypt pacs_run.xcr` で実行するxcrypt スクリプト
- `grompp_bush.sh`
    - `short_run_bush.sh` のための`tpr` ファイルを出力するスクリプト
- `runset.mdp`
    - short MD の設定ファイル; `grompp_bush.sh` で読み込まれる
- `short_run_bush.sh`
    - short MD の実行スクリプト
- `transdist`
    - 各サイクルのベストスコアを抽出,print

- `backtrace_traj.pl`
	- トラジェクトリを結合する準備として，最終サイクルから順に結合するトラジェクトリをトレースする。
- `trajcat.pl`
	- `backtrace_traj.pl` の結果をもとに，トラジェクトリを結合する。

使い方
======

環境
----

- Reedbush-h (GPU計算ノード) で動作を確認しています。
- Reedbush にはxcryptがもともとインストールされていますが、他の計算機でもインストールすることで可能らしいです
- Gromacs-2016-x


ディレクトリ構成(想定)
----------------

- 実行前

```
└── pacsをしたいディレクトリ
    ├── pacs_run.xcr
    ├── grompp_bush.sh
    ├── runset.mdp
    ├── short_run_bush.sh
    └── transdist
```
- 実行後

```
    ├── cyc0
    ├── cyc1
    ├── ... (2 - 100)
    ├── cyc100
    ├── pacs_run.xcr
    ├── dist_trans.dat
    ├── grompp_bush.sh
    ├── inv_watch
    ├── runset.mdp
    ├── short_run_bush.sh
    └── transdist
```

実行手順
--------

- `grompp_bush.sh` のトポロジーファイルやインデックスファイルのパスをお使いの環境に合わせて編集

- ディレクトリに必要なファイルを揃える

```
## reedbushでは
$ module load xcrypt
$ cd your_pacs_dir
## 必ずしも screen である必要はないですが、xcrypt 自体はqsubされているわけではないため、reedbushからログアウトするとスクリプトが停止してしまうことの対策として、nohup なり screen なりを使うのが良いと思います。
## 進行状況が出力されるのでscreenかtmuxがおすすめです。(nohup でどうなるか知りません）
$ screen 
## in screen session
$ xcrypt pacs_run.xcr
job_name <= initialized
job_name <= prepared
job_name <= submitted
...
...
```
- 全てのサイクルを終えたら，トラジェクトリをつなげる作業に入る。
- 好きなディレクトリ（ここでは`cat_all`) を`your_pacs_dir`直下に作成し

```
$ mkdir cat_all
$ cd cat_all/
## やり方は一例です。例えばこんな感じ
## seq 0 120 の部分は適宜変えてください
$ for c in `seq 0 120`;do sort -k3n,3 ../cyc$c/ranking | head -8 | awk '{printf("%s-%d ",$1, $2)} END{print "\n"}' | sed '/^$/d' >> inits.dat ; done
$ less inits.dat
$ perl backtrace_traj.pl inits.dat 120
$ perl backtrace_traj.pl inits.dat  120 > pacs_trace_res.dat
$ less pacs_trace_res.dat
$ cp pacs_trace_res.dat trace1.dat
## PaCS-MDではTHREADSで指定した本数のトラジェクトリが作れるため，１つだけを除いて一旦削除します
## //で区切られているので好きなものを１つ残して，他は消す
$ vim trace1.dat # to delete other than top.
## index.ndx を参照するのでパスを適宜変更してください
## 上手く行けば，full_traj.trr が出力されるはずです
$ perl trajcat.pl trace1.dat
```

xcryptスクリプトの設定
--------------------- 

- パラメータ設定
    - 随時書き足します
    - 計算資源量の部分はややこしいので要注意です
        - 自分が間違えていたら教えてください……
        
- サブルーチン設定
    - サブルーチンは
        1. mesure        : ランキング作成のための測定
        2. merge_xvg     : データの併合
        3. get_next_init : 次サイクルの初期構造の取得 & 生成
    - です。
- xcrypt の使い方 [参考資料](http://www.cc.u-tokyo.ac.jp/support/kosyu/62/shiryou-20160906-5.pdf)


