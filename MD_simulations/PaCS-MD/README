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
    ├── dis_nap-un1_xcrypt.xcr
    ├── dist_trans.dat
    ├── grompp_bush.sh
    ├── inv_watch
    ├── runset.mdp
    ├── short_run_bush.sh
    └── transdist
```

実行手順
--------

- ディレクトリに必要なファイルを揃える

```
# reedbushでは
$ module load xcrypt
$ cd your_pacs_dir
# 必ずしも screen である必要はないですが、xcrypt 自体はqsubされているわけではないため、reedbushからログアウトするとスクリプトが停止してしまうことの対策として、nohup なり screen なりを使うのが良いと思います。
# 進行状況が出力されるのでscreenがおすすめです。(nohup でどうなるか知りません）
$ screen 
# in screen session
$ xcrypt pacs_run.xcr
job_name <= initialized
job_name <= prepared
job_name <= submitted
...
...
```
- tmux からログインしている場合は screen on tmux となりますので prefix キーにお気をつけください
- reedbushにデフォルトでインストールされている screen の prefix は C-a です(おそらく)


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


