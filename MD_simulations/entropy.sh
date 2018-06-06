#!/bin/bash


PROGNAME=$(basename $0)
VERSION="0.1"
HELP_MSG="'$PROGNAME -h' to show help"

# ヘルプメッセージ
usage() {
  echo "Requirements:"
  echo "  python 3.6 ( developed in 3.6.5 )"
  echo "  MDAnalysis (run \` pip install --upgrade MDAnalysis \` )"
  echo "  Julia 0.6 ( developed in 0.6.2 )"
  echo "  PyCall (run julia> \` Pkg.add(\"PyCall\") \` in julia repl)"
  echo "  _entropy.jl (in the same directory)"
  echo "-- It will NOT work on Julia 0.7, 1.0 or later --"
  echo ""
  echo "Usage: $PROGNAME [-h|--version] | -f trr -s tpr [-t temp] [-b]"
  echo 
  echo "options:"
  echo "  -h, --help"
  echo "      --version"
  echo "  -f, <tragectory file>          <required>"
  echo "  -s, <topology file>            <required>"
  echo "  -t, <temperature = 300.0> [K]  system temperature  <optional>"
  echo "  -b     when topology file has no bond information. <optional>"
  echo
  exit 1
}

ARG_T="300.0"
# オプション解析
for OPT in "$@"
do
  case "$OPT" in
    # ヘルプメッセージ
    '-h'|'--help' )
      usage
      exit 1
      ;;
    # バージョンメッセージ
    '--version' )
      echo $VERSION
      exit 1
      ;;
    # オプション-a、--long-a
    '-f' )
      FLG_F=1
      # オプションに引数がなかった場合（必須）
      if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
        echo "$PROGNAME:\"$f\" needs input trajectory file" 1>&2
        exit 1
      fi
      ARG_F="$2"
      shift 2
      ;;
    # オプション-b、--long-b
    '-s' )
      FLG_S=1
      # オプションに引数がなかった場合（必須）
      if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
        echo "$PROGNAME:\"$s\" needs input topology file" 1>&2
        exit 1
      fi
      ARG_S="$2"
      shift 2
      ;;
    '-t' )
      FLG_T=1
      # オプションに引数がなかった場合（必須）
      if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
        echo "$PROGNAME:\"$t\" needs temperature" 1>&2
        exit 1
      fi
      ARG_T="$2"
      shift 2
      ;;
    
    # オプション-c
    '-b' )
      # オプション指定のみなのでフラグだけ設定（引数がないタイプ）
      FLG_B=1
      shift 1
      ;;
    '--'|'-' )
      # 「-」か「--」だけ打った時
      shift 1
      param+=( "$@" )
      break
      ;;
    -*)
      echo "$PROGNAME: \"$(echo $1 | sed 's/^-*//')\" option does not exist" 1>&2
      exit 1
      ;;
    *)
      # コマンド引数（オプション以外のパラメータ）
      if [[ ! -z "$1" ]] && [[ ! "$1" =~ ^-+ ]]; then
        param+=( "$1" )
        shift 1
      fi
      ;;
  esac
done

if [ -z $FLG_F ] || [ -z $FLG_S ]; then
  echo "$PROGNAME: input file required" 1>&2
  echo $HELP_MSG 1>&2
  exit 1
fi

if ($FLG_B); then
    julia _entropy.jl $ARG_S $ARG_F $ARG_T 0
else
    julia _entropy.jl $ARG_S $ARG_F $ARG_T 1
fi