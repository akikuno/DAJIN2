#!/bin/sh

################################################################################
# コメント行
# パイプを縦に整列
# インデントを2つずつにする
################################################################################

cat "$1"                                                    |
  #== Adjust the length of the comment line ================
  #-- H1 ---------------------------------------------------
  awk '$0 ~ /^###/ && length <= 80 {
    seq=""
    for(i=1;i<=80;i++) seq=seq "#"
    print seq
    next}1'                                                 |
  #-- H2 ---------------------------------------------------
  awk '$0 ~ /^#==/ && length <= 80 {
    seq=""
    for(i=length;i<=78;i++) seq=seq "="
    print $0" "seq
    next}1'                                                 |
  #-- H3 ---------------------------------------------------
  awk '$0 ~ /^#--/ && length <= 80 {
    seq=""
    for(i=length;i<=78;i++) seq=seq "-"
    print $0" "seq
    next}1'                                                 |

  #-- H1 within indent -------------------------------------
  awk '$0 ~ /^  ###/ && length <= 60 {
    seq=""
    for(i=length;i<=60;i++) seq=seq "#"
    print $0 seq
    next}1'                                                 |
  #-- H2 within indent -------------------------------------
  awk '$0 ~ /^  #==/ && length <= 60 {
    seq=""
    for(i=length;i<=60-2;i++) seq=seq "="
    print $0" "seq
    next}1'                                                 |
  #-- H3 within indent -------------------------------------
  awk '$0 ~ /^  #--/ && length <= 60 {
    seq=""
    for(i=length;i<=60-2;i++) seq=seq "-"
    print $0" "seq
    next}1'                                                 |
  # #-- Indent
  # sed "s/^    /  /" |
  # sed "s/^      /    /" |
  # sed "s/^        /      /" |
  #-- Align pipes(|) and escape(\) vertically --------------
  awk 'substr($0, length) == "\|" &&  length <= 60 {
    sub(/\|$/, "")
    seq=""
    for(i=length;i<=60-2;i++) seq=seq " "
    print $0" "seq"|"
    next}1'                                                 |
  awk 'substr($0, length) == "\\" &&  length <= 60 {
    sub(/\\$/, "")
    seq=""
    for(i=length;i<=60-2;i++) seq=seq " "
    print $0" "seq"\\"
    next}1'