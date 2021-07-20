#!/bin/sh
# shellcheck disable=SC1090,SC1091,SC2002,SC2086,SC2097,SC2098,SC2016

embed_document() {
  while read -r line; do
    echo 'cat <<'\'EMBED\'' >.DAJIN_temp/'"$line" >/tmp/tmp_header
    echo "EMBED" >/tmp/tmp_footer
    grep ^ "$line" |
      cat /tmp/tmp_header - /tmp/tmp_footer
  done
}

echo '#!/bin/sh' >DAJIN
cat document/VERSION.md |
  sed 's/```/==================================/' |
  sed "s/^/# /" |
  cat >>DAJIN

: 'Initialization' && {
  grep ^ script/00-initialization/initialization.sh
} >>DAJIN

: 'Documentation' && {
  find document/*.md | embed_document | grep -v '```'
} >>DAJIN

: 'util' && {
  find util/*.html | embed_document >>DAJIN
  find util/*.csv | embed_document >>DAJIN
}

: 'script' && {
  find script/ -type f |
    grep -v -e "compile" -e "99-past" -e "00-initialization" |
    sort |
    embed_document >>DAJIN
  echo "chmod -R +x .DAJIN_temp/script/*" >>DAJIN
}

: 'Library' && {
  : 'Generate shell library' && {
    find library/*.sh | embed_document >>DAJIN
  }

  : 'Load shell library' && {
    find library/*.sh |
      while read -r line; do
        echo '. .DAJIN_temp/'"$line"
      done >>DAJIN
  }

  : 'Python library' && {
    find library/*.py | embed_document >>DAJIN
  }
  echo "chmod +x .DAJIN_temp/library/*" >>DAJIN
}

: 'Arguments' && {
  echo 'ARGS="$*" && export ARGS' >>DAJIN
  find script/ -type f |
    grep -v -e "past" -e "init" |
    sort |
    sed "s|^|. .DAJIN_temp/|" >>DAJIN
}

chmod +x DAJIN

./DAJIN -h

rm /tmp/tmp_*
