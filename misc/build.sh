#!/bin/sh
# shellcheck disable=SC1090,SC1091,SC2002,SC2086,SC2097,SC2098,SC2016

embed_document() {
  while read -r line; do
    outdir=.DAJIN_temp/"$(dirname "$line")"
    mkdir -p "$outdir"
    echo 'cat <<'\'EMBED\'' >.DAJIN_temp/'"$line" >/tmp/tmp_header
    echo "EMBED" >/tmp/tmp_footer
    grep ^ "$line" |
      cat /tmp/tmp_header - /tmp/tmp_footer
  done
}

echo '#!/bin/sh' >DAJIN

cat document/version.md |
  sed 's/```/==================================/' |
  sed "s/^/# /" >>DAJIN

# # Make directories
# cat <<EOF >>DAJIN
# find ./ -type d |
#   grep -e "./library" -e "./script" -e "./document" -e "./util" |
#   grep -v -e "past" -e ".DAJIN_temp" |
#   sort -u |
#   sed "s|^./|.DAJIN_temp/|" |
#   xargs mkdir -p
# EOF

# Embed
find document/*.md | embed_document | grep -v '```' >>DAJIN
find util/*.html | embed_document >>DAJIN
find util/*.csv | embed_document >>DAJIN

find library/* -type f | embed_document >>DAJIN
find library/ -name "*.sh" | sed "s|^|. .DAJIN_temp/|" >>DAJIN

echo 'ARGS="$*" && export ARGS' >>DAJIN
find script/ -type f |
  grep -v -e "past" |
  sort |
  embed_document >>DAJIN
find script/ -type f |
  grep -v -e "past" |
  sort |
  sed "s|^|. .DAJIN_temp/|" >>DAJIN

# Execute
chmod +x DAJIN
./DAJIN -h

# Clean up
rm /tmp/tmp_*
