#!/bin/bash
# Bundles index.html + app.js into one self-contained page.
cd "$(dirname "$0")"
awk '{ if ($0 ~ /<script src="app.js"><\/script>/) { print "<script>"; while ((getline line < "app.js") > 0) print line; print "</script>" } else print }' index.html > ../bomb-calorimetry-lab.html
wc -c ../bomb-calorimetry-lab.html
