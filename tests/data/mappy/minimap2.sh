#!/bin/bash

ref="tests/data/mappy/ref.fa"
que="tests/data/mappy/query.fq"

minimap2 -ax map-ont "$ref" "$que" --cs >tests/data/mappy/align_cs.sam
minimap2 -ax map-ont "$ref" "$que" --cs=long >tests/data/mappy/align_csling.sam
