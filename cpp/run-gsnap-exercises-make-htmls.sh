#!/bin/bash

rm deleteme-*

make

src/GsnapParseTest1

java -jar ~gobyweb/reads-for-paper/goby.jar -m att \
  -f html deleteme-seq.entries  > deleteme-seq.html

java -jar ~gobyweb/reads-for-paper/goby.jar -m sort \
  deleteme-seq.entries -o deleteme-seq-sorted

java -jar ~gobyweb/reads-for-paper/goby.jar -m att \
  -f html deleteme-seq-sorted.entries  > deleteme-seq-sorted.html

java -jar ~gobyweb/reads-for-paper/goby.jar -m ca \
  deleteme-seq-sorted.entries -o deleteme-seq-hybrid \
 -x MessageChunksWriter:codec=hybrid-1 \
 -x AlignmentCollectionHandler:enable-domain-optimizations=false \
 -x MessageChunksWriter:template-compression=true \
 -x AlignmentCollectionHandler:ignore-read-origin=true \
 -x AlignmentWriterImpl:permutate-query-indices=true

java -jar ~gobyweb/reads-for-paper/goby.jar -m att \
  -f html deleteme-seq-hybrid.entries  > deleteme-seq-hybrid.html