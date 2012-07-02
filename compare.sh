rm -f times.tsv
for symbolModeling in ORDER_ONE PLUS ORDER_ZERO
    do
        for filename in /data/bam-comparison/EJOYQAZ-hybrid-domain.entries /data/bam-comparison/XAAOBVT-hybrid-domain.entries
        do
            input=`basename ${filename}`
            echo Processing ${filename}
            time -p java -Xmx3g -jar goby.jar --mode ca ${filename} -o discard -x MessageChunksWriter:codec=hybrid-1 \
            -x MessageChunksWriter:chunk-size=100000 -x AlignmentCollectionHandler:debug-level=1 \
            -x AlignmentCollectionHandler:stats-filename=stats.tsv  \
            -x AlignmentCollectionHandler:basename=${input}-${symbolModeling}  \
            -x AlignmentCollectionHandler:symbol-modeling=${symbolModeling} --max-entries  100002 >& /dev/null ;

            echo ${input}-${symbolModeling}

        done
done

#cat times.tsv