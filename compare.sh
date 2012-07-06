rm -f times.tsv
for symbolModeling in PLUS ORDER_ONE  ORDER_ZERO
    do
        for filename in /data/bam-comparison/EJOYQAZ-gzip.header /data/bam-comparison/XAAOBVT-gzip.header
        do
            input=`basename ${filename}`
            echo Processing ${filename} with method ${symbolModeling}
            java -Xmx3g -jar goby.jar --mode ca ${filename} -o discard -x MessageChunksWriter:codec=hybrid-1 \
            -x MessageChunksWriter:chunk-size=100000 -x AlignmentCollectionHandler:debug-level=1 \
            -x AlignmentCollectionHandler:stats-filename=stats.tsv  \
            -x AlignmentCollectionHandler:basename=${input}-${symbolModeling}  \
            -x AlignmentCollectionHandler:symbol-modeling=${symbolModeling} --max-entries  100002 >& /dev/null ;

            time -p java -Xmx3g -jar goby.jar --mode cfs discard.entries >& /dev/null
            echo ${input}-${symbolModeling}

        done
done

#cat times.tsv