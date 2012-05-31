
type alignment {
  string basename;
}
type tsv;
type textfile;

app (textfile t) version_goby() {
   goby "1g" "version" stdout=@filename(t);
}

(textfile result) tr(textfile text, string from, string to) {
    app {
        tr from to stdin=@filename(text) stdout=@filename(result);
    }
}

(string result) replace(string text, string from, string to) {
 textfile f <concurrent_mapper;prefix="string-input", suffix=".txt">;
 f=writeData(text);
 textfile replacedText <concurrent_mapper;prefix="replace-string", suffix=".txt">;
 replacedText=tr(f,from,to);
 result=readData(replacedText);
}

app (tsv stats) alignment_to_annotation_counts(
string groupId1,
string group1_basenames,
string groupId2,
string group2_basenames,
string fullPathEntries[],
string annotationFile,
string useWeights,
string adjustGCBias) {



  goby "3g" "alignment-to-annotation-counts" fullPathEntries
   "--groups" @strcat(groupId1,"=",group1_basenames,"/",groupId2,"=",group2_basenames )
   "--compare" @strcat(groupId1,"/",groupId2)
   "--annotation" annotationFile
   "--use-weights" useWeights
   "--adjust-gc-bias" adjustGCBias
   "--normalization-methods" "aligned-count"
   "--include-annotation-types" "gene" "--write-annotation-counts" "false"
   "--eval" "group-averages"
   "--stats" @stats;

}



(tsv t) call_de(     string groupId1,
                     string group1_basenames,
                     string groupId2,
                     string group2_basenames,
                     string annotationFile,
                     string useWeights,
                     string adjustGCBias) {

  // string spaceSeparatedBasename1 = @regexp(group1_basenames, "[,]"," ");
  string spaceSeparatedBasename1 = replace(group1_basenames, "[,]","[ ]");
 // trace(spaceSeparatedBasename1);
  string spaceSeparatedBasename2 = replace(group2_basenames, "[,]","[ ]");
 // trace(spaceSeparatedBasename2);
  string allBasenames = @strcat(spaceSeparatedBasename1," ",spaceSeparatedBasename2);

  string statsFilename=@strcat(groupId1 ,"-",groupId2 ,"-",useWeights ,"-",adjustGCBias);
  tsv stats  <concurrent_mapper;prefix=statsFilename, suffix=".tsv">;
  trace(@filename(stats));
  // trace(replace("AAABBBCCC","[B]","[Z]"));
  string fullPathEntries[];
  string spaceSeparatedBasenames=allBasenames;
  string entries[]=@strsplit(spaceSeparatedBasenames,"\\s");

  string currentDirectory="/data/helicos-ILM-SOLID/";

  foreach entry,i in entries {
    fullPathEntries[i]=@strcat(currentDirectory,entry);
  }

  stats=  alignment_to_annotation_counts(groupId1=groupId1, group1_basenames,
                                   groupId2=groupId2, group2_basenames=group2_basenames,
                                   fullPathEntries=fullPathEntries,
                                   annotationFile,
                                   useWeights, adjustGCBias);
  
}

                               


string helicos_HBR_Basenames= "UMTVLVQ-helicos-brain.entries";
string bullard_HBR_Basenames= "DLTTEJH-Bullard-HBR-SRR037439.entries,DOWTGPI-Bullard-HBR-SRR037444.entries,HHNVSNR-Bullard-HBR-SRR037440.entries,ORVQUWJ-Bullard-HBR-SRR035678.entries,RRTFBOP-Bullard-HBR-SRR037443.entries,WKQCRQC-Bullard-HBR-SRR037442.entries,ZOBXCNB-Bullard-HBR-SRR037441.entries";
string illumina_HBR_Basenames="ALTZLBT-Illumina-brains_2_sequence.entries,DWFWKHJ-Illumina-brains_4_sequence.entries,FZCSFFY-Illumina-brains_8_sequence.entries,LNZIGRF-Illumina-brains_1_sequence.entries,QNDPONW-Illumina-brains_3_sequence.entries,SYXULGD-Illumina-brains_7_sequence.entries,UUBMNUK-Illumina-brains_6_sequence.entries";
string solid_HBR_Basenames=   "LHFEDQE-solid-HBR.entries";


// textfile out <"output.txt">;
foreach adjustBias in ["false","formula2","formula3","formula4"] {

    foreach useWeights in ["false","gc"]      {

        if (!(adjustBias=="true" && useWeights !="gc")) {

            // adjustBias requires gc weights
            tsv stats;
            stats=call_de(groupId1="Bullard-ILM-HBR",
                       bullard_HBR_Basenames,
                       groupId2="Helicos-HBR",
                       helicos_HBR_Basenames,
                       "/Users/fac2003/IdeaProjects/goby/data/biomart_human_exon_esmbl57genes_NCBI_GRCh37.txt",
                       useWeights,
                       adjustBias );

        }
    }
}