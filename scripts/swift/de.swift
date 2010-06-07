
type alignment {
  string basename;
}
type tsv;
type textfile;

app (textfile t) version_goby() {
   goby "1g" "version" stdout=@filename(t);
}

app (tsv stats) alignment_to_annotation_counts(
string groupId1,
string group1_basenames,
string groupId2,
string group2_basenames,
string spaceSeparatedBasenames,
string annotationFile,
string useWeights,
string adjustGCBias) {


   goby "3g" "alignment-to-annotation-counts" spaceSeparatedBasenames
   "--groups " @strcat(groupId1,"=",group1_basenames,
   "/",
   groupId2,"=",group2_basenames )
   "--compare " @strcat(groupId1,"/",groupId2)
   "--annotation" annotationFile
   "--use-weights" useWeights
   "--adjust-gc-bias" adjustGCBias
   "--parallel --normalization-methods aligned-count"
   "--include-annotation-types gene --write-annotation-counts false"
   "--eval group-averages"
   "--stats" @stats ;

}



(tsv t) call_de(string groupId1,
                     string group1_basenames,
                     string groupId2,
                     string group2_basenames,
                     string annotationFile,
                     string useWeights,
                     string adjustGCBias) {

  string spaceSeparatedBasename1 = @regexp(group1_basenames, ","," ");
  string spaceSeparatedBasename2 = @regexp(group2_basenames, ","," ");

  string allBasenames = @strcat(spaceSeparatedBasename1," ",spaceSeparatedBasename2);
  

  t=  alignment_to_annotation_counts(groupId1="A", group1_basenames,
                                   groupId2="B", group2_basenames=group2_basenames,
                                   spaceSeparatedBasenames=allBasenames,
                                   annotationFile,
                                   useWeights, adjustGCBias);

}



string group1_basenames="basename1,basename2";
string group2_basenames="basename3,basename4";
string groupId1="A";
string groupId2="B";
string useWeights="true";
string adjustGCBias="false";

tsv stats <"output.tsv">;
// textfile out <"output.txt">;

stats =  call_de(groupId1="A",
                       group1_basenames,
                       groupId2="B",
                       group2_basenames,
                       "/Users/fac2003/IdeaProjects/goby/data/biomart_human_exon_esmbl57genes_NCBI_GRCh37.txt",
                       "true",
                       "false");