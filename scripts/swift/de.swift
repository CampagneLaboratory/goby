
type alignment {
  string basename;
}
type tsv;
type textfile;

app (textfile t) version_goby() {
   goby "1g" "version" stdout=@filename(t);
}

(string result) tr(string text, string from, string to) {
    app {
        tr from to stdin=text stdout=result;  
    }
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



(tsv t) call_de(     string groupId1,
                     string group1_basenames,
                     string groupId2,
                     string group2_basenames,
                     string annotationFile,
                     string useWeights,
                     string adjustGCBias) {

  // string spaceSeparatedBasename1 = @regexp(group1_basenames, "[,]"," ");
  string spaceSeparatedBasename1 = tr(group1_basenames, "[,]","[ ]");
  trace(spaceSeparatedBasename1);
  string spaceSeparatedBasename2 = tr(group2_basenames, "[,]","[ ]");
  trace(spaceSeparatedBasename2);
  string allBasenames = @strcat(spaceSeparatedBasename1," ",spaceSeparatedBasename2);

  string statsFilename=@strcat(groupId1 ,"-",groupId2 ,"-",useWeights ,"-",adjustGCBias ,".tsv");
  tsv stats <single_file_mapper;file=statsFilename>;

  trace(tr("AAABBBCCC","[B]","[Z]"));
 
 /* t=  alignment_to_annotation_counts(groupId1="A", group1_basenames,
                                   groupId2="B", group2_basenames=group2_basenames,
                                   spaceSeparatedBasenames=allBasenames,
                                   annotationFile,
                                   useWeights, adjustGCBias);
   */
}


