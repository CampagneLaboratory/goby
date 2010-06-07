type textfile;

type alignment {
  string basename;
}

app (textfile t) version_goby() {
   goby "1g" "version" stdout=@filename(t);
}

app (textfile t) alignment_to_annotation_counts(
string groupId1,
string group1_basenames,
string groupId2,
string group2_basenames,
string spaceSeparatedBasenames,
string annotationFile) {


   goby "3g" "alignment-to-annotation-counts" spaceSeparatedBasenames
   "--groups " @strcat(groupId1,"=",group1_basenames,
   "/",
   groupId2,"=",group2_basenames )
   "--compare " @strcat(groupId1,"/",groupId2)
   "--annotation" annotationFile 
   stdout=@filename(t);

}



(textfile t) call_de(string groupId1,string group1_basenames,
                     string groupId2,string group2_basenames,
                     string annotationFile) {

  string spaceSeparatedBasename1 = @regexp(group1_basenames, ","," ");
  string spaceSeparatedBasename2 = @regexp(group2_basenames, ","," ");

  string allBasenames = @strcat(spaceSeparatedBasename1," ",spaceSeparatedBasename2);

  t=alignment_to_annotation_counts(groupId1="A", group1_basenames,
                                   groupId2="B", group2_basenames=group2_basenames,
                                   spaceSeparatedBasenames=allBasenames,
                                   annotationFile);
}

textfile outfile <"version.txt">;

string group1_basenames="basename1,basename2";
string group2_basenames="basename3,basename4";

    
outfile =  call_de(groupId1="A", group1_basenames,
                  groupId2="B", group2_basenames,
                  "annotations.tsv");