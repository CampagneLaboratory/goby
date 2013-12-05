package edu.cornell.med.icb.goby.alignments;

import it.unimi.dsi.fastutil.objects.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 *     @author Eric Minwei Liu
 *     Date: 11/27/13
 *     Time: 5:56 PM
 */
public class ContigHelper {

    // 0__len__721:203-203	-	chr2	9630375  9630375
    // 1__len__681:107-107	+	chr10	92672164  92672164

    // aa  aa

    private static Object2ObjectMap<String, Boolean> contigMappingStrandHash = new Object2ObjectAVLTreeMap<String, Boolean>();
    private static Object2ObjectMap<String, String> contigMappingChrHash = new Object2ObjectAVLTreeMap<String, String>();
    private static Object2ObjectMap<String, Integer> contigMappingStartHash = new Object2ObjectAVLTreeMap<String, Integer>();
    private static Object2ObjectMap<String, Integer> contigMappingEndHash = new Object2ObjectAVLTreeMap<String, Integer>();


    public ContigHelper(File alnTable) {

        BufferedReader br = null;

        try
        {
            String sCurrentLine;
            br = new BufferedReader(new FileReader(alnTable));

            while ((sCurrentLine = br.readLine()) != null)
            {

                //System.out.println(sCurrentLine);

                String[] splitArray = sCurrentLine.split("\\s+");
                String contigID = splitArray[0];
                String contigMappingStrand = splitArray[1];
                boolean contigMappingStrandBool = contigMappingStrand.equals("+");
                String contigMappingChr = splitArray[2];
                int contigMappingStart = Integer.parseInt(splitArray[3]);
                int contigMappingEnd = Integer.parseInt(splitArray[4]);

                //System.out.printf("%s ## %s\n",contigID,contigMappingStrand);
                //System.out.printf("%s##\n",contigMappingChr);
                //System.out.printf("%d##%d\n",contigMappingStart,contigMappingEnd);

                contigMappingStrandHash.put(contigID,contigMappingStrandBool);
                contigMappingChrHash.put(contigID,contigMappingChr);
                contigMappingStartHash.put(contigID,contigMappingStart);
                contigMappingEndHash.put(contigID,contigMappingEnd);

                //System.out.printf("%s ## %s\n",contigID,contigMappingStrandHash.get(contigID));


            }

        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                if (br != null) br.close();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }




    }


   public String translateContigID(String name) {
          String contigGenomeChr;

          if (contigMappingChrHash.get(name) != null )
          {
              contigGenomeChr = contigMappingChrHash.get(name);
              return (contigGenomeChr);
          } else {
              return (name);


          }
   }

   public int translateContigPosition(String name, int position) {
          int contigGenomePos = 0;

          int mapped_var_start = 0;
          //int mapped_var_end = 0;

     if (contigMappingChrHash.get(name) != null ) {

       if( contigMappingStrandHash.get(name) )
       {

           mapped_var_start = position - 1 + contigMappingStartHash.get(name);
           //mapped_var_end = position -1 + contigMappingStartHash.get(name);
           position = mapped_var_start;
       } else {

           mapped_var_start = contigMappingEndHash.get(name) - position + 1;
           //mapped_var_end = contigMappingEndHash.get(name) - position + 1;
           position = mapped_var_start;
       }

       contigGenomePos = position;
       return(contigGenomePos);

     }else { return (position); }


   }


}
