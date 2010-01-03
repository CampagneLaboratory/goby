package edu.cornell.med.icb.nextgen.datamodel;

import java.io.File;
import java.io.Serializable;

/**
 * A copy of some fileset. The copy may be local to the node or remote.
 * @author Fabien Campagne
 *         Date: Jul 24, 2009
 *         Time: 5:29:35 PM
 */
public interface FileSetCopy extends Serializable {

    /**
     * Returns true when this file set is local to the node.
     * @return
     */
    boolean isLocal();

    /**
     * Delete this file set.
     */
    void delete();

    /**
     * Local file set copies will return the files that correspond to the given file type.
     * Remote file set copies will return null;
     * @param fileType
     * @return
     */
    File[] getFiles(String... fileType);

    /**
     * The size of this file set. This is the sum of file sizes in the file set.
     * @return
     */
    long size();
    long getTimeStamp();

    /**
     * Write the content of the local copy to this resource.
     * @param localCopy
     */
    void write(LocalFileSetCopy localCopy);
}
