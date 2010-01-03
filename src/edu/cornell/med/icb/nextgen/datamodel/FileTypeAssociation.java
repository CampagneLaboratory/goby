package edu.cornell.med.icb.nextgen.datamodel;

import edu.cornell.med.icb.util.LoggingEqualsBuilder;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.log4j.Level;

import java.io.File;
import java.io.Serializable;

/**
 * @author Fabien Campagne
 *         Date: Jul 31, 2009
 *         Time: 4:56:38 PM
 */
class FileTypeAssociation implements Serializable {
    String fileType;
    String filename;

    public FileTypeAssociation(File dataFile, String fileType) {
        this.filename = dataFile.getName();
      //  System.out.format("new file association: %s %s", this.filename, this.fileType);
        this.fileType = fileType;
    }

    /**
     * Empty constructor is required by serialization mechanism.
     */
    public FileTypeAssociation() {
    }


    public File getFile(LocalFileSetCopy localFileSetCopy) {
        final String concatResult = FilenameUtils.concat(localFileSetCopy.location, filename);
        assert concatResult != null : String.format("concat result cannot be null, parent path: %s filename: %s",
                localFileSetCopy.location, filename);
        return new File(concatResult);
    }

    public String toString() {
        return filename + ":" + fileType;
    }

    @Override
    public int hashCode() {
        // you pick a hard-coded, randomly chosen, non-zero, odd number
        // ideally different for each class
        return new HashCodeBuilder(19, 23)
            .append(fileType)
            .append(filename)
            .toHashCode();
    }

    @Override
    public boolean equals(final Object obj) {
        if (obj == null) {
            return false;
        }
        if (obj == this) {
            return true;
        }
        if (obj.getClass() != getClass()) {
            return false;
        }
        final FileTypeAssociation rhs = (FileTypeAssociation) obj;
        return new LoggingEqualsBuilder(this.getClass(), Level.INFO)
            .append("fileType", fileType, rhs.fileType)
            .append("filename", filename, rhs.filename)
            .isEquals();
    }
}
