package edu.cornell.med.icb.nextgen.datamodel;

import java.io.File;
import java.io.Serializable;

/**
 * A resource that a job needs to operate. A set of files or software that must be installed
 * on the local node before some job can start.
 *
 * @author Fabien Campagne
 *         Date: Jul 21, 2009
 *         Time: 4:30:50 PM
 */
public interface Resource extends Serializable {
    String FILE_LIST = "file-list";

    /**
     * Copy or install the resource on the local node. May throw an exception if more node
     * resources than available are required to install the resource locally.
     */
    void makeLocal();

    /**
     * Remove local copies of the resource.
     */
    void clearLocal();

    /**
     * Returns true if the resource is available locally.
     *
     * @return True or false.
     */
    boolean isLocal();

    /**
     * The files that make up this resource on the node filesystem and are of a certain type.
     * File types are resource dependent. Some file types are logical. For instance basename,
     * returns the basename for a set of files contained in a resource, even though no such
     * file actually exists with this name in the resource.
     *
     * @param fileType The type of file whose filenames are sought.
     * @return
     */
    File[] getLocalFiles(String... fileType);

    /**
     * The size of this resource in bytes.
     *
     * @return resource size in bytes.
     */
    long getSize();

    /**
     * Find the size of the resource without making it local.
     * @return
     */
    long peekAtSize();

    /**
     * Write the content of this local copy to the resource.
     *
     * @param localCopy
     */
    void write(LocalFileSetCopy localCopy);

    /**
     * Return the identifier of this resource. The id is unique among all the resources of
     * a job submission.
     *
     * @return the identifier of this resource.
     */
    String getId();

    long getTimeStamp();

    void synchronizeState(Resource remoteResource);
}
