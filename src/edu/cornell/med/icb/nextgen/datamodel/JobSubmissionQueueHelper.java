/*
 * Copyright (C) 2009 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 * WEILL MEDICAL COLLEGE OF CORNELL UNIVERSITY MAKES NO REPRESENTATIONS
 * ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY PURPOSE. IT IS PROVIDED
 * "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY. THE WEILL MEDICAL COLLEGE
 * OF CORNELL UNIVERSITY SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY
 * THE USERS OF THIS SOFTWARE.
 */

package edu.cornell.med.icb.nextgen.datamodel;

import com.xerox.amazonws.sqs2.Message;
import com.xerox.amazonws.sqs2.MessageQueue;
import com.xerox.amazonws.sqs2.SQSException;
import com.xerox.amazonws.sqs2.SQSUtils;
import edu.cornell.med.icb.config.ConfigHelper;
import edu.cornell.med.icb.util.GobyPropertyKeys;
import edu.cornell.med.icb.util.GroovyProperties;
import edu.cornell.med.icb.util.ICBFilenameUtils;
import org.apache.commons.codec.binary.Base64;
import org.apache.commons.lang.SerializationUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.jets3t.service.S3Service;
import org.jets3t.service.S3ServiceException;
import org.jets3t.service.impl.rest.httpclient.RestS3Service;
import org.jets3t.service.model.S3Bucket;
import org.jets3t.service.model.S3Object;
import org.jets3t.service.security.AWSCredentials;
import org.jets3t.service.utils.ServiceUtils;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.Serializable;
import java.security.NoSuchAlgorithmException;
import java.util.HashMap;
import java.util.Map;
import java.util.UUID;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * JobSubmission to / from SQS message body. If a message is too large to fit in the SQS
 * queue message size, this will store the message on S3.
 *
 * @author Kevin Dorff
 */
public class JobSubmissionQueueHelper implements Serializable {

    private static final long serialVersionUID = 8889899335790944014L;

    /**
     * Logger.
     */
    private static final Log LOG = LogFactory.getLog(JobSubmissionQueueHelper.class);

    /**
     * Prefix when a message fits completely in queue so we know it is a queue message.
     */
    private static final String INLINE_PREFIX = "INLINE:";

    /**
     * Default SQS message size (including INLINE_PREFIX).
     */
    public static final int SQS_MESSAGE_SIZE_DEFAULT = 8192;

    /**
     * Maximum SQS message size.
     */
    public static final int SQS_MESSAGE_SIZE_MAXIMUM = 8192;

    /**
     * Minimum SQS message size.
     */
    public static final int SQS_MESSAGE_SIZE_MINIMUM = 250;

    /**
     * Regex pattern for matching S3 URIs (bucket name and path).
     */
    private static final Pattern S3_URI_PATTERN = Pattern.compile("s3://(\\w+)/(.+)");

    /**
     * The S3 service we are using to access S3 files.
     */
    private final S3Service s3Service;

    /**
     * The s3Bucket used to store files when a queue message is too long.
     */
    private final S3Bucket s3Bucket;

    /**
     * The queue to delete messages from.
     */
    private final MessageQueue queue;

    /**
     * The default queue data path to store files that are too large to put in the queue.
     */
    private final String s3QueueDataPath;

    /**
     * The job submission id to the message map so we can delete the message when we're
     * done with it.
     */
    private final Map<String, Message> jobSubmissionIdToMessageMap;

    /**
     * The SQS Message size we are using.
     */
    private int sqsMaxMessageSize;
    private final String queueName;

    /**
     * Create the JobSubmissionHelper. This helps to read, write, and delete
     * JobSubmissions to/from an SQS queue. Note that the data here will Base64 encoded
     * so this will queue.setEncoding(false) so as to not double encode/decode the data.
     *
     * @param s3Service       the s3Service to write data to when a message won't fit in the queue
     * @param s3Bucket        the bucket to write data to when a message won't fit in the queue
     * @param s3QueueDataPath the path in the s3 bucket to write data to when a message won't
     *                        fit in the queue
     * @param queue           the queue to read / write messages from / to
     * @throws IOException error finding the specified bucket
     */
    public JobSubmissionQueueHelper(
            final S3Service s3Service, final S3Bucket s3Bucket,
            final String s3QueueDataPath, final MessageQueue queue) throws IOException {
        this.jobSubmissionIdToMessageMap = new HashMap<String, Message>();
        this.s3Service = s3Service;
        this.queue = queue;
        this.queue.setEncoding(false);
        this.s3Bucket = s3Bucket;
        sqsMaxMessageSize = SQS_MESSAGE_SIZE_DEFAULT;
        if (s3QueueDataPath.endsWith("/")) {
            this.s3QueueDataPath = s3QueueDataPath;
        } else {
            this.s3QueueDataPath = s3QueueDataPath + "/";
        }
        this.queueName = queue.getUrl().getPath();
    }

    /**
     * Create the JobSubmissionHelper. This helps to read, write, and delete
     * JobSubmissions to/from an SQS queue. Note that the data here will Base64 encoded
     * so this will queue.setEncoding(false) so as to not double encode/decode the data.
     *
     * @param s3BucketName    the bucket to write data to when a message won't fit in the queue
     * @param s3QueueDataPath the path in the s3 bucket to write data to when a message won't
     *                        fit in the queue
     * @param queueName       the name of the queue to read / write messages from / to
     * @throws IOException error finding the specified bucket
     */
    public JobSubmissionQueueHelper(final String s3BucketName,
                                    final String s3QueueDataPath, final String queueName) throws IOException {
        this.jobSubmissionIdToMessageMap = new HashMap<String, Message>();

        // From the properties file (LocalConfig.groovy)
        final GroovyProperties properties = ConfigHelper.loadConfiguration();

        final String amazonAccessKey = properties.assertGet(GobyPropertyKeys.AWS_ACCESS_KEY);
        final String amazonSecretKey = properties.assertGet(GobyPropertyKeys.AWS_SECRET_ACCESS_KEY);

        // Connect to S3
        final AWSCredentials awsCredentials = new AWSCredentials(amazonAccessKey, amazonSecretKey);

        final S3Bucket s3Bucket;
        try {
            this.s3Service = new RestS3Service(awsCredentials);
            s3Bucket = this.s3Service.getOrCreateBucket(s3BucketName);
        } catch (S3ServiceException e) {
            throw new IOException("Error connecting to S3 and retrieving a bucket", e);
        }
        try {
            this.queue = SQSUtils.connectToQueue(
                    queueName, amazonAccessKey, amazonSecretKey);
        } catch (SQSException e) {
            throw new IOException("Cannot connect to queue " + queueName, e);
        }
        this.queue.setEncoding(false);
        this.s3Bucket = s3Bucket;
        sqsMaxMessageSize = SQS_MESSAGE_SIZE_DEFAULT;
        if (s3QueueDataPath.endsWith("/")) {
            this.s3QueueDataPath = s3QueueDataPath;
        } else {
            this.s3QueueDataPath = s3QueueDataPath + "/";
        }
        this.queueName = queueName;
    }

    /**
     * Get the SQS maximum message size.
     *
     * @return the maximum message size
     */
    public int getSqsMaxMessageSize() {
        return sqsMaxMessageSize;
    }

    /**
     * Set the SQS maximum message size. Any values > 8192 will be set to 8192.
     * Any value < 250 will be set to 250.
     *
     * @param sqsMaxMessageSize the maximum message size to use with SQS
     */
    public void setSqsMaxMessageSize(final int sqsMaxMessageSize) {
        this.sqsMaxMessageSize = sqsMaxMessageSize;
        if (this.sqsMaxMessageSize < SQS_MESSAGE_SIZE_MINIMUM) {
            this.sqsMaxMessageSize = SQS_MESSAGE_SIZE_MINIMUM;
        } else if (this.sqsMaxMessageSize > SQS_MESSAGE_SIZE_MAXIMUM) {
            this.sqsMaxMessageSize = SQS_MESSAGE_SIZE_MAXIMUM;
        }
    }

    /**
     * Submit a JobSubmission to the queue.
     *
     * @param jobSubmission the job submission to submit
     * @throws IOException              error serializing or creating an S3Object
     * @throws S3ServiceException       error writing to S3
     * @throws SQSException             error writing to SQS
     * @throws NoSuchAlgorithmException error creating an S3Object
     */
    public void addJobSubmissionToSQSQueue(final JobSubmission jobSubmission)
            throws IOException, S3ServiceException, SQSException, NoSuchAlgorithmException {
        if (jobSubmission == null) {
            return;
        }
        final ByteArrayOutputStream serialized = new ByteArrayOutputStream();
        SerializationUtils.serialize(jobSubmission, serialized);
        serialized.close();
        final String b64encoded = new String(Base64.encodeBase64(serialized.toByteArray()));
        final String message;
        LOG.info("Completed job details: " + jobSubmission.toString());
        if (b64encoded.length() + INLINE_PREFIX.length() <= sqsMaxMessageSize) {
            // Message fits in the queue
            LOG.info(String.format("Submitting JobSubmission %s 'inline'", jobSubmission.getId()));
            message = "INLINE:" + b64encoded;
        } else {
            final String s3Filename = s3QueueDataPath + UUID.randomUUID().toString();
            final S3Object s3Object = new S3Object(s3Filename, b64encoded);
            s3Service.putObject(s3Bucket, s3Object);
            message = ICBFilenameUtils.concatPathParts("s3://", s3Bucket.getName(), s3Filename);
            LOG.info(String.format("Saved JobSubmission %s to S3 ", jobSubmission.getId()) + message);
        }
        queue.sendMessage(message);
        LOG.info("Job submitted to queue " + queueName);
    }

    /**
     * Retrieve the next JobSubmission from the queue or null if there wasn't one.
     *
     * @return the JobSubmission or null if there wasn't one waiting in the queue
     * @throws IOException        bad message type found in SQS
     * @throws S3ServiceException error reading from S3
     * @throws SQSException       error reading from SQS
     */
    public JobSubmission retrieveJobSubmissionFromQueue()
            throws IOException, S3ServiceException, SQSException {
        final Message message;
        final int fiveMinutes = 5 * 60;
        // message will time out if we do not acknowledge execution within 5 minutes. It will then reappear on the queue.
        message = queue.receiveMessage(fiveMinutes);

        if (message == null) {
            return null;
        }
        final String b64Encoded;
        final String messageBody = message.getMessageBody();
        if (messageBody.startsWith(INLINE_PREFIX)) {
            LOG.info("JobSubmission read 'inline'");
            b64Encoded = messageBody.substring(INLINE_PREFIX.length());
        } else {
            final Matcher matcher = S3_URI_PATTERN.matcher(messageBody);
            if (!matcher.find()) {
                throw new IOException(
                        "Message prefix unrecognized (" + messageBody.substring(0, 10) + ")");
            }
            final String matchedBucketName = matcher.group(1);
            final String filePath = matcher.group(2);
            LOG.info(String.format(
                    "JobSubmission will be read from S3 s3://%s/%s", matchedBucketName, filePath));
            final S3Bucket matchedS3Bucket;
            if (matchedBucketName.equals(s3Bucket.getName())) {
                matchedS3Bucket = s3Bucket;
            } else {
                matchedS3Bucket = s3Service.getOrCreateBucket(matchedBucketName);
            }
            final S3Object s3Object = s3Service.getObject(matchedS3Bucket, filePath);

            b64Encoded = ServiceUtils.readInputStreamToString(
                    s3Object.getDataInputStream(), "UTF-8");
        }
        final Object o;
        try {
            o = SerializationUtils.deserialize(Base64.decodeBase64(b64Encoded.getBytes()));
        } catch (Exception e) {
            LOG.info("Unable to deserialize a jobSubmission from queue message. "
                    + "Skipping this message. Exception message=" + e.getMessage());
            // we don't acknowledge this message, so it will reappear on the queue after the visibility period,
            // hopefully to be picked up by a compatible listener.
            return null;
        }
        if (o instanceof JobSubmission) {
            final JobSubmission js = (JobSubmission) o;
            jobSubmissionIdToMessageMap.put(js.getId(), message);
            LOG.info("JobSubmission id=" + js.getId() + " retrieved");
            return js;
        } else {
            LOG.error("Object from queue was not a job submission: " + o.toString());
            throw new IOException("Object from queue was not a JobSubmission");
        }
    }

    /**
     * This will delete the queue's message AND if the queue message points to a file
     * in S3, this will ALSO delete the file in S3.
     *
     * @param jobSubmission the job submission that we are done with that came from the queue
     * @throws S3ServiceException error deleting from S3
     * @throws SQSException       error removing message from SQS
     */
    public void deleteJobSubmissionFromQueue(final JobSubmission jobSubmission)
            throws S3ServiceException, SQSException {
        if (jobSubmission == null) {
            return;
        }
        final Message message = jobSubmissionIdToMessageMap.get(jobSubmission.getId());
        if (message == null) {
            return;
        }
        LOG.info(String.format("Going to delete the JobSubmission %s %s",
                jobSubmission.getId(), jobSubmission.getStatus()));
        final Matcher matcher = S3_URI_PATTERN.matcher(message.getMessageBody());

        if (matcher.find()) {
            final String matchedBucketName = matcher.group(1);
            final String filePath = matcher.group(2);
            final S3Bucket matchedS3Bucket;
            if (matchedBucketName.equals(s3Bucket.getName())) {
                matchedS3Bucket = s3Bucket;
            } else {
                matchedS3Bucket = s3Service.getOrCreateBucket(matchedBucketName);
            }
            s3Service.deleteObject(matchedS3Bucket, filePath);
            LOG.info(String.format(
                    "Deleted S3 cache file s3://%s/%s", matchedBucketName, filePath));
        }
        queue.deleteMessage(message);
        LOG.info("Deleted from queue message id=" + message.getMessageId());
        jobSubmissionIdToMessageMap.remove(jobSubmission.getId());
    }


    public String getName() {
        return queueName;
    }

    /**
     * Purge all messages from this queue.
     */
    public void purge() {
        try {
            while (true) {
                final Message message = queue.receiveMessage();
                if (message == null) {
                    break;
                }
                queue.deleteMessage(message);
                LOG.info(String.format("Purged queue=%s, message=%s", queueName, message));
            }

        } catch (SQSException e) {
            throw new RuntimeException(e);
        }
    }
}
