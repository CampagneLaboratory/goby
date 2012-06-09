 export GOBY=.
 export AUTOJAR_SCRIPT_DIR=${GOBY}/scripts/autojar
 java -jar ${AUTOJAR_SCRIPT_DIR}/autojar.jar  -A -a  -p ${GOBY}/lib/dsi*.jar -c  ${GOBY}/lib/fastutil*.jar:${GOBY}/lib/dsi*.jar:${GOBY}/goby-io.jar:${GOBY}/lib/log4j-*.jar:${GOBY}/lib/commons-logging*.jar:${GOBY}/lib/JSAP-*.jar -o ${GOBY}/goby-io-igv.jar \
 @${AUTOJAR_SCRIPT_DIR}/classes.lst

rm -fr META-INF

# extract the META-INF directory from original jar file
 jar xf ${GOBY}/goby-io.jar META-INF/
# Include in reduced auto-jar archive:
 jar uf ${GOBY}/goby-io-igv.jar META-INF/

rm -fr META-INF
