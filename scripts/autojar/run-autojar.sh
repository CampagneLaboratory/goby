 export GOBY=.
 export AUTOJAR_SCRIPT_DIR=${GOBY}/scripts/autojar
 java -jar ${AUTOJAR_SCRIPT_DIR}/autojar.jar  -A -a  -p ${GOBY}/lib/dsi*.jar -c  ${GOBY}/lib/fastutil*.jar:${GOBY}/lib/dsi*.jar:${GOBY}/goby-io.jar:${GOBY}/lib/log4j-*.jar:${GOBY}/lib/commons-logging*.jar:${GOBY}/lib/JSAP-*.jar -o ${GOBY}/goby-io-igv.jar \
 @${AUTOJAR_SCRIPT_DIR}/classes.lst
