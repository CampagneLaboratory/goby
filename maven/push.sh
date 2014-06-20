#!/bin/bash
if [ "$#" -eq "0" ]; then
  echo "Usage: push.sh (version) install"
  exit 1;
fi
VERSION=$1
COMMAND=$2
mvn -X -e -Dreleases.repository.url.deploy=http://repository.campagnelab.org/artifactory/CampagneLab-SNAPSHOT -Dversion=${VERSION} -f pom-goby-build.xml ${COMMAND}
