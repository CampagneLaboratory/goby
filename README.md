Goby is a next-gen data management framework designed to facilitate the implementation of efficient data analysis pipelines. The program is distributed under the GNU General Public License (GPL).

### File formats
Goby provides very efficient file formats to store next-generation sequencing data and intermediary analysis results. The software has been under development and released since 2010. In June 2012, we released Goby 2.0, a version with compression methods that provide state of the art compression of High-Throughput Sequencing (HTS) alignment data.

### Framework
Goby is a framework to help bioinformaticians program efficient analysis tools quickly. The framework was engineered for performance and flexibility. Tools written with Goby often have much better performance and scalability compared to programs developed with other approaches.

### Algorithms
Goby provides efficient algorithms for most computational tasks required when analyzing HTS data. For instance, an ultra-fast local realignment around indels algorithm works directly with Goby HTS alignments and can realign reads on the fly as the alignment is read.

### Authors and Contributors
Goby is currently being developed by the members of the [Campagne laboratory](http://campagnelab.org).

### Source code
Goby source code is now on GitHub.  You can obtain and build the project as follows:
   ```
   git clone git://github.com/CampagneLaboratory/goby.git
   cd goby
   git checkout 2.0 (adjust the version number as needed, or use master for the development branch)
   git submodule update --init   (this will make sure the submodules are fetched in a new repository)
   ```
#### Compilation:
   ```
   ant -f build.xml jar
   ```
#### Testing:
   ```
   ant -f build.xml test
   ```
### Documentation and forums
You will find extensive documentation at [goby.campagnelab.org](http://goby.campagnelab.org).
Questions and feedback should be addressed to the [Goby user forum](https://groups.google.com/forum/?fromgroups#!forum/goby-framework).