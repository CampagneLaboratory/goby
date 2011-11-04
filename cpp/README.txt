1. On UNIX/Linux/Mac systems (not necessary for Cygiwn), assuming you are
    using the BASH shell,  Edit the .bash_profile file so that pkgconfig
    will find libs/includes installed "locally"

      export LOCAL_LIB=${HOME}/local-lib
      export PKG_CONFIG_PATH=/usr/lib/pkgconfig:${LOCAL_LIB}/lib/pkgconfig
      export PATH=${LOCAL_LIB}/bin:${PATH}
      export LD_LIBRARY_PATH=${LOCAL_LIB}/lib:${LD_LIBRARY_PATH}

    Logout and re-login so these environment variables are set in your
    environment.
    
    Make the "local-lib" directories to store local libraries and binaries.

      mkdir -p ${LOCAL_LIB}/lib/pkgconfig/
      mkdir -p ${LOCAL_LIB}/bin/

2. Check your version of autoconf with the command "autoconf --version".
   If you aren't running at least version 2.61, you should update your
   autoconf with the following commands
   
      wget http://ftp.gnu.org/gnu/autoconf/autoconf-2.68.tar.gz
      tar zxvf autoconf-2.68.tar.gz
      cd autoconf-2.68
      ./configure --prefix=${LOCAL_LIB}
      make
      make install

3. You can use the normally distributed version of Protobuf 2.4.1,
   but we strongly recommend you use a version of Protobuf 2.4.1
   that we have patched to better handle large files.

      wget http://campagnelab.org/files/protobuf-2.4.1-icb.tgz
      tar zxvf protobuf-2.4.1-icb
      cd protobuf-2.4.1-icb
      #
      # for root or cygwin, don't use the --prefix option
      #
      ./configure --prefix=${LOCAL_LIB}
      make
      make install

4. Download, build, and install the PCRE (Perl Compatible Regular
   Expressions) library (8.20 or later) from http://pcre.org

      wget ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.20.tar.gz
      tar zxvf pcre-8.20.tar.gz
      cd pcre-8.20
      #
      # for root or cygwin, don't use the --prefix option
      #
      ./configure --prefix=${LOCAL_LIB}
      make
      make install

5. Build the Goby C++ API library, requires the Goby source distribution.
   The following steps install this library:

      wget http://chagall.med.cornell.edu/goby/releases/goby_latest-cpp.zip
      unzip goby_latest-cpp.zip
      cd goby_VERSION/cpp/
      chmod +x autogen.sh
      ./autogen.sh
      #
      # For root or cygwin, don't use the --prefix option.
      #
      ./configure --prefix=${LOCAL_LIB}
      make
      make install

   
