1. On UNIX/Linux (and possibly Mac) systems (not necessary for Cygiwn), assuming you are
    using the BASH shell,  Edit the .bash_profile file so that pkgconfig
    will find libs/includes installed "locally"

      export LOCAL_LIB=${HOME}/local-lib
      export PKG_CONFIG_PATH=/usr/lib/pkgconfig:${LOCAL_LIB}/lib/pkgconfig
      export PATH=${LOCAL_LIB}/bin:${PATH}
      export LD_LIBRARY_PATH=${LOCAL_LIB}/lib:${LD_LIBRARY_PATH}

    ************************************************************************
    ** Logout and re-login so these environment variables are set in your **
    ** environment.                                                       **
    ************************************************************************
    
    Make the "local-lib" directories to store local libraries and binaries.

      mkdir -p ${LOCAL_LIB}/lib/pkgconfig/
      mkdir -p ${LOCAL_LIB}/bin/

2. Check your version of autoconf with the command "autoconf --version".
   If you aren't running _at_least_version 2.61, you should update your
   autoconf with the following commands
   
      wget http://ftp.gnu.org/gnu/autoconf/autoconf-2.68.tar.gz
      tar zxvf autoconf-2.68.tar.gz
      cd autoconf-2.68
      ./configure --prefix=${LOCAL_LIB}
      make
      make install

3. Install Protobuf 2.4.1.

      wget http://protobuf.googlecode.com/files/protobuf-2.4.1.tar.gz
      tar zxvf protobuf-2.4.1.tar.gz
      cd protobuf-2.4.1
      #
      # for root or cygwin, don't use the --prefix option
      #
      ./configure --prefix=${LOCAL_LIB}
      make
      make install

4. Download, build, and install the PCRE (Perl Compatible Regular
   Expressions) library (8.21 or later) from http://pcre.org

      wget ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.21.tar.gz
      tar zxvf pcre-8.21.tar.gz
      cd pcre-8.21
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
