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

3. Download, build, and install Protobuf (2.3.0 or later) from
    http://code.google.com/p/protobuf/

      wget http://protobuf.googlecode.com/files/protobuf-2.3.0.tar.gz
      tar zxvf protobuf-2.3.0.tar.gz
      cd protobuf-2.3.0
      #
      # for root or cygwin, don't use the --prefix option
      #
      ./configure --prefix=${LOCAL_LIB}
      make
      make install

4. Download, build, and install the PCRE (Perl Compatible Regular
   Expressions) library (8.10 or later) from http://pcre.org

      wget ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.10.tar.gz
      tar zxvf pcre-8.10.tar.gz
      cd pcre-8.10
      #
      # for root or cygwin, don't use the --prefix option
      #
      ./configure --prefix=${LOCAL_LIB}
      make
      make install

5. >>OPTIONAL<< The Boost libraries cause problems on some systems/compilers.
   Portions of the Goby C++ API library OPTIONALLY use the Boost
   library. If you choose to use this, first download and install Boost:

      wget http://downloads.sourceforge.net/project/boost/boost/1.44.0/boost_1_44_0.tar.gz
      tar zxvf boost_1_44_0.tar.gz
      cd boost_1_44_0
      #
      # for root or cygwin, don't use the --prefix option
      #
      ./bootstrap.sh --prefix=${LOCAL_LIB}
      ./bjam install


6. Build the Goby C++ API library, requires the Goby source distribution.
   The following steps install this library:

      wget http://chagall.med.cornell.edu/goby/releases/goby_latest-src.zip
      unzip goby_latest-src.zip
      cd goby_1.8/cpp/
      #
      # If you have chosen to use the Boost library from step #4 above,
      # edit the "configure.ac" file to uncomment the AX_BOOST_ lines.
      #
      chmod +x autogen.sh
      ./autogen.sh
      #
      # For root or cygwin, don't use the --prefix option.
      #
      ./configure --prefix=${LOCAL_LIB}
      make
      make install

   