include(ExternalProject)
ExternalProject_Add( blitz++-external
    PREFIX "blitz++-external"
    GIT_REPOSITORY "https://github.com/syntheticpp/blitz.git"
    GIT_TAG "master"

    BINARY_DIR "blitz++-external/src/blitz++-external"

    CONFIGURE_COMMAND pwd
    COMMAND libtoolize --force
    COMMAND aclocal
    COMMAND autoheader
    COMMAND automake --force-missing --add-missing
    COMMAND autoconf
    COMMAND ./configure

    BUILD_COMMAND pwd
    COMMAND make lib



)


