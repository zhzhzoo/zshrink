include_directories(${zshrink_SOURCE_DIR}/zshrink)
link_directories(${zshrink_BINARY_DIR}/zshrink)

set(zshrink_SOURCES
    main.c
    zomb_arithmetic_coder.c
    zomb_arithmetic_coder.h
    zomb_datastream.c
    zomb_datastream.h
    zomb_encrypt.c
    zomb_encrypt.h
    zomb_fileio.c
    zomb_fileio.h
    zshrink.c
    zshrink.h
)

add_executable(zshrink ${zshrink_SOURCES})
