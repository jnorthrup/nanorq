project(oblas)
set(SOURCE_FILES
        oblas_avx.c
        oblas_neon.c
        oblas_classic.c
        oblas_sse.c
        octmat.c)

add_library(oblas ${SOURCE_FILES})
add_executable(tablegen  tablegen.c  tablegen_sse.c ${SOURCE_FILES})
target_include_directories( oblas PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

