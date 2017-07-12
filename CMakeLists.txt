cmake_minimum_required(VERSION 3.8)
project(nanorq)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        oblas/oblas.c
        oblas/oblas_avx.c
        oblas/oblas_neon.c
        oblas/oblas_classic.c
        oblas/oblas_sse.c
        oblas/octmat.c
        bitmask.c
        chooser.c
        graph.c
        io.c
        nanorq.c
        params.c
        precode.c
        rand.c )

add_executable( tablegen oblas/tablegen.c oblas/tablegen_sse.c ${SOURCE_FILES})
add_executable( encode encode.c ${SOURCE_FILES})
add_executable( decode decode.c ${SOURCE_FILES})