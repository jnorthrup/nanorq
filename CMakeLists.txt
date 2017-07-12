cmake_minimum_required(VERSION 3.8)
project(nanorq)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS "-O3 -std=c99 -Wall -mtune=native -march=native -funroll-loops")

#include_directories(after system oblas)
#add_definitions(-DOBLAS_SSE  )
#add_compile_options(-O3)
#add_compile_options(  -funroll-loops   )
#add_compile_options(-march=native   )
#add_compile_options( -mtune=native )


set(CMAKE_VERBOSE_MAKEFILE 1)
set(SOURCE_FILES
        bitmask.c
        chooser.c
        graph.c
        io.c
        nanorq.c
        params.c
        precode.c
        rand.c)

add_subdirectory(oblas )
add_executable(nanorq_encode encode.c ${SOURCE_FILES}  )
add_executable(nanorq_decode decode.c ${SOURCE_FILES}  )
target_link_libraries( nanorq_encode LINK_PUBLIC oblas )
target_link_libraries( nanorq_decode LINK_PUBLIC oblas )