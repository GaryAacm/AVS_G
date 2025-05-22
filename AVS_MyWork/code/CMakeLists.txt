cmake_minimum_required(VERSION 3.6)
project(AVS_api)

set(CMAKE_CXX_STANDARD 11)

set(FLAGS "${FLAGS} -Wall -Wno-unused-function -Wno-unused-variable -Wno-char-subscripts -Wno-sign-compare -g -march=native")
set(FLAGS "${FLAGS} -DHAVE_PTHREAD -DUSE_MALLOC_WRAPPERS -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -DLIBBSC_OPENMP_SUPPORT")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${FLAGS} -O0 ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${FLAGS} -O0 -std=c++11 ")

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${FLAGS} -Ofast -DNDEBUG ")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${FLAGS} -Ofast -DNDEBUG -std=c++11 ")

# find openmp
FIND_PACKAGE(OpenMP REQUIRED)
message(STATUS "OPENMP FOUND")
set(CMAKE_C_FLAGS_DEBUG  "${CMAKE_C_FLAGS_DEBUG} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} ${OpenMP_CXX_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE  "${CMAKE_C_FLAGS_RELEASE} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")

# find boost library, boost version >= 1.69.
# For CentOS, use "sudo yum install boost169-devel.x86_64" to install boost 1.69, and cmake generation command is :
#  cmake3 -DBOOST_INCLUDEDIR=/usr/include/boost169 -DBOOST_LIBRARYDIR=/usr/lib64/boost169 -DCMAKE_BUILD_TYPE=RELEASE ..
find_package(Boost REQUIRED COMPONENTS system filesystem)

# find zlib
find_package(ZLIB REQUIRED)
if (ZLIB_FOUND)
    add_definitions(-DHAVE_ZLIB)
    include_directories(${ZLIB_INCLUDE_DIRS})
endif()

# add lemon library
set(LEMON_ENABLE_GLPK NO CACHE STRING "Disable GLPK solver backend.")
set(LEMON_ENABLE_COIN NO CACHE STRING "Disable ILOG (CPLEX) solver backend.")
set(LEMON_ENABLE_ILOG NO CACHE STRING "Disable COIN solver backend.")
set(LEMON_ENABLE_SOPLEX NO CACHE STRING "Disable SoPlex solver backend.")
add_subdirectory(lemon EXCLUDE_FROM_ALL)

# add FastPFor library
add_subdirectory(FastPFor EXCLUDE_FROM_ALL)

# add libbsc
set(LIBBSC_SOURCE_FILES
        libbsc/bsc.cpp
        libbsc/libbsc/adler32/adler32.cpp
        libbsc/libbsc/bwt/divsufsort/divsufsort.c
        libbsc/libbsc/bwt/bwt.cpp
        libbsc/libbsc/coder/qlfc/qlfc.cpp
        libbsc/libbsc/coder/qlfc/qlfc_model.cpp
        libbsc/libbsc/coder/coder.cpp
        libbsc/libbsc/filters/detectors.cpp
        libbsc/libbsc/filters/preprocessing.cpp
        libbsc/libbsc/libbsc/libbsc.cpp
        libbsc/libbsc/lzp/lzp.cpp
        libbsc/libbsc/platform/platform.cpp)

message(STATUS "** CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "** CMAKE_C_FLAGS_DEBUG: ${CMAKE_C_FLAGS_DEBUG}")
message(STATUS "** CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "** CMAKE_C_FLAGS_RELEASE: ${CMAKE_C_FLAGS_RELEASE}")

set(MAIN_SOURCE_FILES
    main.cpp
    AcoCompress.cpp
    AlignInfo.cpp
    BufPool.cpp
    CommandParse.cpp
    Encap.cpp
    EncodeSEWorker.cpp
    EncodePEWorker.cpp
    EncodeWorker.cpp
    FastqRead.cpp
    NameProcess.cpp
    Param.cpp
    QualProcess.cpp
    SeqLenProcess.cpp
    RefFreeSeqProcess.cpp
    ThreadPool.cpp
    util.cpp
    DecodeWorker.cpp
    DecodeSEWorker.cpp
    DecodePipeWorker.cpp
    DecodePEWorker.cpp
    ArcFile.cpp
    Context.cpp
    BlockParamProcess.cpp
    ExtractWorker.cpp
    NameLenProcess.cpp
    Md5.cpp
    HashRef.cpp
    HashRefIndex.cpp
    CombineCompress.cpp
    COMTXTCompress.cpp
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/scripts/parallel_commands.sh.in"
        "${CMAKE_BINARY_DIR}/bin/parallel_commands.sh"
)

set(BSC_PATH "./bsc")
set(PARALLEL_COMMANDS_PATH "./parallel_commands.sh")
configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/environment.h.in"
        "${CMAKE_BINARY_DIR}/avs-g/environment.h"
)

include_directories(SYSTEM
        ${BOOST_INCLUDEDIR}
        lemon/
        FastPFor/headers/
        "${CMAKE_BINARY_DIR}/lemon" # lemon/config.h
        "${CMAKE_BINARY_DIR}/avs-g" # avs-g/environment.h
        )

add_executable(AVS_api ${MAIN_SOURCE_FILES})
# add_executable(AVS_api_static ${MAIN_SOURCE_FILES}) # 静态链接版本
add_executable(bsc ${LIBBSC_SOURCE_FILES})

if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    target_link_libraries(AVS_api m pthread ${ZLIB_LIBRARIES}  rt  lemon FastPFor  ${Boost_LIBRARIES})
    # target_link_libraries(AVS_api_static m pthread ${ZLIB_LIBRARIES}  rt  lemon FastPFor ${Boost_LIBRARIES} -static-libstdc++)
else()
    target_link_libraries(AVS_api m pthread ${ZLIB_LIBRARIES} lemon FastPFor ${Boost_LIBRARIES})
    # target_link_libraries(AVS_api_static m pthread ${ZLIB_LIBRARIES} lemon FastPFor ${Boost_LIBRARIES} -static-libstdc++)
endif()
