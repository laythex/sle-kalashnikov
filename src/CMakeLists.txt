if (USE_GLOB)
    file(GLOB_RECURSE SOURCES *.hpp *.cpp)
    add_library(${MODULE_NAME} STATIC ${SOURCES})
else ()
    set(HEADERS
        CSRMatrix.hpp
        DenseMatrix.hpp
        TridiagMatrix.hpp
        tools.hpp
        HouseholderQR.hpp
        solvers.hpp
        HessenbergMatrix.hpp
        )

    set(SOURCES
        CSRMatrix.cpp
        DenseMatrix.cpp
        TridiagMatrix.cpp
        tools.cpp
        HouseholderQR.cpp
        solvers.cpp
        HessenbergMatrix.cpp
        )

    source_group("Header Files" FILES ${HEADERS})
    source_group("Source Files" FILES ${SOURCES})

    add_library(${PROJECT_NAME} STATIC ${HEADERS} ${SOURCES})
endif ()

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)