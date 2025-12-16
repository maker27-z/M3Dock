cmake_minimum_required(VERSION 3.31)
project(LSHADE_Adam_final)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 清晰指定 Boost 安装路径
set(BOOST_ROOT "/data2/zsf/Packages/boost_1_70_0_install")
set(Boost_NO_SYSTEM_PATHS ON)

set(BOOST_ROOT /data2/zsf/Packages/boost_1_70_0_install)
set(BOOST_INCLUDEDIR ${BOOST_ROOT}/include)
set(BOOST_LIBRARYDIR ${BOOST_ROOT}/lib)

find_package(Boost 1.70 REQUIRED
    COMPONENTS thread filesystem program_options serialization
)

if(Boost_FOUND)
    message(STATUS "Boost include: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost libs: ${Boost_LIBRARIES}")
else()
    message(FATAL_ERROR "Boost not found")
endif()

# 查找 Boost
find_package(Boost 1.70 REQUIRED COMPONENTS
    serialization
    program_options
    filesystem
    thread
)

if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost not found!")
endif()

# include 目录
include_directories(
    ${Boost_INCLUDE_DIRS}
    /home/zsf/anaconda3/include/python3.9
    /home/zsf/anaconda3/lib/python3.9/site-packages/numpy/core/include
)

# 链接目录
link_directories(
    ${Boost_LIBRARY_DIRS}
    /home/zsf/anaconda3/lib
)

add_executable(LSHADE_Adam_final
        adam_cal.h
        array3d.h
        atom.h
        atom_base.h
        atom_constants.h
        atom_type.h
        bfgs.h
        brick.h
        cache.cpp
        cache.h
        common.h
        conf.h
        convert_substring.h
        coords.cpp
        coords.h
        curl.h
        current_weights.cpp
        current_weights.h
        energy.h
        everything.cpp
        everything.h
        file.h
        grid.cpp
        grid.h
        grid_dim.h
        igrid.h
        incrementable.h
        int_pow.h
        lshade.cpp
        lshade.h
        macros.h
        main.cpp
        manifold.cpp
        manifold.h
        matrix.h
        model.cpp
        model.h
        monte_carlo.cpp
        monte_carlo.h
        mutate.cpp
        mutate.h
        my_pid.cpp
        my_pid.h
        naive_non_cache.cpp
        naive_non_cache.h
        non_cache.cpp
        non_cache.h
        parallel.h
        parallel_mc.cpp
        parallel_mc.h
        parallel_progress.cpp
        parallel_progress.h
        parse_error.h
        parse_pdbqt.cpp
        parse_pdbqt.h
        pdb.cpp
        pdb.h
        precalculate.h
        quasi_newton.cpp
        quasi_newton.h
        quaternion.cpp
        quaternion.h
        random.cpp
        random.h
        recent_history.h
        scoring_function.h
        ssd.cpp
        ssd.h
        statistics.h
        szv_grid.cpp
        szv_grid.h
        tee.h
        terms.cpp
        terms.h
        tree.h
        triangular_matrix_index.h
        weighted_terms.cpp
        weighted_terms.h nsga2.cpp nsga2.h SMPSO.cpp SMPSO.h MOEAD.cpp MOEAD.h MPSOD.cpp MPSOD.h xscore.cpp xscore.h file_error.cpp xtools.h molecule.cpp Ligand.cpp Protein.cpp forcefiled.cpp logp.cpp input.cpp basic.cpp SMoG2016.cpp SMoG2016.h
        ad4cache.h ad4cache.cpp MODEA.cpp MODEA.h
        array3d.h array1.h dfire.h dfire_misc.cc misc.h molecule.cc molecule.h protein.cc protein.h restpl.cc restpl.h restype.cc restype.h rotamer.cc rotamer.h xyzatm.cc matplotlibcpp.h weighted_terms_xscore.cpp weighted_terms_xscore.h weighted_terms_DLIGAND2.cpp weighted_terms_DLIGAND2.h)


target_include_directories(LSHADE_Adam_final PRIVATE
    /data2/zsf/Packages/boost_1_70_0_install/include
)

target_link_directories(LSHADE_Adam_final PRIVATE
    /data2/zsf/Packages/boost_1_70_0_install/lib
)

# 链接库（顺序：Boost → Python → 系统库）
target_link_libraries(LSHADE_Adam_final PRIVATE
    Boost::serialization
    Boost::program_options
    Boost::filesystem
#    Boost::thread
#    Threads::Threads
    /home/zsf/anaconda3/lib/libpython3.9.so
)
#target_link_libraries(LSHADE_Adam_final PRIVATE
#    ${Boost_LIBRARIES}
#    /home/zsf/anaconda3/lib/libpython3.9.so
#    pthread
#    m
#)
target_compile_options(LSHADE_Adam_final PRIVATE -pthread)

# 强制 C++17 和 ABI
add_compile_options(-std=c++17 -D_GLIBCXX_USE_CXX11_ABI=1)