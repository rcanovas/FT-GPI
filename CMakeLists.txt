cmake_minimum_required(VERSION 3.5)
project(fdgpi)

set(CMAKE_BUILD_TYPE "Release")

set(CMAKE_CXX_STANDARD 11)

#make sure that boost is installed correctly and that the path are correct
find_package(Boost COMPONENTS program_options REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
		#link_directories(/usr/lib/x86_64-linux-gnu/)
		link_libraries(boost_iostreams boost_system boost_filesystem)
		set(SOURCE_FILES FT_GPI.cpp)
		add_executable(FT-GPI ${SOURCE_FILES})
		target_link_libraries(FT-GPI ${Boost_LIBRARIES})
endif()


