include_directories(/opt/local/include/eigen3)

add_definitions(--std=c++11)
add_executable(test test.cpp sympp.cpp)

#add_executable(testmat testmat.cpp msympp2.cpp)