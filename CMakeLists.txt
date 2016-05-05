include_directories(/opt/local/include/eigen3)

add_definitions(--std=c++11)
add_executable(test test.cpp sympp.cpp)

#msympp is just an experiment that does not compile
#add_executable(testmat testmat.cpp sympp.cpp msympp.cpp)