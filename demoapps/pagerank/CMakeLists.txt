project(GraphLab)
link_libraries(hiredis)
add_graphlab_executable(pagerank pagerank.cpp)
add_graphlab_executable(advanced_pagerank advanced_pagerank.cpp)
