FILE(GLOB c_files "*.c")
FILE(GLOB header_files "*.h")

if (BUILD_SHARED)
    add_library(tdigest SHARED ${c_files} ${header_files})
    target_link_libraries(tdigest m)
    target_include_directories(tdigest SYSTEM PUBLIC)# ${CMAKE_CURRENT_SOURCE_DIR})
    set_target_properties(tdigest PROPERTIES PUBLIC_HEADER "${header_files}")
    install(TARGETS tdigest DESTINATION lib${LIB_SUFFIX} PUBLIC_HEADER DESTINATION include)
endif(BUILD_SHARED)

if (BUILD_STATIC) 
    add_library(tdigest_static STATIC ${c_files} ${header_files}) 
    target_link_libraries(tdigest_static m)
    target_include_directories(tdigest_static SYSTEM PUBLIC)# ${CMAKE_CURRENT_SOURCE_DIR})
    set_target_properties(tdigest_static PROPERTIES PUBLIC_HEADER "${header_files}")
    install(TARGETS tdigest_static DESTINATION lib${LIB_SUFFIX} PUBLIC_HEADER DESTINATION include)
endif(BUILD_STATIC)
