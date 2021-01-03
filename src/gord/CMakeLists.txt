project(gord)

file(GLOB SOURCES "./*.cpp")
add_library( gord ${SOURCES} )
target_compile_options(gord PUBLIC -fPIC -fpermissive)
