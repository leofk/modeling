set(TARGET_NAME minimeshcore)

##
## SOURCE FILES
##
set(SOURCE_FILES
  ##
  ## UTIL
  ##
  util/assert.hpp
  util/console_colors.hpp
  util/numbers.hpp
  util/foldertools.hpp
  util/foldertools.cpp
  util/macros.hpp
  util/disable_irrelevant_warnings.hpp
  util/enable_irrelevant_warnings.hpp
  util/bit_vector.hpp
  util/get_basefilename.hpp
  util/get_basefilename.cpp
  ##
  ## Manifold oriented half edge
  ##
  mohe/mesh_connectivity.hpp
  mohe/mesh_connectivity.cpp
  mohe/mesh_io.hpp
  mohe/mesh_io.cpp
  mohe/mesh_modifier.hpp
  mohe/mesh_modifier.cpp
  mohe/mesh_subdivision.hpp
  mohe/mesh_subdivision.cpp
  mohe/mesh_simplify.hpp
  mohe/mesh_simplify.cpp
  mohe/mesh_fixed_param.hpp
  mohe/mesh_fixed_param.cpp
  mohe/mesh_free_param.hpp
  mohe/mesh_free_param.cpp
  mohe/mesh_deform.cpp
  mohe/mesh_deform.hpp
  mohe/mesh_mapping.cpp
  mohe/mesh_mapping.hpp
 )

 
# create a library that the viewer can use to link to.
add_library(${TARGET_NAME} ${SOURCE_FILES})

##
## Create a folder structure inside the visual studio project
##
minimesh_folders_lib(${TARGET_NAME}  ${SOURCE_FILES})
