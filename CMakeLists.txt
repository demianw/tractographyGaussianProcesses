project(vtkPythonWeave)

cmake_minimum_required(VERSION 2.4)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# --------------------------------------------------------------------------
# Find Slicer3

#if(NOT Slicer3_SOURCE_DIR)
  #find_package(Slicer3 REQUIRED)
  #include(${Slicer3_USE_FILE})
  #    slicer3_set_default_install_prefix_for_external_projects()
#endif(NOT Slicer3_SOURCE_DIR)


find_package(VTK REQUIRED)
include(${VTK_USE_FILE})


configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in
  ${CMAKE_CURRENT_SOURCE_DIR}/setup.py
)

#
#
## --------------------------------------------------------------------------
## Find Slicer3
#
#if(NOT Slicer3_SOURCE_DIR)
#  find_package(Slicer3 REQUIRED)
#  include(${Slicer3_USE_FILE})
#  slicer3_set_default_install_prefix_for_external_projects()
#endif(NOT Slicer3_SOURCE_DIR)


#if(Slicer3_USE_PYTHON)
#
#configure_file(
#  ${CMAKE_CURRENT_SOURCE_DIR}/TractographyGP/tractographyGP/tractographyGP/setup.py.in
#  ${CMAKE_CURRENT_SOURCE_DIR}/TractographyGP/tractographyGP/tractographyGP/setup.py
#)
#
#  add_custom_command(
#    OUTPUT weave_files.txt
#    COMMAND "/usr/bin/env" "python" "${CMAKE_CURRENT_SOURCE_DIR}/TractographyGP/weave/setup.py" "install" "--record" "${CMAKE_CURRENT_BINARY_DIR}/weave_files.txt" 
#    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/TractographyGP/weave"
#    )
#
#
#  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tractographyGP")
#
#  set(PYTHON_VTK_INCLUDE_DIRS "")
#  foreach(ldir ${VTK_INCLUDE_DIRS})
#    set(PYTHON_VTK_INCLUDE_DIRS "${PYTHON_VTK_INCLUDE_DIRS}:${ldir}")
#  endforeach(ldir)
#
#  set(PYTHON_VTK_LIBRARY_DIRS "")
#  foreach(ldir ${VTK_LIBRARY_DIRS})
#    set(PYTHON_VTK_LIBRARY_DIRS "${PYTHON_VTK_LIBRARY_DIRS}:${ldir}")
#  endforeach(ldir)
#
#
#  add_custom_command(
#    OUTPUT TractographyGPlib
#    COMMAND "/usr/bin/env" "python" "${CMAKE_CURRENT_SOURCE_DIR}/TractographyGP/setup.py" "build_ext" "--build-lib" "${CMAKE_CURRENT_BINARY_DIR}/tractographyGP/build" "--include-dirs" "${PYTHON_VTK_INCLUDE_DIRS}" "--library-dirs" "${PYTHON_VTK_LIBRARY_DIRS}"
#     WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/tractographyGP"
#    )
#
#  add_custom_command(
#    OUTPUT TractographyGPsrc
#    COMMAND "/usr/bin/env" "python" "${CMAKE_CURRENT_SOURCE_DIR}/TractographyGP/setup.py" "build" "--build-base" "${CMAKE_CURRENT_BINARY_DIR}/tractographyGP/build" 
#     WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/tractographyGP"
#    )
#
#
## Compile the reader of VTK PolyData "
#  add_custom_command(
#    OUTPUT readVtkPolyData_ext
#    COMMAND "/usr/bin/env" "python" "${CMAKE_CURRENT_BINARY_DIR}/pythonBuildVTK.py" 
#     WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
#    )
#
#  execute_process(
#    #   COMMAND "/usr/bin/env" "python" "-c" "'import sys; from distutils import sysconfig; sys.stderr.write(sysconfig.get_python_lib()+\"\\n\"); sys.stderr.flush()'" 
#   COMMAND "/usr/bin/env python" "-c" "'print 2; f = open(\"t.txt\"); f.write(2); f.close()'"
#   #COMMAND /bin/echo test2
#   OUTPUT_VARIABLE PYTHON_LIBRARY_DIR
#   ERROR_VARIABLE PYTHON_ERROR
#   )
#
# MESSAGE("test ${PEDRO} ${PYTHON_LIBRARY_DIR} ${PYTHON_ERROR}")
#
# file(GLOB PYTHON_READ_VTK_FILES RELATIVE "${CMAKE_CURRENT_BINARY_DIR}" "readVTKPolyData_ext*")
#
#
#  add_custom_target(
#    weave ALL DEPENDS weave_files.txt
#    )
#
#  add_custom_target(
#    vtkPythonWeave ALL DEPENDS weave_files.txt readVtkPolyData_ext
#    )
#
#  add_custom_target(
#    TractographyGP ALL DEPENDS weave_files.txt TractographyGPsrc TractographyGPlib
#    )
#
#  install(CODE
#    "execute_process(
#      COMMAND \"/usr/bin/env\" \"python\" \"${CMAKE_CURRENT_SOURCE_DIR}/TractographyGP/setup.py\" \"install\" 
#      WORKING_DIRECTORY \"${CMAKE_CURRENT_SOURCE_DIR}/tractographyGP\"
#      )"
#  )
#
# install(
#    FILES ${PYTHON_READ_VTK_FILES}
#    DESTINATION ${Slicer3_INSTALL_PLUGINS_BIN_DIR}
#    )
#
#
#
## --------------------------------------------------------------------------
## Copy over the Python files
#
#
##  execute_process(COMMAND "/usr/bin/env" "python" "${CMAKE_CURRENT_BINARY_DIR}/pythonBuildVTK.py" WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
#                 
#
#
##file(GLOB PYTHONFILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.py")
#
## foreach(file ${PYTHONFILES})
##   configure_file(
##     ${CMAKE_CURRENT_SOURCE_DIR}/${file}
##     ${CMAKE_BINARY_DIR}/${Slicer3_INSTALL_PLUGINS_BIN_DIR}/${file}
##     COPYONLY)
## endforeach(file)
#
##    install(
##    FILES ${PYTHONFILES}
##    DESTINATION ${Slicer3_INSTALL_PLUGINS_BIN_DIR}
##    )
#
#endif(Slicer3_USE_PYTHON)
