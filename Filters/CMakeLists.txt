set(sources vtkH5PartReader.cxx vtkH5PartWriter.cxx)

set(headers vtkH5PartReader.h vtkH5PartWriter.h)

vtk_module_add_module(pvMeshlessFilters SOURCES ${sources} HEADERS ${headers})

vtk_module_link(pvMeshlessFilters PRIVATE H5hut)

vtk_module_include(pvMeshlessFilters PRIVATE
                   "$<BUILD_INTERFACE:${HDF5_INCLUDE_DIRS}>")

vtk_module_definitions(pvMeshlessFilters PRIVATE
                       "$<$<BOOL:${HDF5_IS_PARALLEL}>:PARALLEL_IO>")

paraview_add_server_manager_xmls(XMLS vtkH5PartReader.xml)
