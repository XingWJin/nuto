create_symlink("MeshCompanionGmsh.msh")
create_symlink("MeshCompanionGmshFine.msh")

add_executable(MeshCompanion MeshCompanion.cpp)
target_link_libraries(MeshCompanion Mechanics Math Base
    Boost::unit_test_framework)
if(ENABLE_VISUALIZE)
    target_link_libraries(MeshCompanion Visualize)
endif()
add_test(almost::mechanics::mesh::MeshCompanion
    ${CMAKE_CURRENT_BINARY_DIR}/MeshCompanion)

add_unit_test(Mesh
    mechanics/elements/ElementShapeFunctions)
