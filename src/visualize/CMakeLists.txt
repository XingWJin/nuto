set(VisualizeSources
    Cell.cpp
    ComponentName.cpp
    Point.cpp
    UnstructuredGrid.cpp
    XMLWriter.cpp
    DataArray.cpp
)

create_nuto_module(Visualize "${VisualizeSources}")
