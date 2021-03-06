#include "mechanics/mesh/MeshCompanion.h"

#include "mechanics/structures/unstructured/Structure.h"

#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/nodes/NodeEnum.h"

#include "mechanics/interpolationtypes/InterpolationBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"

#include "base/Exception.h"


struct GmshHeader
{
    double version;
    bool isBinary;
    int double_size;
};

struct GmshNode
{
    int id = 0;
    double Coordinates[3];
};

struct GmshElement
{
    int id = 0;
    int type = 0;
    std::vector<int> tags;
    std::vector<int> nodes;
};

struct TmpGroup
{
    std::vector<int> elementIds;
    int interpolationTypeId = -42;
};


GmshHeader ReadGmshHeader(std::ifstream& fileName)
{
    // ignore first line
    std::string line;
    getline(fileName, line);


    GmshHeader header;
    int binary;
    fileName >> header.version;
    fileName >> binary;
    fileName >> header.double_size;

    header.isBinary = binary == 1;

    std::getline(fileName, line); // endl

    if (static_cast<int>(std::floor(header.version)) != 2)
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Incompatible version. Version 2.x required.");

    return header;
}

std::vector<GmshNode> ReadNodesASCII(std::ifstream& fileName)
{
    std::string line;
    getline(fileName, line);
    if (line != "$EndMeshFormat")
        throw NuTo::Exception(__PRETTY_FUNCTION__, "$EndMeshFormat not found.");

    // begin node section
    getline(fileName, line);
    if (line != "$Nodes")
        throw NuTo::Exception(__PRETTY_FUNCTION__, "$Nodes not found.");

    // read number of nodes
    getline(fileName, line);
    int numNodes = std::stoi(line);

    // read node data
    std::vector<GmshNode> nodes(numNodes);
    for (GmshNode& node : nodes)
    {
        fileName >> node.id;
        fileName >> node.Coordinates[0];
        fileName >> node.Coordinates[1];
        fileName >> node.Coordinates[2];
        getline(fileName, line); // endl
    }

    // end node section
    getline(fileName, line);
    if (line != "$EndNodes")
        throw NuTo::Exception(__PRETTY_FUNCTION__, "$EndNodes not found.");

    return nodes;
}

int GetNumNodesPerElementType(int gmshElementTypeId)
{
    const std::vector<int> numNodesPerElement(
            {0, 2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27, 18, 14, 1, 8, 20, 15, 13, 0, 10, 0, 15});
    return numNodesPerElement[gmshElementTypeId];
}

std::vector<GmshElement> ReadElementsASCII(std::ifstream& fileName)
{
    std::string line;
    getline(fileName, line);
    if (line != "$Elements")
        throw NuTo::Exception(__PRETTY_FUNCTION__, "$Elements not found.");

    // read number of elements
    getline(fileName, line);
    int numElements = std::stoi(line);
    std::vector<GmshElement> elements(numElements);


    // read element data
    for (GmshElement& element : elements)
    {
        int numTags;
        fileName >> element.id;
        fileName >> element.type;
        fileName >> numTags;

        int numNodesInThisElement = GetNumNodesPerElementType(element.type);

        element.tags.resize(numTags);
        for (auto& tag : element.tags)
            fileName >> tag;

        element.nodes.resize(numNodesInThisElement);
        for (auto& node : element.nodes)
            fileName >> node;

        getline(fileName, line); // endl;
    }
    // end element section
    getline(fileName, line);
    if (line != "$EndElements")
        throw NuTo::Exception(__PRETTY_FUNCTION__, "$EndElements not found.");

    return elements;
}

std::vector<GmshNode> ReadNodesBinary(std::ifstream& fileName)
{
    std::string line;

    // read the first two lines again
    getline(fileName, line);
    getline(fileName, line);

    // check size of integer
    int one;
    fileName.read((char*)&one, sizeof(int));
    if (one != 1)
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Invalid binary format.");
    fileName.seekg(1, std::ios::cur);

    getline(fileName, line);
    if (line != "$EndMeshFormat")
        throw NuTo::Exception(__PRETTY_FUNCTION__, "$EndMeshFormat not found.");

    // begin node section
    getline(fileName, line);
    if (line != "$Nodes")
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "$Nodes not found.");
    }

    // read number of nodes
    int numNodes;
    fileName >> numNodes;
    getline(fileName, line); // endl

    // read node data
    std::vector<GmshNode> nodes(numNodes);

    for (auto& node : nodes)
    {
        fileName.read((char*)&node.id, sizeof(int));
        fileName.read((char*)node.Coordinates, 3 * sizeof(double));
    }
    // endl
    getline(fileName, line);

    getline(fileName, line);
    if (line != "$EndNodes")
        throw NuTo::Exception(__PRETTY_FUNCTION__, "$EndNodes not found.");

    // begin element section
    getline(fileName, line);
    if (line != "$Elements")
        throw NuTo::Exception(__PRETTY_FUNCTION__, "$Elements not found.");

    return nodes;
}

std::vector<GmshElement> ReadElementsBinary(std::ifstream& fileName)
{
    std::string line;


    // read number of elements
    int numElements;
    fileName >> numElements;
    getline(fileName, line); // endl

    std::vector<GmshElement> elements(numElements);

    // read element data
    int element_type;
    int num_elm_follow;
    int num_tags;

    for (int elemCount = 0; elemCount < numElements; elemCount++)
    {

        // Read element type
        fileName.read((char*)&element_type, sizeof(int));

        // Read num of Elem with the same header
        fileName.read((char*)&num_elm_follow, sizeof(int));

        // set num_elemt_node
        int cur_num_elm_nodes = GetNumNodesPerElementType(element_type);

        // Read numOfTags
        fileName.read((char*)&num_tags, sizeof(int));

        for (int indexH = 0; indexH < num_elm_follow; indexH++)
        {

            // set element type
            elements[elemCount].type = element_type;

            // read element number
            fileName.read((char*)&elements[elemCount].id, sizeof(int));

            elements[elemCount].tags.resize(num_tags);
            elements[elemCount].nodes.resize(cur_num_elm_nodes);

            // read tags
            for (int tagCount = 0; tagCount < num_tags; tagCount++)
                fileName.read((char*)&elements[elemCount].tags[tagCount], sizeof(int));

            // read nodes
            for (int nodeCount = 0; nodeCount < cur_num_elm_nodes; nodeCount++)
                fileName.read((char*)&elements[elemCount].nodes[nodeCount], sizeof(int));

            elemCount += indexH;
        }
    }
    getline(fileName, line); // endl

    // end element section
    getline(fileName, line);
    if (line != "$EndElements")
        throw NuTo::Exception(__PRETTY_FUNCTION__, "$EndElements not found.");
    return elements;
}


std::map<int, int> CreateNodes(NuTo::Structure& s, const std::vector<GmshNode>& gmshNodes)
{
    // create the nodes
    Eigen::VectorXd coordinates;
    switch (s.GetDimension())
    {
    case 2:
        coordinates.resize(2);
        break;
    case 3:
        coordinates.resize(3);
        break;
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Only implemented for 2D and 3D.");
    }
    std::map<int, int> newNodeNumbers;
    for (const auto& gmshNode : gmshNodes)
    {
        coordinates(0) = gmshNode.Coordinates[0];
        coordinates(1) = gmshNode.Coordinates[1];
        if (s.GetDimension() == 3)
            coordinates(2) = gmshNode.Coordinates[2];
        newNodeNumbers[gmshNode.id] = s.NodeCreate(coordinates);
    }
    return newNodeNumbers;
}

std::vector<std::pair<int, int>> NuTo::MeshCompanion::ImportFromGmsh(Structure& s, const std::string& fileName)
{
    std::ifstream file(fileName.c_str(), std::ios::in);
    if (not file.is_open())
    {
        std::cout << fileName << std::endl;
        throw Exception(__PRETTY_FUNCTION__, "Error opening input file "
                                             "" + fileName +
                                                     ""
                                                     " for read access.");
    }


    GmshHeader header = ReadGmshHeader(file);

    std::map<int, TmpGroup> groups;

    std::vector<GmshNode> nodes;
    std::vector<GmshElement> elements;

    if (header.isBinary)
    {
        if (header.double_size != sizeof(double))
            throw NuTo::Exception(__PRETTY_FUNCTION__, "Invalid size of double.");

        // close fileName and open as binary
        file.close();
        file.open(fileName.c_str(), std::ios::in | std::ios::binary);
        if (file.is_open() == false)
            throw Exception(__PRETTY_FUNCTION__, "Error opening input fileName for read access.");
        nodes = ReadNodesBinary(file);
        elements = ReadElementsBinary(file);
    }
    else
    {
        nodes = ReadNodesASCII(file);
        elements = ReadElementsASCII(file);
    }

    auto newNodeNumber = CreateNodes(s, nodes);

    std::vector<int> nodeNumbers;
    for (auto& element : elements)
    {
        nodeNumbers.resize(element.nodes.size());
        for (unsigned int countNode = 0; countNode < element.nodes.size(); countNode++)
            nodeNumbers[countNode] = newNodeNumber[element.nodes[countNode]];
        Interpolation::eShapeType shapeType;
        Interpolation::eTypeOrder typeOrder;

        switch (element.type)
        {
        case 1: //  2-node line in 2d/3d
            shapeType = Interpolation::eShapeType::TRUSSXD;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT1;
            break;
        case 2: // 3-node triangle.
            shapeType = Interpolation::eShapeType::TRIANGLE2D;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT1;
            break;

        case 3: // 4-node quadrangle.
            shapeType = Interpolation::eShapeType::QUAD2D;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT1;
            break;

        case 4: // 4-node tetrahedron.
            shapeType = Interpolation::eShapeType::TETRAHEDRON3D;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT1;
            break;

        case 5: // 8-node hexahedron.
            shapeType = Interpolation::eShapeType::BRICK3D;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT1;
            break;

        //      case 6: // 6-node prism.

        //      case 7: // 5-node pyramid.

        case 8: // 3-node second order line (2 nodes associated with the vertices and 1 with the edge).
            shapeType = Interpolation::eShapeType::TRUSSXD;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT2;
            // ordering is different than in gmsh
            std::swap(nodeNumbers[1], nodeNumbers[2]);
            break;
        case 9: // 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges).
            shapeType = Interpolation::eShapeType::TRIANGLE2D;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT2;
            break;

        case 10: // 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with
            // the face).
            {
                shapeType = Interpolation::eShapeType::QUAD2D;
                typeOrder = Interpolation::eTypeOrder::LOBATTO2;
                // ordering is different than in gmsh, fix this first
                std::vector<int> nodeNumbersGmsh(nodeNumbers);
                nodeNumbers[0] = nodeNumbersGmsh[0];
                nodeNumbers[1] = nodeNumbersGmsh[4];
                nodeNumbers[2] = nodeNumbersGmsh[1];
                nodeNumbers[3] = nodeNumbersGmsh[7];
                nodeNumbers[4] = nodeNumbersGmsh[8];
                nodeNumbers[5] = nodeNumbersGmsh[5];
                nodeNumbers[6] = nodeNumbersGmsh[3];
                nodeNumbers[7] = nodeNumbersGmsh[6];
                nodeNumbers[8] = nodeNumbersGmsh[2];
                break;
            }
        case 11: // 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges).
            shapeType = Interpolation::eShapeType::TETRAHEDRON3D;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT2;
            break;

        case 12: // 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the
            // faces and 1 with the volume).
            {
                shapeType = Interpolation::eShapeType::BRICK3D;
                typeOrder = Interpolation::eTypeOrder::LOBATTO2;
                // ordering is different than in gmsh, fix this first
                std::vector<int> nodeNumbersGmsh = nodeNumbers;
                nodeNumbers[0] = nodeNumbersGmsh[4];
                nodeNumbers[1] = nodeNumbersGmsh[16];
                nodeNumbers[2] = nodeNumbersGmsh[5];
                nodeNumbers[3] = nodeNumbersGmsh[10];
                nodeNumbers[4] = nodeNumbersGmsh[21];
                nodeNumbers[5] = nodeNumbersGmsh[12];
                nodeNumbers[6] = nodeNumbersGmsh[0];
                nodeNumbers[7] = nodeNumbersGmsh[8];
                nodeNumbers[8] = nodeNumbersGmsh[1];
                nodeNumbers[9] = nodeNumbersGmsh[17];
                nodeNumbers[10] = nodeNumbersGmsh[25];
                nodeNumbers[11] = nodeNumbersGmsh[18];
                nodeNumbers[12] = nodeNumbersGmsh[22];
                nodeNumbers[13] = nodeNumbersGmsh[26];
                nodeNumbers[14] = nodeNumbersGmsh[23];
                nodeNumbers[15] = nodeNumbersGmsh[9];
                nodeNumbers[16] = nodeNumbersGmsh[20];
                nodeNumbers[17] = nodeNumbersGmsh[11];
                nodeNumbers[18] = nodeNumbersGmsh[7];
                nodeNumbers[19] = nodeNumbersGmsh[19];
                nodeNumbers[20] = nodeNumbersGmsh[6];
                nodeNumbers[21] = nodeNumbersGmsh[15];
                nodeNumbers[22] = nodeNumbersGmsh[24];
                nodeNumbers[23] = nodeNumbersGmsh[14];
                nodeNumbers[24] = nodeNumbersGmsh[3];
                nodeNumbers[25] = nodeNumbersGmsh[13];
                nodeNumbers[26] = nodeNumbersGmsh[2];
                break;
            }
        //      case 13: // 18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3
        //    with the quadrangular faces).

        //      case 14: // 14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1
        //    with the quadrangular face).

        //      case 15: // 1-node point.

        case 16: // 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges).
            shapeType = Interpolation::eShapeType::QUAD2D;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT2;
            break;

        case 17: // 20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges).
        {
            shapeType = Interpolation::eShapeType::BRICK3D;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT2;
            // ordering is different than in gmsh, fix this first
            std::vector<int> nodeNumbersGmsh = nodeNumbers;
            // vertices
            nodeNumbers[0] = nodeNumbersGmsh[0];
            nodeNumbers[1] = nodeNumbersGmsh[1];
            nodeNumbers[2] = nodeNumbersGmsh[2];
            nodeNumbers[3] = nodeNumbersGmsh[3];
            nodeNumbers[4] = nodeNumbersGmsh[4];
            nodeNumbers[5] = nodeNumbersGmsh[5];
            nodeNumbers[6] = nodeNumbersGmsh[6];
            nodeNumbers[7] = nodeNumbersGmsh[7];

            // "rear" plane
            nodeNumbers[8] = nodeNumbersGmsh[8];
            nodeNumbers[9] = nodeNumbersGmsh[11];
            nodeNumbers[10] = nodeNumbersGmsh[13];
            nodeNumbers[11] = nodeNumbersGmsh[9];

            // "middle" plane
            nodeNumbers[12] = nodeNumbersGmsh[10];
            nodeNumbers[13] = nodeNumbersGmsh[12];
            nodeNumbers[14] = nodeNumbersGmsh[14];
            nodeNumbers[15] = nodeNumbersGmsh[15];

            // "front" plane
            nodeNumbers[16] = nodeNumbersGmsh[16];
            nodeNumbers[17] = nodeNumbersGmsh[18];
            nodeNumbers[18] = nodeNumbersGmsh[19];
            nodeNumbers[19] = nodeNumbersGmsh[17];
            break;
        }

        //      case 18: // 15-node second order prism (6 nodes associated with the vertices and 9 with the edges).

        //      case 19: // 13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges).

        //      case 20: // 9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the
        //    edges)

        case 21: // 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the
            // face)
            shapeType = Interpolation::eShapeType::TRIANGLE2D;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT3;
            break;

        //      case 22: // 12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the
        //    edges)

        case 23: // 15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the
            // face)
            shapeType = Interpolation::eShapeType::TRIANGLE2D;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT4;
            break;

        //      case 24: // 15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the
        //    edges)

        //      case 25: // 21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the
        //    edges, 6 with the face)

        //      case 26: // 4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge)

        //      case 27: // 5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge)

        //      case 28: // 6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge)

        //      case 29: // 20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4
        //    with the faces)

        //      case 30: // 35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges,
        //    12 with the faces, 1 in the volume)

        //      case 31: // 56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24
        //    with the faces, 4 in the volume)

        //      case 92: // 64-node third order hexahedron (8 nodes associated with the vertices, 24 with the edges, 24
        //    with the faces, 8 in the volume)

        //      case 93: // 125-node fourth order hexahedron (8 nodes associated with the vertices, 36 with the edges,
        //    54 with the faces, 27 in the volume)

        default:
            std::cout << "element type in gmsh " << element.type << std::endl;
            throw Exception(__PRETTY_FUNCTION__, "Element type not implemented in the import routine.");
        }

        // get gmsh group id and create a corresponding nuto group if needed
        int groupId = element.tags[0]; // NuTo groupId == gmsh groupId. // This might cause errors if groups exist
        // before the gmsh import.

        // there is one interpolation type for each group. Else: throw
        TmpGroup& tmpGroup = groups[groupId];
        if (tmpGroup.interpolationTypeId == -42) // does not exist
        {
            // create new
            tmpGroup.interpolationTypeId = s.InterpolationTypeCreate(Interpolation::ShapeTypeToString(shapeType));
            s.InterpolationTypeAdd(tmpGroup.interpolationTypeId, Node::eDof::COORDINATES, typeOrder);
        }
        else
        {
            // check if current element matches the interpolation type
            const InterpolationType& interpolationType = *s.InterpolationTypeGet(tmpGroup.interpolationTypeId);
            Interpolation::eShapeType groupShapeType = interpolationType.GetShapeType();
            Interpolation::eTypeOrder groupTypeOrder = interpolationType.Get(Node::eDof::COORDINATES).GetTypeOrder();
            if (groupShapeType != shapeType or groupTypeOrder != typeOrder)
            {
                throw Exception(
                        __PRETTY_FUNCTION__,
                        "ElementType and InterpolationOrder must be equal for all elements in one physical group.");
            }
        }
        tmpGroup.elementIds.push_back(s.ElementCreate(tmpGroup.interpolationTypeId, nodeNumbers));
    }


    std::vector<std::pair<int, int>> ids;
    for (auto& group : groups)
    {
        int groupId = group.first;
        TmpGroup& tmpGroup = group.second;
        s.GroupCreate(groupId, NuTo::eGroupId::Elements);
        for (int elementId : tmpGroup.elementIds)
            s.GroupAddElement(groupId, elementId);

        ids.push_back({groupId, tmpGroup.interpolationTypeId});
    }

    return ids;
}
