// $Id: VisualizeTruss1D2N.cpp 147 2009-12-02 17:12:54Z eckardt4 $

#include "nuto/math/MathException.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"


int main()
{
	try
	{
		// create structure
		NuTo::Structure myStructure(2);
		myStructure.Info();

		// create nodes
		NuTo::FullMatrix<double> Coordinates(2,9);
		NuTo::FullMatrix<double> Displacements(2,9);

		Coordinates(0,0) = 0; Coordinates(1,0) = 0;
		Coordinates(0,1) = 1; Coordinates(1,1) = 0;
		Coordinates(0,2) = 2; Coordinates(1,2) = 0;
		Coordinates(0,3) = 0; Coordinates(1,3) = 1;
		Coordinates(0,4) = 1; Coordinates(1,4) = 1;
		Coordinates(0,5) = 2; Coordinates(1,5) = 1;
		Coordinates(0,6) = 0; Coordinates(1,6) = 2;
		Coordinates(0,7) = 1; Coordinates(1,7) = 2;
		Coordinates(0,8) = 2; Coordinates(1,8) = 2;

		DBG_POSITION_INFO("Coordinates Matrix")
		Coordinates.Info();

		NuTo::FullMatrix<int> Nodes = myStructure.NodesCreate("displacements", Coordinates);

		// create elements
		NuTo::FullMatrix<int> Incidences(4,4);

		// element1
		Incidences(0,0) = Nodes(0,0);
		Incidences(1,0) = Nodes(1,0);
		Incidences(2,0) = Nodes(4,0);
		Incidences(3,0) = Nodes(3,0);
		// element2
		Incidences(0,1) = Nodes(1,0);
		Incidences(1,1) = Nodes(2,0);
		Incidences(2,1) = Nodes(5,0);
		Incidences(3,1) = Nodes(4,0);
		// element3
		Incidences(0,2) = Nodes(3,0);
		Incidences(1,2) = Nodes(4,0);
		Incidences(2,2) = Nodes(7,0);
		Incidences(3,2) = Nodes(6,0);
		// element4
		Incidences(0,3) = Nodes(4,0);
		Incidences(1,3) = Nodes(5,0);
		Incidences(2,3) = Nodes(8,0);
		Incidences(3,3) = Nodes(7,0);

		DBG_POSITION_INFO("Incidence Matrix")
		Coordinates.Info();

	    NuTo::FullMatrix<int> Elements = myStructure.ElementsCreate("Plane2D4N", Incidences);

	    // create constitutive law
	    myStructure.ConstitutiveLawCreate("myMatLin","LinearElastic");
	    myStructure.ConstitutiveLawSetYoungsModulus("myMatLin",10);
	    myStructure.ConstitutiveLawSetPoissonsRatio("myMatLin",0.1);

	    // create section
	    myStructure.SectionCreate("mySection1","PLANE_STRAIN");
	    myStructure.SectionSetThickness("mySection1",0.01);

	    // assign material, section and integration type
	    myStructure.ElementTotalSetIntegrationType("2D4NGauss4Ip");
	    myStructure.ElementTotalSetConstitutiveLaw("myMatLin");
	    myStructure.ElementTotalSetSection("mySection1");

		// visualize element
		myStructure.ExportVtkDataFile("Plane2D4N.vtk","displacements engineering_strain engineering_stress");
	}
	catch (NuTo::MathException& e)
	{
		std::cout << e.ErrorMessage() << std::endl;
	}
	catch (NuTo::Exception& e)
	{
		std::cout << e.ErrorMessage() << std::endl;
	}
	catch (...)
	{
		std::cout << "Unexpected" << std::endl;
	}

    return 0;
}