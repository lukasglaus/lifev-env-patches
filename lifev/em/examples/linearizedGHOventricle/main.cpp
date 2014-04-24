#include <lifev/core/LifeV.hpp>
#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicMinimalModel.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/em/solver/EMStructuralOperator.hpp>
#include <lifev/em/solver/EMNeoHookeanActivatedMaterial.hpp>
#include <lifev/em/solver/EMGeneralizedActiveHolzapfelOgdenMaterial.hpp>
#include <lifev/em/solver/EMETAFunctors.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/core/interpolation/RBFhtp.hpp>
#include <lifev/core/interpolation/RBFhtpVectorial.hpp>
#include <lifev/core/mesh/MeshLoadingUtility.hpp>
#include <lifev/core/mesh/MeshTransformer.hpp>
#include <lifev/core/interpolation/RBFlocallyRescaledVectorial.hpp>
#include <lifev/core/interpolation/RBFlocallyRescaledScalar.hpp>
#include <lifev/core/interpolation/RBFrescaledVectorial.hpp>
#include <lifev/core/interpolation/RBFrescaledScalar.hpp>
//#include <lifev/core/interpolation/RBFscalar.hpp>
#include <lifev/core/interpolation/RBFvectorial.hpp>

#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>
#include <sys/stat.h>
#include <fstream>

#include <lifev/em/solver/EMSolver.hpp>
#include <lifev/em/solver/EMActiveStrainSolver.hpp>
#include <lifev/core/fem/GradientRecovery.hpp>

using namespace LifeV;

void EpetraPow ( VectorEpetra& vector, const Real p )
{
	Int size = vector.epetraVector().MyLength();

	for(int j(0); j < size; j++ )
	{
	int gid = vector.blockMap().GID(j);
	vector[gid] = std::pow(vector[gid],p);
	}
}
void EpetraSqrt ( VectorEpetra& vector )
{
	Int size = vector.epetraVector().MyLength();

	for(int j(0); j < size; j++ )
	{
	int gid = vector.blockMap().GID(j);
	vector[gid] = std::sqrt(vector[gid]);
	}
}


Real bcZero (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  0.;
}
Real d0(const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  0.;
}

Real initialVlid(const Real& /*t*/, const Real&  X, const Real& /*Y*/, const Real& Z, const ID& /*i*/)
{
	if( X > 0.95 ) return 1.0;
	else return  0.;
}

Real initialV0left(const Real& /*t*/, const Real&  X, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
	if( X == 0 ) return 1.0;
	else return  0.;
}


Real FL(const Real I4f)
{
    	if(I4f > 0.87277 && I4f < 1.334)
    	{
			Real d0 = -4.333618335582119e3;
			Real d1 = 2.570395355352195e3;
			Real e1 = -2.051827278991976e3;
			Real d2 = 1.329536116891330e3;
			Real e2 = 0.302216784558222e3;
			Real d3 = 0.104943770305116e3;
			Real e3 = 0.218375174229422e3;
			Real l0 = 1.95;

			Real Force = d0/2 + d1 * std::sin(I4f * l0)
							  + e1 * std::cos(I4f * l0)
							  + d2 * std::sin(2 * I4f * l0)
							  + e2 * std::cos(2 * I4f * l0)
							  + d3 * std::sin(3 * I4f * l0)
							  + e3 * std::cos(3 * I4f * l0);
			return Force;
    	}
    	else
    		return 0.0;
}


Real fiberRotationRing(const Real& /*t*/, const Real&  X, const Real& Y, const Real&/*Z*/, const ID& i)
{
	Real R = std::sqrt( X * X + Y * Y);
	//Real teta = std::atan( Y / X );
	Real fz = 0.0;
	Real fx =  Y / R;
	Real fy = - X / R;
	Real sx = X / R;
	Real sy = Y / R;
	Real m = -1.9040;
	Real q = 3.5224;
	Real theta = m * R + q;

//	f01a f001*cos(teta)+f001*s01^2*(1-cos(teta))+s01*s02*f002*(1-cos(teta))
//	f02a s01*s02*f001*(1-cos(teta))+f002*cos(teta)+f002*s02^2*(1-cos(teta))
//	f03a s01*f002*sin(teta)-s02*f001*sin(teta)

    switch (i)
    {
        case 0:
            return  fx * std::cos(theta) + fx * sx * sx * ( 1.0  - std::cos(theta) ) + sx * sy * fy * ( 1.0  - std::cos(theta) );
            break;
        case 1:
            return sx * sy * fy *  ( 1.0  - std::cos(theta) ) + fy * std::cos(theta) + fy * sy * sy * ( 1.0  - std::cos(theta) ) ;
            break;
        case 2:
            return sx * fy * std::sin(theta) - sy * fx * std::sin(theta);
            break;
        default:
            ERROR_MSG ("This entry is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }

}

static Real f0fun(const Real&, const Real& x, const Real&, const Real& , const ID& comp)
{
	Real p = 3.14159265358979;
    Real alpha = -2.0 * p/3. * x + p/3.;

    Real compx = 0.0;
    Real compy = std::cos(alpha);
    Real compz = std::sin(alpha);

//    compx = 0.0;
//    compy = 0.0;
//    compz = 1.0;
    if (comp == 0)
        return compx;
    else if (comp == 1)
        return compy;
    else
        return compz;
}

Real rescalingGamma(const Real&, const Real& x, const Real&, const Real& , const ID& /*comp*/)
{
	Real r = (1-x);
	Real p = 3.14159265358979;
	Real alpha = -2.0 * p/3. * x + p/3.;
	Real circumferentialShortening = (0.05-0.2) * r +0.2;
	Real max_gamma =  circumferentialShortening / std::cos(alpha);
	return max_gamma;
}


void createPositionVector(const RegionMesh<LinearTetra>& mesh,  VectorEpetra& positionVector)
{
	Int nLocalDof = positionVector.epetraVector().MyLength();
	Int nComponentLocalDof = nLocalDof / 3;
	for (int k(0); k < nComponentLocalDof; k++)
	{
		UInt iGID = positionVector.blockMap().GID(k);
		UInt jGID = positionVector.blockMap().GID(k + nComponentLocalDof);
		UInt kGID = positionVector.blockMap().GID(k + 2 * nComponentLocalDof);

		positionVector[iGID] = mesh.point(iGID).x();
		positionVector[jGID] = mesh.point(iGID).y();
		positionVector[kGID] = mesh.point(iGID).z();
	}


}

void createPositionXVector(const RegionMesh<LinearTetra>& mesh,  VectorEpetra& positionVector)
{
	Int nLocalDof = positionVector.epetraVector().MyLength();
	Int nComponentLocalDof = nLocalDof / 3;
	for (int k(0); k < nComponentLocalDof; k++)
	{
		UInt iGID = positionVector.blockMap().GID(k);
		positionVector[iGID] = mesh.point(iGID).x();
	}


}

void computeX( VectorEpetra& positionVector )
{
	Int nLocalDof = positionVector.epetraVector().MyLength();
	Int nComponentLocalDof = nLocalDof / 3;
	for (int k(0); k < nComponentLocalDof; k++)
	{
		UInt iGID = positionVector.blockMap().GID(k);
		UInt jGID = positionVector.blockMap().GID(k + nComponentLocalDof);
		UInt kGID = positionVector.blockMap().GID(k + 2 * nComponentLocalDof);

		positionVector[jGID] = 0.0;
		positionVector[kGID] = 0.0;
	}


}


void computeI4( VectorEpetra& i4, VectorEpetra& sx, VectorEpetra& sy, VectorEpetra& sz, VectorEpetra& f )
{
	i4 *= 0.0;
	Int nLocalDof = i4.epetraVector().MyLength();
	//Int nComponentLocalDof = nLocalDof / 3;
	for (int k(0); k < nLocalDof; k++)
	{
		UInt iGID = sx.blockMap().GID(k);
		UInt jGID = sx.blockMap().GID(k + nLocalDof);
		UInt kGID = sx.blockMap().GID(k + 2 * nLocalDof);

		Real fx, fy, fz;
		Real F11 = sx[iGID] + 1.0;
		Real F12 = sy[iGID];
		Real F13 = sz[iGID];
		Real F21 = sx[jGID];
		Real F22 = sy[jGID] + 1.0;
		Real F23 = sz[jGID];
		Real F31 = sx[kGID];
		Real F32 = sy[kGID];
		Real F33 = sz[kGID] + 1.0;
		fx = F11 * f[iGID];
		fx += ( F12 * f[jGID] );
		fx += ( F13 * f[kGID] );
		fy = F21 * f[iGID];
		fy += ( F22 * f[jGID] );
		fy += ( F23 * f[kGID] );
		fz = F31 * f[iGID];
		fz += ( F32 * f[jGID] );
		fz += ( F33 * f[kGID] );
		Real J = F11 * (F22 * F33 - F32 * F23) - F22 * (F21 * F33 - F31 * F23) + F33 * (F21 * F32 - F31 * F22);

		i4[iGID] = fx * fx + fy * fy + fz * fz;
		i4[iGID] *= std::pow(J,-2.0/3.0); 
		
	}


}

template<typename space> Real ComputeVolume(     const boost::shared_ptr<RegionMesh<LinearTetra> > localMesh,
						VectorEpetra positionVector,
						const VectorEpetra& disp,
						const boost::shared_ptr<ETFESpace<RegionMesh<LinearTetra>, MapEpetra, 3, 1> > ETFESpace,
						const boost::shared_ptr< space > dETFESpace,
						int bdFlag,
						boost::shared_ptr<Epetra_Comm>  comm)
{

	Real fluidVolume;

    MatrixSmall<3,3> Id;
    Id(0,0) = 1.; Id(0,1) = 0., Id(0,2) = 0.;
    Id(1,0) = 0.; Id(1,1) = 1., Id(1,2) = 0.;
    Id(2,0) = 0.; Id(2,1) = 0., Id(2,2) = 1.;
    VectorSmall<3> E1;
    E1(0) = 1.;
    E1(1) = 0.;
    E1(2) = 0.;

	Int nLocalDof = positionVector.epetraVector().MyLength();
	for (int k(0); k < nLocalDof; k++)
	{
		UInt iGID = positionVector.blockMap().GID(k);

		positionVector[iGID] += disp[iGID];
	}

	boost::shared_ptr<VectorEpetra> intergral( new VectorEpetra( positionVector.map() ) );

	{
	  	using namespace ExpressionAssembly;

	  	BOOST_AUTO_TPL(I,      value(Id) );
	  	BOOST_AUTO_TPL(vE1,      value(E1) );
	  	BOOST_AUTO_TPL(Grad_u, grad( dETFESpace, disp, 0) );
	    BOOST_AUTO_TPL(x,     value( ETFESpace, positionVector ) );
	  	BOOST_AUTO_TPL(F,      ( Grad_u + I ) );
	  	BOOST_AUTO_TPL(FmT,    minusT(F) );
	  	BOOST_AUTO_TPL(J,       det(F) );
	  	BOOST_AUTO_TPL(x1,     dot(x, vE1) );

    	QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria4pt) );

        *intergral *= 0.0;
        integrate ( boundary ( localMesh, bdFlag),
        		    myBDQR,
        		    ETFESpace,
        		    value(-1.0) * J * dot( vE1, FmT * Nface ) *phi_i
                  ) >> intergral;

        intergral -> globalAssemble();
//        *position = *positionR;
//        for
//        fluidVolume = position ->

        fluidVolume =  positionVector.dot(*intergral);

	if( comm->MyPID() == 0 )
		{
			std::cout << "\nFluid volume: " << fluidVolume << " in processor " << comm->MyPID() << std::endl;
		}
        return fluidVolume;
	}
}

Real ComputeCubicVolume( const RegionMesh<LinearTetra> fullMesh,
							VectorEpetra positionVector,
					  const VectorEpetra& disp,
					  int bdFlag,
					  boost::shared_ptr<Epetra_Comm>  comm)
{
	positionVector += disp;
	Real xMin(0.);
	Real xMax(0.);
	Real yMin(0.);
	Real yMax(0.);
	Real zMin(0.);
	Real zMax(0.);
	Int nLocalDof = positionVector.epetraVector().MyLength();
	Int nComponentLocalDof = nLocalDof / 3;

	for (int k(0); k < nComponentLocalDof; k++)
	{
		UInt iGID = positionVector.blockMap().GID(k);
		UInt jGID = positionVector.blockMap().GID(k + nComponentLocalDof);
		UInt kGID = positionVector.blockMap().GID(k + 2 * nComponentLocalDof);
		if( fullMesh.point ( iGID ).markerID() == bdFlag )
		{
			if(positionVector[iGID] > xMax) xMax = positionVector[iGID];
			if(positionVector[iGID] < xMin) xMin = positionVector[iGID];
			if(positionVector[jGID] > yMax) yMax = positionVector[jGID];
			if(positionVector[jGID] < yMin) yMin = positionVector[jGID];
			if(positionVector[kGID] > zMax) zMax = positionVector[kGID];
			if(positionVector[kGID] < zMin) zMin = positionVector[kGID];

		}
	}


//	Real p = 3.14159265358979;
//	Real Volume;
//	if(xDiameter > yDiameter) Volume = p * xDiameter * xDiameter / 4.0 * (zMax - zMin) / 3.0;
//	else Volume = p * yDiameter * xDiameter / 4.0 * (zMax - zMin) / 3.0;
	int numProc = comm -> NumProc();
	Real xMinGlobal(0.);
	Real xMaxGlobal(0.);
	Real yMinGlobal(0.);
	Real yMaxGlobal(0.);
	Real zMinGlobal(0.);
	Real zMaxGlobal(0.);

	if(numProc == 1)
	{
		Real xDiameter = (xMax - xMin );
		Real yDiameter = (yMax - yMin );
		Real Volume = xDiameter * yDiameter * (zMax - zMin) / 2;
		return Volume;
	}
	else
	{
		comm -> MaxAll(&xMax, &xMaxGlobal,1);
		comm -> MinAll(&xMin, &xMinGlobal,1);
		comm -> MaxAll(&yMax, &yMaxGlobal,1);
		comm -> MinAll(&yMin, &yMinGlobal,1);
		comm -> MaxAll(&zMax, &zMaxGlobal,1);
		comm -> MinAll(&zMin, &zMinGlobal,1);

		Real xDiameter = (xMaxGlobal - xMinGlobal );
		Real yDiameter = (yMaxGlobal - yMinGlobal );
		Real Volume = xDiameter * yDiameter * (zMaxGlobal - zMinGlobal) / 2;
		return Volume;
	}

}


Real evaluatePressure(Real Volume, Real dV, Real pn, Real dp_temporal, Real dV_temporal, Real Cp, Int phase = 0, Real t = 0.0, Real dt = 1.0 )
{
	Real pressure;
        Real R = 150.0 * 1333.22; //mmHg ms / ml
	Real C = 1 / 1333.22; //ml / mmHg

	switch(phase)
	{
		case 0:
			pressure = 1.0 * pn + dV_temporal / Cp;
			break;
		case 1:	
			pressure = pn - ( pn * dt / C / R + dV_temporal / C );
	//		pressue = pn + dt / C ;			
			break;
		case 2:
			pressure = 1.0 * pn + dV_temporal / Cp;
			break;
		case 3:	
			pressure = pn;
//			R = 10.0 * 1333.0;
//			C = 10.0 / 1333.0; 
//			pressure = pn - ( pn * dt / C / R - dV_temporal / C );
//			if(t > 250)pressure = 12000 + 0.35 * (t-250.0) * (t-250.0);
//			else pressure = pn;
			break;
		default:
			std::cout << "\nThis case is not yet available\n I'm not gonna change the pressure.\n";
			pressure = pn;
			break;
	}	
	return pressure;
}

int main (int argc, char** argv)
{

    typedef RegionMesh<LinearTetra>                         mesh_Type;
    typedef boost::shared_ptr<mesh_Type>                    meshPtr_Type;
    typedef boost::function < Real (const Real& /*t*/,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real& /*z*/,
                                    const ID&   /*i*/ ) >   function_Type;
    typedef IonicMinimalModel					ionicModel_Type;
    typedef boost::shared_ptr< ionicModel_Type >  ionicModelPtr_Type;

    typedef EMMonodomainSolver< mesh_Type, ionicModel_Type >        monodomainSolver_Type;
    typedef boost::shared_ptr< monodomainSolver_Type >  monodomainSolverPtr_Type;
    typedef VectorEpetra				vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

	typedef MatrixEpetra<Real> matrix_Type;
	typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;

    typedef BCHandler                                          bc_Type;
    typedef boost::shared_ptr< bc_Type >                       bcPtr_Type;
    typedef StructuralOperator< RegionMesh<LinearTetra> >		physicalSolver_Type;
    typedef BCInterface3D< bc_Type, physicalSolver_Type >              bcInterface_Type;
    typedef boost::shared_ptr< bcInterface_Type >              bcInterfacePtr_Type;



#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
#endif

    boost::shared_ptr<Epetra_Comm>  comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    //*********************************************//
    // creating output folder
    //*********************************************//
    GetPot commandLine ( argc, argv );
    std::string problemFolder = commandLine.follow ( "Output", 2, "-o", "--output" );
    // Create the problem folder
    if ( problemFolder.compare ("./") )
    {
        problemFolder += "/";

        if ( comm->MyPID() == 0 )
        {
            mkdir ( problemFolder.c_str(), 0777 );
        }
    }

  	//===========================================================
  	//===========================================================
  	//				ELECTROPHYSIOLOGY
  	//===========================================================
  	//===========================================================

    if ( comm->MyPID() == 0 )
    {
        cout << "% using MPI" << endl;
    }

    //********************************************//
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//

    if ( comm->MyPID() == 0 )
    {
        std::cout << "Importing parameters list...";
    }
    Teuchos::ParameterList parameterList = * ( Teuchos::getParametersFromXmlFile ( "ParamList.xml" ) );
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }


    //********************************************//
    // In the parameter list we need to specify   //
    // the mesh name and the mesh path.           //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Reading Mesh Name and Path...\n";
    }

    std::string meshName = parameterList.get ("mesh_name", "lid16.mesh");
    std::string meshPath = parameterList.get ("mesh_path", "./");

    Real Cp = parameterList.get ("Cp", 1e-5);
    Real pvtol = parameterList.get ("pvtol", 1e-3);
//    meshPtr_Type mesh ( new mesh_Type ( comm ) );
//    meshPtr_Type fullMesh ( new mesh_Type ( comm ) );
//    MeshUtility::fillWithFullMesh (mesh, fullMesh, meshName, meshPath);
    if ( comm->MyPID() == 0 )
        {
            std::cout << " Done!" << endl;
        }

    //********************************************//
    // We need the GetPot datafile for to setup   //
    // the preconditioner.                        //
    //********************************************//
    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);

    //********************************************//
    // Creates a new model object representing the//
    // model from Aliev and Panfilov 1996.  The   //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Building Constructor for Minimal Model with parameters ... ";
    }
    ionicModelPtr_Type  ionicModel ( new ionicModel_Type() );
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }

    //********************************************//
    // set up the monodomain solver               //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Building Monodomain Solvers... ";
    }

    monodomainSolverPtr_Type monodomain ( new monodomainSolver_Type ( meshName, meshPath, dataFile, ionicModel ) );
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Splitting solver done... ";
    }

    bool load4restart = parameterList.get("load4restart", false);
//    ionicModel -> initialize( monodomain -> globalSolution() );

    if(load4restart)
    {
    	std::string V0filename = parameterList.get("V0filename", "V0");
    	std::string V0fieldname = parameterList.get("V0fieldname", "V0");
    	ElectrophysiologyUtility::importScalarField(monodomain -> globalSolution().at(0),V0filename,V0fieldname,monodomain -> localMeshPtr() );
    	std::string V1filename = parameterList.get("V1filename", "V1");
    	std::string V1fieldname = parameterList.get("V1fieldname", "V1");
    	ElectrophysiologyUtility::importScalarField(monodomain -> globalSolution().at(1),V1filename,V1fieldname,monodomain -> localMeshPtr() );
    	std::string V2filename = parameterList.get("V2filename", "V2");
    	std::string V2fieldname = parameterList.get("V2fieldname", "V2");
    	ElectrophysiologyUtility::importScalarField(monodomain -> globalSolution().at(2),V2filename,V2fieldname,monodomain -> localMeshPtr() );
    	std::string V3filename = parameterList.get("V3filename", "V3");
    	std::string V3fieldname = parameterList.get("V3fieldname", "V3");
    	ElectrophysiologyUtility::importScalarField(monodomain -> globalSolution().at(3),V3filename,V3fieldname,monodomain -> localMeshPtr() );
    }
    else{

		monodomain -> setInitialConditions();

//		for(int i(0); i < ionicModel -> Size(); i++ )
//		{
//		std::cout << "Norm Inf variable " << i  << " = " <<  (  *( monodomain -> globalSolution().at(i) ) ).normInf() << std::endl;
//		}


		ElectrophysiologyUtility::setValueOnBoundary( *(monodomain -> potentialPtr() ), monodomain -> fullMeshPtr(), 1.0, 10 );
		//function_Type Vlid = &initialVlid;
		//monodomain -> setPotentialFromFunction( Vlid );

    }

    std::cout << "Norm Inf potential = " <<  (  *( monodomain -> globalSolution().at(0) ) ).normInf() << std::endl;

    monodomain -> setParameters ( parameterList );

    //********************************************//
    // Creating exporters to save the solution    //
    //********************************************//
    ExporterHDF5< RegionMesh <LinearTetra> > expElectro;

    for(int pid(0); pid < 4 ; pid ++){
    if ( comm->MyPID() == pid )
    {
        cout << "\nExporter setup:  " ;
    }
    }
    monodomain -> setupExporter ( expElectro, parameterList.get ("ElectroOutputFile", "ElectroOutput") );
    expElectro.setPostDir ( problemFolder );
    if ( comm->MyPID() == 0 )
    {
        cout << "\nExport at 0:  " ;
    }

    monodomain -> exportSolution ( expElectro, 0.0 );

    if ( comm->MyPID() == 0 )
    {
        cout << "\nsolve system:  " ;
    }

  	//===========================================================
  	//===========================================================
  	//				SOLID MECHANICS
  	//===========================================================
  	//===========================================================



    if ( comm->MyPID() == 0 )
    {
        std::cout << "monodomain: passed!" << std::endl;
    }

    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >               solidFESpace_Type;
    typedef boost::shared_ptr<solidFESpace_Type>                        solidFESpacePtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 1 >       scalarETFESpace_Type;
    typedef boost::shared_ptr<scalarETFESpace_Type>                      scalarETFESpacePtr_Type;
    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >       solidETFESpace_Type;
    typedef boost::shared_ptr<solidETFESpace_Type>                      solidETFESpacePtr_Type;
    if ( comm->MyPID() == 0 )
    {
        std::cout << "\n\ninitialization bc handler" << std::endl;
    }


    for(int pid(0); pid < 4 ; pid ++){
    if ( comm->MyPID() == pid )
    {
        std::cout << "\nparameters" << std::endl;
    }
    }
    Real rho, poisson, young, bulk, alpha, gammai, mu;
    rho     = dataFile ( "solid/physics/density", 1. );
    young   = dataFile ( "solid/physics/young",   1. );
    poisson = dataFile ( "solid/physics/poisson", 1. );
    bulk    = dataFile ( "solid/physics/bulk",    1. );
    alpha   = dataFile ( "solid/physics/alpha",   1. );
    gammai   = dataFile ( "solid/physics/gamma",   1. );
    mu   = dataFile ( "solid/physics/mu",   1. );
  //  M_gammaf  = dataFile ( "solid/physics/gammaf",  0. );

    if ( comm->MyPID() == 0 )
        {
    std::cout << "density = " << rho     << std::endl
              << "young   = " << young   << std::endl
              << "poisson = " << poisson << std::endl
              << "bulk    = " << bulk    << std::endl
              << "alpha   = " << alpha   << std::endl
              << "gamma   = " << gammai   << std::endl;
        }

    for(int pid(0); pid < 4 ; pid ++){
    if ( comm->MyPID() == pid )
    {
        std::cout << "\ninitialization constitutive law" << std::endl;
    }
    }
    boost::shared_ptr<StructuralConstitutiveLawData> dataStructure (new StructuralConstitutiveLawData( ) );
    dataStructure->setup (dataFile);

    if ( comm->MyPID() == 0 )
    {
        std::cout << "setup spaces" << std::endl;
    }

    meshPtr_Type fullSolidMesh;
    meshPtr_Type localSolidMesh;

    std::string solidMeshName = parameterList.get ("solid_mesh_name", "no_solid_mesh");
    bool usingDifferentMeshes = true;
    if(solidMeshName == "no_solid_mesh" || solidMeshName == meshName )
    	{
    	usingDifferentMeshes = false;
    	}
    else
    {
        if ( comm->MyPID() == 0 )
        {
            std::cout << "\nI am using different meshes: " << meshName << " for the electro and " << solidMeshName << " for the solid\n\n" << std::endl;
        }
    }

    if( usingDifferentMeshes  )
    {
    	fullSolidMesh.reset(new mesh_Type ( comm ) );
    	localSolidMesh.reset(new mesh_Type ( comm ) );
    	MeshUtility::loadMesh (localSolidMesh, fullSolidMesh,  solidMeshName,  parameterList.get ("solid_mesh_path", "") );
    }
    else
    {
    	fullSolidMesh = monodomain -> fullMeshPtr();
    	localSolidMesh = monodomain -> localMeshPtr();
    }

    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (localSolidMesh, dOrder, 3, comm) );
    solidFESpacePtr_Type aFESpace ( new solidFESpace_Type (monodomain -> localMeshPtr(), dOrder, 1, comm) );
    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (localSolidMesh, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), comm) );
    scalarETFESpacePtr_Type aETFESpace ( new scalarETFESpace_Type (monodomain -> localMeshPtr(), & (aFESpace->refFE() ), & (aFESpace->fe().geoMap() ), comm) );
    solidFESpacePtr_Type solidaFESpace ( new solidFESpace_Type (localSolidMesh, "P1", 1, comm) );

    for(int pid(0); pid < 4 ; pid ++){
    if ( comm->MyPID() == pid )
    {
        std::cout << "\nsetup boundary conditions" << std::endl;
    }
    }
    bcInterfacePtr_Type                     solidBC( new bcInterface_Type() );
    solidBC->createHandler();
    solidBC->fillHandler ( data_file_name, "solid" );

    //M_FSIoperator->setSolidBC ( M_solidBC->handler() );
    for(int pid(0); pid < 4 ; pid ++){
    if ( comm->MyPID() == pid )
    {
        std::cout << "\nsetup structural operator" << std::endl;
    }
    }
    //! 1. Constructor of the structuralSolver
     EMStructuralOperator< RegionMesh<LinearTetra> > solid;
     solid.setup (dataStructure,
                  dFESpace,
                  dETFESpace,
                  solidBC -> handler(),
                  comm);
     for(int pid(0); pid < 4 ; pid ++){
     if ( comm->MyPID() == pid )
     {
         std::cout << "\ninitial guess" << std::endl;
     }
     }
     solid.setDataFromGetPot (dataFile);

     if(load4restart)
     {
     	std::string Dfilename = parameterList.get("Gfilename", "G");
     	std::string Dfieldname = parameterList.get("Gfieldname", "G");
     	ElectrophysiologyUtility::importVectorField( solid.displacementPtr(), Dfilename, Dfieldname,  localSolidMesh );
     }
    	//===========================================================
    	//===========================================================
    	//				FIBERS
    	//===========================================================
    	//===========================================================
     for(int pid(0); pid < 4 ; pid ++){
     if ( comm->MyPID() == pid )
     {
         std::cout << "\nreading fibers ... " << std::endl;
     }
     }


     function_Type fibersDirection = &f0fun;
     vectorPtr_Type solidFibers( new vector_Type( dFESpace -> map() ) );
     MPI_Barrier(MPI_COMM_WORLD);
     dFESpace -> interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( fibersDirection ), *solidFibers , 0);
     ElectrophysiologyUtility::normalize(*solidFibers);

     MPI_Barrier(MPI_COMM_WORLD);


     if ( comm->MyPID() == 0 )
     {
         std::cout << "\nread fibers" << std::endl;
     }

     ElectrophysiologyUtility::importVectorField(solidFibers,  parameterList.get ("solid_fiber_file", ""),  parameterList.get ("solid_fiber_field", ""), localSolidMesh );
     //template<typename Mesh> inline void importVectorField(  boost::shared_ptr<VectorEpetra> vector, const std::string& fileName, const std::string& fieldName, boost::shared_ptr< Mesh > localMesh  )

//     std::vector<Real> fvec(3, 0.0);
//     fvec.at(0)  = parameterList.get ("fiber_X", 1.0);
//     fvec.at(1)  = parameterList.get ("fiber_Y", 0.0);
//     fvec.at(2)  = parameterList.get ("fiber_Z", 0.0);
//     HeartUtility::setupFibers(*solidFibers, fvec);
     solid.activeMaterial() -> setFiberVector( *solidFibers );


//     Real sx = 1.0;//parameterList.get ("sheet_X", 1.0);
//     Real sy = 0.0;//parameterList.get ("sheet_Y", 0.0);
//     Real sz = 0.0;//parameterList.get ("sheet_Z", 0.0);;
     vectorPtr_Type solidSheets( new vector_Type( solidFibers -> map() ) );
     ElectrophysiologyUtility::importVectorField(solidSheets,  parameterList.get ("solid_sheets_file", ""),  parameterList.get ("solid_sheets_field", ""), localSolidMesh );

     solid.activeMaterial()->setSheetVector(*solidSheets);

     if ( comm->MyPID() == 0 )
     {
         std::cout << "\nset fibers" << std::endl;
     }

//     solid.activeMaterial() -> setFiberVector( *solidFibers );

//     monodomain -> setupFibers();

     vectorPtr_Type gammaf( new vector_Type( ( monodomain -> globalSolution().at(3) ) -> map() ) );
     vectorPtr_Type gammas( new vector_Type( ( monodomain -> globalSolution().at(3) ) -> map() ) );
     vectorPtr_Type gamman( new vector_Type( ( monodomain -> globalSolution().at(3) ) -> map() ) );
     vectorPtr_Type solidGammaf;
     vectorPtr_Type emDisp;
     solidFESpacePtr_Type electroFiberFESpace;
     solidETFESpacePtr_Type electrodETFESpace;
     if(usingDifferentMeshes)
     {



    	 electroFiberFESpace.reset ( new solidFESpace_Type (monodomain -> localMeshPtr(), "P1", 3, comm) );
    	 electrodETFESpace.reset ( new solidETFESpace_Type (monodomain -> localMeshPtr(), & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), comm) );

         vectorPtr_Type electroFibers( new vector_Type( electroFiberFESpace -> map() ) );
//       HeartUtility::setupFibers(*electroFibers, fvec);
//       vectorPtr_Type fibersRotated( new vector_Type( monodomain -> feSpacePtr() -> map() ) );

         electroFiberFESpace -> interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( fibersDirection ), *electroFibers , 0);
         ElectrophysiologyUtility::normalize(*electroFibers);
//     HeartUtility::importFibers(electroFibers, parameterList.get ("fiber_file", ""), monodomain -> localMeshPtr() );
         monodomain -> setFiberPtr( electroFibers );
    	 emDisp.reset(  new vector_Type( electroFibers -> map() ) );
    	 solidGammaf.reset( new vector_Type( solidaFESpace -> map() ) );

     }
     else
     {
    	solidGammaf = gammaf;
    	 monodomain -> setFiberPtr( solidFibers );
    	 emDisp = solid.displacementPtr();
    	 electroFiberFESpace = dFESpace;
    	 electrodETFESpace = dETFESpace;
     }


     monodomain -> exportFiberDirection(problemFolder);
     boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporterSheets;
     exporterSheets.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, "SheetsDirection" ) );

       //      exporter->setPostDir ( "./" );
             exporterSheets -> setPostDir ( problemFolder );
       exporterSheets->setMeshProcId ( localSolidMesh, comm->MyPID() );
       exporterSheets->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "sheets", dFESpace, solidSheets, UInt (0) );
       exporterSheets->postProcess(0);
       exporterSheets->closeFile();
       //********************************************//
     // Create the global matrix: mass + stiffness in ELECTROPHYSIOLOGY //
     //********************************************//
     if ( comm->MyPID() == 0 )
     {
         cout << "\nSetup operators:  dt = " << monodomain -> timeStep() << "\n" ;
     }

     monodomain -> setDisplacementPtr( emDisp );
     monodomain -> setupLumpedMassMatrix();
     monodomain -> setupStiffnessMatrix();
     monodomain -> setupGlobalMatrix();

     if ( comm->MyPID() == 0 )
     {
         cout << "Done! \n" ;
     }

     //==================================================================//
     //==================================================================//
     //					SETUP Activation								//
     //==================================================================//
     //==================================================================//

     //   vectorPtr_Type gammaf( new vector_Type( monodomain -> globalSolution().at(3) -> map() ) );

     if(load4restart)
     {
     	std::string Gfilename = parameterList.get("Gfilename", "G");
     	std::string Gfieldname = parameterList.get("Gfieldname", "G");
     	ElectrophysiologyUtility::importScalarField(gammaf,Gfilename,Gfieldname,monodomain -> localMeshPtr() );
     }
     else *gammaf *= 0;
        *solidGammaf = parameterList.get("initial_gamma_f",0.0);


     if ( comm->MyPID() == 0 )
     {
         std::cout << "\nset gammaf and fibers" << std::endl;
     }
     solid.activeMaterial() -> setGammaf( *solidGammaf );

     vectorPtr_Type solidGammas( new vector_Type( solidGammaf -> map() ) );
     vectorPtr_Type solidGamman( new vector_Type( solidGammaf -> map() ) );
     Int gcase = parameterList.get ("case", 0);
     if ( comm->MyPID() == 0 )
     {
         std::cout << "\ncase: " << gcase << std::endl;
     }
     if(gcase == 1)
     {
		 Real gfactor = parameterList.get ("gfactor", 3.0);
		 *solidGamman = gfactor * *solidGammaf;
		 solid.activeMaterial() -> setGamman(*solidGamman);
		 *solidGammas = 1.0;
		 *solidGammas /= (1.0 + *solidGammaf);
		 *solidGammas /= (1.0 + *solidGamman);
		 *solidGammas -= 1.0;
		 solid.activeMaterial() -> setGammas(*solidGammas);

		 if(usingDifferentMeshes)
		 {
		 *gamman = gfactor * *gammaf;
		 *gammas = 1.0;
		 *gammas /= (1.0 + *gammaf);
		 *gammas /= (1.0 + *gamman);
		 *gammas -= 1.0;
		 }
		 else
		 {
			 gamman = solidGamman;
			 gammas = solidGammas;
		 }
     }
     else
     {
		 *solidGammas = 1.0;
		 *solidGammas /= (1.0 + *solidGammaf);
		 EpetraSqrt(*solidGammas);
		 *solidGammas -= 1.0;
		 solid.activeMaterial() -> setGamman(*solidGammas);
		 *solidGamman = *solidGammas;
		 solid.activeMaterial() -> setGammas(*solidGamman);

		 if(usingDifferentMeshes)
		 {
		 *gammas = 1.0;
		 *gammas /= (1.0 + *gammaf);
		 EpetraSqrt(*gammas);
		 *gammas -= 1.0;
 		 *gamman = *gammas;
		 }
		 else
		 {
			 gamman = solidGamman;
			 gammas = solidGammas;
		 }
     }

     //==================================================================//
     //==================================================================//
     //					SETUP INTERPOLATION								//
     //==================================================================//
     //==================================================================//

	 typedef RBFInterpolation<mesh_Type>           interpolation_Type;
	 typedef boost::shared_ptr<interpolation_Type> interpolationPtr_Type;

	 //Coarse To Fine ( C2F )
	 interpolationPtr_Type C2F;
	 //Fine To Coarse ( F2C )
	 interpolationPtr_Type F2C;
     if(usingDifferentMeshes)
     {
		    Teuchos::RCP< Teuchos::ParameterList > belosList = Teuchos::rcp ( new Teuchos::ParameterList );
		    belosList = Teuchos::getParametersFromXmlFile ( "ParamList.xml" );

		     if ( comm->MyPID() == 0 )
		     {
		         std::cout << "\nresetting" << std::endl;
		     }
		     std::string c2f = parameterList.get ("c2f", "RBFrescaledVectorial");
		     std::string f2c = parameterList.get ("f2c", "RBFrescaledScalar");

	     C2F.reset ( interpolation_Type::InterpolationFactory::instance().createObject ( c2f ) );
		 F2C.reset ( interpolation_Type::InterpolationFactory::instance().createObject ( f2c ) );
//		 C2F.reset ( interpolation_Type::InterpolationFactory::instance().createObject ( "RBFlocallyRescaledVectorial" ) );
//		 F2C.reset ( interpolation_Type::InterpolationFactory::instance().createObject ( "RBFlocallyRescaledScalar" ) );

			int nFlags = 1;
			std::vector<int> flags (nFlags);
			flags[0] = -1;

		     if ( comm->MyPID() == 0 )
		     {
		         std::cout << "\nintepolation setup" << std::endl;
		     }

			C2F->setup( fullSolidMesh , localSolidMesh, monodomain -> fullMeshPtr(), monodomain -> localMeshPtr(), flags);
			F2C->setup( monodomain -> fullMeshPtr(), monodomain -> localMeshPtr(), fullSolidMesh , localSolidMesh, flags);


		     if ( comm->MyPID() == 0 )
		     {
		         std::cout << "\nC2F: set Radius, ..." << std::endl;
		     }

			C2F->setRadius( 2.0 * (double) MeshUtility::MeshStatistics::computeSize (* (fullSolidMesh) ).maxH );
		     //C2F->setRadius( ( (double) MeshUtility::MeshStatistics::computeSize (* (monodomain -> fullMeshPtr())  ).maxH ) );

		     if ( comm->MyPID() == 0 )
		     {
		         std::cout << "\nC2F: set data..." << std::endl;
		     }
			C2F->setupRBFData ( solid.displacementPtr(), emDisp , dataFile, belosList);

		     if ( comm->MyPID() == 0 )
		     {
		         std::cout << "\nC2F: Build operator..." << std::endl;
		     }
		     if(c2f == "RBFvectorial") C2F->setBasis("TPS");
			C2F->buildOperators();

		     if ( comm->MyPID() == 0 )
		     {
		         std::cout << "\nF2C: set Radius, data, and build operator..." << std::endl;
		     }

			F2C->setRadius( (double) MeshUtility::MeshStatistics::computeSize (* (monodomain -> fullMeshPtr()) ).maxH );
			F2C->setupRBFData ( gammaf, solidGammaf , dataFile, belosList);
			F2C->buildOperators();

		     if ( comm->MyPID() == 0 )
		     {
		         std::cout << "\nC2F: interpolate..." << std::endl;
		     }

			C2F->interpolate();
			C2F->solution (emDisp);

		     if ( comm->MyPID() == 0 )
		     {
		         std::cout << "\nF2C: interpolate..." << std::endl;
		     }
			F2C->interpolate();
			F2C->solution (solidGammaf);

     }

     //     function_Type initialGuess = &d0;
//     vectorPtr_Type initd( new vector_Type( dFESpace -> map() ) );
//     dFESpace -> interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( initialGuess ), *initd , 0);





//     if ( comm->MyPID() == 0 )
//	  {
//		  std::cout << "\nnorm inf gammaf: " << solid.activeMaterial() -> gammaf() -> normInf() << std::endl;
//		  std::cout << "\nnorm inf fiber: " << solid.activeMaterial() -> fiberVector() -> normInf() << std::endl;
//	  }
     if ( comm->MyPID() == 0 )
     {
         std::cout << "\nbuild solid system" << std::endl;
     }

     solid.buildSystem(1.0);
     vectorPtr_Type rhs (new vector_Type (solid.displacement(), Unique) );
     vectorPtr_Type disp (new vector_Type (solid.displacement(), Unique) );
     vectorPtr_Type initialDisplacement (new vector_Type (solid.displacement(), Unique) );
     solid.initialize ( initialDisplacement );


     MPI_Barrier (MPI_COMM_WORLD);

     if ( comm->MyPID() == 0 )
     {
         std::cout << "\nsetup solid exporter" << std::endl;
     }

      boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;
      exporter.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, parameterList.get ("StructureOutputFile", "StructureOutput") ) );

      //      exporter->setPostDir ( "./" );
            exporter -> setPostDir ( problemFolder );
      exporter->setMeshProcId ( localSolidMesh, comm->MyPID() );

      vectorPtr_Type solidDisp ( new vector_Type (solid.displacement(), exporter->mapType() ) );
      exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", dFESpace, solidDisp, UInt (0) );
      exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "solid_gammaf", solidaFESpace, solidGammaf, UInt (0) );





      //================================================================//
      //================================================================//
      //					SETUP COUPLING SOLVER						//
      //																//
      //================================================================//
      //================================================================//
      ExporterHDF5< RegionMesh <LinearTetra> > expGammaf;
	expGammaf.setMeshProcId(monodomain -> localMeshPtr(), comm->MyPID());
	expGammaf.setPrefix(parameterList.get ("ActivationOutputFile", "ActivationOutput"));
	expGammaf.setPostDir ( problemFolder );

//  	expGammaf.addVariable(ExporterData<mesh_Type>::ScalarField, "gammaf",
//  			monodomain -> feSpacePtr(), gammaf, UInt(0));
//    expGammaf.postProcess(0.0);
//    Real min =  0.2;
//    Real max =  0.85;
//
//    Real beta = -0.3;
//
//    HeartUtility::rescaleVector(*gammaf, min, max, beta);


      matrixPtr_Type mass(new matrix_Type( monodomain -> massMatrixPtr() -> map() ) ) ;




     	boost::shared_ptr<FLRelationshipGamma> flg (new FLRelationshipGamma);
     	boost::shared_ptr<FLRelationship> fl (new FLRelationship);

     	boost::shared_ptr<HeavisideFct> H (new HeavisideFct);

     	boost::shared_ptr<Exp> EXP(new Exp);
     	boost::shared_ptr<Exp2> EXP2(new Exp2);
     	boost::shared_ptr<Psi4f> psi4f (new Psi4f);
     	boost::shared_ptr<ShowValue> sv(new ShowValue);

      MatrixSmall<3,3> Id;
      Id(0,0) = 1.; Id(0,1) = 0., Id(0,2) = 0.;
      Id(1,0) = 0.; Id(1,1) = 1., Id(1,2) = 0.;
      Id(2,0) = 0.; Id(2,1) = 0., Id(2,2) = 1.;
      VectorSmall<3> E1;
      E1(0) = 1.;
      E1(1) = 0.;
      E1(2) = 0.;
  	{
  		using namespace ExpressionAssembly;

  		integrate(elements(monodomain -> localMeshPtr() ), monodomain -> feSpacePtr() -> qr(), monodomain -> ETFESpacePtr(),
  				monodomain -> ETFESpacePtr(), phi_i * phi_j) >> mass;

  	}
  	mass -> globalAssemble();


  	vectorPtr_Type rhsActivation( new vector_Type( *gammaf ) );
  	*rhsActivation *= 0;




      if ( comm->MyPID() == 0 )
      {
          std::cout << "\nSolve system" << std::endl;
      }

      //==================================================================//
      //==================================================================//
      //					SETUP LINEAR SOLVER	    						//
      //						ACTIVATION									//
      //==================================================================//
      //==================================================================//

      if ( comm->MyPID() == 0 )
      {
          std::cout << "\nset up linear solver... Does it work???" << std::endl;
      }
  	typedef LinearSolver linearSolver_Type;
  	typedef boost::shared_ptr<LinearSolver> linearSolverPtr_Type;
	typedef LifeV::Preconditioner basePrec_Type;
	typedef boost::shared_ptr<basePrec_Type> basePrecPtr_Type;
	typedef LifeV::PreconditionerIfpack prec_Type;
	typedef boost::shared_ptr<prec_Type> precPtr_Type;


	prec_Type* precRawPtr;
	basePrecPtr_Type precPtr;
	precRawPtr = new prec_Type;
	precRawPtr->setDataFromGetPot(dataFile, "prec");
	precPtr.reset(precRawPtr);

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nprec done!!!!" << std::endl;
    }

	Teuchos::RCP < Teuchos::ParameterList > solverParamList = Teuchos::rcp(
			new Teuchos::ParameterList);

	std::string xmlpath = dataFile("electrophysiology/monodomain_xml_path",
			"./");
	std::string xmlfile = dataFile("electrophysiology/monodomain_xml_file",
			"MonodomainSolverParamList.xml");

	solverParamList = Teuchos::getParametersFromXmlFile(xmlpath + xmlfile);

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nreading file done!!!!" << std::endl;
    }

	linearSolver_Type linearSolver;
    linearSolver.setCommunicator ( comm );
    linearSolver.setParameters ( *solverParamList );
    linearSolver.setPreconditioner ( precPtr );
	linearSolver.setOperator( mass );
//	linearSolver.setOperator( monodomain -> massMatrixPtr() );

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nIt does!!!!" << std::endl;
    }











   	vectorPtr_Type tmpRhsActivation( new vector_Type ( rhsActivation -> map(), Repeated ) );
    solidFESpacePtr_Type emDispFESpace ( new solidFESpace_Type ( monodomain -> localMeshPtr(), "P1", 3, comm) );
  	expGammaf.addVariable(ExporterData<mesh_Type>::ScalarField, "gammaf",
  			monodomain -> feSpacePtr(), gammaf, UInt(0));
  	expGammaf.addVariable(ExporterData<mesh_Type>::VectorField, "interpolated displacement",
  			emDispFESpace, emDisp, UInt(0));
  	expGammaf.addVariable(ExporterData<mesh_Type>::ScalarField, "rhs",
  			monodomain -> feSpacePtr(), rhsActivation, UInt(0));


  	expGammaf.postProcess(0.0);

    //================================================================//
    //================================================================//
    //					SETUP gamma rescaling						  //
    //																  //
    //================================================================//
    //================================================================//
  	vectorPtr_Type rescaling(new vector_Type( solidGammaf -> map() ) );
  	function_Type resc = &rescalingGamma;
    solidaFESpace -> interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( resc ), *rescaling , 0);
  //  vectorPtr_Type emDisp0(new vector_Type( emDisp -> map() ) );




    //===========================================================
    //===========================================================
    // CREATE POSITION VECTOR
    //===========================================================
    //===========================================================
//    vectorPtr_Type referencePosition( new vector_Type( solidDisp -> map() ) );
    vectorPtr_Type referencePositionX( new vector_Type( monodomain -> potentialPtr() -> map() ) );
//    vectorPtr_Type position( new vector_Type( solidDisp -> map() ) );
//    vectorPtr_Type intergralR( new vector_Type( dETFESpace->map() ) );
//    vectorPtr_Type intergralR2( new vector_Type( monodomain -> potentialPtr() -> map() ) );
//    vectorPtr_Type oneVec( new vector_Type( dETFESpace->map() ) );
//    *oneVec = 1.0;
//    createPositionVector(*fullSolidMesh, *referencePosition);
    createPositionXVector(*fullSolidMesh, *referencePositionX);
    Real fluidVolume; //= ComputeVolume(*fullSolidMesh, *referencePosition, *solidDisp, 10, comm);
    fluidVolume =  ComputeVolume( monodomain -> localMeshPtr(), *referencePositionX, *solidDisp, monodomain->ETFESpacePtr(),dETFESpace,10,comm);
    //===========================================================
  	//===========================================================
  	//				Initializing solid
  	//===========================================================
  	//===========================================================
    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporterRamp;

    exporterRamp.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, "RampOutput" ) );
    exporterRamp -> setPostDir ( problemFolder );
    exporterRamp->setMeshProcId ( localSolidMesh, comm->MyPID() );
    exporterRamp->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "ramp displacement", dFESpace, solidDisp, UInt (0) );

    vector_Type endoVec (0.0 * (*solidDisp) , Repeated);
//    vector_Type pressure (aFESpace->map(), Repeated);
    Real pressure = parameterList.get("pressure", 0.0);
    vectorPtr_Type savePressure( new vector_Type( endoVec.map() )  );
    vectorPtr_Type saveVolume( new vector_Type( endoVec.map()) );
	*saveVolume = fluidVolume;
	   *savePressure = -pressure;
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "pressure", aFESpace, savePressure, UInt (0) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "volume", aFESpace, saveVolume, UInt (0) );

    boost::shared_ptr<BCVector> pEndo( new BCVector (endoVec, dFESpace -> dof().numTotalDof(), 1) );

	if(parameterList.get("debug", false)){

    solidBC -> handler() -> addBC("Endocardium", 10, Natural, Full, *pEndo, 3);
//    solidBC -> handler() -> addBC
	}
	if ( comm->MyPID() == 0 )
		{
			std::cout << "\nContractile fraction: " << solid.data() -> contractileFraction() << std::endl;
		}
//	solid.data() -> showMe();
//	solid.activeMaterial() -> showMyParameters();

    if(parameterList.get("pressure_ramp", false) == true){


    	Real ramp_dt = parameterList.get("ramp_timestep", 0.1);
		if ( comm->MyPID() == 0 )
		{
			std::cout << "\nSTARTING PRESSURE RAMP!\n" << ramp_dt << "\n" << std::endl;
		}
    	for(Real pseudot(0); pseudot < 1; ){
    		pseudot += ramp_dt;
    		if ( comm->MyPID() == 0 )
    		{
    			std::cout << "\nPRESSURE RAMP: " << pseudot;
    		}
    		solid.data() -> dataTime() -> setTime(pseudot);
    		if(parameterList.get("debug", false)){
    		endoVec = -(pressure * pseudot);
    		pEndo.reset( ( new BCVector (endoVec, dFESpace -> dof().numTotalDof(), 1) ) );
    	    solidBC -> handler() -> modifyBC(10, *pEndo);
    		}
//    		if ( comm->MyPID() == 0 )
//    		{
//    			std::cout << "\nnorm displacement: " << solid.displacement().norm2();
//    			std::cout << "\nnorm gammaf: " << solid.activeMaterial() -> gammaf() -> norm2();
//    			std::cout << "\nnorm gamman: " << solid.activeMaterial() -> gamman() -> norm2();
//    			std::cout << "\nnorm gammas: " << solid.activeMaterial() -> gammas() -> norm2();
//    		}
//    		solidBC -> handler() -> showMe();
    		solid.iterate ( solidBC -> handler() );

    	    exporterRamp->postProcess(pseudot);
    	    *solidDisp = solid.displacement();

		fluidVolume =  ComputeVolume( monodomain -> localMeshPtr(), *referencePositionX, *solidDisp, monodomain->ETFESpacePtr(),dETFESpace,10,comm);

		}
	      exporterRamp->closeFile();
    }
    *saveVolume = fluidVolume;
    exporter->postProcess ( 0 );


	if(usingDifferentMeshes)
	{
		if(solid.displacementPtr() -> normInf() > 0)
		{
			if ( comm->MyPID() == 0 )
			{
				std::cout << "\nINTERPOLATING FROM COARSE TO FINE!\n" << std::endl;
			}

			C2F -> updateRhs( solid.displacementPtr() );
			C2F -> interpolate();
			C2F -> solution( emDisp );
		}
		else *emDisp *= 0.0;
	}




	vectorPtr_Type emDisp0(new vector_Type( emDisp -> map() ) );
    *emDisp0 = *emDisp;
    cout << "\n\naddress emDisp: " << emDisp << ", emDisp0: " << emDisp0 << "\n\n";
    //return 0.0;

  	bool twoWayCoupling = parameterList.get("two_way", false);

	  if(twoWayCoupling)
	  {
			if ( comm->MyPID() == 0 )
			{
				std::cout << "\nREASSEMBLING STIFFNESS MATRIX FOR TOW WAY COUPLING!\n" << std::endl;
			}

			 monodomain -> setupStiffnessMatrix();
			 monodomain -> setupGlobalMatrix();
	  }


    //===========================================================
  	//===========================================================
  	//				TIME LOOP
  	//===========================================================
  	//===========================================================
      Real emdt = parameterList.get("emdt",1.0);
      int iter((emdt / monodomain -> timeStep()));
      int k(0);
      Real saveStep = parameterList.get("save_step",1.0);
      int saveIter((saveStep / monodomain -> timeStep()));
      Real meth = parameterList.get("meth",1.0);
      int subiter = parameterList.get("subiter",100);

        Real dt_min = 0.01;


		Real Ca_diastolic = dataFile( "solid/physics/Ca_diastolic", 0.02155 );

		Real p0 = 0;
		Real pnm1 = pressure;
		Real dp = pressure - p0;
		Real preload = pressure;
		
		fluidVolume =  ComputeVolume( monodomain -> localMeshPtr(), *referencePositionX, *solidDisp, monodomain->ETFESpacePtr(),dETFESpace,10,comm);
		Real V0 = fluidVolume;
		Real dV = fluidVolume - V0;
		Real dV_temporal(dV);
		Real dp_temporal(dp);

		bool isoPV = true;
		Real res;
		std::string outputPV="/"+problemFolder+"/pv.txt";
		const char * pvOutputFile = outputPV.c_str();
		std::ofstream pv(pvOutputFile);
		if(  comm->MyPID() == 0 ){
		std::string outputPV="/"+problemFolder+"/pv.txt";
		const char * pvOutputFile = outputPV.c_str();
		pv << "Pressure, Volume";
		pv << "\n" << pressure << ", " << fluidVolume;
		}

	Int pv_phase(0);



    vectorPtr_Type sx(new vector_Type( dFESpace -> map() ) );
    vectorPtr_Type sy(new vector_Type( dFESpace -> map() ) );
    vectorPtr_Type sz(new vector_Type( dFESpace -> map() ) );
    vectorPtr_Type I4(new vector_Type( monodomain -> feSpacePtr() -> map() ) );


    *sx = GradientRecovery::ZZGradient(dFESpace, *emDisp, 0);
    *sy = GradientRecovery::ZZGradient(dFESpace, *emDisp, 1);
    *sz = GradientRecovery::ZZGradient(dFESpace, *emDisp, 2);

    computeI4(*I4, *sx, *sy, *sz, *( monodomain -> fiberPtr() ) );


Real active_coefficient = dataFile( "solid/physics/active_coefficient", -2.5 );
Real viscosity = dataFile( "solid/physics/viscosity", 0.0005 );


     for( Real t(0.0); t< monodomain -> endTime(); )
	 {
		  t = t + monodomain -> timeStep();
		  k++;


//		    if ( comm->MyPID() == 0 )
//		    {
//		        std::cout << "\nSolve REACTION step with ROS3P!\n" << std::endl;
//		    }
		    LifeChrono timer;
		    timer.start();
//		    if(meth == 1.0) monodomain -> solveOneReactionStepROS3P(dt_min);
//		    else
//		    	{
		    	 for(int j(0); j<subiter; j++) monodomain -> solveOneReactionStepFE(subiter);
//		    	}

          timer.stop();
          for(int pid(0); pid < 4; pid++)
          {
				if ( comm->MyPID() == pid )
				{
						std::cout << "\nDone in " << timer.diff() << std::endl;
				}
          }

		if ( comm->MyPID() == 0 )
		{
			std::cout << "\nSolve DIFFUSION step with BE!\n" << std::endl;
		}
		timer.reset();
		timer.start();

	    (*monodomain -> rhsPtrUnique()) *= 0.0;
	      monodomain -> updateRhs();
          monodomain -> solveOneDiffusionStepBE();
  		timer.stop();
        for(int pid(0); pid < 4; pid++)
        {
				if ( comm->MyPID() == pid )
				{
						std::cout << "\nDone in " << timer.diff() << std::endl;
				}
        }
		timer.reset();

		  *tmpRhsActivation *= 0;
			if ( comm->MyPID() == 0 )
			{
				std::cout << "\nASSEMBLING ACTIVATION EQUATION!\n" << std::endl;
			}


			  if( monodomain -> globalSolution().at(3) -> normInf() >= Ca_diastolic
					  ||
				   emDisp -> normInf() >=  1e-7)
			  {

//					if ( comm->MyPID() == 0 )
//					{
//						std::cout << "\nI SHALL SOLVE THE ACTIVATION EQUATION!\n" << std::endl;
//					}

					{
						using namespace ExpressionAssembly;


						BOOST_AUTO_TPL(I,      value(Id) );
						BOOST_AUTO_TPL(Grad_u, grad( electrodETFESpace, *emDisp, 0) );
						BOOST_AUTO_TPL(F,      ( Grad_u + I ) );
						//BOOST_AUTO_TPL(FmT,    minusT(F) );
						BOOST_AUTO_TPL(J,       det(F) );
						BOOST_AUTO_TPL(Jm23,    pow(J, -2./3));
						BOOST_AUTO_TPL(I1,     dot(F, F));

						// Fibres
						BOOST_AUTO_TPL(f0,     value( electrodETFESpace, *( monodomain -> fiberPtr() ) ) );
						BOOST_AUTO_TPL(f,      F * f0 );
						BOOST_AUTO_TPL(I4f,    dot(f, f) );


					    BOOST_AUTO_TPL(s0,     value(electrodETFESpace, *(solid.activeMaterial() -> sheetVectorPtr() ) ) );
						BOOST_AUTO_TPL(I4s,    dot(F * s0, F * s0));

						BOOST_AUTO_TPL(I1iso,   Jm23 * I1);
						BOOST_AUTO_TPL(I4fiso,  Jm23 * I4f);
					    BOOST_AUTO_TPL(I4siso,  Jm23 * I4s);


						// Generalised invariants
					    BOOST_AUTO_TPL(gf,  value(aETFESpace, *gammaf));
				    	BOOST_AUTO_TPL(gs,  value(aETFESpace, *gammas));
				    	BOOST_AUTO_TPL(gn,  value(aETFESpace, *gamman));
//					    BOOST_AUTO_TPL(gs,  pow(gf + value(1.0), -0.5) + value(-1.0) );
//				    	BOOST_AUTO_TPL(gn,  gs);


				        BOOST_AUTO_TPL(dI1edI1,   value(1.0) / (   (gn + value(1.0))  *  (gn + value(1.0))  )    );
				        BOOST_AUTO_TPL(dI1edI4f,  value(1.0) / (   (gf + value(1.0))  *  (gf + value(1.0))  ) -  value(1.0) / (   (gn + value(1.0))  *  (gn + value(1.0))  ) );
				        BOOST_AUTO_TPL(dI1edI4s,  value(1.0) / (   (gs + value(1.0))  *  (gs + value(1.0))  ) -  value(1.0) / (   (gn + value(1.0))  *  (gn + value(1.0))  ) );
				        BOOST_AUTO_TPL(dI4fedI4f, value(1.0) / (   (gf + value(1.0)) *   (gf + value(1.0))  )    );

					    BOOST_AUTO_TPL(I1eiso,   dI1edI1 * I1iso + dI1edI4f * I4fiso + dI1edI4s * I4siso);
					    BOOST_AUTO_TPL(I4feiso,  dI4fedI4f * I4fiso);
						BOOST_AUTO_TPL(I4feisom1, ( I4feiso - value(1.0) ) );

						Real A = dataFile( "solid/physics/a", 4960 );
						Real B = dataFile( "solid/physics/b_activation", 0. );
						Real Af = dataFile( "solid/physics/af", 0. );
						Real Bf = dataFile( "solid/physics/bf", 0. );
						//cout << "\n\nparameters: a: " << A << ", af: " << Af << ", b: " << B << ", bf: " << Bf << "\n\n";
						BOOST_AUTO_TPL(a, value( A ) );
						BOOST_AUTO_TPL(b, value(  B ) );
						BOOST_AUTO_TPL(af, value( Af ) );
						BOOST_AUTO_TPL(bf, value(  Bf ) );
						BOOST_AUTO_TPL(psi_iso_e, a *  Jm23 * pow( eval(EXP, ( I1eiso + value(-3.0) ) ), B ) );
						BOOST_AUTO_TPL(psi_f_e,  value(2.0) * af * eval( H, I4feisom1 ) /* ( I4feiso + value(-1.0) )*/  * pow( eval(EXP2, ( I4feiso + value(-1.0) ) ), Bf ) );

						//initial
						BOOST_AUTO_TPL(Grad_u_i, grad( electrodETFESpace, *emDisp0, 0) );
						BOOST_AUTO_TPL(F_i,      ( Grad_u_i + I ) );
						BOOST_AUTO_TPL(J_i,       det(F_i) );
						BOOST_AUTO_TPL(Jm23_i,    pow(J_i, -2./3));
						BOOST_AUTO_TPL(I1_i,     dot(F_i, F_i));
						BOOST_AUTO_TPL(I1iso_i,   Jm23_i * I1_i);
						// Fibres
						BOOST_AUTO_TPL(f_i,      F_i * f0 );
						BOOST_AUTO_TPL(I4f_i,    dot(f_i, f_i) );
						BOOST_AUTO_TPL(I4fiso_i,  Jm23_i * I4f_i);
						BOOST_AUTO_TPL(I1eiso_i,   I1iso_i );
						BOOST_AUTO_TPL(I4feiso_i,  I4fiso_i);
						BOOST_AUTO_TPL(I4feisom1_i, ( I4feiso_i - value(1.0) ) );
						BOOST_AUTO_TPL(psi_iso_e_i, a *  Jm23 * pow( eval(EXP, ( I1eiso_i + value(-3.0) ) ), B ) );
						BOOST_AUTO_TPL(psi_f_e_i,  value(2.0) * af * eval( H, I4feisom1_i ) /* ( I4feiso + value(-1.0) )*/  * pow( eval(EXP2, ( I4feiso_i + value(-1.0) ) ), Bf ) );

						//BOOST_AUTO_TPL(dW01, value(-1.0) * a + eval( H, I4feisom1_i ) *  value(-2.0) * af * ( I4feiso + value(-1.0) ) ) ;
//						BOOST_AUTO_TPL(dW0, value(-1.0) * ( psi_iso_e_i + psi_f_e_i ) * I4fiso_i ) ;
	//					BOOST_AUTO_TPL(dW, value(-1.0) * ( psi_iso_e + psi_f_e ) * I4fiso * pow(gf + value(1.0), -3) );
						BOOST_AUTO_TPL(dW0, value(-2.0) * I4fiso_i) ;
						//BOOST_AUTO_TPL(dW, value(2.0) * I4fiso * ( value(3.0) * gf + value(-6.0) * gf * gf + value(10.0) * gf * gf * gf + value(-15.0) * gf * gf * gf * gf  + value(21.0) * gf * gf * gf * gf * gf) );
						BOOST_AUTO_TPL(dW, value(2.0) * I4f * ( value(3.0) * gf + value(-6.0) * gf * gf + value(10.0) * gf * gf * gf + value(-15.0) * gf * gf * gf * gf  + value(21.0) * gf * gf * gf * gf * gf) );


						BOOST_AUTO_TPL(Ca,    value( aETFESpace, *( monodomain -> globalSolution().at(3)  ) ) );
						BOOST_AUTO_TPL(Ca2, Ca * Ca );

						BOOST_AUTO_TPL(dCa, ( Ca - value(Ca_diastolic) ) );
					//    Real alpha1 = -2.5;
						//Real active_coefficient = dataFile( "solid/physics/active_coefficient", -2.5 );
					//	BOOST_AUTO_TPL(coeff, a * Jm23 * pow( eval(EXP, ( I1iso + value(-3.0) ) ), B )/*+ value(2.0) * af * eval( H, I4feisom1)*/ );
						//BOOST_AUTO_TPL(Pa, value(active_coefficient) * psi_iso_e  * eval(H, dCa) * eval(H, dCa) * eval(fl, I4fiso) + dW0 );
						//BOOST_AUTO_TPL(Pa, value(active_coefficient) * eval(H, dCa) * eval(H, dCa) * eval(fl, I4fiso) );
						//BOOST_AUTO_TPL(Pa, value(active_coefficient) * eval(H, dCa) * eval(H, dCa) * eval(fl, I4f) );
						//BOOST_AUTO_TPL(Pa, value(active_coefficient) * psi_iso_e  * eval(H, dCa) * eval(H, dCa) * eval(fl, I4fiso) + dW0 );

						//BOOST_AUTO_TPL(Pag, value(active_coefficient) * a * eval(H, dCa) * eval(H, dCa) * eval(flg, value( aETFESpace, *gammaf ) ) + dW0 );
					  //  Real delta = 0.001;
						//Real viscosity = dataFile( "solid/physics/viscosity", 0.0005 );
						//BOOST_AUTO_TPL(beta, value(viscosity ) /*/ coeff  /* Jm23 * pow( eval(EXP, ( I1iso + value(-3.0) ) ), -B )*/ );

//						BOOST_AUTO_TPL(dWs, a * pow(g, -1) );
						//BOOST_AUTO_TPL(gamma_dot, beta / ( Ca2 ) * ( Pa - dW )  );


						timer.start();
/*						integrate ( elements ( monodomain -> localMeshPtr() ),
								monodomain -> feSpacePtr() -> qr() ,
								monodomain -> ETFESpacePtr(),
								gamma_dot  * phi_i
						) >> tmpRhsActivation;*/

					}
					timer.stop();
//			        for(int pid(0); pid < 4; pid++)
//			        {
//							if ( comm->MyPID() == pid )
//							{
//									std::cout << "\nDone in " << timer.diff() << std::endl;
//							}
//			        }
					timer.reset();




				/*	*rhsActivation *= 0;
					*rhsActivation = ( *(mass) * ( *gammaf ) );
					*rhsActivation += ( ( monodomain -> timeStep() * *tmpRhsActivation ) );

					linearSolver.setRightHandSide(rhsActivation);
*/

			Int nLocalDof = (*I4).epetraVector().MyLength();
			//Int nComponentLocalDof = nLocalDof / 3;
			bool Pa0 = true;
			bool i40 = true;
			bool dCa0 = true;
			bool fl0 = true;

			for (int ik(0); ik < nLocalDof; ik++)
			{
				int iGID = (*I4).blockMap().GID (ik);
				Real Pa, dW;
				Real i4 = (*I4)[iGID];
				Real Ca = (*( monodomain -> globalSolution().at(3)  ) )[iGID];
				Real dCa = Ca - Ca_diastolic;
				Real fl = FL( i4 );
				if(dCa > 0.0) Pa = active_coefficient * dCa * dCa * fl;
				if(dCa <= 0.0) Pa = 0.0;

				Real g = (*gammaf)[iGID];
				Real g2 = g * g;
				Real g3 = g * g2;
				Real g4 = g * g3;
				Real g5 = g * g4;

				dW = 2.0 * i4 * ( 3.0 * g - 6.0 * g2 + 10.0 * g3 - 15.0 * g4  + 21.0 * g5 );
				Real grhs = viscosity * ( Pa - dW ) / Ca / Ca; 
				grhs *= monodomain -> timeStep();
				(*gammaf)[iGID] += grhs;			
				if(Pa > 0) Pa0 = false;
				if(dCa > 0) dCa0 = false;
				if(fl > 0) fl0 = false;
				if(i4 > 0.8 && i4 < 1.3) i40 = false;

			}

					if ( comm->MyPID() == 0 )
					{
						std::cout << "\nPa = 0 ?: " << Pa0 << ", i4 = 0 ?: " << i40 << ", dCa = 0 ?: " << dCa0 << ", fl = 0 ?: " << fl0 << "\n" << std::endl;
					}
					if ( comm->MyPID() == 0 )
					{
						std::cout << "\nSOLVING ACTIVATION EQUATION!\n" << std::endl;
					}


//					linearSolver.solve(gammaf);
			  }
			  else *gammaf *= 0.0;

//				if( gammaf -> maxValue() > 0.0)
//				{
//					std::cout << "\n*****************************************************";
//					std::cout << "\nCannot be positive: " ;
//					std::cout << "\n*****************************************************";
//					return 0;
//					int d = gammaf -> epetraVector().MyLength();
//					int size =  gammaf -> size();
//					for(int l(0); l < d; l++)
//					{
//						if ( comm->MyPID() == 0 )
//						{
//							std::cout << "\n*****************************************************";
//							std::cout << "\nChanging the gamma back to zero: " ;
//							std::cout << "\n*****************************************************";
//
//						}
//						int m1 = gammaf -> blockMap().GID(l);
//						//cout << m1 << "\t" << size << "\n";
//						if( (*gammaf)[m1] > 0)
//						{
//							(*gammaf)[m1] = 0.0;
//						}
//						std::cout << "\n";
//
//					}
//				}


			  if ( k % iter == 0)
			  {
					if(usingDifferentMeshes)
					{
						if(gammaf -> normInf() > 0)
						{
							if ( comm->MyPID() == 0 )
							{
								std::cout << "\nINTERPOLATING FROM FINE TO COARSE!\n" << std::endl;
							}

							F2C -> updateRhs( gammaf );
							F2C -> interpolate();
							F2C -> solution( solidGammaf );

						}
						else *solidGammaf *= 0.0;
					}


					solid.activeMaterial() -> setGammaf( *solidGammaf );

				     if(gcase == 1)
				     {
						 Real gfactor2 = parameterList.get ("gfactor", 3.0);
						 *solidGamman = gfactor2 * *solidGammaf;
						 solid.activeMaterial() -> setGamman(*solidGamman);
						 *solidGammas = 1.0;
						 *solidGammas /= (1.0 + *solidGammaf);
						 *solidGammas /= (1.0 + *solidGamman);
						 *solidGammas -= 1.0;
						 solid.activeMaterial() -> setGammas(*solidGammas);

						 if(usingDifferentMeshes)
						 {
						 *gamman = gfactor2 * *gammaf;
						 *gammas = 1.0;
						 *gammas /= (1.0 + *gammaf);
						 *gammas /= (1.0 + *gamman);
						 *gammas -= 1.0;
						 }
						 else
						 {
							 gamman = solidGamman;
							 gammas = solidGammas;
						 }
				     }
				     else
				     {
						 *solidGammas = 1.0;
						 *solidGammas /= (1.0 + *solidGammaf);
						 EpetraSqrt(*solidGammas);
						 *solidGammas -= 1.0;
						 solid.activeMaterial() -> setGamman(*solidGammas);
						 *solidGamman = *solidGammas;
						 solid.activeMaterial() -> setGammas(*solidGamman);

						 if(usingDifferentMeshes)
						 {
						 *gammas = 1.0;
						 *gammas /= (1.0 + *gammaf);
						 EpetraSqrt(*gammas);
						 *gammas -= 1.0;
				 		 *gamman = *gammas;
						 }
						 else
						 {
							 gamman = solidGamman;
							 gammas = solidGammas;
						 }
				     }

					if ( comm->MyPID() == 0 )
					{
						std::cout << "\nSOLVING STATIC MECHANICS!\n" << std::endl;
					}



					 res = 2 * fluidVolume;
					 //
					 p0 = pressure;
					 if( pv_phase == 0 && pressure > 100000) 
					{	
						isoPV = false;
						pv_phase = 1;
					}	
					if( pv_phase == 1 && dV_temporal > 0 )					 
					{	
						isoPV = true;
						pv_phase = 2;
						V0 = fluidVolume;
					}
					if( pv_phase == 2 && pressure<preload )					 
					{	
						isoPV = false;
						pv_phase = 3;
						pressure = preload;
					}
					if( pv_phase == 3 && dV_temporal < 0 && pressure > preload )					 
					{	
						isoPV = true;
						pv_phase = 0;
						V0 = fluidVolume;
					}	
						if ( comm->MyPID() == 0 )
						{
						std::cout  << "\nisoPV: " << isoPV << " \n";
						}



					if ( comm->MyPID() == 0 )
					{
						std::cout << "\n*****************************************************";
						std::cout << "\nWE ARE AT TIME: "<< t << ", PV PHASE: " << pv_phase;
						std::cout << "\n*****************************************************";

					}

					 if(isoPV)
					 {
						if ( comm->MyPID() == 0 )
						{
						std::cout  << "\nEntering in the loop (dV / V): " << dV / V0 << " \n";
						}
						int iter(0); 
						while( res > pvtol )
						{
							iter++;

					if ( comm->MyPID() == 0 )
					{
						std::cout << "\n*****************************************************";
						std::cout << "\nITERATION: "<< iter << ", TIME: " << t ;
						std::cout << "\n*****************************************************\n";

					}

							solid.iterate ( solidBC -> handler() );
							*solidDisp = solid.displacement();

							if ( comm->MyPID() == 0 )
							{
								std::cout << "\nisoPV: Initial Volume of Fluid: ... " <<  V0;
							}
							Real Vn = fluidVolume;
							fluidVolume =  ComputeVolume( monodomain -> localMeshPtr(), *referencePositionX, *solidDisp, monodomain->ETFESpacePtr(),dETFESpace,10,comm);
							//Real newFluidVolume = 1;//ComputeVolume(*fullSolidMesh, *referencePosition, *solidDisp, 10, comm);
						
							dV = fluidVolume - V0;
							//fluidVolume = newFluidVolume;
							if ( comm->MyPID() == 0 )
							{
								std::cout  << "\nVariation %: " << dV / V0 << " \n";
							}


							if ( comm->MyPID() == 0 )
							{
								std::cout << "\nComputing pressure: ... ";
							}
							Real newPressure = evaluatePressure(fluidVolume, dV, pressure, dp_temporal, dV, Cp, pv_phase);
							dp = newPressure - pressure;
							pressure = newPressure;
							if ( comm->MyPID() == 0 )
							{
								std::cout <<  newPressure <<", Variation %: " << dp / pressure << " \n";
							}
							endoVec = -newPressure;
							pEndo.reset( ( new BCVector (endoVec, dFESpace -> dof().numTotalDof(), 1) ) );
							solidBC -> handler() -> modifyBC(10, *pEndo);

							res = std::abs(dV) / V0;
							if ( comm->MyPID() == 0 )
							{
								std::cout <<  "\nresidual: " << res <<"\n";
							}

							if( pv_phase == 0 && pressure > 100000) break;
							if( pv_phase == 2 && pressure < preload) break;
						}
						dV_temporal = fluidVolume  - V0;
						dp_temporal = pressure - p0;
					 }
					 else
					 {
							Real Vn = fluidVolume;
							solid.iterate ( solidBC -> handler() );
							*solidDisp = solid.displacement();

							if ( comm->MyPID() == 0 )
							{
								std::cout << "\nno isoPV: Computing Volume of Fluid: ... ";
							}
							fluidVolume =  ComputeVolume( monodomain -> localMeshPtr(), *referencePositionX, *solidDisp, monodomain->ETFESpacePtr(),dETFESpace,10,comm);//Real newFluidVolume = 1.0;//ComputeVolume(*fullSolidMesh, *referencePosition, *solidDisp, 10, comm);
							dV = fluidVolume - Vn;
							//fluidVolume = newFluidVolume;
							if ( comm->MyPID() == 0 )
							{
								std::cout <<  fluidVolume <<", Variation %: " << dV / fluidVolume << " \n";
							}


							if ( comm->MyPID() == 0 )
							{
								std::cout << "\nComputing pressure: ... ";
							}
							Real newPressure = evaluatePressure(fluidVolume, dV, pressure, dp_temporal, dV_temporal, Cp, pv_phase, t);
							dp = newPressure - pressure;
							pressure = newPressure;
							if ( comm->MyPID() == 0 )
							{
								std::cout <<  newPressure <<", Variation %: " << dp / pressure << " \n";
							}
							endoVec = -newPressure;
							pEndo.reset( ( new BCVector (endoVec, dFESpace -> dof().numTotalDof(), 1) ) );
							solidBC -> handler() -> modifyBC(10, *pEndo);

							dV_temporal = fluidVolume  - Vn;
							dp_temporal = pressure - p0;
					 }

					    *savePressure = -pressure;
						*saveVolume = fluidVolume;

					//        timeAdvance->shiftRight ( solid.displacement() );


					/*if(parameterList.get("time_prestretch",false))
					{
						if( monodomain -> globalSolution().at(3)-> maxValue() < Ca_diastolic)
						{
							(*emDisp0) = (*emDisp);
						}
						else if( monodomain -> globalSolution().at(3)-> minValue() < Ca_diastolic)
						{
							int d = monodomain -> globalSolution().at(3) -> epetraVector().MyLength();
//							int size =  monodomain -> globalSolution().at(3) -> size();
							for(int l(0); l < d; l++)
							{
								if ( comm->MyPID() == 0 )
								{
									std::cout << "\n*****************************************************";
									std::cout << "\nChanging the initial displacement: " ;
									std::cout << "\n*****************************************************";

								}
								int m1 = monodomain -> globalSolution().at(3) -> blockMap().GID(l);
								//cout << m1 << "\t" << size << "\n";
								if( (*(monodomain -> globalSolution().at(3)))[m1] <= Ca_diastolic)
								{
									std::cout << "\nchanging at: " << m1 ;
									int m2 = solid.displacementPtr() -> blockMap().GID(l + d);
									int m3 = solid.displacementPtr() -> blockMap().GID(l + 2 * d);
									std::cout << "\n" << (*emDisp0)[m1] << " has become: ";
									(*emDisp0)[m1] = (*emDisp)[m1];
									std::cout << (*emDisp0)[m1] << "\n";
									(*emDisp0)[m2] = (*emDisp)[m2];
									(*emDisp0)[m3] = (*emDisp)[m3];
								}
								std::cout << "\n";

							}
						}
					}*/
					if(usingDifferentMeshes)
					{
						if(solid.displacementPtr() -> normInf() > 0)
						{
							if ( comm->MyPID() == 0 )
							{
								std::cout << "\nINTERPOLATING FROM COARSE TO FINE!\n" << std::endl;
							}

							C2F -> updateRhs( solid.displacementPtr() );
							C2F -> interpolate();
							C2F -> solution( emDisp );
						}
						else *emDisp *= 0.0;
					}



					  if(twoWayCoupling)
					  {
							if ( comm->MyPID() == 0 )
							{
								std::cout << "\nREASSEMBLING STIFFNESS MATRIX FOR TWO WAY COUPLING!\n" << std::endl;
							}

							 monodomain -> setupStiffnessMatrix();
							 monodomain -> setupGlobalMatrix();

					  }
					    *sx *= 0.0;
    					    *sy *= 0.0;
					    *sz *= 0.0;
					    *sx = GradientRecovery::ZZGradient(dFESpace, *emDisp, 0);
    					    *sy = GradientRecovery::ZZGradient(dFESpace, *emDisp, 1);
					    *sz = GradientRecovery::ZZGradient(dFESpace, *emDisp, 2);

					    computeI4(*I4, *sx, *sy, *sz, *( monodomain -> fiberPtr() ) );
			  }

							if ( comm->MyPID() == 0 )
							{
			pv << "\n" << pressure << ", " << fluidVolume;
	}
		  //cout << "\n\n save every " << saveIter << "iteration\n";
		  if ( k % saveIter == 0)
		  {
			  monodomain -> exportSolution(expElectro, t);
			  expGammaf.postProcess(t);
			  exporter->postProcess ( t );
		  }

		  //emDisp -> spy("interpolatedDisplacement");
		  //solid.displacementPtr() -> spy("displacement");

//	      if ( comm->MyPID() == 0 )
//	      {
//	    	  std::cout << "\n===================="
//	    	  std::cout << "\nWE ARE AT TIME: " << t ;
//	    	  std::cout << "\n===================="
//
//	      }


      }
	if ( comm->MyPID() == 0 )
	{
			pv.close();
	}
	
      expElectro.closeFile();
      expGammaf.closeFile();

      exporter -> closeFile();



    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nActive strain example: Passed!" << std::endl;
    }


    bool save4restart = parameterList.get("save4restart", true);
    if(save4restart)
    {
		ExporterHDF5< RegionMesh <LinearTetra> > restartV0;
		ExporterHDF5< RegionMesh <LinearTetra> > restartV1;
		ExporterHDF5< RegionMesh <LinearTetra> > restartV2;
		ExporterHDF5< RegionMesh <LinearTetra> > restartV3;
		ExporterHDF5< RegionMesh <LinearTetra> > restartG;
		ExporterHDF5< RegionMesh <LinearTetra> > restartD;

		restartV0.setMeshProcId(monodomain -> localMeshPtr(), comm->MyPID());
		restartV1.setMeshProcId(monodomain -> localMeshPtr(), comm->MyPID());
		restartV2.setMeshProcId(monodomain -> localMeshPtr(), comm->MyPID());
		restartV3.setMeshProcId(monodomain -> localMeshPtr(), comm->MyPID());
		restartG.setMeshProcId(monodomain -> localMeshPtr(), comm->MyPID());
		restartD.setMeshProcId(localSolidMesh, comm->MyPID());

		restartV0.setPrefix(parameterList.get ("variable0", "V0"));
		restartV1.setPrefix(parameterList.get ("variable1", "V1"));
		restartV2.setPrefix(parameterList.get ("variable2", "V2"));
		restartV3.setPrefix(parameterList.get ("variable3", "V3"));
		restartG.setPrefix(parameterList.get ("gamma_exporter", "A"));
		restartD.setPrefix(parameterList.get ("displacement_exporter", "D"));

		restartV0.setPostDir ( "./save_restart/" );
		restartV1.setPostDir ( "./save_restart/" );
		restartV2.setPostDir ( "./save_restart/" );
		restartV3.setPostDir ( "./save_restart/" );
		restartG.setPostDir ( "./save_restart/" );
		restartD.setPostDir ( "./save_restart/" );

		restartV0.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "V0", monodomain -> feSpacePtr(), monodomain -> globalSolution().at(0), UInt (0) );
		restartV1.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "V1", monodomain -> feSpacePtr(), monodomain -> globalSolution().at(1), UInt (0) );
		restartV2.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "V2", monodomain -> feSpacePtr(), monodomain -> globalSolution().at(2), UInt (0) );
		restartV3.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "V3", monodomain -> feSpacePtr(), monodomain -> globalSolution().at(3), UInt (0) );
		restartG.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "G", monodomain -> feSpacePtr(), gammaf, UInt (0) );
		restartD.addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "D", dFESpace, solidDisp, UInt (0) );

		restartV0.closeFile();
		restartV1.closeFile();
		restartV2.closeFile();
		restartV3.closeFile();
		restartG.closeFile();
		restartD.closeFile();
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}


