//
//  HeartSolver.cpp
//  LifeV
//
//  Created by Thomas Kummer on 03.05.16.
//  Copyright © 2016 Thomas Kummer. All rights reserved.
//

#ifndef _HEARTDATA_H_
#define _HEARTDATA_H_


#include <stdio.h>
#include <lifev/em/solver/EMSolver.hpp>
#include <lifev/em/solver/circulation/Circulation.hpp>

#include <lifev/em/solver/HeartData.hpp>


namespace LifeV
{

#define PI 3.14159265359

    
template <class EmSolver>
class HeartSolver {
   
public:
    
    typedef RegionMesh<LinearTetra>                         mesh_Type;
    typedef boost::shared_ptr<mesh_Type>                    meshPtr_Typ;
    typedef ExporterHDF5<mesh_Type>                         exporter_Type;
    typedef boost::shared_ptr<exporter_Type>                exporterPtr_Type;
    
    
    HeartSolver(EmSolver& emSolver,  Circulation& circulationSolver) :
        M_emSolver          (emSolver),
        M_circulationSolver (circulationSolver),
        M_heartData         (HeartData())
    {}
    
    virtual ~HeartSolver() {}

    EmSolver& emSolver()
    {
        return M_emSolver;
    }
    
    Circulation& circulation()
    {
        return M_heartData;
    }
    
    const HeartData& data() const
    {
        return M_heartData;
    }
    
    void setup(const GetPot& datafile)
    {
        M_heartData.setup(datafile);
    }
    
    template <class lambda>
    void preload(const lambda& modifyFeBC, const std::vector<Real>& bcValues)
    {
        M_emSolver.structuralOperatorPtr() -> data() -> dataTime() -> setTime(0.0);
        
        auto preloadPressure = [] (std::vector<double> p, const int& step, const int& steps)
        {
            for (auto& i : p) {i *= double(step) / double(steps);}
            return p;
        };
        
        LifeChrono chronoSave;
        chronoSave.start();
        
        M_emSolver.saveSolution (-1.0);
        
        if ( 0 == M_emSolver.comm()->MyPID() )
        {
            std::cout << "\n*****************************************************************";
            std::cout << "\nData stored in " << chronoSave.diff() << " s";
            std::cout << "\n*****************************************************************\n";
        }
        
        LifeChrono chronoPreload;
        chronoPreload.start();
        
        for (int i (1); i <= data().preloadSteps(); i++)
        {
            if ( 0 == M_emSolver.comm()->MyPID() )
            {
                std::cout << "\n*****************************************************************";
                std::cout << "\nPreload step: " << i << " / " << data().preloadSteps();
                std::cout << "\n*****************************************************************\n";
            }
            
            // Update pressure b.c.
            modifyFeBC(preloadPressure(bcValues, i, data().preloadSteps() ));
            
            // Solve mechanics
            M_emSolver.bcInterfacePtr() -> updatePhysicalSolverVariables();
            M_emSolver.solveMechanics();
            
            // Safe preload steps
            if ( data().safePreload() ) M_emSolver.saveSolution (i-1);
        }
        
        if ( 0 == M_emSolver.comm()->MyPID() )
        {
            std::cout << "\n*****************************************************************";
            std::cout << "\nPreload done in: " << chronoPreload.diff();
            std::cout << "\n*****************************************************************\n";
        }

    }
    
    
    void restart(std::string& restartInput, const GetPot& command_line, Real& t)
    {
        const std::string restartDir = ""; //command_line.follow (problemFolder.c_str(), 2, "-rd", "--restartDir");
        
        Real dtExport = 10.;
        
        // Set time variable
        const unsigned int restartInputStr = std::stoi(restartInput);
        const unsigned int nIter = (restartInputStr - 1) * dtExport / data().dt_mechanics();
        t = nIter * data().dt_mechanics();
        
        // Set time exporter time index
        M_emSolver.setTimeIndex(restartInputStr + 1);
        //solver.importHdf5();

        // Load restart solutions from output files
        std::string polynomialDegree = data().elementOrder();

        ElectrophysiologyUtility::importVectorField ( M_emSolver.structuralOperatorPtr() -> displacementPtr(), "MechanicalSolution" , "displacement", M_emSolver.localMeshPtr(), restartDir, polynomialDegree, restartInput );
        ElectrophysiologyUtility::importScalarField (M_emSolver.activationModelPtr() -> fiberActivationPtr(), "ActivationSolution" , "Activation", M_emSolver.localMeshPtr(), restartDir, polynomialDegree, restartInput );
        //ElectrophysiologyUtility::importScalarField (solver.activationTimePtr(), "ActivationTimeSolution" , "Activation Time", solver.localMeshPtr(), restartDir, polynomialDegree, restartInput );
        
        for ( unsigned int i = 0; i < M_emSolver.electroSolverPtr()->globalSolution().size() ; ++i )
        {
            ElectrophysiologyUtility::importScalarField (M_emSolver.electroSolverPtr()->globalSolution().at(i), "ElectroSolution" , ("Variable" + std::to_string(i)), M_emSolver.localMeshPtr(), restartDir, polynomialDegree, restartInput );
        }

        if ( 0 == M_emSolver.comm()->MyPID() )
        {
            std::cout << "\nLoad from restart: " << restartInput << ",  nIterCirculation = " << nIter << ",  time = " << t << std::endl;
        }
        
        M_circulationSolver.restartFromFile ( restartDir + "solution.dat" , nIter );
    }
    
    
    static Real Iapp (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
    {
        bool coords ( Y < -7. );
        //bool coords ( Y > 4. ); //( Y > 1.5 && Y < 3 );
        bool time ( fmod(t, 800.) < 4 && fmod(t, 800.) > 2);
        return ( coords && time ? 30 : 0 );
    }
    
    Real sinSquared (const Real& time, const Real& Tmax, const Real& tmax, const Real& tduration) const
    {
        Real timeInPeriod = fmod(time-tmax+0.5*tduration, 800.);
        bool inPeriod ( timeInPeriod < tduration && timeInPeriod > 0);
        Real sinusSquared = std::pow( std::sin(timeInPeriod * PI / tduration) , 2 ) * Tmax;
        return ( inPeriod ? sinusSquared : 0 );
    }
    
    
    template<class bcVectorType>
    void extrapolate4thOrderAdamBashforth(bcVectorType& bcValues, bcVectorType& bcValuesPre, const Real& dpMax)
    {
        VectorSmall<4> ABcoef;
        ABcoef (0) = 55/24; ABcoef (1) = -59/24; ABcoef (2) = 37/24; ABcoef (3) = -3/8;
        
        for ( unsigned int i = ABcoef.size() - 1; i > 0; --i )
        {
            m_ABdplv(i) = m_ABdplv(i-1);
            m_ABdprv(i) = m_ABdprv(i-1);
        }
        
        m_ABdplv(0) = bcValues[0] - bcValuesPre[0];
        m_ABdprv(0) = bcValues[1] - bcValuesPre[1];
        
        bcValuesPre = bcValues;
        
        bcValues[0] += std::min( std::max( ABcoef.dot( m_ABdplv ) , - dpMax ) , dpMax );
        bcValues[1] += std::min( std::max( ABcoef.dot( m_ABdprv ) , - dpMax ) , dpMax );
    }
    
    
    void setupExporter(std::string problemFolder = "./", std::string outputFileName = "humanHeartSolution")
    {
        m_exporter.reset (new exporter_Type());
        setupExporter<mesh_Type>(*m_exporter, M_emSolver.localMeshPtr(), M_emSolver.comm(), outputFileName, problemFolder);

        m_exporter->addVariable (    ExporterData<RegionMesh<LinearTetra> >::VectorField,
                                     "Displacement",
                                     M_emSolver.structuralOperatorPtr()->dispFESpacePtr(),
                                     M_emSolver.structuralOperatorPtr()->displacementPtr(),
                                     UInt (0) );
        
        m_exporter->addVariable (    ExporterData<RegionMesh<LinearTetra> >::ScalarField,
                                     "Von Mises stress",
                                     M_emSolver.electroSolverPtr()->feSpacePtr(),
                                     M_emSolver.tensionEstimator().vonMisesStressPtr(),
                                     UInt (0) );
        
        m_exporter->addVariable (    ExporterData<RegionMesh<LinearTetra> >::VectorField,
                                     "Principal stress",
                                     M_emSolver.structuralOperatorPtr()->dispFESpacePtr(),
                                     M_emSolver.tensionEstimator().principalStressesPtr(),
                                     UInt (0) );
        
        m_exporter->addVariable (    ExporterData<RegionMesh<LinearTetra> >::VectorField,
                                     "X stress",
                                     M_emSolver.structuralOperatorPtr()->dispFESpacePtr(),
                                     M_emSolver.tensionEstimator().sigmaXPtr(),
                                     UInt (0) );
    
        m_exporter->addVariable (    ExporterData<RegionMesh<LinearTetra> >::VectorField,
                                     "Y stress",
                                     M_emSolver.structuralOperatorPtr()->dispFESpacePtr(),
                                     M_emSolver.tensionEstimator().sigmaYPtr(),
                                     UInt (0) );
    
        m_exporter->addVariable (    ExporterData<RegionMesh<LinearTetra> >::VectorField,
                                     "Z stress",
                                     M_emSolver.structuralOperatorPtr()->dispFESpacePtr(),
                                     M_emSolver.tensionEstimator().sigmaZPtr(),
                                     UInt (0) );
        
        m_exporter->addVariable (    ExporterData<RegionMesh<LinearTetra> >::VectorField,
                                     "Patch displacement",
                                     M_emSolver.structuralOperatorPtr()->dispFESpacePtr(),
                                     patchDisplacementSumPtr(),
                                     UInt (0) );
        
        m_exporter->addVariable (    ExporterData<RegionMesh<LinearTetra> >::VectorField,
                                     "Fibers",
                                     M_emSolver.structuralOperatorPtr()->dispFESpacePtr(),
                                     M_emSolver.structuralOperatorPtr()->fPtr(),
                                     // M_emSolver.structuralOperatorPtr()->EMMaterial()->fiberVectorPtr(),
                                     UInt (0) );
        
        m_exporter->addVariable (    ExporterData<RegionMesh<LinearTetra> >::VectorField,
                                     "Sheets",
                                     M_emSolver.structuralOperatorPtr()->dispFESpacePtr(),
                                     M_emSolver.structuralOperatorPtr()->sPtr(),
                                     //M_emSolver.structuralOperatorPtr()->EMMaterial()->sheetVectorPtr(),
                                     UInt (0) );
        
        m_exporter->addVariable (    ExporterData<RegionMesh<LinearTetra> >::ScalarField,
                                     "Activation",
                                     M_emSolver.electroSolverPtr()->feSpacePtr(),
                                     M_emSolver.activationModelPtr()->fiberActivationPtr(),
                                     UInt (0) );
            
        m_exporter -> addVariable (  ExporterData<RegionMesh<LinearTetra> >::ScalarField,
                                     "Activation time",
                                     M_emSolver.electroSolverPtr() -> feSpacePtr(),
                                     M_emSolver.activationTimePtr(),
                                     UInt (0) );
        
        for (int i = 0; i < M_emSolver.electroSolverPtr()->ionicModelPtr()->Size(); ++i)
        {
            std::string variableName = "Ionic Variable " + boost::lexical_cast<std::string> (i);
            m_exporter->addVariable (ExporterData<mesh_Type>::ScalarField,
                                     variableName,
                                     M_emSolver.electroSolverPtr()->feSpacePtr(),
                                     M_emSolver.electroSolverPtr()->globalSolution().at(i),
                                     UInt (0) );
        }

        
    }

    
    template<class Mesh>
    void setupExporter (ExporterHDF5<Mesh>& exporter,
                        boost::shared_ptr<Mesh> localMeshPtr,
                        boost::shared_ptr<Epetra_Comm> commPtr,
                        std::string fileName,
                        std::string folder)
    {
        exporter.setMeshProcId (localMeshPtr, commPtr->MyPID() );
        exporter.setPrefix (fileName);
        exporter.exportPID (localMeshPtr, commPtr);
        exporter.setPostDir (folder);
    }
    
    
    void postProcess(const Real& time)
    {
        std::cout << patchDisplacementSumPtr()->size();

        
        // Compute Von Mises stress, principal stresses and Cauchy stresses
        M_emSolver.tensionEstimator().setDisplacement ( M_emSolver.structuralOperatorPtr()->displacement() );
        M_emSolver.tensionEstimator().analyzeTensionsRecoveryCauchyStresses();
        M_emSolver.tensionEstimator().analyzeTensionsRecoveryVonMisesStress();
        //M_emSolver.tensionEstimator().analyzeTensionsRecoveryEigenvalues();

        // Compute deformed fiber direction
        M_emSolver.computeDeformedFiberDirection (M_emSolver.structuralOperatorPtr()->f(), *M_emSolver.structuralOperatorPtr()->EMMaterial()->fiberVectorPtr(), *M_emSolver.structuralOperatorPtr()->displacementPtr(), M_emSolver.structuralOperatorPtr()->dispFESpacePtr());

        // Compute deformed sheet direction
        M_emSolver.computeDeformedFiberDirection (M_emSolver.structuralOperatorPtr()->s(), *M_emSolver.structuralOperatorPtr()->EMMaterial()->sheetVectorPtr(), *M_emSolver.structuralOperatorPtr()->displacementPtr(), M_emSolver.structuralOperatorPtr()->dispFESpacePtr());
        
        // Write on hdf5 output file
        m_exporter->postProcess(time);
    }
    
    
    exporterPtr_Type exporter()
    {
        return m_exporter;
    }
    
    Real
    externalPower ( const VectorEpetra& dispCurrent,
                   const VectorEpetra& dispPrevious,
                   const boost::shared_ptr<ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 > > dispETFESpace,
                   Real pressure,
                   Real dt,
                   const unsigned int bdFlag) const
    {
        VectorEpetra traction ( dispCurrent.map() );
        VectorEpetra velocity ( (dispCurrent - dispPrevious) / dt );
        
        MatrixSmall<3,3> Id;
        Id(0,0) = 1.; Id(0,1) = 0., Id(0,2) = 0.;
        Id(1,0) = 0.; Id(1,1) = 1., Id(1,2) = 0.;
        Id(2,0) = 0.; Id(2,1) = 0., Id(2,2) = 1.;
        
        {
            using namespace ExpressionAssembly;
            
            auto I = value(Id);
            auto Grad_u = grad( dispETFESpace, dispCurrent, 0);
            auto F =  Grad_u + I;
            auto FmT = minusT(F);
            auto J = det(F);
            auto p = value(pressure);
            
            QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria7pt) );
            
            integrate ( boundary ( dispETFESpace->mesh(), bdFlag),
                       myBDQR,
                       dispETFESpace,
                       value(-1.0) * p * J * dot( FmT * Nface,  phi_i)
                       //p * J * dot( FmT * Nface,  phi_i)
                       //value(-1.0) * J * dot (vE1, FmT * Nface) * phi_i) >> intergral
                       ) >> traction;
            
            traction.globalAssemble();
        }
        
        return traction.dot(velocity);
    }
    
    void setPatchDisplacementSumPtr(vectorPtr_Type patchDisplacementSumPtr)
    {
        m_patchDisplacementSumPtr = patchDisplacementSumPtr;
    }
    
    boost::shared_ptr<VectorEpetra> patchDisplacementSumPtr()
    {
        return m_patchDisplacementSumPtr;
    }
    
    
protected:
    
    EmSolver& M_emSolver;
    Circulation& M_circulationSolver;
    
    HeartData M_heartData;
    
    exporterPtr_Type m_exporter;
    
    VectorSmall<2> M_pressure;
    VectorSmall<2> M_volume;

    VectorSmall<4> m_ABdplv, m_ABdprv;

    boost::shared_ptr<VectorEpetra> m_patchDisplacementSumPtr;
    
    
    
    
    
    std::string pipeToString ( const char* command )
    {
        FILE* file = popen( command, "r" ) ;
        std::ostringstream stm ;
        char line[6] ;
        fgets( line, 6, file );
        stm << line;
        pclose(file) ;
        return stm.str() ;
    };
    


    Real patchDispFunNormal (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& i)
    {
        return (-0.000 - 0.00005*t);// sinSquared(t, 0.1, 50, 100)); // -0.001;// (t * 1e-5);
    }
    
    Real patchFunction (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
    {
        Real disp = std::pow( std::sin(fmod(t, 800.) * 3.14159265359/300) , 2 )*15;
        return disp;
    }
    
    Real potentialMultiplyerFcn (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
    {
        bool time ( fmod(t, 800.) < 4 && fmod(t, 800.) > 2);
        return 1.4 * time; // ( Y < 2.5 && Y > 0.5 ? 1.0 : 0.0 );
    }
    
    Real patchDispFun1 (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& i)
    {
        switch (i)
        {
            case 0:
                return (t);
                break;
            case 1:
                return 0;
                break;
            case 2:
                return (t);
                break;
            default:
                ERROR_MSG ("This entry is not allowed");
                return 0.;
                break;
        }
    }
    
    Real patchDispFun2 (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& i)
    {
        switch (i)
        {
            case 0:
                return (-t);
                break;
            case 1:
                return 0;
                break;
            case 2:
                return (-t);
                break;
            default:
                ERROR_MSG ("This entry is not allowed");
                return 0.;
                break;
        }
    }
    
    Real normalDirection ( const Real& /*t*/, const Real& x , const Real& y, const Real& z, const ID& i)
    {
        Real nnorm = std::sqrt(x * x + y * y + z * z);
        switch (i)
        {
            case 0:
                return 1/std::sqrt(2);
                break;
            case 1:
                return 0;
                break;
            case 2:
                return 1/std::sqrt(2);
                break;
            default:
                ERROR_MSG ("This entry is not allowed");
                return 0.;
                break;
        }
    }

    
};

}

#endif
