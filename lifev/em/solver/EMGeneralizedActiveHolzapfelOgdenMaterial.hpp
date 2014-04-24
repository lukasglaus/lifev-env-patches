//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/**
 *  @file
 *  @brief Holzapfel-Ogden material definition with generalized activation.
 *
 *  @date   05-08-2013
 *  @author Simone Rossi <simone.rossi@epfl.ch>
**/

#ifndef _EMGENERALIZEDACTIVEHOLZAPFELOGDENMATERIAL_H_
#define _EMGENERALIZEDACTIVEHOLZAPFELOGDENMATERIAL_H_

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <lifev/em/solver/EMActiveStructuralConstitutiveLaw.hpp>
#include <lifev/electrophysiology/util/HeartUtility.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <boost/typeof/typeof.hpp>

namespace
{

class PositivePart
{
public:
    typedef LifeV::Real return_Type;

    return_Type operator() (const LifeV::Real& p)
    {
        return (p > 0.0) ? p : 0.0;
    }

    PositivePart() {}
    PositivePart (const PositivePart&) {}
    ~PositivePart() {}
};

// Strain-energy density function
// ------------------------------
namespace StrainEnergyGHO
{



class GAHOShowValue
{
public:
    typedef LifeV::Real return_Type;

    return_Type operator() (const LifeV::Real& a)
    {
        std::cout.precision (15);
        std::cout << "value is: " << a << " \n";
        return 1.0;
    }

    return_Type operator() (const LifeV::VectorSmall<3>& a)
    {
        std::cout << "value is: " << a[0]  << " \n" << a[1]  << " \n" << a[2]  << " \n";
        return 1.0;
    }

    return_Type operator() (const LifeV::MatrixSmall<3, 3>& a)
    {
        std::cout << "value is\n";
        std::cout << a[0][0]  << " \t" << a[0][1]  << " \t" << a[0][2]  << " \n";
        std::cout << a[1][0]  << " \t" << a[1][1]  << " \t" << a[1][2]  << " \n";
        std::cout << a[2][0]  << " \t" << a[2][1]  << " \t" << a[2][2]  << " \n";

        return 1.0;
    }


    GAHOShowValue() {}
    GAHOShowValue (const GAHOShowValue&) {}
    ~GAHOShowValue() {}
};



class Normalize0
{
public:
    typedef LifeV::VectorSmall<3> return_Type;

    return_Type operator() (const LifeV::VectorSmall<3>& a)
    {
        LifeV::VectorSmall<3>   n;

        LifeV::Real norm = std::sqrt (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
        if ( norm != 0 )
        {
            n[0] = a[0] / norm;
            n[1] = a[1] / norm;
            n[2] = a[2] / norm;
        }
        else
        {
            n[0] = 1.0;
            n[1] = 0.0;
            n[2] = 0.0;
        }
        return n;
    }

    Normalize0() {}
    ~Normalize0() {}
};

class Normalize1
{
public:
    typedef LifeV::VectorSmall<3> return_Type;

    return_Type operator() (const LifeV::VectorSmall<3>& a)
    {
        LifeV::VectorSmall<3>   n;

        LifeV::Real norm = std::sqrt (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
        if ( norm != 0 )
        {
            n[0] = a[0] / norm;
            n[1] = a[1] / norm;
            n[2] = a[2] / norm;
        }
        else
        {
            n[0] = 0.0;
            n[1] = 1.0;
            n[2] = 0.0;
        }
        return n;
    }

    Normalize1() {}
    ~Normalize1() {}
};

class Normalize2
{
public:
    typedef LifeV::VectorSmall<3> return_Type;

    return_Type operator() (const LifeV::VectorSmall<3>& a)
    {
        LifeV::VectorSmall<3>   n;

        LifeV::Real norm = std::sqrt (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
        if ( norm != 0 )
        {
            n[0] = a[0] / norm;
            n[1] = a[1] / norm;
            n[2] = a[2] / norm;
        }
        else
        {
            n[0] = 0.0;
            n[1] = 0.0;
            n[2] = 1.0;
        }
        return n;
    }

    Normalize2() {}
    ~Normalize2() {}
};

// Isotropic Part
class dW1
{
public:
    typedef LifeV::Real return_Type;

    return_Type operator() (const LifeV::Real& I1)
    {
        if (I1 >= 3.0 - 1.e-15 && I1 <= 3.0 + 1.e-15)
        {
            return 0.5 * M_a;
        }
        else
        {
            return M_a / 2. * std::exp (M_b * (I1 - 3.0) );
        }
    }

    dW1 (LifeV::Real a, LifeV::Real b)
        : M_a (a), M_b (b)
    {}
    dW1 (const dW1& copy)
        : M_a (copy.M_a), M_b (copy.M_b)
    {}
    ~dW1() {}

private:
    LifeV::Real M_a, M_b;
};

class d2W1
{
public:
    typedef LifeV::Real return_Type;

    return_Type operator() (const LifeV::Real& I1)
    {
        if (I1 >= 3.0 - 1.e-15 && I1 <= 3.0 + 1.e-15)
        {
            return 0.5 * M_a * M_b;
        }
        else
        {
            return M_a * M_b / 2. * std::exp (M_b * (I1 - 3.0) );
        }
    }

    d2W1 (LifeV::Real a, LifeV::Real b)
        : M_a (a), M_b (b)
    {}
    d2W1 (const d2W1& copy)
        : M_a (copy.M_a), M_b (copy.M_b)
    {}
    ~d2W1() {}

private:
    LifeV::Real M_a, M_b;
};

// {Fiber,Sheet}-specific Part
class dW4
{
public:
    typedef LifeV::Real return_Type;

    return_Type operator() (const LifeV::Real& I4)
    {
        if (I4 == 1.0)
        {
            return 0.0;
        }
        else return M_a * (I4 > 1.)
                        * (I4 - 1)
                        * std::exp (M_b * std::pow (I4 - 1, 2.0) );
    }

    dW4 (LifeV::Real a, LifeV::Real b)
        : M_a (a), M_b (b)
    {}
    dW4 (const dW4& copy)
        : M_a (copy.M_a), M_b (copy.M_b)
    {}
    ~dW4() {}

private:
    LifeV::Real M_a, M_b;
};

class d2W4
{
public:
    typedef LifeV::Real return_Type;

    return_Type operator() (const LifeV::Real& I4)
    {
        if (I4 == 1.0)
        {
            return M_a;
        }
        else return M_a * (I4 >= 1.)
                        * (1.0 + 2.0 * M_b * std::pow (I4 - 1, 2.0) )
                        * std::exp (M_b * std::pow (I4 - 1, 2.0) );
    }

    d2W4 (LifeV::Real a, LifeV::Real b)
        : M_a (a), M_b (b)
    {}
    d2W4 (const d2W4& copy)
        : M_a (copy.M_a), M_b (copy.M_b)
    {}
    ~d2W4() {}

private:
    LifeV::Real M_a, M_b;
};

// Cross FiberSheet-specific Part
class dW8
{
public:
    typedef LifeV::Real return_Type;

    return_Type operator() (const LifeV::Real& I8)
    {
        if (I8 == 0)
        {
            return 0.0;
        }
        else
        {
            return M_a * I8 * std::exp (M_b * std::pow (I8, 2.0) );
        }
    }

    dW8 (LifeV::Real a, LifeV::Real b)
        : M_a (a), M_b (b)
    {}
    dW8 (const dW8& copy)
        : M_a (copy.M_a), M_b (copy.M_b)
    {}
    ~dW8() {}

private:
    LifeV::Real M_a, M_b;
};

class d2W8
{
public:
    typedef LifeV::Real return_Type;

    return_Type operator() (const LifeV::Real& I8)
    {
        if (I8 == 0)
        {
            return M_a;
        }
        else return  M_a * (1.0 + 2.0 * M_b * std::pow (I8, 2.0) )
                         * std::exp (M_b * pow (I8, 2.0) );
    }

    d2W8 (LifeV::Real a, LifeV::Real b)
        : M_a (a), M_b (b)
    {}
    d2W8 (const d2W8& copy)
        : M_a (copy.M_a), M_b (copy.M_b)
    {}
    ~d2W8() {}

private:
    LifeV::Real M_a, M_b;
};

// Volumetric part
class dWvol
{
public:
    typedef LifeV::Real return_Type;

    return_Type operator() (const LifeV::Real& J)
    {
        return M_kappa / 2. * (J - 1. / J);
    }

    dWvol (LifeV::Real kappa)
        : M_kappa (kappa)
    {}
    dWvol (const dWvol& copy)
        : M_kappa (copy.M_kappa)
    {}
    ~dWvol() {}

private:
    LifeV::Real M_kappa;
};

class d2Wvol
{
public:
    typedef LifeV::Real return_Type;

    return_Type operator() (const LifeV::Real& J)
    {
        LifeV::Real Jm1 = 1. / J;
        return M_kappa / 2. * (1 + Jm1 * Jm1);
    }

    d2Wvol (LifeV::Real kappa)
        : M_kappa (kappa)
    {}
    d2Wvol (const d2Wvol& copy)
        : M_kappa (copy.M_kappa)
    {}
    ~d2Wvol() {}

private:
    LifeV::Real M_kappa;
};

}

}


namespace LifeV
{

template <typename MeshType>
class EMGeneralizedActiveHolzapfelOgdenMaterial :
    public EMActiveStructuralConstitutiveLaw<MeshType>
{
    //!@name Type definitions
    //@{

public:
    typedef EMActiveStructuralConstitutiveLaw<MeshType>          super;

    typedef StructuralConstitutiveLawData            data_Type;

    typedef typename super::vector_Type              vector_Type;
    typedef typename super::matrix_Type              matrix_Type;

    typedef typename super::matrixPtr_Type           matrixPtr_Type;
    typedef typename super::vectorPtr_Type           vectorPtr_Type;
    typedef typename super::dataPtr_Type             dataPtr_Type;
    typedef typename super::displayerPtr_Type        displayerPtr_Type;

    typedef typename super::mapMarkerVolumesPtr_Type mapMarkerVolumesPtr_Type;
    typedef typename super::mapMarkerVolumes_Type mapMarkerVolumes_Type;
    typedef typename mapMarkerVolumes_Type::const_iterator mapIterator_Type;

    typedef typename super::vectorVolumes_Type       vectorVolumes_Type;
    typedef boost::shared_ptr<vectorVolumes_Type>    vectorVolumesPtr_Type;

    typedef std::vector<UInt>                        vectorIndexes_Type;
    typedef boost::shared_ptr<vectorIndexes_Type>    vectorIndexesPtr_Type;

    typedef typename super::mapMarkerIndexesPtr_Type mapMarkerIndexesPtr_Type;
    typedef typename super::mapMarkerIndexes_Type    mapMarkerIndexes_Type;
    typedef typename mapMarkerIndexes_Type::const_iterator mapIteratorIndex_Type;

    typedef typename super::FESpacePtr_Type          FESpacePtr_Type;
    typedef typename super::ETFESpacePtr_Type        ETFESpacePtr_Type;
    typedef typename super::ETFESpace_Type              ETFESpace_Type;

    typedef MeshType                                        mesh_Type;
    typedef ETFESpace< mesh_Type, MapEpetra, 3, 1 >                        scalarETFESpace_Type;
    typedef boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > >    scalarETFESpacePtr_Type;

    //Vector for vector parameters
    typedef typename super::vectorsParameters_Type       vectorsParameters_Type;
    typedef typename super::vectorsParametersPtr_Type    vectorsParametersPtr_Type;

    //@}

    // Parameters (public, can be freely changed)
    Real M_aiso, M_biso, M_af, M_bf, M_as, M_bs, M_afs, M_bfs;
    Real M_kappa, M_contractileFraction;

    //! @name Constructor &  Destructor
    //@{

    EMGeneralizedActiveHolzapfelOgdenMaterial();

    virtual  ~EMGeneralizedActiveHolzapfelOgdenMaterial();

    //@}

    //!@name Methods
    //@{

    //! Setup the created object of the class StructuralConstitutiveLaw
    /*!
      \param dFespace: the FiniteElement Space
      \param monolithicMap: the MapEpetra
      \param offset: the offset parameter used assembling the matrices
    */
    inline void setup ( const FESpacePtr_Type& dFESpace,
                        const ETFESpacePtr_Type& dETFESpace,
                        const boost::shared_ptr<const MapEpetra>&   monolithicMap,
                        const UInt offset, const dataPtr_Type& dataMaterial,
                        const displayerPtr_Type& displayer  );


    void setup ( const FESpacePtr_Type& dFESpace,
                 const ETFESpacePtr_Type& dETFESpace,
                 const ETFESpacePtr_Type& activationSpace,
                 const boost::shared_ptr<const MapEpetra>&  monolithicMap,
                 const UInt offset, const dataPtr_Type& dataMaterial, const displayerPtr_Type& displayer );

    void setup ( const FESpacePtr_Type& dFESpace,
                 const ETFESpacePtr_Type& dETFESpace,
                 const ETFESpacePtr_Type& activationSpace,
                 const boost::shared_ptr<const MapEpetra>&  monolithicMap,
                 const UInt offset, const dataPtr_Type& dataMaterial, const displayerPtr_Type& displayer,
                 const vector_Type& gammaf );

    void setup ( const FESpacePtr_Type& dFESpace,
                 const ETFESpacePtr_Type& dETFESpace,
                 const ETFESpacePtr_Type& activationSpace,
                 const boost::shared_ptr<const MapEpetra>&  monolithicMap,
                 const UInt offset, const dataPtr_Type& dataMaterial, const displayerPtr_Type& displayer,
                 const vector_Type& gammaf,
                 const vector_Type& fiberVector);

    void setDefaultParams();

    //! Compute the Stiffness matrix in StructuralSolver::buildSystem()
    /*!
      \param dataMaterial the class with Material properties data
    */
    void computeLinearStiff ( dataPtr_Type& /*dataMaterial*/, const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                              const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/ );


    //! Updates the Jacobian matrix in StructualSolver::updateJacobian
    /*!
      \param disp: solution at the k-th iteration of NonLinearRichardson Method
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get
                           the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void updateJacobianMatrix ( const vector_Type& disp,
                                const dataPtr_Type& dataMaterial,
                                const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                const displayerPtr_Type& displayer);


    //! Updates the nonlinear terms in the Jacobian matrix in StructualSolver::updateJacobian
    /*!
      \param stiff: stiffness matrix provided from outside
      \param disp: solution at the k-th iteration of NonLinearRichardson Method
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get
                           the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void updateNonLinearJacobianTerms ( matrixPtr_Type& jacobian,
                                        const vector_Type& disp,
                                        const dataPtr_Type& dataMaterial,
                                        const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                                        const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/,
                                        const displayerPtr_Type& displayer);

    //! Interface method to compute the new Stiffness matrix in StructuralSolver::evalResidual and in
    //! StructuralSolver::updateSystem since the matrix is the expression of the matrix is the same.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get
                           the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void computeStiffness ( const vector_Type& disp, Real factor, const dataPtr_Type& dataMaterial,
                            const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                            const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/,
                            const displayerPtr_Type& displayer );

    void computeResidual ( const vector_Type& disp );

    //! Computes the new Stiffness vector for Neo-Hookean and Exponential materials in
    //! StructuralSolver given a certain displacement field.
    //! This function is used both in StructuralSolver::evalResidual and in StructuralSolver::updateSystem
    //! since the matrix is the expression of the matrix is the same.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get
                           the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    // void computeVector( const vector_Type& sol,
    //                     Real factor,
    //                     const dataPtr_Type& dataMaterial,
    //                     const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
    //                     const displayerPtr_Type& displayer );


    //! Computes the deformation gradient F, the cofactor matrix Cof(F),
    //! the determinant of F (J = det(F)), the trace of right Cauchy-Green tensor tr(C)
    //! This function is used in StructuralConstitutiveLaw::computeStiffness
    /*!
      \param dk_loc: the elemental displacement
    */
    void computeKinematicsVariables ( const VectorElemental& dk_loc ) {}

    //! ShowMe method of the class (saved on a file the stiffness vector and the jacobian)
    void showMe ( std::string const& fileNameVectStiff,
                  std::string const& fileNameJacobain);

    //! ShowMe method of the class (saved on a file the stiffness vector and the jacobian)
    void showMyParameters ();


    //! Compute the First Piola Kirchhoff Tensor
    /*!
       \param firstPiola Epetra_SerialDenseMatrix that has to be filled
       \param tensorF Epetra_SerialDenseMatrix the deformation gradient
       \param cofactorF Epetra_SerialDenseMatrix cofactor of F
       \param invariants std::vector with the invariants of C and the detF
       \param material UInt number to get the material parameteres form the VenantElasticData class
    */
    void computeLocalFirstPiolaKirchhoffTensor ( Epetra_SerialDenseMatrix& firstPiola,
                                                 const Epetra_SerialDenseMatrix& tensorF,
                                                 const Epetra_SerialDenseMatrix& cofactorF,
                                                 const std::vector<Real>& invariants,
                                                 const UInt marker);


    void computeRes ( vectorPtr_Type& res, const vector_Type& disp, const displayerPtr_Type& displayer );


    //@}

    //! @name Get Methods
    //@{

    //! Get the Stiffness matrix
    inline  matrixPtr_Type const stiffMatrix() const
    {
        return super::M_jacobian;
    }

    //! Get the stiffness vector
    inline  vectorPtr_Type const stiffVector() const
    {
        return M_stiff;
    }

    inline  vectorPtr_Type fiberVectorPtr()
    {
        return M_fiberVector;
    }

    inline  vectorPtr_Type sheetVectorPtr()
    {
        return M_sheetVector;
    }

    inline  vectorPtr_Type gammaf()
    {
        return M_Gammaf;
    }
    inline  vectorPtr_Type gammas()
    {
        return M_Gammas;
    }
    inline  vectorPtr_Type gamman()
    {
        return M_Gamman;
    }
    inline  scalarETFESpacePtr_Type activationSpace()
    {
        return M_activationSpace;
    }


    void apply ( const vector_Type& sol, vector_Type& res,
                 const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                 const mapMarkerIndexesPtr_Type mapsMarkerIndexes);

    //@}


    //! @name set Methods
    //@{
    inline void setGammaf ( const vector_Type& gammaf)
    {
        M_Gammaf.reset ( new vector_Type ( gammaf ) );
    }
    inline void setGammas ( const vector_Type& gammas)
    {
        M_Gammas.reset ( new vector_Type ( gammas ) );
    }
    inline void setGamman ( const vector_Type& gamman)
    {
        M_Gamman.reset ( new vector_Type ( gamman ) );
    }

    inline void setFiberVector ( const vector_Type& fiberVector)
    {
        M_fiberVector.reset ( new vector_Type ( fiberVector ) );
        ElectrophysiologyUtility::normalize (*M_fiberVector);
    }

    inline void setSheetVector ( const vector_Type& sheetVector)
    {
        M_sheetVector.reset ( new vector_Type ( sheetVector ) );
        ElectrophysiologyUtility::normalize (*M_sheetVector);
    }

    //@}

    //! @name Methods
    //@{

    inline void setupFiberVector ( std::string& name, boost::shared_ptr<mesh_Type> mesh )
    {
        ElectrophysiologyUtility::importFibers ( M_fiberVector, name, mesh  );
        ElectrophysiologyUtility::normalize (*M_fiberVector);
    }

    inline void setupFiberVector ( std::string& name, std::string& path )
    {
        ElectrophysiologyUtility::importFibers ( M_fiberVector, name, path  );
        ElectrophysiologyUtility::normalize (*M_fiberVector);
    }

    void setupFiberVector ( Real& fx, Real& fy, Real& fz )
    {
        ElectrophysiologyUtility::setupFibers ( *M_fiberVector, fx, fy, fz  );
        ElectrophysiologyUtility::normalize (*M_fiberVector);
    }

    inline void setupSheetVector ( std::string& name, boost::shared_ptr<mesh_Type> mesh )
    {
        ElectrophysiologyUtility::importFibers ( M_sheetVector, name, mesh  );
        ElectrophysiologyUtility::normalize (*M_sheetVector);
    }

    inline void setupSheetVector ( std::string& name, std::string& path )
    {
        ElectrophysiologyUtility::importFibers ( M_sheetVector, name, path  );
        ElectrophysiologyUtility::normalize (*M_sheetVector);
    }

    void setupSheetVector ( Real& sx, Real& sy, Real& sz )
    {
        ElectrophysiologyUtility::setupFibers ( *M_sheetVector, sx, sy, sz);
        ElectrophysiologyUtility::normalize (*M_sheetVector);
    }
    //@}

protected:

    //! construct the vectors for the parameters
    /*!
      \param VOID
      \return VOID
    */
    void setupVectorsParameters ( void );

    //! Vector: stiffness non-linear
    vectorPtr_Type                      M_stiff;
    vectorPtr_Type                      M_fiberVector;
    vectorPtr_Type                      M_sheetVector;

    vectorPtr_Type                      M_Gammaf;
    vectorPtr_Type                      M_Gammas;
    vectorPtr_Type                      M_Gamman;
    scalarETFESpacePtr_Type             M_activationSpace;

};

template <typename MeshType>
EMGeneralizedActiveHolzapfelOgdenMaterial<MeshType>::EMGeneralizedActiveHolzapfelOgdenMaterial() :
    super           ( ),
    M_stiff         ( ),
    M_Gammaf        ( ),
    M_Gammas        ( ),
    M_Gamman        ( ),
    M_fiberVector   ( ),
    M_sheetVector   ( ),
    M_activationSpace( )
{
}





template <typename MeshType>
EMGeneralizedActiveHolzapfelOgdenMaterial<MeshType>::~EMGeneralizedActiveHolzapfelOgdenMaterial()
{}


template <typename MeshType>
void
EMGeneralizedActiveHolzapfelOgdenMaterial<MeshType>::setup ( const FESpacePtr_Type& dFESpace,
                                                           const ETFESpacePtr_Type& dETFESpace,
                                                           const boost::shared_ptr<const MapEpetra>&   monolithicMap,
                                                           const UInt offset, const dataPtr_Type& dataMaterial,
                                                           const displayerPtr_Type& displayer  )
{

    this->M_displayer = displayer;
    this->M_dataMaterial  = dataMaterial;

    //    std::cout<<"I am setting up the Material"<<std::endl;

    //    M_stiff.
    this->M_dispFESpace                 = dFESpace;
    this->M_dispETFESpace               = dETFESpace;
    this->M_localMap                    = monolithicMap;
    this->M_offset                      = offset;
    this->M_dataMaterial                = dataMaterial;
    this->M_displayer                   = displayer;

    M_stiff.reset                   ( new vector_Type (*this->M_localMap) );
    M_fiberVector.reset             ( new vector_Type (*this->M_localMap) );
    M_sheetVector.reset             ( new vector_Type (*this->M_localMap) );
    M_activationSpace.reset         ( new scalarETFESpace_Type ( dETFESpace -> mesh(),
                                                                 &feTetraP1,
                                                                 dFESpace->map().commPtr() ) );
    M_Gammaf.reset                  ( new vector_Type ( M_activationSpace -> map() ) );
    M_Gammas.reset                  ( new vector_Type ( M_activationSpace -> map() ) );
    M_Gamman.reset                  ( new vector_Type ( M_activationSpace -> map() ) );



    //   this->setDefaultParams();

    // The 2 is because the law uses three parameters (mu, bulk).
    // another way would be to set up the number of constitutive parameters of the law
    // in the data file to get the right size. Note the comment below.
    this->M_vectorsParameters.reset ( new vectorsParameters_Type ( 2 ) );

    this->setupVectorsParameters();

}

template <typename MeshType>
void EMGeneralizedActiveHolzapfelOgdenMaterial<MeshType>::setDefaultParams()
{
    // Bulk modulus
    M_kappa = 350.0; // in [kPa]
    // Scalar parameters used in the model (default value, from G\"oktepe, 2011)
    M_aiso = 0.496;
    M_af = 15.193;
    M_as =  3.283;
    M_afs = 0.662; // in [kPa]
    M_biso = 7.209;
    M_bf = 20.417;
    M_bs = 11.176;
    M_bfs = 9.466; // adimensional
    M_contractileFraction = 0.7;
}

template <typename MeshType>
void
EMGeneralizedActiveHolzapfelOgdenMaterial<MeshType>::setup ( const FESpacePtr_Type&                      dFESpace,
                                                           const ETFESpacePtr_Type&                    dETFESpace,
                                                           const ETFESpacePtr_Type& activationSpace,
                                                           const boost::shared_ptr<const MapEpetra>&   monolithicMap,
                                                           const UInt                                  offset,
                                                           const dataPtr_Type&                         dataMaterial,
                                                           const displayerPtr_Type&                    displayer)
{
    this->M_displayer = displayer;
    this->M_dataMaterial  = dataMaterial;

    //    std::cout<<"I am setting up the Material"<<std::endl;

    this->M_dispFESpace                     = dFESpace;
    this->M_dispETFESpace                     = dETFESpace;
    this->M_localMap                    = monolithicMap;
    this->M_offset                      = offset;
    this->M_dataMaterial                = dataMaterial;
    this->M_displayer                   = displayer;
    M_activationSpace               = activationSpace;
    M_stiff.reset                   ( new vector_Type (*this->M_localMap) );
    M_fiberVector.reset             ( new vector_Type (*this->M_localMap) );
    M_sheetVector.reset             ( new vector_Type (*this->M_localMap) );
    M_Gammaf.reset                  ( new vector_Type ( M_activationSpace -> map() ) );
    M_Gammas.reset                  ( new vector_Type ( M_activationSpace -> map() ) );
    M_Gamman.reset                  ( new vector_Type ( M_activationSpace -> map() ) );


    //    this->setDefaultParams();

    // This parameters are not used, but we are keeping them to avoid
    // inconsistencies among different materials
    this->M_vectorsParameters.reset ( new vectorsParameters_Type (2) );
    this->setupVectorsParameters();
}


template <typename MeshType>
void
EMGeneralizedActiveHolzapfelOgdenMaterial<MeshType>::setup ( const FESpacePtr_Type&                      dFESpace,
                                                           const ETFESpacePtr_Type&                    dETFESpace,
                                                           const ETFESpacePtr_Type& activationSpace,
                                                           const boost::shared_ptr<const MapEpetra>&   monolithicMap,
                                                           const UInt                                  offset,
                                                           const dataPtr_Type&                         dataMaterial,
                                                           const displayerPtr_Type&                    displayer,
                                                           const vector_Type&                           gammaf   )
{
    setup ( dFESpace, dETFESpace, activationSpace, monolithicMap, offset, dataMaterial, displayer);
    M_Gammaf = gammaf;
}

template <typename MeshType>
void
EMGeneralizedActiveHolzapfelOgdenMaterial<MeshType>::setup ( const FESpacePtr_Type&                      dFESpace,
                                                           const ETFESpacePtr_Type&                    dETFESpace,
                                                           const ETFESpacePtr_Type& activationSpace,
                                                           const boost::shared_ptr<const MapEpetra>&   monolithicMap,
                                                           const UInt                                  offset,
                                                           const dataPtr_Type&                         dataMaterial,
                                                           const displayerPtr_Type&                    displayer,
                                                           const vector_Type&                           gammaf,
                                                           const vector_Type&                           fiberVector)
{
    setup ( dFESpace, dETFESpace, activationSpace, monolithicMap, offset, dataMaterial, displayer, gammaf );
    M_fiberVector = fiberVector;
}

template <typename MeshType>
void
EMGeneralizedActiveHolzapfelOgdenMaterial<MeshType>::setupVectorsParameters ( void )
{
    // Paolo Tricerri: February, 20th
    // In each class, the name of the parameters has to inserted in the law
    // TODO: move the saving of the material parameters to more abstract objects
    //       such that in the class of the material we do not need to call explicitly
    //       the name of the parameter.

    // Number of volume on the local part of the mesh
    UInt nbElements = this->M_dispFESpace->mesh()->numVolumes();

    // 1. resize the vector in the first element of the vector
    (* (this->M_vectorsParameters) ) [0].resize ( nbElements );
    (* (this->M_vectorsParameters) ) [1].resize ( nbElements );

    for (UInt i (0); i < nbElements; i++ )
    {
        // Extracting the marker
        UInt markerID = this->M_dispFESpace->mesh()->element ( i ).markerID();

        Real mu = this->M_dataMaterial->mu ( markerID );
        Real bulk = this->M_dataMaterial->bulk ( markerID );

        ( (* (this->M_vectorsParameters) ) [0]) [ i ] = mu;
        ( (* (this->M_vectorsParameters) ) [1]) [ i ] = bulk;

        M_aiso  = this->M_dataMaterial->A ( markerID );
        M_af    = this->M_dataMaterial->Af ( markerID );
        M_as    = this->M_dataMaterial->As ( markerID );
        M_afs   = this->M_dataMaterial->Afs ( markerID );
        M_biso  = this->M_dataMaterial->B ( markerID );
        M_bf    = this->M_dataMaterial->Bf ( markerID );
        M_bs    = this->M_dataMaterial->Bs ( markerID );
        M_bfs   = this->M_dataMaterial->Bfs ( markerID );
        M_kappa = this->M_dataMaterial->bulk ( markerID );
        M_contractileFraction = this->M_dataMaterial->contractileFraction ( markerID );
    }
}


template <typename MeshType>
void EMGeneralizedActiveHolzapfelOgdenMaterial<MeshType>::computeLinearStiff (dataPtr_Type& /*dataMaterial*/,
                                                                            const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                                                                            const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/)
{
    //! Empty method for neo-hookean material
}


template <typename MeshType>
void EMGeneralizedActiveHolzapfelOgdenMaterial<MeshType>::updateJacobianMatrix ( const vector_Type&       disp,
                                                                               const dataPtr_Type&      dataMaterial,
                                                                               const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                                               const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                                                               const displayerPtr_Type& displayer )
{
    this->M_jacobian.reset (new matrix_Type (*this->M_localMap) );

    displayer->leaderPrint (" \n*********************************\n  ");
    updateNonLinearJacobianTerms (this->M_jacobian, disp, dataMaterial, mapsMarkerVolumes, mapsMarkerIndexes, displayer);
    displayer->leaderPrint (" \n*********************************\n  ");
    std::cout << std::endl;
}



template <typename MeshType>
void EMGeneralizedActiveHolzapfelOgdenMaterial<MeshType>::updateNonLinearJacobianTerms ( matrixPtr_Type&       jacobian,
        const vector_Type&    disp,
        const dataPtr_Type&   dataMaterial,
        const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
        const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
        const displayerPtr_Type&  displayer )
{
    {
        using namespace ExpressionAssembly;

        displayer->leaderPrint ("   Non-Linear S-  updating non linear terms in the Jacobian Matrix (Generalized Active Holzapfel-Ogden)");

        * (jacobian) *= 0.0;

        // Activation is Fa = gamma (fo x fo) + 1/sqrt(gamma) * (I - fo x fo)
        BOOST_AUTO_TPL (gammaf,  value (this->M_activationSpace, *M_Gammaf) );
        BOOST_AUTO_TPL (gammas,  value (this->M_activationSpace, *M_Gammas) );
        BOOST_AUTO_TPL (gamman,  value (this->M_activationSpace, *M_Gamman) );

        MatrixSmall<3, 3> Id;
        Id (0, 0) = 1.;
        Id (0, 1) = 0., Id (0, 2) = 0.;
        Id (1, 0) = 0.;
        Id (1, 1) = 1., Id (1, 2) = 0.;
        Id (2, 0) = 0.;
        Id (2, 1) = 0., Id (2, 2) = 1.;

        // Strain-energy terms
        boost::shared_ptr<StrainEnergyGHO::dW1> dW1fun (new StrainEnergyGHO::dW1 (M_aiso, M_biso) );
        boost::shared_ptr<StrainEnergyGHO::d2W1> d2W1fun (new StrainEnergyGHO::d2W1 (M_aiso, M_biso) );

        boost::shared_ptr<StrainEnergyGHO::dW4> dW4ffun (new StrainEnergyGHO::dW4 (M_af, M_bf) );
        boost::shared_ptr<StrainEnergyGHO::d2W4> d2W4ffun (new StrainEnergyGHO::d2W4 (M_af, M_bf) );

        boost::shared_ptr<StrainEnergyGHO::dW4> dW4sfun (new StrainEnergyGHO::dW4 (M_as, M_bs) );
        boost::shared_ptr<StrainEnergyGHO::d2W4> d2W4sfun (new StrainEnergyGHO::d2W4 (M_as, M_bs) );

        boost::shared_ptr<StrainEnergyGHO::dW8> dW8fsfun (new StrainEnergyGHO::dW8 (M_afs, M_bfs) );
        boost::shared_ptr<StrainEnergyGHO::d2W8> d2W8fsfun (new StrainEnergyGHO::d2W8 (M_afs, M_bfs) );

        boost::shared_ptr<StrainEnergyGHO::dWvol> dWvolfun (new StrainEnergyGHO::dWvol (M_kappa) );
        boost::shared_ptr<StrainEnergyGHO::d2Wvol> d2Wvolfun (new StrainEnergyGHO::d2Wvol (M_kappa) );

        boost::shared_ptr<StrainEnergyGHO::GAHOShowValue> sv (new StrainEnergyGHO::GAHOShowValue() );

        boost::shared_ptr<StrainEnergyGHO::Normalize0> normalize0 (new StrainEnergyGHO::Normalize0() );
        boost::shared_ptr<StrainEnergyGHO::Normalize1> normalize1 (new StrainEnergyGHO::Normalize1() );
        boost::shared_ptr<StrainEnergyGHO::Normalize2> normalize2 (new StrainEnergyGHO::Normalize2() );

        // Kinematics
        BOOST_AUTO_TPL (I,      value (Id) );
        BOOST_AUTO_TPL (Grad_u, grad (this->M_dispETFESpace, disp, this->M_offset) );
        BOOST_AUTO_TPL (F,      Grad_u + I);
        BOOST_AUTO_TPL (FmT,    minusT (F) );
        BOOST_AUTO_TPL (J,      det (F) );
        // Fibres
        BOOST_AUTO_TPL (f01,    value (this->M_dispETFESpace, *M_fiberVector) );
        //BOOST_AUTO_TPL(f0,     f01 / sqrt( dot(f01, f01) ) );
        BOOST_AUTO_TPL (f0,     eval (normalize0, f01) );

        BOOST_AUTO_TPL (s01,     value (this->M_dispETFESpace, *M_sheetVector) );
        BOOST_AUTO_TPL (s02,     eval (normalize1, s01) );
        BOOST_AUTO_TPL (s0,      eval (normalize2, s02 - dot (f0, s02) * f0) );

        // Invariants
        BOOST_AUTO_TPL (I1,     dot (F, F) );
        BOOST_AUTO_TPL (C,      transpose (F) * F);
        BOOST_AUTO_TPL (I2,    value (0.5) * (I1 * I1 - dot (C, C) ) );


        BOOST_AUTO_TPL (I4f,    dot (F * f0, F * f0) );
        BOOST_AUTO_TPL (I4s,    dot (F * s0, F * s0) );
        BOOST_AUTO_TPL (I5f,    dot ( C * f0, C * f0) );
        BOOST_AUTO_TPL (I5s,    dot ( C * s0, C * s0) );

        BOOST_AUTO_TPL (I8fs,   dot (F * f0, F * s0) );
        //    BOOST_AUTO_TPL(I8fs2,   I2 + I4f * I4s + I5f + I5s - I1 * ( I4f + I4s ) );
        //    BOOST_AUTO_TPL(I8fs,    sqrt ( I8fs2 ) );
        // Reduced invariants
        BOOST_AUTO_TPL (Jm23,    pow (J, -2. / 3) );
        BOOST_AUTO_TPL (I1iso,   Jm23 * I1);
        BOOST_AUTO_TPL (I4fiso,  Jm23 * I4f);
        BOOST_AUTO_TPL (I4siso,  Jm23 * I4s);
        BOOST_AUTO_TPL (I8fsiso, Jm23 * I8fs);
        // Generalised invariants
        //BOOST_AUTO_TPL(gammac, pow(gamma, -2) - gamma);
        BOOST_AUTO_TPL (dI1edI1,   value (1.0) - gamman * ( gamman + value (2.0) ) * pow (gamman + value (1.0), -2.0) );
        BOOST_AUTO_TPL (dI1edI4f,   gamman * ( gamman + value (2.0) ) * pow (gamman + value (1.0), -2.0) - gammaf * ( gammaf + value (2.0) ) * pow (gammaf + value (1.0), -2.0) );
        BOOST_AUTO_TPL (dI1edI4s,   gamman * ( gamman + value (2.0) ) * pow (gamman + value (1.0), -2.0) - gammas * ( gammas + value (2.0) ) * pow (gammas + value (1.0), -2.0) );
        BOOST_AUTO_TPL (dI4fedI4f,   value (1.0) / ( ( gammaf + value (1.0) ) * ( gammaf + value (1.0) ) ) );
        BOOST_AUTO_TPL (dI4sedI4s,   value (1.0) / ( ( gammas + value (1.0) ) * ( gammas + value (1.0) ) ) );
        //    BOOST_AUTO_TPL(gammafp1,  gammaf + value(1.0));
        //    BOOST_AUTO_TPL(gammasp1,  gammaf + value(1.0));
        //    BOOST_AUTO_TPL(dI8fsedI8fs, value(1.0) / ( gammafp1 * gammasp1 ) );
        BOOST_AUTO_TPL (dI8fsedI8fs, value (1.0) / ( ( gammaf + value (1.0) ) * ( gammas + value (1.0) ) ) );


        BOOST_AUTO_TPL (I1eiso,   dI1edI1 * I1iso + dI1edI4f * I4fiso + dI1edI4s * I4siso);
        BOOST_AUTO_TPL (I4feiso,  dI4fedI4f * I4fiso);
        BOOST_AUTO_TPL (I4seiso,  dI4sedI4s * I4siso);
        BOOST_AUTO_TPL (I8fseiso, dI8fsedI8fs * I8fsiso);

        // Strain-energy derivatives
        BOOST_AUTO_TPL (dW1,    value (1.0 - M_contractileFraction) * eval (dW1fun, I1iso) );
        BOOST_AUTO_TPL (d2W1,   value (1.0 - M_contractileFraction) * eval (d2W1fun, I1iso) );
        BOOST_AUTO_TPL (dW4f,   value (1.0 - M_contractileFraction) * eval (dW4ffun, I4fiso) );
        BOOST_AUTO_TPL (d2W4f,  value (1.0 - M_contractileFraction) * eval (d2W4ffun, I4fiso) );
        BOOST_AUTO_TPL (dW4s,   value (1.0 - M_contractileFraction) * eval (dW4sfun, I4siso) );
        BOOST_AUTO_TPL (d2W4s,  value (1.0 - M_contractileFraction) * eval (d2W4sfun, I4siso) );
        BOOST_AUTO_TPL (dW8fs,  value (1.0 - M_contractileFraction) * eval (dW8fsfun, I8fsiso) );
        BOOST_AUTO_TPL (d2W8fs, value (1.0 - M_contractileFraction) * eval (d2W8fsfun, I8fsiso) );

        BOOST_AUTO_TPL (dW1e,    value ( M_contractileFraction ) * eval (dW1fun, I1eiso) );
        BOOST_AUTO_TPL (d2W1e,   value ( M_contractileFraction ) * eval (d2W1fun, I1eiso) );
        BOOST_AUTO_TPL (dW4fe,   /* value( M_contractileFraction ) */ eval (dW4ffun, I4feiso) );
        BOOST_AUTO_TPL (d2W4fe,  /* value( M_contractileFraction ) */ eval (d2W4ffun, I4feiso) );
        BOOST_AUTO_TPL (dW4se,   /* value( M_contractileFraction ) */ eval (dW4sfun, I4seiso) );
        BOOST_AUTO_TPL (d2W4se,  /* value( M_contractileFraction ) */ eval (d2W4sfun, I4seiso) );
        BOOST_AUTO_TPL (dW8fse,  /* value( M_contractileFraction ) */ eval (dW8fsfun, I8fseiso) );
        BOOST_AUTO_TPL (d2W8fse, /* value( M_contractileFraction ) */ eval (d2W8fsfun, I8fseiso) );


        BOOST_AUTO_TPL (dWvol,  eval (dWvolfun, J) );
        BOOST_AUTO_TPL (d2Wvol, eval (d2Wvolfun, J) );

        BOOST_AUTO_TPL (dF,   grad (phi_j) );
        BOOST_AUTO_TPL (dJ,   J * FmT);
        BOOST_AUTO_TPL (dJm23,   value (-2.0 / 3.0) * Jm23 * FmT);
        BOOST_AUTO_TPL (dFmTdF, value (-1.0) * FmT * transpose (dF) * FmT);

        BOOST_AUTO_TPL (dI1,   value (2.0) * F );
        BOOST_AUTO_TPL (dI4f,  value (2.0) * F * outerProduct (f0, f0) );
        BOOST_AUTO_TPL (dI4s,  value (2.0) * F * outerProduct (s0, s0) );
        BOOST_AUTO_TPL (dI8fs, F * ( outerProduct (f0, s0) + outerProduct (s0, f0) ) );

        BOOST_AUTO_TPL (dI1iso,   dJm23 * I1 + Jm23 * dI1 );
        BOOST_AUTO_TPL (dI4fiso,  dJm23 * I4f + Jm23 * dI4f );
        BOOST_AUTO_TPL (dI4siso,  dJm23 * I4s + Jm23 * dI4s );
        BOOST_AUTO_TPL (dI8fsiso, dJm23 * I8fs + Jm23 * dI8fs );

        BOOST_AUTO_TPL (dI1eiso,  dI1edI1 * dI1iso + dI1edI4f * dI4fiso + dI1edI4s * dI4siso );
        BOOST_AUTO_TPL (dI4feiso, dI4fedI4f * dI4fiso );
        BOOST_AUTO_TPL (dI4seiso, dI4sedI4s * dI4siso );
        BOOST_AUTO_TPL (dI8fseiso, dI8fsedI8fs * dI8fsiso );

        BOOST_AUTO_TPL (Piso_1,   dW1 * dI1eiso);
        BOOST_AUTO_TPL (Piso_4f,  dW4f * dI4feiso);
        BOOST_AUTO_TPL (Piso_4s,  dW4s * dI4seiso);
        BOOST_AUTO_TPL (Piso_8fs, dW8fs * dI8fseiso);

        BOOST_AUTO_TPL (dI1dF,   value (2.0) * dot ( F, dF ) );
        BOOST_AUTO_TPL (d2I1dF,  value (2.0) * dF );
        BOOST_AUTO_TPL (dI4fdF,  value (2.0) * dot ( F * outerProduct (f0, f0), dF ) );
        BOOST_AUTO_TPL (d2I4fdF, value (2.0) * dF * outerProduct (f0, f0) );
        BOOST_AUTO_TPL (dI4sdF,  value (2.0) * dot ( F * outerProduct (s0, s0), dF ) );
        BOOST_AUTO_TPL (d2I4sdF, value (2.0) * dF * outerProduct (s0, s0) );
        BOOST_AUTO_TPL (dI8fsdF, dot ( F * ( outerProduct (f0, s0) + outerProduct (s0, f0) ), dF ) );
        BOOST_AUTO_TPL (d2I8fsdF, dF * ( outerProduct (f0, s0) + outerProduct (s0, f0) ) );

        BOOST_AUTO_TPL (dJm23dF,   dot ( dJm23, dF ) );
        BOOST_AUTO_TPL (d2Jm23dF,   value (-2.0 / 3.0) * dJm23dF * FmT + value (-2.0 / 3.0) * Jm23 * dFmTdF );

        BOOST_AUTO_TPL (d2I1isodF,   d2Jm23dF * I1 + dJm23 * dI1dF + dJm23dF * dI1 + Jm23 * d2I1dF );
        BOOST_AUTO_TPL (d2I4fisodF,  d2Jm23dF * I4f + dJm23 * dI4fdF + dJm23dF * dI4f + Jm23 * d2I4fdF );
        BOOST_AUTO_TPL (d2I4sisodF,  d2Jm23dF * I4s + dJm23 * dI4sdF + dJm23dF * dI4s + Jm23 * d2I4sdF );
        BOOST_AUTO_TPL (d2I8fsisodF,  d2Jm23dF * I8fs + dJm23 * dI8fsdF + dJm23dF * dI8fs + Jm23 * d2I8fsdF );

        BOOST_AUTO_TPL (d2I1eisodF,  dI1edI1 * d2I1isodF + dI1edI4f * d2I4fisodF + dI1edI4s * d2I4sisodF);
        BOOST_AUTO_TPL (d2I4feisodF, dI4fedI4f * d2I4fisodF );
        BOOST_AUTO_TPL (d2I4seisodF, dI4sedI4s * d2I4sisodF );
        BOOST_AUTO_TPL (d2I8fseisodF, dI8fsedI8fs * d2I8fsisodF );

        BOOST_AUTO_TPL (dPeiso_1,   dW1e * d2I1eisodF     + d2W1e * dot ( dI1eiso,  dF) * dI1eiso  );
        BOOST_AUTO_TPL (dPeiso_4f,  dW4fe * d2I4feisodF   + d2W4fe * dot ( dI4feiso, dF) * dI4feiso );
        BOOST_AUTO_TPL (dPeiso_4s,  dW4se * d2I4seisodF   + d2W4se * dot ( dI4seiso, dF) * dI4seiso  );
        BOOST_AUTO_TPL (dPeiso_8fs, dW8fse * d2I8fseisodF + d2W8fse * dot ( dI8fseiso, dF) * dI8fseiso  );

        BOOST_AUTO_TPL (dPiso_1,   dW1 * d2I1isodF     + d2W1 * dot ( dI1iso,  dF) * dI1iso  );
        BOOST_AUTO_TPL (dPiso_4f,  dW4f * d2I4fisodF   + d2W4f * dot ( dI4fiso, dF) * dI4fiso );
        BOOST_AUTO_TPL (dPiso_4s,  dW4s * d2I4sisodF   + d2W4s * dot ( dI4siso, dF) * dI4siso  );
        BOOST_AUTO_TPL (dPiso_8fs, dW8fs * d2I8fsisodF + d2W8fs * dot ( dI8fsiso, dF) * dI8fsiso  );

        // Volumetric part
        BOOST_AUTO_TPL (Pvol,  dWvol * dJ);
        BOOST_AUTO_TPL (dPvol, (d2Wvol + J * dWvol) *dot (FmT, dF) *dJ + dWvol * J * dFmTdF);

        // Assembly
        integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                    this->M_dispFESpace->qr(),
                    this->M_dispETFESpace,
                    this->M_dispETFESpace,
                    dot ( dPiso_1 + dPeiso_1 + dPvol, grad (phi_i) )
                  ) >> jacobian;
        if (M_af > 0)
        {
            integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                        this->M_dispFESpace->qr(),
                        this->M_dispETFESpace,
                        this->M_dispETFESpace,
                        dot (/* dPiso_4f +*/ dPeiso_4f , grad (phi_i) )
                      ) >> jacobian;
        }
        if (M_sheetVector)
        {
            if (M_as > 0)
            {
                integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                            this->M_dispFESpace->qr(),
                            this->M_dispETFESpace,
                            this->M_dispETFESpace,
                            dot ( /*dPiso_4s + */dPeiso_4s, grad (phi_i) )
                          ) >> jacobian;
            }
            if (M_afs > 0)
            {
                integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                            this->M_dispFESpace->qr(),
                            this->M_dispETFESpace,
                            this->M_dispETFESpace,
                            dot ( dPeiso_8fs, grad (phi_i) )
                          ) >> jacobian;
            }
        }

    }//end expression assembly


    jacobian->globalAssemble();
}


template <typename MeshType>
void EMGeneralizedActiveHolzapfelOgdenMaterial<MeshType>::apply ( const vector_Type& sol, vector_Type& res,
                                                                const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                                const mapMarkerIndexesPtr_Type mapsMarkerIndexes)
{
    computeStiffness (sol, 0., this->M_dataMaterial, mapsMarkerVolumes, mapsMarkerIndexes, this->M_displayer);
    res += *M_stiff;
}

template <typename MeshType>
void EMGeneralizedActiveHolzapfelOgdenMaterial<MeshType>::computeStiffness ( const vector_Type&       disp,
                                                                           Real                     /*factor*/,
                                                                           const dataPtr_Type&      dataMaterial,
                                                                           const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                                           const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                                                           const displayerPtr_Type& displayer )
{
    displayer->leaderPrint (" \n******************************************************************\n  ");
    displayer->leaderPrint (" Non-Linear S-  Computing the GAHO residual vector"     );
    displayer->leaderPrint (" \n******************************************************************\n  ");

    computeResidual (disp);
}


template <typename MeshType>
void EMGeneralizedActiveHolzapfelOgdenMaterial<MeshType>::computeResidual ( const vector_Type& disp )
{
    using namespace ExpressionAssembly;

    this->M_stiff.reset (new vector_Type (*this->M_localMap) );

    M_stiff.reset (new vector_Type (*this->M_localMap) );
    * (M_stiff) *= 0.0;

    // Activation is Fa = gamma (fo x fo) + 1/sqrt(gamma) * (I - fo x fo)
    BOOST_AUTO_TPL (gammaf,  value (this->M_activationSpace, *M_Gammaf) );
    BOOST_AUTO_TPL (gammas,  value (this->M_activationSpace, *M_Gammas) );
    BOOST_AUTO_TPL (gamman,  value (this->M_activationSpace, *M_Gamman) );

    MatrixSmall<3, 3> Id;
    Id (0, 0) = 1.;
    Id (0, 1) = 0., Id (0, 2) = 0.;
    Id (1, 0) = 0.;
    Id (1, 1) = 1., Id (1, 2) = 0.;
    Id (2, 0) = 0.;
    Id (2, 1) = 0., Id (2, 2) = 1.;

    // Strain-energy terms
    boost::shared_ptr<StrainEnergyGHO::dW1> dW1fun (new StrainEnergyGHO::dW1 (M_aiso, M_biso) );
    boost::shared_ptr<StrainEnergyGHO::d2W1> d2W1fun (new StrainEnergyGHO::d2W1 (M_aiso, M_biso) );

    boost::shared_ptr<StrainEnergyGHO::dW4> dW4ffun (new StrainEnergyGHO::dW4 (M_af, M_bf) );
    boost::shared_ptr<StrainEnergyGHO::d2W4> d2W4ffun (new StrainEnergyGHO::d2W4 (M_af, M_bf) );

    boost::shared_ptr<StrainEnergyGHO::dW4> dW4sfun (new StrainEnergyGHO::dW4 (M_as, M_bs) );
    boost::shared_ptr<StrainEnergyGHO::d2W4> d2W4sfun (new StrainEnergyGHO::d2W4 (M_as, M_bs) );

    boost::shared_ptr<StrainEnergyGHO::dW8> dW8fsfun (new StrainEnergyGHO::dW8 (M_afs, M_bfs) );
    boost::shared_ptr<StrainEnergyGHO::d2W8> d2W8fsfun (new StrainEnergyGHO::d2W8 (M_afs, M_bfs) );

    boost::shared_ptr<StrainEnergyGHO::dWvol> dWvolfun (new StrainEnergyGHO::dWvol (M_kappa) );
    boost::shared_ptr<StrainEnergyGHO::d2Wvol> d2Wvolfun (new StrainEnergyGHO::d2Wvol (M_kappa) );

    boost::shared_ptr<StrainEnergyGHO::GAHOShowValue> sv (new StrainEnergyGHO::GAHOShowValue() );

    boost::shared_ptr<StrainEnergyGHO::Normalize0> normalize0 (new StrainEnergyGHO::Normalize0() );
    boost::shared_ptr<StrainEnergyGHO::Normalize1> normalize1 (new StrainEnergyGHO::Normalize1() );
    boost::shared_ptr<StrainEnergyGHO::Normalize2> normalize2 (new StrainEnergyGHO::Normalize2() );

    // Kinematics
    BOOST_AUTO_TPL (I,      value (Id) );
    BOOST_AUTO_TPL (Grad_u, grad (this->M_dispETFESpace, disp, this->M_offset) );
    BOOST_AUTO_TPL (F,      Grad_u + I);
    BOOST_AUTO_TPL (FmT,    minusT (F) );
    BOOST_AUTO_TPL (J,      det (F) );
    // Fibres
    BOOST_AUTO_TPL (f01,    value (this->M_dispETFESpace, *M_fiberVector) );
    //BOOST_AUTO_TPL(f0,     f01 / sqrt( dot(f01, f01) ) );
    BOOST_AUTO_TPL (f0,     eval (normalize0, f01) );

    BOOST_AUTO_TPL (s01,     value (this->M_dispETFESpace, *M_sheetVector) );
    BOOST_AUTO_TPL (s02,     eval (normalize1, s01) );
    BOOST_AUTO_TPL (s0,      eval (normalize2, s02 - dot (f0, s02) * f0) );

    // Invariants
    BOOST_AUTO_TPL (I1,     dot (F, F) );
    BOOST_AUTO_TPL (C,      transpose (F) * F);
    BOOST_AUTO_TPL (I2,    value (0.5) * (I1 * I1 - dot (C, C) ) );


    BOOST_AUTO_TPL (I4f,    dot (F * f0, F * f0) );
    BOOST_AUTO_TPL (I4s,    dot (F * s0, F * s0) );
    BOOST_AUTO_TPL (I5f,    dot ( C * f0, C * f0) );
    BOOST_AUTO_TPL (I5s,    dot ( C * s0, C * s0) );

    BOOST_AUTO_TPL (I8fs,   dot (F * f0, F * s0) );
    //    BOOST_AUTO_TPL(I8fs2,   I2 + I4f * I4s + I5f + I5s - I1 * ( I4f + I4s ) );
    //    BOOST_AUTO_TPL(I8fs,    sqrt ( I8fs2 ) );
    // Reduced invariants
    BOOST_AUTO_TPL (Jm23,    pow (J, -2. / 3.0) );
    BOOST_AUTO_TPL (I1iso,   Jm23 * I1);
    BOOST_AUTO_TPL (I4fiso,  Jm23 * I4f);
    BOOST_AUTO_TPL (I4siso,  Jm23 * I4s);
    BOOST_AUTO_TPL (I8fsiso, Jm23 * I8fs);
    // Generalised invariants
    //BOOST_AUTO_TPL(gammac, pow(gamma, -2) - gamma);
    BOOST_AUTO_TPL (dI1edI1,   value (1.0) - gamman * ( gamman + value (2.0) ) * pow (gamman + value (1.0), -2.0) );
    BOOST_AUTO_TPL (dI1edI4f,   gamman * ( gamman + value (2.0) ) * pow (gamman + value (1.0), -2.0) - gammaf * ( gammaf + value (2.0) ) * pow (gammaf + value (1.0), -2.0) );
    BOOST_AUTO_TPL (dI1edI4s,   gamman * ( gamman + value (2.0) ) * pow (gamman + value (1.0), -2.0) - gammas * ( gammas + value (2.0) ) * pow (gammas + value (1.0), -2.0) );
    BOOST_AUTO_TPL (dI4fedI4f,   value (1.0) / ( ( gammaf + value (1.0) ) * ( gammaf + value (1.0) ) ) );
    BOOST_AUTO_TPL (dI4sedI4s,   value (1.0) / ( ( gammas + value (1.0) ) * ( gammas + value (1.0) ) ) );
    BOOST_AUTO_TPL (dI8fsedI8fs, value (1.0) / ( ( gammaf + value (1.0) ) * ( gammas + value (1.0) ) ) );


    BOOST_AUTO_TPL (I1eiso,   dI1edI1 * I1iso + dI1edI4f * I4fiso + dI1edI4s * I4siso);
    BOOST_AUTO_TPL (I4feiso,  dI4fedI4f *  I4fiso);
    BOOST_AUTO_TPL (I4seiso,  dI4sedI4s * I4siso);
    BOOST_AUTO_TPL (I8fseiso, dI8fsedI8fs * I8fsiso);

    // Strain-energy derivatives
    BOOST_AUTO_TPL (dW1,    value (1.0 - M_contractileFraction) * eval (dW1fun, I1iso) );
    BOOST_AUTO_TPL (dW4f,   value (1.0 - M_contractileFraction) * eval (dW4ffun, I4fiso) );
    BOOST_AUTO_TPL (dW4s,   value (1.0 - M_contractileFraction) * eval (dW4sfun, I4siso) );
    BOOST_AUTO_TPL (dW8fs,  value (1.0 - M_contractileFraction) * eval (dW8fsfun, I8fsiso) );

    // Strain-energy derivatives
    BOOST_AUTO_TPL (dW1e,    value (M_contractileFraction) * eval (dW1fun, I1eiso) );
    BOOST_AUTO_TPL (dW4fe,   /* value(M_contractileFraction) */ eval (dW4ffun, I4feiso) );
    BOOST_AUTO_TPL (dW4se,   /* value(M_contractileFraction) */ eval (dW4sfun, I4seiso) );
    BOOST_AUTO_TPL (dW8fse,  /* value(M_contractileFraction) */ eval (dW8fsfun, I8fseiso) );

    BOOST_AUTO_TPL (dWvol,  eval (dWvolfun, J) );
    BOOST_AUTO_TPL (d2Wvol, eval (d2Wvolfun, J) );

    BOOST_AUTO_TPL (dF,   grad (phi_j) );
    BOOST_AUTO_TPL (dJ,   J * FmT);
    BOOST_AUTO_TPL (dFmT, value (-1.0) * FmT * transpose (dF) * FmT);

    BOOST_AUTO_TPL (dI1iso,   value (2.0) * Jm23 * (F - I1 * value (1.0 / 3.0) * FmT ) );
    BOOST_AUTO_TPL (dI4fiso,  value (2.0) * Jm23 * (F * outerProduct (f0, f0) - I4f * value (1.0 / 3.0) * FmT ) );
    BOOST_AUTO_TPL (dI4siso,  value (2.0) * Jm23 * (F * outerProduct (s0, s0) - I4s * value (1.0 / 3.0) * FmT ) );
    BOOST_AUTO_TPL (dI8fsiso,  value (2.0) * Jm23 * ( value (0.5) * F * ( outerProduct (f0, s0) + outerProduct (s0, f0) ) - I8fs * value (1.0 / 3.0) * FmT ) );

    BOOST_AUTO_TPL (dI1eiso,  dI1edI1 * dI1iso + dI1edI4f * dI4fiso + dI1edI4s * dI4siso );
    BOOST_AUTO_TPL (dI4feiso, dI4fedI4f * dI4fiso );
    BOOST_AUTO_TPL (dI4seiso, dI4sedI4s * dI4siso );
    BOOST_AUTO_TPL (dI8fseiso, dI8fsedI8fs * dI8fsiso );

    BOOST_AUTO_TPL (Peiso_1,   dW1e * dI1eiso);
    BOOST_AUTO_TPL (Peiso_4f,  dW4fe * dI4feiso);
    BOOST_AUTO_TPL (Peiso_4s,  dW4se * dI4seiso);
    BOOST_AUTO_TPL (Peiso_8fs, dW8fse * dI8fseiso);

    BOOST_AUTO_TPL (Piso_1,   dW1 * dI1iso);
    BOOST_AUTO_TPL (Piso_4f,  dW4f * dI4fiso);
    BOOST_AUTO_TPL (Piso_4s,  dW4s * dI4siso);
    BOOST_AUTO_TPL (Piso_8fs, dW8fs * dI8fsiso);

    // Volumetric part
    BOOST_AUTO_TPL (Pvol,  dWvol * dJ);


    // Assemble residual
    integrate ( elements ( this->M_dispETFESpace->mesh() ),
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                dot ( Piso_1 + Peiso_1 + Pvol, grad (phi_i) )
              ) >> M_stiff;
    // Assemble residual
    if (M_af > 0)
    {
        integrate ( elements ( this->M_dispETFESpace->mesh() ),
                    this->M_dispFESpace->qr(),
                    this->M_dispETFESpace,
                    dot ( /*Piso_4f + */Peiso_4f , grad (phi_i) )
                  ) >> M_stiff;
    }
    // Assemble residual
    if (M_sheetVector)
    {
        if (M_as > 0)
        {

            integrate ( elements ( this->M_dispETFESpace->mesh() ),
                        this->M_dispFESpace->qr(),
                        this->M_dispETFESpace,
                        dot ( /* Piso_4s +*/ Peiso_4s , grad (phi_i) )
                      ) >> M_stiff;
        }
        // Assemble residual
        if (M_afs > 0)
        {
            integrate ( elements ( this->M_dispETFESpace->mesh() ),
                        this->M_dispFESpace->qr(),
                        this->M_dispETFESpace,
                        dot ( /*Piso_8fs +*/ /*eval(sv, I8fs2)  */ /* eval(sv, I2) *//* eval(sv, I8fs2) */ Peiso_8fs, grad (phi_i) )
                      ) >> M_stiff;
        }
    }
    this->M_stiff->globalAssemble();
}

template <typename MeshType>
void EMGeneralizedActiveHolzapfelOgdenMaterial<MeshType>::showMe ( std::string const& fileNameStiff,
                                                                 std::string const& fileNameJacobian)
{
    this->M_stiff->spy (fileNameStiff);
    this->M_jacobian->spy (fileNameJacobian);
}


template <typename MeshType>
void EMGeneralizedActiveHolzapfelOgdenMaterial<MeshType>::showMyParameters ()
{
    std::cout << "\n\n********************************************";
    std::cout << "\na: " << M_aiso;
    std::cout << "\nb: " << M_biso;
    std::cout << "\naf: " << M_af;
    std::cout << "\nbf: " << M_bf;
    std::cout << "\nas: " << M_as;
    std::cout << "\nbs: " << M_bs;
    std::cout << "\nafs: " << M_afs;
    std::cout << "\nbfs: " << M_bfs;
    std::cout << "\ncontractile fraction: " << M_contractileFraction;
    std::cout << "\n********************************************\n\n";
}

template <typename MeshType>
void EMGeneralizedActiveHolzapfelOgdenMaterial<MeshType>::computeLocalFirstPiolaKirchhoffTensor ( Epetra_SerialDenseMatrix& /*firstPiola*/,
        const Epetra_SerialDenseMatrix& /*tensorF*/,
        const Epetra_SerialDenseMatrix& /*cofactorF*/,
        const std::vector<Real>& /*invariants*/,
        const UInt /*marker*/)
{
    assert ("Not implemented yet for Holzapfel-Ogden material!");
}


template <typename MeshType>
inline EMActiveStructuralConstitutiveLaw<MeshType>* createGeneralizedActiveHolzapfelOgdenMaterial()
{
    return new EMGeneralizedActiveHolzapfelOgdenMaterial<MeshType>();
}
namespace
{
static bool registerGHO = EMActiveStructuralConstitutiveLaw<LifeV::RegionMesh<LinearTetra> >::StructureMaterialFactory::instance().registerProduct ("GAHO", &createGeneralizedActiveHolzapfelOgdenMaterial<LifeV::RegionMesh<LinearTetra> > );
}






} //Namespace LifeV

#endif /* __GeneralizedActiveHolzapfelOgdenMaterial_H */
