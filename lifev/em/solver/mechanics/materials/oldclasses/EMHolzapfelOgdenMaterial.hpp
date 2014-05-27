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
 *  @brief Holzapfel-Ogden material definition.
 *
 *  @date   05-08-2013
 *  @author Simone Pezzuto <simone.pezzuto@mail.polimi.it>
**/

#ifndef _HOLZAPFELOGDENMATERIAL_H_
#define _HOLZAPFELOGDENMATERIAL_H_

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <lifev/em/solver/EMActiveStructuralConstitutiveLaw.hpp>
#include <lifev/electrophysiology/util/HeartUtility.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <boost/typeof/typeof.hpp>

namespace
{

class PositivePartFct
{
public:
    typedef LifeV::Real return_Type;

    return_Type operator() (const LifeV::Real& p)
    {
        return (p > 0.0) ? p : 0.0;
    }

    PositivePartFct() {}
    PositivePartFct (const PositivePartFct&) {}
    ~PositivePartFct() {}
};

// Strain-energy density function
// ------------------------------
namespace StrainEnergyHO
{

// Isotropic Part
class dW1
{
public:
    typedef LifeV::Real return_Type;

    return_Type operator() (const LifeV::Real& I1)
    {
        return M_a / 2. * std::exp (M_b * (I1 - 3) );
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
        return M_a * M_b / 2. * std::exp (M_b * (I1 - 3) );
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
        return M_a * (I4 >= 1.)
               * (I4 - 1)
               * std::exp (M_b * pow (I4 - 1, 2) );
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
        return M_a * (I4 >= 1.)
               * (1.0 + 2.0 * M_b * pow (I4 - 1, 2) )
               * std::exp (M_b * pow (I4 - 1, 2) );
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
        return M_a * I8 * std::exp (M_b * pow (I8, 2) );
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
        return  M_a * (1.0 + 2.0 * M_b * pow (I8, 2) )
                * std::exp (M_b * pow (I8, 2) );
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
class EMHolzapfelOgdenMaterial :
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
    Real M_kappa;

    //! @name Constructor &  Destructor
    //@{

    EMHolzapfelOgdenMaterial();

    virtual  ~EMHolzapfelOgdenMaterial();

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

    inline void setFiberVector ( const vector_Type& fiberVector)
    {
        M_fiberVector.reset ( new vector_Type ( fiberVector ) );
    }

    inline void setSheetVector ( const vector_Type& sheetVector)
    {
        M_sheetVector.reset ( new vector_Type ( sheetVector ) );
    }

    //@}

    //! @name Methods
    //@{

    inline void setupFiberVector ( std::string& name, boost::shared_ptr<mesh_Type> mesh )
    {
        ElectrophysiologyUtility::importFibers ( M_fiberVector, name, mesh  );
    }

    inline void setupFiberVector ( std::string& name, std::string& path )
    {
        ElectrophysiologyUtility::importFibers ( M_fiberVector, name, path  );
    }

    void setupFiberVector ( Real& fx, Real& fy, Real& fz )
    {
        ElectrophysiologyUtility::setupFibers ( *M_fiberVector, fx, fy, fz  );
    }

    inline void setupSheetVector ( std::string& name, boost::shared_ptr<mesh_Type> mesh )
    {
        ElectrophysiologyUtility::importFibers ( M_sheetVector, name, mesh  );
    }

    inline void setupSheetVector ( std::string& name, std::string& path )
    {
        ElectrophysiologyUtility::importFibers ( M_sheetVector, name, path  );
    }

    void setupSheetVector ( Real& sx, Real& sy, Real& sz )
    {
        ElectrophysiologyUtility::setupFibers ( *M_sheetVector, sx, sy, sz);
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
    scalarETFESpacePtr_Type             M_activationSpace;

};

template <typename MeshType>
EMHolzapfelOgdenMaterial<MeshType>::EMHolzapfelOgdenMaterial() :
    super           ( ),
    M_stiff         ( ),
    M_Gammaf        ( ),
    M_fiberVector   ( ),
    M_sheetVector   ( ),
    M_activationSpace( )
{
}





template <typename MeshType>
EMHolzapfelOgdenMaterial<MeshType>::~EMHolzapfelOgdenMaterial()
{}


template <typename MeshType>
void
EMHolzapfelOgdenMaterial<MeshType>::setup ( const FESpacePtr_Type& dFESpace,
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



    //   this->setDefaultParams();

    // The 2 is because the law uses three parameters (mu, bulk).
    // another way would be to set up the number of constitutive parameters of the law
    // in the data file to get the right size. Note the comment below.
    this->M_vectorsParameters.reset ( new vectorsParameters_Type ( 2 ) );

    this->setupVectorsParameters();

}

template <typename MeshType>
void EMHolzapfelOgdenMaterial<MeshType>::setDefaultParams()
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
}

template <typename MeshType>
void
EMHolzapfelOgdenMaterial<MeshType>::setup ( const FESpacePtr_Type&                      dFESpace,
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


    //    this->setDefaultParams();

    // This parameters are not used, but we are keeping them to avoid
    // inconsistencies among different materials
    this->M_vectorsParameters.reset ( new vectorsParameters_Type (2) );
    this->setupVectorsParameters();
}


template <typename MeshType>
void
EMHolzapfelOgdenMaterial<MeshType>::setup ( const FESpacePtr_Type&                      dFESpace,
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
EMHolzapfelOgdenMaterial<MeshType>::setup ( const FESpacePtr_Type&                      dFESpace,
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
EMHolzapfelOgdenMaterial<MeshType>::setupVectorsParameters ( void )
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
        M_as    = this->M_dataMaterial->Bs ( markerID );
        M_bfs   = this->M_dataMaterial->Bfs ( markerID );
        M_kappa = this->M_dataMaterial->bulk ( markerID );
    }
}


template <typename MeshType>
void EMHolzapfelOgdenMaterial<MeshType>::computeLinearStiff (dataPtr_Type& /*dataMaterial*/,
                                                             const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                                                             const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/)
{
    //! Empty method for neo-hookean material
}


template <typename MeshType>
void EMHolzapfelOgdenMaterial<MeshType>::updateJacobianMatrix ( const vector_Type&       disp,
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
void EMHolzapfelOgdenMaterial<MeshType>::updateNonLinearJacobianTerms ( matrixPtr_Type&       jacobian,
                                                                        const vector_Type&    disp,
                                                                        const dataPtr_Type&   dataMaterial,
                                                                        const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                                        const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                                                        const displayerPtr_Type&  displayer )
{
    {
        using namespace ExpressionAssembly;

        displayer->leaderPrint ("   Non-Linear S-  updating non linear terms in the Jacobian Matrix (Holzapfel-Ogden)");

        * (jacobian) *= 0.0;

        // Activation is Fa = gamma (fo x fo) + 1/sqrt(gamma) * (I - fo x fo)
        BOOST_AUTO_TPL (gamma,  value (1.0) + value (this->M_activationSpace, *M_Gammaf) );

        MatrixSmall<3, 3> Id;
        Id (0, 0) = 1.;
        Id (0, 1) = 0., Id (0, 2) = 0.;
        Id (1, 0) = 0.;
        Id (1, 1) = 1., Id (1, 2) = 0.;
        Id (2, 0) = 0.;
        Id (2, 1) = 0., Id (2, 2) = 1.;

        // Strain-energy terms
        boost::shared_ptr<StrainEnergyHO::dW1> dW1fun (new StrainEnergyHO::dW1 (M_aiso, M_biso) );
        boost::shared_ptr<StrainEnergyHO::d2W1> d2W1fun (new StrainEnergyHO::d2W1 (M_aiso, M_biso) );

        boost::shared_ptr<StrainEnergyHO::dW4> dW4ffun (new StrainEnergyHO::dW4 (M_af, M_bf) );
        boost::shared_ptr<StrainEnergyHO::d2W4> d2W4ffun (new StrainEnergyHO::d2W4 (M_af, M_bf) );

        boost::shared_ptr<StrainEnergyHO::dW4> dW4sfun (new StrainEnergyHO::dW4 (M_as, M_bs) );
        boost::shared_ptr<StrainEnergyHO::d2W4> d2W4sfun (new StrainEnergyHO::d2W4 (M_as, M_bs) );

        boost::shared_ptr<StrainEnergyHO::dW8> dW8fsfun (new StrainEnergyHO::dW8 (M_afs, M_bfs) );
        boost::shared_ptr<StrainEnergyHO::d2W8> d2W8fsfun (new StrainEnergyHO::d2W8 (M_afs, M_bfs) );

        boost::shared_ptr<StrainEnergyHO::dWvol> dWvolfun (new StrainEnergyHO::dWvol (M_kappa) );
        boost::shared_ptr<StrainEnergyHO::d2Wvol> d2Wvolfun (new StrainEnergyHO::d2Wvol (M_kappa) );

        // Kinematics
        BOOST_AUTO_TPL (I,      value (Id) );
        BOOST_AUTO_TPL (Grad_u, grad (this->M_dispETFESpace, disp, this->M_offset) );
        BOOST_AUTO_TPL (F,      Grad_u + I);
        BOOST_AUTO_TPL (FmT,    minusT (F) );
        BOOST_AUTO_TPL (J,      det (F) );
        // Fibres
        BOOST_AUTO_TPL (f0,     value (this->M_dispETFESpace, *M_fiberVector) );
        BOOST_AUTO_TPL (s0,     value (this->M_dispETFESpace, *M_sheetVector) );
        // Invariants
        BOOST_AUTO_TPL (I1,     dot (F, F) );
        BOOST_AUTO_TPL (I4f,    dot (F * f0, F * f0) );
        BOOST_AUTO_TPL (I4s,    dot (F * s0, F * s0) );
        BOOST_AUTO_TPL (I8fs,   dot (F * f0, F * s0) );
        // Reduced invariants
        BOOST_AUTO_TPL (Jm23,    pow (J, -2. / 3) );
        BOOST_AUTO_TPL (I1iso,   Jm23 * I1);
        BOOST_AUTO_TPL (I4fiso,  Jm23 * I4f);
        BOOST_AUTO_TPL (I4siso,  Jm23 * I4s);
        BOOST_AUTO_TPL (I8fsiso, Jm23 * I8fs);
        // Generalised invariants
        BOOST_AUTO_TPL (gammac, pow (gamma, -2) - gamma);

        BOOST_AUTO_TPL (I1eiso,   gamma * I1iso + gammac * I4fiso);
        BOOST_AUTO_TPL (I4feiso,  pow (gamma, -2) * I4fiso);
        BOOST_AUTO_TPL (I4seiso,  gamma * I4siso);
        BOOST_AUTO_TPL (I8fseiso, pow (gamma, -1. / 2) * I8fsiso);

        // Strain-energy derivatives
        BOOST_AUTO_TPL (dW1,    eval (dW1fun, I1eiso) );
        BOOST_AUTO_TPL (d2W1,   eval (d2W1fun, I1eiso) );
        BOOST_AUTO_TPL (dW4f,   eval (dW4ffun, I4feiso) );
        BOOST_AUTO_TPL (d2W4f,  eval (d2W4ffun, I4feiso) );
        BOOST_AUTO_TPL (dW4s,   eval (dW4sfun, I4seiso) );
        BOOST_AUTO_TPL (d2W4s,  eval (d2W4sfun, I4seiso) );
        BOOST_AUTO_TPL (dW8fs,  eval (dW8fsfun, I8fseiso) );
        BOOST_AUTO_TPL (d2W8fs, eval (d2W8fsfun, I8fseiso) );

        BOOST_AUTO_TPL (dWvol,  eval (dWvolfun, J) );
        BOOST_AUTO_TPL (d2Wvol, eval (d2Wvolfun, J) );

        // Deviatorics
        BOOST_AUTO_TPL (Fdev1,   F - value (1. / 3) *I1 * FmT);
        BOOST_AUTO_TPL (Fdev4f,  outerProduct (F * f0, f0) - value (1. / 3) *I4f * FmT);
        BOOST_AUTO_TPL (Fdev4s,  outerProduct (F * s0, s0) - value (1. / 3) *I4s * FmT);
        BOOST_AUTO_TPL (Fdev8fs, value (0.5) * (outerProduct (F * s0, f0) + outerProduct (F * f0, s0) ) - value (1. / 3) *I8fs * FmT);

        // Derivatives
        BOOST_AUTO_TPL (dF,   grad (phi_j) );
        BOOST_AUTO_TPL (dJ,   J * FmT);
        BOOST_AUTO_TPL (dFmT, value (-1.0) *FmT * transpose (dF) *FmT);

        BOOST_AUTO_TPL (dI1eiso,   value (2.) *Jm23 * (gamma * Fdev1 + gammac * Fdev4f) );
        BOOST_AUTO_TPL (dI4feiso,  value (2.) *Jm23 * pow (gamma, -2) *Fdev4f);
        BOOST_AUTO_TPL (dI4seiso,  value (2.) *Jm23 * gamma * Fdev4s);
        BOOST_AUTO_TPL (dI8fseiso, value (2.) *Jm23 * pow (gamma, -1. / 2) *Fdev8fs);

        BOOST_AUTO_TPL (dFdev1,   dF + value (-1. / 3) * (value (2.0) *dot (F, dF) *FmT + I1 * dFmT) );
        BOOST_AUTO_TPL (dFdev4f,  outerProduct (dF * f0, f0) - value (1. / 3) * (value (2.0) *dot (F * f0, dF * f0) *FmT + I4f * dFmT) );
        BOOST_AUTO_TPL (dFdev4s,  outerProduct (dF * s0, s0) - value (1. / 3) * (value (2.0) *dot (F * s0, dF * s0) *FmT + I4s * dFmT) );
        BOOST_AUTO_TPL (dFdev8fs,   value (0.5) * (outerProduct (dF * s0, f0) + outerProduct (dF * f0, s0) )
                        - value (1. / 3) * ( (dot (F * s0, dF * f0) + dot (F * f0, dF * s0) ) *FmT + I8fs * dFmT) );

        BOOST_AUTO_TPL (d2I1eiso,    value (2.) *Jm23 * (gamma * dFdev1 + gammac * dFdev4f
                                                         - value (2. / 3) *dot (FmT, dF) * (gamma * Fdev1 + gammac * Fdev4f) ) );
        BOOST_AUTO_TPL (d2I4feiso,   value (2.) *Jm23 * pow (gamma, -2) * (dFdev4f
                                                                           - value (2. / 3) *dot (FmT, dF) *Fdev4f) );
        BOOST_AUTO_TPL (d2I4seiso,   value (2.) *Jm23 * gamma * (dFdev4s
                                                                 - value (2. / 3) *dot (FmT, dF) *Fdev4s) );
        BOOST_AUTO_TPL (d2I8fseiso,   value (2.) *Jm23 * pow (gamma, -1. / 2) * (dFdev8fs
                                                                                 - value (2. / 3) *dot (FmT, dF) *Fdev8fs) );

        // Isochoric part
        BOOST_AUTO_TPL (Piso_1,   dW1 * dI1eiso);
        BOOST_AUTO_TPL (Piso_4f,  dW4f * dI4feiso);
        BOOST_AUTO_TPL (Piso_4s,  dW4s * dI4seiso);
        BOOST_AUTO_TPL (Piso_8fs, dW8fs * dI8fseiso);

        BOOST_AUTO_TPL (dPiso_1,   d2W1 * dot (dI1eiso, dF) *dI1eiso + dW1 * d2I1eiso );
        BOOST_AUTO_TPL (dPiso_4f,  d2W4f * dot (dI4feiso, dF) *dI4feiso + dW4f * d2I4feiso );
        BOOST_AUTO_TPL (dPiso_4s,  d2W4s * dot (dI4seiso, dF) *dI4seiso + dW4s * d2I4seiso );
        BOOST_AUTO_TPL (dPiso_8fs, d2W8fs * dot (dI8fseiso, dF) *dI8fseiso + dW8fs * d2I8fseiso );

        // Volumetric part
        BOOST_AUTO_TPL (Pvol,  dWvol * dJ);
        BOOST_AUTO_TPL (dPvol, (d2Wvol + J * dWvol) *dot (FmT, dF) *dJ + dWvol * J * dFmT);

        // Assembly
        integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                    this->M_dispFESpace->qr(),
                    this->M_dispETFESpace,
                    this->M_dispETFESpace,
                    dot (dPiso_1 + dPvol, grad (phi_i) )
                  ) >> jacobian;
        if (M_af > 0)
        {
            integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                        this->M_dispFESpace->qr(),
                        this->M_dispETFESpace,
                        this->M_dispETFESpace,
                        dot ( dPiso_4f , grad (phi_i) )
                      ) >> jacobian;
        }
        if (M_as > 0)
        {
            integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                        this->M_dispFESpace->qr(),
                        this->M_dispETFESpace,
                        this->M_dispETFESpace,
                        dot ( dPiso_4s, grad (phi_i) )
                      ) >> jacobian;
        }
        if (M_afs > 0)
        {
            integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                        this->M_dispFESpace->qr(),
                        this->M_dispETFESpace,
                        this->M_dispETFESpace,
                        dot ( dPiso_8fs, grad (phi_i) )
                      ) >> jacobian;
        }

    }


    jacobian->globalAssemble();
}


template <typename MeshType>
void EMHolzapfelOgdenMaterial<MeshType>::apply ( const vector_Type& sol, vector_Type& res,
                                                 const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                 const mapMarkerIndexesPtr_Type mapsMarkerIndexes)
{
    computeStiffness (sol, 0., this->M_dataMaterial, mapsMarkerVolumes, mapsMarkerIndexes, this->M_displayer);
    res += *M_stiff;
}

template <typename MeshType>
void EMHolzapfelOgdenMaterial<MeshType>::computeStiffness ( const vector_Type&       disp,
                                                            Real                     /*factor*/,
                                                            const dataPtr_Type&      dataMaterial,
                                                            const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                            const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                                            const displayerPtr_Type& displayer )
{
    using namespace ExpressionAssembly;

    this->M_stiff.reset (new vector_Type (*this->M_localMap) );

    displayer->leaderPrint (" \n******************************************************************\n  ");
    displayer->leaderPrint (" Non-Linear S-  Computing the Holzapfel-Ogden nonlinear stiffness vector"     );
    displayer->leaderPrint (" \n******************************************************************\n  ");

    M_stiff.reset (new vector_Type (*this->M_localMap) );
    * (M_stiff) *= 0.0;

    // Activation is Fa = gamma (fo x fo) + 1/sqrt(gamma) * (I - fo x fo)
    BOOST_AUTO_TPL (gamma,  value (1.0) + value (this->M_activationSpace, *M_Gammaf) );

    MatrixSmall<3, 3> Id;
    Id (0, 0) = 1.;
    Id (0, 1) = 0., Id (0, 2) = 0.;
    Id (1, 0) = 0.;
    Id (1, 1) = 1., Id (1, 2) = 0.;
    Id (2, 0) = 0.;
    Id (2, 1) = 0., Id (2, 2) = 1.;

    // Strain-energy terms
    boost::shared_ptr<StrainEnergyHO::dW1> dW1fun (new StrainEnergyHO::dW1 (M_aiso, M_biso) );
    boost::shared_ptr<StrainEnergyHO::d2W1> d2W1fun (new StrainEnergyHO::d2W1 (M_aiso, M_biso) );

    boost::shared_ptr<StrainEnergyHO::dW4> dW4ffun (new StrainEnergyHO::dW4 (M_af, M_bf) );
    boost::shared_ptr<StrainEnergyHO::d2W4> d2W4ffun (new StrainEnergyHO::d2W4 (M_af, M_bf) );

    boost::shared_ptr<StrainEnergyHO::dW4> dW4sfun (new StrainEnergyHO::dW4 (M_as, M_bs) );
    boost::shared_ptr<StrainEnergyHO::d2W4> d2W4sfun (new StrainEnergyHO::d2W4 (M_as, M_bs) );

    boost::shared_ptr<StrainEnergyHO::dW8> dW8fsfun (new StrainEnergyHO::dW8 (M_afs, M_bfs) );
    boost::shared_ptr<StrainEnergyHO::d2W8> d2W8fsfun (new StrainEnergyHO::d2W8 (M_afs, M_bfs) );

    boost::shared_ptr<StrainEnergyHO::dWvol> dWvolfun (new StrainEnergyHO::dWvol (M_kappa) );
    boost::shared_ptr<StrainEnergyHO::d2Wvol> d2Wvolfun (new StrainEnergyHO::d2Wvol (M_kappa) );

    // Kinematics
    BOOST_AUTO_TPL (I,      value (Id) );
    BOOST_AUTO_TPL (Grad_u, grad (this->M_dispETFESpace, disp, this->M_offset) );
    BOOST_AUTO_TPL (F,      Grad_u + I);
    BOOST_AUTO_TPL (FmT,    minusT (F) );
    BOOST_AUTO_TPL (J,      det (F) );
    // Fibres
    BOOST_AUTO_TPL (f0,     value (this->M_dispETFESpace, *M_fiberVector) );
    BOOST_AUTO_TPL (s0,     value (this->M_dispETFESpace, *M_sheetVector) );
    // Invariants
    BOOST_AUTO_TPL (I1,     dot (F, F) );
    BOOST_AUTO_TPL (I4f,    dot (F * f0, F * f0) );
    BOOST_AUTO_TPL (I4s,    dot (F * s0, F * s0) );
    BOOST_AUTO_TPL (I8fs,   dot (F * f0, F * s0) );
    // Reduced invariants
    BOOST_AUTO_TPL (Jm23,    pow (J, -2. / 3) );
    BOOST_AUTO_TPL (I1iso,   Jm23 * I1);
    BOOST_AUTO_TPL (I4fiso,  Jm23 * I4f);
    BOOST_AUTO_TPL (I4siso,  Jm23 * I4s);
    BOOST_AUTO_TPL (I8fsiso, Jm23 * I8fs);
    // Generalised invariants
    BOOST_AUTO_TPL (gammac, pow (gamma, -2) - gamma);

    BOOST_AUTO_TPL (I1eiso,   gamma * I1iso + gammac * I4fiso);
    BOOST_AUTO_TPL (I4feiso,  pow (gamma, -2) * I4fiso);
    BOOST_AUTO_TPL (I4seiso,  gamma * I4siso);
    BOOST_AUTO_TPL (I8fseiso, pow (gamma, -1. / 2) * I8fsiso);

    // Strain-energy derivatives
    BOOST_AUTO_TPL (dW1,    eval (dW1fun, I1eiso) );
    BOOST_AUTO_TPL (d2W1,   eval (d2W1fun, I1eiso) );
    BOOST_AUTO_TPL (dW4f,   eval (dW4ffun, I4feiso) );
    BOOST_AUTO_TPL (d2W4f,  eval (d2W4ffun, I4feiso) );
    BOOST_AUTO_TPL (dW4s,   eval (dW4sfun, I4seiso) );
    BOOST_AUTO_TPL (d2W4s,  eval (d2W4sfun, I4seiso) );
    BOOST_AUTO_TPL (dW8fs,  eval (dW8fsfun, I8fseiso) );
    BOOST_AUTO_TPL (d2W8fs, eval (d2W8fsfun, I8fseiso) );

    BOOST_AUTO_TPL (dWvol,  eval (dWvolfun, J) );
    BOOST_AUTO_TPL (d2Wvol, eval (d2Wvolfun, J) );

    // Deviatorics
    BOOST_AUTO_TPL (Fdev1,   F - value (1. / 3) *I1 * FmT);
    BOOST_AUTO_TPL (Fdev4f,  outerProduct (F * f0, f0) - value (1. / 3) *I4f * FmT);
    BOOST_AUTO_TPL (Fdev4s,  outerProduct (F * s0, s0) - value (1. / 3) *I4s * FmT);
    BOOST_AUTO_TPL (Fdev8fs, value (0.5) * (outerProduct (F * s0, f0) + outerProduct (F * f0, s0) ) - value (1. / 3) *I8fs * FmT);

    // Derivatives
    BOOST_AUTO_TPL (dF,   grad (phi_j) );
    BOOST_AUTO_TPL (dJ,   J * FmT);
    BOOST_AUTO_TPL (dFmT, value (-1.0) *FmT * transpose (dF) *FmT);

    BOOST_AUTO_TPL (dI1eiso,   value (2.) *Jm23 * (gamma * Fdev1 + gammac * Fdev4f) );
    BOOST_AUTO_TPL (dI4feiso,  value (2.) *Jm23 * pow (gamma, -2) *Fdev4f);
    BOOST_AUTO_TPL (dI4seiso,  value (2.) *Jm23 * gamma * Fdev4s);
    BOOST_AUTO_TPL (dI8fseiso, value (2.) *Jm23 * pow (gamma, -1. / 2) *Fdev8fs);

    BOOST_AUTO_TPL (dFdev1,   dF + value (-1. / 3) * (value (2.0) *dot (F, dF) *FmT + I1 * dFmT) );
    BOOST_AUTO_TPL (dFdev4f,  outerProduct (dF * f0, f0) - value (1. / 3) * (value (2.0) *dot (F * f0, dF * f0) *FmT + I4f * dFmT) );
    BOOST_AUTO_TPL (dFdev4s,  outerProduct (dF * s0, s0) - value (1. / 3) * (value (2.0) *dot (F * s0, dF * s0) *FmT + I4s * dFmT) );
    BOOST_AUTO_TPL (dFdev8fs,   value (0.5) * (outerProduct (dF * s0, f0) + outerProduct (dF * f0, s0) )
                    - value (1. / 3) * ( (dot (F * s0, dF * f0) + dot (F * f0, dF * s0) ) *FmT + I8fs * dFmT) );

    BOOST_AUTO_TPL (d2I1eiso,    value (2.) *Jm23 * (gamma * dFdev1 + gammac * dFdev4f
                                                     - value (2. / 3) *dot (FmT, dF) * (gamma * Fdev1 + gammac * Fdev4f) ) );
    BOOST_AUTO_TPL (d2I4feiso,   value (2.) *Jm23 * pow (gamma, -2) * (dFdev4f
                                                                       - value (2. / 3) *dot (FmT, dF) *Fdev4f) );
    BOOST_AUTO_TPL (d2I4seiso,   value (2.) *Jm23 * gamma * (dFdev4s
                                                             - value (2. / 3) *dot (FmT, dF) *Fdev4s) );
    BOOST_AUTO_TPL (d2I8fseiso,   value (2.) *Jm23 * pow (gamma, -1. / 2) * (dFdev8fs
                                                                             - value (2. / 3) *dot (FmT, dF) *Fdev8fs) );

    // Isochoric part
    BOOST_AUTO_TPL (Piso_1,   dW1 * dI1eiso);
    BOOST_AUTO_TPL (Piso_4f,  dW4f * dI4feiso);
    BOOST_AUTO_TPL (Piso_4s,  dW4s * dI4seiso);
    BOOST_AUTO_TPL (Piso_8fs, dW8fs * dI8fseiso);

    BOOST_AUTO_TPL (dPiso_1,   d2W1 * dot (dI1eiso, dF) *dI1eiso + dW1 * d2I1eiso );
    BOOST_AUTO_TPL (dPiso_4f,  d2W4f * dot (dI4feiso, dF) *dI4feiso + dW4f * d2I4feiso );
    BOOST_AUTO_TPL (dPiso_4s,  d2W4s * dot (dI4seiso, dF) *dI4seiso + dW4s * d2I4seiso );
    BOOST_AUTO_TPL (dPiso_8fs, d2W8fs * dot (dI8fseiso, dF) *dI8fseiso + dW8fs * d2I8fseiso );

    // Volumetric part
    BOOST_AUTO_TPL (Pvol,  dWvol * dJ);
    BOOST_AUTO_TPL (dPvol, (d2Wvol + J * dWvol) *dot (FmT, dF) *dJ + dWvol * J * dFmT);

    // Assemble residual
    integrate ( elements ( this->M_dispETFESpace->mesh() ),
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                dot (Piso_1 + Pvol, grad (phi_i) )
              ) >> M_stiff;
    // Assemble residual
    if (M_af > 0)
    {
        integrate ( elements ( this->M_dispETFESpace->mesh() ),
                    this->M_dispFESpace->qr(),
                    this->M_dispETFESpace,
                    dot ( Piso_4f , grad (phi_i) )
                  ) >> M_stiff;
    }
    // Assemble residual
    if (M_as > 0)
    {

        integrate ( elements ( this->M_dispETFESpace->mesh() ),
                    this->M_dispFESpace->qr(),
                    this->M_dispETFESpace,
                    dot ( Piso_4s , grad (phi_i) )
                  ) >> M_stiff;
    }
    // Assemble residual
    if (M_afs > 0)
    {
        integrate ( elements ( this->M_dispETFESpace->mesh() ),
                    this->M_dispFESpace->qr(),
                    this->M_dispETFESpace,
                    dot ( Piso_8fs, grad (phi_i) )
                  ) >> M_stiff;
    }
    this->M_stiff->globalAssemble();
}

template <typename MeshType>
void EMHolzapfelOgdenMaterial<MeshType>::showMe ( std::string const& fileNameStiff,
                                                  std::string const& fileNameJacobian)
{
    this->M_stiff->spy (fileNameStiff);
    this->M_jacobian->spy (fileNameJacobian);
}

template <typename MeshType>
void EMHolzapfelOgdenMaterial<MeshType>::computeLocalFirstPiolaKirchhoffTensor ( Epetra_SerialDenseMatrix& /*firstPiola*/,
        const Epetra_SerialDenseMatrix& /*tensorF*/,
        const Epetra_SerialDenseMatrix& /*cofactorF*/,
        const std::vector<Real>& /*invariants*/,
        const UInt /*marker*/)
{
    assert ("Not implemented yet for Holzapfel-Ogden material!");
}


template <typename MeshType>
inline EMActiveStructuralConstitutiveLaw<MeshType>* createHolzapfelOgdenMaterial()
{
    return new EMHolzapfelOgdenMaterial<MeshType>();
}
namespace
{
static bool registerHO = EMActiveStructuralConstitutiveLaw<LifeV::RegionMesh<LinearTetra> >::StructureMaterialFactory::instance().registerProduct ("HolzapfelOgden", &createHolzapfelOgdenMaterial<LifeV::RegionMesh<LinearTetra> > );
}

} //Namespace LifeV

#endif /* __HOLZAPFELOGDENMATERIAL_H */
