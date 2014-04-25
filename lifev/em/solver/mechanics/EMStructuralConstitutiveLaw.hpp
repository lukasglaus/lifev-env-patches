//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
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

/*!
 *  @file
 *  @brief This file contains an abstract class to implement different kinds of materials for structural dynamic problems (St. Venant-Kirchhoff, Neo-Hookean and Exponential materials right now )
 *
 *  @version 1.0
 *  @date 04-2014
 *  @author Simone Rossi <simone.rossi@epfl.ch>
 */

#ifndef _EMACTIVESTRUCTURALCONSTITUTIVELAW_H_
#define _EMACTIVESTRUCTURALCONSTITUTIVELAW_H_ 1


#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>

namespace LifeV
{
/*!
  \class StructuralConstitutiveLaw
  \brief
  This class is an abstract class to define different type of models for the arterial wall.
  This class has just pure virtual methods. They are implemented in the specific class for one material
*/

template <typename MeshType>
class EMActiveStructuralConstitutiveLaw : public StructuralConstitutiveLaw<MeshType>
{
public:

    //!@name Type definitions
    //@{
    typedef StructuralConstitutiveLawData          data_Type;
    typedef StructuralConstitutiveLaw<MeshType>    super;
    typedef MatrixEpetra<Real>            matrix_Type;
    typedef boost::shared_ptr<matrix_Type>         matrixPtr_Type;
    typedef VectorEpetra           vector_Type;
    typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;

    typedef FactorySingleton<Factory<EMActiveStructuralConstitutiveLaw<MeshType>, std::string> >  StructureMaterialFactory;


    typedef ETFESpace<MeshType, MapEpetra, 3, 3 >         ETFESpace_Type;
    typedef boost::shared_ptr<ETFESpace_Type>             ETFESpacePtr_Type;

    typedef FESpace< MeshType, MapEpetra >                FESpace_Type;
    typedef boost::shared_ptr<FESpace_Type>               FESpacePtr_Type;

    typedef MeshType                                        mesh_Type;
    typedef ETFESpace< mesh_Type, MapEpetra, 3, 1 >                        scalarETFESpace_Type;
    typedef boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > >    scalarETFESpacePtr_Type;

    //@}



    //! @name Constructor &  Deconstructor
    //@{

    EMActiveStructuralConstitutiveLaw() : super() {}

    virtual ~EMActiveStructuralConstitutiveLaw() {}

    //@}


    //SOME METHODS FOR ACTIVATED MATERIALS
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



    inline  vectorPtr_Type const fiberVectorPtr() const
    {
        return M_fiberVector;
    }

    inline  vectorPtr_Type fiberVectorPtr()
    {
        return M_fiberVector;
    }

    inline  vectorPtr_Type const sheetVectorPtr() const
    {
        return M_sheetVector;
    }

    inline  vectorPtr_Type sheetVectorPtr()
    {
        return M_sheetVector;
    }


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

    inline  scalarETFESpacePtr_Type activationSpace()
    {
        return M_activationSpace;
    }


    inline virtual MatrixSmall<3, 3>& identity()
    {
        MatrixSmall<3, 3> I;
        I (0, 0) = 1.0;
        I (0, 1) = 0.0;
        I (0, 2) = 0.0;
        I (1, 0) = 0.0;
        I (1, 1) = 1.0;
        I (1, 2) = 0.0;
        I (2, 0) = 0.0;
        I (2, 1) = 0.0;
        I (2, 2) = 1.0;
        return I;
    }



    inline virtual vectorPtr_Type gammaf()
    {
        vectorPtr_Type k;
        k.reset ( new vector_Type ( this -> M_dispFESpace -> map() ) );
        return k;
    }
    inline virtual vectorPtr_Type gammas()
    {
        vectorPtr_Type k;
        k.reset ( new vector_Type ( this -> M_dispFESpace -> map() ) );
        return k;
    }
    inline virtual vectorPtr_Type gamman()
    {
        vectorPtr_Type k;
        k.reset ( new vector_Type ( this -> M_dispFESpace -> map() ) );
        return k;
    }

    inline virtual void setGammaf (const vector_Type& /*gammaf*/) {}
    inline virtual void setGammas (const vector_Type& /*gammas*/) {}
    inline virtual void setGamman (const vector_Type& /*gamman*/) {}
    //@}

protected:
    vectorPtr_Type                      M_fiberVector;
    vectorPtr_Type                      M_sheetVector;
    scalarETFESpacePtr_Type             M_activationSpace;

};


}
#endif /*_STRUCTURALMATERIAL_H*/
