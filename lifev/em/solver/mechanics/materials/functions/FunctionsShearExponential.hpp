/*
 * EMMaterialFunctions.hpp
 *
 *  Created on: 28/apr/2014
 *      Author: srossi
 */

#ifndef FUNCTIONSSHEAREXPONENTIAL_HPP_
#define FUNCTIONSSHEAREXPONENTIAL_HPP_

#include <lifev/em/util/EMUtility.hpp>

//#include <lifev/em/solver/mechanics/EMETAAssembler.hpp>
#include <lifev/em/solver/mechanics/materials/functions/EMMaterialFunctions.hpp>

//using namespace LifeV;

//namespace MaterialFunctions
namespace LifeV
{

namespace MaterialFunctions
{

typedef VectorEpetra           vector_Type;
typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;

typedef MatrixEpetra<Real>           matrix_Type;
typedef boost::shared_ptr<matrix_Type>         matrixPtr_Type;

////////////////////////////////////////////////////////////////////////
//  SHEAR EXPONENTIAL FUNCTIONS
////////////////////////////////////////////////////////////////////////

template <class Mesh>
class ShearExponential : public virtual EMMaterialFunctions<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef EMData          data_Type;

    void showValue(Real v, std::string name = "variable")
    {
    	std::cout << name << " : " << v << std::endl;
    }

    return_Type operator() (const VectorSmall<3>& f, const VectorSmall<3>& s)
    {
        LifeV::Real I8fs = f.dot (s);
		return M_a * I8fs * std::exp ( M_b * I8fs * I8fs );
    }

    return_Type operator() (const Real& I8)
    {
        return M_a * I8 * std::exp ( M_b * I8 * I8 );
    }


    //    ShearExponential() : M_a(4170), M_b(11.602) {} // 0.33 KPa
    ShearExponential (Real a = 4170, Real b = 11.602) : M_a (a), M_b (b) {} // 0.33 KPa
    ShearExponential (const ShearExponential& ShearExponential)
    {
        M_a = ShearExponential.M_a;
        M_b = ShearExponential.M_b;
    }
    virtual ~ShearExponential() {}

    virtual void computeJacobian ( const vector_Type& disp,
                                  boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                  const vector_Type& fibers,
                                  const vector_Type& sheets,
                                  matrixPtr_Type           jacobianPtr)
    {
        EMAssembler::computeI8JacobianTerms (disp, dispETFESpace, fibers, sheets, jacobianPtr, this->getMe() );
    }

    virtual void computeResidual ( const vector_Type& disp,
                                  boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                  const vector_Type& fibers,
                                  const vector_Type& sheets,
                                  vectorPtr_Type           residualVectorPtr)
    {
        EMAssembler::computeI8ResidualTerms (disp, dispETFESpace, fibers, sheets, residualVectorPtr, this->getMe() );
    }

    void showMe()
    {
        std::cout << "Shear Exponential Function\n";
        std::cout << "Coefficient a: " << M_a;
        std::cout << ", coefficient b: " << M_b << "\n";
    }

    void setParameters (data_Type& data)
    {
		M_a = data.solidParameter<Real>("afs");
		M_b = data.solidParameter<Real>("bfs");
    }

private:
    Real M_a;
    Real M_b;
};



template <class Mesh>
class dShearExponential : public virtual EMMaterialFunctions<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    return_Type operator() (const VectorSmall<3>& f, const VectorSmall<3>& s)
    {
        LifeV::Real I8 = f.dot (s);
        return M_a * std::exp ( M_b * I8 * I8 ) * ( 2.0 * M_b * I8 * I8 + 1.0 );
    }

    return_Type operator() (const Real& I8)
    {
        return M_a * std::exp ( M_b * I8 * I8 ) * ( 2.0 * M_b * I8 * I8 + 1.0 );
    }
    //    ShearExponential() : M_a(4170), M_b(11.602) {} // 0.33 KPa
    dShearExponential (Real a = 4170, Real b = 11.602) : M_a (a), M_b (b) {} // 0.33 KPa
    dShearExponential (const dShearExponential& dShearExponential)
    {
        M_a = dShearExponential.M_a;
        M_b = dShearExponential.M_b;
    }
    virtual ~dShearExponential() {}

    void computeResidual ( const vector_Type& disp,
                                  boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                  const vector_Type& fibers,
                                  const vector_Type& sheets,
                                  vectorPtr_Type           residualVectorPtr)
    {
    }

    void computeJacobian ( const vector_Type& disp,
                                  boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                  const vector_Type& fibers,
                                  const vector_Type& sheets,
                                  matrixPtr_Type           jacobianPtr)
    {
        EMAssembler::computeI8JacobianTermsSecondDerivative (disp, dispETFESpace, fibers, sheets, jacobianPtr, this->getMe() );
    }

    void showMe()
    {
        std::cout << "Derivative Shear Exponential Function\n";
        std::cout << "Coefficient a: " << M_a;
        std::cout << ", coefficient b: " << M_b << "\n";
    }

    typedef EMData data_Type;
    void setParameters (data_Type& data)
    {
		M_a = data.solidParameter<Real>("afs");
		M_b = data.solidParameter<Real>("bfs");
    }

private:
    Real M_a;
    Real M_b;
};




} //EMMaterialFunctions

} //LifeV
#endif /* EMMATERIALFUNCTIONS_HPP_ */
