/*
 * EMMaterialFunctions.hpp
 *
 *  Created on: 28/apr/2014
 *      Author: srossi
 */

#ifndef FUNCTIONSORTHOTROPICFUNG_HPP_
#define FUNCTIONSORTHOTROPICFUNG_HPP_

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
//  ORTHOTROPIC FUNG FUNCTIONS
////////////////////////////////////////////////////////////////////////

namespace FungFunctors
{

template <class Mesh>
class W
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    W (  Real a = 2000,        // 0.23 kPa
         Real bff = 35.7,
         Real bss = 18.9,
         Real bnn = 13.0,
         Real bfs = 10.5,
         Real bfn = 11.4,
         Real bsn = 8.35) :
        M_a (a),
        M_bff (bff),
        M_bss (bss),
        M_bnn (bnn),
        M_bfs (bfs),
        M_bfn (bfn),
        M_bsn (bsn){}

    W (const W& OrthotropicFung)
    {
        M_a = OrthotropicFung.M_a;
        M_bff = OrthotropicFung.M_bff;
        M_bss = OrthotropicFung.M_bss;
        M_bnn = OrthotropicFung.M_bnn;
        M_bfs = OrthotropicFung.M_bfs;
        M_bfn = OrthotropicFung.M_bfn;
        M_bsn = OrthotropicFung.M_bsn;
    }
    virtual ~W() {}

    Real computeQ (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
    	Real I4f = Elasticity::I4(F, f0);
    	Real I4fm1 = I4f - 1;
    	Real I4s = Elasticity::I4(F, s0);
    	Real I4sm1 = I4s - 1;
    	Real I1bar = Elasticity::I1bar(F);
//    	Real I1bar = Elasticity::I1(F);
    	Real I8fs = Elasticity::I8(F, f0,s0);
        VectorSmall<3> n0 = Elasticity::crossProduct(f0, s0);
    	Real I8fn = Elasticity::I8(F, f0, n0);
    	Real I8sn = Elasticity::I8(F, s0, n0);
		Real Q  = 0.25 * M_bff * I4fm1 * I4fm1;
		     Q += 0.25 * M_bss * I4sm1 * I4sm1;
    	     Q += 0.25 * M_bnn * (I1bar - I4f - I4s - 1.0) * (I1bar - I4f - I4s - 1.0);
    	     Q += 0.5 * M_bfs * I8fs * I8fs;
    	     Q += 0.5 * M_bfn * I8fn * I8fn;
    	     Q += 0.5 * M_bsn * I8sn * I8sn;
        	return Q;
    }

    Real computedQdI1 (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
    	Real I4f = Elasticity::I4(F, f0);
    	Real I4s = Elasticity::I4(F, s0);
    	Real I1bar = Elasticity::I1bar(F);
//    	Real I1bar = Elasticity::I1(F);
    	return 0.5 * M_bnn * (I1bar - I4f - I4s - 1.0);
    }
    Real computed2QdI12 (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        return 0.5 * M_bnn;
    }
    Real computed2QdI1dI4f (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        return - 0.5 * M_bnn;
    }
    Real computed2QdI1dI4s (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        return - 0.5 * M_bnn;
    }
    Real computed2QdI1dI8fs (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        return 0.0;
    }
    Real computed2QdI1dI8fn (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        return 0.0;
    }
    Real computed2QdI1dI8sn (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        return 0.0;
    }

    Real computedQdI4f (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
    	Real I4f = Elasticity::I4(F, f0);
    	Real I4fm1 = I4f - 1;
    	Real I4s = Elasticity::I4(F, s0);
    	Real I1bar = Elasticity::I1bar(F);
//    	Real I1bar = Elasticity::I1(F);
    	return 0.5 * ( M_bff * I4fm1  - M_bnn * (I1bar - I4f - I4s - 1.0) );
    }
    Real computed2QdI4f2 (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        return 0.5 * ( M_bff  + M_bnn );
    }
    Real computed2QdI4fdI4s (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        return 0.5 * M_bnn;
    }
    Real computed2QdI4fdI8fs (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        return 0.0;
    }
    Real computed2QdI4fdI8fn (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        return 0.0;
    }
    Real computed2QdI4fdI8sn (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        return 0.0;
    }


    Real computedQdI4s (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
    	Real I4f = Elasticity::I4(F, f0);
    	Real I4s = Elasticity::I4(F, s0);
    	Real I4sm1 = I4s - 1;
    	Real I1bar = Elasticity::I1bar(F);
//    	Real I1bar = Elasticity::I1(F);
    	return 0.5 * ( M_bss * I4sm1  - M_bnn * (I1bar - I4f - I4s - 1.0) );
    }
    Real computed2QdI4s2 (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        return 0.5 * ( M_bss + M_bnn );
    }
    Real computed2QdI4sdI8fs (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        return 0.0;
    }
    Real computed2QdI4sdI8fn (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        return 0.0;
    }
    Real computed2QdI4sdI8sn (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        return 0.0;
    }


    Real computedQdI8fs (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real I8fs = Elasticity::I8 (F, f0, s0);
        return M_bfs * I8fs;
    }
    Real computed2QdI8fs2 (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        return M_bfs;
    }
    Real computed2QdI8fsdI8fn (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        return 0.0;
    }
    Real computed2QdI8fsdI8sn (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        return 0.0;
    }

    Real computedQdI8fn (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        VectorSmall<3> n0 = Elasticity::crossProduct (f0, s0);
        Real I8fn = Elasticity::I8 (F, f0, n0);
        return M_bfn * I8fn;
    }
    Real computed2QdI8fn2 (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        return M_bfn;
    }
    Real computed2QdI8fndI8sn (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        return 0.0;
    }

    Real computedQdI8sn (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        VectorSmall<3> n0 = Elasticity::crossProduct (f0, s0);
        Real I8sn = Elasticity::I8 (F, s0, n0);
        return M_bsn * I8sn;
    }
    Real computed2QdI8sn2 (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        return M_bsn;
    }






protected:
    Real M_a;
    Real M_bff;
    Real M_bss;
    Real M_bnn;
    Real M_bfs;
    Real M_bfn;
    Real M_bsn;
};

template <class Mesh>
class W1 : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    W1 (Real a = 2000,        // 0.23 kPa
        Real bff = 35.7,
        Real bss = 18.9,
        Real bnn = 13.0,
        Real bfs = 10.5,
        Real bfn = 11.4,
        Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    W1 (const W1& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~W1() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI1bar = super::computedQdI1 (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * dQdI1bar;
    }
};

template <class Mesh>
class dW1dI1 : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW1dI1 (Real a = 2000,        // 0.23 kPa
            Real bff = 35.7,
            Real bss = 18.9,
            Real bnn = 13.0,
            Real bfs = 10.5,
            Real bfn = 11.4,
            Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW1dI1 (const dW1dI1& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW1dI1() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI1bar = super::computedQdI1 (F, f0, s0);
        Real d2QdI1bar2 = super::computed2QdI12 (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI1bar * dQdI1bar + d2QdI1bar2 );
    }
};

template <class Mesh>
class dW1dI4f : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW1dI4f (Real a = 2000,        // 0.23 kPa
             Real bff = 35.7,
             Real bss = 18.9,
             Real bnn = 13.0,
             Real bfs = 10.5,
             Real bfn = 11.4,
             Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW1dI4f (const dW1dI4f& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW1dI4f() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI1bar = super::computedQdI1 (F, f0, s0);
        Real dQdI4f = super::computedQdI4f (F, f0, s0);
        Real d2QdI1bardI4f = super::computed2QdI1dI4f (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI1bar * dQdI4f + d2QdI1bardI4f );
    }
};

template <class Mesh>
class dW1dI4s : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW1dI4s (Real a = 2000,        // 0.23 kPa
             Real bff = 35.7,
             Real bss = 18.9,
             Real bnn = 13.0,
             Real bfs = 10.5,
             Real bfn = 11.4,
             Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW1dI4s (const dW1dI4s& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW1dI4s() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI1bar = super::computedQdI1 (F, f0, s0);
        Real dQdI4s = super::computedQdI4s (F, f0, s0);
        Real d2QdI1bardI4s = super::computed2QdI1dI4s (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI1bar * dQdI4s + d2QdI1bardI4s );
    }
};

template <class Mesh>
class dW1dI8fs : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW1dI8fs (Real a = 2000,        // 0.23 kPa
              Real bff = 35.7,
              Real bss = 18.9,
              Real bnn = 13.0,
              Real bfs = 10.5,
              Real bfn = 11.4,
              Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW1dI8fs (const dW1dI8fs& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW1dI8fs() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI1 = super::computedQdI1 (F, f0, s0);
        Real dQdI8fs = super::computedQdI8fs (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI1 * dQdI8fs ); // + d2QdI1dI8fs ( == 0 ) );
    }
};

template <class Mesh>
class dW1dI8fn : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW1dI8fn (Real a = 2000,        // 0.23 kPa
              Real bff = 35.7,
              Real bss = 18.9,
              Real bnn = 13.0,
              Real bfs = 10.5,
              Real bfn = 11.4,
              Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW1dI8fn (const dW1dI8fn& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW1dI8fn() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI1 = super::computedQdI1 (F, f0, s0);
        Real dQdI8fn = super::computedQdI8fn (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI1 * dQdI8fn ); // + d2QdI1dI8fs ( == 0 ) );
    }
};

template <class Mesh>
class dW1dI8sn : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW1dI8sn (Real a = 2000,        // 0.23 kPa
              Real bff = 35.7,
              Real bss = 18.9,
              Real bnn = 13.0,
              Real bfs = 10.5,
              Real bfn = 11.4,
              Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW1dI8sn (const dW1dI8sn& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW1dI8sn() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI1 = super::computedQdI1 (F, f0, s0);
        Real dQdI8sn = super::computedQdI8sn (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI1 * dQdI8sn ); // + d2QdI1dI8fs ( == 0 ) );
    }
};


template <class Mesh>
class W4f : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    W4f (Real a = 2000,        // 0.23 kPa
         Real bff = 35.7,
         Real bss = 18.9,
         Real bnn = 13.0,
         Real bfs = 10.5,
         Real bfn = 11.4,
         Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    W4f (const W4f& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~W4f() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI4f = super::computedQdI4f (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * dQdI4f;

    }
};

template <class Mesh>
class dW4fdI1 : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW4fdI1 (Real a = 2000,        // 0.23 kPa
             Real bff = 35.7,
             Real bss = 18.9,
             Real bnn = 13.0,
             Real bfs = 10.5,
             Real bfn = 11.4,
             Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW4fdI1 (const dW4fdI1& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW4fdI1() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI4f = super::computedQdI4f (F, f0, s0);
        Real dQdI1 = super::computedQdI1 (F, f0, s0);
        Real d2QdI4fdI1 = super::computed2QdI1dI4f (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI4f * dQdI1 + d2QdI4fdI1 );
    }
};

template <class Mesh>
class dW4fdI4f : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW4fdI4f (Real a = 2000,        // 0.23 kPa
              Real bff = 35.7,
              Real bss = 18.9,
              Real bnn = 13.0,
              Real bfs = 10.5,
              Real bfn = 11.4,
              Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW4fdI4f (const dW4fdI4f& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW4fdI4f() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI4f = super::computedQdI4f (F, f0, s0);
        Real d2QdI4f2 = super::computed2QdI4f2 (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI4f * dQdI4f + d2QdI4f2 );
    }
};

template <class Mesh>
class dW4fdI4s : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW4fdI4s (Real a = 2000,        // 0.23 kPa
              Real bff = 35.7,
              Real bss = 18.9,
              Real bnn = 13.0,
              Real bfs = 10.5,
              Real bfn = 11.4,
              Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW4fdI4s (const dW4fdI4s& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW4fdI4s() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI4f = super::computedQdI4f (F, f0, s0);
        Real dQdI4s = super::computedQdI4s (F, f0, s0);
        Real d2QdI4fdI4s = super::computed2QdI4fdI4s (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI4f * dQdI4s + d2QdI4fdI4s );
    }
};

template <class Mesh>
class dW4fdI8fs : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW4fdI8fs (Real a = 2000,        // 0.23 kPa
               Real bff = 35.7,
               Real bss = 18.9,
               Real bnn = 13.0,
               Real bfs = 10.5,
               Real bfn = 11.4,
               Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW4fdI8fs (const dW4fdI8fs& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW4fdI8fs() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI4f = super::computedQdI4f (F, f0, s0);
        Real dQdI8fs = super::computedQdI8fs (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI4f * dQdI8fs ); // + d2QdI4fdI8fs ( == 0 ) );
    }
};

template <class Mesh>
class dW4fdI8fn : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW4fdI8fn (Real a = 2000,        // 0.23 kPa
               Real bff = 35.7,
               Real bss = 18.9,
               Real bnn = 13.0,
               Real bfs = 10.5,
               Real bfn = 11.4,
               Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW4fdI8fn (const dW4fdI8fn& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW4fdI8fn() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI4f = super::computedQdI4f (F, f0, s0);
        Real dQdI8fn = super::computedQdI8fn (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI4f * dQdI8fn ); // + d2QdI4fdI8fn ( == 0 ) );
    }
};

template <class Mesh>
class dW4fdI8sn : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW4fdI8sn (Real a = 2000,        // 0.23 kPa
               Real bff = 35.7,
               Real bss = 18.9,
               Real bnn = 13.0,
               Real bfs = 10.5,
               Real bfn = 11.4,
               Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW4fdI8sn (const dW4fdI8sn& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW4fdI8sn() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI4f = super::computedQdI4f (F, f0, s0);
        Real dQdI8sn = super::computedQdI8sn (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI4f * dQdI8sn ); // + d2QdI4fdI8sn ( == 0 ) );
    }
};


template <class Mesh>
class W4s : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    W4s (Real a = 2000,        // 0.23 kPa
         Real bff = 35.7,
         Real bss = 18.9,
         Real bnn = 13.0,
         Real bfs = 10.5,
         Real bfn = 11.4,
         Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    W4s (const W4s& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~W4s() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {

        Real Q = super::computeQ (F, f0, s0);
        Real dQdI4s =  super::computedQdI4s (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * dQdI4s;
    }
};

template <class Mesh>
class dW4sdI1 : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW4sdI1 (Real a = 2000,        // 0.23 kPa
             Real bff = 35.7,
             Real bss = 18.9,
             Real bnn = 13.0,
             Real bfs = 10.5,
             Real bfn = 11.4,
             Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW4sdI1 (const dW4sdI1& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW4sdI1() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI4s =  super::computedQdI4s (F, f0, s0);
        Real dQdI1  =  super::computedQdI1 (F, f0, s0);
        Real d2QdI4sdI1 = super::computed2QdI1dI4s (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI4s * dQdI1 + d2QdI4sdI1 );
    }
};

template <class Mesh>
class dW4sdI4s : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW4sdI4s (Real a = 2000,        // 0.23 kPa
              Real bff = 35.7,
              Real bss = 18.9,
              Real bnn = 13.0,
              Real bfs = 10.5,
              Real bfn = 11.4,
              Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW4sdI4s (const dW4sdI4s& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW4sdI4s() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI4s =  super::computedQdI4s (F, f0, s0);
        Real d2QdI4s2 = super::computed2QdI4s2 (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI4s * dQdI4s + d2QdI4s2 );
    }
};


template <class Mesh>
class dW4sdI4f : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW4sdI4f (Real a = 2000,        // 0.23 kPa
              Real bff = 35.7,
              Real bss = 18.9,
              Real bnn = 13.0,
              Real bfs = 10.5,
              Real bfn = 11.4,
              Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW4sdI4f (const dW4sdI4f& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW4sdI4f() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI4s =  super::computedQdI4s (F, f0, s0);
        Real dQdI4f =  super::computedQdI4f (F, f0, s0);
        Real d2QdI4sdI4f = super::computed2QdI4fdI4s (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI4f * dQdI4s + d2QdI4sdI4f );
    }
};

template <class Mesh>
class dW4sdI8fs : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW4sdI8fs (Real a = 2000,        // 0.23 kPa
               Real bff = 35.7,
               Real bss = 18.9,
               Real bnn = 13.0,
               Real bfs = 10.5,
               Real bfn = 11.4,
               Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW4sdI8fs (const dW4sdI8fs& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW4sdI8fs() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI4s = super::computedQdI4s (F, f0, s0);
        Real dQdI8fs = super::computedQdI8fs (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI4s * dQdI8fs ); // + d2QdI4sdI8fs ( == 0 ) );
    }
};

template <class Mesh>
class dW4sdI8fn : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW4sdI8fn (Real a = 2000,        // 0.23 kPa
               Real bff = 35.7,
               Real bss = 18.9,
               Real bnn = 13.0,
               Real bfs = 10.5,
               Real bfn = 11.4,
               Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW4sdI8fn (const dW4sdI8fn& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW4sdI8fn() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI4s = super::computedQdI4s (F, f0, s0);
        Real dQdI8fn = super::computedQdI8fn (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI4s * dQdI8fn ); // + d2QdI4sdI8fs ( == 0 ) );
    }
};

template <class Mesh>
class dW4sdI8sn : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW4sdI8sn (Real a = 2000,        // 0.23 kPa
               Real bff = 35.7,
               Real bss = 18.9,
               Real bnn = 13.0,
               Real bfs = 10.5,
               Real bfn = 11.4,
               Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW4sdI8sn (const dW4sdI8sn& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW4sdI8sn() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI4s = super::computedQdI4s (F, f0, s0);
        Real dQdI8sn = super::computedQdI8sn (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI4s * dQdI8sn ); // + d2QdI4sdI8sn ( == 0 ) );
    }
};



template <class Mesh>
class W8fs : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    W8fs (Real a = 2000,        // 0.23 kPa
          Real bff = 35.7,
          Real bss = 18.9,
          Real bnn = 13.0,
          Real bfs = 10.5,
          Real bfn = 11.4,
          Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    W8fs (const W8fs& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~W8fs() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI8fs =  super::computedQdI8fs (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * dQdI8fs;
    }
};

template <class Mesh>
class dW8fsdI8fs : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW8fsdI8fs (Real a = 2000,         // 0.23 kPa
                Real bff = 35.7,
                Real bss = 18.9,
                Real bnn = 13.0,
                Real bfs = 10.5,
                Real bfn = 11.4,
                Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW8fsdI8fs (const dW8fsdI8fs& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW8fsdI8fs () {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI8fs =  super::computedQdI8fs (F, f0, s0);
        Real d2QdI8fs2 = super::computed2QdI8fs2 (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI8fs * dQdI8fs + d2QdI8fs2 );
    }
};

template <class Mesh>
class dW8fsdI8fn : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW8fsdI8fn (Real a = 2000,        // 0.23 kPa
                Real bff = 35.7,
                Real bss = 18.9,
                Real bnn = 13.0,
                Real bfs = 10.5,
                Real bfn = 11.4,
                Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW8fsdI8fn (const dW8fsdI8fn& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW8fsdI8fn() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI8fs = super::computedQdI8fs (F, f0, s0);
        Real dQdI8fn = super::computedQdI8fn (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI8fs * dQdI8fn ); // + d2QdI4sdI8fs ( == 0 ) );
    }
};

template <class Mesh>
class dW8fsdI8sn : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW8fsdI8sn (Real a = 2000,        // 0.23 kPa
                Real bff = 35.7,
                Real bss = 18.9,
                Real bnn = 13.0,
                Real bfs = 10.5,
                Real bfn = 11.4,
                Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW8fsdI8sn (const dW8fsdI8sn& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW8fsdI8sn() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI8fs = super::computedQdI8fs (F, f0, s0);
        Real dQdI8sn = super::computedQdI8sn (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI8fs * dQdI8sn ); // + d2QdI4sdI8fs ( == 0 ) );
    }
};


template <class Mesh>
class W8fn : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    W8fn (Real a = 2000,        // 0.23 kPa
          Real bff = 35.7,
          Real bss = 18.9,
          Real bnn = 13.0,
          Real bfs = 10.5,
          Real bfn = 11.4,
          Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    W8fn (const W8fn& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~W8fn() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI8fn =  super::computedQdI8fn (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * dQdI8fn;
    }
};

template <class Mesh>
class dW8fndI8fn : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW8fndI8fn (Real a = 2000,         // 0.23 kPa
                Real bff = 35.7,
                Real bss = 18.9,
                Real bnn = 13.0,
                Real bfs = 10.5,
                Real bfn = 11.4,
                Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW8fndI8fn (const dW8fndI8fn& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW8fndI8fn () {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI8fn =  super::computedQdI8fn (F, f0, s0);
        Real d2QdI8fn2 = super::computed2QdI8fn2 (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI8fn * dQdI8fn + d2QdI8fn2 );
    }
};

template <class Mesh>
class dW8fndI8sn : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW8fndI8sn (Real a = 2000,        // 0.23 kPa
                Real bff = 35.7,
                Real bss = 18.9,
                Real bnn = 13.0,
                Real bfs = 10.5,
                Real bfn = 11.4,
                Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW8fndI8sn (const dW8fndI8sn& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW8fndI8sn() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI8fn = super::computedQdI8fn (F, f0, s0);
        Real dQdI8sn = super::computedQdI8sn (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI8fn * dQdI8sn ); // + d2QdI8fndI8sn ( == 0 ) );
    }
};


template <class Mesh>
class W8sn : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    W8sn (Real a = 2000,        // 0.23 kPa
          Real bff = 35.7,
          Real bss = 18.9,
          Real bnn = 13.0,
          Real bfs = 10.5,
          Real bfn = 11.4,
          Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    W8sn (const W8sn& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~W8sn() {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI8sn =  super::computedQdI8sn (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * dQdI8sn;
    }
};

template <class Mesh>
class dW8sndI8sn : public virtual W<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef W<Mesh> super;

    dW8sndI8sn (Real a = 2000,         // 0.23 kPa
                Real bff = 35.7,
                Real bss = 18.9,
                Real bnn = 13.0,
                Real bfs = 10.5,
                Real bfn = 11.4,
                Real bsn = 8.35) : super (a, bff, bss, bnn, bfs, bfn, bsn) {}

    dW8sndI8sn (const dW8sndI8sn& OrthotropicFung) : super (OrthotropicFung) {}
    virtual ~dW8sndI8sn () {}

    return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f0, const VectorSmall<3>& s0)
    {
        Real Q = super::computeQ (F, f0, s0);
        Real dQdI8sn =  super::computedQdI8sn (F, f0, s0);
        Real d2QdI8sn2 = super::computed2QdI8sn2 (F, f0, s0);
        return  0.5 * super::M_a * std::exp (Q) * ( dQdI8sn * dQdI8sn + d2QdI8sn2 );
    }
};




}

template <class Mesh>
class OrthotropicFung : public virtual EMMaterialFunctions<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef  boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > ETFESpacePtr_Type;

    //Values taken from Schmid et al. 2008, Myocardial Parameter Estimation Table 2 Exp 5
    OrthotropicFung (Real a = 2000,        // 0.2 kPa
                     Real bff = 35.7,
                     Real bss = 18.9,
                     Real bnn = 13.0,
                     Real bfs = 10.5,
                     Real bfn = 11.4,
                     Real bsn = 8.35) :
        M_a (a),
        M_bff (bff),
        M_bss (bss),
        M_bnn (bnn),
        M_bfs (bfs),
        M_bfn (bfn),
        M_bsn (bsn)
    {
    	std::cout << "C  = " << M_a   << "\n";
    	std::cout << "bff= " << M_bff << "\n";
    	std::cout << "bss= " << M_bss << "\n";
    	std::cout << "bnn= " << M_bnn << "\n";
    	std::cout << "bfs= " << M_bfs << "\n";
    	std::cout << "bfn= " << M_bfn << "\n";
    	std::cout << "bsn= " << M_bsn << "\n";
    }

    OrthotropicFung (const OrthotropicFung& OrthotropicFung)
    {
        M_a = OrthotropicFung.M_a;
        M_bff = OrthotropicFung.M_bff;
        M_bss = OrthotropicFung.M_bss;
        M_bnn = OrthotropicFung.M_bnn;
        M_bfs = OrthotropicFung.M_bfs;
        M_bfn = OrthotropicFung.M_bfn;
        M_bsn = OrthotropicFung.M_bsn;
    }
    virtual ~OrthotropicFung() {}



    void computeJacobian ( const vector_Type& disp,
                           ETFESpacePtr_Type  dispETFESpace,
                           const vector_Type& fibers,
                           const vector_Type& sheets,
                           matrixPtr_Type           jacobianPtr);

    void computeResidual ( const vector_Type& disp,
                           ETFESpacePtr_Type dispETFESpace,
                           const vector_Type& fibers,
                           const vector_Type& sheets,
                           vectorPtr_Type           residualVectorPtr);

    typedef EMData          data_Type;
    void setParameters (data_Type& data)
    {
        M_a   = data.parameter("C");
        M_bff = data.parameter("bff");
        M_bss = data.parameter("bss");
        M_bnn = data.parameter("bnn");
        M_bfs = data.parameter("bfs");
        M_bfn = data.parameter("bfn");
        M_bsn = data.parameter("bsn");
    }

    void showMe()
    {
    	std::cout << "C  = " << M_a   << "\n";
    	std::cout << "bff= " << M_bff << "\n";
    	std::cout << "bss= " << M_bss << "\n";
    	std::cout << "bnn= " << M_bnn << "\n";
    	std::cout << "bfs= " << M_bfs << "\n";
    	std::cout << "bfn= " << M_bfn << "\n";
    	std::cout << "bsn= " << M_bsn << "\n";
    }

protected:
    Real M_a;
    Real M_bff;
    Real M_bss;
    Real M_bnn;
    Real M_bfs;
    Real M_bfn;
    Real M_bsn;

};

template <class Mesh>
void OrthotropicFung<Mesh>::computeResidual (  const vector_Type& disp,
                                               ETFESpacePtr_Type dispETFESpace,
                                               const vector_Type& fibers,
                                               const vector_Type& sheets,
                                               vectorPtr_Type           residualVectorPtr)
{
    using namespace ExpressionAssembly;

    boost::shared_ptr<FungFunctors::W1<Mesh> > W1_bar ( new FungFunctors::W1<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
    boost::shared_ptr<FungFunctors::W4f<Mesh> > W4_f ( new FungFunctors::W4f<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
    boost::shared_ptr<FungFunctors::W4s<Mesh> > W4_s ( new FungFunctors::W4s<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
    boost::shared_ptr<FungFunctors::W8fs<Mesh> > W8_fs ( new FungFunctors::W8fs<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
    boost::shared_ptr<FungFunctors::W8fn<Mesh> > W8_fn ( new FungFunctors::W8fn<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
    boost::shared_ptr<FungFunctors::W8sn<Mesh> > W8_sn ( new FungFunctors::W8sn<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );

    boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
    boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
    boost::shared_ptr<CrossProduct> cross ( new CrossProduct );

    auto F = _F (dispETFESpace, disp, 0);
    auto f_0 = _v0 (dispETFESpace, fibers);
    auto f0 = eval (normalize0, f_0);
    auto s_0 = _v0 (dispETFESpace, sheets);
    auto s_00 = s_0 - dot (f0, s_0) * s_0;
    auto s0 = eval (normalize1, s_00);

    auto n0 = eval (cross, f0, s0);

    auto W1bar = eval(W1_bar, F, f0, s0);
    auto dI1bar = _dI1bar(dispETFESpace, disp, 0);
//    auto dI1bar = _dI1(dispETFESpace, disp, 0);
	auto P1bar = W1bar * dI1bar;

    auto W4F = eval (W4_f, F, f0, s0);
    auto dI4F = value (2.0) * F * outerProduct (f0, f0);
    auto P4F = W4F * dI4F;

    auto W4S = eval (W4_s, F, f0, s0);
    auto dI4S = value (2.0) * F * outerProduct (s0, s0);
    auto P4S = W4S * dI4S;

    auto W8FS = eval (W8_fs, F, f0, s0);
    auto dI8FS = _dI8 (dispETFESpace, disp, 0, f0, s0);
    auto P8FS = W8FS * dI8FS;

    auto W8FN = eval (W8_fn, F, f0, s0);
    auto dI8FN = _dI8 (dispETFESpace, disp, 0, f0, n0);
    auto P8FN = W8FN * dI8FN;

    auto W8SN = eval (W8_sn, F, f0, s0);
    auto dI8SN = _dI8 (dispETFESpace, disp, 0, s0, n0);
    auto P8SN = W8SN * dI8SN;

    auto P = P4F + P4S + P1bar + P8FS + P8FN + P8SN;

    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dot ( P , grad (phi_i) )
              ) >> residualVectorPtr;

}

template <class Mesh>
void OrthotropicFung<Mesh>::computeJacobian (  const vector_Type& disp,
                                               ETFESpacePtr_Type dispETFESpace,
                                               const vector_Type& fibers,
                                               const vector_Type& sheets,
                                               matrixPtr_Type           jacobianPtr)
{
	using namespace ExpressionAssembly;

	boost::shared_ptr<FungFunctors::W1<Mesh> > W1_bar( new FungFunctors::W1<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW1dI1<Mesh> > dW1_bar( new FungFunctors::dW1dI1<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW1dI4f<Mesh> > dW1_dI4f( new FungFunctors::dW1dI4f<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW1dI4s<Mesh> > dW1_dI4s( new FungFunctors::dW1dI4s<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW1dI8fs<Mesh> > dW1_dI8fs( new FungFunctors::dW1dI8fs<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW1dI8fn<Mesh> > dW1_dI8fn( new FungFunctors::dW1dI8fn<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW1dI8sn<Mesh> > dW1_dI8sn( new FungFunctors::dW1dI8sn<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );

	boost::shared_ptr<FungFunctors::W4f<Mesh> > W4_f( new FungFunctors::W4f<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW4fdI1<Mesh> > dW4f_dI1( new FungFunctors::dW4fdI1<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW4fdI4f<Mesh> > dW4f_dI4f( new FungFunctors::dW4fdI4f<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW4fdI4s<Mesh> > dW4f_dI4s( new FungFunctors::dW4fdI4s<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW4fdI8fs<Mesh> > dW4f_dI8fs( new FungFunctors::dW4fdI8fs<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW4fdI8fn<Mesh> > dW4f_dI8fn( new FungFunctors::dW4fdI8fn<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW4fdI8sn<Mesh> > dW4f_dI8sn( new FungFunctors::dW4fdI8sn<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );

	boost::shared_ptr<FungFunctors::W4s<Mesh> > W4_s( new FungFunctors::W4s<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW4sdI1<Mesh> > dW4s_dI1( new FungFunctors::dW4sdI1<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW4sdI4f<Mesh> > dW4s_dI4f( new FungFunctors::dW4sdI4f<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW4sdI4s<Mesh> > dW4s_dI4s( new FungFunctors::dW4sdI4s<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW4sdI8fs<Mesh> > dW4s_dI8fs( new FungFunctors::dW4sdI8fs<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW4sdI8fn<Mesh> > dW4s_dI8fn( new FungFunctors::dW4sdI8fn<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW4sdI8sn<Mesh> > dW4s_dI8sn( new FungFunctors::dW4sdI8sn<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );


	boost::shared_ptr<FungFunctors::W8fs<Mesh> > W8_fs( new FungFunctors::W8fs<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW8fsdI8fs<Mesh> > dW8fs_dI8fs( new FungFunctors::dW8fsdI8fs<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW8fsdI8fn<Mesh> > dW8fs_dI8fn( new FungFunctors::dW8fsdI8fn<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW8fsdI8sn<Mesh> > dW8fs_dI8sn( new FungFunctors::dW8fsdI8sn<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );

	boost::shared_ptr<FungFunctors::W8fn<Mesh> > W8_fn( new FungFunctors::W8fn<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW8fndI8fn<Mesh> > dW8fn_dI8fn( new FungFunctors::dW8fndI8fn<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW8fndI8sn<Mesh> > dW8fn_dI8sn( new FungFunctors::dW8fndI8sn<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );

	boost::shared_ptr<FungFunctors::W8sn<Mesh> > W8_sn( new FungFunctors::W8sn<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );
	boost::shared_ptr<FungFunctors::dW8sndI8sn<Mesh> > dW8sn_dI8sn( new FungFunctors::dW8sndI8sn<Mesh>(M_a, M_bff, M_bss, M_bnn, M_bfs, M_bfn, M_bsn) );

	boost::shared_ptr<orthonormalizeFibers> normalize0(new orthonormalizeFibers);
    boost::shared_ptr<orthonormalizeFibers> normalize1(new orthonormalizeFibers(1) );
	boost::shared_ptr<CrossProduct> cross( new CrossProduct );

	auto F = _F(dispETFESpace, disp, 0);
	auto f_0 = _v0(dispETFESpace, fibers);
	auto f0 = eval(normalize0, f_0);
	auto s_0 = _v0(dispETFESpace, sheets);
    auto s_00 = s_0 - dot(f0, s_0) * s_0;
    auto s0 = eval(normalize1, s_00);

    auto n0 = eval(cross, f0, s0);


    auto W1bar = eval(W1_bar, F, f0, s0);
    auto dW1bardI1= eval(dW1_bar,F, f0, s0);
    auto dW1bardI4F= eval(dW1_dI4f,F, f0, s0);
    auto dW1bardI4S= eval(dW1_dI4s,F, f0, s0);
    auto dW1bardI8FS=eval(dW1_dI8fs,F, f0, s0);
    auto dW1bardI8FN=eval(dW1_dI8fn,F, f0, s0);
    auto dW1bardI8SN=eval(dW1_dI8sn,F, f0, s0);
    auto dI1bar = _dI1bar(dispETFESpace, disp, 0);
    auto dI1bardF = _dI1bardF(dispETFESpace, disp, 0);
	auto d2I1bardF= _d2I1bardF(dispETFESpace, disp, 0);
//    auto dI1bar = _dI1(dispETFESpace, disp, 0);
//    auto dI1bardF = _dI1dF(dispETFESpace, disp, 0);
//	auto d2I1bardF= _d2I1dF;

    auto W4F = eval(W4_f, F, f0, s0);
    auto dW4FdI4F = eval(dW4f_dI4f, F, f0, s0);
    auto dW4FdI1= eval(dW4f_dI1,F, f0, s0);
    auto dW4FdI4S= eval(dW4f_dI4s,F, f0, s0);
    auto dW4FdI8FS= eval(dW4f_dI8fs,F, f0, s0);
    auto dW4FdI8FN= eval(dW4f_dI8fn,F, f0, s0);
    auto dW4FdI8SN= eval(dW4f_dI8sn,F, f0, s0);
    auto dI4F = value(2.0) * F * outerProduct(f0, f0);
    auto dI4FdF = dot( dI4F, _dF );
    auto d2I4FdF = value(2.0) * _dF * outerProduct(f0, f0);

    auto W4S = eval(W4_s, F, f0, s0);
    auto dW4SdI4S = eval(dW4s_dI4s, F, f0, s0);
    auto dW4SdI1= eval(dW4s_dI1,F, f0, s0);
    auto dW4SdI4F= eval(dW4s_dI4f,F, f0, s0);
    auto dW4SdI8FS= eval(dW4s_dI8fs,F, f0, s0);
    auto dW4SdI8FN= eval(dW4s_dI8fn,F, f0, s0);
    auto dW4SdI8SN= eval(dW4s_dI8sn,F, f0, s0);
    auto dI4S = value(2.0) * F * outerProduct(s0, s0);
    auto dI4SdF = dot( dI4S, _dF );
    auto d2I4SdF = value(2.0) * _dF * outerProduct(s0, s0);

	auto W8FS     = eval(W8_fs, F, f0, s0);
    auto dW8FSdI1 = eval(dW1_dI8fs,F, f0, s0);
    auto dW8FSdI4F= eval(dW4f_dI8fs,F, f0, s0);
    auto dW8FSdI4S= eval(dW4s_dI8fs,F, f0, s0);
    auto dW8FSdI8FS=eval(dW8fs_dI8fs,F, f0, s0);
    auto dW8FSdI8FN=eval(dW8fs_dI8fn,F, f0, s0);
    auto dW8FSdI8SN=eval(dW8fs_dI8sn,F, f0, s0);
	auto dI8FS= _dI8(dispETFESpace, disp, 0, f0, s0);
	auto dI8FSdF= _dI8dF(dispETFESpace, disp, 0, f0, s0);
	auto d2I8FSdF= _d2I8dF(dispETFESpace, f0, s0);

	auto W8FN     = eval(W8_fn, F, f0, s0);
    auto dW8FNdI1 = eval(dW1_dI8fn,F, f0, s0);
    auto dW8FNdI4F= eval(dW4f_dI8fn,F, f0, s0);
    auto dW8FNdI4S= eval(dW4s_dI8fn,F, f0, s0);
    auto dW8FNdI8FS=eval(dW8fs_dI8fn,F, f0, s0);
    auto dW8FNdI8FN=eval(dW8fn_dI8fn,F, f0, s0);
    auto dW8FNdI8SN=eval(dW8fn_dI8sn,F, f0, s0);
	auto dI8FN= _dI8(dispETFESpace, disp, 0, f0, n0);
	auto dI8FNdF= _dI8dF(dispETFESpace, disp, 0, f0, n0);
	auto d2I8FNdF= _d2I8dF(dispETFESpace, f0, n0);

	auto W8SN     = eval(W8_sn, F, f0, s0);
    auto dW8SNdI1 = eval(dW1_dI8sn,F, f0, s0);
    auto dW8SNdI4F= eval(dW4f_dI8sn,F, f0, s0);
    auto dW8SNdI4S= eval(dW4s_dI8sn,F, f0, s0);
    auto dW8SNdI8FS=eval(dW8fs_dI8sn,F, f0, s0);
    auto dW8SNdI8FN=eval(dW8fn_dI8sn,F, f0, s0);
    auto dW8SNdI8SN=eval(dW8sn_dI8sn,F, f0, s0);
	auto dI8SN= _dI8(dispETFESpace, disp, 0, s0, n0);
	auto dI8SNdF= _dI8dF(dispETFESpace, disp, 0, s0, n0);
	auto d2I8SNdF= _d2I8dF(dispETFESpace, s0, n0);

	auto dP1bar = W1bar * d2I1bardF
			    + ( dW1bardI1 * dI1bardF
			      + dW1bardI4F * dI4FdF
			      + dW1bardI4S * dI4SdF
			      + dW1bardI8FS * dI8FSdF
			      + dW1bardI8FN * dI8FNdF
			      + dW1bardI8SN * dI8SNdF  ) * dI1bar;

    auto dP4F = W4F * d2I4FdF
                + ( dW4FdI1 * dI1bardF
                    + dW4FdI4F * dI4FdF
                    + dW4FdI4S * dI4SdF
                    + dW4FdI8FS * dI8FSdF
                    + dW4FdI8FN * dI8FNdF
                    + dW4FdI8SN * dI8SNdF  ) * dI4F;

    auto dP4S = W4S * d2I4SdF
                + ( dW4SdI1 * dI1bardF
                    + dW4SdI4F * dI4FdF
                    + dW4SdI4S * dI4SdF
                    + dW4SdI8FS * dI8FSdF
                    + dW4SdI8FN * dI8FNdF
                    + dW4SdI8SN * dI8SNdF  ) * dI4S;

    auto dP8FS = W8FS * d2I8FSdF
                 + ( dW8FSdI1 * dI1bardF
                     + dW8FSdI4F * dI4FdF
                     + dW8FSdI4S * dI4SdF
                     + dW8FSdI8FS * dI8FSdF
                     + dW8FSdI8FN * dI8FNdF
                     + dW8FSdI8SN * dI8SNdF ) * dI8FS;

    auto dP8FN = W8FN * d2I8FNdF
                 + ( dW8FNdI1 * dI1bardF
                     + dW8FNdI4F * dI4FdF
                     + dW8FNdI4S * dI4SdF
                     + dW8FNdI8FS * dI8FSdF
                     + dW8FNdI8FN * dI8FNdF
                     + dW8FNdI8SN * dI8SNdF ) * dI8FN;

    auto dP8SN = W8SN * d2I8SNdF
                 + ( dW8SNdI1 * dI1bardF
                     + dW8SNdI4F * dI4FdF
                     + dW8SNdI4S * dI4SdF
                     + dW8SNdI8FS * dI8FSdF
                     + dW8SNdI8FN * dI8FNdF
                     + dW8SNdI8SN * dI8SNdF ) * dI8SN;

    //    auto dP = dP4F + dP4F;
    auto dP = dP4F + dP4S + dP1bar + dP8FS + dP8FN + dP8SN;
    //    auto dP = W4S * d2I4SdF + dW4SdI1 * dI1bardF * dI4S + dW4SdI4S * dI4SdF * dI4S;
    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dispETFESpace,
                dot ( dP , grad (phi_i) )
              ) >> jacobianPtr;
}



} //EMMaterialFunctions

} //LifeV
#endif /* EMMATERIALFUNCTIONS_HPP_ */
