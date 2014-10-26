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
    @file
    @brief File where the structures for the element-wise multiplication between expressions are defined.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    @date 08-2012
 */

#ifndef EXPRESSION_COFACTOR_HPP
#define EXPRESSION_COFACTOR_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>
#include <lifev/eta/expression/ExpressionMatrix.hpp>

namespace LifeV
{

namespace ExpressionAssembly
{

//! class ExpressionEmult  Class for representing the transpose operation of an expression
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class represents the transpose operation in the expression tree.

  <b> Template parameters </b>

  <i>ExpressionType</i>: The expression to be transposed.

  <b> Template requirements </b>

  <i>ExpressionType</i>: Copiable, static display method


*/
template <typename ExpressionType>
class ExpressionCofactor : public ExpressionBase< ExpressionCofactor<ExpressionType> >
{
public:

    //! @name Public Types
    //@{

    // No real need, just for ease of coding
    typedef ExpressionBase< ExpressionCofactor <ExpressionType> > base_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Full constructor
    ExpressionCofactor (const ExpressionType& expr)
        : base_Type(), M_expr (expr) {}

    //! Copy constructor
    ExpressionCofactor (const ExpressionCofactor<ExpressionType>& expression)
        : base_Type(), M_expr (expression.M_expr) {}

    //! Destructor
    ~ExpressionCofactor() {}

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        out << " transpose ";
        ExpressionType::display (out);
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the expression that we transpose
    const ExpressionType& exprEx() const
    {
        return M_expr;
    }


    //@}

private:

    //! @name Private Methods
    //@{

    //! No default constructor
    ExpressionCofactor();

    //@}

    // Expression that we transpose
    ExpressionType M_expr;

};


// transpose  The generic function for the transpose of an expression.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  Function used in construction of the expression tree. To avoid shadowing
  other functions, it uses the ExpressionBase type to distinguish expressions
  from other types.

  <b> Template parameters </b>

  <i>ExpressionType</i>: The expression to be transposed.

  <b> Template requirements </b>

  <i>ExpressionType</i>: Same as in LifeV::ExpressionCofactor

*/
template< typename ExpressionType >
ExpressionCofactor<ExpressionType>
cofactor (const ExpressionBase<ExpressionType>& expr)
{
    return ExpressionCofactor<ExpressionType> (expr.cast() );
}


// Specialization for the matricial constants
template< UInt Dim1, UInt Dim2>
ExpressionCofactor<ExpressionMatrix<Dim1, Dim2> >
cofactor (const MatrixSmall<Dim1, Dim2>& m)
{
    return ExpressionCofactor<ExpressionMatrix<Dim1, Dim2> > (ExpressionMatrix<Dim1, Dim2> (m) );
}


} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
