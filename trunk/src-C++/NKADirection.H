//   NONLINEAR_KRYLOV_ACCELERATOR - C++/Trilinos implementation

//   Neil N. Carlson <neil.n.carlson@gmail.com>
//   Markus Berndt <markus.berndt@gmail.com>

//   This code implements the nonlinear Krylov accelerator introduced in [1]
//   for inexact Newton's (IN) method, where the correction equation of
//   Newton's method is only approximately solved because the Jacobian matrix
//   is approximated and/or the linear system is not solved exactly.  Placed
//   in the iteration loop, this black-box accelerator listens to the sequence
//   of inexact corrections and replaces them with accelerated corrections;
//   the resulting method is a type of accelerated inexact Newton (AIN) method.
//   Note that an IN iteration is merely a standard fixed point iteration for
//   a preconditioned system, and so this accelerator is more generally
//   applicable to fixed point iterations.

//   This code is a straightforward translation of the original Fortran 95
//   implementation into C++/Trilinos.

//   [1] N.N.Carlson and K.Miller, "Design and application of a gradient-
//       weighted moving finite element code I: in one dimension", SIAM J.
//       Sci. Comput;, 19 (1998), pp. 728-765.  See section 9.

//   ************************************************************************

//   Copyright (c) 2010  Neil N. Carlson, Markus Berndt

//   Permission is hereby granted, free of charge, to any person obtaining a
//   copy of this software and associated documentation files (the "Software"),
//   to deal in the Software without restriction, including without limitation
//   the rights to use, copy, modify, merge, publish, distribute, sublicense,
//   and/or sell copies of the Software, and to permit persons to whom the
//   Software is furnished to do so, subject to the following conditions:

//   The above copyright notice and this permission notice shall be included
//   in all copies or substantial portions of the Software.

//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//   DEALINGS IN THE SOFTWARE.


#ifndef __NKADIRECTION_H__
#define __NKADIRECTION_H__


#include "NOX.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Direction_Generic.H"
#include "NOX_GlobalData.H"

#include "NKA.H"


class NKADirection : public NOX::Direction::Generic {


public:

  NKADirection(const Teuchos::RCP<NOX::GlobalData>&, 
	       Teuchos::ParameterList&, const NOX::Epetra::Vector&);

  ~NKADirection();

  bool reset (const Teuchos::RCP<NOX::GlobalData>&, 
	      Teuchos::ParameterList&);
  
  bool compute (NOX::Abstract::Vector&, NOX::Abstract::Group&,
		const NOX::Solver::Generic&);

  bool compute (NOX::Abstract::Vector&, NOX::Abstract::Group&,
		const NOX::Solver::LineSearchBased&);

private:
  
  void throwError(const string&, const string&);

  // Printing Utils
  Teuchos::RCP<NOX::Utils> utils;

  // Global data pointer.  
  // Keep this so any stored parameter list remains valid.
  Teuchos::RCP<NOX::GlobalData> globalDataPtr;

  Teuchos::ParameterList *paramPtr; 

  nka *state;

  Teuchos::RCP<NOX::Abstract::Vector> tmpVecPtr;

};

#endif
