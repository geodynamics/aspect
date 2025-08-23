/*
  Copyright (C) 2011 - 2025 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/
#ifndef aspect_block_stokes_preconditioner_h
#define aspect_block_stokes_preconditioner_h

namespace aspect
{

  namespace internal
  {
    /**
     * Implement the block Schur preconditioner
     * (A B^T; 0 S)^{-1}.
     */
    template <class AInvOperator, class SInvOperator, class BTOperator,  class VectorType>
    class BlockSchurPreconditioner : public Subscriptor
    {
      public:
        /**
         * @brief Constructor
         * @param A_inverse_operator Approximation of the inverse of the velocity block.
         * @param S_inverse_operator Approximation for the inverse Schur complement.
         * @param BToperator Operator for the B^T block of the Stokes system.
         */
        BlockSchurPreconditioner (
          const AInvOperator                         &A_inverse_operator,
          const SInvOperator                         &S_inverse_operator,
          const BTOperator                           &BT_operator);

        /**
         * Matrix vector product with this preconditioner object.
         */
        void vmult (VectorType       &dst,
                    const VectorType &src) const;

      private:
        /**
         * References to the various operators this preconditioner works with.
         */

        const AInvOperator                     &A_inverse_operator;
        const SInvOperator                     &S_inverse_operator;
        const BTOperator                       &BT_operator;
    };


    template <class AInvOperator, class SInvOperator, class BTOperator,  class VectorType>
    BlockSchurPreconditioner<AInvOperator, SInvOperator, BTOperator, VectorType>::
    BlockSchurPreconditioner (
      const AInvOperator                         &A_inverse_operator,
      const SInvOperator                         &S_inverse_operator,
      const BTOperator                           &BT_operator)
      :
      A_inverse_operator (A_inverse_operator),
      S_inverse_operator (S_inverse_operator),
      BT_operator        (BT_operator)
    {}



    template <class AInvOperator, class SInvOperator, class BTOperator, class VectorType>
    void
    BlockSchurPreconditioner<AInvOperator, SInvOperator, BTOperator, VectorType>::
    vmult (VectorType       &dst,
           const VectorType &src) const
    {
      typename VectorType::BlockType utmp(src.block(0));

      // first apply the Schur Complement inverse operator.
      {
        S_inverse_operator.vmult(dst.block(1),src.block(1));
        dst.block(1) *= -1.0;
      }

      // apply the top right block
      {
        BT_operator.vmult(utmp, dst.block(1)); // B^T or J^{up}
        utmp *= -1.0;
        utmp += src.block(0);
      }

      A_inverse_operator.vmult(dst.block(0), utmp);
    }
  }
}

#endif
