/**
 *  @brief  Implementation of bem matrix with H-matrix preconditioner.
 *  @file  bemmatrix.h
 *  @author  M.Bebendorf
 */

#ifndef BEMMATRIX_H_AHMED
#define BEMMATRIX_H_AHMED

#include "matrix.h" // base class


namespace AHMED {

  // forward declarations:
  class blcluster;
  template<typename> class mblock;

  /**
   *  @brief Class storing a H-matrix and its preconditioner and
   *         in H-matrix format used for BEM
   *  @tparam T1,T2 Data types clustered in the matrix block cluster tree
   */
  template<class T1,class T2>
  struct BEMMatrix : public Matrix<dcomp>
  {
    mblock<dcomp>** blcks;         //!<  The blocks of the stiffness matrix A
    bemblcluster<T1,T2>* blclTree;     //!< the underlying block cluster tree

    mblock<dcomp>** U;   //!<  The blocks of approximate Cholesky decomp U
    blcluster *blU;      //!< the underlying block cluster tree of U

    dcomp* xf;           //!< an auxiliary array required for precond_apply()

    //----------------------------------------------------------------  

    /**
     *  @brief Constructor for initialsation
     *  @param N Dimension of the matrix
     *  @param bl block cluster tree for the H-matrix
     */
  BEMMatrix(const unsigned N, bemblcluster<T1,T2>* bl) 
    : Matrix<dcomp>(N, N), blclTree(bl), U(NULL), blU(NULL) { 
      allocmbls(bl, blcks);
      xf = new dcomp[N]; 
    }

    /// destructor
    virtual ~BEMMatrix() {
      freembls(blclTree, blcks); if(U) freembls(blU, U);
      delete [] xf;
      if (blU) delete blU;
      delete blclTree;
    }

    /**
     *  @brief Adds a scaled matrix-vector multiplication (\f$y \gets y + \alpha A x\f$).
     *  @param d A scalar value.
     *  @param x A pointer of size n (columns).
     *  @param y A pointer of size m (rows).
     */
    void amux(const dcomp& d, const dcomp* x, dcomp* y) const {
      mltaGeHVec(d, blclTree, blcks, x, y);
    }

    /**
     *  @brief Applies preconditioner to a vector (\f$xf \gets x, \, xf \gets
     *                                            Cxf, \, x \gets xf, \, C \approx A^{-1}\f$).
     *  @param T A pointer of size n (columns).
     */
    void precond_apply(dcomp* x) const {
      if (blU) {
	std::copy_n(x, blU->getn1(), xf);
	HCholesky_solve(blU, U, xf);
	std::copy_n(xf, blU->getn1(), x);
      }
    }
  };
} // namespace AHMED

#endif // BEMMATRIX_H_AHMED
