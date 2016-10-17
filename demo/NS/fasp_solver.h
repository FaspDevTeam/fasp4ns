/** \file fasp_solver.h
 *  \brief solver from FASP
 *  \author Kai Yang
 *  \date 09/26/2013
 * 
 */


#include <dolfin.h>
//#include "util.h"
#include <stdio.h>

extern "C"
{
#include "fasp.h"
#include "fasp_functs.h"
#include "fasp4ns.h"
#include "fasp4ns_functs.h"
}


using namespace dolfin;

namespace FASPFSI{

  /**
   * \class FaspSolver
   *
   * \brief fasp solver for the linear system of fsi problem
   *
   */
  class FaspSolver{
    
  public:

    /**
     * \fn FaspSolver()
     *
     * \brief Class Constructor
     *
     */
    FaspSolver();

    /**
     * \fn ~FaspSolver()
     *
     * \brief Class Destructor
     *
     */
    ~FaspSolver();
    
    /**
     * \fn int solve(bool)
     *
     * \brief solve the FSI system based on the F,S block and interface mapping
     * \param bool: true for blk solver, false for umfpack
     *
     */
    int solve(EigenVector &x,EigenVector &b,bool solver_type = false, const int output = 0);

    /**
     * \fn set_operator(GenericMatrix &,GenericMatrix &, GenericVector &, GenericVector &, std::vector<std::vector<size_t> >& )

     * \brief get matrices, RHSs and interface mapping from fenics
     * \param A: matrix A
     * \param b: the rhs
     *
     */
    int set_operator(EigenMatrix &A);
     
    void set_pressure_mass(EigenMatrix &P);
      
    /**
     * \fn assemble_blk()
     * \brief assemble global stiffness matrices(block) and RHSs
     * \param
     */
    int assemble_blk();
     
    void assign_idx(FunctionSpace * W);

//added by Chen Yuyan
    void print_FS(){
        std::cout<<FS.row<<","<<FS.col<<std::endl;
    }

/*    void read(char *fileA,
              char *fileB,
              char *fileC,
              char *filerhs){
        fasp_bdcsr_read("data/A.dat", "data/B.dat", "data/C.dat", "data/rhs.dat", &blkFS, &b_FS);
    }*/

 //   void solve(EigenVector &x, bool solver_type = false);

    void write_A(const char *filename){
        fasp_matrix_write(filename, &blkFS.blocks[0], 1);
    }
    
    void write_B(const char *filename){
		fasp_matrix_write(filename, &blkFS.blocks[2], 1);
    }

    void write_C(const char *filename){
        fasp_matrix_write(filename, &blkFS.blocks[3], 1);
    }

    void write_b(const char *filename){
        fasp_dvec_write(filename, &b_FS);
    }
  protected:
      INT *Iv,*Ip;
      dCSRmat FS;
      dvector b_FS;
      dvector blk_b_FS;
      dvector blk_u_FS;
      dvector u_FS;
      block_dCSRmat blkFS;
      dCSRmat M_P;

    std::vector<dolfin::la_index> gidx_u;
    std::vector<dolfin::la_index> gidx_p;
      
    int print_level;
    int n_FS;
    int n_v,n_p;

    std::vector<int> idx_A;
    std::vector<int> idx_B;
    std::vector<size_t> idx_map;

	input_ns_param     inparam;  // parameters from input files
	itsolver_ns_param  itparam;  // parameters for itsolver
	AMG_ns_param       amgparam; // parameters for AMG
	ILU_param       iluparam; // parameters for ILU
    Schwarz_param   schparam; // parameters for Schwarz
    
  };
    
   

}
