

/** \file   fasp_solver.cpp
 *  \brief  solver by FASP
 *  \author Kai Yang
 *  \date 09/26/2013
 */

#include<iostream>
#include<math.h>
#include "fasp_solver.h"
#include <boost/assign/list_of.hpp>
//#include <boost/numeric/itl/itl.hpp>


using namespace std;
using namespace dolfin;
using namespace FASPFSI;

static std::vector<int> *base_array2;
bool compar2(int a, int b)
{
   return ((*base_array2)[a] < (*base_array2)[b]);
}

/**
 *  \fn FaspSolver()
 *  \brief Class Constructor
 */

FaspSolver::FaspSolver(){
    fasp_dcsr_null(&FS);
    fasp_dvec_null(&b_FS);
    fasp_dvec_null(&u_FS);
    fasp_dvec_null(&blk_u_FS);
    fasp_dvec_null(&blk_b_FS);
  
    blkFS.brow=0;
    blkFS.bcol=0;
    blkFS.blocks=NULL;
    
    char inputfile[] = "ini/ns.dat";
    fasp_ns_param_input(inputfile,&inparam);
    fasp_ns_param_init(&inparam, &itparam, &amgparam, &iluparam, &schparam);
}


/**
 *  \fn ~FaspSolver()
 *  \brief Class Destructor
 */

FaspSolver::~FaspSolver(){
    fasp_dvec_free(&b_FS);
    fasp_dvec_free(&u_FS);
    fasp_dvec_free(&blk_b_FS);
    fasp_dvec_free(&blk_u_FS);
    fasp_dcsr_free(&FS);
  
  if(blkFS.blocks!=NULL){
    //clean up the block matrix
    for(int i=0;i<blkFS.brow;i++)
      for(int j=0;j<blkFS.bcol;j++){
          fasp_dcsr_free(blkFS.blocks[i*blkFS.brow+j]);
          fasp_mem_free(blkFS.blocks[i*blkFS.brow+j]);
      }
    
    fasp_mem_free(blkFS.blocks);

  }
}

/**
 * \fn solve
 * solve linear system
 */
int FaspSolver::solve(EigenVector &x,EigenVector &b,bool solver_type, const int output)
{
    if(b.size() > b_FS.row){
        if(b_FS.val!=NULL&&b_FS.row>0)
            fasp_mem_free(b_FS.val);
        b_FS.row = b.size();
        fasp_dvec_alloc(b_FS.row,&b_FS);
    }
    else{
        b_FS.row = b.size();
    }
    
    double*bval = (double *)b.data();
    for(int i=0;i<b_FS.row;i++){
        b_FS.val[i] = bval[i];
	}

    assemble_blk();
  
  if(output){
  fasp_matrix_write("data/A.dat",blkFS.blocks[0],1);
  fasp_matrix_write("data/Bt.dat",blkFS.blocks[1],1);
  fasp_matrix_write("data/B.dat",blkFS.blocks[2],1);
  fasp_matrix_write("data/C.dat",blkFS.blocks[3],1);
  fasp_dvec_write("data/rhs.dat",&blk_b_FS);
  }
    print_level=10;
    //fasp_dcoo_write("M4FS.mtx", &FS);
    fasp_dvec_alloc(FS.row, &u_FS);
    fasp_dvec_alloc(FS.row, &blk_u_FS);

    if(solver_type){
        //dCSRmat A;
        //dCSRmat ptrA = fasp_format_bdcsr_dcsr(&blkFS);
        //double abs_error,rel_error;
        
        //fasp_dcsr_trans(&ptrA,&A);
        //fasp_solver_umfpack(&A, &blk_b_FS,&blk_u_FS, print_level);
        //fasp_dcsr_free(&ptrA);
        
        //%%%%%%%%%%% fasp-ns solve  %%%%%%%%%%%%%%%
        int flag = fasp_solver_bdcsr_krylov_navier_stokes(&blkFS, &blk_b_FS, &blk_u_FS, &itparam, &amgparam, &iluparam, &schparam);
        //recover solution
        for(int i=0;i<gidx_u.size();i++){
            u_FS.val[gidx_u[i]] = blk_u_FS.val[i];
        }
        for(int i=0;i<gidx_p.size();i++){
            u_FS.val[gidx_p[i]] = blk_u_FS.val[n_v+i];
        }
    }
    else{
        //%%%%%%%%%%% fasp-ns solve with pressure %%%%%%%%%%%%%%%
        int flag = fasp_solver_bdcsr_krylov_navier_stokes_schur_complement_with_pressure_mass(&blkFS, &blk_b_FS, &blk_u_FS, &itparam, &amgparam, &iluparam, &schparam,&M_P);
        //recover solution
        for(int i=0;i<gidx_u.size();i++){
            u_FS.val[gidx_u[i]] = blk_u_FS.val[i];
        }
        for(int i=0;i<gidx_p.size();i++){
            u_FS.val[gidx_p[i]] = blk_u_FS.val[n_v+i];
        }
    //%%%%%%%%%%%%%%%   direct solve  %%%%%%%%%%%%%
        /*
        dCSRmat A;
        dCSRmat *ptrA = &A;
        double abs_error,rel_error;

        fasp_dcsr_trans(&FS,ptrA);
        fasp_solver_umfpack(ptrA, &b_FS,&u_FS, print_level);
        fasp_dcsr_free(ptrA);
         */
    }

    // output the result from u_FS to x
    double*xval = (double *)x.data();
    
    assert(n_FS==x.size());
    
    for(int i = 0;i < n_FS;i++){
        xval[i] = u_FS.val[i];
    }
    fasp_dvec_free(&u_FS);
    return 1;
}



void FaspSolver::assign_idx(FunctionSpace * W)
{
    const dolfin::la_index n0 = W->dofmap()->ownership_range().first;
    const dolfin::la_index n1 = W->dofmap()->ownership_range().second;
    const dolfin::la_index num_dofs = n1 - n0;
    //Assign velocity and pressure indices
    std::vector<std::size_t> component(1);
    component[0] = 0;
    std::shared_ptr<GenericDofMap> dofmap_u =
    W->dofmap()->extract_sub_dofmap(component,
                                    *(W->mesh()));
    component[0] = 1;
    std::shared_ptr<GenericDofMap> dofmap_p =
    W->dofmap()->extract_sub_dofmap(component,
                                    *(W->mesh()));
    gidx_u.clear();
    gidx_p.clear();
    gidx_u = dofmap_u->dofs();
    gidx_p = dofmap_p->dofs();
    std::sort(gidx_u.begin(), gidx_u.end());
    std::sort(gidx_p.begin(), gidx_p.end());
    
}

void FaspSolver::set_pressure_mass(EigenMatrix &P)
{
    int i,row,col,*IA, *JA;
    double *val;
    unsigned int nnz;
    // change format of A
    nnz = P.nnz();
    row = P.size(0);
    col = P.size(1);
    dCSRmat M;
    //Assign row and col
    M.col = col;
    M.IA = (int*)fasp_mem_calloc(row+1,sizeof(int));
    M.row = row;
    
    //Assign nnz
    M.JA = (int*)fasp_mem_calloc(nnz,sizeof(int));
    M.val = (double*)fasp_mem_calloc(nnz,sizeof(double));
    M.nnz = nnz;
    
    IA = M.IA;
    JA = M.JA;
    val = M.val;
    n_FS = row;
    
    auto tp=P.data();
    const int*IA_tmp = std::get<0>(tp);
    for(i = 0; i< row+1; i++){
        IA[i] = IA_tmp[i];
    }
    const int*JA_tmp = std::get<1>(tp);
    for(i = 0; i< nnz;i++){
        JA[i] = JA_tmp[i];
    }
    
    double *val_tmp = (double *)std::get<2>(tp);
    for(i = 0; i< nnz;i++){
        val[i] = val_tmp[i];
    }
    n_p = gidx_p.size();
    Ip = (INT*)fasp_mem_calloc(n_p,sizeof(INT));
    
    for(int i=0;i<gidx_p.size();i++){
        Ip[i] = gidx_p[i];
    }
    
    fasp_dcsr_getblk(&M,Ip,Ip,n_p,n_p,&M_P);
    //fasp_dcoo_write("P.dat",&M_P);
    fasp_mem_free(M.IA);
    fasp_mem_free(M.JA);
    fasp_mem_free(M.val);
    fasp_mem_free(Ip);
}

/**
 *\fn fenicsinput()
 *\brief get matrices, RHSs and interface mapping from fenics
 *\param A :fluid matrix
 *\param B :structure matrix
 */

int FaspSolver::set_operator(EigenMatrix &A)
{
    int i,row,col,*IA, *JA;
    double *val;
    unsigned int nnz;
    // change format of A
    nnz = A.nnz();
    row = A.size(0);
    col = A.size(1);
    
    //Assign row and col
    FS.col = col;
    if(row>FS.row){
        if(FS.IA!=NULL)
            fasp_mem_free(FS.IA);
        FS.IA = (int*)fasp_mem_calloc(row+1,sizeof(int));
      FS.row = row;
    }
    else{
      FS.row = row;
    }
    //Assign nnz
    if(nnz>FS.nnz){
        if((FS.JA!=NULL)&&(FS.row>0))
            fasp_mem_free(FS.JA);
        if(FS.val!=NULL&&(FS.row>0))
            fasp_mem_free(FS.val);
      FS.JA = (int*)fasp_mem_calloc(nnz,sizeof(int));
      FS.val = (double*)fasp_mem_calloc(nnz,sizeof(double));
      FS.nnz = nnz;
    }
    else{
      FS.nnz = nnz;
    }
    
    IA = FS.IA;
    JA = FS.JA;
    val = FS.val;
    n_FS = row;
    
    
    auto tp=A.data();
    const int*IA_tmp = std::get<0>(tp);
    for(i = 0; i< row+1; i++){
      IA[i] = IA_tmp[i];
    }
    const int*JA_tmp = std::get<1>(tp);   
    for(i = 0; i< nnz;i++){
      JA[i] = JA_tmp[i];
    }

    double *val_tmp = (double *)std::get<2>(tp);
    for(i = 0; i< nnz;i++){
      val[i] = val_tmp[i];
    }
    
    return 1;
}


// given FS, we assemble blkFS in block_dCSRmat.
int FaspSolver::assemble_blk(){

  // here we split it into 2 by 2 blocks for velocity and pressure.
  
  //make indices, note that index start from 1 instead of 0

  int count=0,i,count_v=0,count_p=0;
  
  blkFS.brow=2;
  blkFS.bcol=2;
  blkFS.blocks=(dCSRmat**)fasp_mem_calloc(blkFS.brow*blkFS.bcol, sizeof(dCSRmat*));

  dCSRmat *A, *B, *Bt, *C;
  A = (dCSRmat*) fasp_mem_calloc(1,sizeof(dCSRmat));
  B = (dCSRmat*) fasp_mem_calloc(1,sizeof(dCSRmat));
  Bt = (dCSRmat*) fasp_mem_calloc(1,sizeof(dCSRmat));
  C = (dCSRmat*) fasp_mem_calloc(1,sizeof(dCSRmat));

  blkFS.blocks[0] = A;
  blkFS.blocks[1] = Bt;
  blkFS.blocks[2] = B;
  blkFS.blocks[3] = C;
  
  n_p = gidx_p.size();
  n_v = gidx_u.size();
  Iv = (INT*)fasp_mem_calloc(n_v,sizeof(INT));
  Ip = (INT*)fasp_mem_calloc(n_p,sizeof(INT));

  //printf("debug: parameters: Iv:%d,Ip:%d, size of pi: %d, n_FS:%d\n",n_v,n_p,p_indicator.size(),n_FS);
  for(int i=0;i<gidx_u.size();i++){
        Iv[i] = gidx_u[i];
  }
    
    for(int i=0;i<gidx_p.size();i++){
        Ip[i] = gidx_p[i];
    }

  fasp_dcsr_getblk(&FS,Iv,Iv,n_v,n_v,A);
  fasp_dcsr_getblk(&FS,Iv,Ip,n_v,n_p,Bt);
  fasp_dcsr_getblk(&FS,Ip,Iv,n_p,n_v,B);
  fasp_dcsr_getblk(&FS,Ip,Ip,n_p,n_p,C);
    
#ifdef _DEBUG
  printf("info of block matrix assembling:\n");
  printf("size of FS:%d %d\n", FS.row,FS.col);
  printf("size of A:%d %d\n", A->row,A->col);
  printf("size of B:%d %d\n", B->row,B->col);
  printf("size of b':%d %d\n", Bt->row,Bt->col);
  printf("size of C:%d %d\n", C->row,C->col);
#endif
  dCSRmat**FS_blocks=blkFS.blocks;
  FS_blocks[0]=A;
  FS_blocks[1]=Bt;
  FS_blocks[2]=B;
  FS_blocks[3]=C;

  dCSRmat* temp_dcsr;
#ifdef _DEBUG
  printf("n_v:%d,n_p:%d\n",n_v,n_p);
#endif
    
  //reorder the RHS
    blk_b_FS = fasp_dvec_create(n_FS);
    for(int i=0;i<gidx_u.size();i++){
        blk_b_FS.val[i] = b_FS.val[gidx_u[i]];
    }
    for(int i=0;i<gidx_p.size();i++){
        blk_b_FS.val[n_v+i] = b_FS.val[gidx_p[i]];
    }
    fasp_mem_free(Iv);
    fasp_mem_free(Ip);
    return 1;
}
