/* Struct linear solvers header */
#include "HYPRE_struct_ls.h"
//#include "/usr/local/hypre-2.11.1/src/hypre/include/HYPRE_struct_ls.h"
//#include "HYPRE_struct_mv.h"
//#include "_hypre_utilities.h"
//#include "_hypre_utilities.h"
//#include "vis.c"

//const unsigned int Nstens = 3;
//const unsigned int ndim  = 2;

 MPI_Comm                  MPI_CW;
 HYPRE_StructGrid            grid;
 HYPRE_StructStencil      stencil;
 HYPRE_StructMatrix             A;
 HYPRE_StructVector             b; // rhs
 HYPRE_StructVector             x; // initial sol
 HYPRE_StructSolver        solver;
 int ilower[2];
 int iupper[2];


void hypre_solver_( double *bb, double *xx )
{
  /************************************************/
  {
    /* Set the vector coefficients */
    HYPRE_StructVectorSetBoxValues(b, ilower, iupper, bb);
    HYPRE_StructVectorSetBoxValues(x, ilower, iupper, xx);

    /* This is a collective call finalizing the vector assembly. Vectors are now ''ready to be used'' */
    HYPRE_StructVectorAssemble(b);
    HYPRE_StructVectorAssemble(x);
    }

  /* solve the system */
  {

     /*setup solver and solve system */
     HYPRE_StructPFMGSetup(solver, A, b, x);
     HYPRE_StructPFMGSolve(solver, A, b, x);
    // HYPRE_StructPCGSetup(solver, A, b, x);
    // HYPRE_StructPCGSolve(solver, A, b, x);

     HYPRE_StructVectorGetBoxValues(x, ilower, iupper, xx);


  }

}

void hypre_allocate_( MPI_Fint *COMM, int *period, int *topo, double *values){
  /* values : stencil
     ncells : number of cell in the local domain
     topo : sx,sy,ex,ey
     COMM : communicateur
     myid : rang ?
     period : periodicity of the domain
   */

 MPI_CW = MPI_Comm_f2c(*COMM);

 // printf("MPI comm ok\n");

 // local topo = { sx, sy,ex, ey }
 int ilower[2] = {topo[0], topo[1]};
 int iupper[2] = {topo[2], topo[3]};

  /************************************************/

  /* 1. Set up a grid */
   {
   /* Create an empty 2D grid object */
   HYPRE_StructGridCreate(MPI_CW, 2, &grid);

   /* Add boxes to the grid */
   HYPRE_StructGridSetExtents(grid, ilower, iupper);

   /* Periodic boundaries */
   HYPRE_StructGridSetPeriodic(grid,period);

   /* This is a collective call finalizing the grid assembly, grid is now ''ready to be used'' */
   HYPRE_StructGridAssemble(grid);
   }

  /************************************************/

  /* 2. Define the discretization stencil */
   {
   HYPRE_StructStencilCreate( 2, 3, &stencil);

      /* Define the geometry of the stencil. Each represents a 123 relative offset (in the index space). */
      {
      int entry;
      /*                     C      N      E    */
      int offsets[3][2] = {{0,0}, {1,0}, {0,1}};

      /* Assign each of the 3 stencil entries */
      for (entry = 0; entry < 3; entry++)
          HYPRE_StructStencilSetElement(stencil, entry, offsets[entry]);
      }
   }

  /************************************************/

  /* 3. Set up a Struct Matrix */
   {
   /* Create an empty matrix object */
   HYPRE_StructMatrixCreate(MPI_CW, grid, stencil, &A);
   /* Symmetric matrix storage */
   HYPRE_StructMatrixSetSymmetric(A,1);

   /* Indicate that the matrix coefficients are ready to be set */
   HYPRE_StructMatrixInitialize(A);

   /* labels for the stencils entries, corrsponding to the offsets defined above */
   int stencil_indices[3] = {0,1,2};

   /* set boxes - matrix for each CPU */
   HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 3 , stencil_indices, values);

   /* This is a collective call finalizing the matrix assembly. Matrix is now ''ready to be used'' */
   HYPRE_StructMatrixAssemble(A);
   }
  // HYPRE_StructMatrixPrint("Amatrix", A, 0);

  /************************************************/

  /* 4. Set up Struct Vectors for b and x.  Each processor sets the vectors
   corresponding to its boxes. */
   {
   /* Create an empty vector object */
   HYPRE_StructVectorCreate(MPI_CW, grid, &b);
   HYPRE_StructVectorCreate(MPI_CW, grid, &x);

   /* Indicate that the vector coefficients are ready to be set */
   HYPRE_StructVectorInitialize(b);
   HYPRE_StructVectorInitialize(x);
   }
   /* 5. Set up and use a solver to solve the system */
   {
      /* Create an empty PFMG Struct solver */
      HYPRE_StructPFMGCreate(MPI_CW,  &solver);
      //HYPRE_StructPCGCreate(MPI_CW,  &solver);

      /*convergence tolerance */
      HYPRE_StructPFMGSetTol(solver, 1.0e-06);
      /*amount of info printed on screen */
      //HYPRE_StructPCGSetPrintLevel(solver, 0);

      /*diverse shit for PFMG solver */
      HYPRE_StructPFMGSetRAPType(solver, 0);
      HYPRE_StructPFMGSetRelaxType(solver, 1);
      HYPRE_StructPFMGSetNumPreRelax(solver, 1);
      HYPRE_StructPFMGSetNumPostRelax(solver, 1);
      HYPRE_StructPFMGSetSkipRelax(solver, 0);

    }
}


void hypre_deallocate_(){

  /* Destroy the PFMG solver to prevent memory leak - better way has to be found */
  HYPRE_StructPFMGDestroy(solver);
  // HYPRE_StructPCGDestroy(solver);

  HYPRE_StructGridDestroy(grid);
  HYPRE_StructStencilDestroy(stencil);
  HYPRE_StructMatrixDestroy(A);
  HYPRE_StructVectorDestroy(b);
  HYPRE_StructVectorDestroy(x);
}
