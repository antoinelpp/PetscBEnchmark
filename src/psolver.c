static char help[] = "Solves 2D of poisson equations.\n\n";

//#include <stdio.h>
//#include <stdlib.h>
#include <petsc.h>

//#undef __FUNCT__
//#define __FUNCT__ "main"

 DM               da;
 MatStencil row, col[5];
 Mat             A,B;
 Vec               b; // rhs
 Vec               x; // initial sol
 KSP             ksp;
 PetscViewer  viewer, viewer1;
 PC               pc;
 PetscErrorCode ierr;
 PetscInt       i,j,k,mx,my,xm,ym,xs,ys,num,xi,y,Iglob,J,Istart,Iend;
 Vec            b_local, x_local,bb_local, xx_local;
 Vec  x_exact;
 PetscScalar    **arrayb,**arrayx,**arraybb,**arrayxx;
 PetscScalar    **arrayx_exact;
 PetscScalar    stencilValue[5];
 PCType         pc_type_fromArg;
 KSPType        ksp_type_fromArg;
 //char pc_type[10], ksp_type[10];

 PetscRandom rctx;

void petsc_destroy_() {
  ierr = MatDestroy(&A);//CHKERRQ(ierr);
  ierr = VecDestroy(&b);
  ierr = VecDestroy(&x);
  //  ierr = PCDestroy(&pc);
  ierr = KSPDestroy(&ksp);//CHKERRQ(ierr);
  PetscFinalize();
}

void petsc_initialize_(MPI_Fint *COMM, int *M_global, int *N_global,
   int *lx, int *ly, int *M_proc, int *N_proc, char *ksp_type, char *pc_type) {

    //Initialize PETSC Environnement
    PetscInitializeNoArguments();

    //Print arguments
    PetscPrintf(PETSC_COMM_WORLD,"KSP_Type : %s\n", ksp_type);
    PetscPrintf(PETSC_COMM_WORLD,"PC_Type : %s\n", pc_type);

    if ( strcmp(ksp_type, "cg")){
      ksp_type_fromArg = KSPCG;
    } else if ( strcmp(ksp_type, "bcgs") ) {
      ksp_type_fromArg = KSPBCGS;

    }

    if ( strcmp(pc_type, "gamg")){
      pc_type_fromArg = PCGAMG;
    } else if ( strcmp(pc_type, "hypre") ) {
      pc_type_fromArg = PCHYPRE;

    }


    // Creat the DMDA object
    ierr = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC,
      DMDA_STENCIL_STAR,*M_global,*N_global,*M_proc,*N_proc,1,1,lx,ly,&da);
    //   CHKERRQ(ierr);
    //ierr = DMSetFromOptions(da);
    ierr = DMSetUp(da);

    //ierr = DMView(da, PETSC_VIEWER_STDOUT_WORLD);

    ierr = DMSetMatType(da,MATMPIAIJ);
    ierr = DMCreateMatrix(da, &A);
    ierr = DMCreateGlobalVector(da,&b);
    ierr = VecDuplicate(b,&x);

    ierr = DMDAGetInfo(da,0,&mx,&my,0,0,0,0,0,0,0,0,0,0);
    /* mx = Total domain size in x
       my = ..          ..       y
    */
    ierr = DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);//CHKERRQ(ierr);
    /* xs, ys = start position (ie: offset)
       xm, ym = Local (CPU) domain size in x, y
    */

    for (j=ys; j<ys+ym; j++) {
      for (i=xs; i<xs+xm; i++) {
        /*
        The stencil:
              1
          1  -4   1
              1
        */

        row.j = j; row.i = i;

        if (i==0 || i == mx -1   ) { // Lower or upper grounded electrod
           stencilValue[0] = -1;
           MatSetValuesStencil(A,1,&row,1,&row,stencilValue,INSERT_VALUES);
         }
         else {
           stencilValue[0] = 1;              col[0].i = i;   col[0].j = j-1;
           stencilValue[1] = 1;              col[1].i = i-1; col[1].j = j;
           stencilValue[2] =-4;              col[2].i = i;   col[2].j = j;
           stencilValue[3] = 1;              col[3].i = i+1; col[3].j = j;
           stencilValue[4] = 1;              col[4].i = i;   col[4].j = j+1;
           MatSetValuesStencil(A,1,&row,5,col,stencilValue,INSERT_VALUES);
         }

       }
    }

    // Assemble the Maxtrix from the DMDA and Stencil
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

    // ierr = MatView(A, PETSC_VIEWER_DRAW_WORLD);
    ierr = DMCreateGlobalVector(da, &b);
    ierr = DMCreateGlobalVector(da, &x);


    // Create the solver
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
    ierr = KSPSetOperators(ksp,A,A);


    ierr = KSPSetType(ksp,ksp_type_fromArg);

    // Create the Pre-Conditionner

    ierr = KSPGetPC(ksp,&pc);
    ierr = PCSetType(pc,pc_type_fromArg);
    ierr = PCSetUp(pc);

    //ierr = KSPSetFromOptions(ksp);



    ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
    KSPSetUp(ksp);

    //KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);

}

void petsc_solver_(double *bb, double *xx, int *log_view){
  int jj, ii;

if (*log_view == 1) {
    ierr = PetscLogDefaultBegin();
  }

  /** UNPACK the inputs
  to be changes !!
  */
  DMDAVecGetArray(da, b, &arrayb);
  DMDAVecGetArray(da, x, &arrayx);

  for(j=ys; j<ys+ym; j++) {
    for(i=xs; i<xs+xm; i++) {
     k = j-ys+ym*(i-xs);
     arrayb[j][i] = bb[k];
     arrayx[j][i] = xx[k];
   }
  }

  DMDAVecRestoreArray(da, b, &arrayb);
  DMDAVecRestoreArray(da, x, &arrayx);


  // Generate the random sol:
  PetscRandomCreate(PETSC_COMM_WORLD,&rctx);

  DMCreateGlobalVector(da, &x_exact);
  VecSetRandom(x_exact, rctx);
  PetscRandomDestroy(&rctx);

  // calculte the good RhS
  MatMult(A, x_exact, b);
  /****************************************************************/
  /* Solve the system */
  ierr = KSPSolve(ksp,b,x);
  /***************************************************************/

  {
    PetscInt its, reason;
    KSPGetConvergedReason(ksp,&reason);
    if (reason==KSP_DIVERGED_INDEFINITE_PC) {
     PetscPrintf(PETSC_COMM_WORLD,"\nDivergence because of indefinite preconditioner;\n");
     PetscPrintf(PETSC_COMM_WORLD,"Run the executable again but with '-pc_factor_shift_type POSITIVE_DEFINITE' option.\n");
    } else if (reason<0) {
     PetscPrintf(PETSC_COMM_WORLD,"\nOther kind of divergence: this should not happen.\n");
    } else {
     KSPGetIterationNumber(ksp,&its);
     PetscPrintf(PETSC_COMM_WORLD,"\nnumber of Iteration : %d.\n", its);
    }

    {
      Vec          error;
      PetscScalar  errorNorm, xnorm;

      VecDuplicate(x, &error);
      VecCopy(x, error);
      VecAXPY(error, -1, x_exact);
      VecNorm(error, NORM_2, &errorNorm);
      VecNorm(x, NORM_2, &xnorm);
      PetscPrintf(PETSC_COMM_WORLD, "norm (L2) of x : %6.4e \n", xnorm );
      PetscPrintf(PETSC_COMM_WORLD, "Error norm (L2) %6.4e, relative norm %6.4e \n", errorNorm, errorNorm / xnorm  );
      VecNorm(error, NORM_INFINITY, &errorNorm);
      VecNorm(x, NORM_INFINITY, &xnorm);
      PetscPrintf(PETSC_COMM_WORLD, "Error norm (Linf) %6.4e, relative norm %6.4e \n", errorNorm, errorNorm / xnorm  );

    }


  }

  if (*log_view == 1) {
    ierr = PetscLogView(PETSC_VIEWER_STDOUT_WORLD);
  }


  DMDAVecGetArray(da, x, &arrayxx);
  for (j=ys; j<ys+ym; j++) {
   for (i=xs; i<xs+xm; i++) {
       k = j-ys+ym*(i-xs);
       xx[k]=arrayxx[j][i];

     }
  }

  DMDAVecRestoreArray(da, x, &arrayxx);

}
