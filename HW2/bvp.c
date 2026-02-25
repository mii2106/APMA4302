//STARTWHOLE
static char help[] = "Solve a tridiagonal system of arbitrary size.\n"
"Option prefix = tri_.\n";

#include <petsc.h>
#include <petscviewerhdf5.h>

int main(int argc,char **args) {
    Vec         u, usol,f;
    Mat         A;
    KSP         ksp;
    PetscInt    m = 4, i, Istart, Iend, j[3], bc[2];
    PetscReal   v[3], errnorm=0.0, unorm=0.0, x ,u_i ,f_i ;
    PetscViewer viewer;

    // agregado
    PetscReal gamma = 0.0, c=0.0, h = 1;
    PetscInt  k=5;


    PetscCall(PetscInitialize(&argc,&args,NULL,help));

    PetscOptionsBegin(PETSC_COMM_WORLD,"bvp_","options for bvp",NULL);
    PetscCall(PetscOptionsInt("-m","dimension of linear system","bvp.c",m,&m,NULL));
    PetscCall(PetscOptionsReal("-gamma","dimension of linear system","bvp.c",gamma,&gamma,NULL));
    PetscCall(PetscOptionsInt("-k","dimension of linear system","bvp.c",k,&k,NULL));
    PetscCall(PetscOptionsReal("-c","dimension of linear system","bvp.c",c,&c,NULL));
    PetscOptionsEnd();


    PetscCall(VecCreate(PETSC_COMM_WORLD,&u));
    PetscCall(VecSetSizes(u,PETSC_DECIDE,m));
    PetscCall(VecSetFromOptions(u));
    PetscCall(VecDuplicate(u,&f));
    PetscCall(VecDuplicate(u,&usol));

    PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
    PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,m));
    PetscCall(MatSetOptionsPrefix(A,"a_"));
    PetscCall(MatSetFromOptions(A));
    PetscCall(MatSetUp(A));
    PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
    h = 1/(PetscReal)(m-1);
    bc[0] = 0; 
    bc[1] = m-1;
    for (i=Istart; i<Iend; i++) {
        x = i*h;
        u_i = PetscSinReal(k*PETSC_PI*x) + c*(x-0.5)*(x-0.5)*(x-0.5);
        PetscCall(VecSetValues(usol,1,&i,&u_i,INSERT_VALUES));
        if (i == 0 || i == m-1) {
            PetscScalar one = 1.0;
            PetscCall(MatSetValues(A,1,&i,1,&i,&one,INSERT_VALUES));
        } else {
            v[0] = -1.0/(h*h); v[1] = 2.0/(h*h)+gamma; v[2] = -1.0/(h*h);
            j[0] = i-1; j[1] = i; j[2] = i+1;
            PetscCall(MatSetValues(A,1,&i,3,j,v,INSERT_VALUES));
            f_i = (k*k*PETSC_PI*PETSC_PI + gamma)*PetscSinReal(k*PETSC_PI*x)
                + gamma*c*(x-0.5)*(x-0.5)*(x-0.5) - 6*c*(x-0.5);
            PetscCall(VecSetValues(f,1,&i,&f_i,INSERT_VALUES));
        }
    }
    PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
    PetscCall(VecAssemblyBegin(usol));
    PetscCall(VecAssemblyEnd(usol));
    PetscCall(VecAssemblyBegin(f));
    PetscCall(VecAssemblyEnd(f));

    bc[0] =0;
    bc[1] = m-1;
    PetscCall(MatZeroRowsColumns(A, 2, bc, 1.0, usol, f));
    PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp));
    PetscCall(KSPSetOperators(ksp,A,A));
    PetscCall(KSPSetFromOptions(ksp));
    PetscCall(KSPSolve(ksp,f,u));


    // output the solution, rhs, and exact solution to an HDF5 file
    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, "bvp_solution.h5",
    FILE_MODE_WRITE, &viewer));
    PetscCall(PetscObjectSetName((PetscObject) usol, "uexact"));
    PetscCall(PetscObjectSetName((PetscObject) f, "f"));
    PetscCall(PetscObjectSetName((PetscObject) u, "u"));
    PetscCall(VecView(f, viewer));
    PetscCall(VecView(u, viewer));
    PetscCall(VecView(usol, viewer));
    PetscCall(PetscViewerDestroy(&viewer));

    PetscCall(VecAXPY(u,-1.0,usol));
    PetscCall(VecNorm(u,NORM_2,&errnorm));
    PetscCall(VecNorm(usol,NORM_2,&unorm));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
    "realtive error for m = %d system is = %.1e\n",m,errnorm/unorm));

    PetscCall(VecDestroy(&u));
    PetscCall(VecDestroy(&f));
    PetscCall(VecDestroy(&usol));
    PetscCall(MatDestroy(&A));
    PetscCall(KSPDestroy(&ksp));
    PetscCall(PetscFinalize());
    return 0;
}
//ENDWHOLE
