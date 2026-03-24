static char help[] = "2D nonlinear reaction-diffusion problem with DMDA and SNES. Option prefix -rct_.\n\n";

#include <petsc.h>

typedef struct {
    PetscInt   p;
    PetscReal  gamma;
    PetscBool  linear_f;
} AppCtx;

extern PetscReal  ufunction(PetscReal, PetscReal);
extern PetscReal  d2ufunction(PetscReal, PetscReal);
extern PetscReal  fRHS(PetscReal, PetscReal, AppCtx *);

extern PetscErrorCode FormInitial(DM, Vec);
extern PetscErrorCode FormExact(DM, Vec);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo *, PetscReal **, PetscReal **, AppCtx *);
extern PetscErrorCode FormJacobianLocal(DMDALocalInfo *, PetscReal **, Mat, Mat, AppCtx *);

int main(int argc, char **args) {
    DM            da;
    SNES          snes;
    AppCtx        user;
    Vec           u, uexact;
    PetscReal     errnorm, uexactnorm;
    DMDALocalInfo info;
    PetscViewer   viewer;
    PetscCall(PetscInitialize(&argc, &args, NULL, help));

    user.p       = 2;
    user.gamma   = 1.0;
    user.linear_f = PETSC_FALSE;

    PetscOptionsBegin(PETSC_COMM_WORLD, "rct_", "options for reaction2d", "");
    PetscCall(PetscOptionsInt("-p", "exponent p in gamma*u^p","reaction2d.c", user.p, &user.p, NULL));
    PetscCall(PetscOptionsReal("-gamma", "coefficient gamma", "reaction2d.c", user.gamma, &user.gamma, NULL));
    PetscCall(PetscOptionsBool("-linear_f","use linear RHS f = -Laplacian(uexact) instead of full MMS RHS",
                                "reaction2d.c", user.linear_f, &user.linear_f, NULL));
    PetscOptionsEnd();
    // create the grid
    PetscCall(DMDACreate2d(PETSC_COMM_WORLD,
                           DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                           DMDA_STENCIL_STAR,
                           9, 9,
                           PETSC_DECIDE, PETSC_DECIDE,
                           1, 1, NULL, NULL, &da));
    PetscCall(DMSetFromOptions(da));
    PetscCall(DMSetUp(da));
    PetscCall(DMDASetUniformCoordinates(da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0));
    PetscCall(DMSetApplicationContext(da, &user));

    PetscCall(DMCreateGlobalVector(da, &u));
    PetscCall(VecDuplicate(u, &uexact));
    PetscCall(FormInitial(da, u));
    // exact solution
    PetscCall(FormExact(da, uexact));
    // functions for snes
    PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));
    PetscCall(SNESSetDM(snes, da));
    PetscCall(DMDASNESSetFunctionLocal(da, INSERT_VALUES,
              (DMDASNESFunctionFn *)FormFunctionLocal, &user));
    PetscCall(DMDASNESSetJacobianLocal(da,
              (DMDASNESJacobianFn *)FormJacobianLocal, &user));
    PetscCall(SNESSetFromOptions(snes));

    // Solveee!!
    PetscCall(SNESSolve(snes, NULL, u));

    // VTK just solution u and uexact
    PetscCall(PetscViewerVTKOpen(PETSC_COMM_WORLD, "reaction2d.vtr",
                                 FILE_MODE_WRITE, &viewer));
    PetscCall(PetscObjectSetName((PetscObject)u,      "u"));
    PetscCall(PetscObjectSetName((PetscObject)uexact, "uexact"));
    PetscCall(VecView(u,      viewer));
    PetscCall(VecView(uexact, viewer));
    PetscCall(DMView(da,      viewer));
    PetscCall(PetscViewerDestroy(&viewer));

    // For computing error
    PetscCall(VecNorm(uexact, NORM_2, &uexactnorm));
    PetscCall(VecAXPY(u, -1.0, uexact));  
    PetscCall(VecNorm(u, NORM_2, &errnorm));
    PetscCall(DMDAGetLocalInfo(da, &info));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "on %d x %d grid:  rel_error |u-uexact|_2/|uexact|_2 = %g\n",
        info.mx, info.my, errnorm / uexactnorm));

    // --- cleanup ---
    PetscCall(VecDestroy(&u));
    PetscCall(VecDestroy(&uexact));
    PetscCall(SNESDestroy(&snes));
    PetscCall(DMDestroy(&da));
    PetscCall(PetscFinalize());
    return 0;
}

// Exact solution
PetscReal ufunction(PetscReal x, PetscReal y) {
    PetscReal sigma = 0.3, x0 = 0.65, y0 = 0.65, amp = 1.0;
    PetscReal r2 = (x-x0)*(x-x0) + (y-y0)*(y-y0);
    return amp * PetscExpReal(-r2 / (sigma*sigma));
}
// nabla^2 u
PetscReal d2ufunction(PetscReal x, PetscReal y) {
    PetscReal sigma = 0.3, x0 = 0.65, y0 = 0.65, amp = 1.0;
    PetscReal r2 = (x-x0)*(x-x0) + (y-y0)*(y-y0);
    PetscReal expterm = PetscExpReal(-r2 / (sigma*sigma));
    return amp * expterm * 4.0 / (sigma*sigma) * (r2/(sigma*sigma) - 1.0);
}
// funciont f 
PetscReal fRHS(PetscReal x, PetscReal y, AppCtx *user) {
    PetscReal f = -d2ufunction(x, y);
    if (!user->linear_f) {
        PetscReal uex = ufunction(x, y);
        f += user->gamma * PetscPowReal(uex, (PetscReal)user->p);
    }
    return f;
}

// Initial guess: u = 0
PetscErrorCode FormInitial(DM da, Vec u) {
    PetscCall(VecSet(u, 0.0));
    return 0;
}

// Exact sol

PetscErrorCode FormExact(DM da, Vec uexact) {
    DMDALocalInfo  info;
    PetscReal      hx, hy, x, y, **au;

    PetscCall(DMDAGetLocalInfo(da, &info));
    hx = 1.0 / (info.mx - 1);
    hy = 1.0 / (info.my - 1);
    PetscCall(DMDAVecGetArray(da, uexact, &au));
    for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
        y = j * hy;
        for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
            x = i * hx;
            au[j][i] = ufunction(x, y);
        }
    }
    PetscCall(DMDAVecRestoreArray(da, uexact, &au));
    return 0;
}

// NON LINEAR RESIDUAL!!!!!!!!
PetscErrorCode FormFunctionLocal(DMDALocalInfo *info,
                                 PetscReal **u,
                                 PetscReal **FF,
                                 AppCtx    *user) {
    PetscReal h  = 1.0 / (info->mx - 1);  //Assuming dx=dy
    PetscReal h2 = h * h;
    for (PetscInt j = info->ys; j < info->ys + info->ym; j++) {
        PetscReal y = j * h;
        for (PetscInt i = info->xs; i < info->xs + info->xm; i++) {
            PetscReal x = i * h;
            // Dirichlet u = uexact
            if (i == 0 || i == info->mx-1 || j == 0 || j == info->my-1) {
                FF[j][i] = u[j][i] - ufunction(x, y);
            } else {
                // interiors
                PetscReal lap = (4.0*u[j][i] - u[j][i-1] - u[j][i+1]
                                              - u[j-1][i] - u[j+1][i]) / h2;
                PetscReal reaction = user->gamma
                                     * PetscPowReal(u[j][i], (PetscReal)user->p);
                FF[j][i] = lap + reaction - fRHS(x, y, user);
            }
        }
    }
    return 0;
}

// Jacobian  J(u)
PetscErrorCode FormJacobianLocal(DMDALocalInfo *info,
                                 PetscReal    **u,
                                 Mat           J,
                                 Mat           P,
                                 AppCtx       *user) {
    PetscReal    h  = 1.0 / (info->mx - 1);
    PetscReal    h2 = h * h;
    MatStencil   row, col[5];
    PetscReal    v[5];
    PetscInt     ncols;
    for (PetscInt j = info->ys; j < info->ys + info->ym; j++) {
        for (PetscInt i = info->xs; i < info->xs + info->xm; i++) {
            row.j = j;  row.i = i;

            if (i == 0 || i == info->mx-1 || j == 0 || j == info->my-1) {
                //boundary
                col[0].j = j;  col[0].i = i;
                v[0] = 1.0;
                PetscCall(MatSetValuesStencil(P, 1, &row, 1, col, v, INSERT_VALUES));
            } else {
                //interior
                PetscReal diag = 4.0 / h2
                    + user->gamma * (PetscReal)user->p
                      * PetscPowReal(u[j][i], (PetscReal)(user->p - 1));

                ncols = 0;
                // diagonal
                col[ncols].j = j;  col[ncols].i = i;
                v[ncols++] = diag;
                // left 
                col[ncols].j = j;  col[ncols].i = i-1;
                v[ncols++] = -1.0 / h2;
                // right 
                col[ncols].j = j;  col[ncols].i = i+1;
                v[ncols++] = -1.0 / h2;
                // bottom
                col[ncols].j = j-1;  col[ncols].i = i;
                v[ncols++] = -1.0 / h2;
                // top
                col[ncols].j = j+1;  col[ncols].i = i;
                v[ncols++] = -1.0 / h2;
                PetscCall(MatSetValuesStencil(P, 1, &row, ncols, col, v, INSERT_VALUES));
            }
        }
    }

    PetscCall(MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY));
    if (J != P) {
        PetscCall(MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY));
    }
    return 0;
}
