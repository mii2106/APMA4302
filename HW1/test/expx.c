#include <petsc.h>

int main(int argc, char **argv)
{
  PetscMPIInt rank, size;
  PetscInt    N = 10,i, base, rem,start, count,end   ;
  PetscReal   x = 1.0, xs, localsum = 0.0, globalsum = 0.0, B, cn, shift, n, fact;
  PetscReal approx, exact, relerr;
  PetscBool   invert = PETSC_FALSE;

  PetscCall(PetscInitialize(&argc, &argv, NULL,
    "Compute exp(x) with an N-term Taylor polynomial in parallel.\n"));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));

  PetscOptionsBegin(PETSC_COMM_WORLD, "", "options for expx", "");
  PetscCall(PetscOptionsReal("-x", "input to exp(x) function", NULL, x, &x, NULL));
  PetscCall(PetscOptionsInt ("-N", "number of Taylor terms (>0)", NULL, N, &N, NULL));
  PetscOptionsEnd();
 //check if N is valid
  if (N < 1) {
    if (rank == 0) PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Error!!!! N must be > 0\n"));
    PetscCall(PetscFinalize());
    return 0;
  }
  // for negative x
  xs = x;
  if (x < 0.0) { xs = -x; invert = PETSC_TRUE; }
  // split terms, and 
  base  = N/size;
  rem = N%size;
  start = rank*base+ (rank < rem ? rank : rem);
  count = base+(rank<rem ? 1 :0);
  end = start+count;

//PetscCall(PetscPrintf(PETSC_COMM_SELF,
 // "rank %d: start = %d, end = %d\n", (int)rank, (int)start, (int)end));
  
  if (count > 0) {
    // apply horners method for every bracket
    B = 0.0;
    for (n = end - 1; n >= start; --n) {
      PetscCall(PetscDTFactorial(n, &fact));  
      cn = 1.0 / fact;  
      B = cn + xs * B;
    }
   //multiply by the firsts x**start
    if (xs == 0.0) {
      localsum = (start == 0) ? 1.0 : 0.0;
    } else {
      shift = 1.0;
        for (i = 0; i < start; ++i)
            shift *= xs;
      localsum = shift * B;
    }
  }

  //Make the sum of the blocks
  PetscCallMPI(MPI_Reduce(&localsum, &globalsum, 1, MPIU_REAL, MPIU_SUM, 0, PETSC_COMM_WORLD));

  if (rank == 0) {
     approx = globalsum;
    if (invert) approx = 1.0 / approx;
    exact = exp(x);
    relerr = (exact != 0.0) ? (PetscReal)fabs((approx - exact) / exact)
                                      : (PetscReal)fabs(approx - exact);
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
      "exp(% .15e) = % .15e\n"
      "rel error = %.3e (%.3e * machine epsilon)\n",
      (double)x, (double)approx,
      (double)relerr, (double)(relerr/PETSC_MACHINE_EPSILON)));
  }
  PetscCall(PetscFinalize());
  return 0;
}
