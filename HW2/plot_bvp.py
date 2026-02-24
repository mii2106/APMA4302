import numpy as np
from petsc4py import PETSc

def read_hdf5_vec(filename, vec_name):
    """
    Read PETSc HDF5 viewer output and convert to numpy arrays.
    
    Parameters:
    filename: str - path to the HDF5 file
    vec_name: str - name of the vector to read  
    
    Returns:
    numpy array containing the data
    """
    # Create a viewer for reading HDF5 files
    viewer = PETSc.Viewer().createHDF5(filename, 'r')   
    
    # Create a Vec to load the data
    vec = PETSc.Vec().create(comm=PETSc.COMM_WORLD)
    vec.setName(vec_name)
    vec.load(viewer)
   
    # Convert to numpy array
    array = vec.getArray()
    
    # Clean up
    vec.destroy()
    viewer.destroy()
    
    return array.copy()


def plot_bvp_solution(x, u_numeric, u_exact):
    """
    Plot the numerical and exact solutions of the BVP.
    
    Parameters:
    x: numpy array - grid points
    u_numeric: numpy array - numerical solution
    u_exact: numpy array - exact solution
    """
    import matplotlib.pyplot as plt
    
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax2 = ax1.twinx()

    ax1.plot(x, u_numeric, 'b-', label='Numerical Solution', linewidth=2)
    ax1.plot(x, u_exact, 'r--', label='Exact Solution', linewidth=2)
    ax1.set_xlabel('x', fontsize=14)
    ax1.set_ylabel('u(x)', fontsize=14)
    ax1.set_title('BVP Numerical vs Exact Solution', fontsize=16)
    ax1.legend(fontsize=12)
    ax1 .grid(True)

    ax2.plot(x, u_numeric - u_exact, 'g--', label='Error', linewidth=1)
    ax2.set_ylabel('Error', fontsize=14)
    ax2.legend(loc='lower right', fontsize=12)
    ax2.set_ylim(-np.max(np.abs(u_numeric - u_exact)) * 3., np.max(np.abs(u_numeric - u_exact)) * 3.)

    plt.show()

if __name__ == "__main__":
    # Example usage
    # read numerical and exact solutions from HDF5 files
    h5_filename = 'bvp_solution.h5'  # Update with your actual filename
    u = read_hdf5_vec(h5_filename, 'u') 
    u_exact = read_hdf5_vec(h5_filename, 'uexact')

    x = np.linspace(0, 1, len(u))  
    
    plot_bvp_solution(x, u, u_exact)