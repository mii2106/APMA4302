import numpy as np
import subprocess
import pytest
import time

class TestExpx:
    """Test suite for single and multiprocessor PETSc expx code"""
    
    def setup_method(self):
        """Setup test parameters"""
        self.x_values = [-709., -500.,-100., -50., -10.0, -1., -0.1, 0. , 0.1, 1., 10., 50., 100., 500., 709.]
        self.N = lambda x: int(4 * np.abs(x)  + 30)
        self.num_processors = [1, 2, 4, 8]
        self.tolerance = 1e-6
    
    def run_expx_command(self, x, np_proc=1):
        """Run the expx PETSc executable with given parameters"""
        if np_proc == 1:
            cmd = ['./expx', '-x', str(x), '-N', str(self.N(x))]
        else:
            cmd = ['mpirun', '-np', str(np_proc), './expx', '-x', str(x), '-N', str(self.N(x))]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            if result.returncode != 0:
                raise RuntimeError(f"Command failed: {result.stderr}")
            return result.stdout
        except subprocess.TimeoutExpired:
            raise RuntimeError(f"Command timed out for x={x}, np={np_proc}")
    
    def parse_output(self, output):
        """Parse the output from expx executable to extract result"""
        # Assuming the output contains a line like "Result: 2.718281828"
        for line in output.strip().split('\n'):
            if 'Result:' in line or 'exp(' in line:
                return float(line.split()[-1])
        raise ValueError("Could not parse result from output")
    
    
if __name__ == "__main__":

    # check compilation
    print("Compiling expx code...")
    try:
        subprocess.run(['make', 'clean'], check=True)
        subprocess.run(['make', 'expx'], check=True)
        print("Compilation successful!")
    except subprocess.CalledProcessError as e:
        print(f"Compilation failed: {e}")
        exit(1)

    # Run tests directly
    test_suite = TestExpx()
    test_suite.setup_method()
    
    print("Running basic functionality tests...")
    try:
        for x in test_suite.x_values:   
            print(f"Testing with x = {x}, N = {test_suite.N(x)}")
            for nprocs in test_suite.num_processors:
                result = test_suite.run_expx_command(x, np_proc=nprocs)  
                print(f"\tnp={nprocs}, result={result}")  
        print("All tests passed!")
    except Exception as e:
        print(f"Test failed: {e}")