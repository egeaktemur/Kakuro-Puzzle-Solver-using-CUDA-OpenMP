# Kakuro-Puzzle-Solver-using-CUDA-OpenMP

## OpenMP
### Implementation (Algorithm Overview and Heuristics)
Implementation includes a backtracking algorithm that attempts all possible numbers (1 to 9) in a cell and moves on to the next cell if the current partial solution is still valid. If the partial solution becomes invalid, it backtracks to the previous cell and tries the next partially valid number. This process continues until a solution is found or all possibilities are exhausted.
The verifyRemainingNumbers() function checks whether the remaining numbers needed for the sum
can be achieved using the remaining empty cells, considering the constraints of the puzzle. It does so by calculating the minimum and maximum possible sums with the given constraints.
The wrongSolution() function checks whether the current partial solution is incorrect, considering the current sum, the direction of the sum, and the cells involved.
The solve kakuro() function is the main backtracking algorithm implementation. It recursively attempts all possible numbers in a cell and moves to the next cell if the partial solution is still valid.
The solution() function initializes the data structures, assigns cell values, and calls the solve kakuro() function to start the solving process.
### Execution times
Small board sizes: As the number of threads increases, the execution time tends to increase as well, except for some cases. This is likely due to the overhead of creating and managing threads.
Medium board sizes: The execution times for varying numbers of threads do not follow a clear pattern, with some cases improving as the thread count increases, while others show the opposite trend. This suggests that the optimal number of threads may depend on the specific puzzle.
Large board sizes: In general, the execution times tend to decrease as the number of threads increases, with a few exceptions. This suggests that parallelism becomes more beneficial as the problem size grows.
It’s worth noting that the optimal number of threads is not consistent across different problem sizes and specific puzzles which might be due to various factors such as thread management overhead or specific characteristics of the puzzle.
### Efficiency Tricks
#### Spatial and temporal locality
##### Bitset for spatial locality
A std::bitset is used to store the ”exists” variable in the verifyRemainingNumbers and wrongSolution functions. This data structure provides spatial locality by packing boolean values into a contiguous block of memory, allowing for faster access times and cache performance. 
##### DeepCopy for temporal locality
The deepCopy function is used to create a copy of the sol mat matrix before modifying it. This ensures that the original matrix remains unchanged, allowing for a more efficient backtracking algorithm.
##### Exploit spatial locality in loops
Loops in the verifyRemainingNumbers and wrongSolution functions have been designed in such a way that they access elements in a contiguous memory region. This helps to improve cache performance and overall efficiency.
#### Preprocessing
Preprocessing which sums belong to which cell: The which sums cell exists data structure is a 3D vector that maps each cell to the sum objects it belongs to. This preprocessing step helps to speed up the solving process by allowing the solver to quickly identify which sums are associated with a given cell.
#### Search Space Optimizations
##### Pruning
In the solve kakuro function, the search space is pruned by checking whether the current cell’s value is valid (i.e., it doesn’t cause any conflict with the given sums). This helps to reduce the total number of possibilities that need to be explored during the solving process.
##### Backtracking
The backtracking algorithm is used to explore the solution space. By incrementally building the solution and undoing any changes when a conflict is found, the algorithm can search for a solution more efficiently than a brute-force approach.
#### Parallelism
The code uses the OpenMP to parallelize the computation of valid numbers in solve kakuro function. To elaborate on that more while in the serial algorithm code it calculates the next number if the current one is not valid, in the parallel algorithm it calculates all of the numbers’ validness at the start in a parallel manner. While this precalculation seems efficient in large boards high thread numbers compared to 1 thread execution, it is still more inefficient than the serial algorithm since it creates overhead does unnecessary
calculations.
### Execution Instruction
Code is compiled using: ” g++ kakuro solver hw2.cpp -O3 -fopenmp”
It can be run by just ”./a.out” but one my want to edit the main for testing purposes.
## CUDA
### General Explanation of the Implementation and Used Algorithm
The primary algorithm used in this solver is Backtracking with Pruning. Backtracking is a common technique for solving constraint satisfaction problems such as this, where the goal is to find a solution that meets a series of constraints, and it involves trying out potential solutions and then undoing them if they don’t work out (i.e., ”backtracking”). Pruning is a technique that involves eliminating branches in the solution space that cannot possibly contain a valid solution. It is used here to reduce the number of potential solutions that need to be examined.
### Code Flow
An initial empty solution matrix (d sol mat) is created and placed in the parent variable (parent).
The code then loops over each cell of the matrix. For each cell, it creates 9 child solutions, each representing one of the nine possible values for the cell. These children are stored in the childs array.
The solve kakuro function is launched as a CUDA kernel to evaluate the validity of each child solution in parallel.
It checks each associated sum for the current cell if it is still valid after inserting the new number. If a solution is found to be invalid, it’s pruned (set to NULL in childs).
After all child solutions have been evaluated, the countpruned kernel is run to count the number of pruned solutions.
The eliminatepruned kernel is then run to remove the pruned solutions from the list of child solutions, moving the valid solutions to the parent array for the next iteration.
The process repeats until all cells have been filled, at which point all remaining solutions in the parent array are valid solutions to the Kakuro puzzle.
### Heuristics
The code uses a simple heuristic to speed up the process of finding solutions which is Pruning. By immediately eliminating any child solutions that do not meet the game’s rules, the code is able to significantly reduce the number of potential solutions that it needs to check. This greatly reduces the total computational cost of the solver.
### Efficiency Tricks
#### Parallelism
By implementing the algorithm on a GPU using CUDA, the code can evaluate multiple potential solutions simultaneously. This massively parallel approach greatly reduces the overall running time of the solver.
#### Memory Management
The code uses CUDA memory management functions to allocate and deallocate memory on the device. Using cudaMallocManaged allows the code to allocate memory that is accessible from both the host and the device, simplifying memory management.
#### Reusing Memory
The parent and childs pointers are reused at each level of the recursion, reducing the amount of memory allocation and deallocation required.
#### Coalesced Memory Access
The childs array is accessed in a coalesced manner, meaning that the threads in a warp access consecutive memory locations. This leads to efficient memory access and maximum memory bandwidth utilization.
#### Pruning
It uses a pruning mechanism to eliminate the wrong solutions from the computation process early (by the function wrongSolution). It eliminates branches that lead to an incorrect solution, reducing the number of possibilities to be checked.
#### Concurrency
The code also uses the function cudaDeviceSynchronize to synchronize the execution of threads, ensuring that all preceding commands have completed before proceeding. This is important in situations where the sequence of execution matters.
### Execution Example
nvcc kakuro solver backup.cu -o kakuro.out
./kakuro.out board3 1.kakuro
### Results
board3 1.kakuro: Number of solutions: 1 Time to generate: 3.2 ms
board3 2.kakuro: Number of solutions: 1 Time to generate: 2.9 ms
board3 3.kakuro: Number of solutions: 1 Time to generate: 2.8 ms
board4 1.kakuro: Number of solutions: 44 Time to generate: 35.2 ms
board4 2.kakuro: Number of solutions: 1 Time to generate: 11.3 ms
board5 1.kakuro: Number of solutions: 12 Time to generate: 70.4 ms
board5 2.kakuro: Number of solutions: 1 Time to generate: 12.8 ms
board20 1.kakuro: Number of solutions: 1 Time to generate: 1031.1 ms
### Further Improvements
In the version that I am uploading I am allocating memory for successor in parent and copying each parent to child. This method is slow compared to each child copying parent in their execution. However this implementation was giving me crashes in 20 1 so I will upload the first version.
