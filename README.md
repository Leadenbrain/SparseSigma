# Sparse Signature Matrix Implementation

## Summary:

In the evaluation of Differential Algebraic Equations (DAEs), one of the first steps is to determine the highest order of which every independent variable appears in each given equation. This produces a NxN matrix for a system of N equations. This technique is referred to as structural analysis. The technique used in my simulations is the Sigma Matrix technique for structural analysis. For systems with high interdependency of variables, this matrix is relative dense. However, if there is little coupling between differential equations, the final signature matrix will be mostly empty. Working with one such equation, for molecular diffusion, motivated the need for a sparse implementation of this sigma matrix. This reduces our matrix in memory from an N<sup>2</sup> problem to one that is only N, allowing significantly larger problems to be handled.

I successfully implemented such a technique, and demonstrate, in my code, how this can be applied to handle very large systems with little impact on performance and memory. In my code, the DAE functions are defined as a void function (which is cast into a custom type, SigmaMatrixFCN). This allows a signature matrix to be created as an object, which contains a function to compute and print the matrix.

## Technicals

The code is split into several files:

    hello.cc : The main file
    sparseSigma.cc/h : Class implementation of a sparse vector (sigma map) for individual equations
    sparseSigMat.cc/h : Class Implementation of a sparse matrix (sigma matrix) for whole DAE
    compSpSigma.cc : File implementing the computation of the sparse signature matrix

In `hello.cc`, the DAEs are encoded as void functions casting a template allows for overloading of our types further in our code, and plays into a huge part of the AD that happens later, but is not so relevant here. I further define a quick function that shows the code running as it should be, this will be explained further in the Results Section. In my main loop, I set two numbers to a variable (1 mil and 10) to show that the sparse matrix works, and that the results are correct; respectively. Next I call the function to test our code, twice on the molecular diffusion problem (to verify results and sparsity) and once on the Layne Watson mechanism, to ensure that dense problems are still solved correctly.

In `sparseSigma.cc/h` I created a class for handling the individual equations of the DAE in a SigmaVector class. These objects will be later joined together to make the entire signature matrix. There looks to be a lot going on in this class, but it's all basically the same, and relatively simple once examined closely. I create multiple overloaded constructors, these will allow me to efficiently enable a std::move constructor, which further increases the ability to get quick results with minimal memory overhead. For this SigmaVector class, I overload the individual arithmetic and mathematical operators. This allows a sigmavector object to be cast directly into the DAE and to compute. Instead of computing the resulting value of the DAE, it merely returns a 0 if an independent variable is present in an arithmetic operation. The only value that changes this from being 0 to anything else is when differentiation occurs. Thus differntiation will have a seperate overloaded function to increase the value of the entry by the order of the differentiation. Therefore if an array of SigmaVectors is fed into the DAE, we can determine the resulting signature matrix without any additional implementation. This will further allow this code to implement an automatic differentiation capacity later if need be.

In `sparseSigMat.cc/h` I created a class for the sparse signature matrix. I put the actual implementation of the computation in `compSpSigma.cc`, so the `sparseSigMat` class files really just hold the constructor, variables, and destructor for my signature matrix.

In `compSpSigma.cc` I implement both the computation of the sparse signature matrix, and the function to print out the resulting matrix that is produced. The function for computing the matrix is relatively simple and straightforward. It starts with the declaration of two pointers of our SigmaVector class, one will be used for handling the independent varialbs (y) and one will be used for our sparse vectors (f). The DAE is then executed, with the SigmaVectors, y and f, linked in to incorporate the values they get from the operator overloading functions. Provided there is no error, the sigma_matrix is then initialized with all values from each equation in f, by looping over all f's and inserting their maps into the array entry. The print function is very basic, I just loop over the numbers both ways and print out the entry if present, or 0 if not. This function is just for the project, and I personally restrict it to 10. I would generally want to loop those automatically only over the values, rather than looping over n_ twice, but this was just quick and easy to implement, and will be removed as soon as the project is handed in.

## Results

I set up three different test cases in this project:
    
    1) A simulation of a molecular diffusion problem (a sparse signature matrix) at N=1 million. This will showcase the ability to handle large sparse problems. We expect to see 3000001 non-zero entries
    2) A simulation of a molecular diffusion problem but this time a N=10. This is just to show that the result is correct. We expect tridiagonality with 1 in the diagonal, 0 in the off diagonal, and no entry in the rest of the matrix. The final row will be the same as the second last row, but missing the entry of '1'. (This is due to the boundary condition of the problem).
    3) A simulation of the Layne Watson mechanism, a DAE benchmark. This is a dense problem, so ensure the results are correct validates that our sparse matrix is truly representing the actual dense structure. We expect 1 in every entry of the first row and 0 in every other entry.

When I compile and run my main function (with g++ -o main *.cc -O2 [note: additionally you can just run through .vscode, but make sure to include all source files, and the -O2 flag helps with speed.]) I obtained the following results:

    1)Computed Sparse; number non-zero: 3000001
    Time to run compSparseSigma w/ 1000000 mol: 0.674716 s

    2)Time to run compSparseSigma w/ 10 mol: 0.000010 s
    1  0  -  -  -  -  -  -  -  -  - 
    0  1  0  -  -  -  -  -  -  -  - 
    -  0  1  0  -  -  -  -  -  -  - 
    -  -  0  1  0  -  -  -  -  -  - 
    -  -  -  0  1  0  -  -  -  -  - 
    -  -  -  -  0  1  0  -  -  -  - 
    -  -  -  -  -  0  1  0  -  -  - 
    -  -  -  -  -  -  0  1  0  -  - 
    -  -  -  -  -  -  -  0  1  0  - 
    -  -  -  -  -  -  -  -  0  1  0 
    -  -  -  -  -  -  -  -  0  -  0

    3) Computed Sparse; number non-zero: 121
    Time to run compSparseSigma w/ 10 mol: 0.000043 s
    1  1  1  1  1  1  1  1  1  1  1 
    0  0  0  0  0  0  0  0  0  0  0 
    0  0  0  0  0  0  0  0  0  0  0 
    0  0  0  0  0  0  0  0  0  0  0 
    0  0  0  0  0  0  0  0  0  0  0 
    0  0  0  0  0  0  0  0  0  0  0 
    0  0  0  0  0  0  0  0  0  0  0 
    0  0  0  0  0  0  0  0  0  0  0 
    0  0  0  0  0  0  0  0  0  0  0 
    0  0  0  0  0  0  0  0  0  0  0 
    0  0  0  0  0  0  0  0  0  0  0 

These results confirm that the code is producing the expected results. Furthermore, it indicates that even large problems can be solved very quickly, ensuring that the implementation in the DAE library will be sufficiently quick.