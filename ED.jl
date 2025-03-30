# Exact diagonalization of SU(4) Heisengerg model

using LinearAlgebra
using SparseArrays
using KrylovKit

L   = 4
dim = 6
global jobID = round(Int,time()*1e4)

## Define all the generators

H1 = sparse([0 0 0 0 0 0;
             0 1/2 0 0 0 0; 
             0 0 1/2 0 0 0;
             0 0 0 -1/2 0 0;
             0 0 0 0 -1/2 0;
             0 0 0 0 0 0]) # H1 matrix

H2 = sparse([2/3 0 0 0 0 0;
             0 -1/3 0 0 0 0; 
             0 0  1/3 0 0 0;
             0 0 0 -1/3 0 0;
             0 0 0 0 1/3 0;
             0 0 0 0 0 -2/3]) # H2 matrix


H3 = sparse([1 0 0 0 0 0;
             0 1 0 0 0 0;
             0 0 -1 0 0 0; 
             0 0 0 1 0 0;
             0 0 0 0 -1 0;
             0 0 0 0 0 -1]) # H3 matrix

Ep1 = sparse([0 0 0 0 0 0;
              0 0 0 1 0 0; 
              0 0 0 0 1 0;
              0 0 0 0 0 0;
              0 0 0 0 0 0;
              0 0 0 0 0 0]) # Ep1 matrix


Ep2 = sparse([0 0 0 -1 0 0;
              0 0 0 0 0 0; 
              0 0 0 0 0 1;
              0 0 0 0 0 0;
              0 0 0 0 0 0;
              0 0 0 0 0 0]) # Ep2 matrix

Ep3 = sparse([0 1 0 0 0 0;
              0 0 0 0 0 0; 
              0 0 0 0 0 0;
              0 0 0 0 0 0;
              0 0 0 0 0 1;
              0 0 0 0 0 0]) # Ep3 matrix


Ep4 = sparse([0 0 0 0 -1 0;
              0 0 0 0 0 -1; 
              0 0 0 0 0 0;
              0 0 0 0 0 0;
              0 0 0 0 0 0;
              0 0 0 0 0 0]) # Ep4 matrix              

Ep5 = sparse([0 0 1 0 0 0;
              0 0 0 0 0 0; 
              0 0 0 0 0 0;
              0 0 0 0 0 -1;
              0 0 0 0 0 0;
              0 0 0 0 0 0]) # Ep5 matrix

Ep6 = sparse([0 0 0 0 0 0;
              0 0 1 0 0 0; 
              0 0 0 0 0 0;
              0 0 0 0 1 0;
              0 0 0 0 0 0;
              0 0 0 0 0 0]) # Ep6 matrix

Em1 = Ep1'
Em2 = Ep2'
Em3 = Ep3'
Em4 = Ep4'
Em5 = Ep5'
Em6 = Ep6'


# Define the Hamiltonian
J1 = 1.0;
J2 = 1.0;
J3 = 1.0;
J4 = 1.0;
J5 = 1.0;
J6 = 1.0;
J7 = 1.0;
J8 = 1.0;
J9 = 1.0;

Coupling = J1*kron(H1, H1)+J2*kron(H2, H2)+J3*kron(H3, H3)+
           0.5*J4*kron(Ep1, Em1)+0.5*J5*kron(Ep2, Em2)+0.5*J6*kron(Ep3, Em3)+
           0.5*J7*kron(Ep4, Em4)+0.5*J8*kron(Ep5, Em5)+0.5*J9*kron(Ep6, Em6)+
           0.5*J4*kron(Em1, Ep1)+0.5*J5*kron(Em2, Ep2)+0.5*J6*kron(Em3, Ep3)+
           0.5*J7*kron(Em4, Ep4)+0.5*J8*kron(Em5, Ep5)+0.5*J9*kron(Em6, Ep6)

let
    Hamiltonian = spzeros(dim^L, dim^L)

    for i in 1:L-1
        Hamiltonian += kron(sparse(I,dim^(i-1),dim^(i-1)),kron(Coupling,sparse(I,dim^(L-i-1),dim^(L-i-1))))
    end


    ## Diagonalization of the Hamiltonian
    n_lowest = 5
    which = :SR
    Abstract_Array(x) = Hamiltonian*x
    dimension = size(Hamiltonian)[1]
    elapsed_time = @elapsed begin
        (vals, vecs, res) = eigsolve(Abstract_Array, dimension, n_lowest, which, tol = 1e-6)
    end
    # Save the eigenvalues to a text file
    open("result.txt", "a") do file
        println(file, "\n\n=======================")
        println(file, "Exact Diagonalization of SU(4) Heisengerg model")
        println(file, "=======================")
        println(file, "jobID: ", jobID)
        println(file, "\nHamiltonian Information:")
        if ishermitian(Hamiltonian)
            println(file, "\nThe Hamiltonian is Hermitian")
        else
            println(file, "\nThe Hamiltonian is NOT Hermitian.")
        end
        println(file, "Number of sites: ", L)
        println(file, "Dimension of the Hamiltonian: ", dimension)
        println(file, "Number of lowest eigenvalues computed: ", n_lowest)
        println(file, "\nEigenvalues:")
        for val in vals
            println(file, val)
        end
        println(file, "\ntimes taken to diagonalize the Hamiltonian: ", elapsed_time)
        println(file, " =============DONE=============")
    end
    println("DONE")
    
    #@show vals
end