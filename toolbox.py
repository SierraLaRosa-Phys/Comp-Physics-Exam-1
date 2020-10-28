import numpy as np

def Sieve(n):

    # We need an array of n true entries. entries 0 and 1 will be marked false as they are not primes
    arr = [ True for i in range(n+1) ]
    arr[0] = False
    arr[1] = False

    p = 2
    q = 0

    primes = []
    while ( q <= n ):

        #  If arr[p] is not changed we say that it is prime.
        if arr[p]:
            primes.append(p)

            # now we will update the array arr for multiples of p
            for inst in range( p*2, n+1, p):
                arr[inst] = False

        for jnst in range(n+1):
            if arr[jnst] and jnst > p:
                q = jnst
                p = q
                break
            elif jnst > p:
                q = n + 1


    return primes

#---define functions for iv problems

#define the differential equation for Euler's method
def iv_diff_eq(xi, b):
    px_i, py_i, vx_i, vy_i = xi[0], xi[1], xi[2], xi[3]
    ax = - b*vx_i
    ay = -1 - b*vy_i
    return np.array([vx_i, vy_i, ax, ay], float)

#define Euler's method
def euler(xi, dt, b):
    xi += dt*iv_diff_eq(xi, b)
    return xi

#---define definitions for boundary condition solvers---

#first to define the gaussian elimination
def exchange_rows(A, I, J):
    for inst in range(len(a)):
        temp = A[I][inst]
        A[I][inst] = A[J][inst]
        A[J][inst] = temp
    return A

# Define a function to turn a matrix into a Upper triangle for Gaussian eliminations
def gaussian_UT(matrix, vector):
    matrix_new = np.column_stack((matrix,vector))
    for inst in range(len(matrix_new)):
        # Switch rows if one of the diaganal elements in zero
        if matrix_new[inst][inst] == 0:
            exchange_rows(matrix_new,inst,inst+1)
        else:
            diag = matrix_new[inst][inst]
            for jnst in range(len(matrix_new) + 1):
                matrix_new[inst][jnst] = float(matrix_new[inst][jnst])/float(diag)

            knst = inst + 1
            while knst < len(matrix_new):
                temp = matrix_new[knst][inst]
                for jnst in range(len(matrix_new)+1):
                    matrix_new[knst][jnst] = matrix_new[knst][jnst] - temp*matrix_new[inst][jnst]
                knst +=1
    return matrix_new

#For boundaries I am just using an array instead of a dict
def back_sub_bc(matrix, boundaries):
    row = np.zeros(len(matrix))

    for jnst in boundaries.keys():
        row[jnst] = boundaries[jnst]

    inst = len(row) - 2
    while inst >= 0:
        holder = 0.0
        jnst = inst + 1
        while jnst < len(row):
            holder += matrix[inst][jnst]*row[jnst]
            jnst += 1
        row[inst] = matrix[inst][len(row)] - holder
        inst -= 1
    return row

def Elimination_bc(Matrix, Vector, Boundaries):
    New_Mat = gaussian_UT(Matrix, Vector)
    result = back_sub_bc(New_Mat, Boundaries)
    return result

def Jacobi_bc(Matrix, Vector, Boundaries, epsilon=1e-10, max_iter=5000):
    Diagonal = np.diag(np.diag(Matrix))
    Lower_Upper = Matrix - Diagonal
    Result = np.zeros(len(Vector))

    Inverse_Diagonal = np.diag(1/np.diag(Diagonal))

    for inst in range(max_iter):
        temp = np.dot(Inverse_Diagonal, Vector - np.dot(Lower_Upper, Result))

        for jnst in Boundaries.keys():
            temp[jnst] = Boundaries[jnst]

        if np.linalg.norm(temp - Result) < epsilon:
            return temp

        Result = temp
    return Result

def bc_diff_eq(coeff, rhs, a, b, Boundaries, n, solver):
    delta = (b-a)/n
    N = n + 1
    Matrix =  np.zeros((N,N))
    Delta = delta*delta

    if len(coeff) < 3:
        raise Exception("Not enough coefficient functions given.")

    x = lambda i: a+i*delta

    inst = 0
    Matrix[inst][inst] = -2.0*coeff[2](x(inst))/Delta
    Matrix[inst][inst + 1] = 1.0*coeff[2](x(inst))/Delta
    for inst in range(1, N-1):
        Matrix[inst][inst - 1] = 1.0*coeff[2](x(inst))/Delta
        Matrix[inst][inst] = -2.0*coeff[2](x(inst))/Delta
        Matrix[inst][inst + 1] = 1.0*coeff[2](x(inst))/Delta
    Matrix[N-1][N-1]=-2.0*coeff[2](x(inst))/Delta
    Matrix[N-1][N-2]=1.0*coeff[2](x(inst))/Delta

    for inst in range(0, N-1):
        Matrix[inst][inst] = Matrix[inst][inst] - 1.0*coeff[1](x(inst))/delta
        Matrix[inst][inst+1] = Matrix[inst][inst+1] + 1.0*coeff[1](x(inst))/delta
    Matrix[N-1][N-1] = Matrix[N-1][N-1] - 1.0*coeff[1](x(inst))/delta

    for inst in range(0, N):
        Matrix[inst][inst] = Matrix[inst][inst] + coeff[0](x(inst))

    equals = np.zeros(N)

    for inst in range(0,N):
        equals[inst] = rhs(x(inst))

    return solver(Matrix, equals, Boundaries)

def equation(l=2, n=500, solver=Elimination_bc ):
    p = lambda x: 1.0 - x*x
    q = lambda x: 0
    r = lambda x: l*(l-1)

    bc={n:1.0}

    if solver == Elimination_bc:
        a=-1.0
        b=1
        Delta = (b-a)/(n+1)
    else:
        a=-1
        b=1
        Delta = (b-a)/(n+1)

    x = np.arange(a,b,Delta)

    coeff = [r,q,p]

    rhs = lambda x: 0
    """
    def rhs(x):
        if x == 0:
            return -1/(Delta**3)
        else:
            return 0
    """

    return x,bc_diff_eq(coeff, rhs, a, b, bc, n, solver)

if __name__ == '__main__':
    n = 100
    print (f"The test is {n} integers and the prime numbers are\n{Sieve(n)}")


#    x = legendre(solver=np.linalg)
    x,y = equation(l=5, solver=Elimination_bc)
    #x,z = equation(solver=Jacobi_bc)
    import matplotlib.pyplot as mp

    x *= 2.0e-9

#    mp.plot(x)
    mp.plot(x[:len(y)],y)
    #mp.plot(x[:len(z)],z)

    mp.show()
