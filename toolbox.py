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
                print(p)
                q = jnst
                p = q
                break
            elif jnst > p:
                q = n + 1


    return primes

if __name__ == '__main__':
    n = 100
    print (f"The test is {n} integers and the prime numbers are\n{Sieve(n)}")
