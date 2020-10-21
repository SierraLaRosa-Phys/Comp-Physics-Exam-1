from time import time
from toolbox import Sieve

print("-"*40, "Part A", "-"*40)
print("Pleas see toolbox.py for the algorithm of this question. It is called Sieve")

print("-"*40, "Part B", "-"*40)
print( "The primes in the range of [0,30] are:")
print( "2,3,5,7,11,13,17,19,23,29" )
print( f"My algorithm for the Sieve of Eratosthenes gives:\n{Sieve(30)}" )

print("-"*40, "Part C", "-"*40)
n = 100000000
print( f"The first {n} primes are:\n{Sieve(n)}" )

