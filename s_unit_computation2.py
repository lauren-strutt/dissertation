K.<x > = CyclotomicField (8)
G = K.galois_group () #The Galois group used in the Log - embedding
alpha =-18*x^3-15*x^2-12*x-6
S_generators=[1+ x + x ^(-1), x ^2+ x -1,x ^3 + x ^2 +1,2 - x, 2 - x ^3, 2+ x ^3,2+ x] # s-unit group generators
def Log(f,sgens): #logarithmic embedding function
LogEmbedding=[float ( log (abs ( f * conjugate(f) ) ) ), float ( log (abs ( G [3](f) * conjugate ( G [3](f) ) ) ) )] #units under log embedding
for i in range(1,len(sgens)):
    LogEmbedding.append( float ( log ( sgens[i] . norm () ^( -( K . valuation ( sgens[i] ) (f) ) ) ) )) #adding places
return LogEmbedding

print(f"log(alpha) is {Log(alpha,S_generators)}")
def Loglattice(units):
loglattice=[]
for s in units:
    loglattice.append(Log(s,S_generators))
loglattice.append([1 ,1 ,1 ,1  ,1 ,1 ,1 ,1]) #after we have created our s-unit lattice we add a row of 1s to make it square so that we can perform later matrix calculations
return loglattice
M=matrix(Loglattice([1+ x + x ^(-1),x ^2+ x -1, x ^3 + x ^2 +1, 2 - x, 2 - x ^3, 2+ x ^3,2+ x])) # we turn our log s-unit lattice into a matrix
t= M .solve_left(vector(Log(alpha,S_generators))) #we calculate log(alpha)M^{-1} to calculate the coefficients of ti in Babai's algorithm
print(t)
v=alpha/(S_generators[1]*S_generators[2]*S_generators[6]) #round the vector t and use the values as powers of our prime ideals to divide by alpha
print(v)
print(alpha.norm())
print(v.norm()) #by computing norms we see that the vector has become shorter
