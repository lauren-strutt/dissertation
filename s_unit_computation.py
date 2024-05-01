K.<x > = CyclotomicField (16)
G = K.galois_group () #The Galois group used in the Log - embedding
OK= K . ring_of_integers ()
t=vector ([ 0, 0 , 0, 0 ])
def Log(f): #logarithmic embedding function
    LogEmbedding=[float ( log (abs ( f * conjugate(f) ) ) )] #units under log embedding
    for i in [3,1,2]:
        LogEmbedding.append( float ( log (abs ( G [int(i)]( f ) * conjugate ( G [int(i)]( f ) ) ) ) ))
    return vector(LogEmbedding)
def Loglattice(units):
    loglattice=[]
    for u in units:
        loglattice.append(Log(u))
    loglattice.append([1 ,1 ,1 ,1]) #after we have created our s-unit lattice we add a row of 1s to make it square so that we can perform later matrix calculations
    return loglattice

while(norm(t[:3])<2): # so that we divide by units and get a nontrivial example
    alpha = OK.random_element()
    latticematrix=matrix(Loglattice([1+x+x^(-1),1+x^3+x^(-3),1+x^5+x^(-5)]))
    t = latticematrix . solve_left (Log(alpha))
print(alpha)
print(t)
