RBASIS3 = [[1.0,0.0,0.0],[1/2,sqrt(3)/2,0.0],[0.,0.,1.]]
KBASIS2 = kbasis2(RBASIS3)
RBASIS2 = map(x -> x[1:2],RBASIS3[1:2])
FBZpoint = [[0,0],[0,1/2],[1/2,1/2],[1/2,0]]

