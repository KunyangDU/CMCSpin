RBASIS3 = [[1/2,sqrt(3)/2,0.],[-1/2,sqrt(3)/2,0.],[0.,0.,1.]]
KBASIS2 = kbasis2(RBASIS3)
RBASIS2 = map(x -> x[1:2],RBASIS3[1:2])
KitaevBasis2 = [[1/2,sqrt(3)/2],[1/2,-sqrt(3)/2]]
FBZpoint = [[1/3,2/3],[2/3,1/3],[1/3,-1/3],[-1/3,-2/3],[-2/3,-1/3],[-1/3,1/3]]
BASIS2 = [[1.,0.],[0.,1.]]
Triapoint = [[-sqrt(3)/2,-1/2],[sqrt(3)/2,-1/2],[0.,1]]
JBASIS2 = Triapoint
KitaevBDpoiont = let 
    a,b,c = JBASIS2
    [(a+b)/2,(b+c)/2,(a+c)/2]
end
MFBZpoint = [[1,0],[1,1],[0,1],[-1,0],[-1,-1],[0,-1]]
dZZFBZpoint = map(x -> x/2,FBZpoint)

P = [
    1/sqrt(2) 1/sqrt(6) 1/sqrt(3);
    -1/sqrt(2) 1/sqrt(6) 1/sqrt(3);
    0 -2/sqrt(6) 1/sqrt(3)
]

function CubicHamiltonian(Latt::AbstractLattice;
    J,K,Γ,Γ′,Hv)
    As = [
    [J+K Γ′    Γ′;
    Γ′   J     Γ;
    Γ′   Γ     J],
    [J    Γ′    Γ;
    Γ′   J+K  Γ′;
    Γ    Γ′    J],
    [J    Γ     Γ′;
    Γ    J     Γ′;
    Γ′   Γ′    J+K]
    ]

    bonds = getxyzbonds(Latt)

    H = Hamiltonian()
    map(1:3) do i
        addIntr2!(H,Latt,"J$(i)",As[i],Tuple(bonds[i]))
    end
    addIntr1!(H,Latt,"H",- Hv, Tuple(1:length(Latt)))
    setgroup!(H,Latt.group)
    return H
end

function XCPH1(Lx::Int64 , Ly::Int64)
    Latt = XCPeriHoneycomb(Lx,Ly)
    groups = group(Latt,[1,],2)
    Latt = XCPeriHoneycomb(Lx,Ly,reverse_map(Tuple(vcat(groups...))))
    groups = map(g -> Tuple(Latt.sitemap[g]),groups)
    Latt.group = Tuple(groups)
    return Latt
end