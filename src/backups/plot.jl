function snapshot(ax::Axis,Latt::AbstractLattice;kwargs...)

    for bond in neighbor(Latt)
        x =  map(bond) do i
            coordinate(Latt,i)[1]
        end
        y =  map(bond) do i
            coordinate(Latt,i)[2]
        end
        lines!(ax,collect.((x,y))...;color = :grey,linewidth = 2)
    end

    for i in 1:length(Latt)
        co = coordinate(Latt,i)
        scatter!(ax,co;color = :grey,markersize = 14)
        text!(ax,(co .+ 0.1)...;text = "$(i)")
    end


end

function snapshot(Latt::AbstractLattice;kwargs...)
    Lx,Ly = size(Latt)
    fig = Figure()
    ax = Axis(fig[1,1],autolimitaspect = true;
    xticks = 0:Lx,yticks = 0:Ly)
    snapshot(ax,Latt)
    resize_to_layout!(fig)
    display(fig)
    return fig,ax
end
