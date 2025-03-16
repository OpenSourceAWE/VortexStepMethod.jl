using Pkg

try
    Pkg.instantiate()
catch
    try
        Pkg.resolve()
    catch
        Pkg.update()
    end
end