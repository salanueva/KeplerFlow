using Documenter, KeplerFlow

makedocs(modules=[KeplerFlow],
        doctest=true)
 
deploydocs(deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/salanueva/KeplerFlow.jl.git",
    julia  = "1.0",
    osname = "linux")