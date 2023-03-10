10/03/23

Problems with the DiferentialEquations packages that do not let run the code of the main.jl fail

The lines to test are from 174 to 330.

The problem appears in the line 300 with this comment:

> ERROR: Failed to precompile DifferentialEquations [0c46a032-eb83-5123-abaf-570d42b7fbaa] to /home/melihan/.julia/compiled/v1.7/DifferentialEquations/jl_LFWZoY.
> Stacktrace:
>  [1] error(s::String)
>    @ Base ./error.jl:33
>  [2] compilecache(pkg::Base.PkgId, path::String, internal_stderr::IO, internal_stdout::IO, ignore_loaded_modules::Bool)
>    @ Base ./loading.jl:1466
>  [3] compilecache(pkg::Base.PkgId, path::String)
>    @ Base ./loading.jl:1410
>  [4] _require(pkg::Base.PkgId)
>    @ Base ./loading.jl:1120
>  [5] require(uuidkey::Base.PkgId)
>    @ Base ./loading.jl:1013
>  [6] require(into::Module, mod::Symbol)
>    @ Base ./loading.jl:997
