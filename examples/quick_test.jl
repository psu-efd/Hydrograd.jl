#println("Current working directory: ", pwd())
#println("LOAD_PATH: ", LOAD_PATH)

#push!(LOAD_PATH, dirname(dirname(@__FILE__)))
import Hydrograd

# Now try to use print_banner
Hydrograd.print_banner()