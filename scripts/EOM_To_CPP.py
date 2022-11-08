import sys
import itertools

n_fields = 8

#internal
filename_i="EOMs_I.txt"
with open(filename_i,"r") as infile:
    variables = ["mu","mu1","Gx", "Gy", "ax","ay", "x","y", "z", "nperiodsx", "nperiodsy", "phasex", "phasey", "B", "c1"]

    types = ["const double" for _ in range(len(variables))]

    field = "const double f[][10]"

    var_comb = ", ".join((" ".join((typ,var)) for typ,var in zip(types,variables)))

    lines = [line.strip() for line in infile.readlines() if len(line.strip()) > 0]

    j = 0
    with open("output/HeaderEOMs_I.h","w") as header:
        header.write("#include <petscsys.h>\n") 
        header.write("#include <petscmath.h>\n") 
        for eom in range(n_fields):
            with open(f"output/EOMs_I_{eom}.c", "w") as output:

                output.write("#include <petscsys.h>\n") 
                output.write("#include <petscmath.h>\n") 

                output.write("#include \"math.h\"\n")
                output.write(f" double eom_{eom}("
                        f"{field}, {var_comb}){{\n return {lines[j]};\n}}\n")
                header.write(f" double eom_{eom}("
                        f"{field},{var_comb});\n")
                j += 1 



for bdy in [2]:
    for end in [0,1]:
        filename_bdy = f"EOMs_{bdy}_{end}.txt"
        with open(filename_bdy,"r") as infile:
                variables = ["mu","mu1","Gx", "Gy", "ax","ay", "x","y", "z", "nperiodsx", "nperiodsy", "phasex", "phasey", "B", "c1"]
        
                types = ["const double" for _ in range(len(variables))]
        
                field = "const double f[][10]"
        
                var_comb = ", ".join((" ".join((typ,var)) for typ,var in zip(types,variables)))
        
                lines = [line.strip() for line in infile.readlines() if len(line.strip()) > 0]
        
                j = 0
                with open(f"output/HeaderEOMs_{bdy}_{end}.h","w") as header:
                    header.write("#include <petscsys.h>\n") 
                    header.write("#include <petscmath.h>\n") 
                    for eom in range(n_fields):
                        with open(f"output/EOMs_{bdy}_{end}_{eom}.c", "w") as output:

                            output.write("#include <petscsys.h>\n") 
                            output.write("#include <petscmath.h>\n") 

                            output.write("#include \"math.h\"\n")
                            output.write(f" double eom_{bdy}_{end}_{eom}("
                                    f"{field}, {var_comb}){{\n return {lines[j]};\n}}\n")
                            header.write(f" double eom_{bdy}_{end}_{eom}("
                                    f"{field},{var_comb});\n")
                            j += 1 
