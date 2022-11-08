import sys
import itertools

n_fields = 9
nderivs = 10
n_eom=n_fields

#internal
filename_i="Coefs_I.txt"
with open(filename_i,"r") as infile:
    variables = ["mu","Q","Gx", "Gy", "ax","ay", "x","y", "z", "nperiodsx", "nperiodsy", "phasex", "phasey", "B", "c1"]

    types = ["const double" for _ in range(len(variables))]

    field = "const double f[][10]"

    var_comb = ", ".join((" ".join((typ,var)) for typ,var in zip(types,variables)))

    lines = [line.strip() for line in infile.readlines() if len(line.strip()) > 0]

    j = 0
    with open("output_inlined/HeaderCoefs_I.h","w") as header:
        header.write("#include <petscsys.h>\n")
        header.write("#include <petscmath.h>\n")

        for nfield in range(n_fields):
            for deriv in range(nderivs):
                for eom in range(n_eom):
                    # Inline the trivial ones
                    if len(lines[j]) <= 2:
                        header.write(f"inline static double Coefs_I_{nfield}_{deriv}_{eom}("
                                f"{field},{var_comb}){{\n return {lines[j]};\n}}\n")
                    else:

                        with open(f"output_inlined/Coefs_I_{nfield}_{deriv}_{eom}.c", "w") as output:
                            output.write("#include <petscsys.h>\n")
                            output.write("#include <petscmath.h>\n")
                            output.write("#include \"math.h\"\n")
                            output.write(f"double Coefs_I_{nfield}_{deriv}_{eom}("
                                    f"{field}, {var_comb}){{\n return {lines[j]};\n}}\n")
                            header.write(f"double Coefs_I_{nfield}_{deriv}_{eom}("
                                    f"{field},{var_comb});\n")
                    j += 1 



for bdy in [2]:
    for end in [0,1]:
        filename_bdy = f"Coefs_{bdy}_{end}.txt"
        with open(filename_bdy,"r") as infile:
                variables = ["mu","Q","Gx", "Gy", "ax","ay", "x","y", "z", "nperiodsx", "nperiodsy", "phasex", "phasey", "B", "c1"]
        
                types = ["const double" for _ in range(len(variables))]
        
                field = "const double f[][10]"
        
                var_comb = ", ".join((" ".join((typ,var)) for typ,var in zip(types,variables)))
        
                lines = [line.strip() for line in infile.readlines() if len(line.strip()) > 0]
        
                j = 0
                with open(f"output_inlined/HeaderCoefs_{bdy}_{end}.h","w") as header:
                    header.write("#include <petscsys.h>\n")
                    header.write("#include <petscmath.h>\n")
                    for nfield in range(n_fields):
                        for deriv in range(nderivs):
                            for eom in range(n_eom):
                                if len(lines[j]) <=2:
                                    header.write(f"inline static double Coefs_{bdy}_{end}_{nfield}_{deriv}_{eom}("
                                    f"{field},{var_comb}){{\n return {lines[j]};\n}}\n")
                                else:
                                    with open(f"output_inlined/Coefs_{bdy}_{end}_{nfield}_{deriv}_{eom}.c", "w") as output:
                                        output.write("#include <petscsys.h>\n")
                                        output.write("#include <petscmath.h>\n")
                                        output.write("#include \"math.h\"\n")
                                        output.write(f"double Coefs_{bdy}_{end}_{nfield}_{deriv}_{eom}("
                                            f"{field}, {var_comb}){{\n return {lines[j]};\n}}\n")
                                        header.write(f"double Coefs_{bdy}_{end}_{nfield}_{deriv}_{eom}("
                                            f"{field},{var_comb});\n")
                                j += 1 
