#!/usr/bin/env python3
# usage: dynamically generate C code according to the number of variables
# to keep
import sys

if len(sys.argv) < 2:
    sys.exit('usage: <script> <number of vars to keep>')


k = int(sys.argv[1])
decl_f = 'gc_decl_lsys.def'
cpy_f = 'gc_copy_lsys.def'
check_f = 'gc_check_lsys.def'
solve_f = 'gc_solve_lsys.def'
extract_f = 'gc_extract_sol.def'
dep_f = 'gc_dep_lsys.def'
gauss_f = "gc_gauss.def"

with open(decl_f, "w") as fp:  # declare lsys
    for i in range(0, k+1):
        fp.write("register uint32_t lsys%d;\n" % (i))


with open(cpy_f, "w") as fp:  # copy lsys from gc data structure
    for i in range(0, k+1):
        fp.write("lsys%d = eterm_lsys(wdata, warp_tid, %d);\n" % (i, i))


with open(check_f, "w") as fp:  # check lsys
    fp.write("uint32_t rmask = ~0x0U;\n")

    for i in range(0, k):
        fp.write("{\n")
        fp.write("\tuint32_t tmp = lsys%d & rmask;\n" % (i))
        fp.write("\tuint32_t sf = (!tmp) ? 0x0U : ~0x0U;\n")
        fp.write("\tuint32_t piv = ctz(tmp);\n")
        fp.write("\trmask ^= (0x1U << piv) & sf;\n")
        fp.write("\tuint32_t mask = (lsys%d ^ (0x1U << piv) ) & sf;\n" % (i))
        fp.write("\tlsys%d ^= mask;\n" % (i))

        for j in range(i+1, k+1):
            fp.write("\t\tlsys%d ^= mask & (((lsys%d >> piv) & 0x1U) ? ~0x0U : 0x0U);\n" % (j, j))
        fp.write("}\n")

    fp.write("solvable = !(lsys%d & rmask);\n" % (k))


with open(extract_f, "w") as fp:
    for i in range(0, k):
        fp.write("sol |= ( (lsys%d >> ctz(lsys%d)) & 0x1U ) << %d;\n" % (k, i, i))


with open(gauss_f, "w") as fp:
    for i in range(0, k+1):
        fp.write("if( lsys%d >> %d ) {\n" % (i, i))
        fp.write("\tuint32_t piv = ctz(lsys%d >> %d) + %d;\n" % (i, i, i))
        fp.write("\tuint32_t mask = lsys%d ^ (0x1U << piv);\n" % (i))
        fp.write("\tlsys%d ^= mask;\n" % (i))

        for j in range(0, i+1):
            fp.write("\t{\n")
            fp.write("\t\tuint32_t tmp = ((lsys%d >> %d) ^ (lsys%d >> piv)) & 0x1U;\n" % (j, i, j))
            fp.write("\t\tlsys%d ^= (tmp << %d) | (tmp << piv);\n" % (j, i))
            fp.write("\t}\n")

        for j in range(i+1, k+1):
            fp.write("\tif((lsys%d >> piv) & 0x1U) {\n" % (j))
            fp.write("\t\tlsys%d ^= mask;\n" % (j))
            fp.write("\t}\n")
            fp.write("\t{\n")
            fp.write("\t\tuint32_t tmp = ((lsys%d >> %d) ^ (lsys%d >> piv)) & 0x1U;\n" % (j, i, j))
            fp.write("\t\tlsys%d ^= (tmp << %d) | (tmp << piv);\n" % (j, i))
            fp.write("\t}\n")

        fp.write("}\n")

    fp.write("solvable = !((lsys%d >> %d) & 0x1U);\n" % (k, k))

with open(dep_f, "w") as fp:
    for i in  range(0, k):
        fp.write("if(!lsys%d) {\n" % (i))
        fp.write("\tdep = 0x1U;\n")
        fp.write("}\n")
