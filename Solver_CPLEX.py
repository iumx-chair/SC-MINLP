import cplex 
import math
from wrapperc import ASLPY

def getSolver(asl : ASLPY) :
    n_var = asl.get_n_var()
    n_con = asl.get_n_con()

    varLB = asl.get_LUv()
    cplexvarLB = [varLB[i] for i in range(n_var)]
    varUB = asl.get_Uvx()
    cplexvarUB = [varUB[i] for i in range(n_var)]
    varInt = asl.get_intvar()

    objtype = asl.get_objtype()
    og = asl.get_Ograd(0)

    A_vals = asl.get_A_vals()
    A_rownos = asl.get_A_rownos()
    A_colstarts = asl.get_A_colstarts()

    UB = asl.get_Urhsx()
    LB = asl.get_LUrhs()

    cpxObjCoeff = [0]*n_var
    for i in range(n_var) :
        if og :
            j = og.contents.varno
            coef = og.contents.coef
            if math.fabs(coef) > 1e-30 :
                cpxObjCoeff[j] = coef
            og = og.contents.next
        else :
            continue
    
    cpx = cplex.Cplex()

    if (objtype[0] == b'\x00') :
        cpx.objective.set_sense(cpx.objective.sense.minimize)
    else :
        cpx.objective.set_sense(cpx.objective.sense.maximize)

    const = asl.get_objcons(0)
    cpx.objective.set_offset(const)

    varType = ['C']*n_var
    for j in range(n_var) :
        if varInt[j] :
            if cplexvarLB[j] == 0 and cplexvarUB[j] == 1 :
                varType[j] = 'B'
            else :
                varType[j] = 'I'

    varX = cpx.variables.add(obj=cpxObjCoeff,
                            lb=cplexvarLB, ub=cplexvarUB, 
                            types=varType, names=['x%d' % (i+1) for i in range(n_var)])

    my_rhs = [0]*n_con
    my_sense = ['E']*n_con

    for j in range(n_con) :
        my_rhs[j] = UB[j]
        if UB[j] >  cplex.infinity:
            my_sense[j] = 'G'
            my_rhs[j] = LB[j]
        elif LB[j] <  -cplex.infinity:
            my_sense[j] = 'L'


    cpx.linear_constraints.add(rhs = my_rhs, senses = my_sense)

    rows = []
    cols = []
    vals = []

    for j in range(n_var) :
        for i in range(A_colstarts[j],A_colstarts[j+1]) :
            if ([A_rownos[i]] and math.fabs(A_vals[i]) > 1e-30) :
                rows.append(A_rownos[i])
                cols.append(j)
                vals.append(A_vals[i])
    
    cpx.linear_constraints.set_coefficients(zip(rows, cols, vals))

    nlCons = asl.get_nlinearCons()
    nlCount = asl.get_nl_count()

    delNl = [i for i in range(n_con) if nlCons[i]]
    cpx.linear_constraints.delete(delNl)
    cpx.linear_constraints.set_names([(i, 'c%d' % (i+1)) for i in range(n_con-nlCount)])

    return cpx