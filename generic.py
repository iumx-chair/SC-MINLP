#!/usr/bin/env python3

import json
import cplex

from sympy import *
from math import *

from sys import argv
import os


class NLProblem():
    PC_EPSCHCKLOOP = 1e-9
    PC_RELAXCUTLIM = 3000
    PC_PARAM_EPSILON = 1e-6
    PC_THRSHLD = 1e-6
    CPLEX_MIPGAP = 1e-7  # First 1e-7
    CPLEX_ABSMIPGAP = 1e-7

    PC_REF_DISTPOINTS = 0  # First 1E-7
    PC_REF_THRSHLD = 1E-5

    GLOBAL_GAPTOL = 1E-4

    def __init__(self):
        self.poolconst = []
        self.poolsense = []
        self.poolrhs = []
        self.local_lb = []
        self.local_ub = []

        self.prevxbar = {}
        self.relaxcutcount = 0
        self.cutscount = 0

        # Refinement
        self.iterations = 0
        self.maxIterations = -1
        self.totalBreakpointsAdded = 0

        # Options
        self.display = -1
        self.relax = -1
        self.reformulationType = -1
        self.timelimit = -1
        self.totaltime = 0

    def LoadBreakpoints(self, bpcalc):
        self.varxname = bpcalc.varxname
        self.varyname = bpcalc.varyname

        self.varxid = bpcalc.varxid
        self.varyid = bpcalc.varyid
        self.varzid = bpcalc.varzid

        self.varx_lb = bpcalc.varxlb
        self.varx_ub = bpcalc.varxub

        self.breakpoints = bpcalc.setBreakpoints
        self.estimation = bpcalc.setEstimation

        self.functions = bpcalc.functions1
        self.funcdiff2 = bpcalc.functions2
        self.ineq = bpcalc.ineq
        self.funcSense = [('G' if i == 1 else 'L') for i in self.ineq]

        self.nFunctions = len(self.functions)

        # Change from json original to newOne
        self.Nitem = len(self.functions)
        self.Nb = [len(self.breakpoints[i]) for i in range(self.Nitem)]
        # self.Nb = []
        # for i in range(self.Nitem ) :
        #     self.Nb.append(len(self.breakpoints[i]))
        self.Lparam = self.breakpoints
        self.Concave = self.estimation

        self.CompileFunctions()

    def LoadFileJSON(self, FileName="dataNew.json"):
        f = open(FileName)
        jsonList = json.load(f)

        self.varxname = jsonList["varx_name"]
        self.varyname = jsonList["vary_name"]

        self.varxid = jsonList["varx_id"]
        self.varyid = jsonList["vary_id"]
        self.varzid = jsonList["varz_id"]

        self.varx_lb = jsonList["varx_lb"]
        self.varx_ub = jsonList["varx_ub"]

        self.breakpoints = jsonList["breakpoints"]
        self.estimation = jsonList["estimation"]

        self.functions = jsonList["functions"]
        self.funcdiff2 = jsonList["funcdiff2"]

        if "ineq" in jsonList:
            self.ineq = jsonList["ineq"]
        else:
            self.ineq = [2]*len(self.functions)
            print(
                "Warning: Inequalities sense are missing on .json file. Assuming >= as default")

        self.funcSense = [('G' if i == 1 else 'L') for i in self.ineq]

        self.nFunctions = len(self.functions)

        # Change from json original to newOne
        self.Nitem = len(self.functions)
        self.Nb = [len(self.breakpoints[i]) for i in range(self.Nitem)]

        # for i in range(self.Nitem ) :
        #     self.Nb.append(len(self.breakpoints[i]))

        self.Lparam = self.breakpoints
        self.Concave = self.estimation

        self.CompileFunctions()

    def CompileFunctions(self):
        from types import FunctionType

        self.NLfunc = [FunctionType(compile(f"def myNLFunction{str(i)} (x):\n\ty = {self.functions[i]} \n\treturn y", '<string>', 'exec').co_consts[0], globals())
                       for i in range(self.nFunctions)]

        self.NLfuncDiff = [FunctionType(compile(f"def myNLFunctionDiff{str(i)} (x):\n\ty = {self.funcdiff2[i]} \n\treturn -y", '<string>', 'exec').co_consts[0], globals())
                           for i in range(self.nFunctions)]

    def LoadOptions(self, options):
        # KW((char *)"display",     L_val, &optionsreturn[0],  (char *)"Frequency of information display; default 1."),
        # KW((char *)"pc_refmode",  L_val, &optionsreturn[1],  (char *)"PC reformulation mode: 1 Incremental, 2 Multiple Choice; default = 1."),
        # KW((char *)"relax",       L_val, &optionsreturn[2],  (char *)"Ignore integrality; default 0."),
        # KW((char *)"time"  ,      L_val, &optionsreturn[3],  (char *)"Time limit in seconds; default = 1e75."),
        # KW((char *)"timing",      L_val, &optionsreturn[4],  (char *)"display timings for the run."),

        # display
        if (options[0] >= 0):
            self.display = options[0]
            pass

        # pc_refmode
        if (options[1] >= 0):
            self.reformulationType = options[1]

        # relax
        if (options[2] >= 0):
            self.relax = options[2]

        # time
        if (options[3] >= 0):
            self.timelimit = options[3]

        # timing

    def CreateProblemSolv(self, solver : cplex, mode=1):
        'Chose the reformulation model.\n option = 1 for Incremental \n option = 2 for Multiple Choice.'

        if (self.reformulationType >= 0):
            mode = self.reformulationType

        self.p_mode = mode

        cpx = cplex.Cplex(solver)
        self.cpx = cpx
        self.Problem_Sense = self.cpx.objective.sense[self.cpx.objective.get_sense()]

        self.varXOri = self.varxid
        self.varYOri = self.varyid
        self.varZOri = self.varzid

        self.originalConsNum = cpx.linear_constraints.get_num()

        self.varX = [i-1 for i in self.varxid]
        self.varY = [i-1 for i in self.varyid]

        # Map original indexes from ampl
        self.originalVarNum = cpx.variables.get_num()

        # self.varnamesOri = []
        # for i in range(self.originalVarNum):
        #     self.varnamesOri.append("x" + str(i+1))

        # self.newVarIndex = cpx.variables.get_indices(self.varnamesOri)
        self.newVarIndex = range(self.originalVarNum)

        self.cpxCopy = cplex.Cplex(cpx)
        self.cpxCopy.parameters.paramdisplay.set(0)
        self.cpxCopy.parameters.simplex.tolerances.optimality.set(1e-9)
        self.cpxCopy.parameters.simplex.tolerances.feasibility.set(1e-9)
        self.cpxCopy.parameters.read.datacheck.set(0)
        self.cpxCopy.parameters.mip.tolerances.mipgap.set(self.CPLEX_MIPGAP)
        self.cpxCopy.parameters.mip.tolerances.absmipgap.set(
            self.CPLEX_ABSMIPGAP)

        if self.iterations == 0 and self.Problem_Sense == 'minimize':
            self.objValGlobalLB = -cplex.infinity
            self.objValGlobalUB = cplex.infinity
            self.varValGlobalUB = [0]*self.cpxCopy.variables.get_num()

        elif self.iterations == 0 and self.Problem_Sense == 'maximize':
            self.objValGlobalLB = cplex.infinity
            self.objValGlobalUB = -cplex.infinity
            self.varValGlobalUB = [0]*self.cpxCopy.variables.get_num()

        if mode == 1:
            print("\nStarting Incremental Model...")
            self.CreateProblem_Incremental()
        if mode == 2:
            print("\nStarting Multiple Choice Model...")
            self.CreateProblem_MultipleChoiceTest()
        if mode == 3:
            print("\nStarting Multiple Choice Model - Incremental...")
            self.CreateProblem_MultipleChoiceInc()

        self.listcplexvarnum = list(range(cpx.variables.get_num()))
        self.listcplexvartypes = cpx.variables.get_types(self.listcplexvarnum)
        self.listcplexvarInt = [
            i for i in self.listcplexvarnum if self.listcplexvartypes[i] != 'C']

        self.prevSolVal = [-1] * cpx.variables.get_num()

    def CreateProblemMPS(self, mpsfile="problem.mps", mode=1):
        'Chose the reformulation model.\n option = 1 for Incremental \n option = 2 for Multiple Choice.'

        if (self.reformulationType >= 0):
            mode = self.reformulationType

        self.p_mode = mode

        cpx = cplex.Cplex(mpsfile)
        self.cpx = cpx
        self.Problem_Sense = self.cpx.objective.sense[self.cpx.objective.get_sense(
        )]

        self.varXOri = self.varxid
        self.varYOri = self.varyid
        self.varZOri = self.varzid

        self.originalConsNum = cpx.linear_constraints.get_num()

        self.varX = cpx.variables.get_indices(self.varxname)
        self.varY = cpx.variables.get_indices(self.varyname)

        # Map original indexes from ampl
        self.originalVarNum = cpx.variables.get_num()

        self.varnamesOri = []
        for i in range(self.originalVarNum):
            self.varnamesOri.append("x" + str(i+1))

        self.newVarIndex = cpx.variables.get_indices(self.varnamesOri)

        self.cpxCopy = cplex.Cplex(cpx)
        self.cpxCopy.parameters.paramdisplay.set(0)
        self.cpxCopy.parameters.simplex.tolerances.optimality.set(1e-9)
        self.cpxCopy.parameters.simplex.tolerances.feasibility.set(1e-9)
        self.cpxCopy.parameters.read.datacheck.set(0)
        self.cpxCopy.parameters.mip.tolerances.mipgap.set(self.CPLEX_MIPGAP)
        self.cpxCopy.parameters.mip.tolerances.absmipgap.set(
            self.CPLEX_ABSMIPGAP)

        if self.iterations == 0 and self.Problem_Sense == 'minimize':
            self.objValGlobalLB = -cplex.infinity
            self.objValGlobalUB = cplex.infinity
            self.varValGlobalUB = [0]*self.cpxCopy.variables.get_num()

        elif self.iterations == 0 and self.Problem_Sense == 'maximize':
            self.objValGlobalLB = cplex.infinity
            self.objValGlobalUB = -cplex.infinity
            self.varValGlobalUB = [0]*self.cpxCopy.variables.get_num()

        if mode == 1:
            print("\nStarting Incremental Model...")
            self.CreateProblem_Incremental()
        if mode == 2:
            print("\nStarting Multiple Choice Model...")
            self.CreateProblem_MultipleChoiceTest()

        self.listcplexvarnum = list(range(cpx.variables.get_num()))
        self.listcplexvartypes = cpx.variables.get_types(self.listcplexvarnum)
        self.listcplexvarInt = [
            i for i in self.listcplexvarnum if self.listcplexvartypes[i] != 'C']

        self.prevSolVal = [-1] * cpx.variables.get_num()

    def CreateProblem(self, mode=1):
        'Chose the reformulation model.\n option = 1 for Incremental \n option = 2 for Multiple Choice.'
        self.p_mode = mode

        Nitem = self.Nitem
        Umax = self.Umax
        Cap = self.Cap
        weight = self.weight
        UBvar = self.UBvar
        Nb = self.Nb
        Lparam = self.Lparam
        Concave = self.Concave

        setItem = list(range(Nitem))

        cpx = cplex.Cplex()
        self.cpx = cpx

        cpx.objective.set_sense(cpx.objective.sense.minimize)

        varX = cpx.variables.add(obj=[0] * Nitem,
                                 lb=[0] * Nitem, ub=UBvar,
                                 types=['C'] * Nitem,
                                 names=['x(%d)' % (j+1) for j in setItem])

        varY = cpx.variables.add(obj=[-1] * Nitem,
                                 lb=[-cplex.infinity] * Nitem, ub=[cplex.infinity] * Nitem,
                                 types=['C'] * Nitem,
                                 names=['y(%d)' % (j+1) for j in setItem])

        cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(varX, weight)],
            senses=['L'],
            rhs=[Cap])

        self.varX = varX
        self.varY = varY

        if mode == 1:
            self.CreateProblem_Incremental()
        if mode == 2:
            self.CreateProblem_MultipleChoice()
        if mode == 3:
            self.CreateProblem_MultipleChoiceInc()

        self.listcplexvarnum = list(range(cpx.variables.get_num()))
        self.listcplexvartypes = cpx.variables.get_types(self.listcplexvarnum)

    # Incremental reformulation
    def CreateProblem_Incremental(self):
        Nitem = self.Nitem
        # Umax = self.Umax
        # Cap = self.Cap
        # weight = self.weight
        # UBvar = self.UBvar
        Nb = self.Nb
        Lparam = self.Lparam
        Concave = self.Concave
        Ineq = self.ineq

        cpx = self.cpx
        varX = self.varX
        varY = self.varY

        # Sets
        setVars = list(range(Nitem))
        setS = [list(range(Nb[j]-1)) for j in setVars]

        setSConcave = [
            (j, s)
            for j in setVars for s in setS[j]
            if int(Concave[j][s]) == 1
        ]

        setSConvex = [
            (j, s)
            for j in setVars for s in setS[j]
            if int(Concave[j][s]) == -1
        ]

        setSComplet = [
            (j, s)
            for j in setVars for s in setS[j]
        ]

        self.setVars = setVars
        self.setS = setS
        self.setSConcave = setSConcave
        self.setSConvex = setSConvex
        self.setSComplet = setSComplet

        # Variables
        x_ref = [cpx.variables.add(obj=[0] * len(setS[j]),
                                   lb=[0] * len(setS[j]),
                                   ub=[Lparam[j][s+1]-Lparam[j][s]
                                       for s in setS[j]],
                                   # ub=[cplex.infinity] * len(setS[j]),
                                   types=['C'] * len(setS[j]),
                                   names=['x_ref(%d)(%d)' % (j+1, s+1) for s in setS[j]])
                 for j in setVars]

        y_ref = [cpx.variables.add(obj=[0] * len(setS[j]),
                                   lb=[0] * len(setS[j]),
                                   ub=[1] * len(setS[j]),
                                   types=['B'] * len(setS[j]),
                                   names=['y_ref(%d)(%d)' % (j+1, s+1) for s in setS[j]])
                 for j in setVars]

        for j in setVars:
            cpx.variables.set_lower_bounds(y_ref[j][0], 1)

        self.setVars_WithZ = [j for j in setVars if self.varzid[j] > 0]
        self.setVars_NoZ = [j for j in setVars if self.varzid[j] <= 0]

        if len(self.setVars_WithZ) > 0:
            x_ref_bar = [-1] * len(setVars)

            for j in self.setVars_WithZ:
                x_ref_bar[j] = cpx.variables.add(obj=[0],
                                                 lb=[0],
                                                 ub=[Lparam[j][len(
                                                     setS[j])] - Lparam[j][0]],
                                                 #ub=[cplex.infinity] ,
                                                 types=['C'],
                                                 names=['x_ref_bar(%d)' % (j+1)])[0]

        for j in self.setVars_WithZ:
            y_ref[j] = list(y_ref[j])
            cpx.variables.set_types(y_ref[j][0], 'C')
            cpx.variables.set_upper_bounds(y_ref[j][0], 0)
            cpx.variables.set_lower_bounds(y_ref[j][0], 0)
            y_ref[j][0] = self.newVarIndex[self.varzid[j]-1]
            cpx.order.set(
                [(y_ref[j][0], 10000, cpx.order.branch_direction.up)])
        z_ref = [cpx.variables.add(obj=[0] * len(setS[j]),
                                   lb=[-cplex.infinity] * len(setS[j]),
                                   ub=[cplex.infinity] * len(setS[j]),
                                   types=['C'] * len(setS[j]),
                                   names=['z_ref(%d)(%d)' % (j+1, s+1) for s in setS[j]])
                 for j in setVars]

        self.x_ref = x_ref
        self.y_ref = y_ref
        self.z_ref = z_ref

        # Constraints
        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [varX[j]] + [x_ref_bar[j]] + [x_ref[j][k] for k in setS[j]],
                [1] + [-1] + [-1] * len([x_ref[j][k] for k in setS[j]]))],
            senses=['E'],
            rhs=[Lparam[j][0]])
         for j in self.setVars_WithZ]

        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [varX[j]] + [x_ref[j][k] for k in setS[j]],
                [1] + [-1] * len([x_ref[j][k] for k in setS[j]]))],
            senses=['E'],
            rhs=[Lparam[j][0]])
         for j in self.setVars_NoZ]

        #  #######################################

        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [x_ref_bar[j]] + [y_ref[j][0]],
                [1] + [(Lparam[j][len(setS[j])] - Lparam[j][0])])],
            senses=['L'],
            rhs=[(Lparam[j][len(setS[j])] - Lparam[j][0])])
         for j in self.setVars_WithZ]

        #  #######################################

        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [x_ref[j][s]] + [y_ref[j][s]],
                [1] + [-(Lparam[j][s+1] - Lparam[j][s])])],
            senses=['L'],
            rhs=[0])
         for j in setVars for s in setS[j]]

        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [x_ref[j][s]] + [y_ref[j][s+1]],
                [1] + [-(Lparam[j][s+1] - Lparam[j][s])])],
            senses=['G'],
            rhs=[0])
         for j in setVars for s in range(Nb[j]-2)]

        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [z_ref[j][s]] + [x_ref[j][s]],
                [1]+[(self.NLfunc[j](Lparam[j][s+1]) - self.NLfunc[j](Lparam[j][s]))/(Lparam[j][s+1]-Lparam[j][s])])],
            senses=['E'],
            rhs=[0])
         for j, s in setSConcave]

        for j in setVars:
            if (Ineq[j] == 1):
                newsense = 'L'
            else:
                newsense = 'G'
            cpx.linear_constraints.add(
                lin_expr=[cplex.SparsePair(
                    [varY[j]] + [z_ref[j][s] for s in setS[j]] + [y_ref[j][0]],
                    [1] + [-1]*len(setS[j]) + [self.NLfunc[j](Lparam[j][0])])],
                senses=[newsense],
                rhs=[0])

        # Creating Cuts
        nz_ref = {}

        # Clone variables
        for j, s in setSConvex:
            nz_ref[(j, s)] = cpx.variables.add(obj=[0],
                                               lb=[-cplex.infinity],
                                               ub=[cplex.infinity],
                                               types=['C'],
                                               names=['nz_ref(%d)(%d)' % (j+1, s+1)])[0]

        self.nz_ref = nz_ref

        # Clone variables constraints
        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [z_ref[j][s]] + [nz_ref[(j, s)]],
                [1] + [1])],
            senses=['E'],
            rhs=[0])
         for j, s in setSConvex]

        # Add initial constraints
        for j, s in setSConvex:
            for newValue in [0]+[Lparam[j][s+1] - Lparam[j][s]]:
                Lparamn = Lparam[j][s]
                if (Ineq[j] == 1):
                    newsense = 'G'
                else:
                    newsense = 'L'
                newrhs = 0
                valx = -self.NLfuncDiff[j](newValue+Lparamn)
                valy = -(self.NLfunc[j](newValue+Lparamn) - self.NLfunc[j]
                         (0+Lparamn) - self.NLfuncDiff[j](newValue+Lparamn)*newValue)
                cpx.linear_constraints.add(
                    lin_expr=[cplex.SparsePair(
                        [nz_ref[(j, s)]] + [x_ref[j][s]] + [y_ref[j][s]],
                        [1] + [valx] + [valy])],
                    senses=[newsense],
                    rhs=[newrhs])

        # Lparam
        # Mode: 1 Incremental, 2 Multiple choice
        self.lParamMode = [[self.Lparam[j][s]
                            for s in range(len(self.Lparam[j]))] for j in setVars]
        self.lParamFunc = [[self.NLfunc[j](self.Lparam[j][s]) for s in range(
            len(self.Lparam[j]))] for j in setVars]
        self.lParamDiff = [[self.NLfuncDiff[j](self.Lparam[j][s]) for s in range(
            len(self.Lparam[j]))] for j in setVars]

        # # Add initial constraints
        # for j,s in setSConvex :
        #     for newValue in [0]+[ Lparam[j][s+1] - Lparam[j][s] ] :
        #         valx = -self.NLfuncDiff[j](newValue+Lparam[j][s])
        #         valy = -(self.NLfunc[j](newValue+Lparam[j][s]) - self.NLfunc[j](0+Lparam[j][s]) - self.NLfuncDiff[j](newValue+Lparam[j][s])*newValue )
        #         cpx.linear_constraints.add(
        #                             lin_expr=[ cplex.SparsePair(
        #                                 [nz_ref[(j,s)]] + [ x_ref[j][s] ] + [ y_ref[j][s] ] ,
        #                                 [1] + [valx] + [valy] )],
        #                             senses=['G'],
        #                             rhs=[0] )

# MultipleChoice reformulation
    def CreateProblem_MultipleChoice(self):
        Nitem = self.Nitem
        # Umax = self.Umax
        # Cap = self.Cap
        # weight = self.weight
        # UBvar = self.UBvar
        Nb = self.Nb
        Lparam = self.Lparam
        Concave = self.Concave
        Ineq = self.ineq

        cpx = self.cpx
        varX = self.varX
        varY = self.varY

        # Sets
        setVars = list(range(Nitem))
        setS = [list(range(Nb[j]-1)) for j in setVars]

        setSConcave = [
            (j, s)
            for j in setVars for s in setS[j]
            if int(Concave[j][s]) == 1
        ]

        setSConvex = [
            (j, s)
            for j in setVars for s in setS[j]
            if int(Concave[j][s]) == -1
        ]

        setSComplet = [
            (j, s)
            for j in setVars for s in setS[j]
        ]

        self.setVars = setVars
        self.setS = setS
        self.setSConcave = setSConcave
        self.setSConvex = setSConvex
        self.setSComplet = setSComplet

        # Variables
        x_ref = [cpx.variables.add(obj=[0] * len(setS[j]),
                                   lb=[Lparam[j][0]] * len(setS[j]),
                                   # ub=[Lparam[j][s+1] for s in setS[j] ],
                                   ub=[cplex.infinity] * len(setS[j]),
                                   types=['C'] * len(setS[j]),
                                   names=['x_ref(%d)(%d)' % (j+1, s+1) for s in setS[j]])
                 for j in setVars]

        y_ref = [cpx.variables.add(obj=[0] * len(setS[j]),
                                   lb=[0] * len(setS[j]),
                                   ub=[1] * len(setS[j]),
                                   types=['B'] * len(setS[j]),
                                   names=['y_ref(%d)(%d)' % (j+1, s+1) for s in setS[j]])
                 for j in setVars]

        z_ref = [cpx.variables.add(obj=[0] * len(setS[j]),
                                   lb=[-cplex.infinity] * len(setS[j]),
                                   ub=[cplex.infinity] * len(setS[j]),
                                   types=['C'] * len(setS[j]),
                                   names=['z_ref(%d)(%d)' % (j+1, s+1) for s in setS[j]])
                 for j in setVars]

        self.x_ref = x_ref
        self.y_ref = y_ref
        self.z_ref = z_ref

        # Constraints
        # (17)
        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [y_ref[j][s] for s in setS[j]],
                [1] * len(setS[j]))],
            senses=['E'],
            rhs=[1])
         for j in setVars]

        # (15)
        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [varX[j]] + [x_ref[j][k] for k in setS[j]],
                [1] + [-1] * len([x_ref[j][k] for k in setS[j]]))],
            senses=['E'],
            rhs=[0])
         for j in setVars]

        # (16)
        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [x_ref[j][s]] + [y_ref[j][s]],
                [1] + [-Lparam[j][s+1]])],
            senses=['L'],
            rhs=[0])
         for j in setVars for s in setS[j]]

        # (16)
        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [x_ref[j][s]] + [y_ref[j][s]],
                [1] + [-Lparam[j][s]])],
            senses=['G'],
            rhs=[0])
         for j in setVars for s in setS[j]]

        # z_ref concave
        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [z_ref[j][s]] + [x_ref[j][s]] + [y_ref[j][s]],
                [1]+[(self.NLfunc[j](Lparam[j][s+1]) - self.NLfunc[j](Lparam[j][s]))/(Lparam[j][s+1]-Lparam[j][s])] + [self.NLfunc[j](Lparam[j][s]) - Lparam[j][s]*(self.NLfunc[j](Lparam[j][s+1]) - self.NLfunc[j](Lparam[j][s]))/(Lparam[j][s+1]-Lparam[j][s])])],
            senses=['E'],
            rhs=[0])
         for j, s in setSConcave]

        # (13)
        for j in setVars:
            if (Ineq[j] == 1):
                newsense = 'L'
            else:
                newsense = 'G'
            cpx.linear_constraints.add(
                lin_expr=[cplex.SparsePair(
                    [varY[j]] + [z_ref[j][s] for s in setS[j]] + [y_ref[k][s]
                                                                  for k, s in setSConvex if k == j],
                    [1] + [-1]*len(setS[j]) + [self.NLfunc[j](Lparam[j][0])] * len([y_ref[k][s] for k, s in setSConvex if k == j]))],
                senses=[newsense],
                rhs=[0])

        # Creating Cuts
        nz_ref = {}

        # Clone variables
        for j, s in setSConvex:
            nz_ref[(j, s)] = cpx.variables.add(obj=[0],
                                               lb=[-cplex.infinity],
                                               ub=[cplex.infinity],
                                               types=['C'],
                                               names=['nz_ref(%d)(%d)' % (j+1, s+1)])[0]

        self.nz_ref = nz_ref

        # Clone variables constraints
        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [z_ref[j][s]] + [nz_ref[(j, s)]],
                [1] + [1])],
            senses=['E'],
            rhs=[0])
         for j, s in setSConvex]

        # Add initial constraints
        for j, s in setSConvex:
            for newValue in [Lparam[j][s]]+[Lparam[j][s+1]]:
                Lparamn = Lparam[j][s]
                Lparamn = 0
                if (Ineq[j] == 1):
                    newsense = 'G'
                else:
                    newsense = 'L'
                newrhs = 0
                valx = -self.NLfuncDiff[j](newValue+Lparamn)
                valy = -(self.NLfunc[j](newValue+Lparamn) - self.NLfunc[j]
                         (0+Lparamn) - self.NLfuncDiff[j](newValue+Lparamn)*newValue)
                cpx.linear_constraints.add(
                    lin_expr=[cplex.SparsePair(
                        [nz_ref[(j, s)]] + [x_ref[j][s]] + [y_ref[j][s]],
                        [1] + [valx] + [valy])],
                    senses=[newsense],
                    rhs=[newrhs])

        # Lparam
        # Mode: 1 Incremental, 2 Multiple choice
        self.lParamMode = [
            [0 for s in range(len(self.Lparam[j]))] for j in setVars]
        self.lParamFunc = [[self.NLfunc[j](0) for s in range(
            len(self.Lparam[j]))] for j in setVars]

# MultipleChoice reformulation
    def CreateProblem_MultipleChoiceTest(self):
        Nitem = self.Nitem
        # Umax = self.Umax
        # Cap = self.Cap
        # weight = self.weight
        # UBvar = self.UBvar
        Nb = self.Nb
        Lparam = self.Lparam
        Concave = self.Concave
        Ineq = self.ineq

        cpx = self.cpx
        varX = self.varX
        varY = self.varY

        # Sets
        setVars = list(range(Nitem))
        setS = [list(range(Nb[j]-1)) for j in setVars]

        setSConcave = [
            (j, s)
            for j in setVars for s in setS[j]
            if int(Concave[j][s]) == 1
        ]

        setSConvex = [
            (j, s)
            for j in setVars for s in setS[j]
            if int(Concave[j][s]) == -1
        ]

        setSComplet = [
            (j, s)
            for j in setVars for s in setS[j]
        ]

        self.setVars = setVars
        self.setS = setS
        self.setSConcave = setSConcave
        self.setSConvex = setSConvex
        self.setSComplet = setSComplet

        # Variables
        x_ref = [cpx.variables.add(obj=[0] * len(setS[j]),
                                   # lb=[0] * len(setS[j]), Lparam[j][s+1]
                                   lb=[-cplex.infinity] * len(setS[j]),
                                   # ub=[Lparam[j][s+1] for s in setS[j] ],
                                   ub=[cplex.infinity] * len(setS[j]),
                                   types=['C'] * len(setS[j]),
                                   names=['x_ref(%d)(%d)' % (j+1, s+1) for s in setS[j]])
                 for j in setVars]

        y_ref = [cpx.variables.add(obj=[0] * len(setS[j]),
                                   lb=[0] * len(setS[j]),
                                   ub=[1] * len(setS[j]),
                                   types=['B'] * len(setS[j]),
                                   names=['y_ref(%d)(%d)' % (j+1, s+1) for s in setS[j]])
                 for j in setVars]

        z_ref = [cpx.variables.add(obj=[0] * len(setS[j]),
                                   lb=[-cplex.infinity] * len(setS[j]),
                                   ub=[cplex.infinity] * len(setS[j]),
                                   types=['C'] * len(setS[j]),
                                   names=['z_ref(%d)(%d)' % (j+1, s+1) for s in setS[j]])
                 for j in setVars]

        self.x_ref = x_ref
        self.y_ref = y_ref
        self.z_ref = z_ref

        self.setVars_WithZ = [j for j in setVars if self.varzid[j] > 0]
        self.setVars_NoZ = [j for j in setVars if self.varzid[j] <= 0]

        # self.newVarIndex[self.varzid[j]
        
        # Constraints
        # (17)
        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [y_ref[j][s] for s in setS[j]],
                [1] * len(setS[j]))],
            senses=['E'],
            rhs=[1])
         for j in self.setVars_NoZ]

        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [y_ref[j][s] for s in setS[j]] +
                [self.newVarIndex[self.varzid[j]-1]],
                [1] * len(setS[j]) + [-1])],
            senses=['E'],
            rhs=[0])
         for j in self.setVars_WithZ]

        # (15)
        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [varX[j]] + [x_ref[j][k] for k in setS[j]],
                [1] + [-1] * len([x_ref[j][k] for k in setS[j]]))],
            senses=['E'],
            rhs=[0])
         for j in setVars]

        # (16)
        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [x_ref[j][s]] + [y_ref[j][s]],
                [1] + [-Lparam[j][s+1]])],
            senses=['L'],
            rhs=[0])
         for j in setVars for s in setS[j]]

        # (16)
        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [x_ref[j][s]] + [y_ref[j][s]],
                [1] + [-Lparam[j][s]])],
            senses=['G'],
            rhs=[0])
         for j in setVars for s in setS[j]]

        # z_ref concave
        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [z_ref[j][s]] + [x_ref[j][s]] + [y_ref[j][s]],
                [1] +
                [(self.NLfunc[j](Lparam[j][s+1]) - self.NLfunc[j](Lparam[j][s]))/(Lparam[j][s+1]-Lparam[j][s])] +
                [self.NLfunc[j](Lparam[j][s]) - Lparam[j][s]*((self.NLfunc[j](
                    Lparam[j][s+1]) - self.NLfunc[j](Lparam[j][s]))/(Lparam[j][s+1]-Lparam[j][s]))]
            )],
            senses=['E'],
            rhs=[0])
         for j, s in setSConcave]

        # (13)
        for j in setVars:
            if (Ineq[j] == 1):
                newsense = 'L'
            else:
                newsense = 'G'
            cpx.linear_constraints.add(
                lin_expr=[cplex.SparsePair(
                    [varY[j]] + [z_ref[j][s] for s in setS[j]],
                    [1] + [-1]*len(setS[j]))],
                senses=[newsense],
                rhs=[0])

        # Creating Cuts
        nz_ref = {}

        # Clone variables
        for j, s in setSConvex:
            nz_ref[(j, s)] = cpx.variables.add(obj=[0],
                                               lb=[-cplex.infinity],
                                               ub=[cplex.infinity],
                                               types=['C'],
                                               names=['nz_ref(%d)(%d)' % (j+1, s+1)])[0]

        self.nz_ref = nz_ref

        # Clone variables constraints
        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [z_ref[j][s]] + [nz_ref[(j, s)]] + [y_ref[j][s]],
                [1] + [1] + [self.NLfunc[j](0.00)])],
            senses=['E'],
            rhs=[0])
         for j, s in setSConvex]

        # Add initial constraints
        for j, s in setSConvex:
            for newValue in [Lparam[j][s]]+[Lparam[j][s+1]]:
                Lparamn = Lparam[j][s]
                Lparamn = 0
                if (Ineq[j] == 1):
                    newsense = 'G'
                else:
                    newsense = 'L'
                newrhs = 0
                valx = -self.NLfuncDiff[j](newValue+Lparamn)
                valy = -(self.NLfunc[j](newValue+Lparamn) - self.NLfunc[j]
                         (0+Lparamn) - self.NLfuncDiff[j](newValue+Lparamn)*newValue)
                cpx.linear_constraints.add(
                    lin_expr=[cplex.SparsePair(
                        [nz_ref[(j, s)]] + [x_ref[j][s]] + [y_ref[j][s]],
                        [1] + [valx] + [valy])],
                    senses=[newsense],
                    rhs=[newrhs])

        # Lparam
        # Mode: 1 Incremental, 2 Multiple choice
        self.lParamMode = [
            [0 for s in range(len(self.Lparam[j]))] for j in setVars]
        self.lParamFunc = [[self.NLfunc[j](0) for s in range(
            len(self.Lparam[j]))] for j in setVars]

#MultipleChoice-Incremental reformulation
    def CreateProblem_MultipleChoiceInc(self):
        Nitem = self.Nitem
        # Umax = self.Umax
        # Cap = self.Cap
        # weight = self.weight
        # UBvar = self.UBvar
        Nb = self.Nb
        Lparam = self.Lparam
        Concave = self.Concave
        Ineq = self.ineq

        cpx = self.cpx
        varX = self.varX
        varY = self.varY

        # Sets
        setVars = list(range(Nitem))
        setS = [list(range(Nb[j]-1)) for j in setVars]

        setSConcave = [
            (j, s)
            for j in setVars for s in setS[j]
            if int(Concave[j][s]) == 1
        ]

        setSConvex = [
            (j, s)
            for j in setVars for s in setS[j]
            if int(Concave[j][s]) == -1
        ]

        setSComplet = [
            (j, s)
            for j in setVars for s in setS[j]
        ]

        self.setVars = setVars
        self.setS = setS
        self.setSConcave = setSConcave
        self.setSConvex = setSConvex
        self.setSComplet = setSComplet

        # Variables
        x_ref = [cpx.variables.add(obj=[0] * len(setS[j]),
                                   # lb=[0] * len(setS[j]), Lparam[j][s+1]
                                   lb=[-cplex.infinity] * len(setS[j]),
                                   # ub=[Lparam[j][s+1] for s in setS[j] ],
                                   ub=[cplex.infinity] * len(setS[j]),
                                   types=['C'] * len(setS[j]),
                                   names=['x_ref(%d)(%d)' % (j+1, s+1) for s in setS[j]])
                 for j in setVars]

        y_ref = [cpx.variables.add(obj=[0] * len(setS[j]),
                                   lb=[0] * len(setS[j]),
                                   ub=[1] * len(setS[j]),
                                   types=['C'] * len(setS[j]),
                                   names=['y_ref(%d)(%d)' % (j+1, s+1) for s in setS[j]])
                 for j in setVars]

        self.setVars_WithZ = [j for j in setVars if self.varzid[j] > 0]
        self.setVars_NoZ = [j for j in setVars if self.varzid[j] <= 0]

        if len(self.setVars_WithZ) > 0:
            x_ref_bar = [-1] * len(setVars)

            for j in self.setVars_WithZ:
                x_ref_bar[j] = cpx.variables.add(obj=[0],
                                                 lb=[0],
                                                 ub=[Lparam[j][len(
                                                     setS[j])] - Lparam[j][0]],
                                                 #ub=[cplex.infinity] ,
                                                 types=['C'],
                                                 names=['x_ref_bar(%d)' % (j+1)])[0]

        for j in self.setVars_WithZ:
            y_ref[j] = list(y_ref[j])
            cpx.variables.set_types(y_ref[j][0], 'C')
            cpx.variables.set_upper_bounds(y_ref[j][0], 0)
            cpx.variables.set_lower_bounds(y_ref[j][0], 0)
            y_ref[j][0] = self.newVarIndex[self.varzid[j]-1]
            cpx.order.set(
                [(y_ref[j][0], 10000, cpx.order.branch_direction.up)])

        y2_ref = [cpx.variables.add(obj=[0] * len(range(Nb[j])),
                                   lb=[0] * len(range(Nb[j])),
                                   ub=[1] * len(range(Nb[j])),
                                   types=['B'] * len(range(Nb[j])),
                                   names=['y2_ref(%d)(%d)' % (j+1, s+1) for s in range(Nb[j]) ])
                 for j in setVars]
        
        for j in setVars:
            cpx.variables.set_lower_bounds(y2_ref[j][0], 1)
            cpx.variables.set_upper_bounds(y2_ref[j][-1], 0)

        z_ref = [cpx.variables.add(obj=[0] * len(setS[j]),
                                   lb=[-cplex.infinity] * len(setS[j]),
                                   ub=[cplex.infinity] * len(setS[j]),
                                   types=['C'] * len(setS[j]),
                                   names=['z_ref(%d)(%d)' % (j+1, s+1) for s in setS[j]])
                 for j in setVars]

        self.x_ref = x_ref
        self.y_ref = y_ref
        self.z_ref = z_ref

        # (14) {i in VARS, k in Intervals[i]}:
        #   y_ref[i,k] = y2_ref[i,k] - y2_ref[i,k+1];
        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [y_ref[j][s]] + [y2_ref[j][s]] + [y2_ref[j][s+1]],
                [1] + [-1] + [1])],
            senses=['E'],
            rhs=[0])
         for j in setVars for s in range(Nb[j]-1)]

        # Constraints
        # (15){i in VARS} : 
        #   x[i] == sum{ k in Intervals[i]} x_ref[i,k];
        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [varX[j]] + [x_ref[j][k] for k in setS[j]],
                [1] + [-1] * len([x_ref[j][k] for k in setS[j]]))],
            senses=['E'],
            rhs=[0])
         for j in setVars]
        
        # (14) {i in VARS, k in Intervals[i]}:
        #   y_ref[i,k] >= y_ref[i,k+1];
        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [y2_ref[j][s]] + [y2_ref[j][s+1]],
                [1] + [-1])],
            senses=['G'],
            rhs=[0])
         for j in setVars for s in range(Nb[j]-2)]

        # (16) {i in VARS, k in Intervals[i]}:
	    #   x_ref[i,k] <= y_ref[i,k]*breakpointsL[i,k+1] - y_ref[i,k+1]*breakpointsL[i,k+1]
        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [x_ref[j][s]] + [y2_ref[j][s]] + [y2_ref[j][s+1]],
                [1] + [-Lparam[j][s+1]] + [Lparam[j][s+1]])],
            senses=['L'],
            rhs=[0])
         for j in setVars for s in setS[j]]

        # (16) {i in VARS, k in Intervals[i]}:
        #   x_ref[i,k] >= y_ref[i,k]*breakpointsL[i,k] - y_ref[i,k+1]*breakpointsL[i,k]
        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [x_ref[j][s]] + [y2_ref[j][s]] + [y2_ref[j][s+1]],
                [1] + [-Lparam[j][s]] + [Lparam[j][s]])],
            senses=['G'],
            rhs=[0])
         for j in setVars for s in setS[j]]

        # z_ref concave
        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [z_ref[j][s]] + [x_ref[j][s]] + [y2_ref[j][s]] + [y2_ref[j][s+1]],
                [1] +
                [(self.NLfunc[j](Lparam[j][s+1]) - self.NLfunc[j](Lparam[j][s]))/(Lparam[j][s+1]-Lparam[j][s])] +
                [self.NLfunc[j](Lparam[j][s]) - Lparam[j][s]*((self.NLfunc[j](
                    Lparam[j][s+1]) - self.NLfunc[j](Lparam[j][s]))/(Lparam[j][s+1]-Lparam[j][s]))] +
                [-(self.NLfunc[j](Lparam[j][s]) - Lparam[j][s]*((self.NLfunc[j](
                    Lparam[j][s+1]) - self.NLfunc[j](Lparam[j][s]))/(Lparam[j][s+1]-Lparam[j][s])))]
            )],
            senses=['E'],
            rhs=[0])
         for j, s in setSConcave]

        # (13)
        for j in setVars:
            if (Ineq[j] == 1):
                newsense = 'L'
            else:
                newsense = 'G'
            cpx.linear_constraints.add(
                lin_expr=[cplex.SparsePair(
                    [varY[j]] + [z_ref[j][s] for s in setS[j]],
                    [1] + [-1]*len(setS[j]))],
                senses=[newsense],
                rhs=[0])

        # Creating Cuts
        nz_ref = {}

        # Clone variables
        for j, s in setSConvex:
            nz_ref[(j, s)] = cpx.variables.add(obj=[0],
                                               lb=[-cplex.infinity],
                                               ub=[cplex.infinity],
                                               types=['C'],
                                               names=['nz_ref(%d)(%d)' % (j+1, s+1)])[0]

        self.nz_ref = nz_ref

        # Clone variables constraints
        [cpx.linear_constraints.add(
            lin_expr=[cplex.SparsePair(
                [z_ref[j][s]] + [nz_ref[(j, s)]] + [y_ref[j][s]],
                [1] + [1] + [self.NLfunc[j](0.00)])],
            senses=['E'],
            rhs=[0])
         for j, s in setSConvex]

        # Add initial constraints
        for j, s in setSConvex:
            for newValue in [Lparam[j][s]]+[Lparam[j][s+1]]:
                Lparamn = Lparam[j][s]
                Lparamn = 0
                if (Ineq[j] == 1):
                    newsense = 'G'
                else:
                    newsense = 'L'
                newrhs = 0
                valx = -self.NLfuncDiff[j](newValue+Lparamn)
                valy = -(self.NLfunc[j](newValue+Lparamn) - self.NLfunc[j]
                         (0+Lparamn) - self.NLfuncDiff[j](newValue+Lparamn)*newValue)
                cpx.linear_constraints.add(
                    lin_expr=[cplex.SparsePair(
                        [nz_ref[(j, s)]] + [x_ref[j][s]] + [y_ref[j][s]],
                        [1] + [valx] + [valy])],
                    senses=[newsense],
                    rhs=[newrhs])

        # Lparam
        # Mode: 1 Incremental, 2 Multiple choice
        self.lParamMode = [
            [0 for s in range(len(self.Lparam[j]))] for j in setVars]
        self.lParamFunc = [[self.NLfunc[j](0) for s in range(
            len(self.Lparam[j]))] for j in setVars]


    def PrintProblem(self):
        lp = self.cpx.write_as_string("lp")
        print(lp)

    def ExportFinalLP(self):
        cpxCopy = cplex.Cplex(self.cpx)

        # Relax integrality
        for i in range(cpxCopy.variables.get_num()):
            if cpxCopy.variables.get_types(i) == 'B':
                cpxCopy.variables.set_types(
                    i, cpxCopy.variables.type.continuous)

        for i in range(cpxCopy.variables.get_num()):
            cpxCopy.variables.set_lower_bounds(i, self.local_lb[i])
            cpxCopy.variables.set_upper_bounds(i, self.local_ub[i])

        # Add cuts from the pool
        cpxCopy.linear_constraints.add(
            lin_expr=self.poolconst,
            senses=self.poolsense,
            rhs=self.poolrhs)

        # lp = cpxCopy.write_as_string("lp")
        cpxCopy.write("pippo.lp")

    def getSolution(self):
        # solutionList = self.cpx.solution.get_values(self.newVarIndex)
        solutionList = [self.varValGlobalUB[i] for i in self.newVarIndex]
        return solutionList

    def getSolutionLB(self):
        solutionList = self.cpx.solution.get_values(self.newVarIndex)
        # solutionList = [self.varValGlobalUB[i] for i in self.newVarIndex]
        return solutionList

    def ExportSol(self):

        solutionList = self.cpx.solution.get_values(self.newVarIndex)

        f = open("tmp/solution.run", "w")
        for i in range(0, self.originalVarNum):
            # f.write("let _var[" + str(i+1) + "] := " + str(self.cpx.solution.get_values( i )) + ";\n")
            f.write("let _var[" + str(i+1) + "] := " +
                    str(solutionList[i]) + ";\n")
        # for k in range(len(self.varyid)) :
        #     i = self.varyid[k]
        #     j = self.varxid[k]
        #     XValue = self.cpx.solution.get_values( j )
        #     newYValue = -self.NLfunc[k](XValue)
        #     f.write("let _var[" + str(i+1) + "] := " + str(newYValue) + ";\n")
        f.close()

        nconstraints = self.nFunctions + self.originalConsNum

        f = open("tmp/solution.sol", "w")
        f.write("Solved\n\n\nOptions\n3\n0\n1\n0\n")
        f.write(str(nconstraints) + "\n")
        f.write("0\n")
        f.write(str(self.originalVarNum) + "\n")
        f.write(str(self.originalVarNum) + "\n")
        for i in range(0, self.originalVarNum):
            f.write(str(solutionList[i]) + "\n")
        f.close()

    def ExportSolIpopt(self):
        f = open("tmp/solution.run", "w")
        for i in range(0, self.originalVarNum):
            f.write("let _var[" + str(i+1) + "] := " +
                    str(self.varsolipopt[i]) + ";\n")
        for k in range(len(self.varyid)):
            i = self.varyid[k]
            j = self.varxid[k]
            XValue = self.varsolipopt[j]
            newYValue = -self.NLfunc[k](XValue)
            f.write("let _var[" + str(i+1) + "] := " + str(newYValue) + ";\n")
        f.close()

    def IpoptCall(self, filename_in="tmp/prob.nl", filename_out="tmp/prob2.nl"):
        from inputnl import nlFile

        nlfile = nlFile()
        nlfile.loadnl(filename_in)

        x = []
        for i in range(0, self.originalVarNum):
            x.append(self.cpx.solution.get_values(i))

        for k in range(len(self.varyid)):
            i = self.varyid[k]
            j = self.varxid[k]
            XValue = self.cpx.solution.get_values(j)
            newYValue = -self.NLfunc[k](XValue)
            x[i] = newYValue

        nlfile.exportFile(x, filename_out)

        print("\n\n******************************************************************************")
        print("START IPOPT\n")
        # stream = os.system("ipopt " + filename_out + " wantsol=10 print_level=3")
        stream = os.popen("ipopt " + filename_out +
                          " wantsol=10 print_level=3")
        output = stream.read()
        firstDivision = output.split("EXIT:")
        # print(output)

        # Solution
        sol = firstDivision[1].split("value\n")[1].splitlines()

        for i in range(len(sol)):
            pairValues = sol[i].split()
            x[i] = float(pairValues[1])

        # Solution info
        info = firstDivision[0].split("Objective...............:")[
            1].splitlines()
        obj = float(info[0].split()[0])
        # dualinf = float(info[1].split("Dual infeasibility......:")[1].split()[1])
        # consviol = float(info[2].split("Constraint violation....:")[1].split()[1])
        # varviol = float(info[3].split("Variable bound violation:")[1].split()[1])
        # comp = float(info[4].split("Complementarity.........:")[1].split()[1])
        NLPerror = float(info[5].split(
            "Overall NLP error.......:")[1].split()[1])

        Status = firstDivision[1].splitlines()[0].lstrip()

        totaltime = float(firstDivision[0].split(
            "Total seconds in IPOPT")[1].split()[1])

        self.varsolipopt = x

        print(Status)
        print("Optimal value:                      " + str(obj))
        print("Solution time:                      " + str(totaltime))
        print()

        print("END IPOPT")
        print("******************************************************************************\n\n")

    def RelaxIntegrality(self):
        try:
            # Parameter to limit the number of nodes
            self.cpx.parameters.mip.limits.nodes.set(0)  # Important

            # Converting all variables to continuous
            for i in range(self.cpx.variables.get_num()):
                if self.cpx.variables.get_types(i) != 'C':
                    self.cpx.variables.set_types(
                        i, self.cpx.variables.type.continuous)
        except cplex.exceptions.errors.CplexSolverError as err:
            print("Relaxed problem error!\n")

    def ConfigParametersTotal(self):
        cpx = self.cpx

        cpx.parameters.mip.strategy.search.set(1)
        cpx.parameters.mip.strategy.lbheur.set(0)

        cpx.parameters.preprocessing.repeatpresolve.set(0)
        cpx.parameters.preprocessing.linear.set(1)
        cpx.parameters.preprocessing.aggregator.set(0)
        cpx.parameters.preprocessing.presolve.set(0)  # Important
        cpx.parameters.preprocessing.numpass.set(0)
        cpx.parameters.preprocessing.reduce.set(0)
        cpx.parameters.preprocessing.linear.set(0)
        cpx.parameters.preprocessing.relax.set(0)

        cpx.parameters.preprocessing.dual.set(-1)
        cpx.parameters.simplex.dgradient.set(1)

        cpx.parameters.mip.cuts.disjunctive.set(-1)
        cpx.parameters.mip.cuts.gubcovers.set(-1)
        cpx.parameters.mip.cuts.mcfcut.set(-1)
        cpx.parameters.mip.cuts.localimplied.set(-1)
        cpx.parameters.mip.cuts.rlt.set(-1)

        cpx.parameters.mip.limits.cutsfactor.set(0)

        cpx.parameters.mip.display.set(5)
        cpx.parameters.mip.display.set(2)

        cpx.parameters.simplex.tolerances.optimality.set(1e-9)
        cpx.parameters.simplex.tolerances.feasibility.set(1e-9)

        # OFF ###################################################################
        # cpx.parameters.mip.tolerances.integrality.set(1e-6)       #Default 1e-5
        # cpx.parameters.mip.internal.set(1)
        # CPXPARAM_MIP_Interval                            1

        # cpx.parameters.mip.limits.nodes.set(1)            #Important
        # cpx.parameters.mip.cuts.nodecuts.set(3)

        # cpx.parameters.mip.strategy.callbackreducedlp.set(0)
        # CPXPARAM_MIP_Strategy_CallbackReducedLP          0

    def setParametersPaperCode(self):
        cpx = self.cpx
        cpx.parameters.reset()
        #  CPLEX Parameter File Version 12.10.0.0
        #  CPXPARAM_Preprocessing_Reduce                    1
        cpx.parameters.preprocessing.reduce.set(1)
        #  CPXPARAM_Preprocessing_Linear                    0
        cpx.parameters.preprocessing.linear.set(0)
        #  CPXPARAM_Threads                                 1
        cpx.parameters.threads.set(4)
        #  CPXPARAM_MIP_Cuts_FlowCovers                     -1
        cpx.parameters.mip.cuts.flowcovers.set(-1)
        #  CPXPARAM_MIP_Strategy_CallbackReducedLP          0
        #  CPXPARAM_MIP_Interval                            1
        # cpx.parameters.mip.interval.set(1)
        #  CPXPARAM_Simplex_Tolerances_Optimality           1.0000000000000001e-09
        cpx.parameters.simplex.tolerances.optimality.set(1e-9)
        #  CPXPARAM_Simplex_Tolerances_Feasibility          1.0000000000000001e-09
        cpx.parameters.simplex.tolerances.feasibility.set(1e-9)

        # New parameters
        cpx.parameters.read.datacheck.set(0)
        cpx.parameters.mip.tolerances.mipgap.set(self.CPLEX_MIPGAP)
        cpx.parameters.mip.tolerances.absmipgap.set(self.CPLEX_ABSMIPGAP)

        # cpx.parameters.simplex.tolerances.markowitz.set(0.99999)
        # cpx.parameters.simplex.perturbation.constant.set(1e-8)
        # cpx.parameters.mip.strategy.probe.set(3)
        version = cpx.get_version()
        if version == '22.1.0.0' or version == '20.1.0.0':
            cpx.parameters.preprocessing.reformulations.set(2)
        cpx.parameters.preprocessing.repeatpresolve.set(0)
        # cpx.parameters.mip.strategy.variableselect.set(2)

        # cpx.parameters.mip.strategy.heuristicfreq.set(-1)
        # cpx.parameters.mip.strategy.search.set(1)

        if (self.iterations != 0):
            self.cpx.parameters.paramdisplay.set(0)

    def SolveIterative(self):
        # self.cpx.parameters.mip.display.set(0)
        # self.cpx.set_log_stream(None)
        # self.cpx.set_error_stream(None)
        # self.cpx.set_warning_stream(None)
        # self.cpx.set_results_stream(None)

        addedcuts = True

        poolconst_temp = []
        poolsense_temp = []
        poolrhs_temp = []

        ncuts = 0

        while addedcuts:
            # self.cpx.write("tmp/prob.lp")
            self.cpx.solve()

            addedcuts = False

            for j, s in self.setSConvex:
                varX = self.cpx.solution.get_values(self.x_ref[j][s])
                varY = self.cpx.solution.get_values(self.y_ref[j][s])
                VarZ = self.cpx.solution.get_values(self.nz_ref[(j, s)])

                if self.p_mode == 1:
                    lParam = self.Lparam[j][s]
                else:
                    if self.p_mode == 2:
                        lParam = 0

                if varY > 0:
                    # lParam = self.Lparam[j][s]
                    xbar = varX/varY
                    # xbar = varX

                    # Calculating violation
                    violation = varY * \
                        self.NLfunc[j](xbar+lParam) - varY * \
                        self.NLfunc[j](0+lParam) - VarZ
                    # violation = self.NLfunc[j](xbar+lParam) - self.NLfunc[j](0+lParam) - VarZ

                    if (violation > self.PC_THRSHLD):
                        xCoeff = -self.NLfuncDiff[j](xbar+lParam)
                        yCoeff = -(self.NLfunc[j](xbar+lParam) - self.NLfunc[j]
                                   (0+lParam) - self.NLfuncDiff[j](xbar+lParam)*xbar)
                        # yCoeff = 0

                        newconst = cplex.SparsePair([self.x_ref[j][s]] + [self.y_ref[j][s]] + [self.nz_ref[(j, s)]],
                                                    [xCoeff]+[yCoeff]+[1])
                        newsense = 'G'
                        newrhs = 0
                        # newrhs = (self.NLfunc[j](xbar+lParam) - self.NLfunc[j](0+lParam) - self.NLfuncDiff[j](xbar+lParam)*xbar)

                        poolconst_temp.append(newconst)
                        poolsense_temp.append(newsense)
                        poolrhs_temp.append(newrhs)

                        ncuts += 1
                        addedcuts = True

            if addedcuts == True:
                self.cpx.linear_constraints.add(
                    lin_expr=poolconst_temp,
                    senses=poolsense_temp,
                    rhs=poolrhs_temp)

            self.poolconst.append(poolconst_temp)
            self.poolsense.append(poolsense_temp)
            self.poolrhs.append(poolrhs_temp)
        print("Total cuts:", ncuts)

    def Solve(self, option=1):
        'Solve the problem already defined.\n option = 1 for new implementation \n option = 2 for legacy.'
        cpx = self.cpx

        self.setParametersPaperCode()

        if self.display >= 0:
            # Display options
            self.cpx.parameters.mip.display.set(self.display)
            if self.display == 0:
                self.cpx.set_log_stream(None)
                self.cpx.set_error_stream(None)
                self.cpx.set_warning_stream(None)
                self.cpx.set_results_stream(None)

        if self.relax >= 1:
            self.RelaxIntegrality()

        if self.timelimit >= 0:
            if (self.timelimit - self.totaltime > 0):
                cpx.parameters.timelimit.set(self.timelimit - self.totaltime)

        # Callback
        # cpx.parameters.paramdisplay.set(0)
        cpx.parameters.write_file("current_param.prm")

        # Cleaning callbacks registration
        cpx.set_callback(None)
        cpx.unregister_callback(NLProblem_Callback_LazyConstraint)
        cpx.unregister_callback(NLProblem_Callback_UserCut)

        # cpx.unregister_callback(NLProblem_Callback_SolveCallback)
        # cpx.unregister_callback(NLProblem_Callback_BranchCallback)
        # cpx.unregister_callback(NLProblem_Callback_NodeCallback)
        # cpx.unregister_callback(NLProblem_Callback_Continuous)
        # cpx.unregister_callback(NLProblem_Callback_Heuristic)
        # cpx.unregister_callback(NLProblem_Callback_Incumbent)

        for j, s in self.setSConvex:
            self.prevxbar[(j, s)] = -1

        # self.cutscount = 0

        totaltime = 0
        if (option == 0):
            timestart = self.cpx.get_time()
            self.SolveIterative()
            totaltime = self.cpx.get_time() - timestart

        if (option == 1):
            contextmask = 0
            contextmask |= cplex.callbacks.Context.id.candidate
            contextmask |= cplex.callbacks.Context.id.relaxation

            # contextmask |= cplex.callbacks.Context.id.branching
            # contextmask |= cplex.callbacks.Context.id.local_progress
            # contextmask |= cplex.callbacks.Context.id.global_progress

            cpx.set_callback(self, contextmask)

            timestart = self.cpx.get_time()
            cpx.solve()
            self.totaltime = self.cpx.get_time() - timestart

        if (option == 2):
            cpx.register_callback(NLProblem_Callback_LazyConstraint)
            cpx.register_callback(NLProblem_Callback_UserCut)

            # cpx.register_callback(NLProblem_Callback_SolveCallback)
            # cpx.register_callback(NLProblem_Callback_BranchCallback)
            # cpx.register_callback(NLProblem_Callback_NodeCallback)
            # cpx.register_callback(NLProblem_Callback_Continuous)
            # cpx.register_callback(NLProblem_Callback_Heuristic)
            # cpx.register_callback(NLProblem_Callback_Incumbent)

            NLProblem_Callback_LazyConstraint.NLProblemData = self
            NLProblem_Callback_UserCut.NLProblemData = self
            # NLProblem_Callback_Incumbent.NLProblemData = self

            timestart = self.cpx.get_time()
            cpx.solve()
            self.totaltime += self.cpx.get_time() - timestart

        currentBestBound = self.cpx.solution.MIP.get_best_objective()

        if self.cpx.objective.sense[self.cpx.objective.get_sense()] == 'minimize':
            if self.objValGlobalLB < currentBestBound:
                self.objValGlobalLB = currentBestBound

        elif self.cpx.objective.sense[self.cpx.objective.get_sense()] == 'maximize':
            if self.objValGlobalLB > currentBestBound:
                self.objValGlobalLB = currentBestBound

        print('\n*** Iteration ' + str(self.iterations) +
              ', current Lower Bounding problem solved ***')
        print('Solution status:                   ',
              self.cpx.solution.get_status_string())
        print('Best solution:                     ',
              cpx.solution.get_objective_value())
        print('Best bound:                        ',
              cpx.solution.MIP.get_best_objective())
        print('Gap CPLEX:                          %.2f %%' %
              self.cpx.solution.MIP.get_mip_relative_gap())
        print('Solution time (total):              %.2f s' % self.totaltime)
        print('Cuts:                              ', self.cutscount)
        # print('Iteration:                         ', self.iterations)

        # self.plot()

    def get_Status(self) -> int:
        if (self.cpx.solution.get_status() == cplex._internal._constants.CPXMIP_OPTIMAL or
                self.cpx.solution.get_status() == cplex._internal._constants.CPXMIP_OPTIMAL_TOL):
            return 1
        if (self.cpx.solution.get_status() == cplex._internal._constants.CPXMIP_ABORT_FEAS or
                self.cpx.solution.get_status() == cplex._internal._constants.CPXMIP_TIME_LIM_FEAS):
            return 0
        return 0

    def __update_sets(self):
        # Sets
        setVars = list(range(self.Nitem))
        setS = [list(range(self.Nb[j]-1)) for j in setVars]

        setSConcave = [
            (j, s)
            for j in setVars for s in setS[j]
            if int(self.Concave[j][s]) == 1
        ]

        setSConvex = [
            (j, s)
            for j in setVars for s in setS[j]
            if int(self.Concave[j][s]) == -1
        ]

        setSComplet = [
            (j, s)
            for j in setVars for s in setS[j]
        ]

        self.setVars = setVars
        self.setS = setS
        self.setSConcave = setSConcave
        self.setSConvex = setSConvex
        self.setSComplet = setSComplet

    def initialSolutionFase(self):
        # Add points from global optimum
        for j in self.setVars:
            for s in self.setS[j]:
                self.Concave[j][s] = 1

        for b in range(20):
            valuesX = [b for j in self.setVars]
            valuesY = [self.NLfunc[j](b) for j in self.setVars]

            self.__refinement(valuesX, valuesY, -1, 1)
            self.__update_sets()

    def refinementFase(self):
        self.PC_REF_X_LIMITS = 1E-9
        self.PC_REF_X_INTVL_LIM = 1E-9
        self.PC_REF_Y_VIOLATION = 1E-6

        print("\n*** Refinement *** ")
        added = 0

        if (self.maxIterations > -1):
            if (self.maxIterations <= self.iterations):
                return added

        if (self.timelimit > -1):
            if (self.timelimit <= self.totaltime):
                return added

        REF_TOL = 0

        # Add points from global optimum
        valuesX = [self.varValGlobalUB[self.varX[j]] +
                   REF_TOL for j in self.setVars]
        valuesY = [self.varValGlobalUB[self.varY[j]] for j in self.setVars]

        print("From IPOPT")
        added += self.__refinement(valuesX, valuesY, -1, 1)
        self.__update_sets()

        # valuesX = [self.varValGlobalUB[self.varX[j]] - REF_TOL for j in self.setVars ]
        # valuesY = [self.varValGlobalUB[self.varY[j]] for j in self.setVars ]

        # added += self.__refinement(valuesX, valuesY, -1, 1)
        # self.__update_sets()

        # Add points from lower bounding
        valuesX = self.cpx.solution.get_values(
            [self.varX[j] for j in self.setVars])
        valuesY = self.cpx.solution.get_values(
            [self.varY[j] for j in self.setVars])

        # valuesX = [valuesX[j] + REF_TOL for j in self.setVars]

        print("From relaxation problem")
        added += self.__refinement(valuesX, valuesY,
                                   self.PC_REF_Y_VIOLATION, 0)
        self.__update_sets()

        # valuesX = self.cpx.solution.get_values([self.varX[j] for j in self.setVars])
        # valuesY = self.cpx.solution.get_values([self.varY[j] for j in self.setVars])

        # valuesX = [valuesX[j] - REF_TOL for j in self.setVars]

        # added += self.__refinement(valuesX, valuesY, self.PC_THRSHLD, 0)
        # self.__update_sets()

        self.iterations += 1
        self.totalBreakpointsAdded += added

        return added

    def __refinement(self, valuesX, valuesY, tol, mode=0):
        added = 0

        setVars = self.setVars

        lindex_add = []
        lvaluex_add = []

        for j in setVars:
            valueX = valuesX[j]

            # First test : Respect the limits
            if valueX < self.Lparam[j][0] + self.PC_REF_X_LIMITS or valueX > self.Lparam[j][-1] - self.PC_REF_X_LIMITS:
                continue

            valueYreal = -self.NLfunc[j](valueX)

            # # if 1 : <=
            # # if 2 : >=
            # Second test: Function Violation test
            if (mode == 0):
                if (self.ineq[j] == 1):
                    if (valueYreal - valuesY[j] > tol):
                        continue
                elif (valuesY[j] - valueYreal > tol):
                    continue

            # if (abs(valueYreal - valuesY[j]) > tol):
            setS_new = [s for (jn, s) in self.setSConcave
                        if jn == j and
                        valueX > self.Lparam[j][s] and
                        valueX < self.Lparam[j][s+1]]

            if (len(setS_new) == 1):
                s = setS_new[0]
                if (valueX > self.Lparam[j][s] + self.PC_REF_X_INTVL_LIM
                        and valueX < self.Lparam[j][s+1] - self.PC_REF_X_INTVL_LIM):
                    lindex_add.append((j, s))
                    lvaluex_add.append(valueX)

        # for i in range(len(lvaluex_add)):
        for i in range(len(lvaluex_add)-1, -1, -1):
            (j, s) = lindex_add[i]
            valuex = lvaluex_add[i]

            valueF = self.Lparam[j][s]
            valueS = self.Lparam[j][s+1]

            # if (abs(valuex - valueF) > self.PC_REF_THRSHLD and abs(valueS - valuex) > self.PC_REF_THRSHLD) :
            if ((valuex - valueF) > self.PC_REF_X_INTVL_LIM and
                    (valueS - valuex) > self.PC_REF_X_INTVL_LIM):
                self.Lparam[j].insert(s+1, valuex)
                self.Concave[j].insert(s, 1)
                self.Nb[j] += 1
                added += 1

        print("New breakpoints:                    " + str(added))

        return added

    def refinement(self):
        added = 0
        if (self.maxIterations > -1):
            if (self.maxIterations <= self.iterations):
                return added

        if (self.timelimit > -1):
            if (self.timelimit <= self.totaltime):
                return added

        lindex_add = []
        lvaluex_add = []
        for (j, s) in self.setSConcave:
            valuey = self.cpx.solution.get_values(self.y_ref[j][s])
            if (valuey > 0.5):
                valuex = self.cpx.solution.get_values(
                    self.x_ref[j][s]) + self.lParamMode[j][s]
                lindex_add.append((j, s))
                lvaluex_add.append(valuex - self.PC_REF_DISTPOINTS)
                if (self.PC_REF_DISTPOINTS > 0):
                    lindex_add.append((j, s))
                    lvaluex_add.append(valuex + self.PC_REF_DISTPOINTS)

        # for i in range(len(lvaluex_add)):
        for i in range(len(lvaluex_add)-1, -1, -1):
            (j, s) = lindex_add[i]
            valuex = lvaluex_add[i]

            valueF = self.Lparam[j][s]
            valueS = self.Lparam[j][s+1]

            # if (abs(valuex - valueF) > self.PC_REF_THRSHLD and abs(valueS - valuex) > self.PC_REF_THRSHLD) :
            if ((valuex - valueF) > self.PC_REF_THRSHLD and (valueS - valuex) > self.PC_REF_THRSHLD):
                self.Lparam[j].insert(s+1, valuex)
                self.Concave[j].insert(s, 1)
                self.Nb[j] += 1
                added += 1

        print("\n*** Refinement *** ")
        print("New breakpoints:                    " + str(added))
        print('')

        self.iterations += 1
        self.totalBreakpointsAdded += added
        return added

    def setIPOPT(self, ipopt):
        self.ipopt = ipopt

    def runIPOPT(self):
        varValues = self.cpx.solution.get_values()

        # Update solution with non-relaxed values
        for j in self.setVars:
            valueX = varValues[self.varX[j]]
            # valueY = varValues[self.varY[j]]

            valueZ = [varValues[self.y_ref[j][i]]
                      for i in range(len(self.y_ref[j]))]
            valueZ = sum(valueZ)

            if (valueZ > 0.1):
                valueY2 = -self.NLfunc[j](valueX)
            else:
                valueY2 = 0.0

            varValues[self.varY[j]] = valueY2

        # Mapping
        solutionList = [varValues[i] for i in self.newVarIndex]
        indexList = list(range(len(solutionList)))
        self.ipopt.resetValues(indexList, solutionList)

        cpxCopy = self.cpxCopy
        j = 0
        integetVar = cpxCopy.variables.get_num_binary(
        ) or cpxCopy.variables.get_num_integer()
        for i in self.newVarIndex:
            if integetVar and cpxCopy.variables.get_types(i) != 'C':
                self.ipopt.fixVariable(j, solutionList[j])
            j += 1

        self.ipopt.optimize()
        # print(self.ipopt.messageSolution)

        if self.ipopt.optimalFound():

            currentFixVarSolution_temp = self.ipopt.getSolutionVar()
            currentFixSolution = self.ipopt.getSolutionObj()

            # REVIEW
            if cpxCopy.objective.sense[cpxCopy.objective.get_sense()] == 'maximize':
                currentFixSolution = - currentFixSolution

            currentFixVarSolution = [0]*len(currentFixVarSolution_temp)
            for i in range(len(self.newVarIndex)):
                j = self.newVarIndex[i]
                currentFixVarSolution[j] = currentFixVarSolution_temp[i]

            print(
                f'Ipopt solution:                     {currentFixSolution}\n')

            if cpxCopy.objective.sense[cpxCopy.objective.get_sense()] == 'minimize':
                if self.objValGlobalUB > currentFixSolution:
                    self.objValGlobalUB = currentFixSolution
                    self.varValGlobalUB = currentFixVarSolution

                    # self.refinement2(1)

            elif cpxCopy.objective.sense[cpxCopy.objective.get_sense()] == 'maximize':
                if self.objValGlobalUB < currentFixSolution:
                    self.objValGlobalUB = currentFixSolution
                    self.varValGlobalUB = currentFixVarSolution

                    # self.refinement2(1)
        else:
            print("IPOPT: " + self.ipopt.messageStatus + "\n")

    def fixingLocal(self):
        print('*** Upper Bounding problem ***')
        print("Checking incumbent solution...")
        cpxCopy = self.cpxCopy

        cpxCopy.parameters.simplex.tolerances.optimality.set(1e-8)
        cpxCopy.parameters.simplex.tolerances.feasibility.set(1e-8)
        cpxCopy.parameters.paramdisplay.set(0)
        cpxCopy.parameters.mip.display.set(0)
        cpxCopy.set_log_stream(None)
        cpxCopy.set_error_stream(None)
        cpxCopy.set_warning_stream(None)
        cpxCopy.set_results_stream(None)

        try:
            varValues = self.cpx.solution.get_values()

            if cpxCopy.get_problem_type() != 0:
                for i in range(cpxCopy.variables.get_num()):
                    if cpxCopy.variables.get_types(i) != 'C':
                        # valueInt = varValues[i]
                        # cpxCopy.variables.set_types(i, cpxCopy.variables.type.continuous)
                        cpxCopy.variables.set_lower_bounds(i, varValues[i])
                        cpxCopy.variables.set_upper_bounds(i, varValues[i])

            for j in self.setVars:
                valueX = varValues[self.varX[j]]
                valueY = varValues[self.varY[j]]
                valueZ = [varValues[self.y_ref[j][i]]
                          for i in range(len(self.y_ref[j]))]
                valueZ = sum(valueZ)

                if (valueZ > 0.1):
                    valueY2 = -self.NLfunc[j](valueX)
                else:
                    valueY2 = 0.0

                cpxCopy.variables.set_lower_bounds(self.varX[j], valueX)
                cpxCopy.variables.set_upper_bounds(self.varX[j], valueX)

                if self.funcSense[j] == 'L':
                    if cpxCopy.variables.get_lower_bounds(self.varY[j]) < valueY2:
                        cpxCopy.variables.set_lower_bounds(
                            self.varY[j], valueY2)
                else:
                    if cpxCopy.variables.get_upper_bounds(self.varY[j]) > valueY2:
                        cpxCopy.variables.set_upper_bounds(
                            self.varY[j], valueY2)

            cpxCopy.solve()

            currentBestBound = self.cpx.solution.MIP.get_best_objective()
            currentFixSolution = cpxCopy.solution.get_objective_value()

            print(f'Incumbent solution:                 {currentFixSolution}')

            if self.Problem_Sense == 'minimize':
                if self.objValGlobalLB < currentBestBound:
                    self.objValGlobalLB = currentBestBound

                if self.objValGlobalUB > currentFixSolution:
                    self.objValGlobalUB = currentFixSolution
                    self.varValGlobalUB = cpxCopy.solution.get_values()

            elif self.Problem_Sense == 'maximize':
                if self.objValGlobalLB > currentBestBound:
                    self.objValGlobalLB = currentBestBound

                if self.objValGlobalUB < currentFixSolution:
                    self.objValGlobalUB = currentFixSolution
                    self.varValGlobalUB = cpxCopy.solution.get_values()

        except cplex.exceptions.errors.CplexSolverError as err:
            print("Not possible to find incubent solution.")

        self.runIPOPT()

        print('*** Current solution ***')
        if self.Problem_Sense == 'minimize':
            print(f'Upper bounding problem solution:    {self.objValGlobalUB}')
            print(f'Lower bounding problem solution:    {self.objValGlobalLB}')
        elif self.Problem_Sense == 'maximize':
            print(f'Lower bounding problem solution:    {self.objValGlobalUB}')
            print(f'Upper bounding problem solution:    {self.objValGlobalLB}')

        print(
            f'Gap:                                {self.getGlobalGap()*100:.3f} %\n')

    def getGlobalGap(self):
        # |bestbound-bestinteger|/(1e-10+|bestinteger|)
        return abs(self.objValGlobalLB - self.objValGlobalUB)/(1e-10+abs(self.objValGlobalUB))

    def PrintSolution(self):
        for i in range(self.cpx.variables.get_num()):
            name = self.cpx.variables.get_names(i)
            value = self.cpx.solution.get_values(i)
            print(name + '\t' + str(value))

    def addConstraintFromPool(self):
        self.cpx.linear_constraints.add(
            lin_expr=self.poolconst,
            senses=self.poolsense,
            rhs=self.poolrhs)

    # Defines the callback function
    def invoke(self, context):
        if context.in_candidate():
            print(
                "Custom set callback <================================================== Candidate")

            vars = context.get_candidate_point(self.listcplexvarnum)

            for j, s in self.setSConvex:
                varX = vars[self.x_ref[j][s]]
                varY = vars[self.y_ref[j][s]]
                VarZ = vars[self.nz_ref[(j, s)]]

                if self.p_mode == 1:
                    lParam = self.Lparam[j][s]
                else:
                    if self.p_mode == 2:
                        lParam = 0

                # self.local_lb = context.get_local_lower_bounds()
                # self.local_ub = context.get_local_upper_bounds()

                if varY > 0:
                    xbar = varX/varY

                    # Check increase
                    increase = false
                    part1 = varY * \
                        self.NLfunc[j](xbar+lParam) - self.NLfunc[j](0+lParam)
                    increase = part1 - VarZ
                    if increase >= self.PC_THRSHLD:
                        increase = true
                    else:
                        increase = false

                    # Check Loop
                    prevxbar = self.prevxbar[(j, s)]
                    checkloop = true
                    if (prevxbar > self.PC_THRSHLD):
                        ratio = (xbar - prevxbar)/prevxbar
                        if ((ratio < self.PC_EPSCHCKLOOP) and (ratio > -self.PC_EPSCHCKLOOP)):
                            checkloop = false
                    else:
                        diff = xbar - prevxbar
                        if ((diff < self.PC_EPSCHCKLOOP) and (diff > -self.PC_EPSCHCKLOOP)):
                            checkloop = false

                    # Calculating violation
                    violation = varY * \
                        self.NLfunc[j](xbar+lParam) - varY * \
                        self.NLfunc[j](0+lParam) - VarZ

                    if (violation > self.PC_THRSHLD and checkloop and increase):
                        self.prevxbar[(j, s)] = xbar

                        xCoeff = -self.NLfuncDiff[j](xbar+lParam)
                        yCoeff = -(self.NLfunc[j](xbar+lParam) - self.NLfunc[j]
                                   (0+lParam) - self.NLfuncDiff[j](xbar+lParam)*xbar)
                        # yCoeff = 0

                        newconst = cplex.SparsePair([self.x_ref[j][s]] + [self.y_ref[j][s]] + [self.nz_ref[(j, s)]],
                                                    [xCoeff]+[yCoeff]+[1])
                        newsense = 'G'
                        newrhs = 0
                        # newrhs = (self.NLfunc[j](xbar+lParam) - self.NLfunc[j](0+lParam) - self.NLfuncDiff[j](xbar+lParam)*xbar)

                        context.reject_candidate(
                            constraints=[newconst],
                            senses=[newsense],
                            rhs=[newrhs])

                        # self.poolconst.append(newconst)
                        # self.poolsense.append(newsense)
                        # self.poolrhs.append(newrhs)

        if (context.in_relaxation() and self.relaxcutcount < self.PC_RELAXCUTLIM):
            # print ("Custom set callback <================================================== Relaxation")
            vars = context.get_relaxation_point(self.listcplexvarnum)

            for j, s in self.setSConvex:
                varX = vars[self.x_ref[j][s]]
                varY = vars[self.y_ref[j][s]]
                VarZ = vars[self.nz_ref[(j, s)]]

                Ineq = self.ineq[j]

                if self.p_mode == 1:
                    lParam = self.Lparam[j][s]
                else:
                    if self.p_mode == 2:
                        lParam = 0

                # self.local_lb = context.get_local_lower_bounds()
                # self.local_ub = context.get_local_upper_bounds()

                if varY > 0:
                    xbar = varX/varY

                    # Check increase
                    increase = false
                    part1 = varY * \
                        self.NLfunc[j](xbar+lParam) - self.NLfunc[j](0+lParam)
                    inc = part1 - VarZ
                    if inc >= self.PC_THRSHLD:
                        increase = true
                    else:
                        increase = false

                    # Check Loop
                    prevxbar = self.prevxbar[(j, s)]
                    checkloop = true
                    if (prevxbar > self.PC_THRSHLD):
                        ratio = (xbar - prevxbar)/prevxbar
                        if ((ratio < self.PC_EPSCHCKLOOP) and (ratio > -self.PC_EPSCHCKLOOP)):
                            checkloop = false
                    else:
                        diff = xbar - prevxbar
                        if ((diff < self.PC_EPSCHCKLOOP) and (diff > -self.PC_EPSCHCKLOOP)):
                            checkloop = false

                    # Calculating violation
                    violation = varY * \
                        self.NLfunc[j](xbar+lParam) - varY * \
                        self.NLfunc[j](0+lParam) - VarZ

                    if (violation > self.PC_THRSHLD and checkloop and increase):
                        self.prevxbar[(j, s)] = xbar

                        xCoeff = -self.NLfuncDiff[j](xbar+lParam)
                        yCoeff = -(self.NLfunc[j](xbar+lParam) - self.NLfunc[j]
                                   (0+lParam) - self.NLfuncDiff[j](xbar+lParam)*xbar)

                        cutmanagement = cplex.callbacks.UserCutCallback.use_cut.force

                        newconst = cplex.SparsePair([self.x_ref[j][s]] + [self.y_ref[j][s]] + [self.nz_ref[(j, s)]],
                                                    [xCoeff]+[yCoeff]+[1])
                        if (Ineq == 1):
                            newsense = 'G'
                        else:
                            newsense = 'L'
                        newrhs = 0

                        self.relaxcutcount += 1
                        context.add_user_cut(
                            cut=newconst,
                            sense=newsense,
                            rhs=newrhs,
                            cutmanagement=cutmanagement,
                            local=False)

                        # self.poolconst.append(newconst)
                        # self.poolsense.append(newsense)
                        # self.poolrhs.append(newrhs)

        if context.in_branching():
            print(
                "Custom set callback <================================================== Branching")
            # context.prune_current_node()
            # context.abort()
        if context.in_local_progress():
            print(
                "Custom set callback <================================================== Local Progress")
            # context.exit_cut_loop()
            # context.prune_current_node()
        if context.in_global_progress():
            print(
                "Custom set callback <================================================== Global Progress")


class NLProblem_Callback_Incumbent(cplex.callbacks.IncumbentCallback):

    # NLProblemData = (NLProblem)

    def __call__(self):
        print("Custom reg callback <================================================== Incumbent")

        # obj = self.get_best_objective_value()


class NLProblem_Callback_Heuristic(cplex.callbacks.HeuristicCallback):
    def __call__(self):
        print("Custom reg callback <================================================== Heuristic")
        print("\t\t   " + str(self.get_best_objective_value()))

        pass


class NLProblem_Callback_Continuous(cplex.callbacks.ContinuousCallback):
    def __call__(self):
        print("Custom reg callback <================================================== Continuous")
        pass


class NLProblem_Callback_SolveCallback(cplex.callbacks.SolveCallback):
    def __call__(self):
        print(
            "Custom reg callback <================================================== Solve")
        print("\t\t   " + str(self.get_best_objective_value()))


class NLProblem_Callback_BranchCallback(cplex.callbacks.BranchCallback):
    def __call__(self):
        print(
            "Custom reg callback <================================================== Branch")
        # self.prune()


class NLProblem_Callback_NodeCallback(cplex.callbacks.NodeCallback):
    def __call__(self):
        print(
            "Custom reg callback <================================================== Node")


class NLProblem_Callback_LazyConstraint(cplex.callbacks.LazyConstraintCallback):
    'Defines the CPLEX Callback for Lazy Constraints call.'

    def __call__(self):
        # print ("Custom reg callback <================================================== Lazy Const.")
        CutsSeparation(self, self.use_constraint.force)


class NLProblem_Callback_UserCut(cplex.callbacks.UserCutCallback):

    def __call__(self):
        # print ("Custom reg callback <================================================== User cut " )
        CutsSeparation(self, self.use_cut.force)


def CutsSeparation(self: cplex.callbacks.LazyConstraintCallback, useMode=0):
    self2: NLProblem = self.NLProblemData
    # self2 = self.NLProblemData
    # setSConvex = self2.setS

    vars = self.get_values()

    # startCutscount = self2.cutscount

    for j, s in self2.setSConvex:
        varX = vars[self2.x_ref[j][s]]
        varY = vars[self2.y_ref[j][s]]
        VarZ = vars[self2.nz_ref[(j, s)]]

        if (self2.Lparam[j][0] > 0):
            if abs(varX - 0) < self2.PC_PARAM_EPSILON:
                continue
        elif abs(varX - self2.Lparam[j][0]) < self2.PC_PARAM_EPSILON:
            continue

        # Mode: 1 Incremental, 2 Multiple choice
        lParam = self2.lParamMode[j][s]

        Ineq = self2.ineq[j]

        if varY > self2.PC_PARAM_EPSILON:  # 9.9999999999999995E-7
            xbar = varX/varY

            NLfuncX = self2.NLfunc[j](xbar+lParam)
            # NLfunc0 = self2.NLfunc[j](lParam)
            NLfunc0 = self2.lParamFunc[j][s]

            # NLfuncX = self2.myNLFunction(j,xbar+lParam)
            # NLfunc0 = self2.myNLFunction(j,lParam)

            # Check increase
            part1 = varY * NLfuncX - varY * NLfunc0
            # part1 = NLfuncX - NLfunc0

            if (Ineq == 1):
                inc = part1 - VarZ
            else:
                inc = -part1 + VarZ
            if (inc < self2.PC_THRSHLD):  # 9.9999999999999995E-7
                continue

            # Check Loop
            prevxbar = self2.prevSolVal[self2.x_ref[j][s]]
            if prevxbar > -1:
                if prevxbar >= self2.PC_PARAM_EPSILON:
                    ratio = abs((xbar - prevxbar)/prevxbar)
                    if ((ratio < self2.PC_EPSCHCKLOOP)):
                        continue
                else:
                    diff = abs(xbar - prevxbar)
                    if ((diff < self2.PC_EPSCHCKLOOP)):
                        continue

            self2.prevSolVal[self2.x_ref[j][s]] = xbar

            # xCoeff = -self2.NLfuncDiff[j](xbar+lParam)
            # yCoeff = -(self2.NLfunc[j](xbar+lParam) - self2.NLfunc[j](lParam) - self2.NLfuncDiff[j](xbar+lParam)*xbar)

            NLDiffX = self2.NLfuncDiff[j](xbar+lParam)
            # NLDiffX = self2.myNLFunctionDiff(j,xbar+lParam)

            xCoeff = -NLDiffX
            yCoeff = -(NLfuncX - NLfunc0 - NLDiffX*xbar)
            # yCoeff = 0

            newconst = cplex.SparsePair([self2.x_ref[j][s]] + [self2.y_ref[j][s]] + [self2.nz_ref[(j, s)]],
                                        [xCoeff]+[yCoeff]+[1])

            # self2.relaxcutcount += 1
            self2.cutscount += 1
            newsense = self2.funcSense[j]
            newrhs = 0.0
            # newrhs = (NLfuncX - NLfunc0 - NLDiffX*xbar)

            self.add(
                newconst,
                sense=newsense,
                rhs=newrhs,
                use=useMode
            )

            # # Print cut
            # print ( "#####" + str(-self2.NLfunc[j](0)) + " + " + str(xCoeff) + "*x + " + str(yCoeff) + " " + newsense + "= " + str(newrhs))
            # print (str(lParam) + ", " + str(lParam) )
            # print("ncuts=",str(self2.cutscount-startCutscount))


def test_IntegerMPS_Relax():
    mpsNameFile = "tmp/prob.mps"
    jsonNameFile = 'tmp/data.json'
    reformulationType = 2

    if len(argv) >= 4:
        mpsNameFile = argv[1]
        jsonNameFile = argv[2]
        reformulationType = int(argv[3])

    prob = NLProblem()
    prob.LoadFileJSON(jsonNameFile)
    prob.CreateProblemMPS(mpsNameFile, reformulationType)
    prob.RelaxIntegrality()
    # prob.ConfigParametersTotal()
    prob.Solve(0)
    prob.ExportSol()
    # prob.IpoptCall()
    # prob.ExportSolIpopt()
    # prob.IpoptCall("wb.nl", "wb3.nl")


def test_IntegerMPS():
    mpsNameFile = "tmp/prob.lp"
    jsonNameFile = 'tmp/dataNew.json'
    reformulationType = 2

    if len(argv) >= 4:
        mpsNameFile = argv[1]
        jsonNameFile = argv[2]
        reformulationType = int(argv[3])

    prob = NLProblem()
    prob.LoadFileJSON(jsonNameFile)
    prob.CreateProblemMPS(mpsNameFile, reformulationType)
    # prob.RelaxIntegrality()
    # prob.ConfigParametersTotal()
    prob.Solve(2)
    prob.ExportSol()

    # prob.IpoptCall()
    # prob.ExportSolIpopt()

    # prob.IpoptCall("wb.nl", "wb3.nl")


def test_Integer():
    prob = NLProblem()
    prob.LoadFileJSON('data.json')
    prob.CreateProblem(1)
    # prob.ConfigParametersTotal()
    prob.Solve(0)


def test_Relax():
    prob = NLProblem()
    prob.LoadFileJSON('data.json')
    prob.CreateProblem(2)
    prob.RelaxIntegrality()
    prob.ConfigParametersTotal()
    prob.Solve(0)


def main():
    if len(argv) >= 5:
        if (int(argv[4]) == 1):
            test_IntegerMPS_Relax()
        if (int(argv[4]) == 0):
            test_IntegerMPS()
    else:
        test_IntegerMPS()
    # test_Integer()
    # test_Relax()

    # prob = NLProblem()
    # prob.CreateProblem(1)
    # prob.RelaxIntegrality()
    # prob.PrintProblem()
    # prob.ConfigParametersTotal()
    # prob.Solve(0)
    # prob.ExportFinalLP()
    # prob.PrintSolution()


if __name__ == "__main__":
    main()
