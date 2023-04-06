#!/usr/local/bin/python3

import sympy as syp
from mpmath import findroot
from math import *
import types
import time
from sys import argv

class BPcalc(object):
    def __init__(self):
        self.samplen = 1000

        # Import
        self.counterlines = 0
        self.indexx = []
        self.indexy = []
        self.varxlb = []
        self.varxub = []
        self.consrhslb = []
        self.consrhsub = []
        self.nlfunc = []

        # Export
        self.varxname = []
        self.varxid = []
        self.varxub = []
        self.varxlb = []
        self.ineq = []

        self.varyname = []
        self.varyid = []

        self.varzid = []

        self.functions1 = []
        self.functions2 = []

        self.setBreakpoints = []
        self.setEstimation = []

        self.execute_time = 0

    def exportJSON(self, outNameFile='tmp/dataNew.json') : 
        import json
        exportJson = {}

        exportJson["varx_name"] = self.varxname
        exportJson["varx_id"] = self.varxid
        exportJson["varx_ub"] = self.varxub
        exportJson["varx_lb"] = self.varxlb


        exportJson["vary_name"] = self.varyname
        exportJson["vary_id"] = self.varyid

        exportJson["varz_id"] = self.varzid

        exportJson["functions"] = self.functions1
        exportJson["funcdiff2"] = self.functions2

        exportJson["breakpoints"] = self.setBreakpoints
        exportJson["estimation"] = self.setEstimation

        exportJson["ineq"] = self.ineq

        with open(outNameFile, 'w') as f:
            json.dump(exportJson, f)
    
        print("File exported to: " + outNameFile)

    def set_asl(self, asl) :
        self.counterlines = asl.get_nl_count()
        self.indexx = asl.get_nl_varx(range(self.counterlines))
        self.indexy = asl.get_nl_vary(range(self.counterlines))
        self.indexz = asl.get_nl_varz(range(self.counterlines))
        varxlb = asl.get_nl_varxlb(range(self.counterlines))
        varxub = asl.get_nl_varxub(range(self.counterlines))
        self.consrhslb = asl.get_nl_consub(range(self.counterlines))
        self.consrhsub = asl.get_nl_conslb(range(self.counterlines))
        self.nlfunc = asl.get_nl_functionslist(range(self.counterlines))

        x, y, z = syp.symbols("x,y,z")
        for i in range(0, self.counterlines):
            idx = self.indexx[i]
            idy = self.indexy[i]
            idz = self.indexz[i]

            varxnam = "x"+str(idx+1)
            varynam = "x"+str(idy+1)

            # x0, y0 = symbols(varxnam + "," + varynam)
            fun = syp.sympify(self.nlfunc[i])
            # fun = fun.replace(x0, x)
            fun = fun.replace(y, 0)
            fun = fun.replace(z, 1)

            if (self.consrhslb[i] < inf) :
                self.functions1.append(str(  fun-self.consrhslb[i]  ) )
                self.ineq.append(1)
                self.varxid.append(idx+1)
                self.varyid.append(idy+1)
                self.varzid.append(idz+1)
                self.varxname.append(varxnam)
                self.varyname.append(varynam)
                self.varxlb.append(varxlb[idx])
                self.varxub.append(varxub[idx])

            if (self.consrhsub[i] > -inf) :
                self.functions1.append(str( (fun-self.consrhsub[i]) ) )
                self.ineq.append(-1)
                self.varxid.append(idx+1)
                self.varyid.append(idy+1)
                self.varzid.append(idz+1)
                self.varxname.append(varxnam)
                self.varyname.append(varynam)
                self.varxlb.append(varxlb[idx])
                self.varxub.append(varxub[idx])

    def execute(self):
        start_time = time.time()
        x = syp.symbols("x")
        for j in range(len(self.functions1)):
            lb = self.varxlb[j]
            ub = self.varxub[j]


            a = syp.sympify("-(" + self.functions1[j] + ")")  
            b = syp.diff(a, x)
            c = syp.diff(b, x)

            new_func = 0
            funcDefStr = "def myNLFunction" + str(0) + "(x):\n\ty = " + str(a) + "\n\treturn y"
            code_obj = compile(funcDefStr, '<string>', 'exec')
            new_func = types.FunctionType(code_obj.co_consts[0], globals())

            new_func_diff2 = 0
            funcDefStr = "def myNLFunction" + str(0) + "(x):\n\ty = " + str(c) + "\n\treturn y"
            code_obj = compile(funcDefStr, '<string>', 'exec')
            new_func_diff2 = types.FunctionType(code_obj.co_consts[0], globals())

            self.functions2.append(str(b))

            breakpoints = []
            breakpoints.append(lb)

            # Sample the interval in $self.samplen points
            for i in range(self.samplen):
                q = lb + i*((ub-lb)/self.samplen)

                extreme1 = new_func_diff2(q)
                extreme2 = new_func_diff2(q + ((ub-lb)/self.samplen))

                signal1 = (extreme1 > 0)
                signal2 = (extreme2 > 0)
                if (signal1 != signal2) :
                    # r = nsolve(c,x,(q, q + ((ub-lb)/self.samplen) ), verbose=false, solver='anderson')
                    # r = nsolve(c,x,(q, q + ((ub-lb)/self.samplen) ), solver='anderson')
                    # r = findroot(lambda x:-0.0462637706301063*sin(0.0785398163397448*x - 0.785398163397448) + 0.0231318853150531*cos(0.0392699081698724*x - 0.392699081698724) , (q, q + ((ub-lb)/100) ) , verify=False , verbose=False, solver='anderson')
                    r = findroot(lambda x:new_func_diff2(x) , (q, q + ((ub-lb)/self.samplen) ) , verify=False , verbose=False, solver='anderson')

                    # print("Zero found: " + str('% 4.18f'%c.evalf(subs={x: r})) + " [" + str('% 4.10f'%q) + ", " + str('% 4.10f'%(q + ((ub-lb)/self.samplen))) + "] = " + str('% 4.16f'% r) )
                    # print("Zero found: " + str('% 4.18f'%new_func_diff2(r) ) + " [" + str('% 4.10f'%q) + ", " + str('% 4.10f'%(q + ((ub-lb)/self.samplen))) + "] = " + str('% 4.16f'% r) )

                    if abs(lb - float(r)) >= 1e-10 :
                        breakpoints.append(float(r))
                    else :
                        print("Different")

            breakpoints.append(ub)
            self.setBreakpoints.append(breakpoints)

            estimation = []
            # Concave and convex information
            for b in range(0,len(breakpoints)-1) :
                valuex1 = breakpoints[b]
                valuex2 = breakpoints[b+1]

                valuexm = breakpoints[b] + (breakpoints[b+1] - breakpoints[b])/2

                valuey1 = new_func(valuex1)
                valuey2 = new_func(valuex2)

                alpha = (valuey2 - valuey1)/(valuex2 - valuex1)

                valueym = new_func(valuexm)
                valuepwm = (valuexm)*(alpha) + (valuey1-(valuex1*alpha))

                if(self.ineq[j] == 1) :
                    if(valuepwm < valueym) : 
                        estimation.append(-1)
                    else :
                        estimation.append(1)
                else :
                    if(valuepwm < valueym) : 
                        estimation.append(1)
                    else :
                        estimation.append(-1)
            
            self.setEstimation.append(estimation)
        self.execute_time = time.time() - start_time

def main(argv):
    print ("Start")

    from wrapperc import ASLPY
    asl = ASLPY()
    if(len(argv) > 1) :
        asl.set_argv(argv)
    else :
        asl.set_argv_fileonly("/Users/renanspencertrindade/Desktop/zGenericSolver/tmp/prob2.nl")
        asl.start()

    bpcalc = BPcalc()
    bpcalc.set_asl(asl)
    bpcalc.execute()
    bpcalc.exportJSON()
    print("Execution time:" + str(bpcalc.execute_time) )
    print ("End")


if __name__ == '__main__':
    main(argv)