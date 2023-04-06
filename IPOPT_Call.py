import struct

def read_Char(f):
    byte = f.read(1)
    # byteT = bytes(byte).decode("utf-8") 
    # print (byteT)
    return byte

def read_Str(f):
    byte = f.read(4)
    byteT = int.from_bytes(byte,byteorder='little')
    byte2 = f.read(byteT)
    # byteT = bytes(byte2).decode("utf-8") 
    # print (byteT)
    return byte+byte2

def read_Int(f):
    byte = f.read(4)
    # byteT = int.from_bytes(byte,byteorder='little')
    # print (byteT)
    return byte

def read_Float(f):
    byte = f.read(8)
    # byteT = struct.unpack("d",byte)[0]
    # print (byteT)
    return byte

# suffix values
def read_Suffix(f):
    Kind = read_Int(f)

    string = b'S'+Kind

    Kind = int.from_bytes(Kind,byteorder='little')
    integer = Kind & 4 == 0
    Kind = Kind & 3

    numCons = read_Int(f)
    string += numCons

    string += read_Str(f)

    numCons = int.from_bytes(numCons,byteorder='little')
    for i in range(numCons) :
        string += read_Int(f)

        if (integer) :
            string += read_Int(f)
        else :
            string += read_Float(f)

    return string

# bounds on variable
def read_Bounds(f, numVar):
    list = Nl_VarBoundList()
    for i in range(numVar) :
        ref = read_Char(f)
        lb = ""
        ub = ""
        if ref == b'0' : # >= <=
            lb = read_Float(f)
            ub = read_Float(f)
        if ref == b'1' : # <=
            ub = read_Float(f)
        if ref == b'2' : # >=
            lb = read_Float(f)
        if ref == b'4' : # ==
            lb = read_Float(f)
            ub = lb

        list.add(ref, lb, ub)
    return list

# Primal initial guess
def read_Xvalues(f):
    numVar = read_Int(f)
    list = Nl_VarInitialValues(numVar)

    numVar = int.from_bytes(numVar,byteorder='little')

    for i in range(numVar) :
        index = read_Int(f)
        value = read_Float(f)
        list.add(index, value)

    return list   

class Nl_VarInitialValues :
    def __init__(self, numVar=b'\x00\x00\x00\x00') :
        self.numVar = numVar
        self.indexList = []
        self.valueList = []

    def add(self, index, value):
        self.indexList.append(index)
        self.valueList.append(value)

    def resetValues(self, indexList, valueList):
        self.indexList = [int.to_bytes(intval,length=4,byteorder='little') for intval in indexList]
        self.valueList = [struct.pack("d", floatval) for floatval in valueList]

    def __getnumvar(self) :
        numVar = len(self.indexList)
        return int.to_bytes(numVar,length=4,byteorder='little')

    def getStringB(self) : 
        if len(self.indexList) == 0:
            return b''

        # string = b'x'+self.numVar
        string = b'x'+self.__getnumvar()
        for i in range(len(self.indexList)) :
            string += self.indexList[i]+self.valueList[i]
        return string

class Nl_VarBoundList : 
    def __init__(self) :
        from typing import List
        self.list : List[Nl_VarBound] = list()

    def add(self, code, lb, ub ) :
        self.list.append(Nl_VarBound(code, lb, ub))

    def printStr(self) : 
        print("b")
        for bd in self.list :
            bd.printStr()

    def getStringB(self) : 
        string = b'b'
        for bd in self.list :
            string += bd.getStringB()
        return string

    def setValue(self, index, code, lb, ub) :
        code = str.encode(str(code))
        lb = struct.pack("d", lb)
        ub = struct.pack("d", ub)
        self.list[index].code = code
        self.list[index].lb = lb
        self.list[index].ub = ub


class Nl_VarBound :
    def __init__(self, code, lb, ub) :
        self.code = code
        self.lb = lb
        self.ub = ub
    
    def printStr(self) :
        code = bytes(self.code).decode("utf-8") 
        lb = self.lb
        ub = self.ub

        if lb != "" :
            lb = struct.unpack("d",lb)[0]
            if float.is_integer(lb) :
                lb = int(lb)
    
        if ub != "" :
            ub = struct.unpack("d",ub)[0]
            if float.is_integer(ub) :
                ub = int(ub)

        if code == '0' : # >= <=
            print (f'{code} {lb} {ub}')
        if code == '1' : # <=
            print (f'{code} {ub}')
        if code == '2' : # >=
            print (f'{code} {lb}')
        if code == '3' : # free
            print (f'{code}')
        if code == '4' : # ==
            print (f'{code} {lb}')

    def getStringB(self) :
        code = self.code
        lb = self.lb
        ub = self.ub

        if code == b'0' : # >= <=
            return code+lb+ub
        if code == b'1' : # <=
            return code+ub
        if code == b'2' : # >=
            return code+lb
        if code == b'3' : # free
            return code
        if code == b'4' : # ==
            return code+lb


class IPOPT_Call :
    def __init__(self) :
        self.headFile = b''
        self.suffixFile = b''
        self.boundFile = b''
        self.initialXFile = b''
        self.finalPartFile = b''

        self.boundsList : Nl_VarBoundList = Nl_VarBoundList()
        self.xValuesList : Nl_VarInitialValues = Nl_VarInitialValues()

        self.varValues = []
        self.objective = 0.0

        self.messageSolution = ""
        self.messageStatus = ""

    def getSolutionObj(self) :
        return self.objective

    def getSolutionVar(self) :
        return self.varValues

    def fixVariable(self, index, value) :
        self.boundsList.setValue(index, 4, value, value)

    def resetValues(self, indexList = [], valueList = []) :
        self.xValuesList.resetValues(indexList, valueList)

    def exportFileB(self, nameFileExport:str="problen_Created.nl") :
        self.boundFile = self.boundsList.getStringB()
        self.initialXFile = self.xValuesList.getStringB()

        with open(nameFileExport, "wb") as f:
            f.write(self.headFile+self.suffixFile+self.boundFile+self.initialXFile+self.finalPartFile)
        f.close()

    def readFileB(self, nameFileImport:str="problen.nl") :
        try:
            with open(nameFileImport, "rb") as f:
                # First Line
                line = f.readline()
                self.headFile = line
                StrLine = line[:-1].decode("utf-8")
                # print(StrLine)
                binary = False
                char = StrLine[0]
                if char == "b" :
                    binary = True

                # Second Line
                line = f.readline()
                self.headFile += line
                StrLine = line[:-1].decode("utf-8")
                numVar = int(StrLine.split()[0])

                for i in range(8) : 
                    line = f.readline()
                    self.headFile += line
                    # print(line[:-1].decode("utf-8") )

                ref = read_Char(f)

                while ref:
                    if (ref == b'S') : 
                        self.suffixFile += read_Suffix(f)
                    elif (ref == b'b') : 
                        self.boundsList = read_Bounds(f, numVar)
                        self.boundFile = self.boundsList.getStringB()
                        # self.boundsList.printStr()
                    elif (ref == b'x') : 
                        self.xValuesList = read_Xvalues(f)
                        self.initialXFile = self.xValuesList.getStringB()
                    else :
                        self.finalPartFile += ref
                        ref = 0
                        byte = f.read(1)
                        while byte:
                            self.finalPartFile += byte
                            byte = f.read(1)
                        continue

                    ref = read_Char(f)

        except IOError:
            print('Error While Opening the file! : ' + nameFileImport)  


    def optimalFound(self) :
        if self.messageStatus == "Optimal Solution Found." :
            return True
        return False

    def optimize(self) :
        import os

        filename_out = "tmp/ipopt.nl"

        self.exportFileB(filename_out)

        # stream = os.system("ipopt " + filename_out + " wantsol=10 print_level=3")
        stream = os.popen("ipopt " + filename_out + " wantsol=10 print_level=3")
        output = stream.read()
        firstDivision = output.split("EXIT:")
        # print(output)
        self.messageSolution = output

        # Converged to a point of local infeasibility. Problem may be infeasible.
        # Optimal Solution Found.
        self.messageStatus = firstDivision[1].split("\n")[0].lstrip()

        # Solution
        sol = firstDivision[1].split("value\n")[1].splitlines()

        self.varValues = []
        for i in range(len(sol)) :
            pairValues = sol[i].split()
            self.varValues.append(float(pairValues[1]))

        # # Solution info
        info = firstDivision[0].split("Objective...............:")[1].splitlines()
        self.objective = float(info[0].split()[1])
        # dualinf = float(info[1].split("Dual infeasibility......:")[1].split()[1])
        # # consviol = float(info[2].split("Constraint violation....:")[1].split()[1])
        # # varviol = float(info[3].split("Variable bound violation:")[1].split()[1])
        # # comp = float(info[4].split("Complementarity.........:")[1].split()[1])
        # NLPerror = float(info[5].split("Overall NLP error.......:")[1].split()[1])

        # Status = firstDivision[1].splitlines()[0].lstrip()

        # totaltime = float(firstDivision[0].split("Total seconds in IPOPT")[1].split()[1])

        # self.varsolipopt = x

        # print(Status)
        # print("Optimal value:                      " + str(obj))
        # print("Solution time:                      " + str(totaltime))
        # print()

        # print ("END IPOPT")
        # print ("******************************************************************************\n\n")


# nameFileImport = "prob 2.nl"
# ipopt = IPOPT_Call()
# ipopt.readFileB(nameFileImport)
# ipopt.exportFileB()

# ipopt.fixVariable(6,2)

# ipopt.optimize()

# print(ipopt.messageSolution)
# print(ipopt.getSolutionObj())


# f1 = open(nameFileImport, "rb")
# f2 = open("problen_Created.nl", "rb")

# f1str = f1.read()
# f2str = f2.read()

# if (f1str == f2str) :
#     print ("Files are identical")