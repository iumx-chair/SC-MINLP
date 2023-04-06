#!/usr/local/bin/python3

import ctypes
from ctypes import cdll, Structure, Union, c_double, c_int, c_void_p, c_bool ,POINTER
# from re import A

from sys import argv

class ograd(Structure):
    pass

class expr_n(Structure):
    pass

class ei(Union):
    pass

class expr(Structure):
    pass

ograd._fields_=[("coef",c_double),
               ("next",POINTER(ograd)),
               ("varno",c_int)
               ]

expr_n._fields_=[("op",POINTER(c_void_p)),
               ("v",c_double)
               ]

ei._fields_=[("e",POINTER(expr) ),
               ("ep",POINTER(POINTER(expr)) ),
               ("eif",POINTER(c_void_p) ),
               ("en",POINTER(expr_n) ),
               ("i",c_int),
               ("p",POINTER(c_void_p) ),
               ("d",POINTER(c_void_p) ),
               ("rp",POINTER(c_double)),
               ("D",POINTER(c_void_p) ),
               ("ce",POINTER(c_void_p) )
               ]

expr._fields_=[("op",POINTER(c_void_p)),
               ("a",c_int),
               ("dL",c_double),
               ("L",ei ),
               ("R",ei ),
               ("dR",c_double)
               ]

# lib = cdll.LoadLibrary('./ASL_Lib/Rose_Part/libaslpyc.so')
# lib = cdll.LoadLibrary('/Applications/amplide.macosx64/libaslpyc.so')

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "aslpy/lib"))

class ASLPY(object):
   def __init__(self):
      import os
      dll_name = "aslpy/lib/libaslpyc.so"
      dllabspath = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + dll_name
      self.lib = cdll.LoadLibrary(dllabspath)

      # Return types
      self.lib.ASLPYCgetX0.restype = ctypes.c_double
      self.lib.ASLPYCgetnvar.restype = ctypes.c_int
      self.lib.ASLPYCgetncon.restype = ctypes.c_int
      self.lib.ASLPYCgetnl.restype = ctypes.c_char_p
      self.lib.ASLPYCgetlp.restype = ctypes.c_char_p
      self.lib.ASLPYCgetvarlb.restype = ctypes.c_double
      self.lib.ASLPYCgetvarub.restype = ctypes.c_double
      self.lib.ASLPYCgetnlfunc.restype = ctypes.c_char_p
      self.lib.ASLPYCgetconsub.restype = ctypes.c_double
      self.lib.ASLPYCgetconslb.restype = ctypes.c_double
      self.lib.ASLPYCgetoptionchar.restype = ctypes.c_char_p

      self.lib.ASLPYCsetX0.argtypes = [ctypes.c_int, ctypes.c_double]

      self.lib.ASLPYgetexpr.argtypes = [ctypes.c_int]
      self.lib.ASLPYgetexpr.restype = POINTER(expr)

      self.lib.ASLPYgetexpr_e.argtypes = [POINTER(expr)]
      self.lib.ASLPYgetexpr_e.restype = ctypes.c_double

      self.lib.ASLPYgetOgrad.argtypes = [ctypes.c_int]
      self.lib.ASLPYgetOgrad.restype = POINTER(ograd)

      self.lib.ASLPYgetA_vals_.restype = POINTER(ctypes.c_double)
      self.lib.ASLPYgetA_rownos_.restype = POINTER(ctypes.c_int)
      self.lib.ASLPYgetA_colstarts_.restype = POINTER(ctypes.c_int)

      self.lib.ASLPYCgetLUv.restype = POINTER(ctypes.c_double)
      self.lib.ASLPYCgetUvx.restype = POINTER(ctypes.c_double)
      self.lib.ASLPYCgetUrhsx.restype = POINTER(ctypes.c_double)
      self.lib.ASLPYCgetLUrhs.restype = POINTER(ctypes.c_double)

      self.lib.ASLPYgetintvar.restype = POINTER(ctypes.c_bool)

      self.lib.ASLPYgetobjtype.restype = POINTER(ctypes.c_char)

      self.lib.ASLPYgetobjconst.argtypes = [ctypes.c_int]
      self.lib.ASLPYgetobjconst.restype = ctypes.c_double

      self.lib.ASLPYgetnlinearCons.restype = POINTER(ctypes.c_bool)

      # self.lib.ASLPYCwrite_solf = [ctypes.c_char_p]
      # self.obj = self.lib.ASLPY_new()

   def set_argv(self, argv):
      self.lib.ASLPYCsetargc(len(argv))
      for i in range(len(argv)):
         arg = str(argv[i]).encode('utf-8')
         self.lib.ASLPYCsetargv(i, arg )

   def set_argv_fileonly(self, fname=""):
      self.lib.ASLPYCsetargc(2)
      arg = str(fname).encode('utf-8')
      self.lib.ASLPYCsetargv(1, arg )

   def start(self):
      self.lib.ASLPYCstart()

   def get_option(self, i) :
      return self.lib.ASLPYCgetoptionint(i)

   def get_optionC(self, i) :
      return str(self.lib.ASLPYCgetoptionchar(i), 'utf-8')

   def write_sol(self):
      self.lib.ASLPYCwrite_sol()

   def write_sol(self, newsolution=[]):
      for i in range(self.get_n_var()) :
         self.lib.ASLPYCsetX0(i, newsolution[i])
      self.lib.ASLPYCwrite_sol()

   def write_solf(self, filename=""):
      file = filename.encode('utf-8')
      self.lib.ASLPYCwrite_solf(file)

   def get_n_var(self):
      return self.lib.ASLPYCgetnvar()
   
   def get_n_con(self):
      return self.lib.ASLPYCgetncon()

   def get_X0(self, indexes=[]):
      lst = []

      if(len(indexes) == 0):
         indexes = range(self.get_n_var())

      for i in indexes :
         lst.append(self.lib.ASLPYCgetX0(i))
      return lst

   def get_file_nl(self):
      return str(self.lib.ASLPYCgetnl(), 'utf-8')

   def get_file_lp(self):
      str1 = self.lib.ASLPYCgetlp()
      return str(str1, 'utf-8')

   def write_file_lp(self, filename):
      self.lib.ASLPYCwritelp(filename.encode('utf-8') )

   def get_nl_count(self):
      return self.lib.ASLPYCgetnlcount()

   def get_nl_varx(self, indexes=[]):
      lst = []

      if(len(indexes) == 0):
         indexes = range(self.get_nl_count())

      for i in indexes :
         lst.append(self.lib.ASLPYCgetnlvarx(i))
      return lst

   def get_nl_vary(self, indexes=[]):
      lst = []

      if(len(indexes) == 0):
         indexes = range(self.get_nl_count())

      for i in indexes :
         lst.append(self.lib.ASLPYCgetnlvary(i))
      return lst

   def get_nl_varz(self, indexes=[]):
      lst = []

      if(len(indexes) == 0):
         indexes = range(self.get_nl_count())

      for i in indexes :
         lst.append(self.lib.ASLPYCgetnlvarz(i))
      return lst

   def get_nl_varxlb(self, indexes=[]):
      lst = []

      if(len(indexes) == 0):
         indexes = self.get_nl_varx()     #TODO: Solve overcharge

      for i in indexes :
         lst.append(self.lib.ASLPYCgetvarlb(i))
      return lst

   def get_nl_varxub(self, indexes=[]):
      lst = []

      if(len(indexes) == 0):
         indexes = self.get_nl_varx()     #TODO: Solve overcharge

      for i in indexes :
         lst.append(self.lib.ASLPYCgetvarub(i))
      return lst

   def get_nl_functionslist(self, indexes=[]):
      lst = []

      if(len(indexes) == 0):
         indexes = range(self.get_nl_count())

      for i in indexes :
         lst.append(str(self.lib.ASLPYCgetnlfunc(i), 'utf-8') )
      return lst

   def get_nl_consub(self, indexes=[]):
      lst = []

      if(len(indexes) == 0):
         indexes = range(self.get_nl_count())

      for i in indexes :
         lst.append( self.lib.ASLPYCgetconsub(i) )
      return lst

   def get_nl_conslb(self, indexes=[]):
      lst = []

      if(len(indexes) == 0):
         indexes = range(self.get_nl_count())

      for i in indexes :
         lst.append( self.lib.ASLPYCgetconslb(i) )
      return lst
   
   def get_expr(self, index) :
      return self.lib.ASLPYgetexpr(index)
   
   def get_expr_n(self, expr) :
      return self.lib.ASLPYgetexpr_e(expr)
   
   def get_Ograd(self, index) :
      return self.lib.ASLPYgetOgrad(index)
   
   def get_LUv(self) :
      return self.lib.ASLPYCgetLUv()

   def get_Uvx(self) :
      return self.lib.ASLPYCgetUvx()
   
   def get_Urhsx(self) :
      return self.lib.ASLPYCgetUrhsx()
   
   def get_LUrhs(self) :
      return self.lib.ASLPYCgetLUrhs()
   
   def get_A_vals(self) :
      return self.lib.ASLPYgetA_vals_()
   
   def get_A_rownos(self) :
      return self.lib.ASLPYgetA_rownos_()
   
   def get_A_colstarts(self) :
      return self.lib.ASLPYgetA_colstarts_()
   
   def get_intvar(self) :
      return self.lib.ASLPYgetintvar()
   
   def get_objtype(self) :
      return self.lib.ASLPYgetobjtype()
   
   def get_objcons(self, index) :
      return self.lib.ASLPYgetobjconst(index)
   
   def get_nlinearCons(self) :
      return self.lib.ASLPYgetnlinearCons()


# asl = ASLPY()
# if(len(argv) > 1) :
#    asl.set_argv(argv)
# else :
#    asl.set_argv_fileonly("/Users/renanspencertrindade/Desktop/zGenericSolver/tmp/prob2.nl")
# asl.start()
# asl.write_sol()
# # asl.write_solf("tmp/nlFuncNew4.sol")

# n_var = asl.get_n_var()
# x0 = asl.get_X0()
# print(x0)
# print(n_var)

# print( asl.get_file_nl() )
# print( "Number of NL func: " + str(asl.get_nl_count() ) )
# varx = asl.get_nl_varx()
# vary = asl.get_nl_vary()
# print( varx )
# print( vary )

# varxlb = asl.get_nl_varxlb()
# varxub = asl.get_nl_varxub()
# print(varxlb)
# print(varxub)

# listnl = asl.get_nl_functionslist()
# print(listnl)

# print("Final test")

# print(argv)