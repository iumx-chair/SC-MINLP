#include <iostream>
#include <cstdio>
#include <vector>
#include <map>
#include <string>

#include <sstream>
#include <iomanip>

#include "asl.h"
#include "nlp.h"
#include "getstub.h"

#include "r_opn.hd" /* for N_OPS */
#include "opcode.hd"

#define EPSILONTOLERANCE 1e-30

using namespace std;

bool is_expr_zero(expr *e)
{
  efunc *op = e->op;
  int opnum = Intcast op;
  if (opnum == OPNUM && fabs(e->dL) < EPSILONTOLERANCE)
  {
    return true;
  }
  else
  {
    return false;
  }
}

void display_expr(std::ostringstream &ss, expr *e, vector<bool> &integrality, int *index_X, int *index_Z)
{

  using namespace std;

  efunc *op;
  expr **ep;
  int opnum;
  float number;
  float y, z, k, ind;
  char buffer[50];
  char *pEnd;
  double t;

  op = e->op;
  opnum = Intcast op;

  /* ss << ("op %d  optype %d  ", opnum, optype[opnum]);*/

  switch (opnum)
  {

  case OPPLUS:
    ss << "(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << "+";
    display_expr(ss,e->R.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OPMINUS:
    ss << "(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << "-";
    display_expr(ss,e->R.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OPMULT:
    ss << "(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << "*";
    display_expr(ss,e->R.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OPDIV:
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << "/";
    display_expr(ss,e->R.e, integrality, index_X, index_Z);
    return;

  case OPREM:
    ss << "remainder\n";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    display_expr(ss,e->R.e, integrality, index_X, index_Z);
    return;

  case OPPOW:
    ss << "(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << ")^(";
    display_expr(ss,e->R.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OPLESS:
    ss << "less\n";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    display_expr(ss,e->R.e, integrality, index_X, index_Z);
    return;

  case MINLIST:
    ss << "min -- not implemented\n";
    exit(1);

  case MAXLIST:
    ss << "max -- not implemented\n";
    exit(1);

  case FLOOR:
    ss << "floor -- not implemented\n";
    exit(1);

  case CEIL:
    ss << "ceil -- not implemented\n";
    exit(1);

  case ABS:
    ss << "abs";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    return;

  case OPUMINUS:
    ss << "(-(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << "))";
    return;

  case OPIFnl:
    ss << "if nl -- not implemented\n";
    exit(1);

  case OP_tanh:
    ss << "tanh(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OP_tan:
    ss << "tan(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OP_sqrt:
    ss << "sqrt(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OP_sinh:
    ss << "sinh(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OP_sin:
    ss << "sin(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OP_log10:
    ss << "log10(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OP_log:
    ss << "log(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OP_exp:
    ss << "exp(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OP_cosh:
    ss << "cosh(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OP_cos:
    ss << "cos(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OP_atanh:
    ss << "atanh(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OP_atan2:
    ss << "atan2(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OP_atan:
    ss << "atan(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OP_asinh:
    ss << "asinh(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OP_asin:
    ss << "asin(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OP_acosh:
    ss << "acosh(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OP_acos:
    ss << "acos(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OPSUMLIST:
    ss << "(";
    for (ep = e->L.ep; ep < e->R.ep; *ep++)
    {
      display_expr(ss,*ep, integrality, index_X, index_Z);
      if (ep < e->R.ep - 1)
      {
        ss << "+";
      }
    }
    ss << ")";
    return;

  case OPintDIV:
    ss << "int division";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    display_expr(ss,e->R.e, integrality, index_X, index_Z);

  case OPprecision:
    ss << "precision";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    return;

  case OPround:
    ss << "round";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    return;

  case OPtrunc:
    ss << "trunc";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    return;

  case OP1POW:
    ss << "(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    t = e->R.en->v;
    if (t >= 0)
    {
      ss << ")^" << e->R.en->v;
    }
    else
    {
      ss << ")^(" << e->R.en->v << ")";
    }
    return;

  case OP2POW:
    ss << "(";
    display_expr(ss,e->L.e, integrality, index_X, index_Z);
    ss << ")^2";
    return;

  case OPCPOW:
    t = e->L.en->v;
    if (t >= 0)
    {
      ss << e->L.en->v << "^(";
    }
    else
    {
      ss << "(" << e->L.en->v << ")^(";
    }
    display_expr(ss,e->R.e, integrality, index_X, index_Z);
    ss << ")";
    return;

  case OPFUNCALL:
    ss << "function call -- not implemented\n";
    exit(1);

  case OPNUM:
    t = ((expr_n *)e)->v;
    if (t >= 0)
    {
      ss << std::setprecision(30) << t;
    }
    else
    {
      ss << std::setprecision(30) << std::showpos << "(" << t << ")";
    }
    return;

  case OPPLTERM:
    ss << "pl term -- not implemented\n";
    exit(1);

  case OPIFSYM:
    ss << "if sym -- not implemented\n";
    exit(1);

  case OPHOL:
    ss << "string argument -- not implemented\n";
    exit(1);

  case OPVARVAL:
    // ss << ("(x%d)", e->a + 1);

    // Identifying variables
    if (integrality[e->a]== false){
      if (*index_X == -1){
        *index_X = e->a;
      } else {
        if (*index_X != e->a) {
          cout << "Error! Constraint with two or more continuous variables.";
          exit(1);
        }
      }
      ss << "x";
      // ss << std::noshowpos << e->a + 1;
      ss << " ";
    } else {
      if (*index_Z == -1){
        *index_Z = e->a;
      } else {
        if (*index_Z != e->a) {
          cout << "Error! Constraint with two or more binary variables.";
          exit(1);
        }
      }
      ss << "z";
      // ss << std::noshowpos << e->a + 1;
      ss << " ";
    }

    // cout << ((integrality[*index_X] == true) ? "T" : "F") << ", " << *index_X << "\n";
    // if (integrality[*index_X] == true) {
    //   ss << "T";
    //   // ss << std::noshowpos << e->a + 1;
    //   ss << " ";
    // } else {
    //   ss << "x";
    //   // ss << std::noshowpos << e->a + 1;
    //   ss << " ";
    // }
    break;

  default:
    cout << "other -- not implemented\n";
    // exit(1);
  }
}

int returnIndexVar(expr *e)
{

  using namespace std;

  efunc *op;
  expr **ep;
  int opnum;
  float number;
  float y, z, k, ind;
  char buffer[50];
  char *pEnd;
  double t;

  op = e->op;
  opnum = Intcast op;

  if (opnum == OPVARVAL) {
    return e->a;
  }

  int indexL = -1;
  int indexR = -1;

  switch (opnum)
  {

  case OPPLUS:
    indexL = returnIndexVar(e->L.e);
    indexR = returnIndexVar(e->R.e);
    break;
    
  case OPMINUS:
    indexL = returnIndexVar(e->L.e);
    indexR = returnIndexVar(e->R.e);
    break;
    
  case OPMULT:
    indexL = returnIndexVar(e->L.e);
    indexR = returnIndexVar(e->R.e);
    break;
    
  case OPDIV:
    indexL = returnIndexVar(e->L.e);
    indexR = returnIndexVar(e->R.e);
    break;
    
  case OPREM:
    indexL = returnIndexVar(e->L.e);
    indexR = returnIndexVar(e->R.e);
    break;
    
  case OPPOW:
    indexL = returnIndexVar(e->L.e);
    indexR = returnIndexVar(e->R.e);
    break;
    
  case OPLESS:
    indexL = returnIndexVar(e->L.e);
    indexR = returnIndexVar(e->R.e);
    break;
    
  case MINLIST:
    exit(1);
    break;
    
  case MAXLIST:
    exit(1);
    break;
    
  case FLOOR:
    exit(1);
    break;
    
  case CEIL:
    exit(1);
    break;
    
  case ABS:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OPUMINUS:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OPIFnl:
    exit(1);
    break;
    
  case OP_tanh:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OP_tan:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OP_sqrt:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OP_sinh:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OP_sin:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OP_log10:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OP_log:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OP_exp:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OP_cosh:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OP_cos:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OP_atanh:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OP_atan2:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OP_atan:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OP_asinh:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OP_asin:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OP_acosh:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OP_acos:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OPSUMLIST:
    for (ep = e->L.ep; ep < e->R.ep; *ep++)
    {
      indexL = returnIndexVar(*ep);
    }
    break;
    
  case OPintDIV:
    indexL = returnIndexVar(e->L.e);
    indexR = returnIndexVar(e->R.e);
    break;
    
  case OPprecision:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OPround:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OPtrunc:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OP1POW:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OP2POW:
    indexL = returnIndexVar(e->L.e);
    break;
    
  case OPCPOW:
    indexR = returnIndexVar(e->R.e);
  }

  if (indexL >= 0){
    if (indexR >= 0){ 
      if (indexL != indexR){
        return -2;
      }
    }
    return indexL;
  }
  if (indexR >= 0){ 
    return indexR;
  }

  return -1;
}