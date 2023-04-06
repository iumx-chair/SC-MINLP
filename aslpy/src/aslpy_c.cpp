#include <iostream>
#include <cstdio>
#include <vector>
#include <map>
#include <string>
#include <vector>

#include <sstream>
#include <fstream>
#include <iomanip>

#include "asl.h"
#include "nlp.h"
#include "getstub.h"
#include "r_opn.hd" /* for N_OPS */
#include "opcode.hd"

#include "display.h"

#define R_OPS ((ASL_fg *)asl)->I.r_ops_
#define OBJ_DE ((ASL_fg *)asl)->I.obj_de_
#define CON_DE ((ASL_fg *)asl)->I.con_de_

#define EPSILONTOLERANCE 1e-30

FILE *nl;
ASL *asl;
char *stub, *getenvErr;
char **argvN = 0;
int argcN = 0;

using namespace std;

std::ostringstream nlFuncFile, linearFuncFile;
vector<char *> nlFuncList;
char *nlFuncFileChar = 0;
char *linearFuncFileChar = 0;
vector<int> varx_indexes;
vector<int> vary_indexes;
vector<int> varz_indexes;
int countLin = 0;
int countNLin = 0;

fint optionsreturn[5];
char *coptionsreturn;

bool *intvar;
bool *nlinearCons;

// fint pc_refmode = 1;
// fint timelimit = (fint)1e75;
// fint timing = 0;

static keyword keywds[] = {
    /* must be alphabetical */
    KW((char *)"display", L_val, &optionsreturn[0], (char *)"Frequency of information display; default 1."),
    KW((char *)"pc_fileph1", *C_val, &coptionsreturn, (char *)"PC reformulation mode: 1 Incremental, 2 Multiple Choice; default = 1."),
    KW((char *)"pc_refmode", L_val, &optionsreturn[1], (char *)"PC reformulation mode: 1 Incremental, 2 Multiple Choice; default = 1."),
    KW((char *)"relax", L_val, &optionsreturn[2], (char *)"Ignore integrality; default 0."),
    KW((char *)"time", L_val, &optionsreturn[3], (char *)"Time limit in seconds; default = 1e75."),
    KW((char *)"timing", L_val, &optionsreturn[4], (char *)"display timings for the run."),
};

char *usage[] = {"Usage", "Test"};

static Option_Info Oinfo = {
    (char *)"SCMINLP",                          /* invocation name of solver */
    (char *)"SCMINLP solver\nv0.01 20230406\n", /* solver name in startup "banner" */
    (char *)"SCMINLP_options",                  /* name of solver_options environment var */
    keywds, nkeywds,                            /* key words */
    0,                                          /* whether funcadd will be called, etc.: */
    (char *)"V0.01 20230406",                   /* for -v and Ver_key_ASL() */
    (char **)usage,                             /* solver-specific usage message */
};

SufDecl suftab[] = {
    {(char *)"_scmnlp", 0, ASL_Sufkind_con}
    //  {"_scmnlp_x", 0, ASL_Sufkind_var},
    //  {"_scmnlp_y", 0, ASL_Sufkind_var}
    // { "integer", 0, ASL_Sufkind_var | ASL_Sufkind_output },
    // { "linear", 0, ASL_Sufkind_con | ASL_Sufkind_output },
    // { "feasible", 0, ASL_Sufkind_con | ASL_Sufkind_output },
    // { "evconvex", 0, ASL_Sufkind_con | ASL_Sufkind_output },
    // { "evconcave", 0, ASL_Sufkind_con | ASL_Sufkind_output },
    // { "amplsolverindex", 0, ASL_Sufkind_var | ASL_Sufkind_output },
    // { "linear", 0, ASL_Sufkind_obj | ASL_Sufkind_output },
    // { "infeasquant", 0, ASL_Sufkind_con | ASL_Sufkind_real | ASL_Sufkind_output}
};

void ASLPYC_start()
{
   optionsreturn[0] = -1;
   optionsreturn[1] = -1;
   optionsreturn[2] = -1;
   optionsreturn[3] = -1;
   optionsreturn[4] = -1;

   coptionsreturn = (char *)malloc(50 * sizeof(char));
   strcpy(coptionsreturn, "\0");

   // Oinfo.keywds[1]

   // char *str1 = (char *)malloc(200 * sizeof(char));
   // strcpy(str1, "/Users/renanspencertrindade/Desktop/zGenericSolver/wrapper");
   // char *str2 = (char *)malloc(200 * sizeof(char));
   // strcpy(str2, "/Users/renanspencertrindade/Desktop/zGenericSolver/tmp/prob.nl");
   // // strcpy(str2,"example/prob.nl");
   // char *str3 = NULL;

   // char **argv2 = (char **)malloc(3 * sizeof(char *));
   // argv2[0] = str1;
   // argv2[1] = str2;
   // argv2[2] = NULL;

   asl = ASL_alloc(ASL_read_fg);
   stub = getstub_ASL(asl, &argvN, &Oinfo);
   int a = (fint)strlen(stub);
   // printf("%s %d", stub, a);
   nl = jac0dim(stub, (fint)strlen(stub));

   // Allocate vectors
   X0 = (real *)Malloc(n_var * sizeof(real));  // Initial guess (if nonzero)
   LUv = (real *)Malloc(n_var * sizeof(real)); // variable lower (and, if Uvx == 0, upper) bounds
   Uvx = (real *)Malloc(n_var * sizeof(real)); // Variable upper bounds (if nonzero)

   LUrhs = (real *)Malloc(n_con * sizeof(real)); // Constraint lower (and, if Urhsx == 0, upper) bounds
   Urhsx = (real *)Malloc(n_con * sizeof(real)); // Constraint upper bounds (if nonzero)
   A_vals = (real *)Malloc(nzc * sizeof(real));  /* If nonzero, store constant Jacobian values */
                                                 /* (possibly 0 when nonlinearities are involved) */
                                                 /* in A_vals, A_rownos, and A_colstarts, */
                                                 /* rather than in Cgrad_. */

   // Declare suffixes
   suf_declare_ASL(asl, suftab, sizeof(suftab) / sizeof(SufDecl));

   // Options
   efunc *r_ops_int[N_OPS];

   for (int i = 0; i < N_OPS; i++)
      r_ops_int[i] = (efunc *)(unsigned long)i;
   R_OPS = r_ops_int;
   want_derivs = 0;
   fg_read(nl, 0);
   R_OPS = 0;

   /*** Get and record options ***/

   if (getopts(argvN, &Oinfo))
      exit(1);

   // int n_badvals = 0;

   // if (timelimit >= 0){
   //    cout << "Option value " << timelimit << " for directive time" << endl;
   // } else {
   //    cout << "Invalid value " << timelimit << " for directive time" << endl;
   //    cout << "Time limit in seconds; default = 1e75." << endl;
   //    timelimit = (fint)1e75;
   //    n_badvals++;
   // }

   // switch (timing)
   // {
   // case 0:
   //    break;
   // case 1:
   //    break;
   // default:
   //    cout << "Invalid value " << timing
   //          << " for directive timing" << endl;
   //    n_badvals++;
   // }

   // // Get tree
   // efunc *r_ops_int[N_OPS];

   // for (int i = 0; i < N_OPS; i++)
   // {
   //    r_ops_int[i] = (efunc *)(unsigned long)i;
   // }
   // R_OPS = r_ops_int;
   // want_derivs = 0;

   // // Read nl file
   // fg_read(nl, 0);

   // Display initial guess x0
   // for (int i = 0; i < n_var; i++)
   // {
   //    printf("%s = %f, >= %f <= %f \n", var_name(i), X0[i], LUv[i], Uvx[i]);
   // }

   // printf("Started");
   // Reading suffix
   // SufDesc *dp = suf_get("_scmnlp_x", ASL_Sufkind_var);
   // SufDesc *dp2 = suf_get("_scmnlp_y", ASL_Sufkind_var);
   SufDesc *dp3 = suf_get("_scmnlp", ASL_Sufkind_con);
   // for(int i = 0; i < n_var ; i++){
   //     printf("Var %d = %d, %d\n", i, dp->u.i[i], dp2->u.i[i]);
   // }

   // Transform linear jacobian in sparse form
   map<pair<int, int>, double> A;
   pair<int, int> p;
   for (int j = 0; j < n_var; j++)
   {
      p.second = j;
      for (int i = A_colstarts[j]; i < A_colstarts[j + 1]; i++)
      {
         p.first = A_rownos[i];
         if (fabs(A_vals[i]) < EPSILONTOLERANCE)
         {
            A_vals[i] = 0;
         }
         A[p] = A_vals[i];
         // printf("(%d, %d) = %f \n", p.first, p.second, A_vals[i]);
      }
   }

   // integrality of variables
   int thenlv = max(nlvc, nlvo);
   int theother = n_var - thenlv - nwv - nbv - niv;
   vector<bool> integrality(n_var, false);
   intvar = new bool[n_var];
   bool theresInteger = false;
   int i = 0;
   for (; i < nlvb - nlvbi; i++)
   {
      integrality[i] = false;
      intvar[i] = false;
   }
   for (; i < nlvb; i++)
   {
      theresInteger = true;
      integrality[i] = true;
      intvar[i] = true;
   }
   for (; i < nlvc - nlvci; i++)
   {
      integrality[i] = false;
      intvar[i] = false;
   }
   for (; i < nlvc; i++)
   {
      theresInteger = true;
      integrality[i] = true;
      intvar[i] = true;
   }
   for (; i < nlvo - nlvoi; i++)
   {
      integrality[i] = false;
      intvar[i] = false;
   }
   for (; i < nlvo; i++)
   {
      theresInteger = true;
      integrality[i] = true;
      intvar[i] = true;
   }
   for (; i < thenlv + nwv; i++)
   {
      integrality[i] = false;
      intvar[i] = false;
   }
   for (; i < thenlv + nwv + theother; i++)
   {
      integrality[i] = false;
      intvar[i] = false;
   }
   for (; i < n_var - niv; i++)
   {
      theresInteger = true;
      integrality[i] = true;
      intvar[i] = true;
      LUv[i] = 0;
      Uvx[i] = 1;
   }
   for (; i < n_var; i++)
   {
      theresInteger = true;
      integrality[i] = true;
   }

   map<pair<int, int>, double>::iterator mit;

   // objective function
   if (n_obj > 1)
   {
      cout << "nl2ros: cannot deal with more than 1 objective function\n";
      exit(2);
   }
   // for (i = 0; i < n_obj; i++)
   // {
   //    if (objtype[0] == 0)
   //    {
   //       linearFuncFile << "Minimize" << endl;
   //    }
   //    else
   //    {
   //       linearFuncFile << "Maximize" << endl;
   //    }
   //    for (ograd *og = Ograd[i]; og; og = og->next)
   //    {
   //       int j = og->varno;

   //       linearFuncFile << std::setprecision(30) << std::showpos << og->coef;
   //       linearFuncFile << " x" + to_string(j + 1) + " ";
   //    }
   // }
   // linearFuncFile << std::setprecision(30) << std::showpos << objconst(0);
   // linearFuncFile << endl;

   // Constraints
   bool nlinear;
   nlinearCons = new bool[n_con];
   // linearFuncFile << "Subject To" << endl;
   for (int i = 0; i < n_con; i++)
   {
      // Non-linearity test
      nlinear = !is_expr_zero((CON_DE + i)->e);
      nlinearCons[i] = nlinear;
      if (nlinear)
      {
         countNLin++;
      }

      std::ostringstream ss;
      p.first = i;

      // Load Non linear constraints
      int indexi = -1; // Non-linear variable
      int index_Z = -1;
      if (nlinear)
      {
         ss << "+";
         display_expr(ss, (CON_DE + i)->e, integrality, &indexi, &index_Z);
         // indexi = returnIndexVar((CON_DE + i)->e);
      }

      // Load Linear constraints
      int indexj = -1; // Linear variable
      for (int j = 0; j < n_var; j++)
      {
         p.second = j;
         mit = A.find(p);
         if (mit != A.end() && (mit->second > EPSILONTOLERANCE || mit->second < -EPSILONTOLERANCE))
         {
            ss << std::setprecision(30) << std::showpos << mit->second;
            if (nlinear)
            {
               if (j != indexi) {
                  ss << "*y ";
                  // ss << "*x" + to_string(j + 1) + " ";
                  // countNLin++;
                  indexj = j;
               } else {
                  ss << "*x ";
                  countLin++;
                  indexj = j; 
               }
            }
            else
            {
               ss << " x" + to_string(j + 1) + " ";
               countLin++;
               indexj = j;
            }
         }
      }

      // Export functions
      if (dp3->u.i[i] > 0)
      {
         nlFuncFile << "x" + to_string(indexi + 1) + " " + to_string(indexi + 1) << " ";
         nlFuncFile << to_string(LUv[indexi]) + " " + to_string(Uvx[indexi]) + "\n";
         nlFuncFile << "x" + to_string(indexj + 1) + " " + to_string(indexj + 1) + " \n";
         nlFuncFile << "subject to ";
         nlFuncFile << con_name(i);
         nlFuncFile << " :\n";
         nlFuncFile << std::setprecision(30) << "\t" << ss.str() << " <= " << Urhsx[i] << ";" << endl
                    << endl;

         // Store data
         varx_indexes.push_back(indexi);
         vary_indexes.push_back(indexj);
         varz_indexes.push_back(index_Z);

         std::ostringstream streamtmp;
         streamtmp << ss.str() << "\0";
         // char *nlTemp = (char *) malloc (streamtmp.str().size()*sizeof(char));
         char *nlTemp = new char[streamtmp.str().size() + 1];
         string temp = streamtmp.str();
         strcpy(nlTemp, temp.c_str());
         nlFuncList.push_back(nlTemp);
         // cout << ss.str() << endl;
      }

      // if (!nlinear)
      // {
      //    linearFuncFile << " C" << to_string(i) << "  : ";
      //    // linearFuncFile << con_name(i);

      //    if (fabs(LUrhs[i] - Urhsx[i]) > EPSILONTOLERANCE)
      //    {
      //       if (Urhsx[i] < Infinity)
      //       {
      //          linearFuncFile << std::setprecision(30) << ss.str() << " <= " << Urhsx[i] << endl;
      //       }
      //       else
      //       {
      //          linearFuncFile << std::setprecision(30) << ss.str() << " >= " << LUrhs[i] << endl;
      //       }
      //    }
      //    else
      //    {
      //       linearFuncFile << std::setprecision(30) << ss.str() << " = " << Urhsx[i] << endl;
      //    }
      // }
   }

   // // Include bounds LP File
   // linearFuncFile << "Bounds" << endl;
   // for (int i = 0; i < n_var; i++)
   // {
   //    linearFuncFile << std::setprecision(30) << LUv[i];
   //    linearFuncFile << " <= ";
   //    linearFuncFile << " x" + to_string(i + 1) + " ";
   //    linearFuncFile << " <= ";
   //    linearFuncFile << std::setprecision(30) << Uvx[i];
   //    linearFuncFile << " " << endl;
   // }

   // // General and Binaries
   // bool binexists = false;
   // bool firstgeneral = true;
   // if (theresInteger)
   // {
   //    for (int i = 0; i < n_var; i++)
   //    {
   //       if (integrality[i]) 
   //       {
   //          bool binCheck = LUv[i] >= 0.0 - EPSILONTOLERANCE && LUv[i] <= 0.0 + EPSILONTOLERANCE &&
   //                         Uvx[i] >= 1.0 - EPSILONTOLERANCE && Uvx[i] <= 1.0 + EPSILONTOLERANCE;
   //          if (!binCheck)
   //          {
   //             if (firstgeneral)
   //             {
   //                linearFuncFile << "General" << endl;
   //                firstgeneral = false;
   //             }
   //             linearFuncFile << " x" + to_string(i + 1) << endl;
   //          }
   //          else
   //          {
   //             binexists = true;
   //          }
   //       }
   //    }

   //    if (binexists)
   //    {
   //       linearFuncFile << "Binaries" << endl;
   //       for (int i = 0; i < n_var; i++)
   //       {
   //          if (integrality[i]){
   //             bool binCheck = LUv[i] >= 0.0 - EPSILONTOLERANCE && LUv[i] <= 0.0 + EPSILONTOLERANCE &&
   //                            Uvx[i] >= 1.0 - EPSILONTOLERANCE && Uvx[i] <= 1.0 + EPSILONTOLERANCE;
   //             if (binCheck)
   //             {
   //                linearFuncFile << " x" + to_string(i + 1) << endl;
   //             }
   //          }
   //       }
   //    }
   // }

   // End
   // linearFuncFile << "End";

   // std::ofstream fileNLfile;
   // fileNLfile.open("tmp/nlFunc2.txt");
   // fileNLfile << to_string(countNLin) << endl;
   // fileNLfile << nlFuncFile.str();
   // fileNLfile.close();

   // std::ofstream LPfile;
   // LPfile.open("tmp/nlFunc2.lp");
   // LPfile << linearFuncFile.str();
   // LPfile.close();

   // Change initial guess x0
   // for (int i = 0; i < n_var; i++)
   // // {
   //    X0[i] = Uvx[i];
   // // }

   // write_sol((char *)"ASL call is finished\n", X0, pi0, &Oinfo);
   // write_solf_ASL(asl, "Solved", X0, pi0, &Oinfo, "tmp/nlFunc2.sol");

   // printf("Filename:\n%s\n", asl->i.filename_);
}

extern "C"
{
   void ASLPYCstart() { ASLPYC_start(); }
   void ASLPYCwrite_solf(const char *fnamestr) { write_solf_ASL(asl, "Solved", X0, pi0, &Oinfo, fnamestr); }
   void ASLPYCwrite_sol(const char *fnamestr) { write_sol_ASL(asl, (char *)"ASL call is finished\n", X0, pi0, &Oinfo); }
   double ASLPYCgetX0(const int i) { return asl->i.X0_[i]; }
   void ASLPYCsetX0(const int i, double value) { asl->i.X0_[i] = value; }
   int ASLPYCgetnvar() { return asl->i.n_var_; }
   int ASLPYCgetncon() { return asl->i.n_con_; }
   void ASLPYCsetargc(int argc)
   {
      argvN = (char **)malloc((argc + 1) * sizeof(char *));
      argvN[argc] = NULL;
      argcN = argc;
   }
   void ASLPYCsetargv(int i, char *str)
   {
      char *str1 = (char *)malloc(strlen(str) * sizeof(char));
      strcpy(str1, str);
      argvN[i] = str1;
   }
   void ASLPYCwritelp(char *fnamestr)
   {
      std::ofstream LPfile;
      LPfile.open(fnamestr);
      LPfile << linearFuncFile.str();
      LPfile.close();
   }
   void ASLPYCwritenl(char *fnamestr)
   {
      std::ofstream LPfile;
      LPfile.open(fnamestr);
      LPfile << nlFuncFile.str();
      LPfile.close();
   }
   const char *ASLPYCgetlp()
   {
      if (!linearFuncFileChar)
      {
         free(linearFuncFileChar);
      }
      linearFuncFileChar = (char *)malloc(linearFuncFile.tellp() * sizeof(char));
      string temp = linearFuncFile.str();
      strcpy(linearFuncFileChar, temp.c_str());
      return (const char *)linearFuncFileChar;
   }

   // Return file .nls as a string.
   const char *ASLPYCgetnl()
   {
      if (!nlFuncFileChar)
      {
         free(nlFuncFileChar);
      }
      nlFuncFileChar = (char *)malloc(nlFuncFile.tellp() * sizeof(char));
      string temp = nlFuncFile.str();
      strcpy(nlFuncFileChar, temp.c_str());
      return (const char *)nlFuncFileChar;
   }

   // Return the nl function i.
   const char *ASLPYCgetnlfunc(int i)
   {
      return nlFuncList[i];
   }

   int ASLPYCgetnlcount()
   {
      return countNLin;
   }

   int ASLPYCgetnlvarx(int i)
   {
      return varx_indexes[i];
   }

   int ASLPYCgetnlvary(int i)
   {
      return vary_indexes[i];
   }

   int ASLPYCgetnlvarz(int i)
   {
      return varz_indexes[i];
   }

   double ASLPYCgetvarlb(int i)
   {
      return LUv[i];
   }

   double ASLPYCgetvarub(int i)
   {
      return Uvx[i];
   }

   double *ASLPYCgetLUv()
   {
      return LUv;
   }

   double *ASLPYCgetUvx()
   {
      return Uvx;
   }


   double *ASLPYCgetUrhsx()
   {
      return Urhsx;
   }

   double *ASLPYCgetLUrhs()
   {
      return LUrhs;
   }


   double ASLPYCgetconsub(int i)
   {
      return Urhsx[i];
   }

   double ASLPYCgetconslb(int i)
   {
      return LUrhs[i];
   }

   int ASLPYCgetoptionint(int i)
   {
      return optionsreturn[i];
   }

   char *ASLPYCgetoptionchar(int i)
   {
      return coptionsreturn;
   }

   expr *ASLPYgetexpr(int i)
   {
      return (CON_DE + i)->e;
   }

   ograd *ASLPYgetOgrad(int i)
   {
      return Ograd[i];
   }

   double ASLPYgetexpr_e(expr *e)
   {
      return ((expr_n *)e)->v;
   }

   double *ASLPYgetA_vals_()
   {
      return A_vals;
   }

   int *ASLPYgetA_rownos_()
   {
      return A_rownos;
   }

   int *ASLPYgetA_colstarts_()
   {
      return A_colstarts;
   }

   bool *ASLPYgetintvar(){
      return intvar;
   }

   char * ASLPYgetobjtype(){
      return objtype;
   }

   double ASLPYgetobjconst(int i){
      return objconst(i);
   }

   bool * ASLPYgetnlinearCons(){
      return nlinearCons;
   }
}

int main(int argc, char **argv)
{
   ASLPYCsetargc(2);
   ASLPYCsetargv(0, (char *)"/Users/renanspencertrindade/Desktop/zGenericSolver 2/wrapper");
   ASLPYCsetargv(1, (char *)"/Users/renanspencertrindade/Desktop/zGenericSolver 2/tmp/prob.nl");
   ASLPYC_start();
   // usage_ASL(&Oinfo, 1);
   ASLPYCwritelp( (char *)"/Users/renanspencertrindade/Desktop/zGenericSolver 2/ASL_Lib/Rose_Part/tmp/prob2.lp" );
}