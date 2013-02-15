
#include "ngraph/ngraph.hpp"
#include "mergesort.hpp"
//#include "complem_knap.hpp"
#include "LP_relax.hpp"
#include "conflict_problem.hpp"
#include "scip/pub_var.h"          // SCIPvarGetObj
#include "scip/pub_cons.h"         // SCIPconsGetHdlr    

#define _INF std::numeric_limits<double>::infinity();

using namespace std;
using namespace NGraph;




static
SCIP_RETCODE cutgen(SCIP* scip)
{
  ostringstream namebuf;
  SCIP_CALL(SCIPsetMessagehdlr(scip,NULL)); //suppress messages

  Graph Gr; //Gr = conflict Graph of scip

  ///////////////////
  // get scip data //
  ///////////////////
  SCIP_CONS** cons=SCIPgetConss(scip); // cons = constraints of scip
  SCIP_VAR** vars=SCIPgetVars(scip); //vars = variables of scip
  int numbercons=SCIPgetNConss(scip); //numbercons = number of constraints
  int n=SCIPgetNVars(scip); // n = number of items/variables

  
  ////////////////////
  //   presolving   //
  ////////////////////

  SCIP_CALL(SCIPpresolve(scip));


  //////////////////////////////////////////
  // compute x_star=solution of LP relax. //
  //////////////////////////////////////////

  double *x_star;
  x_star = new double[n];
  SCIP_CALL(LP_relax(scip,x_star)); 
  
  ////////////////////////////////
  //    create conflict graph   //
  ////////////////////////////////

  for (int i=0; i<numbercons; i++)
  {
    if(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons[i])), "SOS1") == 0)
    {
      int numvarsSOS1=SCIPgetNVarsSOS1(scip,cons[i]);
      if(numvarsSOS1==2)
      {
        SCIP_VAR** getvarssos=SCIPgetVarsSOS1(scip,cons[i]);
        int save_getvar1=-1;
        int save_getvar2=-1;
        for (int j=0; j<n; j++)
        {    
          if(getvarssos[0]==vars[j])
          {
            save_getvar1=j;
          }
          if(getvarssos[1]==vars[j])
          {
            save_getvar2=j;
          }
          if(save_getvar1!=-1 && save_getvar2!=-1)
          {
            Gr.insert_undirected_edge(save_getvar1,save_getvar2);
            break;
          }
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////
  //   for every linear constraint: try to generate a cut  //
  ///////////////////////////////////////////////////////////
  
  for (int d=0; d<numbercons; d++)
  {
    if(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons[d])), "linear") == 0)
    {
      ////////////////////////////////////////
      // create and solve conflict problem  //
      ////////////////////////////////////////

      Graph Gr_conf;//conflict graph of conflict problem
      int n_vars_cons_d=SCIPgetNVarsLinear(scip,cons[d]);//number of variables of constraint d 
      double OBJVAL_conf;//variable to store objective value if conflict problem 
      double *w;//pointer to store solution of conflict problem
      w = new double[n_vars_cons_d];
      SCIP_CALL(conflict_problem(scip,n,numbercons,d,vars,cons,x_star,w,&OBJVAL_conf,&Gr_conf));


      delete [] w;      //delete solution of conflict problem
    }
  }
  


  delete [] x_star; //solution of LP relax
  //delete [] Gr;   //conflict Graph
  //delete [] Gr_conf;   //conflict Graph
  return SCIP_OKAY;
}
