// 1.                   get scip data
// 2.                   presolving 
// 3.                   compute x_star=solution of LP relax.
// 4.                   create conflict graph
// 5.                   for every linear constraint d: try to generate a cut
//    5.1               get data of constraint d
//    5.2               create and solve conflict problem
//    5.3               check if a cover is at hand
//       5.3.1          compute cover
//       5.3.2          determine inequality \sum_{i \in C without N(nue)} a_i^(d) + a_nue^(d)
//                      ( N(nue)=neighbors of nue in conflict graph Gr_conf )
//            5.3.2.1   get data for inequality
//            5.3.2.2   determine nue 



#include "ngraph/ngraph.hpp"
#include "mergesort.hpp"
//#include "complem_knap.hpp"
#include "LP_relax.hpp"
#include "conflict_problem.hpp"
#include "scip/pub_var.h"          // SCIPvarGetObj
#include "scip/pub_cons.h"         // SCIPconsGetHdlr    

#define _INF std::numeric_limits<double>::infinity();
#define delta 0.0001; //tolerance

using namespace std;
using namespace NGraph;





static
SCIP_RETCODE cutgen(SCIP* scip)
{
  ostringstream namebuf;
  SCIP_CALL(SCIPsetMessagehdlr(scip,NULL)); //suppress messages

  double tol=delta;
  Graph Gr; //Gr = conflict Graph of scip
 
  /////////////////////////
  // 1.   get scip data  //
  /////////////////////////
  SCIP_CONS** cons=SCIPgetConss(scip); // cons = constraints of scip
  SCIP_VAR** vars=SCIPgetVars(scip); //vars = variables of scip
  int numbercons=SCIPgetNConss(scip); //numbercons = number of constraints
  int n=SCIPgetNVars(scip); // n = number of items/variables

  
  //////////////////////////
  // 2.     presolving   //
  //////////////////////////

  SCIP_CALL(SCIPpresolve(scip));


  /////////////////////////////////////////////////
  // 3.    compute x_star=solution of LP relax. //
  /////////////////////////////////////////////////

  double *x_star;
  x_star = new double[n];
  SCIP_CALL(LP_relax(scip,x_star)); 
  
  /////////////////////////////////////
  // 4.     create conflict graph   //
  /////////////////////////////////////

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

  /////////////////////////////////////////////////////////////
  // 5.  for every linear constraint d: try to generate a cut  //
  /////////////////////////////////////////////////////////////
  
  for (int d=0; d<numbercons; d++)
  {
    if(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons[d])), "linear") == 0)
    {

      //////////////////////////////////////
      // 5.1    get data of constraint d  //
      //////////////////////////////////////

      const char* name_cons_d=SCIPconsGetName(cons[d]);        //name of constraint d
      SCIP_VAR** vars_cons_d=SCIPgetVarsLinear(scip,cons[d]);  //vars of constraint d
      // xxx varslin
      SCIP_Real b= SCIPgetRhsLinear(scip,cons[d]);             //rhs of constraint d
      SCIP_Real* a= SCIPgetValsLinear(scip,cons[d]);           //values of constraint d
      int n_vars_cons_d=SCIPgetNVarsLinear(scip,cons[d]);      //number of variables of constraint d    
      // xxx numvarslinear

      // vector to store position of variables of constraint d
      int posvars_cons_d[n_vars_cons_d];
      // xxx posvarslin

      int vars_it=0; //iterator
      for (int i=0; i<n; i++)
      {
        if(vars[i]==vars_cons_d[vars_it])
        {
          posvars_cons_d[vars_it]=i;
          vars_it+=1;
        }
      }


      ////////////////////////////////////////////////
      // 5.2     create and solve conflict problem  //
      ////////////////////////////////////////////////

      Graph Gr_conf;//conflict graph of conflict problem
      double OBJVAL_conf;//variable to store objective value if conflict problem 
      double *w;//pointer to store solution of conflict problem
      w = new double[n_vars_cons_d];
      int null_norm_w;//variable to store the null (quasi) norm of w 
      SCIP_CALL(conflict_problem(scip,n,numbercons,d,n_vars_cons_d,vars_cons_d,a,vars,cons,x_star,w,&OBJVAL_conf,&null_norm_w,&Gr_conf));


      /////////////////////////////////////////
      // 5.3    check if a cover is at hand  //
      /////////////////////////////////////////

      if(OBJVAL_conf>b+tol) // cover detected
      {
        cout << endl << "cover detected" << endl << endl;

        ///////////////////////////////
        // 5.3.1    compute cover    //
        ///////////////////////////////

        int size_cover=null_norm_w; // size of cover
        int cover_xvars[size_cover]; // cover (in indices of x_vars)
        int cover_wvars[size_cover]; // cover (in indices of w_vars)
        int cover_it=0; //iterator
        for (int i=0; i<n_vars_cons_d; i++)
        {
          if(w[i]>0.5) //(w[i]==1)
	  {
            cover_xvars[cover_it]=posvars_cons_d[i];
            cover_wvars[cover_it]=i;
            cover_it+=1;
	  }
        }


        ///////////////////////////////////////////////////////////////////////////////////////
        // 5.3.2    determine inequality\sum_{i \in C without N(nue)} a_i^(d) + a_nue^(d)    //
        //          ( N(nue)=neighbors of nue in conflict graph Gr_conf )                    //
        ///////////////////////////////////////////////////////////////////////////////////////


        /////////////////////////////////////////
        // 5.3.2.1    get data for inequality  //
	/////////////////////////////////////////

        double sum_C=OBJVAL_conf;                        // sum_C := \sum_{i \in C } a_i^(d)
        double sum_L;                                    // sum_L := \sum_{i \in C without N(nue)} a_i^(d) + a_nue^(d)

        int n_wnull=n_vars_cons_d-size_cover;            // number of w[i] with w[i]==0
        double a_wnull[n_wnull];                         // set {a_i^{d} : w[i]==0}
        double a_wnull_sort[n_wnull];                    // a_wnull sorted (maximum value - first position)
        int pos_w_a_wnull[n_wnull];                      // vector to store position of a_i^{d} with w[i]==0 in w_vars
        int pos_w_a_wnull_sort[n_wnull];                 // pos_w_a_wnull in same order as a_wnull_sorted
        //int pos_x_a_wnull[n_vars_cons_d-size_cover];   // vector to store position of a_i^{d} with w[i]==0 in x_vars
        
        vars_it=0; //iterator
        for (int i=0; i<n_vars_cons_d; i++)
        {
          if(w[i]<0.5) //(w[i]==0)
          {
            a_wnull[vars_it]=a[i];
            a_wnull_sort[vars_it]=a[i];
            pos_w_a_wnull[vars_it]=i;
            pos_w_a_wnull_sort[vars_it]=i;
            //pos_x_a_wnull[vars_it]=posvars_cons_d[i];
            vars_it+=1;
          }
        }
        mergesort(a_wnull_sort,pos_w_a_wnull_sort,0,n_wnull-1);  // sort a_wnull
        
        ///////////////////////////////
        // 5.3.2.2    determine nue  //
	///////////////////////////////
        
        double a_nue;      //a_nue
        int nue_w;         //index nue in indices of w_vars
        int nue_x;         //index nue in indeces of x_vars
        int success=0;     //success=0 -> nue does not exist  

        for (int i=0; i<n_wnull; i++)
        {
          sum_L=sum_C;
          a_nue=a_wnull_sort[i];          //a_nue     
	  nue_w=pos_w_a_wnull_sort[i];    //index nue in indices of w_vars
          nue_x=posvars_cons_d[nue_w];    //index nue in indices of x_vars
          sum_L+=a_nue;
          // determine neighbors of nue in conflict graph of the conflict problem
          Graph::vertex_set S = Gr_conf.out_neighbors(nue_w); 
          for (Graph::vertex_set::const_iterator t = S.begin(); t !=S.end(); t++)
          {
            if(w[*t]>0.5)  // (w[*t]==1)
	    {
              sum_L-=a[*t];
	    }
	  }
          if(sum_L<b-tol)
	  {
            success=1; //nue detected
            break;
	  }     
        }

        if(success==1)
	{

	} //end if success==1 (nue detected)
        



      } //end check if a cover is at hand




      delete [] w;      //delete solution of conflict problem
    } //end if cons[d] ^= linear
  } //end for every linear constraint: try to generate a cut
  


  delete [] x_star; //solution of LP relax
  //delete [] Gr;   //conflict Graph
  //delete [] Gr_conf;   //conflict Graph
  return SCIP_OKAY;
}
