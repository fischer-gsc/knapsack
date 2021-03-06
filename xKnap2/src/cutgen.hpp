// 1.                   presolving
// 2.                   get scip data
// 3.                   compute x_star=solution of LP relax.
// 4.                   create conflict graph
// 5.                   for every linear constraint d: try to generate a cut
//    5.1               get data of constraint d
//    5.2               create and solve conflict problem
//    5.3               check if a cover is at hand
//       5.3.1          compute cover
//       5.3.2          determine inequality \sum_{i \in C without N(nue)} a_i^(d) + a_nue^(d) < b^(d)
//                      ( N(nue)=neighbors of nue in conflict graph Gr_conf )
//            5.3.2.1   get data for inequality
//            5.3.2.2   determine nue
//            5.3.2.3   determine L and size of L ( L=C without N(nue) with nue )
//       5.3.3          lifting
//       5.3.4          include lifted inequality in scip





struct SCIP_VarData
{
  int index;
};



//#define SCIP_DEBUG
#include "ngraph/ngraph.hpp"
#include "mergesort.hpp"
#include "complem_knap.hpp"
#include "LP_relax.hpp"
#include "conflict_problem.hpp"
#include "lifting.hpp"
#include "scip/pub_var.h"          // SCIPvarGetObj   // SCIPvarSetData
#include "scip/pub_cons.h"         // SCIPconsGetHdlr    

#define _INF std::numeric_limits<double>::infinity();
#define delta 0.0001; //tolerance

using namespace std;
using namespace NGraph;




static
SCIP_RETCODE cutgen(SCIP* scip)
{
  ostringstream namebuf;

  int cut_it=0; //count number of generated cuts
  vector<double> ZFW;

  for (int iter=0; iter<10; iter++)
  {

    double tol=delta; 
    Graph Gr; //Gr = conflict Graph of scip


    //////////////////////////
    // 1.     presolving   //
    //////////////////////////

    //SCIP_CALL(SCIPpresolve(scip));

    SCIP_CALL(SCIPsetMessagehdlr(scip,NULL)); //suppress messages

    /////////////////////////
    // 2.   get scip data  //
    /////////////////////////

    SCIP_CONS** cons=SCIPgetConss(scip); // cons = constraints of scip
    SCIP_VAR** vars=SCIPgetVars(scip); //vars = variables of scip
    int numbercons=SCIPgetNConss(scip); //numbercons = number of constraints
    int n=SCIPgetNVars(scip); // n = number of items/variables


    // define indices of variables
    SCIP_VARDATA vardata[n];
    for (int i=0; i<n; i++)
    {
      (vardata[i]).index=i;
      SCIPvarSetData(vars[i],&vardata[i]);
    }


    /////////////////////////////////////////////////
    // 3.    compute x_star=solution of LP relax. //
    /////////////////////////////////////////////////

    SCIP_Real objval_LP;
    double *x_star;
    x_star = new double[n];
    SCIP_CALL(LP_relax(scip,cons,vars,numbercons,n,&vardata[n],x_star,&objval_LP));
    ZFW.push_back(objval_LP);


    /////////////////////////////////////
    // 4.     create conflict graph   //
    /////////////////////////////////////

    int numvarsSOS1;
    SCIP_VAR** getvarssos;
    SCIP_VARDATA *vardataout0;
    SCIP_VARDATA *vardataout1;

    for (int i=0; i<numbercons; i++)
    {
      if(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons[i])), "SOS1") == 0)
      {
	numvarsSOS1=SCIPgetNVarsSOS1(scip,cons[i]);
	if(numvarsSOS1==2)
	{
	  getvarssos=SCIPgetVarsSOS1(scip,cons[i]);
	  vardataout0=SCIPvarGetData(getvarssos[0]);
	  vardataout1=SCIPvarGetData(getvarssos[1]);
	  Gr.insert_edge((*vardataout0).index,(*vardataout1).index);
	  Gr.insert_edge((*vardataout1).index,(*vardataout0).index);
	}
      }
    }

    ///////////////////////////////////////////////////////////////
    // 5.  for every linear constraint d: try to generate a cut  //
    ///////////////////////////////////////////////////////////////

    for (int d=0; d<numbercons; d++)
    {
      if(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons[d])), "linear") == 0)
      {

	//////////////////////////////////////
	// 5.1    get data of constraint d  //
	//////////////////////////////////////

	SCIP_VAR** vars_cons_d=SCIPgetVarsLinear(scip,cons[d]);  //vars of constraint d
	// xxx varslin
	SCIP_Real b= SCIPgetRhsLinear(scip,cons[d]);             //rhs of constraint d
	SCIP_Real* a= SCIPgetValsLinear(scip,cons[d]);           //values of constraint d
	int n_vars_cons_d=SCIPgetNVarsLinear(scip,cons[d]);      //number of variables of constraint d    
	// xxx numvarslinear


	int posvars_cons_d[n_vars_cons_d]; // vector to store pos. of variables of constraint d in obj-func
	int* posvars_obj_d; // vector to store pos. of variables of obj-func in constraint d
			    // (-1 if vatiable does not exist in constraint d)
	posvars_obj_d = new int[n];
	// xxx posvarslin

	int vars_it=0; //iterator
	for (int i=0; i<n; i++)
	{
	  if(vars[i]==vars_cons_d[vars_it])
	  {
	    posvars_cons_d[vars_it]=i;
	    posvars_obj_d[i]=vars_it;
	    vars_it+=1;
	  }
	  else
	  {
	    posvars_obj_d[i]=-1;
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
	SCIP_CALL(conflict_problem(scip,n,numbercons,d,n_vars_cons_d,vars_cons_d,a,posvars_obj_d,&vardata[n],vars,cons,x_star,w,&OBJVAL_conf,&null_norm_w,&Gr_conf));


	/////////////////////////////////////////
	// 5.3    check if a cover is at hand  //
	/////////////////////////////////////////

	if(OBJVAL_conf>b+tol) // cover detected
	{
	  cout << endl << "cover detected" << endl << endl;

	  ///////////////////////////////
	  // 5.3.1    compute cover    //
	  ///////////////////////////////

	  int size_cover=null_norm_w;  // size of cover
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


	  //////////////////////////////////////////////////////////////////////////////////////////////
	  // 5.3.2    determine inequality \sum_{i \in C without N(nue)} a_i^(d) + a_nue^(d) < b^(d)  //
	  //          ( N(nue)=neighbors of nue in conflict graph Gr_conf )                           //
	  //////////////////////////////////////////////////////////////////////////////////////////////


	  /////////////////////////////////////////
	  // 5.3.2.1    get data for inequality  //
	  /////////////////////////////////////////

	  double sum_C=OBJVAL_conf;                        // sum_C := \sum_{i \in C } a_i^(d)   (C=cover_wvars)
	  double sum_L;                                    // sum_L := \sum_{i \in C without N(nue)} a_i^(d) + a_nue^(d)

	  int n_wnull=0;                                   // number of w[i] with w[i]==0 and x_star[posvars_cons_d[i]]!=0
	  vector<double> a_wnull;                          // set {a_i^{d} : w[i]==0 and x_star[posvars_cons_d[i]]!=0}
	  vector<double> a_wnull_sort;                     // a_wnull sorted (minimum value - first position)
	  vector<int> pos_w_a_wnull;                       // vector to store position of a_i^{d} with w[i]==0 and 
							   // x_star[posvars_cons_d[i]]!=0  in w_vars
	  vector<int> pos_w_a_wnull_sort;                  // pos_w_a_wnull in same order as a_wnull_sorted
	  //int pos_x_a_wnull[n_vars_cons_d-size_cover];   // vector to store position of a_i^{d} with w[i]==0 in x_vars

	  for (int i=0; i<n_vars_cons_d; i++)
	  {
	    if((w[i]<0.5)&&(x_star[posvars_cons_d[i]]>tol))       //(w[i]<0.5 ^= w[i]==0)
	    {
	      a_wnull.push_back(a[i]);
	      a_wnull_sort.push_back(a[i]);
	      pos_w_a_wnull.push_back(i);
	      pos_w_a_wnull_sort.push_back(i);
	      //pos_x_a_wnull[vars_it]=posvars_cons_d[i];
	      n_wnull+=1;
	    }
	  }

	  if(n_wnull>1)
	  {
	    mergesort(&a_wnull_sort,&pos_w_a_wnull_sort,0,n_wnull-1);  // sort a_wnull
	  }
	  else if(n_wnull==0)
	  {
            continue;  //continue 5. : for every linear constraint d: try to generate a cut 
	  }


	  //////////////////////////////
	  // 5.3.2.2    determine nue //
	  //////////////////////////////

	  double a_nue;              //a_nue
	  int nue_w;                 //index nue in indices of w_vars
	  int nue_x;                 //index nue in indeces of x_vars
	  int success=0;             //success=0 -> nue does not exist

	  Graph::vertex_set N_nue;   //neighbors of nue
	  for (int i=n_wnull-1; i>=0; i--)
	  {
	    sum_L=sum_C;
	    a_nue=a_wnull_sort.at(i);          //a_nue     
	    nue_w=pos_w_a_wnull_sort.at(i);    //index nue in indices of w_vars
	    nue_x=posvars_cons_d[nue_w];    //index nue in indices of x_vars
	    sum_L+=a_nue;
	    // determine neighbors of nue in conflict graph of the conflict problem
	    N_nue = Gr_conf.out_neighbors(nue_w); 
	    for (Graph::vertex_set::const_iterator t = N_nue.begin(); t !=N_nue.end(); t++)
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
	    //////////////////////////////////////////////////////////////////////////
	    // 5.3.2.3   determine L and size of L ( L=C with nue )  //
	    //////////////////////////////////////////////////////////////////////////

	    int size_L=size_cover+1; //size of cover + index nue
	    vector<int> cover_longnull(n_vars_cons_d,0); // cover_longnull[i]=0 , if i not in C
							 // cover_longnull[i]=one of the cover indices in xvars , else
	    for (int i=0; i<size_cover; i++)
	    {
	      cover_longnull.at(cover_wvars[i])=1;
	    }

	    /*
	    for (Graph::vertex_set::const_iterator t = N_nue.begin(); t !=N_nue.end(); t++)
	    {
	      if(w[*t]>0.5)  // (w[*t]>0.5 -> w[*t]==1 -> *t \in C)
	      {
		cover_longnull.at(*t)=0;;
		size_L-=1;
	      }
	    }
	    */

	    vector<int> L;
	    int L_it=0; //iterator
	    for (int i=0; i<n_vars_cons_d; i++)
	    {
	      if(cover_longnull.at(i)==1)
	      {
		L.push_back(i); //L[L_it]=i;
		L_it+=1;
	      }
	    }
	    L.push_back(nue_w); //remark: L sorted without regard to last entry
	    //sort L
	    for (int i=0; i<size_L-1; i++)
	    {
	      if(L.at(i)>nue_w)
	      {
		for (int j=size_L-2; j>=i; j--)
		{
		  L.at(j+1)=L.at(j);
		}
		L.at(i)=nue_w;
		break;
	      }
	    }

	    //////////////////////
	    // 5.3.3   lifting  //
	    //////////////////////

	    // determine alpha_nue = b^(d)-\sum_{i \in C without N(nue)} a_i^(d)
	    double alpha_nue = b - (sum_L-a_nue);


	    // determine I=complement of L (I in w_vars)
	    int size_I=n_vars_cons_d-size_L;
	    vector<int> I;

	    L_it=0;     //L iterator 
	    for (int i=0; i<n_vars_cons_d; i++)
	    {
	      if(L.at(L_it)==i)
	      {
		if(L_it<size_L-1)
		{
		  L_it+=1;
		}
	      }
	      else
	      {
		I.push_back(i);
	      }
	    }

	    //sort I in increasing order of values of x_star restricted to I_in_x_vars
	    vector<double> x_star_restricted; //x_star restricted to I_in_x_vars
	    for (int i=0; i<size_I; i++)
	    {
	      x_star_restricted.push_back(x_star[posvars_cons_d[I.at(i)]]);
	    }

	    if(size_I>1)
	    {
	      //sort I in increasing order of values of x_star_restricted
	      mergesort(&x_star_restricted,&I,0,size_I-1);
	    }

	    //define pointer to store values of lifted inequality
	    double *alpha;
	    alpha = new double[n_vars_cons_d];

	    //lift inequality in 5.3.2
	    SCIP_CALL(lifting(n_vars_cons_d,alpha,a,b,nue_w,alpha_nue,&Gr_conf,&L,&I,size_L,size_I));


	    /////////////////////////////////////////////////////////////////
	    // 5.3.4   when inequality is beneficial - include it in scip  //
	    /////////////////////////////////////////////////////////////////

	    //does the inequality violate the cuttent LP-relax. ?
	    double sum=0;
	    for (int i=0; i<n_vars_cons_d; i++)
	    {
	      sum+=alpha[i]*x_star[posvars_cons_d[i]];
	    }
	    if(sum > b+tol)   // inequality is beneficial
	    {
	      SCIP_CONS * cons_cut;
	      namebuf.str("");
	      namebuf<<"cut_"<<cut_it;
	      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons_cut, namebuf.str().c_str(), 0, NULL, NULL, 0.0 , b ) );

	      for (int i=0; i<n_vars_cons_d; i++)
	      {
		SCIP_CALL( SCIPaddCoefLinear(scip, cons_cut, vars[posvars_cons_d[i]], alpha[i] ) );
	      }
	      SCIP_CALL( SCIPaddCons(scip, cons_cut) );
	      cut_it+=1;
	      break; //no more cuts needed in this iteration
	    }


	    delete [] alpha;
	  } //end if success==1 (nue detected)
	} //end check if a cover is at hand




	delete [] w;      //delete solution of conflict problem
      } //end if cons[d] ^= linear

    } //end for every linear constraint: try to generate a cut



    delete [] x_star; //solution of LP relax
    //delete [] Gr;   //conflict Graph
    //delete [] Gr_conf;   //conflict Graph



    // for (int i=0; i<n; i++)
    //{
    //   cout << endl << x_star[i];
    //}

  }

  for (int i=0; i<static_cast<int> (ZFW.size ()); i++)
  {
    cout << endl << ZFW.at(i);
  }





  return SCIP_OKAY;
}
