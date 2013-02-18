//////////#include "scip/pub_var.h"

using namespace std;
using namespace NGraph;

static
SCIP_RETCODE lifting(SCIP* scip,             //SCIP data structure 
		     int n,                  //number of variables
		     int numbercons,         //number of constraints
		     SCIP_VAR** vars,        //variables of scip
                     SCIP_CONS** cons,       //constraints of scip
                     SCIP_Real* a,           //values of constraint d
                     SCIP_Real b,            //right hand side of constraint d
		     Graph* Gr_conf,         //pointer to store conflict graph of conflict problem
		     int* L,                 //set L=C without N(nue) with nue
		     int size_L              //size of L
		     )
{
  ostringstream namebuf;


  int k=1; //(lifting variable)


    /////////////////////////////
    // create lifting problem  //
    /////////////////////////////

    SCIP* lift = NULL;
    SCIP_CALL( SCIPcreate(&lift) );
    SCIP_CALL( SCIPincludeDefaultPlugins(lift) );
    SCIP_CALL( SCIPcreateProbBasic(lift, "lift") );
    SCIP_CALL(SCIPsetObjsense(lift, SCIP_OBJSENSE_MAXIMIZE));


    /////////////////////////////////////////////////////////
    // get neigbors of k and sort them in increasing order //
    /////////////////////////////////////////////////////////

    int deg_k = (*Gr_conf).degree(k);                       //degree of node k
    Graph::vertex_set N_k= (*Gr_conf).out_neighbors(k);     //neighbors of node k
    vector<double> N_k_copy(deg_k);                         //copy of N_k
    vector<int> nullvec(deg_k,0);                           //no need for mergesort  
    int it=0;                                               //iterator
    for (Graph::vertex_set::const_iterator t = N_k.begin(); t !=N_k.end(); t++)
    {
      N_k_copy.at(it)=*t;
      it+=1;
    }

    if(deg_k>1)
    {
      mergesort(&N_k_copy,&nullvec,0,deg_k-1);  // sort N_k_copy
    }
    nullvec.clear();

    ///////////////////////////////////////////////////
    // for all i\in L without N(k) create a variable //
    ///////////////////////////////////////////////////

    vector<int> L_wo_N_k; // L without N(k)
    int size_L_wo_N_k;    // size of L without N(k)

    int vars_it=0;  //vars iterator
    it=0;           //neigbor of k iterator
    for (int i = 0; i < size_L; ++i)
    {
      while(L[i]>N_k_copy.at(it))
      {
        if(it<deg_k-1)
	{
          it+=1;
	}
        else
	{
	  break;
        }
      }
      if(it==deg_k-1)
      {
        if(L[i]!=N_k_copy.at(it))
	{
          SCIP_VAR * var;
          namebuf.str("");
          namebuf << "x_lift#" << vars_it;

          // create the SCIP_VAR object
          SCIP_CALL( SCIPcreateVarBasic(lift, &var, namebuf.str().c_str(), 0.0, 1.0, a[L[i]], SCIP_VARTYPE_CONTINUOUS) );

          //add L[i] to L without N(k)
          L_wo_N_k.push_back(L[i]);

          // add the SCIP_VAR object to the scip problem
          SCIP_CALL( SCIPaddVar(lift, var) );
          vars_it+=1;
	}
      }
      else
      {
        if(L[i]<N_k_copy.at(it))
	{
          SCIP_VAR * var;
          namebuf.str("");
          namebuf << "x_lift#" << vars_it;

          // create the SCIP_VAR object
          SCIP_CALL( SCIPcreateVarBasic(lift, &var, namebuf.str().c_str(), 0.0, 1.0, a[L[i]], SCIP_VARTYPE_CONTINUOUS) );

          //add L[i] to L without N(k)
          L_wo_N_k.push_back(L[i]);

          // add the SCIP_VAR object to the scip problem
          SCIP_CALL( SCIPaddVar(lift, var) );
          vars_it+=1;
        }
      }
    }
 
    size_L_wo_N_k=vars_it;  // size of L without N(k)
          
    //get variables
    SCIP_VAR** vars_lift=SCIPgetVars(lift); //vars_lift = variables of the lifting problem
    

    /////////////////////////////////
    //create the linear constraint //
    /////////////////////////////////

    SCIP_CONS * cons_lift;
    namebuf.str("");
    namebuf<<"row_"<<1;

    // create SCIP_CONS object
    SCIP_CALL( SCIPcreateConsBasicLinear(lift, &cons_lift, namebuf.str().c_str(), 0, NULL, NULL, -SCIPinfinity(lift), b) );
          
    // add the vars to the constraint
    for (int j = 0; j < size_L_wo_N_k; ++j)
    {
      SCIP_CALL( SCIPaddCoefLinear(lift, cons_lift, vars_lift[j],a[L_wo_N_k.at(j)]) );
    }

    /////////////////////////////////
    // create the SOS1 constraints //
    /////////////////////////////////

    Graph::vertex_set S;
    for (int j = 0; j < size_L_wo_N_k; ++j)
    {
      S.insert(L_wo_N_k.at(j));
    }
    
    Graph sub_Gr=(*Gr_conf).subgraph(S);
    

  return SCIP_OKAY;
}


