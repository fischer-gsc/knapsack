
#define _INF std::numeric_limits<double>::infinity();

using namespace std;
using namespace NGraph;

static
SCIP_RETCODE lifting(
		     int n_vars_cons_d,      //number of variables of constraint d
                     SCIP_Real* a,           //values of constraint d
                     SCIP_Real b,            //right hand side of constraint d
		     Graph* Gr_conf,         //pointer to store conflict graph of conflict problem
		     int* L,                 //set L=C without N(nue) with nue
                     int* I,                 //set I=complement of L
		     int size_L              //size of L
		     )
{
  ostringstream namebuf;


  int k=I[0];            //(lifting variable)
  //k=1;
    double a_k=a[k];  //coefficient of lifting variable


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

    ///////////////////////////////
    // create set L without N(k) //
    ///////////////////////////////

    vector<int> L_wo_N_k;   // L without N(k)
    int size_L_wo_N_k;      // size of L without N(k)

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
          //add L[i] to L without N(k)
          L_wo_N_k.push_back(L[i]);
	}
      }
      else
      {
        if(L[i]<N_k_copy.at(it))
	{
          //add L[i] to L without N(k)
          L_wo_N_k.push_back(L[i]);

          vars_it+=1;
        }
      }
    }
 
    size_L_wo_N_k=vars_it;  // size of L without N(k)
          
    

    ///////////////////////////////////////////////////
    //create vector of weights and vector of profits //
    ///////////////////////////////////////////////////

    vector<double> w;     // vector of weights
    vector<int> p;        // vector of profits

    for (int j = 0; j < size_L_wo_N_k; ++j)
    {
      p.push_back(a[L_wo_N_k.at(j)]);
      w.push_back(a[L_wo_N_k.at(j)]);
    }

    
    //////////////////////////////////////////////////////
    // create the conflict Graph of the lifting problem //
    //////////////////////////////////////////////////////

    // remark:
    // conflict graph of lifting problem ^= subgraph of conflict graph
    // of conflict problem

    Graph::vertex_set S;
    for (int j = 0; j < size_L_wo_N_k; ++j)
    {
      Graph::vertex vert=L_wo_N_k.at(j);
      S.insert(vert);
    }
    
    Graph sub_Gr=(*Gr_conf).subgraph(S);
    


    //we have to change indices of vertex nodes
    
    //Example: {2,5,6} -> {0,1,2}
    //define ind_ind_d with
    //ind_ind_d[2]=0
    //ind_ind_d[5]=1
    //ind_ind_d[6]=2
    //ind_ind_d[i]=-1 if i not \in {2,5,6}
    
    int *ind_in_d;
    ind_in_d = new int[n_vars_cons_d];
    it=0;
    for (int j = 0; j < n_vars_cons_d; ++j)
    {
      if(j==L_wo_N_k[it])
      {
        ind_in_d[j]=it;
        if(it<size_L_wo_N_k-1)
	{
          it+=1;
	}
      }
      else
      {
        ind_in_d[j]=-1;
      }
    }

    //create conflict graph with changed vertex indices 

    Graph compGraph;
    for (int i = 0;  i< size_L_wo_N_k; ++i)
    {
      int node=L_wo_N_k.at(i);
      if(sub_Gr.degree(node)!=0)
      {
        Graph::vertex_set Sout = sub_Gr.out_neighbors(node);
        for (Graph::vertex_set::const_iterator t = Sout.begin(); t !=Sout.end(); t++)
        {
          cout << endl << node << "   xxx   " << *t;
          compGraph.insert_edge(i,ind_in_d[*t]);
        }
      }
    }


    double OBJVAL_knap;          //objective value of knapsack problem
    int items=size_L_wo_N_k;     //number of items
    double *x_knap;              //solution vector of knapsack problem
    x_knap = new double[items];
    vector<double> lhs_values;   //vector to store lhs-values of linear 
                                 //constraint for every rhs<=c 
    vector<double> obj_values;   //vector to store obj.-values for every rhs<=c


    //solve knapsack problem s.t. x \in {0,1}^n
    //in addition, get all solutions for smaller rhs
    SCIP_CALL(complementarity_knapsack(&compGraph,&w,&p,x_knap,&OBJVAL_knap,&lhs_values,&obj_values,b,a_k,items));



    int size_lhs_values= static_cast<int>(lhs_values.size ());  // size of vector lhs_values
    double y;                    // y=a_k * x[k]
    double f_y;                  // f(y)=b-z(y)     ( z(y)=sol of knapsack problem with rhs=b-a_k*y ) 
    double alpha;                // alpha=lifting coefficient of x[k]
    double alpha_save=_INF;
    cout << endl << size_lhs_values;
    for(int j=0;j<size_lhs_values;j++)
    {
      y=(b-lhs_values[j])/a_k;
      if(y>1)
      {
        y=1;
      }
      f_y=b-obj_values[j];
      cout << endl << y << endl <<"  "<< f_y;
      alpha=f_y/y;              // y != 0 ??????
      if(alpha<alpha_save)
      {
        alpha_save=alpha;
      }
    }
    alpha=alpha_save;

    cout << endl << alpha;



    delete [] x_knap;
    delete [] ind_in_d;
    lhs_values.clear ();
    obj_values.clear ();


  return SCIP_OKAY;
}


