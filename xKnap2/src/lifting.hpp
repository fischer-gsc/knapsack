#include <set>
#define _INF std::numeric_limits<double>::infinity();

using namespace std;
using namespace NGraph;

static
SCIP_RETCODE lifting(
		     int n_vars_cons_d,      //number of variables of constraint d
                     SCIP_Real* a,           //values of constraint d
                     SCIP_Real b,            //right hand side of constraint d
		     Graph* Gr_conf,         //pointer to store conflict graph of conflict problem
		     vector<int>* L,         //set L=C without N(nue) with nue
                     int* I,                 //set I=complement of L
		     int size_L,             //size of L
		     int size_I              //size of I
		     )
{
  ostringstream namebuf;


  //get copy of a
  double alpha[n_vars_cons_d];
  for (int j = 0; j < n_vars_cons_d; ++j)
  {
    alpha[j]=a[j];
  }


  //for every k \in I :  lift
  for (int f = 0; f < size_I; ++f)
  {
    int k=I[f];        //k=lifting variable
    double a_k=a[k];  //coefficient of lifting variable
    

    /////////////////////////////////////////////////////////
    // get neigbors of k and sort them in increasing order //
    /////////////////////////////////////////////////////////

    int deg_k = (*Gr_conf).degree(k);                       //degree of node k
    Graph::vertex_set N_k= (*Gr_conf).out_neighbors(k);     //neighbors of node k
    vector<double> N_k_sort(deg_k);                         //copy of N_k
    vector<int> nullvec(deg_k,0);                           //no need for mergesort  
    int it=0;                                               //iterator
    for (Graph::vertex_set::const_iterator t = N_k.begin(); t !=N_k.end(); t++)
    {
      N_k_sort.at(it)=*t;
      it+=1;
    }

    if(deg_k>1)
    {
      mergesort(&N_k_sort,&nullvec,0,deg_k-1);  // sort N_k_sort
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
      while((*L).at(i)>N_k_sort.at(it))
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
        if((*L).at(i)!=N_k_sort.at(it))
	{
          //add L[i] to L without N(k)
          L_wo_N_k.push_back((*L).at(i));
          vars_it+=1;
	}
      }
      else
      {
        if((*L).at(i)<N_k_sort.at(it))
	{
          //add L[i] to L without N(k)
          L_wo_N_k.push_back((*L).at(i));
          vars_it+=1;
        }
      }
    }
 
    size_L_wo_N_k=vars_it;  // size of L without N(k)


    ////////////////////////////////////////////////////
    // create vector of weights and vector of profits //
    ////////////////////////////////////////////////////

    vector<double> w;     // vector of weights
    vector<int> p;        // vector of profits

    for (int j = 0; j < size_L_wo_N_k; ++j)
    {
      p.push_back(alpha[L_wo_N_k.at(j)]);
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
    


    //we have to change the indices of the nodes
    
    //Example: {2,5,6} --> {0,1,2}
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
          compGraph.insert_edge(i,ind_in_d[*t]);
        }
      }
    }


    /////////////////////////////
    // get lifting coefficient //
    /////////////////////////////


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
    double alpha_k;              // alpha_k=lifting coefficient of x[k]
    double alpha_k_save=_INF;

    for(int j=0;j<size_lhs_values;j++)
    {
        cout << endl << "lhs_values " << lhs_values[j];
        cout << endl << "obj_values " << obj_values[j];
      y=(b-lhs_values[j])/a_k;  // remark: a_k>0 always fulfilled
      if(y>1)
      {
        y=1;
      }
      if(y>0)
      {
        f_y=b-obj_values[j];
        cout << endl << "y= " << y;
	cout << endl << "f_y= " << f_y;
        alpha_k=f_y/y;          // y != 0
        if(alpha_k<alpha_k_save)
        {
          alpha_k_save=alpha_k;
        }
      }
    }
    alpha_k=alpha_k_save;


    cout << endl << "L\N(k)= " <<endl;
    for(int j=0;j<size_L_wo_N_k;j++)
    {
      cout << endl << L_wo_N_k.at(j);
    }
    cout << endl;

    cout << endl << "k=  "<< k << "  " << "alpha_k=  " << alpha_k << endl;


    /////////////////////////
    // update L and size_L //
    /////////////////////////

    size_L+=1;
    (*L).push_back(k); //remark: L sorted without regard to last entry
    //sort L
    for (int i=0; i<size_L-1; i++)
    {
      if((*L).at(i)>k)
      {
        for (int j=size_L-2; j>=i; j--)
        {
          (*L).at(j+1)=(*L).at(j);
	}
        (*L).at(i)=k;
        break;
      }
    }


    //////////////////
    // update alpha //
    //////////////////

    alpha[k]=alpha_k;


    ///////////////
    //clear data //
    ///////////////

    delete [] x_knap;
    delete [] ind_in_d;
    lhs_values.clear ();
    obj_values.clear ();
}

  return SCIP_OKAY;
}


