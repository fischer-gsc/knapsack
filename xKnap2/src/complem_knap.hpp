#include <cstdlib>
#include <sstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <limits>
#include <ctime>
#include <vector>
#include "ngraph/ngraph.hpp"
#include "scip/pub_var.h"          // SCIPvarGetObj   

#define _INF std::numeric_limits<double>::infinity();

using namespace std;
using namespace NGraph;


void algtr(int,int,int,int,int*,vector<int>*,vector<vector<double> >*,vector<vector<double> >*,vector<vector<double> >*,vector<vector<double> >*,vector<int>*,vector<double>*,int,Graph* compGraph);
int sol(int,vector<vector<double> >*,int,double);

static
SCIP_RETCODE complementarity_knapsack(
				       Graph* compGraph,              //conflict graph (forest consisting of tree subgraphs)
				       vector<double>* w,             //vector of weights
				       vector<int>* p,                //vector of profits
                                       double* x_sol,                 //pointer to store solution
                                       double* OBJVAL_sol,            //pointer to store objective value
                                       vector<double> * lhs_values,   //pointer to store lhs-values of linear constraint
                                                                      //for every rhs<=c  
                                       vector<double> * obj_values,   //pointer to store obj.-values for every rhs<=c
                                       SCIP_Real c,                   //right hand side of constraint d
                                       double a_k,                    //coefficient of lifting index
                                       int items                      //number of items
                                     )
{
  time_t tstart, tend;
  tstart = time(0);


  int n=items+1;
	
  //compute upper bound for objective value
  int P=0;
  for (int i=0; i<items; i++)
  {
    P+=(*p).at(i);
  }

  //compute maximum left hand side of constraint
  double W=0;
  for (int i=0; i<items; i++)
  {
    W+=(*w).at(i);
  }

  //define root node
  int r=items;
  (*p).push_back(0);
  (*w).push_back(W+1);



  //create matrix Z
  //Z[i][j]=solution with minimal weight found in the subtree T(j) that leads to a profit of i
  //with i included into the knapsack solution
  vector<vector<double> > Z;  
  Z = vector<vector<double> > (P+1);
  for(int i = 0; i < P+1; i++)
  {
    Z[i] = vector<double>(n,W+1);
  }

  //create matrix Y
  //Y[i][j]=solution with minimal weight found in the subtree T(j) that leads to a profit of i
  //with i excluded from the the knapsack solution
  vector<vector<double> > Y;  
  Y = vector<vector<double> > (P+1);
  for(int i = 0; i < P+1; i++)
  {
    Y[i] = vector<double>(n,W+1);
  }
  
  //create matrix XZ
  //XY[i][j]=to Y[i][j] associated solution vector (binary encoding)
  vector<vector<double> > XZ;  
  XZ = vector<vector<double> > (P+1);
  for(int i = 0; i < P+1; i++)
  {
    XZ[i] = vector<double>(n,0);
  }
  
  //create matrix XY
  //XY[i][j]=to Z[i][j] associated solution vector (binary encoding)
  vector<vector<double> > XY;  
  XY = vector<vector<double> > (P+1);
  for(int i = 0; i < P+1; i++)
  {
    XY[i] = vector<double>(n,0);
  }


  vector<int> mark_node(n,0); // store which node is marked
  mark_node.at(r)=1; // mark root node
  
  int parentnode=-1;    // root node has no parent node
  int mark_iterator=0;  // number of visited nodes (root node excluded)

  //compute Z,Y,XZ,XY
  algtr(n,r,r,parentnode,&mark_iterator,&mark_node,&Z,&Y,&XZ,&XY,p,w,P,compGraph);

  int OBJVAL;
  double d;

  //compute objective value
  OBJVAL=sol(r,&Y,P,c);

  //compute solution vector
  double X;
  X=XY.at(OBJVAL).at(r);
  d=Y.at(OBJVAL).at(r);

  int x[n];
  for(int i=n-1;i>=0;i--)
  {
    double a;
    a=X-pow(2.0,i);
    if(a>=0)
    {
      X=a;
      x[i]=1;
    }
    else
    {
      x[i]=0;
    }
  }

  *OBJVAL_sol=OBJVAL;
  for(int i=0;i<items;i++)
  {
    x_sol[i]=x[i];
  }

  
  //for each upper right hand side of linear constraint compute objective value
  vector<int> Data_rhs;

  for(int i=0;i<=OBJVAL;i++)
  {
    double dw; //left hand side of linear constraint
    dw=Y.at(i).at(r);
    if((dw<=c)&&(dw>=c-a_k))    // c-1*a_k <= lhs <=c- 0*a_k
    {
      (*lhs_values).push_back(dw);
      (*obj_values).push_back(i);
    }
  }
  

  Z.clear ();
  Y.clear ();
  XZ.clear ();
  XY.clear ();
  

  tend = time(0); 
  cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< endl;

  return SCIP_OKAY;
}




void algtr(int n,int r,int j,int parentnode,int* mark_iterator,vector<int>* mark_node,vector<vector<double> >* Z,vector<vector<double> >* Y,vector<vector<double> >* XZ,vector<vector<double> >* XY,vector<int>* p,vector<double>* w,int P,Graph* compGraph)
{  
  (*Z).at((*p).at(j)).at(j)=(*w).at(j);
  (*XZ).at((*p).at(j)).at(j)=pow(2.0,j);
  (*Y).at(0).at(j)=0;
  
  if(j==r) // j=root node
  {
    for (int i=*mark_iterator; i<n-1; i++)  //(n-1=items)
    {
      if((*mark_node).at(i)==0)
      {
        (*mark_node).at(i)=1;
        *mark_iterator+=1;
        (*compGraph).insert_edge(r,i);
        (*compGraph).insert_edge(i,r);
        algtr(n,r,i,j,mark_iterator,mark_node,Z,Y,XZ,XY,p,w,P,compGraph);
        //break;
      }
    }
  }
  
  Graph::vertex_set S = (*compGraph).out_neighbors(j);
  for (Graph::vertex_set::const_iterator t = S.begin(); t !=S.end(); t++)
  {
    if(*t!=parentnode)
    {
      (*mark_node).at(*t)=1;
      parentnode=j;
      if(j!=r)
      {
        algtr(n,r,*t,parentnode,mark_iterator,mark_node,Z,Y,XZ,XY,p,w,P,compGraph);
      }
      vector<double> save_vector_1(P+1);
      vector<double> save_vector_2(P+1);
      vector<double> save_vectorXY(P+1);
      vector<double> save_vectorXZ(P+1);
      for (int nue=0; nue<P+1; nue++)
      {
        save_vector_1.at(nue)=(*Y).at(nue).at(j);
        save_vector_2.at(nue)=(*Z).at(nue).at(j);
        save_vectorXY.at(nue)=(*XY).at(nue).at(j);
        save_vectorXZ.at(nue)=(*XZ).at(nue).at(j);
      }

      for (int d=0; d<P+1; d++)
      {
        int save_k=0;
        int _isZ=0;
        double v=_INF;
        for (int k=0; k<d+1; k++)
        {
          double u;
          int isZ;
          if((*Y).at(k).at(*t)<=(*Z).at(k).at(*t))
          {
            u=(*Y).at(k).at(*t);
            isZ=0;
          }
          else
          {
            u=(*Z).at(k).at(*t);
            isZ=1;
          }
          u+=save_vector_1.at(d-k);
          if(v>u)
          {
            v=u;
            save_k=k;
            _isZ=isZ;
          }
        }
        (*Y).at(d).at(j)=v;
        if(_isZ==0)
        {
          (*XY).at(d).at(j)=save_vectorXY.at(d-save_k)+(*XY).at(save_k).at(*t);
        }
        else
        {
          (*XY).at(d).at(j)=save_vectorXY.at(d-save_k)+(*XZ).at(save_k).at(*t);
        }
      }
      for (int d=(*p).at(j); d<P+1; d++)
      {
        int save_k=0;
        double v=_INF;
        for (int k=0; k<d-(*p).at(j)+1; k++)
        {
          double u;
          u=(*Y).at(k).at(*t);
          u+=save_vector_2.at(d-k);
          if(v>u)
          {
            v=u;
            save_k=k;
          }
        }
        (*Z).at(d).at(j)=v;
        (*XZ).at(d).at(j)=save_vectorXZ.at(d-save_k)+(*XY).at(save_k).at(*t);
      }
    }
  }
}


int sol(int j,                         // root node j=0
        vector<vector<double> >* Y,    // 
        int P,                         // P ^= 1-norm of objective vector 
        double c                       // right hand side of the linear constraint
        )
{
  int OBJVAL=0;
  for(int i=0;i<P+1;i++)
  {
    double u;
    u=(*Y).at(i).at(j);
    if(u<=c)
    {
      OBJVAL=i;
    }
  }
  return OBJVAL;
}
