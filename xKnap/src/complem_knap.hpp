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

Graph compGraph;


void algtr(int,int,int,int,int*,vector<int>*,vector<vector<double> >*,vector<vector<double> >*,vector<vector<double> >*,vector<vector<double> >*,vector<int>*,vector<double>*,int);
int sol(int,vector<vector<double> >*,int,double);

static
SCIP_RETCODE complementarity_knapsack(SCIP* scip,double* x_sol,double* OBJVAL_sol,vector<double> * lhs_values,vector<double> * obj_values,double a_k)
{
  time_t tstart, tend;
  tstart = time(0);

  //SCIP_CALL(SCIPpresolve(scip)); 
  
  SCIP_CONS** cons=SCIPgetConss(scip);
  int numbercons=SCIPgetNConss(scip); //number of constraints
  SCIP_VAR** vars=SCIPgetVars(scip); 
  int indexlinearcons;
  for (int i=0; i<numbercons; i++)
  {
    if(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons[i])), "linear") == 0)
    {
      indexlinearcons=i;
      break;
    }
  }

  SCIP_VAR** varslin=SCIPgetVarsLinear(scip,cons[indexlinearcons]);
  int items=SCIPgetNVars(scip); // number of items
  int n=items+1; //(items + root node)
  SCIP_Real c=SCIPgetRhsLinear(scip,cons[indexlinearcons]);  //c=capacity



	
  //create vector of weights and vector of profits
  int linvarsnr=0;
  vector<double> w(n);
  vector<int> p(n);
  
  SCIP_Real* getvalslin= SCIPgetValsLinear(scip,cons[indexlinearcons]);
  for (int i=0; i<items; i++)
  {
    if(vars[i]==varslin[linvarsnr])
    {
      SCIP_Real vargetobj=SCIPvarGetObj(vars[i]);
      p.at(i+1)=-vargetobj;
      w.at(i+1)=getvalslin[linvarsnr];
      linvarsnr+=1;     
    }
    else
    {
      SCIP_Real vargetobj=SCIPvarGetObj(vars[i]);
      p.at(i+1)=-vargetobj;
      w.at(i+1)=0;
    }
  }


  //compute upper bound for objective value
  int P=0;
  for (int i=1; i<n; i++)
  {
    P+=p.at(i);
  }

  //compute maximum left hand side of constraint
  double W=0;
  for (int i=1; i<n; i++)
  {
    W+=w.at(i);
  }

  //define root node
  int r=0;
  p.at(0)=0;
  w.at(0)=W+1;


  //create conflict tree 
  for (int i=0; i<numbercons; i++)
  {
    if(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons[i])), "SOS1") == 0)
    {
      int numvarsSOS1=SCIPgetNVarsSOS1(scip,cons[i]);
      if(numvarsSOS1>2)
      {
        std::cout << "no tree structure in complementarity constraints" << std::endl;
        return SCIP_ERROR;
      }
      else if(numvarsSOS1==2)
      {
        SCIP_VAR** getvarssos=SCIPgetVarsSOS1(scip,cons[i]);
        int save_getvar1=-1;
        int save_getvar2=-1;
        for (int j=0; j<items; j++)
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
            compGraph.insert_undirected_edge(save_getvar1+1,save_getvar2+1);
            break;
          }
	}
      }
    }
  }

  /*
  double epsilon=0.01;
  int p[n];
  for (int j=0; j<n; j++)
  {
    p[j]=K*floor(p_double[j]/K);
  }
  */



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


  int mark_iterator;
  vector<int> mark_node(n,0); // store which node is marked
  mark_node.at(0)=1; // mark root node
  mark_iterator=1;
  
  int parentnode=-1; // root node has no parent node

  //compute Z,Y,XZ,XY
  algtr(n,r,r,parentnode,&mark_iterator,&mark_node,&Z,&Y,&XZ,&XY,&p,&w,P);

  int OBJVAL;
  double d;

  //compute objective value s.t. x\in {0,1}^n
  OBJVAL=sol(r,&Y,P,c);


  //compute solution vector s.t. x\in {0,1}^n
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
  for(int i=0;i<n-1;i++)
  {
    x_sol[i]=x[i+1];
  }

  
  //for each upper right hand side of linear constraint compute objective value
  vector<int> Data_rhs;

  for(int i=0;i<OBJVAL;i++)
  {
    double dw;
    dw=Y.at(i).at(r);
    if((dw<=c)&&(dw>=dw-a_k))
    {
      (*lhs_values).push_back(dw);
      (*obj_values).push_back(i);
    }
  }
  



  // delete [] Z;
  //delete [] Y;
  //delete [] XZ;
  //delete [] XY;
  

  tend = time(0); 
  cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< endl;
  // system("PAUSE");

  return SCIP_OKAY;
}




void algtr(int n,int r,int j,int parentnode,int* mark_iterator,vector<int>* mark_node,vector<vector<double> >* Z,vector<vector<double> >* Y,vector<vector<double> >* XZ,vector<vector<double> >* XY,vector<int>* p,vector<double>* w,int P)
{  
  (*Z).at((*p).at(j)).at(j)=(*w).at(j);
  (*XZ).at((*p).at(j)).at(j)=pow(2.0,j);
  (*Y).at(0).at(j)=0;
  
  if(j==r) // j=root node
  {
    for (int i=*mark_iterator; i<n; i++)
    {
      if((*mark_node).at(i)==0)
      {
        (*mark_node).at(i)=1;
        *mark_iterator+=1;
        compGraph.insert_undirected_edge(r,i);
        algtr(n,r,i,j,mark_iterator,mark_node,Z,Y,XZ,XY,p,w,P);
        //break;
      }
    }
  }
  
  Graph::vertex_set S = compGraph.out_neighbors(j);
  for (Graph::vertex_set::const_iterator t = S.begin(); t !=S.end(); t++)
  {
    if(*t!=parentnode)
    {
      (*mark_node).at(*t)=1;
      parentnode=j;
      if(j!=r)
      {
        algtr(n,r,*t,parentnode,mark_iterator,mark_node,Z,Y,XZ,XY,p,w,P);
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


int sol(int j,vector<vector<double> >* Y,int P,double c)
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
