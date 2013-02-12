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
#include "mergesort.hpp"
#include "scip/pub_var.h"          // SCIPvarGetObj
#include "scip/pub_cons.h"         // SCIPconsGetHdlr    

#define _INF std::numeric_limits<double>::infinity();
#define negative 0.0;

using namespace std;
using namespace NGraph;

Graph Gr;

static
SCIP_RETCODE complementarity_knapsack(SCIP* scip)
{


  //SCIP_CALL(SCIPpresolve(scip)); 
  
  SCIP_CONS** cons=SCIPgetConss(scip);
  SCIP_VAR** vars=SCIPgetVars(scip); 
  int numbercons=SCIPgetNConss(scip); //number of constraints
  int n=SCIPgetNVars(scip); // number of items

  double x_star[9]={0.333333,1,0,1,1,1,0,1,1};


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


  /*
  // get objective values
  SCIP_Real c[n];
  for (int i=0; i<n; i++)
  {
      c[i]=-SCIPvarGetObj(vars[i]);
  }
  */



  ostringstream namebuf;
  /*
  for (int d=0; d<numbercons; i++)
  {
    if(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons[i])), "linear") == 0)
    {

    }
  }
  */
  int d=0;
  const char* name_cons_d=SCIPconsGetName(cons[d]); 
 
    // create (4) - Tree conflict problem
    SCIP_VAR** varslin=SCIPgetVarsLinear(scip,cons[d]);
    SCIP_Real b= SCIPgetRhsLinear(scip,cons[d]);
    SCIP_Real* a= SCIPgetValsLinear(scip,cons[d]);
    int numvarslinear=SCIPgetNVarsLinear(scip,cons[d]);
    int posvarslin[numvarslinear];

    SCIP* conf = NULL;
    SCIP_CALL( SCIPcreate(&conf) );
    SCIP_CALL( SCIPincludeDefaultPlugins(conf) );
    SCIP_CALL( SCIPcreateProbBasic(conf, "conf") );
    SCIP_CALL(SCIPsetObjsense(conf, SCIP_OBJSENSE_MAXIMIZE));


    //create variables
    SCIP_VAR *  vars_conf[numvarslinear];
    int vars_it=0;
    for (int i=0; i<numvarslinear; i++)
    {
       SCIP_VAR * var_conf;
       namebuf.str("");
       namebuf << "x_conf#" << i;

       // create the SCIP_VAR object
       
       for (int j=vars_it; j<n; j++)
       {
         if(vars[vars_it]==varslin[i])
	 {
           posvarslin[i]=vars_it;
           if(x_star[vars_it]>=0.0001)
	   {
             SCIP_CALL( SCIPcreateVarBasic(conf, &var_conf, namebuf.str().c_str(), 0.0, 1.0 ,a[i], SCIP_VARTYPE_BINARY) );
	   }
           else       
	   {
             SCIP_Real y=negative;
             SCIP_CALL( SCIPcreateVarBasic(conf, &var_conf, namebuf.str().c_str(), 0.0, 1.0 ,y, SCIP_VARTYPE_BINARY) );
	   }
           vars_it+=1;
           break;
	 }
         vars_it+=1;
       }
       // add the SCIP_VAR object to the scip problem
       SCIP_CALL( SCIPaddVar(conf, var_conf) );

       // storing for later access
       vars_conf[i] = var_conf;
    }


    int r=0;
    Graph A;
    //create SOS1 constraints
    for (int i=0; i<numbercons; i++)
    {
      if(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons[i])), "SOS1") == 0)
      {
        int numvarsSOS1=SCIPgetNVarsSOS1(scip,cons[i]);
        SCIP_VAR** getvarssos=SCIPgetVarsSOS1(scip,cons[i]);
        int save_getvar1=-1;
        int save_getvar2=-1;
        if(numvarsSOS1==2)
	{
          for (int j=0; j<numvarslinear; j++)
          {    
            if(getvarssos[0]==varslin[j])
            {
              save_getvar1=j;
            }
            if(getvarssos[1]==varslin[j])
            {
              save_getvar2=j;
            }
            if(save_getvar1!=-1 && save_getvar2!=-1)
            {
              A.insert_undirected_edge(save_getvar1,save_getvar2);
              SCIP_CONS * sos1;
              namebuf.str("");
              namebuf<<"sosrow_"<<r; 
              r+=1;
              // create SCIP_CONS object
              SCIP_CALL( SCIPcreateConsBasicSOS1( conf, &sos1, namebuf.str().c_str(), 0, NULL, NULL ) ); 
              // add variables to the SOS1 constraint
              SCIP_CALL( SCIPappendVarSOS1(conf, sos1, vars_conf[save_getvar1]) );
              SCIP_CALL( SCIPappendVarSOS1(conf, sos1, vars_conf[save_getvar2]) );
              // add the constraint to scip
              SCIP_CALL( SCIPaddCons(conf, sos1) );
              break;
            }
	  }
        }
      }
    }
    
    /*
  SCIP_CONS** consk=SCIPgetConss(conf);
  int numberconsk=SCIPgetNConss(conf); //number of constraints
    
          for (int k=0; k<numberconsk; k++)
          {    
    int numvarsSOS1k=SCIPgetNVarsSOS1(conf,consk[k]);
    SCIP_VAR** getvarssos=SCIPgetVarsSOS1(conf,consk[k]);
    cout << endl << numberconsk;
    cout << endl << numvarsSOS1k;
    cout << endl << getvarssos[0] << "  "<< getvarssos[1];
	  } 
    */

    

    

    SCIP_CALL( SCIPsolve(conf) );
    SCIP_CALL( SCIPprintBestSol(conf, NULL, FALSE) );
    
    SCIP_SOL* sol_conf = SCIPgetBestSol(conf);
    SCIP_Real objval_conf= SCIPgetSolOrigObj(conf,sol_conf);


    if(objval_conf>b) // cover detected
    {
      cout << endl << "cover detected" << endl << endl;
      SCIP_Real x[numvarslinear];
      SCIP_Real x2[numvarslinear];
      SCIP_Real c[numvarslinear];
      SCIP_Real c2[numvarslinear];
      int save_ind[numvarslinear];
      double sum=0;
      int size_cover=0;
      int change;
      for (int i=0; i<numvarslinear; i++)
      {
        x[i]=SCIPgetSolVal(conf,sol_conf,vars_conf[i]);
        x2[i]=x[i];
        if(x[i]>0.5) //(x[i]==1)
	{
          size_cover+=1;
	}
        c[i]=SCIPvarGetObj(vars_conf[i]);
        c2[i]=c[i];
        sum+=x[i]*c[i];
        save_ind[i]=i;
      }
      double save_sum=sum;
      mergesort(c2,save_ind,0,numvarslinear-1);

      for (int i=numvarslinear-1; i>=0; i--)
      {
        if((x[save_ind[i]]<=0.5)&&(c[save_ind[i]]!=0.0)) //(x[save_ind[i]]==1)
	{
          x[save_ind[i]]=1;
          sum+=c[save_ind[i]];
          Graph::vertex_set S = A.out_neighbors(save_ind[i]);
          for (Graph::vertex_set::const_iterator t = S.begin(); t !=S.end(); t++)
          {
            if(x[*t]==1)
	    {
              x[*t]=0;
              sum-=c[*t];
	    }
          }  
	  //  const char* SCIPvarGetName();
          if(sum<b)
	  {
            change=save_ind[i]; //nue=posvarslin[change]
            break;
	  }
          else
	  {
            sum=save_sum;
	  }
	}
      }
      SCIP_Bool valid;
      valid=FALSE;
      SCIP* lift = NULL;
      SCIP_CALL( SCIPcreate(&lift) );
      SCIP_CALL( SCIPcopy(scip,lift,NULL,NULL,"lift",TRUE,FALSE,FALSE,&valid) );  
      SCIP_CALL(SCIPsetObjsense(lift, SCIP_OBJSENSE_MAXIMIZE));
      SCIP_CONS** cons_lift=SCIPgetConss(lift);
      SCIP_VAR** vars_lift=SCIPgetVars(lift);
      double cover[size_cover];
      // delete all the linear constraints except constraint d
      for (int i=0; i<numbercons; i++)
      {
        if((strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons_lift[i])), "linear") == 0)&&(strcmp(SCIPconsGetName(cons_lift[i]), name_cons_d) != 0))
        {
	  SCIP_CALL(SCIPdelCons(lift,cons_lift[i]));
        }
      }
      // define objective function
      SCIP_CALL( SCIPchgVarObj(lift,vars_lift[posvarslin[change]],sum-c[change])); //(posvarslin[change]=nue)
      int xy=0;
      int xy2=0;
      for (int i=0; i<n; i++)
      {
        if(i==posvarslin[xy])
	{
          if((x2[xy]>0.5)&&(c[xy]>0))
	  {
          SCIP_CALL( SCIPchgVarObj(lift,vars_lift[i],c[xy]));
          cover[xy2]=i;
          xy2+=1;
          }
          else
	  {
            SCIP_CALL( SCIPchgVarObj(lift,vars_lift[i],0.0));
	  }
          xy+=1;
	}
        else
	{
          SCIP_CALL( SCIPchgVarObj(lift,vars_lift[i],0.0));
	}
      }
      // set x_i=0 for all i \in {1,...,n}\cover
      int cover_iterator=0;
      for (int i=0; i<n; i++)
      {
        if(i!=cover[cover_iterator])
	{
          SCIP_CALL( SCIPchgVarUb(lift,vars_lift[i],0.0) );
	}
        else
	{
          cover_iterator+=1;
	}
      }
      int k=6;/////// lifting variable
      double dd=0.1;
      // set k==dd
      SCIP_CALL( SCIPchgVarUb(lift,vars_lift[k],dd) );
      SCIP_CALL( SCIPchgVarLb(lift,vars_lift[k],dd) );
      // set neighbors of k null
      Graph::vertex_set S1 = Gr.out_neighbors(k);
      for (Graph::vertex_set::const_iterator t = S1.begin(); t !=S1.end(); t++)
      {
	SCIP_CALL( SCIPchgVarUb(lift,vars_lift[*t],0.0) );
      }     
      SCIP_CALL( SCIPsolve(lift) );
      // SCIP_CALL( SCIPprintBestSol(lift, NULL, FALSE) );
      SCIP_SOL* sol_lift = SCIPgetBestSol(lift);
      SCIP_Real objval_lift= SCIPgetSolOrigObj(lift,sol_lift);
      cout << endl << "f[y]= " << b-objval_lift;

      
    } //end  if(objval_conf>b) 
   




  /*
  //for (int k=0; k<nlincons; k++)
  int k=0;
  int d=numbercons-1-k;
  
    // create (4) - Tree conflict problem
    SCIP_Bool valid;
    valid=FALSE;
    SCIP* conf = NULL;
    SCIP_CALL( SCIPcreate(&conf) );
    SCIP_CALL( SCIPcopy(scip,conf,NULL,NULL,"confx",TRUE,FALSE,FALSE,&valid) );
    SCIP_CALL(SCIPsetObjsense(conf, SCIP_OBJSENSE_MAXIMIZE));



    // get variables & get constraints
    SCIP_VAR** varsk=SCIPgetVars(conf);  
    SCIP_CONS** consk=SCIPgetConss(conf);


    // get values, rhs, vars and number of vars of linear constraint d
    SCIP_VAR** varslin=SCIPgetVarsLinear(conf,consk[d]);
    SCIP_Real b= SCIPgetRhsLinear(conf,consk[d]);
    SCIP_Real* a= SCIPgetValsLinear(conf,consk[d]);
    int numvarslinear=SCIPgetNVarsLinear(conf,consk[d]);



    // replace objective function with left hand side of constraint d
    int i_lin=0;
    for (int i=0; i<n; i++)
    {
      if(varsk[i]==varslin[i_lin])
      {
        SCIP_CALL( SCIPchgVarObj(conf,varsk[i],a[i_lin]));
        i_lin+=1;
      }
      else
      {
        SCIP_Bool deleted;
        SCIP_CALL(SCIPdelVar(conf,varsk[i],&deleted));
      }
    }



    
    SCIP_Real c_neu[n];
    for (int i=0; i<numvarslinear; i++)
    {
      c_neu[i]=-SCIPvarGetObj(varslin[i]);
      cout << endl << c_neu[i];
    } 
    


    //delete some SOS1 constraints
    for (int i=0; i<numbercons; i++)
    {
      if(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(consk[i])), "SOS1") == 0)
      {
        int numvarsSOS1=SCIPgetNVarsSOS1(conf,consk[i]);
        SCIP_VAR** getvarssos=SCIPgetVarsSOS1(conf,consk[i]);
        int save_getvar1=-1;
        int save_getvar2=-1;
        if(numvarsSOS1==2)
	{
          for (int j=0; j<numvarslinear; j++)
          {    
            if(getvarssos[0]==varslin[j])
            {
              save_getvar1=j;
            }
            if(getvarssos[1]==varslin[j])
            {
              save_getvar2=j;
            }
            if(save_getvar1!=-1 && save_getvar2!=-1)
            {
              A.insert_undirected_edge(save_getvar1,save_getvar2);
            }
            if(j==numvarslinear-1)
	    {
              if(save_getvar1==-1 && save_getvar2==-1)
              {
                SCIP_CALL(SCIPdelCons(conf,consk[i]));
              }
	    }
	  }
        }
      }
      else if(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(consk[i])), "linear") == 0)
      {
        SCIP_CALL(SCIPdelCons(conf,consk[i]));
      }
    }






    SCIP_VAR** varsk_new=SCIPgetVars(conf); 
    SCIP_CALL( SCIPsolve(conf) );
    SCIP_CALL( SCIPprintBestSol(conf, NULL, FALSE) );
    SCIP_SOL* sol_conf = SCIPgetBestSol(conf);
    SCIP_Real objval_conf= SCIPgetSolOrigObj(conf,sol_conf);
    SCIP_Real x_conf[numvarslinear];


    if(objval_conf>b) // cover detected
    {
      cout << endl << "cover detected";
      for (int i=0; i<n; i++)
      {
        x_conf[i]=SCIPgetSolVal(conf,sol_conf,varsk_new[i]);
	cout << endl << i << "  "<< x_conf[i];
      } 
    }
*/

   return SCIP_OKAY;
}
