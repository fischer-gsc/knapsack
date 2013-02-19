#include "scip/pub_var.h"

#define delta 0.001; //error tol

using namespace std;
using namespace NGraph;


static
SCIP_RETCODE conflict_problem(SCIP* scip,             //SCIP data structure 
			      int n,                  //number of variables
			      int numbercons,         //number of constraints
                              int d,                  //constraint being considered
                              int n_vars_cons_d,      //number of variables of constraint d
                              SCIP_VAR** vars_cons_d, //vars of constraint d
                              SCIP_Real* a,           //values of constraint d
                              int* posvars_obj_d,     //position of variables of obj in cons d
                                                      //( -1 if variable does not exist )
                              SCIP_VARDATA* vardata,  //vector with indices of variables of obj
			      SCIP_VAR** vars,        //variables of scip
                              SCIP_CONS** cons,       //constraints of scip
                              double* x_star,         //solution of LP relax.
                              double* w,              //pointer to store solution of conflict problem
                              double* OBJVAL,         //pointer to store objective value of conflict problem
                              int* null_norm_w,       //pointer to store the null (quasi) norm of w   
                              Graph* A                //pointer to store conflict graph of conflict problem
			      )
{
  ostringstream namebuf;


  //////////////////////////////
  // create conflict problem  //
  //////////////////////////////

  SCIP* conf = NULL;
  SCIP_CALL( SCIPcreate(&conf) );
  SCIP_CALL( SCIPincludeDefaultPlugins(conf) );
  SCIP_CALL( SCIPcreateProbBasic(conf, "conf") );
  SCIP_CALL(SCIPsetObjsense(conf, SCIP_OBJSENSE_MAXIMIZE));

  //////////////////////////////////////////////////
  // create the variables of the conflict problem //
  //////////////////////////////////////////////////

  SCIP_VAR* vars_conf[n_vars_cons_d]; //vector to store variables
  int vars_it=0;     //vars iterator
  int cons_it=0;     //cons iterator
  for (int i=0; i<n_vars_cons_d; i++)
  {
     SCIP_VAR * var_conf;
     namebuf.str("");
     namebuf << "w_conf#" << i;

     // create the SCIP_VAR object
    for (int j=vars_it; j<n; j++)
    {
      if(vars[vars_it]==vars_cons_d[i])
      {
        SCIP_CALL( SCIPcreateVarBasic(conf, &var_conf, namebuf.str().c_str(), 0.0, 1.0 ,a[i], SCIP_VARTYPE_BINARY) );
        double error_tol=delta;
        if(x_star[vars_it]<error_tol)
        {
          /////////////////////////////////////////
	  // add constraint var_conf[vars_it]==0 //
          /////////////////////////////////////////
          SCIP_CONS * cons_conf;
          namebuf.str("");
          namebuf<<"row_"<<cons_it;
          SCIP_CALL( SCIPcreateConsBasicLinear(conf, &cons_conf, namebuf.str().c_str(), 0, NULL, NULL, 0.0 , 0.0 ) );
          SCIP_CALL( SCIPaddCoefLinear(conf, cons_conf, var_conf, 1.0 ) );
          SCIP_CALL( SCIPaddCons(conf, cons_conf) );
	  cons_it+=1;
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

  ////////////////////////////////////////////////////////////////////////////////
  // create the SOS1 constraints and the conflict graph of the conflict problem //
  ////////////////////////////////////////////////////////////////////////////////

  int numvarsSOS1;
  SCIP_VAR** getvarssos;
  SCIP_VARDATA *vardataout0;
  SCIP_VARDATA *vardataout1;
  int node1;
  int node2;

  cons_it=0; // cons iterator
  //create SOS1 constraints
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
        node1=posvars_obj_d[(*vardataout0).index];
        node2=posvars_obj_d[(*vardataout1).index];
        if((node1!=-1)&&(node2!=-1))
	{
        //add edge to conflict graph A
        (*A).insert_edge(node1,node2);
        (*A).insert_edge(node2,node1);
        //create constraint
        SCIP_CONS * sos1;
        namebuf.str("");
        namebuf<<"sosrow_"<<cons_it; 
        // create SCIP_CONS object
        SCIP_CALL( SCIPcreateConsBasicSOS1( conf, &sos1, namebuf.str().c_str(), 0, NULL, NULL ) ); 
        // add variables to the SOS1 constraint
	SCIP_CALL( SCIPaddVarSOS1(conf, sos1, vars_conf[node1],0.5) );
        SCIP_CALL( SCIPaddVarSOS1(conf, sos1, vars_conf[node2],0.5) );
        // add the constraint to scip
        SCIP_CALL( SCIPaddCons(conf, sos1) );
        cons_it+=1;
        }
      }
    }
  }

  /////////////////////////////
  // solve conflict problem  //
  /////////////////////////////

  SCIP_CALL( SCIPsolve(conf) );  
  SCIP_SOL* sol_conf = SCIPgetBestSol(conf); 
  *OBJVAL= SCIPgetSolOrigObj(conf,sol_conf);
  
  *null_norm_w=0;
  for (int i=0; i<n_vars_cons_d; i++)
  {
    w[i]=SCIPgetSolVal(conf,sol_conf,vars_conf[i]);
    if(w[i]>0.5)  //(w[i]==1)
    {
      *null_norm_w+=1;
    }
  }


  /////////////////////////////
  //  free conflict problem  //
  /////////////////////////////

  SCIP_CONS** cons_conf=SCIPgetConss(conf);
  int numbercons_conf=SCIPgetNConss(conf); //number of constraints
  for (int i=0; i<numbercons_conf; i++)
  {
    SCIP_CALL( SCIPreleaseCons(conf,&cons_conf[i] ) );
  }
  for (int i=0; i<n_vars_cons_d; i++)
  {
    SCIP_CALL( SCIPreleaseVar(conf,&vars_conf[i] ) );
  }
  SCIP_CALL(SCIPfree(&conf) );

  return SCIP_OKAY;
}


