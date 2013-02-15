

#define delta 0.001; //error tol

using namespace NGraph;

static
SCIP_RETCODE conflict_problem(SCIP* scip,         //SCIP data structure 
			      int n,              //number of variables
			      int numbercons,     //number of constraints
                              int d,              //constraint being considered
			      SCIP_VAR** vars,    //variables of scip
			      SCIP_CONS** cons,   //constraints of scip
                              double* x_star,     //solution of LP relax.
                              double* w,          //pointer to store solution of conflict problem
                              double* OBJVAL,     //pointer to store objective value of conflict problem
                              Graph* A            //pointer to store conflict graph of conflict problem
			      )
{
  ostringstream namebuf;

  ///////////////////////////////
  // get data of constraint d  //
  ///////////////////////////////

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
  SCIP_VAR *  vars_conf[n_vars_cons_d]; 
  vars_it=0;     //vars iterator
  int cons_it=0; //cons iterator
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

  cons_it=0; // cons iterator
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
        for (int j=0; j<n_vars_cons_d; j++)
        {    
          if(getvarssos[0]==vars_cons_d[j])
          {
            save_getvar1=j;
          }
          if(getvarssos[1]==vars_cons_d[j])
          {
            save_getvar2=j;
          }
          if(save_getvar1!=-1 && save_getvar2!=-1)
          {
            (*A).insert_undirected_edge(save_getvar1,save_getvar2);
            SCIP_CONS * sos1;
            namebuf.str("");
            namebuf<<"sosrow_"<<cons_it; 
            // create SCIP_CONS object
            SCIP_CALL( SCIPcreateConsBasicSOS1( conf, &sos1, namebuf.str().c_str(), 0, NULL, NULL ) ); 
            // add variables to the SOS1 constraint
            SCIP_CALL( SCIPappendVarSOS1(conf, sos1, vars_conf[save_getvar1]) );
            SCIP_CALL( SCIPappendVarSOS1(conf, sos1, vars_conf[save_getvar2]) );
            // add the constraint to scip
            SCIP_CALL( SCIPaddCons(conf, sos1) );
            cons_it+=1;
            break;
          }
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
  for (int i=0; i<n_vars_cons_d; i++)
  {
    w[i]=SCIPgetSolVal(conf,sol_conf,vars_conf[i]);
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


