
static
SCIP_RETCODE LP_relax(SCIP* scip,             //SCIP data structure
		      SCIP_CONS** cons,       // cons = constraints of scip
		      SCIP_VAR** vars,        //vars = variables of scip
		      int numbercons,         //numbercons = number of constraints
                      int n,                  //n=number of variables
                      SCIP_VARDATA* vardata,  //vector with indices of variables of obj
		      double* x_star          //pointer to store solution of LP rel.
		      )
{


  // create LP relax.
  SCIP* LP = NULL;
  SCIP_CALL( SCIPcreate(&LP) );
  SCIP_CALL( SCIPincludeDefaultPlugins(LP) );
  SCIP_CALL( SCIPcreateProbBasic(LP, "LP") );
  SCIP_CALL(SCIPsetObjsense(LP, SCIP_OBJSENSE_MAXIMIZE));

  SCIP_VAR* vars_LP[n]; //vector to store variables of LP rel.


  // create variables
  ostringstream namebuf;
  for (int i = 0; i < n; ++i)
  {
     SCIP_VAR * xvar;
     SCIP_Real c;
     namebuf.str("");
     namebuf << "x#" << i;

     // create the SCIP_VAR object
     c=SCIPvarGetObj(vars[i]);
     SCIP_CALL( SCIPcreateVarBasic(LP, &xvar, namebuf.str().c_str(), 0.0, 1.0, c, SCIP_VARTYPE_CONTINUOUS) );

     // add the SCIP_VAR object to the scip problem
     SCIP_CALL( SCIPaddVar(LP, xvar) );

     // storing for later access
     vars_LP[i] = xvar;
   }



  // create linear constraints
  int cons_it=0; //iterator
  for (int i=0; i<numbercons; i++)
  {
    if(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons[i])), "linear") == 0)
    {  
      SCIP_Real b= SCIPgetRhsLinear(scip,cons[i]);             //rhs of constraint i
      SCIP_Real* a= SCIPgetValsLinear(scip,cons[i]);           //values of constraint i
      SCIP_VAR** vars_cons_i=SCIPgetVarsLinear(scip,cons[i]);  //vars of constraint i
      int n_vars_cons_i=SCIPgetNVarsLinear(scip,cons[i]);      //number of variables of constraint i 
      SCIP_VARDATA *vardataout;

      SCIP_CONS * constr;
      namebuf.str("");
      namebuf<<"row_"<<cons_it;

      // create SCIP_CONS object
      SCIP_CALL( SCIPcreateConsBasicLinear(LP, &constr, namebuf.str().c_str(), 0, NULL, NULL, -SCIPinfinity(LP), b) );

      // add the vars belonging to field in this row to the constraint
      for (int j = 0; j < n_vars_cons_i; ++j)
      {
        vardataout=SCIPvarGetData(vars_cons_i[j]);
        SCIP_CALL( SCIPaddCoefLinear(LP, constr, vars_LP[(*vardataout).index], a[j]) );
      }

      // add the constraint to scip
      SCIP_CALL( SCIPaddCons(LP, constr) );

      cons_it+=1;

    }
  }

  /*
  // set var-type continuous
  for (int i=0; i<n_LP; i++)
  {
    SCIP_Bool infeasible; 
    SCIP_CALL(SCIPchgVarType(LP,vars_LP[i],SCIP_VARTYPE_CONTINUOUS,&infeasible));
  }
  */

  //get solution
  SCIP_CALL( SCIPsolve(LP) );
  SCIP_SOL* sol_LP = SCIPgetBestSol(LP);
  //SCIP_Real objval_LP= SCIPgetSolOrigObj(LP,sol_LP);
  
  for (int i=0; i<n; i++)
  {
    x_star[i]=SCIPgetSolVal(LP,sol_LP,vars_LP[i]);
  }


  //release LP
  /*
  SCIP_CONS** cons_LP=SCIPgetConss(LP);
  SCIP_VAR** varss_LP=SCIPgetVars(LP);
  int numbercons_LP=SCIPgetNConss(LP); //number of constraints
  for (int i=0; i<numbercons_LP; i++)
  {
    SCIP_CALL( SCIPreleaseCons(LP,&cons_LP[i] ) );
  }
  for (int i=0; i<n; i++)
  {
    SCIP_CALL( SCIPreleaseVar(LP,&varss_LP[i] ) );
  }
  */
  SCIP_CALL(SCIPfree(&LP) );
  
  return SCIP_OKAY;
}


