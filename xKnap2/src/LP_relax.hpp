
static
SCIP_RETCODE LP_relax(SCIP* scip,         //SCIP data structure
		      double* x_star      //pointer to store solution of LP rel.
		      )
{

  // create LP relax.
  SCIP_Bool valid;
  valid=FALSE;
  SCIP* LP = NULL;
  SCIP_CALL( SCIPcreate(&LP) );
  SCIP_CALL( SCIPcopy(scip,LP,NULL,NULL,"LP",TRUE,FALSE,FALSE,&valid) );
  SCIP_CALL(SCIPsetMessagehdlr(scip,NULL)); //suppress messages

  // get LP relax. data
  SCIP_CONS** cons_LP=SCIPgetConss(LP); // cons = constraints of scip
  SCIP_VAR** vars_LP=SCIPgetVars(LP); //vars = variables of scip
  int numbercons_LP=SCIPgetNConss(LP); //numbercons_LP = number of constraints
  int n_LP=SCIPgetNVars(LP); // n = number of variables

  // delete all the SOS1 constraints
  for (int i=0; i<numbercons_LP; i++)
  {
    if(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons_LP[i])), "SOS1") == 0)
    {
      SCIP_CALL(SCIPdelCons(LP,cons_LP[i]));
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
  
  for (int i=0; i<n_LP; i++)
  {
    x_star[i]=SCIPgetSolVal(LP,sol_LP,vars_LP[i]);
  }

  SCIP_CALL(SCIPfree(&LP) );

  return SCIP_OKAY;
}


