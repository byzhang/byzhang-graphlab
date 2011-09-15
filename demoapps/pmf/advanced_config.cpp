#include "pmf.h"
#include "../gabp/advanced_config.h"
extern problem_setup ps;
extern advanced_config ac;

void advanced_config::init_command_line_options(graphlab::command_line_options & clopts){

  clopts.attach_option("datafile", &datafile, datafile, "Input matrix/tensor");
  clopts.add_positional("datafile");

  clopts.attach_option("algorithm", &algorithm, algorithm, "ps.algorithm");
  clopts.add_positional("algorithm");
  
  clopts.attach_option("debug", &debug, debug, "Display debug output. (optional)");
  clopts.attach_option("float", &FLOAT, FLOAT, "is data in float format?");
  clopts.attach_option("D", &D, D, "number of features (dimension of computed weight vector)");
  clopts.attach_option("zero", &zero, zero, "support zero edges");  
  clopts.attach_option("scaling", &scaling, scaling, "time scaling factor (optional)");  
  clopts.attach_option("truncating", &truncating, truncating, "time truncation factor (optional)");  
  clopts.attach_option("savegraph", &savegraph, savegraph, "save graphs to file");  
  clopts.attach_option("loadgraph", &loadgraph, loadgraph, "load graphs to file");  
  clopts.attach_option("savefactors", &savefactors, savefactors, "save factors to file");  
  clopts.attach_option("loadfactors", &loadfactors, loadfactors, "load factors from file");  
  clopts.attach_option("stats", &stats, stats, "compute graph statistics");  
  clopts.attach_option("binaryoutput", &binaryoutput, binaryoutput, "export U,V,T to a binary file"); 
  clopts.attach_option("matrixmarket", &matrixmarket, matrixmarket, "give input in matrix market format"); 
 
  //BPTF related switches
  clopts.attach_option("bptf_alpha", &bptf_alpha, bptf_alpha, "BPTF alpha (noise parameter)");  
  clopts.attach_option("bptf_burn_in", &bptf_burn_in, bptf_burn_in, "BPTF burn-in period");
  clopts.attach_option("bptf_delay_alpha", &bptf_delay_alpha, bptf_delay_alpha, "BPTF start sampling alpha (noise level) the bptf_delay_alpha round ");  
  
  //ALS related switches
  clopts.attach_option("regnormal", &regnormal, regnormal, "ALS - use identical normalization for each variable? (default is weighted regularization by the number of edges");  
  clopts.attach_option("lambda", &als_lambda, als_lambda, "ALS regularization weight");  
  
  clopts.attach_option("scalerating", &scalerating, scalerating, "scale rating value ");  
  clopts.attach_option("aggregatevalidation", &aggregatevalidation, aggregatevalidation, "aggregate training and validation into one dataset ");  
  clopts.attach_option("maxval", &maxval, maxval, "maximal allowed value in matrix/tensor");
  clopts.attach_option("minval", &minval, minval, "minimal allowed value in matrix/tensor");
  clopts.attach_option("outputvalidation", &outputvalidation, outputvalidation, "output prediction on vadlidation data in kdd format");
  clopts.attach_option("unittest", &unittest, unittest, "unit testing. ");

  //SVD++ related switches
  clopts.attach_option("svdpp_step_dec", &svdpp_step_dec, svdpp_step_dec, "SVD++ step decrement ");
 
  //SGD related switches
  clopts.attach_option("sgd_lambda", &sgd_lambda, sgd_lambda, "SGD step size");
  clopts.attach_option("sgd_gamma", &sgd_gamma, sgd_gamma, "SGD starting step size");
  clopts.attach_option("sgd_step_dec", &sgd_step_dec, sgd_step_dec, "SGD step decrement");
 
  //SVD related switches
  clopts.attach_option("svd_iter", &svd_iter, svd_iter, "SVD iteration number"); 
  clopts.attach_option("printhighprecision", &printhighprecision, printhighprecision, "print RMSE output with high precision");

  //implicit rating options (see reference 10 in pmf.h)
  clopts.attach_option("implicitratingtype", &implicitratingtype, implicitratingtype, "type can be: user/item/uniform");
  clopts.attach_option("implicitratingpercentage", &implicitratingpercentage, implicitratingpercentage, " precentage of implicit added edges (0-100)");
  clopts.attach_option("implicitratingvalue", &implicitratingvalue, implicitratingvalue, "value for implicit negative ratings");
  clopts.attach_option("implicitratingweight", &implicitratingweight, implicitratingweight, "weight/time for implicit negative ratings");

  //sparsity enforcing priors (see reference 11 in pmf.h)
  clopts.attach_option("user_sparsity", &user_sparsity, user_sparsity, "user sparsity [0.5->1) (for L1 regularization for sparsity enforcing priors - run modes 10,11");
  clopts.attach_option("movie_sparsity", &movie_sparsity, movie_sparsity, "movie sparsity [0.5->1) (for L1 regularization for sparsity enforcing priors - run modes 11,12");
  clopts.attach_option("lasso_max_iter", &lasso_max_iter, lasso_max_iter, "max iter for lasso sparsity (run modes 10-12)");


}

void problem_setup::verify_setup(){
  switch(algorithm){
  // iterative matrix factorization using alternating least squares
  // or SVD ++
  case ALS_MATRIX:
  case ALS_SPARSE_USR_FACTOR:
  case ALS_SPARSE_USR_MOVIE_FACTORS:
  case ALS_SPARSE_MOVIE_FACTOR:
  case WEIGHTED_ALS:
  case SVD_PLUS_PLUS:
  case STOCHASTIC_GRADIENT_DESCENT:
  case LANCZOS:
  case NMF:
    tensor = false; BPTF = false;
    break;

    // MCMC tensor factorization
  case BPTF_TENSOR:
    // tensor factorization , allow for multiple edges between user and movie in different times
  case BPTF_TENSOR_MULT:
    tensor = true; BPTF = true;
   break;
    //MCMC matrix factorization
  case BPTF_MATRIX:
    tensor = false; BPTF = true;
    break;
   // tensor factorization
  case ALS_TENSOR_MULT:
    tensor = true; BPTF = false;
    break;
  default:
    assert(0);
  }


//INPUT SANITY CHECKS
#ifdef GL_NO_MCMC
  if (BPTF){
    logstream(LOG_ERROR) << "Can not run MCMC method with GL_NO_MCMC flag. Please comment flag on pmf.h and recompile\n";
    exit(1); 
 }
#endif

#ifdef GL_NO_MULT_EDGES
  if (algorithm == ALS_TENSOR_MULT || algorithm == BPTF_TENSOR_MULT){
    logstream(LOG_ERROR) << "Can not have support for multiple edges with GL_NO_MULT_EDGES flag. Please comment flag on pmf.h and recompile\n";
   exit(1);
  }
#endif  
#ifdef GL_SVD_PP
  if (algorithm != SVD_PLUS_PLUS){
    logstream(LOG_ERROR) << "Can not run required algorithm with GL_SVD_PP flag. Please comment flag on pmf.h and recompile\n";
    exit(1);
  }
#else
  if (algorithm == SVD_PLUS_PLUS){
    logstream(LOG_ERROR) << "Can not run required algorithm without GL_SVD_PP flag. Please define flag on pmf.h and recompile\n";
    exit(1);
  }
#endif

#ifndef GL_NO_MULT_EDGES
#ifdef GL_SVD_PP
   logstream(LOG_ERROR) << "When working with SVD++ (GL_SVD_PP is defined in pmf.h) you should also define GL_NO_MULT_EDGES\n";
   exit(1); 
#endif
#endif

  if (ac.bptf_delay_alpha != 0 && (algorithm != BPTF_TENSOR_MULT && algorithm != BPTF_TENSOR))
	logstream(LOG_WARNING) << "Delaying alpha (sampling of noise level) is ignored in non-MCMC methods" << std::endl;

  if (ac.bptf_burn_in != 10 && (algorithm != BPTF_TENSOR_MULT && algorithm != BPTF_TENSOR && algorithm != BPTF_MATRIX))
	logstream(LOG_WARNING) << "Markov chain burn in period is ignored in non-MCMC methods" << std::endl;

  if (ac.user_sparsity < 0.5 || ac.user_sparsity >= 1){
	logstream(LOG_ERROR) << "user_sparsity of factor matrix has to be in the range [0.5 1)" << std::endl;
        exit(1);
  }
  if (ac.movie_sparsity < 0.5 || ac.movie_sparsity >= 1){
	logstream(LOG_ERROR) << "movie_sparsity of factor matrix has to be in the range [0.5 1)" << std::endl;
        exit(1);
  }

}
