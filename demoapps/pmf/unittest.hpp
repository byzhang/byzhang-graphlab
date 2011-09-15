#ifndef _UNITTEST_HPP
#define _UNITTEST_HPP

extern advanced_config ac;

void test_cosamp();

/**
 * UNIT TESTING
 */
void verify_result(double obj, double train_rmse, double validation_rmse){
   assert(ac.unittest > 0);
   switch(ac.unittest){
      case 1: //ALS: Final result. Obj=0.0114447, TRAIN RMSE= 0.0033 VALIDATION RMSE= 1.1005.
	 assert(pow(obj - 0.012,2) < 1e-3);
         assert(pow(train_rmse - 0.0033,2) < 1e-3);
	 assert(pow(validation_rmse - 1.1005,2) < 1e-2);
	 break;

      case 91: //WEIGHTED_ALS: -Iter100... UV. objective=0.0207271, RMSE=0.0043/0.6344. Time to finish=0.00hr.
         assert(pow(obj -  0.0207271,2)<1e-5);
         assert(pow(train_rmse - 0.0043,2)<1e-5);
         assert(pow(validation_rmse - 0.6344,2)<1e-3);
         break;

   }
}


//UNIT TESTING OF VARIOUS METHODS
void unit_testing(int unittest, command_line_options& clopts){

   if (unittest == 1){ //ALTERNATING LEAST SQUARES
      ac.datafile = "als"; ps.algorithm = ALS_MATRIX; ac.algorithm = ALS_MATRIX; ac.FLOAT = true; ac.als_lambda = 0.001;
      clopts.set_scheduler_type("round_robin(max_iterations=100,block_size=1)");
   }

   else if (unittest == 91){ //WEIGHTED ALTERNATING LEAST SQUARES
      ac.datafile = "wals"; ac.algorithm = WEIGHTED_ALS; ps.algorithm = WEIGHTED_ALS; ac.FLOAT = true; ac.als_lambda = 0.001;
      clopts.set_scheduler_type("round_robin(max_iterations=100,block_size=1)");
   }
   else if (unittest == 101){
      test_cosamp();
      exit(0);
   } 
   else {
      logstream(LOG_ERROR) << " Unit test mode " << unittest << " is not supported yet! " << std::endl;
      exit(1);
   }
}



#endif
