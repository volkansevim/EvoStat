#define T0    RCONST(0.0)      /* initial time           */

//*************************************************************************************************
double DiffSq(VecofDs timeseries, double desired, unsigned todiscard){
	unsigned size = timeseries.size();
	double diffsq = 0.;

	if(todiscard > size) MyErr("todiscard > size");
	for(unsigned i=todiscard; i<size; i++){
		diffsq += SQR(timeseries[i] - desired);
	}
	return diffsq/(size-todiscard);
}
//*************************************************************************************************
// Construction of the state change matrix. Each row i contains the changes in each reactant j 
// after firing a reaction i.
void ConstructStateChMx(Reactor &indiv, MatrixInt & StateChMx){
	size_t ReactantCount = indiv.Reactants.size();
	VecofInts row;
	row.resize(ReactantCount);
	StateChMx.clear();
	
	for(size_t i=0; i < indiv.Reactions.size(); i++) {
		int a = indiv.Reactions[i].IN[0];
		int b = indiv.Reactions[i].IN[1];
		int c = indiv.Reactions[i].OUT[0];
		int d = indiv.Reactions[i].OUT[1];
		foreach(int &a, row) a=0; //zero row out.
		
		if(a>=0) row[a]-=1;
		if(b>=0) row[b]-=1;
		if(c>=0) row[c]+=1;
		if(d>=0) row[d]+=1;
		StateChMx.push_back(row);
	}
	//----------------------------------------
	//foreach(VecofInts &v, StateChMx){
	//	cout << endl;
	//	foreach(int &a, v)
	//		cout << a << ", ";
}
//*************************************************************************************************
void RecordTimeSeries(Reactor &indiv, MatrixDouble &timeseries, double t){
	if(timeseries.size() != p.ReactantsToRecord.size()+1) MyErr("Need to resize timeseries");
	unsigned i=1;
	
	timeseries[0].push_back(t);
	foreach(int &r, p.ReactantsToRecord)
		timeseries.at(i++).push_back(indiv.Reactants[r].conc);
}
//*************************************************************************************************
void RunReactorGillespie(Reactor & indiv, MatrixDouble & timeseries, double runtime, double t_init = 0.0){
	size_t ReactantCount = indiv.Reactants.size();
	size_t ReactionCount = indiv.Reactions.size();
	double t=0, samplepoint, asum, sum, atarget, tau;
	const double asum_min = 1e-20;
	
	timeseries.clear();
	timeseries.resize(p.ReactantsToRecord.size()+1); //Need a vector for each reactant and the time
	//vector<double> time, A, B;
	unsigned m;
	double *a; //Propensity array
	
	t = t_init;
	a = (double *) malloc(ReactionCount*sizeof(double));
	MatrixInt StateChMx;
	ConstructStateChMx(indiv, StateChMx);
	
	RecordTimeSeries(indiv, timeseries, t);
	samplepoint = t + p.stepsize;
	
	while(t<p.runtime){
		asum = 0;
		//cout << t << ", ";
		
		for(size_t i=0; i < ReactionCount; i++){
			int in1 = indiv.Reactions[i].IN[0];
			int in2 = indiv.Reactions[i].IN[1];
			int ri = indiv.Reactions[i].RateIndex;
			a[i] = indiv.Rates[ri];
			
			if(in1<0 && in2<0) MyErr("in1 & in2 <=0 "); //Reactant number -1 means no reactant.
			if(in1 >= 0) a[i] *= indiv.Reactants[in1].conc;
			if(in2 >= 0) { //second reactant is present...
				if(in1 != in2) //not a dimerization reaction
					a[i] *= indiv.Reactants[in2].conc;
				else //dimerization reaction
					a[i] *= (indiv.Reactants[in2].conc - 1)/2;
			}
			asum += a[i];			
			
			if(in1 >= 0) //I rewrote this part to ensure that -1 is not inserted as reactant number.
				if(indiv.Reactants[in1].conc < 0) 
					MyErr("Conc < 0");
			if(in2 >= 0)
				if(indiv.Reactants[in2].conc < 0) 
					MyErr("Conc < 0");
			//cout << i << ", " << indiv.Reactants[in1].conc << ", " << indiv.Reactants[in2].conc << ", " << indiv.Reactions[i].Rate << ", " << a[i] << endl;
 			//cout << a[i] << ", ";
		}
			
		if(asum < asum_min) {
			ostringstream oss;
			oss << "ASUM TO SMALL!!!    asum: " << asum << endl; 				
			MyErr(oss.str().c_str());
		}
				
		do{ 
			tau= -log(RNDF())/asum;  //time to the next reaction
		}while(isinf(tau));

				
		atarget = RNDF() * asum; //pick the firing reaction
		sum = a[0];
		m=0;
		
		while(sum < atarget) sum+=a[++m]; //locate the firing reaction
		//cout << "asum=" << asum << ", tau=" << tau << ", Fire: " << m << endl;
		for(size_t i=0; i < ReactantCount; i++)  //fire reaction# m		
			indiv.Reactants[i].conc += StateChMx[m][i];		
		t+=tau;
		//cout << t << ", " << tau << endl;
		if(t >= samplepoint){ //Sample data points with p.stepsize
			if(t-samplepoint > p.stepsize) {
				//cout << t-samplepoint << endl;
				//MyErr("Stepsize too small for Gillespie.");
			}
			samplepoint += p.stepsize;
			//time.push_back(t);
			//A.push_back(indiv.Reactants[1].conc); B.push_back(indiv.Reactants[3].conc);
			RecordTimeSeries(indiv, timeseries, t);
		}		
		
	}

// 	timeseries.clear();
// 	timeseries.push_back(time);
// 	timeseries.push_back(A);
// 	timeseries.push_back(B);
}
//*************************************************************************************************

// //*************************************************************************************************
// void RunReactorOld(Reactor & indiv, MatrixDouble & timeseries, double runtime, double t_init = 0){
// 	size_t SYSSIZE = indiv.Reactants.size();
// 	const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
//  	vector<double> time, a, A, b, B;
//  	 	
// 	gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, SYSSIZE);
// 	gsl_odeiv_control * c = gsl_odeiv_control_y_new (1e-6, 0.0);
// 	gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (SYSSIZE);
// 	gsl_odeiv_system sys = {func, NULL, SYSSIZE, &indiv};
//  
// 	double t = t_init + 0.0, t1 = t_init + runtime;
// 	double h = 1e-6;
// 	double *y;
// 	y = (double *) malloc(SYSSIZE*sizeof(double));
// 	//________________________________________________________________________________________
// 	while (t < t1) {
// 		indiv.GetConc(y);
// 
// 		//printf ("%.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", t, y[0], y[1], y[2], y[3], y[4], y[5]);
// 		time.push_back(t);
// 		a.push_back(y[0]); A.push_back(y[1]); b.push_back(y[2]); B.push_back(y[3]);
// 
// 		int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, t1, &h, y);
// 		if (status != GSL_SUCCESS) { cout << "ERR"; break; }
// 	
// 		indiv.SetConc(y);
//    	}
// 	//________________________________________________________________________________________
// 	gsl_odeiv_evolve_free (e);
//    	gsl_odeiv_control_free (c);
//    	gsl_odeiv_step_free (s); 		
// 
// 	timeseries.clear();
// 	timeseries.push_back(time);
// 	timeseries.push_back(A);
// 	timeseries.push_back(B);
// }
// 
// //*************************************************************************************************
// void RunReactor(Reactor & indiv, MatrixDouble & timeseries, double runtime, double t_init = 0.0){
// 	size_t SYSSIZE = indiv.Reactants.size();
// 	const gsl_odeiv_step_type * T = gsl_odeiv_step_rk2; //gsl_odeiv_step_rk8pd;
//  	//vector<double> time, a, A, b, B;
// 	timeseries.clear();
// 	timeseries.resize(p.ReactantsToRecord.size()+1); //Need a vector for each reactant and the time
// 	
// 	gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, SYSSIZE);
// 	gsl_odeiv_control * c = gsl_odeiv_control_y_new (1e-6, 0.0);
// 	gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (SYSSIZE);
// 	gsl_odeiv_system sys = {func, NULL, SYSSIZE, &indiv};
//  
// 	double t = t_init + 0.0; //t1 = t_init + runtime;
// 	double h = 1e-6;
// 	double *y;
// 	y = (double *) malloc(SYSSIZE*sizeof(double));
// 	//Record initial data points
// 	indiv.GetConc(y);
// 	RecordTimeSeries(indiv, timeseries, t);
// 	//time.push_back(t);a.push_back(y[0]); A.push_back(y[1]); b.push_back(y[2]); B.push_back(y[3]);
// 	//________________________________________________________________________________________
// 	unsigned totalsteps = (unsigned) (p.runtime/p.stepsize);
// 	for(unsigned i=1; i<=totalsteps; i++) {
// 		double ti = t_init + i * p.stepsize;
// 		indiv.GetConc(y);		
// 		while (t < ti) {
// 			int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, ti, &h, y);
// 			if (status != GSL_SUCCESS) { cout << "ERR"; break; }		
// 		}
// 		//time.push_back(t);a.push_back(y[0]); A.push_back(y[1]); b.push_back(y[2]); B.push_back(y[3]);
//    		indiv.SetConc(y);
// 		RecordTimeSeries(indiv, timeseries, t);
// 	}
// 	//________________________________________________________________________________________
// 	gsl_odeiv_evolve_free (e);
//    	gsl_odeiv_control_free (c);
//    	gsl_odeiv_step_free (s); 		
// 
// // 	timeseries.clear();
// // 	timeseries.push_back(time);
// // 	timeseries.push_back(A);
// // 	timeseries.push_back(B);
// }
//***************************************************************************************************
static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
//*************************************************************************************************
void RunReactorCVode(Reactor & indiv, MatrixDouble & timeseries, realtype runtime, realtype t_init = 0.0){
	//----------------- Initialize ---------------------------------------------------------------/
	size_t SYSSIZE = indiv.Reactants.size();
	realtype reltol, abstol, t, tout;
	N_Vector y;
	void *cvode_mem;
	int flag, iout;
// 	double TMULT = 2;   	 /* output time factor     */
// 	double NOUT = 10;        /* number of output times */

	timeseries.clear();
	timeseries.resize(p.ReactantsToRecord.size()+1); //Need a vector for each reactant and the time

	y =  NULL;
	cvode_mem = NULL;
	/* Set the tolerances */
	reltol = 1.0e-8;
	abstol = 1.0e-8;
	tout = (realtype) p.stepsize; iout = 0;

	/* Create serial vector of length NEQ for I.C. and abstol */
	y = N_VNew_Serial(SYSSIZE);
	if (check_flag((void *)y, "N_VNew_Serial", 0)) MyErr("out");
	indiv.GetConc(NV_DATA_S(y));

	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) MyErr("out");

	flag = CVodeInit(cvode_mem, f, t_init, y);
	if (check_flag(&flag, "CVodeInit", 1)) MyErr("out");

	flag = CVodeSStolerances(cvode_mem, reltol, abstol);
	if (check_flag(&flag, "CVodeSStolerances", 1)) MyErr("out");

	flag = CVodeSetUserData(cvode_mem, &indiv);
	if (check_flag((void *)cvode_mem, "CVodeSetUserData", 0)) MyErr("out");
	
	flag = CVodeSetMaxNumSteps(cvode_mem, 2000);
	if (check_flag((void *)cvode_mem, "CVodeSetMaxNumSteps", 0)) MyErr("out");
	
	/* Call CVDense to specify the CVDENSE dense linear solver */
	flag = CVDense(cvode_mem, SYSSIZE);
	if (check_flag(&flag, "CVDense", 1)) MyErr("out");
	//----------------- End Initialize --------------------------------------------------------------/
	
	RecordTimeSeries(indiv, timeseries, t_init);
	while(1) {    
		indiv.GetConc(NV_DATA_S(y));		
		flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
		//cout << t << ", " <<  NV_Ith_S(y,0) <<  ", " <<  NV_Ith_S(y,1) << ", " <<  NV_Ith_S(y,2) << endl;
		indiv.SetConc(NV_DATA_S(y));
		RecordTimeSeries(indiv, timeseries, t);

		if (check_flag(&flag, "CVode", 1)) break;
		if (flag == CV_SUCCESS) {
			//iout++;
			tout += p.stepsize;
		}
		if (tout > p.runtime) {
			//RecordTimeSeries(indiv, timeseries, t);
			break;
		}
			
	}
			
/* Free y and abstol vectors */
	N_VDestroy_Serial(y);
	/* Free integrator memory */
	CVodeFree(&cvode_mem);
	if (check_flag(&flag, "CVode", 1)) {
		cout << "\n ******* CRASH ******** " << endl;
		cout << indiv.PrintReactions();
		MyErr("Cvode error.");
	}
}
#define CLBACTBARRIERLOW 100
#define CLBACTBARRIERHIGH 250
#define TIME 0
#define CLN 1
#define CLB 2
#define SIC1 3
//****************************************************************************************************
/*
void ScoreForTimelyandSharpTransition(Reactor & indiv, Reactor initcond, MatrixDouble & timeseries, int step=0){
	// I use fitnessparam to specify desired timing, and fitnessparam2 to specify the fitness function.
	// fitnessparam2=0 means use timing+sharpness, 1 means use only timing, 2 use only shaprness.
	int activatedlow = false, activatedhigh = false, clnactivated=false;
	unsigned int lowindex=-1, highindex=-1;
	const double DESIREDCLNLEVEL = p.FitnessParam;	
	double DESIREDTIME=0;
	double timing=0, sharpness=0, ActivationTime=0;
	
	timeseries.clear();	
	//cout << "Fitnesses: " ;
	double fit = 0;
	
	for(unsigned j=0; j < p.FITNESSAVGOVER; j++){
		activatedlow = false, activatedhigh = false;
		lowindex=-1, highindex=-1;
		double f = 0;
		unsigned int i = 0;
		
		//indiv.SetInitCond(initcond); ///IC is now set before the reaction set passed to this function.
		RunReactorCVode(indiv, timeseries, (realtype) p.runtime);
		//cout << "End run." <<endl;				
		while(i < timeseries[TIME].size()) {// [clb]>100??
			if(timeseries[CLB][i] >= CLBACTBARRIERLOW) {
				activatedlow = true; 
				lowindex=i; 
				break; 				
			}
			i++;
		}
		//using last value of i!
		while(i < timeseries[TIME].size()) {// clb>>500??
			if(timeseries[CLB][i] >= CLBACTBARRIERHIGH) { 
				activatedhigh = true; 
				highindex=i; 
				break; 				
			}
			i++;
		}
		
		for(i=0; i < timeseries[TIME].size(); i++) {// cln>TRIGGER POINT?
			if(timeseries[CLN][i] >= DESIREDCLNLEVEL) { // Translate trigger point to time
				DESIREDTIME=timeseries[TIME][i]; 
				clnactivated=true; //Cln crossed the trigger point
				break; 				
			}
		}
		
		//.cout << "desired time: " << DESIREDTIME << " "; 	
		if (activatedlow && clnactivated) {
			ActivationTime = timeseries[TIME][lowindex];
			timing = pow(ActivationTime - DESIREDTIME, 2);	
			
			if(p.FitnessParam2 < 2 ) f = timing; 
			
			if (activatedhigh) {
				double ActivationHighTime = timeseries[TIME][highindex];
				sharpness = abs(ActivationHighTime - ActivationTime);
				if(p.FitnessParam2==0 || p.FitnessParam2==2) f += sharpness;
				//.cout << "act time " << ActivationTime << ", " << ", high: " << ActivationHighTime << endl;
			}
			else f = 5000.;
		}
		else f = 10000.;
		//cout << f << ", ";
		fit += f;
	}
	
	indiv.FitnessScore = fit / p.FITNESSAVGOVER;
	
	//cout << "activated:" << activatedlow << " at " << lowindex << " , t=" << ActivationTime \
	//	<<" score is " << indiv.FitnessScore << ", timing is" << timing << ", sharpness is" << sharpness << endl;
	 
	return;
}
 */
//*****************************************************************************************************
void SmoothTimeSeries(VecofDs & series, const int avrange){
	int size = series.size();
	if(avrange > size) MyErr("SmoothTimeSeries: avrange>size");
	
	for(int i=0; i<size-avrange; i++){
		double avg=0;
		for(int j=0; j<avrange; j++){
			avg+=series.at(i+j);
		}
		avg/=avrange;
		//cout << series.at(i) << ", " << avg << endl; 
		series.at(i)=avg;		
	} 
}
//*****************************************************************************************************
double CalcTimingScore(const MatrixDouble & timeseries, const VecofDs & ALLSIC1, const VecofDs & ALLCLB){
	/// Timing score is (T_act - T_desired)^2
	/// Returns -1 if not activated.
	unsigned int rows, cols;
	rows = timeseries[TIME].size();
	cols = timeseries.size();	
	/// Locate start
	double diffmax=0, START=-1, score=-1;	
	
	if(p.startiszero) 
	    START=0;
	else {
	    for(unsigned int i=1; i < rows; i++){
		    double diff = timeseries[CLN][i] - timeseries[CLN][i-1];
		    if(diff > diffmax){
			    diffmax = diff;
			    START = timeseries[TIME][i];
		    }	
	    }
	}
	  
	/// Locate activation
	bool activated=false;
	double ActivationTime=-1;
	for(unsigned int i=0; i<rows; i++) {
		if(ALLSIC1[i] < ALLCLB[i]) {
			activated = true;
			ActivationTime=timeseries[TIME][i];
			break;
		}		
	}
	
	if(activated){		
		double timing = ActivationTime - (START + p.FitnessParam);
		score = pow(timing,2);
	}
	
	return score;
}

//*****************************************************************************************************
double CalcSharpnessScoreUsingfClb(const MatrixDouble timeseries, const VecofDs & ALLSIC1, const VecofDs & ALLCLB){
	double score=0;
	double high, low, mid, deltay, deltat;
	int i;
	const int last = timeseries[TIME].size()-1;
	
	low =timeseries[CLB][0];
	high=timeseries[CLB][last];
	mid =low + (high-low)/2.0;	
	
	for(i=0; i<=last; i++){
		/// locate the midpoint;
		if(timeseries[CLB][i] > mid) break;		
	}
	if(i==last) 
		return score;
	
	deltay= timeseries[CLB][i] - timeseries[CLB][i+1];
	deltat= timeseries[TIME][i] -timeseries[TIME][i+1];
	//cout << high << ", " << mid << ", " << timeseries[TIME][i] << ", " << deltay << ", " << deltat << endl;
	score = -deltay/deltat;
	if(ALLCLB[last] < ALLSIC1[last]) /// Use Sic1 level to score if Clb is not activated.
	  score=ALLSIC1[last]/ALLSIC1[0]; /// The lower Sic1 the fitter.
	
	return score;
}
//*****************************************************************************************************
int LocateMinMaxSlope(const VecofDs timeseries, bool FindMin, double &deltay) {
	int i, rows=timeseries.size();
	double diffmax=0, diffmin=0;	
	int location = 0;
	
	if(!FindMin) { ///Locate Maximum
		for(i=1; i < rows; i++){
			double diff = timeseries[i] - timeseries[i-1];
			if(diff > diffmax){
				diffmax = diff;
				location = i;
			}	
		}	
	} 
	else { ///Locate minumum
		for(i=1; i < rows; i++){
			double diff = timeseries[i] - timeseries[i-1];
			if(diff < diffmin){
				diffmin = diff;
				location = i;
			}	
		}	
	}
	deltay = diffmin;
	return location;
}
//*******************************************************************************************************
void CalcScoresUsingSic1(const MatrixDouble & timeseries, const VecofDs & ALLSIC1, double &tscore, double &sscore){
	const bool MIN=true, MAX=false;
	int i, rows, cols;
	double deltay=0, deltat=0, START, ACTIVATIONTIME;
	double DESIREDTIME = p.FitnessParam;
		
	rows = timeseries[TIME].size();
	cols = timeseries.size();	
	
	/// Locate start
	if(p.startiszero) 
	    START=0;
	else {
	    i = LocateMinMaxSlope(timeseries[CLN], MAX, deltay); 
	    START = timeseries[TIME][i];
	}
	
	/// Locate activation by fastest Sic1 degredation	
	i = LocateMinMaxSlope(ALLSIC1, MIN, deltay); 
	ACTIVATIONTIME = timeseries[TIME][i];
	deltat = timeseries[TIME][i] - timeseries[TIME][i-1];
	
	tscore = (ACTIVATIONTIME-START) - DESIREDTIME;
	tscore = pow(tscore, 2);
	sscore = deltay/deltat;
	//cout << "START:" << START << ", ACTIVATIONTIME:" << ACTIVATIONTIME << ", tscore" << tscore << ", sscore:" << sscore << endl;
}
//*******************************************************************************************************
double CalcSharpnessScoreUsingSic1(const MatrixDouble timeseries, const VecofDs & ALLSIC1, const VecofDs & ALLCLB){
	double score=0;
	double high, low, mid, deltay, deltat;
	int i;
	const int last = timeseries[TIME].size()-1;
	
	high=ALLSIC1[0];
	low =ALLSIC1[last];
	mid =low + (high-low)/2.0;
	
	for(i=0; i<=last; i++){
		/// locate the midpoint;
		if(ALLSIC1[i] < mid) break;		
	}
	if(i==last) 
		return score;
	
	deltay= ALLSIC1[i] - ALLSIC1[i+1];
	deltat= timeseries[TIME][i] -timeseries[TIME][i+1];
	//cout << mid << ", " << deltay << ", " << deltat << endl;
	score = deltay/deltat;
	return score;
}
//*****************************************************************************************************
double CalcSharpnessScore(const MatrixDouble timeseries, const VecofDs & ALLSIC1, const VecofDs & ALLCLB){
	
	bool activatedhigh = false;
	unsigned int lowindex=-1, highindex=-1, i;
	double sharpness=-1;
	//const double CLBMAX=1500;
	const int last = timeseries[TIME].size()-1;
	
	///include sic1 level in the score 

	i=0;
	while(i < timeseries[TIME].size()) {// [clb]>100??
		if(timeseries[CLB][i] >= CLBACTBARRIERLOW) {
			lowindex=i; 
			break; 				
		}
		i++;
	}
	
	//using last value of i!
	while(i < timeseries[TIME].size()) {// clb>>500??
		if(timeseries[CLB][i] >= CLBACTBARRIERHIGH) { 
			activatedhigh = true; 
			highindex=i; 
			break; 				
		}
		i++;
	}
	
	if(activatedhigh){
		double ClbLowTime = timeseries[TIME][lowindex];
		double ClbHighTime = timeseries[TIME][highindex];
		sharpness = abs(ClbHighTime - ClbLowTime);
	}
	/*Can't measure shaprness. Calculate score based on max fClb level*/
	if(sharpness==-1) 
		sharpness =  ALLSIC1[last];
	//else
	//	sharpness +=  ALLSIC1[last]/20; 
	
	return sharpness;

}
//****************************************************************************************************
void ScoreTransition(Reactor & indiv, Reactor initcond, MatrixDouble & timeseries, int step=0){
	// I use fitnessparam to specify desired timing, and fitnessparam2 to specify the fitness function.
	// fitnessparam2=0 means use timing+sharpness, 1 means use only timing, 2 use only shaprness.

	VecofDs ALLSIC1, ALLCLB;
	unsigned int rows, cols;
	double tscore, sscore;
	
	//indiv.PrintReactants();
	//cout << "-------------------------------" << endl;
	
	timeseries.clear();	
	if(p.stochastic) 
		RunReactorGillespie(indiv, timeseries, (realtype) p.runtime);	
	else 
		RunReactorCVode(indiv, timeseries, (realtype) p.runtime);
	
	if(p.noscoring) {
	  indiv.FitnessScore=0;
	  return;
	}
	
	rows = timeseries[TIME].size();
	cols = timeseries.size();	
	//cout  <<"rows: " << rows << ", cols: " << cols << endl;
	//ALLSIC1 = (double *) malloc(rows*sizeof(double));
	//ALLCLB  = (double *) malloc(rows*sizeof(double));
	ALLSIC1.resize(rows);
	ALLCLB.resize(rows);
	
	/// Calculate total sic1
	for(unsigned int i=0; i<rows; i++){
		ALLSIC1[i]=0;
		ALLCLB[i]=0;
/*		for(unsigned int j=SIC1; j<cols; j++){
			cout << j << ", ";
			//ALLSIC1[i]+=timeseries[j][i];
		}*/
		foreach(int &j, p.SIC1cols){
			ALLSIC1[i]+=timeseries[j][i];
		}
		foreach(int &j, p.CSPcols){
			ALLSIC1[i]+=timeseries[j][i];
			ALLCLB[i]+=timeseries[j][i];
		}
		foreach(int &j, p.fCLBcol){
			ALLCLB[i]+=timeseries[j][i];
		}
		
		//cout << timeseries[TIME][i] << ", " << ALLSIC1[i] << endl;		
	}

	/// Smooth the time series before calculating scores.
	if(p.stochastic) { 
		int range = p.smoothingrange;
		SmoothTimeSeries(ALLSIC1,range);
		SmoothTimeSeries(ALLCLB,range);
		SmoothTimeSeries(timeseries[CLN],range);
		SmoothTimeSeries(timeseries[CLB],range);
		SmoothTimeSeries(timeseries[SIC1],range);
		int newsize = rows-range;
		ALLCLB.resize(newsize); //Last range# of points can't be averaged. Remove them.
		ALLSIC1.resize(newsize);
		foreach(MatrixDouble::value_type &v, timeseries){
			v.resize(newsize);
		}
	}
	if(p.usesic1forfitness) { 
	    /// Use Sic1 sharpest degradation point for timing and sharpness measurements.
	    CalcScoresUsingSic1(timeseries, ALLSIC1, tscore, sscore);
	}
	else {
	    /// Use Clb activation point (Clb>Sic1) for timing. 
	    sscore=CalcSharpnessScoreUsingfClb(timeseries,ALLSIC1, ALLCLB);
	    tscore=CalcTimingScore(timeseries, ALLSIC1, ALLCLB);
	}
	
	 ///only timing
	if(p.FitnessParam2==1) {
		if(tscore==-1) tscore=1000000;
		indiv.FitnessScore=tscore;
		return;
	}
	
	///only sharpness		
	if(p.FitnessParam2==2) { 
		if(sscore==-1) {
			sscore=1000000;
			if(tscore>0) sscore=500000; //it's not as bad if activated.
		}
		indiv.FitnessScore=sscore;
		return;
	}
	
	///timing and sharpness
	if(tscore==-1) // couldn't measure timing. this is the worst.
		tscore = 500000;
	if (sscore==-1) //measured timing but not sharpness. this is not as bad.  
		sscore = 500000;
	
	indiv.FitnessScore = (p.SharpnessWeight*sscore) + (p.TimingWeight*tscore);
		
	//cout << "Sscore:" << sscore << ", Tscore:" << tscore << endl;
 	//cout << "act " << ActivationTime << ", START " << START << ", FitParam " << p.FitnessParam << ", timing" 
 	//	<< ActivationTime - (START + p.FitnessParam) << ", Fitness "<< indiv.FitnessScore << endl;
}

//****************************************************************************************************
void ScoreForTimelyTransitionAtClnThreshold(Reactor & indiv, Reactor initcond, MatrixDouble & timeseries, int step=0){

	double *ALLSIC1;
	unsigned int rows, cols;
	
	timeseries.clear();	
	RunReactorCVode(indiv, timeseries, (realtype) p.runtime);	
	
	rows = timeseries[TIME].size();
	cols = timeseries.size();	
	//cout  <<"rows: " << rows << ", cols: " << cols << endl;
	ALLSIC1 = (double *) malloc(rows*sizeof(double));
		
	/// Calculate total sic1
	for(unsigned int i=0; i<rows; i++){
		ALLSIC1[i]=0;
		for(unsigned int j=SIC1; j<cols; j++){
			ALLSIC1[i]+=timeseries[j][i];
		}
		//cout << timeseries[TIME][i] << ", " << ALLSIC1[i] << endl;
	}
	
	/// Locate activation
	bool activated=false;
	double ActivationTime=-1, ClnActLevel=-1;
	for(unsigned int i=0; i<rows; i++) {
		if(ALLSIC1[i] < timeseries[CLB][i]) {
			activated = true;
			ActivationTime=timeseries[TIME][i];
			ClnActLevel=timeseries[CLN][i];
			break;
		}		
	}
	
	if(activated){		
		double distance = ClnActLevel - p.FitnessParam;
		indiv.FitnessScore = pow(distance,2);
	}
	else indiv.FitnessScore = 1000000;

 	//cout << "act " << ActivationTime << ", FitParam " << p.FitnessParam << ", ClnActLevel " 
	//<< ClnActLevel << ", " << p.FitnessParam << ", Fitness "<< indiv.FitnessScore << endl;
}

