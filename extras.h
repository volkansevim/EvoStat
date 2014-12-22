//--------------------------------------------------------------------------
inline double GRand(double Sigma) {
	return Sigma * sqrt(-2*log(RGauss())) * cos(2*M_PI*RGauss());
}
//--------------------------------------------------------------------------
double Gaussian(double Sigma){ 
	float x1, x2, w, y1, y2;

	if(Sigma < 10e-10) return 0.;
	do {
		x1 = (2.0 * RNDF()) - 1.0;
		x2 = (2.0 * RNDF()) - 1.0;
		w = (x1 * x1) + (x2 * x2);		
	} while ( w >= 1.0 );

	w = sqrt((-2.0 * log(w))/w);
	y1 = x1 * w;
	y2 = x2 * w;
	
	return Sigma * y1;
}
//--------------------------------------------------------------------------
double ExtNoiseGaussian(double sigma){
	double gr = Gaussian(sigma);
	const double CUTOFF = 0.4; 
	while(abs(gr) > CUTOFF ) {
		gr = Gaussian(sigma); 	//Assert %noise < %CUTOFF		
	}
	return 1. + gr; 	// %extrinsic noise level 
	
}

//--------------------------------------------------------------------------
double ExtNoise(double sigma){
	double gr = (RNDF(2*sigma) - sigma);

	return 1. + gr; 	// %extrinsic noise level 	
}

//--------------------------------------------------------------------------
class Parameters {
public:
	psimag::KISS::seed_type SEED;
	unsigned GENERATIONS;
	unsigned POPSIZE;
	unsigned FITNESSAVGOVER;
	double MuReacRate, MuNewTFDNA, MuNewDimer;	
	double runtime; 
	double FitnessParam;	// I implemented this for timing. No specific purpose.
	int    FitnessParam2; 	// Another parameter to specify the fitness function. 
	double SharpnessWeight, TimingWeight;
				// Refer to the definition in scoring.h for usage.
	double MaxRate, MinRate, MaxRndRate, MinRndRate;
	
	double stepsize;
	double lognoisemax;		// lognoise is produced as 10^r. r is a rnd number on [-max, max]
	double sigma; 			// Stdev of the random Gaussian distr for extrinsic noise. 
					      // sigma>0 implies cell to cell variability
	bool stochastic; 		// ODE or Gillespie? (Intrinsic noise)
	bool extrinsicnoise;    	// Vary rates of transcription? (works with both ODE and Gillespie)
	bool detailedlogs;		// Write out mutatables after every time step?
	bool usesic1forfitness;		// Use Sic1 sharpest degr. for timing and sharpness. False uses Clb.
	bool startiszero;		// Take t=0 as START.
	int randomize;			// randomize: 0=none, 1=only mutatable reaction rates, 2=all rates
	unsigned smoothingrange;// #of elements to use in moving average.
	bool noscoring;
	
	string FileNameTemplate;
	vector <int> ReactantsToRecord, CSPcols, SIC1cols, CLNcol, fCLBcol;
	vector <int> RatesToMutate;
	vector <string> ReactantsToRecordNames;
	string FitnessFunc, parameterstrlong, parameterstrshort;	
	
	map< string, int > ReactantTypes;
	double MaxSynthRate, MaxDNATFAssocDissocRate;
	char reactionfile[300], reactantfile[300];
	char exec[200];
	
	Parameters(){
		SEED=12345678;
		GENERATIONS = 1;
		POPSIZE = 2;
		FITNESSAVGOVER = 1;
		FitnessParam = 35;  // Desired timing in Clb activation
		FitnessParam2= 0;   // Use both timing and sharpness
		MuReacRate = 0.3; 	// Mutation rate of reaction rates
		MuNewTFDNA = 0.3; MuNewDimer = 0.3;
		MaxSynthRate = 0.1;
		runtime = 50;  	//Integration duration of ODEs
		stepsize = .25; 		//Report concentrations in large timesteps
		MaxDNATFAssocDissocRate = 1;
		strcpy(reactantfile, "reactants.dat");
		strcpy(reactionfile, "reactions.dat");
		FileNameTemplate = "output";
		usesic1forfitness = false;
		stochastic = false;
		extrinsicnoise = false;
		detailedlogs = false;
		randomize = 0;
		sigma = 0.0;
		lognoisemax = 1.0;
		TimingWeight=1.;
		SharpnessWeight=1.;
		startiszero=false;
		MaxRate=MaxRndRate=MAXRATE;
		MinRate=MinRndRate=MINRATE;
		
		smoothingrange=10;
		noscoring=false;
		
		ReactantTypes["GENE"] = GENE; 
		ReactantTypes["PROTEIN"] = PROTEIN;
		ReactantTypes["PCOMPLEX"] = PCOMPLEX; 
		ReactantTypes["PGCOMPLEX"] = PGCOMPLEX;		
		ReactantTypes["NOISYPROTEIN"] = NOISYPROTEIN;		
	}
	
	double RandomSynthRate() {
		return MaxSynthRate*RNDF();
	}
	
	double RandomDNATFAssocDissocRate() {
		return MaxDNATFAssocDissocRate*RNDF();
	}
	
	double RandomReactionRate() {		
		double range, exp, result;
		//log sampling
		range = log10(MaxRndRate) - log10(MinRndRate);
		exp = RNDF()*range;	//Sample between 10^-6 and 1
		result = MaxRndRate*pow(10, -exp);
		//cout << "range: " <<  range << ", exp:" << exp << ", result" << result << endl;
		if(result>MaxRate) result=MaxRate;
		if(result<MinRate) result=MinRate;
		return result;
	}

	void Assign(string variable, string value) {
		boost::to_lower(variable);
		if(variable == "generations")
			GENERATIONS = boost::lexical_cast<unsigned>(value);		
		else if(variable == "fitnessavgover")
			FITNESSAVGOVER = boost::lexical_cast<unsigned>(value);		
		else if(variable == "popsize")
			POPSIZE = boost::lexical_cast<unsigned>(value);		
		else if(variable == "seed")
			SEED =  (psimag::KISS::seed_type) boost::lexical_cast<int>(value);		
		else if(variable == "runtime")
			runtime = boost::lexical_cast<double>(value);		
		else if(variable == "stochastic")
			stochastic = boost::lexical_cast<bool>(value);		
		else if(variable == "extrinsicnoise")
			extrinsicnoise = boost::lexical_cast<bool>(value);				
		else if(variable == "detailedlogs")
			 detailedlogs = boost::lexical_cast<bool>(value);									
		else if(variable == "sigma")
			sigma = boost::lexical_cast<double>(value);		
		else if(variable == "randomize")
			randomize = boost::lexical_cast<unsigned>(value);
		else if(variable == "fitnessparam")
			FitnessParam = boost::lexical_cast<double>(value);	
		else if(variable == "fitnessparam2")
			FitnessParam2 = boost::lexical_cast<int>(value);	
		else if(variable == "mureacrate"){
			MuReacRate = boost::lexical_cast<double>(value);		
			assert(MuReacRate >= 0 && MuReacRate <= 1);
		}
		else if(variable == "lognoisemax"){
			lognoisemax = boost::lexical_cast<double>(value);		
			assert(lognoisemax >= 0);	
		}
		else if(variable == "sharpnessweight"){
			SharpnessWeight = boost::lexical_cast<double>(value);		
			assert(SharpnessWeight >= 0.);
		}
		else if(variable == "timingweight"){
			TimingWeight = boost::lexical_cast<double>(value);		
			assert(TimingWeight >= 0.);
		}
		else if(variable == "usesic1forfitness"){
			usesic1forfitness = boost::lexical_cast<bool>(value);					
		}		
		else if(variable == "startiszero"){
			startiszero = boost::lexical_cast<bool>(value);					
		}		
		else if(variable == "maxrate")
			MaxRate = boost::lexical_cast<double>(value);	
		else if(variable == "minrate")
			MinRate = boost::lexical_cast<double>(value);	
		else if(variable == "maxrndrate")
			MaxRndRate = boost::lexical_cast<double>(value);	
		else if(variable == "minrndrate")
			MinRndRate = boost::lexical_cast<double>(value);	
		else if(variable == "noscoring")
			noscoring = boost::lexical_cast<bool>(value);		
		else if(variable == "smoothingrange")
			smoothingrange = boost::lexical_cast<unsigned>(value);		

		else MyErr("Argument not recognized");

	}
};
Parameters p;


//********************************************************************************************
void WriteParameters() {
	char filename[150];
	
	sprintf(filename, "%s-parameters.prm", p.FileNameTemplate.c_str());
	ofstream file (filename);
	ostringstream oss;
	time_t rawtime;
	time ( &rawtime );
		
	oss << "# Date & Time = " <<  ctime (&rawtime);
	oss << "# Executable = " << p.exec << endl;
	oss << "# Reactant File = " << p.reactantfile << endl;
	oss << "# Reaction File = " << p.reactionfile << endl;
	oss << "# SEED = " << p.SEED << endl;
	oss << "# Runtime = " << p.runtime << endl;
	oss << "# GENERATIONS = " << p.GENERATIONS << endl;
	oss << "# POPSIZE = " << p.POPSIZE << endl;
	oss << "# stepsize = " << p.stepsize << endl;
	oss << "# stochastic? " << p.stochastic << endl;
	oss << "# randomize = " << p.randomize << endl;
	oss << "# extrinsic noise? " << p.extrinsicnoise << endl;
	oss << "# Sigma (for extr. noise) = " << p.sigma << endl;
	oss << "# lognoisemax (for reaction class 2) = " << p.lognoisemax << endl;
	oss << "# \n";
	oss << "# Fitness Func. = " << p.FitnessFunc << endl;
	oss << "# Fitness Param = " << p.FitnessParam << endl;
	oss << "# Fitness Param2 = " << p.FitnessParam2 << endl;
	oss << "# Use Sic1 Sharpest Degr. for Fitness = " << p.usesic1forfitness << endl;
	oss << "# SharpnessWeight = " << p.SharpnessWeight << endl;
	oss << "# TimingWeight = " << p.TimingWeight << endl;
	oss << "# MaxRate = " << p.MaxRate << endl;
	oss << "# MinRate = " << p.MinRate << endl;
	oss << "# smoothingrange = " << p.smoothingrange << endl;
	oss << "# MuReacRate = " << p.MuReacRate <<endl;
	oss << "# MuNewTFDNA = " << p.MuNewTFDNA << endl; 
	oss << "# MuNewDimer = " << p.MuNewDimer << endl;
	oss << "# MaxSynthRate = " << p.MaxSynthRate << endl;
	oss << "# MaxDNATFAssocDissocRate = " << p.MaxDNATFAssocDissocRate << endl;	
	p.parameterstrlong = oss.str();
	//cout << oss.str();
	file << oss.str();
	file.close();
}
//********************************************************************************************
string ShortParameterString() {
	ostringstream oss;
		
	//Construct the short version here:	
	oss << "Fitness Func. = " << p.FitnessFunc << endl;
	oss << "GENERATIONS = " << p.GENERATIONS << endl;
	oss << "POPSIZE = " << p.POPSIZE << endl;
	oss << "MuReacRate = " << p.MuReacRate <<endl;
	oss << "MuNewTFDNA = " << p.MuNewTFDNA << endl; 
	oss << "MuNewDimer = " << p.MuNewDimer << endl;
	oss << "Sigma = " << p.sigma << endl;
	oss << "Extrinsic noise? " << p.extrinsicnoise << endl;
	
	p.parameterstrshort = oss.str();	
	return p.parameterstrshort;
}
//**************************************************************************************************
void CmdLineParser(int argc, char *argv[]) {
	char msg[]=
	"Need at least 2 arguments: reactantfile reactionfile VARIABLE=VALUE\n \
	Allowed variables are GENERATIONS, POPSIZE, SEED, sigma, extrinsicnoise, mureacrate, lognoisemax, randomize";
	boost::regex pat("([\\w]+)[=:]([\\w\\d\\.]+)"); 
	
	if(argc < 3) MyErr(msg);
	strcpy(p.reactantfile, argv[1]);
	strcpy(p.reactionfile, argv[2]);
	
	boost::regex pat2("(.+?)(\\.[^.]*$|$)"); //match filename and extension
	boost::smatch matches;
	string f(p.reactionfile);
	if (boost::regex_match(f, matches, pat2)) { //remove the extension in the filename
		p.FileNameTemplate = matches[1];		
	}
	
	strcpy(p.exec, argv[0]);
	for(int i = 3; i < argc; i++){
		boost::smatch matches;
		string argument = argv[i];
		if(boost::regex_match(argument, matches, pat)) {
			p.Assign(matches[1], matches[2]); //assign variable the value
		}
	}	

	WriteParameters();
}

//**************************************************************************************************
void GetReactantNamesForSummation() {
	int count=1; //Start from 1 because the column 0 is time.
	boost::smatch matches;
	boost::regex patCSP("CSP.+");
	boost::regex patSIC(".+SIC1.+");
	boost::regex patCLN("CLN");
	boost::regex patfCLB("fCLB");
		
	p.CSPcols.clear();
	foreach(string &n, p.ReactantsToRecordNames) {
		if(boost::regex_match(n, matches, patCSP)) 
			p.CSPcols.push_back(count);

		if(boost::regex_match(n, matches, patSIC)) 
			p.SIC1cols.push_back(count);

		if(boost::regex_match(n, matches, patCLN)) 
			p.CLNcol.push_back(count);
	
		if(boost::regex_match(n, matches, patfCLB)) 
			p.fCLBcol.push_back(count);
		
		count++;
	}
	
}
