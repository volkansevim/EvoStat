# define SQR(x)((x)*(x)) 
#include "mpi.h"
#include <iostream>
#include <algorithm>
#include "my_gnuplot_i.hpp" //Gnuplot class handles POSIX-Pipe-communication with Gnuplot
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iterator>
#include <map>
#include <deque>
#include <set>
#include <fstream>
#include <iomanip>
#include <string>
#include <assert.h>
#ifdef __LP64__
  #define  KISS_INT 1
#endif 
#include <boost/regex.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
//#include <omp.h>
#include "kiss64bit.h"
	

#include <gsl/gsl_errno.h>
// #include <gsl/gsl_matrix.h>
// #include <gsl/gsl_odeiv.h>

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */

#define foreach BOOST_FOREACH
#define MAXRATE 1
#define MINRATE 1E-10
#define BLANK -1 
#define GENE 0
#define PROTEIN 1
#define PCOMPLEX 2
#define PGCOMPLEX 3
#define NOISYPROTEIN 4
#define NOTNOISY 0
#define NOISYLIN 1
#define NOISYLOG 2
#define NOISYLINCORR 3
#define ICNOISYLIN 4
#define ICNOISYLOG 5
#define ICNOISYLINCORR 6
#define HILLKINETICS 7
#define MASSACTION 0
#define MASTER 0
#define WORKREQTAG 0
#define DEFTAG 1
#define RESULTTAG 2
#define TERMINATETAG 3
#define DETAILSSIZE 4
#define ORDER 0
#define REALIZATIONNO 1
#define GEN 2
#define SLAVESEED 3
#define WRITEOUT 1
#define MUTATABLE 99
// noisylin and log are reaction classes. lin and log specifiy wheter noise is generated 
// on a linear or log scale.
#define __CLUSTER__ 1
#define __PARALLEL__ 1

using namespace std;
using psimag::KISS;
KISS RNDF, RGauss, RSlaveSeeder;

#include "err.h"
#include "extras.h"
#include "reactor.h"
#include "scoring.h"
#include "importer.h"

//*************************************************************************************************
void WriteTimeSeries(MatrixDouble & timeseries, char *name, string extrastuff = " "){
	unsigned rows = timeseries[0].size();
	unsigned cols = timeseries.size();
	char buf[100];

	ofstream file (name);
	file.precision(5);
	file << p.parameterstrlong; //Write the parameters first.
	file << "#-------------------------------------\n";
	file << extrastuff << endl;
	
	if (!file.is_open()) MyErr("Can't open timeseries file.");
	file << "#t"; 
	foreach (string &n, p.ReactantsToRecordNames)
		file << "\t" << n;
	file << endl;
		
	for(unsigned i=0; i < rows; i++) {
		for(unsigned j=0; j < cols; j++) {
			sprintf(buf, "%10.4f", timeseries[j][i]);
			file <<  setprecision(8) << timeseries[j][i]; 
			//file << buf;
			if (j < cols-1) 
				file << "\t";
		}
		file << endl;
	}
			
	file.close();
}
//*************************************************************************************************
int Init(Reactor & indiv){
	ImportReactants(p.reactantfile, indiv);
	ImportReactions(p.reactionfile, indiv);
//	ImportReactants("reactants.dat", indiv);
//	ImportReactions("reactions.dat", indiv);

	return 0;
}
//*************************************************************************************************
void InitColony(ReactorMap &col, int size, int randomize, int verbose=0) {
	Reactor dummy;

	ImportReactants(p.reactantfile, dummy, verbose);
	ImportReactions(p.reactionfile, dummy, verbose);
	GetReactantNamesForSummation();
	
/*	cout << "CSPcols: ";
	foreach(int i, p.CSPcols) cout << i << " ,"; 
	cout << endl;	

	cout << "SIC1cols: ";
	foreach(int i, p.SIC1cols) cout << i << " ,"; 
	cout << endl;	

	cout << "cln cols:" << p.CLNcol.at(0) << endl;
	cout << "fclb col:" << p.fCLBcol.at(0) << endl;		*/																			
	
	col[0] = dummy; //put the originial in slot 0
	unsigned reactioncount = dummy.Rates.size();
	
	if(verbose){
		if(randomize==0) cout << "Rates will be used as provided. No randomization.\n";
		else if(randomize==1) cout << "Randomizing mutatable rates.\n";
		else if(randomize==2) cout << "Randomizing all rates.\n";
	}
	
	for(int i=0; i<size; i++){ //Randomize the rest of the colony
		col[i] = dummy;		
		switch(randomize){
			case 1:  //This randomizes ONLY MUTATABLE rates.
				for(unsigned int j=0; j < p.RatesToMutate.size(); j++) { 
					int reactionindex = p.RatesToMutate.at(j);
					col[i].Rates.at(reactionindex) = p.RandomReactionRate();						
				}				
				break;
				
			case 2: //This randomizes all rates.
				for(unsigned j=0; j<reactioncount; j++)
					col[i].Rates[j] = p.RandomReactionRate();
				break;
		
		}
	}
	//dummy.PrintReactants();
}
//*************************************************************************************************
void PlotTimeseriesPS(MatrixDouble & timeseries, char *filename, Reactor &indiv){
    Gnuplot g1;
    g1.reset_plot();
    g1.set_grid();
    g1.set_style("lines");
	g1.savetops("dummy.ps");
	
	string label = "set label \"" + indiv.PrintReactions() + "\n" + ShortParameterString() + "\" at graph 0.05, graph 0.97 font\"Helvetica,10\"";
	boost::regex newline("\n",boost::regex_constants::icase|boost::regex_constants::perl); 
	string rep("\\\\n");
	label  = boost::regex_replace(label, newline, rep); //replace all newlines with \n

	//cout << label;
	g1.cmd(label);
	//g1.savetops(filename);
	
	for(unsigned i=1; i<timeseries.size(); i++) {
		string n = p.ReactantsToRecordNames[i-1]; //names vector don't contain time
		g1.plot_xy(timeseries[0], timeseries[i], n.c_str());
	}	
	g1.savetops(filename);
	g1.replot();
}
//*************************************************************************************************
void PlotTimeseries(MatrixDouble & timeseries){
  	Gnuplot g1;
    g1.reset_plot();
    g1.set_grid();
    g1.set_style("lines");
	for(unsigned i=1; i<timeseries.size(); i++) {
		string n = p.ReactantsToRecordNames[i-1]; //names vector don't contain time
		g1.plot_xy(timeseries[0], timeseries[i], n.c_str());
	}
	cin.get();
}
//*************************************************************************************************
void ExperimentalRun(Reactor &indiv, Reactor initcond) {
	MatrixDouble timeseries;
	char txtfile[150], psfile[150];
	
	indiv.PrintReactants();
	for(unsigned i=0; i < p.GENERATIONS; i++) {
		indiv.SetInitCond(initcond);
		indiv.PrintReactants();
		cout << indiv.PrintReactions();
		cout << "Experimental run mode.\n";
		if(p.stochastic){ 
			cout << "Running stochastic.\n";
			RunReactorGillespie(indiv, timeseries, p.runtime);
			sprintf(txtfile, "%s-%03d-Gillespie.txt", p.FileNameTemplate.c_str(), i);
			sprintf(psfile, "%s-%03d-Gillespie", p.FileNameTemplate.c_str(), i);
		}
		else {
			RunReactorCVode(indiv, timeseries, p.runtime);
			sprintf(txtfile, "%s-%03d-ODE.txt", p.FileNameTemplate.c_str(), i);
			sprintf(psfile, "%s-%03d-ODE", p.FileNameTemplate.c_str(), i);
		}
		
		if(i==0) { 
			PlotTimeseries(timeseries); 
			PlotTimeseriesPS(timeseries, psfile, indiv);	
		}
		WriteTimeSeries(timeseries, txtfile, indiv.PrintReactions(true));
		if(p.stochastic == 0 && !p.extrinsicnoise) break; // No need for multiple deterministic runs.
		//indiv.RandomizeRates();
	}
}
//*************************************************************************************************
void master(const int processno, const int numprocesses){		
	Reactor dummy, InitCond;
	//time_t t1, t2;	
	ReactorMap Colony; // map<int, Reactor>
	MatrixDouble timeseries;
	typedef multimap <double, int> DImultimap;
	ofstream file;
	int details[DETAILSSIZE];
	
	InitColony(Colony, p.POPSIZE, p.randomize, 1); ///Randomize mutatable rates
	InitCond = Colony[0];		
	cout << p.parameterstrlong;	
	if(Colony.size() % 2 > 0)  MyErr("Colony size is odd!");

	///-- OPEN OUTPUT---------------------------------------------------------------------------------------
	//Final reaction rates and fitnesses.
	char name[100];
	sprintf(name, "proc%03i.txt", processno);
	file.open(name);
	if (!file.is_open()) MyErr("Can't open output file.");

	file << "____ Beginning " << " ____\n";
	for(unsigned j=0; j<Colony.size(); j++){
		file << "\nReactor " << j << " --- Process " << processno << " ____________________\n";
		//Colony[j].PrintReactants();
		file << Colony[j].PrintReactions(false, false);
//		file << "Fitness Score: " << Colony[j].FitnessScore << "\n";
	}
	file << " _____________________ END BEGINNING ____________________\n";
	
	///-- RUN COLONY --------------------------------------------------------------------------------
	unsigned half = p.POPSIZE/2;
	unsigned ratescount, conccount;
	double fitness, rcvbuf;
	int source, tag, tagfilter, dummyint, error;
	MPI_Status status;		
	MPI_Request pending;
	ratescount = InitCond.Rates.size();
	conccount = InitCond.Reactants.size();
	int *processes=(int *) calloc(numprocesses, sizeof(int));	
	double *rates = (double *) calloc(ratescount, sizeof(double));
	double *concentrations = (double *) calloc(conccount, sizeof(double));
	
	for(unsigned gens=0; gens<p.GENERATIONS; gens++) {
		///-- DISTRIBUTE ---------------------------------------------------------------------------
		unsigned assigncount=0;
		unsigned resultcount=0;
		tagfilter = MPI_ANY_TAG;	
		
		cout << "Generation: " << gens << endl;		
		do {
			if(assigncount == p.POPSIZE)
				tagfilter=RESULTTAG;
			//cout << "Master recv wait\n" ;
			MPI_Recv(&rcvbuf, 1, MPI_DOUBLE, MPI_ANY_SOURCE, tagfilter, MPI_COMM_WORLD, &status);
			source = status.MPI_SOURCE;					
			tag = status.MPI_TAG;
			error = status.MPI_ERROR;
			//cout << "Master RECV. Source:" << source << ", Tag:" << tag << endl;
			
			if(tag==RESULTTAG) { 
				fitness = rcvbuf;
				//cout << "Master received result from " << source << endl;
				int indivno = processes[source];
				Colony[indivno].FitnessScore = fitness;
				resultcount++;
			}
			else if(tag==WORKREQTAG) {
				//Assign a new job to the slave.
				//cout << "Master assigning job to " << source << endl;
				Colony[assigncount].SetInitCond(InitCond);
				Colony[assigncount].GetRates(rates);				
				Colony[assigncount].GetConc(concentrations);				
				
				if(gens==p.GENERATIONS-1 && (assigncount<half || p.GENERATIONS==1)) {
					details[ORDER] = WRITEOUT;
				}
				else details[ORDER] = 0;
				details[REALIZATIONNO] = assigncount;
				details[GEN]=gens;
				details[SLAVESEED]=RSlaveSeeder(4294967295UL);
				MPI_Send(rates, ratescount, MPI_DOUBLE, source, DEFTAG, MPI_COMM_WORLD);
				MPI_Send(concentrations, conccount, MPI_DOUBLE, source, DEFTAG, MPI_COMM_WORLD);
				MPI_Send(details, DETAILSSIZE, MPI_INT, source, DEFTAG, MPI_COMM_WORLD);
				//cout << "Master: assignment complete at #" << source << endl;
				processes[source] = assigncount;
				assigncount++;
			}
		} while(resultcount<p.POPSIZE);				
		//.cout << "END GENERATION " << gens << endl;
		///---------------------------------------------------------------------------------------------
		double meanfitness=0.;	

		///-- WRITE OUTPUT -------------------------------------------------------------------------
		//Final reaction rates and fitnesses.
		meanfitness = 0.;	
		file << "***********\t\tGENERATION " << gens << "\t\t**************\n";
		for(unsigned j=0; j<half; j++){ 
			// Write out only last generation if detailed logging is not requested.
			if (p.detailedlogs || gens==p.GENERATIONS-1) { 
				file << "\nReactor " << j << endl;			
				file << Colony[j].PrintReactions(true);
				file << "Fitness Score: " << Colony[j].FitnessScore << "\n";
			}
			meanfitness+=Colony[j].FitnessScore;
		}
		file << "MeanFitness=" << meanfitness/half  << endl << endl;		

		///-- MUTATE AND SELECT -----------------------------------------------------------------------			
		DImultimap scores;
		ReactorMap dummycolony;
		
		scores.clear(); 
		dummycolony.clear();
		
		cout << "Fitnesses_:";
		for(unsigned j=0; j<Colony.size(); j++){
			//Fill the map with scores.
			double fit=Colony[j].FitnessScore;
			scores.insert(DImultimap::value_type(fit, j));		
			cout << fit << ", ";
		}
		cout << endl;
				
		///-- SELECT ---------------------------------------------------------------------------------			
		unsigned i=0;
		DImultimap::const_iterator it = scores.begin();	
		while(it != scores.end()) {	
			//Selection. Get top most fit half.
			int index = it->second;
			dummycolony[i]=Colony[index];
			dummycolony[i+half] = Colony[index];
			i++; it++;
			if(i>=half) break;
		}
		///-- MUTATE --------------------------------------------------------------------------------
		Colony = dummycolony;		
		for(unsigned k=0; k<Colony.size(); k++){ 
			//Mutate the upper half.
			if(k >= half && gens<p.GENERATIONS-1) {
				Colony[k].MutateRates();	
			}
		}	
		meanfitness/=Colony.size();	
		//.cout << "end loop\n";
	}
	//cout << "outside loop\n";
	for(int i=1; i < numprocesses; i++) {
		//cout << "Master terminating " << i << endl;
		MPI_Isend(&dummyint, 1, MPI_INT, i, TERMINATETAG, MPI_COMM_WORLD, &pending);
	}
	cout << "master complete\n";
	file.close();	
}

///-------------------------------------------------------------------------------------------------------
void slave(const int processno, const int numprocesses) { 
	Reactor dummy, InitCond;	
	MPI_Status stat;
	ReactorMap Colony; // map<int, Reactor>
	MatrixDouble timeseries;
	unsigned ratescount, conccount;
	double fitness, rcvbuf;
	int tag, dummyint;
	int details[DETAILSSIZE];
	MPI_Status status;		
	MPI_Request pending;

	InitColony(Colony, 1, 0); ///No randomization.
	InitCond = Colony[0];	
	ratescount = InitCond.Rates.size();
	conccount = InitCond.Reactants.size();
	double *rates = (double *) calloc(ratescount, sizeof(double));
	double *concentrations = (double *) calloc(conccount, sizeof(double));
	time_t rawtime;
	
	while(true) {
		//Request work
		//cout << "Slave #" << processno << " requesting work." << endl;
		MPI_Isend(&rcvbuf, 1, MPI_DOUBLE, MASTER, WORKREQTAG, MPI_COMM_WORLD, &pending); 			
		MPI_Probe(MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &status);				
		tag = status.MPI_TAG;
		//cout << "Slave#" << processno << " Message in buffer." << ", Tag:" << tag << endl;
		
		if(tag==TERMINATETAG) {
			MPI_Recv(&dummyint, 1, MPI_INT, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			//cout << "Slave #" << processno << " terminated." << endl;
			break;					
		}
	
		MPI_Recv(rates, ratescount, MPI_DOUBLE, MASTER, DEFTAG, MPI_COMM_WORLD, &stat);
		MPI_Recv(concentrations, conccount, MPI_DOUBLE, MASTER, DEFTAG, MPI_COMM_WORLD, &stat);
		MPI_Recv(details, DETAILSSIZE, MPI_INT, MASTER, DEFTAG, MPI_COMM_WORLD, &stat);
		
		Colony[0].PutRates(rates);
		Colony[0].SetConc(concentrations);
		RNDF.seed(details[SLAVESEED]);
		
/*		for(unsigned i=0; i<conccount; i++) cout << concentrations[i] << ", ";
		cout << endl;
		Colony[0].PrintReactants();
		cout << "--------------------------------" << endl;*/
		
		//time ( &rawtime );
		//cout << "Slave#" << processno << " received parameters:" << ctime(&rawtime) << endl;
		ScoreTransition(Colony[0], InitCond, timeseries);
		if(details[ORDER]==WRITEOUT) {
			char psfile[100];
			sprintf(psfile, "gen%03d_no%03d.txt", details[GEN], details[REALIZATIONNO]);
			ostringstream message;
			message << Colony[0].PrintReactions(false, true) << "\n" << "#Fitness: " \
					<< Colony[0].FitnessScore << endl; ;
			WriteTimeSeries(timeseries, psfile, message.str());
		}		
		
		time ( &rawtime );
			/*cout << "Slave#" << processno << " finished calculation. FITNESS: " 
						<< Colony[0].FitnessScore << ctime(&rawtime) << endl;*/
		
		fitness=Colony[0].FitnessScore;
		MPI_Send(&fitness, 1, MPI_DOUBLE, MASTER, RESULTTAG, MPI_COMM_WORLD);		
		
		//time ( &rawtime );
		//cout << "Slave#" << processno << " sent parameters." << ctime(&rawtime ) << endl;
	}		
}		
//*************************************************************************************************
int main (int argc, char *argv[]) {
	int numprocesses, processno;
	CmdLineParser(argc, argv);
	RNDF.seed(p.SEED);
	RGauss(p.SEED + 112233);
	RSlaveSeeder(p.SEED);
		
	setenv("DISPLAY", "0:0", 0);
	#if defined(__CLUSTER__)
		setenv("GNUTERM", "dumb", 1);
	#endif	
		
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &processno);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocesses);
	if(numprocesses<2) 
		MyErr("Need at least 2 processes.");

	if(processno==MASTER)
		master(processno, numprocesses);
	else
		slave(processno, numprocesses);
		
	MPI_Finalize();	
	return 0;
}
