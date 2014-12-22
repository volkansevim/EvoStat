#define HILLk 0
#define HILLK 1
#define HILLn 2
typedef map <int, double> IDMap;

typedef struct {
	int kindex, Kindex, nindex;
} HillFuncPrm;

class Reaction {
	// Rates are held in a seperate map. Each reaction has a rate index. 
	// This makes easier to track sister reactions.
	public:
		int IN[2];
		int OUT[2];
		//double Rate; 
		int ReactionClass;
		int RateIndex; 
		HillFuncPrm Hill;
		
		double CalcReacRate(double inputlevel, IDMap &Rates){
			//Calc rate for Hill kinetics: k x^n / (K^n + x^n)
			double inputn, k, n, Kn, r;
			
			k = Rates[Hill.kindex] ;
			n = Rates[Hill.nindex] ;
			Kn = pow(Rates[Hill.Kindex], n);
			inputn = pow(inputlevel, n);
			
			r = k*inputn/(Kn+inputn);
			//cout << "rate: " << r << ", " << inputlevel << ", n= " << k << ", Kn=" << Kn << ", n=" << n << ", " << endl;
			//cout << "kindex: " << Hill.kindex << " ," << Hill.Kindex << ", " << Hill.nindex<< endl;
			/*for(int i=0; i<Rates.size(); i++){
				cout << i << ", " << Rates[i]<< endl;
			}*/
			
			return r;
		}		
};

typedef struct {
	double conc;
	string name;
	int type; 
	int produces;
	bool recorded;
} Reactant;

typedef deque< Reaction > ReactionDeque;
typedef map< int, Reactant, less<int> > ReactantsMap;
typedef vector <double> VecofDs;
typedef vector <int> VecofInts;
typedef vector <VecofDs> MatrixDouble;
typedef vector <VecofInts> MatrixInt;
typedef map <double, double> DDMap;
//*******************************************
// A safety feature for the code, IntegrityCheck():
// This function should go thru the reactions to check that
// 1- a gene does not produce more than one protein in multiple reactions.
// 2- each gene produces at least one protein
// 3- all proteins degrade
// 4- all proteins are used in at least one reaction
// 5- each  reaction has at least one input reactant
//*******************************************

//-------------------------------------------------------------------------------------------------
class Reactor{
	public:
	ReactionDeque Reactions;
	ReactantsMap Reactants;
	/// All indexed by reaction rate number. 
	/// I didn't put them all in a single object for backward compatibility.
	IDMap Rates, NoiseLevels, NoiseTypes; 
	
	int ReactantCount;
	double FitnessScore;
	vector <int> Genes, Proteins, PGComplexes, PComplexes;
	//--------------------------------------------------------------------------------------------	
	void SetDefaultInitCond(){
		foreach(ReactantsMap::value_type &r, Reactants) {
			r.second.conc = (r.second.type==GENE ? 1.:0.);
		}		
	}
	//--------------------------------------------------------------------------------------------	
	void SetInitCond(Reactor IC){
		//Reactions = IC.Reactions;
		Reactants = IC.Reactants;
		int reactantcount=Reactants.size();
		
		int sz = Rates.size();
	
		if(p.extrinsicnoise) {
			for(int i=0; i<sz; i++) {				
				int noisetype = IC.NoiseTypes[i];
				//cout << "Set IC: " << i << ", " << noisetype << endl;
				if (noisetype == 0) continue;
				
				double noiselevel = IC.NoiseLevels[i];
				if(noisetype == NOISYLIN) 
					Rates[i] = IC.Rates[i] * ExtNoise(noiselevel);					
				if(noisetype  == NOISYLOG)
					Rates[i] = IC.Rates[i] * pow(10., noiselevel * 2 *(RNDF()-0.5));
					/// Equals rate*= 10^RND(-max,+max)
			}
			
			for(int i=0; i<reactantcount; i++){
			    if(Reactants[i].type == NOISYPROTEIN) {
				Reactants[i].conc *= ExtNoise(p.sigma);
			    }			    
			}
			
		}
	}		
	//--------------------------------------------------------------------------------------------	
	void GetConc(double *y){
		int i=0;
		foreach(ReactantsMap::value_type &r, Reactants) 
				y[i++] = r.second.conc; 		
	}
	//--------------------------------------------------------------------------------------------
	void SetConc(double *y){
		int i=0;
		foreach(ReactantsMap::value_type &r, Reactants) { 
			//check if concentration is negative.
			if(y[i]< -0.005) { 
			  cout << i << ", " << y[i] << endl;
			  MyErr("SetConc: negative concentration!");	 
			}
			// check if concentration is too small to integrate and decreasing.
			// this happens when a reactant is not produced but only degrades.
			// ** I decided not to do this because it could cause unstability.
			// ** Cvode suggests keeping the
			//if(y[i] < 1e-18 && r.second.conc > y[i]) y[i] = 0.;   
			r.second.conc = y[i++]; 		
		}
	}
	//--------------------------------------------------------------------------------------------
	const Reactor &operator=(const Reactor &right){
		if(&right != this) {
			Reactions = right.Reactions;
			Reactants = right.Reactants;
			Rates 	  = right.Rates;
			NoiseLevels = right.NoiseLevels;
			NoiseTypes = right.NoiseTypes;
			ReactantCount = right.ReactantCount;
			FitnessScore = right.FitnessScore;
			Genes = right.Genes;
			Proteins = right.Proteins;
			PGComplexes = right.PGComplexes;
			PComplexes = right.PComplexes;
		}
		return *this;
	}
	//--------------------------------------------------------------------------------------------	
	int ReactantNo(string reactantname){
	//Returns number of a given reactant name. -1 if NULL.
		if(reactantname=="") return BLANK;
		map< int, Reactant, less<int> >::const_iterator iter;
		for(iter=Reactants.begin(); iter!=Reactants.end(); ++iter)
			if(iter->second.name == reactantname) return iter->first;
		
		cout << "Not found:" << reactantname << endl;		
		MyErr("Reactant Name Not Found");		
		return 0;
	}
	//--------------------------------------------------------------------------------------------	
	string PrintReactions(bool noisyonly=false, bool remark=false){
		ostringstream buffer;
		string s="+";
		
		for(unsigned int i = 0; i < Reactions.size(); i++){
				int a = Reactions[i].IN[0];
				int b = Reactions[i].IN[1];
				int c = Reactions[i].OUT[0];
				int d = Reactions[i].OUT[1];
				int ri = Reactions[i].RateIndex;
				double r = Rates[ri];
				int type=NoiseTypes[ri];
				int rc = Reactions[i].ReactionClass;				
				
				if(noisyonly && (rc == NOTNOISY && type != MUTATABLE)) // For runs w/ extr. noise, write only noisy reactions.
					continue;
				
				if(rc==HILLKINETICS) {
					double k = Rates[Reactions[i].Hill.kindex];
					double K = Rates[Reactions[i].Hill.Kindex];
					double n = Rates[Reactions[i].Hill.nindex];
									
					buffer << (remark?"# ":"") << (a>=0?Reactants[a].name:"") << 
						(b>=0?s + Reactants[b].name:"") << " > "<<
						(c>=0?Reactants[c].name:"") << 
						(d>=0?s + Reactants[d].name:"") << ",\t Hill(" << k << "*" << ", " 
														<< K <<  ","
														<< n <<  "), "
														<< rc << ",ri=" << ri <<  endl;	
				}
				else {
					buffer << (remark?"# ":"") << (a>=0?Reactants[a].name:"") << 
						(b>=0?s + Reactants[b].name:"") << " > "<<
						(c>=0?Reactants[c].name:"") << 
						(d>=0?s + Reactants[d].name:"") << ",\t" << r << ", " << rc << ",ri=" << ri <<  endl;	
				}
		}
		return buffer.str();
	}
	//--------------------------------------------------------------------------------------------
	void PrintReactants(bool nonzeroonly=false) {
		cout << "Reactants" << endl;
		foreach(ReactantsMap::value_type &v, Reactants) {
			Reactant r = v.second;
			int reactantno = v.first;

			if((nonzeroonly && r.conc <= 0) || r.type == GENE) continue; //Skip if concentration is zero.			
			cout << reactantno << ",\t\t" << r.name << ",\t\t" << r.type << ",\t\t" 
				<< r.produces << ",\t\tConc:" << r.conc << ",\t\tRecorded:" << r.recorded << endl;
		}
	}
	//--------------------------------------------------------------------------------------------
	void MutateRates(){
		map <unsigned int, double> sisters; // Store rates for sister reactions. (SisterNo, Rate)

		/*	foreach(Reaction &r, Reactions){
			if(RNDF() < p.MuReacRate)			
				r.Rate *= RNDF(2.);
		}*/
		for(unsigned int i=0; i < p.RatesToMutate.size(); i++) { 
			int reactionindex = p.RatesToMutate.at(i);
			
			if(RNDF() < p.MuReacRate){			
				//cout << "mutating " << reactionindex << endl;				
				double r = Rates.at(reactionindex);
				r *= RNDF(2.);	
				//Set limits for reaction rates.
				if(r > p.MaxRate) r=p.MaxRate;
				if(r < p.MinRate) r=p.MinRate;
				Rates.at(reactionindex)=r;
			}
		}
	}
	//--------------------------------------------------------------------------------------------
	int GetRates(double *r) { //Put the Rates map in an array.
		unsigned s = Rates.size();
		
		for(unsigned i=0; i<s; i++)  {
			r[i]=Rates[i];
			//cout << Rates[i] << ", ";
		}
	
		return s;
	}
	//--------------------------------------------------------------------------------------------
	void PutRates(double *r) { //Get the rates from an array, and put them in the Rates map.
		for(unsigned i=0; i<Rates.size(); i++) 
			Rates[i]=r[i];
	}
	//--------------------------------------------------------------------------------------------
	/*void AddNewTF_DNA_Interaction() {
		int tf, gene, product = BLANK;
		int trials = 5; 
		bool tryagain = true;
			
		Genes.insert(Genes.begin(), PGComplexes.begin(), PGComplexes.end());
		Proteins.insert(Proteins.begin(), PComplexes.begin(), PComplexes.end()); 
		int gcount = Genes.size();
		int pcount = Proteins.size();

		//@@@
		cout << "Genes: ";
		ostream_iterator<int> out_it (cout,", ");
		copy (Genes.begin(), Genes.end(), out_it );	
		cout << endl;
		//@@@
		
		while(tryagain && trials) {		
			tryagain = 0;	
			tf = Proteins[RNDF(pcount)];
			gene = Genes[RNDF(gcount)];
			cout << "RND pairs: tf=" << tf << ", gene=" << gene << endl;
			//Pick a random gene and protein. If they already interact pick another.
			foreach(Reaction &r, Reactions){
				if((r.IN[0]==tf && r.IN[1]==gene) || (r.IN[1]==tf && r.IN[0]==gene))
					{ tryagain = 1; trials--;}					
			}
		}	

		product = Reactants[gene].produces;	
		cout << "Pair:" << tf << ", " << gene << endl; 		
		cout << "Gene:" << gene << " Product: " << product << endl; //cin.get();
		if(!trials) { cout << "Can't add new reactions\n"; return; }
		if (product==BLANK) MyErr("product==BLANK");

		//Add 3 reactions to the reactor: 1) a+B->aB 2) aB->a+B 3) aB->aB+B
		//First add the new reactant aB to the reactant list.
		Reactant newreactant; string colon=":";
		int newreactantindex = Reactants.size();
		newreactant.name = Reactants[gene].name + colon + Reactants[tf].name;
		newreactant.type = PGCOMPLEX;
		newreactant.conc = 0;
		newreactant.produces = product;
		newreactant.recorded = 0;
		Reactants[newreactantindex] = newreactant;

		//Then add three reactions to the reactor:
		//First: a + B -> a:B
		Reaction dummyr;
		dummyr.IN[0]  = gene;
		dummyr.IN[1]  = tf;
		dummyr.OUT[0] = newreactantindex;
		dummyr.OUT[1] = BLANK;
		dummyr.Rate   = p.RandomDNATFAssocDissocRate(); 
		dummyr.ReactionClass = 0;
		Reactions.push_back(dummyr);			
		//Second: reverse reaction a:B -> a + B
		dummyr.IN[0] = newreactantindex;
		dummyr.IN[1] = BLANK;
		dummyr.OUT[0]  = gene;
		dummyr.OUT[1]  = tf;
		dummyr.Rate   = p.RandomDNATFAssocDissocRate(); 
		dummyr.ReactionClass = 0;
		Reactions.push_back(dummyr);			
		//Third: synthesis reaction a:B -> a:B + A
		dummyr.IN[0] = newreactantindex;
		dummyr.IN[1] = BLANK;
		dummyr.OUT[0]  = newreactantindex;
		dummyr.OUT[1]  = product;
		dummyr.Rate   = p.RandomSynthRate(); 
		dummyr.ReactionClass = 0;		
		Reactions.push_back(dummyr);			
		
		UpdateReactantClasses();
		cout << "New reactant: " << newreactant.name << endl;
		//PrintReactants();
		//PrintReactions();
		//cin.get();
	}*/
	//--------------------------------------------------------------------------------------------
	/*void AddNewDimerization() {
		int p1, p2;
		int trials = 5; 
		bool tryagain = true;
			
		Proteins.insert(Proteins.begin(), PComplexes.begin(), PComplexes.end()); 
		int pcount = Proteins.size();

		//@@@
		cout << "Proteins&PComplexes: ";
		ostream_iterator<int> out_it (cout,", ");
		copy (Proteins.begin(), Proteins.end(), out_it );	
		cout << endl;
		//@@@
		
		while(tryagain && trials) {		
			tryagain = 0;	
			p1 = Proteins[RNDF(pcount)];
			p2 = Proteins[RNDF(pcount)];
			cout << "RND pairs: p1=" << p1 << ", p2=" << p2 << endl;
			//Pick two random proteins. If they already interact pick another.
			foreach(Reaction &r, Reactions){
				if((r.IN[0]==p1 && r.IN[1]==p2) || (r.IN[1]==p1 && r.IN[0]==p2))
					{ tryagain = 1; trials--;}					
			}
		}	
		
		cout << "Pair:" << p1 << ", " << p2 << endl; 		
		if(!trials) { cout << "Can't add new reactions\n"; return; }

		//Add 4 new reactions to the reactor: 1) A+B->A:B   2) A:B->A+B
		//First add the new reactant A:B to the reactant list.
		Reactant newreactant; string colon=":";
		int newreactantindex = Reactants.size();
		newreactant.name = Reactants[p1].name + colon + Reactants[p2].name;
		
		newreactant.type = PCOMPLEX;
		newreactant.conc = 0;
		newreactant.produces = -1;
		newreactant.recorded = 1;
		Reactants[newreactantindex] = newreactant;

		//Then add FOUR NEW reactions to the reactor:
		//First, association: A + B -> A:B
		Reaction dummyr;
		dummyr.IN[0]  = p1;
		dummyr.IN[1]  = p2;
		dummyr.OUT[0] = newreactantindex;
		dummyr.OUT[1] = BLANK;
		dummyr.Rate   = p.RandomReactionRate(); 
		dummyr.ReactionClass = 0;		
		Reactions.push_back(dummyr);			
		
		//Second, dissociation: A:B -> A + B
		dummyr.IN[0] = newreactantindex;
		dummyr.IN[1] = BLANK;
		dummyr.OUT[0]  = p1;
		dummyr.OUT[1]  = p2;
		dummyr.Rate   = p.RandomReactionRate(); 
		dummyr.ReactionClass = 0;		
		Reactions.push_back(dummyr);	
		
		//Third, degradation A:B -> *
		dummyr.IN[0] = newreactantindex;
		dummyr.IN[1] = BLANK;
		dummyr.OUT[0]  = BLANK;
		dummyr.OUT[1]  = BLANK;
		dummyr.Rate   = p.RandomReactionRate(); 
		dummyr.ReactionClass = 0;		
		Reactions.push_back(dummyr);	
		
		UpdateReactantClasses();
		cout << " reactant: " << newreactant.name << endl;
		//PrintReactants();
		//PrintReactions();
		//cin.get();
	}*/
	//--------------------------------------------------------------------------------------------

	void UpdateReactantClasses () {
		Genes.clear(); Proteins.clear(); PComplexes.clear(); PGComplexes.clear(); 

		foreach(ReactantsMap::value_type &v, Reactants) {
			Reactant r = v.second;
			int reactantno = v.first;
			if(r.type == GENE)
				Genes.push_back(reactantno);
			else if(r.type == PROTEIN || r.type == NOISYPROTEIN)
				Proteins.push_back(reactantno);
			else if(r.type == PGCOMPLEX)
				PGComplexes.push_back(reactantno);
			else if(r.type == PCOMPLEX)
				PComplexes.push_back(reactantno);
			else MyErr("Reactant type unknown;");	
		}		
	
/*		ostream_iterator<int> out_it (cout,", ");	
		copy (Genes.begin(), Genes.end(), out_it);				
		cout << endl;
		copy (Proteins.begin(), Proteins.end(), out_it);				
		cout << endl;
		copy (PGComplexes.begin(), PGComplexes.end(), out_it);	
		cout << "--" <<endl;	*/		
	}
	
/*	void RandomizeRates() {
		unsigned size = Reactions.size();
		int oldclass = -1;
		double r;
				
		for(unsigned i=0; i < size; i++) {
			int newclass = Reactions[i].ReactionClass;
			if(oldclass != newclass) {
				oldclass = newclass;
				r = p.RandomReactionRate(); //Generate a new random rate for each reaction class
			}
			
			if(newclass != 0) {//Class 0 will not be randomized!!!
				int ri = Reactions[i].RateIndex; 
				Rates[ri] = r;	
			}
		}			
	} */
	
};

typedef	map <int, Reactor> ReactorMap;

//

//-------------------------------------------------------------------------------------------------
//int func (double t, const double y[], double f[], void *prm) {
static int f(realtype t, N_Vector y, N_Vector f, void *prm){
  Reactor ind = *(Reactor *) prm;
  double c1, c2;

  for(int i=0; i < (int) ind.Reactants.size() ;i++){
	  NV_Ith_S(f,i) = 0.;	
  }			
//````````````````
//   for(int i=0; i < (int) ind.Reactants.size(); i++) {
//      cout << "in f:  y[" << i << "]=" << NV_Ith_S(y,i) <<  ", " ;
//    }
//    cout << endl;
//````````````````

  for(int i=0; i < (int) ind.Reactants.size() ;i++){
    for(int j=0; j < (int) ind.Reactions.size() ;j++){
      double k;
      //cout << i << ", #k: " << k << endl;
	  int a = ind.Reactions[j].IN[0];
	  int b = ind.Reactions[j].IN[1];
	  int c = ind.Reactions[j].OUT[0];
	  int d = ind.Reactions[j].OUT[1];	  
	  c1 = (a >= 0. ? NV_Ith_S(y,a) : 1.);
	  c2 = (b >= 0. ? NV_Ith_S(y,b) : 1.);

		if(ind.Reactions[j].ReactionClass == HILLKINETICS){		
			///Hill kinetics
			k = ind.Reactions[j].CalcReacRate(c1, ind.Rates); 
			if(a == i) 	NV_Ith_S(f,i) -= k;
			if(b == i) 	NV_Ith_S(f,i) -= k;
			if(c == i)	NV_Ith_S(f,i) += k;				
			if(d == i)	NV_Ith_S(f,i) += k;				
		}		
		else{
			///Mass action
			k = ind.Rates[ind.Reactions[j].RateIndex];			
			//If you want to number reactants starting from 0, use y[a-1] and y[b-1].
			if(a == i) 	NV_Ith_S(f,i) -= k*c1*c2;
			if(b == i) 	NV_Ith_S(f,i) -= k*c1*c2;
			if(c == i)	NV_Ith_S(f,i) += k*c1*c2;				
			if(d == i)	NV_Ith_S(f,i) += k*c1*c2;				
		}
	  //printf("i:%i, reaction:%i, a:%i, b:%i, c:%i, d:%i, c1:%.5e, c2:%.5e, f[i]:%.5e \n", i, j, a, b, c, d, c1, c2, NV_Ith_S(f,i));
	  //cout << "i:" << i << ", reaction:" << j 
	  //<< ", c1: " << c1 << ", c2:" << c2 <<endl; 			
    }	
  }	
  return 0;
}
//-------------------------------------------------------------------------------------------------
int func (double t, const double y[], double f[], void *prm) {
	Reactor ind = *(Reactor *) prm;
	double c1, c2;
	for(int i=0; i < (int) ind.Reactants.size() ;i++){
		f[i]=0.;	
	}			

	for(int i=0; i < (int) ind.Reactants.size() ;i++){
		for(int j=0; j < (int) ind.Reactions.size() ;j++){
			int ri = ind.Reactions[j].RateIndex;
			double k = ind.Rates[ri];
			//cout << i << ", #k: " << k << endl;
			int a = ind.Reactions[j].IN[0];
			int b = ind.Reactions[j].IN[1];
			int c = ind.Reactions[j].OUT[0];
			int d = ind.Reactions[j].OUT[1];
			c1 = (a >= 0. ? y[a] : 1.);
			c2 = (b >= 0. ? y[b] : 1.);
			//If you want to number reactants starting from 0, use y[a-1] and y[b-1].
			if(a == i) 	f[i] -= k*c1*c2;
			if(b == i) 	f[i] -= k*c1*c2;
			if(c == i)	f[i] += k*c1*c2;				
			if(d == i)	f[i] += k*c1*c2;				
			//printf("i:%i, reaction:%i, a:%i, b:%i, c:%i, d:%i, c1:%.5e, c2:%.5e, f[i]:%.5e \n", i, j, a, b, c, d, c1, c2, f[i]);
			//cout << "i:" << i << ", reaction:" << j 
						//<< ", c1: " << c1 << ", c2:" << c2 <<endl; 			
		}	
	}	
				
	return GSL_SUCCESS;
}

