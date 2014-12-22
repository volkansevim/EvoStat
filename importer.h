boost::regex remark("^#");
//*************************************************************************************************
int ImportReactants(const char *filename, Reactor &individual, int verbose=0){
	bool recording = false;
	std::string line, a, b, c, d, e;
	boost::regex pat6(" *([\\w:]+) *, *(\\w+) *, *([0-9]*\\.?[0-9]+) *"); 
	//Matches "a:B, PGCOMPLEX, 0.1". 
	//This is the old format of reactants.dat. ReactantName, ComplexType, InitialConcentration
	boost::regex pat7(" *([\\w:]+) *, *(\\w+) *, *([0-9]*\\.?[0-9]+) *, *(\\w) *"); 
	//Matches "a:B, PGCOMPLEX, 0.1". 
	//This is the new format of reactants.dat. Name, Type, InitialConcentration, Recorded?

	ifstream file (filename);
	Reactant dummy;
	p.ReactantsToRecord.clear();
	p.ReactantsToRecordNames.clear();
	
	if (!file.is_open()) MyErr("Can't open reactant file.");
	int counter=0;
	while (file.good()) {
	    std::getline(file, line);
	    if (line.length() == 0) continue; 
	    boost::smatch matches;

		if (boost::regex_search(line, remark)) continue; //skip remarks in the file
	    if (boost::regex_match(line, matches, pat6) || boost::regex_match(line, matches, pat7)) {
			//I support both old and new formats of the reactant.dat for compatibility.
	        a = matches[1]; b = matches[2]; c=matches[3]; 
		
			dummy.conc = atof(c.c_str());; 
			dummy.name = a; 
			dummy.type = p.ReactantTypes[b]; //
			dummy.produces = BLANK;		
			dummy.recorded = 0;
			
			if(p.stochastic && fmod(dummy.conc,1)>0) 
				MyErr("Initial concentration has to be an integer for stochastic runs.");
			
			if (boost::regex_match(line, matches, pat7)) {// New format?
				d=matches[4];
				boost::to_lower(d);
				if(d == "r") { //Record the time series for this reactant. 
					p.ReactantsToRecord.push_back(counter); //vector containing recorded reactant #s
					p.ReactantsToRecordNames.push_back(dummy.name);
					dummy.recorded = true;
					recording = true;
				}
			}
			individual.Reactants.insert(pair<int, Reactant>(counter, dummy));			
			//cout << "\t" << dummy.name << " ; " << "(" << b << ")" << dummy.type  << " ; " << endl;
			counter++;
		}
		else { 
			cout << "Not matched: " << line << endl;
			getline(cin, line) ;
		}
	}
	if(!recording) MyErr("Need to record at least one reactant!"); 

	if(verbose){
		cout << "ReactansToRecord: ";
		foreach(int r, p.ReactantsToRecord) 
			cout << individual.Reactants[r].name << ", ";
		cout << endl;
	}
	file.close();
	individual.UpdateReactantClasses();
	return 0;
}
//#define MUTATABLE 99
//*************************************************************************************************
int ImportReactions(const char *filename, Reactor &individual, int verbose=0){
	std::string line, a, b, c, d, e;
	double rate=0, noiselevel=0;
	int noisetype = 0;
	//map <int, int> sisters;
	//map <int, int>::iterator it;
	map <string, int> variablenames;
	map <string, int>::iterator it;	
	unsigned int rindex=0; //count the number of non-sister reactions
	
	string word = "([\\w:]+)", floatn = "([0-9]*\\.?[0-9]+)", comma=",";
	//string threefloats= floatn+comma+floatn+comma+floatn;
	string threewords = word+comma+word+comma+word;
	string Hill = "(.*),Hill\\(";
	Hill += threewords;
	Hill += "\\)";
	
// 	string dum="funkyvar";
// 	variablenames[dum]=0;
// 	dum="otherfunkyvar";
// 	variablenames[dum]=1;
// 	rindex=2;
// 	individual.Rates[0]=0.001;
// 	individual.Rates[1]=100;
	
	Reaction dummyr;
    boost::regex pat1("([\\w:]+)\\+([\\w:]+)>([\\w:]+)\\+([\\w:]+)" ); 
	///Matches "reactant1+reactant2>reactant3+reactant4,rate". 
    boost::regex pat2("([\\w:]+)\\+([\\w:]+)>([\\w:]+)" ); 
	///Matches "reactant1+reactant2>reactant3,rate". 
    boost::regex pat3("([\\w:]+)>([\\w:]+)\\+([\\w:]+)" ); 
	///Matches "reactant1>reactant2+reactant3,rate". 
    boost::regex pat4("([\\w:]+)>([\\w:]+)" ); 
	///Matches "reactant1>reactant2,rate". 
    boost::regex pat5("([\\w:]+)>" ); 
	///Matches "reactant1>,rate". (Degradation of reactant1)
	boost::regex pat6("(.*);(\\d+)@([0-9]*\\.?[0-9]+)" ); 
	///Semicolon seperates reaction class specification.
	boost::regex pat7("(.*);(\\d+),(\\d+)" ); 
	///Semicolon seperates reaction class and sister reaction specifications.
	boost::regex pat8(Hill); 	
	///Matches "*, Hill(k,K,n)"
	boost::regex pat9("^variable(.*)" ); 
	///Matches a variable declaration
	boost::regex pat10("(.*),(.*)" ); 
	///Matches a reaction declaration
	boost::regex pat11("([\\w:]+)=([0-9]*\\.?[0-9]+)" ); 
	///Matches "variablename=FLOAT"
	boost::regex pat12("([\\w:]+)=([0-9]*\\.?[0-9]+),mutatable" ); 
	///Matches "variablename=FLOAT, mutatable"	
	boost::regex pat13("([\\w:]+)=([0-9]*\\.?[0-9]+),([\\w:]+),([0-9]*\\.?[0-9]+)"); 
	///Matches "variablename=FLOAT, linear, 0.2"	
	boost::regex patword(word); 	
	///Matches a word
	boost::regex patfloat(floatn); 
	///Matches a float
	
	ifstream myfile (filename);
	if (!myfile.is_open()) MyErr("Can't open reaction file.");
	p.RatesToMutate.clear();
	
    while (myfile.good()) {
        string name, namek, nameK, namen;
		getline(myfile, line);
		boost::regex wspace("\\s+",boost::regex_constants::icase|boost::regex_constants::perl); 
		string clear("");		
		line  = boost::regex_replace(line, wspace, clear);		
		//line is now free of whitespace. 		
		
        if (line.length() == 0) continue; 
		if (boost::regex_search(line, remark)) continue; //skip remarks in the file
        boost::smatch matches;
		a=b=c=d=e="";
		int rclass = MASSACTION;		

		///-------- VARIABLE DECLARATION
		if (boost::regex_match(line, matches, pat9)) { ///Variable declaration or reaction?
			line=matches[1];
			if (boost::regex_match(line, matches, pat11)) { ///Match "rate1=2.02"
				a = matches[1]; b = matches[2]; 
				//cout << "Variable: " << a << "=" << b<< endl;
				name=a;
				rate=atof(b.c_str());
				noisetype=0;
				noiselevel=0;
			}
			else if (boost::regex_match(line, matches, pat12)) { ///Match "rate1=0.1, mutatable"
				a = matches[1]; b = matches[2];
				//cout << "Mutatable: " << a << "=" << b<< endl;
				name=a;
				rate=atof(b.c_str());
				noisetype=MUTATABLE;
				noiselevel=0;				
			}
			else if (boost::regex_match(line, matches, pat13)) { ///Match "rate1=0.1, linear ,0.2
				a = matches[1]; b = matches[2]; c=matches[3]; d=matches[4]; 	
				//cout << "in match: " << a << "|" << b << "|" << c << "|" << d << endl;
				string strlin("linear"); 
				string strloglin("loglinear");
				name=a;
				rate=atof(b.c_str());
				if (c.compare(strlin)==0) { /// Linear noise
					noisetype=NOISYLIN;
					noiselevel=atof(d.c_str());
					if(noiselevel>1 || noiselevel<0) MyErr("Linear NoiseLevel is out of bounds.");
				}
				else if (c.compare(strloglin)==0) { /// Log-linear noise
					noisetype=NOISYLOG;
					noiselevel=atof(d.c_str());									
					if(noiselevel<0) MyErr("Log-Linear noiselevel is out of bounds.");
				}
				else MyErr("Can't match noisy variable declaration.");
				
			}
			else {
				cout << line << endl;
				MyErr("Not matched! Faulty variable declaration.");
			}		
			
			variablenames[a]=rindex;
			individual.Rates[rindex]=rate;	
			individual.NoiseTypes[rindex]=noisetype;
			individual.NoiseLevels[rindex]=noiselevel;				
			rindex++;
			continue;									
		}
		///-------- REACTION DECLARATION	
		else if (boost::regex_match(line, matches, pat8)){ ///Hill Kinetics
			a = matches[1]; b = matches[2]; c=matches[3]; d=matches[4];
			rclass = HILLKINETICS;
			namek = b;
			nameK = c;
			namen = d;
			line = a; ///strip the Hill func part from the line
			//line += ",0";///React. w/ Hill kinetics don't have a specified rate. Concatanate ",0" for compatibility.
			cout << "Hill kinetics: " << namek << ", " << nameK << ", " << namen << endl;
			a=b=c=d=e="";
		}
		else if (boost::regex_match(line, matches, pat10)) { ///Mass Action: split at the comma.
			a = matches[1]; b = matches[2];
			line = a;
			std::string aftercomma=b;
			
			if(boost::regex_match(aftercomma, matches, patfloat)){ ///Is the rate a number?
				a = matches[1];
				rate = atof(a.c_str());
			}
			else if (boost::regex_match(aftercomma, matches, patword)) { /// Is it a variable name?
				a = matches[1];				
				name=a;
			}
			else { /// Faulty declaration
				cout << line << endl;
				MyErr("Not matched! Faulty reaction rate decleration.");
			}
		}
		else { /// Not a reaction or a variable declaration.
			cout << line << endl;
			MyErr("Can't match line!", line);
		}

		///Get the part before comma
		a=b=c=d=e="";
		if (boost::regex_match(line, matches, pat1)) {
			a = matches[1]; b = matches[2]; c=matches[3]; d=matches[4]; 				
			//cout << "\t" << a << " + " << b << " -> " << c << " + " << d << " : " << e << endl;
		}
		else if (boost::regex_match(line, matches, pat2)) {
			a = matches[1]; b = matches[2]; c=matches[3]; 				
			//cout << "\t" << a << " + " << b << " -> " << c << " + " << d << " : " << e << endl;
		}
		else if (boost::regex_match(line, matches, pat3)) {
			a = matches[1]; c = matches[2]; d = matches[3]; 				
			//cout << "\t" << a << " + " << b << " -> " << c << " + " << d << " : " << e << endl;
		}
		else if (boost::regex_match(line, matches, pat4)) {
			a = matches[1]; c = matches[2]; 				
			//cout << "\t" << a << " + " << b << " -> " << c << " + " << d << " : " << e << endl;
		}
		else if (boost::regex_match(line, matches, pat5)) {
			a = matches[1]; 				
			//cout << "\t" << a << " + " << b << " -> " << c << " + " << d << " : " << e << endl;
		}
		else { /// Faulty 
			cout << line << endl;
			MyErr("Not matched! Faulty reaction.");
		}
		
		
		//CC variability?			
/*		if(p.sigma > 0) 
			if(rclass == 1) 
				rate *= p.ccvariability; //Add extrinsic noise to transcription (Class 1). */
		
		int in0, in1, out0, out1;
		dummyr.IN[0]  = in0 = individual.ReactantNo(a); 
		dummyr.IN[1]  = in1 = individual.ReactantNo(b); 
		dummyr.OUT[0] = out0 = individual.ReactantNo(c); 
		dummyr.OUT[1] = out1 = individual.ReactantNo(d); 
		dummyr.ReactionClass = rclass;

		/// Is it a Hill kinetics reaction?
		if(rclass==HILLKINETICS){										
			it = variablenames.find(namek);					
			if(it != variablenames.end()) 
				dummyr.Hill.kindex = variablenames[namek];	
			else { cout << "\n" <<  namek << endl; MyErr("Variable not defined."); }
			
			it = variablenames.find(nameK);					
			if(it != variablenames.end()) 
				dummyr.Hill.Kindex = variablenames[nameK];	
			else { cout << "\n" << nameK << endl; MyErr("Variable not defined."); }

			it = variablenames.find(namen);					
			if(it != variablenames.end()) 
				dummyr.Hill.nindex = variablenames[namen];	
			else { cout << "\n" << namen << endl; MyErr("Variable not defined."); }
			
			cout << "Imported Hill\n";
			cout << individual.Rates[dummyr.Hill.kindex] << endl;
			cout << individual.Rates[dummyr.Hill.Kindex] << endl;
			cout << individual.Rates[dummyr.Hill.nindex] << endl;
			
		}
		else {	/// Mass action
			///Is the reaction rate a variable?
			if(!name.empty()) { 
				it = variablenames.find(name);					
				if(it != variablenames.end()) { //Is the name defined?
					dummyr.RateIndex = variablenames[name]; //Yes. Use its index.
					//dummyr.ReactionClass = individual.NoiseTypes[dummyr.RateIndex];
				}
				else {
					cout << name << endl;
					MyErr("Variable not defined.");
				}			
			}
			///The rate is given as a float.
			else { 
				dummyr.RateIndex = rindex;
				dummyr.ReactionClass = 0;
				individual.Rates[rindex]=rate; //also record the rate	
				individual.NoiseTypes[rindex]=0;
				individual.NoiseLevels[rindex]=0;				
				rindex++;
			}
		}
				
		individual.Reactions.push_back(dummyr);			
		//if this is a synthesis reaction, mark the gene or PGcomplex
		if (boost::regex_match(line, matches, pat3)) {//check if reaction looks like a->b+c
			if(in0 == out0 || in0 == out1) //check if it looks like a->a+B or a->B+a
				if( (individual.Reactants[in0].type == GENE) ||
									 (individual.Reactants[in0].type == PGCOMPLEX) ){
					int product = out1; 
					if(product == in0) product=out0; //in case reaction looks like a->B+a
					individual.Reactants[in0].produces = product;
				}
		}					
		
    }
	
	//cout << "NOISE TYPES\n";
	IDMap::iterator ita;
	for(ita=individual.NoiseTypes.begin(); ita!=individual.NoiseTypes.end(); ita++) {	
		//cout << ita->first << " " << ita->second << endl;
		if(ita->second==MUTATABLE)
			p.RatesToMutate.push_back(ita->first);
	}
	
	if(verbose){
		cout << "Mutatables:\n";
		for(unsigned int i=0; i<p.RatesToMutate.size(); i++){ cout << p.RatesToMutate.at(i) << ", ";}
		cout << endl;

		cout << "Rates:\n";
		for(unsigned int i=0; i<individual.Rates.size(); i++) { cout << individual.Rates[i] <<  ", "; }
		cout << endl;
		
		cout << "NoiseTypes & Levels:\n";
		for(unsigned int i=0; i<individual.NoiseTypes.size(); i++) { 
			cout << i << ":" << individual.NoiseTypes[i] << "@" << individual.NoiseLevels[i] <<  ",   ";
		}
		cout << endl;

	}
    myfile.close();	
	return 0;
}
