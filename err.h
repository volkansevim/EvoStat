inline void MyErr(const char *s, const string message="") {
	cout << "E R R O R  !\n"<< s << endl << message << endl ;
    #if defined(__PARALLEL__) 
		MPI_Finalize();
	#endif
	exit(-1);
}
