#g++ -Wall -ggdb3 -I/usr/local/include/sundials/  -c  $1.cpp 
mpicxx -g -O3 -g -Wall -O3 -I/usr/local/include/sundials/  -c  $1.cpp 
#g++ -L/usr/local/lib -static -pthread $1.o libsundials_cvode.a libsundials_nvecserial.a -lboost_regex -lgsl -lgslcblas -lm 
###!!! mpicxx -g -O3 -g -Wall -O3 -Wl,-Bsymbolic-functions -Wall -L/usr/local/lib $1.o -lsundials_cvode -lsundials_nvecserial -lboost_regex 
mpicxx -g -O3 -g -Wall -O3 -Wl,-Bsymbolic-functions -Wall $1.o -lsundials_cvode -lsundials_nvecserial -lboost_regex -lgsl -llapack  
## I removed -static -pthread
rm $1.o
