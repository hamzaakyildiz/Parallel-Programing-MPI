/*
Student Name: HAMZA AKYILDIZ
Student Number: 2019400249
Compile Status: Compiling
Program Status: Working
Notes:  In mac even if it gives the correct outputs, it gives an some kind of system call errors.
		but "export OMPI_MCA_btl=self,tcp" execution of this line before running the program solves the problem.
		but it may not be the case in linux environment.
*/
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <mpi.h>
#include <cmath>
#include <iterator>
#include <set>
using namespace std;
//allocates 2D array dynamically
double **dynamic2DArray(int rows, int cols);

//calculates distance between two instances
double manhattanDistance(double* row1, double* row2, int n);

// it has a memory with minA and maxA updated by every call with different row such that
// at the end of the firs iteration we have the maxA and minA that keeps the values
void min_max(double* row, double* minA, double* maxA, int n);

//sorts by selection sort
void sort(int *arr, int size);


int main(int argc, char *argv[])
{
	int rank, numprocs;
	int N,A,M,T;
	int dataPart;
	double** dataSet;
	int* classTag;

	//initializes the MPI environment
	MPI_Init(&argc,&argv); 
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(rank==0){
	//master part
		fstream fin;
		fin.open(argv[1], ios::in);
		string line;
		getline(fin,line); // first line
		int numworker=stoi(line)-1; //number of slaves
		// reding fin line into line string
		getline(fin,line);
		stringstream s(line);

		//reading second line
		for(int i=0;getline(s,line,'\t');i++){
			//stoi is a string to integer converter
			if(i==0)
				N=stoi(line);
			else if(i==1)
				A=stoi(line);
			else if(i==2)
				M=stoi(line);
			else
				T=stoi(line);
		}

		dataPart=N/numworker; // how many lines each slave will process
		dataSet=dynamic2DArray(N,A); //memory space for instances features
		classTag=new int[N]; // memory space for class tags

		for(int i=0;i<N;i++)//reading input file into the programs relevant memory spaces
		{	
			getline(fin,line);
			stringstream s(line);
			string str;
			for (int j = 0; j < A; ++j)
			{
				getline(s,str,'\t');
			 	dataSet[i][j]=stod(str);//features
			}
			getline(s,str,'\t');
			classTag[i]=stoi(str);//class tags

		}
		// broadcasting the values that slaves need
		// key part that slaves also need to call MPI_Bcast since int he function source is declared by rank 0, master process and 
		// every body gets the value that master has
		MPI_Bcast(&dataPart, 1, MPI_INT, 0, MPI_COMM_WORLD); 
		MPI_Bcast(&A, 1, MPI_INT, 0, MPI_COMM_WORLD); 
		MPI_Bcast(&T, 1, MPI_INT, 0, MPI_COMM_WORLD); 
		MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);


		for(int dest =1; dest<numprocs;dest++){//sending the data part and relevant class tags to slaves with offset
			//these are blocking functions so until slave gets the message master stops
			MPI_Send(&dataSet[(dest-1)*(dataPart)][0],A*dataPart, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
			MPI_Send(&classTag[(dest-1)*(dataPart)],dataPart, MPI_INT,dest,1,MPI_COMM_WORLD);
		}
		//set is used for unites the outputs since it does not have same value twice
		set<int, greater<int> > output;
		int* index=new int[T];
		for (int recv = 1; recv < numprocs; recv++)// recieving outputs from slaves and putting them into a set
		{
			// this is also blocking recieve such that we dont need extra synchronization since master will wait until slave numprocs-1 sends the message
			MPI_Recv(&index[0],T,MPI_INT,recv,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for (int i = 0; i < T; ++i)
			{
				output.insert(index[i]);
			}
		}
		//reverse_iterator is used to arranging the printing order. Not reverse iterator iterateS the values descending order but output should be ascending order.
		set<int, greater<int> >::reverse_iterator itr;
		cout<<"Master P0 :";
		for (itr = output.rbegin();itr != output.rend(); ++itr) 
    	{

        	cout <<" "<< *itr;
    	}
    	cout << endl;
    	// closing the file that was opened
		fin.close();
	}else{
		//slaves call MPI_Bcast to gets the broadcasted value
		MPI_Bcast(&dataPart, 1, MPI_INT, 0, MPI_COMM_WORLD); 
		MPI_Bcast(&A, 1, MPI_INT, 0, MPI_COMM_WORLD); 
		MPI_Bcast(&T, 1, MPI_INT, 0, MPI_COMM_WORLD); 
		MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
		//allocating 2D array 
		dataSet=dynamic2DArray(dataPart,A);
		classTag=new int[dataPart];

		//recieves the data that will be processed in the slaves. Every slaves gets different part of the data but they seem to put into same address. 
		//However, they all have different memory spaces such that it is not the same memory address at all 
		MPI_Recv(&(dataSet[0][0]), A*dataPart, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&classTag[0],dataPart,MPI_INT,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		

		double* W = new double[A];
		double* minA = new double[A];
		double* maxA=new double[A];
		for (int i = 0; i < A; ++i)
		{

			W[i]=0;
			minA[i]=dataSet[0][i];
			maxA[i]=dataSet[0][i];
		}

		for(int i=0;i<M;i++)
		{
			double* row1=dataSet[i];//target instance
			double* row2;// instance that will be used to calculate distance from
			double nearestHit=100000000,nearestMiss=100000000,distance;// to find the minimums, at first i initialize the distance very far
			int hitIndex=0,missIndex=0;
			for (int j = 0; j < dataPart; j++)
			{
				row2=dataSet[j];// second instance
				if (i==0)//this is sufficinet to run this kine only for first iteration
					min_max(row2,minA,maxA,A);//updates the minA and maxA according to row2
				
				if(i==j)// it does not calculate distance from itself
					continue;
				
				distance=manhattanDistance(row1,row2,A);//distance calculations
				if(classTag[i]==classTag[j])
				{
					if(nearestHit>distance){//checking which is nearest if it is a hit
						nearestHit=distance;
						hitIndex=j;
					}

				}else if(nearestMiss>distance){//checking which is nearest if it is a miss
					nearestMiss=distance;
					missIndex=j;
				}
			}
			double numeratorHit,numeratorMiss,denominator;
			for (int i = 0; i < A; ++i)
			{
				numeratorHit=abs(row1[i]-dataSet[hitIndex][i]);//numerator of diff function for hit
				numeratorMiss=abs(row1[i]-dataSet[missIndex][i]);//numerator of diff function for mis
				denominator=(maxA[i]-minA[i])*M;//(maxA[i]-minA[i]) is the denominator of diff function and M is defined in this way for W updates
				W[i]=W[i]-(numeratorHit/denominator)+(numeratorMiss/denominator);
			}
		}

		double* output = new double[T];
		int* outputIndex = new int[T];
		for (int i = 0; i < T; ++i)
		{
			output[i]=-10;
			outputIndex[i]=-10;
		}
		for (int i = 0; i < T; ++i)//finds the biggest T weights
		{
			for(int j = 0 ; j<A;j++){
				if(W[j]>output[i]){//biggest weight is put the first index of the output if so, the index of the feature is put into outputindex
					output[i]=W[j];
					outputIndex[i]=j;
				}
			}
			W[outputIndex[i]]=-1000000;//if biggest is put into output, not to chosee it twice it is declared ver small number
		}
		sort(outputIndex,T);//for the sake of the printing order, sorts the output indexes of the features
		cout<<"Slave P"<<rank<<" :";
		for (int i = 0; i < T; ++i)
		{
			cout<<" "<<outputIndex[i];
		}
		cout<<" "<<endl;
		//slaves send results that are selected feature in their partitions to the master process
		MPI_Send(&outputIndex[0],T,MPI_INT,0,1,MPI_COMM_WORLD);
	}
	//finalize the MPI environment
	MPI_Finalize();
	return 0;
}
//distance calculation
double manhattanDistance(double* row1, double* row2, int n){
	double tot=0;
	for(int i=0; i <n;i++){
		tot+=abs(row1[i]-row2[i]);
	}
	return tot;
}

//updates minA and maxA. if the row has smaller value than minA in the same index, minA will be updated.
//if the row has greater value than maxA in the same index, maxA will be updated.
void min_max(double* row, double* minA, double* maxA, int n){
	for (int i = 0; i < n; ++i)
	{
		if (minA[i]>row[i])
		{
			minA[i]=row[i];
		}
		if (maxA[i]<row[i])
		{
			maxA[i]=row[i];
		}
	}
}

//it allocates 1D array with size rows*cols and assigns it addreses into 2D array
double** dynamic2DArray(int rows, int cols) {
  int i;
  double *data = (double *)malloc(rows*cols*sizeof(double));
  double **array= (double **)malloc(rows*sizeof(double*));
  for (i=0; i<rows; i++)
    array[i] = &(data[cols*i]);
  return array;
}

//selection sort
void sort(int *arr,int size){

	for (int i = 0; i < size; ++i)
	{
		double target;
		for (int j = i; j < size; ++j)
		{
			if (arr[j]<arr[i])
			{
				swap(arr[j],arr[i]);
			}
		}
	}
}


