/* Tri Minh Ngo
Johns Hopkins University
10/07/2010
*/

/*
Mex compilation command
the -g flag is for debugging
mex fastmtpolyhedralFTcellpoly.cpp -g -I/usr/local/include/InsightToolkit/Utilities/vxl/core -I/usr/local/include/InsightToolkit/Utilities/vxl/vcl -I/usr/local/include/InsightToolkit/ -I/usr/local/include/InsightToolkit/Common -L/usr/local/lib/InsightToolkit/ -lITKCommon -lvcl -lvnl
mex fastmtpolyhedralFTcellpoly.cpp -g -I/usr/include/Utilities/vxl/core -I/usr/include/InsightToolkit/Utilities/vxl/vcl -I/usr/include/InsightToolkit/ -I/usr/include/InsightToolkit/Common -L/usr/lib/InsightToolkit/ -lITKCommon -lvcl -lvnl
*on the box
*mex fastmtpolyhedralFTcellpoly.cpp -I/usr/local/include/InsightToolkit/Utilities/vxl/core -I/usr/local/include/InsightToolkit/Utilities/vxl/vcl -I/usr/local/include/InsightToolkit/ -I/usr/local/include/InsightToolkit/Common -L/usr/local/lib/InsightToolkit/ -lITKCommon -litkvcl -litkvnl
*in windows
 *mex fastmtpolyhedralFTcellpoly.cpp -I"C:\Program Files (x86)\ITK\include\InsightToolkit\Utilities\vxl\core" -I"C:\Program Files (x86)\ITK\include\InsightToolkit\Utilities\vxl\vcl"
 * How to use the mex command
*fastmtpolyhedralFT(vertices,faces,pointsvector)
*/

#include "mex.h"
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_transpose.h>
#include <vcl_complex.h>
#include <vnl/vnl_cross.h>
#include <itkMultiThreader.h>
#include <vector>
#include <vnl/vnl_math.h>
#include <vnl/vnl_complex.h>
#include <cmath>
#include <string>

using namespace itk;
using namespace std;

double sinc(double x){
	return (x==0)?1:sin(vnl_math::pi*x)/(vnl_math::pi*x);
}

struct ThreadStruct{
	vector< vnl_matrix<double> >* r_v_comps;
	vector< vector< vector<unsigned int> > >* face_node_X_comps;
	vnl_matrix<double>* kvec;
	vector<double>* values;
	vnl_matrix< vcl_complex<double> >* S_star_partial;
};

struct SumStruct{
	double* datar;
	double* datai;
	vnl_matrix< vcl_complex<double> >* S_star_partial;
	vnl_matrix<double>* kvec;  
};

/*
class stringex: public exception{
    public:
    string mymsg;
    stringex(string msg):mymsg(msg){};
		const char* what() const throw() { return mymsg.c_str(); }
};
*/

/* 
struct ThreadStruct{
	vnl_matrix<double>* r;
	vnl_matrix<double>* k;
	vnl_matrix<double>* l;
	vnl_matrix<double>* r_nc;
	double* datai;
	double* datar;
};
*/
/* 
Calculate the Polygonal FT for a subset of k space samples in k matrix
[Inputs]
r: matrix holding vertices of the polygon.  Should be ordered counterclockwise. (3xN)
k: matrix holding sampling locations in kspace. (3xP)
startk: the index of the first k location to be calculated.
numksamps: the total number of k locations to be calculated.
[Outputs]
datar: the real part of the output
datai: the imaginary part of the output
*/

/*
Calculate face contribution to polygonal FT for a subset of faces
Assumption made that all faces have same number of vertices
[Inputs]
r_v: matrix holding vertices of the polyhedron.  There is no particular ordering (3xE)
kvec: matrix holding samplig locations in kspace (3xP)
face_node_X: matrix holding the indices of the vertices composing each face (3xF)
[Outputs]
S_star: matrix holding the contribution of each face for each k value (PxF)
*/

void calcFaceContribution(vector< vnl_matrix<double> >* r_v_comps, vector< vector< vector<unsigned int> > >* face_node_X_comps, vector<double>* values,
	vnl_matrix<double>* kvec, vnl_matrix< vcl_complex<double> >* S_star_partial,
	int threadCount, int threadID){
	vnl_vector_fixed<double,3> accumv;
	vcl_complex<double> j(0,1);	//imaginary number j
	vcl_complex<double> accum, solution;
// 	std::cout<<"Thread started\n";
	//number of components
	int C=r_v_comps->size();
	
	//per component pointers
	vector< vnl_matrix<double> >::iterator r_v;
	vector< vector< vector<unsigned int> > >::iterator face_node_X;
	
	//assume that all partial sums have been initialized to zero
	/*
	calculate associated parameters for this face
	N_f: matrix holding the unit normal vector of each face (3xF)
	n_f: a vector of matrices holding the unit normal to each edge composing a face (vector(F)x(3xE))
	t: a vector of matrices holding the directional unit vectors along each edge of a face (vector(F)x(3xE))
	r_c: a vector of matrices holding the position vector of the midpoint of each edge in a face (vector(F)x(3xE))
	*/
	for(int c=0; c<C; c++){
		face_node_X=face_node_X_comps->begin()+c;
		r_v=r_v_comps->begin()+c;
		//vector< vector< unsigned int> >& face_node_X=(*face_node_X_comps)[c];
		//vnl_matrix<double>& r_v=(*r_v_comps)[c];
		//divvy up faces to different threads
		int F = face_node_X->size();
		int chunksize = F/threadCount;
		int startface = threadID*chunksize;
		//last thread must take care of the rest, thread numbering begins at zero
		int numfaces = (threadID<threadCount-1)? chunksize:F-threadID*chunksize;
		
		//vector to hold kspace contribution for a component
		vnl_vector< vcl_complex<double> > kdata_comp(kvec->columns());
		
		kdata_comp.fill(0.0);
		for(int f=startface; f<startface+numfaces; f++){
			unsigned int E=(*face_node_X)[f].size()-1;
			//length of the edges of the face
			vnl_vector_fixed<double,3> L_hat; //says hat but not really unit length, should say L_bar or L_vec
			vector<double> L;
			vector< vnl_vector_fixed<double,3> > t;
			vector< vnl_vector_fixed<double,3> > n;
			vector< vnl_vector_fixed<double,3> > r_c;
			vnl_vector_fixed<double,3> N_f;
			//std::cout<<"face:"<<f<<std::endl;
			for(int e=0; e<E;e++){
				L_hat=r_v->get_column((*face_node_X)[f][e+1])-r_v->get_column((*face_node_X)[f][e]);
				L.push_back(L_hat.magnitude());
				t.push_back(L_hat/L[e]);
			}
			N_f=vnl_cross_3d(t[1],t[2]).normalize();
			for(int e=0; e<E;e++){	//this needs to be changed to allow varying number of edges per face.
				n.push_back(vnl_cross_3d(t[e],N_f));
				r_c.push_back(r_v->get_column((*face_node_X)[f][e])+t[e]*(L[e]/2));
			}
			for(int i=0;i < kvec->columns();i++){
				const vnl_vector<double>& k=kvec->get_column(i);
				//if(k.is_zero()){
				if(k.squared_magnitude()<1e-16){
					//calculate face contribution to volume of polyhedron
					accumv.fill(0);
					for(int e=0; e<E;e++)
						accumv+=vnl_cross_3d(r_v->get_column((*face_node_X)[f][e]),
							r_v->get_column((*face_node_X)[f][e+1]));
					kdata_comp(i)+=dot_product(r_v->get_column((*face_node_X)[f][0]),N_f)*
						abs(dot_product(N_f,accumv));
					//if(!isfinite(real(kdata_comp(i))) || !isfinite(imag(kdata_comp(i)))){
					//	std::cerr<< "nan or inf detected!\n";
					//}
				}
				else{
					//calculate appropriate face contribution
					if(vnl_cross_3d(k,N_f).squared_magnitude()<1e-16){
						//k is perpendicular to plane of face
						accumv.fill(0);
						for(int e=0; e<E;e++)
							accumv+=vnl_cross_3d(r_v->get_column((*face_node_X)[f][e]),
								r_v->get_column((*face_node_X)[f][e+1]));
                        //note below that factor of two is gone because expression for P_f contains (1/2) which cancels out factor of 2
						kdata_comp(i)+=-j*vnl_math::pi*dot_product(k,N_f)*
							abs(dot_product(N_f,accumv))*vcl_exp(-vnl_math::pi*2*j*dot_product(k,r_v->get_column((*face_node_X)[f][0])));
						//if(!isfinite(real(kdata_comp(i))) || !isfinite(imag(kdata_comp(i)))){
						//	std::cerr<< "nan or inf detected!\n";
						//}
					}
					else{
						//the usual contribution
						accum=0.0;
						for(int e=0; e<E;e++){
							accum+=L[e]*dot_product(k,n[e])*sinc(dot_product(k,t[e])*L[e])*
								vcl_exp(-vnl_math::pi*2*j*dot_product(k,r_c[e]));
								/*
								if(!isfinite(real(accum)) || !isfinite(imag(accum))){
									std::cerr << "nan or inf detected!\n";
									std::cerr << "L\n " << L[e];
									std::cerr << "k\n " << k;
									std::cerr << "n\n " << n[e];
									std::cerr << "r_c\n " << r_c[e];
								}
								*/
						}
						kdata_comp(i)+=(accum*dot_product(k,N_f))/(k.squared_magnitude()-pow(dot_product(k,N_f),2));
					}
				}
			}
		}
		//multiply partial sum by value of component
		(*S_star_partial).set_column(threadID,(*S_star_partial).get_column(threadID)+(kdata_comp*((*values)[c])));
	}
//   std::cout<<"Thread finished:"<<threadID <<std::endl;
}

ITK_THREAD_RETURN_TYPE ThreaderCallback( void *arg ){
	
	int threadId = ((MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
	int threadCount = ((MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;
	ThreadStruct* str = (ThreadStruct *)(((MultiThreader::ThreadInfoStruct *)(arg))->UserData);
	
//   std::cout<<"Thread started:"<<threadId <<std::endl;
	calcFaceContribution(str->r_v_comps, str->face_node_X_comps, str->values,
		str->kvec, str->S_star_partial,threadCount, threadId);
}
ITK_THREAD_RETURN_TYPE SumCallback( void *arg ){
	
	int threadID = ((MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
	int threadCount = ((MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;
	SumStruct* str = (SumStruct *)(((MultiThreader::ThreadInfoStruct *)(arg))->UserData);
	
//divvy up rows to different threads (turn this into a function later)
int N= (str->kvec)->columns();
int chunksize = N/threadCount;
int startcolumn = threadID*chunksize;
//last thread must take care of the rest, thread numbering begins at zero
int numcolumns = (threadID<threadCount-1)? chunksize:N-threadID*chunksize;  

	//calculate final values
	vcl_complex<double> S;
	
	for(int i=startcolumn;i<startcolumn+numcolumns;i++){
		S=(str->kvec)->get_column(i).is_zero()? abs((str->S_star_partial)->get_row(i).sum())/6:
			-(str->S_star_partial)->get_row(i).sum()/pow((2*vnl_math::pi*(str->kvec)->get_column(i).magnitude()),2);
		str->datar[i]=real(S);
		str->datai[i]=imag(S);
	}
}

void mexFunction(int nlhs, mxArray *plhs[], 
	int nrhs, const mxArray *prhs[])
{
	//prhs[0]: r_v_comps: cell array where components are rows
	//prhs[1]: face_node_X: cell array containing a cell array of face nodes (don't need the last node to repeat
	//prhs[2]: component values
	//prhs[3]: kvec
	//get the number of components
	int C = mxGetNumberOfElements(prhs[0]);  //get the total number of componenets
	vector< vnl_matrix<double> > r_v_comps(C);

	//this needs to be changed to a vector of vectors of vnl vectors.
	//vector< vnl_matrix<unsigned int> > face_node_X_comps(C);
	vector< vector< vector<unsigned int> > > face_node_X_comps(C);
	vector<double> values(C);
	
	double* temp=NULL;
	
	for(int c=0;c<C;c++){
		int N = mxGetN(mxGetCell(prhs[0],c));	//total number of vertices in polyhedron
		//int E = mxGetM(mxGetCell(prhs[1],c));  //vertices per polygonal face
		int F = mxGetNumberOfElements(mxGetCell(prhs[1],c));  //total number of faces in polyhedron
		//vcl_complex<double> j(0,1);	//imaginary number j
		
		//allocate storage
		//copy matlab data to vnl matrices and convert to 3D vectors
		//NOTE: Try to optimize the below
		r_v_comps[c].set_size(3,N);
		vnl_matrix<double>& r_v=r_v_comps[c];
		temp = mxGetPr(mxGetCell(prhs[0],c));
		for(int i=0;i<N;i++){
				r_v(0,i)=temp[i*3];
				r_v(1,i)=temp[i*3+1];
				r_v(2,i)=temp[i*3+2];
		}
		
		vector< vector<unsigned int> >& face_node_X=face_node_X_comps[c];
		face_node_X.reserve(F); //this line is optional?
		for(int f=0;f<F;f++){
			mxArray* cface=mxGetCell(mxGetCell(prhs[1],c),f);
			int E=mxGetNumberOfElements(cface);
			face_node_X.push_back(vector<unsigned int>());
			//face_node_X[f].reserve(E);  //this line is optional
			for(int e=0;e<E;e++){
				temp=mxGetPr(cface);
				face_node_X[f].push_back(temp[e]-1);
			}
			face_node_X[f].push_back(temp[0]-1);
		}

		values[c] = mxGetPr(prhs[2])[c];
	}
	int P = mxGetN(prhs[3]);	//total number of k space samples
	vnl_matrix<double> kvec(3,P);

	temp=mxGetPr(prhs[3]);
	for(int i=0;i<P;i++){
		kvec(0,i)=temp[i*3];
		kvec(1,i)=temp[i*3+1];
		kvec(2,i)=temp[i*3+2];
	}
	
	//setup multithreading
	MultiThreader::Pointer multithreader = MultiThreader::New();
	//multithreader->SetNumberOfThreads(16);
	vnl_matrix< vcl_complex<double> > 
		S_star_partial(P,multithreader->GetNumberOfThreads(),0.0); //partial sums for each thread
	cout<<multithreader->GetNumberOfThreads()<<endl;
	
	ThreadStruct iodata;
	iodata.r_v_comps=&r_v_comps;
	iodata.face_node_X_comps=&face_node_X_comps;
	iodata.values=&values;
	iodata.kvec=&kvec;
	iodata.S_star_partial=&S_star_partial;
	
	multithreader->SetSingleMethod(ThreaderCallback, &iodata);
	multithreader->SingleMethodExecute();
// 	std::cout << "done computing\n";
	//create the output matlab matrix
	plhs[0] = mxCreateDoubleMatrix(P, 1 , mxCOMPLEX);
	double* datar = mxGetPr(plhs[0]);
	double* datai = mxGetPi(plhs[0]);
	//multithread the column sum
//   SumStruct sumdata;
//   sumdata.datar=mxGetPr(plhs[0]);
//   sumdata.datai=mxGetPr(plhs[0]);
//   sumdata.S_star_partial=&S_star_partial;
//   sumdata.kvec=&kvec;
//   
//   multithreader->SetSingleMethod(SumCallback, &sumdata);
//   multithreader->SingleMethodExecute();
		//calculate final values
//   std::cout << "computing final sum\n";
	vcl_complex<double> S;
	for(int i=0;i<P;i++){
		S=kvec.get_column(i).is_zero()? abs(S_star_partial.get_row(i).sum())/6:
			-S_star_partial.get_row(i).sum()/pow((2*vnl_math::pi*kvec.get_column(i).magnitude()),2);
		datar[i]=real(S);
		datai[i]=imag(S);
	}
}
