#include<algorithm>
#include<vector>
#include<iostream>
#include<Eigen/Eigenvalues>
#include<fstream>
#include<math.h>

using namespace Eigen;
using namespace std;

typedef struct Point3D
{
	float x, y, z;
};

typedef struct Rotation
{
	float x1, x2, x3, y1, y2, y3, z1, z2, z3;
};

float sq(float x) { return x*x; }

void print(vector<Point3D> &P)
{
    vector<Point3D>::iterator it;
    for (it=P.begin(); it!=P.end(); it++)
    {
        cout<<it->x<<" "<<it->y<<" "<<it->z<<endl;
    }
}

void CalculateMeanPoint3D(vector<Point3D> &P, Point3D &mean)
{
	vector<Point3D>::iterator it;
	mean.x = 0;
	mean.y = 0;
	mean.z = 0;
	for (it = P.begin(); it != P.end(); it++)
	{
		mean.x += it->x;
		mean.y += it->y;
		mean.z += it->z;
	}
	mean.x = mean.x / P.size();
	mean.y = mean.y / P.size();
	mean.z = mean.z / P.size();
}

void MovingPointSet(vector<Point3D> &P, Point3D &T)
{
	vector<Point3D>::iterator it;
	for (it = P.begin(); it != P.end(); it++)
	{
		it->x += T.x;
		it->y += T.y;
		it->z += T.z;
	}
}

void FindCorrespondingPoint(vector<Point3D> &P, vector<Point3D> &Q, vector<Point3D> &X, float &E) {

	vector<Point3D>::iterator itp, itq;
	E = 0;
	for (itp = P.begin(); itp != P.end(); itp++) {
		float mindist = -1;int qcount=0; int pos=0;
		for (itq = Q.begin(); itq != Q.end(); itq++) {
            qcount++;
			float dist = sq(itq->x - itp->x) + sq(itq->y - itp->y) + sq(itq->z - itp->z);
			if (mindist<0 || mindist>dist) { mindist = dist; pos=qcount;}
		}
		E += mindist / P.size();//cout<<E<<" "<<(E>100)<<endl;
		qcount=1;itq=Q.begin();
		while(qcount!=pos)
		{
		    itq++; qcount++;
		}
		X.push_back(*itq);
		Q.erase(itq);
	}

}

void CalculateRotation(vector<Point3D> &P, vector<Point3D> &X, Rotation &R, Point3D P_mean, Point3D X_mean) {
	//build covariance matrix between P and X
	float cov[9] = { 0,0,0,0,0,0,0,0,0 };
	vector<Point3D>::iterator itx, itp;	for (itp = P.begin(), itx = X.begin(); itp != P.end(); itp++, itx++) {
		cov[0] += (itp->x-P_mean.x)*(itx->x-X_mean.x) / P.size();//(1,1)
		cov[1] += (itp->x-P_mean.x)*(itx->y-X_mean.y) / P.size();//(1,2)
		cov[2] += (itp->x-P_mean.x)*(itx->z-X_mean.z) / P.size();//(1,3)
		cov[3] += (itp->y-P_mean.y)*(itx->x-X_mean.x) / P.size();//(2,1)
		cov[4] += (itp->y-P_mean.y)*(itx->y-X_mean.y) / P.size();//(2,2)
		cov[5] += (itp->y-P_mean.y)*(itx->z-X_mean.z) / P.size();//(2,3)
		cov[6] += (itp->z-P_mean.z)*(itx->x-X_mean.x) / P.size();//(3,1)
		cov[7] += (itp->z-P_mean.z)*(itx->y-X_mean.y) / P.size();//(3,2)
		cov[8] += (itp->z-P_mean.z)*(itx->z-X_mean.z) / P.size();//(3,3)
	}
        cout<<"cov: "<<endl;
        Matrix3f B;
        B  << cov[0], cov[1], cov[2], cov[3], cov[4], cov[5], cov[6], cov[7], cov[8];
        cout<<B<<endl;
	//build 4*4 symetric matrix
	float Q[16];
	//first line
	Q[0] = cov[0] + cov[4] + cov[8];
	Q[1] = cov[5] - cov[7];
	Q[2] = cov[6] - cov[2];
	Q[3] = cov[1] - cov[3];
	//second line
	Q[4] = cov[5] - cov[7];
	Q[5] = cov[0] - cov[4] - cov[8];
	Q[6] = cov[1] + cov[3];
	Q[7] = cov[6] + cov[2];
	//third line
	Q[8] = cov[6] - cov[2];
	Q[9] = cov[1] + cov[3];
	Q[10] = cov[4] - cov[0] - cov[8];
	Q[11] = cov[5] + cov[7];
	//last line
	Q[12] = cov[1] - cov[3];
	Q[13] = cov[2] + cov[6];
	Q[14] = cov[5] + cov[7];
	Q[15] = cov[8] - cov[0] - cov[4];

	Matrix4f A;
	A << Q[0], Q[1], Q[2], Q[3], Q[4], Q[5], Q[6], Q[7], Q[8], Q[9], Q[10], Q[11], Q[12], Q[13], Q[14], Q[15];
	cout<<A<<endl;
	EigenSolver<Matrix4f> es(A);

	Matrix4f D = es.pseudoEigenvalueMatrix();
	Matrix4f V = es.pseudoEigenvectors();

	float biggestValue=A(0,0);
	int pos=0;
	for (int i=1; i<4; i++){
        if (biggestValue<D(i,i)){biggestValue=D(i,i);pos=i;}
	}
    float q[]={V(0,pos), V(1,pos), V(2,pos), V(3,pos)};

    cout<<"D: "<<endl;cout<<D<<endl;
    cout<<"V: "<<endl;cout<<V<<endl;
    cout<<"unit quaternion:"<<endl;
    cout<<q[0]<<" "<<q[1]<<" "<<q[2]<<" "<<q[3]<<endl;

	//calculate rotation matrix with unit quaternion

	R.x1 = q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];//(1,1)

	R.y1 = 2*(q[1]*q[2]-q[0]*q[3]);//(1,2)

	R.z1 = 2*(q[1]*q[3]+q[0]*q[2]);//(1,3)

	R.x2 = 2*(q[1]*q[2]+q[0]*q[3]);//(2,1)

	R.y2 = q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];//(2,2)

	R.z2 = 2*(q[2]*q[3]-q[0]*q[1]);//(2,3)

	R.x3 = 2*(q[1]*q[3]-q[0]*q[2]);//(3,1)

	R.y3 = 2*(q[2]*q[3]+q[0]*q[1]);//(3,2)

	R.z3 = q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];//(3,3)
}

void Rotate(vector<Point3D> &P, Rotation &R){
    vector<Point3D>::iterator itp;
    int itcount=0;
    while(itcount<P.size())
    {
        itp=P.begin();
        Point3D res;
        res.x = (itp->x)*R.x1 + (itp->y)*R.y1 + (itp->z)*R.z1;
        res.y = (itp->x)*R.x2 + (itp->y)*R.y2 + (itp->z)*R.z2;
        res.z = (itp->x)*R.x3 + (itp->y)*R.y3 + (itp->z)*R.z3;
        P.erase(itp);
        P.push_back(res);
        itcount++;
    }
}

int main()
{
    FILE *fp=NULL;

	float f;

	Point3D xyz,T;

	Rotation R;

	vector<Point3D> P,X,Q;
	int ncount=0;



	char* F_PATH = "C:\\Users\\LIJ\\Desktop\\icp0\\data\\bunny.asc";
    //char* F_PATH = "C:\\Users\\LIJ\\Desktop\\icp0\\data\\1.txt";
	//int iLength ;

	//iLength = WideCharToMultiByte(CP_ACP, 0, argv[1], -1, NULL, 0, NULL, NULL);

	//WideCharToMultiByte(CP_ACP, 0, argv[1], -1, F_PATH, iLength, NULL, NULL);

	fp=fopen(F_PATH,"r");

	while(fscanf(fp,"%f",&f)!=EOF){

		xyz.x=f;//printf("%f",f);

		fscanf(fp,"%f",&f);

		xyz.y=f;//printf("%f",f);

		fscanf(fp,"%f",&f);

		xyz.z=f;//printf("%f   ",f);

		P.push_back(xyz); //ncount++;cout<<ncount<<endl;

	}

	if (fclose(fp)) {printf("error, file fails to close.\n");}

	fp=NULL;

	F_PATH = "C:\\Users\\LIJ\\Desktop\\icp0\\data\\bunny_perturbed.asc";
    //F_PATH = "C:\\Users\\LIJ\\Desktop\\icp0\\data\\2.txt";
	fp=fopen(F_PATH,"r");

	while(fscanf(fp,"%f",&f)!=EOF){

		xyz.x=f;
		fscanf(fp,"%f",&f);

		xyz.y=f;

		fscanf(fp,"%f",&f);

		xyz.z=f;

		Q.push_back(xyz);

	}

	if (fclose(fp)) {printf("error, file fails to close.\n");}

	fp=NULL;

    float E=1;
    Point3D P_mean, Q_mean,X_mean;
    //initialzing P,Q
    if(P.size() > Q.size()){P.swap(Q);}//make sure that P always has less elements
    //CalculateMeanPoint3D(P, P_mean);
    //CalculateMeanPoint3D(Q, Q_mean);
    //T.x = Q_mean.x - P_mean.x;
    //T.y = Q_mean.y - P_mean.y;
    //T.z = Q_mean.z - P_mean.z;
    //MovingPointSet(P, T);

    //begin iterating
    int itcount=0;
    while(E>0.0001 && itcount<10){
        FindCorrespondingPoint(P, Q, X, E);
        CalculateMeanPoint3D(P, P_mean);
        CalculateMeanPoint3D(X, X_mean);
        CalculateRotation(P, X, R, P_mean, X_mean);
        Rotate(P, R);

        Q.insert(Q.end(),X.begin(),X.end());
        cout<<"Q: "<<endl; print(Q);
        CalculateMeanPoint3D(Q, Q_mean);
        T.x = Q_mean.x - P_mean.x;
        T.y = Q_mean.y - P_mean.y;
        T.z = Q_mean.z - P_mean.z;
        MovingPointSet(P, T);cout<<"P: "<<endl; print(P);
        X.clear();
        itcount++;
        cout<<itcount<<" "<<E<<endl;
        cout<<"T: "<<endl;
        cout<<T.x<<" "<<T.y<<" "<<T.z<<endl;
        cout<<"R: "<<endl;
        cout<<R.x1<<" "<<R.y1<<" "<<R.z1<<endl;
        cout<<R.x2<<" "<<R.y2<<" "<<R.z2<<endl;
        cout<<R.x3<<" "<<R.y3<<" "<<R.z3<<endl;

    }
    //cout<<E<<endl;

    return 0;
}
