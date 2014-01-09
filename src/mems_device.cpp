#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <vector>
#include <math.h>

using namespace std;

struct holder_atoms{
	
	int i,j,k;
	float q,x,y,z;

	holder_atoms(int I, int J, int K, float Q, float X, float Y, float Z){
		i=I; j=J; k=K; x=X; y=Y; z=Z; q=Q;
	}

	friend ostream& operator<< (ostream& stream, const holder_atoms& in){

		stream << in.i << "\t " << in.j << "\t " << in.k << "\t " << in.q << "\t " \
			<< in.x << "\t " << in.y << "\t " << in.z << endl; 
		return stream;

	}

#if 0
	friend ostream& operator<< (ostream& stream, const holder_atoms& in){

		stream << in.x << "\t " << in.y << "\t " << in.z << endl; 
		return stream;

	}

#endif	

};
	
struct holder_bonds{
	
	int i,j,k,l;

	holder_bonds(int I, int J, int K, int L){
		i=I; j=J; k=K; l=L;
	}

	friend ostream& operator<< (ostream& stream, const holder_bonds& in){

		stream << in.i << "\t " << in.j << "\t " << in.k << "\t " << in.l << endl;
		return stream;

	}
};
		
struct holder_angles{
	
	int i,j,k,l,m;

	holder_angles(int I, int J, int K, int L, int M){
		i=I; j=J; k=K; l=L, m=M;
	}

	friend ostream& operator<< (ostream& stream, const holder_angles& in){

		stream << in.i << "\t " << in.j << "\t " << in.k << "\t " << in.l << "\t " << in.m << endl;
		return stream;

	}
};
	
int main(int argc, char *argv[]){

	   if (argc != 9){
                cerr << "USAGE : mems_device.x <cube side n> <lattice parameter> <limit in z> <outer limit in r=x^2+y^2, metal wall> <inner limit in r=x^2+y^2, metal wall> <helix pitch b (radius a == inner limit> <helix width w (along radial direction)> <helix thickness (fraction lattice parameter, down along z dirn)>\n"<< endl;
                exit(1);
        }



	int n 		= atoi(argv[1]);
	float a 	= atof(argv[2]);
	float z 	= atof(argv[3]);
	float ro 	= atof(argv[4]);
	float ri 	= atof(argv[5]);
	float b 	= atof(argv[6]);
	float w 	= atof(argv[7]);
	float th 	= atof(argv[8]);



	float hbond = 0.942;
	float angle = 104.52/180.0 * 3.14159265359;
	float m_charge = 0.0;

	float * vecA = new float[3*n*n*n];
	float * vecB = new float[3*n*(n-1)*(n-1)];
	float * vecC = new float[3*n*(n-1)*(n-1)];
	float * vecD = new float[3*n*(n-1)*(n-1)];


	//put x-y origin in middle; offset n/2
	float off = n/2.0;

	//cube corners
	for (int i=0; i<n; i++)
		for (int j=0; j<n; j++)
			for (int k=0; k<n; k++){

				int index = k+j*(n-1) + i*(n-1)*(n-1);
				vecA[3*index+0]=((float) i - off)*a;
				vecA[3*index+1]=((float) j - off)*a;
				vecA[3*index+2]=((float) k - 0.0)*a;
			}

	//top and bottom faces
	for (int i=0; i<n-1; i++)
		for (int j=0; j<n-1; j++)
			for (int k=0; k<n; k++){

				int index = k+j*(n-1) + i*(n-1)*(n-1);
				vecB[3*index+0]=((float) i + 0.5 - off)*a;
				vecB[3*index+1]=((float) j + 0.5 - off)*a;
				vecB[3*index+2]=((float) k - 0.0)*a;
			}

	//sides in x plane
	for (int i=0; i<n-1; i++)
		for (int j=0; j<n; j++)
			for (int k=0; k<n-1; k++){

				int index = k+j*(n-1) + i*(n-1)*(n-1);
				vecC[3*index+0]=((float) i + 0.5 - off)*a;
				vecC[3*index+1]=((float) j - off)*a;
				vecC[3*index+2]=((float) k + 0.5 - 0.0)*a;
			}

	//sides in y plane
	for (int i=0; i<n; i++)
		for (int j=0; j<n-1; j++)
			for (int k=0; k<n-1; k++){

				int index = k+j*(n-1) + i*(n-1)*(n-1);
				vecD[3*index+0]=((float) i - off)*a;
				vecD[3*index+1]=((float) j + 0.5 - off)*a;
				vecD[3*index+2]=((float) k + 0.5 - 0.0)*a;
			}


	//filter for interior water
	vector<holder_atoms> out;

	float minx=1e24,miny=1e24,minz=1e24,maxx=-1e24,maxy=-1e24,maxz=-1e24;

	int k=1,l=1;	
	for (int i=0; i< n*n*n; i++){

		// cylinder region
		float rad = (sqrt(vecA[3*i+0]*vecA[3*i+0] + vecA[3*i+1]*vecA[3*i+1]));
		bool val = (rad < ri);

		if (val){

			//solve for parameter t
			float t 	= vecA[3*i+2] / b;

			// angle center of helix section
			float tt = t - (float) ((int) (t/(2*M_PI))) * 2*M_PI;

			//angle extents
			float a_min = tt - b/(2*a*th);
			float a_max = tt + b/(2*a*th);


			//the angle of the atom in question
			float pt_angle = atan2(vecA[3*i+1], vecA[3*i+0]);

			pt_angle = (pt_angle < 0) ? 2*M_PI + pt_angle : pt_angle;

			bool val2 = ((pt_angle > a_min) && (pt_angle < a_max));
			val2 &= (rad > (ri-w));


			std::cerr << a_min << " " << a_max << " " << pt_angle << " " << val2 << endl;

			if (!val2){

				float h1z = vecA[3*i+2] + hbond;	
				float h2y = vecA[3*i+1] + cos(angle)*hbond;
				float h2z = vecA[3*i+2] - sin(angle)*hbond;

				out.push_back(holder_atoms(k,l,2,-0.834,vecA[3*i+0],vecA[3*i+1],vecA[3*i+2]));
				k++;
				out.push_back(holder_atoms(k,l,1,0.417,vecA[3*i+0],vecA[3*i+1],h1z));
				k++;
				out.push_back(holder_atoms(k,l,1,0.417,vecA[3*i+0],h2y,h2z));
				k++;
				l++;
			}

		}
	}

	for (int i=0; i< n*(n-1)*(n-1); i++){

		// cylinder region
		float rad = (sqrt(vecB[3*i+0]*vecB[3*i+0] + vecB[3*i+1]*vecB[3*i+1]));
		bool val = (rad < ri);

		if (val){

			//solve for parameter t
			float t 	= vecA[3*i+2] / b;

			// angle center of helix section
			float tt = t - (float) ((int) (t/(2*M_PI))) * 2*M_PI;

			//angle extents
			float a_min = tt - b/(2*a*th);
			float a_max = tt + b/(2*a*th);


			//the angle of the atom in question
			float pt_angle = atan2(vecB[3*i+1], vecB[3*i+0]);
			pt_angle = (pt_angle < 0) ? 2*M_PI + pt_angle : pt_angle;


			bool val2 = ((pt_angle > a_min) && (pt_angle < a_max));
			val2 &= (rad > (ri-w));

			std::cerr << a_min << " " << a_max << " " << pt_angle << " " << val2 << endl;

			if (!val2){

				float h1z = vecB[3*i+2] + hbond;	
				float h2y = vecB[3*i+1] + cos(angle)*hbond;
				float h2z = vecB[3*i+2] - sin(angle)*hbond;

				out.push_back(holder_atoms(k,l,2,-0.834,vecB[3*i+0],vecB[3*i+1],vecB[3*i+2]));
				k++;
				out.push_back(holder_atoms(k,l,1,0.417,vecB[3*i+0],vecB[3*i+1],h1z));
				k++;
				out.push_back(holder_atoms(k,l,1,0.417,vecB[3*i+0],h2y,h2z));
				k++;
				l++;
			}


		}
	}


	for (int i=0; i< n*(n-1)*(n-1); i++){
		// cylinder region
		float rad = (sqrt(vecC[3*i+0]*vecC[3*i+0] + vecC[3*i+1]*vecC[3*i+1]));
		bool val = (rad < ri);

		if (val){

			//solve for parameter t
			float t 	= vecA[3*i+2] / b;

			// angle center of helix section
			float tt = t - (float) ((int) (t/(2*M_PI))) * 2*M_PI;

			//angle extents
			float a_min = tt - b/(2*a*th);
			float a_max = tt + b/(2*a*th);


			//the angle of the atom in question
			float pt_angle = atan2(vecC[3*i+1], vecC[3*i+0]);

			pt_angle = (pt_angle < 0) ? 2*M_PI + pt_angle : pt_angle;



			bool val2 = ((pt_angle > a_min) && (pt_angle < a_max));
			val2 &= (rad > (ri-w));


			std::cerr << a_min << " " << a_max << " " << pt_angle << " " << val2 << endl;


			if (!val2){

				float h1z = vecC[3*i+2] + hbond;	
				float h2y = vecC[3*i+1] + cos(angle)*hbond;
				float h2z = vecC[3*i+2] - sin(angle)*hbond;

				out.push_back(holder_atoms(k,l,2,-0.834,vecC[3*i+0],vecC[3*i+1],vecC[3*i+2]));
				k++;
				out.push_back(holder_atoms(k,l,1,0.417,vecC[3*i+0],vecC[3*i+1],h1z));
				k++;
				out.push_back(holder_atoms(k,l,1,0.417,vecC[3*i+0],h2y,h2z));
				k++;
				l++;
			}

		}
	}

	for (int i=0; i< n*(n-1)*(n-1); i++){
		// cylinder region
		float rad = (sqrt(vecD[3*i+0]*vecD[3*i+0] + vecD[3*i+1]*vecD[3*i+1]));
		bool val = (rad < ri);

		if (val){

			//solve for parameter t
			float t 	= vecA[3*i+2] / b;

			// angle center of helix section
			float tt = t - (float) ((int) (t/(2*M_PI))) * 2*M_PI;

			//angle extents
			float a_min = tt - b/(2*a*th);
			float a_max = tt + b/(2*a*th);


			//the angle of the atom in question
			float pt_angle = atan2(vecD[3*i+1], vecD[3*i+0]);

			pt_angle = (pt_angle < 0) ? 2*M_PI + pt_angle : pt_angle;




			bool val2 = ((pt_angle > a_min) && (pt_angle < a_max));
			val2 &= (rad > (ri-w));

			std::cerr << a_min << " " << a_max << " " << pt_angle << " " << val2 << endl;


			if (!val2){

				float h1z = vecD[3*i+2] + hbond;	
				float h2y = vecD[3*i+1] + cos(angle)*hbond;
				float h2z = vecD[3*i+2] - sin(angle)*hbond;

				out.push_back(holder_atoms(k,l,2,-0.834,vecD[3*i+0],vecD[3*i+1],vecD[3*i+2]));
				k++;
				out.push_back(holder_atoms(k,l,1,0.417,vecD[3*i+0],vecD[3*i+1],h1z));
				k++;
				out.push_back(holder_atoms(k,l,1,0.417,vecD[3*i+0],h2y,h2z));
				k++;
				l++;
			}

		}
	}


	vector<holder_bonds> out2;
	//bonds
	for (int i=0; i<out.size()/3; i++){

		out2.push_back(holder_bonds(2*i+1,1,3*i+1,3*i+2));
		out2.push_back(holder_bonds(2*i+2,1,3*i+1,3*i+3));

	}

	vector<holder_angles> out3;
	//angles
	for (int i=0; i<out.size()/3; i++)
		out3.push_back(holder_angles(i+1,1,3*i+2,3*i+1,3*i+3));


	//filter for metal container

	vector<holder_atoms> out4;

	for (int i=0; i< n*n*n; i++){

		// cylinder region
		float rad = (sqrt(vecA[3*i+0]*vecA[3*i+0] + vecA[3*i+1]*vecA[3*i+1]));
		bool val = (rad > ri) && (rad < ro);


			//solve for parameter t
			float t 	= vecA[3*i+2] / b;

			// angle center of helix section
			float tt = t - (float) ((int) (t/(2*M_PI))) * 2*M_PI;

			//angle extents
			float a_min = tt - b/(2*a*th);
			float a_max = tt + b/(2*a*th);


			//the angle of the atom in question
			float pt_angle = atan2(vecA[3*i+1], vecA[3*i+0]);

			pt_angle = (pt_angle < 0) ? 2*M_PI + pt_angle : pt_angle;

			bool val2 = ((pt_angle > a_min) && (pt_angle < a_max));
			val2 &= ((rad > (ri-w)) && (rad < ri));


			std::cerr << a_min << " " << a_max << " " << pt_angle << " " << val2 << endl;

			if (val || val2){

				float h1z = vecA[3*i+2] + hbond;	
				float h2y = vecA[3*i+1] + cos(angle)*hbond;
				float h2z = vecA[3*i+2] - sin(angle)*hbond;

				out4.push_back(holder_atoms(k,l,3,m_charge,vecA[3*i+0],vecA[3*i+1],vecA[3*i+2]));
				k++;
				l++;
			}

	}

	for (int i=0; i< n*(n-1)*(n-1); i++){

		// cylinder region
		float rad = (sqrt(vecB[3*i+0]*vecB[3*i+0] + vecB[3*i+1]*vecB[3*i+1]));
		bool val = (rad > ri) && (rad < ro);


			//solve for parameter t
			float t 	= vecA[3*i+2] / b;

			// angle center of helix section
			float tt = t - (float) ((int) (t/(2*M_PI))) * 2*M_PI;

			//angle extents
			float a_min = tt - b/(2*a*th);
			float a_max = tt + b/(2*a*th);


			//the angle of the atom in question
			float pt_angle = atan2(vecB[3*i+1], vecB[3*i+0]);
			pt_angle = (pt_angle < 0) ? 2*M_PI + pt_angle : pt_angle;


			bool val2 = ((pt_angle > a_min) && (pt_angle < a_max));
			val2 &= ((rad > (ri-w)) && (rad < ri));

			std::cerr << a_min << " " << a_max << " " << pt_angle << " " << val2 << endl;

			if (val || val2){

				float h1z = vecB[3*i+2] + hbond;	
				float h2y = vecB[3*i+1] + cos(angle)*hbond;
				float h2z = vecB[3*i+2] - sin(angle)*hbond;

				out4.push_back(holder_atoms(k,l,3,m_charge,vecB[3*i+0],vecB[3*i+1],vecB[3*i+2]));
				k++;
				l++;
			}


	}


	for (int i=0; i< n*(n-1)*(n-1); i++){
		// cylinder region
		float rad = (sqrt(vecC[3*i+0]*vecC[3*i+0] + vecC[3*i+1]*vecC[3*i+1]));
		bool val = (rad > ri) && (rad < ro);


			//solve for parameter t
			float t 	= vecA[3*i+2] / b;

			// angle center of helix section
			float tt = t - (float) ((int) (t/(2*M_PI))) * 2*M_PI;

			//angle extents
			float a_min = tt - b/(2*a*th);
			float a_max = tt + b/(2*a*th);


			//the angle of the atom in question
			float pt_angle = atan2(vecC[3*i+1], vecC[3*i+0]);

			pt_angle = (pt_angle < 0) ? 2*M_PI + pt_angle : pt_angle;



			bool val2 = ((pt_angle > a_min) && (pt_angle < a_max));
			val2 &= ((rad > (ri-w)) && (rad < ri));


			std::cerr << a_min << " " << a_max << " " << pt_angle << " " << val2 << endl;


			if (val || val2){

				float h1z = vecC[3*i+2] + hbond;	
				float h2y = vecC[3*i+1] + cos(angle)*hbond;
				float h2z = vecC[3*i+2] - sin(angle)*hbond;

				out4.push_back(holder_atoms(k,l,3,m_charge,vecC[3*i+0],vecC[3*i+1],vecC[3*i+2]));
				k++;
				l++;
			}

	}

	for (int i=0; i< n*(n-1)*(n-1); i++){
		// cylinder region
		float rad = (sqrt(vecD[3*i+0]*vecD[3*i+0] + vecD[3*i+1]*vecD[3*i+1]));
		bool val = (rad > ri) && (rad < ro);


			//solve for parameter t
			float t 	= vecA[3*i+2] / b;

			// angle center of helix section
			float tt = t - (float) ((int) (t/(2*M_PI))) * 2*M_PI;

			//angle extents
			float a_min = tt - b/(2*a*th);
			float a_max = tt + b/(2*a*th);


			//the angle of the atom in question
			float pt_angle = atan2(vecD[3*i+1], vecD[3*i+0]);
			pt_angle = (pt_angle < 0) ? 2*M_PI + pt_angle : pt_angle;




			bool val2 = ((pt_angle > a_min) && (pt_angle < a_max));
			val2 &= ((rad > (ri-w)) && (rad < ri));

			std::cerr << a_min << " " << a_max << " " << pt_angle << " " << val2 << endl;


			if (val2 || val){

				float h1z = vecD[3*i+2] + hbond;	
				float h2y = vecD[3*i+1] + cos(angle)*hbond;
				float h2z = vecD[3*i+2] - sin(angle)*hbond;

				out4.push_back(holder_atoms(k,l,3,m_charge,vecD[3*i+0],vecD[3*i+1],vecD[3*i+2]));
				k++;
				l++;
			}

	}



	//for (std::vector<holder_atoms>::iterator it = out4.begin() ; it != out4.end(); ++it)
	//	cout << *it ;

	for (int i=0; i<out.size(); i++){
	
		if (out[i].x < minx)
			minx=out[i].x;
		if (out[i].x > maxx)
			maxx=out[i].x;
		if (out[i].y < miny)
			miny=out[i].y;
		if (out[i].y > maxy)
			maxy=out[i].y;
		if (out[i].z < minz)
			minz=out[i].z;
		if (out[i].z > maxz)
			maxz=out[i].z;

	}

	for (int i=0; i<out4.size(); i++){
	
		if (out4[i].x < minx)
			minx=out[i].x;
		if (out4[i].x > maxx)
			maxx=out[i].x;
		if (out4[i].y < miny)
			miny=out[i].y;
		if (out4[i].y > maxy)
			maxy=out[i].y;
		if (out4[i].z < minz)
			minz=out[i].z;
		if (out4[i].z > maxz)
			maxz=out[i].z;

	}


	cout << "#mems device atoms file" << endl; 
	cout << endl;
	cout << out.size()+out4.size() << " atoms" << endl; 
	cout << 2*out.size()/3 << " bonds" << endl; 
	cout << out.size()/3 << " angles" << endl; 
	cout << "0 dihedrals" << endl;
	cout << "0 impropers" << endl;
	cout << endl;	
	cout << "3 atom types" << endl;
	cout << "1 bond types" << endl;
	cout << "1 angle types" << endl;
	cout << "0 dihedral types" << endl;
	cout << "0 improper types" << endl;
	cout << endl;
	cout << minx << " " << maxx << " " << "xlo xhi" << endl;
	cout << miny << " " << maxy << " " << "ylo yhi" << endl;
	cout << minz << " " << maxz << " " << "zlo zhi" << endl;
	cout << endl;
	cout << "Masses" << endl;
	cout << endl;
	cout << "1 1.008" << endl;
	cout << "2 15.994" << endl;
	cout << "3 12.011" << endl;
	cout << endl;
	//cout << "Pair Coeffs" << endl;
	//cout << endl;
	//cout << "1 0.046 0.4000135 0.046 0.4000135" << endl;
	//cout << "2 0.1521 3.150574 0.1521 3.150574" << endl;
	//cout << endl;
	cout << "Atoms" << endl;
	cout << endl;
	for (std::vector<holder_atoms>::iterator it = out.begin() ; it != out.end(); ++it)
		cout << *it ;
	for (std::vector<holder_atoms>::iterator it = out4.begin() ; it != out4.end(); ++it)
		cout << *it ;
	cout << endl;
	//cout << "Bond Coeffs" << endl;
	//cout << endl;
	//cout << " 1 450 0.9572" << endl;
	//cout << endl;
	cout << "Bonds" << endl;
	cout << endl;
	for (std::vector<holder_bonds>::iterator it = out2.begin() ; it != out2.end(); ++it)
		cout << *it ;
	cout << endl;
	//cout << "Angle Coeffs" << endl;
	//cout << endl;
	//cout << "1 55 104.52 0 0 " << endl;
	cout << endl;
	cout << "Angles" << endl;
	cout << endl;
	for (std::vector<holder_angles>::iterator it = out3.begin() ; it != out3.end(); ++it)
		cout << *it ;
	
	

#if 0

	for (std::vector<holder_atoms>::iterator it = out.begin() ; it != out.end(); ++it)
		cout << *it ;

	for (std::vector<holder_bonds>::iterator it = out2.begin() ; it != out2.end(); ++it)
		cout << *it ;

	for (std::vector<holder_angles>::iterator it = out3.begin() ; it != out3.end(); ++it)
		cout << *it ;
#endif


	delete vecA;
	delete vecB;
	delete vecC;
	delete vecD;

	return 0;

}
