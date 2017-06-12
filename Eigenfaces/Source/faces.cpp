
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: faces.cpp                                                                            //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////



#include "stdafx.h"

#include "jacob.h"

Faces::Faces()
:
Array<Face>(),
width(0),
height(0),
vector_size(0)
{
	//empty
}

Faces::Faces(int count, int width, int height)
:
Array<Face>(count),
width(width),
height(height),
vector_size(width*height)
{
	for (int i=0; i<getSize(); i++) {
		(*this)[i].resize(width, height, 1);
        //printf("Face %d\n", (*this)[i]);
        //printf("width %d\n", width);
        //printf("height %d\n", height);
        
	}
    //printf("count %d\n", count);
}

void Faces::load(BinaryFileReader& file)
{
	resize(file.readInt());
	width=file.readInt();
	height=file.readInt();
    //printf("width %d\n", width);
    //printf("height %d\n", height);
	vector_size=width*height;
	for (int i=0; i<getSize(); i++) {
		(*this)[i].load(file);
	}
	average_face.load(file);
	//std::cout << "Loaded faces from '" << file.getFilename() << "'" << std::endl;

}

void Faces::load(std::string filename)
{
	BinaryFileReader file(filename);
	load(file);
}

void Faces::save(std::string filename) const
{
	BinaryFileWriter file(filename);
	save(file);
}

void Faces::save(BinaryFileWriter& file) const
{
	file.write(getSize());
	file.write(width);
	file.write(height);
	for (int i=0; i<getSize(); i++) {
		(*this)[i].save(file);
	}
	average_face.save(file);
	std::cout << "Saved faces to '" << file.getFilename() << "'" << std::endl;
}

void Faces::output(std::string filepattern) const
{
	for (int i=0; i<getSize(); i++) {
		// normalize for output
		Image out_image;
		(*this)[i].normalize(0.0, 255.0, out_image);
		std::string filename=Functions::filenameNumber(filepattern, i, getSize()-1);
		out_image.saveTarga(filename);
	}
}

void Faces::eigenFaces(EigFaces& results, int n) const
{
	// size the results vector
	results.resize(n);
	results.setHeight(height);
	results.setWidth(width);
    
	// allocate matrices
	double **matrix = Jacobi::matrix(1, vector_size, 1, vector_size);
	double **eigmatrix = Jacobi::matrix(1, vector_size, 1, vector_size);
	double *eigenvec = Jacobi::vector(1, vector_size);
        
	// --------- TODO #1: fill in your code to prepare a matrix whose eigenvalues and eigenvectors are to be computed.
	// Also be sure you store the average face in results.average_face (A "set" method is provided for this).
    Face grayFace;
    std::vector<Vector> faces;
    Face avgFace (width, height);
    
    for(int x = 0; x < width; x++) {
        for(int y=0; y < height; y++){
            double sum=0;
            double avg=0;
            for (int i=0; i<getSize(); i++) {
                grayFace = (*this)[i];
                //printf("face[%d].pixel[%d, %d] = %f\n", i, x, y, (*this)[i].pixel(x,y,0));
                sum+= grayFace.pixel(x,y,0);
            }
            avgFace.pixel(x,y,0) = sum/getSize();
        }
    }
    results.setAverage(avgFace);
    
    //put all images into n x 1 vectors (n = height*width)
    for (int i=0; i<getSize(); i++) {
        Vector face (vector_size);
        Vector faceMinusAvg (vector_size);
        //printf("i %d\n", i);
        (*this)[i].grayscale(grayFace);
        
        //int width1 = grayFace.getWidth();
        //printf("width %d\n", width1);   output 25
        int count = 0; //used to convert image into one column of pixels
        double sum = 0;
        double average = 0;
        for(int x = 0; x < width; x++) {
            for(int y=0; y < height; y++){
                face[count] = grayFace.pixel(y,x,0);
                //printf("pixel %f\n", grayFace.pixel(x,y,0));
                //printf("face[%d]:  %f\n", count, face[count]);
                sum += face[count];
                count++;
            }
        }
        faces.push_back(face);
        //printf("\n\n\n\n\n\n\n\n\n\n\n");
    }
    
    Vector sumFace (vector_size);
    Vector avgFaceV (vector_size);
    std::vector<Vector> normFaces;
    
    //sum faces
    for(int j=1; j< faces.size(); j++){
        faces[0].add(faces[j], sumFace);
    }
    //avg faces
    for(int i=0; i<vector_size; i++){
        avgFaceV[i] = sumFace[i]/getSize();
    }

    //normalize faces
    for(int i=0; i<faces.size(); i++){
        Vector v (vector_size);
        faces[i].sub(avgFaceV, v);
        normFaces.push_back(v);
    }
    
    //initialize to zero
    for(int i=1; i<=vector_size; i++){
        for(int j=1; j<=vector_size; j++){
            matrix[i][j] = 0;
        }
    }
    
    //get covariance matrix
    for(int i=0; i<getSize(); i++){
        for(int j=0; j< vector_size; j++){
            for(int k=0; k< vector_size; k++){
                //printf("\n hi again \n");
                //printf("normFaces[i][j] = %f\n", normFaces[i][j]);
                //printf("normFaces[i][k] = %f\n", normFaces[i][k]);
                matrix[j+1][k+1] += normFaces[i][k]*normFaces[i][j];
                //printf("matrix[%d][%d] =  %f\n", j, k, matrix[j][k]);
            }
        }
    }
    
    //Computes all eigenvalues and eigenvectors of a real symmetric matrix a[1..n][1..n]. On output, elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a. v[1..n][1..n] is a matrix whose columns contain, on output, the normalized eigenvectors of a. nrot returns the number of Jacobi rotations that were required.
    //void jacobi(double **a, int n, double d[], double **v, int *nrot);


	// find eigenvectors
	int nrot;
	Jacobi::jacobi(matrix, vector_size, eigenvec, eigmatrix, &nrot);
	// sort eigenvectors
	Array<int> ordering;
	sortEigenvalues(eigenvec, ordering);
	for (int i=0; i<n; i++) {
		for (int k=0; k<vector_size; k++) {
			results[i][k] = eigmatrix[k+1][ordering[i]+1];
		}
	}
	// free matrices
	Jacobi::free_matrix(matrix, 1, vector_size, 1, vector_size);
	Jacobi::free_matrix(eigmatrix, 1, vector_size, 1, vector_size);
	Jacobi::free_vector(eigenvec, 1, vector_size);
}



int Faces::getWidth() const
{
	return width;
}

int Faces::getHeight() const
{
	return height;
}

void Faces::setWidth(int width)
{
	width=width;
	vector_size=width*height;
}

void Faces::setHeight(int height)
{
	height=height;
	vector_size=width*height;
}

void Faces::sortEigenvalues(double *eigenvec, Array<int>& ordering) const
{
	// for now use simple bubble sort
	ordering.resize(vector_size);
	std::list<EigenVectorIndex> list;
	for (int i=0; i<vector_size; i++) {
		EigenVectorIndex e;
		e.eigenvalue=eigenvec[i+1];
		e.index=i;
		list.push_back(e);
	}
	bool change=true;
	list.sort();
	std::list<EigenVectorIndex>::iterator it=list.begin();
	int n=0;
	while (it!=list.end()) {
		ordering[n] = (*it).index;
		it++;
		n++;
	}
}

