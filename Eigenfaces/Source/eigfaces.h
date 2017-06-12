

/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: eigfaces.h                                                                              //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////

// for sorting best face positions
struct FacePosition
{
	FacePosition() : error(DBL_MAX) {}
	int x, y;
	double scale;
	double error;
    
    //std::list::remove uses the comparison operator 
    bool operator==(const FacePosition& p) const
    {
        //printf("== operator working\n");
        return (x == p.x && y == p.y && error == p.error && scale == p.scale);
    }
};

// for sorting matches
struct Match
{
	int index;
	double mse;
	bool operator<(const Match& rhs) { return mse<rhs.mse; }
};

class EigFaces : public Faces
{
public:
	EigFaces();
	EigFaces(int count, int width, int height);
	void projectFace(const Face& face, Vector& coefficients) const;
	void constructFace(const Vector& coefficients, Face& result) const;
	bool isFace(const Face& face, double max_reconstructed_mse, double& mse) const;
	bool verifyFace(const Face& face, const Vector& user_coefficients, double max_coefficients_mse, double& mse) const; 
	void recognizeFace(const Face& face, Users& users) const;
	void findFace(const Image& img, double min_scale, double max_scale, double step, int n, bool crop, Image& result) const;
	void morphFaces(const Face& face1, const Face& face2, double distance, Face& result) const;
	const Face& getAverage() const;
	void setAverage(const Face& average);
private:

};
