
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: eigfaces.cpp                                                                         //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////



#include "stdafx.h"
static bool overlap(FacePosition, double, FacePosition, double, int, int);
static void mark(Image &, int, int, int, int, double);



EigFaces::EigFaces()
:
Faces()
{
	//empty
}

EigFaces::EigFaces(int count, int width, int height)
:
Faces(count, width, height)
{
	//empty
}

void EigFaces::projectFace(const Face& face, Vector& coefficients) const
{
	if (face.getWidth()!=width || face.getHeight()!=height) {
		throw Error("Project: Face to project has different dimensions");
	}

	coefficients.resize(getSize());
	// ----------- TODO #2: compute the coefficients for the face and store in coefficients.
    
    Face norm;
    face.sub(getAverage(), norm);
    
    for (int i=0; i<getSize(); i++) {
        Face eigFace = (*this)[i];
        //eigFace.normalize(0.0, 255.0, eigFace);
        coefficients[i] = norm.dot(eigFace);
    }
    
    for (int i=0; i<getSize(); i++) {
        //coefficients[i] /= coefficients.mag();
        //printf("coefficients[%d] = %f\n", i, coefficients[i]);
    }

    
}

void EigFaces::constructFace(const Vector& coefficients, Face& result) const
{	
	// ----------- TODO #3: construct a face given the coefficients
    //Vector faceV (width*height);
    //Vector faceAvg (width*height);
    Face f;

    result.resize(width, height,1);
    
    f = getAverage();
    //f = (*this)[0];
    
    for(int i=0; i< getSize(); i++){
        Face eigFace = (*this)[i];
        eigFace *= coefficients[i];
       /*
        for(int x = 0; x < width; x++) {
            for(int y=0; y< height; y++){
                eigFace.pixel(x,y,0)*= coefficients[i];
                //printf("result.pixel[%d, %d] = %f\n", y, x, result.pixel(y,x,0));
                //printf("face[%d].pixel[%d, %d] = %f\n", i, x, y, (*this)[i].pixel(x,y,0)*coefficients[i]);
            }
        }
        */
        f.add(eigFace, f);
    }
    
    result = f;
    
    /*
    for(int x = 0; x < width; x++) {
        for(int y=0; y < height; y++){
            //printf("result.pixel[%d, %d] = %f\n", x, y, result.pixel(x,y,0));
            result.pixel(x,y,0) += getAverage().pixel(x,y,0);
            //printf("result.pixel[%d, %d] = %f\n", x, y, result.pixel(x,y,0));
        }
    }
     

    result.normalize(0.0, 255.0, result);
    
    //test shows (*this) is an eigenface but needs to be normalized
    //(*this)[1].normalize(0.0, 255.0, result);
    //result = (*this)[0];
    */

}

bool EigFaces::isFace(const Face& face, double max_reconstructed_mse, double& mse) const
{
	// ----------- TODO #4: Determine if an image is a face and return true if it is. Return the actual
	// MSE you calculated for the determination in mse
	// Be sure to test this method out with some face images and some non face images
	// to verify it is working correctly.
    Vector coefficients;
    Face projFace = face;
    projectFace(face, coefficients);
    constructFace(coefficients, projFace);
    mse = projFace.mse(face);
    
    if (mse < max_reconstructed_mse) {
        return true;
    }
    else {
        return false;
    }
    
}

bool EigFaces::verifyFace(const Face& face, const Vector& user_coefficients, double max_coefficients_mse, double& mse) const
{
	// ----------- TODO #5 : Determine if face is the same user give the user's coefficients.
	// return the MSE you calculated for the determination in mse. (need to return a boolean though?)
    Vector coefficients;
    Face projFace = face;
    projectFace(projFace, coefficients);
    mse = coefficients.mse(user_coefficients);
    //printf("mse = %f\n", mse);
    
    if (mse > max_coefficients_mse) {
        return false;
    }
    else {
        return true;
    }
}

void EigFaces::recognizeFace(const Face& face, Users& users) const
{
	// ----------- TODO #6: Sort the users by closeness of match to the face
    
    Vector coefficients;
    Face projFace = face;
    projectFace(projFace,coefficients);
    
    for (int i = 0; i < users.getSize(); i++) {
        double mse = users[i].mse(coefficients);
        //printf("mse = %f\n", mse);
        users[i].setMse(mse);
    }
    users.sort();
    printf("\n\n");

}

void EigFaces::findFace(const Image& img, double min_scale, double max_scale, double step, int n, bool crop, Image& result) const
{
	// ----------- TODO #7: Find the faces in Image. Search image scales from min_scale to max_scale inclusive,
	// stepping by step in between. Find the best n faces that do not overlap each other. If crop is true,
	// n is one and you should return the cropped original img in result. The result must be identical
	// to the original besides being cropped. It cannot be scaled and it must be full color. If crop is
	// false, draw green boxes (use r=100, g=255, b=100) around the n faces found. The result must be
	// identical to the original image except for the addition of the boxes.
   
    
    std::list<FacePosition> bestMatches;

    // iterate through scales first as per the hint
    for (double s = min_scale; s <= max_scale; s += step)
    {
        //scale and resample image, s is scale, sImage is scaled image
        Image sImage(s*img.getWidth(), s*img.getHeight(), img.getColors());
        img.resample(sImage);
        
        // iterate through windows in image
        for (int x = 0; x < sImage.getWidth(); x++)
        {
            for (int y = 0; y < sImage.getHeight(); y++)
            {
                double mse;
                
                // create subimage, or window, and find the mse of window & its projection on face space
                Face window((*this).width, (*this).height);
                window.subimage(x, x+(*this).width-1, y, y+(*this).height-1, sImage, false);
                isFace(window, 200, mse);
                
                //FacePosition struct of the window
                FacePosition cur_Face;
                cur_Face.x = x, cur_Face.y = y, cur_Face.scale = s;
                cur_Face.error = mse * window.mse(getAverage())/window.var(); //suggested formula from project description
                
                FacePosition fpWorst; //worst match out of n bestMatches
                fpWorst.error = 0; //dummy value to tell if has a real face assignment
                bool moveOn = false; // will set to true if the current face overlap a better face
                std::list<FacePosition> temp(0); //to temporarily store the worst match for easier manipulation of bestMatches
                
                for (std::list<FacePosition>::iterator iterator = bestMatches.begin(), end = bestMatches.end(); iterator != end; ++iterator)
                {
                    FacePosition p = *iterator;
                    //if current face is a worse match than face in bestMatch
                    if (cur_Face.error > p.error){
                        //move on from this face if it's an overlap, continue with another face
                        if (overlap(cur_Face, cur_Face.scale, p, p.scale, (*this).width, (*this).height)){
                            moveOn = true;
                            break;
                        }
                    }
                    //if current face is better match than face in bestMatch
                    else{
                        if (overlap(cur_Face, cur_Face.scale, p, p.scale, (*this).width, (*this).height)){
                            temp.push_back(p);
                        }
                        else{
                            //if no overlap, its the worst fp
                            if (fpWorst.error < p.error)
                                fpWorst = p;
                        }
                    }
                }
                
                if (moveOn == true) continue;
                //remove the overlapping faces if they're not better than the current face
                for (std::list<FacePosition>::iterator iterator = temp.begin(); iterator != temp.end(); iterator++){
                    FacePosition p = (*iterator);
                    bestMatches.remove(p);
                }
                //if the current face is better than the worst in bestmatches
                if (bestMatches.size() == n && fpWorst.error != 0){
                    bestMatches.remove(fpWorst);
                }
                //if the list size < n, add the current face
                if (bestMatches.size() < n){
                    bestMatches.push_back(cur_Face);
                }
            }
        }
    }
    
    //crop
    if (crop == false)
    {
        result.resize(img.getWidth(), img.getHeight(), img.getColors());
        for (int i = 0; i < img.getSize(); i++){
            result[i] = img[i];
        }
        
        //draw borders about faces
        for (std::list<FacePosition>::iterator iterator = bestMatches.begin(); iterator != bestMatches.end(); iterator++){
            FacePosition p = *iterator;
            mark(result, p.x, p.y, (*this).width, (*this).height, p.scale);
        }
    }
    
    //mark with borders
    else
    {
        FacePosition final = *(bestMatches.begin());
        double fs = final.scale;
        double fx = final.x, fy = final.y;
        result.resize(((*this).width/fs), ((*this).height/fs), img.getColors());
        for (int i=0; i<((*this).height/fs); i++){
            for (int j=0; j<((*this).width/fs); j++){
                for (int k=0; k<img.getColors(); k++){
                    result.pixel(j, i, k) = img.pixel((fx/fs) + j, (fy/fs) + i, k);
                }
            }
        }
    }


}


static bool overlap(FacePosition p1, double scale1, FacePosition p2, double scale2, int width, int height)
{
    // FacePosition 1 (leftmost, rightmost, top, and bottom coordinates)
    int l1 = (double) p1.x/scale1, r1 = (double) (p1.x+width-1)/scale1;
    int t1 = (double) p1.y/scale1, b1 = (double) (p1.y+height-1)/scale1;
    
    // FacePosition 2
    int l2 = (double) p2.x/scale2, r2 = (double) (p2.x+width-1)/scale2;
    int t2 = (double) p2.y/scale2, b2 = (double) (p2.y+height-1)/scale2;
    
    //if ((overlap on horizontal plane) && (overlap on vertical plane)), faces overlap!!!
    /*
    if (((l1<l2 && r1>l2) || (l1>l2, l1<r2)) && ((t1>t2 && t2>b1) || (t1<t2 && t1>b2))){
        return true;
    }
    else {
        return false;
    }
     */
    return (l1 < r2 && r1 > l2 && b1 > t2 && t1 < b2);
}


static void mark(Image &img, int x, int y, int width, int height, double scale)
{
    int l = x/scale, r = (x + width - 1)/scale;
    int t = y/scale, b = (y + height - 1)/scale;
    
    // draw top, bottom
    for (int i = l; i <= r; i++)
    {
        img.pixel(i, t, 0) = 100, img.pixel(i, t, 1) = 255, img.pixel(i, t, 2) = 100;
        img.pixel(i, b, 0) = 100, img.pixel(i, b, 1) = 255, img.pixel(i, b, 2) = 100;
    }
    
    // draw left, right
    for (int i = t; i <= b; i++)
    {
        img.pixel(l, i, 0) = 100, img.pixel(l, i, 1) = 255, img.pixel(l, i, 2) = 100;
        img.pixel(r, i, 0) = 100, img.pixel(r, i, 1) = 255, img.pixel(r, i, 2) = 100;
    }
}

void EigFaces::morphFaces(const Face& face1, const Face& face2, double distance, Face& result) const
{
	// TODO (extra credit): MORPH along *distance* fraction of the vector from face1 to face2 by
	// interpolating between the coefficients for the two faces and reconstructing the result.
	// For example, distance 0.0 will approximate the first, while distance 1.0 will approximate the second.
	// Negative distances are ok two.

}

const Face& EigFaces::getAverage() const
{
	return average_face;
}

void EigFaces::setAverage(const Face& average)
{
	average_face=average;
}



