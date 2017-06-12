#include <assert.h>
#include <math.h>
#include <FL/Fl.H>
#include <FL/Fl_Image.H>
#include "features.h"
#include "ImageLib/FileIO.h"
#include "ImageLib/Convolve.h"

#define PI 3.14159265358979323846

//checks for potential bound errors in pixel operations
double imagePixel(int x, int y, CFloatImage img, int band)
{
    int h = img.Shape().height;
    int w = img.Shape().width;
    int numBands = img.Shape().nBands;
    
    if(x<0 || y<0 || x>=w || y>=h || band>=numBands)
        return 0;
    
    else
        return img.Pixel(x, y, band);
}

// Compute features of an image.
bool computeFeatures(CFloatImage &image, FeatureSet &features, int featureType) {
	// TODO: Instead of calling dummyComputeFeatures, write your own
	// feature computation routines and call them here.

 //   printf("hi computeFeatures");
	switch (featureType) {
	case 1:
		dummyComputeFeatures(image, features);
		break;
	case 2:
		ComputeHarrisFeatures(image, features);
		break;
	default:
		return false;
	}

    
	// This is just to make sure the IDs are assigned in order, because
	// the ID gets used to index into the feature array.
	for (unsigned int i=0; i<features.size(); i++) {
		features[i].id = i+1;
	}

	return true;
}

// Perform a query on the database.  This simply runs matchFeatures on
// each image in the database, and returns the feature set of the best
// matching image.
bool performQuery(const FeatureSet &f, const ImageDatabase &db, int &bestIndex, vector<FeatureMatch> &bestMatches, double &bestScore, int matchType) {
	// Here's a nice low number.
	bestScore = -1e100;

    //printf("hi performQuery");
	vector<FeatureMatch> tempMatches;
	double tempScore;

	for (unsigned int i=0; i<db.size(); i++) {
		if (!matchFeatures(f, db[i].features, tempMatches, tempScore, matchType)) {
			return false;
		}

		if (tempScore > bestScore) {
			bestIndex = i;
			bestScore = tempScore;
			bestMatches = tempMatches;
		}
	}

	return true;
}

// Match one feature set with another.
bool matchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore, int matchType) {
	// TODO: We have given you the ssd matching function, you must write your own
	// feature matching function for the ratio test.
	
	printf("\nMatching features.......\n");

	switch (matchType) {
	case 1:
		ssdMatchFeatures(f1, f2, matches, totalScore);
		return true;
	case 2:
		ratioMatchFeatures(f1, f2, matches, totalScore);
		return true;
	default:
		return false;
	}
}

// Evaluate a match using a ground truth homography.  This computes the
// average SSD distance between the matched feature points and
// the actual transformed positions.
double evaluateMatch(const FeatureSet &f1, const FeatureSet &f2, const vector<FeatureMatch> &matches, double h[9]) {
	double d = 0;
	int n = 0;

	double xNew;
	double yNew;

    unsigned int num_matches = matches.size();
	for (unsigned int i=0; i<num_matches; i++) {
		int id1 = matches[i].id1;
        int id2 = matches[i].id2;
        applyHomography(f1[id1-1].x, f1[id1-1].y, xNew, yNew, h);
		d += sqrt(pow(xNew-f2[id2-1].x,2)+pow(yNew-f2[id2-1].y,2));
		n++;
	}	

	return d / n;
}

void addRocData(const FeatureSet &f1, const FeatureSet &f2, const vector<FeatureMatch> &matches, double h[9],vector<bool> &isMatch,double threshold,double &maxD) {
	double d = 0;

	double xNew;
	double yNew;

    unsigned int num_matches = matches.size();
	for (unsigned int i=0; i<num_matches; i++) {
		int id1 = matches[i].id1;
        int id2 = matches[i].id2;
		applyHomography(f1[id1-1].x, f1[id1-1].y, xNew, yNew, h);

		// Ignore unmatched points.  There might be a better way to
		// handle this.
		d = sqrt(pow(xNew-f2[id2-1].x,2)+pow(yNew-f2[id2-1].y,2));
		if (d<=threshold)
		{
			isMatch.push_back(1);
		}
		else
		{
			isMatch.push_back(0);
		}

		if (matches[i].score>maxD)
			maxD=matches[i].score;
	}	
}

vector<ROCPoint> computeRocCurve(vector<FeatureMatch> &matches,vector<bool> &isMatch,vector<double> &thresholds)
{
	vector<ROCPoint> dataPoints;

	for (int i=0; i < (int)thresholds.size();i++)
	{
		//printf("Checking threshold: %lf.\r\n",thresholds[i]);
		int tp=0;
		int actualCorrect=0;
		int fp=0;
		int actualError=0;
		int total=0;

        int num_matches = (int) matches.size();
		for (int j=0;j < num_matches;j++)
		{
			if (isMatch[j])
			{
				actualCorrect++;
				if (matches[j].score<thresholds[i])
				{
					tp++;
				}
			}
			else
			{
				actualError++;
				if (matches[j].score<thresholds[i])
				{
					fp++;
				}
            }
			
			total++;
		}

		ROCPoint newPoint;
		//printf("newPoints: %lf,%lf",newPoint.trueRate,newPoint.falseRate);
		newPoint.trueRate=(double(tp)/actualCorrect);
		newPoint.falseRate=(double(fp)/actualError);
		//printf("newPoints: %lf,%lf",newPoint.trueRate,newPoint.falseRate);

		dataPoints.push_back(newPoint);
	}

	return dataPoints;
}


// Compute silly example features.  This doesn't do anything
// meaningful.
void dummyComputeFeatures(CFloatImage &image, FeatureSet &features) {
	CShape sh = image.Shape();
	Feature f;

	for (int y=0; y<sh.height; y++) {
		for (int x=0; x<sh.width; x++) {
			double r = image.Pixel(x,y,0);
			double g = image.Pixel(x,y,1);
			double b = image.Pixel(x,y,2);

			if ((int)(255*(r+g+b)+0.5) % 100  == 1) {
				// If the pixel satisfies this meaningless criterion,
				// make it a feature.
				
				f.type = 1;
				f.id += 1;
				f.x = x;
				f.y = y;

				f.data.resize(1);
				f.data[0] = r + g + b;

				features.push_back(f);
			}
		}
	}
}


double thresh = 1.0;
int recur = 0;

void ComputeHarrisFeatures(CFloatImage &image, FeatureSet &features)
{
	//Create grayscale image used for Harris detection
	CFloatImage grayImage=ConvertToGray(image);

	//Create image to store Harris values
	CFloatImage harrisImage(image.Shape().width,image.Shape().height,1);

	//Create image to store local maximum harris values as 1, other pixels 0
	CByteImage harrisMaxImage(image.Shape().width,image.Shape().height,1);
    
    CFloatImage orientationImage(image.Shape().width, image.Shape().height, 1);
	
	//compute Harris values puts harris values at each pixel position in harrisImage. 
	//You'll need to implement this function.
    computeHarrisValues(grayImage, harrisImage, orientationImage);
	
	// Threshold the harris image and compute local maxima.  You'll need to implement this function.

    computeLocalMaxima(harrisImage,harrisMaxImage, thresh);
    
    // Prints out the harris image for debugging purposes
	CByteImage tmp(harrisImage.Shape());
	convertToByteImage(harrisImage, tmp);
    WriteFile(tmp, "harris.tga");
    

	// TO DO--------------------------------------------------------------------
	//Loop through feature points in harrisMaxImage and create feature descriptor 
	//for each point above a threshold

    
    ///// Simple Feature Descriptor /////
    /*
     for (int y=0;y<harrisMaxImage.Shape().height;y++) {
         for (int x=0;x<harrisMaxImage.Shape().width;x++) {
     
             // Skip over non-maxima
             if (harrisMaxImage.Pixel(x, y, 0) == 0)
                 continue;
     
             Feature f;
             int id = 0;
     
             f.type = 2;
             f.id = id;
             f.x = x;
             f.y = y;
        
             f.data.resize(5 * 5);
             int range = 5 / 2;
             int index = 1;
     
             for (int row = -range; row <= range; row++ ){
                 for (int col = -range; col <= range; col++){
                     f.data[index]= imagePixel(x+col, y+row, grayImage, 0);
                     index++;
                 }
             }
     
             // Add the feature to the list of features
             features.push_back(f);
             id++;
     
        }
     }
     */
    
    
    ///// Circular Feature Descriptor /////
    
    int id = 0;
    int w = image.Shape().width;
    int h = image.Shape().height;
    
    for (int y=0;y<harrisMaxImage.Shape().height;y++) {
        for (int x=0;x<harrisMaxImage.Shape().width;x++) {
            
            // Skip over non-maxima
            if (harrisMaxImage.Pixel(x, y, 0) == 0)
                continue;
            
            Feature f;
            f.type = 2;
            f.id=id;
        
            f.x = x;
            f.y = y;
            f.angleRadians=orientationImage.Pixel(x,y,0);
            
            double angle = f.angleRadians;
            
            //printf("\n Angle radians: %f\n", angle);
            double stepAngle = 20.0/180.0 * PI;
        
            //feature descriptor is a vector wih (#steps, 18)*(#radius lengths, 3) = 54 elements
            f.data.resize(54);
            
            // r: radius of circle about feature
            double r=6.0;
            
            // 3 radius lengths: r, 2r, 3r
            for (int i=1; i<=3; i++)
            {
                //since (pi/9)*18= 2pi, 18 steps used per rotation
                for (int step=0; step<18; step++)
                {
                    // xRad, yRad are positions on the circumference of the scanning circle
                    int xr = x + floor(r * (double)i * cos(angle - step*stepAngle) + .5);
                    int yr = y + floor(r * (double)i * sin(angle - step*stepAngle) + .5);
                    if (!image.Shape().InBounds(xr, yr))
                    {
                        f.data[step + 18*(i-1)] = 0;
                        continue;
                    }
                    f.data[step + 18*(i-1)] = grayImage.Pixel(xr, yr, 0);
                }
            }

            
            // Normalize intensity for illuminosity invariance //
            
            // average
            double sum = 0.0;
            for (int i = 0; i < 54; i++)
                sum += f.data[i];
            double avg = sum/54.0;
            
            // standard deviation
            double sum2 = 0.0;
            for (int i = 0; i < 54; i++)
                sum2 += pow(f.data[i]-avg, 2);
            double std = sqrt(sum2/54);
            
            // fill in descriptor data
            for (int i = 0; i < 54; i++){
                f.data[i] = (f.data[i]-avg)/(std+.0000001);
            }
            
            features.push_back(f);
            id++;
       
        }
    }
    
 
}


//TO DO---------------------------------------------------------------------
//Loop through the image to compute the harris corner values as described in class
// srcImage:  grayscale of original image
// harrisImage:  populate the harris values per pixel in this image
void computeHarrisValues(CFloatImage &srcImage, CFloatImage &harrisImage, CFloatImage &orientationImage)
{
    
	int w = srcImage.Shape().width;
    int h = srcImage.Shape().height;
    
    /* Sobel filters */
    static float k_SobelX[9] = { -1, 0, 1,
        -2, 0, 2,
        -1, 0, 1 };
    
    static float k_SobelY[9] = { -1, -2, -1,
        0,  0,  0,
        1,  2,  1 };
    
    CFloatImage ConvolveKernel_SobelX(3, 3, 1);
    CFloatImage ConvolveKernel_SobelY(3, 3, 1);
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ConvolveKernel_SobelX.Pixel(i,j,0) = k_SobelX[i * 3 + j];
            ConvolveKernel_SobelY.Pixel(i,j,0) = k_SobelY[i * 3 + j];
        }
    }
    
    CFloatImage Ix(srcImage.Shape());
    CFloatImage Iy(srcImage.Shape());
    CFloatImage IxIy(srcImage.Shape());
    CFloatImage Ix2(srcImage.Shape());
    CFloatImage Iy2(srcImage.Shape());
    
    Convolve(srcImage, Ix, ConvolveKernel_SobelX);
    Convolve(srcImage, Iy, ConvolveKernel_SobelY);

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            
            // TODO:  Compute the harris score for 'srcImage' at this pixel and store in 'harrisImage'.  See the project
            //   page for pointers on how to do this
            
            double ix = Ix.Pixel(x,y,0);
            //printf("%f\n", ix);
            double iy = Iy.Pixel(x,y,0);
            //printf("%f\n", iy);
            
            double ix2 = pow(ix, 2);
            double iy2 = pow(iy, 2);
            double ixiy = ix*iy;
            double iyix = iy*ix;
            
            Ix2.Pixel(x,y,0) = ix2;
            Iy2.Pixel(x,y,0) = iy2;
            IxIy.Pixel(x,y,0) = ixiy;
            
        }
    }
    
    //5x5 Gaussian
    const double gaussian5x5[25] = {0.003663, 0.014652,  0.025641,  0.014652,  0.003663,
        0.014652, 0.0586081, 0.0952381, 0.0586081, 0.014652,
        0.025641, 0.0952381, 0.150183,  0.0952381, 0.025641,
        0.014652, 0.0586081, 0.0952381, 0.0586081, 0.014652,
        0.003663, 0.014652,  0.025641,  0.014652,  0.003663 };
    
    CFloatImage GaussianKernel(5, 5, 1);
    CFloatImage wIx2(srcImage.Shape());
    CFloatImage wIy2(srcImage.Shape());
    CFloatImage wIxIy(srcImage.Shape());
    
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            GaussianKernel.Pixel(i,j,0) = gaussian5x5[i * 5 + j];
        }
    }
    
    Convolve(Ix2, wIx2, GaussianKernel);
    Convolve(Iy2, wIy2, GaussianKernel);
    Convolve(IxIy, wIxIy, GaussianKernel);
    
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            double wix2 = wIx2.Pixel(x,y,0);
            double wiy2 = wIy2.Pixel(x,y,0);
            double wixiy =  wIxIy.Pixel(x,y,0);
            
            double minEig = .5*((wix2 + wiy2) - sqrt(4 * wixiy * wixiy + pow(wix2-wiy2, 2)));
            //printf("%f\n", ix2+iy2);
            //printf("%f\n", sqrt(4 * wixiy * wixiy + pow(wix2-wiy2, 2)));
            double maxEig = .5*((wix2 + wiy2) + sqrt(4 * wixiy * wixiy + pow(wix2-wiy2, 2)));
            //printf("%f\n", maxEig);
    
            double c = 2*(minEig * maxEig)/(minEig + maxEig + 0.00000001); //corner strength function
    
            harrisImage.Pixel(x,y,0) = c;
            //printf("%f\n", c);
            //printf("hi");
            
            //compute dominant orientation for use in feature description
            orientationImage.Pixel(x,y,0) = (Ix.Pixel(x,y,0) > 0) ? atan((maxEig-wix2)/wixiy) : atan((maxEig-wix2/wixiy)) + PI;
        }
    }
}



// TO DO---------------------------------------------------------------------
// Loop through the harrisImage to threshold and compute the local maxima in a neighborhood
// srcImage:  image with Harris values
// destImage: Assign 1 to a pixel if it is above a threshold and is the local maximum in 3x3 window, 0 otherwise.
//    You'll need to find a good threshold to use.
void computeLocalMaxima(CFloatImage &srcImage,CByteImage &destImage, double thresh)
{

    
    int w = srcImage.Shape().width;
    int h = srcImage.Shape().height;
    
    double max = 0;
    double min = 0;
    for (int y = 1; y < h; y++) {
        for (int x = 1; x < w; x++) {
            if(srcImage.Pixel(x, y, 0) > max){
                max = srcImage.Pixel(x, y, 0);
            }
            
            if(srcImage.Pixel(x, y, 0) < min){
                min = srcImage.Pixel(x, y, 0);
            }
        }
    }
    printf("%f\n", max);
    printf("%f\n", min);
    double th = (min + .2 * (max-min))*thresh;
    double numFeatures = 0;
    
    //printf("compute local maxima\n");
    for(int y = 0 ; y < h ; y++){
        for(int x = 0 ; x < w ; x++){
            
            destImage.Pixel(x,y,0) = 0;
            double score = srcImage.Pixel(x,y,0);
            bool Max = false;
            if(score > th){
                Max = true;
            }
            
            for(int i = -2 ; i < 3 ; i++){
                for(int j = -2 ; j < 3 ; j++){
                    double locScore = imagePixel(x+j, y+i, srcImage, 0);
                    if(locScore > score){
                        Max = false;
                    }
                }
            }
            
            if(Max){
                destImage.Pixel(x,y,0) = 1;
                numFeatures++;
            }
        }
    }
    
    // make sure number of features is in between 500 and 1000– this produces the best results
    // for images given and allows for consistent number of features when testing multiple images
   // while (recur < 10){
    
        //printf("numFeatures: %d\n", recur);
        printf("%f\n", numFeatures);
        
        if (numFeatures > 1000){
            computeLocalMaxima(srcImage, destImage, thresh*1.5);
        }
        
        if (numFeatures < 500){
            computeLocalMaxima(srcImage, destImage, thresh*.5);
        }
                              
      //  recur++;
 //   }
}

// Perform simple feature matching.  This just uses the SSD
// distance between two feature vectors, and matches a feature in the
// first image with the closest feature in the second image.  It can
// match multiple features in the first image to the same feature in
// the second image.
void ssdMatchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore) {
	int m = f1.size();
	int n = f2.size();

	matches.resize(m);
	totalScore = 0;

	double d;
	double dBest;
	int idBest;

	for (int i=0; i<m; i++) {
		dBest = 1e100;
		idBest = 0;

		for (int j=0; j<n; j++) {
			d = distanceSSD(f1[i].data, f2[j].data);

			if (d < dBest) {
				dBest = d;
				idBest = f2[j].id;
			}
		}

        matches[i].id1 = f1[i].id;
		matches[i].id2 = idBest;
		matches[i].score = dBest;
		totalScore += matches[i].score;
	}
}

// TODO: Write this function to perform ratio feature matching.  
// This just uses the ratio of the SSD distance of the two best matches as the score
// and matches a feature in the first image with the closest feature in the second image.
// It can match multiple features in the first image to the same feature in
// the second image.  (See class notes for more information, and the sshMatchFeatures function above as a reference)
void ratioMatchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore) 
{
    int m = f1.size();
    int n = f2.size();
    
    double d;
    double dBest = -1;
    double dNextBest = -1;
    
    int idBest = 0;
    int idNextBest = 0;
    
    // Iterate through features in FeatureSet 1
    for(int i=0; i<m; i++)
    {
        FeatureMatch newMatch;
        
        //reset distances, id's
        idBest = 0;
        idNextBest = 0;
        dBest = -1;
        dNextBest = -1;
        
        // Iterate through features in FeatureSet 2
        for(int j=0; j<n; j++)
        {
            d = distanceSSD(f1[i].data, f2[j].data);
            
            // if d is the first distance or shortest distance, set nextBest to best
            if(dBest<0 || d<dBest)
            {
                idNextBest = idBest;
                dNextBest = dBest;
                
                idBest = f2[j].id;
                dBest = d;
            }
            
            // Else check if this is the second best distance
            else if(dNextBest<0 || d<dNextBest)
            {
                idNextBest = f2[j].id;
                dNextBest = d;
            }
        }
        
        // Check if best match and second-best match distances pass ratio test
        newMatch.id1 = f1[i].id;
        newMatch.id2 = idBest;
        newMatch.score = (dBest/dNextBest);
        matches.push_back(newMatch);
    }
}


// Convert Fl_Image to CFloatImage.
bool convertImage(const Fl_Image *image, CFloatImage &convertedImage) {
	if (image == NULL) {
		return false;
	}

	// Let's not handle indexed color images.
	if (image->count() != 1) {
		return false;
	}

	int w = image->w();
	int h = image->h();
	int d = image->d();

	// Get the image data.
	const char *const *data = image->data();

	int index = 0;

	for (int y=0; y<h; y++) {
		for (int x=0; x<w; x++) {
			if (d < 3) {
				// If there are fewer than 3 channels, just use the
				// first one for all colors.
				convertedImage.Pixel(x,y,0) = ((uchar) data[0][index]) / 255.0f;
				convertedImage.Pixel(x,y,1) = ((uchar) data[0][index]) / 255.0f;
				convertedImage.Pixel(x,y,2) = ((uchar) data[0][index]) / 255.0f;
			}
			else {
				// Otherwise, use the first 3.
				convertedImage.Pixel(x,y,0) = ((uchar) data[0][index]) / 255.0f;
				convertedImage.Pixel(x,y,1) = ((uchar) data[0][index+1]) / 255.0f;
				convertedImage.Pixel(x,y,2) = ((uchar) data[0][index+2]) / 255.0f;
			}

			index += d;
		}
	}
	
	return true;
}

// Convert CFloatImage to CByteImage.
void convertToByteImage(CFloatImage &floatImage, CByteImage &byteImage) {
	CShape sh = floatImage.Shape();

    assert(floatImage.Shape().nBands == byteImage.Shape().nBands);
	for (int y=0; y<sh.height; y++) {
		for (int x=0; x<sh.width; x++) {
			for (int c=0; c<sh.nBands; c++) {
				float value = floor(255*floatImage.Pixel(x,y,c) + 0.5f);

				if (value < byteImage.MinVal()) {
					value = byteImage.MinVal();
				}
				else if (value > byteImage.MaxVal()) {
					value = byteImage.MaxVal();
				}

				// We have to flip the image and reverse the color
				// channels to get it to come out right.  How silly!
				byteImage.Pixel(x,sh.height-y-1,sh.nBands-c-1) = (uchar) value;
			}
		}
	}
}

// Compute SSD distance between two vectors.
double distanceSSD(const vector<double> &v1, const vector<double> &v2) {
	int m = v1.size();
	int n = v2.size();

	if (m != n) {
		// Here's a big number.
		return 1e100;
	}

	double dist = 0;

	for (int i=0; i<m; i++) {
		dist += pow(v1[i]-v2[i], 2);
	}

	
	return sqrt(dist);
}

// Transform point by homography.
void applyHomography(double x, double y, double &xNew, double &yNew, double h[9]) {
	double d = h[6]*x + h[7]*y + h[8];

	xNew = (h[0]*x + h[1]*y + h[2]) / d;
	yNew = (h[3]*x + h[4]*y + h[5]) / d;
}

// Compute AUC given a ROC curve
double computeAUC(vector<ROCPoint> &results)
{
	double auc=0;
	double xdiff,ydiff;
	for (int i = 1; i < (int) results.size(); i++)
    {
        //fprintf(stream,"%lf\t%lf\t%lf\n",thresholdList[i],results[i].falseRate,results[i].trueRate);
		xdiff=(results[i].falseRate-results[i-1].falseRate);
		ydiff=(results[i].trueRate-results[i-1].trueRate);
		auc=auc+xdiff*results[i-1].trueRate+xdiff*ydiff/2;
    	    
    }
	return auc;
}

