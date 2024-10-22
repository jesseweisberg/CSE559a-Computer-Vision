///////////////////////////////////////////////////////////////////////////
//
// NAME
//  FeatureAlign.h -- image registration using feature matching
//
// SEE ALSO
//  FeatureAlign.h      longer description
//
// Copyright ?Richard Szeliski, 2001.  See Copyright.h for more details
// (modified for CSE576 Spring 2005)
//
///////////////////////////////////////////////////////////////////////////

#include "ImageLib/ImageLib.h"
#include "FeatureAlign.h"
#include <math.h>

/******************* TO DO *********************
* alignPair:
*	INPUT:
*		f1, f2: source feature sets
*		matches: correspondences between f1 and f2
*               *NOTE* Each match in 'matches' contains two feature ids of matching features, id1 (in f1) and id2 (in f2).
*               These ids are 1-based indices into the feature arrays,
*               so you access the appropriate features as f1[id1-1] and f2[id2-1].
*		m: motion model
*		f: focal length
*		nRANSAC: number of RANSAC iterations
*		RANSACthresh: RANSAC distance threshold
*		M: transformation matrix (output)
*	OUTPUT:
*		repeat for nRANSAC iterations:
*			choose a minimal set of feature matches
*			estimate the transformation implied by these matches
*			count the number of inliers
*		for the transformation with the maximum number of inliers,
*		compute the least squares motion estimate using the inliers,
*		and store it in M
*/
int alignPair(const FeatureSet &f1, const FeatureSet &f2,
              const vector<FeatureMatch> &matches, MotionModel m, float f,
              int nRANSAC, double RANSACthresh, CTransform3x3& M)
{
    // BEGIN TODO
    // write this entire method
    int bestInliersCount = 0;
    
    //fill setMatches with nRANSAC random matches
    for (int i=0; i<nRANSAC; i++){
        
        vector<int> sampledMatch;
        int r = rand() % matches.size();
        sampledMatch.push_back(r);
        //printf("\nRandom Sample %d\n", r);
        //printf("\nSample Match %d\n", matches[r]);
        //printf("\nTotal Matches %d\n", matches.size());
        
        /* Generate the transformation from the sampledMatches
         use leastSquaresFit rather than duplicate the code */
        CTransform3x3 tempM;
        leastSquaresFit(f1, f2, matches, m, f, sampledMatch, tempM);
        
        /* Count the inliers*/
        vector<int> inliers;
        countInliers(f1,f2,matches,m, f, tempM, RANSACthresh, inliers);
        
        /*
        for(int i=0; i<inliers.size(); i++){
            printf("\ninlier %d\n", inliers[i]);
        }
        */
        
        /* If the count is better, recalculate using all the inliers
         and save the transformation */
        if (inliers.size() > bestInliersCount){
            bestInliersCount = inliers.size();
            leastSquaresFit(f1,f2,matches, m, f, inliers, M);
        }
        
    }
    
    // END TODO

    return 0;
}

/******************* TO DO *********************
* countInliers:
*	INPUT:
*		f1, f2: source feature sets
*		matches: correspondences between f1 and f2
*               *NOTE* Each match contains two feature ids of matching features, id1 (in f1) and id2 (in f2).
*               These ids are 1-based indices into the feature arrays,
*               so you access the appropriate features as f1[id1-1] and f2[id2-1].
*		m: motion model
*		f: focal length
*		M: transformation matrix
*		RANSACthresh: RANSAC distance threshold
*		inliers: inlier feature IDs
*	OUTPUT:
*		transform the matched features in f1 by M
*
*		count the number of matching features for which the transformed
*		feature f1[id1-1] is within SSD distance RANSACthresh of its match
*		f2[id2-1]
*
*		store the indices of these matches in inliers
*
*		
*/
int countInliers(const FeatureSet &f1, const FeatureSet &f2,
                 const vector<FeatureMatch> &matches, MotionModel m, float f,
                 CTransform3x3 M, double RANSACthresh, vector<int> &inliers)
{
    inliers.clear();
    int count = 0;

    for (unsigned int i=0; i<(int) matches.size(); i++) {
        // BEGIN TODO
        // determine if the ith matched feature f1[id1-1], when transformed by M,
        // is within RANSACthresh of its match in f2
        //
        // if so, increment count and append i to inliers
        //
        // *NOTE* Each match contains two feature ids of matching features, id1 and id2.
        //        These ids are 1-based indices into the feature arrays,
        //        so you access the appropriate features as f1[id1-1] and f2[id2-1].
        
        //fp1, fp2 designate feature points 1, 2
        CVector3 fp1, fp2;
        
        fp1[0] = f1[matches[i].id1-1].x;
        fp1[1] = f1[matches[i].id1-1].y;
        fp1[2] = 1;
    
        fp2[0] = f2[matches[i].id2-1].x;
        fp2[1] = f2[matches[i].id2-1].y;
        fp2[2] = 1;
        
        //transform p1 to the plane of p2 (translate fp1.x,y by u,v and normalize to z=1 plane)
        fp1[0] += M[0][2];
        fp1[1] += M[1][2];
        
        //fp1[0] = fp1[0]/fp1[2];
        //fp1[1] = fp1[1]/fp1[2];
        
        double dist = sqrt(pow(fp1[0] - fp2[0], 2) + pow(fp1[1] - fp2[1], 2));
        //printf("\n dist: %f\n", dist);
        
        if (dist < RANSACthresh){
            count+=1;
            inliers.push_back(i);
            //printf("\n** ** Inlier %d: %f ** **", i, dist);
        }
        
        //else {printf("\nOutlier %d: %f", i, dist);}
        

        // END TODO
    }
    //printf("\n\nInlier Count: %d\n", count);
    
    return count;
}

/******************* TO DO *********************
* leastSquaresFit:
*	INPUT:
*		f1, f2: source feature sets
*		matches: correspondences between f1 and f2
*		m: motion model
*		f: focal length
*		inliers: inlier match indices (indexes into 'matches' array)
*		M: transformation matrix (output)
*	OUTPUT:
*		compute the transformation from f1 to f2 using only the inliers
*		and return it in M
*/
int leastSquaresFit(const FeatureSet &f1, const FeatureSet &f2,
                    const vector<FeatureMatch> &matches, MotionModel m, float f,
                    const vector<int> &inliers, CTransform3x3& M)
{
    // for project 2, the transformation is a translation and
    // only has two degrees of freedom
    //
    // therefore, we simply compute the average translation vector
    // between the feature in f1 and its match in f2 for all inliers
    double u = 0;
    double v = 0;

    for (int i=0; i<inliers.size(); i++) {
        double xTrans, yTrans;

        // BEGIN TODO
        // compute the translation implied by the ith inlier match
        // and store it in (xTrans,yTrans)
        
        FeatureMatch match = matches[inliers[i]];
        
        xTrans = f2[match.id2-1].x - f1[match.id1-1].x;
        yTrans = f2[match.id2-1].y - f1[match.id1-1].y;

        // END TODO
        
        u += xTrans;
        v += yTrans;
        
    }

    u /= inliers.size();
    v /= inliers.size();

    M[0][0] = 1;
    M[0][1] = 0;
    M[0][2] = u;
    M[1][0] = 0;
    M[1][1] = 1;
    M[1][2] = v;
    M[2][0] = 0;
    M[2][1] = 0;
    M[2][2] = 1;

    return 0;
}
