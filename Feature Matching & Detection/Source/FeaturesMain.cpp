#include <assert.h>

#include <fstream>
#include <FL/Fl.H>
#include <FL/Fl_Shared_Image.H>
#include "features.h"
#include "FeaturesUI.h"
#include "FeaturesDoc.h"

FeaturesUI *ui;
FeaturesDoc *doc;

void saveRocFile(const char* filename,vector<double> &thresholdList,vector<ROCPoint> &results);
void saveAUCFile(const char* filename,vector<double> &thresholdList,vector<ROCPoint> &results);

void convertToFloatImage(CByteImage &byteImage, CFloatImage &floatImage) {
	CShape sh = byteImage.Shape();

    assert(floatImage.Shape().nBands == min(byteImage.Shape().nBands, 3));
	for (int y=0; y<sh.height; y++) {
		for (int x=0; x<sh.width; x++) {
			for (int c=0; c<min(3,sh.nBands); c++) {
				float value = byteImage.Pixel(x,y,c) / 255.0f;

				if (value < floatImage.MinVal()) {
					value = floatImage.MinVal();
				}
				else if (value > floatImage.MaxVal()) {
					value = floatImage.MaxVal();
				}

				// We have to flip the image and reverse the color
				// channels to get it to come out right.  How silly!
				floatImage.Pixel(x,sh.height-y-1,min(3,sh.nBands)-c-1) = value;
			}
		}
	}
}

bool LoadImageFile(const char *filename, CFloatImage &image) 
{
	// Load the query image.
	Fl_Shared_Image *fl_image = Fl_Shared_Image::get(filename);

	if (fl_image == NULL) {
		// printf("couldn't load query image\n");
        CByteImage byteImage;
		ReadFile(byteImage, filename);

        CShape sh = byteImage.Shape();
        sh.nBands = 3;

        image = CFloatImage(sh);
        convertToFloatImage(byteImage, image);
        
        return true;
    } else {
	    CShape sh(fl_image->w(), fl_image->h(), 3);
	    image = CFloatImage(sh);

	    // Convert the image to the CImage format.
	    if (!convertImage(fl_image, image)) {
		    printf("couldn't convert image to RGB format\n");
		    return false;
	    }
    
        return true;
    }
}

// Compute the features for a single image.
int mainComputeFeatures(int argc, char **argv) {
	if ((argc < 4) || (argc > 5)) {
		printf("usage: %s computeFeatures imagefile featurefile [featuretype]\n", argv[0]);

		return -1;
	}

	// Use feature type 1 as default.
	int type = 1;

	if (argc > 4) {
		type = atoi(argv[4]);
	}

#if 0
	// Load the query image.
	Fl_Shared_Image *queryImage = Fl_Shared_Image::get(argv[2]);

	if (queryImage == NULL) {
		printf("couldn't load query image\n");
		return -1;
	}

    CShape sh(queryImage->w(), queryImage->h(), 3);
	CFloatImage floatQueryImage(sh);

	// Convert the image to the CImage format.
    printf("converting image to CImage format\n");
	if (!convertImage(queryImage, floatQueryImage)) {
		printf("couldn't convert image to RGB format\n");
		return -1;
	}
#else
    CFloatImage floatQueryImage;
    bool success = LoadImageFile(argv[2], floatQueryImage);

    if (!success) {
        printf("couldn't load query image\n");
        return -1;
    }
#endif

	// Compute the image features.
	FeatureSet features;
    printf("computing featureSet\n");
    
	computeFeatures(floatQueryImage, features, type);
    printf("features computed\n");
    
	// Save the image features.
	features.save(argv[3]);
    printf("features saved\n");
    
	return 0;
}

// Match the features of one image to another, the output file matches to a file
int mainMatchFeatures(int argc, char **argv) {
	if ((argc < 6) || (argc > 7)) {
		printf("usage: %s matchFeatures featurefile1 featurefile2 threshold matchfile [matchtype]\n", argv[0]);
		return -1;
	}

	// Use match type 1 as default.
	int type = 1;

	if (argc > 6) {
		type = atoi(argv[6]);
	}

	FeatureSet f1;
	FeatureSet f2;

	if (!f1.load(argv[2])) {
		printf("couldn't load feature file %s\n", argv[2]);
		return -1;
	}

	if (!f2.load(argv[3])) {
		printf("couldn't load feature file %s\n", argv[3]);
		return -1;
	}

    double threshold = atof(argv[4]);

	vector<FeatureMatch> matches;
	double totalScore;

	// Compute the match.
	if (!matchFeatures(f1, f2, matches, totalScore, type)) {
		printf("matching failed, probably due to invalid match type\n");
		return -1;
	}

    // Output the matches 
    const char *matchFile = argv[5];
    FILE *f = fopen(matchFile, "w");

    // Count number of matches

    int num_matches = matches.size();
    int num_good_matches = 0;
    for (int i = 0; i < num_matches; i++) {
        if (matches[i].score < threshold) {
            num_good_matches++;
        }
    }

    fprintf(f, "%d\n", num_good_matches);
    for (int i = 0; i < num_matches; i++) {
        if (matches[i].score < threshold) {
            fprintf(f, "%d %d %lf\n", matches[i].id1, matches[i].id2, matches[i].score);
        }
    }

    fclose(f);

	return 0;
}

// Match the features of one image to another, the output file matches to a file
int mainMatchSIFTFeatures(int argc, char **argv) {
	if ((argc < 6) || (argc > 7)) {
		printf("usage: %s matchFeatures featurefile1 featurefile2 threshold matchfile [matchtype]\n", argv[0]);
		return -1;
	}

	// Use match type 1 as default.
	int type = 1;

	if (argc > 6) {
		type = atoi(argv[6]);
	}

	FeatureSet f1;
	FeatureSet f2;

	if (!f1.load_sift(argv[2])) {
		printf("couldn't load feature file %s\n", argv[2]);
		return -1;
	}

	if (!f2.load_sift(argv[3])) {
		printf("couldn't load feature file %s\n", argv[3]);
		return -1;
	}

    double threshold = atof(argv[4]);

	vector<FeatureMatch> matches;
	double totalScore;

	// Compute the match.
	if (!matchFeatures(f1, f2, matches, totalScore, type)) {
		printf("matching failed, probably due to invalid match type\n");
		return -1;
	}

    // Output the matches 
    const char *matchFile = argv[5];
    FILE *f = fopen(matchFile, "w");

    // Count number of matches

    int num_matches = matches.size();
    int num_good_matches = 0;
    for (int i = 0; i < num_matches; i++) {
        if (matches[i].score < threshold) {
            num_good_matches++;
        }
    }

    fprintf(f, "%d\n", num_good_matches);
    for (int i = 0; i < num_matches; i++) {
        if (matches[i].score < threshold) {
            fprintf(f, "%d %d %lf\n", matches[i].id1, matches[i].id2, matches[i].score);
        }
    }

    fclose(f);

	return 0;
}

// Match the features of one image to another, then compare the match
// with a ground truth homography.
int mainTestMatch(int argc, char **argv) {
	if ((argc < 5) || (argc > 6)) {
		printf("usage: %s testMatch featurefile1 featurefile2 homographyfile [matchtype]\n", argv[0]);

		return -1;
	}

	// Use feature type 1 as default.
	int type = 1;

	if (argc > 5) {
		type = atoi(argv[5]);
	}

	FeatureSet f1;
	FeatureSet f2;

	if (!f1.load(argv[2])) {
		printf("couldn't load feature file %s\n", argv[2]);
		return -1;
	}

	if (!f2.load(argv[3])) {
		printf("couldn't load feature file %s\n", argv[3]);
		return -1;
	}

	double h[9];

	ifstream is(argv[4]);

	if (!is.is_open()) {
		printf("couldn't open homography file %s\n", argv[4]);
		return -1;
	}

	// Read in the homography matrix.
	is >> h[0] >> h[1] >> h[2];
	is >> h[3] >> h[4] >> h[5];
	is >> h[6] >> h[7] >> h[8];

	vector<FeatureMatch> matches;
	double totalScore;

	// Compute the match.
	if (!matchFeatures(f1, f2, matches, totalScore, type)) {
		printf("matching failed, probably due to invalid match type\n");
		return -1;
	}

	double d = evaluateMatch(f1, f2, matches, h);

	// The total error is the average SSD distance between a
	// (correctly) transformed feature point in the first image and its
	// matched feature point in the second image.
	printf("%f\n", d);

	return 0;
}

// Match the features of one image to another, then compare the match
// with a ground truth homography. Compute the ROC points for various thresholds.
int mainRocTestMatch(int argc, char **argv) {
	if ((argc < 7) || (argc > 8)) {
		printf("usage: %s roc featurefile1 featurefile2 homographyfile [matchtype] rocfile aucfile\n", argv[0]);

		return -1;
	}

	// Use feature type 1 as default.
	int type = 1;

	const char* filename;
	const char* aucfilename;

	if (argc == 8) {
		type = atoi(argv[5]);
		filename=argv[6];
		aucfilename=argv[7];
	}

	if (argc==7)
	{
		filename=argv[5];
		aucfilename=argv[6];
	}


	FeatureSet f1;
	FeatureSet f2;

	if (!f1.load(argv[2])) {
		printf("couldn't load feature file %s\n", argv[2]);
		return -1;
	}

	if (!f2.load(argv[3])) {
		printf("couldn't load feature file %s\n", argv[3]);
		return -1;
	}

	double h[9];

	ifstream is(argv[4]);

	if (!is.is_open()) {
		printf("couldn't open homography file %s\n", argv[4]);
		return -1;
	}

	// Read in the homography matrix.
	is >> h[0] >> h[1] >> h[2];
	is >> h[3] >> h[4] >> h[5];
	is >> h[6] >> h[7] >> h[8];

	vector<FeatureMatch> matches;
	vector<bool> isMatch;
	double totalScore;
	double maxDistance=0;

	// Compute the match.
	if (!matchFeatures(f1, f2, matches, totalScore, type)) {
		printf("matching failed, probably due to invalid match type\n");
		return -1;
	}

	//double d = evaluateMatch(f1, f2, matches, h);
	addRocData(f1,f2,matches,h,isMatch,5,maxDistance);

	vector<double> thresholdList;

	for (int i=0;i<102;i++)
	{
		thresholdList.push_back(maxDistance/100.0*i);
	}

	vector<ROCPoint> results=computeRocCurve(matches,isMatch,thresholdList);

	// Write the ROC data to a vile
    saveRocFile(filename,thresholdList,results);

	printf("\nroc file complete.\n");

	// Write the AUC data to a file
    saveAUCFile(aucfilename,thresholdList,results);

	printf("\nauc file complete.\n");
	
	return 0;
}
// Match the SIFT features of one image to another, then compare the
// match with a ground truth homography.
int mainTestSIFTMatch(int argc, char **argv) {
	if ((argc < 5) || (argc > 6)) {
		printf("usage: %s testSIFTMatch featurefile1 featurefile2 homographyfile [matchtype]\n", argv[0]);

		return -1;
	}

	// Use feature type 1 as default.
	int type = 1;

	if (argc > 5) {
		type = atoi(argv[5]);
	}

	FeatureSet f1;
	FeatureSet f2;

	if (!f1.load_sift(argv[2])) {
		printf("couldn't load feature file %s\n", argv[2]);
		return -1;
	}

	if (!f2.load_sift(argv[3])) {
		printf("couldn't load SIFT feature file %s\n", argv[3]);
		return -1;
	}

	double h[9];

	ifstream is(argv[4]);

	if (!is.is_open()) {
		printf("couldn't open homography file %s\n", argv[4]);
		return -1;
	}

	// Read in the homography matrix.
	is >> h[0] >> h[1] >> h[2];
	is >> h[3] >> h[4] >> h[5];
	is >> h[6] >> h[7] >> h[8];

	vector<FeatureMatch> matches;
	double totalScore;

	// Compute the match.
	if (!matchFeatures(f1, f2, matches, totalScore, type)) {
		printf("matching failed, probably due to invalid match type\n");
		return -1;
	}

	double d = evaluateMatch(f1, f2, matches, h);

	// The total error is the average SSD distance between a
	// (correctly) transformed feature point in the first image and its
	// matched feature point in the second image.
	printf("%f\n", d);

	return 0;
}

// Match the SIFT features of one image to another, then compare the
// match with a ground truth homography.  Compute the ROC points for various thresholds.
int mainRocTestSIFTMatch(int argc, char **argv) {
	if ((argc < 7) || (argc > 8)) {
		printf("usage: %s rocSIFT featurefile1 featurefile2 homographyfile [matchtype] rocfile aucfile\n", argv[0]);

		return -1;
	}

	// Use feature type 1 as default.
	int type = 1;
	const char* filename;
	const char* aucfilename;

	if (argc == 8) {
		type = atoi(argv[5]);
		filename=argv[6];
		aucfilename=argv[7];
	}

	if (argc==7)
	{
		filename=argv[5];
		aucfilename=argv[6];
	}

	FeatureSet f1;
	FeatureSet f2;

	if (!f1.load_sift(argv[2])) {
		printf("couldn't load feature file %s\n", argv[2]);
		return -1;
	}

	if (!f2.load_sift(argv[3])) {
		printf("couldn't load SIFT feature file %s\n", argv[3]);
		return -1;
	}

	double h[9];

	ifstream is(argv[4]);

	if (!is.is_open()) {
		printf("couldn't open homography file %s\n", argv[4]);
		return -1;
	}

	// Read in the homography matrix.
	is >> h[0] >> h[1] >> h[2];
	is >> h[3] >> h[4] >> h[5];
	is >> h[6] >> h[7] >> h[8];

	vector<FeatureMatch> matches;
	vector<bool> isMatch;
	double maxDistance=0.0;
	double totalScore;

	// Compute the match.
	if (!matchFeatures(f1, f2, matches, totalScore, type)) {
		printf("matching failed, probably due to invalid match type\n");
		return -1;
	}

	//double d = evaluateMatch(f1, f2, matches, h);
	addRocData(f1,f2,matches,h,isMatch,5,maxDistance);

	
	
	vector<double> thresholdList;

	for (int i=0;i<102;i++)
	{
		thresholdList.push_back(maxDistance/100.0*i);
	}


	vector<ROCPoint> results=computeRocCurve(matches,isMatch,thresholdList);

	// Write the ROC data to a vile
    saveRocFile(filename,thresholdList,results);

	printf("\nSIFT roc file complete.\n");

	// Write the AUC data to a file
    saveAUCFile(aucfilename,thresholdList,results);

	printf("\nSIFT auc file complete.\n");
	

	return 0;
}

// Compute the features of all the images in one of the benchmark sets,
// then match the first image in the set with all of the others,
// comparing the resulting match with the ground truth homography.
int mainBenchmark(int argc, char **argv) {
	if ((argc != 3) && (argc != 5)) {
		printf("usage: %s benchmark imagedir [featuretype matchtype]\n", argv[0]);
		return -1;
	}

	int featureType = 1;
	int matchType = 1;

	if (argc == 5) {
		featureType = atoi(argv[3]);
		matchType = atoi(argv[4]);
	}

	// Get the directory containing the images.
	string imageDir(argv[2]);

	if ((imageDir[imageDir.size()-1] != '/') && (imageDir[imageDir.size()-1] != '\\')){
		imageDir += '/';
	}

	string imageFile;
	FeatureSet features[6];

	// Compute the features for each of the six images in the set.
	for (int i=0; i<6; i++) {
		imageFile = imageDir + "img" + (char)('1'+i) + ".ppm";

		// Load the query image.
		Fl_Shared_Image *tempImage = Fl_Shared_Image::get(imageFile.c_str());

		if (tempImage == NULL) {
			printf("couldn't load image %d\n", i+1);
			return -1;
		}

		CShape sh(tempImage->w(), tempImage->h(), 3);
		CFloatImage floatImage(sh);
 
		// Convert the image to the CImage format.
		if (!convertImage(tempImage, floatImage)) {
			printf("couldn't convert image %d to RGB format\n", i);
			return -1;
		}

		// Compute the image features.
		printf("computing features for image %d\n", i+1);
		computeFeatures(floatImage, features[i], featureType);
	}

	string homographyFile;
	double h[9];
	vector<FeatureMatch> matches;

	double d = 0;

	//AUC parameter
	double auc=0;
	vector<bool> isMatch;
	double maxDistance=0.0;

	// Match the first image with each of the other images, and compare
	// the results to the ground truth homography.
	for (int i=1; i<6; i++) {
		// Open the homography file.
		homographyFile = imageDir + "H1to" + (char)('1'+i) + "p";

		ifstream is(homographyFile.c_str());

		if (!is.is_open()) {
			printf("couldn't open homography file %s\n", homographyFile.c_str());
			return -1;
		}

		// Read in the homography matrix.
		is >> h[0] >> h[1] >> h[2];
		is >> h[3] >> h[4] >> h[5];
		is >> h[6] >> h[7] >> h[8];

		is.close();

		matches.clear();
		double totalScore;

		// Compute the match.
		printf("matching image 1 with image %d\n", i+1);
		if (!matchFeatures(features[0], features[i], matches, totalScore, matchType)) {
			printf("matching failed, probably due to invalid match type\n");
			return -1;
		}

		d += evaluateMatch(features[0], features[i], matches, h);

		//compute AUC
		isMatch.clear();
		maxDistance=0.0;
		addRocData(features[0],features[i],matches,h,isMatch,5,maxDistance);

		vector<double> thresholdList;

		for (int i=0;i<102;i++)
		{
			thresholdList.push_back(maxDistance/100.0*i);
		}

		vector<ROCPoint> results=computeRocCurve(matches,isMatch,thresholdList);
        printf("\n auc of image 1 with image %d\n: %f\n", i+1, computeAUC(results));
		auc=auc+computeAUC(results);
	}

	printf("\naverage error: %f pixels\n", d/5);
	printf("\naverage AUC: %f\n", auc/5);

    return 0;
}


void saveRocFile(const char* filename,vector<double> &thresholdList,vector<ROCPoint> &results)
{
	FILE *stream = fopen(filename, "wt");
    if (stream == 0)
        throw CError("saveRoc: could not open %s", filename);
   
    // Write out the numbers

	fprintf(stream,"Threshold\tFP Rate\tTP Rate\n");
	for (int i = 0; i < (int) results.size(); i++)
    {
        fprintf(stream,"%lf\t%lf\t%lf\n",thresholdList[i],results[i].falseRate,results[i].trueRate);
    	    
    }

    if (fclose(stream))
        throw CError("saveRoc: error closing file", filename);
}

void saveAUCFile(const char* filename,vector<double> &thresholdList,vector<ROCPoint> &results)
{
	double auc;
	FILE *stream = fopen(filename, "wt");
    if (stream == 0)
        throw CError("saveAUC: could not open %s", filename);
   
	auc=computeAUC(results);
	fprintf(stream,"%lf\n",auc);
    if (fclose(stream))
        throw CError("saveAUC: error closing file", filename);
}

int main(int argc, char **argv) {
	// This lets us load various image formats.
	fl_register_images();

	if (argc > 1) {
		if (strcmp(argv[1], "computeFeatures") == 0) {
			return mainComputeFeatures(argc, argv);
		}
		else if (strcmp(argv[1], "matchFeatures") == 0) {
			return mainMatchFeatures(argc, argv);
		}
		else if (strcmp(argv[1], "matchSIFTFeatures") == 0) {
			return mainMatchSIFTFeatures(argc, argv);
		}
		else if (strcmp(argv[1], "testMatch") == 0) {
			return mainTestMatch(argc, argv);
		}
		else if (strcmp(argv[1], "testSIFTMatch") == 0) {
			return mainTestSIFTMatch(argc, argv);
		}
		else if (strcmp(argv[1], "benchmark") == 0) {
			return mainBenchmark(argc, argv);
		}
		else if (strcmp(argv[1], "rocSIFT") == 0)
		{
			//return saveRoc(argc,argv);
			mainRocTestSIFTMatch(argc,argv);
		}
		else if (strcmp(argv[1], "roc") == 0)
		{
			//return saveRoc(argc,argv);
			mainRocTestMatch(argc,argv);
		}
		
		else {
			printf("usage:\n");
			printf("\t%s\n", argv[0]);
			printf("\t%s computeFeatures imagefile featurefile [featuretype]\n", argv[0]);
			printf("\t%s matchFeatures featurefile1 featurefile2 matchfile [matchtype]\n", argv[0]);
			printf("\t%s matchSIFTFeatures featurefile1 featurefile2 matchfile [matchtype]\n", argv[0]);
			// printf("\t%s testMatch featurefile1 featurefile2 homographyfile [matchtype]\n", argv[0]);
			// printf("\t%s testSIFTMatch featurefile1 featurefile2 homographyfile [matchtype]\n", argv[0]);
			// printf("\t%s benchmark imagedir [featuretype matchtype]\n", argv[0]);
			printf("\t%s rocSIFT featurefile1 featurefile2 homographyfile [matchtype] filename\n", argv[0]);
			printf("\t%s roc featurefile1 featurefile2 homographyfile [matchtype] filename\n", argv[0]);

			return -1;
		}
	}
	else {
		// Use the GUI.
		doc = new FeaturesDoc();
		ui = new FeaturesUI();

		ui->set_document(doc);
		doc->set_ui(ui);

		Fl::visual(FL_DOUBLE|FL_INDEX);

		ui->show();

		return Fl::run();
	}
}