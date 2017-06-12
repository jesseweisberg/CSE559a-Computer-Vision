
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: main.cpp                                                                             //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////




#include "stdafx.h"





Main::Main(int argc, char* argv[])
:
args(argc, argv),
programName(argv[0])
{
	commandFunctions[int(INVALID_COMMAND)]=&Main::invalidCommand;
	commandFunctions[int(COMPUTE_EIGENFACES)]=&Main::computeEigenfaces;
	commandFunctions[int(PROJECT_FACE)]=&Main::projectFace;
	commandFunctions[int(CONSTRUCT_USERBASE)]=&Main::constructUserbase;
	commandFunctions[int(IS_FACE)]=&Main::isFace;
	commandFunctions[int(VERIFY_FACE)]=&Main::verifyFace;
	commandFunctions[int(RECOGNIZE_FACE)]=&Main::recognizeFace;
	commandFunctions[int(FIND_FACE)]=&Main::findFace;
	commandFunctions[int(MORPH_FACES)]=&Main::morphFaces;

}

void Main::usage() const
{
	std::cout << programName                                                               << std::endl;
	std::cout << "has the following options:"                                              << std::endl;
	std::cout << "   Compute p Eigenfaces:"                                                << std::endl;
	std::cout << "    --eigenfaces <p> <width> <height> <textfile image list>"             << std::endl;
	std::cout << "      <outfile for eigenfaces>"                                          << std::endl;
	std::cout                                                                              << std::endl;
	std::cout << "   Project a face to the Eigenfaces:"                                    << std::endl;
	std::cout << "    --projectface <targa filename> <eigenfaces filename>"                << std::endl;
	std::cout                                                                              << std::endl;
	std::cout << "   Construct user database (coefficients vector+name for each user):"    << std::endl;
	std::cout << "    --constructuserbase <eigenfaces filename> <textfile image list>"     << std::endl;
	std::cout << "      <database output filename>"                                        << std::endl;
	std::cout                                                                              << std::endl;
	std::cout << "   Determine if an image is a face:"                                     << std::endl;
	std::cout << "    --isface <targa filename> <eigenfaces filename>"                     << std::endl;
	std::cout << "      <maximum face reconstruction MSE>"                                 << std::endl;
	std::cout                                                                              << std::endl;
	std::cout << "   Verify an identity (determine if a face is person x's):"              << std::endl;
	std::cout << "    --verifyface <targa filename> <userbase file> <person name>"         << std::endl;
	std::cout << "      <eigenfaces filename> <maximum v coefficients MSE>"                << std::endl;
	std::cout                                                                              << std::endl;
	std::cout << "   Find user given face:"                                                << std::endl;
	std::cout << "    --recognizeface <targa filename> <userbase file> <eigenfaces file>"  << std::endl;
	std::cout << "    <matches>"                                                           << std::endl;
	std::cout                                                                              << std::endl;
	std::cout << "   Find face in image:"                                                  << std::endl;
	std::cout << "    --findface <targa filename> <eigenfaces file> <min scale>"           << std::endl;
	std::cout << "    <max scale> <scale step> crop|mark <number to find> <result file>"   << std::endl;
	std::cout                                                                              << std::endl;
	std::cout << "   Morph two faces"                                                      << std::endl;
	std::cout << "    --morphfaces <targa filename1> <targa filename2>"                    << std::endl;
	std::cout << "      <eigenfaces filename>"                                             << std::endl;
	std::cout << "      <distance (0.0=first, 1.0=second, in between is a morph)>"         << std::endl;
	std::cout << "      <outfile name>"                                                    << std::endl;
	//std::cout << "      <ect....>"                                                       << std::endl;
}

void Main::run()
{
	if (args.moreArgs()) {
		Command c = args.getCommand();
		// call the function from the pointer to member function found in the array at location c
		(this->*commandFunctions[int(c)])();
	}
	else {
		usage();
	}
}

// Compute p Eigenfaces:
// --eigenfaces <p> <width> <height> <textfile image list> <outfile for eigenfaces>
void Main::computeEigenfaces()
{
	int count= args.nextPositiveInt();
	int width = args.nextPositiveInt();
	int height= args.nextPositiveInt();
	int n=0;
	std::string filelist = args.nextArg();
	std::string facefile=args.nextArg();

	{
		TextFileReader file(filelist);
		std::string image_name=file.readString();
		while (image_name!="") {
			n++;
			image_name=file.readString();
		}
	}
	Faces faces(n, width, height);
	EigFaces efaces(count, width, height);
	TextFileReader file(filelist);
	std::string image_name;
	for (int i=0; i<n; i++) {
		std::string filename=file.readString();
		Face face(width, height);
		face.loadTarga(filename);
		faces[i]=face;
	}
	faces.eigenFaces(efaces, count);
	efaces.output("eigen_face_%%.tga");
	Image average;
	efaces.getAverage().normalize(0.0, 255.0, average);
	average.saveTarga("average_face.tga");
	efaces.save(facefile);

}

// Project a face to the Eigenfaces:
// --projectface <targa filename> <eigenfaces filename>
void Main::projectFace()
{
	std::string filename=args.nextArg();
	std::string facesfile=args.nextArg();
	EigFaces eigfaces;
	eigfaces.load(facesfile);
	Face face(eigfaces.getWidth(), eigfaces.getHeight());
	face.loadTarga(filename);
	Vector coefficients;
	eigfaces.projectFace(face, coefficients);
	Face result;
	eigfaces.constructFace(coefficients, result);
	eigfaces.getAverage().saveTarga("average_face.tga");
	result.saveTarga("reconstruction.tga");
	face.saveTarga("original.tga");
}

// Construct user database (coefficients vector+name for each user):
// --constructuserbase <eigenfaces filename> <textfile image list> <database output filename>
void Main::constructUserbase()
{
	std::string eigenfacesname=args.nextArg();
	std::string imagelistname=args.nextArg();
	std::string databasename=args.nextArg();
	EigFaces eigfaces;
	eigfaces.load(eigenfacesname);
	int n=0;
	{
		TextFileReader file(imagelistname);
		std::string image_name=file.readString();
		while (image_name!="") {
			n++;
			image_name=file.readString();
		}
	}
	Users userbase;
	TextFileReader file(imagelistname);
	for (int i=0; i<n; i++) {
		std::string image_name=file.readString();
		Face face(eigfaces.getWidth(), eigfaces.getHeight());
		face.loadTarga(image_name);
		User user(image_name.substr(0, image_name.length()-4), eigfaces.getHeight()*eigfaces.getWidth());
		eigfaces.projectFace(face, user);
		userbase.addUser(user);
	}
	userbase.save(databasename);
}

// Determine if an image is a face:
// --isface <targa filename> <eigenfaces filename> <maximum face reconstruction MSE>
void Main::isFace()
{
	std::string filename = args.nextArg();
	std::string eigfacesname = args.nextArg();
	double max_reconstruction_mse=args.nextPositiveFloat();
	EigFaces eigenfaces;
	eigenfaces.load(eigfacesname);
	Face face(eigenfaces.getWidth(), eigenfaces.getHeight());
	face.loadTarga(filename);
	double mse;
	bool isface=eigenfaces.isFace(face, max_reconstruction_mse, mse);
	std::cout << "Image '" << filename << "' is ";
	if (isface) {
		std::cout << "a face; ";
	}
	else {
		std::cout << "not a face; ";
	}
	std::cout << "MSE: " << mse << std::endl;
}

// Verify an identity (determine if a face is person x's):
// --verifyface <targa filename> <userbase file> <person name> <eigenfaces filename> <maximum v coefficients MSE>
double fp;
double fn;

void Main::verifyFace()
{
	std::string filename = args.nextArg();
	std::string userbasename = args.nextArg();
	std::string username = args.nextArg();
	std::string eigfilename = args.nextArg();
	double max_coefficients_mse = args.nextPositiveFloat();
	EigFaces eigenfaces;
	eigenfaces.load(eigfilename);
	Users userbase;
	userbase.load(userbasename);
	Face face(eigenfaces.getWidth(), eigenfaces.getHeight());
	face.loadTarga(filename);
	const User& user= userbase[username];
	double mse;
	bool verified=eigenfaces.verifyFace(face, user, max_coefficients_mse, mse);
    
    
    if (verified){
        std::cout << "true\n";
    }
    else{
        std::cout << "false\n";
    }
    
	//std::cout << "MSE: " << mse << std::endl;
	//std::cout << "Image '" << filename << "' is";
	//std::cout << (verified?" ":" not ");
	//std::cout << "a picture of " << username << "'s face." << std::endl;
}

// Find user given face
// --recognizeface <targa filename> <userbase file> <eigenfaces file> <matches>
void Main::recognizeFace()
{
	std::string filename=args.nextArg();
	std::string userbasename = args.nextArg();
	std::string eigfacesname = args.nextArg();
	int matches=args.nextPositiveInt();
	EigFaces eigenfaces;
	eigenfaces.load(eigfacesname);
	Users userbase;
	userbase.load(userbasename);
	Face face(eigenfaces.getWidth(), eigenfaces.getHeight());
	face.loadTarga(filename);
	std::string name;
	eigenfaces.recognizeFace(face, userbase);
	std::cout << "Face '" << filename << "' ";
	std::cout << "recognized as being closest too:" << std::endl;
	for (int i=0; i< matches; i++) {
		const User& user=userbase[i];
		std::cout << i << ": " << user.getName() << "; MSE: " << user.getMse() << std::endl;
	}
}

// Find face in image:
// --findface <targa filename> <eigenfaces file> <min scale> <max scale> <scale step> crop|mark <number to find> <result file>
void Main::findFace()
{
	std::string filename=args.nextArg();
	std::string eigfilename=args.nextArg();
	double min_scale = args.nextPositiveFloat();
	double max_scale = args.nextPositiveFloat();
	double step = args.nextPositiveFloat();
	bool crop = (args.nextExpect("crop", "mark")=="crop");
	int n = args.nextPositiveInt();
	std::string resultfile=args.nextArg();
	Image img;
	img.loadTarga(filename);
	EigFaces eigfaces;
	eigfaces.load(eigfilename);
	Image result;
	eigfaces.findFace(img, min_scale, max_scale, step, n, crop, result);
	result.saveTarga(resultfile);
}

// Morph two faces
// --morphfaces <targa filename1> <targa filename2> <eigenfaces filename> <distance (0.0=first, 1.0=second, in between is a morph)> <outfile name>
void Main::morphFaces()
{
	std::string filename1=args.nextArg();
	std::string filename2=args.nextArg();
	std::string eigfilename=args.nextArg();
	double distance = args.nextFloat();
	std::string outfilename=args.nextArg();
	EigFaces eigenfaces;
	eigenfaces.load(eigfilename);
	Face face1(eigenfaces.getWidth(), eigenfaces.getHeight());
	Face face2(eigenfaces.getWidth(), eigenfaces.getHeight());
	Face result;
	face1.loadTarga(filename1);
	face2.loadTarga(filename2);
	eigenfaces.morphFaces(face1, face2, distance, result);
	result.saveTarga(outfilename);
}

void Main::invalidCommand()
{
	throw Error("Command '%%' is not recognized", args.getArg());
}


