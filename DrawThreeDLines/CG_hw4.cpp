/**
Course: CS536 - Computer Graphics
Student: Zhichao Cao
Email: zc77@drexel.edu
Title: CS536 - Homework 4
*/

#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
using namespace std ;

string readFile(string strPath);
double calSlope(int qx, int qy, int rx, int ry);
void bresenham(int qx, int qy, int rx, int ry, int sign, int index, double scal, int rota, int xTran, int yTran);
void writePixel(int x, int y, int sign, int index, int length);
void transform(double& x, double& y, int sign, int index, double scal, int rota, int xTran, int yTran, int umin, int vmin, int umax, int vmax, int xl, int yl, int xu, int yu);
int calXpoint(int qx, int rx);
void drawPic(string & strOutput, int countIndex, int xLowerBound, int xUpperBound, int yLowerBound, int yUpperBound);
void scanLineFill(int viewXLowerBound, int viewYLowerBound, int viewXUpperBound, int viewYUpperBound, int xyLnIndex);
void matrixMul(double a[4][4], double b[4][4]);
void matrixMulVer(double a[4][4], double b[4][1]);
void matrixSet(double a[4][4], double b[4][4]);

double xyOp[10000][4];  //Set of pairs of points
double xyOpSlope[10000];  //Set of slope of lines
double xOpPath[10000][2000];  //Set x-coop of points on line
double yOpPath[10000][2000];  //Set y-coop of points on line
double typeOpPath[10000];  //Type of slope of lines
double lenOpPath[10000];  //Number of points on line
int realWorldCoor[10000][2000];  //For cooperation of real world
int viewWorldCoor[10000][2000];  //For viewport
int scanLine[501][501];
int scanCount[501];
int scanLineCount[501];
int scanLineSet[501][501];
int sortedInter[501][501];
int yMax = -1;
int yMin = 2001;  //For scan-line filling algorithm
double xyThreeDV[10000][3];
int xyThreeDF[10000][3];
double xyThreeViewD[10000][2];
int xyThreeRealD[10000][2];
double tempMat[4][4];
double tempMatBack[4][4];
double tempMatVer[4][4];
double tempMatVerBack[4][4];
double nPar[4][4];
double nPer[4][4];
double paraProj[4][4];
double persProj[4][4];

const double PI = 3.14159265358979323846;

int main(int argc, char* argv[])
{
	string strOutput = "";  //String for output xpm file
	//string strOutputFile = "out.xpm";  //Output file name
	string strInputFile = "bound-lo-sphere.smf";  //-f
	double scalingFactor = 1.0;  //-s
	int cClockwiseRotation = 0;  //-r
	int xTranslation = 0;  //-m
	int yTranslation = 0;  //-n
	int xLowerBound = 0;  //-a
	int yLowerBound = 0;  //-b
	int xUpperBound = 500;  //-c
	int yUpperBound = 500;  //-d
	int viewXLowerBound = 0;  //-j
	int viewYLowerBound = 0;  //-k
	int viewXUpperBound = 500;  //-o
	int viewYUpperBound = 500;  //-p
	double xPRP = 0.0;  //-x
	double yPRP = 0.0;  //-y
	double zPRP = 1.0;  //-z
	double xVRP = 0.0;  //-X
	double yVRP = 0.0;  //-Y
	double zVRP = 0.0;  //-Z
	double xVPN = 0.0;  //-q
	double yVPN = 0.0;  //-r
	double zVPN = 1.0;  //-w
	double xVUP = 0.0;  //-Q
	double yVUP = 1.0;  //-R
	double zVUP = 0.0;  //-W
	double uMinVRC = -0.7;  //-u
	double vMinVRC = -0.7;  //-v
	double uMaxVRC = 0.7;  //-U
	double vMaxVRC = 0.7;  //-V
	int parallelProj = 0;  //-P  If this flag is not presented (be 0), use perspective projection, else (be 1), use parallel projection
	double backFace = -0.6;  //-B
	double frontFace = 0.6;  //-F
	string strPS = "";  //Content string in the .ps file
	string strSetPS[10000];  //Split strPS string
	int stIndex = 0;  
	int edIndex = 0;
	int stIndexOp = 0;
	string tempStrOp = "";
	int xyLnIndex = 0;  //Number of lines
	int xyThreeLnIndex = 0;
	int xyThreeLnIndexV = 0;
	int xyThreeLnIndexF = 0;
	int xyOpIndex = 0;  //Index for x,y-coop, for line command
	int negNum = 0;		//If it is negtive number of the input
	int fstPointIndex = 0;  //Temp index for first point of moveto and lineto command
	int sndPointIndex = 0;  //Temp index for second point of moveto and lineto command
	int countIndex = 0;  //Number of lines in .ps file
	string strLength = "";  //Length of pic
	string strWidth = "";  //Width of pic
	int lineType = 0;  //Type of command in .ps file. 1 for line, 2 for moveto, 3 for lineto

	double mn[3];
	double mvp[3];
	double mu[3];
	double mv[3];

	if(argv[1] == ">" && argv[2] == "out.xpm")
	//if(argv[1] == ">")
	{
		//./CG_hw1 -f hw1.ps -a 0 -b 0 -c 499 -d 499 -s 1.0 -m 0 -n 0 -r 0 > out.xpm
		strInputFile = "bound-lo-sphere.smf";
		scalingFactor = 1.0;
		cClockwiseRotation = 0;
		xTranslation = 0;
		yTranslation = 0;
		xLowerBound = 0;
		yLowerBound = 0;
		xUpperBound = 500;
		yUpperBound = 500;
		viewXLowerBound = 0;
		viewYLowerBound = 0;
		viewXUpperBound = 500;
		viewYUpperBound = 500;
		//strOutputFile = "out.xpm";
		xPRP = 0.0;
		yPRP = 0.0;
		zPRP = 1.0;
		xVRP = 0.0;
		yVRP = 0.0;
		zVRP = 0.0;
		xVPN = 0.0; 
	    yVPN = 0.0;
		zVPN = 1.0;
		xVUP = 0.0;
		yVUP = 1.0;
		zVUP = 0.0;
		uMinVRC = -0.7;
		vMinVRC = -0.7;
		uMaxVRC = 0.7;
		vMaxVRC = 0.7;
		parallelProj = 0;
		backFace = -0.6;
		frontFace = 0.6;
	}
	else
	{
		for(int i = 0; i < argc; i++)
		{
			string tempStr = "";
			tempStr = argv[i];
			if(tempStr == "-f")
			{
				strInputFile = argv[i+1];
			}
			else if(tempStr == "-s")
			{
				scalingFactor = atof(argv[i+1]);
			}
			/*else if(tempStr == "-r")
			{
				cClockwiseRotation = atoi(argv[i+1]);
			}*/
			else if(tempStr == "-m")
			{
				xTranslation = atoi(argv[i+1]);
			}
			else if(tempStr == "-n")
			{
				yTranslation = atoi(argv[i+1]);
			}
			else if(tempStr == "-a")
			{
				xLowerBound = atoi(argv[i+1]);
			}
			else if(tempStr == "-b")
			{
				yLowerBound = atoi(argv[i+1]);
			}
			else if(tempStr == "-c")
			{
				xUpperBound = atoi(argv[i+1]);
			}
			else if(tempStr == "-d")
			{
				yUpperBound = atoi(argv[i+1]);
			}
			else if(tempStr == "-j")
			{
				viewXLowerBound = atoi(argv[i+1]);
			}
			else if(tempStr == "-k")
			{
				viewYLowerBound = atoi(argv[i+1]);
			}
			else if(tempStr == "-o")
			{
				viewXUpperBound = atoi(argv[i+1]);
			}
			else if(tempStr == "-p")
			{
				viewYUpperBound = atoi(argv[i+1]);
			}
			else if(tempStr == ">")
			{
				//strOutputFile = argv[i+1];
			}
			else if(tempStr == "-x")
			{
				xPRP = atof(argv[i+1]);
			}
			else if(tempStr == "-y")
			{
				yPRP = atof(argv[i+1]);
			}
			else if(tempStr == "-z")
			{
				zPRP = atof(argv[i+1]);
			}
			else if(tempStr == "-X")
			{
				xVRP = atof(argv[i+1]);
			}
			else if(tempStr == "-Y")
			{
				yVRP = atof(argv[i+1]);
			}
			else if(tempStr == "-Z")
			{
				zVRP = atof(argv[i+1]);
			}
			else if(tempStr == "-q")
			{
				xVPN = atof(argv[i+1]);
			}
			else if(tempStr == "-r")
			{
				yVPN = atof(argv[i+1]);
			}
			else if(tempStr == "-w")
			{
				zVPN = atof(argv[i+1]);
			}
			else if(tempStr == "-Q")
			{
				xVUP = atof(argv[i+1]);
			}
			else if(tempStr == "-R")
			{
				yVUP = atof(argv[i+1]);
			}
			else if(tempStr == "-W")
			{
				zVUP = atof(argv[i+1]);
			}
			else if(tempStr == "-u")
			{
				uMinVRC = atof(argv[i+1]);
			}
			else if(tempStr == "-v")
			{
				vMinVRC = atof(argv[i+1]);
			}
			else if(tempStr == "-U")
			{
				uMaxVRC = atof(argv[i+1]);
			}
			else if(tempStr == "-V")
			{
				vMaxVRC = atof(argv[i+1]);
			}
			else if(tempStr == "-P")
			{
				parallelProj = 1;
			}
			else if(tempStr == "-B")
			{
				backFace = atoi(argv[i+1]);
			}
			else if(tempStr == "-F")
			{
				frontFace = atoi(argv[i+1]);
			}
		}
		//cout << strInputFile << " " << scalingFactor << " "  << cClockwiseRotation << " "  << xTranslation << " "  << yTranslation << " "  << xLowerBound << " "  << yLowerBound << " "  << xUpperBound << " "  << yUpperBound << " "  << endl;
	}
	double kx = static_cast<double>(viewXUpperBound - viewXLowerBound) / (2);
	double ky = static_cast<double>(viewYUpperBound - viewYLowerBound) / (2);
	double mTVRP[4][4] = {
				{1, 0, 0, -xVRP},
				{0, 1, 0, -yVRP},
				{0, 0, 1, -zVRP},
				{0, 0, 0, 1}
				};
	double mTPRP[4][4] = {
				{1, 0, 0, -xPRP},
				{0, 1, 0, -yPRP},
				{0, 0, 1, -zPRP},
				{0, 0, 0, 1}
				};
	double sqVPN = sqrt(xVPN * xVPN + yVPN * yVPN + zVPN * zVPN);
	double sqVUP = sqrt(xVUP * xVUP + yVUP * yVUP + zVUP * zVUP);
	mn[0] = xVPN / sqVPN;
	mn[1] = yVPN / sqVPN;
	mn[2] = zVPN / sqVPN;
	/*mvp[0] = xVUP / sqVUP;
	mvp[1] = yVUP / sqVUP;
	mvp[2] = zVUP / sqVUP;*/
	mvp[0] = xVUP;
	mvp[1] = yVUP;
	mvp[2] = zVUP;
	double ta = mn[1] * mvp[2] - mn[2] * mvp[1];
	double tb = mn[2] * mvp[0] - mn[0] * mvp[2];
	double tc = mn[0] * mvp[1] - mn[1] * mvp[0];
	double sqt = sqrt(ta * ta + tb * tb + tc * tc);
	mu[0] = ta / sqt;
	mu[1] = tb / sqt;
	mu[2] = tc / sqt;
	mv[0] = mn[1] * mu[2] - mn[2] * mu[1];
	mv[1] = mn[2] * mu[0] - mn[0] * mu[2];
	mv[2] = mn[0] * mu[1] - mn[1] * mu[0];
	double mR[4][4] = {
				{mu[0], mu[1], mu[2], 0},
				{mv[0], mv[1], mv[2], 0},
				{mn[0], mn[1], mn[2], 0},
				{0, 0, 0, 1}
				};
	double mSHpar[4][4] = {
				{1, 0, ((uMaxVRC + uMinVRC) / 2 - xPRP) / zPRP, 0},
				{0, 1, ((vMaxVRC + vMinVRC) / 2 - yPRP) / zPRP, 0},
				{0, 0, 1, 0},
				{0, 0, 0, 1}
				};
	double mTpar[4][4] = {
				{1, 0, 0, -(uMaxVRC + uMinVRC) / 2},
				{0, 1, 0, -(vMaxVRC + vMinVRC) / 2},
				{0, 0, 1, -frontFace},
				{0, 0, 0, 1}
				};
	double mSpar[4][4] = {
				{2 / (uMaxVRC - uMinVRC), 0, 0, 0},
				{0, 2 / (vMaxVRC - vMinVRC), 0, 0},
				{0, 0, 1 / (frontFace - backFace), 0},
				{0, 0, 0, 1}
				};
	// parallel
	double mMort[4][4] = {
			{1, 0, 0, 0},
			{0, 1, 0, 0},
			{0, 0, 0, 0},
			{0, 0, 0, 1}
			};
	// perspective
	double mMper[4][4] = {
			{1, 0, 0, 0},
			{0, 1, 0, 0},
			{0, 0, 1, 0},
			{0, 0, -1, 0} // to back plane
			};
	double negzPRP = -zPRP;
	double mSper[4][4] = {
				{2 * negzPRP / ((uMaxVRC - uMinVRC) * (negzPRP + backFace)), 0, 0, 0},
				{0, 2 * negzPRP / ((vMaxVRC - vMinVRC) * (negzPRP + backFace)), 0, 0},
				{0, 0, -1 / (negzPRP + backFace), 0},
				{0, 0, 0, 1}
				};
	if(parallelProj == 0)
	{//Nper = Sper.mTimes(SHpar.mTimes(TPRP.mTimes(R.mTimes(TVRP))));
		matrixMul(mR, mTVRP);
		matrixSet(tempMatBack, tempMat);
		matrixMul(mTPRP, tempMatBack);
		matrixSet(tempMatBack, tempMat);
		matrixMul(mSHpar, tempMatBack);
		matrixSet(tempMatBack, tempMat);
		matrixMul(mSper, tempMatBack);
		matrixSet(tempMatBack, tempMat);
		matrixSet(nPer, tempMat);
		matrixMul(mMper, tempMatBack);
		matrixSet(persProj, tempMat);
	}
	else
	{//Npar = Spar.mTimes(Tpar.mTimes(SHpar.mTimes(R.mTimes(TVRP))));
		matrixMul(mR, mTVRP);
		matrixSet(tempMatBack, tempMat);
		matrixMul(mSHpar, tempMatBack);
		matrixSet(tempMatBack, tempMat);
		matrixMul(mTpar, tempMatBack);
		matrixSet(tempMatBack, tempMat);
		matrixMul(mSpar, tempMatBack);
		matrixSet(tempMatBack, tempMat);
		matrixSet(nPar, tempMat);
		matrixMul(mMort, tempMatBack);
		matrixSet(paraProj, tempMat);
	}

	strPS = readFile(strInputFile);
	//cout << strPS << endl;

	for(int i = 0; i < strPS.length(); i++)
	{
		if(strPS[i] == '\n')
		{
			strSetPS[countIndex] = strPS.substr(stIndex, i - stIndex);
			stIndex = i + 1;
			countIndex++;
		}
	}
	/*for(int i = 0; i < countIndex; i++)
	{
		cout<<strSetPS[i]<<endl;
	}*/
	stIndex = 0;
	for(int i = 0; i < countIndex; i++) 
	{
		if(strSetPS[i].empty() != true && strSetPS[i] != "")
		{
			std::size_t posV = strSetPS[i].find("v");
			std::size_t posF = strSetPS[i].find("f");
			if(posV != std::string::npos)
			{
				lineType = 1;
			}
			else if(posF != std::string::npos)
			{
				lineType = 2;
			}
			for (int j = 0; j < strSetPS[i].length(); j++)
			{
				if(lineType == 1)
				{
					if(strSetPS[i][j] == ' ' && strSetPS[i][j - 1] == 'v')
					{
						stIndexOp = j + 1;
					}
					else if(strSetPS[i][j] == ' ' && strSetPS[i][stIndexOp] != '-')
					{
						tempStrOp = strSetPS[i].substr(stIndexOp, j - stIndexOp);
						xyThreeDV[xyThreeLnIndexV][xyOpIndex] = atof(tempStrOp.c_str());
						stIndexOp = j + 1;
						xyOpIndex++;
					}
					else if(strSetPS[i][j] == ' ' && strSetPS[i][stIndexOp] == '-')
					{
						stIndexOp = stIndexOp + 1;
						tempStrOp = strSetPS[i].substr(stIndexOp + 1, j - stIndexOp);
						xyThreeDV[xyThreeLnIndexV][xyOpIndex] = atof(tempStrOp.c_str());
						xyThreeDV[xyThreeLnIndexV][xyOpIndex] = 0 - xyThreeDV[xyThreeLnIndexV][xyOpIndex];
						stIndexOp = j + 1;
						xyOpIndex++;
					}
					else if(xyOpIndex == 2 && strSetPS[i][stIndex] != '-')
					{
						tempStrOp = strSetPS[i].substr(stIndexOp, strSetPS[i].length() - 1);
						xyThreeDV[xyThreeLnIndexV][xyOpIndex] = atof(tempStrOp.c_str());
						xyOpIndex++;
					}
					else if(xyOpIndex == 2 && strSetPS[i][stIndex] == '-')
					{
						tempStrOp = strSetPS[i].substr(stIndexOp + 1, strSetPS[i].length() - 1);
						xyThreeDV[xyThreeLnIndexV][xyOpIndex] = atof(tempStrOp.c_str());
						xyThreeDV[xyThreeLnIndexV][xyOpIndex] = 0 - xyThreeDV[xyThreeLnIndexV][xyOpIndex];
						xyOpIndex++;
					}
				}
				else if(lineType == 2)
				{
					if(strSetPS[i][j] == ' ' && strSetPS[i][j - 1] == 'f')
					{
						stIndexOp = j + 1;
					}
					else if(strSetPS[i][j] == ' ' && strSetPS[i][stIndex] != '-')
					{
						tempStrOp = strSetPS[i].substr(stIndexOp, j - stIndexOp);
						xyThreeDF[xyThreeLnIndexF][xyOpIndex] = atoi(tempStrOp.c_str());
						stIndexOp = j + 1;
						xyOpIndex++;
					}
					else if(strSetPS[i][j] == ' ' && strSetPS[i][stIndex] == '-')
					{
						stIndexOp = stIndexOp + 1;
						tempStrOp = strSetPS[i].substr(stIndexOp + 1, j - stIndexOp);
						xyThreeDF[xyThreeLnIndexF][xyOpIndex] = atoi(tempStrOp.c_str());
						xyThreeDF[xyThreeLnIndexF][xyOpIndex] = 0 - xyThreeDF[xyThreeLnIndexF][xyOpIndex];
						stIndexOp = j + 1;
						xyOpIndex++;
					}
					else if(xyOpIndex == 2 && strSetPS[i][stIndex] != '-')
					{
						tempStrOp = strSetPS[i].substr(stIndexOp, strSetPS[i].length() - 1);
						xyThreeDF[xyThreeLnIndexF][xyOpIndex] = atoi(tempStrOp.c_str());
						xyOpIndex++;
					}
					else if(xyOpIndex == 2 && strSetPS[i][stIndex] == '-')
					{
						tempStrOp = strSetPS[i].substr(stIndexOp + 1, strSetPS[i].length() - 1);
						xyThreeDF[xyThreeLnIndexF][xyOpIndex] = atoi(tempStrOp.c_str());
						xyThreeDF[xyThreeLnIndexF][xyOpIndex] = 0 - xyThreeDF[xyThreeLnIndexF][xyOpIndex];
						xyOpIndex++;
					}
				}
			}
			xyOpIndex = 0;
			stIndexOp = 0;
			if(lineType == 1)
			{
				xyThreeLnIndexV++;
			}
			else if(lineType == 2)
			{
				xyThreeLnIndexF++;
			}
			xyLnIndex++;
		}
	}

	/*for(int i = 0; i < xyLnIndex; i++)
	{
		cout << xyOp[i][0] << " " << xyOp[i][1] << " "  << xyOp[i][2] << " "  << xyOp[i][3] << endl;
	}*/
	/*for(int i = 0; i < xyThreeLnIndexV; i++)
	{
		cout << xyThreeDV[i][0] << " " << xyThreeDV[i][1] << " "  << xyThreeDV[i][2] << " " << xyThreeDV[i][3] << endl;
	}
	for(int i = 0; i < xyThreeLnIndexF; i++)
	{
		cout << xyThreeDF[i][0] << " " << xyThreeDF[i][1] << " "  << xyThreeDF[i][2] << endl;
	}*/

	for(int i = 0; i < xyThreeLnIndexV; i++)
	{
		double xabcd = xyThreeDV[i][0];
		double yabcd = xyThreeDV[i][1];
		double zabcd = xyThreeDV[i][2];
		double tP[4][1] = {{xabcd},{yabcd},{zabcd},{1}};
		double x1 = 0;
		double y1 = 0;
		if (parallelProj == 0) 
		{
			matrixMulVer(persProj, tP);
			double W = tempMatVer[3][0];
			x1 = tempMatVer[0][0] / W;
			y1 = tempMatVer[1][0] / W;
		}
		else
		{
			//Matrix pPrime = parallelProj.mTimes(P);
			matrixMulVer(paraProj, tP);
			x1 = tempMatVer[0][0];
			y1 = tempMatVer[1][0];
		}
		/*xyThreeRealD[i][0] = x1;
		xyThreeRealD[i][1] = y1;*/
		int x = (int)((double)((x1 + 1) * kx) + viewXLowerBound);
		int y = (int)((double)((y1 + 1) * ky) + viewYLowerBound);
		xyThreeRealD[i][0] = x;
		xyThreeRealD[i][1] = y;
	}

	//if(parallelProj == 0)
	//{
	//	//Perspective projection
	//	for(int i = 0; i < xyThreeLnIndexV; i++)
	//	{
	//		double zDis = zPRP - zVRP;
	//		xyThreeViewD[i][0] = xyThreeDV[i][0] / (1 - (xyThreeDV[i][2] / zDis));
	//		xyThreeViewD[i][1] = xyThreeDV[i][1] / (1 - (xyThreeDV[i][2] / zDis));
	//	}
	//}
	//else
	//{
	//	//Parallel projection
	//	for(int i = 0; i < xyThreeLnIndexV; i++)
	//	{
	//		xyThreeViewD[i][0] = xyThreeDV[i][0];
	//		xyThreeViewD[i][1] = xyThreeDV[i][1];
	//	}
	//}
	//for(int i = 0; i < xyThreeLnIndexV; i++)
	//{
	//	xyThreeRealD[i][0] = 250 + xyThreeViewD[i][0] / uMaxVRC * 250;
	//	xyThreeRealD[i][1] = 250 + xyThreeViewD[i][1] / vMaxVRC * 250;
	//	/*if(xyThreeViewD[i][0] < 0)
	//	{
	//		xyThreeRealD[i][0] = 250 + xyThreeViewD[i][0] / uMinVRC * 250;
	//	}
	//	else
	//	{
	//		xyThreeRealD[i][0] = 250 - xyThreeViewD[i][0] / uMaxVRC * 250;
	//	}
	//	if(xyThreeViewD[i][1] < 0)
	//	{
	//		xyThreeRealD[i][1] = 250 - xyThreeViewD[i][1] / vMinVRC * 250;
	//	}
	//	else
	//	{
	//		xyThreeRealD[i][1] = 250 + xyThreeViewD[i][1] / vMaxVRC * 250;
	//	}*/
	//}
	xyLnIndex = 0;
	for(int i = 0; i < xyThreeLnIndexF; i++)
	{
		//xyLnIndex
		int p1 = 0;
		int p2 = 0;
		int p3 = 0;
		p1 = xyThreeDF[i][0] - 1;
		p2 = xyThreeDF[i][1] - 1;
		p3 = xyThreeDF[i][2] - 1;
		xyOp[xyLnIndex][0] = xyThreeRealD[p1][0];
		xyOp[xyLnIndex][1] = xyThreeRealD[p1][1];
		xyOp[xyLnIndex][2] = xyThreeRealD[p2][0];
		xyOp[xyLnIndex][3] = xyThreeRealD[p2][1];
		xyLnIndex++;
		xyOp[xyLnIndex][0] = xyThreeRealD[p2][0];
		xyOp[xyLnIndex][1] = xyThreeRealD[p2][1];
		xyOp[xyLnIndex][2] = xyThreeRealD[p3][0];
		xyOp[xyLnIndex][3] = xyThreeRealD[p3][1];
		xyLnIndex++;
		xyOp[xyLnIndex][0] = xyThreeRealD[p3][0];
		xyOp[xyLnIndex][1] = xyThreeRealD[p3][1];
		xyOp[xyLnIndex][2] = xyThreeRealD[p1][0];
		xyOp[xyLnIndex][3] = xyThreeRealD[p1][1];
		xyLnIndex++;
	}

	/*for(int i = 0; i < xyThreeLnIndexV; i++)
	{
		cout << xyThreeViewD[i][0] << " " << xyThreeViewD[i][1] << endl;
	}
	for(int i = 0; i < xyThreeLnIndexV; i++)
	{
		cout << xyThreeRealD[i][0] << " " << xyThreeRealD[i][1] << endl;
	}*/

	/*double xScale = static_cast<double>(viewXUpperBound - viewXLowerBound) / (xUpperBound - xLowerBound);
	double yScale = static_cast<double>(viewYUpperBound - viewYLowerBound) / (yUpperBound - yLowerBound);
	cout<<viewXLowerBound<<" "<<viewXUpperBound<<" "<<viewYLowerBound<<" "<<viewYUpperBound<<endl;
	cout<<xLowerBound<<" "<<xUpperBound<<" "<<yLowerBound<<" "<<yUpperBound<<endl;
	cout<<xScale<<" "<<yScale<<endl;*/

	/*for (int i = 0; i < xyLnIndex; i++)
	{
		cout<<xyOp[i][0] << " " << xyOp[i][1] << " " << xyOp[i][2] << " "<< xyOp[i][3] << " " << i << endl;
	}*/

	for (int i = 0; i < xyLnIndex; i++)
	{
		double tempX, tempY = 0;			
		transform(xyOp[i][0], xyOp[i][1], 0, i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation, viewXLowerBound, viewYLowerBound, viewXUpperBound, viewYUpperBound, xLowerBound, yLowerBound, xUpperBound, yUpperBound);
		transform(xyOp[i][2], xyOp[i][3], 0, i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation, viewXLowerBound, viewYLowerBound, viewXUpperBound, viewYUpperBound, xLowerBound, yLowerBound, xUpperBound, yUpperBound);
		if(calXpoint(xyOp[i][0], xyOp[i][2]) == -1)
		{
			//xyOpSlope[i] = calSlope(xyOp[i][2], xyOp[i][3], xyOp[i][0], xyOp[i][1]);
			tempX = xyOp[i][0];
			tempY = xyOp[i][1];
			xyOp[i][0] = xyOp[i][2];
			xyOp[i][1] = xyOp[i][3];
			xyOp[i][2] = tempX;
			xyOp[i][3] = tempY;
		}
		//else
		//{
		//	xyOpSlope[i] = calSlope(xyOp[i][0], xyOp[i][1], xyOp[i][2], xyOp[i][3]);
		//}
		xyOpSlope[i] = calSlope(xyOp[i][0], xyOp[i][1], xyOp[i][2], xyOp[i][3]);
		
		//cout << xyOp[i][0] << " " << xyOp[i][1] << " "  << xyOp[i][2] << " "  << xyOp[i][3] << " " << xyOpSlope[i] << endl;
		
		if(xyOpSlope[i] >= 0 && xyOpSlope[i] <= 1)
		{
			typeOpPath[i] = 0;
			//No need for modification
			bresenham(xyOp[i][0], xyOp[i][1], xyOp[i][2], xyOp[i][3], typeOpPath[i], i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation);
		}
		else if(xyOpSlope[i] >= -1 && xyOpSlope[i] < 0)
		{
			typeOpPath[i] = 1;
			//Switch y1 and y2 to make 0<=slope<=1
			bresenham(xyOp[i][0], -xyOp[i][1], xyOp[i][2], -xyOp[i][3], typeOpPath[i], i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation);
		}
		else if(xyOpSlope[i] > 1)
		{
			typeOpPath[i] = 2;
			//(y1, x1)(y2, x2) to make 0<=slope<=1
			bresenham(xyOp[i][1], xyOp[i][0], xyOp[i][3], xyOp[i][2], typeOpPath[i], i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation);
		}
		else if(xyOpSlope[i] < -1)
		{
			typeOpPath[i] = 3;
			//(y2, -x2) and (y1, -x1) to make 0<=slope<=1
			bresenham(xyOp[i][3], -xyOp[i][2], xyOp[i][1], -xyOp[i][0], typeOpPath[i], i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation);
		}
		//else if(xyOpSlope[i] == 1000)
		//{
		//	for (int j = xyOp[i][1]; j <= xyOp[i][3]; j++)
		//	{
		//		writePixel(xyOp[i][0], j, 4, i, (j - xyOp[i][1]));
		//	}
		//	//bresenham(xyOp[i][0], xyOp[i][1], xyOp[i][2], xyOp[i][3], 4, i);
		//}
		//else if(xyOpSlope[i] == 1000)
		//{
		//	for (int j = xyOp[i][3]; j <= xyOp[i][1]; j++)
		//	{
		//		writePixel(xyOp[i][0], j, 5, i, (j - xyOp[i][3]));
		//	}
		//	//bresenham(xyOp[i][0], xyOp[i][1], xyOp[i][2], xyOp[i][3], 5, i);
		//}
	}
	for(int k = 0; k < xyLnIndex; k++)
	{
		int xCoop = 0;
		int yCoop = 0;
		for (int l = 0; l < lenOpPath[k]; l++)
		{					
			xCoop = xOpPath[k][l];
			yCoop = yOpPath[k][l];
			realWorldCoor[xCoop][yCoop] = 1;
			if((xCoop > viewXLowerBound && xCoop <= viewXUpperBound) && (yCoop > viewYLowerBound && yCoop <= viewYUpperBound))
			{
				viewWorldCoor[xCoop][yCoop] = 1;
				if(yCoop > yMax)
				{
					yMax = yCoop;
				}
				else if(yCoop < yMin)
				{
					yMin = yCoop;
				}
			}
		}
	}

	//scanLineFill(viewXLowerBound, viewYLowerBound, viewXUpperBound, viewYUpperBound, xyLnIndex);

	/*std::stringstream sl;
	sl << (xUpperBound - xLowerBound + 1);
	strLength = sl.str();
	std::stringstream sw;
	sw << (yUpperBound - yLowerBound + 1);
	strWidth = sw.str();*/


	strOutput = "/* XPM */\n";
	strOutput += "static char *quad_bw[] = {\n";
	strOutput += "/* columns rows colors chars-per-pixel */\n";
	//strOutput += "\"" + strLength + " " + strWidth + " 2 1\",\n";
	strOutput += "\" 501 501 2 1\",\n";
	strOutput += "/* pixels */\n";
	strOutput += "\"@ c #000000\",\n";
	strOutput += "\"  c #FFFFFF\",\n";

	drawPic(strOutput, xyLnIndex, xLowerBound, xUpperBound, yLowerBound, yUpperBound);

	cout << strOutput;

	return 0;
}

string readFile(string strPath)
{
	//Read the ps file
	string strPS = "";
	char chStrPath[20];
	strcpy(chStrPath, strPath.c_str());
	char buffer[256];
	ifstream hsfile;
	hsfile.open(chStrPath);
	if(!hsfile){
        cout << "Unable to open the file";
        exit(1);  // terminate with error
	}
	while (!hsfile.eof())
	{
		//strPS = hsfile.get();
		hsfile.getline(buffer,50);
		//cout << buffer << endl;
		strPS += buffer;
		strPS += "\n";
	}
	hsfile.close();
	return strPS;
}

void bresenham(int qx, int qy, int rx, int ry, int sign, int index, double scal, int rota, int xTran, int yTran)
{
	//Bresenham's algorithm
	int length_line = 0;
	int dx, dy, D, x, y = 0;
	//cout<<qx<<","<<qy<<" "<<rx<<","<<ry<<endl;
	dx = rx - qx;
	dy = ry - qy;
	D = 2 * dy - dx;
	y = qy;
	for (x = qx; x <= rx; x++)
	{
		writePixel(x, y, sign, index, length_line);
		//cout << x << "," << y << endl;
		if(D <= 0)
		{
			D += 2 * dy;
		}
		else
		{
			D += 2 * (dy -dx) ;
			y++;
		}
		length_line++;
	}
	lenOpPath[index] = length_line;
}

void writePixel(int x, int y, int sign, int index, int length)
{
	//Generating mapping point of the line
	switch (sign)
	{
		case 0:
			xOpPath[index][length] = x;
			yOpPath[index][length] = y;
			break;
		case 1:
			xOpPath[index][length] = x;
			yOpPath[index][length] = -y;
			break;
		case 2:	
			xOpPath[index][length] = y;
			yOpPath[index][length] = x;
			break;
		case 3:
			xOpPath[index][length] = -y;
			yOpPath[index][length] = x;
			break;
		/*case 4:
		case 5:
			xOpPath[index][length] = x;
			yOpPath[index][length] = y;
			break;*/
		default:
			xOpPath[index][length] = 0;
			yOpPath[index][length] = 0;
			cout << "Points fail" << endl;
			break;
	}
}

double calSlope(int qx, int qy, int rx, int ry)
{
	//Calculate the slope of the line
	double slope = 0.0;
	if(rx - qx != 0)
		slope = static_cast<double> (ry - qy) / (rx - qx);
	else
	{
		if(ry - qy < 0)
			slope = -1000;
		else if(ry - qy > 0)
			slope = 1000;
		else if (ry - qy == 0)
			slope = 0;
	}
	return slope;
}

int calXpoint(int qx, int rx)
{
	//Calculate the distance of points of x-coop
	if(rx - qx < 0)
	{
		return -1;
	}
	else if(rx - qx == 0)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

void drawPic(string & strOutput, int countIndex, int xLowerBound, int xUpperBound, int yLowerBound, int yUpperBound)
{
	//Building output string for xpm file
	int pExist = 0;
	int xCoop = 0;
	int yCoop = 0;
	yUpperBound = 500;
	yLowerBound = 0;
	xUpperBound = 500;
	xLowerBound = 0;
	for(int i = 0; i <= (yUpperBound - yLowerBound); i++)
	{
		//Clipping line by Sutherland-Hodgman Algorithm 
		strOutput += "\"";
		for (int j = xLowerBound; j <= xUpperBound; j++)
		{
			if(realWorldCoor[j][yUpperBound -i] == 1 && viewWorldCoor[j][yUpperBound -i] == 1)
			//if(realWorldCoor[j][yUpperBound -i] == 1)
			{
				strOutput += "@";
			}
			else
			{
				strOutput += " ";
			}
		}
		if(i == yUpperBound)
		{
			strOutput += "\"\n};";
		}
		else
		{
			strOutput += "\",\n";
		}
	}
}

void transform(double& x, double& y, int sign, int index, double scal, int rota, int xTran, int yTran, int umin, int vmin, int umax, int vmax, int xl, int yl, int xu, int yu)
{
	//Transformation of lines
	double xScale = static_cast<double>(umax - umin) / (xu - xl);
	double yScale = static_cast<double>(vmax - vmin) / (yu - yl);
	double tempX = x * scal;
	double tempY = y * scal;;
	double cosValue = cos((rota*PI)/180.0);
	double sinValue = sin((rota*PI)/180.0);
	x = tempX * cosValue - tempY * sinValue + xTran;
	y = tempX * sinValue + tempY * cosValue + yTran;
	/*x = (x - xl) * xScale + umin;
	y = (y - yl) * yScale + vmin;*/
}

void scanLineFill(int viewXLowerBound, int viewYLowerBound, int viewXUpperBound, int viewYUpperBound, int xyLnIndex)
{
	
	for(int i = 0; i < 501; i++)
	{
		//Initialization
		scanCount[i] = 0;
		scanLineCount[i] = 0;
		for(int j = 0; j < 501; j++)
		{
			scanLine[i][j] = 0;
			scanLineSet[i][j] = -1;
			sortedInter[i][j] = -1;
		}
	}
	for(int i = yMax; i >= yMin; i--)
	{
		int count = 0;
		for (int j = 0; j < xyLnIndex; j++)
		{
			if(xyOpSlope[j] == 0)
			{
				//Ignore
			}
			else if(xyOp[j][3] > xyOp[j][1] && i == xyOp[j][3])
			{
				//Ignore
			}
			else if(xyOp[j][3] < xyOp[j][1] && i == xyOp[j][1])
			{
				//Ignore
			}
			else if(xyOp[j][3] > xyOp[j][1] && (i < xyOp[j][3] && i >= xyOp[j][1]))
			{
				//Add to scan line set
				scanLine[count][i] = j;
				count++;
				scanLineCount[i]++;
			}
			else if(xyOp[j][3] < xyOp[j][1] && (i >= xyOp[j][3] && i < xyOp[j][1]))
			{
				//Add to scan line set
				scanLine[count][i] = j;
				count++;
				scanLineCount[i]++;
			}
		}
	}
	for(int i = yMax; i >= yMin; i--)
	{
		//Calculate intersections
		int count = 0;
		for(int j = 0; j < scanLineCount[i]; j++)
		{
			int lIndex = scanLine[j][i];
			for(int k = 0; k < 2000; k++)
			{
				if(yOpPath[lIndex][k] == i)
				{
					sortedInter[count][i] = xOpPath[lIndex][k];
					count++;
					scanCount[i]++;
					break;
				}
			}
		}
	}
	for(int i = yMax; i >= yMin; i--)
	{
		//Sort the intersections in x
		for(int j = 0; j < scanLineCount[i]; j++)
		{
			for(int k = 0; k < scanLineCount[i] - 1; k++)
			{
				int temp = 0;
				if(sortedInter[k][i] > sortedInter[k + 1][i])
				{
					temp = sortedInter[k][i];
					sortedInter[k][i] = sortedInter[k + 1][i];
					sortedInter[k + 1][i] = temp;
				}
			}
		}
	}
	/*for(int i = yMax; i >= yMin; i--)
	{
		cout<<i<<endl;
		for(int j = 0; j < scanLineCount[i]; j++)
		{
			cout<<sortedInter[j][i]<<" ";
		}
		cout<<scanLineCount[i]<<endl;
	}*/
	for(int i = yMax; i >= yMin; i--)
	{
		int count = 0;
		for (int j = 0; j < 501; j++)
		{
			if(realWorldCoor[j][i] == 0 && j > sortedInter[count][i] && j < sortedInter[count + 1][i])
			{
				realWorldCoor[j][i] = 1;
				if((j >= viewXLowerBound && j <= viewXUpperBound) && (i >= viewYLowerBound && i <= viewYUpperBound))
				{
					viewWorldCoor[j][i] = 1;
				}
			}
			else if(j == sortedInter[count + 1][i])
			{
				if(sortedInter[count][i] == sortedInter[count + 1][i])
				{
					j = j - 1;
				}
				if(scanLineCount[i] > (count + 2) || scanCount[i] > (count + 2))
				{					
					count = count + 2;
				}
				else
				{
					break;
				}
			}
		}
	}
}

void matrixMul(double a[4][4], double b[4][4])
{
	for(int i = 0;i < 4;i++) 
	{
		for(int j = 0;j < 4;j++)
		{ 
			for(int k = 0;k < 4;k++) 
			{
				tempMat[i][j] += a[i][k] * b[k][j];
			}
		}
	}
}

void matrixMulVer(double a[4][4], double b[4][1])
{
	tempMatVer[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0] + a[0][3] + b[3][0];
	tempMatVer[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0] + a[1][3] + b[3][0];
	tempMatVer[2][0] = a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0] + a[2][3] + b[3][0];
	tempMatVer[3][0] = a[3][0] * b[0][0] + a[3][1] * b[1][0] + a[3][2] * b[2][0] + a[3][3] + b[3][0];
}

void matrixSet(double a[4][4], double b[4][4])
{
	for(int i = 0;i < 4;i++) 
	{
		for(int j = 0;j < 4;j++)
		{ 
			a[i][j] = b[i][j];
		}
	}
}