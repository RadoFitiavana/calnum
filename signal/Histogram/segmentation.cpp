/***************************************
Build:
g++ -Wall -o "%e" "%f" -IC:\opencv\build\include -LC:\opencv\build\x86\mingw\lib -lopencv_core310 -lopencv_imgcodecs310 -lopencv_highgui310 -lopencv_imgproc310
***************************************/
#include <iostream>
#include <fstream>
#include <cmath>
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"

using namespace cv;
using namespace std;
void histogram(Mat img, Mat histImage, int hist_w, int hist_h);
void threshold1(int &estim, int &m1, int &m2, Mat img);

int main(){
    cout << "Segmentation" << endl;
    Mat img;

/// Get the original image
    img = imread("docphys.jpg", IMREAD_GRAYSCALE);
    imshow("Original", img);

/// Compute histogram
    int bins=256, histo[bins];
    for(int i=0; i<bins; i++) histo[i] = 0;
    for(int i=0; i<img.rows*img.cols; i++){
        histo[img.data[i]]++;
    }
/// Plot histogram with gnuplot
    FILE *GP = popen("c:\\maxima-5.46.0\\gnuplot\\bin\\gnuplot -persist","w");
    if (GP) {
		fprintf(GP, "set term wxt size 400,300\n");
        fprintf(GP, "set term wxt title 'Histogram'\n");
		fprintf(GP, "set xrange [0:255]\n");
		fprintf(GP, "set nokey\n");    // No legend
        fprintf(GP, "$histo << EOD\n");
		for(int i=0; i<bins; i++){
			fprintf(GP,"%d %d\n", i, histo[i]);
		}
		fprintf(GP, "EOD\n");
		fprintf(GP, "plot $histo with lines\n");
		fflush(GP);
	}						// pipe closing delayed till end of program
	else cout << "gnuplot not found..." << endl;

/// Iterated threshold calculation for bimodal image
	int thres1(200), thres2(255), m1(0), m2(0), xtop(0), vtop(0);

/// Find top peak of histogram
	for(int i=0; i<bins; i++){
		if(histo[i]>vtop){
			xtop = i;
			vtop = histo[i];
		}
	}
	thres1 = xtop;
	cout << "\nInit thres = " << xtop << endl;
	while(abs(thres1-thres2)>0){
		thres2 = thres1;
		threshold1(thres1, m1, m2, img);
		cout << "\nThreshold = " << thres1 << endl;
		cout << "Means of the 2 groups: " << m1 << " : " << m2 << endl;
	}

/// Threshold segmentation
	Mat imgout = img;

	for(int i=0; i<img.rows*img.cols; i++){
		if(imgout.data[i] < thres1)
			imgout.data[i]= 0;
		else
			imgout.data[i]=255;
	}


	threshold( img, imgout, thres1, 255,0 );
	imshow("Segmented", imgout);

	Mat imgout1 = img;
	GaussianBlur( imgout, imgout1, Size(3,3), 0, 0, BORDER_DEFAULT);
	imshow("Smoothed", imgout1);

	waitKey();

	if(GP)	pclose(GP);

    return 0;
}

void threshold1(int &estim, int &m1, int &m2, Mat img){
	long t0(estim), t1(255);	// assume the other pick is to the left of top pick
    long i1(0), i2(0);
    while((t0 - t1)< 2){
        m1 = 0; m2 = 0;
        i1 = 1; i2 = 1;
        for (int i=0; i< img.rows*img.cols; i++){
            if(img.data[i] < t0){
                m1 += img.data[i];
                i1++;
            }
            else{
                m2 += img.data[i];
                i2++;
            }
        }
        m1 /= i1;
        m2 /= i2;
        t0 = t1;
        t1 = (m1+m2)/2;
    }
    estim = t1;
}

void histogram(Mat img, Mat histImage, int hist_w, int hist_h){
    Mat hist;
    int channels[]={0};
    int histSize[] ={hist_w};  /// Number of bins
    float range[] = { 0, 256 };
    const float* ranges[] = { range };

    int bin_w  = cvRound( (double) hist_w/histSize[0] );
    calcHist( &img, 1, channels, Mat(),
              hist, 1, histSize, ranges,
              true, false );

/// Normalize the result to [ 0, histImage.rows ]
    normalize(hist, hist, 0, histImage.rows, NORM_MINMAX, -1, Mat() );
    for( int i = 1; i < histSize[0]; i++ )
        rectangle(histImage, Point( bin_w*(i-1), hist_h ) ,
                  Point( bin_w*(i), hist_h - cvRound(hist.at<float>(i)) ) ,
                  Scalar( 255, 0, 0), -1, 8, 0);
}
