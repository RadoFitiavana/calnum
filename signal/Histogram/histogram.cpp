/**********************************
Build:
g++ -Wall -o "%e" "%f" -IC:\opencv\include -LC:\opencv\mybuild\lib -lopencv_core460 -lopencv_imgcodecs460 -lopencv_highgui460 -lopencv_imgproc460

***********************************/
#include <iostream>
#include <fstream>
#include <cmath>
#include "opencv2/opencv.hpp"

using namespace cv;
using namespace std;
void histogram(Mat img, Mat histImage, int hist_w, int hist_h);
void histogramc(Mat img, Mat histImage, int hist_w, int hist_h);
void CallBackFunc(int event, int x, int y, int flags, void* userdata);
void gpHistogram(Mat img);

int main(){
    cout << "Histogramme(s)" << endl;
    Mat img;
/// Get an image
    img = imread("empreinte.png",IMREAD_COLOR);
//	img = imread("lena-soderberg.png",0);
    imshow("Lena", img);

/// Compute histogram with openCV
    int hist_w(256), hist_h(256);
    Mat histImage ( hist_h, hist_w, CV_8UC3, Scalar( 255,255,255) );
    histogramc(img, histImage, hist_w, hist_h);		// image en couleur
//    histogram(img, histImage, hist_w, hist_h);			// image en niveaux de gris

/// Display histogram
	namedWindow("Histogram", 1); // pour identifier l'image d'histogramme
	setMouseCallback("Histogram", CallBackFunc, NULL);
    imshow("Histogram", histImage );
    waitKey();

/// Display histogram with gnuplot 5.x

	vector<Mat> bgrPlanes;
    split(img, bgrPlanes);
	gpHistogram(bgrPlanes[1]);

//	gpHistogram(img);

    return 0;
}

void histogram(Mat img, Mat histImage, int hist_w, int hist_h){
    Mat hist; 							// to store the histogram(s)
    int histSize[] ={hist_w};   		// Number of bins
    float range[] = { 0, 256 }; 		// Range of values to be measured
    const float* ranges[] = { range };

    int bin_w  = cvRound( (double) hist_w/histSize[0] );
    calcHist( &img, 1, 0, Mat(),
              hist, 1, histSize, ranges,
              true, false );

/// Normalize the result to [ 0, histImage.rows ]
    normalize(hist, hist, 0, histImage.rows, NORM_MINMAX, -1, Mat() );
    for( int i = 1; i < histSize[0]; i++ )
        rectangle(histImage, Point( bin_w*(i-1), hist_h ) ,
                  Point( bin_w*(i), hist_h - cvRound(hist.at<float>(i)) ) ,
                  Scalar( 255, 0, 0), -1, 8, 0);
}

void histogramc(Mat img, Mat histImage, int hist_w, int hist_h){
    Mat bhist, ghist, rhist; 			// to store the histograms
    int histSize[] = {hist_w};   		// Number of bins
    float range[] = { 0, 256 }; 		// Range of values to be measured
    const float* ranges[] = { range };

    vector<Mat> bgrPlanes;
    split(img, bgrPlanes);

    int bin_w  = cvRound( (double) hist_w/histSize[0] );
    calcHist( &bgrPlanes[0], 1, 0, Mat(),
              bhist, 1, histSize, ranges,
              true, false );
    calcHist( &bgrPlanes[1], 1, 0, Mat(),
              ghist, 1, histSize, ranges,
              true, false );
    calcHist( &bgrPlanes[2], 1, 0, Mat(),
              rhist, 1, histSize, ranges,
              true, false );

/// Normalize the result to [ 0, histImage.rows ]
    normalize(bhist, bhist, 0, histImage.rows, NORM_MINMAX, -1, Mat() );
    normalize(ghist, ghist, 0, histImage.rows, NORM_MINMAX, -1, Mat() );
    normalize(rhist, rhist, 0, histImage.rows, NORM_MINMAX, -1, Mat() );
    for( int i = 1; i < histSize[0]; i++ ){
        line(histImage, Point( bin_w*(i-1), hist_h - cvRound(bhist.at<float>(i-1))) ,
                  Point( bin_w*(i), hist_h - cvRound(bhist.at<float>(i)) ) ,
                  Scalar( 0, 0, 255), 2, 8, 0);
        line(histImage, Point( bin_w*(i-1), hist_h - cvRound(ghist.at<float>(i-1))) ,
                  Point( bin_w*(i), hist_h - cvRound(ghist.at<float>(i)) ) ,
                  Scalar( 0, 255, 0), 2, 8, 0);
        line(histImage, Point( bin_w*(i-1), hist_h - cvRound(rhist.at<float>(i-1))) ,
                  Point( bin_w*(i), hist_h - cvRound(rhist.at<float>(i)) ) ,
                  Scalar( 255, 0, 0), 2, 8, 0);
    }
}

void gpHistogram(Mat img){
    int bins=256, histo[bins];
    for(int i=0; i<bins; i++) histo[i] = 0;
    for(int i=0; i<img.rows*img.cols; i++){
        histo[img.data[i]]++;
    }
    FILE *gnuplotPipe = popen("c:\\maxima-5.47.0\\gnuplot\\bin\\gnuplot -persist","w");
    cout << "Display histogram using gnuplot: press q to continue..." << endl;
    if (gnuplotPipe) {  // If gnuplot is found
        fprintf(gnuplotPipe, "set title \"Histogram\"\n");
        fprintf(gnuplotPipe, "set terminal wxt size 500, 500\n");
        fprintf(gnuplotPipe, "set style data boxes\n");
		fprintf(gnuplotPipe, "set xrange [0:255]\n");
		fprintf(gnuplotPipe, "set nokey\n");    // No legend
        fprintf(gnuplotPipe, "$toto << EOD\n");
		for(int i=0; i<bins; i++){
			fprintf(gnuplotPipe,"%d %d\n", i, histo[i]);
		}
		fprintf(gnuplotPipe, "EOD\n");
		fprintf(gnuplotPipe, "plot $toto with lines\n");
		fflush(gnuplotPipe);
		pclose(gnuplotPipe);
	}
	else cout << "gnuplot not found..." << endl;
}

void CallBackFunc(int event, int x, int y, int flags, void* userdata){
	switch(event){
		case(EVENT_LBUTTONDOWN):
			cout << "(" << x << ", " << y << ")" << endl;
			break;
		case(EVENT_RBUTTONDOWN):
			cout << "(" << x << ", " << y << ")" << endl;
			break;
		case(EVENT_MBUTTONDOWN):
			cout << "(" << x << ", " << y << ")" << endl;
			break;
		case(EVENT_MOUSEMOVE):			// lecture de l'histogramme
			cout << "(" << x << ", " << y << ")" << endl;
			break;
		default:
			break;
    }
}
