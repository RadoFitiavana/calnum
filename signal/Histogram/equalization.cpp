/***************************************
Build:
g++ -Wall -o "%e" "%f"  -IC:\opencv\include -LC:\opencv\mybuild\lib -lopencv_core400 -lopencv_imgcodecs400 -lopencv_highgui400 -lopencv_imgproc400
***************************************/
#include <opencv2/opencv.hpp>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
using namespace cv;

void histogram(Mat img, Mat histImage, int hist_w, int hist_h);

int main(){
    cout << "Equalisation et filtrage..." << endl;
    Mat img;
/// Get an image
    img = imread("empreinte.png", 0);
    imshow("Original", img);

/// Compute histogram with openCV
    int hist_w(256), hist_h(256);
    Mat histImage ( hist_h, hist_w, CV_8UC3, Scalar( 255,255,255) );
    histogram(img, histImage, hist_w, hist_h);

/// Display histogram
    imshow("Histogram", histImage );

/// Histogram equalization by openCV
    Mat imgEq = img.clone();
    Mat histEq ( hist_h, hist_w, CV_8UC3, Scalar( 255,255,255) );
    equalizeHist( img, imgEq );
    histogram(imgEq, histEq, hist_w, hist_h);
    imshow("openCV Equalized histogram", histEq);
    imshow("openCV Equalized", imgEq);
    waitKey();

/// Compute histogram from definition
    Mat imgEqh = img.clone();
    Mat histEqh ( hist_h, hist_w, CV_8UC3, Scalar( 255,255,255) );
    int bins=256, histo[bins];
    for(int i=0; i<bins; i++) histo[i] = 0;
    for(int i=0; i<img.rows*img.cols; i++){
        histo[img.data[i]]++;
    }

/// Histogram equalization
    int N=img.rows*img.cols, i=0, k=0;
    float pk[bins], s=0;
    for(k=0; k<bins; k++)
        pk[k] = float(histo[k])/N;
    for(i=0; i< N; i++){
        for(k=0, s=0; k< img.data[i]; k++)
            s += pk[k];
        imgEqh.data[i] = floor(bins*s);
    }

    histogram(imgEqh, histEqh, hist_w, hist_h);
    imshow("Formula Equalized histogram", histEqh);
    imshow("Formula Equalized", imgEqh);
    waitKey();
    return 0;
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
