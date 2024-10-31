#include <iostream>
#include "linalg_img.h"
#include <opencv2/opencv.hpp>

using namespace std ;
using namespace linalg ;
using namespace cv ;

int main(){
    Mat img = imread("lena.png",IMREAD_GRAYSCALE) ;
    img.convertTo(img, CV_32F) ;
    float ratio ;
    cout << "ratio: " ; cin >> ratio ;
    Mat imgC = compressGray(img,ratio) ;
    imgC.convertTo(imgC,CV_8U) ;
    img.convertTo(img,CV_8U) ;
    imshow("lena",img) ;
    imshow("lena compressed",imgC) ;
    waitKey() ;

    return 0 ;
}