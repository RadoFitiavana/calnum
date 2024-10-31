#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>

using namespace std ;
using namespace cv ;

int main() {
    cout << "Hello image" << endl ;
    Mat img ;
    vector<int> compression ;

    // Read and display 
    img = imread("lena.png",IMREAD_COLOR) ; // IMREAD_GRAYSCALE
    imshow("Lena",img) ;

    // Save: the saving format is according to the extension of the file
    compression.push_back(IMWRITE_JPEG_QUALITY) ; // 1
    compression.push_back(100) ; //default is 95
    imwrite("lena1.png",img,compression) ;
    waitKey() ;


    return 0 ;
}