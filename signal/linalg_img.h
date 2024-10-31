#ifndef LINALG_IMG_H
#define LINALG_IMG_H

#include <opencv2/opencv.hpp>

namespace linalg {

static uint MIN = 0 ;
static uint AVG = 1 ;
static uint MAX = 2 ; 

cv::Mat compressGray(const cv::Mat&, float) ;
cv::Mat compressColor(const cv::Mat&, float, const uint) ;

}

#endif