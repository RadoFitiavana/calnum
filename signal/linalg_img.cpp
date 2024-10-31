#include "linalg_img.h"
#include <opencv2/opencv.hpp>
#include <vector>
#include <iostream>
#include <cmath>

using namespace cv ;
using namespace std ;

namespace linalg {

Mat compressGray(const Mat& img, float ratio) {
    Mat w, u, vt, diag ;
    vector<float> sv ;
    SVD::compute(img,w,u,vt) ;
    sv.assign((float*)w.data, (float*)w.data + w.rows) ;
    float T = 0 ;
    for (uint i=0; i<sv.size(); i++) {
        T += abs(sv[i]) ;
    }
    uint n = 1 ;
    float r = abs(sv[0])/T ;
    for (uint i=0; i<sv.size(); i++){
        if (r >= ratio){
            break ;
        }
        r += sv[i+1]/T ;
        n += 1 ;
    }
    sv.clear() ;
    sv.assign((float*)w.data, (float*)w.data + n);
    diag = Mat::zeros(n, n, CV_32F) ;
    for (uint i = 0; i < n; i++) {
        diag.at<float>(i, i) = sv[i];
    }

    cout << n << endl ;
    return u.colRange(0, n) * diag * vt.rowRange(0, n) ;
}

}