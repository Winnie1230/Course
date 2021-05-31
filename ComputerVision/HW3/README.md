# Homework3: Projective Geometry
## Part1. Homography Estimation
Use DLT(Direct Linear Transformation) algorithm to calculate homography between two images.

**Objective**:
Given n>=4 2D to 2D point correspondence xi <=> xi' determine the 2D homography matrix H such that xi'=Hxi.

**Algorithm**:
Step1. For each correspondence compute Ai, and pick first two columns of Ai. 
Step2. Assemble n 2x9 matrices Ai into a single 2nx9 matrix A. 
Step3. Obtain SVD of A. Solution for h is the **last column of V**.
Step4. Determine H from h.(reshape 1x9 to 3x3)

## Part2. Marker-Based Planar AR
### Warping
**Forward warping**: Given homography between source image and destination image, transform source pixel coordinate directly into destination coordinate. 
One can use nearest neighbor or bilinear interpolation to get the destination pixel coordinate.
> Problem: 
> There may be holes in the destination image since possibly not every destination locations is mapped(ex. the destination region is larger than source image etc.) 

**Backward warping**: Given homography between source image and destination image, inverse transform destination coordinate to source to get the corresponding pixel value. 
There should be not holes since every destination locations are mapped. 
> the transformation must not be singular(i.e. the inverse should exist)

### ArUco
[Detecting ArUco markers with OpenCV and Python](https://www.pyimagesearch.com/2020/12/21/detecting-aruco-markers-with-opencv-and-python/)
Use `aruco.detectMarkers()` to get the detected corner

## Part3. Unwarp the Secret
Unwarp QR code and see if one can get the correct QR code link from both images.
Use backward warping.

## Part4. Panorama
Step1. Feature detection and matching
> Use opencv built-in ORB detector for keypoint detection and opencv brute force matcher for feature matching. 

Step2. Apply RANSAC to choose best H 
Step3. Chain the homographies
> In order to map frame2 and frame3 to frame1, you determine the homography between frame1, frame2(H12) and frame2, frame3(H23) and the homography between frame1 and frame3 is H12xH23. 

Step4. Apply warping


