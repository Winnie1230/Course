import numpy as np
import cv2
import matplotlib.pyplot as plt


class Harris_corner_detector(object):
    def __init__(self, threshold):
        self.threshold = threshold

    def detect_harris_corners(self, img):
        ### TODO ####
        # Step 1: Smooth the image by Gaussian kernel
        # - Function: cv2.GaussianBlur (kernel = 3, sigma = 1.5)
        img_blur = cv2.GaussianBlur(img, (3, 3), 1.5)

        # Step 2: Calculate Ix, Iy (1st derivative of image along x and y axis)
        # - Function: cv2.filter2D (kernel = [[1.,0.,-1.]] for Ix or [[1.],[0.],[-1.]] for Iy)
        kernel_x = np.array(([[1., 0., -1.]]), dtype="float64")
        kernel_y = np.array(([1.], [0.], [-1.]), dtype="float64")

        Ix = cv2.filter2D(img_blur, -1, kernel_x)
        Iy = cv2.filter2D(img_blur, -1, kernel_y)

        # Step 3: Compute Ixx, Ixy, Iyy (Ixx = Ix*Ix, ...)
        Ixx = Ix*Ix
        Ixy = Ix*Iy
        Iyy = Iy*Iy

        # Step 4: Compute Sxx, Sxy, Syy (weighted summation of Ixx, Ixy, Iyy in neighbor pixels)
        # - Function: cv2.GaussianBlur (kernel = 3, sigma = 1.)
        Sxx = cv2.GaussianBlur(Ixx, (3, 3), 1.).flatten()
        Sxy = cv2.GaussianBlur(Ixy, (3, 3), 1.).flatten()
        Syy = cv2.GaussianBlur(Iyy, (3, 3), 1.).flatten()

        # Step 5: Compute the det and trace of matrix M (M = [[Sxx, Sxy], [Sxy, Syy]])
        det_M = Sxx*Syy - Sxy*Sxy
        trace_M = Sxx + Syy

        # Step 6: Compute the response of the detector by det/(trace+1e-12)
        response = (det_M/(trace_M + 1e-12)).reshape(img.shape)

        return response

    def post_processing(self, response):
        ### TODO ###
        # Step 1: Thresholding
        response[response < self.threshold] = 0
        candidate_list = np.array(np.argwhere(
            response > self.threshold).tolist())

        # Step 2: Find local maximum
        # ----- zero padding -----
        response_padding = cv2.copyMakeBorder(
            response, 2, 2, 2, 2, cv2.BORDER_CONSTANT, value=0)

        window_size = 5
        # check if candidate point is window's local maximum
        for index in candidate_list:
            sgn_value = np.sign(
                response_padding[index[0]:index[0]+window_size, index[1]:index[1]+window_size] - response_padding[index[0]+window_size//2, index[1]+window_size//2])

            if (1 in sgn_value):
                response[index[0], index[1]] = 0

        local_max = np.argwhere(response > 0).tolist()

        return local_max
