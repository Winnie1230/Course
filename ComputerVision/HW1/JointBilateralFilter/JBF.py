import numpy as np
import cv2


class Joint_bilateral_filter(object):
    def __init__(self, sigma_s, sigma_r):
        self.sigma_r = sigma_r  # range sigma
        self.sigma_s = sigma_s  # spatial sigma
        self.wndw_size = 6*sigma_s+1
        self.pad_w = 3*sigma_s

    def getGaussianKernel(self):
        x, y = np.mgrid[-(self.wndw_size - self.wndw_size//2-1):(self.wndw_size - self.wndw_size // 2),
                        -(self.wndw_size - self.wndw_size//2-1):(self.wndw_size - self.wndw_size//2)]
        gaussian_kernel = np.exp(-(x**2+y**2)/(2*self.sigma_s**2))
        gaussian_kernel = gaussian_kernel / gaussian_kernel.sum()  # normalize
        return gaussian_kernel

    def getRangeKernel(self, guidance):
        # normalize pixel value to [0,1] to construct range kernel
        if (guidance.ndim == 2):  # one channel
            sort = self.SortConvImg(guidance[:, :])
            sort = sort / 255
            diff = (sort - sort[:, sort.shape[1]//2].reshape(-1, 1))**2
        else:  # three channel
            img_h, img_w, ch = guidance.shape
            diff = np.zeros(
                ((img_h-2*self.pad_w)*(img_w-2*self.pad_w), (self.wndw_size**2)))
            for i in range(ch):
                sort = self.SortConvImg(guidance[:, :, i])
                sort = sort / 255
                diff += (sort - sort[:, sort.shape[1]//2].reshape(-1, 1))**2

        range_kernel = np.exp(-(diff) / (2*self.sigma_r**2))
        return range_kernel

    def SortConvImg(self, img):
        img_h, img_w = img.shape
        # first index of conv_img
        i0 = np.repeat(np.arange(img_h-2*self.pad_w), img_w-2*self.pad_w)
        i1 = np.repeat(np.arange(self.wndw_size), self.wndw_size)
        i = i0.reshape(-1, 1) + i1.reshape(1, -1)

        # second index of conv_img
        j0 = np.tile(np.arange(img_w-2*self.pad_w), img_h-2*self.pad_w)
        j1 = np.tile(np.arange(self.wndw_size), self.wndw_size)
        j = j0.reshape(-1, 1)+j1.reshape(1, -1)

        sorted_img = img[i, j]
        return sorted_img

    def joint_bilateral_filter(self, img, guidance):
        BORDER_TYPE = cv2.BORDER_REFLECT
        padded_img = cv2.copyMakeBorder(
            img, self.pad_w, self.pad_w, self.pad_w, self.pad_w, BORDER_TYPE)
        padded_guidance = cv2.copyMakeBorder(
            guidance, self.pad_w, self.pad_w, self.pad_w, self.pad_w, BORDER_TYPE)

        ### TODO ###
        # Step1. get spatial kernel(Gaussian kernel)
        spatial_kernel = self.getGaussianKernel().flatten()

        # Step2. Sort convolution matrix
        # input img
        sorted_r = self.SortConvImg(padded_img[:, :, 0])
        sorted_g = self.SortConvImg(padded_img[:, :, 1])
        sorted_b = self.SortConvImg(padded_img[:, :, 2])

        # Step3. get Range kernel(use guidance img as reference)
        range_kernel = self.getRangeKernel(padded_guidance)

        # Step4. Convolution
        kernel = spatial_kernel*range_kernel
        # normalize kernel (let kernel sum = 1)
        kernel = kernel / (kernel.sum(axis=1).reshape(-1, 1))

        img_h, img_w, ch = img.shape  # origin img size
        jbl_r = (sorted_r*kernel).sum(axis=1).reshape(img_h, img_w)
        jbl_g = (sorted_g*kernel).sum(axis=1).reshape(img_h, img_w)
        jbl_b = (sorted_b*kernel).sum(axis=1).reshape(img_h, img_w)

        output = np.zeros((img.shape), dtype=np.uint8)
        output[:, :, 0] = jbl_r
        output[:, :, 1] = jbl_g
        output[:, :, 2] = jbl_b

        return np.clip(output, 0, 255).astype(np.uint8)
