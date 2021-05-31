from types import coroutine
import numpy as np
import cv2
import random
from numpy.core.defchararray import rpartition
from numpy.testing._private.utils import print_assert_equal
from tqdm import tnrange, tqdm
from utils import solve_homography, warping

random.seed(999)


def panorama(imgs):
    """
    Image stitching with estimated homograpy between consecutive
    :param imgs: list of images to be stitched
    :return: stitched panorama
    """
    h_max = max([x.shape[0] for x in imgs])
    w_max = sum([x.shape[1] for x in imgs])

    # create the final stitched canvas
    dst = np.zeros((h_max, w_max, imgs[0].shape[2]), dtype=np.uint8)
    dst[:imgs[0].shape[0], :imgs[0].shape[1]] = imgs[0]
    last_best_H = np.eye(3)
    out = None

    # for all images to be stitched:
    for idx in range(len(imgs) - 1):
        im1 = imgs[idx]
        im2 = imgs[idx + 1]

        # TODO: 1.feature detection & matching
        # Initialize the ORB detector algorithm
        orb = cv2.ORB_create()
        Keypoints1, Descriptors1 = orb.detectAndCompute(im1, None)
        Keypoints2, Descriptors2 = orb.detectAndCompute(im2, None)

        # Initialize the Matcher for matching the keypoints and then match the keypoints
        matcher = cv2.BFMatcher()
        matches = matcher.match(Descriptors1, Descriptors2)
        coor_kp1 = [Keypoints1[mat.queryIdx].pt for mat in matches]
        coor_kp2 = [Keypoints2[mat.trainIdx].pt for mat in matches]
        coor_kp1 = np.asarray(coor_kp1)
        coor_kp2 = np.asarray(coor_kp2)

        # TODO: 2. apply RANSAC to choose best H
        H = RANSAC(coor_kp2, coor_kp1)

        # TODO: 3. chain the homographies
        last_best_H = np.dot(last_best_H, H)

        # TODO: 4. apply warping
        h, w, c = im2.shape

        # img_temp = warping(im1, temp, H, ymin, ymax, xmin, xmax, 'b')
        out = warping(im2, dst, last_best_H, 0, h, 0, w*(idx+2), 'b')

    return out


def RANSAC(u, v):
    '''
    # sample point: s
    selected sample point: selected_s
    # of coordinate that is smaller than distance threshold: m
    # of all matches: M
    distance threshold: dis_theshold
    threshold of terminate: T
    # of trial: N
    '''
    s = 4
    dis_threshold = 5
    match_num = u.shape[0]
    thres_ratio = 0.5
    N = 1000
    best_H = np.zeros((3, 3))
    best_inlier_num = 0  # most inlier number

    for i in range(N):
        selected_s = np.random.randint(match_num, size=s)
        selected_u = u[selected_s, :]
        selected_v = v[selected_s, :]
        H = solve_homography(selected_u, selected_v)

        # convert coor_kp1 to homogeneous coordinate
        homo_u = np.ones((3, u.shape[0]))
        homo_u[0:-1, :] = u.T

        trans_u = np.dot(H, homo_u)

        np.seterr(divide='ignore', invalid='ignore')  # TODO
        # nearest neighbor
        trans_u = np.round((trans_u/trans_u[-1, :])[0:-1, :]).astype('int64')

        distance = np.sum((trans_u-v.T)**2, axis=0)
        inlier = distance[distance < dis_threshold]
        inlier_num = inlier.shape[0]

        # ----- check if H is singular -----
        try:
            H_inv = np.linalg.inv(H)
        except:
            # print("MATRIX IS SINGULAR")
            continue
        else:
            if ((inlier_num/match_num) >= thres_ratio and best_inlier_num > 0):
                break
            else:
                # ----- update best H -----
                if (inlier_num > best_inlier_num or best_inlier_num == 0):
                    best_inlier_num = inlier_num
                    best_H = H
                pass

    return best_H


if __name__ == "__main__":

    # ================== Part 4: Panorama ========================
    # TODO: change the number of frames to be stitched
    FRAME_NUM = 3
    imgs = [cv2.imread('../resource/frame{:d}.jpg'.format(x))
            for x in range(1, FRAME_NUM + 1)]

    output4 = panorama(imgs)
    cv2.imwrite('output4.png', output4)
