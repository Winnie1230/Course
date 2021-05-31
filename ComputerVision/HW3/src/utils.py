import numpy as np


def solve_homography(u, v):
    """
    This function should return a 3-by-3 homography matrix,
    u, v are N-by-2 matrices, representing N corresponding points for v = T(u)
    :param u: N-by-2 source pixel location matrices
    :param v: N-by-2 destination pixel location matrices
    :return:
    """
    N = u.shape[0]
    H = None

    if v.shape[0] is not N:
        print('u and v should have the same size')
        return None
    if N < 4:
        print('At least 4 points should be given')

    homo_u = np.ones((u.shape[0], u.shape[1]+1))
    homo_v = np.ones((v.shape[0], v.shape[1]+1))
    homo_u[:, :-1] = u
    homo_v[:, :-1] = v
    A = np.zeros((u.shape[0]*2, 9))

    # TODO: 1.forming A
    for i in range(u.shape[0]):
        A[i*2, :] = np.concatenate((np.zeros((3)),
                                    np.dot(-homo_v[i][2], homo_u[i]), np.dot(homo_v[i][1], homo_u[i])))
        A[i*2+1, :] = np.concatenate((np.dot(homo_v[i][2], homo_u[i]),
                                      np.zeros(3), np.dot(-homo_v[i][0], homo_u[i])))

    # TODO: 2.solve H with A
    u, s, vh = np.linalg.svd(A, full_matrices=True)
    H = (vh[-1, :]).reshape(3, 3)
    # nomalize
    H = H/H.sum()
    return H


def warping(src, dst, H, ymin, ymax, xmin, xmax, direction='b'):
    """
    Perform forward/backward warpping without for loops. i.e.
    for all pixels in src(xmin~xmax, ymin~ymax),  warp to destination
          (xmin=0,ymin=0)  source                       destination
                         |--------|              |------------------------|
                         |        |              |                        |
                         |        |     warp     |                        |
    forward warp         |        |  --------->  |                        |
                         |        |              |                        |
                         |--------|              |------------------------|
                                 (xmax=w,ymax=h)

    for all pixels in dst(xmin~xmax, ymin~ymax),  sample from source
                            source                       destination
                         |--------|              |------------------------|
                         |        |              | (xmin,ymin)            |
                         |        |     warp     |           |--|         |
    backward warp        |        |  <---------  |           |__|         |
                         |        |              |             (xmax,ymax)|
                         |--------|              |------------------------|

    :param src: source image
    :param dst: destination output image
    :param H:
    :param ymin: lower vertical bound of the destination(source, if forward warp) pixel coordinate
    :param ymax: upper vertical bound of the destination(source, if forward warp) pixel coordinate
    :param xmin: lower horizontal bound of the destination(source, if forward warp) pixel coordinate
    :param xmax: upper horizontal bound of the destination(source, if forward warp) pixel coordinate
    :param direction: indicates backward warping or forward warping
    :return: destination output image
    """

    h_src, w_src, ch = src.shape
    h_dst, w_dst, ch = dst.shape
    H_inv = np.linalg.inv(H)

    # TODO: 1.meshgrid the (x,y) coordinate pairs
    x = np.tile(np.arange(xmin, xmax, 1), ymax-ymin)
    y = np.repeat(np.arange(ymin, ymax, 1), xmax-xmin)

    # TODO: 2.reshape the destination pixels as N x 3 homogeneous coordinate
    u = np.ones((3, x.shape[0])).astype('int64')
    u[0, :] = x
    u[1, :] = y

    if direction == 'b':
        # TODO: 3.apply H_inv to the destination pixels and retrieve (u,v) pixels, then reshape to (ymax-ymin),(xmax-xmin)
        # u: destination pixel coordinate(in homogeneous coordinate)
        # v: source pixel coordinate(in homogeneous coordinate)
        H_inv = np.linalg.inv(H)
        v = np.dot(H_inv, u)
        v = np.round((v/v[-1, :])[0:-1, :]).astype('int64')  # nearest neighbor

        # TODO: 4.calculate the mask of the transformed coordinate (should not exceed the boundaries of source image)
        mask_y = np.logical_and(v[0, :] < w_src, v[0, :] >= 0)
        mask_x = np.logical_and(v[1, :] < h_src, v[1, :] >= 0)
        mask = np.logical_and(mask_x, mask_y)
        mask = np.tile(mask, 2)
        u = u[0:2, :]  # convert homogeneous coordinate to cartesian coordinate
        v = (v.flatten())[mask]
        u = (u.flatten())[mask]
        v = v.reshape(2, -1)
        u = u.reshape(2, -1)

        # TODO: 5.sample the source image with the masked and reshaped transformed coordinates

        # TODO: 6. assign to destination image with proper masking
        dst[u[1, :], u[0, :]] = src[v[1, :], v[0, :]]
        pass

    elif direction == 'f':
        # TODO: 3.apply H to the source pixels and retrieve (u,v) pixels, then reshape to (ymax-ymin),(xmax-xmin)
        # u: source pixel coordinate(in homogeneous coordinate)
        # v: destination pixel coordinate(homogeneous coordinate)
        #u = np.ones((3, x.shape[0])).astype('int64')
        #u[0, :] = x
        #u[1, :] = y
        #u[0, :] = np.tile(np.arange(w_src), h_src)
        #u[1, :] = np.repeat(np.arange(h_src), w_src)
        # print(u.shape)
        v = np.dot(H, u)
        v = np.round((v/v[-1, :])[0:-1, :]).astype('int64')

        # TODO: 4.calculate the mask of the transformed coordinate (should not exceed the boundaries of destination image)
        mask_y = np.logical_and(v[0, :] < w_dst, v[0, :] >= 0)
        mask_x = np.logical_and(v[1, :] < h_dst, v[1, :] >= 0)
        mask = np.logical_and(mask_x, mask_y)
        mask = np.tile(mask, 2)

        # TODO: 5.filter the valid coordinates using previous obtained mask
        u = u[0:2, :]  # convert homogeneous coordinate to cartesian coordinate
        v = (v.flatten())[mask]
        u = (u.flatten())[mask]
        v = v.reshape(2, -1)
        u = u.reshape(2, -1)

        # TODO: 6. assign to destination image using advanced array indicing
        #print(v[0, 0], v[1, 0])
        #print(u[0, 0], u[1, 0])

        '''
        canvas_corners1 = np.array(
            [[749, 521], [883, 525], [883, 750], [750, 750]])
        for i in canvas_corners1:
            cv2.circle(dst, (i[0], i[1]), radius=1,
                       color=(0, 0, 255), thickness=10)
        plt.imshow(cv2.cvtColor(dst, cv2.COLOR_BGR2RGB))
        plt.show()
        '''
        '''
        print(v.shape)
        for i in range(v.shape[1]):
            #print(v[0, i], v[1, i])
            cv2.circle(dst, (v[0, i], v[1, i]), radius=1,
                       color=(0, 0, 255), thickness=1)
        '''
        '''
        tmp = np.zeros((dst.shape))
        tmp[v[0, :], v[1, :]] = src[u[0, :], u[1, :]]
        plt.imshow(cv2.cvtColor(tmp.astype(np.uint8), cv2.COLOR_BGR2RGB))
        plt.show()
        '''
        dst[v[1, :], v[0, :]] = src[u[1, :], u[0, :]]
        pass

    return dst
