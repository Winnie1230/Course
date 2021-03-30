import numpy as np
import cv2
import argparse
from matplotlib import pyplot as plt
from HCD import Harris_corner_detector


def main():
    parser = argparse.ArgumentParser(
        description='main function of Harris corner detector')
    parser.add_argument('--threshold', default=100., type=float,
                        help='threshold value to determine corner')
    parser.add_argument(
        '--image_path', default='./testdata/1.png', help='path to input image')
    args = parser.parse_args()

    print('Processing %s ...' % args.image_path)
    img = cv2.imread(args.image_path)
    img_gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY).astype(np.float64)

    ### TODO ###
    # create HCD class
    HCD = Harris_corner_detector(args.threshold)
    response = HCD.detect_harris_corners(img_gray)
    result = HCD.post_processing(response)
    print(len(result))

    detected_img = np.array(img)
    for i in result:
        cv2.circle(detected_img, (i[1], i[0]), radius=1,
                   color=(0, 0, 255), thickness=1)

    plt.subplot(121), plt.imshow(cv2.cvtColor(
        img, cv2.COLOR_BGR2RGB)), plt.title('Original')
    plt.xticks([]), plt.yticks([])
    plt.subplot(122), plt.imshow(cv2.cvtColor(
        detected_img, cv2.COLOR_BGR2RGB)), plt.title('Harris Corner Detection')
    plt.xticks([]), plt.yticks([])
    plt.show()

    '''
    plt.imshow(cv2.cvtColor(detected_img, cv2.COLOR_BGR2RGB)), plt.title(
        'Harris Corner Detection (threshold = '+str(args.threshold)+')')
    plt.xticks([]), plt.yticks([])
    plt.show()
    '''


if __name__ == '__main__':
    main()
