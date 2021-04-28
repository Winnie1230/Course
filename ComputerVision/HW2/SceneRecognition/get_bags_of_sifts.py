from PIL import Image
import numpy as np
from scipy.spatial import distance
import pickle
import scipy.spatial.distance as distance
from cyvlfeat.sift.dsift import dsift
from time import time


def get_bags_of_sifts(image_paths, vocab_file):
    ############################################################################
    # TODO:                                                                    #
    # This function assumes that 'vocab.pkl' exists and contains an N x 128    #
    # matrix 'vocab' where each row is a kmeans centroid or visual word. This  #
    # matrix is saved to disk rather than passed in a parameter to avoid       #
    # recomputing the vocabulary every time at significant expense.            #
    #                                                                          #
    # image_feats is an N x d matrix, where d is the dimensionality of the     #
    # feature representation. In this case, d will equal the number of clusters#
    # or equivalently the number of entries in each image's histogram.         #
    #                                                                          #
    # You will construct SIFT features here in the same way you did in         #
    # build_vocabulary (except for possibly changing the sampling rate)        #
    # and then assign each local feature to its nearest cluster center         #
    # and build a histogram indicating how many times each cluster was used.   #
    # Don't forget to normalize the histogram, or else a larger image with more#
    # SIFT features will look very different from a smaller version of the same#
    # image.                                                                   #
    ############################################################################
    '''
    Input : 
        image_paths : a list(N) of training images
    Output : 
        image_feats : (N, d) feature, each row represent a feature of an image
    '''

    image_feats = []
    vocab_cen = pickle.load(open(vocab_file, 'rb'))
    for i in image_paths:
        img = np.asarray(Image.open(i), dtype='float32')
        frames, descriptors = dsift(img, step=[1, 1], fast=True)
        distance_matrix = distance.cdist(descriptors, vocab_cen, 'euclidean')
        index = np.argmin(distance_matrix, axis=1)
        values, counts = np.unique(index, return_counts=True)
        counter = dict(zip(values, counts))

        histogram = np.zeros(vocab_cen.shape[0])
        for idx, count in counter.items():
            histogram[idx] = count

        # normalize histogram
        histogram = histogram / histogram.sum()
        image_feats.append(histogram)

    image_feats = np.asarray(image_feats)

    #############################################################################
    #                                END OF YOUR CODE                           #
    #############################################################################
    return image_feats
