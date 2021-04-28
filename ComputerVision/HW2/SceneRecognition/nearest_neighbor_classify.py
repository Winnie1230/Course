from __future__ import print_function

import numpy as np
import scipy.spatial.distance as distance
from scipy.stats import mode


def nearest_neighbor_classify(train_image_feats, train_labels, test_image_feats):
    ###########################################################################
    # TODO:                                                                   #
    # This function will predict the category for every test image by finding #
    # the training image with most similar features. Instead of 1 nearest     #
    # neighbor, you can vote based on k nearest neighbors which will increase #
    # performance (although you need to pick a reasonable value for k).       #
    ###########################################################################
    ###########################################################################
    # NOTE: Some useful functions                                             #
    # distance.cdist :                                                        #
    #   This function will calculate the distance between two list of features#
    #       e.g. distance.cdist(? ?)                                          #
    ###########################################################################
    '''
    Input : 
        train_image_feats : 
            image_feats is an (N, d) matrix, where d is the 
            dimensionality of the feature representation.

        train_labels : 
            image_feats is a list of string, each string
            indicate the ground truth category for each training image. 

        test_image_feats : 
            image_feats is an (M, d) matrix, where d is the 
            dimensionality of the feature representation.
    Output :
        test_predicts : 
            a list(M) of string, each string indicate the predict
            category for each testing image.
    '''

    test_predicts = []
    train_labels = np.asarray(train_labels)
    distance_matrix = distance.cdist(test_image_feats, train_image_feats)
    k = 6
    kmin_index = np.argsort(distance_matrix, axis=1)[:, :k]
    for i in range(kmin_index.shape[0]):
        # get the label with a majority vote
        index2label = [train_labels[j] for j in kmin_index[i, :]]
        maxlabel = max(index2label, key=index2label.count)
        test_predicts.append(maxlabel)

    #############################################################################
    #                                END OF YOUR CODE                           #
    #############################################################################
    return test_predicts
