# Bag-of-words Scene Recognition
Use two methods to implement Scene Recognition.
1. Tiny images representation and nearest neighbor classifier
2. Bag of SIFT representation and nearest neighbor classifier

Model Performance: 
Accuracy should be above the baseline score to get points
- Tiny images + kNN: 0.2
- Bag of SIFT + kNN : 0.55 (soft baseline)

### Tiny image
Resize image to 16x16 image and flatten it (1x256 dimension), and use this as the feature of the image.

### Bag of SIFT
Use SIFT algorithm to get the feature of the image. To get the baseline score, I have tried many combinations of parameters. The result is record in bag_of_sift.txt
The best accuracy by using different combination of parameters can only get to 0.564.
changable parameters:
1. kNN's k
2. vocabulary size
3. build vocabulary dsift step
4. get_bag_of_sift dsiftt step

### Reference
[Project 3: Scene recognition with bag of words](https://github.com/hschao/cv-hw4)

