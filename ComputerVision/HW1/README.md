# Joint Bilateral Filter
> The difference between Joint Bilateral Filter and Bilateral Filter is the range kernel. Joint Bilateral Filter uses guidance image as reference to generate range kernel 
while Bilateral Filter uses the same image to construct range kernel.

> Project Constrains: Cannot use cv2.filter2D, cv2.GaussianBlur in implementing Joint Bilateral Filter

## Implement Steps (JBL.py)
Step1. Add padding to input image and guidance image
Step2. Get spatial kernel (Gaussian Kernel)
Step3. Vetorize image to make the convolution operation just becomes a matrix product(without using for loop)
Step4. Get range kernel (use guidance image as reference)
> Before calculating range kernel, pixel values should be normalized to [0,1] to construct range kernel
> Divide all channel with 255

Step5. Do convolution on each channel
Step6. Output filtered image

## Check if Joint Bilateral Filter algorithm is correct and its execution time for one image (eval.py)
> Use eval.py to evaluate your JBL.py
```python
python eval.py --image_path './testdata/ex.png' --gt_bf_path './testdata/ex_gt_bf.png' --gt_jbf_paht './testdata/ex_gt_jbf.png'
```
> The error of bilateral filter and joint bilateral filter should both be 0.
```python
[Error] Bilateral: 0
[Error] Joint bilateral: 0
```

## Compare the perceptual similarity using different gray conversion image as reference (main.py)
> In the setting file gives sigma_s, sigma_r and five kinds of gray conversion parameters.
> Need to use those five and also original cv2 gray conversions as guidance to run Joint Bilateral Filter and compute L1-norm as our cost function.
