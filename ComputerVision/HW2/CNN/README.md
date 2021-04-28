# Handwritten Digits Classification using CNN
In this assignment, you will need to perform image classification which predicts labels to each image with CNN models. 
input: monochrome image 
output: classification label 

### ConvNet (LeNet-5)
#### Network architecture and number of parameters
ConvNet(
  (conv1): Conv2d(1, 6, kernel_size=(5, 5), stride=(1, 1))
  (conv2): Conv2d(6, 16, kernel_size=(5, 5), stride=(1, 1))
  (fc1): Linear(in_features=256, out_features=120, bias=True)
  (fc2): Linear(in_features=120, out_features=84, bias=True)
  (fc3): Linear(in_features=84, out_features=10, bias=True)
)
----------------------------------------------------------------
        Layer (type)               Output Shape         Param #
================================================================
            Conv2d-1            [-1, 6, 24, 24]             156
            Conv2d-2             [-1, 16, 8, 8]           2,416
            Linear-3                  [-1, 120]          30,840
            Linear-4                   [-1, 84]          10,164
            Linear-5                   [-1, 10]             850
================================================================
Total params: 44,426
Trainable params: 44,426
Non-trainable params: 0
----------------------------------------------------------------
Input size (MB): 0.00
Forward/backward pass size (MB): 0.04
Params size (MB): 0.17
----------------------------------------------------------------

#### Learning curve
![conv_train](https://user-images.githubusercontent.com/40762011/116354174-9b263f80-a82a-11eb-802a-bab022cfd6e6.png)

![conv_val](https://user-images.githubusercontent.com/40762011/116354205-a9745b80-a82a-11eb-91e0-1481b576030e.png)

### MyNet
#### Network architecture and number of parameters
MyNet(
  (conv1): Conv2d(1, 20, kernel_size=(5, 5), stride=(1, 1))
  (conv2): Conv2d(20, 50, kernel_size=(5, 5), stride=(1, 1))
  (fc1): Linear(in_features=800, out_features=500, bias=True)
  (fc2): Linear(in_features=500, out_features=10, bias=True)
)
----------------------------------------------------------------
        Layer (type)               Output Shape         Param #
================================================================
            Conv2d-1           [-1, 20, 24, 24]             520
            Conv2d-2             [-1, 50, 8, 8]          25,050
            Linear-3                  [-1, 500]         400,500
            Linear-4                   [-1, 10]           5,010
================================================================
Total params: 431,080
Trainable params: 431,080
Non-trainable params: 0
----------------------------------------------------------------
Input size (MB): 0.00
Forward/backward pass size (MB): 0.12
Params size (MB): 1.64
Estimated Total Size (MB): 1.76
----------------------------------------------------------------

#### Learning curve
![mynet_train](https://user-images.githubusercontent.com/40762011/116354108-85b11580-a82a-11eb-8462-b3c891f40e2e.png)

![mynet_val](https://user-images.githubusercontent.com/40762011/116354151-92356e00-a82a-11eb-8990-f06c370ad075.png)

