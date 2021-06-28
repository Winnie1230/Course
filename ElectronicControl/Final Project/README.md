# Final Project
The topic of our final open project is to implement SLAM and generate trajectory using obstacle point cloud map scanned by IR sensor. 
Our project can divide into two parts: SLAM, and trajectory planning.

## SLAM
Set three IR sensor on the car and save the three IR values and car position. After traveling, one can get the point cloud of obstacles.

## Trajectory planning
Try three methods to do trajectory planning 
1. electric field method(`field.m`)
2. conventional potential field method(`cpf.m`)
- rectangular obstacle map
3. modified conventional potential field method
- rectangular obtacle map(`mod_cpf_recmap.m`)
- actual SLAM point cloud(`mod_cpf_map_kp.m`)
