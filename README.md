# Ultrasound Joint Morphology MQP
MATLAB and Arduino code used for the development of a system that can adapt ultrasound (US) for imaging joint morphology. The final design consisted of a US transducer mounted on a frame-rail fixture that would traverse across the suprapatellar region of the knee. The orientation of the transducer was collected with individial B-scans and fed into a reconstruction algorithm to better visualize the underlying anatomy.

## Ultrasound Image Acquisition
MATLAB Image Acquisition Toolbox used to acquire individual B-scans (frames) from ultrasound system.

## Inertial Measurement Unit Data Acquisition
Measurements for rotational displacement acquired through MPU 6050 inertial measurement unit (IMU).

## Image Reconstruction
Custom polar panoramic image reconstruction algorithm used to reconstruct multiple B-scans into single image.
