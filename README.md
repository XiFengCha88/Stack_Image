# Stack_Image
This is my python script for staking .fits files.

## Description
- __init__.py: An empty file
- maincode.py: Start the file: "StackImg.py"
- StackImg.py: the main script
  - So far, the section in the code
    - creating_direc_BDF(): Create directories BIAS, DARK, FLAT (BDF).
    - moving_BDF(): Move Image types of BDF files to corresponded directories.
    - search_header(): Search the header of .fits file.
    - stacking_BDF(): Stack BDF files.
    - calibrate_img(): Calibrate target images and improve.

- Update time stamp: 2023.Dec.8th 12:31 
